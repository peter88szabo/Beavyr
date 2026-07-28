use bevy::camera::visibility::RenderLayers;
use bevy::camera::RenderTarget;
use bevy::ecs::observer::On;
use bevy::prelude::*;
use bevy::render::render_resource::{
    Extent3d, TextureDescriptor, TextureDimension, TextureFormat, TextureUsages,
};
use bevy::render::view::screenshot::{Captured, Screenshot, ScreenshotCaptured};
use rfd::FileDialog;
use std::path::{Path, PathBuf};

use crate::scene::LAYER_MAIN;

#[derive(Component)]
pub(crate) struct ExportCamera;

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExportFormat {
    Png,
    Jpg,
}

impl ExportFormat {
    pub(crate) fn ext(self) -> &'static str {
        match self {
            ExportFormat::Png => "png",
            ExportFormat::Jpg => "jpg",
        }
    }
}

#[derive(Resource)]
pub(crate) struct ExportSettings {
    pub(crate) format: ExportFormat,
    pub(crate) canvas_width: u32,
    pub(crate) canvas_height: u32,
    pub(crate) dpi: u32,
    pub(crate) last_dir: PathBuf,
}

impl Default for ExportSettings {
    fn default() -> Self {
        Self {
            format: ExportFormat::Png,
            canvas_width: 1920,
            canvas_height: 1080,
            dpi: 300,
            last_dir: Path::new(env!("CARGO_MANIFEST_DIR")).join("exports"),
        }
    }
}

#[derive(Resource, Default)]
pub(crate) struct ExportCounter {
    pub(crate) next_id: u32,
}

#[derive(Resource, Default, Clone)]
pub(crate) struct ExportSelectArea {
    pub(crate) active: bool,
    pub(crate) dragging: bool,
    pub(crate) start: Option<Vec2>,
    pub(crate) end: Option<Vec2>,
    pub(crate) has_area: bool,
}

impl ExportSelectArea {
    pub(crate) fn bounds(&self) -> Option<(Vec2, Vec2)> {
        let (Some(a), Some(b)) = (self.start, self.end) else {
            return None;
        };
        let min = Vec2::new(a.x.min(b.x), a.y.min(b.y));
        let max = Vec2::new(a.x.max(b.x), a.y.max(b.y));
        if (max.x - min.x).abs() < 1e-3 || (max.y - min.y).abs() < 1e-3 {
            None
        } else {
            Some((min, max))
        }
    }
}

#[derive(Component)]
pub(crate) struct ExportCleanup {
    pub(crate) camera: Entity,
    pub(crate) image: Handle<Image>,
}

#[derive(Component)]
pub(crate) struct PendingCanvasCapture {
    pub(crate) frames_until_capture: u8,
    pub(crate) path: PathBuf,
    pub(crate) camera: Entity,
    pub(crate) image: Handle<Image>,
    pub(crate) crop_min: UVec2,
    pub(crate) crop_size: UVec2,
}

pub(crate) fn prompt_export_save_path(
    export_settings: &mut ExportSettings,
    suggested_stem: &str,
    format: ExportFormat,
) -> Option<PathBuf> {
    let mut dialog = FileDialog::new().set_file_name(&format!("{suggested_stem}.{}", format.ext()));
    if export_settings.last_dir.is_dir() {
        dialog = dialog.set_directory(&export_settings.last_dir);
    }
    dialog = match format {
        ExportFormat::Png => dialog.add_filter("PNG image", &["png", "PNG"]),
        ExportFormat::Jpg => dialog.add_filter("JPEG image", &["jpg", "jpeg", "JPG", "JPEG"]),
    };

    let mut path = dialog.save_file()?;
    let ext_ok = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase())
        .map(|e| match format {
            ExportFormat::Png => e == "png",
            ExportFormat::Jpg => e == "jpg" || e == "jpeg",
        })
        .unwrap_or(false);
    if !ext_ok {
        path.set_extension(format.ext());
    }
    if let Some(parent) = path.parent() {
        export_settings.last_dir = parent.to_path_buf();
    }
    Some(path)
}

pub(crate) fn update_export_select_area(
    mut area: ResMut<ExportSelectArea>,
    buttons: Res<ButtonInput<MouseButton>>,
    windows: Query<&Window>,
    egui_wants_input: Res<bevy_egui::input::EguiWantsInput>,
) {
    if !area.active {
        return;
    }
    let Ok(window) = windows.single() else {
        return;
    };
    let pointer_over_ui = egui_wants_input.is_using_pointer() || egui_wants_input.is_popup_open();

    if buttons.just_pressed(MouseButton::Left) && !pointer_over_ui {
        if let Some(pos) = window.cursor_position() {
            area.start = Some(pos);
            area.end = Some(pos);
            area.dragging = true;
            area.has_area = false;
        }
    }

    if buttons.pressed(MouseButton::Left) && area.dragging {
        if let Some(pos) = window.cursor_position() {
            area.end = Some(pos);
        }
    }

    if buttons.just_released(MouseButton::Left) && area.dragging {
        if let Some(pos) = window.cursor_position() {
            area.end = Some(pos);
        }
        area.dragging = false;
        area.active = false;
        area.has_area = area.bounds().is_some();
    }
}

pub(crate) fn start_export_capture(
    commands: &mut Commands,
    images: &mut Assets<Image>,
    export_settings: &mut ExportSettings,
    export_counter: &mut ExportCounter,
    area: &ExportSelectArea,
    window: &Window,
    main_proj: &Projection,
    main_cam_xform: &GlobalTransform,
) -> bool {
    let Some((min, max)) = area.bounds() else {
        return false;
    };

    let dpi_scale = export_settings.dpi as f32 / 300.0;
    let window_w = window.physical_width() as f32;
    let window_h = window.physical_height() as f32;
    if window_w < 1.0 || window_h < 1.0 {
        return false;
    }

    let max_w = (export_settings.canvas_width as f32 * dpi_scale).max(1.0);
    let max_h = (export_settings.canvas_height as f32 * dpi_scale).max(1.0);
    let mut render_scale = dpi_scale.max(0.01);
    let target_w = window_w * render_scale;
    let target_h = window_h * render_scale;
    let limit_scale = (max_w / target_w).min(max_h / target_h).min(1.0);
    render_scale *= limit_scale.max(0.01);

    let render_w = (window_w * render_scale).ceil().max(1.0) as u32;
    let render_h = (window_h * render_scale).ceil().max(1.0) as u32;

    let scale_factor = window.scale_factor() as f32;
    let crop_min = Vec2::new(
        min.x.max(0.0).min(window_w / scale_factor),
        min.y.max(0.0).min(window_h / scale_factor),
    ) * scale_factor;
    let crop_max = Vec2::new(
        max.x.max(0.0).min(window_w / scale_factor),
        max.y.max(0.0).min(window_h / scale_factor),
    ) * scale_factor;

    let crop_min = crop_min * render_scale;
    let crop_max = crop_max * render_scale;

    let crop_min_u = UVec2::new(
        crop_min.x.floor().max(0.0) as u32,
        crop_min.y.floor().max(0.0) as u32,
    );
    let crop_max_u = UVec2::new(
        crop_max.x.ceil().max(0.0) as u32,
        crop_max.y.ceil().max(0.0) as u32,
    );
    let crop_w = crop_max_u.x.saturating_sub(crop_min_u.x).max(1);
    let crop_h = crop_max_u.y.saturating_sub(crop_min_u.y).max(1);
    let crop_size = UVec2::new(crop_w, crop_h);

    let id = export_counter.next_id;
    let format = export_settings.format;
    let Some(path) =
        prompt_export_save_path(export_settings, &format!("beavyr-canvas-{id:04}"), format)
    else {
        return false;
    };
    if let Some(parent) = path.parent() {
        let _ = std::fs::create_dir_all(parent);
    }
    export_counter.next_id += 1;

    let extent = Extent3d {
        width: render_w,
        height: render_h,
        depth_or_array_layers: 1,
    };
    let mut image = Image {
        texture_descriptor: TextureDescriptor {
            label: None,
            size: extent,
            dimension: TextureDimension::D2,
            format: TextureFormat::Bgra8UnormSrgb,
            mip_level_count: 1,
            sample_count: 1,
            usage: TextureUsages::TEXTURE_BINDING
                | TextureUsages::COPY_SRC
                | TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        },
        ..default()
    };
    image.resize(extent);
    let image_handle = images.add(image);

    let cam = Camera::default();
    let render_layers = RenderLayers::layer(LAYER_MAIN);
    let transform = Transform::from_matrix(main_cam_xform.compute_transform().to_matrix());

    let mut cloned_projection = main_proj.clone();
    if let Projection::Perspective(p) = &mut cloned_projection {
        p.aspect_ratio = render_w as f32 / render_h as f32;
    } else if let Projection::Orthographic(o) = &mut cloned_projection {
        let half_height = o.area.height() * 0.5;
        let half_width = half_height * (render_w as f32 / render_h as f32);
        o.area.min = Vec2::new(-half_width, -half_height);
        o.area.max = Vec2::new(half_width, half_height);
    }

    let camera_entity = commands
        .spawn((
            Camera3d::default(),
            ExportCamera,
            cam,
            RenderTarget::Image(image_handle.clone().into()),
            cloned_projection,
            transform,
            GlobalTransform::from(transform),
            render_layers,
        ))
        .id();
    commands.spawn(PendingCanvasCapture {
        frames_until_capture: 1,
        path,
        camera: camera_entity,
        image: image_handle,
        crop_min: crop_min_u,
        crop_size,
    });
    true
}

pub(crate) fn process_pending_canvas_captures(
    mut commands: Commands,
    mut pending: Query<(Entity, &mut PendingCanvasCapture)>,
) {
    for (entity, mut pending) in &mut pending {
        if pending.frames_until_capture > 0 {
            pending.frames_until_capture -= 1;
            continue;
        }
        let path = pending.path.clone();
        let crop_min = pending.crop_min;
        let crop_size = pending.crop_size;
        let image = pending.image.clone();
        commands
            .entity(entity)
            .remove::<PendingCanvasCapture>()
            .insert((
                Screenshot::image(image.clone()),
                ExportCleanup {
                    camera: pending.camera,
                    image,
                },
            ))
            .observe(save_cropped_to_disk(path, crop_min, crop_size));
    }
}

pub(crate) fn cleanup_export_captures(
    mut commands: Commands,
    mut images: ResMut<Assets<Image>>,
    captured: Query<(Entity, &ExportCleanup), With<Captured>>,
) {
    for (entity, cleanup) in &captured {
        commands.entity(cleanup.camera).despawn();
        images.remove(cleanup.image.id());
        commands.entity(entity).despawn();
    }
}

fn save_cropped_to_disk(
    path: PathBuf,
    crop_min: UVec2,
    crop_size: UVec2,
) -> impl FnMut(On<ScreenshotCaptured>) {
    move |screenshot_captured| {
        let img = screenshot_captured.image.clone();
        let dyn_img = match img.try_into_dynamic() {
            Ok(img) => img,
            Err(err) => {
                eprintln!("Failed to read screenshot image: {err}");
                return;
            }
        };
        let img_w = dyn_img.width();
        let img_h = dyn_img.height();
        if img_w == 0 || img_h == 0 {
            eprintln!("Empty screenshot image.");
            return;
        }
        let x = crop_min.x.min(img_w.saturating_sub(1));
        let y = crop_min.y.min(img_h.saturating_sub(1));
        let w = crop_size.x.min(img_w.saturating_sub(x)).max(1);
        let h = crop_size.y.min(img_h.saturating_sub(y)).max(1);
        let cropped = dyn_img.crop_imm(x, y, w, h);
        if let Err(err) = cropped.save(&path) {
            eprintln!("Failed to save screenshot: {err}");
        }
    }
}
