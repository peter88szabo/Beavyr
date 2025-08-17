use bevy::prelude::*;
use bevy::render::camera::{Projection, RenderTarget};
use bevy::render::render_resource::*;
use bevy::render::view::RenderLayers;
use bevy_image_export::{ImageExport, ImageExportPlugin, ImageExportSettings, ImageExportSource};
use std::path::PathBuf;

use crate::scene::{MainCamera, LAYER_MAIN};

/// Public plugin you’ll add in `main.rs`.
pub struct ExportPlugin;

impl Plugin for ExportPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(ImageExportPlugin::default())
            .init_resource::<ExportUiState>()
            .init_resource::<ExportRequestQueue>()
            .add_systems(
                Update,
                (
                    consume_export_requests,
                    tick_and_cleanup_temp_cameras,
                    tick_and_trigger_delayed_exports, // spawn exporter after a few frames
                    tick_and_cleanup_exporters,        // exporter lives 1 frame -> exactly one image
                ),
            );
    }
}

/// UI state for export panel.
#[derive(Resource)]
pub struct ExportUiState {
    pub open: bool,
    pub dpi: u32,        // 300 / 600 / 1200
    pub width_mm: f32,   // physical size
    pub height_mm: f32,  // physical size
    pub last_dir: Option<PathBuf>,
}
impl Default for ExportUiState {
    fn default() -> Self {
        Self {
            open: false,
            dpi: 300,
            width_mm: 100.0,
            height_mm: 100.0,
            last_dir: None,
        }
    }
}

/// Parameters of a pending export, pushed by the UI (pixels + folder).
#[derive(Clone)]
pub struct ExportParams {
    pub output_dir: PathBuf,
    pub extension: String, // "png" or "exr"
    pub width: u32,        // px
    pub height: u32,       // px
}

/// Queue for pending exports.
#[derive(Resource, Default)]
pub struct ExportRequestQueue {
    pub queue: Vec<ExportParams>,
}

/// Temp camera kept alive briefly to ensure a rendered frame exists.
#[derive(Component)]
struct TempExportCamera {
    frames_waited: u8,
    wait_frames: u8,
}

/// Spawns the exporter after a small delay, so exactly one non-empty frame is written.
#[derive(Component)]
struct DelayedExportArm {
    frames_left: u8,
    source: Handle<ImageExportSource>,
    output_dir: PathBuf,
    extension: String,
}

/// Despawn the exporter after exactly one frame (so only one image is produced).
#[derive(Component)]
struct OneShotExport {
    frames_waited: u8,
    wait_frames: u8,
}

fn consume_export_requests(
    mut commands: Commands,
    mut images: ResMut<Assets<Image>>,
    mut export_sources: ResMut<Assets<ImageExportSource>>,
    mut queue: ResMut<ExportRequestQueue>,
    main_cam_q: Query<(&Camera, &Projection, &GlobalTransform), With<MainCamera>>,
) {
    if queue.queue.is_empty() {
        return;
    }

    let (_main_cam, main_proj, main_cam_xform) = match main_cam_q.get_single() {
        Ok(t) => t,
        Err(_) => return,
    };

    for req in queue.queue.drain(..) {
        let w = req.width.max(16);
        let h = req.height.max(16);
        let aspect = w as f32 / h as f32;

        // Output texture
        let size = Extent3d { width: w, height: h, depth_or_array_layers: 1 };
        let mut export_image = Image {
            texture_descriptor: TextureDescriptor {
                label: Some("export_texture"),
                size,
                dimension: TextureDimension::D2,
                format: TextureFormat::Rgba8UnormSrgb,
                mip_level_count: 1,
                sample_count: 1,
                usage: TextureUsages::COPY_DST
                    | TextureUsages::COPY_SRC
                    | TextureUsages::RENDER_ATTACHMENT,
                view_formats: &[],
            },
            ..default()
        };
        export_image.resize(size);
        let output_image = images.add(export_image);

        // Temp camera rendering to that texture (main layer only, no axes)
        let cam3d = Camera3d::default();

        let mut cam = Camera::default();
        cam.target = RenderTarget::Image(output_image.clone().into());
        let render_layers = RenderLayers::layer(LAYER_MAIN);

        let transform =
            Transform::from_matrix(main_cam_xform.compute_transform().compute_matrix());

        let mut cloned_projection = main_proj.clone();
        if let Projection::Perspective(p) = &mut cloned_projection {
            p.aspect_ratio = aspect;
        } else if let Projection::Orthographic(o) = &mut cloned_projection {
            let half_height = o.area.height() * 0.5;
            let half_width = half_height * aspect;
            o.area.min = Vec2::new(-half_width, -half_height);
            o.area.max = Vec2::new( half_width,  half_height);
        }

        // Keep camera alive long enough for warm-up + one rendered frame
        commands.spawn((
            TempExportCamera { frames_waited: 0, wait_frames: 4 },
            cam3d, cam, cloned_projection,
            transform, GlobalTransform::from(transform),
            render_layers,
        ));

        // Arm the export to trigger AFTER a short delay (no immediate exporting).
        let source = export_sources.add(output_image);
        commands.spawn(DelayedExportArm {
            frames_left: 2, // wait a couple frames, then spawn exporter for exactly one frame
            source,
            output_dir: req.output_dir,
            extension: req.extension,
        });
    }
}

fn tick_and_cleanup_temp_cameras(
    mut commands: Commands,
    mut q: Query<(Entity, &mut TempExportCamera)>,
) {
    for (e, mut m) in &mut q {
        m.frames_waited = m.frames_waited.saturating_add(1);
        if m.frames_waited >= m.wait_frames {
            commands.entity(e).despawn_recursive();
        }
    }
}

fn tick_and_trigger_delayed_exports(
    mut commands: Commands,
    mut q: Query<(Entity, &mut DelayedExportArm)>,
) {
    for (e, mut arm) in &mut q {
        if arm.frames_left > 0 {
            arm.frames_left -= 1;
            continue;
        }

        // Spawn exporter now, for exactly one frame.
        commands.spawn((
            ImageExport(arm.source.clone()),
            ImageExportSettings {
                output_dir: arm
                    .output_dir
                    .display()
                    .to_string()
                    .into(),
                extension: if arm.extension.eq_ignore_ascii_case("exr") { "exr".into() } else { "png".into() },
            },
            OneShotExport { frames_waited: 0, wait_frames: 1 },
        ));

        // Remove the arming entity so it doesn't trigger again.
        commands.entity(e).despawn_recursive();
    }
}

fn tick_and_cleanup_exporters(
    mut commands: Commands,
    mut q: Query<(Entity, &mut OneShotExport), With<ImageExport>>,
) {
    for (e, mut m) in &mut q {
        m.frames_waited = m.frames_waited.saturating_add(1);
        if m.frames_waited >= m.wait_frames {
            commands.entity(e).despawn_recursive(); // stop after exactly one frame
        }
    }
}

