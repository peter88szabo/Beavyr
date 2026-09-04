// src/camera.rs
use bevy::input::mouse::{MouseMotion, MouseWheel};
use bevy::prelude::*;
use bevy::window::{CursorGrabMode, CursorOptions};
use bevy_egui::input::EguiWantsInput;

use crate::scene::MainCamera;

/// Orbit camera with correct screen-space panning (eye translates with target).
#[derive(Resource)]
pub struct OrbitCamera {
    pub radius: f32,
    /// Free orbit orientation. Quaternions avoid Euler-angle limits and allow
    /// continuous turns around both axes.
    pub orientation: Quat,
    pub target: Vec3, // world position we look at
    pub rotate_sensitivity: f32,
    pub zoom_sensitivity: f32,
    pub pan_sensitivity: f32,
    pub rotating: bool,
    pub panning: bool,
}

impl Default for OrbitCamera {
    fn default() -> Self {
        let radius: f32 = 12.0;
        let theta: f32 = 0.7;
        let phi: f32 = 0.9;
        let offset = Vec3::new(theta.cos() * phi.sin(), phi.cos(), theta.sin() * phi.sin());
        let orientation = Quat::from_rotation_arc(Vec3::Z, offset.normalize());
        OrbitCamera {
            radius,
            target: Vec3::ZERO,
            orientation,
            rotate_sensitivity: 0.005,
            zoom_sensitivity: 1.0,
            pan_sensitivity: 0.004,
            rotating: false,
            panning: false,
        }
    }
}

pub fn orbit_camera_system(
    mut q_cam: Query<&mut Transform, (With<Camera3d>, With<MainCamera>)>,
    mut cam: ResMut<OrbitCamera>,
    mut mouse_evr: MessageReader<MouseMotion>,
    buttons: Res<ButtonInput<MouseButton>>,
    mut scroll_evr: MessageReader<MouseWheel>, // <-- fixed: no extra '>'
    keys: Res<ButtonInput<KeyCode>>,
    time: Res<Time>,
    mut windows: Query<(&mut Window, &mut CursorOptions)>,
    mut settings: ResMut<crate::settings::MolSettings>,
    egui_wants_input: Res<EguiWantsInput>,
    panel_regions: Res<crate::ui::UiPanelRegions>,
    export_select_area: Res<crate::export_image::ExportSelectArea>,
) {
    let Ok((window, mut cursor_options)) = windows.single_mut() else {
        return;
    };
    // A pointer over any egui panel must never start or continue a viewport
    // interaction.
    //
    // `EguiWantsInput` alone is not enough here.  The control panels are built
    // into a hand-made root `Ui`, so egui's own hover test never learns that
    // the right-hand strip is occupied and reports the pointer as being over
    // the 3D view.  It only turns true once a widget is actively being dragged,
    // which is why pressing a slider was handled but hovering and scrolling a
    // list was not -- the wheel reached the camera and zoomed the molecule.
    // `UiPanelRegions` carries the geometry egui cannot supply.
    let cursor = window.cursor_position();
    let pointer_over_panel = cursor.is_some_and(|p| panel_regions.covers(p));
    let pointer_over_ui = egui_wants_input.wants_any_pointer_input() || pointer_over_panel;
    if pointer_over_ui {
        cam.rotating = false;
        cam.panning = false;
        cursor_options.grab_mode = CursorGrabMode::None;
        cursor_options.visible = true;
    }
    if export_select_area.active || export_select_area.dragging {
        cam.rotating = false;
        cam.panning = false;
        cursor_options.grab_mode = CursorGrabMode::None;
        cursor_options.visible = true;
        return;
    }

    // --- begin/end interaction modes (mutually exclusive) ---
    if window.cursor_position().is_some() {
        if buttons.just_pressed(MouseButton::Left) && !pointer_over_ui {
            cam.rotating = true;
            cam.panning = false;
            cursor_options.grab_mode = CursorGrabMode::Locked;
            cursor_options.visible = false;
        }
        if buttons.just_pressed(MouseButton::Right) && !pointer_over_ui {
            cam.panning = true;
            cam.rotating = false;
            cursor_options.grab_mode = CursorGrabMode::Locked;
            cursor_options.visible = false;
        }
    }
    if buttons.just_released(MouseButton::Left) {
        cam.rotating = false;
    }
    if buttons.just_released(MouseButton::Right) {
        cam.panning = false;
    }
    if !buttons.pressed(MouseButton::Left) && !buttons.pressed(MouseButton::Right) {
        // Derive the idle state from the buttons themselves rather than relying
        // on a release message.  Egui's input absorbing can clear a release that
        // happens over a panel, which would otherwise leave the camera stuck in
        // rotate or pan mode.
        cam.rotating = false;
        cam.panning = false;
        cursor_options.grab_mode = CursorGrabMode::None;
        cursor_options.visible = true;
    }

    // collect mouse delta
    let mut mouse_delta = Vec2::ZERO;
    for ev in mouse_evr.read() {
        mouse_delta += ev.delta;
    }

    // Rotate (LMB) with no angular limits.
    if cam.rotating {
        let yaw = Quat::from_axis_angle(Vec3::Y, -mouse_delta.x * cam.rotate_sensitivity);
        let right_axis = (cam.orientation * Vec3::X).normalize();
        let pitch = Quat::from_axis_angle(right_axis, -mouse_delta.y * cam.rotate_sensitivity);
        cam.orientation = (yaw * pitch * cam.orientation).normalize();
    }

    // Use the updated basis for panning and for the final camera transform.
    let right = cam.orientation * Vec3::X;
    let up = cam.orientation * Vec3::Y;
    // pan (RMB): move target AND eye together in camera plane
    if cam.panning {
        let pan_scale = cam.pan_sensitivity * cam.radius;
        // drag right -> move +right; drag up (negative y) -> move +up
        let delta_world = right * (mouse_delta.x * pan_scale) + up * (-mouse_delta.y * pan_scale);
        cam.target += delta_world;
    }

    // keyboard panning uses the same basis
    let dt = time.delta_secs();
    let pan_key_speed = cam.pan_sensitivity * cam.radius * 240.0 * dt;
    let mut key_delta = Vec3::ZERO;
    if keys.pressed(KeyCode::ArrowLeft) {
        key_delta -= right * pan_key_speed;
    }
    if keys.pressed(KeyCode::ArrowRight) {
        key_delta += right * pan_key_speed;
    }
    if keys.pressed(KeyCode::ArrowUp) {
        key_delta += up * pan_key_speed;
    }
    if keys.pressed(KeyCode::ArrowDown) {
        key_delta -= up * pan_key_speed;
    }
    if key_delta != Vec3::ZERO {
        cam.target += key_delta;
    }

    // Zoom.
    //
    // Always drain the reader, even when the pointer is over a panel.  Skipping
    // the read leaves the messages queued, and Bevy keeps them alive for two
    // frames, so a scroll spent on an egui list would still reach the camera a
    // frame later.  Egui's input absorbing normally clears these first; draining
    // here keeps the behaviour correct if that is ever turned off.
    let mut scroll_delta = 0.0;
    for ev in scroll_evr.read() {
        scroll_delta += ev.y;
    }
    if !pointer_over_ui && scroll_delta != 0.0 {
        cam.radius = (cam.radius - scroll_delta * cam.zoom_sensitivity).clamp(2.0, 200.0);
    }

    let final_eye = cam.target + (cam.orientation * Vec3::Z) * cam.radius;

    if let Ok(mut tf) = q_cam.single_mut() {
        *tf = Transform::from_translation(final_eye).looking_at(cam.target, up);
    }

    // quick BG toggles (unchanged)
    if keys.just_pressed(KeyCode::KeyB) {
        settings.bg_color = Color::srgb(0.0, 0.0, 0.0);
        settings.lighting_dirty = true;
    }
    if keys.just_pressed(KeyCode::KeyW) {
        settings.bg_color = Color::srgb(1.0, 1.0, 1.0);
        settings.lighting_dirty = true;
    }
}
