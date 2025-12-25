// src/camera.rs
use bevy::input::mouse::{MouseMotion, MouseWheel};
use bevy::prelude::*;
use bevy::window::CursorGrabMode;
use bevy_egui::EguiContexts;

use crate::scene::MainCamera;

/// Orbit camera with correct screen-space panning (eye translates with target).
#[derive(Resource)]
pub struct OrbitCamera {
    pub radius: f32,
    #[allow(dead_code)]
    pub theta: f32,
    #[allow(dead_code)]
    pub phi: f32,
    pub target: Vec3, // world position we look at
    pub orientation: Quat,
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
        let offset = Vec3::new(
            theta.cos() * phi.sin(),
            phi.cos(),
            theta.sin() * phi.sin(),
        );
        let orientation = Quat::from_rotation_arc(Vec3::Z, offset.normalize());
        OrbitCamera {
            radius,
            theta,
            phi,
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
    mut mouse_evr: EventReader<MouseMotion>,
    buttons: Res<ButtonInput<MouseButton>>,
    mut scroll_evr: EventReader<MouseWheel>,  // <-- fixed: no extra '>'
    keys: Res<ButtonInput<KeyCode>>,
    time: Res<Time>,
    mut windows: Query<&mut Window>,
    mut settings: ResMut<crate::settings::MolSettings>,
    mut contexts: EguiContexts,
) {
    let Ok(mut window) = windows.single_mut() else { return; };
    let pointer_over_ui = contexts
        .ctx_mut()
        .ok()
        .map(|ctx| ctx.is_pointer_over_area())
        .unwrap_or(false);

    // --- begin/end interaction modes (mutually exclusive) ---
    if window.cursor_position().is_some() {
        if buttons.just_pressed(MouseButton::Left) && !pointer_over_ui {
            cam.rotating = true;
            cam.panning = false;
            window.cursor_options.grab_mode = CursorGrabMode::Locked;
            window.cursor_options.visible = false;
        }
        if buttons.just_pressed(MouseButton::Right) && !pointer_over_ui {
            cam.panning = true;
            cam.rotating = false;
            window.cursor_options.grab_mode = CursorGrabMode::Locked;
            window.cursor_options.visible = false;
        }
    }
    if buttons.just_released(MouseButton::Left)  { cam.rotating = false; }
    if buttons.just_released(MouseButton::Right) { cam.panning  = false; }
    if !buttons.pressed(MouseButton::Left) && !buttons.pressed(MouseButton::Right) {
        window.cursor_options.grab_mode = CursorGrabMode::None;
        window.cursor_options.visible = true;
    }

    // collect mouse delta
    let mut mouse_delta = Vec2::ZERO;
    for ev in mouse_evr.read() {
        mouse_delta += ev.delta;
    }

    // offset from target (using quaternion orientation)
    let offset = (cam.orientation * Vec3::Z) * cam.radius;

    // eye = target + offset
    let mut eye = cam.target + offset;

    // view basis from current orientation
    let right = cam.orientation * Vec3::X;
    let up = cam.orientation * Vec3::Y;

    // rotate (LMB)
    if cam.rotating {
        let yaw = -mouse_delta.x * cam.rotate_sensitivity;
        let pitch = -mouse_delta.y * cam.rotate_sensitivity;
        let yaw_q = Quat::from_axis_angle(Vec3::Y, yaw);
        let right_axis = (cam.orientation * Vec3::X).normalize();
        let pitch_q = Quat::from_axis_angle(right_axis, pitch);
        cam.orientation = (yaw_q * pitch_q) * cam.orientation;
        cam.orientation = cam.orientation.normalize();
    }
    // pan (RMB): move target AND eye together in camera plane
    else if cam.panning {
        let pan_scale = cam.pan_sensitivity * cam.radius;
        // drag right -> move +right; drag up (negative y) -> move +up
        let delta_world = right * (mouse_delta.x * pan_scale)
                        + up    * (-mouse_delta.y * pan_scale);
        cam.target += delta_world;
        eye        += delta_world;
    }

    // keyboard panning uses the same basis
    let dt = time.delta_secs();
    let pan_key_speed = cam.pan_sensitivity * cam.radius * 240.0 * dt;
    let mut key_delta = Vec3::ZERO;
    if keys.pressed(KeyCode::ArrowLeft)  { key_delta -= right * pan_key_speed; }
    if keys.pressed(KeyCode::ArrowRight) { key_delta += right * pan_key_speed; }
    if keys.pressed(KeyCode::ArrowUp)    { key_delta += up    * pan_key_speed; }
    if keys.pressed(KeyCode::ArrowDown)  { key_delta -= up    * pan_key_speed; }
    if key_delta != Vec3::ZERO {
        cam.target += key_delta;
        eye        += key_delta;
    }

    // zoom
    if !pointer_over_ui {
        for ev in scroll_evr.read() {
            cam.radius = (cam.radius - ev.y * cam.zoom_sensitivity).clamp(2.0, 200.0);
        }
    }

    // recompute final eye from (theta,phi,radius) & updated target
    let final_eye = cam.target + (cam.orientation * Vec3::Z) * cam.radius;

    if let Ok(mut tf) = q_cam.single_mut() {
        *tf = Transform::from_translation(final_eye).looking_at(cam.target, up);
    }

    // quick BG toggles (unchanged)
    if keys.just_pressed(KeyCode::KeyB) {
        settings.bg_color = Color::srgb(0.0, 0.0, 0.0);
        settings.dirty = true;
    }
    if keys.just_pressed(KeyCode::KeyW) {
        settings.bg_color = Color::srgb(1.0, 1.0, 1.0);
        settings.dirty = true;
    }
}
