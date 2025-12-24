// src/camera.rs
use bevy::input::mouse::{MouseMotion, MouseWheel};
use bevy::prelude::*;
use bevy_egui::EguiContexts;

use crate::scene::MainCamera;

/// Orbit camera with correct screen-space panning (eye translates with target).
#[derive(Resource)]
pub struct OrbitCamera {
    pub radius: f32,
    pub theta: f32,
    pub phi: f32,
    pub target: Vec3, // world position we look at
    pub rotate_sensitivity: f32,
    pub zoom_sensitivity: f32,
    pub pan_sensitivity: f32,
    pub rotating: bool,
    pub panning: bool,
}

impl Default for OrbitCamera {
    fn default() -> Self {
        OrbitCamera {
            radius: 12.0,
            theta: 0.7,
            phi: 0.9,
            target: Vec3::ZERO,
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
    windows: Query<&Window>,
    mut settings: ResMut<crate::settings::MolSettings>,
    mut contexts: EguiContexts,
) {
    let Ok(window) = windows.single() else { return; };
    let pointer_over_ui = contexts
        .ctx_mut()
        .ok()
        .map(|ctx| ctx.is_pointer_over_area())
        .unwrap_or(false);

    // --- begin/end interaction modes (mutually exclusive) ---
    if window.cursor_position().is_some() {
        if buttons.just_pressed(MouseButton::Left) {
            cam.rotating = true;
            cam.panning = false;
        }
        if buttons.just_pressed(MouseButton::Right) {
            cam.panning = true;
            cam.rotating = false;
        }
    }
    if buttons.just_released(MouseButton::Left)  { cam.rotating = false; }
    if buttons.just_released(MouseButton::Right) { cam.panning  = false; }

    // collect mouse delta
    let mut mouse_delta = Vec2::ZERO;
    for ev in mouse_evr.read() {
        mouse_delta += ev.delta;
    }

    // spherical offset from target
    let offset = Vec3::new(
        cam.theta.cos() * cam.phi.sin(),
        cam.phi.cos(),
        cam.theta.sin() * cam.phi.sin(),
    ) * cam.radius;

    // eye = target + offset
    let mut eye = cam.target + offset;

    // view basis from current view direction
    let forward = (cam.target - eye).normalize();
    let world_up = Vec3::Y;
    let mut right = forward.cross(world_up);
    if right.length_squared() < 1e-6 {
        right = Vec3::X;
    } else {
        right = right.normalize();
    }
    let up = right.cross(forward).normalize();

    // rotate (LMB)
    if cam.rotating {
        cam.theta -= mouse_delta.x * cam.rotate_sensitivity;
        cam.phi   -= mouse_delta.y * cam.rotate_sensitivity;
        cam.phi = cam.phi.clamp(0.01, std::f32::consts::PI - 0.01);
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
    let final_offset = Vec3::new(
        cam.theta.cos() * cam.phi.sin(),
        cam.phi.cos(),
        cam.theta.sin() * cam.phi.sin(),
    ) * cam.radius;
    let final_eye = cam.target + final_offset;

    if let Ok(mut tf) = q_cam.single_mut() {
        *tf = Transform::from_translation(final_eye).looking_at(cam.target, Vec3::Y);
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
