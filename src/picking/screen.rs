use bevy::prelude::*;
use bevy_egui::egui;

/// Screen-space nearest-atom helper (PHYSICAL px API).
///
/// Returns (atom_index, mouse_px_pos) if within `radius_px`.
pub fn find_nearest_atom_screen_space(
    cam: &Camera,
    cam_xform: &GlobalTransform,
    atoms: &[Vec3],
    mouse_px: Vec2,
    radius_px: f32,
) -> Option<(usize, Vec2)> {
    let mut best: Option<(usize, f32, Vec2)> = None;

    for (i, &p_world) in atoms.iter().enumerate() {
        if let Ok(screen_px) = cam.world_to_viewport(cam_xform, p_world) {
            let d = screen_px.distance(mouse_px);
            if d <= radius_px {
                match best {
                    Some((_, best_d, _)) if d >= best_d => {}
                    _ => best = Some((i, d, screen_px)),
                }
            }
        }
    }

    best.map(|(i, _d, s)| (i, s))
}

/// Convenience wrapper using egui logical points (egui::Pos2) and returning the same.
/// Useful from panels/overlays that already operate in egui coordinates.
pub fn find_nearest_atom_screen_space_egui(
    cam: &Camera,
    cam_xform: &GlobalTransform,
    atoms: &[Vec3],
    mouse_px: egui::Pos2,
    radius_px: f32,
) -> Option<(usize, egui::Pos2)> {
    let mouse_v = Vec2::new(mouse_px.x, mouse_px.y);
    find_nearest_atom_screen_space(cam, cam_xform, atoms, mouse_v, radius_px)
        .map(|(i, v)| (i, egui::pos2(v.x, v.y)))
}
