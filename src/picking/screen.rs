use bevy::prelude::*;

pub fn find_nearest_visible_atom(
    cam: &Camera,
    cam_xform: &GlobalTransform,
    atoms: &[Vec3],
    mouse_px: Vec2,
    radius_px: f32,
    visible: impl Fn(usize) -> bool,
) -> Option<(usize, Vec2)> {
    let mut best: Option<(usize, f32, Vec2)> = None;

    for (i, &p_world) in atoms.iter().enumerate() {
        if !visible(i) {
            continue;
        }
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
