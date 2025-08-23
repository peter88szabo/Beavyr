use bevy::prelude::*;
use bevy_egui::EguiContexts;

use crate::molecule::{covalent_radius_angstrom, Molecule};
use crate::settings::MolSettings;

#[derive(Default, Reflect, GizmoConfigGroup)]
pub struct MeasurementGizmos;

#[derive(Debug, Clone)]
pub struct MeasurePair {
    pub id: u32,
    pub a: usize,
    pub b: usize,
    pub distance: f32, // Å
    pub visible: bool,
    pub label_on: bool,
}

#[derive(Resource, Debug, Clone)]
pub struct Measurements {
    pub pairs: Vec<MeasurePair>,
    pub is_active: bool,
    pub pending: Option<usize>, // first picked atom (if any)

    // Global style (uniform for all pairs)
    pub color: Color,
    pub line_width: f32,
    pub dash_gap_scale: f32,
    pub dash_line_scale: f32,
    pub show_labels: bool,

    // Picking settings
    /// Extra slack added to the atom sphere radius when ray-testing (Å). (Kept for future use.)
    pub pick_radius: f32,
    /// Screen-space picking radius in pixels (used when measurement mode is active)
    pub pick_radius_px: f32,

    next_id: u32,
}

impl Default for Measurements {
    fn default() -> Self {
        Self {
            pairs: Vec::new(),
            is_active: false,
            pending: None,
            color: Color::srgba(0.0, 0.9, 0.9, 1.0),
            line_width: 20.0,
            dash_gap_scale: 1.5,
            dash_line_scale: 2.0,
            show_labels: true,
            pick_radius: 0.2,
            pick_radius_px: 18.0,
            next_id: 1,
        }
    }
}

impl Measurements {
    pub fn new_pair(&mut self, a: usize, b: usize, dist: f32) {
        let (ai, bi) = if a <= b { (a, b) } else { (b, a) };
        if self.pairs.iter().any(|p| p.a == ai && p.b == bi) {
            return;
        }
        let id = self.next_id;
        self.next_id += 1;
        self.pairs.push(MeasurePair {
            id,
            a: ai,
            b: bi,
            distance: dist,
            visible: true,
            label_on: true,
        });
    }

    pub fn remove_pair(&mut self, id: u32) {
        self.pairs.retain(|p| p.id != id);
    }

    pub fn clear(&mut self) {
        self.pairs.clear();
        self.pending = None;
    }
}

/// Configure dashed-line gizmo style (global, uniform for all pairs)
pub fn configure_measurement_gizmos(
    mut cfg_store: ResMut<GizmoConfigStore>,
    measurements: Res<Measurements>,
) {
    let (cfg, _) = cfg_store.config_mut::<MeasurementGizmos>();
    cfg.enabled = true;
    cfg.line.width = measurements.line_width.max(1.0).min(50.0);
    cfg.line.perspective = true;
    cfg.line.style = GizmoLineStyle::Dashed {
        gap_scale: measurements.dash_gap_scale.max(0.05),
        line_scale: measurements.dash_line_scale.max(0.05),
    };
}

/// Screen-space picking:
/// - Active only when measurement tool is toggled on
/// - Ignores clicks when mouse is over egui
/// - Left-click near an atom selects it; picking two atoms creates a new pair
pub fn handle_measurement_picking(
    mut contexts: EguiContexts,
    buttons: Res<ButtonInput<MouseButton>>,
    windows: Query<&Window>,
    // Match ui.rs: main scene camera
    q_cam: Query<(&Camera, &GlobalTransform), (With<Camera3d>, With<crate::scene::MainCamera>)>,
    _settings: Res<MolSettings>,
    mol: Option<Res<Molecule>>,
    mut measurements: ResMut<Measurements>,
) {
    let Some(mol) = mol else { return; };
    if !measurements.is_active { return; }
    if !buttons.just_pressed(MouseButton::Left) { return; }

    // Skip if pointer is over any egui area
    let Ok(ctx) = contexts.ctx_mut() else { return; };
    if ctx.is_pointer_over_area() {
        return;
    }

    // Camera + window + cursor (physical px)
    let Ok(window) = windows.single() else { return; };
    let Ok((cam, cam_xform)) = q_cam.single() else { return; };
    let Some(mouse_px) = window.cursor_position() else { return; };

    // Find nearest atom within pixel radius
    let hit = find_nearest_atom_screen_space(
        cam, cam_xform, &mol.pos, mouse_px, measurements.pick_radius_px,
    );

    // If we didn't hit anything, keep pending highlight untouched
    let Some((hit_idx, _)) = hit else { return; };

    // Only clear pending after a successful pair creation or explicit restart
    if let Some(first) = measurements.pending {
        if first == hit_idx {
            // Same atom clicked again: keep waiting for a different second atom
            measurements.pending = Some(first);
        } else {
            // Valid pair -> create it and clear pending
            let dist = mol.pos[first].distance(mol.pos[hit_idx]);
            measurements.new_pair(first, hit_idx, dist);
            measurements.pending = None;
        }
    } else {
        // Start a new pair
        measurements.pending = Some(hit_idx);
    }
}

/// Update distances every frame from Molecule.pos (Å)
pub fn update_measurement_distances(
    mol: Option<Res<Molecule>>,
    mut measurements: ResMut<Measurements>,
) {
    let Some(mol) = mol else { return; };
    let n = mol.pos.len();
    for p in &mut measurements.pairs {
        if p.a < n && p.b < n {
            p.distance = mol.pos[p.a].distance(mol.pos[p.b]);
        }
    }
}

/// Draw:
/// - the in-progress selection highlight (if `pending` is Some)
/// - dashed lines for each visible measured pair (uniform style)
pub fn draw_measurement_lines(
    mut gizmos: Gizmos<MeasurementGizmos>,
    measurements: Res<Measurements>,
    mol: Option<Res<Molecule>>,
    settings: Res<MolSettings>,
) {
    let Some(mol) = mol else { return; };
    let n = mol.pos.len();
    if n == 0 { return; }

    // 1) Pending atom highlight: visible magenta sphere around the first picked atom
    if let Some(i) = measurements.pending {
        if i < n {
            let center = mol.pos[i];
            // Slightly larger than the atom so it reads clearly
            let base_r = covalent_radius_angstrom(&mol.atoms[i]) * settings.atom_scale;
            let r = (base_r * 1.25).max(0.15);
            let highlight = Color::srgba(1.0, 0.0, 0.8, 1.0); // magenta
            gizmos.sphere(center, r, highlight);
        }
    }

    // 2) Draw all measured pairs as dashed lines
    for pair in measurements.pairs.iter().filter(|p| p.visible) {
        if pair.a >= n || pair.b >= n { continue; }

        let p0 = mol.pos[pair.a];
        let p1 = mol.pos[pair.b];

        let d = p1 - p0;
        let len = d.length();
        if len <= 1e-4 { continue; }
        let dn = d / len;

        // inset start/end just inside atom spheres
        let ri = covalent_radius_angstrom(&mol.atoms[pair.a]) * settings.atom_scale;
        let rj = covalent_radius_angstrom(&mol.atoms[pair.b]) * settings.atom_scale;

        let surface_inset_frac = 0.12;
        let inset_i = (ri * surface_inset_frac).min(ri * 0.9);
        let inset_j = (rj * surface_inset_frac).min(rj * 0.9);

        let a = p0 + dn * (ri - inset_i);
        let b = p1 - dn * (rj - inset_j);
        if (b - a).length() <= 1e-4 { continue; }

        gizmos.line(a, b, measurements.color);
    }
}

// ---------- helpers ----------

/// Return (index, mouse_px_pos) of nearest atom within `radius_px` in PHYSICAL px.
/// Mirrors your ui.rs helper.
fn find_nearest_atom_screen_space(
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
                if let Some((_, best_d, _)) = best {
                    if d < best_d {
                        best = Some((i, d, screen_px));
                    }
                } else {
                    best = Some((i, d, screen_px));
                }
            }
        }
    }

    best.map(|(i, _d, s)| (i, s))
}

