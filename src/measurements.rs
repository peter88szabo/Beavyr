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

#[derive(Debug, Clone)]
pub struct AngleMeasure {
    pub id: u32,
    pub a: usize,
    pub b: usize,
    pub c: usize,
    pub degrees: f32,
}

#[derive(Debug, Clone)]
pub struct DihedralMeasure {
    pub id: u32,
    pub a: usize,
    pub b: usize,
    pub c: usize,
    pub d: usize,
    pub degrees: f32,
}

#[derive(Resource, Debug, Clone)]
pub struct Measurements {
    // ---- Distances ----
    pub pairs: Vec<MeasurePair>,
    pub is_active: bool,           // distance picking active
    pub pending: Option<usize>,    // first picked atom (distance)

    // ---- Angles ----
    pub angles: Vec<AngleMeasure>,
    pub angle_active: bool,
    pub pending_angle: Vec<usize>, // will hold up to 2 indices

    // ---- Dihedrals ----
    pub dihedrals: Vec<DihedralMeasure>,
    pub dihedral_active: bool,
    pub pending_dihedral: Vec<usize>, // will hold up to 3 indices

    // Global style for distance lines
    pub color: Color,
    pub line_width: f32,
    pub dash_gap_scale: f32,
    pub dash_line_scale: f32,
    pub show_labels: bool,

    // Picking settings
    pub pick_radius: f32,
    pub pick_radius_px: f32,

    next_id: u32,
}

impl Default for Measurements {
    fn default() -> Self {
        Self {
            pairs: Vec::new(),
            is_active: false,
            pending: None,

            angles: Vec::new(),
            angle_active: false,
            pending_angle: Vec::new(),

            dihedrals: Vec::new(),
            dihedral_active: false,
            pending_dihedral: Vec::new(),

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
    // ---- Distance API ----
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

    // ---- Angle API ----
    pub fn add_angle(&mut self, a: usize, b: usize, c: usize, deg: f32) {
        if self.angles.iter().any(|x| x.a == a && x.b == b && x.c == c) {
            return;
        }
        let id = self.next_id;
        self.next_id += 1;
        self.angles.push(AngleMeasure { id, a, b, c, degrees: deg });
    }
    pub fn clear_angles(&mut self) {
        self.angles.clear();
        self.pending_angle.clear();
    }

    // ---- Dihedral API ----
    pub fn add_dihedral(&mut self, a: usize, b: usize, c: usize, d: usize, deg: f32) {
        if self
            .dihedrals
            .iter()
            .any(|x| x.a == a && x.b == b && x.c == c && x.d == d)
        {
            return;
        }
        let id = self.next_id;
        self.next_id += 1;
        self.dihedrals.push(DihedralMeasure { id, a, b, c, d, degrees: deg });
    }
    pub fn clear_dihedrals(&mut self) {
        self.dihedrals.clear();
        self.pending_dihedral.clear();
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

/// Screen-space picking for all three modes.
/// Priority if multiple toggles are on: dihedral > angle > distance.
pub fn handle_measurement_picking(
    mut contexts: EguiContexts,
    buttons: Res<ButtonInput<MouseButton>>,
    windows: Query<&Window>,
    q_cam: Query<(&Camera, &GlobalTransform), (With<Camera3d>, With<crate::scene::MainCamera>)>,
    _settings: Res<MolSettings>,
    mol: Option<Res<Molecule>>,
    mut measurements: ResMut<Measurements>,
) {
    let Some(mol) = mol else { return; };

    let any_active = measurements.dihedral_active || measurements.angle_active || measurements.is_active;
    if !any_active { return; }
    if !buttons.just_pressed(MouseButton::Left) { return; }

    let Ok(ctx) = contexts.ctx_mut() else { return; };
    if ctx.is_pointer_over_area() {
        return;
    }

    let Ok(window) = windows.single() else { return; };
    let Ok((cam, cam_xform)) = q_cam.single() else { return; };
    let Some(mouse_px) = window.cursor_position() else { return; };

    let hit = find_nearest_atom_screen_space(
        cam, cam_xform, &mol.pos, mouse_px, measurements.pick_radius_px,
    );
    let Some((hit_idx, _)) = hit else { return; };

    // --- Priority: dihedral > angle > distance
    if measurements.dihedral_active {
        measurements.pending_dihedral.push(hit_idx);
        if measurements.pending_dihedral.len() == 4 {
            let a = measurements.pending_dihedral[0];
            let b = measurements.pending_dihedral[1];
            let c = measurements.pending_dihedral[2];
            let d = measurements.pending_dihedral[3];
            let deg = dihedral_degrees(mol.pos[a], mol.pos[b], mol.pos[c], mol.pos[d]);
            measurements.add_dihedral(a, b, c, d, deg);
            measurements.pending_dihedral.clear();
            measurements.dihedral_active = false; // auto-off when complete
        }
        return;
    }

    if measurements.angle_active {
        measurements.pending_angle.push(hit_idx);
        if measurements.pending_angle.len() == 3 {
            let a = measurements.pending_angle[0];
            let b = measurements.pending_angle[1];
            let c = measurements.pending_angle[2];
            let deg = angle_degrees(mol.pos[a], mol.pos[b], mol.pos[c]);
            measurements.add_angle(a, b, c, deg);
            measurements.pending_angle.clear();
            measurements.angle_active = false; // auto-off when complete
        }
        return;
    }

    // Distance
    if measurements.is_active {
        if let Some(first) = measurements.pending {
            if first == hit_idx {
                measurements.pending = Some(first);
            } else {
                let dist = mol.pos[first].distance(mol.pos[hit_idx]);
                measurements.new_pair(first, hit_idx, dist);
                measurements.pending = None;
                measurements.is_active = false; // auto-off when complete
            }
        } else {
            measurements.pending = Some(hit_idx);
        }
    }
}

/// Update measured values every frame
pub fn update_measurement_distances(
    mol: Option<Res<Molecule>>,
    mut measurements: ResMut<Measurements>,
) {
    let Some(mol) = mol else { return; };
    let n = mol.pos.len();

    // distances
    for p in &mut measurements.pairs {
        if p.a < n && p.b < n {
            p.distance = mol.pos[p.a].distance(mol.pos[p.b]);
        }
    }

    // angles
    for a in &mut measurements.angles {
        if a.a < n && a.b < n && a.c < n {
            a.degrees = angle_degrees(mol.pos[a.a], mol.pos[a.b], mol.pos[a.c]);
        }
    }

    // dihedrals
    for d in &mut measurements.dihedrals {
        if d.a < n && d.b < n && d.c < n && d.d < n {
            d.degrees = dihedral_degrees(
                mol.pos[d.a], mol.pos[d.b], mol.pos[d.c], mol.pos[d.d]
            );
        }
    }
}

/// Draw gizmos (distance lines + pending highlights for distance/angle/dihedral)
pub fn draw_measurement_lines(
    mut gizmos: Gizmos<MeasurementGizmos>,
    measurements: Res<Measurements>,
    mol: Option<Res<Molecule>>,
    settings: Res<MolSettings>,
) {
    let Some(mol) = mol else { return; };
    let n = mol.pos.len();
    if n == 0 { return; }

    // ---------- Pending highlights ----------
    // Common radius calculator
    let mut draw_highlight = |idx: usize, color: Color| {
        if idx < n {
            let center = mol.pos[idx];
            let base_r = covalent_radius_angstrom(&mol.atoms[idx]) * settings.atom_scale;
            let r = (base_r * 1.25).max(0.15);
            gizmos.sphere(center, r, color);
        }
    };

    // Colors
    let magenta = Color::srgba(1.0, 0.0, 0.8, 1.0);
    let cyan_blue = Color::srgba(0.0, 0.8, 1.0, 1.0);
    let neon_green = Color::srgba(0.2, 1.0, 0.2, 1.0);

    // Distance pending (first atom)
    if let Some(i) = measurements.pending {
        draw_highlight(i, magenta);
    }

    // Angle pending: first = magenta, second = cyan-blue
    if !measurements.pending_angle.is_empty() {
        if let Some(&i0) = measurements.pending_angle.get(0) {
            draw_highlight(i0, magenta);
        }
        if let Some(&i1) = measurements.pending_angle.get(1) {
            draw_highlight(i1, cyan_blue);
        }
    }

    // Dihedral pending: first = magenta, second = cyan-blue, third = neon-green
    if !measurements.pending_dihedral.is_empty() {
        if let Some(&i0) = measurements.pending_dihedral.get(0) {
            draw_highlight(i0, magenta);
        }
        if let Some(&i1) = measurements.pending_dihedral.get(1) {
            draw_highlight(i1, cyan_blue);
        }
        if let Some(&i2) = measurements.pending_dihedral.get(2) {
            draw_highlight(i2, neon_green);
        }
    }

    // ---------- Distance lines (as before) ----------
    for pair in measurements.pairs.iter().filter(|p| p.visible) {
        if pair.a >= n || pair.b >= n { continue; }

        let p0 = mol.pos[pair.a];
        let p1 = mol.pos[pair.b];

        let d = p1 - p0;
        let len = d.length();
        if len <= 1e-4 { continue; }
        let dn = d / len;

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

fn angle_degrees(p1: Vec3, p2: Vec3, p3: Vec3) -> f32 {
    let v1 = (p1 - p2);
    let v2 = (p3 - p2);
    let n1 = v1.length();
    let n2 = v2.length();
    if n1 <= 1e-6 || n2 <= 1e-6 {
        return 0.0;
    }
    let c = (v1.dot(v2) / (n1 * n2)).clamp(-1.0, 1.0);
    c.acos().to_degrees()
}

fn dihedral_degrees(p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3) -> f32 {
    // Standard atan2 formulation for torsion angle
    let b1 = p2 - p1;
    let b2 = p3 - p2;
    let b3 = p4 - p3;

    let n1 = b1.cross(b2);
    let n2 = b2.cross(b3);

    let n1_len = n1.length();
    let n2_len = n2.length();
    if n1_len <= 1e-6 || n2_len <= 1e-6 {
        return 0.0;
    }

    let m1 = n1.cross(b2.normalize());
    let x = n1.dot(n2) / (n1_len * n2_len);
    let y = m1.dot(n2) / (n1_len * n2_len);

    y.atan2(x).to_degrees()
}

