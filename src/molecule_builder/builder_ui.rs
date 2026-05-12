use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};

use crate::events::{AtomPicked, ToolKind};
use crate::molecule::{Molecule, covalent_radius_angstrom};
use crate::settings::MolSettings;
use crate::color_schemes::color_for;
use super::rotator::{bend_side, rotate_side, split_sides, translate_side, RotateSide};
use super::zmat2xyz::{self, ZAtom};
use crate::molecule::parse_xyz_angstrom;
use super::fragments;

const ANGLE_SINGLE_DEG: f64 = 109.47;
const ANGLE_DOUBLE_DEG: f64 = 120.0;
const ANGLE_TRIPLE_DEG: f64 = 179.95;
const ANGLE_MOLDEN_SINGLE_DEG: f64 = 109.471;
const ANGLE_MOLDEN_NO2_DEG: f64 = 117.38;
const ANGLE_MOLDEN_CHCH_DEG: f64 = 90.0;

const DIHEDRAL_SINGLE_DEG: f64 = 0.0;
const DIHEDRAL_DOUBLE_DEG: f64 = 120.0;
const DIHEDRAL_TRIPLE_DEG: f64 = 180.0;
const DIHEDRAL_MOLDEN_DEG: f64 = 180.0;
const DIHEDRAL_MOLDEN_CYCLOPENTANE_DEG: f64 = -90.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
}

impl BondOrder {
    fn label(self) -> &'static str {
        match self {
            BondOrder::Single => "Single",
            BondOrder::Double => "Double",
            BondOrder::Triple => "Triple",
        }
    }

    fn default_angle(self) -> f64 {
        match self {
            BondOrder::Single => ANGLE_SINGLE_DEG,
            BondOrder::Double => ANGLE_DOUBLE_DEG,
            BondOrder::Triple => ANGLE_TRIPLE_DEG,
        }
    }

    fn default_dihedral(self) -> f64 {
        match self {
            BondOrder::Single => DIHEDRAL_SINGLE_DEG,
            BondOrder::Double => DIHEDRAL_DOUBLE_DEG,
            BondOrder::Triple => DIHEDRAL_TRIPLE_DEG,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FragmentInsertMode {
    Connect,
    Replace,
}

fn single_bond_radius(sym: &str) -> Option<f64> {
    match sym {
        "H" => Some(0.31),
        "C" => Some(0.76),
        "N" => Some(0.71),
        "O" => Some(0.66),
        "F" => Some(0.57),
        "Cl" => Some(1.02),
        "Br" => Some(1.20),
        "I" => Some(1.39),
        "P" => Some(1.07),
        "S" => Some(1.05),
        "Al" => Some(1.21),
        "Si" => Some(1.11),
        "Fe" => Some(1.24),
        _ => None,
    }
}

fn double_bond_radius(sym: &str) -> Option<f64> {
    match sym {
        "C" => Some(0.67),
        "N" => Some(0.60),
        "O" => Some(0.57),
        "S" => Some(0.94),
        "P" => Some(1.00),
        "Si" => Some(1.02),
        "Fe" => Some(1.16),
        _ => None,
    }
}

fn triple_bond_radius(sym: &str) -> Option<f64> {
    match sym {
        "C" => Some(0.60),
        "N" => Some(0.54),
        "P" => Some(0.94),
        _ => None,
    }
}

fn bond_length_for(
    sym_a: &str,
    sym_b: &str,
    order: BondOrder,
    fallback_a: f32,
    fallback_b: f32,
) -> f64 {
    let (a, b) = if sym_a <= sym_b { (sym_a, sym_b) } else { (sym_b, sym_a) };
    match order {
        BondOrder::Triple => {
            if a == "C" && b == "C" { return 1.20; }
            if a == "C" && b == "N" { return 1.16; }
            if a == "C" && b == "P" { return 1.55; }
            if a == "N" && b == "P" { return 1.55; }
        }
        BondOrder::Double => {}
        BondOrder::Single => {}
    }

    match order {
        BondOrder::Single => {
            if let (Some(ra), Some(rb)) = (single_bond_radius(a), single_bond_radius(b)) {
                return ra + rb;
            }
        }
        BondOrder::Double => {
            if let (Some(ra), Some(rb)) = (double_bond_radius(a), double_bond_radius(b)) {
                return ra + rb;
            }
        }
        BondOrder::Triple => {
            if let (Some(ra), Some(rb)) = (triple_bond_radius(a), triple_bond_radius(b)) {
                return ra + rb;
            }
        }
    }

    let sum = (fallback_a + fallback_b) as f64;
    match order {
        BondOrder::Single => sum,
        BondOrder::Double => sum * 0.90,
        BondOrder::Triple => sum * 0.80,
    }
}

fn clamp_cos(c: f64) -> f64 {
    if c > 1.0 {
        1.0
    } else if c < -1.0 {
        -1.0
    } else {
        c
    }
}

fn calc_angle_deg(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let ab = b - a;
    let cb = b - c;
    let denom = ab.length() * cb.length();
    if denom <= 1.0e-12 {
        return 0.0;
    }
    let cos_angle = clamp_cos((ab.dot(cb) / denom) as f64);
    cos_angle.acos().to_degrees()
}

fn pick_perpendicular(v: Vec3) -> Vec3 {
    let axis = if v.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    v.cross(axis).normalize_or_zero()
}

fn calc_dihedral_deg(p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3) -> f64 {
    let e1 = (p1 - p2).normalize_or_zero();
    let e2 = (p2 - p3).normalize_or_zero();
    let mut n = e1.cross(e2);
    if n.length() < 1.0e-6 {
        n = pick_perpendicular(e1);
    } else {
        n = n.normalize();
    }
    let m = n.cross(e1);
    let v = p4 - p1;
    (v.dot(n)).atan2(v.dot(m)).to_degrees() as f64
}

fn centroid_excluding(pos: &[Vec3], exclude: usize) -> Option<Vec3> {
    if pos.len() <= 1 {
        return None;
    }
    let mut sum = Vec3::ZERO;
    let mut count = 0;
    for (idx, p) in pos.iter().enumerate() {
        if idx == exclude {
            continue;
        }
        sum += *p;
        count += 1;
    }
    if count == 0 {
        return None;
    }
    Some(sum / count as f32)
}

fn fragment_adjacency(symbols: &[String], coords: &[Vec3]) -> Vec<Vec<usize>> {
    let n = symbols.len();
    let mut adj = vec![Vec::new(); n];
    let scale = 1.3_f32;
    for i in 0..n {
        for j in (i + 1)..n {
            let ri = covalent_radius_angstrom(&symbols[i]);
            let rj = covalent_radius_angstrom(&symbols[j]);
            let cutoff = (ri + rj) * scale;
            let d = coords[i].distance(coords[j]);
            if d > 0.01 && d <= cutoff {
                adj[i].push(j);
                adj[j].push(i);
            }
        }
    }
    adj
}

fn fragment_atom_order(symbols: &[String], coords: &[Vec3]) -> Vec<usize> {
    let n = symbols.len();
    if n == 0 {
        return Vec::new();
    }
    let adj = fragment_adjacency(symbols, coords);
    let mut order = Vec::with_capacity(n);
    let mut visited = vec![false; n];

    order.push(0);
    visited[0] = true;

    let mut remaining_heavy: Vec<usize> = (0..n)
        .filter(|&i| symbols[i] != "H" && i != 0)
        .collect();

    while !remaining_heavy.is_empty() {
        let mut best: Option<(usize, f32)> = None;
        for &idx in &remaining_heavy {
            let mut min_dist = f32::MAX;
            for &v in &order {
                if adj[idx].contains(&v) {
                    let d = coords[idx].distance(coords[v]);
                    if d < min_dist {
                        min_dist = d;
                    }
                }
            }
            if min_dist < f32::MAX {
                if best.map_or(true, |(_, bd)| min_dist < bd) {
                    best = Some((idx, min_dist));
                }
            }
        }
        let next = if let Some((idx, _)) = best {
            idx
        } else {
            remaining_heavy[0]
        };
        order.push(next);
        visited[next] = true;
        remaining_heavy.retain(|&i| i != next);
    }

    for &heavy_idx in order.clone().iter() {
        let mut hs: Vec<usize> = adj[heavy_idx]
            .iter()
            .copied()
            .filter(|&i| symbols[i] == "H" && !visited[i])
            .collect();
        hs.sort_by(|&a, &b| {
            coords[a]
                .distance(coords[heavy_idx])
                .partial_cmp(&coords[b].distance(coords[heavy_idx]))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for h in hs {
            order.push(h);
            visited[h] = true;
        }
    }

    for i in 0..n {
        if !visited[i] {
            order.push(i);
        }
    }

    order
}

fn build_fragment_zmat(symbols: &[String], coords: &[Vec3]) -> Vec<ZAtom> {
    let n = symbols.len();
    if n == 0 {
        return Vec::new();
    }
    let adj = fragment_adjacency(symbols, coords);
    let order = fragment_atom_order(symbols, coords);
    let mut index_of = vec![usize::MAX; n];
    for (pos, &idx) in order.iter().enumerate() {
        index_of[idx] = pos;
    }

    let mut zmat = Vec::with_capacity(n);
    for (pos, &idx) in order.iter().enumerate() {
        let symbol = symbols[idx].clone();
        if pos == 0 {
            zmat.push(ZAtom {
                symbol,
                bond_ref: None,
                bond_len: 0.0,
                angle_ref: None,
                angle_deg: 0.0,
                dihedral_ref: None,
                dihedral_deg: 0.0,
            });
            continue;
        }

        let mut bond_ref_idx = None;
        let mut bond_dist = f32::MAX;
        for &nb in &adj[idx] {
            let nb_pos = index_of[nb];
            if nb_pos < pos {
                let d = coords[idx].distance(coords[nb]);
                if d < bond_dist {
                    bond_dist = d;
                    bond_ref_idx = Some(nb);
                }
            }
        }
        let bond_ref_idx = bond_ref_idx.unwrap_or_else(|| order[pos - 1]);
        let bond_len = coords[idx].distance(coords[bond_ref_idx]) as f64;

        let mut angle_ref_idx = None;
        for &nb in &adj[bond_ref_idx] {
            let nb_pos = index_of[nb];
            if nb_pos < pos && nb != bond_ref_idx && nb != idx {
                angle_ref_idx = Some(nb);
                break;
            }
        }
        if angle_ref_idx.is_none() {
            angle_ref_idx = order[..pos]
                .iter()
                .copied()
                .find(|&j| j != bond_ref_idx && j != idx);
        }
        let angle_deg = angle_ref_idx
            .map(|a| calc_angle_deg(coords[idx], coords[bond_ref_idx], coords[a]))
            .unwrap_or(0.0);

        let mut dihedral_ref_idx = None;
        if let Some(angle_idx) = angle_ref_idx {
            for &nb in &adj[angle_idx] {
                let nb_pos = index_of[nb];
                if nb_pos < pos && nb != bond_ref_idx && nb != angle_idx && nb != idx {
                    dihedral_ref_idx = Some(nb);
                    break;
                }
            }
            if dihedral_ref_idx.is_none() {
                dihedral_ref_idx = order[..pos]
                    .iter()
                    .copied()
                    .find(|&j| j != bond_ref_idx && j != angle_idx && j != idx);
            }
        }
        let dihedral_deg = dihedral_ref_idx
            .zip(angle_ref_idx)
            .map(|(d, a)| calc_dihedral_deg(coords[bond_ref_idx], coords[a], coords[d], coords[idx]))
            .unwrap_or(0.0);

        zmat.push(ZAtom {
            symbol,
            bond_ref: Some(bond_ref_idx + 1),
            bond_len,
            angle_ref: angle_ref_idx.map(|v| v + 1),
            angle_deg,
            dihedral_ref: dihedral_ref_idx.map(|v| v + 1),
            dihedral_deg,
        });
    }

    zmat
}

fn build_bond_adjacency(n: usize, bonds: &[(usize, usize, f32)]) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); n];
    for &(i, j, _) in bonds {
        if i < n && j < n {
            adj[i].push(j);
            adj[j].push(i);
        }
    }
    adj
}

fn pick_neighbor(adj: &[Vec<usize>], center: usize, exclude: &[usize]) -> Option<usize> {
    adj.get(center).and_then(|nei| {
        nei.iter()
            .copied()
            .find(|idx| !exclude.contains(idx))
    })
}

fn pick_any_except(n: usize, exclude: &[usize]) -> Option<usize> {
    (0..n).find(|idx| !exclude.contains(idx))
}

fn distance_to_axis(p: Vec3, a: Vec3, b: Vec3) -> f32 {
    let axis = b - a;
    let denom = axis.length_squared();
    if denom <= 1.0e-8 {
        return (p - a).length();
    }
    let t = (p - a).dot(axis) / denom;
    let proj = a + axis * t;
    p.distance(proj)
}

fn rotate_indices_around_axis(
    coords: &[Vec3],
    a: Vec3,
    b: Vec3,
    indices: &[usize],
    angle_deg: f32,
) -> Vec<Vec3> {
    let axis = (b - a).normalize_or_zero();
    if axis.length() <= 1.0e-6 {
        return coords.to_vec();
    }
    let angle = angle_deg.to_radians();
    let mut out = coords.to_vec();
    for &i in indices {
        let v = coords[i] - a;
        let v_rot = v * angle.cos()
            + axis.cross(v) * angle.sin()
            + axis * (axis.dot(v)) * (1.0 - angle.cos());
        out[i] = a + v_rot;
    }
    out
}

fn bend_indices_towards_angle(
    coords: &[Vec3],
    a: Vec3,
    b: Vec3,
    indices: &[usize],
    picked_idx: usize,
    target_angle_deg: f32,
) -> Vec<Vec3> {
    let axis_vec = (b - a).normalize_or_zero();
    if axis_vec.length() <= 1.0e-6 {
        return coords.to_vec();
    }
    let center = a;
    let v = coords[picked_idx] - center;
    let v_len = v.length();
    if v_len <= 1.0e-6 {
        return coords.to_vec();
    }
    let v_n = v / v_len;

    let mut dot = axis_vec.dot(v_n);
    dot = dot.clamp(-1.0, 1.0);
    let current_angle = dot.acos();
    let target_angle = target_angle_deg.to_radians();
    let delta = target_angle - current_angle;
    if delta.abs() <= 1.0e-6 {
        return coords.to_vec();
    }

    let mut rot_axis = axis_vec.cross(v_n);
    if rot_axis.length() <= 1.0e-6 {
        rot_axis = axis_vec.cross(Vec3::X);
        if rot_axis.length() <= 1.0e-6 {
            rot_axis = axis_vec.cross(Vec3::Y);
        }
        if rot_axis.length() <= 1.0e-6 {
            return coords.to_vec();
        }
    }
    let rot_axis = rot_axis.normalize();

    let mut out = coords.to_vec();
    let angle = delta;
    for &i in indices {
        let p = coords[i] - center;
        let p_rot = p * angle.cos()
            + rot_axis.cross(p) * angle.sin()
            + rot_axis * (rot_axis.dot(p)) * (1.0 - angle.cos());
        out[i] = center + p_rot;
    }
    out
}

fn min_clearance_between(
    coords: &[Vec3],
    frag_indices: &[usize],
    exclude_other: &[usize],
) -> f32 {
    if coords.is_empty() || frag_indices.is_empty() {
        return 0.0;
    }
    let mut is_frag = vec![false; coords.len()];
    for &idx in frag_indices {
        if idx < is_frag.len() {
            is_frag[idx] = true;
        }
    }
    let mut min_d = f32::MAX;
    for &fi in frag_indices {
        if fi >= coords.len() || fi == frag_indices[0] {
            continue;
        }
        let fpos = coords[fi];
        for (mi, mpos) in coords.iter().enumerate() {
            if is_frag[mi] {
                continue;
            }
            if exclude_other.contains(&mi) {
                continue;
            }
            let d = fpos.distance(*mpos);
            if d < min_d {
                min_d = d;
            }
        }
    }
    min_d
}

#[derive(Default, Reflect, GizmoConfigGroup)]
pub struct BuilderGizmos;

#[derive(Resource, Default, Clone)]
pub struct EditorRotateState {
    pub active: bool,                     // capture 3D clicks when true
    pub picks: Vec<usize>,                // [A, B, side-indicator (any atom on that side)]
    pub chosen_side: Option<RotateSide>,  // inferred from 3rd pick
    pub angle_deg: f32,                   // degrees
    pub axis_len_target: f32,             // Å, absolute axis length target
    pub axis_pair: Option<(usize, usize)>,
    pub bend_angle_deg: f32,              // degrees
    pub bend_ref: Option<(usize, usize, usize)>,
    pub last_rotate_snapshot: Option<Vec<Vec3>>,    // undo buffer
    pub last_translate_snapshot: Option<Vec<Vec3>>, // undo buffer
    pub last_bend_snapshot: Option<Vec<Vec3>>,      // undo buffer
    pub bond_th_hx: f32,                  // Å
    pub bond_th_xx: f32,                  // Å
}

#[derive(Resource, Clone)]
pub struct ZMatrixBuilderState {
    pub zmat: Vec<ZAtom>,
    pub new_symbol: String,
    pub new_bond_order: BondOrder,
    pub new_angle_deg: f64,
    pub new_dihedral_deg: f64,
    pub frag_name: String,
    pub frag_bond_order: BondOrder,
    pub frag_angle_deg: f64,
    pub frag_dihedral_deg: f64,
    pub frag_pick_active: bool,
    pub frag_pick_indices: Vec<usize>,
    pub frag_pick_hint: Option<String>,
    pub frag_scan_score: Option<f32>,
    pub frag_mode: FragmentInsertMode,
    pub last_frag_snapshot: Option<Vec<ZAtom>>,
    pub frag_undo_visible: bool,
    pub pick_active: bool,
    pub pick_indices: Vec<usize>,
    pub pick_hint: Option<String>,
    pub edit_rows: Vec<ZAtomEditRow>,
    pub edit_refresh: bool,
    pub edit_preview: Vec<usize>,
    pub selected_index: Option<usize>,
    pub selected_symbol: Option<String>,
    pub remove_popup_open: bool,
    pub remove_popup_pos: Option<egui::Pos2>,
    pub last_remove_snapshot: Option<Vec<ZAtom>>,
    pub redo_remove_visible: bool,
    pub original_atoms: Option<Vec<String>>,
    pub original_pos: Option<Vec<Vec3>>,
    pub last_error: Option<String>,
}

#[derive(Clone, Default)]
pub struct ZAtomEditRow {
    pub symbol: String,
    pub bond_ref: String,
    pub bond_len: String,
    pub angle_ref: String,
    pub angle_deg: String,
    pub dihedral_ref: String,
    pub dihedral_deg: String,
}

impl ZAtomEditRow {
    fn from_atom(atom: &ZAtom, index: usize) -> Self {
        Self {
            symbol: atom.symbol.clone(),
            bond_ref: if index >= 1 {
                atom.bond_ref.unwrap_or(1).to_string()
            } else {
                String::new()
            },
            bond_len: if index >= 1 {
                format!("{:.3}", atom.bond_len)
            } else {
                String::new()
            },
            angle_ref: if index >= 2 {
                atom.angle_ref.unwrap_or(1).to_string()
            } else {
                String::new()
            },
            angle_deg: if index >= 2 {
                format!("{:.2}", atom.angle_deg)
            } else {
                String::new()
            },
            dihedral_ref: if index >= 3 {
                atom.dihedral_ref.unwrap_or(1).to_string()
            } else {
                String::new()
            },
            dihedral_deg: if index >= 3 {
                format!("{:.2}", atom.dihedral_deg)
            } else {
                String::new()
            },
        }
    }
}

impl Default for ZMatrixBuilderState {
    fn default() -> Self {
        Self {
            zmat: Vec::new(),
            new_symbol: "H".to_string(),
            new_bond_order: BondOrder::Single,
            new_angle_deg: ANGLE_SINGLE_DEG,
            new_dihedral_deg: DIHEDRAL_SINGLE_DEG,
            frag_name: fragments::FRAGMENTS
                .first()
                .map(|f| f.name.to_string())
                .unwrap_or_else(|| "CH3".to_string()),
            frag_bond_order: BondOrder::Single,
            frag_angle_deg: ANGLE_SINGLE_DEG,
            frag_dihedral_deg: DIHEDRAL_SINGLE_DEG,
            frag_pick_active: false,
            frag_pick_indices: Vec::new(),
            frag_pick_hint: None,
            frag_scan_score: None,
            frag_mode: FragmentInsertMode::Connect,
            last_frag_snapshot: None,
            frag_undo_visible: false,
            pick_active: false,
            pick_indices: Vec::new(),
            pick_hint: None,
            edit_rows: Vec::new(),
            edit_refresh: false,
            edit_preview: Vec::new(),
            selected_index: None,
            selected_symbol: None,
            remove_popup_open: false,
            remove_popup_pos: None,
            last_remove_snapshot: None,
            redo_remove_visible: false,
            original_atoms: None,
            original_pos: None,
            last_error: None,
        }
    }
}

fn color_to_egui(c: Color) -> egui::Color32 {
    let s = c.to_srgba();
    egui::Color32::from_rgba_premultiplied(
        (s.red * 255.0).round() as u8,
        (s.green * 255.0).round() as u8,
        (s.blue * 255.0).round() as u8,
        (s.alpha * 255.0).round() as u8,
    )
}

fn can_remove_zmat_index(zmat: &[ZAtom], idx: usize) -> bool {
    let target = idx + 1;
    for (i, atom) in zmat.iter().enumerate() {
        if i <= idx {
            continue;
        }
        if atom.bond_ref == Some(target)
            || atom.angle_ref == Some(target)
            || atom.dihedral_ref == Some(target)
        {
            return false;
        }
    }
    true
}

fn remove_zmat_index(zmat: &mut Vec<ZAtom>, idx: usize) {
    if idx >= zmat.len() {
        return;
    }
    zmat.remove(idx);
    let removed = idx + 1;
    for atom in zmat.iter_mut() {
        if let Some(r) = atom.bond_ref {
            if r > removed {
                atom.bond_ref = Some(r - 1);
            }
        }
        if let Some(r) = atom.angle_ref {
            if r > removed {
                atom.angle_ref = Some(r - 1);
            }
        }
        if let Some(r) = atom.dihedral_ref {
            if r > removed {
                atom.dihedral_ref = Some(r - 1);
            }
        }
    }
}

fn add_fragment_to_zmat(
    zmat: &mut Vec<ZAtom>,
    fragment_xyz: &str,
    connector_refs: Option<(usize, Option<usize>, Option<usize>)>,
    bond_len: f64,
    angle_deg: f64,
    dihedral_deg: f64,
) {
    let (_n, symbols, coords) = parse_xyz_angstrom(fragment_xyz);
    let frag_coords: Vec<Vec3> = coords
        .into_iter()
        .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
        .collect();
    let frag_atoms = build_fragment_zmat(&symbols, &frag_coords);
    if frag_atoms.is_empty() {
        return;
    }

    let offset = zmat.len();
    for (i, mut atom) in frag_atoms.into_iter().enumerate() {
        if i == 0 {
            if let Some((bond_ref, angle_ref, dihedral_ref)) = connector_refs {
                atom.bond_ref = Some(bond_ref + 1);
                atom.bond_len = bond_len;
                atom.angle_ref = angle_ref.map(|v| v + 1);
                atom.angle_deg = angle_deg;
                atom.dihedral_ref = dihedral_ref.map(|v| v + 1);
                atom.dihedral_deg = dihedral_deg;
            }
        } else {
            if let Some(r) = atom.bond_ref {
                atom.bond_ref = Some(r + offset);
            }
            if let Some(r) = atom.angle_ref {
                atom.angle_ref = Some(r + offset);
            }
            if let Some(r) = atom.dihedral_ref {
                atom.dihedral_ref = Some(r + offset);
            }
        }
        zmat.push(atom);
    }
}

fn find_fragment(name: &str) -> Option<&'static fragments::FragmentDef> {
    fragments::FRAGMENTS.iter().find(|f| f.name == name)
}

fn fragment_molden_defaults(name: &str) -> (f64, f64) {
    match name {
        "-CH3" => (ANGLE_MOLDEN_SINGLE_DEG, DIHEDRAL_MOLDEN_DEG),
        "-CH=CH2" => (ANGLE_DOUBLE_DEG, DIHEDRAL_MOLDEN_DEG),
        "-CH=O" => (ANGLE_DOUBLE_DEG, DIHEDRAL_MOLDEN_DEG),
        "-COOH" => (ANGLE_DOUBLE_DEG, DIHEDRAL_MOLDEN_DEG),
        "-NH2" => (ANGLE_MOLDEN_SINGLE_DEG, DIHEDRAL_MOLDEN_DEG),
        "-OH" => (ANGLE_MOLDEN_SINGLE_DEG, DIHEDRAL_MOLDEN_DEG),
        "-CCH" => (ANGLE_MOLDEN_CHCH_DEG, DIHEDRAL_MOLDEN_DEG),
        "-Phenyl" => (ANGLE_DOUBLE_DEG, DIHEDRAL_MOLDEN_DEG),
        "-Pyrrole" => (ANGLE_DOUBLE_DEG, DIHEDRAL_MOLDEN_DEG),
        "-OCH3" => (ANGLE_MOLDEN_SINGLE_DEG, DIHEDRAL_MOLDEN_DEG),
        "-NO2" => (ANGLE_MOLDEN_NO2_DEG, DIHEDRAL_MOLDEN_DEG),
        "-OOH" => (ANGLE_MOLDEN_SINGLE_DEG, DIHEDRAL_MOLDEN_DEG),
        "-CycloPentane" => (ANGLE_MOLDEN_SINGLE_DEG, DIHEDRAL_MOLDEN_CYCLOPENTANE_DEG),
        "-CycloHexane" => (ANGLE_MOLDEN_SINGLE_DEG, DIHEDRAL_MOLDEN_DEG),
        _ => (ANGLE_MOLDEN_SINGLE_DEG, DIHEDRAL_MOLDEN_DEG),
    }
}

fn replace_atom_with_fragment_zmat(
    zmat: &mut Vec<ZAtom>,
    frag_atoms: Vec<ZAtom>,
    replace_idx: usize,
    bond_len: f64,
    angle_deg: f64,
    dihedral_deg: f64,
) -> Vec<usize> {
    if replace_idx >= zmat.len() || frag_atoms.is_empty() {
        return Vec::new();
    }
    let connector_symbol = frag_atoms[0].symbol.clone();
    let mut target = zmat[replace_idx].clone();
    target.symbol = connector_symbol;
    if target.bond_ref.is_some() {
        target.bond_len = bond_len;
        target.angle_deg = angle_deg;
        target.dihedral_deg = dihedral_deg;
    }
    zmat[replace_idx] = target;

    let base_len = zmat.len();
    let mut map: Vec<usize> = Vec::with_capacity(frag_atoms.len());
    map.push(replace_idx);
    for i in 1..frag_atoms.len() {
        map.push(base_len + (i - 1));
    }

    for (_i, mut atom) in frag_atoms.into_iter().enumerate().skip(1) {
        if let Some(r) = atom.bond_ref {
            atom.bond_ref = Some(map[r - 1] + 1);
        }
        if let Some(r) = atom.angle_ref {
            atom.angle_ref = Some(map[r - 1] + 1);
        }
        if let Some(r) = atom.dihedral_ref {
            atom.dihedral_ref = Some(map[r - 1] + 1);
        }
        zmat.push(atom);
    }
    map
}

fn ensure_zmat_edit_buffers(zmat_state: &mut ZMatrixBuilderState) {
    if zmat_state.edit_refresh || zmat_state.edit_rows.len() != zmat_state.zmat.len() {
        zmat_state.edit_rows = zmat_state
            .zmat
            .iter()
            .enumerate()
            .map(|(i, atom)| ZAtomEditRow::from_atom(atom, i))
            .collect();
        zmat_state.edit_refresh = false;
    }
}

fn apply_zmat_edits(
    zmat_state: &mut ZMatrixBuilderState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
) {
    fn parse_usize(s: &str, label: &str, row_idx: usize) -> Result<usize, String> {
        match s.trim().parse::<usize>() {
            Ok(v) if v >= 1 => Ok(v),
            _ => Err(format!("Row {}: invalid {} reference.", row_idx + 1, label)),
        }
    }

    fn parse_f64(s: &str, label: &str, row_idx: usize) -> Result<f64, String> {
        match s.trim().parse::<f64>() {
            Ok(v) if v.is_finite() => Ok(v),
            _ => Err(format!("Row {}: invalid {} value.", row_idx + 1, label)),
        }
    }

    let mut next = Vec::with_capacity(zmat_state.edit_rows.len());
    for (i, row) in zmat_state.edit_rows.iter().enumerate() {
        let symbol = row.symbol.trim();
        if symbol.is_empty() {
            zmat_state.last_error = Some(format!("Row {}: missing atom symbol.", i + 1));
            return;
        }

        let clamp_ref = |v: usize, max_ref: usize| -> usize {
            v.max(1).min(max_ref)
        };

        let atom = if i == 0 {
            ZAtom {
                symbol: symbol.to_string(),
                bond_ref: None,
                bond_len: 0.0,
                angle_ref: None,
                angle_deg: 0.0,
                dihedral_ref: None,
                dihedral_deg: 0.0,
            }
        } else if i == 1 {
            let bond_ref = match parse_usize(&row.bond_ref, "bond", i) {
                Ok(v) => v,
                Err(err) => { zmat_state.last_error = Some(err); return; }
            };
            let bond_len = match parse_f64(&row.bond_len, "bond length", i) {
                Ok(v) => v,
                Err(err) => { zmat_state.last_error = Some(err); return; }
            };
            let bond_ref = clamp_ref(bond_ref, i);
            ZAtom {
                symbol: symbol.to_string(),
                bond_ref: Some(bond_ref),
                bond_len,
                angle_ref: None,
                angle_deg: 0.0,
                dihedral_ref: None,
                dihedral_deg: 0.0,
            }
        } else if i == 2 {
            let bond_ref = match parse_usize(&row.bond_ref, "bond", i) {
                Ok(v) => v,
                Err(err) => { zmat_state.last_error = Some(err); return; }
            };
            let bond_len = match parse_f64(&row.bond_len, "bond length", i) {
                Ok(v) => v,
                Err(err) => { zmat_state.last_error = Some(err); return; }
            };
            let angle_ref = match parse_usize(&row.angle_ref, "angle", i) {
                Ok(v) => v,
                Err(err) => { zmat_state.last_error = Some(err); return; }
            };
            let angle_deg = match parse_f64(&row.angle_deg, "angle", i) {
                Ok(v) => v,
                Err(err) => { zmat_state.last_error = Some(err); return; }
            };
            let bond_ref = clamp_ref(bond_ref, i);
            let angle_ref = clamp_ref(angle_ref, i);
            ZAtom {
                symbol: symbol.to_string(),
                bond_ref: Some(bond_ref),
                bond_len,
                angle_ref: Some(angle_ref),
                angle_deg,
                dihedral_ref: None,
                dihedral_deg: 0.0,
            }
        } else {
            let bond_ref = match parse_usize(&row.bond_ref, "bond", i) {
                Ok(v) => v,
                Err(err) => { zmat_state.last_error = Some(err); return; }
            };
            let bond_len = match parse_f64(&row.bond_len, "bond length", i) {
                Ok(v) => v,
                Err(err) => { zmat_state.last_error = Some(err); return; }
            };
            let angle_ref = match parse_usize(&row.angle_ref, "angle", i) {
                Ok(v) => v,
                Err(err) => { zmat_state.last_error = Some(err); return; }
            };
            let angle_deg = match parse_f64(&row.angle_deg, "angle", i) {
                Ok(v) => v,
                Err(err) => { zmat_state.last_error = Some(err); return; }
            };
            let dihedral_ref = match parse_usize(&row.dihedral_ref, "dihedral", i) {
                Ok(v) => v,
                Err(err) => { zmat_state.last_error = Some(err); return; }
            };
            let dihedral_deg = match parse_f64(&row.dihedral_deg, "dihedral", i) {
                Ok(v) => v,
                Err(err) => { zmat_state.last_error = Some(err); return; }
            };
            let bond_ref = clamp_ref(bond_ref, i);
            let angle_ref = clamp_ref(angle_ref, i);
            let dihedral_ref = clamp_ref(dihedral_ref, i);
            ZAtom {
                symbol: symbol.to_string(),
                bond_ref: Some(bond_ref),
                bond_len,
                angle_ref: Some(angle_ref),
                angle_deg,
                dihedral_ref: Some(dihedral_ref),
                dihedral_deg,
            }
        };
        next.push(atom);
    }

    zmat_state.zmat = next;
    let coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
    mol.atoms = zmat_state
        .zmat
        .iter()
        .map(|atom| atom.symbol.clone())
        .collect();
    mol.pos = coords;
    mol.recompute_bonds(2.0, 3.0);
    settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
    zmat_state.edit_refresh = true;
    zmat_state.last_error = None;
}

#[derive(Clone, Copy)]
enum PreviewKind {
    Bond,
    Angle,
    Dihedral,
}

fn build_preview_indices(
    row_idx: usize,
    row: &ZAtomEditRow,
    kind: PreviewKind,
    total_atoms: usize,
) -> Vec<usize> {
    let mut out = Vec::new();
    let mut push_unique = |idx: usize| {
        if idx < total_atoms && !out.contains(&idx) {
            out.push(idx);
        }
    };

    let parse_ref = |s: &str, max_ref: usize| -> Option<usize> {
        match s.trim().parse::<usize>() {
            Ok(v) if v >= 1 && v <= max_ref => Some(v - 1),
            _ => None,
        }
    };

    push_unique(row_idx);
    let bond_ref = parse_ref(&row.bond_ref, row_idx);
    if let Some(b) = bond_ref {
        push_unique(b);
    }

    if matches!(kind, PreviewKind::Angle | PreviewKind::Dihedral) {
        if let Some(a) = parse_ref(&row.angle_ref, row_idx) {
            push_unique(a);
        }
    }

    if matches!(kind, PreviewKind::Dihedral) {
        if let Some(d) = parse_ref(&row.dihedral_ref, row_idx) {
            push_unique(d);
        }
    }

    out
}

impl EditorRotateState {
    pub fn reset(&mut self) {
        *self = Self {
            active: false,
            picks: Vec::new(),
            chosen_side: None,
            angle_deg: 0.0,
            axis_len_target: 0.0,
            axis_pair: None,
            bend_angle_deg: 0.0,
            bend_ref: None,
            last_rotate_snapshot: None,
            last_translate_snapshot: None,
            last_bend_snapshot: None,
            bond_th_hx: 1.4,
            bond_th_xx: 1.8,
        };
    }
}

/// Separate gizmo styling (so we don't conflict with measurements).
pub fn configure_builder_gizmos(mut cfg_store: ResMut<GizmoConfigStore>) {
    let (cfg, _) = cfg_store.config_mut::<BuilderGizmos>();
    cfg.enabled = true;
    cfg.line.width = 16.0;
    cfg.line.perspective = true;
}

/// Consume AtomPicked events (from shared picker) for the builder.
pub fn handle_builder_atom_picked(
    mut ev: MessageReader<AtomPicked>,
    mol: Option<ResMut<Molecule>>,
    mut state: ResMut<EditorRotateState>,
    mut zmat_state: ResMut<ZMatrixBuilderState>,
    mut settings: ResMut<MolSettings>,
) {
    let Some(mut mol) = mol else { return; };

    if zmat_state.pick_active {
        for AtomPicked { index: hit_idx, tool } in ev.read().copied() {
            if tool != ToolKind::Builder { continue; }

            if zmat_state.zmat.len() != mol.atoms.len() {
                zmat_state.last_error = Some(
                    "Z-matrix is out of sync with molecule. Sync first.".to_string(),
                );
                zmat_state.pick_active = false;
                zmat_state.pick_indices.clear();
                zmat_state.pick_hint = None;
                zmat_state.edit_refresh = true;
                break;
            }
            if zmat_state.zmat.is_empty() {
                zmat_state.last_error = Some(
                    "Add the first atom before using pick-add.".to_string(),
                );
                zmat_state.pick_active = false;
                zmat_state.pick_indices.clear();
                zmat_state.pick_hint = None;
                zmat_state.edit_refresh = true;
                break;
            }

            if !zmat_state.pick_indices.contains(&hit_idx) {
                zmat_state.pick_indices.push(hit_idx);
            }

            let required_picks = zmat_state.zmat.len().min(3);
            if zmat_state.pick_indices.len() >= required_picks {
                let bond_ref = zmat_state.pick_indices[0];
                let angle_ref = zmat_state.pick_indices.get(1).copied();
                let dihedral_ref = zmat_state.pick_indices.get(2).copied();

                let symbol = zmat_state.new_symbol.trim().to_string();
                if symbol.is_empty() {
                    zmat_state.last_error = Some("Atom symbol cannot be empty.".to_string());
                } else {
                    let r_new = covalent_radius_angstrom(&symbol);
                    let r_ref = covalent_radius_angstrom(&mol.atoms[bond_ref]);
                    let bond_len = bond_length_for(
                        &symbol,
                        &mol.atoms[bond_ref],
                        zmat_state.new_bond_order,
                        r_new,
                        r_ref,
                    );
                    let angle_deg = zmat_state.new_angle_deg;
                    let dihedral_deg = zmat_state.new_dihedral_deg;

                    let entry = match required_picks {
                        1 => ZAtom {
                            symbol,
                            bond_ref: Some(bond_ref + 1),
                            bond_len,
                            angle_ref: None,
                            angle_deg: 0.0,
                            dihedral_ref: None,
                            dihedral_deg: 0.0,
                        },
                        2 => ZAtom {
                            symbol,
                            bond_ref: Some(bond_ref + 1),
                            bond_len,
                            angle_ref: angle_ref.map(|v| v + 1),
                            angle_deg,
                            dihedral_ref: None,
                            dihedral_deg: 0.0,
                        },
                        _ => ZAtom {
                            symbol,
                            bond_ref: Some(bond_ref + 1),
                            bond_len,
                            angle_ref: angle_ref.map(|v| v + 1),
                            angle_deg,
                            dihedral_ref: dihedral_ref.map(|v| v + 1),
                            dihedral_deg,
                        },
                    };
                    zmat_state.zmat.push(entry);

                    let coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
                    mol.atoms = zmat_state
                        .zmat
                        .iter()
                        .map(|atom| atom.symbol.clone())
                        .collect();
                    mol.pos = coords;
                    mol.recompute_bonds(2.0, 3.0);
                    settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                    zmat_state.last_error = None;
                }

                zmat_state.pick_active = false;
                zmat_state.pick_indices.clear();
                zmat_state.pick_hint = None;
                zmat_state.edit_refresh = true;
            }
        }
        return;
    }

    if zmat_state.frag_pick_active {
        for AtomPicked { index: hit_idx, tool } in ev.read().copied() {
            if tool != ToolKind::Builder { continue; }

            if zmat_state.zmat.len() != mol.atoms.len() {
                zmat_state.last_error = Some(
                    "Z-matrix is out of sync with molecule. Sync first.".to_string(),
                );
                zmat_state.frag_pick_active = false;
                zmat_state.frag_pick_indices.clear();
                zmat_state.frag_pick_hint = None;
                zmat_state.edit_refresh = true;
                break;
            }
            if zmat_state.zmat.is_empty() {
                zmat_state.last_error = Some(
                    "Add the first atom before using fragment pick-add.".to_string(),
                );
                zmat_state.frag_pick_active = false;
                zmat_state.frag_pick_indices.clear();
                zmat_state.frag_pick_hint = None;
                zmat_state.edit_refresh = true;
                break;
            }

            if !zmat_state.frag_pick_indices.contains(&hit_idx) {
                zmat_state.frag_pick_indices.push(hit_idx);
            }

            let required_picks = match zmat_state.frag_mode {
                FragmentInsertMode::Connect => zmat_state.zmat.len().min(3),
                FragmentInsertMode::Replace => 1,
            };
            if zmat_state.frag_pick_indices.len() >= required_picks {
                let Some(frag) = find_fragment(&zmat_state.frag_name) else {
                    zmat_state.last_error = Some("Fragment not found.".to_string());
                    break;
                };
                match zmat_state.frag_mode {
                    FragmentInsertMode::Connect => {
                        zmat_state.last_frag_snapshot = Some(zmat_state.zmat.clone());
                        zmat_state.frag_undo_visible = true;
                        let bond_ref = zmat_state.frag_pick_indices[0];
                        let angle_ref = zmat_state.frag_pick_indices.get(1).copied();
                        let dihedral_ref = zmat_state.frag_pick_indices.get(2).copied();

                        let (_n, symbols, _coords) = parse_xyz_angstrom(frag.xyz);
                        let connector_symbol = symbols.get(0).cloned().unwrap_or_else(|| "C".to_string());
                        let r_new = covalent_radius_angstrom(&connector_symbol);
                        let r_ref = covalent_radius_angstrom(&mol.atoms[bond_ref]);
                        let bond_len = bond_length_for(
                            &connector_symbol,
                            &mol.atoms[bond_ref],
                            zmat_state.frag_bond_order,
                            r_new,
                            r_ref,
                        );

                        let frag_angle_deg = zmat_state.frag_angle_deg;
                        let frag_dihedral_deg = zmat_state.frag_dihedral_deg;
                        add_fragment_to_zmat(
                            &mut zmat_state.zmat,
                            frag.xyz,
                            Some((bond_ref, angle_ref, dihedral_ref)),
                            bond_len,
                            frag_angle_deg,
                            frag_dihedral_deg,
                        );

                        let coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
                        mol.atoms = zmat_state
                            .zmat
                            .iter()
                            .map(|atom| atom.symbol.clone())
                            .collect();
                        mol.pos = coords;
                        mol.recompute_bonds(2.0, 3.0);
                        settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                        zmat_state.last_error = None;
                    }
                    FragmentInsertMode::Replace => {
                        let replace_idx = zmat_state.frag_pick_indices[0];
                        let (_n, symbols, coords) = parse_xyz_angstrom(frag.xyz);
                        if symbols.is_empty() {
                            zmat_state.last_error = Some("Fragment is empty.".to_string());
                            break;
                        }
                        let frag_coords: Vec<Vec3> = coords
                            .into_iter()
                            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
                            .collect();

                        zmat_state.last_frag_snapshot = Some(zmat_state.zmat.clone());
                        zmat_state.frag_undo_visible = true;

                        let connector_symbol = symbols[0].clone();
                        let Some(target) = zmat_state.zmat.get(replace_idx) else {
                            zmat_state.last_error = Some("Invalid replacement index.".to_string());
                            break;
                        };
                        let target_bond_ref = target.bond_ref;
                        let mut target_angle_ref = target.angle_ref;
                        let mut target_dihedral_ref = target.dihedral_ref;
                        let mut bond_len = 0.0;
                        let mut angle_deg = 0.0;
                        let mut dihedral_deg = 0.0;
                        let Some(bond_ref) = target_bond_ref else {
                            zmat_state.last_error = Some(
                                "Selected atom has no bond reference to replace.".to_string(),
                            );
                            break;
                        };
                        let ref_idx = bond_ref - 1;
                        let r_new = covalent_radius_angstrom(&connector_symbol);
                        let r_ref = covalent_radius_angstrom(&mol.atoms[ref_idx]);
                        bond_len = bond_length_for(
                            &connector_symbol,
                            &mol.atoms[ref_idx],
                            BondOrder::Single,
                            r_new,
                            r_ref,
                        );

                        let frag_atoms = build_fragment_zmat(&symbols, &frag_coords);
                        let (default_angle, default_dihedral) =
                            fragment_molden_defaults(&zmat_state.frag_name);
                        angle_deg = default_angle;
                        dihedral_deg = default_dihedral;

                        if target_angle_ref.is_none() {
                            angle_deg = 0.0;
                        }
                        if target_dihedral_ref.is_none() {
                            dihedral_deg = 0.0;
                        }
                        let frag_indices = replace_atom_with_fragment_zmat(
                            &mut zmat_state.zmat,
                            frag_atoms,
                            replace_idx,
                            bond_len,
                            angle_deg,
                            dihedral_deg,
                        );

                        if !frag_indices.is_empty() {
                            if target_angle_ref.is_none() || target_dihedral_ref.is_none() {
                                let adj = build_bond_adjacency(mol.atoms.len(), &mol.bonds);
                                if target_angle_ref.is_none() {
                                    target_angle_ref = pick_neighbor(&adj, ref_idx, &[replace_idx])
                                        .map(|v| v + 1)
                                        .or_else(|| pick_any_except(mol.atoms.len(), &[ref_idx, replace_idx]).map(|v| v + 1));
                                }
                                if target_dihedral_ref.is_none() {
                                    let angle_idx = target_angle_ref.map(|v| v - 1);
                                    if let Some(aidx) = angle_idx {
                                        target_dihedral_ref = pick_neighbor(&adj, aidx, &[ref_idx, replace_idx])
                                            .map(|v| v + 1)
                                            .or_else(|| {
                                                pick_any_except(mol.atoms.len(), &[ref_idx, replace_idx, aidx])
                                                    .map(|v| v + 1)
                                            });
                                    }
                                }
                                zmat_state.zmat[replace_idx].angle_ref = target_angle_ref;
                                zmat_state.zmat[replace_idx].dihedral_ref = target_dihedral_ref;
                            }
                            let coords_base = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
                            let atoms: Vec<String> = zmat_state
                                .zmat
                                .iter()
                                .map(|atom| atom.symbol.clone())
                                .collect();
                            let (side_a, side_b) = split_sides(
                                &atoms,
                                &coords_base,
                                ref_idx,
                                replace_idx,
                                state.bond_th_hx,
                                state.bond_th_xx,
                            );
                            let side_indices = if side_b.contains(&replace_idx) {
                                side_b
                            } else {
                                side_a
                            };
                            let mut picked_idx = replace_idx;
                            let mut best_axis_dist = -1.0_f32;
                            for &idx in &side_indices {
                                if idx == replace_idx {
                                    continue;
                                }
                                if !frag_indices.contains(&idx) {
                                    continue;
                                }
                                let d = distance_to_axis(
                                    coords_base[idx],
                                    coords_base[ref_idx],
                                    coords_base[replace_idx],
                                );
                                if d > best_axis_dist {
                                    best_axis_dist = d;
                                    picked_idx = idx;
                                }
                            }

                            let angle_candidates: Vec<f32> = (12..=36).map(|i| i as f32 * 5.0).collect();
                            let dihedral_candidates: Vec<f32> = (0..36).map(|i| i as f32 * 10.0).collect();
                            let mut best_angle = angle_deg;
                            let mut best_dihedral = dihedral_deg;
                            let mut best_clear = -1.0_f32;
                            for &cand_angle in &angle_candidates {
                                let bent = bend_indices_towards_angle(
                                    &coords_base,
                                    coords_base[ref_idx],
                                    coords_base[replace_idx],
                                    &side_indices,
                                    picked_idx,
                                    cand_angle,
                                );
                                for &cand_dihedral in &dihedral_candidates {
                                    let rotated = rotate_indices_around_axis(
                                        &bent,
                                        coords_base[ref_idx],
                                        coords_base[replace_idx],
                                        &side_indices,
                                        cand_dihedral,
                                    );
                                    let clear = min_clearance_between(
                                        &rotated,
                                        &frag_indices,
                                        &[ref_idx, replace_idx],
                                    );
                                    if clear > best_clear {
                                        best_clear = clear;
                                        best_angle = cand_angle as f64;
                                        best_dihedral = cand_dihedral as f64;
                                    }
                                }
                            }
                            zmat_state.zmat[replace_idx].angle_deg = best_angle;
                            zmat_state.zmat[replace_idx].dihedral_deg = best_dihedral;
                            zmat_state.frag_scan_score = if best_clear >= 0.0 {
                                Some(best_clear)
                            } else {
                                None
                            };
                        }

                        let coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
                        mol.atoms = zmat_state
                            .zmat
                            .iter()
                            .map(|atom| atom.symbol.clone())
                            .collect();
                        mol.pos = coords;
                        mol.recompute_bonds(2.0, 3.0);
                        settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                        zmat_state.edit_refresh = true;
                        zmat_state.last_error = None;
                    }
                }

                zmat_state.frag_pick_active = false;
                zmat_state.frag_pick_indices.clear();
                zmat_state.frag_pick_hint = None;
                zmat_state.edit_refresh = true;
            }
        }
        return;
    }

    if !state.active { return; }

    for AtomPicked { index: hit_idx, tool } in ev.read().copied() {
        if tool != ToolKind::Builder { continue; }

        if state.picks.len() < 2 {
            if !state.picks.contains(&hit_idx) {
                state.picks.push(hit_idx);
            }
        } else if state.picks.len() == 2 {
            if hit_idx != state.picks[0] && hit_idx != state.picks[1] {
                let (a, b) = (state.picks[0], state.picks[1]);
                let (side_a, side_b) = split_sides(&mol.atoms, &mol.pos, a, b, state.bond_th_hx, state.bond_th_xx);

                // If the clicked atom is clearly in A or B side, use that.
                if side_a.contains(&hit_idx) {
                    state.picks.push(hit_idx);
                    state.chosen_side = Some(RotateSide::A);
                } else if side_b.contains(&hit_idx) {
                    state.picks.push(hit_idx);
                    state.chosen_side = Some(RotateSide::B);
                } else {
                    // Robust fallback: infer side by reachability from the clicked atom.
                    // If we can reach A without passing through B → side A; else if we can
                    // reach B without passing through A → side B.
                    let reach_a = reaches_without(&mol.atoms, &mol.pos, hit_idx, a, b, state.bond_th_hx, state.bond_th_xx);
                    let reach_b = reaches_without(&mol.atoms, &mol.pos, hit_idx, b, a, state.bond_th_hx, state.bond_th_xx);

                    if reach_a {
                        state.picks.push(hit_idx);
                        state.chosen_side = Some(RotateSide::A);
                    } else if reach_b {
                        state.picks.push(hit_idx);
                        state.chosen_side = Some(RotateSide::B);
                    } else {
                        // Final fallback: pick the closer endpoint so the UI never dead-ends.
                        let da = mol.pos[hit_idx].distance(mol.pos[a]);
                        let db = mol.pos[hit_idx].distance(mol.pos[b]);
                        state.picks.push(hit_idx);
                        state.chosen_side = if da <= db { Some(RotateSide::A) } else { Some(RotateSide::B) };
                    }
                    // else: ignore click (not connected to either side by our thresholds)
                }
            }
        } else {
            // had A,B,side already → restart with new A
            state.picks.clear();
            state.chosen_side = None;
            state.picks.push(hit_idx);
        }
    }
}

/// Draw pending highlights and axis line using the same color/radius style as measurements.
pub fn draw_builder_highlights(
    mut gizmos: Gizmos<BuilderGizmos>,
    state: Res<EditorRotateState>,
    zmat_state: Res<ZMatrixBuilderState>,
    mol: Option<Res<Molecule>>,
    settings: Res<MolSettings>,
) {
    let Some(mol) = mol else { return; };
    if mol.pos.is_empty() { return; }

    let magenta = Color::srgba(1.0, 0.0, 0.8, 1.0);
    let cyan    = Color::srgba(0.0, 0.8, 1.0, 1.0);
    let green   = Color::srgba(0.2, 1.0, 0.2, 1.0);

    let draw_hl = |gizmos: &mut Gizmos<BuilderGizmos>, idx: usize, color: Color| {
        if idx < mol.pos.len() {
            let base_r = covalent_radius_angstrom(&mol.atoms[idx]) * settings.atom_scale;
            let r = (base_r * 1.25).max(0.15);
            gizmos.sphere(mol.pos[idx], r, color);
        }
    };

    if state.active {
        if let Some(&i0) = state.picks.get(0) {
            draw_hl(&mut gizmos, i0, magenta);
        }
        if let Some(&i1) = state.picks.get(1) {
            draw_hl(&mut gizmos, i1, cyan);
            if let Some(&i0) = state.picks.get(0) {
                gizmos.line(mol.pos[i0], mol.pos[i1], Color::WHITE);
            }
        }
        if let Some(&i2) = state.picks.get(2) {
            draw_hl(&mut gizmos, i2, green);
        }
    }

    if zmat_state.pick_active {
        if let Some(&i0) = zmat_state.pick_indices.get(0) {
            draw_hl(&mut gizmos, i0, magenta);
        }
        if let Some(&i1) = zmat_state.pick_indices.get(1) {
            draw_hl(&mut gizmos, i1, cyan);
            if let Some(&i0) = zmat_state.pick_indices.get(0) {
                gizmos.line(mol.pos[i0], mol.pos[i1], Color::WHITE);
            }
        }
        if let Some(&i2) = zmat_state.pick_indices.get(2) {
            draw_hl(&mut gizmos, i2, green);
        }
    }
    if zmat_state.frag_pick_active {
        if let Some(&i0) = zmat_state.frag_pick_indices.get(0) {
            draw_hl(&mut gizmos, i0, magenta);
        }
        if let Some(&i1) = zmat_state.frag_pick_indices.get(1) {
            draw_hl(&mut gizmos, i1, cyan);
            if let Some(&i0) = zmat_state.frag_pick_indices.get(0) {
                gizmos.line(mol.pos[i0], mol.pos[i1], Color::WHITE);
            }
        }
        if let Some(&i2) = zmat_state.frag_pick_indices.get(2) {
            draw_hl(&mut gizmos, i2, green);
        }
    }

    if !zmat_state.edit_preview.is_empty() {
        let orange = Color::srgba(1.0, 0.6, 0.1, 1.0);
        let colors = [magenta, cyan, green, orange];
        for (i, idx) in zmat_state.edit_preview.iter().copied().enumerate() {
            let color = colors.get(i).copied().unwrap_or(orange);
            draw_hl(&mut gizmos, idx, color);
        }
    }
}

/// Left-side, resizable SidePanel with a collapsible “Molecule Editor” section.
/// NOTE: We now do a clean jump: no immediate atom transform sync; instead we
/// update `mol.pos` + `mol.bonds` and flag geometry/topology dirty so the scene
/// rebuilds atoms & bonds together next frame (no mismatched frame).
pub fn builder_ui_panel(
    mut contexts: EguiContexts,
    mut state: ResMut<EditorRotateState>,
    mut zmat_state: ResMut<ZMatrixBuilderState>,
    mut mol: ResMut<Molecule>,
    mut settings: ResMut<MolSettings>, // mark scene dirty to rebuild
) {
    let Ok(ctx) = contexts.ctx_mut() else { return; };

    egui::SidePanel::left("molecule_editor_panel")
        .resizable(true)
        .show(&ctx, |ui| {
            if zmat_state.original_atoms.is_none() && !mol.atoms.is_empty() {
                zmat_state.original_atoms = Some(mol.atoms.clone());
                zmat_state.original_pos = Some(mol.pos.clone());
            }
            ui.add_space(8.0);
            ui.collapsing("Molecule Builder", |ui| {
                ui.horizontal(|ui| {
                    if ui.button("Sync from Molecule").clicked() {
                        zmat_state.zmat = zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos);
                        zmat_state.pick_active = false;
                        zmat_state.pick_indices.clear();
                        zmat_state.pick_hint = None;
                        zmat_state.edit_refresh = true;
                        zmat_state.last_error = None;
                    }
                    if ui.button("Clear Display").clicked() {
                        mol.atoms.clear();
                        mol.pos.clear();
                        mol.bonds.clear();
                        mol.hydrogen_bonds.clear();
                        settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                        zmat_state.zmat.clear();
                        zmat_state.pick_indices.clear();
                        zmat_state.pick_active = false;
                        zmat_state.pick_hint = None;
                        zmat_state.edit_refresh = true;
                        zmat_state.last_error = None;
                    }
                    if ui.button("Load Original").clicked() {
                        if let (Some(atoms), Some(pos)) = (
                            zmat_state.original_atoms.clone(),
                            zmat_state.original_pos.clone(),
                        ) {
                            mol.atoms = atoms;
                            mol.pos = pos;
                            mol.recompute_bonds(2.0, 3.0);
                            settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                            zmat_state.zmat = zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos);
                            zmat_state.pick_indices.clear();
                            zmat_state.pick_active = false;
                            zmat_state.pick_hint = None;
                            zmat_state.edit_refresh = true;
                            zmat_state.last_error = None;
                        } else {
                            zmat_state.last_error = Some("Original structure not captured.".to_string());
                        }
                    }
                });

                ui.add_space(6.0);
                ui.separator();
                ui.add_space(6.0);

                let row_count = zmat_state.zmat.len();
                ui.label(format!("Z-matrix rows: {}", row_count));
                ui.weak("References are 1-based, targeting previous rows.");
                if row_count > 0 {
                    ensure_zmat_edit_buffers(&mut zmat_state);
                    zmat_state.edit_preview.clear();
                    let mut preview_indices: Option<Vec<usize>> = None;
                    let mut grid_rect: Option<egui::Rect> = None;
                    let mut suppress_clear = false;
                    let col_idx = [20.0, 28.0, 28.0, 68.0, 28.0, 60.0, 28.0, 60.0];
                    let row_height = 22.0;
                    let max_height = row_height * 15.0;
                    egui::ScrollArea::vertical()
                        .max_height(max_height)
                        .show(ui, |ui| {
                            let grid_resp = egui::Grid::new("zmat_rows")
                                .striped(true)
                                .spacing(egui::vec2(4.0, 4.0))
                                .show(ui, |ui| {
                                    ui.add_sized([col_idx[0], 0.0], egui::Label::new("#"));
                                    ui.add_sized([col_idx[1], 0.0], egui::Label::new("Atom"));
                                    ui.add_sized([col_idx[2], 0.0], egui::Label::new(""));
                                    ui.add_sized([col_idx[3], 0.0], egui::Label::new("Bond"));
                                    ui.add_sized([col_idx[4], 0.0], egui::Label::new(""));
                                    ui.add_sized([col_idx[5], 0.0], egui::Label::new("Angle"));
                                    ui.add_sized([col_idx[6], 0.0], egui::Label::new(""));
                                    ui.add_sized([col_idx[7], 0.0], egui::Label::new("Dihedral"));
                                    ui.end_row();

                                    let selected_idx = zmat_state.selected_index;
                                    let selected_sym = zmat_state.selected_symbol.clone();
                                    #[derive(Copy, Clone, Debug, PartialEq, Eq)]
                                    enum SelectionAction {
                                        Clear,
                                    }
                                    let mut pending_select: Option<(usize, String)> = None;
                                    let mut pending_popup_pos: Option<egui::Pos2> = None;
                                    let mut pending_open_popup = false;
                                    let mut pending_action: Option<SelectionAction> = None;

                                    for (i, row) in zmat_state.edit_rows.iter_mut().enumerate() {
                                        let is_selected = selected_idx == Some(i);
                                        let same_element = selected_sym
                                            .as_deref()
                                            .map(|s| s == row.symbol.as_str())
                                            .unwrap_or(false);
                                        let sym_color = if !row.symbol.trim().is_empty() {
                                            settings
                                                .element_colors
                                                .get(row.symbol.trim())
                                                .copied()
                                                .unwrap_or_else(|| color_for(row.symbol.trim(), settings.scheme))
                                        } else {
                                            Color::WHITE
                                        };
                                        let mut hl_bg: Option<egui::Color32> = None;
                                        if same_element && is_selected {
                                            let c = color_to_egui(sym_color);
                                            hl_bg = Some(egui::Color32::from_rgba_premultiplied(
                                                c.r(),
                                                c.g(),
                                                c.b(),
                                                if is_selected { 90 } else { 40 },
                                            ));
                                        } else if is_selected {
                                            hl_bg = Some(egui::Color32::from_rgba_premultiplied(200, 200, 200, 60));
                                        }
                                        let idx_label = format!("{}", i + 1);
                                        let idx_resp = ui.add_sized(
                                            [col_idx[0], 0.0],
                                            egui::Button::new(idx_label).selected(is_selected),
                                        );
                                        if idx_resp.clicked() {
                                            if selected_idx == Some(i) {
                                                pending_action = Some(SelectionAction::Clear);
                                            } else {
                                                pending_select = Some((i, row.symbol.clone()));
                                            }
                                        }
                                        if idx_resp.secondary_clicked() {
                                            pending_select = Some((i, row.symbol.clone()));
                                            pending_popup_pos = Some(idx_resp.rect.right_top());
                                            pending_open_popup = true;
                                        }

                                        let sym_edit = egui::TextEdit::singleline(&mut row.symbol)
                                            .desired_width(col_idx[1] - 4.0);
                                        let mut sym_changed = false;
                                        let sym_resp = if let Some(bg) = hl_bg {
                                            let mut resp = None;
                                            egui::Frame::NONE.fill(bg).show(ui, |ui| {
                                                let r = ui.add_sized([col_idx[1], 0.0], sym_edit);
                                                if r.changed() {
                                                    sym_changed = true;
                                                }
                                                resp = Some(r);
                                            });
                                            resp.unwrap()
                                        } else {
                                            let r = ui.add_sized([col_idx[1], 0.0], sym_edit);
                                            if r.changed() {
                                                sym_changed = true;
                                            }
                                            r
                                        };
                                        if sym_resp.clicked() {
                                            if selected_idx == Some(i) {
                                                pending_action = Some(SelectionAction::Clear);
                                            } else {
                                                pending_select = Some((i, row.symbol.clone()));
                                            }
                                        }
                                        if sym_resp.secondary_clicked() {
                                            pending_select = Some((i, row.symbol.clone()));
                                            pending_popup_pos = Some(sym_resp.rect.right_top());
                                            pending_open_popup = true;
                                        }
                                        if sym_changed {
                                            let trimmed = row.symbol.trim();
                                            if !trimmed.is_empty() {
                                                row.symbol = trimmed.to_string();
                                            }
                                        }
                                        if selected_idx == Some(i)
                                            && pending_action.is_none()
                                            && pending_select.is_none()
                                            && sym_changed
                                        {
                                            pending_select = Some((i, row.symbol.clone()));
                                        }
                                        if i >= 1 {
                                            let bond_ref_resp = ui
                                                .add_sized([col_idx[2], 0.0], egui::TextEdit::singleline(&mut row.bond_ref))
                                                ;
                                            if bond_ref_resp.changed() {
                                                row.bond_ref = row.bond_ref.trim().to_string();
                                            }
                                            if bond_ref_resp.has_focus() {
                                                preview_indices = Some(build_preview_indices(
                                                    i,
                                                    row,
                                                    PreviewKind::Bond,
                                                    row_count,
                                                ));
                                            }

                                            let bond_len_resp = ui
                                                .add_sized([col_idx[3], 0.0], egui::TextEdit::singleline(&mut row.bond_len))
                                                ;
                                            if bond_len_resp.changed() {
                                                row.bond_len = row.bond_len.trim().to_string();
                                            }
                                            if bond_len_resp.has_focus() {
                                                preview_indices = Some(build_preview_indices(
                                                    i,
                                                    row,
                                                    PreviewKind::Bond,
                                                    row_count,
                                                ));
                                            }
                                        } else {
                                            ui.add_sized([col_idx[2], 0.0], egui::Label::new("—"));
                                            ui.add_sized([col_idx[3], 0.0], egui::Label::new("—"));
                                        }
                                        if i >= 2 {
                                            let angle_ref_resp = ui
                                                .add_sized([col_idx[4], 0.0], egui::TextEdit::singleline(&mut row.angle_ref))
                                                ;
                                            if angle_ref_resp.changed() {
                                                row.angle_ref = row.angle_ref.trim().to_string();
                                            }
                                            if angle_ref_resp.has_focus() {
                                                preview_indices = Some(build_preview_indices(
                                                    i,
                                                    row,
                                                    PreviewKind::Angle,
                                                    row_count,
                                                ));
                                            }

                                            let angle_deg_resp = ui
                                                .add_sized([col_idx[5], 0.0], egui::TextEdit::singleline(&mut row.angle_deg))
                                                ;
                                            if angle_deg_resp.changed() {
                                                row.angle_deg = row.angle_deg.trim().to_string();
                                            }
                                            if angle_deg_resp.has_focus() {
                                                preview_indices = Some(build_preview_indices(
                                                    i,
                                                    row,
                                                    PreviewKind::Angle,
                                                    row_count,
                                                ));
                                            }
                                        } else {
                                            ui.add_sized([col_idx[4], 0.0], egui::Label::new("—"));
                                            ui.add_sized([col_idx[5], 0.0], egui::Label::new("—"));
                                        }
                                        if i >= 3 {
                                            let dihedral_ref_resp = ui
                                                .add_sized([col_idx[6], 0.0], egui::TextEdit::singleline(&mut row.dihedral_ref))
                                                ;
                                            if dihedral_ref_resp.changed() {
                                                row.dihedral_ref = row.dihedral_ref.trim().to_string();
                                            }
                                            if dihedral_ref_resp.has_focus() {
                                                preview_indices = Some(build_preview_indices(
                                                    i,
                                                    row,
                                                    PreviewKind::Dihedral,
                                                    row_count,
                                                ));
                                            }

                                            let dihedral_deg_resp = ui
                                                .add_sized([col_idx[7], 0.0], egui::TextEdit::singleline(&mut row.dihedral_deg))
                                                ;
                                            if dihedral_deg_resp.changed() {
                                                row.dihedral_deg = row.dihedral_deg.trim().to_string();
                                            }
                                            if dihedral_deg_resp.has_focus() {
                                                preview_indices = Some(build_preview_indices(
                                                    i,
                                                    row,
                                                    PreviewKind::Dihedral,
                                                    row_count,
                                                ));
                                            }
                                        } else {
                                            ui.add_sized([col_idx[6], 0.0], egui::Label::new("—"));
                                            ui.add_sized([col_idx[7], 0.0], egui::Label::new("—"));
                                        }
                                        ui.end_row();
                                    }

                                    if let Some(SelectionAction::Clear) = pending_action {
                                        zmat_state.selected_index = None;
                                        zmat_state.selected_symbol = None;
                                        zmat_state.edit_preview.clear();
                                        zmat_state.remove_popup_open = false;
                                        zmat_state.remove_popup_pos = None;
                                        zmat_state.redo_remove_visible = false;
                                        zmat_state.frag_undo_visible = false;
                                        zmat_state.last_error = None;
                                    } else if let Some((idx, sym)) = pending_select {
                                        zmat_state.selected_index = Some(idx);
                                        zmat_state.selected_symbol = Some(sym);
                                        zmat_state.edit_preview = vec![idx];
                                        if pending_open_popup {
                                            zmat_state.remove_popup_open = true;
                                            zmat_state.remove_popup_pos = pending_popup_pos;
                                        }
                                        zmat_state.redo_remove_visible = false;
                                        zmat_state.frag_undo_visible = false;
                                        zmat_state.last_error = None;
                                    }
                                });
                            grid_rect = Some(grid_resp.response.rect);
                        });
                    if let Some(preview) = preview_indices {
                        zmat_state.edit_preview = preview;
                    } else if let Some(sel) = zmat_state.selected_index {
                        zmat_state.edit_preview = vec![sel];
                    }
                    if ui.button("Apply Changes").clicked() {
                        apply_zmat_edits(&mut zmat_state, &mut mol, &mut settings);
                        zmat_state.last_error = None;
                        suppress_clear = true;
                    }
                    if let Some(idx) = zmat_state.selected_index {
                        if can_remove_zmat_index(&zmat_state.zmat, idx) {
                            if ui.button("Remove Atom").clicked() {
                                suppress_clear = true;
                                zmat_state.last_remove_snapshot = Some(zmat_state.zmat.clone());
                                remove_zmat_index(&mut zmat_state.zmat, idx);
                                let coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
                                mol.atoms = zmat_state
                                    .zmat
                                    .iter()
                                    .map(|atom| atom.symbol.clone())
                                    .collect();
                                mol.pos = coords;
                                mol.recompute_bonds(2.0, 3.0);
                                settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                                zmat_state.edit_refresh = true;
                                zmat_state.selected_index = None;
                                zmat_state.selected_symbol = None;
                                zmat_state.edit_preview.clear();
                                zmat_state.redo_remove_visible = true;
                                zmat_state.frag_undo_visible = false;
                                zmat_state.last_error = None;
                            }
                        } else {
                            zmat_state.last_error =
                                Some("Cannot remove: referenced by later rows.".to_string());
                        }
                    }
                    if zmat_state.redo_remove_visible {
                        if ui.button("Undo Remove Atom").clicked() {
                            suppress_clear = true;
                            if let Some(prev) = zmat_state.last_remove_snapshot.take() {
                                zmat_state.zmat = prev;
                                let coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
                                mol.atoms = zmat_state
                                    .zmat
                                    .iter()
                                    .map(|atom| atom.symbol.clone())
                                    .collect();
                                mol.pos = coords;
                                mol.recompute_bonds(2.0, 3.0);
                                settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                                zmat_state.edit_refresh = true;
                            }
                            zmat_state.redo_remove_visible = false;
                            zmat_state.frag_undo_visible = false;
                        }
                    }
                    if let Some(err) = &zmat_state.last_error {
                        ui.colored_label(egui::Color32::LIGHT_RED, err);
                    }

                    if let Some(rect) = grid_rect {
                        if !zmat_state.remove_popup_open
                            && !suppress_clear
                            && ctx.input(|i| {
                                i.pointer.any_click()
                                    && i.pointer.latest_pos().is_some_and(|p| !rect.contains(p))
                            })
                        {
                            zmat_state.selected_index = None;
                            zmat_state.selected_symbol = None;
                            zmat_state.edit_preview.clear();
                            zmat_state.redo_remove_visible = false;
                            zmat_state.frag_undo_visible = false;
                            zmat_state.last_error = None;
                        }
                    }
                }

                ui.add_space(6.0);
                ui.separator();
                ui.add_space(6.0);

                ui.label("Add atom (pick references)");
                ui.horizontal(|ui| {
                    let row_h = ui.spacing().interact_size.y;
                    ui.add_sized([0.0, row_h], egui::Label::new("Symbol"));
                    ui.add_sized(
                        [32.0, row_h],
                        egui::TextEdit::singleline(&mut zmat_state.new_symbol),
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Bond");
                    egui::ComboBox::from_id_salt("zmat_new_bond_order")
                        .selected_text(zmat_state.new_bond_order.label())
                        .show_ui(ui, |ui| {
                            for &bo in &[BondOrder::Single, BondOrder::Double, BondOrder::Triple] {
                                if ui
                                    .selectable_value(
                                        &mut zmat_state.new_bond_order,
                                        bo,
                                        bo.label(),
                                    )
                                    .clicked()
                                {
                                    zmat_state.new_angle_deg = bo.default_angle();
                                    zmat_state.new_dihedral_deg = bo.default_dihedral();
                                }
                            }
                        });
                });
                ui.horizontal(|ui| {
                    ui.label("Angle (deg)");
                    ui.add(egui::DragValue::new(&mut zmat_state.new_angle_deg).speed(0.1));
                    ui.add_space(12.0);
                    ui.label("Dihedral (deg)");
                    ui.add(egui::DragValue::new(&mut zmat_state.new_dihedral_deg).speed(0.1));
                });

                ui.horizontal(|ui| {
                    let pick_on = zmat_state.pick_active;
                    let sel = ui.visuals().selection.bg_fill;
                    let dim = ui.visuals().widgets.inactive.bg_fill;
                    if ui
                        .add(egui::Button::new("Add Atom").fill(if pick_on { sel } else { dim }))
                        .clicked()
                    {
                        if zmat_state.zmat.is_empty() {
                            let symbol = zmat_state.new_symbol.trim().to_string();
                            if symbol.is_empty() {
                                zmat_state.last_error = Some("Atom symbol cannot be empty.".to_string());
                            } else {
                                zmat_state.zmat.push(ZAtom {
                                    symbol,
                                    bond_ref: None,
                                    bond_len: 0.0,
                                    angle_ref: None,
                                    angle_deg: 0.0,
                                    dihedral_ref: None,
                                    dihedral_deg: 0.0,
                                });
                                let coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
                                mol.atoms = zmat_state
                                    .zmat
                                    .iter()
                                    .map(|atom| atom.symbol.clone())
                                    .collect();
                                mol.pos = coords;
                                mol.recompute_bonds(2.0, 3.0);
                                settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                                zmat_state.edit_refresh = true;
                                zmat_state.last_error = None;
                            }
                        } else {
                            zmat_state.pick_active = true;
                            zmat_state.pick_indices.clear();
                            zmat_state.frag_pick_active = false;
                            zmat_state.frag_pick_indices.clear();
                            let required_picks = zmat_state.zmat.len().min(3);
                            let hint = match required_picks {
                                1 => "Pick 1 reference atom for the bond.",
                                2 => "Pick 2 reference atoms: bond then angle.",
                                _ => "Pick 3 reference atoms: bond, angle, dihedral.",
                            };
                            zmat_state.pick_hint = Some(hint.to_string());
                            zmat_state.last_error = None;
                        }
                    }
                    if zmat_state.pick_active {
                        if ui.button("Cancel Picking").clicked() {
                            zmat_state.pick_active = false;
                            zmat_state.pick_indices.clear();
                            zmat_state.pick_hint = None;
                            zmat_state.last_error = None;
                        }
                    }
                    if ui.button("Remove Last").clicked() {
                        if zmat_state.zmat.pop().is_some() {
                            let coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
                            mol.atoms = zmat_state
                                .zmat
                                .iter()
                                .map(|atom| atom.symbol.clone())
                                .collect();
                            mol.pos = coords;
                            mol.recompute_bonds(2.0, 3.0);
                            settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                            zmat_state.edit_refresh = true;
                        }
                    }
                });

                ui.add_space(8.0);
                ui.separator();
                ui.add_space(6.0);

                ui.label("Add fragment");
                ui.horizontal(|ui| {
                    ui.label("Mode");
                    if ui
                        .selectable_value(
                            &mut zmat_state.frag_mode,
                            FragmentInsertMode::Connect,
                            "Connect",
                        )
                        .clicked()
                    {
                        zmat_state.frag_pick_active = false;
                        zmat_state.frag_pick_indices.clear();
                        zmat_state.frag_pick_hint = None;
                        zmat_state.last_error = None;
                    }
                    if ui
                        .selectable_value(
                            &mut zmat_state.frag_mode,
                            FragmentInsertMode::Replace,
                            "Replace",
                        )
                        .clicked()
                    {
                        zmat_state.frag_pick_active = false;
                        zmat_state.frag_pick_indices.clear();
                        zmat_state.frag_pick_hint = None;
                        zmat_state.last_error = None;
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("Fragment");
                    egui::ComboBox::from_id_salt("frag_combo")
                        .selected_text(zmat_state.frag_name.clone())
                        .show_ui(ui, |ui| {
                            for frag in fragments::FRAGMENTS {
                                ui.selectable_value(
                                    &mut zmat_state.frag_name,
                                    frag.name.to_string(),
                                    frag.name,
                                );
                            }
                        });
                });
                if matches!(zmat_state.frag_mode, FragmentInsertMode::Connect) {
                    ui.horizontal(|ui| {
                        ui.label("Bond");
                        egui::ComboBox::from_id_salt("frag_bond_order")
                            .selected_text(zmat_state.frag_bond_order.label())
                            .show_ui(ui, |ui| {
                                for &bo in &[BondOrder::Single, BondOrder::Double, BondOrder::Triple] {
                                    if ui
                                        .selectable_value(
                                            &mut zmat_state.frag_bond_order,
                                            bo,
                                            bo.label(),
                                        )
                                        .clicked()
                                    {
                                        zmat_state.frag_angle_deg = bo.default_angle();
                                        zmat_state.frag_dihedral_deg = bo.default_dihedral();
                                    }
                                }
                            });
                    });
                    ui.horizontal(|ui| {
                        ui.label("Angle (deg)");
                        ui.add(egui::DragValue::new(&mut zmat_state.frag_angle_deg).speed(0.1));
                        ui.add_space(12.0);
                        ui.label("Dihedral (deg)");
                        ui.add(egui::DragValue::new(&mut zmat_state.frag_dihedral_deg).speed(0.1));
                    });
                }

                ui.horizontal(|ui| {
                    let pick_on = zmat_state.frag_pick_active;
                    let sel = ui.visuals().selection.bg_fill;
                    let dim = ui.visuals().widgets.inactive.bg_fill;
                    if ui
                        .add(egui::Button::new("Add Fragment").fill(if pick_on { sel } else { dim }))
                        .clicked()
                    {
                        if zmat_state.zmat.is_empty() {
                            if let Some(frag) = find_fragment(&zmat_state.frag_name) {
                                zmat_state.last_frag_snapshot = Some(zmat_state.zmat.clone());
                                zmat_state.frag_undo_visible = true;
                                let frag_angle_deg = zmat_state.frag_angle_deg;
                                let frag_dihedral_deg = zmat_state.frag_dihedral_deg;
                                add_fragment_to_zmat(
                                    &mut zmat_state.zmat,
                                    frag.xyz,
                                    None,
                                    0.0,
                                    frag_angle_deg,
                                    frag_dihedral_deg,
                                );
                                let coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
                                mol.atoms = zmat_state
                                    .zmat
                                    .iter()
                                    .map(|atom| atom.symbol.clone())
                                    .collect();
                                mol.pos = coords;
                                mol.recompute_bonds(2.0, 3.0);
                                settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                                zmat_state.edit_refresh = true;
                                zmat_state.last_error = None;
                            } else {
                                zmat_state.last_error = Some("Fragment not found.".to_string());
                            }
                        } else {
                            zmat_state.frag_pick_active = true;
                            zmat_state.frag_pick_indices.clear();
                            zmat_state.pick_active = false;
                            zmat_state.pick_indices.clear();
                            let required_picks = match zmat_state.frag_mode {
                                FragmentInsertMode::Connect => zmat_state.zmat.len().min(3),
                                FragmentInsertMode::Replace => 1,
                            };
                            let hint = match required_picks {
                                1 if matches!(zmat_state.frag_mode, FragmentInsertMode::Replace) => {
                                    "Pick 1 atom to replace."
                                }
                                1 => "Pick 1 reference atom for the bond.",
                                2 => "Pick 2 reference atoms: bond then angle.",
                                _ => "Pick 3 reference atoms: bond, angle, dihedral.",
                            };
                            zmat_state.frag_pick_hint = Some(hint.to_string());
                            zmat_state.last_error = None;
                        }
                    }
                    if zmat_state.frag_pick_active {
                        if ui.button("Cancel Fragment Picking").clicked() {
                            zmat_state.frag_pick_active = false;
                            zmat_state.frag_pick_indices.clear();
                            zmat_state.frag_pick_hint = None;
                            zmat_state.last_error = None;
                        }
                    }
                });

                if let Some(hint) = &zmat_state.frag_pick_hint {
                    ui.weak(hint);
                }
                if let Some(score) = zmat_state.frag_scan_score {
                    ui.weak(format!("Fragment scan min distance: {:.3} Å", score));
                }
                if zmat_state.frag_undo_visible {
                    if ui.button("Undo Fragment").clicked() {
                        if let Some(snapshot) = zmat_state.last_frag_snapshot.take() {
                            zmat_state.zmat = snapshot;
                            let coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
                            mol.atoms = zmat_state
                                .zmat
                                .iter()
                                .map(|atom| atom.symbol.clone())
                                .collect();
                            mol.pos = coords;
                            mol.recompute_bonds(2.0, 3.0);
                            settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                            zmat_state.edit_refresh = true;
                        }
                        zmat_state.frag_pick_active = false;
                        zmat_state.frag_pick_indices.clear();
                        zmat_state.frag_pick_hint = None;
                        zmat_state.frag_undo_visible = false;
                        zmat_state.last_error = None;
                    }
                }

                if let Some(hint) = &zmat_state.pick_hint {
                    ui.weak(hint);
                }

                if zmat_state.pick_active {
                    let required_picks = zmat_state.zmat.len().min(3);
                    let hint = match required_picks {
                        1 => "Pick bond reference atom.",
                        2 => "Pick bond reference, then angle reference.",
                        _ => "Pick bond, angle, then dihedral reference atoms.",
                    };
                    ui.weak(hint);
                }

            });

    if zmat_state.remove_popup_open {
        let popup_pos = zmat_state
            .remove_popup_pos
            .unwrap_or(egui::pos2(40.0, 40.0));
        let popup_id = egui::Id::new("zmat_remove_popup");
        egui::Area::new(popup_id)
            .fixed_pos(popup_pos)
            .order(egui::Order::Foreground)
            .show(&ctx, |ui| {
                egui::Frame::popup(&ctx.style()).show(ui, |ui| {
                    ui.label("Remove selected atom?");
                    let Some(idx) = zmat_state.selected_index else {
                        if ui.button("Close").clicked() {
                            zmat_state.remove_popup_open = false;
                        }
                        return;
                    };

                    let can_remove = can_remove_zmat_index(&zmat_state.zmat, idx);
                    if !can_remove {
                        ui.colored_label(
                            egui::Color32::LIGHT_RED,
                            "Cannot remove: referenced by later rows.",
                        );
                    }
                    if ui.add_enabled(can_remove, egui::Button::new("Remove")).clicked() {
                        remove_zmat_index(&mut zmat_state.zmat, idx);
                        let coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
                        mol.atoms = zmat_state
                            .zmat
                            .iter()
                            .map(|atom| atom.symbol.clone())
                            .collect();
                        mol.pos = coords;
                        mol.recompute_bonds(2.0, 3.0);
                        settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                        zmat_state.edit_refresh = true;
                        zmat_state.selected_index = None;
                        zmat_state.selected_symbol = None;
                        zmat_state.edit_preview.clear();
                        zmat_state.remove_popup_open = false;
                        zmat_state.remove_popup_pos = None;
                    }
                    if ui.button("Close").clicked() {
                        zmat_state.remove_popup_open = false;
                    }
                });
            });
    }

            ui.add_space(8.0);
            ui.separator();
            ui.add_space(8.0);
            ui.collapsing("Molecule Editor", |ui| {
                ui.horizontal(|ui| {
                    if ui.button(if state.active { "Deactivate" } else { "Activate" }).clicked() {
                        state.active = !state.active;
                        if !state.active {
                            state.picks.clear();
                            state.chosen_side = None;
                        }
                    }
                    if ui.button("Reset").clicked() {
                        state.reset();
                    }
                });

                ui.separator();
                ui.label("Bond-axis Fragment Rotation");
                ui.monospace(format!("Picks: {:?}", state.picks));

                if state.picks.len() < 2 {
                    ui.weak("Pick two atoms in 3D to define the bond axis (A then B).");
                } else if state.picks.len() == 2 {
                    ui.weak("Pick a third atom on the side you want to rotate.");
                } else {
                    match state.chosen_side {
                        Some(RotateSide::A) => ui.colored_label(egui::Color32::LIGHT_GREEN, "Chosen side: A"),
                        Some(RotateSide::B) => ui.colored_label(egui::Color32::LIGHT_GREEN, "Chosen side: B"),
                        None => ui.weak("Chosen side: —"),
                    };
                }

                ui.add_space(6.0);
                ui.add(egui::Slider::new(&mut state.angle_deg, -180.0..=180.0).text("Angle (°)"));
                ui.horizontal(|ui| {
                    if ui.button("Rotate").clicked() {
                        if state.picks.len() >= 3 {
                            let a = state.picks[0];
                            let b = state.picks[1];
                            if let Some(side) = state.chosen_side {
                                // Save undo snapshot
                                state.last_rotate_snapshot = Some(mol.pos.clone());

                                // Rotate fragment coordinates
                                mol.pos = rotate_side(
                                    &mol.atoms, &mol.pos, a, b,
                                    side, state.angle_deg,
                                    state.bond_th_hx, state.bond_th_xx,
                                );

                                settings.coords_dirty = true;
                            }
                        }
                    }
                    if ui.button("Undo").clicked() {
                        if let Some(snapshot) = &mut state.last_rotate_snapshot {
                            mol.pos = snapshot.clone();
                            settings.coords_dirty = true;
                        }
                    }
                });

                ui.add_space(6.0);
                ui.separator();
                ui.add_space(6.0);

                ui.label("Translate along axis");

                let axis_len = if state.picks.len() >= 2 {
                    let a = state.picks[0];
                    let b = state.picks[1];
                    mol.pos[a].distance(mol.pos[b])
                } else {
                    0.0
                };
                if state.picks.len() >= 2 {
                    let pair = (state.picks[0], state.picks[1]);
                    if state.axis_pair != Some(pair) {
                        state.axis_pair = Some(pair);
                        state.axis_len_target = axis_len;
                    }
                } else {
                    state.axis_pair = None;
                    state.axis_len_target = 0.0;
                }

                let max_len = (axis_len * 2.0).max(0.5);
                ui.add_enabled(
                    state.picks.len() >= 2,
                    egui::Slider::new(&mut state.axis_len_target, 0.0..=max_len)
                        .text("Axis length (Å)"),
                );

                ui.horizontal(|ui| {
                    if ui.button("Set Distance").clicked() {
                        if state.picks.len() >= 3 {
                            let a = state.picks[0];
                            let b = state.picks[1];
                            if let Some(side) = state.chosen_side {
                                let current_len = mol.pos[a].distance(mol.pos[b]);
                                let delta = state.axis_len_target - current_len;
                                if delta.abs() > 1e-6 {
                                    state.last_translate_snapshot = Some(mol.pos.clone());
                                    mol.pos = translate_side(
                                        &mol.atoms, &mol.pos, a, b,
                                        side, delta,
                                        state.bond_th_hx, state.bond_th_xx,
                                    );
                                    settings.coords_dirty = true;
                                }
                            }
                        }
                    }
                    if ui.button("Undo").clicked() {
                        if let Some(snapshot) = &mut state.last_translate_snapshot {
                            mol.pos = snapshot.clone();
                            settings.coords_dirty = true;
                        }
                    }
                });

                ui.add_space(6.0);
                ui.separator();
                ui.add_space(6.0);

                ui.label("Bend relative to axis");

                let bend_ready = state.picks.len() >= 3 && state.chosen_side.is_some();
                if bend_ready {
                    let a = state.picks[0];
                    let b = state.picks[1];
                    let p = state.picks[2];
                    let ref_key = (a, b, p);
                    if state.bend_ref != Some(ref_key) {
                        if let Some(side) = state.chosen_side {
                            let center = if side == RotateSide::A { a } else { b };
                            let axis_vec = mol.pos[b] - mol.pos[a];
                            let axis_len = axis_vec.length();
                            let v = mol.pos[p] - mol.pos[center];
                            if axis_len > 1e-6 && v.length() > 1e-6 {
                                let axis = axis_vec / axis_len;
                                let v_n = v.normalize();
                                let mut dot = axis.dot(v_n);
                                dot = dot.clamp(-1.0, 1.0);
                                state.bend_angle_deg = dot.acos().to_degrees();
                            } else {
                                state.bend_angle_deg = 0.0;
                            }
                            state.bend_ref = Some(ref_key);
                        }
                    }
                } else {
                    state.bend_ref = None;
                    state.bend_angle_deg = 0.0;
                }

                ui.add_enabled(
                    bend_ready,
                    egui::Slider::new(&mut state.bend_angle_deg, 0.0..=180.0)
                        .text("Axis angle (°)"),
                );

                ui.horizontal(|ui| {
                    if ui.button("Set Angle").clicked() {
                        if state.picks.len() >= 3 {
                            let a = state.picks[0];
                            let b = state.picks[1];
                            let p = state.picks[2];
                            if let Some(side) = state.chosen_side {
                                state.last_bend_snapshot = Some(mol.pos.clone());
                                mol.pos = bend_side(
                                    &mol.atoms, &mol.pos, a, b,
                                    side, p, state.bend_angle_deg,
                                    state.bond_th_hx, state.bond_th_xx,
                                );
                                settings.coords_dirty = true;
                            }
                        }
                    }
                    if ui.button("Undo").clicked() {
                        if let Some(snapshot) = &mut state.last_bend_snapshot {
                            mol.pos = snapshot.clone();
                            settings.coords_dirty = true;
                        }
                    }
                    if ui.button("Cancel").clicked() {
                        state.picks.clear();
                        state.chosen_side = None;
                        state.angle_deg = 0.0;
                        state.axis_len_target = 0.0;
                        state.axis_pair = None;
                        state.bend_angle_deg = 0.0;
                        state.bend_ref = None;
                        state.last_rotate_snapshot = None;
                        state.last_translate_snapshot = None;
                        state.last_bend_snapshot = None;
                        state.active = false;
                    }
                });

                ui.collapsing("Bond thresholds (Å)", |ui| {
                    ui.add(egui::Slider::new(&mut state.bond_th_hx, 0.6..=2.0).text("H–X (Å)"));
                    ui.add(egui::Slider::new(&mut state.bond_th_xx, 1.0..=3.0).text("X–X (Å)"));
                });
            });
        });
}

// ---------- local reachability helper ----------

/// True if `start` can reach `target` without passing through `exclude`.
fn reaches_without(
    atoms: &[String],
    coords: &[Vec3],
    start: usize,
    target: usize,
    exclude: usize,
    hx: f32,
    xx: f32,
) -> bool {
    use std::collections::{HashMap, VecDeque};
    // Build adjacency with thresholds
    let n = atoms.len();
    let mut adj: HashMap<usize, Vec<usize>> = (0..n).map(|i| (i, Vec::new())).collect();
    for i in 0..n {
        for j in (i + 1)..n {
            let hi = atoms[i] == "H";
            let hj = atoms[j] == "H";
            if hi && hj { continue; }
            let thr = if !hi && !hj { xx } else { hx };
            if coords[i].distance(coords[j]) <= thr {
                adj.get_mut(&i).unwrap().push(j);
                adj.get_mut(&j).unwrap().push(i);
            }
        }
    }

    // BFS from start, skipping `exclude`
    let mut vis = vec![false; n];
    let mut q = VecDeque::new();
    q.push_back(start);
    while let Some(v) = q.pop_front() {
        if v == exclude || vis[v] { continue; }
        if v == target { return true; }
        vis[v] = true;
        if let Some(nei) = adj.get(&v) {
            for &u in nei {
                if !vis[u] && u != exclude {
                    q.push_back(u);
                }
            }
        }
    }
    false
}
