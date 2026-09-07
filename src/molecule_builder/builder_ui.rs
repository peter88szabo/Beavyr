use bevy::prelude::*;
use bevy_egui::egui;

use super::attach;
use super::fragments;
use super::rotator::{
    bend_side_with_bonds, rotate_side_with_bonds, split_sides_from_bonds,
    translate_fragment_containing_atom, RotateSide,
};
use super::zmat2xyz::{self, ZAtom};
use crate::color_schemes::color_for;
use crate::events::{
    AtomPicked, MoleculeChangeReason, MoleculeChanged, ToolKind, ViewportClicked,
};
use crate::molecule::parse_xyz_angstrom;
use crate::molecule::{covalent_radius_angstrom, Molecule};
use crate::settings::MolSettings;

const ANGLE_SINGLE_DEG: f64 = 109.47;
const ANGLE_DOUBLE_DEG: f64 = 120.0;
const ANGLE_TRIPLE_DEG: f64 = 179.95;

const DIHEDRAL_SINGLE_DEG: f64 = 0.0;
const DIHEDRAL_DOUBLE_DEG: f64 = 120.0;
const DIHEDRAL_TRIPLE_DEG: f64 = 180.0;

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
    let (a, b) = if sym_a <= sym_b {
        (sym_a, sym_b)
    } else {
        (sym_b, sym_a)
    };
    match order {
        BondOrder::Triple => {
            if a == "C" && b == "C" {
                return 1.20;
            }
            if a == "C" && b == "N" {
                return 1.16;
            }
            if a == "C" && b == "P" {
                return 1.55;
            }
            if a == "N" && b == "P" {
                return 1.55;
            }
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

    let mut remaining_heavy: Vec<usize> = (0..n).filter(|&i| symbols[i] != "H" && i != 0).collect();

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
            .map(|(d, a)| {
                calc_dihedral_deg(coords[bond_ref_idx], coords[a], coords[d], coords[idx])
            })
            .unwrap_or(0.0);

        // `bond_ref_idx`, `angle_ref_idx` and `dihedral_ref_idx` are indices
        // into the *original* `symbols`/`coords` arrays, found by walking
        // `adj`.  The Z-matrix being built here is in `order` (position)
        // space, so each one must go through `index_of` before it is a valid
        // 1-based reference.  Storing the raw original index instead (as this
        // used to do) is only correct by coincidence for the root atom, whose
        // original index is always 0 -- the same as its position -- and is
        // wrong for almost everything else: `order` reorders heavy atoms by
        // adjacency and appends hydrogens last, so an atom's original index
        // and its position in the built Z-matrix are generally different.
        // The wrong index frequently lands on the row's own position, so the
        // measured angle or dihedral describes a degenerate frame and the
        // rebuilt atom collapses onto one of its neighbours -- which is what
        // made even a plain methyl group come out wrong.
        zmat.push(ZAtom {
            symbol,
            bond_ref: Some(index_of[bond_ref_idx] + 1),
            bond_len,
            angle_ref: angle_ref_idx.map(|v| index_of[v] + 1),
            angle_deg,
            dihedral_ref: dihedral_ref_idx.map(|v| index_of[v] + 1),
            dihedral_deg,
        });
    }

    zmat
}

#[derive(Default, Reflect, GizmoConfigGroup)]
pub struct BuilderGizmos;

#[derive(Default, Reflect, GizmoConfigGroup)]
pub struct BuilderHighlightGizmos;

#[derive(Resource, Default, Clone)]
pub struct EditorRotateState {
    pub active: bool,                    // capture 3D clicks when true
    pub picks: Vec<usize>,               // [A, B, side-indicator (any atom on that side)]
    pub chosen_side: Option<RotateSide>, // inferred from 3rd pick
    pub status: Option<String>,          // selection or topology feedback
    pub angle_deg: f32,                  // degrees
    pub axis_len_target: f32,            // Å, absolute axis length target
    pub axis_pair: Option<(usize, usize)>,
    pub bend_angle_deg: f32, // degrees
    pub bend_ref: Option<(usize, usize, usize)>,
    pub last_rotate_snapshot: Option<Vec<Vec3>>, // undo buffer
    pub last_translate_snapshot: Option<Vec<Vec3>>, // undo buffer
    pub last_bend_snapshot: Option<Vec<Vec3>>,   // undo buffer
    pub bond_th_hx: f32,                         // Å
    pub bond_th_xx: f32,                         // Å
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
    pub frag_scan_score: Option<f32>,
    pub frag_mode: FragmentInsertMode,
    pub last_frag_snapshot: Option<Vec<ZAtom>>,
    pub frag_undo_visible: bool,
    pub edit_rows: Vec<ZAtomEditRow>,
    pub edit_refresh: bool,
    pub edit_preview: Vec<usize>,
    /// The one selected atom, shared by every part of the builder: a
    /// left-click in the 3D view sets it, a click on empty background clears
    /// it, and clicking a row in the Z-matrix table sets it too. "Add Atom"
    /// and "Add Fragment" both act on whatever is selected, so there are no
    /// separate picking modes to arm.
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
            frag_scan_score: None,
            frag_mode: FragmentInsertMode::Connect,
            last_frag_snapshot: None,
            frag_undo_visible: false,
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

/// Splice a fragment onto `zmat`, returning the combined indices it occupies
/// and the clearance achieved.
///
/// The fragment is positioned in Cartesian space and its internal coordinates
/// are then *derived* from that placement.  Filling in only the references and
/// keeping the fragment's own numbers does not work: a self-contained fragment
/// carries placeholder zeros for coordinates it had no reference to measure,
/// and a zero angle places an atom straight back along the bond it arrived by,
/// through the middle of the host.
fn add_fragment_to_zmat(
    zmat: &mut Vec<ZAtom>,
    fragment_xyz: &str,
    connector_refs: Option<(usize, Option<usize>, Option<usize>)>,
    bond_len: f64,
    angle_deg: f64,
    _dihedral_deg: f64,
) -> (Vec<usize>, Option<f32>) {
    let (_n, symbols, coords) = parse_xyz_angstrom(fragment_xyz);
    let frag_coords: Vec<Vec3> = coords
        .into_iter()
        .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
        .collect();
    let frag_atoms = build_fragment_zmat(&symbols, &frag_coords);
    if frag_atoms.is_empty() {
        return (Vec::new(), None);
    }

    // `build_fragment_zmat` reorders atoms, so the local coordinates must be
    // permuted the same way or every derived value refers to the wrong atom.
    let order = fragment_atom_order(&symbols, &frag_coords);
    let local: Vec<Vec3> = order.iter().map(|&i| frag_coords[i]).collect();

    let offset = zmat.len();
    let positions: Vec<usize> = (0..frag_atoms.len()).map(|i| offset + i).collect();
    let mut atoms = frag_atoms;

    let Some((attach, angle_ref, dihedral_ref)) = connector_refs else {
        // Standalone fragment: only its internal references need shifting.
        for atom in atoms.iter_mut().skip(1) {
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
        zmat.extend(atoms);
        return (positions, None);
    };

    attach::assign_fragment_references(
        &mut atoms,
        &positions,
        zmat,
        attach,
        (angle_ref, dihedral_ref),
    );
    zmat.extend(atoms);

    let clearance =
        place_and_derive_fragment(zmat, &positions, attach, &local, bond_len, angle_deg);
    (positions, Some(clearance))
}

/// Search for a non-overlapping rigid placement of a fragment already spliced
/// into `zmat` (references filled by `assign_fragment_references`, connector
/// at `positions[0]`), and bake the result into every fragment row's
/// bond_len/angle_deg/dihedral_deg.
///
/// `local` is the fragment's own Cartesian shape (same order as `positions`,
/// connector first). `attach` is the host atom the fragment bonds to;
/// `bond_len` is the fixed connector-to-host distance; `base_angle_deg` seeds
/// the angle search.
///
/// This is the routine that actually keeps a fragment from overlapping the
/// rest of the molecule.  Filling in *references* (which atoms an angle or
/// dihedral is measured against) is necessary but not sufficient: an atom
/// that is the first or second atom of an isolated fragment has no preceding
/// atoms to measure a real angle or dihedral from, so `build_fragment_zmat`
/// stores a meaningless placeholder (0.0) for those fields.  Pairing that
/// placeholder with a newly filled-in *real* reference is what previously
/// collapsed a fragment atom onto its own reference atom instead of
/// extending outward -- reusing this rigid-body search rather than trusting
/// the fragment's original angle/dihedral numbers is what avoids that.
fn place_and_derive_fragment(
    zmat: &mut [ZAtom],
    positions: &[usize],
    attach: usize,
    local: &[Vec3],
    bond_len: f64,
    base_angle_deg: f64,
) -> f32 {
    let host_coords = zmat2xyz::zmat_to_xyz(zmat);
    // The fragment's own rows are present in `zmat` for both callers (the
    // Replace path splices before calling this; the Connect path now does
    // too, so the two share one code path). Their positions in `host_coords`
    // are whatever the fragment's original, meaningless angle/dihedral
    // placeholders produced -- exactly the values this function exists to
    // replace -- so they must never be treated as real obstacles the search
    // is trying to clear.
    let mut ignore = attach::attachment_partners(zmat, attach);
    ignore.extend_from_slice(positions);

    let connector = positions[0];
    let attach_pos = host_coords.get(attach).copied().unwrap_or(Vec3::ZERO);

    // The connector's own "missing" direction, from its actual bonded local
    // neighbours -- not from every atom in the fragment, which for a
    // multi-atom fragment like -CH=CH2 would wrongly include the terminal
    // =CH2's own hydrogens (two bonds from the connector) as if they were its
    // direct substituents. Computed once: it depends only on the fragment's
    // fixed local shape, not on anything being searched over below.
    let local_symbols: Vec<&str> = positions[1..]
        .iter()
        .map(|&i| zmat[i].symbol.as_str())
        .collect();
    let open = attach::fragment_open_valence(local, &zmat[connector].symbol, &local_symbols);

    // What is the attachment atom's own valence already telling us?  Two
    // existing bonded neighbours at ~120 degrees, or three at ~109.47,
    // determine the completing direction exactly -- there is nothing left to
    // search for the DIRECTION itself, only (optionally) the fragment's spin
    // about it.  This is what makes a conjugated system like butadiene come
    // out flat instead of merely clash-free: clearance alone has no opinion
    // about planarity, since plenty of non-planar orientations have plenty of
    // room.
    let real_neighbors = attach::host_bonded_neighbors(zmat, &host_coords, attach, positions);
    let neighbor_positions: Vec<Vec3> = real_neighbors
        .iter()
        .filter_map(|&i| host_coords.get(i).copied())
        .collect();

    // A 2-neighbour trigonal match is only trustworthy as sp2 when at least
    // one of the existing bonds is short enough to be a real double or
    // aromatic bond. Angle alone cannot tell these apart: an ordinary sp3
    // chain atom caught with only two of its eventual four substituents
    // placed (e.g. the middle carbon of a propane-like chain, before its
    // third and fourth substituents exist) measures ~109.47 degrees, only
    // ~10.5 degrees short of a genuine sp2 atom's ~120 -- comfortably inside
    // any angle tolerance generous enough to admit real-world sp2 variation.
    // Bond length is the actual chemical signal, so it -- not the angle
    // tolerance -- is what has to carry this distinction.
    let bond_ratio = |neighbor: usize, position: Vec3| -> Option<f32> {
        let single_bond_length = covalent_radius_angstrom(&zmat[attach].symbol)
            + covalent_radius_angstrom(&zmat[neighbor].symbol);
        (single_bond_length > 0.0).then(|| attach_pos.distance(position) / single_bond_length)
    };
    const MULTIPLE_BOND_RATIO: f32 = 0.95; // catches C=C (~0.89) and aromatic C-C (~0.93)

    let sp_pattern = attach::ideal_direction_from_neighbors(attach_pos, &neighbor_positions);
    let valence_direction = match (sp_pattern, neighbor_positions.len()) {
        (Some(direction), 2) => {
            let has_multiple_bond = real_neighbors
                .iter()
                .zip(&neighbor_positions)
                .any(|(&n, &pos)| bond_ratio(n, pos).is_some_and(|r| r < MULTIPLE_BOND_RATIO));
            has_multiple_bond.then_some(direction)
        }
        // Three existing neighbours awaiting a fourth is unambiguous: sp2
        // never has room for a fourth substituent, so this can only be sp3
        // regardless of any bond length.
        (Some(direction), _) => Some(direction),
        (None, _) => None,
    }
    .or_else(|| {
        // Only a genuinely single existing neighbour is eligible for the
        // linear/multiple-bond check -- with two or more, checking just the
        // first would silently ignore the rest.
        if real_neighbors.len() != 1 {
            return None;
        }
        let (&single, &only_neighbor) = (&real_neighbors[0], &neighbor_positions[0]);
        let single_bond_length = covalent_radius_angstrom(&zmat[attach].symbol)
            + covalent_radius_angstrom(&zmat[single].symbol);
        attach::linear_direction_for_multiple_bond(attach_pos, only_neighbor, single_bond_length)
    });

    let (clearance, placed) = if let Some(direction) = valence_direction {
        let origin = attach_pos + direction * bond_len as f32;

        // Conjugation: the attachment side is trigonal (exactly two real
        // neighbours) and the fragment's own connector is too (exactly two
        // bonded neighbours in its local, host-less shape -- e.g. vinyl's
        // "=CH-" carbon, locally bonded to its double-bond partner and one
        // H).  Snap the fragment's own plane onto the attachment's plane
        // rather than trusting free rotation to land there by chance.
        let attach_plane_normal = (neighbor_positions.len() == 2).then(|| {
            (neighbor_positions[0] - attach_pos)
                .cross(neighbor_positions[1] - attach_pos)
                .normalize_or_zero()
        });
        let local_others: Vec<(Vec3, &str)> = local
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, &p)| (p, zmat[positions[i]].symbol.as_str()))
            .collect();
        let local_neighbors = attach::geometric_bonded_neighbors(
            local[0],
            &zmat[connector].symbol,
            &local_others,
        );
        let fragment_plane_normal = (local_neighbors.len() == 2).then(|| {
            (local_neighbors[0] - local[0])
                .cross(local_neighbors[1] - local[0])
                .normalize_or_zero()
        });
        let mut best: Option<(f32, f32, Vec<Vec3>)> = None; // (score, clearance, placed)
        let spin_steps = match (attach_plane_normal, fragment_plane_normal) {
            (Some(_), Some(_)) => 72, // 5-degree resolution: precise enough to find s-cis/s-trans
            _ => 24,
        };
        for spin_step in 0..spin_steps {
            let spin = spin_step as f32 * (360.0 / spin_steps as f32);
            let placed = attach::place_fragment(local, open, origin, direction, spin);
            let clearance = attach::clearance_against(&placed, &host_coords, &ignore);

            let score = match (attach_plane_normal, fragment_plane_normal) {
                (Some(host_normal), Some(frag_local_normal)) => {
                    // Carry the fragment's own local plane normal through the
                    // identical rotation `place_fragment` just applied to its
                    // atoms -- not a separately re-derived one, which would
                    // generally disagree with it.
                    let rotation = attach::fragment_rotation(local, open, direction, spin);
                    let rotated_normal = rotation * frag_local_normal;
                    let coplanarity = rotated_normal.normalize_or_zero().dot(host_normal).powi(2);
                    // Planarity dominates; clearance only breaks ties between
                    // the (normally two) orientations that are equally flat --
                    // which is exactly the s-cis vs. s-trans choice, and
                    // favours the roomier s-trans exactly as real dienes do.
                    coplanarity * 1000.0 + clearance
                }
                _ => clearance,
            };

            if best.as_ref().is_none_or(|(best_score, _, _)| score > *best_score) {
                best = Some((score, clearance, placed));
            }
        }
        let (_, clearance, placed) = best.expect("at least one spin is always tried");
        (clearance, placed)
    } else {
        // No unambiguous valence pattern: fall back to the clearance-driven
        // search over the connector's own angle and dihedral relative to the
        // host frame, exactly as before this function knew about valence.
        let frame = zmat[connector]
            .angle_ref
            .zip(zmat[connector].dihedral_ref)
            .and_then(|(b, c)| {
                Some((
                    attach_pos,
                    *host_coords.get(b - 1)?,
                    *host_coords.get(c - 1)?,
                ))
            });

        let mut best: Option<(f32, Vec<Vec3>)> = None;
        let angle_candidates: Vec<f64> = [-20.0, -10.0, 0.0, 10.0, 20.0]
            .iter()
            .map(|delta| (base_angle_deg + delta).clamp(60.0, 180.0))
            .collect();

        match frame {
            Some((a, b, c)) => {
                for &angle in &angle_candidates {
                    for dihedral_step in 0..24 {
                        let dihedral = dihedral_step as f64 * 15.0;
                        let direction = attach::direction_from_internal(a, b, c, angle, dihedral);
                        let origin = a + direction * bond_len as f32;
                        for spin_step in 0..12 {
                            let spin = spin_step as f32 * 30.0;
                            let placed = attach::place_fragment(local, open, origin, direction, spin);
                            let clearance =
                                attach::clearance_against(&placed, &host_coords, &ignore);
                            if best.as_ref().is_none_or(|(score, _)| clearance > *score) {
                                best = Some((clearance, placed));
                            }
                        }
                    }
                }
            }
            None => {
                // Too few host atoms to build a frame; aim away from the centroid.
                let centroid = if host_coords.is_empty() {
                    Vec3::ZERO
                } else {
                    host_coords.iter().copied().sum::<Vec3>() / host_coords.len() as f32
                };
                let direction = (attach_pos - centroid).normalize_or_zero();
                let direction = if direction.length() < 0.5 {
                    Vec3::X
                } else {
                    direction
                };
                let origin = attach_pos + direction * bond_len as f32;
                let placed = attach::place_fragment(local, open, origin, direction, 0.0);
                let clearance = attach::clearance_against(&placed, &host_coords, &ignore);
                best = Some((clearance, placed));
            }
        }
        best.expect("at least one placement is always tried")
    };

    // Measure the internal coordinates off the chosen placement.  `combined`
    // must end up with `placed[k]` at index `positions[k]` regardless of
    // whether those slots already existed in `host_coords` (the Replace path
    // splices the fragment into `zmat` before this runs, so they do, holding
    // the meaningless placeholder positions described above) or not yet
    // existed (the Connect path, where they are new rows) -- so the fragment
    // slots are written after resizing rather than blindly appended.
    let mut combined = host_coords;
    let required_len = positions.iter().copied().max().map_or(0, |m| m + 1);
    if combined.len() < required_len {
        combined.resize(required_len, Vec3::ZERO);
    }
    for (&index, &p) in positions.iter().zip(&placed) {
        combined[index] = p;
    }
    attach::derive_internals(&combined, zmat, positions);
    zmat[connector].bond_len = bond_len;

    clearance
}

fn find_fragment(name: &str) -> Option<&'static fragments::FragmentDef> {
    fragments::FRAGMENTS.iter().find(|f| f.name == name)
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

        let clamp_ref = |v: usize, max_ref: usize| -> usize { v.max(1).min(max_ref) };

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
                Err(err) => {
                    zmat_state.last_error = Some(err);
                    return;
                }
            };
            let bond_len = match parse_f64(&row.bond_len, "bond length", i) {
                Ok(v) => v,
                Err(err) => {
                    zmat_state.last_error = Some(err);
                    return;
                }
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
                Err(err) => {
                    zmat_state.last_error = Some(err);
                    return;
                }
            };
            let bond_len = match parse_f64(&row.bond_len, "bond length", i) {
                Ok(v) => v,
                Err(err) => {
                    zmat_state.last_error = Some(err);
                    return;
                }
            };
            let angle_ref = match parse_usize(&row.angle_ref, "angle", i) {
                Ok(v) => v,
                Err(err) => {
                    zmat_state.last_error = Some(err);
                    return;
                }
            };
            let angle_deg = match parse_f64(&row.angle_deg, "angle", i) {
                Ok(v) => v,
                Err(err) => {
                    zmat_state.last_error = Some(err);
                    return;
                }
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
                Err(err) => {
                    zmat_state.last_error = Some(err);
                    return;
                }
            };
            let bond_len = match parse_f64(&row.bond_len, "bond length", i) {
                Ok(v) => v,
                Err(err) => {
                    zmat_state.last_error = Some(err);
                    return;
                }
            };
            let angle_ref = match parse_usize(&row.angle_ref, "angle", i) {
                Ok(v) => v,
                Err(err) => {
                    zmat_state.last_error = Some(err);
                    return;
                }
            };
            let angle_deg = match parse_f64(&row.angle_deg, "angle", i) {
                Ok(v) => v,
                Err(err) => {
                    zmat_state.last_error = Some(err);
                    return;
                }
            };
            let dihedral_ref = match parse_usize(&row.dihedral_ref, "dihedral", i) {
                Ok(v) => v,
                Err(err) => {
                    zmat_state.last_error = Some(err);
                    return;
                }
            };
            let dihedral_deg = match parse_f64(&row.dihedral_deg, "dihedral", i) {
                Ok(v) => v,
                Err(err) => {
                    zmat_state.last_error = Some(err);
                    return;
                }
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
            status: None,
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
    {
        let (cfg, _) = cfg_store.config_mut::<BuilderGizmos>();
        cfg.enabled = true;
        cfg.line.width = 16.0;
        cfg.line.perspective = true;
    }
    {
        let (cfg, _) = cfg_store.config_mut::<BuilderHighlightGizmos>();
        cfg.enabled = true;
        cfg.line.width = 2.0;
        cfg.line.perspective = true;
        cfg.line.style = GizmoLineStyle::Dashed {
            gap_scale: 2.0,
            line_scale: 1.25,
        };
    }
}

/// Put the structure back as it was before the last fragment was added.
///
/// Both fragment rows call this, so the very first fragment dropped into an
/// empty editor is as undoable as one spliced onto an existing atom.
fn undo_last_fragment(
    zmat_state: &mut ZMatrixBuilderState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
) {
    if let Some(snapshot) = zmat_state.last_frag_snapshot.take() {
        zmat_state.zmat = snapshot;
        mol.atoms = zmat_state
            .zmat
            .iter()
            .map(|atom| atom.symbol.clone())
            .collect();
        mol.pos = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
        mol.recompute_bonds(2.0, 3.0);
        settings.geometry_dirty = true;
        settings.bond_topology_dirty = true;
        zmat_state.edit_refresh = true;
    }
    zmat_state.frag_undo_visible = false;
    // The atom the fragment was hung on may not exist any more (an undone
    // Replace restores the atom that was swapped out), so start clean.
    set_selected_atom(zmat_state, &mol.atoms, None);
    zmat_state.last_error = None;
}

/// Whether a click should clear the Z-matrix table's selection.
///
/// Clicking the panel away from the grid deselects, which is what makes the
/// table feel like a list. A click in the 3D view must NOT: that click is
/// itself a selection (`handle_viewport_click`), and clearing here would undo
/// it in the same frame -- which is why picking an atom stopped working as
/// soon as the Z-matrix section was open, since the grid rect only exists
/// while that section is drawn.
fn click_clears_table_selection(
    clicked: bool,
    pointer_over_egui: bool,
    pointer_pos: Option<egui::Pos2>,
    grid: egui::Rect,
) -> bool {
    clicked
        && pointer_over_egui
        && pointer_pos.is_some_and(|p| !grid.contains(p))
}

/// Select an atom, or clear the selection with `None`.
///
/// The symbol is cached alongside the index because the remove-atom popup and
/// the edit rows want to name the atom without re-deriving it, and because an
/// out-of-range index is rejected here rather than at every use site.
pub fn set_selected_atom(
    zmat_state: &mut ZMatrixBuilderState,
    atoms: &[String],
    index: Option<usize>,
) {
    match index.filter(|&i| i < atoms.len()) {
        Some(i) => {
            zmat_state.selected_index = Some(i);
            zmat_state.selected_symbol = Some(atoms[i].clone());
        }
        None => {
            zmat_state.selected_index = None;
            zmat_state.selected_symbol = None;
        }
    }
    zmat_state.edit_refresh = true;
}

/// Splice the fragment named in `zmat_state.frag_name` onto the selected atom.
/// The placement search fills in every reference and direction that a
/// hand-picked angle/dihedral atom used to supply, via `host_reference_chain`
/// and the valence-aware direction rules.
fn commit_fragment_connect(
    zmat_state: &mut ZMatrixBuilderState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
) {
    let Some(bond_ref) = zmat_state.selected_index else {
        zmat_state.last_error = Some("Click an atom to connect to.".to_string());
        return;
    };
    if !zmat_matches_molecule(zmat_state, mol) {
        return;
    }
    let Some(frag) = find_fragment(&zmat_state.frag_name) else {
        zmat_state.last_error = Some("Fragment not found.".to_string());
        return;
    };

    zmat_state.last_frag_snapshot = Some(zmat_state.zmat.clone());
    zmat_state.frag_undo_visible = true;
    let angle_ref = None;
    let dihedral_ref = None;

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
    // Placement and the clearance search now happen inside
    // the splice, in Cartesian space, so the reported score
    // describes the geometry that was actually adopted.
    let (_fragment_indices, clearance) = add_fragment_to_zmat(
        &mut zmat_state.zmat,
        frag.xyz,
        Some((bond_ref, angle_ref, dihedral_ref)),
        bond_len,
        frag_angle_deg,
        frag_dihedral_deg,
    );
    zmat_state.frag_scan_score = clearance;

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

    // The selection survives: Connect only appends rows, so the selected
    // index still names the same atom and the user can hang a second
    // substituent off it without re-clicking.
    zmat_state.edit_refresh = true;
}

/// Replace the selected atom with the fragment named in
/// `zmat_state.frag_name`.
///
/// Substitution is what it is chemically: take the atom away, then hang the
/// fragment off the atom that was actually holding it. Both steps work on the
/// Cartesian structure and the Z-matrix is rebuilt from the result, rather
/// than the old row being edited in place.
///
/// Editing in place is what made this unusable. A Z-matrix row's `bond_ref` is
/// only "the row this one was measured against" -- for a structure that came
/// from `xyz_to_zmat` that is simply the previous line of the file. Treating it
/// as "the atom this one is bonded to" attached the fragment to whatever
/// happened to be listed before the atom being replaced, which for an
/// aromatic hydrogen is usually *another hydrogen*. Perceiving the real bond
/// from the coordinates and rebuilding afterwards removes the whole class of
/// problem, including having to keep every later row's references intact.
fn commit_fragment_replace(
    zmat_state: &mut ZMatrixBuilderState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
    _editor_state: &EditorRotateState,
) {
    let Some(replace_idx) = zmat_state.selected_index else {
        zmat_state.last_error = Some("Click an atom to replace.".to_string());
        return;
    };
    if !zmat_matches_molecule(zmat_state, mol) {
        return;
    }
    if find_fragment(&zmat_state.frag_name).is_none() {
        zmat_state.last_error = Some("Fragment not found.".to_string());
        return;
    }
    if replace_idx >= mol.atoms.len() {
        zmat_state.last_error = Some("Invalid replacement index.".to_string());
        return;
    }

    // The atom that actually holds the one being replaced, by covalent-radius
    // distance -- the same test that draws the bonds on screen. Nearest first,
    // so a bridging atom hands the fragment to its closest partner.
    let symbols: Vec<&str> = mol.atoms.iter().map(|s| s.as_str()).collect();
    let Some(&host) = attach::bonded_indices(&symbols, &mol.pos, replace_idx, &[]).first() else {
        zmat_state.last_error =
            Some("That atom is not bonded to anything, so there is nothing to replace it on."
                .to_string());
        return;
    };

    // Undo has to restore the structure from before the atom was taken away,
    // so this snapshot must survive the Connect step's own bookkeeping.
    let snapshot = zmat_state.zmat.clone();

    let mut atoms = mol.atoms.clone();
    let mut pos = mol.pos.clone();
    atoms.remove(replace_idx);
    pos.remove(replace_idx);
    let host = if host > replace_idx { host - 1 } else { host };

    mol.atoms = atoms;
    mol.pos = pos;
    mol.recompute_bonds(2.0, 3.0);
    // Rebuilt from the geometry that is left, so no inherited row can carry a
    // reference to the atom that just went away.
    zmat_state.zmat = zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos);

    // From here it is an ordinary Connect onto the host, which places the
    // fragment with the valence-aware rules: with the hydrogen gone, an
    // aromatic carbon has two ring neighbours at 120 degrees, so the sp2 rule
    // points the new group exactly where the hydrogen pointed.
    zmat_state.selected_index = Some(host);
    commit_fragment_connect(zmat_state, mol, settings);
    if zmat_state.last_error.is_some() {
        // Put the structure back rather than leave the atom deleted.
        zmat_state.zmat = snapshot;
        mol.atoms = zmat_state
            .zmat
            .iter()
            .map(|atom| atom.symbol.clone())
            .collect();
        mol.pos = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
        mol.recompute_bonds(2.0, 3.0);
        settings.geometry_dirty = true;
        settings.bond_topology_dirty = true;
        set_selected_atom(zmat_state, &mol.atoms, None);
        return;
    }

    zmat_state.last_frag_snapshot = Some(snapshot);
    zmat_state.frag_undo_visible = true;
    // The atom the user picked no longer exists, so nothing is selected.
    set_selected_atom(zmat_state, &mol.atoms, None);
}

/// Add a new, single atom bonded to the selected atom, using the same
/// valence-aware placement as
/// fragment attachment: `host_reference_chain` fills in the angle/dihedral
/// references, and `ideal_direction_from_neighbors` /
/// `linear_direction_for_multiple_bond` fix the exact direction when the
/// attachment atom's existing bonded neighbours unambiguously imply one
/// (sp2 trigonal, sp3 tetrahedral, or sp linear). Only when no such pattern
/// applies does the new atom fall back to the angle/dihedral the user typed
/// in directly.
fn commit_add_atom(zmat_state: &mut ZMatrixBuilderState, mol: &mut Molecule, settings: &mut MolSettings) {
    let Some(bond_ref) = zmat_state.selected_index else {
        zmat_state.last_error = Some("Click an atom to bond to.".to_string());
        return;
    };
    if !zmat_matches_molecule(zmat_state, mol) {
        return;
    }
    let symbol = zmat_state.new_symbol.trim().to_string();
    if symbol.is_empty() {
        zmat_state.last_error = Some("Atom symbol cannot be empty.".to_string());
        return;
    }

    let r_new = covalent_radius_angstrom(&symbol);
    let r_ref = covalent_radius_angstrom(&mol.atoms[bond_ref]);
    let bond_len = bond_length_for(
        &symbol,
        &mol.atoms[bond_ref],
        zmat_state.new_bond_order,
        r_new,
        r_ref,
    );

    let (chain_angle, chain_dihedral) = attach::host_reference_chain(&zmat_state.zmat, bond_ref);
    let mut entry = ZAtom {
        symbol,
        bond_ref: Some(bond_ref + 1),
        bond_len,
        angle_ref: chain_angle.map(|v| v + 1),
        angle_deg: zmat_state.new_angle_deg,
        dihedral_ref: chain_dihedral.map(|v| v + 1),
        dihedral_deg: zmat_state.new_dihedral_deg,
    };

    // Refine the angle/dihedral against the attachment atom's real bonded
    // neighbours, exactly as a fragment's connector would be placed. A plain
    // added atom has no fragment shape to search a spin over, so this only
    // ever fixes the direction itself, never a rotation around it.
    let host_coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
    if let Some(&attach_pos) = host_coords.get(bond_ref) {
        let real_neighbors =
            attach::host_bonded_neighbors(&zmat_state.zmat, &host_coords, bond_ref, &[]);
        let neighbor_positions: Vec<Vec3> = real_neighbors
            .iter()
            .filter_map(|&i| host_coords.get(i).copied())
            .collect();
        let bond_ratio = |neighbor: usize, position: Vec3| -> Option<f32> {
            let single_bond_length = covalent_radius_angstrom(&zmat_state.zmat[bond_ref].symbol)
                + covalent_radius_angstrom(&zmat_state.zmat[neighbor].symbol);
            (single_bond_length > 0.0).then(|| attach_pos.distance(position) / single_bond_length)
        };
        const MULTIPLE_BOND_RATIO: f32 = 0.95;

        let sp_pattern = attach::ideal_direction_from_neighbors(attach_pos, &neighbor_positions);
        let valence_direction = match (sp_pattern, neighbor_positions.len()) {
            (Some(direction), 2) => real_neighbors
                .iter()
                .zip(&neighbor_positions)
                .any(|(&n, &pos)| bond_ratio(n, pos).is_some_and(|r| r < MULTIPLE_BOND_RATIO))
                .then_some(direction),
            (Some(direction), _) => Some(direction),
            (None, _) => None,
        }
        .or_else(|| {
            if real_neighbors.len() != 1 {
                return None;
            }
            let single_bond_length = covalent_radius_angstrom(&zmat_state.zmat[bond_ref].symbol)
                + covalent_radius_angstrom(&zmat_state.zmat[real_neighbors[0]].symbol);
            attach::linear_direction_for_multiple_bond(
                attach_pos,
                neighbor_positions[0],
                single_bond_length,
            )
        });

        if let (Some(direction), Some(b), Some(c)) =
            (valence_direction, entry.angle_ref, entry.dihedral_ref)
        {
            if let (Some(&pb), Some(&pc)) = (host_coords.get(b - 1), host_coords.get(c - 1)) {
                let (angle, dihedral) =
                    attach::internal_from_direction(attach_pos, pb, pc, direction);
                entry.angle_deg = angle;
                entry.dihedral_deg = dihedral;
            }
        }
    }

    zmat_state.zmat.push(entry);

    let coords = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);
    mol.atoms = zmat_state.zmat.iter().map(|atom| atom.symbol.clone()).collect();
    mol.pos = coords;
    mol.recompute_bonds(2.0, 3.0);
    settings.geometry_dirty = true;
    settings.bond_topology_dirty = true;
    zmat_state.last_error = None;

    // Adding only appends a row, so the selection still names the same atom
    // and can take another substituent straight away.
    zmat_state.edit_refresh = true;
}

/// Turn a left-click in the 3D view into the builder's selection: an atom
/// selects it, empty background clears it. The picker only reports clicks that
/// did not turn into a camera orbit, so rotating the view leaves the selection
/// alone.
pub fn handle_viewport_click(
    mut ev: MessageReader<ViewportClicked>,
    mol: Option<Res<Molecule>>,
    mut zmat_state: ResMut<ZMatrixBuilderState>,
) {
    let Some(mol) = mol else {
        return;
    };
    for ViewportClicked { hit } in ev.read().copied() {
        set_selected_atom(&mut zmat_state, &mol.atoms, hit);
    }
}

/// Bring the builder in line with a newly loaded structure: rebuild the
/// Z-matrix from it and drop the selection.
///
/// Both matter. A stale selection index usually stays in range and would
/// quietly highlight -- and be built onto -- some unrelated atom of the new
/// molecule. And the Z-matrix used to be rebuilt only by the "Sync from
/// Molecule" button, so every builder action after a file load was operating
/// on a description of whatever was in the editor beforehand.
pub fn sync_builder_on_structure_load(
    mut evr: MessageReader<MoleculeChanged>,
    mol: Option<Res<Molecule>>,
    mut zmat_state: ResMut<ZMatrixBuilderState>,
) {
    let structure_replaced = evr
        .read()
        .any(|ev| matches!(ev.reason, MoleculeChangeReason::ParseXyz { .. }));
    if !structure_replaced {
        return;
    }
    let atoms: &[String] = match &mol {
        Some(mol) => {
            zmat_state.zmat = zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos);
            &mol.atoms
        }
        None => &[],
    };
    let atoms = atoms.to_vec();
    set_selected_atom(&mut zmat_state, &atoms, None);
    zmat_state.last_frag_snapshot = None;
    zmat_state.frag_undo_visible = false;
    zmat_state.last_remove_snapshot = None;
    zmat_state.redo_remove_visible = false;
    zmat_state.last_error = None;
}

/// Refuse to build when the Z-matrix does not describe the molecule on screen.
///
/// Every commit rewrites `mol` from the Z-matrix, so acting on a stale one
/// replaces the displayed structure with a different molecule entirely. The
/// builder normally keeps the two in step (`sync_builder_on_structure_load`),
/// and this is the backstop for anything that slips through.
fn zmat_matches_molecule(zmat_state: &mut ZMatrixBuilderState, mol: &Molecule) -> bool {
    if zmat_state.zmat.len() == mol.atoms.len() && !zmat_state.zmat.is_empty() {
        return true;
    }
    zmat_state.last_error = Some(
        "Z-matrix is out of sync with the displayed molecule. Press \"Sync from Molecule\"."
            .to_string(),
    );
    false
}

pub fn handle_builder_atom_picked(
    mut ev: MessageReader<AtomPicked>,
    mol: Option<ResMut<Molecule>>,
    mut state: ResMut<EditorRotateState>,
) {
    // The Z-matrix builder's own selection is handled by
    // `handle_viewport_click`; what is left here is the bond-rotate tool,
    // which needs three atoms picked in a set order and so stays explicitly
    // armed.
    let Some(mol) = mol else {
        return;
    };
    if !state.active {
        return;
    }

    for AtomPicked {
        index: hit_idx,
        tool,
    } in ev.read().copied()
    {
        if tool != ToolKind::Builder {
            continue;
        }

        match state.picks.len() {
            0 => {
                state.picks.push(hit_idx);
                state.status = None;
            }
            1 => {
                if hit_idx == state.picks[0] {
                    continue;
                }
                let a = state.picks[0];
                if split_sides_from_bonds(mol.atoms.len(), &mol.bonds, a, hit_idx).is_some() {
                    state.picks.push(hit_idx);
                    state.status = None;
                } else {
                    state.status = Some(
                        "Select a displayed covalent bond that is not part of a ring.".to_string(),
                    );
                }
            }
            2 => {
                let (a, b) = (state.picks[0], state.picks[1]);
                let Some((side_a, side_b)) =
                    split_sides_from_bonds(mol.atoms.len(), &mol.bonds, a, b)
                else {
                    state.picks.truncate(1);
                    state.chosen_side = None;
                    state.status = Some(
                        "That bond no longer separates two fragments. Pick another bond."
                            .to_string(),
                    );
                    continue;
                };

                if side_a.contains(&hit_idx) {
                    state.picks.push(hit_idx);
                    state.chosen_side = Some(RotateSide::A);
                    state.status = None;
                } else if side_b.contains(&hit_idx) {
                    state.picks.push(hit_idx);
                    state.chosen_side = Some(RotateSide::B);
                    state.status = None;
                } else {
                    state.status = Some(
                        "Pick an atom on one of the two fragments attached to this bond."
                            .to_string(),
                    );
                }
            }
            _ => {
                // A completed selection starts over with a new A endpoint.
                state.picks.clear();
                state.chosen_side = None;
                state.status = None;
                state.picks.push(hit_idx);
            }
        }
    }
}

/// Draw pending highlights and axis line using the same color/radius style as measurements.
pub fn draw_builder_highlights(
    mut gizmos: Gizmos<BuilderGizmos>,
    mut highlights: Gizmos<BuilderHighlightGizmos>,
    state: Res<EditorRotateState>,
    zmat_state: Res<ZMatrixBuilderState>,
    mol: Option<Res<Molecule>>,
    settings: Res<MolSettings>,
) {
    let Some(mol) = mol else {
        return;
    };
    if mol.pos.is_empty() {
        return;
    }

    let magenta = Color::srgba(1.0, 0.0, 0.8, 1.0);
    let cyan = Color::srgba(0.0, 0.8, 1.0, 1.0);
    let green = Color::srgba(0.2, 1.0, 0.2, 1.0);

    let draw_hl = |gizmos: &mut Gizmos<BuilderHighlightGizmos>, idx: usize, color: Color| {
        if idx < mol.pos.len() {
            let base_r = covalent_radius_angstrom(&mol.atoms[idx]) * settings.atom_scale;
            let r = (base_r * 1.25).max(0.15);
            gizmos.sphere(mol.pos[idx], r, color);
        }
    };

    if state.active {
        if let Some(&i0) = state.picks.get(0) {
            draw_hl(&mut highlights, i0, magenta);
        }
        if let Some(&i1) = state.picks.get(1) {
            draw_hl(&mut highlights, i1, cyan);
            if let Some(&i0) = state.picks.get(0) {
                gizmos.line(mol.pos[i0], mol.pos[i1], Color::WHITE);
            }
        }
        if let Some(&i2) = state.picks.get(2) {
            draw_hl(&mut highlights, i2, green);
        }
    }

    // The one selected atom, however it was selected -- clicked in the 3D
    // view or picked by row in the Z-matrix table. Both write the same field,
    // so the highlight follows either.
    if let Some(selected) = zmat_state.selected_index.filter(|&i| i < mol.pos.len()) {
        draw_hl(&mut highlights, selected, magenta);
    }

    if !zmat_state.edit_preview.is_empty() {
        let orange = Color::srgba(1.0, 0.6, 0.1, 1.0);
        let colors = [magenta, cyan, green, orange];
        for (i, idx) in zmat_state.edit_preview.iter().copied().enumerate() {
            let color = colors.get(i).copied().unwrap_or(orange);
            draw_hl(&mut highlights, idx, color);
        }
    }
}

/// Left-side, resizable panel with a collapsible “Molecule Editor” section.
/// NOTE: We now do a clean jump: no immediate atom transform sync; instead we
/// update `mol.pos` + `mol.bonds` and flag geometry/topology dirty so the scene
/// rebuilds atoms & bonds together next frame (no mismatched frame).
/// The molecule editor's contents, independent of the container it is drawn
/// in: the classic layout docks this in a permanent left panel, while the
/// tabbed layout puts it in a floating window opened from a small round
/// button, so it no longer occupies a whole screen edge.
pub fn builder_ui_contents(
    ui: &mut egui::Ui,
    state: &mut EditorRotateState,
    mut zmat_state: &mut ZMatrixBuilderState,
    mut mol: &mut Molecule,
    mut settings: &mut MolSettings,
) {
    let ctx = ui.ctx().clone();
            if zmat_state.original_atoms.is_none() && !mol.atoms.is_empty() {
                zmat_state.original_atoms = Some(mol.atoms.clone());
                zmat_state.original_pos = Some(mol.pos.clone());
            }
            ui.add_space(8.0);
            ui.collapsing("Molecule Builder", |ui| {
                ui.horizontal(|ui| {
                    if ui.button("Sync from Molecule").clicked() {
                        zmat_state.zmat = zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos);
                        set_selected_atom(&mut zmat_state, &mol.atoms, None);
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
                        set_selected_atom(&mut zmat_state, &mol.atoms, None);
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
                            set_selected_atom(&mut zmat_state, &mol.atoms, None);
                            zmat_state.last_error = None;
                        } else {
                            zmat_state.last_error =
                                Some("Original structure not captured.".to_string());
                        }
                    }
                });

                ui.add_space(6.0);
                // One line saying what "Add Atom" and "Add Fragment" will act
                // on. It replaces the old per-section picking hints: there is
                // only one selection now, however it was made.
                match (zmat_state.selected_index, &zmat_state.selected_symbol) {
                    (Some(idx), Some(symbol)) => {
                        ui.label(format!("Selected: {symbol}{}", idx + 1));
                    }
                    _ => {
                        ui.weak("Click an atom to select it; click the background to clear.");
                    }
                }
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
                                                .unwrap_or_else(|| {
                                                    color_for(row.symbol.trim(), settings.scheme)
                                                })
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
                                            hl_bg = Some(egui::Color32::from_rgba_premultiplied(
                                                200, 200, 200, 60,
                                            ));
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
                                            let bond_ref_resp = ui.add_sized(
                                                [col_idx[2], 0.0],
                                                egui::TextEdit::singleline(&mut row.bond_ref),
                                            );
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

                                            let bond_len_resp = ui.add_sized(
                                                [col_idx[3], 0.0],
                                                egui::TextEdit::singleline(&mut row.bond_len),
                                            );
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
                                            let angle_ref_resp = ui.add_sized(
                                                [col_idx[4], 0.0],
                                                egui::TextEdit::singleline(&mut row.angle_ref),
                                            );
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

                                            let angle_deg_resp = ui.add_sized(
                                                [col_idx[5], 0.0],
                                                egui::TextEdit::singleline(&mut row.angle_deg),
                                            );
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
                                            let dihedral_ref_resp = ui.add_sized(
                                                [col_idx[6], 0.0],
                                                egui::TextEdit::singleline(&mut row.dihedral_ref),
                                            );
                                            if dihedral_ref_resp.changed() {
                                                row.dihedral_ref =
                                                    row.dihedral_ref.trim().to_string();
                                            }
                                            if dihedral_ref_resp.has_focus() {
                                                preview_indices = Some(build_preview_indices(
                                                    i,
                                                    row,
                                                    PreviewKind::Dihedral,
                                                    row_count,
                                                ));
                                            }

                                            let dihedral_deg_resp = ui.add_sized(
                                                [col_idx[7], 0.0],
                                                egui::TextEdit::singleline(&mut row.dihedral_deg),
                                            );
                                            if dihedral_deg_resp.changed() {
                                                row.dihedral_deg =
                                                    row.dihedral_deg.trim().to_string();
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
                        let over_egui = ctx.is_pointer_over_egui();
                        let (clicked, pointer_pos) =
                            ctx.input(|i| (i.pointer.any_click(), i.pointer.latest_pos()));
                        if !zmat_state.remove_popup_open
                            && !suppress_clear
                            && click_clears_table_selection(clicked, over_egui, pointer_pos, rect)
                        {
                            zmat_state.selected_index = None;
                            zmat_state.selected_symbol = None;
                            zmat_state.edit_preview.clear();
                            zmat_state.redo_remove_visible = false;
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
                    if zmat_state.zmat.is_empty() {
                        // Nothing to bond to yet, so there is no atom to pick:
                        // the very first atom is just dropped in directly.
                        if ui.button("Add Atom").clicked() {
                            let symbol = zmat_state.new_symbol.trim().to_string();
                            if symbol.is_empty() {
                                zmat_state.last_error =
                                    Some("Atom symbol cannot be empty.".to_string());
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
                        }
                    } else {
                        // No picking mode to arm: the atom to bond to is
                        // whichever one is selected, and selecting is just a
                        // left-click in the 3D view (or a click on a row in
                        // the Z-matrix table).
                        let can_commit = zmat_state.selected_index.is_some();
                        if ui
                            .add_enabled(can_commit, egui::Button::new("Add Atom"))
                            .clicked()
                        {
                            commit_add_atom(&mut zmat_state, &mut mol, &mut settings);
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
                                for &bo in
                                    &[BondOrder::Single, BondOrder::Double, BondOrder::Triple]
                                {
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
                    if zmat_state.zmat.is_empty() {
                        // Nothing exists yet to attach to, so there is no atom
                        // to select: the first fragment in an empty editor is
                        // just dropped in directly.
                        if ui.button("Add Fragment").clicked() {
                            if let Some(frag) = find_fragment(&zmat_state.frag_name) {
                                zmat_state.last_frag_snapshot = Some(zmat_state.zmat.clone());
                                zmat_state.frag_undo_visible = true;
                                let frag_angle_deg = zmat_state.frag_angle_deg;
                                let frag_dihedral_deg = zmat_state.frag_dihedral_deg;
                                // First fragment in an empty editor: nothing to
                                // collide with, so no relaxation.
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
                        }
                        if zmat_state.frag_undo_visible
                            && ui.button("Undo Fragment").clicked()
                        {
                            undo_last_fragment(zmat_state, mol, settings);
                        }
                        return;
                    }

                    // Nothing to arm: the fragment goes onto whichever atom
                    // is selected, and selecting is a left-click in the 3D
                    // view (or a click on a Z-matrix row). Connect needs only
                    // that one atom -- the placement search derives the
                    // angle/dihedral references and the valence-aware
                    // direction itself.
                    let can_commit = zmat_state.selected_index.is_some();
                    if ui
                        .add_enabled(can_commit, egui::Button::new("Add Fragment"))
                        .clicked()
                    {
                        match zmat_state.frag_mode {
                            FragmentInsertMode::Connect => {
                                commit_fragment_connect(zmat_state, mol, settings)
                            }
                            FragmentInsertMode::Replace => {
                                commit_fragment_replace(zmat_state, mol, settings, state)
                            }
                        }
                    }

                    // Right next to "Add Fragment", so the way out of a
                    // placement you did not want is where you just looked.
                    // It stays until it is used or until a structural edit
                    // makes the snapshot stale -- clicking elsewhere no
                    // longer takes it away.
                    if zmat_state.frag_undo_visible
                        && ui.button("Undo Fragment").clicked()
                    {
                        undo_last_fragment(zmat_state, mol, settings);
                    }
                });

                // Why a press did nothing, said where the press happened. The
                // panel's other error line sits above the Z-matrix table,
                // which in a scrolling window is off-screen by the time these
                // controls are in view -- so a failure looked like silence.
                if let Some(err) = &zmat_state.last_error {
                    ui.colored_label(egui::Color32::LIGHT_RED, err);
                }
                ui.weak(format!(
                    "Z-matrix {} rows \u{2022} molecule {} atoms \u{2022} selected {}",
                    zmat_state.zmat.len(),
                    mol.atoms.len(),
                    match (zmat_state.selected_index, &zmat_state.selected_symbol) {
                        (Some(i), Some(sym)) => format!("{sym}{}", i + 1),
                        _ => "nothing".to_string(),
                    }
                ));

                if let Some(score) = zmat_state.frag_scan_score {
                    ui.weak(format!("Fragment scan min distance: {:.3} Å", score));
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
                        egui::Frame::popup(ui.style()).show(ui, |ui| {
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
                            if ui
                                .add_enabled(can_remove, egui::Button::new("Remove"))
                                .clicked()
                            {
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
                    if ui
                        .button(if state.active {
                            "Deactivate"
                        } else {
                            "Activate"
                        })
                        .clicked()
                    {
                        state.active = !state.active;
                        if !state.active {
                            state.picks.clear();
                            state.chosen_side = None;
                            state.status = None;
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
                        Some(RotateSide::A) => {
                            ui.colored_label(egui::Color32::LIGHT_GREEN, "Chosen side: A")
                        }
                        Some(RotateSide::B) => {
                            ui.colored_label(egui::Color32::LIGHT_GREEN, "Chosen side: B")
                        }
                        None => ui.weak("Chosen side: —"),
                    };
                }

                if let Some(status) = &state.status {
                    ui.colored_label(egui::Color32::LIGHT_RED, status);
                }

                ui.add_space(6.0);
                ui.add(egui::Slider::new(&mut state.angle_deg, -180.0..=180.0).text("Angle (°)"));
                ui.horizontal(|ui| {
                    if ui.button("Rotate").clicked() && state.picks.len() >= 3 {
                        let a = state.picks[0];
                        let b = state.picks[1];
                        if let Some(side) = state.chosen_side {
                            if let Some(new_pos) = rotate_side_with_bonds(
                                &mol.pos,
                                &mol.bonds,
                                a,
                                b,
                                side,
                                state.angle_deg,
                            ) {
                                state.last_rotate_snapshot = Some(mol.pos.clone());
                                mol.pos = new_pos;
                                settings.coords_dirty = true;
                                state.status = None;
                            } else {
                                state.status = Some(
                                    "The selected bond no longer separates two fragments."
                                        .to_string(),
                                );
                            }
                        }
                    }
                    if ui.button("Undo").clicked() {
                        if let Some(snapshot) = &state.last_rotate_snapshot {
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
                    if ui.button("Set Distance").clicked() && state.picks.len() >= 3 {
                        let a = state.picks[0];
                        let b = state.picks[1];
                        let selected_atom = state.picks[2];
                        let delta = state.axis_len_target - mol.pos[a].distance(mol.pos[b]);
                        if delta.abs() > 1.0e-6 {
                            if let Some(new_pos) = translate_fragment_containing_atom(
                                &mol.pos,
                                &mol.bonds,
                                a,
                                b,
                                selected_atom,
                                delta,
                            ) {
                                state.last_translate_snapshot = Some(mol.pos.clone());
                                mol.pos = new_pos;
                                settings.coords_dirty = true;
                                state.status = None;
                            } else {
                                state.status = Some(
                                    "The third-picked atom is not on a movable bond fragment."
                                        .to_string(),
                                );
                            }
                        }
                    }
                    if ui.button("Undo").clicked() {
                        if let Some(snapshot) = &state.last_translate_snapshot {
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
                            if axis_len > 1.0e-6 && v.length() > 1.0e-6 {
                                state.bend_angle_deg = (axis_vec / axis_len)
                                    .dot(v.normalize())
                                    .clamp(-1.0, 1.0)
                                    .acos()
                                    .to_degrees();
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
                    if ui.button("Set Angle").clicked() && state.picks.len() >= 3 {
                        let a = state.picks[0];
                        let b = state.picks[1];
                        let p = state.picks[2];
                        if let Some(side) = state.chosen_side {
                            if let Some(new_pos) = bend_side_with_bonds(
                                &mol.pos,
                                &mol.bonds,
                                a,
                                b,
                                side,
                                p,
                                state.bend_angle_deg,
                            ) {
                                state.last_bend_snapshot = Some(mol.pos.clone());
                                mol.pos = new_pos;
                                settings.coords_dirty = true;
                                state.status = None;
                            } else {
                                state.status = Some(
                                    "The selected bond no longer separates two fragments."
                                        .to_string(),
                                );
                            }
                        }
                    }
                    if ui.button("Undo").clicked() {
                        if let Some(snapshot) = &state.last_bend_snapshot {
                            mol.pos = snapshot.clone();
                            settings.coords_dirty = true;
                        }
                    }
                    if ui.button("Cancel").clicked() {
                        state.picks.clear();
                        state.chosen_side = None;
                        state.status = None;
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
}

/// The classic docked-left-panel presentation of the molecule editor.
pub fn builder_ui_panel(
    ui: &mut egui::Ui,
    state: &mut EditorRotateState,
    zmat_state: &mut ZMatrixBuilderState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
) {
    egui::Panel::left("molecule_editor_panel")
        .resizable(true)
        .show(ui, |ui| {
            builder_ui_contents(ui, state, zmat_state, mol, settings);
        });
}


#[cfg(test)]
mod fragment_placement_tests {
    use super::*;
    use crate::molecule_builder::attach;
    use crate::molecule_builder::zmat2xyz::xyz_to_zmat;

    /// Three carbons with a free coordination site on the middle one, so a
    /// fragment attached there has somewhere correct to go.
    const HOST: &str = "\
C   0.000   0.000   0.000
C   1.540   0.000   0.000
C   2.050   1.450   0.000
";

    fn host_zmat() -> Vec<ZAtom> {
        let (_n, symbols, coords) = parse_xyz_angstrom(HOST);
        let positions: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        xyz_to_zmat(&symbols, &positions)
    }

    fn splice_methyl(attach_index: usize) -> (Vec<ZAtom>, Vec<usize>) {
        let mut zmat = host_zmat();
        let fragment = find_fragment("-CH3").expect("-CH3 is a built-in fragment");
        let (indices, _clearance) = add_fragment_to_zmat(
            &mut zmat,
            fragment.xyz,
            Some((attach_index, None, None)),
            1.54,
            109.471,
            180.0,
        );
        (zmat, indices)
    }

    /// Splicing must leave a fully referenced Z-matrix.  A single `None` past
    /// the third row silently resolves against a distant atom.
    #[test]
    fn a_spliced_fragment_is_fully_referenced() {
        let (zmat, indices) = splice_methyl(1);
        assert!(!indices.is_empty(), "nothing was spliced");

        for (i, atom) in zmat.iter().enumerate().skip(3) {
            assert!(atom.bond_ref.is_some(), "row {i} has no bond reference");
            assert!(atom.angle_ref.is_some(), "row {i} has no angle reference");
            assert!(
                atom.dihedral_ref.is_some(),
                "row {i} has no dihedral reference"
            );
        }
    }

    /// The reported bug: an attached fragment overlapped the rest of the
    /// molecule.  After placement no fragment atom may sit on top of a host
    /// atom it is not bonded to.
    #[test]
    fn an_attached_fragment_does_not_overlap_the_host() {
        let attach_index = 1;
        let (zmat, indices) = splice_methyl(attach_index);

        let coords = zmat2xyz::zmat_to_xyz(&zmat);
        let ignore = attach::attachment_partners(&zmat, attach_index);
        let clearance = attach::min_clearance(&coords, &indices, &ignore);

        assert!(
            clearance > 1.2,
            "fragment sits {clearance:.2} A from the host; anything under a bond \
             length means it was planted inside the molecule"
        );
    }

    /// Demonstrates the defect the fix removes.  Pointing the fragment's angle
    /// and dihedral references at atom 1 of the whole molecule -- which is what
    /// an unfilled `None` used to resolve to -- rebuilds its geometry in an
    /// unrelated frame and measurably worsens the placement.
    #[test]
    fn resolving_references_against_atom_one_is_worse() {
        let attach_index = 1;

        let (fixed, indices) = splice_methyl(attach_index);
        let ignore = attach::attachment_partners(&fixed, attach_index);
        let good = attach::min_clearance(&zmat2xyz::zmat_to_xyz(&fixed), &indices, &ignore);

        let (mut broken, indices) = splice_methyl(attach_index);
        for &i in indices.iter().skip(1) {
            broken[i].angle_ref = Some(1);
            broken[i].dihedral_ref = Some(1);
        }
        let bad = attach::min_clearance(&zmat2xyz::zmat_to_xyz(&broken), &indices, &ignore);

        assert!(
            good > bad,
            "inherited references gave {good:.2} A, atom-1 references {bad:.2} A; \
             the fix should be an improvement"
        );
    }

    /// The splice reports the clearance it actually achieved, so the panel's
    /// score is not describing a different geometry from the one on screen.
    #[test]
    fn the_reported_clearance_matches_the_placement() {
        let attach_index = 1;
        let mut zmat = host_zmat();
        let fragment = find_fragment("-CH3").unwrap();
        let (indices, reported) = add_fragment_to_zmat(
            &mut zmat,
            fragment.xyz,
            Some((attach_index, None, None)),
            1.54,
            109.471,
            180.0,
        );
        let reported = reported.expect("a clearance should have been reported");

        let ignore = attach::attachment_partners(&zmat, attach_index);
        let measured = attach::min_clearance(&zmat2xyz::zmat_to_xyz(&zmat), &indices, &ignore);
        assert!(
            (reported - measured).abs() < 0.05,
            "reported {reported:.3} A but the Z-matrix rebuilds to {measured:.3} A"
        );
    }

    /// The connector must keep the bond length it was given; only the angle and
    /// dihedral are free to move.
    #[test]
    fn relaxation_leaves_the_bond_length_alone() {
        let attach_index = 1;
        let (zmat, indices) = splice_methyl(attach_index);
        let connector = indices[0];
        let bond_len = zmat[connector].bond_len;
        assert!(
            (bond_len - 1.54).abs() < 1.0e-9,
            "bond length was not honoured"
        );

        let coords = zmat2xyz::zmat_to_xyz(&zmat);
        let measured = coords[connector].distance(coords[attach_index]) as f64;
        assert!(
            (measured - bond_len).abs() < 1.0e-3,
            "connector ended up {measured:.3} A away, expected {bond_len:.3}"
        );
    }
}

#[cfg(test)]
mod fragment_geometry_tests {
    use super::*;
    use crate::molecule_builder::zmat2xyz::xyz_to_zmat;

    const HOST: &str = "\
C   0.000   0.000   0.000
C   1.540   0.000   0.000
C   2.050   1.450   0.000
";

    /// A spliced methyl must actually look like a methyl: four bonds near
    /// 1.07 A and six H-C-H angles near 109.5 degrees.
    #[test]
    fn a_spliced_methyl_has_tetrahedral_geometry() {
        let (_n, symbols, coords) = parse_xyz_angstrom(HOST);
        let positions: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let mut zmat = xyz_to_zmat(&symbols, &positions);

        let attach = 1usize;
        let fragment = find_fragment("-CH3").unwrap();
        let (indices, _clearance) = add_fragment_to_zmat(
            &mut zmat,
            fragment.xyz,
            Some((attach, None, None)),
            1.54,
            109.471,
            180.0,
        );

        let out = zmat2xyz::zmat_to_xyz(&zmat);
        let carbon = indices[0];
        let hydrogens = &indices[1..];

        for &h in hydrogens {
            let d = out[carbon].distance(out[h]);
            assert!(
                (d - 1.089).abs() < 0.02,
                "C-H came out {d:.3} A, expected 1.089"
            );
        }

        for (a, &h1) in hydrogens.iter().enumerate() {
            for &h2 in &hydrogens[a + 1..] {
                let angle = calc_angle_deg(out[h1], out[carbon], out[h2]);
                assert!(
                    (angle - 109.5).abs() < 5.0,
                    "H-C-H came out {angle:.1} deg, expected 109.5"
                );
            }
        }

        // And the methyl must point away from the host, not back into it.
        let to_host = (out[attach] - out[carbon]).normalize();
        for &h in hydrogens {
            let to_h = (out[h] - out[carbon]).normalize();
            let angle = to_host.dot(to_h).clamp(-1.0, 1.0).acos().to_degrees();
            assert!(
                angle > 90.0,
                "an H sits {angle:.1} deg from the host bond, i.e. folded back \
                 towards the molecule"
            );
        }
    }

    /// Every built-in fragment must produce a self-consistent Z-matrix: each
    /// reference must point at a strictly earlier row.  The original bug
    /// stored an *original atom index* where a *Z-matrix position* was
    /// expected; for several of these fragments that landed exactly on the
    /// row's own position, which silently degenerates the frame instead of
    /// failing to parse.  Checking every built-in fragment rather than just
    /// the one that was reported catches the same defect in the others
    /// (`-Phenyl`, `-CycloHexane`, ... reorder heavy atoms exactly the way
    /// `-CH3` did) before a user finds it fragment by fragment.
    #[test]
    fn every_builtin_fragment_has_a_self_consistent_zmat() {
        for fragment in fragments::FRAGMENTS {
            let (_n, symbols, coords) = parse_xyz_angstrom(fragment.xyz);
            let frag_coords: Vec<Vec3> = coords
                .into_iter()
                .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
                .collect();
            let zmat = build_fragment_zmat(&symbols, &frag_coords);

            for (pos, atom) in zmat.iter().enumerate() {
                for (label, reference) in [
                    ("bond", atom.bond_ref),
                    ("angle", atom.angle_ref),
                    ("dihedral", atom.dihedral_ref),
                ] {
                    if let Some(r) = reference {
                        assert!(
                            r >= 1 && r - 1 < pos,
                            "{}: row {pos} {label}_ref = {r} does not point to an                              earlier row",
                            fragment.name
                        );
                    }
                }
            }
        }
    }

    /// Splicing a fragment must not collapse two of its atoms onto (nearly)
    /// the same point.  This is the shape the original defect actually took:
    /// references were valid array indices, just the wrong ones, so nothing
    /// panicked -- two hydrogens simply landed on top of each other.
    #[test]
    fn every_builtin_fragment_splices_without_atoms_colliding() {
        let host = "\
C   0.000   0.000   0.000
C   1.540   0.000   0.000
C   2.050   1.450   0.000
";
        for fragment in fragments::FRAGMENTS {
            let (_n, symbols, coords) = parse_xyz_angstrom(host);
            let positions: Vec<Vec3> = coords
                .into_iter()
                .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
                .collect();
            let mut zmat = xyz_to_zmat(&symbols, &positions);

            let attach = 1usize;
            let (indices, _clearance) = add_fragment_to_zmat(
                &mut zmat,
                fragment.xyz,
                Some((attach, None, None)),
                1.5,
                109.471,
                180.0,
            );
            if indices.is_empty() {
                continue;
            }

            let out = zmat2xyz::zmat_to_xyz(&zmat);
            for a in 0..indices.len() {
                for b in (a + 1)..indices.len() {
                    let d = out[indices[a]].distance(out[indices[b]]);
                    assert!(
                        d > 0.3,
                        "{}: fragment atoms {a} and {b} are only {d:.3} A apart",
                        fragment.name
                    );
                }
            }
        }
    }

    /// -NH2's geometry was corrected against Molden's canonical values
    /// (N-H = 1.010 A, H-N-H = 120 deg). Measuring it directly from the
    /// stored fragment text, rather than after a splice, isolates the
    /// fixture itself from any splice-time relaxation.
    #[test]
    fn nh2_fragment_matches_moldens_corrected_geometry() {
        let frag = find_fragment("-NH2").expect("-NH2 is a built-in fragment");
        let (_n, symbols, coords) = parse_xyz_angstrom(frag.xyz);
        let positions: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        assert_eq!(symbols, vec!["N", "H", "H"]);

        let d1 = positions[0].distance(positions[1]);
        let d2 = positions[0].distance(positions[2]);
        assert!((d1 - 1.010).abs() < 1e-4, "N-H(1) = {d1}, expected 1.010");
        assert!((d2 - 1.010).abs() < 1e-4, "N-H(2) = {d2}, expected 1.010");

        let angle = calc_angle_deg(positions[1], positions[0], positions[2]);
        assert!(
            (angle - 120.0).abs() < 1e-2,
            "H-N-H = {angle}, expected 120"
        );
    }

    /// -CH3's C-H bond length was rescaled to Molden's value (1.089 A) by
    /// uniformly scaling the existing tetrahedral vectors, specifically to
    /// avoid re-deriving the shape from Molden's raw Z-matrix (which, as
    /// `every_builtin_fragment_has_a_self_consistent_zmat` and the splice
    /// tests already guard, is not a safe transplant here -- see the doc
    /// comment on `FRAG_CH3`). This test is the check that the scaling
    /// preserved perfect tetrahedral symmetry.
    #[test]
    fn ch3_fragment_is_tetrahedral_at_moldens_bond_length() {
        let frag = find_fragment("-CH3").expect("-CH3 is a built-in fragment");
        let (_n, symbols, coords) = parse_xyz_angstrom(frag.xyz);
        let positions: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        assert_eq!(symbols, vec!["C", "H", "H", "H"]);

        for &h in &positions[1..] {
            let d = positions[0].distance(h);
            assert!((d - 1.089).abs() < 1e-4, "C-H = {d}, expected 1.089");
        }
        for a in 1..4 {
            for b in (a + 1)..4 {
                let angle = calc_angle_deg(positions[a], positions[0], positions[b]);
                assert!(
                    (angle - 109.471).abs() < 1e-2,
                    "H{a}-C-H{b} = {angle}, expected 109.471"
                );
            }
        }
    }

    /// The four fragments newly added from Molden's library (Cl, Br, I: single
    /// terminal atoms; SH: two atoms, no angle or dihedral involved) must
    /// parse correctly and splice onto a host without panicking, at a
    /// reasonable bond length.
    #[test]
    fn newly_added_molden_fragments_parse_and_splice() {
        let host = "\
C   0.000   0.000   0.000
C   1.540   0.000   0.000
C   2.050   1.450   0.000
";
        for (name, expected_atoms) in [
            ("-Cl", vec!["Cl"]),
            ("-Br", vec!["Br"]),
            ("-I", vec!["I"]),
            ("-SH", vec!["S", "H"]),
        ] {
            let frag = find_fragment(name).unwrap_or_else(|| panic!("{name} should be built in"));
            let (_n, symbols, _coords) = parse_xyz_angstrom(frag.xyz);
            assert_eq!(symbols, expected_atoms, "{name} atom list");

            let (_n2, symbols2, coords2) = parse_xyz_angstrom(host);
            let positions: Vec<Vec3> = coords2
                .into_iter()
                .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
                .collect();
            let mut zmat = xyz_to_zmat(&symbols2, &positions);

            let attach = 1usize;
            let (indices, _clearance) = add_fragment_to_zmat(
                &mut zmat,
                frag.xyz,
                Some((attach, None, None)),
                1.5,
                109.471,
                180.0,
            );
            assert!(!indices.is_empty(), "{name} produced no spliced atoms");

            let out = zmat2xyz::zmat_to_xyz(&zmat);
            let connector = indices[0];
            let d = out[attach].distance(out[connector]);
            assert!(
                (d - 1.5).abs() < 1e-2,
                "{name}: host-connector bond came out {d:.3} A, expected 1.5"
            );
        }
    }
}

#[cfg(test)]
mod toluene_replace_diagnostic {
    use super::*;
    use crate::molecule_builder::zmat2xyz::xyz_to_zmat;

    fn min_clearance_all(coords: &[Vec3], group: &[usize]) -> f32 {
        let mut in_group = vec![false; coords.len()];
        for &i in group {
            in_group[i] = true;
        }
        let mut smallest = f32::MAX;
        for &i in group {
            for (j, &p) in coords.iter().enumerate() {
                if in_group[j] {
                    continue;
                }
                smallest = smallest.min(coords[i].distance(p));
            }
        }
        smallest
    }

    #[test]
    fn replacing_a_methyl_hydrogen_on_toluene_with_phenyl_does_not_overlap() {
        // Build toluene: benzene ring + -CH3 via the (already-fixed) Connect path.
        let benzene = "\
C   1.400   0.000   0.000
C   0.700   1.212   0.000
C  -0.700   1.212   0.000
C  -1.400   0.000   0.000
C  -0.700  -1.212   0.000
C   0.700  -1.212   0.000
";
        let (_n, symbols, coords) = parse_xyz_angstrom(benzene);
        let positions: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let mut zmat = xyz_to_zmat(&symbols, &positions);

        let ch3 = find_fragment("-CH3").unwrap();
        let (methyl_indices, _clearance) = add_fragment_to_zmat(
            &mut zmat,
            ch3.xyz,
            Some((0, None, None)),
            1.51,
            109.471,
            180.0,
        );
        let methyl_hydrogens = &methyl_indices[1..];

        // Sanity: this is the already-fixed Connect path; toluene must be clean.
        let coords_before = zmat2xyz::zmat_to_xyz(&zmat);
        let toluene_clearance = min_clearance_all(&coords_before, methyl_indices.as_slice());
        assert!(
            toluene_clearance > 0.9,
            "toluene itself already overlaps ({toluene_clearance:.3} A); \
             the bug must be reproduced starting from a clean structure"
        );

        // Now replace one methyl hydrogen with -Phenyl via the Replace path.
        let replace_idx = methyl_hydrogens[0];
        let mut zmat_state = ZMatrixBuilderState::default();
        zmat_state.zmat = zmat;
        zmat_state.frag_name = "-Phenyl".to_string();
        zmat_state.selected_index = Some(replace_idx);

        let mut mol = Molecule {
            atoms: zmat_state.zmat.iter().map(|a| a.symbol.clone()).collect(),
            pos: zmat2xyz::zmat_to_xyz(&zmat_state.zmat),
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        mol.recompute_bonds(1.2, 2.5);
        let mut settings = MolSettings::default();
        let editor_state = EditorRotateState::default();

        commit_fragment_replace(&mut zmat_state, &mut mol, &mut settings, &editor_state);
        assert!(
            zmat_state.last_error.is_none(),
            "replace failed: {:?}",
            zmat_state.last_error
        );

        let out = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);

        // Replace takes the hydrogen away and appends the fragment, so the new
        // phenyl ring is the tail of the Z-matrix: everything from the
        // pre-replace count minus that one removed atom.
        let phenyl_indices: Vec<usize> =
            (coords_before.len() - 1..zmat_state.zmat.len()).collect();

        let clearance = min_clearance_all(&out, &phenyl_indices);
        assert!(
            clearance > 1.2,
            "new phenyl group sits {clearance:.3} A from the rest of toluene; \
             expected it to be pushed away along the new bond, not overlapping"
        );
    }
}

#[cfg(test)]
mod conjugation_tests {
    use super::*;
    use crate::molecule_builder::zmat2xyz::xyz_to_zmat;

    /// Connect a vinyl group to a lone host carbon, then replace the vinyl's
    /// own remaining hydrogen (on its "=CH-" connector, the sp2 carbon) with a
    /// second vinyl group.  This is the reported scenario: two -CH=CH2
    /// fragments joined into a butadiene-like conjugated system.  The
    /// attachment carbon ends up with exactly two real neighbours -- the host
    /// and its double-bond partner -- at exactly 120 degrees (by construction
    /// of the first splice), which is exactly the case the sp2 valence rule
    /// exists to recognise.
    #[test]
    fn joining_two_vinyl_groups_produces_a_planar_conjugated_system() {
        let host = "C   0.000   0.000   0.000\n";
        let (_n, symbols, coords) = parse_xyz_angstrom(host);
        let positions: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let mut zmat = xyz_to_zmat(&symbols, &positions);

        let vinyl = find_fragment("-CH=CH2").expect("-CH=CH2 is a built-in fragment");
        let (first_indices, _clearance) =
            add_fragment_to_zmat(&mut zmat, vinyl.xyz, Some((0, None, None)), 1.51, 120.0, 180.0);

        // Identify C1 (the connector, bonded to the host) and C2 (its double
        // bond partner) by role rather than assuming a fixed index order.
        let c1 = first_indices[0];
        let c2 = *first_indices[1..]
            .iter()
            .find(|&&i| zmat[i].symbol == "C")
            .expect("the vinyl fragment has a second carbon");
        // The hydrogen bonded directly to C1 (not to C2) is the vacancy the
        // user replaces to form the conjugated chain.
        let c1_hydrogen = *first_indices
            .iter()
            .find(|&&i| zmat[i].symbol == "H" && zmat[i].bond_ref == Some(c1 + 1))
            .expect("C1 has exactly one hydrogen of its own");

        let mut zmat_state = ZMatrixBuilderState::default();
        let atoms_before_replace = zmat.len();
        zmat_state.zmat = zmat;
        zmat_state.frag_name = "-CH=CH2".to_string();
        zmat_state.selected_index = Some(c1_hydrogen);

        let mut mol = Molecule {
            atoms: zmat_state.zmat.iter().map(|a| a.symbol.clone()).collect(),
            pos: zmat2xyz::zmat_to_xyz(&zmat_state.zmat),
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        mol.recompute_bonds(1.2, 2.5);
        let mut settings = MolSettings::default();
        let editor_state = EditorRotateState::default();

        commit_fragment_replace(&mut zmat_state, &mut mol, &mut settings, &editor_state);

        // Replace deletes the selected atom, so every later index shifts; the
        // selection must not be left pointing at whatever moved into its slot.
        assert_eq!(zmat_state.selected_index, None);
        assert!(
            zmat_state.last_error.is_none(),
            "replace failed: {:?}",
            zmat_state.last_error
        );

        let out = zmat2xyz::zmat_to_xyz(&zmat_state.zmat);

        // The replaced hydrogen is gone and the second vinyl is appended, so
        // its connector (C1\') is the first of the appended rows and its own
        // double-bond partner (C2\') is the other new carbon.
        let c1_prime = atoms_before_replace - 1;
        assert_eq!(zmat_state.zmat[c1_prime].symbol, "C");
        let c2_prime = (c1_prime + 1..zmat_state.zmat.len())
            .find(|&i| {
                zmat_state.zmat[i].symbol == "C" && zmat_state.zmat[i].bond_ref == Some(c1_prime + 1)
            })
            .expect("the second vinyl's own double-bond carbon");

        // The whole conjugated system (C2=C1-C1'=C2') must be planar: the
        // dihedral around the central C1-C1' bond must land near 0 (s-cis) or
        // 180 (s-trans) degrees, not at some arbitrary clash-avoiding angle.
        let dihedral = calc_dihedral_deg(out[c2], out[c1], out[c1_prime], out[c2_prime]);
        let distance_from_planar = dihedral.abs().min((180.0 - dihedral.abs()).abs());
        assert!(
            distance_from_planar < 5.0,
            "C=C-C=C dihedral was {dihedral:.1} degrees; expected close to 0 or 180 \
             (planar conjugation), not an arbitrary clash-avoiding twist"
        );

        // And the two carbon-carbon-carbon angles at the shared bond must both
        // be close to the ideal 120 degrees, not merely "clash-free".
        let angle_at_c1 = calc_angle_deg(out[c2], out[c1], out[c1_prime]);
        let angle_at_c1_prime = calc_angle_deg(out[c1], out[c1_prime], out[c2_prime]);
        assert!(
            (angle_at_c1 - 120.0).abs() < 2.0,
            "C2-C1-C1' angle was {angle_at_c1:.1} degrees, expected ~120"
        );
        assert!(
            (angle_at_c1_prime - 120.0).abs() < 2.0,
            "C1-C1'-C2' angle was {angle_at_c1_prime:.1} degrees, expected ~120"
        );
    }

    /// A plain sp3 chain carbon that merely happens to have two of its four
    /// eventual substituents in place (an ordinary alkyl chain, not a double
    /// bond) must NOT be misclassified as sp2 -- its ~109.47-degree angle
    /// sits only ~10.5 degrees short of a genuine sp2 angle, well inside any
    /// angle tolerance generous enough to admit real sp2 variation, so this
    /// has to be decided by bond length (a real double/aromatic bond) rather
    /// than angle alone.
    #[test]
    fn an_ordinary_sp3_chain_atom_is_not_treated_as_sp2() {
        let host = "\
C   0.000   0.000   0.000
C   1.540   0.000   0.000
C   2.053   1.450   0.000
";
        let (_n, symbols, coords) = parse_xyz_angstrom(host);
        let positions: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let mut zmat = xyz_to_zmat(&symbols, &positions);

        let methyl = find_fragment("-CH3").unwrap();
        let attach = 1usize; // the middle carbon: two ordinary single bonds already
        let real_neighbors =
            attach::host_bonded_neighbors(&zmat, &zmat2xyz::zmat_to_xyz(&zmat), attach, &[]);
        assert_eq!(real_neighbors.len(), 2, "fixture sanity: two existing neighbours");

        let (indices, _clearance) =
            add_fragment_to_zmat(&mut zmat, methyl.xyz, Some((attach, None, None)), 1.54, 109.471, 180.0);
        let out = zmat2xyz::zmat_to_xyz(&zmat);
        let connector = indices[0];

        // If this were wrongly treated as sp2, the new C-C-C angle would land
        // near 120 degrees, coplanar with the two existing bonds.  Correctly
        // treated as an ordinary sp3 chain atom, the placement search instead
        // uses the fragment's own tetrahedral default (~109.47) with only a
        // small clash-avoidance nudge.
        let angle = calc_angle_deg(out[0], out[attach], out[connector]);
        assert!(
            (angle - 120.0).abs() > 5.0,
            "middle chain carbon's new substituent landed at {angle:.1} degrees, \
             suspiciously close to a wrongly-applied sp2 120-degree rule"
        );
    }
}

#[cfg(test)]
mod add_atom_tests {
    use super::*;
    use crate::molecule_builder::zmat2xyz::xyz_to_zmat;

    fn host() -> (Molecule, ZMatrixBuilderState) {
        let xyz = "\
C   0.000   0.000   0.000
C   1.540   0.000   0.000
";
        let (_n, symbols, coords) = parse_xyz_angstrom(xyz);
        let positions: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let zmat = xyz_to_zmat(&symbols, &positions);
        let mol = Molecule {
            atoms: zmat.iter().map(|a| a.symbol.clone()).collect(),
            pos: zmat2xyz::zmat_to_xyz(&zmat),
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        let mut state = ZMatrixBuilderState::default();
        state.zmat = zmat;
        state.new_symbol = "H".to_string();
        (mol, state)
    }

    /// The selected atom is the whole input: one click in the 3D view (or one
    /// Z-matrix row) is all "Add Atom" needs.
    #[test]
    fn add_atom_commits_from_the_selected_atom() {
        let (mut mol, mut state) = host();
        state.selected_index = Some(0);

        commit_add_atom(&mut state, &mut mol, &mut MolSettings::default());

        assert!(state.last_error.is_none(), "commit failed: {:?}", state.last_error);
        assert_eq!(state.zmat.len(), 3, "the new atom should have been appended");
        assert_eq!(state.zmat[2].symbol, "H");
        assert_eq!(state.zmat[2].bond_ref, Some(1));
        // The selection survives: adding only appends, so index 0 still names
        // the same atom and can take a second substituent right away.
        assert_eq!(state.selected_index, Some(0));
    }

    /// With nothing selected, committing must fail loudly rather than guess.
    #[test]
    fn add_atom_without_a_selection_is_an_error() {
        let (mut mol, mut state) = host();
        commit_add_atom(&mut state, &mut mol, &mut MolSettings::default());
        assert!(state.last_error.is_some());
        assert_eq!(state.zmat.len(), 2, "nothing should have been added");
    }

    /// Bonding to an atom that already has two real neighbours at a
    /// double-bond-confirmed 120 degrees should place the new atom using the
    /// same exact sp2 rule fragments use, not merely the user's typed angle.
    #[test]
    fn add_atom_uses_the_valence_direction_when_available() {
        // A vinyl-like host: C1=C2, plus one H already on C1 at ~120 degrees.
        let xyz = "\
C   0.000   0.000   0.000
C   1.335   0.000   0.000
H  -0.544   0.943   0.000
";
        let (_n, symbols, coords) = parse_xyz_angstrom(xyz);
        let positions: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let zmat = xyz_to_zmat(&symbols, &positions);
        let mut mol = Molecule {
            atoms: zmat.iter().map(|a| a.symbol.clone()).collect(),
            pos: zmat2xyz::zmat_to_xyz(&zmat),
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        mol.recompute_bonds(1.2, 2.5);
        let mut state = ZMatrixBuilderState::default();
        state.zmat = zmat;
        state.new_symbol = "H".to_string();
        // Deliberately a bad angle: if the valence rule fires, this is ignored.
        state.new_angle_deg = 60.0;
        // C1 (index 0) already has two real neighbours (C2, H) at ~120 degrees.
        state.selected_index = Some(0);

        commit_add_atom(&mut state, &mut mol, &mut MolSettings::default());
        assert!(state.last_error.is_none());

        let out = zmat2xyz::zmat_to_xyz(&state.zmat);
        let new_atom = state.zmat.len() - 1;
        let angle_to_c2 = calc_angle_deg(out[1], out[0], out[new_atom]);
        assert!(
            (angle_to_c2 - 120.0).abs() < 2.0,
            "new H landed at {angle_to_c2:.1} degrees from C2; expected ~120 \
             from the sp2 valence rule, not the requested 60"
        );
    }
}

#[cfg(test)]
mod selection_tests {
    use super::*;

    fn state_with(atoms: &[&str]) -> (Vec<String>, ZMatrixBuilderState) {
        (
            atoms.iter().map(|s| s.to_string()).collect(),
            ZMatrixBuilderState::default(),
        )
    }

    #[test]
    fn selecting_an_atom_records_its_index_and_symbol() {
        let (atoms, mut state) = state_with(&["C", "C", "H"]);
        set_selected_atom(&mut state, &atoms, Some(2));
        assert_eq!(state.selected_index, Some(2));
        assert_eq!(state.selected_symbol.as_deref(), Some("H"));
    }

    #[test]
    fn a_background_click_clears_the_selection() {
        let (atoms, mut state) = state_with(&["C", "H"]);
        set_selected_atom(&mut state, &atoms, Some(1));
        set_selected_atom(&mut state, &atoms, None);
        assert_eq!(state.selected_index, None);
        assert_eq!(state.selected_symbol, None);
    }

    /// A stale index -- from a structure that has since shrunk -- must clear
    /// the selection rather than be stored and indexed with later on.
    #[test]
    fn an_out_of_range_index_clears_instead_of_being_stored() {
        let (atoms, mut state) = state_with(&["C", "H"]);
        set_selected_atom(&mut state, &atoms, Some(1));
        set_selected_atom(&mut state, &atoms, Some(7));
        assert_eq!(state.selected_index, None);
        assert_eq!(state.selected_symbol, None);
    }

    /// Loading a different structure must not leave a selection behind: the
    /// old index would still be in range for the new molecule.
    #[test]
    fn loading_a_structure_clears_the_selection() {
        let mut app = App::new();
        app.add_message::<MoleculeChanged>()
            .init_resource::<ZMatrixBuilderState>()
            .add_systems(Update, sync_builder_on_structure_load);
        {
            let mut state = app.world_mut().resource_mut::<ZMatrixBuilderState>();
            state.selected_index = Some(3);
            state.selected_symbol = Some("H".to_string());
        }

        app.world_mut()
            .write_message(MoleculeChanged::parse_xyz(true));
        app.update();

        let state = app.world().resource::<ZMatrixBuilderState>();
        assert_eq!(state.selected_index, None);
        assert_eq!(state.selected_symbol, None);
    }

    /// An unrelated molecule change -- a fragment rotation, say -- leaves the
    /// selection alone; the atom the user picked is still that atom.
    #[test]
    fn a_mere_geometry_change_keeps_the_selection() {
        let mut app = App::new();
        app.add_message::<MoleculeChanged>()
            .init_resource::<ZMatrixBuilderState>()
            .add_systems(Update, sync_builder_on_structure_load);
        {
            let mut state = app.world_mut().resource_mut::<ZMatrixBuilderState>();
            state.selected_index = Some(3);
            state.selected_symbol = Some("H".to_string());
        }

        app.world_mut().write_message(MoleculeChanged::set_pos());
        app.update();

        assert_eq!(
            app.world().resource::<ZMatrixBuilderState>().selected_index,
            Some(3)
        );
    }

    /// Committing an Add Atom leaves the selection in place, so a second
    /// substituent can go on the same centre without re-clicking it.
    #[test]
    fn adding_an_atom_keeps_the_selection() {
        let xyz = "\
C   0.000   0.000   0.000
C   1.540   0.000   0.000
";
        let (_n, symbols, coords) = parse_xyz_angstrom(xyz);
        let positions: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let zmat = crate::molecule_builder::zmat2xyz::xyz_to_zmat(&symbols, &positions);
        let mut mol = Molecule {
            atoms: zmat.iter().map(|a| a.symbol.clone()).collect(),
            pos: zmat2xyz::zmat_to_xyz(&zmat),
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        let mut state = ZMatrixBuilderState::default();
        state.zmat = zmat;
        state.new_symbol = "H".to_string();
        state.selected_index = Some(0);

        commit_add_atom(&mut state, &mut mol, &mut MolSettings::default());

        assert!(state.last_error.is_none(), "{:?}", state.last_error);
        assert_eq!(state.selected_index, Some(0));
    }
}

#[cfg(test)]
mod table_selection_tests {
    use super::*;

    fn grid() -> egui::Rect {
        egui::Rect::from_min_max(egui::pos2(100.0, 100.0), egui::pos2(400.0, 300.0))
    }

    /// The regression that broke picking: with the Z-matrix section open, a
    /// click in the 3D view sits outside the grid rect, and the old rule
    /// cleared the selection the viewport click had just made.
    #[test]
    fn a_click_in_the_3d_view_does_not_clear_the_selection() {
        assert!(!click_clears_table_selection(
            true,
            false, // not over egui: this is the 3D view
            Some(egui::pos2(700.0, 500.0)),
            grid(),
        ));
    }

    /// Clicking the panel away from the grid still deselects, which is what
    /// makes the table behave like a list.
    #[test]
    fn a_click_on_the_panel_away_from_the_grid_clears() {
        assert!(click_clears_table_selection(
            true,
            true,
            Some(egui::pos2(420.0, 500.0)),
            grid(),
        ));
    }

    #[test]
    fn a_click_inside_the_grid_is_left_to_the_rows_themselves() {
        assert!(!click_clears_table_selection(
            true,
            true,
            Some(egui::pos2(200.0, 200.0)),
            grid(),
        ));
    }

    #[test]
    fn merely_hovering_clears_nothing() {
        assert!(!click_clears_table_selection(
            false,
            true,
            Some(egui::pos2(420.0, 500.0)),
            grid(),
        ));
    }

    #[test]
    fn a_click_with_no_known_position_clears_nothing() {
        assert!(!click_clears_table_selection(true, true, None, grid()));
    }
}

#[cfg(test)]
mod fragment_undo_tests {
    use super::*;

    fn host() -> (Molecule, ZMatrixBuilderState) {
        let xyz = "\
C   0.000   0.000   0.000
H   1.090   0.000   0.000
";
        let (_n, symbols, coords) = parse_xyz_angstrom(xyz);
        let positions: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let zmat = crate::molecule_builder::zmat2xyz::xyz_to_zmat(&symbols, &positions);
        let mut mol = Molecule {
            atoms: zmat.iter().map(|a| a.symbol.clone()).collect(),
            pos: zmat2xyz::zmat_to_xyz(&zmat),
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        mol.recompute_bonds(2.0, 3.0);
        let mut state = ZMatrixBuilderState::default();
        state.zmat = zmat;
        state.frag_name = "-CH3".to_string();
        (mol, state)
    }

    /// Adding a fragment must leave an undo available, and taking it must put
    /// the structure back exactly as it was.
    #[test]
    fn undo_restores_the_structure_from_before_the_fragment() {
        let (mut mol, mut state) = host();
        let before = state.zmat.len();
        let atoms_before = mol.atoms.clone();
        state.selected_index = Some(0);

        commit_fragment_connect(&mut state, &mut mol, &mut MolSettings::default());
        assert!(state.last_error.is_none(), "{:?}", state.last_error);
        assert!(state.zmat.len() > before, "the fragment should have landed");
        assert!(
            state.frag_undo_visible,
            "an undo must be offered right after adding a fragment"
        );

        undo_last_fragment(&mut state, &mut mol, &mut MolSettings::default());

        assert_eq!(state.zmat.len(), before);
        assert_eq!(mol.atoms, atoms_before);
        assert!(!state.frag_undo_visible, "the undo is spent");
        assert_eq!(state.selected_index, None);
    }

    /// Undo must also work for a Replace, where the swapped-out atom has to
    /// come back.
    #[test]
    fn undo_brings_back_an_atom_that_replace_swapped_out() {
        let (mut mol, mut state) = host();
        let atoms_before = mol.atoms.clone();
        let zmat_before = state.zmat.clone();
        state.selected_index = Some(1); // the hydrogen
        state.frag_mode = FragmentInsertMode::Replace;

        commit_fragment_replace(
            &mut state,
            &mut mol,
            &mut MolSettings::default(),
            &EditorRotateState::default(),
        );
        assert!(state.last_error.is_none(), "{:?}", state.last_error);

        undo_last_fragment(&mut state, &mut mol, &mut MolSettings::default());

        assert_eq!(mol.atoms, atoms_before, "the hydrogen must be back");
        assert_eq!(state.zmat.len(), zmat_before.len());
    }
}

/// Read one of the sample structures in `examples/`.
///
/// Resolved against the crate root rather than the working directory, so the
/// tests do not depend on where cargo happens to be run from -- and not
/// against anyone's Desktop, which is where these fixtures used to live.
#[cfg(test)]
fn read_example(name: &str) -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join(name);
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read example {}: {e}", path.display()))
}

#[cfg(test)]
mod real_file_replace_tests {
    use super::*;

    fn load(name: &str) -> (Vec<String>, Vec<Vec3>) {
        let text = read_example(name);
        let (_n, symbols, coords) = parse_xyz_angstrom(&text);
        (
            symbols,
            coords
                .into_iter()
                .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
                .collect(),
        )
    }

    /// Every pairwise distance -- the shape, independent of how the molecule
    /// happens to be oriented or centred. `zmat_to_xyz` always rebuilds a
    /// structure in its own frame (first atom at the origin, second along the
    /// placement axis), so absolute coordinates are the wrong thing to compare.
    fn distance_matrix(pos: &[Vec3], skip: Option<usize>) -> Vec<f32> {
        let keep: Vec<usize> = (0..pos.len()).filter(|i| Some(*i) != skip).collect();
        let mut out = Vec::new();
        for (a, &i) in keep.iter().enumerate() {
            for &j in &keep[a + 1..] {
                out.push(pos[i].distance(pos[j]));
            }
        }
        out
    }

    /// The round trip may re-orient the molecule, but it must not deform it.
    #[test]
    fn the_zmat_round_trip_preserves_the_shape() {
        let (symbols, positions) = load("phenyl-OCH3.xyz");
        let zmat = crate::molecule_builder::zmat2xyz::xyz_to_zmat(&symbols, &positions);
        let rebuilt = zmat2xyz::zmat_to_xyz(&zmat);

        let before = distance_matrix(&positions, None);
        let after = distance_matrix(&rebuilt, None);
        let worst = before
            .iter()
            .zip(&after)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f32, f32::max);
        assert!(
            worst < 1.0e-2,
            "the xyz -> zmat -> xyz round trip deformed the molecule by {worst:.4} A"
        );
    }

    /// The realistic path: load a file, sync the Z-matrix, then replace a ring
    /// hydrogen with a methyl. The rest of the molecule must keep its exact
    /// geometry; only the substituted hydrogen goes away.
    #[test]
    fn replacing_a_ring_hydrogen_keeps_the_rest_of_the_shape() {
        let (symbols, positions) = load("phenyl-OCH3.xyz");
        let zmat = crate::molecule_builder::zmat2xyz::xyz_to_zmat(&symbols, &positions);
        let rebuilt = zmat2xyz::zmat_to_xyz(&zmat);
        let mut mol = Molecule {
            atoms: symbols.clone(),
            pos: rebuilt.clone(),
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        mol.recompute_bonds(2.0, 3.0);

        let replaced = 6usize;
        assert_eq!(mol.atoms[replaced], "H");

        let mut state = ZMatrixBuilderState::default();
        state.zmat = zmat;
        state.frag_name = "-CH3".to_string();
        state.frag_mode = FragmentInsertMode::Replace;
        state.selected_index = Some(replaced);

        commit_fragment_replace(
            &mut state,
            &mut mol,
            &mut MolSettings::default(),
            &EditorRotateState::default(),
        );
        assert!(state.last_error.is_none(), "{:?}", state.last_error);
        assert_eq!(
            state.zmat.len(),
            mol.atoms.len(),
            "the rebuilt Z-matrix must describe the molecule it produced"
        );

        // Every surviving original atom, in its new position, must sit exactly
        // where it did before.
        let mut worst = (0usize, 0.0f32);
        for k in 0..positions.len() {
            if k == replaced {
                continue;
            }
            let now = if k < replaced { k } else { k - 1 };
            let moved = rebuilt[k].distance(mol.pos[now]);
            if moved > worst.1 {
                worst = (k, moved);
            }
        }
        assert!(
            worst.1 < 1.0e-2,
            "original atom {} ({}) moved {:.4} A; replacing one hydrogen must \
             not disturb the rest of the structure",
            worst.0,
            symbols[worst.0],
            worst.1
        );
    }

    /// The user's exact report: pick atom 8 in the panel's 1-based numbering
    /// (index 7, an aromatic hydrogen) and replace it with a methyl. The
    /// methyl must land on the ring carbon that actually held that hydrogen --
    /// not on whichever atom happened to precede it in the file.
    #[test]
    fn replacing_ring_hydrogen_number_eight_attaches_to_its_own_carbon() {
        let (symbols, positions) = load("phenyl-OCH3.xyz");
        let zmat = crate::molecule_builder::zmat2xyz::xyz_to_zmat(&symbols, &positions);
        let mut mol = Molecule {
            atoms: symbols.clone(),
            pos: zmat2xyz::zmat_to_xyz(&zmat),
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        mol.recompute_bonds(2.0, 3.0);
        let before = mol.pos.clone();

        let replaced = 7usize; // atom 8 in the panel
        let ring_carbon = 2usize;
        assert_eq!(mol.atoms[replaced], "H");
        assert!(
            (before[replaced].distance(before[ring_carbon]) - 1.089).abs() < 0.01,
            "fixture check: atom {replaced} should be the hydrogen on atom {ring_carbon}"
        );
        // The old code used this row's Z-matrix bond reference instead, which
        // for a file-order Z-matrix is simply the previous line -- another
        // hydrogen. That is what put the methyl on a hydrogen.
        assert_eq!(zmat[replaced].bond_ref.map(|r| r - 1), Some(6));
        assert_eq!(mol.atoms[6], "H");

        let mut state = ZMatrixBuilderState::default();
        state.zmat = zmat;
        state.frag_name = "-CH3".to_string();
        state.frag_mode = FragmentInsertMode::Replace;
        state.selected_index = Some(replaced);

        commit_fragment_replace(
            &mut state,
            &mut mol,
            &mut MolSettings::default(),
            &EditorRotateState::default(),
        );
        assert!(state.last_error.is_none(), "{:?}", state.last_error);

        // One atom removed, four added.
        assert_eq!(mol.atoms.len(), symbols.len() - 1 + 4);

        // The methyl carbon is the first appended row.
        let methyl_c = symbols.len() - 1;
        assert_eq!(mol.atoms[methyl_c], "C");
        let d = mol.pos[methyl_c].distance(mol.pos[ring_carbon]);
        assert!(
            (d - 1.51).abs() < 0.15,
            "the methyl carbon sits {d:.3} A from the ring carbon that held the \
             hydrogen; it should be a single bond away"
        );

        // It must not be bonded to a hydrogen -- that was the bug.
        for h in 0..methyl_c {
            if mol.atoms[h] != "H" {
                continue;
            }
            assert!(
                mol.pos[methyl_c].distance(mol.pos[h]) > 1.3,
                "the methyl carbon is {:.3} A from hydrogen {h}: it has been \
                 bonded to a hydrogen instead of to the ring",
                mol.pos[methyl_c].distance(mol.pos[h])
            );
        }

        // Every other ring carbon keeps its own hydrogen; the substituted one
        // no longer has one.
        let hydrogens_on = |mol: &Molecule, c: usize| -> usize {
            (0..mol.atoms.len())
                .filter(|&j| mol.atoms[j] == "H" && mol.pos[c].distance(mol.pos[j]) < 1.3)
                .count()
        };
        for c in 1..6 {
            let expected = usize::from(c != ring_carbon);
            assert_eq!(
                hydrogens_on(&mol, c),
                expected,
                "ring carbon {c} should have {expected} hydrogen(s)"
            );
        }

        // And nothing else moved.
        for k in 0..positions.len() {
            if k == replaced {
                continue;
            }
            let now = if k < replaced { k } else { k - 1 };
            assert!(
                before[k].distance(mol.pos[now]) < 1.0e-3,
                "original atom {k} ({}) moved {:.4} A",
                symbols[k],
                before[k].distance(mol.pos[now])
            );
        }
    }

    /// Loading a structure must leave the builder ready to work on it, with
    /// no manual "Sync from Molecule" step.
    #[test]
    fn loading_a_structure_syncs_the_zmat() {
        let (symbols, positions) = load("phenyl-OCH3.xyz");
        let mut mol = Molecule {
            atoms: symbols.clone(),
            pos: positions.clone(),
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        mol.recompute_bonds(2.0, 3.0);

        let mut app = App::new();
        app.add_message::<MoleculeChanged>()
            .insert_resource(mol)
            .init_resource::<ZMatrixBuilderState>()
            .add_systems(Update, sync_builder_on_structure_load);
        app.world_mut()
            .write_message(MoleculeChanged::parse_xyz(true));
        app.update();

        let state = app.world().resource::<ZMatrixBuilderState>();
        assert_eq!(
            state.zmat.len(),
            symbols.len(),
            "the Z-matrix should describe the molecule that was just loaded"
        );
        assert_eq!(state.selected_index, None);
        assert!(!state.frag_undo_visible, "no undo carries over a load");
    }

    /// An unsynced Z-matrix must never be allowed to rewrite the displayed
    /// molecule: it describes something else entirely.
    #[test]
    fn replacing_without_a_synced_zmat_is_refused() {
        let (symbols, positions) = load("phenyl-OCH3.xyz");
        let mut mol = Molecule {
            atoms: symbols.clone(),
            pos: positions.clone(),
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        mol.recompute_bonds(2.0, 3.0);

        // Deliberately no Z-matrix: an empty or stale one describes a
        // different molecule, and committing would rewrite the display from it.
        let mut state = ZMatrixBuilderState::default();
        state.frag_name = "-CH3".to_string();
        state.frag_mode = FragmentInsertMode::Replace;
        state.selected_index = Some(6);

        commit_fragment_replace(
            &mut state,
            &mut mol,
            &mut MolSettings::default(),
            &EditorRotateState::default(),
        );

        assert!(
            state.last_error.is_some(),
            "an out-of-sync Z-matrix must be refused, not acted on"
        );
        assert_eq!(mol.atoms, symbols, "the molecule must be left untouched");
        assert_eq!(mol.pos, positions);
    }
}

#[cfg(test)]
mod convention_probe {
    use super::*;

    /// `derive_internals` measures angle/dihedral with
    /// `internal_from_direction`; `zmat_to_xyz` rebuilds positions with its
    /// own formula. If those two disagree, re-measuring a row from correct
    /// coordinates silently moves the atom.
    #[test]
    fn derive_internals_agrees_with_zmat_to_xyz() {
        let text = read_example("phenyl-OCH3.xyz");
        let (_n, symbols, coords) = parse_xyz_angstrom(&text);
        let positions: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let zmat = crate::molecule_builder::zmat2xyz::xyz_to_zmat(&symbols, &positions);
        let rebuilt = zmat2xyz::zmat_to_xyz(&zmat);

        // Re-measure every row from its own rebuilt coordinates, changing no
        // references at all. This must be a no-op.
        let mut remeasured = zmat.clone();
        let all: Vec<usize> = (0..zmat.len()).collect();
        attach::derive_internals(&rebuilt, &mut remeasured, &all);

        let after_coords = zmat2xyz::zmat_to_xyz(&remeasured);
        let worst = rebuilt
            .iter()
            .zip(&after_coords)
            .map(|(a, b)| a.distance(*b))
            .fold(0.0f32, f32::max);
        assert!(
            worst < 1.0e-2,
            "re-measuring the internals from the very coordinates they describe \
             moved atoms by {worst:.4} A: the two conventions disagree"
        );
    }
}

#[cfg(test)]
mod editor_flow_tests {
    use super::*;

    fn load(name: &str) -> Molecule {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join(name);
        let mut mol = Molecule::from_xyz(&std::fs::read_to_string(&path).unwrap());
        mol.recompute_bonds(1.2, 2.5);
        mol
    }

    /// The exact sequence the UI performs: load a structure, sync the
    /// builder, click an atom, press Add Fragment.
    #[test]
    fn connecting_a_fragment_to_a_clicked_atom_grows_the_molecule() {
        let mut mol = load("phenyl-OCH3.xyz");
        let mut zmat_state = ZMatrixBuilderState::default();
        let mut settings = MolSettings::default();

        // What `sync_builder_on_structure_load` does.
        zmat_state.zmat = zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos);
        // What a viewport click does: pick a ring hydrogen.
        let h = mol
            .atoms
            .iter()
            .position(|s| s == "H")
            .expect("the ring has hydrogens");
        set_selected_atom(&mut zmat_state, &mol.atoms, Some(h));
        assert!(zmat_state.selected_index.is_some(), "the click must select");

        let before = mol.atoms.len();
        zmat_state.frag_mode = FragmentInsertMode::Connect;
        zmat_state.frag_name = "-CH3".to_string();
        commit_fragment_connect(&mut zmat_state, &mut mol, &mut settings);

        assert!(
            zmat_state.last_error.is_none(),
            "connect reported: {:?}",
            zmat_state.last_error
        );
        assert!(
            mol.atoms.len() > before,
            "the molecule did not grow: {before} -> {}",
            mol.atoms.len()
        );
    }
}
