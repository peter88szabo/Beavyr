use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};

use crate::events::{AtomPicked, ToolKind};
use crate::molecule::{Molecule, covalent_radius_angstrom};
use crate::settings::MolSettings;
use super::rotator::{bend_side, rotate_side, split_sides, translate_side, RotateSide};
use super::zmat2xyz::{self, ZAtom};

const DEFAULT_ANGLE_DEG: f64 = 109.45;
const DEFAULT_DIHEDRAL_DEG: f64 = 120.0;

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
    pub pick_active: bool,
    pub pick_indices: Vec<usize>,
    pub pick_hint: Option<String>,
    pub edit_rows: Vec<ZAtomEditRow>,
    pub edit_refresh: bool,
    pub edit_preview: Vec<usize>,
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
            new_symbol: "C".to_string(),
            pick_active: false,
            pick_indices: Vec::new(),
            pick_hint: None,
            edit_rows: Vec::new(),
            edit_refresh: false,
            edit_preview: Vec::new(),
            original_atoms: None,
            original_pos: None,
            last_error: None,
        }
    }
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
    settings.dirty = true;
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
    mut ev: EventReader<AtomPicked>,
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
                    let bond_len = (r_new + r_ref) as f64;

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
                            angle_deg: DEFAULT_ANGLE_DEG,
                            dihedral_ref: None,
                            dihedral_deg: 0.0,
                        },
                        _ => ZAtom {
                            symbol,
                            bond_ref: Some(bond_ref + 1),
                            bond_len,
                            angle_ref: angle_ref.map(|v| v + 1),
                            angle_deg: DEFAULT_ANGLE_DEG,
                            dihedral_ref: dihedral_ref.map(|v| v + 1),
                            dihedral_deg: DEFAULT_DIHEDRAL_DEG,
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
                    settings.dirty = true;
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
/// update `mol.pos` + `mol.bonds` and flag `settings.dirty = true` so the scene
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
                        settings.dirty = true;
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
                            settings.dirty = true;
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
                    let col_idx = [20.0, 28.0, 28.0, 68.0, 28.0, 60.0, 28.0, 60.0];
                    let row_height = 22.0;
                    let max_height = row_height * 15.0;
                    egui::ScrollArea::vertical()
                        .max_height(max_height)
                        .show(ui, |ui| {
                            egui::Grid::new("zmat_rows")
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

                                    for (i, row) in zmat_state.edit_rows.iter_mut().enumerate() {
                                        let idx_label = format!("{}", i + 1);
                                        ui.add_sized([col_idx[0], 0.0], egui::Label::new(idx_label));
                                        let sym_edit = egui::TextEdit::singleline(&mut row.symbol)
                                            .desired_width(col_idx[1] - 4.0);
                                        if ui.add_sized([col_idx[1], 0.0], sym_edit).changed() {
                                            let trimmed = row.symbol.trim();
                                            if !trimmed.is_empty() {
                                                row.symbol = trimmed.to_string();
                                            }
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
                                });
                        });
                    if let Some(preview) = preview_indices {
                        zmat_state.edit_preview = preview;
                    }
                    if ui.button("Apply Changes").clicked() {
                        apply_zmat_edits(&mut zmat_state, &mut mol, &mut settings);
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
                                settings.dirty = true;
                                zmat_state.edit_refresh = true;
                                zmat_state.last_error = None;
                            }
                        } else {
                            zmat_state.pick_active = true;
                            zmat_state.pick_indices.clear();
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
                        zmat_state.zmat.pop();
                        zmat_state.edit_refresh = true;
                    }
                });

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

                if let Some(err) = &zmat_state.last_error {
                    ui.colored_label(egui::Color32::LIGHT_RED, err);
                }
            });

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

                                // Update bond list to new geometry
                                mol.recompute_bonds(2.0, 3.0);

                                // Jump-style: rebuild whole scene (atoms & bonds) next frame
                                settings.dirty = true;
                            }
                        }
                    }
                    if ui.button("Undo").clicked() {
                        if let Some(snapshot) = &mut state.last_rotate_snapshot {
                            mol.pos = snapshot.clone();
                            mol.recompute_bonds(2.0, 3.0);
                            settings.dirty = true; // rebuild to keep atoms+bonds in sync
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
                                    mol.recompute_bonds(2.0, 3.0);
                                    settings.dirty = true;
                                }
                            }
                        }
                    }
                    if ui.button("Undo").clicked() {
                        if let Some(snapshot) = &mut state.last_translate_snapshot {
                            mol.pos = snapshot.clone();
                            mol.recompute_bonds(2.0, 3.0);
                            settings.dirty = true;
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
                                mol.recompute_bonds(2.0, 3.0);
                                settings.dirty = true;
                            }
                        }
                    }
                    if ui.button("Undo").clicked() {
                        if let Some(snapshot) = &mut state.last_bend_snapshot {
                            mol.pos = snapshot.clone();
                            mol.recompute_bonds(2.0, 3.0);
                            settings.dirty = true;
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
