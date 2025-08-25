use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};

use crate::events::{AtomPicked, ToolKind};
use crate::molecule::{Molecule, covalent_radius_angstrom};
use crate::settings::MolSettings;
use super::rotator::{rotate_side, split_sides, RotateSide};

#[derive(Default, Reflect, GizmoConfigGroup)]
pub struct BuilderGizmos;

#[derive(Resource, Default, Clone)]
pub struct EditorRotateState {
    pub active: bool,                     // capture 3D clicks when true
    pub picks: Vec<usize>,                // [A, B, side-indicator (any atom on that side)]
    pub chosen_side: Option<RotateSide>,  // inferred from 3rd pick
    pub angle_deg: f32,                   // degrees
    pub last_snapshot: Option<Vec<Vec3>>, // undo buffer
    pub bond_th_hx: f32,                  // Å
    pub bond_th_xx: f32,                  // Å
}

impl EditorRotateState {
    pub fn reset(&mut self) {
        *self = Self {
            active: false,
            picks: Vec::new(),
            chosen_side: None,
            angle_deg: 0.0,
            last_snapshot: None,
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
    mol: Option<Res<Molecule>>,
    mut state: ResMut<EditorRotateState>,
) {
    let Some(mol) = mol else { return; };
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
    mol: Option<Res<Molecule>>,
    settings: Res<MolSettings>,
) {
    let Some(mol) = mol else { return; };
    if !state.active { return; }
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

/// Left-side, resizable SidePanel with a collapsible “Molecule Editor” section.
/// NOTE: We now do a clean jump: no immediate atom transform sync; instead we
/// update `mol.pos` + `mol.bonds` and flag `settings.dirty = true` so the scene
/// rebuilds atoms & bonds together next frame (no mismatched frame).
pub fn builder_ui_panel(
    mut contexts: EguiContexts,
    mut state: ResMut<EditorRotateState>,
    mut mol: ResMut<Molecule>,
    mut settings: ResMut<MolSettings>, // mark scene dirty to rebuild
) {
    let Ok(ctx) = contexts.ctx_mut() else { return; };

    egui::SidePanel::left("molecule_editor_panel")
        .resizable(true)
        .show(&ctx, |ui| {
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
                    if ui.button("Apply").clicked() {
                        if state.picks.len() >= 3 {
                            let a = state.picks[0];
                            let b = state.picks[1];
                            if let Some(side) = state.chosen_side {
                                // Save undo snapshot
                                state.last_snapshot = Some(mol.pos.clone());

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
                        if let Some(snapshot) = &mut state.last_snapshot {
                            mol.pos = snapshot.clone();
                            mol.recompute_bonds(2.0, 3.0);
                            settings.dirty = true; // rebuild to keep atoms+bonds in sync
                        }
                    }
                    if ui.button("Cancel").clicked() {
                        state.picks.clear();
                        state.chosen_side = None;
                        state.angle_deg = 0.0;
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

