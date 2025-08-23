use bevy::prelude::*;
use bevy_egui::egui;

use crate::measurements::{MeasurePair, Measurements};
use crate::molecule::{Molecule, covalent_radius_angstrom};
use crate::settings::MolSettings;

/// Call this from your existing right-side panel in `ui.rs`.
pub fn measurements_panel(
    ui: &mut egui::Ui,
    measurements: &mut Measurements,
    mol: &Molecule,
    settings: &MolSettings,
) {
    // =======================
    // Distance panel (existing)
    // =======================
    ui.collapsing("Distance", |ui| {
        // Controls
        ui.horizontal(|ui| {
            if ui
                .toggle_value(&mut measurements.is_active, "Activate picking")
                .clicked()
            {
                // When activating distance, it's okay if angle/dihedral stay on,
                // picking system gives priority to dihedral > angle > distance.
            }
            if ui.button("New distance").clicked() {
                measurements.is_active = true;
                measurements.pending = None;
            }
            if ui.button("Clear all").clicked() {
                measurements.clear();
            }
        });

        // Pending first atom indicator
        if let Some(i) = measurements.pending {
            let name = mol.atoms.get(i).map(|s| s.as_str()).unwrap_or("?");
            ui.horizontal(|ui| {
                ui.colored_label(
                    egui::Color32::from_rgb(255, 64, 200),
                    format!("First atom: {}({})", name, i + 1),
                );
                if ui.button("Cancel").clicked() {
                    measurements.pending = None;
                }
            });
        }

        // Global style
        ui.collapsing("Change Distance Style", |ui| {
            ui.separator();
            ui.label("Global style");
            ui.add(
                egui::Slider::new(&mut measurements.line_width, 1.0..=50.0)
                    .text("Line width"),
            );
            ui.add(
                egui::Slider::new(&mut measurements.dash_line_scale, 0.05..=4.0)
                    .text("Dash line scale"),
            );
            ui.add(
                egui::Slider::new(&mut measurements.dash_gap_scale, 0.05..=4.0)
                    .text("Dash gap scale"),
            );
            ui.checkbox(&mut measurements.show_labels, "Show labels next to line");

            // Global color picker
            {
                let s = measurements.color.to_srgba();
                let mut eg = egui::Color32::from_rgba_premultiplied(
                    (s.red * 255.0).round() as u8,
                    (s.green * 255.0).round() as u8,
                    (s.blue * 255.0).round() as u8,
                    (s.alpha * 255.0).round() as u8,
                );
                if ui.color_edit_button_srgba(&mut eg).on_hover_text("Global color").changed() {
                    measurements.color = Color::srgba(
                        eg.r() as f32 / 255.0,
                        eg.g() as f32 / 255.0,
                        eg.b() as f32 / 255.0,
                        eg.a() as f32 / 255.0,
                    );
                }
            }

            // Picking configuration (Å slack + pixel radius)
            ui.separator();
            ui.label("Picking");
            ui.add(
                egui::Slider::new(&mut measurements.pick_radius, 0.0..=2.0)
                    .text("Sphere slack (Å)"),
            );
            ui.add(
                egui::Slider::new(&mut measurements.pick_radius_px, 6.0..=40.0)
                    .text("Pixel pick radius (px)"),
            );
            ui.small(format!(
                "Atom scale: {:.2}  (effective world slack ~ {:.2} Å). Pixel radius affects screen-space picking.",
                settings.atom_scale,
                measurements.pick_radius * settings.atom_scale
            ));
        }); // end Change Distance Style

        // Pairs list
        ui.separator();
        if measurements.pairs.is_empty() {
            ui.weak("No measurement pairs yet. Activate picking, then click two atoms.");
        } else {
            egui::ScrollArea::vertical()
                .id_salt("measure_pairs_scroll")
                .auto_shrink([false, false])
                .max_height(240.0)
                .show(ui, |ui| {
                    let ids: Vec<u32> = measurements.pairs.iter().map(|p| p.id).collect();
                    let mut to_delete: Vec<u32> = Vec::new();

                    for id in ids {
                        if let Some(idx) = measurements.pairs.iter().position(|p| p.id == id) {
                            if pair_row(ui, &mut measurements.pairs[idx], mol) {
                                to_delete.push(id); // delete requested
                            }
                            ui.separator();
                        }
                    }

                    if !to_delete.is_empty() {
                        measurements.pairs.retain(|p| !to_delete.contains(&p.id));
                    }
                });
        }
    });

    // =======================
    // Angle panel (NEW)
    // =======================
    ui.collapsing("Angle", |ui| {
        ui.horizontal(|ui| {
            if ui.toggle_value(&mut measurements.angle_active, "Activate picking").clicked() {
                // When activating angles, we leave others as-is; priority is handled in system.
            }
            if ui.button("New angle").clicked() {
                measurements.angle_active = true;
                measurements.pending_angle.clear();
            }
            if ui.button("Clear all").clicked() {
                measurements.clear_angles();
            }
        });

        if !measurements.pending_angle.is_empty() {
            let mut labels: Vec<String> = Vec::new();
            for &idx in &measurements.pending_angle {
                let name = mol.atoms.get(idx).map(|s| s.as_str()).unwrap_or("?");
                labels.push(format!("{}({})", name, idx + 1));
            }
            ui.colored_label(
                egui::Color32::from_rgb(200, 200, 80),
                format!("Picked: {}", labels.join(" - ")),
            );
        }

        ui.separator();
        if measurements.angles.is_empty() {
            ui.weak("No angles yet. Activate picking, then click three atoms (A-B-C).");
        } else {
            egui::ScrollArea::vertical()
                .id_salt("angle_list_scroll")
                .max_height(220.0)
                .show(ui, |ui| {
                    for a in &measurements.angles {
                        let na = mol.atoms.get(a.a).map(|s| s.as_str()).unwrap_or("?");
                        let nb = mol.atoms.get(a.b).map(|s| s.as_str()).unwrap_or("?");
                        let nc = mol.atoms.get(a.c).map(|s| s.as_str()).unwrap_or("?");
                        ui.monospace(format!(
                            "{}({}) - {}({}) - {}({})  |  {:.1}°",
                            na, a.a + 1, nb, a.b + 1, nc, a.c + 1, a.degrees
                        ));
                    }
                });
        }
    });

    // =======================
    // Dihedral panel (NEW)
    // =======================
    ui.collapsing("Dihedral", |ui| {
        ui.horizontal(|ui| {
            if ui
                .toggle_value(&mut measurements.dihedral_active, "Activate picking")
                .clicked()
            {
                // Priority handled in picking system.
            }
            if ui.button("New dihedral").clicked() {
                measurements.dihedral_active = true;
                measurements.pending_dihedral.clear();
            }
            if ui.button("Clear all").clicked() {
                measurements.clear_dihedrals();
            }
        });

        if !measurements.pending_dihedral.is_empty() {
            let mut labels: Vec<String> = Vec::new();
            for &idx in &measurements.pending_dihedral {
                let name = mol.atoms.get(idx).map(|s| s.as_str()).unwrap_or("?");
                labels.push(format!("{}({})", name, idx + 1));
            }
            ui.colored_label(
                egui::Color32::from_rgb(80, 200, 200),
                format!("Picked: {}", labels.join(" - ")),
            );
        }

        ui.separator();
        if measurements.dihedrals.is_empty() {
            ui.weak("No dihedrals yet. Activate picking, then click four atoms (A-B-C-D).");
        } else {
            egui::ScrollArea::vertical()
                .id_salt("dihedral_list_scroll")
                .max_height(220.0)
                .show(ui, |ui| {
                    for d in &measurements.dihedrals {
                        let na = mol.atoms.get(d.a).map(|s| s.as_str()).unwrap_or("?");
                        let nb = mol.atoms.get(d.b).map(|s| s.as_str()).unwrap_or("?");
                        let nc = mol.atoms.get(d.c).map(|s| s.as_str()).unwrap_or("?");
                        let nd = mol.atoms.get(d.d).map(|s| s.as_str()).unwrap_or("?");
                        ui.monospace(format!(
                            "{}({}) - {}({}) - {}({}) - {}({})  |  {:+.1}°",
                            na, d.a + 1, nb, d.b + 1, nc, d.c + 1, nd, d.d + 1, d.degrees
                        ));
                    }
                });
        }
    });
}

/// A Bevy system that paints distance labels near the measured lines in a 2D egui overlay.
pub fn distance_labels_overlay(
    mut contexts: bevy_egui::EguiContexts,
    q_cam: Query<(&Camera, &GlobalTransform), (With<Camera3d>, With<crate::scene::MainCamera>)>,
    windows: Query<&Window>,
    mol: Option<Res<Molecule>>,
    settings: Option<Res<MolSettings>>,
    measurements: Option<Res<Measurements>>,
) {
    let (Some(mol), Some(settings), Some(measurements)) = (mol, settings, measurements) else { return; };
    if !measurements.show_labels || measurements.pairs.is_empty() {
        return;
    }

    // Need egui ctx + camera + window to project 3D -> screen
    let Ok(ctx) = contexts.ctx_mut() else { return; };
    let Ok((cam, cam_xform)) = q_cam.single() else { return; };
    let Ok(window) = windows.single() else { return; };

    egui::Area::new(egui::Id::new("measurement_labels_overlay"))
        .movable(false)
        .interactable(false)
        .order(egui::Order::Tooltip)
        .show(&ctx, |ui| {
            let painter = ui.painter();
            let scale = window.scale_factor() as f32;

            for pair in measurements.pairs.iter().filter(|p| p.visible) {
                let n = mol.pos.len();
                if pair.a >= n || pair.b >= n {
                    continue;
                }

                // Recompute the visible endpoints the same way as gizmo lines
                let p0 = mol.pos[pair.a];
                let p1 = mol.pos[pair.b];
                let d = p1 - p0;
                let len = d.length();
                if len <= 1e-4 {
                    continue;
                }
                let dn = d / len;

                let ri = covalent_radius_angstrom(&mol.atoms[pair.a]) * settings.atom_scale;
                let rj = covalent_radius_angstrom(&mol.atoms[pair.b]) * settings.atom_scale;

                let surface_inset_frac = 0.12;
                let inset_i = (ri * surface_inset_frac).min(ri * 0.9);
                let inset_j = (rj * surface_inset_frac).min(rj * 0.9);

                let a = p0 + dn * (ri - inset_i);
                let b = p1 - dn * (rj - inset_j);
                if (b - a).length() <= 1e-4 {
                    continue;
                }

                // Midpoint in world
                let mid = (a + b) * 0.5;

                if let Ok(screen_px) = cam.world_to_viewport(cam_xform, mid) {
                    let pos = egui::pos2(screen_px.x / scale, screen_px.y / scale);

                    // Text: distance in Å with three decimals (no unit in overlay per your ask)
                    let text = format!("{:.3}", pair.distance);
                    let galley = ctx.fonts(|f| {
                        f.layout_no_wrap(
                            text.clone(),
                            egui::FontId::proportional(16.0),
                            egui::Color32::WHITE,
                        )
                    });

                    // Slight offset so label doesn't sit exactly on the dashed line
                    let offset = egui::vec2(6.0, -6.0);
                    painter.galley(pos + offset, galley, egui::Color32::WHITE);
                }
            }
        });
}

/// Row UI for a distance pair; returns true if user clicked delete.
fn pair_row(ui: &mut egui::Ui, pair: &mut MeasurePair, mol: &Molecule) -> bool {
    let name_a = mol.atoms.get(pair.a).map(|s| s.as_str()).unwrap_or("?");
    let name_b = mol.atoms.get(pair.b).map(|s| s.as_str()).unwrap_or("?");
    let mut delete_me = false;

    ui.horizontal(|ui| {
        // Show/Hide (affects both line and overlay label)
        let eye = if pair.visible { "Hide" } else { "Show" };
        if ui.button(eye).clicked() {
            pair.visible = !pair.visible;
        }

        // Pair info + distance (1-based indices in UI)
        ui.label(format!(
            "{}({}) — {}({})",
            name_a,
            pair.a + 1,
            name_b,
            pair.b + 1
        ));
        ui.separator();
        ui.monospace(format!("{:.3} Å", pair.distance));

        ui.separator();

        // Delete
        if ui.button("Delete").clicked() {
            delete_me = true;
        }
    });

    delete_me
}

