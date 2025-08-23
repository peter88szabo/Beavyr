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
    // Distance panel
    // =======================
    ui.collapsing("Distance", |ui| {
        // Controls
        ui.horizontal(|ui| {
            // Highlighted status (not a toggle)
            active_status_pill(ui, measurements.is_active);

            // Start exclusive picking for distance
            if ui.button("New distance").clicked() {
                measurements.is_active = true;
                measurements.pending = None;

                // Exclusivity
                measurements.angle_active = false;
                measurements.pending_angle.clear();
                measurements.dihedral_active = false;
                measurements.pending_dihedral.clear();
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
                    measurements.is_active = false; // cancel picking
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
            /*
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
            */
        }); // end Change Distance Style

        // Pairs list
        ui.separator();
        if measurements.pairs.is_empty() {
            ui.weak("No measurement pairs yet. Activate picking, then click two atoms.");
        } else {
            egui::ScrollArea::vertical()
                .id_salt("measure_pairs_scroll")
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
    // Angle panel
    // =======================
    ui.add_space(4.0);
    ui.collapsing("Angle", |ui| {
        ui.horizontal(|ui| {
            active_status_pill(ui, measurements.angle_active);

            if ui.button("New angle").clicked() {
                measurements.angle_active = true;
                measurements.pending_angle.clear();

                measurements.is_active = false;
                measurements.pending = None;
                measurements.dihedral_active = false;
                measurements.pending_dihedral.clear();
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
            ui.horizontal(|ui| {
                ui.colored_label(
                    egui::Color32::from_rgb(200, 200, 80),
                    format!("Picked: {}", labels.join(" - ")),
                );
                if ui.button("Cancel").clicked() {
                    measurements.pending_angle.clear();
                    measurements.angle_active = false;
                }
            });
        }

        ui.separator();
        if measurements.angles.is_empty() {
            ui.weak("No angles yet. Activate picking, then click three atoms (A-B-C).");
        } else {
            egui::ScrollArea::vertical()
                .id_salt("angle_list_scroll")
                .max_height(220.0)
                .show(ui, |ui| {
                    let ids: Vec<u32> = measurements.angles.iter().map(|a| a.id).collect();
                    let mut to_delete: Vec<u32> = Vec::new();
                    let mut req_set: Option<Vec<usize>> = None;
                    let mut req_clear: bool = false;

                    for id in ids {
                        if let Some(idx) = measurements.angles.iter().position(|a| a.id == id) {
                            let a = &measurements.angles[idx];
                            let na = mol.atoms.get(a.a).map(|s| s.as_str()).unwrap_or("?");
                            let nb = mol.atoms.get(a.b).map(|s| s.as_str()).unwrap_or("?");
                            let nc = mol.atoms.get(a.c).map(|s| s.as_str()).unwrap_or("?");

                            let is_hl = is_highlighted_for(
                                &measurements.preview_highlight,
                                &[a.a, a.b, a.c],
                            );

                            ui.horizontal(|ui| {
                                // Highlight toggle at the BEGINNING of the row
                                let label = if is_hl { "Hide" } else { "Show" };
                                let mut btn = egui::Button::new(label);
                                if is_hl {
                                    btn = btn.fill(egui::Color32::from_rgb(60, 120, 200));
                                }
                                if ui.add(btn).clicked() {
                                    if is_hl {
                                        req_clear = true;
                                    } else {
                                        req_set = Some(vec![a.a, a.b, a.c]);
                                    }
                                }

                                ui.add_space(8.0);
                                ui.monospace(format!(
                                    "{}({})-{}({})-{}({}) | {:.1}°",
                                    na, a.a + 1, nb, a.b + 1, nc, a.c + 1, a.degrees
                                ));
                                ui.add_space(8.0);
                                if ui.button("Delete").clicked() {
                                    to_delete.push(id);
                                }
                            });
                            ui.separator();
                        }
                    }

                    if !to_delete.is_empty() {
                        measurements.angles.retain(|a| !to_delete.contains(&a.id));
                    }
                    if req_clear {
                        measurements.preview_highlight = None;
                    } else if let Some(v) = req_set {
                        measurements.preview_highlight = Some(v);
                    }
                });
        }
    });

    // =======================
    // Dihedral panel
    // =======================
    ui.add_space(4.0);
    ui.collapsing("Dihedral", |ui| {
        ui.horizontal(|ui| {
            active_status_pill(ui, measurements.dihedral_active);

            if ui.button("New dihedral").clicked() {
                measurements.dihedral_active = true;
                measurements.pending_dihedral.clear();

                measurements.is_active = false;
                measurements.pending = None;
                measurements.angle_active = false;
                measurements.pending_angle.clear();
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
            ui.horizontal(|ui| {
                ui.colored_label(
                    egui::Color32::from_rgb(80, 200, 200),
                    format!("Picked: {}", labels.join(" - ")),
                );
                if ui.button("Cancel").clicked() {
                    measurements.pending_dihedral.clear();
                    measurements.dihedral_active = false;
                }
            });
        }

        ui.separator();
        if measurements.dihedrals.is_empty() {
            ui.weak("No dihedrals yet. Activate picking, then click four atoms (A-B-C-D).");
        } else {
            egui::ScrollArea::vertical()
                .id_salt("dihedral_list_scroll")
                .max_height(220.0)
                .show(ui, |ui| {
                    let ids: Vec<u32> = measurements.dihedrals.iter().map(|d| d.id).collect();
                    let mut to_delete: Vec<u32> = Vec::new();
                    let mut req_set: Option<Vec<usize>> = None;
                    let mut req_clear: bool = false;

                    for id in ids {
                        if let Some(idx) = measurements.dihedrals.iter().position(|d| d.id == id) {
                            let dmeas = &measurements.dihedrals[idx];
                            let na = mol.atoms.get(dmeas.a).map(|s| s.as_str()).unwrap_or("?");
                            let nb = mol.atoms.get(dmeas.b).map(|s| s.as_str()).unwrap_or("?");
                            let nc = mol.atoms.get(dmeas.c).map(|s| s.as_str()).unwrap_or("?");
                            let nd = mol.atoms.get(dmeas.d).map(|s| s.as_str()).unwrap_or("?");

                            let is_hl = is_highlighted_for(
                                &measurements.preview_highlight,
                                &[dmeas.a, dmeas.b, dmeas.c, dmeas.d],
                            );

                            ui.horizontal(|ui| {
                                // Highlight toggle at the BEGINNING of the row
                                let label = if is_hl { "Hide" } else { "Show" };
                                let mut btn = egui::Button::new(label);
                                if is_hl {
                                    btn = btn.fill(egui::Color32::from_rgb(60, 120, 200));
                                }
                                if ui.add(btn).clicked() {
                                    if is_hl {
                                        req_clear = true;
                                    } else {
                                        req_set = Some(vec![dmeas.a, dmeas.b, dmeas.c, dmeas.d]);
                                    }
                                }

                                ui.add_space(8.0);
                                ui.monospace(format!(
                                    "{}({})-{}({})-{}({})-{}({}) | {:+.1}°",
                                    na, dmeas.a + 1, nb, dmeas.b + 1, nc, dmeas.c + 1, nd, dmeas.d + 1, dmeas.degrees
                                ));
                                ui.add_space(8.0);
                                if ui.button("Delete").clicked() {
                                    to_delete.push(id);
                                }
                            });
                            ui.separator();
                        }
                    }

                    if !to_delete.is_empty() {
                        measurements.dihedrals.retain(|d| !to_delete.contains(&d.id));
                    }
                    if req_clear {
                        measurements.preview_highlight = None;
                    } else if let Some(v) = req_set {
                        measurements.preview_highlight = Some(v);
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

                let mid = (a + b) * 0.5;

                if let Ok(screen_px) = cam.world_to_viewport(cam_xform, mid) {
                    let pos = egui::pos2(screen_px.x / scale, screen_px.y / scale);

                    let text = format!("{:.3}", pair.distance);
                    let galley = ctx.fonts(|f| {
                        f.layout_no_wrap(
                            text.clone(),
                            egui::FontId::proportional(16.0),
                            egui::Color32::WHITE,
                        )
                    });

                    let offset = egui::vec2(6.0, -6.0);
                    painter.galley(pos + offset, galley, egui::Color32::WHITE);
                }
            }
        });
}

fn pair_row(ui: &mut egui::Ui, pair: &mut MeasurePair, mol: &Molecule) -> bool {
    let name_a = mol.atoms.get(pair.a).map(|s| s.as_str()).unwrap_or("?");
    let name_b = mol.atoms.get(pair.b).map(|s| s.as_str()).unwrap_or("?");
    let mut delete_me = false;

    ui.horizontal(|ui| {
        let eye = if pair.visible { "Hide" } else { "Show" };
        if ui.button(eye).clicked() {
            pair.visible = !pair.visible;
        }

        ui.label(format!(
            "{}({}) — {}({})",
            name_a,
            pair.a + 1,
            name_b,
            pair.b + 1
        ));
        ui.separator();
        ui.monospace(format!("{:.3} Å", pair.distance));

        ui.add_space(8.0);

        ui.separator();

        if ui.button("Delete").clicked() {
            delete_me = true;
        }
    });

    delete_me
}

fn active_status_pill(ui: &mut egui::Ui, active: bool) {
    if active {
        egui::Frame::none()
            .fill(egui::Color32::from_rgb(40, 180, 120))
            .rounding(egui::Rounding::same(6))        // u8
            .inner_margin(egui::Margin::symmetric(8, 4)) // i8
            .show(ui, |ui| {
                ui.label(
                    egui::RichText::new("Active picking")
                        .strong()
                        .color(egui::Color32::WHITE),
                );
            });
    } else {
        ui.label("Inactive Picking");
    }
}

// Helper: whether the current preview highlight matches exactly this candidate set (ordered).
fn is_highlighted_for(current: &Option<Vec<usize>>, candidate: &[usize]) -> bool {
    match current {
        Some(v) => v.as_slice() == candidate,
        None => false,
    }
}

