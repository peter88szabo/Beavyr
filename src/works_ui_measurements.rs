// ==============================
// FILE: src/ui_measurements.rs
// ==============================
use bevy::prelude::*;
use bevy_egui::egui;

use crate::measurements::{MeasurePair, Measurements};
use crate::molecule::Molecule;
use crate::settings::MolSettings;

/// Call this from your existing right-side panel in `ui.rs`.
pub fn measurements_panel(
    ui: &mut egui::Ui,
    measurements: &mut Measurements,
    mol: &Molecule,
    settings: &MolSettings,
) {
    //ui.add_space(8.0);
    //ui.separator();
    ui.collapsing("Distance", |ui| {
        // Controls
        ui.horizontal(|ui| {
            ui.toggle_value(&mut measurements.is_active, "Activate picking");
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
        ui.checkbox(&mut measurements.show_labels, "Show labels (UI/in-scene later)");

        // Global color picker (same style as elsewhere in your UI)
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
        }); //end of Distance Sytle collapsible

        // Pairs list
        ui.separator();
        if measurements.pairs.is_empty() {
            ui.weak("No measurement pairs yet. Activate picking, then click two atoms.");
            return;
        }

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
        });
}

/// Row UI; returns true if user clicked delete.
fn pair_row(ui: &mut egui::Ui, pair: &mut MeasurePair, mol: &Molecule) -> bool {
    let name_a = mol.atoms.get(pair.a).map(|s| s.as_str()).unwrap_or("?");
    let name_b = mol.atoms.get(pair.b).map(|s| s.as_str()).unwrap_or("?");
    let mut delete_me = false;

    ui.horizontal(|ui| {
        // Show/Hide
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
        // Label toggle (reserved for future in-scene labels)
        ui.toggle_value(&mut pair.label_on, "Label");

        // Delete
        if ui.button("Delete").clicked() {
            delete_me = true;
        }
    });

    delete_me
}
