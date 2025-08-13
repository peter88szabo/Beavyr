use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};

use crate::color_schemes::color_scheme_map;
use crate::molecule::{parse_xyz_angstrom, Molecule};
use crate::settings::{ColorScheme, LightingMode, MolSettings};

#[derive(Resource, Clone)]
pub struct XyzBuffer {
    pub text: String,
}

/// Same behavior as before, but:
/// - XYZ loader is collapsible
/// - Appearance split into 3 collapsibles: Geometry, Color & Materials, Lighting
/// - Emissive controls removed from UI
pub fn ui_panel(
    mut contexts: EguiContexts,
    mut settings: ResMut<MolSettings>,
    mut mol: ResMut<Molecule>,
    mut xyz_buf: ResMut<XyzBuffer>,
) {
    // bevy_egui 0.36: ctx_mut() returns Result; if it fails, skip this frame
    let Ok(ctx) = contexts.ctx_mut() else { return; };

    egui::SidePanel::right("controls")
        .default_width(380.0)
        .resizable(true)
        .show(&ctx, |ui| {
            // ===========================
            // 1) Structure (XYZ) — collapsible
            // ===========================
            ui.collapsing("Structure (XYZ)", |ui| {
                ui.label("Paste or edit XYZ. Input is assumed in Å (Angstrom).");

                ui.add(
                    egui::TextEdit::multiline(&mut xyz_buf.text)
                        .desired_rows(10)
                        .code_editor()
                        .lock_focus(true),
                );

                if ui.button("Load XYZ (Å)").clicked() {
                    let (_n, atoms, qxyz) = parse_xyz_angstrom(&xyz_buf.text);
                    mol.atoms = atoms;
                    mol.pos = qxyz
                        .into_iter()
                        .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
                        .collect();
                    mol.recompute_bonds(settings.bond_thresh_scale);
                    settings.dirty = true;
                }
            });

            ui.add_space(8.0);
            ui.separator();

            // ===========================
            // 2) Geometry — collapsible
            // ===========================
            ui.collapsing("Geometry Scaling", |ui| {
                let mut changed = false;
                changed |= ui
                    .add(egui::Slider::new(&mut settings.atom_scale, 0.3..=3.0).text("Atom scale"))
                    .changed();
                changed |= ui
                    .add(
                        egui::Slider::new(&mut settings.bond_thresh_scale, 0.8..=3.0)
                            .text("Bond cutoff × covalent radii"),
                    )
                    .changed();
                changed |= ui
                    .add(
                        egui::Slider::new(&mut settings.bond_radius_pct, 0.05..=1.0)
                            .text("Bond radius (% of smaller atom)"),
                    )
                    .changed();

                if changed {
                    // recompute bonds for threshold; radius is handled at spawn
                    mol.recompute_bonds(settings.bond_thresh_scale);
                    settings.dirty = true;
                }
            });

            ui.add_space(8.0);

            // ===========================
            // 3) Color & Materials — collapsible
            // ===========================
            ui.collapsing("Color & Materials", |ui| {
                // Background quick picks + picker
                ui.horizontal(|ui| {
                    ui.label("Background:");
                    if ui.button("Black").clicked() {
                        settings.bg_color = Color::srgb(0.0, 0.0, 0.0);
                        settings.dirty = true;
                    }
                    if ui.button("White").clicked() {
                        settings.bg_color = Color::srgb(1.0, 1.0, 1.0);
                        settings.dirty = true;
                    }
                    ui.label("or pick:");
                    let mut eg = color_to_egui(settings.bg_color);
                    if ui.color_edit_button_srgba(&mut eg).changed() {
                        settings.bg_color = egui_to_color(eg);
                        settings.dirty = true;
                    }
                });

                ui.add_space(6.0);
                ui.label("PBR material (atoms + bonds):");
                let mut mat_changed = false;
                mat_changed |= ui
                    .add(egui::Slider::new(&mut settings.metallic, 0.0..=1.0).text("Metallic"))
                    .changed();
                mat_changed |= ui
                    .add(egui::Slider::new(&mut settings.roughness, 0.02..=1.0).text("Roughness"))
                    .changed();
                mat_changed |= ui
                    .add(egui::Slider::new(&mut settings.reflectance, 0.0..=1.0).text("Reflectance"))
                    .changed();
                if mat_changed {
                    settings.dirty = true;
                }

                ui.add_space(6.0);
                ui.label("Colors:");
                ui.horizontal(|ui| {
                    ui.label("Scheme:");
                    egui::ComboBox::from_id_source("scheme_combo")
                        .selected_text(format!("{:?}", settings.scheme))
                        .show_ui(ui, |ui| {
                            for &sc in &[
                                ColorScheme::Custom,
                                ColorScheme::CPK,
                                ColorScheme::Jmol,
                                ColorScheme::VMD,
                                ColorScheme::Molden,
                                ColorScheme::Molden0,
                            ] {
                                if ui
                                    .selectable_label(settings.scheme == sc, format!("{:?}", sc))
                                    .clicked()
                                {
                                    settings.scheme = sc;
                                    if sc != ColorScheme::Custom {
                                        settings.element_colors = color_scheme_map(sc);
                                    }
                                    settings.dirty = true;
                                }
                            }
                        });
                });

                ui.collapsing("Per-element overrides (active in Custom)", |ui| {
                    let keys: Vec<String> = settings.element_colors.keys().cloned().collect();
                    for k in keys {
                        if let Some(old) = settings.element_colors.get(&k).copied() {
                            let mut col = color_to_egui(old);
                            ui.horizontal(|ui| {
                                ui.label(k.as_str());
                                if ui.color_edit_button_srgba(&mut col).changed() {
                                    settings
                                        .element_colors
                                        .insert(k.clone(), egui_to_color(col));
                                    settings.scheme = ColorScheme::Custom;
                                    settings.dirty = true;
                                }
                            });
                        }
                    }
                });
            });

            ui.add_space(8.0);

            // ===========================
            // 4) Lighting — collapsible
            // ===========================
            ui.collapsing("Lighting", |ui| {
                let mut changed = false;

                // Ambient
                changed |= ui.checkbox(&mut settings.use_ambient, "Ambient (omnidirectional fill)").changed();
                ui.horizontal(|ui| {
                    ui.label("Ambient color");
                    let mut ac = color_to_egui(settings.ambient_color);
                    if ui.color_edit_button_srgba(&mut ac).changed() {
                        settings.ambient_color = egui_to_color(ac);
                        settings.dirty = true;
                    }
                });
                changed |= ui
                    .add(egui::Slider::new(&mut settings.ambient_brightness, 0.0..=3000.0).text("Ambient brightness"))
                    .changed();

                ui.separator();

                // Lighting rig
                ui.horizontal(|ui| {
                    ui.label("Rig:");
                    egui::ComboBox::from_id_source("lighting_mode_combo")
                        .selected_text(match settings.lighting_mode {
                            LightingMode::SinglePoint => "Single point",
                            LightingMode::ThreePoint  => "Three-point",
                        })
                        .show_ui(ui, |ui| {
                            if ui.selectable_label(matches!(settings.lighting_mode, LightingMode::SinglePoint), "Single point").clicked() {
                                settings.lighting_mode = LightingMode::SinglePoint;
                                settings.dirty = true;
                            }
                            if ui.selectable_label(matches!(settings.lighting_mode, LightingMode::ThreePoint), "Three-point").clicked() {
                                settings.lighting_mode = LightingMode::ThreePoint;
                                settings.dirty = true;
                            }
                        });
                });

                // Key
                changed |= ui
                    .add(egui::Slider::new(&mut settings.light_intensity, 50_000.0..=8_000_000.0).text("Key intensity"))
                    .changed();
                changed |= ui
                    .add(egui::Slider::new(&mut settings.light_distance, 3.0..=40.0).text("Key distance"))
                    .changed();

                // Fill / Rim (only in three-point)
                if matches!(settings.lighting_mode, LightingMode::ThreePoint) {
                    ui.separator();
                    ui.label("Fill & Rim (three-point)");

                    changed |= ui
                        .add(egui::Slider::new(&mut settings.fill_intensity, 0.0..=6_000_000.0).text("Fill intensity"))
                        .changed();
                    changed |= ui
                        .add(egui::Slider::new(&mut settings.fill_distance, 3.0..=50.0).text("Fill distance"))
                        .changed();

                    changed |= ui
                        .add(egui::Slider::new(&mut settings.rim_intensity, 0.0..=6_000_000.0).text("Rim intensity"))
                        .changed();
                    changed |= ui
                        .add(egui::Slider::new(&mut settings.rim_distance, 3.0..=50.0).text("Rim distance"))
                        .changed();
                }

                if changed { settings.dirty = true; }
            });
        });
}

// --- color helpers (identical behavior) ---
fn color_to_egui(c: Color) -> egui::Color32 {
    let s = c.to_srgba();
    egui::Color32::from_rgba_premultiplied(
        (s.red * 255.0).round() as u8,
        (s.green * 255.0).round() as u8,
        (s.blue * 255.0).round() as u8,
        (s.alpha * 255.0).round() as u8,
    )
}
fn egui_to_color(c: egui::Color32) -> Color {
    Color::srgba(
        c.r() as f32 / 255.0,
        c.g() as f32 / 255.0,
        c.b() as f32 / 255.0,
        c.a() as f32 / 255.0,
    )
}

