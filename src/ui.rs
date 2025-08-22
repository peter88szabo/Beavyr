use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};

use crate::color_schemes::color_scheme_map;
use crate::molecule::{parse_xyz_angstrom, Molecule};
use crate::settings::{BondColorMode, ColorScheme, LightingMode, MolSettings};

// NEW: import inline export section + resources
use crate::export_ui;
use crate::export::{ExportRequestQueue, ExportUiState};

#[derive(Resource, Clone)]
pub struct XyzBuffer {
    pub text: String,
}

/// Same behavior as before +
/// - Two toggles in Structure block: show atom Type, show atom Index
/// - Labels drawn in the 3D view (overlay) using egui painter
/// - "Center molecule" button (and auto-center after load)
/// - Right-click near atom -> popup to use as rotation center
pub fn ui_panel(
    mut contexts: EguiContexts,
    mut settings: ResMut<MolSettings>,
    mut mol: ResMut<Molecule>,
    mut xyz_buf: ResMut<XyzBuffer>,
    // control orbit target
    mut cam: ResMut<crate::camera::OrbitCamera>,
    // need camera/projection to project atoms to screen for picking + labels
    q_cam: Query<(&Camera, &GlobalTransform), (With<Camera3d>, With<crate::scene::MainCamera>)>,
    windows: Query<&Window>,

    // NEW: export resources for the inline export section
    mut export_ui_state: ResMut<ExportUiState>,
    mut export_queue: ResMut<ExportRequestQueue>,
) {
    // bevy_egui 0.36: ctx_mut() returns Result; if it fails, skip this frame
    let Ok(ctx) = contexts.ctx_mut() else { return; };

    // -------------------------
    // Right-click picking: detect click in empty space (not over egui)
    // -------------------------
    let (latest_pos_opt, right_clicked) =
        ctx.input(|i| (i.pointer.latest_pos(), i.pointer.secondary_clicked()));
    let mouse_over_ui = ctx.is_pointer_over_area();
    if let (Some(pos_points), true) = (latest_pos_opt, right_clicked) {
        if !mouse_over_ui {
            if let (Ok((cam_comp, cam_xform)), Ok(win)) = (q_cam.single(), windows.single()) {
                // Egui uses logical points; world_to_viewport returns PHYSICAL px. Convert:
                let scale = win.scale_factor() as f32;
                let mouse_px = pos_points * scale;

                if let Some((hit_idx, hit_px_pos)) =
                    find_nearest_atom_screen_space(cam_comp, cam_xform, &mol.pos, mouse_px, 18.0)
                {
                    // Store the picked atom + popup location in egui memory (persist across frames)
                    let pick_id = egui::Id::new("picked_atom_idx");
                    let popup_id = egui::Id::new("rc_popup_open");
                    let popup_pos_id = egui::Id::new("rc_popup_pos_px");

                    ctx.data_mut(|d| {
                        d.insert_persisted(pick_id, hit_idx as i32);
                        d.insert_persisted(popup_id, true);
                        d.insert_persisted(popup_pos_id, hit_px_pos / scale); // store as logical points
                    });
                }
            }
        }
    }

    // If popup is open, draw a tiny context menu at stored mouse position
    {
        let popup_id = egui::Id::new("rc_popup_open");
        let popup_pos_id = egui::Id::new("rc_popup_pos_px");
        let pick_id = egui::Id::new("picked_atom_idx");

        let (mut open, pos_pts, picked) = ctx.data_mut(|d| {
            (
                d.get_persisted::<bool>(popup_id).unwrap_or(false),
                d.get_persisted::<egui::Pos2>(popup_pos_id)
                    .unwrap_or(egui::pos2(20.0, 20.0)),
                d.get_persisted::<i32>(pick_id).unwrap_or(-1),
            )
        });

        if open && picked >= 0 {
            egui::Area::new(egui::Id::new("rc_atom_popup"))
                .fixed_pos(pos_pts)
                .order(egui::Order::Foreground)
                .show(&ctx, |ui| {
                    egui::Frame::popup(&ctx.style()).show(ui, |ui| {
                        ui.vertical(|ui| {
                            ui.label(format!("Atom #{} (right-click)", picked));
                            if ui.button("Use as rotation center").clicked() {
                                let idx = picked as usize;
                                if idx < mol.pos.len() {
                                    cam.target = mol.pos[idx];
                                }
                                open = false; // close popup
                            }
                            if ui.button("Cancel").clicked() {
                                open = false;
                            }
                        });
                    });
                });

            // write back popup state
            ctx.data_mut(|d| d.insert_persisted(popup_id, open));
        }
    }

    // -------------------------
    // Side panel
    // -------------------------
    egui::SidePanel::right("controls")
        .default_width(350.0)
        .resizable(true)
        .show(&ctx, |ui| {
            // ===========================
            // 1) Structure (XYZ) — collapsible
            // ===========================
            ui.collapsing("Structure (XYZ)", |ui| {
                ui.label("Paste or edit XYZ. Input is assumed in Å (Angstrom).");

                // Two toggles: show type / show index
                let show_type_id = egui::Id::new("show_atom_type");
                let show_index_id = egui::Id::new("show_atom_index");

                // read current state
                let (mut show_type, mut show_index) = ctx.data_mut(|d| {
                    (
                        d.get_persisted::<bool>(show_type_id).unwrap_or(false),
                        d.get_persisted::<bool>(show_index_id).unwrap_or(false),
                    )
                });

                ui.horizontal(|ui| {
                    let sel = ui.visuals().selection.bg_fill;
                    let dim = ui.visuals().widgets.inactive.bg_fill;

                    // Toggle: show element symbol
                    if ui
                        .add(
                            egui::Button::new("Show Atom Type")
                                .fill(if show_type { sel } else { dim }),
                        )
                        .clicked()
                    {
                        show_type = !show_type;
                        ctx.data_mut(|d| d.insert_persisted(show_type_id, show_type));
                    }

                    // Toggle: show 1-based index
                    if ui
                        .add(
                            egui::Button::new("Show Atom Index")
                                .fill(if show_index { sel } else { dim }),
                        )
                        .clicked()
                    {
                        show_index = !show_index;
                        ctx.data_mut(|d| d.insert_persisted(show_index_id, show_index));
                    }
                });

                ui.add_space(8.0);

                ui.add(
                    egui::TextEdit::multiline(&mut xyz_buf.text)
                        //.desired_rows(14)
                        .desired_width(f32::INFINITY)
                        .code_editor()
                        .lock_focus(true),
                );

                ui.horizontal(|ui| {
                    if ui.button("Load XYZ (Å)").clicked() {
                        let (_n, atoms, qxyz) = parse_xyz_angstrom(&xyz_buf.text);
                        mol.atoms = atoms;
                        mol.pos = qxyz
                            .into_iter()
                            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
                            .collect();
                        mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                        // auto-center after loading
                        if let Some(c) = compute_centroid(&mol.pos) {
                            cam.target = c;
                        }
                        settings.dirty = true;
                    }

                    if ui.button("Center molecule").clicked() {
                        if let Some(c) = compute_centroid(&mol.pos) {
                            cam.target = c;
                        }
                    }
                });

                ui.small("Tip: right-click near an atom in the 3D view to set it as the rotation center.");
            });

            ui.add_space(8.0);
            ui.separator();

            // ===========================
            // 2) Geometry — collapsible
            // ===========================
            ui.collapsing("Atom & Bond Scaling", |ui| {
                let mut changed = false;
                changed |= ui
                    .add(egui::Slider::new(&mut settings.atom_scale, 0.1..=3.0).text("Atom scale"))
                    .changed();
                changed |= ui
                    .add(egui::Slider::new(&mut settings.atom_resolution, 0..=10).text("Atom resolution"))
                    .changed();

                ui.add_space(8.0);
                ui.separator();
                ui.add_space(8.0);

                changed |= ui
                    .add(
                        egui::Slider::new(&mut settings.bond_radius_pct, 0.05..=1.0)
                            .text("Bond radius (% of smaller atom)"),
                    )
                    .changed();

                changed |= ui
                    .add(
                        egui::Slider::new(&mut settings.bond_thresh_scale, 0.8..=3.0)
                            .text("Bond cutoff × covalent radii"),
                    )
                    .changed();

                ui.add_space(8.0);
                ui.separator();
                ui.add_space(8.0);

                changed |= ui
                    .add(
                        egui::Slider::new(&mut settings.hbond_cutoff, 1.2..=5.0)
                            .text("Hydrogen bond cutoff in Å"),
                    )
                    .changed();

                // ---- NEW: H-bond thickness & dash gaps ----
                changed |= ui
                    .add(
                        egui::Slider::new(&mut settings.hbond_thickness, 1.0..=80.0)
                            .text("H-bond thickness"),
                    )
                    .changed();
                changed |= ui
                    .add(
                        egui::Slider::new(&mut settings.hbond_gap_scale, 0.2..=5.0)
                            .text("H-bond dash gap scale"),
                    )
                    .changed();

                if changed {
                    mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
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
                    .add(egui::Slider::new(&mut settings.roughness, 0.005..=1.0).text("Roughness"))
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
                    egui::ComboBox::from_id_salt("scheme_combo")
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

                ui.add_space(8.0);
                ui.separator();

                // ---------------------------
                // NEW: Bond appearance block
                // ---------------------------
                ui.label("Bond appearance:");
                ui.horizontal(|ui| {
                    let sel = ui.visuals().selection.bg_fill;
                    let dim = ui.visuals().widgets.inactive.bg_fill;

                    let is_uniform = matches!(settings.bond_color_mode, BondColorMode::Uniform);
                    let is_split   = matches!(settings.bond_color_mode, BondColorMode::AtomSplit);

                    if ui
                        .add(egui::Button::new("Uniform bonds").fill(if is_uniform { sel } else { dim }))
                        .clicked()
                    {
                        settings.bond_color_mode = BondColorMode::Uniform;
                        settings.dirty = true;
                    }

                    if ui
                        .add(egui::Button::new("Atom-split bonds").fill(if is_split { sel } else { dim }))
                        .clicked()
                    {
                        settings.bond_color_mode = BondColorMode::AtomSplit;
                        settings.dirty = true;
                    }
                });

                // Uniform bond color picker (enabled for both modes; useful to preselect)
                ui.horizontal(|ui| {
                    ui.label("Uniform bond color:");
                    let mut bc = color_to_egui(settings.uniform_bond_color);
                    if ui.color_edit_button_srgba(&mut bc).changed() {
                        settings.uniform_bond_color = egui_to_color(bc);
                        settings.dirty = true;
                    }
                });

                // ---- NEW: Hydrogen-bond color ----
                ui.horizontal(|ui| {
                    ui.label("H-bond color:");
                    let mut hc = color_to_egui(settings.hbond_color);
                    if ui.color_edit_button_srgba(&mut hc).changed() {
                        settings.hbond_color = egui_to_color(hc);
                        settings.dirty = true;
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
                    egui::ComboBox::from_id_salt("lighting_mode_combo")
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

            // ===========================
            // 5) Export Image — collapsible (NEW)
            // ===========================
            ui.add_space(8.0);
            ui.separator();
            ui.collapsing("Export Image", |ui| {
                export_ui::export_section(ui, &mut export_ui_state, &mut export_queue);
            });
        });

    // -------------------------
    // Overlay: draw labels next to atoms (types / indices) if toggled
    // -------------------------
    let show_type = ctx.data_mut(|d| d.get_persisted::<bool>(egui::Id::new("show_atom_type")).unwrap_or(false));
    let show_index = ctx.data_mut(|d| d.get_persisted::<bool>(egui::Id::new("show_atom_index")).unwrap_or(false));

    if (show_type || show_index) && mol.pos.len() == mol.atoms.len() {
        if let (Ok((cam_comp, cam_xform)), Ok(win)) = (q_cam.single(), windows.single()) {
            // paint on a full-screen, non-interactive egui area
            egui::Area::new(egui::Id::new("atom_labels_overlay"))
                .movable(false)
                .interactable(false)
                .order(egui::Order::Tooltip) // on top but below popups
                .show(&ctx, |ui| {
                    let painter = ui.painter();
                    let scale = win.scale_factor() as f32;

                    for (i, (&p_world, sym)) in mol.pos.iter().zip(mol.atoms.iter()).enumerate() {
                        if let Ok(screen_px) = cam_comp.world_to_viewport(cam_xform, p_world) {
                            // Convert to egui logical points for painting
                            let pos = egui::pos2(screen_px.x / scale, screen_px.y / scale);

                            // Compose label text
                            let mut text = String::new();
                            if show_type {
                                text.push_str(sym);
                            }
                            if show_index {
                                if !text.is_empty() {
                                    text.push('(');
                                    text.push_str(&(i + 1).to_string());
                                    text.push(')');
                                } else {
                                    text.push_str(&(i + 1).to_string());
                                }
                            }
                            if text.is_empty() {
                                continue;
                            }

                            // Layout text (measure size) with egui fonts
                            let galley = ctx.fonts(|f| {
                                f.layout_no_wrap(
                                    text.clone(),
                                    egui::FontId::proportional(13.0),
                                    egui::Color32::WHITE,
                                )
                            });
                            // Draw the text galley
                            painter.galley(pos, galley, egui::Color32::WHITE);
                        }
                    }
                });
        }
    }
}

// ---- helpers ----

fn compute_centroid(points: &[Vec3]) -> Option<Vec3> {
    if points.is_empty() { return None; }
    let mut acc = Vec3::ZERO;
    for &p in points {
        acc += p;
    }
    Some(acc / (points.len() as f32))
}

/// Return (index, mouse_px_pos) of nearest atom within `radius_px` of `mouse_px`.
fn find_nearest_atom_screen_space(
    cam: &Camera,
    cam_xform: &GlobalTransform,
    atoms: &[Vec3],
    mouse_px: egui::Pos2,
    radius_px: f32,
) -> Option<(usize, egui::Pos2)> {
    // Convert mouse to Vec2 for distance math
    let mouse_v = Vec2::new(mouse_px.x, mouse_px.y);
    let mut best: Option<(usize, f32, Vec2)> = None;

    for (i, &p_world) in atoms.iter().enumerate() {
        // Bevy 0.16: world_to_viewport -> Result<Vec2, ViewportConversionError>
        if let Ok(screen_px) = cam.world_to_viewport(cam_xform, p_world) {
            let d = screen_px.distance(mouse_v);
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

    best.map(|(i, _d, s)| (i, egui::pos2(s.x, s.y)))
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

