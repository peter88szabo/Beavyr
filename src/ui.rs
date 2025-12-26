// src/ui.rs
use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts};
use std::collections::{BTreeSet, HashSet};
use std::fs;
use std::path::PathBuf;

use crate::color_schemes::{color_scheme_map, ELEMENT_SYMBOLS};
use crate::events::MoleculeChanged;
use crate::molecule::{parse_xyz_first_frame_angstrom, Molecule};
use crate::settings::{BondColorMode, ColorScheme, LightingMode, MolSettings, RepresentationMode};

// Export UI + resources
use crate::export_ui;
use crate::export::{ExportRequestQueue, ExportUiState};

// Measurements UI + resource
use crate::measurements::Measurements;
use crate::ui_measurements;
use crate::trajectory;

// NEW: shared picking helper (egui-friendly shim)
use crate::picking::screen::find_nearest_atom_screen_space_egui;

const ELEMENT_FILTER_S_BLOCK: [&str; 14] = [
    "H", "He", "Li", "Be", "Na", "Mg", "K", "Ca", "Rb", "Sr", "Cs", "Ba", "Fr", "Ra",
];
const ELEMENT_FILTER_D_BLOCK: [&str; 29] = [
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
];
const ELEMENT_FILTER_P_BLOCK: [&str; 29] = [
    "B", "C", "N", "O", "F", "Ne",
    "Al", "Si", "P", "S", "Cl", "Ar",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Tl", "Pb", "Bi", "Po", "At",
];
const ELEMENT_FILTER_F_BLOCK: [&str; 30] = [
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
];

#[derive(Resource, Clone)]
pub struct XyzBuffer {
    pub text: String,
    pub last_dir: Option<PathBuf>,
    pub warning: Option<String>,
}

/// Right-side control panel + overlays.
///
/// Key features:
/// - Show atom type / index labels in 3D overlay
/// - Right-click pick atom → popup (set camera target)
/// - "Structure (XYZ)" block:
///     * View mode: live read-only XYZ from Molecule.pos
///     * Edit mode: multiline text editor; Apply → parse → write Molecule + emit MoleculeChanged
///     * Copy XYZ → copies current live geometry
pub fn ui_panel(
    mut contexts: EguiContexts,
    mut settings: ResMut<MolSettings>,
    mut mol: ResMut<Molecule>,
    mut xyz_buf: ResMut<XyzBuffer>,
    mut ev_changed: EventWriter<MoleculeChanged>,
    // orbit target control
    mut cam: ResMut<crate::camera::OrbitCamera>,
    // camera for picking + label projection
    q_cam: Query<(&Camera, &GlobalTransform), (With<Camera3d>, With<crate::scene::MainCamera>)>,
    windows: Query<&Window>,

    // export section
    mut export_ui_state: ResMut<ExportUiState>,
    mut export_queue: ResMut<ExportRequestQueue>,

    // measurements panel
    mut measurements: ResMut<Measurements>,
    // trajectory
    mut traj: ResMut<trajectory::TrajectoryState>,
) {
    // bevy_egui 0.36: ctx_mut() returns Result; if it fails, skip this frame
    let Ok(ctx) = contexts.ctx_mut() else { return; };
    ctx.style_mut(|style| {
        style.text_styles = [
            (egui::TextStyle::Heading, egui::FontId::proportional(22.0)),
            (egui::TextStyle::Body, egui::FontId::proportional(16.0)),
            (egui::TextStyle::Monospace, egui::FontId::monospace(15.0)),
            (egui::TextStyle::Button, egui::FontId::proportional(16.0)),
            (egui::TextStyle::Small, egui::FontId::proportional(13.0)),
        ]
        .into();
    });

    // -------------------------
    // Right-click picking (not over egui)
    // -------------------------
    let (latest_pos_opt, right_clicked) =
        ctx.input(|i| (i.pointer.latest_pos(), i.pointer.secondary_clicked()));
    let mouse_over_ui = ctx.is_pointer_over_area();
    if let (Some(pos_points), true) = (latest_pos_opt, right_clicked) {
        if !mouse_over_ui {
            if let (Ok((cam_comp, cam_xform)), Ok(win)) = (q_cam.single(), windows.single()) {
                // Convert egui logical points → PHYSICAL px
                let scale = win.scale_factor() as f32;
                let mouse_px = egui::pos2(pos_points.x * scale, pos_points.y * scale);

                let (picked_idx, popup_pos) = if let Some((hit_idx, hit_px_pos)) =
                    find_nearest_atom_screen_space_egui(cam_comp, cam_xform, &mol.pos, mouse_px, 18.0)
                {
                    (hit_idx as i32, egui::pos2(hit_px_pos.x / scale, hit_px_pos.y / scale))
                } else {
                    (-1, pos_points)
                };

                // Persist picked atom + popup position (store in logical points)
                let pick_id = egui::Id::new("picked_atom_idx");
                let popup_id = egui::Id::new("rc_popup_open");
                let popup_pos_id = egui::Id::new("rc_popup_pos_px");

                ctx.data_mut(|d| {
                    d.insert_persisted(pick_id, picked_idx);
                    d.insert_persisted(popup_id, true);
                    d.insert_persisted(popup_pos_id, popup_pos);
                });
            }
        }
    }

    // Right-click popup at stored position
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

        if open {
            egui::Area::new(egui::Id::new("rc_atom_popup"))
                .fixed_pos(pos_pts)
                .order(egui::Order::Foreground)
                .show(&ctx, |ui| {
                    egui::Frame::popup(&ctx.style()).show(ui, |ui| {
                        if picked >= 0 {
                            ui.set_min_width(220.0);
                        }
                        ui.vertical(|ui| {
                            ui.label("Background:");
                            ui.horizontal(|ui| {
                                if ui.button("Black").clicked() {
                                    settings.bg_color = Color::srgb(0.0, 0.0, 0.0);
                                    settings.dirty = true;
                                }
                                if ui.button("White").clicked() {
                                    settings.bg_color = Color::srgb(1.0, 1.0, 1.0);
                                    settings.dirty = true;
                                }
                            });
                            if ui.button("Center molecule").clicked() {
                                if let Some(c) = compute_centroid(&mol.pos) {
                                    cam.target = c;
                                }
                                open = false;
                            }
                            if picked >= 0 {
                                ui.separator();
                                ui.label(format!("Atom #{} (right-click)", picked));
                                if ui.button("Use as rotation center").clicked() {
                                    let idx = picked as usize;
                                    if idx < mol.pos.len() {
                                        cam.target = mol.pos[idx];
                                    }
                                    open = false;
                                }
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
    // Right side panel
    // -------------------------
    egui::SidePanel::right("controls")
        .default_width(320.0)
        .min_width(240.0)
        .max_width(520.0)
        .resizable(true)
        .show(&ctx, |ui| {
            // ===========================
            // 1) Structure (XYZ)
            // ===========================
            ui.add_space(8.0);
            ui.collapsing("Structure (XYZ)", |ui| {
                let show_type_id = egui::Id::new("show_atom_type");
                let show_index_id = egui::Id::new("show_atom_index");
                let edit_mode_id = egui::Id::new("xyz_edit_mode");

                // load flags
                let (mut show_type, mut show_index, mut edit_mode) = ctx.data_mut(|d| {
                    (
                        d.get_persisted::<bool>(show_type_id).unwrap_or(false),
                        d.get_persisted::<bool>(show_index_id).unwrap_or(false),
                        d.get_persisted::<bool>(edit_mode_id).unwrap_or(false),
                    )
                });

                ui.horizontal(|ui| {
                    let sel = ui.visuals().selection.bg_fill;
                    let dim = ui.visuals().widgets.inactive.bg_fill;

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
                ui.horizontal(|ui| {
                    let sel = ui.visuals().selection.bg_fill;
                    let dim = ui.visuals().widgets.inactive.bg_fill;

                    if ui.button("Load XYZ file").clicked() {
                        let mut dlg = rfd::FileDialog::new().add_filter("XYZ", &["xyz"]);
                        if let Some(dir) = &xyz_buf.last_dir {
                            dlg = dlg.set_directory(dir);
                        }
                        if let Some(path) = dlg.pick_file() {
                            xyz_buf.last_dir = path.parent().map(|p| p.to_path_buf());
                            if let Ok(text) = fs::read_to_string(&path) {
                                xyz_buf.text = text;
                                xyz_buf.warning = apply_xyz_text(
                                    &xyz_buf.text,
                                    &mut mol,
                                    &mut cam,
                                    &mut ev_changed,
                                );
                            }
                        }
                    }

                    if ui
                        .add(
                            egui::Button::new(if edit_mode { "View mode" } else { "Edit as text" })
                                .fill(if edit_mode { sel } else { dim }),
                        )
                        .on_hover_text("Toggle free-form XYZ editing")
                        .clicked()
                    {
                        edit_mode = !edit_mode;
                        if edit_mode {
                            // entering edit mode: prefill from live geometry
                            xyz_buf.text = format_xyz_from_molecule(&mol);
                        }
                        ctx.data_mut(|d| d.insert_persisted(edit_mode_id, edit_mode));
                    }

                    if ui.button("Copy XYZ").clicked() {
                        let live = format_xyz_from_molecule(&mol);
                        ctx.copy_text(live);
                    }
                });
                ui.add_space(8.0);

                if let Some(warn) = &xyz_buf.warning {
                    ui.colored_label(
                        egui::Color32::from_rgb(255, 170, 0),
                        format!("Warning: {}", warn),
                    );
                }

                if edit_mode {
                    // ---- EDIT MODE ----
                    ui.label("Paste or edit XYZ. Input is assumed in Å (Angstrom).");
                    egui::ScrollArea::vertical()
                        .max_height(260.0)
                        .auto_shrink([false; 2])
                        .show(ui, |ui| {
                            ui.add(
                                egui::TextEdit::multiline(&mut xyz_buf.text)
                                    .desired_width(f32::INFINITY)
                                    .code_editor()
                                    .lock_focus(true),
                            );
                        });

                    ui.horizontal(|ui| {
                        if ui.button("Apply (Parse XYZ)").clicked() {
                            xyz_buf.warning = apply_xyz_text(
                                &xyz_buf.text,
                                &mut mol,
                                &mut cam,
                                &mut ev_changed,
                            );
                        }

                        if ui.button("Clear text").clicked() {
                            xyz_buf.text.clear();
                        }

                        if ui.button("Revert from live").clicked() {
                            xyz_buf.text = format_xyz_from_molecule(&mol);
                        }
                    });

                    ui.add_space(8.0);
                    ui.weak("Tip: Use ‘Copy XYZ’ to copy the current live geometry without leaving edit mode.");
                } else {
                    // ---- VIEW MODE (LIVE XYZ) ----
                    ui.label("Live XYZ (read-only):");
                    egui::ScrollArea::vertical()
                        .max_height(260.0)
                        .auto_shrink([false; 2])
                        .show(ui, |ui| {
                            ui.style_mut().override_text_style = Some(egui::TextStyle::Monospace);
                            for (i, (sym, p)) in mol.atoms.iter().zip(mol.pos.iter()).enumerate() {
                                let label = format!(
                                    "{:>3} {:<2} {:>7.3} {:>7.3} {:>7.3}",
                                    i + 1,
                                    sym,
                                    p.x,
                                    p.y,
                                    p.z
                                );
                                ui.monospace(label);
                            }
                        });

                    ui.add_space(6.0);
                    ui.horizontal(|ui| {
                        if ui.button("Center molecule").clicked() {
                            if let Some(c) = compute_centroid(&mol.pos) {
                                cam.target = c;
                            }
                        }
                        let hbond_label = if settings.show_hbonds {
                            "H-bonds Off"
                        } else {
                            "H-bonds On"
                        };
                        if ui.button(hbond_label).clicked() {
                            settings.show_hbonds = !settings.show_hbonds;
                        }
                    });
                }
            });

            ui.add_space(8.0);
            ui.separator();

            // ===========================
            // 1b) Trajectory
            // ===========================
            ui.add_space(8.0);
            ui.collapsing("Trajectory", |ui| {
                if ui.button("Load trajectory XYZ").clicked() {
                    let mut dlg = rfd::FileDialog::new().add_filter("XYZ", &["xyz"]);
                    if let Some(dir) = &traj.last_dir {
                        dlg = dlg.set_directory(dir);
                    }
                    if let Some(path) = dlg.pick_file() {
                        traj.last_dir = path.parent().map(|p| p.to_path_buf());
                        xyz_buf.warning = None;
                        if let Ok(text) = fs::read_to_string(&path) {
                            match trajectory::parse_multi_xyz(&text) {
                                Ok((atoms, frames)) => {
                                    traj.atoms = atoms;
                                    traj.frames = frames;
                                    traj.current_frame = 0;
                                    traj.playing = false;
                                    traj.accum = 0.0;
                                    traj.last_applied = None;
                                    traj.overlay_dirty = true;
                                    trajectory::apply_current_frame(
                                        &mut traj,
                                        &mut mol,
                                        &mut settings,
                                        &mut ev_changed,
                                        true,
                                        Some(&mut cam),
                                    );
                                }
                                Err(_) => {}
                            }
                        }
                    }
                }

                if traj.frames.is_empty() {
                    ui.weak("No trajectory loaded.");
                } else {
                    let total = traj.frames.len();
                    ui.label(format!("Frame {}/{}", traj.current_frame + 1, total));

                    let mut frame_display = (traj.current_frame + 1) as i32;
                    if ui
                        .add(egui::Slider::new(&mut frame_display, 1..=total as i32).text("Frame"))
                        .changed()
                    {
                        traj.current_frame = (frame_display - 1) as usize;
                        traj.playing = false;
                        trajectory::apply_current_frame(
                            &mut traj,
                            &mut mol,
                            &mut settings,
                            &mut ev_changed,
                            false,
                            None,
                        );
                    }

                    ui.horizontal(|ui| {
                        if ui.button("|◀").clicked() {
                            traj.current_frame = 0;
                            traj.playing = false;
                            trajectory::apply_current_frame(
                                &mut traj,
                                &mut mol,
                                &mut settings,
                                &mut ev_changed,
                                false,
                                None,
                            );
                        }
                        if ui.button("Play ◀").clicked() {
                            traj.direction = -1;
                            traj.playing = true;
                        }
                        if ui.button("Pause").clicked() {
                            traj.playing = false;
                        }
                        if ui.button("Play ▶").clicked() {
                            traj.direction = 1;
                            traj.playing = true;
                        }
                        if ui.button("▶|").clicked() {
                            traj.current_frame = total - 1;
                            traj.playing = false;
                            trajectory::apply_current_frame(
                                &mut traj,
                                &mut mol,
                                &mut settings,
                                &mut ev_changed,
                                false,
                                None,
                            );
                        }
                    });

                    ui.horizontal(|ui| {
                        if ui.button("Step ◀").clicked() {
                            if traj.current_frame > 0 {
                                traj.current_frame -= 1;
                                traj.playing = false;
                                trajectory::apply_current_frame(
                                    &mut traj,
                                    &mut mol,
                                    &mut settings,
                                    &mut ev_changed,
                                    false,
                                    None,
                                );
                            }
                        }
                        if ui.button("Step ▶").clicked() {
                            if traj.current_frame + 1 < total {
                                traj.current_frame += 1;
                                traj.playing = false;
                                trajectory::apply_current_frame(
                                    &mut traj,
                                    &mut mol,
                                    &mut settings,
                                    &mut ev_changed,
                                    false,
                                    None,
                                );
                            }
                        }
                    });

                    ui.horizontal(|ui| {
                        ui.label("Speed (fps)");
                        ui.add(egui::Slider::new(&mut traj.fps, 1.0..=100.0).show_value(true));
                    });

                    ui.horizontal(|ui| {
                        let sel = ui.visuals().selection.bg_fill;
                        let dim = ui.visuals().widgets.inactive.bg_fill;

                        let is_loop = matches!(traj.mode, trajectory::PlaybackMode::Loop);
                        let is_once = matches!(traj.mode, trajectory::PlaybackMode::Once);
                        let is_ping = matches!(traj.mode, trajectory::PlaybackMode::PingPong);

                        if ui
                            .add(egui::Button::new("Loop").fill(if is_loop { sel } else { dim }))
                            .clicked()
                        {
                            traj.mode = trajectory::PlaybackMode::Loop;
                        }
                        if ui
                            .add(egui::Button::new("Once").fill(if is_once { sel } else { dim }))
                            .clicked()
                        {
                            traj.mode = trajectory::PlaybackMode::Once;
                        }
                        if ui
                            .add(egui::Button::new("Ping-pong").fill(if is_ping { sel } else { dim }))
                            .clicked()
                        {
                            traj.mode = trajectory::PlaybackMode::PingPong;
                        }
                    });

                    ui.add_space(8.0);
                    ui.separator();
                    ui.add_space(8.0);

                    ui.horizontal(|ui| {
                        let sel = ui.visuals().selection.bg_fill;
                        let dim = ui.visuals().widgets.inactive.bg_fill;
                        ui.add_space(4.0);
                        if ui
                            .add(egui::Button::new("Overlay").fill(if traj.overlay_enabled { sel } else { dim }))
                            .clicked()
                        {
                            traj.overlay_enabled = !traj.overlay_enabled;
                            traj.overlay_dirty = true;
                        }
                    });

                    ui.horizontal(|ui| {
                        ui.label("Full overlay stride");
                        let mut stride_val = traj.overlay_full_stride as i32;
                        if ui
                            .add(
                                egui::Slider::new(&mut stride_val, 0..=500)
                                    .show_value(true),
                            )
                            .changed()
                        {
                            traj.overlay_full_stride = stride_val.max(0) as usize;
                            traj.overlay_dirty = true;
                        }
                    });

                    ui.horizontal(|ui| {
                        ui.label("Ghost overlay stride");
                        let mut stride_val = traj.overlay_ghost_stride as i32;
                        if ui
                            .add(
                                egui::Slider::new(&mut stride_val, 0..=500)
                                    .show_value(true),
                            )
                            .changed()
                        {
                            traj.overlay_ghost_stride = stride_val.max(0) as usize;
                            traj.overlay_dirty = true;
                        }
                    });

                    ui.horizontal(|ui| {
                        let sel = ui.visuals().selection.bg_fill;
                        let dim = ui.visuals().widgets.inactive.bg_fill;
                        if ui
                            .add(
                                egui::Button::new("Full-last")
                                    .fill(if traj.overlay_full_last { sel } else { dim }),
                            )
                            .clicked()
                        {
                            traj.overlay_full_last = !traj.overlay_full_last;
                            traj.overlay_dirty = true;
                        }
                    });

                    ui.horizontal(|ui| {
                        ui.label("Ghost transparency");
                        if ui
                            .add(egui::Slider::new(&mut traj.ghost_alpha, 0.05..=0.8))
                            .changed()
                        {
                            traj.overlay_dirty = true;
                        }
                    });
                }
            });
            // ===========================
            // 2) Measurements
            // ===========================
            ui.add_space(8.0);
            ui.separator();
            ui.collapsing("Measurements", |ui| {
                ui_measurements::measurements_panel(
                    ui,
                    &mut measurements,
                    &mol,
                    &settings,
                );
            });

            // ===========================
            // 3) Appearance
            // ===========================
            ui.add_space(8.0);
            ui.separator();
            ui.collapsing("Appearance", |ui| {
                ui.collapsing("Representation", |ui| {
                    let mut changed = false;

                    changed |= ui
                        .radio_value(
                            &mut settings.representation,
                            RepresentationMode::BallAndStick,
                            "Ball & stick (best <2k atoms)",
                        )
                        .clicked();
                    changed |= ui
                        .radio_value(
                            &mut settings.representation,
                            RepresentationMode::SpaceFilling,
                            "Space-filling (CPK)",
                        )
                        .clicked();
                    changed |= ui
                        .radio_value(
                            &mut settings.representation,
                            RepresentationMode::SticksRounded,
                            "Sticks only (rounded joints)",
                        )
                        .clicked();
                    changed |= ui
                        .radio_value(
                            &mut settings.representation,
                            RepresentationMode::LowResBallsAndLines,
                            "Low-res balls + lines",
                        )
                        .clicked();
                    changed |= ui
                        .radio_value(
                            &mut settings.representation,
                            RepresentationMode::LinesOnly,
                            "Lines only",
                        )
                        .clicked();
                    changed |= ui
                        .radio_value(
                            &mut settings.representation,
                            RepresentationMode::BackboneTrace,
                            "Backbone/trace (non-H lines)",
                        )
                        .clicked();

                    if matches!(settings.representation, RepresentationMode::LowResBallsAndLines) {
                        ui.add_space(6.0);
                        if ui
                            .add(
                                egui::Slider::new(&mut settings.low_res_atom_resolution, 0..=6)
                                    .text("Low-res atom resolution"),
                            )
                            .changed()
                        {
                            changed = true;
                        }
                        if ui
                            .add(
                                egui::Slider::new(&mut settings.low_res_atom_scale, 0.2..=1.5)
                                    .text("Low-res atom scale"),
                            )
                            .changed()
                        {
                            changed = true;
                        }
                    }

                    if matches!(settings.representation, RepresentationMode::SpaceFilling) {
                        ui.add_space(6.0);
                        if ui
                            .add(
                                egui::Slider::new(&mut settings.cpk_atom_resolution, 0..=10)
                                    .text("CPK atom resolution"),
                            )
                            .changed()
                        {
                            changed = true;
                        }
                        if ui
                            .add(
                                egui::Slider::new(&mut settings.cpk_atom_scale, 0.6..=2.6)
                                    .text("CPK atom scale"),
                            )
                            .changed()
                        {
                            changed = true;
                        }
                    }

                    if matches!(settings.representation, RepresentationMode::SticksRounded) {
                        ui.add_space(6.0);
                        if ui
                            .add(
                                egui::Slider::new(&mut settings.stick_radius, 0.02..=0.40)
                                    .text("Stick thickness"),
                            )
                            .changed()
                        {
                            changed = true;
                        }
                    }

                    if matches!(
                        settings.representation,
                        RepresentationMode::LowResBallsAndLines
                            | RepresentationMode::LinesOnly
                            | RepresentationMode::BackboneTrace
                    ) {
                        ui.add_space(6.0);
                        if ui
                            .add(
                                egui::Slider::new(&mut settings.line_bond_thickness, 1.0..=50.0)
                                    .text("Line thickness"),
                            )
                            .changed()
                        {
                            changed = true;
                        }
                    }

                    if matches!(settings.representation, RepresentationMode::BackboneTrace) {
                        ui.add_space(4.0);
                        if ui
                            .checkbox(&mut settings.trace_show_atoms, "Show trace atoms")
                            .changed()
                        {
                            changed = true;
                        }
                        if settings.trace_show_atoms {
                            if ui
                                .add(
                                    egui::Slider::new(&mut settings.low_res_atom_resolution, 0..=6)
                                        .text("Trace atom resolution"),
                                )
                                .changed()
                            {
                                changed = true;
                            }
                            if ui
                                .add(
                                    egui::Slider::new(&mut settings.trace_atom_scale, 0.2..=1.5)
                                        .text("Trace atom scale"),
                                )
                                .changed()
                            {
                                changed = true;
                            }
                        }
                    }

                    if changed {
                        settings.dirty = true;
                    }
                });

                ui.add_space(8.0);

                ui.collapsing("Element Filter", |ui| {
                    let mut elements: Vec<String> = Vec::new();
                    elements.extend(ELEMENT_FILTER_S_BLOCK.iter().map(|s| s.to_string()));
                    elements.extend(ELEMENT_FILTER_D_BLOCK.iter().map(|s| s.to_string()));
                    elements.extend(ELEMENT_FILTER_P_BLOCK.iter().map(|s| s.to_string()));
                    elements.extend(ELEMENT_FILTER_F_BLOCK.iter().map(|s| s.to_string()));

                    let mut default_set: HashSet<&str> = HashSet::new();
                    default_set.extend(ELEMENT_FILTER_S_BLOCK);
                    default_set.extend(ELEMENT_FILTER_D_BLOCK);
                    default_set.extend(ELEMENT_FILTER_P_BLOCK);
                    default_set.extend(ELEMENT_FILTER_F_BLOCK);

                    let mut extras: BTreeSet<String> = BTreeSet::new();
                    for sym in &mol.atoms {
                        if !default_set.contains(sym.as_str()) {
                            extras.insert(sym.clone());
                        }
                    }
                    elements.extend(extras.iter().cloned());

                    if elements.is_empty() {
                        ui.weak("No atoms loaded.");
                        return;
                    }

                    let mut changed = false;
                    ui.horizontal(|ui| {
                        if ui.button("Show all").clicked() {
                            for sym in &elements {
                                settings.element_visibility.insert(sym.clone(), true);
                            }
                            changed = true;
                        }
                        if ui.button("Hide all").clicked() {
                            for sym in &elements {
                                settings.element_visibility.insert(sym.clone(), false);
                            }
                            changed = true;
                        }
                    });

                    ui.add_space(4.0);
                    ui.spacing_mut().item_spacing = egui::vec2(6.0, 4.0);
                    egui::ScrollArea::vertical()
                        .max_height(220.0)
                        .show(ui, |ui| {
                            let render_block = |ui: &mut egui::Ui,
                                                block: &[&str],
                                                settings: &mut MolSettings,
                                                changed: &mut bool| {
                                let cols = 4;
                                let rows = (block.len() + cols - 1) / cols;
                                for r in 0..rows {
                                    ui.horizontal(|ui| {
                                        for c in 0..cols {
                                            let idx = r * cols + c;
                                            if let Some(sym) = block.get(idx) {
                                                let entry = settings
                                                    .element_visibility
                                                    .entry((*sym).to_string())
                                                    .or_insert(true);
                                                if ui.checkbox(entry, *sym).changed() {
                                                    *changed = true;
                                                }
                                            } else {
                                                ui.add_space(52.0);
                                            }
                                        }
                                    });
                                }
                            };

                            ui.label("Hydrogen");
                            let h_block = ["H"];
                            render_block(ui, &h_block, &mut settings, &mut changed);
                            ui.separator();
                            ui.label("P-block");
                            render_block(ui, &ELEMENT_FILTER_P_BLOCK, &mut settings, &mut changed);
                            ui.separator();
                            ui.label("S-block");
                            render_block(ui, &ELEMENT_FILTER_S_BLOCK, &mut settings, &mut changed);
                            ui.separator();
                            ui.label("D-block");
                            render_block(ui, &ELEMENT_FILTER_D_BLOCK, &mut settings, &mut changed);
                            ui.separator();
                            ui.label("F-block");
                            render_block(ui, &ELEMENT_FILTER_F_BLOCK, &mut settings, &mut changed);

                            if !extras.is_empty() {
                                ui.separator();
                                ui.label("Other");
                                let other: Vec<&str> = extras.iter().map(|s| s.as_str()).collect();
                                render_block(ui, &other, &mut settings, &mut changed);
                            }
                        });

                    if changed {
                        settings.dirty = true;
                    }
                });

                ui.add_space(8.0);

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

                    // H-bond style
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
                        .add(egui::Slider::new(&mut settings.metallic, 0.0..=0.5).text("Metallic"))
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
                        ui.spacing_mut().item_spacing = egui::vec2(6.0, 4.0);
                        egui::ScrollArea::vertical()
                            .max_height(240.0)
                            .show(ui, |ui| {
                                let per_row = 6;
                                for chunk in ELEMENT_SYMBOLS.chunks(per_row) {
                                    ui.horizontal(|ui| {
                                        for &sym in chunk {
                                            if let Some(old) = settings.element_colors.get(sym).copied() {
                                                let key = sym.to_string();
                                                let mut col = color_to_egui(old);
                                                ui.allocate_ui_with_layout(
                                                    egui::vec2(70.0, 0.0),
                                                    egui::Layout::left_to_right(egui::Align::Center),
                                                    |ui| {
                                                        ui.add_sized([22.0, 0.0], egui::Label::new(sym));
                                                        if ui.color_edit_button_srgba(&mut col).changed() {
                                                            settings
                                                                .element_colors
                                                                .insert(key.clone(), egui_to_color(col));
                                                            settings.scheme = ColorScheme::Custom;
                                                            settings.dirty = true;
                                                        }
                                                    },
                                                );
                                            }
                                        }
                                    });
                                }
                            });
                    });

                    ui.add_space(8.0);
                    ui.separator();

                    // Bond appearance
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

                    ui.horizontal(|ui| {
                        ui.label("Uniform bond color:");
                        let mut bc = color_to_egui(settings.uniform_bond_color);
                        if ui.color_edit_button_srgba(&mut bc).changed() {
                            settings.uniform_bond_color = egui_to_color(bc);
                            settings.dirty = true;
                        }
                    });

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
                            .add(egui::Slider::new(&mut settings.fill_intensity, 0.0..=10_000_000.0).text("Fill intensity"))
                            .changed();
                        changed |= ui
                            .add(egui::Slider::new(&mut settings.fill_distance, 3.0..=50.0).text("Fill distance"))
                            .changed();

                        changed |= ui
                            .add(egui::Slider::new(&mut settings.rim_intensity, 0.0..=10_000_000.0).text("Rim intensity"))
                            .changed();
                        changed |= ui
                            .add(egui::Slider::new(&mut settings.rim_distance, 3.0..=50.0).text("Rim distance"))
                            .changed();
                    }

                    if changed { settings.dirty = true; }
                });
            });

            // ===========================
            // 5) Export Image
            // ===========================
            ui.add_space(8.0);
            ui.separator();
            ui.collapsing("Export Image", |ui| {
                export_ui::export_section(ui, &mut export_ui_state, &mut export_queue);
            });
        });

    // -------------------------
    // Overlay: atom labels (type/index)
    // -------------------------
    let show_type = ctx.data_mut(|d| d.get_persisted::<bool>(egui::Id::new("show_atom_type")).unwrap_or(false));
    let show_index = ctx.data_mut(|d| d.get_persisted::<bool>(egui::Id::new("show_atom_index")).unwrap_or(false));

    if (show_type || show_index) && mol.pos.len() == mol.atoms.len() {
        if let (Ok((cam_comp, cam_xform)), Ok(win)) = (q_cam.single(), windows.single()) {
            egui::Area::new(egui::Id::new("atom_labels_overlay"))
                .movable(false)
                .interactable(false)
                .order(egui::Order::Tooltip)
                .show(&ctx, |ui| {
                    let painter = ui.painter();
                    let scale = win.scale_factor() as f32;

                    for (i, (&p_world, sym)) in mol.pos.iter().zip(mol.atoms.iter()).enumerate() {
                        if let Ok(screen_px) = cam_comp.world_to_viewport(cam_xform, p_world) {
                            let pos = egui::pos2(screen_px.x / scale, screen_px.y / scale);

                            let mut text = String::new();
                            if show_type { text.push_str(sym); }
                            if show_index {
                                if !text.is_empty() {
                                    text.push('(');
                                    text.push_str(&(i + 1).to_string());
                                    text.push(')');
                                } else {
                                    text.push_str(&(i + 1).to_string());
                                }
                            }
                            if text.is_empty() { continue; }

                            let galley = ctx.fonts(|f| {
                                f.layout_no_wrap(
                                    text.clone(),
                                    egui::FontId::proportional(13.0),
                                    egui::Color32::WHITE,
                                )
                            });
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

fn apply_xyz_text(
    text: &str,
    mol: &mut Molecule,
    cam: &mut crate::camera::OrbitCamera,
    ev_changed: &mut EventWriter<MoleculeChanged>,
) -> Option<String> {
    let (atoms, qxyz, extra_frames) = parse_xyz_first_frame_angstrom(text);
    mol.atoms = atoms;
    mol.set_pos(
        qxyz
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect(),
    );

    // Emit change event; parsing a new XYZ likely changes topology
    ev_changed.write(MoleculeChanged::parse_xyz(true));

    // Auto-center after load
    if let Some(c) = compute_centroid(&mol.pos) {
        cam.target = c;
    }

    if extra_frames {
        Some("Multiple XYZ frames detected; loaded only the first frame.".to_string())
    } else {
        None
    }
}

/// Format current `Molecule` as simple XYZ lines: "Sym  x  y  z"
fn format_xyz_from_molecule(mol: &Molecule) -> String {
    let mut out = String::new();
    for (sym, p) in mol.atoms.iter().zip(mol.pos.iter()) {
        out.push_str(&format!(
            "{:<2} {:>14.8} {:>14.8} {:>14.8}\n",
            sym,
            p.x as f64,
            p.y as f64,
            p.z as f64
        ));
    }
    out
}
