// src/ui.rs
use bevy::prelude::*;
use bevy_egui::{egui, input::EguiWantsInput, EguiContexts};
use std::collections::BTreeSet;
use std::fs;
use std::path::PathBuf;
use std::sync::Arc;

use crate::color_schemes::{color_scheme_map, ELEMENT_SYMBOLS};
use crate::diagnostics::{self, DiagnosticSeverity};
use crate::events::MoleculeChanged;
use crate::molecule::{parse_xyz_first_frame_angstrom, Molecule};
use crate::settings::{BondColorMode, ColorScheme, LightingMode, MolSettings, RepresentationMode};

// Export UI + resources
use crate::export_image::{ExportCounter, ExportFormat, ExportSelectArea, ExportSettings};

// Measurements UI + resource
use crate::measurements::Measurements;
use crate::molecule_builder::builder_ui::{
    builder_ui_contents, builder_ui_panel, EditorRotateState, ZMatrixBuilderState,
};
use crate::qchem_interfaces::xtb_freq::xtb_frequency_panel;
use crate::rmsd::ui::rmsd_panel;
use crate::ui_layout::{icon_rail, section, Tab, UiLayout};
use crate::qchem_interfaces::xtb_optimize::xtb_optimization_panel;
use crate::trajectory;
use crate::ui_measurements;
use crate::uvvis::{ui::uvvis_panel, UvVisState};

/// Point size of the X/Y/Z letters on the corner axis gizmo. Big enough to
/// read at a glance without crowding the little arrows they label.
const AXIS_LABEL_FONT_SIZE: f32 = 19.0;

// NEW: shared picking helper (egui-friendly shim)
use crate::picking::screen::find_nearest_atom_screen_space_egui;

const ELEMENT_FILTER_S_BLOCK: [&str; 14] = [
    "H", "He", "Li", "Be", "Na", "Mg", "K", "Ca", "Rb", "Sr", "Cs", "Ba", "Fr", "Ra",
];
const ELEMENT_FILTER_D_BLOCK: [&str; 29] = [
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
];
const ELEMENT_FILTER_P_BLOCK: [&str; 29] = [
    "B", "C", "N", "O", "F", "Ne", "Al", "Si", "P", "S", "Cl", "Ar", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "In", "Sn", "Sb", "Te", "I", "Xe", "Tl", "Pb", "Bi", "Po", "At",
];
const ELEMENT_FILTER_F_BLOCK: [&str; 30] = [
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Ac",
    "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
];

/// Screen regions covered by the control panels, in egui points.
///
/// The panels are drawn into a hand-built root `egui::Ui` rather than through
/// the `Context`, so `Context::is_pointer_over_egui` cannot see them: for a
/// background layer it only tests the pointer against `root_ui_available_rect`,
/// which these panels never shrink.  Everything derived from `EguiWantsInput`
/// is therefore blind to a pointer hovering the panel, which let scroll wheel
/// input reach the orbit camera while scrolling a list.  Publishing the
/// geometry here gives the camera a test that does not depend on egui's
/// bookkeeping.
#[derive(Resource, Default, Clone)]
pub struct UiPanelRegions {
    /// The docked control panel (right in the classic layout, the icon rail
    /// in the floating-window one).
    pub controls: Option<Rect>,
    /// Whatever viewport space is left once every panel has been laid out.
    pub free_viewport: Option<Rect>,
    /// Floating tool windows. Unlike docked panels these sit *inside* the
    /// leftover viewport rect, so without listing them a click on a window
    /// would fall through and also pick the atom behind it.
    pub windows: Vec<Rect>,
}

impl UiPanelRegions {
    /// True when `point` (egui points) is over panel chrome rather than the 3D view.
    pub fn covers(&self, point: Vec2) -> bool {
        if self.controls.is_some_and(|r| r.contains(point)) {
            return true;
        }
        if self.windows.iter().any(|r| r.contains(point)) {
            return true;
        }
        // Also treat anything outside the leftover viewport as covered, which
        // catches panels that do not report a rect of their own.
        self.free_viewport.is_some_and(|r| !r.contains(point))
    }
}

#[derive(Resource, Clone)]
pub struct XyzBuffer {
    pub text: String,
    pub last_dir: Option<PathBuf>,
    pub current_file: Option<PathBuf>,
    pub warning: Option<String>,
    // Cached read-only XYZ rows; refreshed only when atoms or coordinates change.
    pub live_xyz_atoms: Vec<String>,
    pub live_xyz_pos: Vec<Vec3>,
    pub live_xyz_lines: Vec<String>,
    // Cached egui layout for atom labels; rebuilt only when atom symbols change.
    pub label_atoms: Vec<String>,
    pub type_label_galleys: Vec<Arc<egui::Galley>>,
    pub index_label_galleys: Vec<Arc<egui::Galley>>,
    pub type_index_label_galleys: Vec<Arc<egui::Galley>>,
    pub element_filter_atoms: Vec<String>,
    pub element_filter_extras: Vec<String>,
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
    mut commands: Commands,
    mut contexts: EguiContexts,
    egui_wants_input: Res<EguiWantsInput>,
    mut images: ResMut<Assets<Image>>,
    mut settings: ResMut<MolSettings>,
    mut mol: ResMut<Molecule>,
    mut xyz_buf: ResMut<XyzBuffer>,
    mut ev_changed: MessageWriter<MoleculeChanged>,
    // orbit target control
    mut cam: ResMut<crate::camera::OrbitCamera>,
    // camera for picking + label projection
    q_cam: Query<
        (&Camera, &Projection, &GlobalTransform),
        (With<Camera3d>, With<crate::scene::MainCamera>),
    >,
    windows: Query<&Window>,

    // export section
    export_resources: (
        ResMut<ExportSettings>,
        ResMut<ExportCounter>,
        ResMut<ExportSelectArea>,
    ),

    // measurements panel
    mut measurements: ResMut<Measurements>,
    mut diagnostics_cache: ResMut<diagnostics::DiagnosticsCache>,
    // trajectory
    mut traj: ResMut<trajectory::TrajectoryState>,
    // `ui_panel` sits at Bevy's 16-parameter system limit, so late additions
    // ride along in this tuple rather than becoming parameters of their own.
    builder_resources: (
        ResMut<EditorRotateState>,
        ResMut<ZMatrixBuilderState>,
        Local<bool>,
        ResMut<crate::orbitals::OrbitalState>,
        ResMut<UiPanelRegions>,
        ResMut<crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask>,
        ResMut<crate::qchem_interfaces::xtb_optimize::XtbPanelState>,
        ResMut<crate::qchem_interfaces::xtb_freq::XtbFrequencyTask>,
        ResMut<crate::qchem_interfaces::xtb_freq::XtbFreqPanelState>,
        ResMut<crate::rmsd::RmsdState>,
        ResMut<UiLayout>,
        ResMut<UvVisState>,
        ResMut<crate::uvvis::run::SpectrumTask>,
    ),
) {
    // bevy_egui 0.41: ctx_mut() returns Result; if it fails, skip this frame
    let ctx = contexts.ctx_mut().expect("missing primary egui context");
    let (mut export_settings, mut export_counter, mut export_select_area) = export_resources;
    let (
        mut editor_rotate_state,
        mut zmat_state,
        mut style_initialized,
        mut orbital_state,
        mut panel_regions,
        mut xtb_task,
        mut xtb_panel_state,
        mut xtb_freq_task,
        mut xtb_freq_panel_state,
        mut rmsd_state,
        mut ui_layout,
        mut uvvis_state,
        mut uvvis_task,
    ) = builder_resources;
    if !*style_initialized {
        ctx.style_mut_of(ctx.theme(), |style| {
            style.text_styles = [
                (egui::TextStyle::Heading, egui::FontId::proportional(22.0)),
                (egui::TextStyle::Body, egui::FontId::proportional(16.0)),
                (egui::TextStyle::Monospace, egui::FontId::monospace(15.0)),
                (egui::TextStyle::Button, egui::FontId::proportional(16.0)),
                (egui::TextStyle::Small, egui::FontId::proportional(13.0)),
            ]
            .into();

            // Tooltips carry the tool names on the icon rail, so they have to
            // appear as soon as the pointer is over an icon. By default egui
            // waits for the pointer to stop moving and then delays further,
            // which makes a rail of unlabelled icons feel unresponsive.
            style.interaction.show_tooltips_only_when_still = false;
            style.interaction.tooltip_delay = 0.0;
        });
        *style_initialized = true;
    }

    // -------------------------
    // Right-click picking (not over egui)
    // -------------------------
    let (latest_pos_opt, right_clicked) =
        ctx.input(|i| (i.pointer.latest_pos(), i.pointer.secondary_clicked()));
    let mouse_over_ui = egui_wants_input.is_using_pointer() || egui_wants_input.is_popup_open();
    if let (Some(pos_points), true) = (latest_pos_opt, right_clicked) {
        if !mouse_over_ui {
            if let (Ok((cam_comp, _proj, cam_xform)), Ok(win)) = (q_cam.single(), windows.single())
            {
                // Convert egui logical points → PHYSICAL px
                let scale = win.scale_factor() as f32;
                let mouse_px = egui::pos2(pos_points.x * scale, pos_points.y * scale);

                let (picked_idx, popup_pos) = if let Some((hit_idx, hit_px_pos)) =
                    find_nearest_atom_screen_space_egui(
                        cam_comp, cam_xform, &mol.pos, mouse_px, 18.0,
                    ) {
                    (
                        hit_idx as i32,
                        egui::pos2(hit_px_pos.x / scale, hit_px_pos.y / scale),
                    )
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
                    egui::Frame::popup(ui.style()).show(ui, |ui| {
                        if picked >= 0 {
                            ui.set_min_width(220.0);
                        }
                        ui.vertical(|ui| {
                            ui.label("Background:");
                            ui.horizontal(|ui| {
                                if ui.button("Black").clicked() {
                                    settings.bg_color = Color::srgb(0.0, 0.0, 0.0);
                                    settings.lighting_dirty = true;
                                }
                                if ui.button("White").clicked() {
                                    settings.bg_color = Color::srgb(1.0, 1.0, 1.0);
                                    settings.lighting_dirty = true;
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
    let mut viewport_ui = egui::Ui::new(
        ctx.clone(),
        "controls_viewport".into(),
        egui::UiBuilder::new()
            .layer_id(egui::LayerId::background())
            .max_rect(ctx.viewport_rect()),
    );
    // Both layouts draw exactly the same sections; only the container
    // differs, so the bodies live in one closure rather than being copied.
    let windowed = ui_layout.windowed;
    let mut open_windows = ui_layout.open;
    // Rects of everything egui drew as a floating window this frame, so the
    // 3D picker knows not to select atoms behind them.
    let mut window_rects: Vec<egui::Rect> = Vec::new();

    let controls_response = {
        let mut draw_sections = |ui: &mut egui::Ui,
                                 windowed: bool,
                                 open: &mut [bool; Tab::COUNT],
                                 rects: &mut Vec<egui::Rect>| {
            // ===========================
            // 1) Structure (XYZ)
            // ===========================
            ui.add_space(8.0);
            section(
                ui,
                windowed,
                &mut open[Tab::Structure.index()],
                rects,
                Tab::Structure.default_size(),
                "Structure (XYZ)",
                |ui| {
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
                        let dlg = crate::recent_dir::open().add_filter("XYZ", &["xyz"]);
                        if let Some(path) = crate::recent_dir::pick_file(dlg) {
                            xyz_buf.last_dir = path.parent().map(|p| p.to_path_buf());
                            xyz_buf.current_file = Some(path.clone());
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
                let xyz_file_label = match &xyz_buf.current_file {
                    Some(path) => format!("File: {}", path.display()),
                    None => "File: (none)".to_string(),
                };
                ui.collapsing("Show Filename", |ui| {
                    ui.weak(xyz_file_label);
                });

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
                    refresh_live_xyz_cache(&mut xyz_buf, &mol);
                    ui.label("Live XYZ (read-only):");
                    egui::ScrollArea::vertical()
                        .max_height(260.0)
                        .auto_shrink([false; 2])
                        .show(ui, |ui| {
                            ui.style_mut().override_text_style = Some(egui::TextStyle::Monospace);
                            for label in &xyz_buf.live_xyz_lines {
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
            // 1b) Geometry Optimization (xTB)
            // ===========================
            ui.add_space(8.0);
            section(
                ui,
                windowed,
                &mut open[Tab::Optimize.index()],
                rects,
                Tab::Optimize.default_size(),
                "Geometry Optimization (xTB)",
                |ui| {
                xtb_optimization_panel(
                    ui,
                    &mut xtb_panel_state,
                    &mut xtb_task,
                    &mol,
                    traj.playing,
                );
            });
            // ===========================
            // 1c) Frequency Analysis
            // ===========================
            ui.add_space(8.0);
            section(
                ui,
                windowed,
                &mut open[Tab::Vibrations.index()],
                rects,
                Tab::Vibrations.default_size(),
                "Frequency Analysis",
                |ui| {
                xtb_frequency_panel(
                    ui,
                    &mut xtb_freq_panel_state,
                    &mut xtb_freq_task,
                    &xtb_panel_state,
                    &mol,
                    &mut traj,
                );
                // A Hessian loaded from a file arrives with its own geometry.
                // Putting it on screen is what makes the mode list mean
                // anything -- the program starts with an empty viewport, so
                // otherwise there would be nothing for the modes to move.
                if let Some((atoms, pos)) = xtb_freq_task.pending_geometry.take() {
                    mol.atoms = atoms;
                    mol.pos = pos;
                    mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                    ev_changed.write(MoleculeChanged::parse_xyz(true));
                }
            });
            // ===========================
            // 1c-bis) UV-Vis / TD-DFT
            // ===========================
            ui.add_space(8.0);
            section(
                ui,
                windowed,
                &mut open[Tab::UvVis.index()],
                rects,
                Tab::UvVis.default_size(),
                Tab::UvVis.title(),
                |ui| {
                    // The Behemoth path is the optimizer panel's, so there
                    // is one place to set it rather than three.
                    let behemoth_path = xtb_panel_state
                        .paths
                        .iter()
                        .find(|(p, _)| {
                            *p == crate::qchem_interfaces::program::QcProgram::Behemoth
                        })
                        .map(|(_, path)| path.clone())
                        .unwrap_or_default();
                    uvvis_panel(
                        ui,
                        &mut uvvis_state,
                        &mut uvvis_task,
                        &mol,
                        &behemoth_path,
                        traj.playing,
                    );
                },
            );
            // ===========================
            // 1d) Trajectory
            // ===========================
            ui.add_space(8.0);
            section(
                ui,
                windowed,
                &mut open[Tab::Trajectory.index()],
                rects,
                Tab::Trajectory.default_size(),
                "Trajectory",
                |ui| {
                if ui.button("Load trajectory XYZ").clicked() {
                    let dlg = crate::recent_dir::open().add_filter("XYZ", &["xyz"]);
                    if let Some(path) = crate::recent_dir::pick_file(dlg) {
                        traj.last_dir = path.parent().map(|p| p.to_path_buf());
                        traj.current_file = Some(path.clone());
                        xyz_buf.warning = None;
                        if let Ok(text) = fs::read_to_string(&path) {
                            match trajectory::parse_multi_xyz(&text) {
                                Ok(frames) => {
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
                                    // `apply_current_frame` only reports a
                                    // coordinate update — it runs for every
                                    // playback frame too.  Loading a file
                                    // replaces the structure, so say so, or
                                    // stale measurements carry over.
                                    // `false`: the call above already centred
                                    // the camera.
                                    ev_changed.write(MoleculeChanged::parse_xyz(false));
                                }
                                Err(_) => {}
                            }
                        }
                    }
                }

                let traj_file_label = match &traj.current_file {
                    Some(path) => format!("File: {}", path.display()),
                    None => "File: (none)".to_string(),
                };
                ui.collapsing("Show Filename", |ui| {
                    ui.weak(traj_file_label);
                });

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
                        if ui
                            .add(
                                egui::Button::new("Fixed bonds")
                                    .fill(if traj.fixed_bonds { sel } else { dim }),
                            )
                            .on_hover_text(
                                "Keep the current bond topology during playback instead of \
                                 perceiving it per frame. Off by default -- a reactive \
                                 trajectory needs bonds that follow the reaction.",
                            )
                            .clicked()
                        {
                            traj.fixed_bonds = !traj.fixed_bonds;
                            traj.last_applied = None;
                            if traj.fixed_bonds {
                                traj.overlay_dirty = true;
                            }
                        }
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
            // 1e) Structure Comparison (RMSD / Kabsch)
            // ===========================
            ui.add_space(8.0);
            section(
                ui,
                windowed,
                &mut open[Tab::Compare.index()],
                rects,
                Tab::Compare.default_size(),
                "Structure Comparison (RMSD)",
                |ui| {
                rmsd_panel(ui, &mut rmsd_state, &mut traj, &mut mol);
            });
            // ===========================
            // 2) Measurements
            // ===========================
            ui.add_space(8.0);
            ui.separator();
            section(
                ui,
                windowed,
                &mut open[Tab::Measure.index()],
                rects,
                Tab::Measure.default_size(),
                "Measurements",
                |ui| {
                ui_measurements::measurements_panel(
                    ui,
                    &mut measurements,
                    &mol,
                    &settings,
                );
            });

            // ===========================
            // 2b) Surface tools
            // ===========================
            ui.add_space(8.0);
            ui.separator();
            section(
                ui,
                windowed,
                &mut open[Tab::Surface.index()],
                rects,
                Tab::Surface.default_size(),
                "Surface Tools",
                |ui| {
                crate::orbitals::ui::orbital_panel(
                    ui,
                    &mut orbital_state,
                    &mut mol,
                    &mut settings,
                    &mut ev_changed,
                );
            });

            // ===========================
            // 3) Diagnostics
            // ===========================
            ui.add_space(8.0);
            ui.separator();
            // The analysis behind this report is O(n) with neighbour queries but
            // still far too heavy to run per frame during playback, so tell the
            // refresh system whether anyone is actually looking at it.
            let diagnostics_visible = section(
                ui,
                windowed,
                &mut open[Tab::Diagnostics.index()],
                rects,
                Tab::Diagnostics.default_size(),
                "Diagnostics",
                |ui| {
                let report = &diagnostics_cache.report;
                let warnings = report.warning_count();
                let info = report.info_count();

                if report.items.is_empty() {
                    ui.weak("No valence, charge, isolation, or close-contact issues detected.");
                } else {
                    ui.horizontal(|ui| {
                        ui.colored_label(
                            egui::Color32::from_rgb(255, 170, 80),
                            format!("Warnings: {}", warnings),
                        );
                        ui.weak(format!("Info: {}", info));
                    });
                    ui.add_space(4.0);
                    egui::ScrollArea::vertical()
                        .max_height(220.0)
                        .auto_shrink([false; 2])
                        .show(ui, |ui| {
                            for item in &report.items {
                                let color = match item.severity {
                                    DiagnosticSeverity::Warning => {
                                        egui::Color32::from_rgb(255, 170, 80)
                                    }
                                    DiagnosticSeverity::Info => egui::Color32::LIGHT_BLUE,
                                };
                                let prefix = match item.atom {
                                    Some(idx) => format!("#{} ", idx + 1),
                                    None => String::new(),
                                };
                                ui.colored_label(color, format!("{}{}", prefix, item.message));
                            }
                        });
                    ui.add_space(4.0);
                    ui.weak("Diagnostics use perceived single-bond connectivity; multiple bonds and formal charges may need chemical review.");
                }
            });
            diagnostics_cache.wanted = diagnostics_visible;

            // ===========================
            // 4) Appearance
            // ===========================
            ui.add_space(8.0);
            ui.separator();
            section(
                ui,
                windowed,
                &mut open[Tab::Representation.index()],
                rects,
                Tab::Representation.default_size(),
                "Representation",
                |ui| {
                // `changed` = the entity set must be rebuilt (mode switch, atom
                // visibility).  `mesh_changed` = same entities, different mesh
                // handles, which is dramatically cheaper.
                let mut changed = false;
                let mut mesh_changed = false;

                if mol.atoms.len() > 1000 {
                    ui.weak("Large molecule: heavy mesh representations auto-switch to low-res.");
                    ui.add_space(4.0);
                }

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
                        mesh_changed = true;
                    }
                    if ui
                        .add(
                            egui::Slider::new(&mut settings.low_res_atom_scale, 0.2..=1.5)
                                .text("Low-res atom scale"),
                        )
                        .changed()
                    {
                        mesh_changed = true;
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
                        mesh_changed = true;
                    }
                    if ui
                        .add(
                            egui::Slider::new(&mut settings.cpk_atom_scale, 0.6..=2.6)
                                .text("CPK atom scale"),
                        )
                        .changed()
                    {
                        mesh_changed = true;
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
                        mesh_changed = true;
                    }
                }

                if matches!(
                    settings.representation,
                    RepresentationMode::LowResBallsAndLines
                        | RepresentationMode::LinesOnly
                        | RepresentationMode::BackboneTrace
                ) {
                    ui.add_space(6.0);
                    ui.add(
                        egui::Slider::new(&mut settings.line_bond_thickness, 1.0..=50.0)
                            .text("Line thickness"),
                    );
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
                            mesh_changed = true;
                        }
                        if ui
                            .add(
                                egui::Slider::new(&mut settings.trace_atom_scale, 0.2..=1.5)
                                    .text("Trace atom scale"),
                            )
                            .changed()
                        {
                            mesh_changed = true;
                        }
                    }
                }

                if changed {
                    settings.geometry_dirty = true;
                }
                if mesh_changed {
                    settings.meshes_dirty = true;
                }
                },
            );
            ui.add_space(8.0);
            ui.separator();
            section(
                ui,
                windowed,
                &mut open[Tab::Appearance.index()],
                rects,
                Tab::Appearance.default_size(),
                "Appearance",
                |ui| {

                ui.add_space(8.0);

                ui.collapsing("Element Filter", |ui| {
                    refresh_element_filter_cache(&mut xyz_buf, &mol);
                    let extras = &xyz_buf.element_filter_extras;
                    if mol.atoms.is_empty() {
                        ui.weak("No atoms loaded.");
                        return;
                    }

                    let mut changed = false;
                    ui.horizontal(|ui| {
                        if ui.button("Show all").clicked() {
                            for block in [
                                &ELEMENT_FILTER_S_BLOCK[..],
                                &ELEMENT_FILTER_D_BLOCK[..],
                                &ELEMENT_FILTER_P_BLOCK[..],
                                &ELEMENT_FILTER_F_BLOCK[..],
                            ] {
                                for sym in block {
                                    settings.element_visibility.insert((*sym).to_string(), true);
                                }
                            }
                            for sym in extras {
                                settings.element_visibility.insert(sym.clone(), true);
                            }
                            changed = true;
                        }
                        if ui.button("Hide all").clicked() {
                            for block in [
                                &ELEMENT_FILTER_S_BLOCK[..],
                                &ELEMENT_FILTER_D_BLOCK[..],
                                &ELEMENT_FILTER_P_BLOCK[..],
                                &ELEMENT_FILTER_F_BLOCK[..],
                            ] {
                                for sym in block {
                                    settings.element_visibility.insert((*sym).to_string(), false);
                                }
                            }
                            for sym in extras {
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
                        settings.geometry_dirty = true;
                    }
                });

                ui.add_space(8.0);

                ui.collapsing("Atom & Bond Scaling", |ui| {
                    // Sizes only re-point mesh handles; the cutoff sliders are the
                    // only ones here that change connectivity.
                    let mut mesh_changed = false;
                    let mut bond_topology_changed = false;
                    let mut hbond_topology_changed = false;

                    mesh_changed |= ui
                        .add(egui::Slider::new(&mut settings.atom_scale, 0.1..=3.0).text("Atom scale"))
                        .changed();
                    mesh_changed |= ui
                        .add(egui::Slider::new(&mut settings.atom_resolution, 0..=10).text("Atom resolution"))
                        .changed();

                    ui.add_space(8.0);
                    ui.separator();
                    ui.add_space(8.0);

                    mesh_changed |= ui
                        .add(
                            egui::Slider::new(&mut settings.bond_radius_pct, 0.05..=1.0)
                                .text("Bond radius (% of smaller atom)"),
                        )
                        .changed();

                    bond_topology_changed |= ui
                        .add(
                            egui::Slider::new(&mut settings.bond_thresh_scale, 0.8..=3.0)
                                .text("Bond cutoff × covalent radii"),
                        )
                        .changed();

                    ui.add_space(8.0);
                    ui.separator();
                    ui.add_space(8.0);

                    hbond_topology_changed |= ui
                        .add(
                            egui::Slider::new(&mut settings.hbond_cutoff, 1.2..=5.0)
                                .text("Hydrogen bond cutoff in Å"),
                        )
                        .changed();

                    // H-bond style
                    ui.add(
                        egui::Slider::new(&mut settings.hbond_thickness, 1.0..=80.0)
                            .text("H-bond thickness"),
                    );
                    ui.add(
                        egui::Slider::new(&mut settings.hbond_gap_scale, 0.2..=5.0)
                            .text("H-bond dash gap scale"),
                    );
                    ui.add(
                        egui::Slider::new(&mut settings.hbond_line_scale, 0.2..=5.0)
                            .text("H-bond dash line scale"),
                    );

                    if bond_topology_changed {
                        // `rebuild_if_dirty` recomputes connectivity itself; doing it
                        // here as well meant two full bond passes in the same frame.
                        settings.geometry_dirty = true;
                        settings.bond_topology_dirty = true;
                    } else if hbond_topology_changed {
                        // H-bonds are drawn straight from `mol.hydrogen_bonds` as
                        // gizmos, so no entity rebuild is needed — just the data.
                        mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                    }
                    if mesh_changed {
                        settings.meshes_dirty = true;
                    }
                });

                ui.add_space(8.0);

                ui.collapsing("Color & Materials", |ui| {
                    // Background quick picks + picker
                    ui.horizontal(|ui| {
                        ui.label("Background:");
                        if ui.button("Black").clicked() {
                            settings.bg_color = Color::srgb(0.0, 0.0, 0.0);
                            settings.lighting_dirty = true;
                        }
                        if ui.button("White").clicked() {
                            settings.bg_color = Color::srgb(1.0, 1.0, 1.0);
                            settings.lighting_dirty = true;
                        }
                        ui.label("or pick:");
                        let mut eg = color_to_egui(settings.bg_color);
                        if ui.color_edit_button_srgba(&mut eg).changed() {
                            settings.bg_color = egui_to_color(eg);
                            settings.lighting_dirty = true;
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
                        settings.materials_dirty = true;
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
                                            settings.custom_scheme_name.clear();
                                        }
                                        settings.materials_dirty = true;
                                    }
                                }

                                // The user's own schemes, alongside the
                                // built-ins: a saved scheme is a Custom one
                                // with a name, so selecting it loads its
                                // colours and adopts its name.
                                let saved = crate::custom_schemes::list_saved();
                                if !saved.is_empty() {
                                    ui.separator();
                                    for name in saved {
                                        let selected = settings.scheme == ColorScheme::Custom
                                            && settings.custom_scheme_name == name;
                                        if ui.selectable_label(selected, &name).clicked() {
                                            if let Some(colors) =
                                                crate::custom_schemes::load(&name)
                                            {
                                                settings.element_colors = colors;
                                                settings.scheme = ColorScheme::Custom;
                                                settings.custom_scheme_name = name.clone();
                                                settings.materials_dirty = true;
                                            }
                                        }
                                    }
                                }
                            });
                    });

                    // Naming and saving. Editing any element colour above
                    // already switches the scheme to Custom, so this is where
                    // that edit becomes something the user keeps.
                    ui.horizontal(|ui| {
                        ui.label("Name:");
                        ui.add(
                            egui::TextEdit::singleline(&mut settings.custom_scheme_name)
                                .desired_width(140.0)
                                .hint_text("my scheme"),
                        );
                    });
                    ui.horizontal(|ui| {
                        let name = settings.custom_scheme_name.trim().to_string();
                        let can_save = !name.is_empty();
                        let save_now = ui
                            .add_enabled(can_save, egui::Button::new("Save Color Scheme"))
                            .on_disabled_hover_text("Give the scheme a name first.")
                            .clicked();
                        let save_default = ui
                            .add_enabled(can_save, egui::Button::new("Save & Make Default"))
                            .on_disabled_hover_text("Give the scheme a name first.")
                            .on_hover_text("Also load this scheme when Beavyr starts.")
                            .clicked();

                        if save_now || save_default {
                            match crate::custom_schemes::save(&name, &settings.element_colors) {
                                Ok(saved_as) => {
                                    settings.custom_scheme_name = saved_as.clone();
                                    settings.scheme = ColorScheme::Custom;
                                    let mut message = format!("Saved \u{201c}{saved_as}\u{201d}.");
                                    if save_default {
                                        match crate::custom_schemes::set_default(&saved_as) {
                                            Ok(()) => message
                                                .push_str(" It will load when Beavyr starts."),
                                            Err(err) => {
                                                message = err;
                                            }
                                        }
                                    }
                                    settings.scheme_message = Some(message);
                                }
                                Err(err) => settings.scheme_message = Some(err),
                            }
                        }

                        // Deleting matters because the names are typed: a
                        // mistyped scheme would otherwise sit in the list for
                        // good.
                        let exists = crate::custom_schemes::list_saved().contains(&name);
                        if ui
                            .add_enabled(exists, egui::Button::new("Delete"))
                            .on_disabled_hover_text("No saved scheme by that name.")
                            .clicked()
                        {
                            settings.scheme_message = match crate::custom_schemes::delete(&name) {
                                Ok(()) => Some(format!("Deleted \u{201c}{name}\u{201d}.")),
                                Err(err) => Some(err),
                            };
                        }
                    });
                    if let Some(message) = &settings.scheme_message {
                        ui.weak(message);
                    }

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
                                                            settings.materials_dirty = true;
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
                            settings.geometry_dirty = true;
                        }

                        if ui
                            .add(egui::Button::new("Atom-split bonds").fill(if is_split { sel } else { dim }))
                            .clicked()
                        {
                            settings.bond_color_mode = BondColorMode::AtomSplit;
                            settings.geometry_dirty = true;
                        }
                    });

                    ui.horizontal(|ui| {
                        ui.label("Uniform bond color:");
                        let mut bc = color_to_egui(settings.uniform_bond_color);
                        if ui.color_edit_button_srgba(&mut bc).changed() {
                            settings.uniform_bond_color = egui_to_color(bc);
                            settings.materials_dirty = true;
                        }
                    });

                    ui.horizontal(|ui| {
                        ui.label("H-bond color:");
                        let mut hc = color_to_egui(settings.hbond_color);
                        if ui.color_edit_button_srgba(&mut hc).changed() {
                            settings.hbond_color = egui_to_color(hc);
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
                            settings.lighting_dirty = true;
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
                                    settings.lighting_dirty = true;
                                }
                                if ui.selectable_label(matches!(settings.lighting_mode, LightingMode::ThreePoint), "Three-point").clicked() {
                                    settings.lighting_mode = LightingMode::ThreePoint;
                                    settings.lighting_dirty = true;
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

                    if changed { settings.lighting_dirty = true; }
                });
            });

            // ===========================
            // 5) Export Image
            // ===========================
            ui.add_space(8.0);
            ui.separator();
            section(
                ui,
                windowed,
                &mut open[Tab::Export.index()],
                rects,
                Tab::Export.default_size(),
                "Export Image",
                |ui| {
                let warn_id = egui::Id::new("export_no_selection_warn");
                let mut warn = ctx.data_mut(|d| d.get_persisted::<bool>(warn_id).unwrap_or(false));
                if export_select_area.has_area {
                    warn = false;
                }

                ui.horizontal(|ui| {
                    ui.label("Format");
                    egui::ComboBox::from_id_salt("export_format")
                        .selected_text(match export_settings.format {
                            ExportFormat::Png => "PNG",
                            ExportFormat::Jpg => "JPG",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut export_settings.format, ExportFormat::Png, "PNG");
                            ui.selectable_value(&mut export_settings.format, ExportFormat::Jpg, "JPG");
                        });
                });
                ui.horizontal(|ui| {
                    ui.label("DPI");
                    egui::ComboBox::from_id_salt("export_dpi")
                        .selected_text(format!("{} DPI", export_settings.dpi))
                        .show_ui(ui, |ui| {
                            for dpi in [150_u32, 200, 300, 600, 1200] {
                                ui.selectable_value(
                                    &mut export_settings.dpi,
                                    dpi,
                                    format!("{dpi} DPI"),
                                );
                            }
                        });
                });
                ui.horizontal(|ui| {
                    ui.label("Canvas");
                    ui.add(egui::DragValue::new(&mut export_settings.canvas_width).range(256..=8192));
                    ui.label("x");
                    ui.add(egui::DragValue::new(&mut export_settings.canvas_height).range(256..=8192));
                });
                ui.horizontal(|ui| {
                    let label = if export_select_area.active {
                        "Selecting..."
                    } else {
                        "Select Area"
                    };
                    if ui.add(egui::Button::new(label)).clicked() {
                        warn = false;
                        export_select_area.active = true;
                        export_select_area.dragging = false;
                        export_select_area.start = None;
                        export_select_area.end = None;
                        export_select_area.has_area = false;
                    }
                    if export_select_area.has_area && ui.button("Clear Area").clicked() {
                        export_select_area.active = false;
                        export_select_area.dragging = false;
                        export_select_area.start = None;
                        export_select_area.end = None;
                        export_select_area.has_area = false;
                    }
                });
                ui.horizontal(|ui| {
                    let button_size = egui::vec2(120.0, 24.0);
                    if ui
                        .add_sized(button_size, egui::Button::new("Export Image"))
                        .clicked()
                    {
                        if !export_select_area.has_area {
                            warn = true;
                            return;
                        }
                        let (Ok((_cam, proj, cam_xform)), Ok(win)) = (q_cam.single(), windows.single()) else {
                            return;
                        };
                        let ok = crate::export_image::start_export_capture(
                            &mut commands,
                            &mut images,
                            &mut export_settings,
                            &mut export_counter,
                            &export_select_area,
                            &win,
                            proj,
                            cam_xform,
                        );
                        if ok {
                            export_select_area.active = false;
                            export_select_area.dragging = false;
                            export_select_area.start = None;
                            export_select_area.end = None;
                            export_select_area.has_area = false;
                        }
                    }
                });
                if warn {
                    ui.add_space(4.0);
                    ui.colored_label(
                        egui::Color32::from_rgb(220, 80, 80),
                        "Select an area before exporting.",
                    );
                }
                ctx.data_mut(|d| d.insert_persisted(warn_id, warn));
            });
        };

        if windowed {
            // A permanent, narrow rail on the left; every tool lives in its
            // own floating window that the rail toggles, so any number can be
            // open at once and dragged where the user wants them.
            let mut switch_layout = false;
            let rail = egui::Panel::left("tool_rail")
                .exact_size(50.0)
                .resizable(false)
                .show(&mut viewport_ui, |ui| {
                    icon_rail(ui, &mut open_windows);
                    // The layout switch lives at the foot of the rail, the
                    // usual home for it, and away from the orientation gizmo
                    // now in the bottom-right corner.
                    ui.with_layout(egui::Layout::bottom_up(egui::Align::Center), |ui| {
                        ui.add_space(8.0);
                        if ui
                            .add(
                                egui::Button::new(egui::RichText::new("▤").size(15.0))
                                    .min_size(egui::vec2(30.0, 26.0)),
                            )
                            .on_hover_text("Switch to the classic panel layout")
                            .clicked()
                        {
                            switch_layout = true;
                        }
                    });
                });
            // The windows themselves are drawn straight onto the context,
            // which is what lets them float. The `Ui` handed to
            // `draw_sections` is therefore only a carrier for the context --
            // and it is deliberately invisible, because the section bodies
            // also emit separators and spacing between sections that made
            // sense inside the old right panel. Drawn into the full-viewport
            // background those became lines straight across the screen.
            let mut scratch = egui::Ui::new(
                ctx.clone(),
                "windowed_sections".into(),
                egui::UiBuilder::new()
                    .layer_id(egui::LayerId::background())
                    // A finite, empty rect: `Rect::NOTHING` is built from
                    // infinities and makes egui's layout arithmetic produce
                    // NaN, which it asserts on.
                    .max_rect(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::Vec2::ZERO))
                    .invisible(),
            );
            draw_sections(&mut scratch, true, &mut open_windows, &mut window_rects);
            if switch_layout {
                ui_layout.windowed = false;
            }
            rail
        } else {
            egui::Panel::right("controls")
                .default_size(320.0)
                .min_size(240.0)
                .max_size(520.0)
                .resizable(true)
                .show(&mut viewport_ui, |ui| {
                    draw_sections(ui, false, &mut open_windows, &mut window_rects);
                })
        }
    };

    ui_layout.open = open_windows;

    // In the classic layout the switch floats bottom-left; the
    // floating-window layout puts it at the foot of its rail instead.
    if !windowed {
    egui::Area::new(egui::Id::new("layout_switch"))
        .anchor(egui::Align2::LEFT_BOTTOM, egui::vec2(12.0, -12.0))
        .show(&ctx, |ui| {
            egui::Frame::popup(ui.style())
                .corner_radius(egui::CornerRadius::same(14))
                .show(ui, |ui| {
                    if ui
                        .add(
                            egui::Button::new(egui::RichText::new("▦").size(16.0))
                                .min_size(egui::vec2(26.0, 26.0)),
                        )
                        .on_hover_text("Switch to the floating-window layout")
                        .clicked()
                    {
                        ui_layout.windowed = true;
                    }
                });
        });
    }

    // A small rounded floating button opens the molecule editor as a window
    // in the tabbed layout, instead of it permanently occupying the left
    // edge the way the classic layout docks it.
    if windowed {
        egui::Area::new(egui::Id::new("builder_launcher"))
            .anchor(egui::Align2::RIGHT_TOP, egui::vec2(-12.0, 12.0))
            .show(&ctx, |ui| {
                egui::Frame::popup(ui.style())
                    .corner_radius(egui::CornerRadius::same(18))
                    .show(ui, |ui| {
                        let label = egui::RichText::new("✏").size(18.0);
                        let button = egui::Button::new(label)
                            .min_size(egui::vec2(30.0, 30.0))
                            .selected(ui_layout.builder_open);
                        if ui
                            .add(button)
                            .on_hover_text("Molecule editor / Z-matrix")
                            .clicked()
                        {
                            ui_layout.builder_open = !ui_layout.builder_open;
                        }
                    });
            });

        let mut builder_open = ui_layout.builder_open;
        // Opens tall and on the right: the editor is a long column of
        // controls with the Z-matrix table under them, and the button that
        // opens it lives on the right edge, so appearing on the left meant
        // crossing the whole window to reach it. `default_*` only apply the
        // first time -- egui remembers wherever the user drags it afterwards.
        const EDITOR_WIDTH: f32 = 440.0;
        const EDITOR_MARGIN: f32 = 16.0;
        let viewport = ctx.viewport_rect();
        let editor_height = (viewport.height() - 2.0 * EDITOR_MARGIN).clamp(420.0, 980.0);
        let builder_window = egui::Window::new("Molecule Editor")
            .open(&mut builder_open)
            .resizable(true)
            .default_size([EDITOR_WIDTH, editor_height])
            .default_pos([
                (viewport.right() - EDITOR_WIDTH - EDITOR_MARGIN).max(viewport.left()),
                viewport.top() + EDITOR_MARGIN,
            ])
            .vscroll(true)
            .show(&ctx, |ui| {
                builder_ui_contents(
                    ui,
                    &mut editor_rotate_state,
                    &mut zmat_state,
                    &mut mol,
                    &mut settings,
                );
            });
        if let Some(window) = &builder_window {
            window_rects.push(window.response.rect);
        }
        ui_layout.builder_open = builder_open;
    }

    if !windowed {
        builder_ui_panel(
            &mut viewport_ui,
            &mut editor_rotate_state,
            &mut zmat_state,
            &mut mol,
            &mut settings,
        );
    }

    // Record what the panels ended up covering, now that they have all been
    // laid out, so the camera can tell panel from viewport next frame.
    let to_rect = |r: egui::Rect| Rect::new(r.min.x, r.min.y, r.max.x, r.max.y);
    panel_regions.controls = Some(to_rect(controls_response.response.rect));
    panel_regions.free_viewport = Some(to_rect(viewport_ui.available_rect_before_wrap()));
    // Includes the molecule-editor window, which is drawn separately below
    // the tool windows but covers the viewport just the same.
    panel_regions.windows = window_rects.iter().map(|r| to_rect(*r)).collect();

    // -------------------------
    // Overlay: X/Y/Z letters on the orientation gizmo
    // -------------------------
    // The gizmo shows which way the *original* axes now point after the
    // molecule has been rotated, which is only readable if the arrows are
    // named.
    if let (Ok((_cam_comp, _proj, cam_xform)), Ok(window)) = (q_cam.single(), windows.single()) {
        let physical = UVec2::new(window.physical_width(), window.physical_height());
        if physical.x > 0 && physical.y > 0 {
            let painter = ctx.layer_painter(egui::LayerId::new(
                egui::Order::Foreground,
                egui::Id::new("axis_gizmo_labels"),
            ));
            let labels = crate::scene::axis_gizmo_label_positions(
                physical,
                window.scale_factor() as f32,
                cam_xform,
            );
            for (name, pos, color) in labels {
                let rgba = color.to_srgba();
                painter.text(
                    egui::pos2(pos.x, pos.y),
                    egui::Align2::CENTER_CENTER,
                    name,
                    egui::FontId::proportional(AXIS_LABEL_FONT_SIZE),
                    egui::Color32::from_rgb(
                        (rgba.red * 255.0) as u8,
                        (rgba.green * 255.0) as u8,
                        (rgba.blue * 255.0) as u8,
                    ),
                );
            }
        }
    }

    // -------------------------
    // Overlay: atom labels (type/index)
    // -------------------------
    let show_type = ctx.data_mut(|d| {
        d.get_persisted::<bool>(egui::Id::new("show_atom_type"))
            .unwrap_or(false)
    });
    let show_index = ctx.data_mut(|d| {
        d.get_persisted::<bool>(egui::Id::new("show_atom_index"))
            .unwrap_or(false)
    });

    if (show_type || show_index) && mol.pos.len() == mol.atoms.len() {
        refresh_atom_label_cache(&mut xyz_buf, &mol, &ctx);
        if let (Ok((cam_comp, _proj, cam_xform)), Ok(win)) = (q_cam.single(), windows.single()) {
            egui::Area::new(egui::Id::new("atom_labels_overlay"))
                .movable(false)
                .interactable(false)
                .order(egui::Order::Tooltip)
                .show(&ctx, |ui| {
                    let painter = ui.painter();
                    let scale = win.scale_factor() as f32;

                    // Pre-laid-out galleys from `refresh_atom_label_cache`.  Laying
                    // text out per atom per frame here was the actual cost of having
                    // labels on for a large structure.
                    let galleys = match (show_type, show_index) {
                        (true, true) => &xyz_buf.type_index_label_galleys,
                        (true, false) => &xyz_buf.type_label_galleys,
                        (false, true) => &xyz_buf.index_label_galleys,
                        (false, false) => return,
                    };

                    for (i, &p_world) in mol.pos.iter().enumerate() {
                        let Some(galley) = galleys.get(i) else {
                            continue;
                        };
                        if let Ok(screen_px) = cam_comp.world_to_viewport(cam_xform, p_world) {
                            let pos = egui::pos2(screen_px.x / scale, screen_px.y / scale);
                            painter.galley(pos, galley.clone(), egui::Color32::WHITE);
                        }
                    }
                });
        }
    }

    if export_select_area.dragging || export_select_area.has_area {
        egui::Area::new(egui::Id::new("export_selection_overlay"))
            .movable(false)
            .interactable(false)
            .order(egui::Order::Foreground)
            .show(&ctx, |ui| {
                if let Some((min, max)) = export_select_area.bounds() {
                    let rect = egui::Rect::from_min_max(
                        egui::pos2(min.x, min.y),
                        egui::pos2(max.x, max.y),
                    );
                    let stroke = egui::Stroke::new(1.5, egui::Color32::from_rgb(70, 190, 90));
                    let fill = egui::Color32::from_rgba_premultiplied(70, 190, 90, 30);
                    ui.painter()
                        .rect(rect, 0.0, fill, stroke, egui::StrokeKind::Inside);
                }
            });
    }
}

// ---- helpers ----

fn is_standard_filter_element(symbol: &str) -> bool {
    ELEMENT_FILTER_S_BLOCK.contains(&symbol)
        || ELEMENT_FILTER_D_BLOCK.contains(&symbol)
        || ELEMENT_FILTER_P_BLOCK.contains(&symbol)
        || ELEMENT_FILTER_F_BLOCK.contains(&symbol)
}

fn refresh_element_filter_cache(cache: &mut XyzBuffer, mol: &Molecule) {
    if cache.element_filter_atoms == mol.atoms {
        return;
    }
    cache.element_filter_atoms.clone_from(&mol.atoms);
    let mut extras = BTreeSet::new();
    for symbol in &mol.atoms {
        if !is_standard_filter_element(symbol) {
            extras.insert(symbol.clone());
        }
    }
    cache.element_filter_extras.clear();
    cache.element_filter_extras.extend(extras);
}

fn refresh_atom_label_cache(cache: &mut XyzBuffer, mol: &Molecule, ctx: &egui::Context) {
    if cache.label_atoms == mol.atoms {
        return;
    }

    cache.label_atoms.clone_from(&mol.atoms);
    cache.type_label_galleys.clear();
    cache.index_label_galleys.clear();
    cache.type_index_label_galleys.clear();
    ctx.fonts_mut(|fonts| {
        for (index, symbol) in mol.atoms.iter().enumerate() {
            let font = egui::FontId::proportional(13.0);
            let color = egui::Color32::WHITE;
            cache.type_label_galleys.push(fonts.layout_no_wrap(
                symbol.clone(),
                font.clone(),
                color,
            ));
            cache.index_label_galleys.push(fonts.layout_no_wrap(
                (index + 1).to_string(),
                font.clone(),
                color,
            ));
            cache.type_index_label_galleys.push(fonts.layout_no_wrap(
                format!("{}({})", symbol, index + 1),
                font,
                color,
            ));
        }
    });
}

fn refresh_live_xyz_cache(cache: &mut XyzBuffer, mol: &Molecule) {
    if cache.live_xyz_atoms == mol.atoms && cache.live_xyz_pos == mol.pos {
        return;
    }

    cache.live_xyz_atoms.clone_from(&mol.atoms);
    cache.live_xyz_pos.clone_from(&mol.pos);
    cache.live_xyz_lines.clear();
    cache.live_xyz_lines.reserve(mol.atoms.len());
    for (index, (symbol, position)) in mol.atoms.iter().zip(&mol.pos).enumerate() {
        cache.live_xyz_lines.push(format!(
            "{:>3} {:<2} {:>7.3} {:>7.3} {:>7.3}",
            index + 1,
            symbol,
            position.x,
            position.y,
            position.z
        ));
    }
}

fn compute_centroid(points: &[Vec3]) -> Option<Vec3> {
    if points.is_empty() {
        return None;
    }
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
    ev_changed: &mut MessageWriter<MoleculeChanged>,
) -> Option<String> {
    let (atoms, qxyz, extra_frames) = parse_xyz_first_frame_angstrom(text);
    mol.atoms = atoms;
    mol.set_pos(
        qxyz.into_iter()
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
            sym, p.x as f64, p.y as f64, p.z as f64
        ));
    }
    out
}

#[cfg(test)]
mod panel_region_tests {
    use super::*;

    fn regions() -> UiPanelRegions {
        UiPanelRegions {
            // A left rail, as the floating-window layout uses.
            controls: Some(Rect::new(0.0, 0.0, 50.0, 800.0)),
            free_viewport: Some(Rect::new(50.0, 0.0, 1200.0, 800.0)),
            windows: vec![Rect::new(300.0, 200.0, 700.0, 600.0)],
        }
    }

    #[test]
    fn the_docked_rail_counts_as_covered() {
        assert!(regions().covers(Vec2::new(25.0, 400.0)));
    }

    #[test]
    fn empty_viewport_space_is_not_covered() {
        assert!(!regions().covers(Vec2::new(1000.0, 100.0)));
    }

    /// The regression this field exists for: a floating window sits *inside*
    /// the leftover viewport rect, so without listing it a click on the
    /// window would fall through and also pick the atom behind it.
    #[test]
    fn a_floating_window_counts_as_covered() {
        assert!(regions().covers(Vec2::new(500.0, 400.0)));
    }

    #[test]
    fn just_outside_a_window_is_not_covered() {
        assert!(!regions().covers(Vec2::new(705.0, 400.0)));
    }

    #[test]
    fn with_no_windows_open_the_viewport_is_clear() {
        let mut r = regions();
        r.windows.clear();
        assert!(!r.covers(Vec2::new(500.0, 400.0)));
    }
}
