// src/ui.rs
use bevy::prelude::*;
use bevy_egui::{egui, input::EguiWantsInput, EguiContexts};
use std::collections::BTreeSet;
use std::fs;
use std::path::PathBuf;
use std::sync::Arc;

use crate::color_schemes::{color_scheme_map, ELEMENT_SYMBOLS};
use crate::diagnostics;
use crate::events::MoleculeChanged;
use crate::molecule::{parse_xyz_first_frame_angstrom, Molecule};
use crate::settings::{BondColorMode, ColorScheme, LightingMode, MolSettings};

// Export UI + resources
use crate::export_image::{ExportCounter, ExportSelectArea, ExportSettings};

// Measurements UI + resource
use crate::measurements::Measurements;
use crate::molecule_builder::builder_ui::{
    builder_ui_contents, EditorRotateState, ZMatrixBuilderState,
};
use crate::qchem_interfaces::xtb_freq::xtb_frequency_panel;
use crate::rmsd::ui::rmsd_panel;
use crate::ui_layout::{icon_rail, section, Tab, UiLayout};
use crate::ui_style::ResponseExt;
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

#[derive(Resource, Clone, Default)]
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
        Res<crate::cli::StartupLoadReport>,
        ResMut<crate::conformer::ConformerRun>,
        (
            ResMut<crate::ts_generation::TsGeneration>,
            ResMut<crate::structure_history::StructureHistory>,
            ResMut<crate::qchem_panel::QcPanelState>,
            ResMut<crate::qchem_panel::run::QcRunTask>,
            ResMut<crate::project::ProjectState>,
        ),
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
        startup_report,
        mut conformer_run,
        (mut ts_generation, mut structure_history, mut qc_panel, mut qc_task, mut project),
    ) = builder_resources;
    crate::project::restore_ui(ctx, &mut project);
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
                    egui::Frame::popup(ui.style()).corner_radius(10).inner_margin(12).show(ui, |ui| {
                        if picked >= 0 {
                            ui.set_min_width(220.0);
                        }
                        ui.vertical(|ui| {
                            if let Some(action) = crate::ui_background_menu::contents(
                                ui, &mut settings, !mol.atoms.is_empty(),
                            ) {
                                use crate::ui_background_menu::Action;
                                let positions = match action {
                                    Action::Center => {
                                        if let Some(center) = compute_centroid(&mol.pos) {
                                            cam.target = center;
                                        }
                                        None
                                    }
                                    Action::Orient(plane) => Some(crate::orientation::orient(
                                        &mol.atoms, &mol.pos, plane,
                                    )),
                                    Action::Mirror(plane) => Some(crate::orientation::mirror(
                                        &mol.atoms, &mol.pos, plane,
                                    )),
                                    Action::Flip(axis) => Some(crate::orientation::flip(
                                        &mol.atoms, &mol.pos, axis,
                                    )),
                                    Action::Close => None,
                                };
                                if let Some(positions) = positions {
                                    mol.set_pos(positions);
                                    ev_changed.write(MoleculeChanged::set_pos());
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
                                if ui
                                    .add_enabled((picked as usize) < mol.atoms.len(), egui::Button::new("Delete atom"))
                                    .on_hover_text("Removes this atom and recreates the Z-matrix. All remaining atoms keep their current coordinates.")
                                    .clicked()
                                    && crate::molecule_builder::builder_ui::delete_atom_cartesian(
                                        &mut mol, &mut settings, &mut zmat_state,
                                        &mut editor_rotate_state, picked as usize,
                                    )
                                {
                                    measurements.reset_all();
                                    traj.clear_for_structure();
                                    xtb_freq_panel_state.selected_mode = None;
                                    xyz_buf.current_file = None;
                                    ev_changed.write(MoleculeChanged::parse_xyz(false));
                                    ctx.data_mut(|d| d.insert_persisted(pick_id, -1_i32));
                                    open = false;
                                }
                            }
                            if ui.button("Close").clicked() {
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
            crate::ui_layout::editor_section(
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

                    ui.horizontal_wrapped(|ui| {
                        let sel = ui.visuals().selection.bg_fill;
                        let dim = ui.visuals().widgets.inactive.bg_fill;

                        crate::project::menu(ui, &mut project);
                        if ui.add(crate::ui_style::primary("Open XYZ…")).clicked_once() {
                            let dlg = crate::recent_dir::open().add_filter("XYZ", &["xyz"]);
                            if let Some(path) = crate::recent_dir::pick_file(dlg) {
                                let previous = structure_history.before(&mol, &traj);
                                match load_xyz_from_path(
                                    &path,
                                    &mut xyz_buf,
                                    &mut traj,
                                    &mut mol,
                                    &mut settings,
                                    &mut cam,
                                    &mut ev_changed,
                                ) {
                                    Ok(_) => {
                                        xtb_freq_panel_state.selected_mode = None;
                                        let name =
                                            path.file_name().unwrap_or_default().to_string_lossy();
                                        structure_history.replaced(previous, &mol, &traj, &name);
                                    }
                                    Err(error) => xyz_buf.warning = Some(error),
                                }
                            }
                        }

                        if ui
                            .add(
                                egui::Button::new(if edit_mode {
                                    "View mode"
                                } else {
                                    "Edit as text"
                                })
                                .fill(if edit_mode {
                                    sel
                                } else {
                                    dim
                                }),
                            )
                            .on_hover_text("Toggle free-form XYZ editing")
                            .clicked_once()
                        {
                            edit_mode = !edit_mode;
                            if edit_mode {
                                // entering edit mode: prefill from live geometry
                                xyz_buf.text = format_xyz_from_molecule(&mol);
                            }
                            ctx.data_mut(|d| d.insert_persisted(edit_mode_id, edit_mode));
                        }

                        if ui
                            .button("Copy XYZ")
                            .on_hover_text(
                                "Copy the current molecule, including while editing a draft.",
                            )
                            .clicked_once()
                        {
                            let live = format_xyz_from_molecule(&mol);
                            ctx.copy_text(live);
                        }

                        // The structure as it is on screen, to a file. Named after
                        // the file it came from when there is one, so saving an
                        // edited copy starts from a sensible name.
                        if ui
                            .add_enabled(!mol.atoms.is_empty(), egui::Button::new("Save XYZ"))
                            .on_disabled_hover_text("There is no structure to save.")
                            .clicked_once()
                        {
                            let suggested = xyz_buf
                                .current_file
                                .as_ref()
                                .and_then(|path| path.file_name())
                                .map(|name| name.to_string_lossy().to_string())
                                .unwrap_or_else(|| "structure.xyz".to_string());
                            let dlg =
                                crate::recent_dir::save(&suggested).add_filter("XYZ", &["xyz"]);
                            if let Some(path) = crate::recent_dir::save_file(dlg) {
                                match std::fs::write(&path, format_xyz_from_molecule(&mol)) {
                                    Ok(()) => xyz_buf.warning = None,
                                    Err(error) => {
                                        xyz_buf.warning = Some(format!(
                                            "Could not save {}: {error}",
                                            path.display()
                                        ))
                                    }
                                }
                            }
                        }
                    });

                    ui.horizontal(|ui| {
                        let summary_width = (ui.available_width() - 72.0).max(60.0);
                        ui.allocate_ui_with_layout(
                            egui::vec2(summary_width, 30.0),
                            egui::Layout::left_to_right(egui::Align::Center),
                            |ui| {
                                crate::ui_structure::summary(
                                    ui,
                                    &crate::molecule::composition(&mol.atoms),
                                    xyz_buf.current_file.as_deref(),
                                );
                            },
                        );
                        ui.menu_button("View", |ui| {
                            if ui.checkbox(&mut show_type, "Element symbols").changed() {
                                ctx.data_mut(|d| d.insert_persisted(show_type_id, show_type));
                            }
                            if ui.checkbox(&mut show_index, "Atom numbers").changed() {
                                ctx.data_mut(|d| d.insert_persisted(show_index_id, show_index));
                            }
                            ui.checkbox(&mut settings.show_hbonds, "Hydrogen bonds");
                            if ui.button("Center molecule").clicked_once() {
                                if let Some(c) = compute_centroid(&mol.pos) {
                                    cam.target = c;
                                }
                                ui.close();
                            }
                            if startup_report.messages.iter().any(|(_, error)| !error) {
                                ui.separator();
                                for (message, is_error) in &startup_report.messages {
                                    if !is_error {
                                        ui.small(message);
                                    }
                                }
                            }
                        });
                    });

                    // Failed command-line loads and invalid drafts must remain
                    // visible; successful file-loading details live in View.
                    for (message, is_error) in &startup_report.messages {
                        if *is_error {
                            ui.colored_label(egui::Color32::from_rgb(220, 110, 60), message);
                        }
                    }
                    if let Some(warn) = &xyz_buf.warning {
                        ui.colored_label(
                            egui::Color32::from_rgb(255, 170, 0),
                            format!("Warning: {}", warn),
                        );
                    }
                    if !edit_mode {
                        refresh_live_xyz_cache(&mut xyz_buf, &mol);
                    }
                    let draft = &mut *xyz_buf;
                    match crate::ui_structure::coordinates(
                        ui,
                        windowed,
                        edit_mode,
                        &mut draft.text,
                        &draft.live_xyz_lines,
                    ) {
                        crate::ui_structure::CoordinateAction::Apply => {
                            traj.clear_for_structure();
                            xtb_freq_panel_state.selected_mode = None;
                            xyz_buf.warning =
                                apply_xyz_text(&xyz_buf.text, &mut mol, &mut cam, &mut ev_changed);
                        }
                        crate::ui_structure::CoordinateAction::Clear => xyz_buf.text.clear(),
                        crate::ui_structure::CoordinateAction::Reset => {
                            xyz_buf.text = format_xyz_from_molecule(&mol);
                        }
                        crate::ui_structure::CoordinateAction::None => {}
                    }
                },
            );

            ui.add_space(8.0);
            ui.separator();

            // ===========================
            // Recent structures, shared between sessions and viewer windows.
            // ===========================
            section(
                ui,
                windowed,
                &mut open[Tab::RecentStructures.index()],
                rects,
                Tab::RecentStructures.default_size(),
                Tab::RecentStructures.title(),
                |ui| {
                    if crate::structure_history::panel(
                        ui,
                        &mut structure_history,
                        &mut mol,
                        &mut traj,
                    ) {
                        xtb_freq_panel_state.selected_mode = None;
                        xyz_buf.current_file = None;
                        xyz_buf.text = format_xyz_from_molecule(&mol);
                        mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                        ev_changed.write(MoleculeChanged::parse_xyz(true));
                    }
                },
            );
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
                },
            );
            // ===========================
            // TS guess generation (RDA / Poor Man's NEB)
            // ===========================
            let mut show_ts_trajectory = false;
            section(
                ui,
                windowed,
                &mut open[Tab::TsGeneration.index()],
                rects,
                Tab::TsGeneration.default_size(),
                Tab::TsGeneration.title(),
                |ui| {
                    let previous = structure_history.before(&mol, &traj);
                    if crate::ts_generation::ui::panel(ui, &mut ts_generation, &mut mol, &mut traj)
                    {
                        structure_history.replaced(previous, &mol, &traj, "TS tool structure");
                        // An explicit TS/endpoint display replaces any vibration animation.
                        xtb_freq_panel_state.selected_mode = None;
                        mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                        ev_changed.write(MoleculeChanged::parse_xyz(true));
                        show_ts_trajectory = traj.frames.len() > 1;
                    }
                },
            );
            if show_ts_trajectory {
                open[Tab::Trajectory.index()] = true;
            }
            // ===========================
            // Conformer Search
            // ===========================
            ui.add_space(8.0);
            section(
                ui,
                windowed,
                &mut open[Tab::Conformers.index()],
                rects,
                Tab::Conformers.default_size(),
                "Conformer Search",
                |ui| {
                    let previous = structure_history.before(&mol, &traj);
                    let loaded = crate::conformer::ui::conformer_panel(
                        ui,
                        &mut conformer_run,
                        &mut mol,
                        &mut traj,
                        &settings,
                    );
                    // Showing a conformer, or loading the set as a trajectory, replaces the
                    // geometry -- so the renderer, measurements and bond lists have to be told.
                    if loaded {
                        xtb_freq_panel_state.selected_mode = None;
                        structure_history.replaced(previous, &mol, &traj, "Conformer");
                    }
                    crate::conformer::ui::announce(loaded, &mut ev_changed);
                },
            );
            // ===========================
            // 1c) Frequency Analysis
            // ===========================
            ui.add_space(8.0);
            // Closing the window stops any mode animation it started and puts
            // the equilibrium geometry back, so the user is not left hunting a
            // long mode list for the row whose Stop button is live.
            crate::qchem_interfaces::xtb_freq::stop_animation_when_panel_closes(
                &mut xtb_freq_panel_state,
                &xtb_freq_task,
                &mut traj,
                open[Tab::Vibrations.index()],
            );
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
                        let previous = structure_history.before(&mol, &traj);
                        traj.clear_for_structure();
                        mol.atoms = atoms;
                        mol.pos = pos;
                        structure_history.replaced(
                            previous,
                            &mol,
                            &traj,
                            "Imported frequency structure",
                        );
                        mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                        ev_changed.write(MoleculeChanged::parse_xyz(true));
                    }
                },
            );
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
                        .find(|(p, _)| *p == crate::qchem_interfaces::program::QcProgram::Behemoth)
                        .map(|(_, path)| path.clone())
                        .unwrap_or_default();
                    uvvis_panel(
                        ui,
                        &mut uvvis_state,
                        &mut uvvis_task,
                        &mol,
                        &behemoth_path,
                        traj.playing,
                        &mut orbital_state,
                    );
                },
            );
            if let Some(path) = uvvis_state.orbital_link.pending_molden.take() {
                let previous = structure_history.before(&mol, &traj);
                if crate::orbitals::ui::load_molden_from_path(&path, &mut orbital_state, &mut mol, &mut settings, &mut ev_changed) {
                    traj.clear_for_structure();
                    structure_history.replaced(previous, &mol, &traj, "Imported spectrum orbitals");
                    xtb_freq_panel_state.selected_mode = None;
                    uvvis_state.orbital_link.bind(&orbital_state);
                    uvvis_state.load_error = None;
                } else {
                    uvvis_state.load_error = orbital_state.warning.clone();
                }
            }
            if std::mem::take(&mut uvvis_state.orbital_link.reveal) {
                open[Tab::Surface.index()] = true;
            }
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
                    crate::ui_style::heading(
                        ui,
                        "Explore a trajectory",
                        "Browse frames, play motion and compare overlays.",
                    );
                    crate::ui_style::group(ui, "Trajectory file", |ui| {
                        if ui
                            .add(crate::ui_style::primary("Open trajectory XYZ…"))
                            .clicked_once()
                        {
                            let dlg = crate::recent_dir::open().add_filter("XYZ", &["xyz"]);
                            if let Some(path) = crate::recent_dir::pick_file(dlg) {
                                traj.last_dir = path.parent().map(|p| p.to_path_buf());
                                traj.current_file = Some(path.clone());
                                xyz_buf.warning = None;
                                if let Ok(text) = fs::read_to_string(&path) {
                                    match trajectory::parse_multi_xyz(&text) {
                                        Ok(frames) => {
                                            let previous = structure_history.before(&mol, &traj);
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
                                            let name = path
                                                .file_name()
                                                .unwrap_or_default()
                                                .to_string_lossy();
                                            structure_history
                                                .replaced(previous, &mol, &traj, &name);
                                            xtb_freq_panel_state.selected_mode = None;
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
                        ui.collapsing("File location", |ui| {
                            ui.weak(traj_file_label);
                        });
                    });

                    if traj.frames.is_empty() {
                        ui.weak("No trajectory loaded.");
                    } else {
                        let total = traj.frames.len();
                        crate::ui_style::group(ui, "Playback", |ui| {
                            ui.label(format!("Frame {}/{}", traj.current_frame + 1, total));

                            let mut frame_display = (traj.current_frame + 1) as i32;
                            if ui
                                .add(
                                    egui::Slider::new(&mut frame_display, 1..=total as i32)
                                        .text("Frame"),
                                )
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

                            ui.horizontal_wrapped(|ui| {
                                if ui.button("First").clicked_once() {
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
                                if ui.button("Play backward").clicked_once() {
                                    traj.direction = -1;
                                    traj.playing = true;
                                }
                                if ui.button("Pause").clicked_once() {
                                    traj.playing = false;
                                }
                                if ui.add(crate::ui_style::primary("Play")).clicked_once() {
                                    traj.direction = 1;
                                    traj.playing = true;
                                }
                                if ui.button("Last").clicked_once() {
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

                            ui.horizontal_wrapped(|ui| {
                                if ui.button("Previous").clicked_once() {
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
                                if ui.button("Next").clicked_once() {
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

                            ui.horizontal_wrapped(|ui| {
                                ui.label("Speed (fps)");
                                ui.add(
                                    egui::Slider::new(&mut traj.fps, 1.0..=100.0).show_value(true),
                                );
                            });

                            ui.horizontal_wrapped(|ui| {
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
                            .clicked_once()
                        {
                            traj.fixed_bonds = !traj.fixed_bonds;
                            traj.last_applied = None;
                            if traj.fixed_bonds {
                                traj.overlay_dirty = true;
                            }
                        }
                            });

                            ui.horizontal_wrapped(|ui| {
                                let sel = ui.visuals().selection.bg_fill;
                                let dim = ui.visuals().widgets.inactive.bg_fill;

                                let is_loop = matches!(traj.mode, trajectory::PlaybackMode::Loop);
                                let is_once = matches!(traj.mode, trajectory::PlaybackMode::Once);
                                let is_ping =
                                    matches!(traj.mode, trajectory::PlaybackMode::PingPong);

                                if ui
                                    .add(egui::Button::new("Loop").fill(if is_loop {
                                        sel
                                    } else {
                                        dim
                                    }))
                                    .clicked_once()
                                {
                                    traj.mode = trajectory::PlaybackMode::Loop;
                                }
                                if ui
                                    .add(egui::Button::new("Once").fill(if is_once {
                                        sel
                                    } else {
                                        dim
                                    }))
                                    .clicked_once()
                                {
                                    traj.mode = trajectory::PlaybackMode::Once;
                                }
                                if ui
                                    .add(egui::Button::new("Ping-pong").fill(if is_ping {
                                        sel
                                    } else {
                                        dim
                                    }))
                                    .clicked_once()
                                {
                                    traj.mode = trajectory::PlaybackMode::PingPong;
                                }
                            });
                        });
                        crate::ui_style::group(ui, "Frame overlays", |ui| {
                            ui.horizontal_wrapped(|ui| {
                                let sel = ui.visuals().selection.bg_fill;
                                let dim = ui.visuals().widgets.inactive.bg_fill;
                                ui.add_space(4.0);
                                if ui
                                    .add(
                                        egui::Button::new("Show overlays")
                                            .fill(if traj.overlay_enabled { sel } else { dim }),
                                    )
                                    .clicked_once()
                                {
                                    traj.overlay_enabled = !traj.overlay_enabled;
                                    traj.overlay_dirty = true;
                                }
                            });

                            ui.horizontal_wrapped(|ui| {
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

                            ui.horizontal_wrapped(|ui| {
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

                            ui.horizontal_wrapped(|ui| {
                                let sel = ui.visuals().selection.bg_fill;
                                let dim = ui.visuals().widgets.inactive.bg_fill;
                                if ui
                                    .add(
                                        egui::Button::new("Keep last frame solid")
                                            .fill(if traj.overlay_full_last { sel } else { dim }),
                                    )
                                    .clicked_once()
                                {
                                    traj.overlay_full_last = !traj.overlay_full_last;
                                    traj.overlay_dirty = true;
                                }
                            });

                            ui.horizontal_wrapped(|ui| {
                                ui.label("Ghost transparency");
                                if ui
                                    .add(egui::Slider::new(&mut traj.ghost_alpha, 0.05..=0.8))
                                    .changed()
                                {
                                    traj.overlay_dirty = true;
                                }
                            });
                        });
                    }
                },
            );
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
                },
            );
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
                    ui_measurements::measurements_panel(ui, &mut measurements, &mol, &settings);
                },
            );

            // ===========================
            // 2b) Surface tools
            // ===========================
            ui.add_space(8.0);
            ui.separator();
            crate::ui_layout::editor_section(
                ui,
                windowed,
                &mut open[Tab::Surface.index()],
                rects,
                Tab::Surface.default_size(),
                "Surface Tools",
                |ui| {
                    let previous = structure_history.before(&mol, &traj);
                    if crate::orbitals::ui::orbital_panel(
                        ui,
                        windowed,
                        &mut orbital_state,
                        &mut mol,
                        &mut settings,
                        &mut ev_changed,
                    ) {
                        traj.clear_for_structure();
                        xtb_freq_panel_state.selected_mode = None;
                        let name = orbital_state
                            .file_name
                            .clone()
                            .unwrap_or_else(|| "Molden structure".into());
                        structure_history.replaced(previous, &mol, &traj, &name);
                    }
                },
            );

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
                    crate::ui_view_panels::diagnostics_panel(
                        ui,
                        &diagnostics_cache.report,
                        !mol.atoms.is_empty(),
                    );
                },
            );
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
                    crate::ui_view_panels::representation_panel(ui, &mut settings, mol.atoms.len());
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

                    crate::ui_style::heading(
                        ui,
                        "Style the molecular view",
                        "Colors, lighting, sizes and element visibility.",
                    );
                    crate::ui_style::collapsible_group(ui, "Colors and materials", true, |ui| {
                        // Background quick picks + picker
                        ui.horizontal_wrapped(|ui| {
                            ui.label("Background:");
                            if ui.button("Black").clicked_once() {
                                settings.bg_color = Color::srgb(0.0, 0.0, 0.0);
                                settings.lighting_dirty = true;
                            }
                            if ui.button("White").clicked_once() {
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
                            .add(
                                egui::Slider::new(&mut settings.metallic, 0.0..=0.5)
                                    .text("Metallic"),
                            )
                            .changed();
                        mat_changed |= ui
                            .add(
                                egui::Slider::new(&mut settings.roughness, 0.005..=1.0)
                                    .text("Roughness"),
                            )
                            .changed();
                        mat_changed |= ui
                            .add(
                                egui::Slider::new(&mut settings.reflectance, 0.0..=1.0)
                                    .text("Reflectance"),
                            )
                            .changed();
                        if mat_changed {
                            settings.materials_dirty = true;
                        }

                        ui.add_space(6.0);
                        ui.label("Colors:");
                        ui.horizontal_wrapped(|ui| {
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
                                            .selectable_label(
                                                settings.scheme == sc,
                                                format!("{:?}", sc),
                                            )
                                            .clicked_once()
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
                                            if ui.selectable_label(selected, &name).clicked_once() {
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
                        ui.horizontal_wrapped(|ui| {
                            ui.label("Name:");
                            ui.add(
                                egui::TextEdit::singleline(&mut settings.custom_scheme_name)
                                    .desired_width(140.0)
                                    .hint_text("my scheme"),
                            );
                        });
                        ui.horizontal_wrapped(|ui| {
                            let name = settings.custom_scheme_name.trim().to_string();
                            let can_save = !name.is_empty();
                            let save_now = ui
                                .add_enabled(can_save, egui::Button::new("Save scheme"))
                                .on_disabled_hover_text("Give the scheme a name first.")
                                .clicked_once();
                            let save_default = ui
                                .add_enabled(can_save, egui::Button::new("Save as default"))
                                .on_disabled_hover_text("Give the scheme a name first.")
                                .on_hover_text("Also load this scheme when Beavyr starts.")
                                .clicked_once();

                            if save_now || save_default {
                                match crate::custom_schemes::save(&name, &settings.element_colors) {
                                    Ok(saved_as) => {
                                        settings.custom_scheme_name = saved_as.clone();
                                        settings.scheme = ColorScheme::Custom;
                                        let mut message =
                                            format!("Saved \u{201c}{saved_as}\u{201d}.");
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
                                .clicked_once()
                            {
                                settings.scheme_message = match crate::custom_schemes::delete(&name)
                                {
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
                                    let per_row =
                                        ((ui.available_width() / 78.0) as usize).clamp(1, 6);
                                    for chunk in ELEMENT_SYMBOLS.chunks(per_row) {
                                        ui.horizontal_wrapped(|ui| {
                                            for &sym in chunk {
                                                if let Some(old) =
                                                    settings.element_colors.get(sym).copied()
                                                {
                                                    let key = sym.to_string();
                                                    let mut col = color_to_egui(old);
                                                    ui.allocate_ui_with_layout(
                                                        egui::vec2(70.0, 0.0),
                                                        egui::Layout::left_to_right(
                                                            egui::Align::Center,
                                                        ),
                                                        |ui| {
                                                            ui.add_sized(
                                                                [22.0, 0.0],
                                                                egui::Label::new(sym),
                                                            );
                                                            if ui
                                                                .color_edit_button_srgba(&mut col)
                                                                .changed()
                                                            {
                                                                settings.element_colors.insert(
                                                                    key.clone(),
                                                                    egui_to_color(col),
                                                                );
                                                                settings.scheme =
                                                                    ColorScheme::Custom;
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
                        ui.horizontal_wrapped(|ui| {
                            let sel = ui.visuals().selection.bg_fill;
                            let dim = ui.visuals().widgets.inactive.bg_fill;

                            let is_uniform =
                                matches!(settings.bond_color_mode, BondColorMode::Uniform);
                            let is_split =
                                matches!(settings.bond_color_mode, BondColorMode::AtomSplit);

                            if ui
                                .add(egui::Button::new("Uniform bonds").fill(if is_uniform {
                                    sel
                                } else {
                                    dim
                                }))
                                .clicked_once()
                            {
                                settings.bond_color_mode = BondColorMode::Uniform;
                                settings.geometry_dirty = true;
                            }

                            if ui
                                .add(egui::Button::new("Atom-split bonds").fill(if is_split {
                                    sel
                                } else {
                                    dim
                                }))
                                .clicked_once()
                            {
                                settings.bond_color_mode = BondColorMode::AtomSplit;
                                settings.geometry_dirty = true;
                            }
                        });

                        ui.horizontal_wrapped(|ui| {
                            ui.label("Uniform bond color:");
                            let mut bc = color_to_egui(settings.uniform_bond_color);
                            if ui.color_edit_button_srgba(&mut bc).changed() {
                                settings.uniform_bond_color = egui_to_color(bc);
                                settings.materials_dirty = true;
                            }
                        });

                        ui.horizontal_wrapped(|ui| {
                            ui.label("H-bond color:");
                            let mut hc = color_to_egui(settings.hbond_color);
                            if ui.color_edit_button_srgba(&mut hc).changed() {
                                settings.hbond_color = egui_to_color(hc);
                            }
                        });
                    });

                    ui.add_space(8.0);

                    crate::ui_style::collapsible_group(ui, "Atom and bond size", false, |ui| {
                        // Sizes only re-point mesh handles; the cutoff sliders are the
                        // only ones here that change connectivity.
                        let mut mesh_changed = false;
                        let mut bond_topology_changed = false;
                        let mut hbond_topology_changed = false;

                        mesh_changed |= ui
                            .add(
                                egui::Slider::new(&mut settings.atom_scale, 0.1..=3.0)
                                    .text("Atom scale"),
                            )
                            .changed();
                        mesh_changed |= ui
                            .add(
                                egui::Slider::new(&mut settings.atom_resolution, 0..=10)
                                    .text("Atom resolution"),
                            )
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

                    crate::ui_style::collapsible_group(ui, "Lighting", false, |ui| {
                        let mut changed = false;

                        // Ambient
                        changed |= ui
                            .checkbox(&mut settings.use_ambient, "Ambient (omnidirectional fill)")
                            .changed();
                        ui.horizontal_wrapped(|ui| {
                            ui.label("Ambient color");
                            let mut ac = color_to_egui(settings.ambient_color);
                            if ui.color_edit_button_srgba(&mut ac).changed() {
                                settings.ambient_color = egui_to_color(ac);
                                settings.lighting_dirty = true;
                            }
                        });
                        changed |= ui
                            .add(
                                egui::Slider::new(&mut settings.ambient_brightness, 0.0..=3000.0)
                                    .text("Ambient brightness"),
                            )
                            .changed();

                        ui.separator();

                        // Lighting rig
                        ui.horizontal_wrapped(|ui| {
                            ui.label("Rig:");
                            egui::ComboBox::from_id_salt("lighting_mode_combo")
                                .selected_text(match settings.lighting_mode {
                                    LightingMode::SinglePoint => "Single point",
                                    LightingMode::ThreePoint => "Three-point",
                                })
                                .show_ui(ui, |ui| {
                                    if ui
                                        .selectable_label(
                                            matches!(
                                                settings.lighting_mode,
                                                LightingMode::SinglePoint
                                            ),
                                            "Single point",
                                        )
                                        .clicked_once()
                                    {
                                        settings.lighting_mode = LightingMode::SinglePoint;
                                        settings.lighting_dirty = true;
                                    }
                                    if ui
                                        .selectable_label(
                                            matches!(
                                                settings.lighting_mode,
                                                LightingMode::ThreePoint
                                            ),
                                            "Three-point",
                                        )
                                        .clicked_once()
                                    {
                                        settings.lighting_mode = LightingMode::ThreePoint;
                                        settings.lighting_dirty = true;
                                    }
                                });
                        });

                        // Key
                        changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut settings.light_intensity,
                                    50_000.0..=8_000_000.0,
                                )
                                .text("Key intensity"),
                            )
                            .changed();
                        changed |= ui
                            .add(
                                egui::Slider::new(&mut settings.light_distance, 3.0..=40.0)
                                    .text("Key distance"),
                            )
                            .changed();

                        // Fill / Rim (only in three-point)
                        if matches!(settings.lighting_mode, LightingMode::ThreePoint) {
                            ui.separator();
                            ui.label("Fill & Rim (three-point)");

                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut settings.fill_intensity,
                                        0.0..=10_000_000.0,
                                    )
                                    .text("Fill intensity"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(&mut settings.fill_distance, 3.0..=50.0)
                                        .text("Fill distance"),
                                )
                                .changed();

                            changed |= ui
                                .add(
                                    egui::Slider::new(
                                        &mut settings.rim_intensity,
                                        0.0..=10_000_000.0,
                                    )
                                    .text("Rim intensity"),
                                )
                                .changed();
                            changed |= ui
                                .add(
                                    egui::Slider::new(&mut settings.rim_distance, 3.0..=50.0)
                                        .text("Rim distance"),
                                )
                                .changed();
                        }

                        if changed {
                            settings.lighting_dirty = true;
                        }
                    });
                    crate::ui_style::collapsible_group(ui, "Element visibility", false, |ui| {
                        refresh_element_filter_cache(&mut xyz_buf, &mol);
                        let extras = &xyz_buf.element_filter_extras;
                        if mol.atoms.is_empty() {
                            ui.weak("No atoms loaded.");
                            return;
                        }

                        let mut changed = false;
                        ui.horizontal_wrapped(|ui| {
                            if ui.button("Show all").clicked_once() {
                                for block in [
                                    &ELEMENT_FILTER_S_BLOCK[..],
                                    &ELEMENT_FILTER_D_BLOCK[..],
                                    &ELEMENT_FILTER_P_BLOCK[..],
                                    &ELEMENT_FILTER_F_BLOCK[..],
                                ] {
                                    for sym in block {
                                        settings
                                            .element_visibility
                                            .insert((*sym).to_string(), true);
                                    }
                                }
                                for sym in extras {
                                    settings.element_visibility.insert(sym.clone(), true);
                                }
                                changed = true;
                            }
                            if ui.button("Hide all").clicked_once() {
                                for block in [
                                    &ELEMENT_FILTER_S_BLOCK[..],
                                    &ELEMENT_FILTER_D_BLOCK[..],
                                    &ELEMENT_FILTER_P_BLOCK[..],
                                    &ELEMENT_FILTER_F_BLOCK[..],
                                ] {
                                    for sym in block {
                                        settings
                                            .element_visibility
                                            .insert((*sym).to_string(), false);
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
                                let render_block =
                                    |ui: &mut egui::Ui,
                                     block: &[&str],
                                     settings: &mut MolSettings,
                                     changed: &mut bool| {
                                        let cols = 4;
                                        let rows = (block.len() + cols - 1) / cols;
                                        for r in 0..rows {
                                            ui.horizontal_wrapped(|ui| {
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
                                render_block(
                                    ui,
                                    &ELEMENT_FILTER_P_BLOCK,
                                    &mut settings,
                                    &mut changed,
                                );
                                ui.separator();
                                ui.label("S-block");
                                render_block(
                                    ui,
                                    &ELEMENT_FILTER_S_BLOCK,
                                    &mut settings,
                                    &mut changed,
                                );
                                ui.separator();
                                ui.label("D-block");
                                render_block(
                                    ui,
                                    &ELEMENT_FILTER_D_BLOCK,
                                    &mut settings,
                                    &mut changed,
                                );
                                ui.separator();
                                ui.label("F-block");
                                render_block(
                                    ui,
                                    &ELEMENT_FILTER_F_BLOCK,
                                    &mut settings,
                                    &mut changed,
                                );

                                if !extras.is_empty() {
                                    ui.separator();
                                    ui.label("Other");
                                    let other: Vec<&str> =
                                        extras.iter().map(|s| s.as_str()).collect();
                                    render_block(ui, &other, &mut settings, &mut changed);
                                }
                            });

                        if changed {
                            settings.geometry_dirty = true;
                        }
                    });

                    ui.add_space(8.0);
                },
            );

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
                    if crate::export_image::ui::panel(
                        ui,
                        &mut export_settings,
                        &mut export_select_area,
                    ) {
                        let (Ok((_cam, proj, cam_xform)), Ok(win)) =
                            (q_cam.single(), windows.single())
                        else {
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
                },
            );
        };

        if windowed {
            // A permanent, narrow rail on the left; every tool lives in its
            // own floating window that the rail toggles, so any number can be
            // open at once and dragged where the user wants them.
            let rail = egui::Panel::left("tool_rail")
                .exact_size(50.0)
                .resizable(false)
                .show(&mut viewport_ui, |ui| {
                    icon_rail(ui, &mut open_windows);
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

    // The editor opens floating in either application layout. Docking is an
    // explicit choice in the editor, independent of the other tool windows.
    {
        let mut editor_launcher_offset = if ui_layout.builder_open
            && crate::molecule_builder::builder_ui::editor_docked(&ctx)
        {
            ctx.data_mut(|data| {
                data.get_temp::<f32>(egui::Id::new("modern_editor_docked_width"))
                    .unwrap_or(420.0)
            })
        } else {
            0.0
        };
        if qc_panel.open && crate::qchem_panel::ui::is_docked(&ctx) {
            editor_launcher_offset += crate::qchem_panel::ui::docked_width(&ctx);
        }
        egui::Area::new(egui::Id::new("builder_launcher"))
            .anchor(
                egui::Align2::RIGHT_TOP,
                egui::vec2(-12.0 - editor_launcher_offset, 12.0),
            )
            .show(&ctx, |ui| {
                egui::Frame::popup(ui.style())
                    .corner_radius(egui::CornerRadius::same(18))
                    .show(ui, |ui| {
                        // Two launchers in one group, the general
                        // quantum-chemistry panel beside the editor. They sit
                        // together because both are entry points to doing
                        // something *to* the structure, rather than tools for
                        // looking at one, which is what the left rail is.
                        ui.horizontal(|ui| {
                            let button = crate::ui_icons::IconButton::new(
                                crate::ui_icons::Icon::QuantumChemistry,
                            )
                            .size(30.0)
                            .selected(qc_panel.open);
                            if ui
                                .add(button)
                                .on_hover_text(
                                    "Quantum chemistry — run any job on this structure",
                                )
                                .clicked_once()
                            {
                                qc_panel.open = !qc_panel.open;
                            }

                            let button = crate::ui_icons::IconButton::new(
                                crate::ui_icons::Icon::MoleculeEditor,
                            )
                            .size(30.0)
                            .selected(ui_layout.builder_open);
                            if ui
                                .add(button)
                                .on_hover_text("Molecule editor / Z-matrix")
                                .clicked_once()
                            {
                                ui_layout.builder_open = !ui_layout.builder_open;
                            }
                        });
                    });
            });

        let mut builder_open = ui_layout.builder_open;
        let mut draw_editor = |ui: &mut egui::Ui| {
            let change = builder_ui_contents(
                ui, &mut editor_rotate_state, &mut zmat_state, &mut mol,
                &mut settings, &mut structure_history, &mut traj,
            );
            if change.happened() {
                xtb_freq_panel_state.selected_mode = None;
                xyz_buf.current_file = None;
                // Only a genuinely different structure re-frames the camera.
                // An edit in place leaves the view alone.
                ev_changed.write(MoleculeChanged::parse_xyz(change.recenter()));
            }
        };
        if builder_open && crate::molecule_builder::builder_ui::editor_docked(&ctx) {
            let panel = egui::Panel::right("modern_molecular_editor")
                .default_size(420.0)
                .min_size(380.0)
                .max_size(620.0)
                .resizable(true)
                .show(&mut viewport_ui, &mut draw_editor);
            ctx.data_mut(|data| {
                data.insert_temp(egui::Id::new("modern_editor_docked_width"), panel.response.rect.width());
            });
            window_rects.push(panel.response.rect);
        } else {
            let builder_window = molecule_editor_window(&ctx, &mut builder_open, draw_editor);
            if let Some(window) = &builder_window {
                window_rects.push(window.response.rect);
            }
        }
        if crate::molecule_builder::builder_ui::take_close_request(&ctx) {
            builder_open = false;
        }
        ui_layout.builder_open = builder_open;
    }

    // The calculator opens from the launcher beside the editor. Its modern
    // form can float or dock; finished results use their own summary window.
    if qc_panel.open {
        if let Some(rect) = crate::qchem_panel::ui::qc_panel(
            ctx,
            &mut viewport_ui,
            &mut qc_panel,
            &mut qc_task,
            &mut xtb_panel_state,
            &mol,
        ) {
            window_rects.push(rect);
        }
    }
    if let Some(rect) = crate::qchem_panel::ui::input_editor_window(
        ctx, &mut qc_panel, &mut qc_task, &xtb_panel_state, &mol,
    ) {
        window_rects.push(rect);
    }
    crate::qchem_panel::ui::result_windows(ctx, &mut xtb_panel_state, &xtb_task);
    if qc_panel.summary_open {
        if let Some(output) = &qc_task.last {
            let duration = qc_task.last_duration;
            let mut open = qc_panel.summary_open;
            crate::qchem_panel::summary::summary_window(
                ctx,
                &mut open,
                &mut qc_panel.summary_confirm_close,
                output,
                duration,
                &mut xtb_panel_state,
                &xtb_task,
            );
            qc_panel.summary_open = open;
        } else {
            qc_panel.summary_open = false;
        }
    }

    if let Some(rect) = crate::ts_generation::help::window(ctx, &mut ts_generation.help_open) {
        window_rects.push(rect);
    }
    if let Some(rect) = crate::ts_generation::plot::window(ctx, &mut ts_generation) {
        window_rects.push(rect);
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
        let physical = crate::scene::axis_gizmo_canvas_size(
            window, panel_regions.free_viewport, &panel_regions.windows,
        );
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

pub(crate) fn compute_centroid(points: &[Vec3]) -> Option<Vec3> {
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

/// The Molecule Editor window.
///
/// The editor owns its scrolling so its header, selection and footer stay visible.
pub(crate) fn molecule_editor_window(
    ctx: &egui::Context,
    open: &mut bool,
    add_contents: impl FnOnce(&mut egui::Ui),
) -> Option<egui::InnerResponse<Option<()>>> {
    const EDITOR_WIDTH: f32 = 420.0;
    const EDITOR_HEIGHT: f32 = 760.0;
    const EDITOR_MARGIN: f32 = 16.0;
    let viewport = ctx.viewport_rect();
    let height = EDITOR_HEIGHT.min(viewport.height() - 2.0 * EDITOR_MARGIN);
    egui::Window::new("Molecule Editor")
        .id(egui::Id::new("molecule_editor_modern_v1"))
        .open(open)
        .resizable(true)
        .default_size([EDITOR_WIDTH, height])
        .min_width(380.0)
        .min_height(440.0)
        .default_pos([
            (viewport.right() - EDITOR_WIDTH - EDITOR_MARGIN).max(viewport.left()),
            viewport.top() + EDITOR_MARGIN,
        ])
        .show(ctx, add_contents)
}

pub(crate) fn apply_xyz_text(
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

/// Read an XYZ file from a path we already have, single-frame or trajectory.
///
/// The Structure panel's button, the Trajectory panel's button and a filename
/// given on the command line all end up here, so the multi-frame case is
/// decided in one place: more than one frame becomes a trajectory rather than
/// a warning that the other frames were thrown away.
///
/// Returns how many frames were loaded.
pub(crate) fn load_xyz_from_path(
    path: &std::path::Path,
    xyz_buf: &mut XyzBuffer,
    traj: &mut trajectory::TrajectoryState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
    cam: &mut crate::camera::OrbitCamera,
    ev_changed: &mut MessageWriter<MoleculeChanged>,
) -> Result<usize, String> {
    let text = fs::read_to_string(path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;

    let frames = trajectory::parse_multi_xyz(&text)
        .map_err(|e| format!("{}: {e}", path.display()))?;
    if frames.is_empty() {
        return Err(format!("{}: no atoms in this file", path.display()));
    }

    xyz_buf.last_dir = path.parent().map(|p| p.to_path_buf());
    xyz_buf.current_file = Some(path.to_path_buf());
    xyz_buf.text = text.clone();
    xyz_buf.warning = None;

    if frames.len() > 1 {
        traj.last_dir = path.parent().map(|p| p.to_path_buf());
        traj.current_file = Some(path.to_path_buf());
        traj.frames = frames;
        traj.current_frame = 0;
        traj.playing = false;
        traj.accum = 0.0;
        traj.last_applied = None;
        traj.overlay_dirty = true;
        let n = traj.frames.len();
        trajectory::apply_current_frame(traj, mol, settings, ev_changed, true, Some(cam));
        // `apply_current_frame` only reports a coordinate update -- it runs for
        // every playback frame too. Loading a file replaces the structure, so
        // say so, or stale measurements carry over. `false`: the call above
        // already centred the camera.
        ev_changed.write(MoleculeChanged::parse_xyz(false));
        Ok(n)
    } else {
        traj.clear_for_structure();
        apply_xyz_text(&text, mol, cam, ev_changed);
        Ok(1)
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
