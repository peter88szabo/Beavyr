//! Runs `xtb --opt` in the background and hands back its optimization
//! trajectory. Ported from `Tautomer`'s `qchem_interface/xtbrun.rs`, which
//! already solved cancellation and non-blocking process management; adapted
//! here to Beavyr's own `AsyncComputeTaskPool` idiom (see
//! `src/orbitals/mod.rs`) rather than importing Tautomer's plugin shape.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use bevy::prelude::*;
use bevy::tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task};

use crate::events::MoleculeChanged;
use crate::molecule::Molecule;
use crate::settings::MolSettings;
use crate::trajectory::{self, TrajectoryState};

use bevy_egui::egui;

use super::config;
use super::program::QcProgram;
use super::valence::validate_electronic_state;

/// Charge/multiplicity inputs and the xTB path field, plus the last
/// validation result -- everything the panel needs to redraw itself and
/// nothing that belongs to the optimization run itself (that lives in
/// `XtbOptimizationTask`).
#[derive(Resource)]
pub struct XtbPanelState {
    pub charge: i32,
    pub multiplicity: i32,
    /// Which external program runs the optimization.
    pub program: QcProgram,
    /// One path per program, so switching programs switches paths rather than
    /// overwriting the one the user set for the other.
    pub paths: Vec<(QcProgram, String)>,
    /// Whether the energy-vs-iteration plot window is open. Kept here rather
    /// than in `XtbOptimizationTask` because it is purely a UI toggle, not
    /// part of the run itself.
    pub energy_plot_open: bool,
    /// Valence warnings are only shown once the user has actually tried to
    /// optimize with the current charge/multiplicity -- not the instant the
    /// panel is drawn for whatever structure happens to be on screen, which
    /// read as unsolicited noise for a perfectly ordinary molecule (an
    /// aromatic ring, before the bond-order fix, or any atom the heuristic
    /// is simply unsure about). A hard `ElectronicStateError` is different:
    /// it explains why the Optimize button is disabled, so it stays visible
    /// immediately.
    pub show_warnings: bool,
}

impl XtbPanelState {
    /// The path currently entered for the selected program.
    pub fn path(&self) -> &str {
        self.paths
            .iter()
            .find(|(p, _)| *p == self.program)
            .map(|(_, path)| path.as_str())
            .unwrap_or("")
    }

    /// Editable access to the selected program's path, so the one text field
    /// always writes to the program it is showing.
    pub fn path_mut(&mut self) -> &mut String {
        let program = self.program;
        if !self.paths.iter().any(|(p, _)| *p == program) {
            self.paths.push((program, String::new()));
        }
        self.paths
            .iter_mut()
            .find(|(p, _)| *p == program)
            .map(|(_, path)| path)
            .expect("just inserted if missing")
    }
}

impl Default for XtbPanelState {
    fn default() -> Self {
        Self {
            program: QcProgram::default(),
            // A saved path per program; xTB falls back to a sensible guess so
            // a first run needs no configuration on a normal installation.
            paths: QcProgram::ALL
                .iter()
                .map(|&p| {
                    let saved = config::load_path(p);
                    let path = match (saved.is_empty(), p) {
                        (true, QcProgram::Xtb) => default_xtb_executable(),
                        _ => saved,
                    };
                    (p, path)
                })
                .collect(),
            charge: 0,
            multiplicity: 1,
            energy_plot_open: false,
            show_warnings: false,
        }
    }
}

static XTB_RUN_SEQUENCE: AtomicU64 = AtomicU64::new(0);

/// One optimizer iteration's energy, for the "how did this run behave" plot.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyHistoryPoint {
    pub iteration: u32,
    pub energy_hartree: f64,
}

/// A successful `--opt` run: the whole trajectory, its final energy, and the
/// per-iteration energy history (from xTB's own stdout) for the energy plot.
pub struct XtbOptimizationOutput {
    pub trajectory_text: String,
    pub final_energy_hartree: Option<f64>,
    pub energy_history: Vec<EnergyHistoryPoint>,
}

/// A run that did not finish successfully. It still carries whatever energy
/// history xTB had printed before it stopped, precisely so a failed run can
/// still be inspected on the same plot -- how far did it get, and was it
/// still descending or already diverging?
pub struct XtbOptimizationFailure {
    pub message: String,
    pub energy_history: Vec<EnergyHistoryPoint>,
}

/// Background xTB optimization, plus everything needed to cancel it cleanly.
#[derive(Resource, Default)]
pub struct XtbOptimizationTask {
    task: Option<Task<Result<XtbOptimizationOutput, XtbOptimizationFailure>>>,
    cancel: Option<Arc<AtomicBool>>,
    child: Option<Arc<Mutex<Option<Child>>>>,
    run_dir: Option<PathBuf>,
    started_at: Option<Instant>,
    /// Set once the task finishes, success or failure, for the panel to show.
    pub last_message: Option<String>,
    pub last_is_error: bool,
    /// The finished run's per-iteration energies, success or failure alike,
    /// kept until the next run starts so the plot button stays available.
    pub last_energy_history: Vec<EnergyHistoryPoint>,
}

impl XtbOptimizationTask {
    pub fn is_running(&self) -> bool {
        self.task.is_some()
    }

    /// How long the current run has been going, for a live "is this stuck?"
    /// readout -- `None` once nothing is running.
    pub fn elapsed(&self) -> Option<Duration> {
        self.started_at.map(|t| t.elapsed())
    }

    /// Where the run is writing its scratch files, so the panel can read
    /// xTB's own live stdout for progress. `None` once nothing is running.
    pub fn run_dir(&self) -> Option<&Path> {
        self.run_dir.as_deref()
    }

    pub fn cancel_now(&self) {
        if let Some(cancel) = &self.cancel {
            cancel.store(true, Ordering::Relaxed);
        }
        if let Some(child) = &self.child {
            if let Ok(mut child) = child.lock() {
                if let Some(mut process) = child.take() {
                    let _ = process.kill();
                    let _ = process.wait();
                }
            }
        }
    }

    /// Start optimizing `atoms`/`pos` with xTB. `charge` and `uhf` (unpaired
    /// electrons, i.e. `multiplicity - 1`) are assumed already validated by
    /// `qchem_interfaces::valence::validate_electronic_state` -- this
    /// function does not re-check them.
    pub fn start(
        &mut self,
        program: QcProgram,
        binary: &Path,
        atoms: &[String],
        pos: &[Vec3],
        charge: i32,
        uhf: i32,
        multiplicity: i32,
    ) {
        if self.is_running() {
            return;
        }
        let input_xyz = write_xyz_string(atoms, pos);
        let binary = binary.to_path_buf();
        let cancel = Arc::new(AtomicBool::new(false));
        let child_slot = Arc::new(Mutex::new(None));
        let run_dir = match create_xtb_run_dir(&xtb_scratch_dir()) {
            Ok(dir) => dir,
            Err(err) => {
                self.last_message = Some(err);
                self.last_is_error = true;
                return;
            }
        };

        let task_cancel = cancel.clone();
        let task_child = child_slot.clone();
        let task_dir = run_dir.clone();
        let task = AsyncComputeTaskPool::get().spawn(async move {
            run_optimize_cancellable(
                program,
                &binary,
                &task_dir,
                &input_xyz,
                charge,
                uhf,
                multiplicity,
                &task_cancel,
                &task_child,
            )
            .map_err(|message| XtbOptimizationFailure {
                energy_history: energy_history_for(program, &task_dir),
                message,
            })
        });

        self.task = Some(task);
        self.cancel = Some(cancel);
        self.child = Some(child_slot);
        self.run_dir = Some(run_dir);
        self.started_at = Some(Instant::now());
        self.last_message = None;
        self.last_energy_history.clear();
    }
}

/// Poll the running task and, once it finishes, apply the result: load a
/// successful trajectory (auto-playing it, per the user's own choice) or
/// surface the failure message. The displayed structure is left untouched on
/// failure.
pub fn poll_xtb_optimization(
    mut xtb_task: ResMut<XtbOptimizationTask>,
    mut traj: ResMut<TrajectoryState>,
    mut mol: ResMut<Molecule>,
    mut settings: ResMut<MolSettings>,
    mut cam: Query<&mut crate::camera::OrbitCamera>,
    mut ev_changed: MessageWriter<MoleculeChanged>,
) {
    let Some(task) = xtb_task.task.as_mut() else {
        return;
    };
    let Some(result) = block_on(future::poll_once(task)) else {
        return;
    };
    xtb_task.task = None;
    xtb_task.cancel = None;
    xtb_task.child = None;
    xtb_task.started_at = None;
    let run_dir = xtb_task.run_dir.take();

    match result {
        Ok(output) => {
            if let Some(dir) = &run_dir {
                discard_xtb_run_dir(dir);
            }
            match trajectory::parse_multi_xyz(&output.trajectory_text) {
                Ok(frames) => {
                    let frame_count = frames.len();
                    traj.frames = frames;
                    traj.current_frame = 0;
                    traj.accum = 0.0;
                    traj.last_applied = None;
                    traj.overlay_dirty = true;
                    let mut cam = cam.iter_mut().next();
                    trajectory::apply_current_frame(
                        &mut traj,
                        &mut mol,
                        &mut settings,
                        &mut ev_changed,
                        true,
                        cam.as_deref_mut(),
                    );
                    ev_changed.write(MoleculeChanged::parse_xyz(false));
                    // Auto-play once, per the user's own choice: switch to
                    // the Trajectory tool and press play for them. Forced to
                    // `Once` regardless of whatever mode was left over from a
                    // previous trajectory -- watching an optimization loop or
                    // ping-pong back and forth makes no physical sense.
                    traj.mode = trajectory::PlaybackMode::Once;
                    traj.direction = 1;
                    traj.playing = true;

                    let energy_text = output
                        .final_energy_hartree
                        .map(|e| format!(", final energy {e:.6} Eh"))
                        .unwrap_or_default();
                    xtb_task.last_message = Some(format!(
                        "xTB optimization succeeded ({frame_count} frames{energy_text}). \
                         Open the Trajectory tool to watch or replay it."
                    ));
                    xtb_task.last_is_error = false;
                }
                Err(err) => {
                    xtb_task.last_message =
                        Some(format!("xTB finished, but its trajectory could not be read: {err}"));
                    xtb_task.last_is_error = true;
                }
            }
            xtb_task.last_energy_history = output.energy_history;
        }
        Err(failure) => {
            // A failed run's scratch directory is kept (not discarded) so
            // `xtb.stdout`/`xtb.stderr` remain inspectable -- matching
            // Tautomer's own rule (confirmed at its `ui.rs:2561` call site).
            //
            // The run may still have written a partial `xtbopt.log` before
            // it stopped. That is loaded into the trajectory player too --
            // paused, not auto-played, since it never converged -- so the
            // user can scrub to any attempted step and use it as a new
            // starting geometry for another try, exactly as they would with
            // a successful run's trajectory.
            let mut partial_note = String::new();
            if let Some(dir) = &run_dir {
                if let Ok(text) = fs::read_to_string(dir.join("xtbopt.log")) {
                    if let Ok(frames) = trajectory::parse_multi_xyz(&text) {
                        if !frames.is_empty() {
                            let frame_count = frames.len();
                            traj.frames = frames;
                            traj.current_frame = frame_count - 1;
                            traj.mode = trajectory::PlaybackMode::Once;
                            traj.direction = 1;
                            traj.playing = false;
                            traj.accum = 0.0;
                            traj.last_applied = None;
                            traj.overlay_dirty = true;
                            let mut cam = cam.iter_mut().next();
                            trajectory::apply_current_frame(
                                &mut traj,
                                &mut mol,
                                &mut settings,
                                &mut ev_changed,
                                true,
                                cam.as_deref_mut(),
                            );
                            ev_changed.write(MoleculeChanged::parse_xyz(false));
                            partial_note = format!(
                                " A partial trajectory ({frame_count} step(s)) was loaded \
                                 into the Trajectory tool -- scrub to any step and use it \
                                 as a new starting geometry."
                            );
                        }
                    }
                }
            }
            xtb_task.last_message = Some(match &run_dir {
                Some(dir) => format!("{}{partial_note} (see {})", failure.message, dir.display()),
                None => format!("{}{partial_note}", failure.message),
            });
            xtb_task.last_is_error = true;
            xtb_task.last_energy_history = failure.energy_history;
        }
    }
}

/// Draws the geometry-optimization section: charge/multiplicity, the xTB
/// path, validation feedback, and the Optimize/Cancel controls. Called from
/// `ui.rs`, matching how `molecule_builder::builder_ui_panel` is a free
/// function the main panel invokes rather than a system of its own.
pub fn xtb_optimization_panel(
    ui: &mut egui::Ui,
    panel_state: &mut XtbPanelState,
    task: &mut XtbOptimizationTask,
    mol: &Molecule,
    // Whether a trajectory or mode animation is playing. While one is, the
    // molecule on screen is a frame of it, not the structure the user means
    // to optimize.
    animating: bool,
) {
    let running = task.is_running();

    ui.horizontal(|ui| {
        ui.label("Charge");
        if ui
            .add_enabled(
                !running,
                egui::DragValue::new(&mut panel_state.charge).range(-10..=10),
            )
            .changed()
        {
            panel_state.show_warnings = false;
        }
        ui.label("Multiplicity");
        if ui
            .add_enabled(
                !running,
                egui::DragValue::new(&mut panel_state.multiplicity).range(1..=10),
            )
            .changed()
        {
            panel_state.show_warnings = false;
        }
        if ui
            .add_enabled(!running, egui::Button::new("Reset to defaults"))
            .clicked()
        {
            panel_state.charge = 0;
            panel_state.multiplicity = 1;
            panel_state.show_warnings = false;
        }
    });

    ui.horizontal(|ui| {
        ui.label("Program");
        egui::ComboBox::from_id_salt("opt_program")
            .selected_text(panel_state.program.label())
            .show_ui(ui, |ui| {
                for program in QcProgram::ALL {
                    ui.selectable_value(&mut panel_state.program, program, program.label());
                }
            });
    });
    ui.horizontal(|ui| {
        let program = panel_state.program;
        ui.label(format!("{} path", program.label()));
        let response = ui.add_enabled(
            !running,
            egui::TextEdit::singleline(panel_state.path_mut())
                .desired_width(220.0)
                .hint_text(program.binary_hint()),
        );
        if response.lost_focus() {
            config::save_path(program, panel_state.path());
        }
        if ui
            .add_enabled(!running, egui::Button::new("Browse…"))
            .clicked()
        {
            let mut dlg = rfd::FileDialog::new();
            if let Some(dir) = std::path::Path::new(panel_state.path())
                .parent()
                .filter(|p| p.is_dir())
            {
                dlg = dlg.set_directory(dir);
            }
            if let Some(path) = dlg.pick_file() {
                *panel_state.path_mut() = path.display().to_string();
                config::save_path(program, panel_state.path());
            }
        }
    });

    // A playing frame is a distorted geometry, so its bond lengths -- and
    // every valence conclusion drawn from them -- are about that frame, not
    // about the molecule. Checking it would produce warnings that change from
    // frame to frame and mean nothing.
    let validation = validate_electronic_state(mol, panel_state.charge, panel_state.multiplicity);
    let (uhf, warnings) = match &validation {
        Ok((uhf, warnings)) => (Some(*uhf), warnings.as_slice()),
        Err(_) => (None, &[][..]),
    };

    if animating {
        ui.weak("Playback is running \u{2014} stop it to optimize this structure.");
    } else if let Err(err) = &validation {
        // A hard error explains why Optimize is disabled, so it is shown
        // immediately rather than waiting for an attempt that the disabled
        // button cannot register anyway.
        ui.colored_label(egui::Color32::from_rgb(220, 80, 80), err.to_string());
    } else if panel_state.show_warnings {
        for warning in warnings {
            ui.colored_label(
                egui::Color32::from_rgb(210, 160, 40),
                format!("Atom {}: {}", warning.atom_index + 1, warning.message),
            );
        }
    }

    ui.horizontal(|ui| {
        let can_optimize = !running && !mol.atoms.is_empty() && uhf.is_some() && !animating;
        if ui
            .add_enabled(can_optimize, egui::Button::new("Optimize"))
            .on_disabled_hover_text(if animating {
                "Stop the animation first: the structure on screen is a frame of it."
            } else {
                "Nothing to optimize with the current structure and electronic state."
            })
            .clicked()
        {
            panel_state.show_warnings = true;
            let program = panel_state.program;
            let entered = panel_state.path().to_string();
            match (uhf, resolve_program_executable(program, &entered)) {
                (Some(uhf), Ok(binary)) => {
                    config::save_path(program, &entered);
                    task.start(
                        program,
                        &binary,
                        &mol.atoms,
                        &mol.pos,
                        panel_state.charge,
                        uhf,
                        panel_state.multiplicity,
                    );
                }
                (_, Err(err)) => {
                    task.last_message = Some(err);
                    task.last_is_error = true;
                }
                _ => {}
            }
        }
        if running && ui.button("Cancel").clicked() {
            task.cancel_now();
        }
    });

    if running {
        ui.horizontal(|ui| {
            let elapsed = task.elapsed().unwrap_or_default();
            ui.weak(format!("Running… {}", format_elapsed(elapsed)));
            let progress = task
                .run_dir()
                .and_then(|dir| std::fs::read_to_string(dir.join("xtb.stdout")).ok())
                .and_then(|text| parse_xtb_progress(&text));
            if let Some(progress) = progress {
                let energy = progress
                    .energy_hartree
                    .map(|e| format!(", energy {e:.6} Eh"))
                    .unwrap_or_default();
                let gnorm = progress
                    .gradient_norm
                    .map(|g| format!(", gradient norm {g:.5} Eh/\u{03b1}"))
                    .unwrap_or_default();
                ui.weak(format!("— iteration {}{energy}{gnorm}", progress.iteration));
            } else {
                ui.weak("— starting…");
            }
        });
    }

    if let Some(message) = &task.last_message {
        let color = if task.last_is_error {
            egui::Color32::from_rgb(220, 80, 80)
        } else {
            egui::Color32::from_rgb(80, 180, 100)
        };
        ui.colored_label(color, message);
    }

    if !task.last_energy_history.is_empty() {
        let label = if panel_state.energy_plot_open {
            "Hide Energy Plot"
        } else {
            "Show Energy Plot"
        };
        if ui.button(label).clicked() {
            panel_state.energy_plot_open = !panel_state.energy_plot_open;
        }
    }

    let ctx = ui.ctx().clone();
    energy_history_window(
        &ctx,
        &mut panel_state.energy_plot_open,
        &task.last_energy_history,
        task.last_is_error,
    );
}

/// A simple line plot of energy vs. iteration, in its own window -- similar
/// in spirit to Molden's own optimization-energy plot. Hand-drawn on
/// `egui::Painter` rather than pulling in a plotting crate: a straight
/// polyline through a handful of points does not need one.
fn energy_history_window(
    ctx: &egui::Context,
    open: &mut bool,
    history: &[EnergyHistoryPoint],
    run_failed: bool,
) {
    if !*open || history.is_empty() {
        return;
    }
    egui::Window::new("xTB Optimization Energy")
        .open(open)
        .resizable(true)
        .default_size([420.0, 300.0])
        .show(ctx, |ui| {
            let status = if run_failed {
                "This run did not converge; the plot shows the energy up to                  where it stopped."
            } else {
                "Energy at each optimizer iteration."
            };
            ui.weak(status);

            let min_e = history
                .iter()
                .map(|p| p.energy_hartree)
                .fold(f64::INFINITY, f64::min);
            let max_e = history
                .iter()
                .map(|p| p.energy_hartree)
                .fold(f64::NEG_INFINITY, f64::max);
            let min_it = history.iter().map(|p| p.iteration).min().unwrap_or(0);
            let max_it = history.iter().map(|p| p.iteration).max().unwrap_or(1);
            let e_span = (max_e - min_e).max(1.0e-9);
            let it_span = (max_it - min_it).max(1) as f64;

            const AXIS_LABEL_SIZE: f32 = 13.0;
            const AXIS_TITLE_SIZE: f32 = 16.0;
            const TICK_COUNT: u32 = 5;

            // Both dimensions must track the window's own size, not just
            // width -- a fixed height here is exactly why the window could be
            // resized horizontally but not vertically: the plot itself never
            // grew to fill the extra space.
            let desired_size = egui::vec2(
                ui.available_width().max(200.0),
                ui.available_height().max(120.0),
            );
            let (rect, response) = ui.allocate_exact_size(desired_size, egui::Sense::hover());
            let painter = ui.painter_at(rect);
            let plot_color = ui.visuals().text_color();
            let grid_color = plot_color.gamma_multiply(0.25);
            let line_color = if run_failed {
                egui::Color32::from_rgb(220, 100, 60)
            } else {
                egui::Color32::from_rgb(90, 160, 220)
            };

            // Margins wide enough for tick labels plus a rotated y-axis
            // title on the left, and an x-axis title along the bottom.
            let plot_rect = egui::Rect::from_min_max(
                rect.min + egui::vec2(78.0, 10.0),
                rect.max - egui::vec2(10.0, 46.0),
            );

            let to_screen = |iteration: u32, energy: f64| -> egui::Pos2 {
                let x = plot_rect.left()
                    + ((iteration - min_it) as f64 / it_span) as f32 * plot_rect.width();
                // Screen y grows downward, so the highest energy is at the top.
                let y = plot_rect.bottom()
                    - ((energy - min_e) / e_span) as f32 * plot_rect.height();
                egui::pos2(x, y)
            };

            // Gridlines with value labels, evenly spaced across each axis --
            // this is what turns a bare outline into a readable set of axes.
            for step in 0..=TICK_COUNT {
                let t = step as f64 / TICK_COUNT as f64;

                let energy = max_e - t * e_span; // top of the axis is the highest energy
                let y = plot_rect.top() + t as f32 * plot_rect.height();
                painter.line_segment(
                    [egui::pos2(plot_rect.left(), y), egui::pos2(plot_rect.right(), y)],
                    egui::Stroke::new(1.0, grid_color),
                );
                painter.text(
                    egui::pos2(plot_rect.left() - 6.0, y),
                    egui::Align2::RIGHT_CENTER,
                    format!("{energy:.6}"),
                    egui::FontId::monospace(AXIS_LABEL_SIZE),
                    plot_color,
                );

                let iteration = min_it as f64 + t * it_span;
                let x = plot_rect.left() + t as f32 * plot_rect.width();
                painter.line_segment(
                    [egui::pos2(x, plot_rect.top()), egui::pos2(x, plot_rect.bottom())],
                    egui::Stroke::new(1.0, grid_color),
                );
                painter.text(
                    egui::pos2(x, plot_rect.bottom() + 6.0),
                    egui::Align2::CENTER_TOP,
                    format!("{iteration:.0}"),
                    egui::FontId::monospace(AXIS_LABEL_SIZE),
                    plot_color,
                );
            }

            painter.rect_stroke(
                plot_rect,
                0.0,
                egui::Stroke::new(1.5, plot_color.gamma_multiply(0.6)),
                egui::StrokeKind::Outside,
            );

            let points: Vec<egui::Pos2> = history
                .iter()
                .map(|p| to_screen(p.iteration, p.energy_hartree))
                .collect();
            if points.len() >= 2 {
                painter.add(egui::Shape::line(points.clone(), egui::Stroke::new(2.0, line_color)));
            }
            for &p in &points {
                painter.circle_filled(p, 2.5, line_color);
            }

            // Highlight and label whichever point the cursor is nearest to,
            // as long as it is actually close -- so hovering empty space
            // inside the plot shows nothing.
            if let Some(hover_pos) = response.hover_pos() {
                if let Some((idx, dist)) = points
                    .iter()
                    .enumerate()
                    .map(|(i, p)| (i, p.distance(hover_pos)))
                    .min_by(|a, b| a.1.total_cmp(&b.1))
                {
                    const HOVER_RADIUS: f32 = 14.0;
                    if dist <= HOVER_RADIUS {
                        let point = history[idx];
                        let highlight_color = ui.visuals().strong_text_color();
                        painter.circle_stroke(
                            points[idx],
                            5.0,
                            egui::Stroke::new(2.0, highlight_color),
                        );

                        // Drawn directly at the cursor rather than via
                        // `on_hover_text_at_pointer`: egui's built-in tooltip
                        // only appears once the pointer has been *still* for
                        // its configured delay, which never fires while
                        // sweeping across the curve to compare points -- the
                        // exact way this is meant to be used.
                        let label = format!(
                            "iteration {}\nenergy {:.8} Eh",
                            point.iteration, point.energy_hartree
                        );
                        let galley = painter.layout_no_wrap(
                            label,
                            egui::FontId::proportional(AXIS_LABEL_SIZE),
                            highlight_color,
                        );
                        let padding = egui::vec2(6.0, 4.0);
                        let box_size = galley.size() + padding * 2.0;
                        let mut box_pos = points[idx] + egui::vec2(12.0, -box_size.y - 12.0);
                        // Keep the label inside the window even near an edge.
                        box_pos.x = box_pos.x.clamp(rect.left(), rect.right() - box_size.x);
                        box_pos.y = box_pos.y.clamp(rect.top(), rect.bottom() - box_size.y);
                        let box_rect = egui::Rect::from_min_size(box_pos, box_size);

                        painter.rect_filled(box_rect, 4.0, ui.visuals().extreme_bg_color);
                        painter.rect_stroke(
                            box_rect,
                            4.0,
                            egui::Stroke::new(1.0, plot_color.gamma_multiply(0.5)),
                            egui::StrokeKind::Outside,
                        );
                        painter.galley(box_rect.min + padding, galley, highlight_color);
                    }
                }
            }

            // Axis titles: "Iteration" along the bottom, "Energy (Eh)"
            // rotated upright along the left edge.
            painter.text(
                egui::pos2(plot_rect.center().x, rect.bottom() - 4.0),
                egui::Align2::CENTER_BOTTOM,
                "Iteration",
                egui::FontId::proportional(AXIS_TITLE_SIZE),
                plot_color,
            );
            let y_title_galley = painter.layout_no_wrap(
                "Energy (Eh)".to_string(),
                egui::FontId::proportional(AXIS_TITLE_SIZE),
                plot_color,
            );
            let y_title_pos = egui::pos2(
                rect.left() + y_title_galley.size().y,
                plot_rect.center().y + y_title_galley.size().x / 2.0,
            );
            painter.add(egui::epaint::TextShape {
                angle: -std::f32::consts::FRAC_PI_2,
                ..egui::epaint::TextShape::new(y_title_pos, y_title_galley, plot_color)
            });
        });
}

/// `mm:ss`, or `h:mm:ss` past an hour -- plain enough to read at a glance
/// while deciding whether a long run is still alive.
pub(super) fn format_elapsed(elapsed: Duration) -> String {
    let total_secs = elapsed.as_secs();
    let (hours, rem) = (total_secs / 3600, total_secs % 3600);
    let (minutes, seconds) = (rem / 60, rem % 60);
    if hours > 0 {
        format!("{hours}:{minutes:02}:{seconds:02}")
    } else {
        format!("{minutes:02}:{seconds:02}")
    }
}

pub(super) fn write_xyz_string(atoms: &[String], pos: &[Vec3]) -> String {
    let mut out = format!("{}\n\n", atoms.len());
    for (symbol, p) in atoms.iter().zip(pos) {
        out.push_str(&format!(
            "{symbol:<3} {:>15.8} {:>15.8} {:>15.8}\n",
            p.x, p.y, p.z
        ));
    }
    out
}

pub fn default_xtb_executable() -> String {
    std::env::var("BEAVYR_XTB_PATH").unwrap_or_else(|_| "xtb".to_string())
}

/// Turns whatever the user typed into an executable to run: an absolute or
/// relative path is checked directly, a bare name is looked up on `PATH`.
///
/// The program is named in every message, because with more than one backend
/// "executable was not found" on its own does not say which one to fix.
pub fn resolve_program_executable(
    program: QcProgram,
    configured: &str,
) -> Result<PathBuf, String> {
    let configured = configured.trim();
    if configured.is_empty() {
        return Err(format!(
            "Set the {} executable path first.",
            program.label()
        ));
    }
    let candidate = PathBuf::from(configured);
    if candidate.components().count() > 1 || candidate.is_absolute() {
        return candidate.is_file().then_some(candidate).ok_or_else(|| {
            format!(
                "{} executable was not found: {configured}",
                program.label()
            )
        });
    }
    let path = std::env::var_os("PATH")
        .ok_or_else(|| format!("Could not find {configured:?} because PATH is not available"))?;
    std::env::split_paths(&path)
        .map(|dir| dir.join(configured))
        .find(|path| path.is_file())
        .ok_or_else(|| format!("Could not find {configured:?} on PATH"))
}

/// Directory that holds xTB working directories. Deliberately not derived
/// from `CARGO_MANIFEST_DIR` (a compile-time, machine-specific path) --
/// `std::env::temp_dir()` works for any installed build, on any machine.
pub(super) fn xtb_scratch_dir() -> PathBuf {
    std::env::temp_dir().join("beavyr_xtb")
}

/// Deletes a finished run directory. Failures are ignored: leaving scratch
/// files behind is untidy but never a reason to fail the user's operation.
pub(super) fn discard_xtb_run_dir(dir: &Path) {
    let _ = fs::remove_dir_all(dir);
}

pub(super) fn create_xtb_run_dir(base_dir: &Path) -> Result<PathBuf, String> {
    fs::create_dir_all(base_dir).map_err(|e| format!("Failed to create xTB directory: {e}"))?;
    for _ in 0..128 {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|duration| duration.as_millis())
            .unwrap_or_default();
        let serial = XTB_RUN_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        let dir = base_dir.join(format!("xtb_{stamp}_{}_{}", std::process::id(), serial));
        match fs::create_dir(&dir) {
            Ok(()) => return Ok(dir),
            Err(err) if err.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(err) => return Err(format!("Failed to create xTB run directory: {err}")),
        }
    }
    Err("Could not allocate a unique xTB run directory".to_string())
}

/// A snapshot of where an in-progress `--opt` run currently stands, read
/// from xTB's own live stdout -- so a long run reads as "working, iteration 12"
/// rather than an unexplained silence the user has to just trust.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XtbProgress {
    pub iteration: u32,
    pub energy_hartree: Option<f64>,
    pub gradient_norm: Option<f64>,
}

/// Parses xTB's optimizer log format, e.g.:
/// ```text
/// .............................. CYCLE   12 ..............................
/// ...
///  * total energy  :    -9.0838683 Eh     change       -0.123E-05 Eh
///    gradient norm :     0.0088123 Eh/α   predicted     ...
/// ```
/// Returns the *latest* iteration reached and, if xTB has printed them yet
/// for it, its energy and gradient norm. `None` if no iteration marker has
/// appeared at all (xTB is still in its initial single-point setup).
fn parse_xtb_progress(stdout_text: &str) -> Option<XtbProgress> {
    let mut iteration = None;
    let mut energy_hartree = None;
    let mut gradient_norm = None;

    for line in stdout_text.lines() {
        let trimmed = line.trim();
        // The marker sits inside a row of dots (".... CYCLE   12 ...."),
        // never at the start of the line, so the digits have to be pulled
        // out of the middle rather than stripped as a prefix.
        if let Some(pos) = trimmed.find("CYCLE") {
            let digits: String = trimmed[pos + "CYCLE".len()..]
                .chars()
                .skip_while(|c| c.is_whitespace())
                .take_while(|c| c.is_ascii_digit())
                .collect();
            if let Ok(n) = digits.parse::<u32>() {
                // A new cycle starting invalidates the previous cycle's
                // energy/gradient readout until this one prints its own.
                iteration = Some(n);
                energy_hartree = None;
                gradient_norm = None;
            }
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix("* total energy") {
            energy_hartree = rest
                .trim_start_matches([':', ' '])
                .split_whitespace()
                .next()
                .and_then(|tok| tok.parse::<f64>().ok());
        } else if let Some(rest) = trimmed.strip_prefix("gradient norm") {
            gradient_norm = rest
                .trim_start_matches([':', ' '])
                .split_whitespace()
                .next()
                .and_then(|tok| tok.parse::<f64>().ok());
        }
    }

    iteration.map(|iteration| XtbProgress {
        iteration,
        energy_hartree,
        gradient_norm,
    })
}

/// Every iteration's energy from a full xTB `--opt` stdout capture, in
/// order -- unlike `parse_xtb_progress`, which only reports the latest one.
/// This is what feeds the energy-vs-iteration plot, and it is read the same
/// way for a converged run and one that stopped partway: whatever iterations
/// printed their energy line before the run ended are included.
fn parse_full_energy_history(stdout_text: &str) -> Vec<EnergyHistoryPoint> {
    let mut history = Vec::new();
    let mut current_iteration = None;

    for line in stdout_text.lines() {
        let trimmed = line.trim();
        if let Some(pos) = trimmed.find("CYCLE") {
            let digits: String = trimmed[pos + "CYCLE".len()..]
                .chars()
                .skip_while(|c| c.is_whitespace())
                .take_while(|c| c.is_ascii_digit())
                .collect();
            current_iteration = digits.parse::<u32>().ok();
            continue;
        }
        let Some(iteration) = current_iteration else {
            continue;
        };
        if let Some(rest) = trimmed.strip_prefix("* total energy") {
            if let Some(energy_hartree) = rest
                .trim_start_matches([':', ' '])
                .split_whitespace()
                .next()
                .and_then(|tok| tok.parse::<f64>().ok())
            {
                history.push(EnergyHistoryPoint {
                    iteration,
                    energy_hartree,
                });
                // Each iteration prints its energy line once; move on so a
                // later, unrelated line starting the same way (there is
                // none in practice, but this keeps the function robust to
                // xTB printing the total energy more than once per
                // iteration) cannot record a duplicate point.
                current_iteration = None;
            }
        }
    }

    history
}

/// Best-effort read of a run directory's xTB stdout, for the energy history.
/// Never errors: a missing or unreadable file just means no history to show.
/// The per-iteration energies to plot, from wherever the chosen program puts
/// them: xTB prints them to its log, Behemoth writes them into the
/// trajectory's comment lines.
///
/// Read separately from the trajectory so a run that failed or was cancelled
/// still plots how far it got.
fn energy_history_for(program: QcProgram, workdir: &Path) -> Vec<EnergyHistoryPoint> {
    match program {
        QcProgram::Xtb => read_energy_history(workdir),
        QcProgram::Behemoth => {
            let text = fs::read_to_string(workdir.join(super::behemoth::TRAJECTORY_FILE))
                .unwrap_or_default();
            super::behemoth::parse_trajectory_energies(&text)
                .into_iter()
                .enumerate()
                .map(|(i, energy)| EnergyHistoryPoint {
                    iteration: i as u32 + 1,
                    energy_hartree: energy,
                })
                .collect()
        }
    }
}

fn read_energy_history(workdir: &Path) -> Vec<EnergyHistoryPoint> {
    fs::read_to_string(workdir.join("xtb.stdout"))
        .map(|text| parse_full_energy_history(&text))
        .unwrap_or_default()
}

/// The last frame of a multi-frame xTB log, parsed for its
/// `energy: <value> ...` comment line.
fn parse_final_energy(trajectory_text: &str) -> Option<f64> {
    trajectory_text
        .lines()
        .filter(|line| line.trim_start().starts_with("energy:"))
        .next_back()
        .and_then(|line| line.split_whitespace().nth(1))
        .and_then(|tok| tok.parse::<f64>().ok())
}

/// Runs a geometry optimization with whichever program was chosen.
///
/// The process handling -- spawn, poll, cancel, kill -- is the same whatever
/// the program; only the command line and the files to read back differ, and
/// both of those are the only things that branch below.
#[allow(clippy::too_many_arguments)]
fn run_optimize_cancellable(
    program: QcProgram,
    binary: &Path,
    workdir: &Path,
    input_xyz: &str,
    charge: i32,
    uhf: i32,
    multiplicity: i32,
    cancel: &AtomicBool,
    child_slot: &Mutex<Option<Child>>,
) -> Result<XtbOptimizationOutput, String> {
    if binary.as_os_str().is_empty() {
        return Err(format!("{} path is empty", program.label()));
    }
    fs::create_dir_all(workdir).map_err(|e| format!("Failed to create workdir: {e}"))?;
    fs::write(workdir.join("input.xyz"), input_xyz)
        .map_err(|e| format!("Failed to write input.xyz: {e}"))?;

    let stdout = fs::File::create(workdir.join("xtb.stdout"))
        .map_err(|e| format!("Failed to create xtb.stdout: {e}"))?;
    let stderr = fs::File::create(workdir.join("xtb.stderr"))
        .map_err(|e| format!("Failed to create xtb.stderr: {e}"))?;
    let mut cmd = match program {
        QcProgram::Xtb => {
            let mut cmd = Command::new(binary);
            cmd.current_dir(workdir)
                .arg("input.xyz")
                .arg("--opt")
                .arg("--chrg")
                .arg(charge.to_string());
            if uhf > 0 {
                cmd.arg("--uhf").arg(uhf.to_string());
            }
            cmd
        }
        QcProgram::Behemoth => super::behemoth::optimize_command(
            binary,
            workdir,
            "input.xyz",
            charge,
            multiplicity,
        ),
    };
    cmd.stdout(Stdio::from(stdout)).stderr(Stdio::from(stderr));

    let child = cmd
        .spawn()
        .map_err(|e| format!("Failed to run {}: {e}", program.label()))?;
    *child_slot
        .lock()
        .map_err(|_| "xTB child-process lock was poisoned".to_string())? = Some(child);
    let status = loop {
        let mut child = child_slot
            .lock()
            .map_err(|_| "xTB child-process lock was poisoned".to_string())?;
        if cancel.load(Ordering::Relaxed) {
            if let Some(mut process) = child.take() {
                let _ = process.kill();
                let _ = process.wait();
            }
            return Err("xTB optimization cancelled".to_string());
        }
        let Some(process) = child.as_mut() else {
            return Err("xTB optimization cancelled".to_string());
        };
        if let Some(status) = process
            .try_wait()
            .map_err(|e| format!("Failed while waiting for xTB: {e}"))?
        {
            child.take();
            break status;
        }
        drop(child);
        std::thread::sleep(Duration::from_millis(50));
    };
    if !status.success() {
        return Err(format!(
            "{} optimization failed (see xtb.stderr in the run directory)",
            program.label()
        ));
    }

    let (trajectory_text, final_energy_hartree) = match program {
        QcProgram::Xtb => {
            let text = fs::read_to_string(workdir.join("xtbopt.log"))
                .or_else(|_| fs::read_to_string(workdir.join("xtbopt.xyz")))
                .map_err(|_| "No optimized geometry found".to_string())?;
            let energy = parse_final_energy(&text);
            (text, energy)
        }
        QcProgram::Behemoth => {
            let text = fs::read_to_string(workdir.join(super::behemoth::TRAJECTORY_FILE))
                .or_else(|_| {
                    fs::read_to_string(workdir.join(super::behemoth::OPTIMIZED_GEOMETRY_FILE))
                })
                .map_err(|_| "No optimized geometry found".to_string())?;
            // The summary is more precise than the last trajectory comment,
            // which is the energy *before* the final step.
            let log = fs::read_to_string(workdir.join("xtb.stdout")).unwrap_or_default();
            let energy = super::behemoth::parse_summary(&log).final_energy_hartree;
            (text, energy)
        }
    };
    let energy_history = energy_history_for(program, workdir);

    Ok(XtbOptimizationOutput {
        trajectory_text,
        final_energy_hartree,
        energy_history,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn write_xyz_string_round_trips_through_parse_xyz_angstrom() {
        let atoms = vec!["O".to_string(), "H".to_string()];
        let pos = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.96)];
        let text = write_xyz_string(&atoms, &pos);
        let (n, parsed_atoms, parsed_pos) = crate::molecule::parse_xyz_angstrom(&text);
        assert_eq!(n, 2);
        assert_eq!(parsed_atoms, atoms);
        assert!((parsed_pos[1][2] - 0.96).abs() < 1.0e-6);
    }

    /// A trimmed but faithful excerpt of real xTB `--opt` stdout (the exact
    /// line shapes confirmed by running xtb 6.7.1 against `examples/h2o2.xyz`
    /// during the design work for this feature): two cycles, the second only
    /// partway printed, matching a run caught mid-cycle.
    const SAMPLE_STDOUT_MID_RUN: &str = "\
........................................................................
.............................. CYCLE    1 ..............................
........................................................................
 * total energy  :    -9.0320735 Eh     change       -0.301E-13 Eh
   gradient norm :     0.2183805 Eh/\u{03b1}   predicted     0.0000000E+00 (-100.00%)
........................................................................
.............................. CYCLE    2 ..............................
........................................................................
 iter      E             dE          RMSdq      gap      omega  full diag
   1     -9.0837968 -0.908380E+01  0.899E-01    9.38       0.0  T
";

    #[test]
    fn parse_xtb_progress_reads_the_latest_cycle_only() {
        let progress = parse_xtb_progress(SAMPLE_STDOUT_MID_RUN).unwrap();
        assert_eq!(progress.iteration, 2, "cycle 2 has started, so it is now current");
        assert_eq!(
            progress.energy_hartree, None,
            "cycle 2 has not printed its own energy line yet"
        );
        assert_eq!(progress.gradient_norm, None);
    }

    #[test]
    fn parse_xtb_progress_reads_a_completed_cycles_energy_and_gradient() {
        let text = "\
.............................. CYCLE    1 ..............................
 * total energy  :    -9.0320735 Eh     change       -0.301E-13 Eh
   gradient norm :     0.2183805 Eh/\u{03b1}   predicted     0.0000000E+00 (-100.00%)
";
        let progress = parse_xtb_progress(text).unwrap();
        assert_eq!(progress.iteration, 1);
        assert!((progress.energy_hartree.unwrap() - -9.0320735).abs() < 1.0e-6);
        assert!((progress.gradient_norm.unwrap() - 0.2183805).abs() < 1.0e-6);
    }

    #[test]
    fn parse_xtb_progress_is_none_before_the_first_cycle() {
        assert_eq!(
            parse_xtb_progress("generating ANC from exact Hessian ...\n"),
            None
        );
    }

    #[test]
    fn format_elapsed_matches_common_readouts() {
        assert_eq!(format_elapsed(Duration::from_secs(7)), "00:07");
        assert_eq!(format_elapsed(Duration::from_secs(65)), "01:05");
        assert_eq!(format_elapsed(Duration::from_secs(3725)), "1:02:05");
    }

    #[test]
    fn parse_full_energy_history_collects_every_completed_iteration_in_order() {
        let history = parse_full_energy_history(SAMPLE_STDOUT_MID_RUN);
        // Iteration 1 completed and printed its energy; iteration 2 has only
        // just started (its own SCF table, no "* total energy" line yet), so
        // only one point should be collected.
        assert_eq!(history.len(), 1);
        assert_eq!(history[0].iteration, 1);
        assert!((history[0].energy_hartree - -9.0320735).abs() < 1.0e-6);
    }

    #[test]
    fn parse_full_energy_history_is_empty_for_text_with_no_iterations() {
        assert!(parse_full_energy_history("nothing relevant here\n").is_empty());
    }

    #[test]
    fn read_energy_history_returns_empty_for_a_missing_file_rather_than_erroring() {
        let dir = std::env::temp_dir().join(format!(
            "beavyr_xtb_missing_stdout_{}_{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        assert!(read_energy_history(&dir).is_empty());
    }

    /// A charge that makes the electron count impossible for xTB's own
    /// internal accounting is a reliable, fast way to make a real xTB
    /// process actually fail (as opposed to `--uhf` on a stable species,
    /// which xTB accepts silently -- see the design spec). This exercises
    /// the failure path end-to-end against the real binary: a failure must
    /// still surface whatever iteration history xTB printed before it
    /// stopped, since that partial run is exactly what the user asked to be
    /// able to inspect.
    #[test]
    fn a_failing_real_xtb_run_still_reports_whatever_history_it_produced() {
        let configured = default_xtb_executable();
        let Ok(xtb_path) = resolve_program_executable(QcProgram::Xtb, &configured) else {
            eprintln!("skipping: no xTB executable found (set BEAVYR_XTB_PATH to run this test)");
            return;
        };

        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("h2o2.xyz");
        let text = std::fs::read_to_string(&path).expect("examples/h2o2.xyz");
        let (_n, atoms, coords) = crate::molecule::parse_xyz_angstrom(&text);
        let pos: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();

        let cancel = AtomicBool::new(false);
        let child_slot = Mutex::new(None);
        let workdir = xtb_scratch_dir().join(format!("test_fail_{}", std::process::id()));
        let input_xyz = write_xyz_string(&atoms, &pos);

        // An enormous, impossible charge: xTB rejects this outright rather
        // than silently proceeding, unlike the `--uhf` case.
        let result = run_optimize_cancellable(
            QcProgram::Xtb,
            &xtb_path,
            &workdir,
            &input_xyz,
            99,
            0,
            1,
            &cancel,
            &child_slot,
        );
        assert!(result.is_err(), "an impossible charge must make xTB fail");

        // Mirror exactly what `start()`'s spawned closure does on failure:
        // read whatever iteration history xTB printed to stdout before it
        // stopped. Measured behaviour: xTB reaches and prints one geometry
        // iteration's energy even for this broken a charge before its
        // internal checks abort the run, so there is a real, non-empty
        // history to surface here -- exactly the case this feature exists
        // for. The assertion is `>=1` rather than an exact count, since that
        // exact number is xTB's own implementation detail, not this code's.
        let history = read_energy_history(&workdir);
        assert!(
            !history.is_empty(),
            "a failed run that reached at least one iteration must still \
             report that iteration's history"
        );

        let _ = fs::remove_dir_all(&workdir);
    }

    #[test]
    fn parse_final_energy_reads_the_last_frames_comment_line() {
        let text = "\
2
 energy: -1.0 gnorm: 0.5 xtb: 6.7.1
H 0.0 0.0 0.0
H 0.0 0.0 0.7
2
 energy: -1.234567 gnorm: 0.0001 xtb: 6.7.1
H 0.0 0.0 0.0
H 0.0 0.0 0.74
";
        assert_eq!(parse_final_energy(text), Some(-1.234567));
    }

    #[test]
    fn parse_final_energy_is_none_for_text_without_an_energy_comment() {
        assert_eq!(parse_final_energy("2\n\nH 0 0 0\nH 0 0 1\n"), None);
    }

    #[test]
    fn create_xtb_run_dir_allocates_a_fresh_unique_directory() {
        let base = std::env::temp_dir().join(format!(
            "beavyr_xtb_run_dir_test_{}_{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let a = create_xtb_run_dir(&base).unwrap();
        let b = create_xtb_run_dir(&base).unwrap();
        assert_ne!(a, b);
        assert!(a.is_dir());
        assert!(b.is_dir());
        let _ = fs::remove_dir_all(&base);
    }

    /// Reproduces the exact scenario the design spec was written against:
    /// runs the real `xtb` binary (if one can be resolved) on hydrogen
    /// peroxide and checks the trajectory it writes is what
    /// `trajectory::parse_multi_xyz` expects.
    #[test]
    fn a_real_xtb_run_produces_a_trajectory_the_player_can_load() {
        let configured = default_xtb_executable();
        let Ok(xtb_path) = resolve_program_executable(QcProgram::Xtb, &configured) else {
            eprintln!("skipping: no xTB executable found (set BEAVYR_XTB_PATH to run this test)");
            return;
        };

        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("h2o2.xyz");
        let text = std::fs::read_to_string(&path).expect("examples/h2o2.xyz");
        let (_n, atoms, coords) = crate::molecule::parse_xyz_angstrom(&text);
        let pos: Vec<Vec3> = coords
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();

        let cancel = AtomicBool::new(false);
        let child_slot = Mutex::new(None);
        let workdir = xtb_scratch_dir().join(format!("test_{}", std::process::id()));
        let input_xyz = write_xyz_string(&atoms, &pos);

        let result = run_optimize_cancellable(
            QcProgram::Xtb,
            &xtb_path,
            &workdir,
            &input_xyz,
            0,
            0,
            1,
            &cancel,
            &child_slot,
        );
        let output = result.expect("a valid closed-shell run should succeed");
        assert!(output.final_energy_hartree.is_some());

        let frames = trajectory::parse_multi_xyz(&output.trajectory_text)
            .expect("xtb's own trajectory format must parse");
        assert!(!frames.is_empty());
        assert_eq!(frames[0].atoms.len(), 4);

        let _ = fs::remove_dir_all(&workdir);
    }
}
