//! What a finished run has to say for itself.
//!
//! Every job gets this window. Which iterations it lists depends on what the
//! job actually did: a geometry job has cycles, a job that never moved the
//! nuclei has only its self-consistent field history. Both end with the
//! energies.
//!
//! The geometry cycles are read back out of the trajectory rather than parsed
//! from the program's log a second time. The trajectory already carries one
//! energy per cycle, in the format every backend writes, so there is nothing
//! here that a second parser would add.

use bevy_egui::egui;

use super::run::QcRunOutput;
use crate::qchem_interfaces::behemoth::parse_trajectory_energies;
use crate::qchem_interfaces::job::SummaryKind;

/// Hartree to kJ/mol, for showing a step's energy change in units a chemist
/// reads without converting.
const HARTREE_TO_KJ_MOL: f64 = 2625.499_639_479_1;

/// Draws the run summary window.
///
/// `opt_panel` and `opt_task` are the optimizer panel's own state. After a
/// geometry optimisation this window ends with buttons into that panel's
/// energy plot and optimisation summary, so everything a finished optimisation
/// has to say is reachable from the one window that reported it.
#[allow(clippy::too_many_arguments)]
pub fn summary_window(
    ctx: &egui::Context,
    open: &mut bool,
    confirm_close: &mut bool,
    output: &QcRunOutput,
    duration: Option<std::time::Duration>,
    opt_panel: &mut crate::qchem_interfaces::xtb_optimize::XtbPanelState,
    opt_task: &crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask,
) {
    // The close button is intercepted rather than obeyed. egui closes a window
    // the moment its X is clicked, and this one sits beside the Save buttons
    // for files that go with the scratch directory -- so a mis-click after a
    // long calculation loses the only view of what came back.
    let mut still_open = *open;
    egui::Window::new("Calculation summary")
        .open(&mut still_open)
        .default_size([470.0, 520.0])
        .vscroll(true)
        .show(ctx, |ui| {
            header(ui, output, duration);
            ui.separator();

            match output.job.summary_kind() {
                SummaryKind::GeometryCycles => {
                    geometry_cycles(ui, output);
                    // A geometry job ran a self-consistent field at every cycle;
                    // the last one's history is still worth having, folded away.
                    if !output.scf.is_empty() {
                        ui.separator();
                        egui::CollapsingHeader::new(format!(
                            "Self-consistent field, last cycle ({} iterations)",
                            output.scf.len()
                        ))
                        .default_open(false)
                        .show(ui, |ui| scf_table(ui, output));
                    }
                }
                SummaryKind::ScfOnly => scf_table(ui, output),
            }

            ui.separator();
            energetics(ui, output);

            if !output.saveable.is_empty() {
                ui.separator();
                saveable_files(ui, output);
            }

            // A geometry job has more to show, and it is the optimizer panel
            // that shows it. These open that panel's own windows rather than
            // second copies, so an optimisation run from here looks exactly
            // like one run from there.
            if output.job.produces_trajectory() {
                ui.separator();
                optimization_buttons(ui, opt_panel, opt_task);
            }
        });

    // The X was clicked. Ask rather than obey.
    if !still_open && *open {
        *confirm_close = true;
    }

    if *confirm_close {
        confirm_close_window(ctx, open, confirm_close, output);
    }
}

/// "Do you really want to close it?"
///
/// Says what closing actually costs, which is the point of asking: the numbers
/// go, and any file not yet saved goes with the scratch directory. A prompt
/// that only asked would be an obstacle; one that says what is at stake is a
/// warning.
fn confirm_close_window(
    ctx: &egui::Context,
    open: &mut bool,
    confirm_close: &mut bool,
    output: &QcRunOutput,
) {
    let unsaved = output
        .saveable
        .iter()
        .filter(|(_, path)| path.is_file())
        .count();

    egui::Window::new("Close the calculation summary?")
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
        .show(ctx, |ui| {
            ui.label("This is the only record of the run in Beavyr. Closing it discards it.");
            if unsaved > 0 {
                ui.label(
                    egui::RichText::new(format!(
                        "{unsaved} file{} from this run {} still only in the scratch directory.                          Save anything you want to keep first.",
                        if unsaved == 1 { "" } else { "s" },
                        if unsaved == 1 { "is" } else { "are" },
                    ))
                    .color(egui::Color32::from_rgb(230, 170, 90)),
                );
            }
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                // The safe choice first and on the left, so the reflexive
                // click is the one that keeps the results.
                if ui.button("Keep it open").clicked() {
                    *confirm_close = false;
                }
                if ui.button("Close").clicked() {
                    *confirm_close = false;
                    *open = false;
                }
            });
        });
}

/// The way into the optimizer panel's energy plot and optimisation summary.
fn optimization_buttons(
    ui: &mut egui::Ui,
    opt_panel: &mut crate::qchem_interfaces::xtb_optimize::XtbPanelState,
    opt_task: &crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask,
) {
    use crate::qchem_interfaces::xtb_optimize::OptimizationReport;

    ui.strong("Optimization");
    ui.horizontal(|ui| {
        let has_history = !opt_task.last_energy_history.is_empty();
        if ui
            .add_enabled(
                has_history,
                egui::Button::new(if opt_panel.energy_plot_open {
                    "Hide Energy Plot"
                } else {
                    "Show Energy Plot"
                }),
            )
            .on_disabled_hover_text("This run recorded no per-cycle energies.")
            .clicked()
        {
            opt_panel.energy_plot_open = !opt_panel.energy_plot_open;
        }

        // The convergence table only exists where the program printed one this
        // reader understands. Shown disabled with the reason rather than
        // hidden, so a missing button is never a mystery.
        let has_report = opt_task
            .last_log
            .as_deref()
            .map(|log| OptimizationReport::from_log(log, opt_task.last_duration))
            .is_some_and(|report| report.has_content());
        if ui
            .add_enabled(
                has_report,
                egui::Button::new(if opt_panel.summary_open {
                    "Hide Optimization Summary"
                } else {
                    "Show Optimization Summary"
                }),
            )
            .on_disabled_hover_text(
                "This program printed no per-cycle convergence table Beavyr can read. The \
                 energy plot and the cycle table above come from the trajectory instead.",
            )
            .clicked()
        {
            opt_panel.summary_open = !opt_panel.summary_open;
        }
    });
}

/// The files this run produced, each with a Save button.
///
/// The scratch directory is discarded once the panel is done with it, so these
/// are offered while they still exist. A gradient or a Hessian is often the
/// whole point of running the job, and losing it to a cleanup would be a poor
/// reward for waiting.
fn saveable_files(ui: &mut egui::Ui, output: &QcRunOutput) {
    ui.strong("Files");
    for (what, path) in &output.saveable {
        ui.horizontal(|ui| {
            ui.label(*what);
            let exists = path.is_file();
            if ui
                .add_enabled(exists, egui::Button::new("Save\u{2026}"))
                .on_disabled_hover_text("This file is no longer in the run directory.")
                .clicked()
            {
                let suggested = path
                    .file_name()
                    .map(|name| name.to_string_lossy().to_string())
                    .unwrap_or_else(|| "orca_output".to_string());
                if let Some(target) = rfd::FileDialog::new()
                    .set_file_name(suggested)
                    .save_file()
                {
                    // A copy, not a move: the panel may still need the original,
                    // and a second Save should work as well as the first.
                    if let Err(err) = std::fs::copy(path, &target) {
                        ui.colored_label(
                            egui::Color32::from_rgb(230, 120, 90),
                            format!("Could not save: {err}"),
                        );
                    }
                }
            }
        });
    }
}

/// What ran, and whether it got there.
fn header(ui: &mut egui::Ui, output: &QcRunOutput, duration: Option<std::time::Duration>) {
    ui.horizontal(|ui| {
        ui.strong(output.job.label());
        ui.weak(format!("\u{2014} {}", output.program.label()));
        if let Some(duration) = duration {
            ui.weak(crate::qchem_interfaces::xtb_optimize::format_elapsed(duration));
        }
    });

    match output.converged {
        Some(true) => {
            ui.colored_label(egui::Color32::from_rgb(120, 190, 120), "Converged.");
        }
        Some(false) => {
            // Not an error. An unconverged geometry is still where the
            // optimiser got to, and is worth looking at -- but it is not an
            // answer, and saying so plainly is the point.
            ui.colored_label(
                egui::Color32::from_rgb(230, 170, 90),
                "\u{26a0} Did not converge. The structure shown is where it stopped, not a \
                 stationary point.",
            );
        }
        None => {}
    }
}

/// One row per geometry cycle: the energy, and the change from the cycle before.
fn geometry_cycles(ui: &mut egui::Ui, output: &QcRunOutput) {
    let Some(text) = &output.trajectory_text else {
        ui.weak("No optimization path was recorded.");
        return;
    };
    let energies = parse_trajectory_energies(text);
    if energies.is_empty() {
        ui.weak("No optimization path was recorded.");
        return;
    }

    ui.strong(format!("Geometry cycles ({})", energies.len()));
    egui::Grid::new("qc_geometry_cycles")
        .striped(true)
        .num_columns(3)
        .show(ui, |ui| {
            ui.strong("cycle");
            ui.strong("E [Eh]");
            ui.strong("dE [kJ/mol]");
            ui.end_row();

            for (index, energy) in energies.iter().enumerate() {
                ui.label(format!("{}", index + 1));
                ui.label(format!("{energy:.10}"));
                if index == 0 {
                    ui.label("\u{2014}");
                } else {
                    let change = (energy - energies[index - 1]) * HARTREE_TO_KJ_MOL;
                    ui.label(format!("{change:+.4}"));
                }
                ui.end_row();
            }
        });

    // The total drop, which is the number worth knowing at a glance.
    if energies.len() > 1 {
        let drop = (energies[energies.len() - 1] - energies[0]) * HARTREE_TO_KJ_MOL;
        ui.weak(format!("Total change: {drop:+.4} kJ/mol"));
    }
}

/// One row per self-consistent field iteration.
fn scf_table(ui: &mut egui::Ui, output: &QcRunOutput) {
    if output.scf.is_empty() {
        ui.weak("No self-consistent field history was printed.");
        return;
    }

    ui.strong(format!(
        "Self-consistent field ({} iterations)",
        output.scf.len()
    ));
    egui::Grid::new("qc_scf_iterations")
        .striped(true)
        .num_columns(3)
        .show(ui, |ui| {
            ui.strong("iter");
            ui.strong("E [Eh]");
            ui.strong("dE [Eh]");
            ui.end_row();

            for step in &output.scf {
                ui.label(format!("{}", step.iteration));
                ui.label(format!("{:.10}", step.energy_hartree));
                ui.label(format!("{:+.3e}", step.delta_e));
                ui.end_row();
            }
        });
}

/// The energies at the end, with the excited-state trap called out.
fn energetics(ui: &mut egui::Ui, output: &QcRunOutput) {
    let e = &output.energetics;
    if !e.has_content() {
        ui.weak("No energies were found in the output.");
        return;
    }

    ui.strong("Energetics");
    egui::Grid::new("qc_energetics")
        .striped(true)
        .num_columns(2)
        .show(ui, |ui| {
            if let Some(value) = e.scf_energy_hartree {
                ui.label("Self-consistent field");
                ui.label(format!("{value:.10} Eh"));
                ui.end_row();
            }
            if let Some(value) = e.correlation_energy_hartree {
                ui.label("Correlation");
                ui.label(format!("{value:.10} Eh"));
                ui.end_row();
            }
            if let Some(value) = e.dispersion_correction_hartree {
                ui.label("Dispersion correction");
                ui.label(format!("{value:.10} Eh"));
                ui.end_row();
            }
            if let Some(value) = e.ground_state_hartree() {
                ui.strong("Ground state");
                ui.strong(format!("{value:.10} Eh"));
                ui.end_row();
            }
        });

    // The one place a plausible-looking number would be the wrong one.
    if e.final_includes_excitation {
        if let Some(value) = e.final_energy_hartree {
            ui.label(
                egui::RichText::new(format!(
                    "ORCA's own \u{201C}final single point energy\u{201D} for this run is \
                     {value:.10} Eh. That is the first excited state, not the ground state: it \
                     is the self-consistent field energy plus the first root's excitation. The \
                     ground-state figure above is the one to quote."
                ))
                .small()
                .weak(),
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::qchem_interfaces::behemoth::parse_trajectory_energies;
    use crate::qchem_interfaces::job::{JobType, SummaryKind};

    /// The cycle table is read from the trajectory, so a trajectory in the
    /// shared format is all the summary needs -- no second parser, and nothing
    /// that could disagree with the energy plot.
    #[test]
    fn the_cycle_table_comes_from_the_trajectory() {
        let text = concat!(
            "3\ncycle 1 E = -76.1000000000 Eh\nO 0 0 0\nH 0 0 0\nH 0 0 0\n",
            "3\ncycle 2 E = -76.2000000000 Eh\nO 0 0 0\nH 0 0 0\nH 0 0 0\n",
            "3\ncycle 3 E = -76.2100000000 Eh\nO 0 0 0\nH 0 0 0\nH 0 0 0\n"
        );
        let energies = parse_trajectory_energies(text);
        assert_eq!(energies.len(), 3);
        // Downhill, as an optimisation should be.
        assert!(energies[2] < energies[0]);

        let drop = (energies[2] - energies[0]) * super::HARTREE_TO_KJ_MOL;
        assert!((drop + 288.8).abs() < 1.0, "about -289 kJ/mol, got {drop}");
    }

    /// A job that never moved the nuclei shows its self-consistent field
    /// history instead of an empty cycle table.
    #[test]
    fn a_static_job_shows_the_scf_history() {
        assert_eq!(JobType::SinglePoint.summary_kind(), SummaryKind::ScfOnly);
        assert_eq!(JobType::Optimize.summary_kind(), SummaryKind::GeometryCycles);
    }
}
