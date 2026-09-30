//! Drawing the general quantum-chemistry panel.
//!
//! Every rule this panel follows lives in `qchem_interfaces::job` and the
//! per-program catalogues, which have no GUI dependency and carry the tests.
//! This file is only the drawing.

use bevy::prelude::*;
use bevy_egui::egui;

use super::run::QcRunTask;
use super::QcPanelState;
use crate::molecule::Molecule;
use crate::qchem_interfaces::config;
use crate::qchem_interfaces::job::{self, JobType, QcMethod};
use crate::qchem_interfaces::method_ui;
use crate::qchem_interfaces::orca_run;
use crate::qchem_interfaces::program::QcProgram;
use crate::qchem_interfaces::valence::validate_electronic_state;
use crate::qchem_interfaces::xtb_optimize::{
    python_environment_row, resolve_program_executable, XtbPanelState,
};

/// Draws the panel and starts a run when the button is pressed.
#[allow(clippy::too_many_arguments)]
pub fn qc_panel(
    ctx: &egui::Context,
    panel: &mut QcPanelState,
    task: &mut QcRunTask,
    opt_panel: &mut XtbPanelState,
    opt_task: &crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask,
    mol: &Molecule,
) {
    let mut open = panel.open;
    egui::Window::new("Quantum Chemistry")
        .open(&mut open)
        .default_size([460.0, 640.0])
        .vscroll(true)
        .show(ctx, |ui| {
            body(ui, panel, task, opt_panel, mol);
        });
    panel.open = open;

    // The optimizer panel's own two windows, driven from its own flags. An
    // optimisation started here and one started there are the same
    // calculation, so they get the same summary and the same energy plot
    // rather than a second pair built to look almost alike.
    let report = opt_task
        .last_log
        .as_deref()
        .filter(|_| {
            opt_task.last_program == Some(QcProgram::Behemoth)
                || opt_task.last_program == Some(QcProgram::Dreiding)
        })
        .map(|log| {
            crate::qchem_interfaces::xtb_optimize::OptimizationReport::from_log(
                log,
                opt_task.last_duration,
            )
        });
    crate::qchem_interfaces::xtb_optimize::energy_history_window(
        ctx,
        &mut opt_panel.energy_plot_open,
        &opt_task.last_energy_history,
        opt_task.last_is_error,
    );
    if let Some(report) = &report {
        crate::qchem_interfaces::xtb_optimize::optimization_summary_window(
            ctx,
            &mut opt_panel.summary_open,
            report,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn body(
    ui: &mut egui::Ui,
    panel: &mut QcPanelState,
    task: &mut QcRunTask,
    opt_panel: &mut XtbPanelState,
    mol: &Molecule,
) {
    let running = task.is_running();

    ui.label(egui::RichText::new(panel.summary_line()).strong());
    ui.separator();

    program_row(ui, panel, opt_panel, running);
    charge_and_multiplicity(ui, panel, running);

    // For a Python backend, whether this environment can import it. Answered
    // before the run rather than minutes into one.
    let python_ready = python_environment_row(ui, panel.program, running);

    ui.separator();
    let level_ok = level_of_theory(ui, panel, running);

    ui.separator();
    job_row(ui, panel, running);
    if panel.needs_excited_state_options() {
        excited_state_row(ui, panel, running);
    }

    wavefunction_row(ui, panel, running);

    ui.separator();
    resources_row(ui, panel, running);

    // The charge and multiplicity have to describe a real electronic state.
    let electronic = validate_electronic_state(mol, panel.charge, panel.multiplicity);
    if let Err(error) = &electronic {
        ui.colored_label(
            egui::Color32::from_rgb(230, 120, 90),
            format!("\u{26a0} {error}"),
        );
    } else if panel.show_warnings {
        if let Ok((_, warnings)) = &electronic {
            for warning in warnings {
                ui.label(egui::RichText::new(format!("\u{26a0} {}", warning.message)).weak());
            }
        }
    }

    ui.separator();
    input_row(ui, panel, mol, running);

    ui.separator();
    run_row(ui, panel, task, opt_panel, mol, running, python_ready, level_ok, electronic.is_ok());
    results_row(ui, panel, task);
}

/// The input file: save it, or write it by hand.
///
/// The controls above cover what most runs need. Every one of these programs
/// can do things no set of dropdowns will reach -- a solvent model, an
/// embedding, a scan written in Python -- and the honest answer to that is a
/// box to type in, not another row of widgets.
fn input_row(ui: &mut egui::Ui, panel: &mut QcPanelState, mol: &Molecule, running: bool) {
    if !panel.has_input_file() {
        // xTB, Behemoth and Sparrow are driven by settings rather than by a
        // file, and DREIDING runs in this process, so there is nothing to edit
        // or save for any of them.
        return;
    }

    let positions: Vec<f64> = mol
        .pos
        .iter()
        .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
        .collect();

    ui.horizontal(|ui| {
        let was = panel.advanced;
        ui.add_enabled(
            !running,
            egui::Checkbox::new(&mut panel.advanced, "Edit the input directly"),
        );
        // Opening the box fills it from what the controls would have written,
        // so editing starts from something that runs rather than a blank page.
        if panel.advanced && !was && panel.advanced_text.trim().is_empty() {
            panel.advanced_text = panel.generated_input(&mol.atoms, &positions);
        }

        if ui
            .add_enabled(!running, egui::Button::new("Save input\u{2026}"))
            .on_hover_text("Write this input to a file without running it.")
            .clicked()
        {
            let text = if panel.advanced {
                panel.advanced_text.clone()
            } else {
                panel.generated_input(&mol.atoms, &positions)
            };
            if let Some(target) = rfd::FileDialog::new()
                .set_file_name(panel.input_file_name())
                .save_file()
            {
                let _ = std::fs::write(target, text);
            }
        }
    });

    if !panel.advanced {
        return;
    }

    ui.horizontal(|ui| {
        if ui
            .add_enabled(!running, egui::Button::new("Rebuild from the settings above"))
            .on_hover_text("Replaces what is in the box. Anything typed here is lost.")
            .clicked()
        {
            panel.advanced_text = panel.generated_input(&mol.atoms, &positions);
        }
        ui.label(
            egui::RichText::new(format!("written as {}", panel.input_file_name()))
                .small()
                .weak(),
        );
    });

    // What the controls above no longer decide, said plainly. A box that
    // silently overrode the dropdowns beside it would be the worst of both.
    ui.label(
        egui::RichText::new(
            "This text is sent exactly as it stands. The method, basis and job controls above \
             no longer apply, including the geometry, so a structure changed on screen will \
             not reach the calculation until you rebuild.",
        )
        .small()
        .color(egui::Color32::from_rgb(230, 170, 90)),
    );

    egui::ScrollArea::vertical()
        .max_height(260.0)
        .id_salt("qc_panel_advanced_input")
        .show(ui, |ui| {
            ui.add_enabled(
                !running,
                egui::TextEdit::multiline(&mut panel.advanced_text)
                    .font(egui::TextStyle::Monospace)
                    .desired_width(f32::INFINITY)
                    .desired_rows(14)
                    .code_editor(),
            );
        });
}

/// The one button that opens the calculation summary.
///
/// Everything else a finished run offers -- the energy plot, the optimisation
/// summary, the files to save -- lives inside that window, at the end, rather
/// than scattered across the form.
fn results_row(ui: &mut egui::Ui, panel: &mut QcPanelState, task: &QcRunTask) {
    if task.last.is_none() {
        return;
    }
    ui.horizontal(|ui| {
        let label = if panel.summary_open {
            "Hide Calculation Summary"
        } else {
            "Show Calculation Summary"
        };
        if ui.button(label).clicked() {
            panel.summary_open = !panel.summary_open;
        }
    });
}

/// The program selector. Changing it resets the level of theory, because the
/// catalogues do not translate between codes.
fn program_row(
    ui: &mut egui::Ui,
    panel: &mut QcPanelState,
    opt_panel: &mut XtbPanelState,
    running: bool,
) {
    ui.horizontal(|ui| {
        ui.label("Program");
        let mut chosen = panel.program;
        ui.add_enabled_ui(!running, |ui| {
            egui::ComboBox::from_id_salt("qc_panel_program")
                .selected_text(panel.program.label())
                .show_ui(ui, |ui| {
                    for program in QcProgram::ALL {
                        ui.selectable_value(&mut chosen, program, program.label());
                    }
                });
        });
        panel.set_program(chosen);
    });

    // The path, shared with the optimizer panel so a program pointed at once is
    // pointed at everywhere.
    if panel.program.needs_binary_path() {
        ui.horizontal(|ui| {
            let program = panel.program;
            ui.label(format!("{} path", program.label()));
            // The optimizer panel owns the path per program; borrow its field so
            // the two cannot disagree about where a program lives.
            let previous = opt_panel.program;
            opt_panel.program = program;
            let response = ui.add_enabled(
                !running,
                egui::TextEdit::singleline(opt_panel.path_mut())
                    .desired_width(210.0)
                    .hint_text(program.binary_hint()),
            );
            if response.lost_focus() {
                config::save_path(program, opt_panel.path());
            }
            if ui
                .add_enabled(!running, egui::Button::new("Browse\u{2026}"))
                .clicked()
            {
                let mut dialog = rfd::FileDialog::new();
                if let Some(dir) = std::path::Path::new(opt_panel.path())
                    .parent()
                    .filter(|p| p.is_dir())
                {
                    dialog = dialog.set_directory(dir);
                }
                if let Some(path) = dialog.pick_file() {
                    *opt_panel.path_mut() = path.display().to_string();
                    config::save_path(program, opt_panel.path());
                }
            }
            opt_panel.program = previous;
        });
    }
}

fn charge_and_multiplicity(ui: &mut egui::Ui, panel: &mut QcPanelState, running: bool) {
    ui.horizontal(|ui| {
        ui.label("Charge");
        ui.add_enabled(
            !running,
            egui::DragValue::new(&mut panel.charge).speed(1.0).range(-10..=10),
        );
        ui.label("Multiplicity");
        ui.add_enabled(
            !running,
            egui::DragValue::new(&mut panel.multiplicity)
                .speed(1.0)
                .range(1..=11),
        );
    });
}

/// The method, functional and basis rows. Returns whether the level of theory
/// can actually run.
fn level_of_theory(ui: &mut egui::Ui, panel: &mut QcPanelState, running: bool) -> bool {
    let methods = job::methods_for(panel.program);
    if !panel.has_level_of_theory() || methods.is_empty() {
        ui.label(
            egui::RichText::new(
                "A generic force field with no level of theory to set: its parameters follow \
                 from the atom types, which are perceived from the structure.",
            )
            .small()
            .weak(),
        );
        return true;
    }

    ui.horizontal(|ui| {
        ui.label("Method");
        ui.add_enabled_ui(!running, |ui| {
            let selected = panel.method.map_or("\u{2014}", |m| m.label());
            egui::ComboBox::from_id_salt("qc_panel_method")
                .selected_text(selected)
                .show_ui(ui, |ui| {
                    for &method in methods {
                        let mut chosen = panel.method;
                        if ui
                            .selectable_value(&mut chosen, Some(method), method.label())
                            .clicked()
                        {
                            panel.method = Some(method);
                            adopt_method_defaults(panel, method);
                        }
                    }
                });
        });
        if panel.method.is_some_and(|m| m.is_expensive()) {
            ui.label(
                egui::RichText::new("expensive")
                    .small()
                    .color(egui::Color32::from_rgb(230, 170, 90)),
            );
        }
    });

    if panel.needs_functional() {
        match panel.program {
            QcProgram::Orca => {
                method_ui::orca_functional_row(ui, &mut panel.config, running, "qc_panel")
            }
            QcProgram::PySCF => {
                method_ui::pyscf_functional_row(ui, &mut panel.config, running, "qc_panel")
            }
            _ => {
                ui.horizontal(|ui| {
                    ui.label("Functional");
                    ui.add_enabled(
                        !running,
                        egui::TextEdit::singleline(&mut panel.config.functional)
                            .desired_width(200.0),
                    );
                });
            }
        }
    }

    if panel.needs_basis() {
        // An explicitly correlated method needs one of its own purpose-built
        // basis sets; an ordinary one will not run. So the list is different
        // rather than the same list with a warning attached.
        if panel.method.is_some_and(|m| m.is_f12()) {
            f12_basis_row(ui, panel, running);
        } else {
            match panel.program {
                QcProgram::Orca => {
                    method_ui::orca_basis_row(ui, &mut panel.config, running, "qc_panel")
                }
                QcProgram::PySCF => {
                    method_ui::pyscf_basis_row(ui, &mut panel.config, running, "qc_panel")
                }
                _ => {
                    ui.horizontal(|ui| {
                        ui.label("Basis set");
                        ui.add_enabled(
                            !running,
                            egui::TextEdit::singleline(&mut panel.config.basis)
                                .desired_width(200.0),
                        );
                    });
                }
            }
        }

        // What else is going onto the keyword line, said out loud, because it
        // changes the input and nobody should have to work out where an extra
        // keyword came from.
        if panel.program == QcProgram::Orca {
            if let Some(method) = panel.method {
                let keywords =
                    orca_run::basis_keywords(method, &panel.config.functional, &panel.config.basis);
                if keywords.len() > 1 {
                    ui.label(
                        egui::RichText::new(format!(
                            "This method needs a fitting basis; {} added automatically.",
                            keywords[1..].join(" and ")
                        ))
                        .small()
                        .weak(),
                    );
                }
            }
        }
    }

    // A dispersion correction, for a density functional that does not already
    // carry one.
    if panel.program == QcProgram::Orca && panel.needs_functional() {
        method_ui::orca_dispersion_row(ui, &mut panel.config, running, "qc_panel");
    }

    // Anything the dropdowns do not cover.
    if panel.program == QcProgram::Orca {
        ui.horizontal(|ui| {
            ui.label("Extra keywords");
            ui.add_enabled(
                !running,
                egui::TextEdit::singleline(&mut panel.config.orca_extra)
                    .desired_width(210.0)
                    .hint_text("e.g. D4 TightSCF CPCM(water)"),
            );
        });
    }

    true
}

/// The purpose-built basis sets an F12 method needs.
fn f12_basis_row(ui: &mut egui::Ui, panel: &mut QcPanelState, running: bool) {
    ui.horizontal(|ui| {
        ui.label("Basis set");
        ui.add_enabled_ui(!running, |ui| {
            egui::ComboBox::from_id_salt("qc_panel_f12_basis")
                .selected_text(&panel.config.basis)
                .show_ui(ui, |ui| {
                    for (label, keyword) in orca_run::F12_BASIS_SETS {
                        let mut chosen = panel.config.basis.clone();
                        if ui
                            .selectable_value(&mut chosen, keyword.to_string(), label)
                            .clicked()
                        {
                            panel.config.basis = keyword.to_string();
                        }
                    }
                });
        });
        ui.label(egui::RichText::new("F12 sets only").small().weak());
    });
}

/// Moving to a method whose basis the old one cannot serve.
///
/// An explicitly correlated method needs an F12 basis, and an ordinary one will
/// not run; moving away from F12 leaves an F12 basis behind that an ordinary
/// method has no fitting set for. Either way the basis is switched to something
/// the new method can use rather than left to fail at run time.
fn adopt_method_defaults(panel: &mut QcPanelState, method: QcMethod) {
    let is_f12_basis = panel.config.basis.to_ascii_uppercase().contains("-F12");
    if method.is_f12() && !is_f12_basis {
        panel.config.basis = orca_run::DEFAULT_F12_BASIS.to_string();
    } else if !method.is_f12() && is_f12_basis {
        let (_, basis) = crate::qchem_interfaces::method::defaults_for(panel.program);
        panel.config.basis = basis;
    }
}

/// The job selector. Jobs the program cannot do are shown disabled with the
/// reason, not hidden: someone choosing between codes should be able to see
/// that ORCA follows a reaction coordinate and xTB does not.
fn job_row(ui: &mut egui::Ui, panel: &mut QcPanelState, running: bool) {
    ui.horizontal(|ui| {
        ui.label("Job");
        ui.add_enabled_ui(!running, |ui| {
            egui::ComboBox::from_id_salt("qc_panel_job")
                .selected_text(panel.job.label())
                .width(240.0)
                .show_ui(ui, |ui| {
                    for job in JobType::ALL {
                        match job::unsupported_reason(panel.program, job) {
                            None => {
                                ui.selectable_value(&mut panel.job, job, job.label());
                            }
                            Some(reason) => {
                                ui.add_enabled_ui(false, |ui| {
                                    let _ = ui.selectable_label(false, job.label());
                                })
                                .response
                                .on_disabled_hover_text(reason);
                            }
                        }
                    }
                });
        });
    });
    ui.label(egui::RichText::new(panel.job.description()).small().weak());
}

fn excited_state_row(ui: &mut egui::Ui, panel: &mut QcPanelState, running: bool) {
    ui.horizontal(|ui| {
        ui.label("States");
        ui.add_enabled(
            !running,
            egui::DragValue::new(&mut panel.excited.states)
                .speed(1.0)
                .range(1..=200),
        );
        ui.add_enabled(
            !running,
            egui::Checkbox::new(&mut panel.excited.tda, "Tamm-Dancoff"),
        );
    });
    ui.label(
        egui::RichText::new(if panel.excited.tda {
            "Tamm-Dancoff on: cheaper and more robust, and what both programs default to."
        } else {
            "Tamm-Dancoff off: the full linear response, costlier and more prone to \
             instabilities near a conical intersection."
        })
        .small()
        .weak(),
    );
}

/// The toggle that also writes the orbitals.
///
/// Every program can produce a Molden file, but each by its own route: ORCA
/// through a converter run afterwards, the others from inside the generated
/// script. What they have in common is that none does it unasked, which is why
/// this is a toggle rather than something always on -- for ORCA it is a second
/// process, and for the correlated methods the file is large.
fn wavefunction_row(ui: &mut egui::Ui, panel: &mut QcPanelState, running: bool) {
    ui.horizontal(|ui| {
        ui.add_enabled(
            !running,
            egui::Checkbox::new(&mut panel.wavefunction, "Write orbitals (Molden)"),
        );
    });
    if panel.wavefunction {
        ui.label(
            egui::RichText::new(
                "Loaded into Surface Tools when the run finishes, and offered to save.",
            )
            .small()
            .weak(),
        );
    }
}

fn resources_row(ui: &mut egui::Ui, panel: &mut QcPanelState, running: bool) {
    ui.horizontal(|ui| {
        ui.label("Memory (MB)");
        ui.add_enabled(
            !running,
            egui::DragValue::new(&mut panel.config.memory_mb)
                .speed(64.0)
                .range(256..=1_048_576),
        );
        ui.label("nproc");
        ui.add_enabled(
            !running,
            egui::DragValue::new(&mut panel.config.nproc).speed(1.0).range(1..=1024),
        );
    });
    if panel.program == QcProgram::Orca && panel.config.nproc > 1 {
        let per_core = (panel.config.memory_mb / panel.config.nproc).max(256);
        ui.label(
            egui::RichText::new(format!(
                "ORCA takes memory per core, so this is written as {per_core} MB each."
            ))
            .small()
            .weak(),
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn run_row(
    ui: &mut egui::Ui,
    panel: &mut QcPanelState,
    task: &mut QcRunTask,
    opt_panel: &XtbPanelState,
    mol: &Molecule,
    running: bool,
    python_ready: bool,
    level_ok: bool,
    electronic_ok: bool,
) {
    ui.horizontal(|ui| {
        let can_run =
            !running && !mol.atoms.is_empty() && python_ready && level_ok && electronic_ok;
        let button = ui.add_enabled(can_run, egui::Button::new("Run"));
        let button = if mol.atoms.is_empty() {
            button.on_disabled_hover_text("There is no structure on screen to calculate.")
        } else if !python_ready {
            button.on_disabled_hover_text(
                "This backend needs a Python that can import it -- see the note above.",
            )
        } else if !electronic_ok {
            button.on_disabled_hover_text(
                "This charge and multiplicity do not describe a possible electronic state.",
            )
        } else {
            button
        };

        if button.clicked() {
            panel.show_warnings = true;
            let program = panel.program;
            // The path the optimizer panel remembers for this program.
            let entered = opt_panel
                .paths
                .iter()
                .find(|(p, _)| *p == program)
                .map(|(_, path)| path.clone())
                .unwrap_or_default();
            match resolve_program_executable(program, &entered) {
                Ok(binary) => {
                    let positions: Vec<f64> = mol
                        .pos
                        .iter()
                        .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
                        .collect();
                    task.start(
                        program,
                        panel.method,
                        panel.job,
                        binary,
                        mol.atoms.clone(),
                        positions,
                        panel.charge,
                        panel.multiplicity,
                        panel.config.clone(),
                        panel.excited,
                        panel.wavefunction,
                        panel.advanced.then(|| panel.advanced_text.clone()),
                    );
                }
                Err(err) => {
                    task.last_message = Some(err);
                    task.last_is_error = true;
                }
            }
        }

        if running {
            if ui.button("Cancel").clicked() {
                task.cancel();
            }
            if let Some(elapsed) = task.elapsed() {
                ui.label(crate::qchem_interfaces::xtb_optimize::format_elapsed(elapsed));
            }
            ui.spinner();
        }

    });

    if let Some(message) = &task.last_message {
        if task.last_is_error {
            ui.colored_label(egui::Color32::from_rgb(230, 120, 90), message);
        } else {
            ui.label(egui::RichText::new(message).weak());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qchem_interfaces::job::QcMethod;

    fn panel() -> QcPanelState {
        QcPanelState::default()
    }

    /// Moving to an explicitly correlated method switches the basis to one it
    /// can actually use. An ordinary basis simply will not run F12, so leaving
    /// it would guarantee a failed job.
    #[test]
    fn choosing_f12_switches_to_an_f12_basis() {
        let mut state = panel();
        state.config.basis = "def2-SVP".to_string();
        adopt_method_defaults(&mut state, QcMethod::Mp2F12);
        assert_eq!(state.config.basis, orca_run::DEFAULT_F12_BASIS);
    }

    /// And moving away from one puts an ordinary basis back, rather than
    /// leaving an F12 set behind that an ordinary method has no fitting basis
    /// for.
    #[test]
    fn leaving_f12_restores_an_ordinary_basis() {
        let mut state = panel();
        state.config.basis = "cc-pVDZ-F12".to_string();
        adopt_method_defaults(&mut state, QcMethod::Dft);
        assert!(!state.config.basis.contains("F12"), "{}", state.config.basis);
        assert_eq!(state.config.basis, "def2-SVP");
    }

    /// An F12 basis already chosen is left alone when another F12 method is
    /// picked, rather than being reset to the default each time.
    #[test]
    fn an_existing_f12_basis_survives_a_change_between_f12_methods() {
        let mut state = panel();
        state.config.basis = "cc-pVTZ-F12".to_string();
        adopt_method_defaults(&mut state, QcMethod::CcsdTF12);
        assert_eq!(state.config.basis, "cc-pVTZ-F12");
    }

    /// The editor shows the input that would actually be sent. An editor
    /// showing something else would be worse than no editor, because what you
    /// edited and what ran would silently differ.
    #[test]
    fn the_editor_starts_from_what_would_be_sent() {
        let mut state = panel();
        state.program = QcProgram::Orca;
        state.method = Some(QcMethod::Dft);
        state.config.functional = "M062X".to_string();
        state.config.basis = "def2-TZVP".to_string();
        state.job = JobType::Optimize;

        let atoms: Vec<String> = ["O", "H", "H"].iter().map(|s| s.to_string()).collect();
        let water = [0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692];
        let text = state.generated_input(&atoms, &water);

        // Everything the controls said, in ORCA's own spelling.
        assert!(text.contains("! M062X def2-TZVP Opt"), "{text}");
        assert!(text.contains("* xyz 0 1"), "{text}");
        assert!(text.contains("0.1173"), "the geometry is in it: {text}");
    }

    /// Only the programs driven by an input file offer one to edit or save.
    /// xTB and Behemoth take everything on their command lines and DREIDING
    /// runs in this process, so for those there is nothing to show.
    #[test]
    fn only_file_driven_programs_offer_an_input() {
        let mut state = panel();
        for program in [
            QcProgram::Orca,
            QcProgram::Psi4,
            QcProgram::Psi4Py,
            QcProgram::PySCF,
        ] {
            state.set_program(program);
            assert!(state.has_input_file(), "{program:?}");
            assert!(!state.input_file_name().is_empty(), "{program:?}");
        }
        // Driven by settings, not by a file. Sparrow belongs here despite
        // Beavyr generating a script for it: that script is this program's
        // scaffolding, not an input format Sparrow has.
        for program in [
            QcProgram::Xtb,
            QcProgram::Behemoth,
            QcProgram::Dreiding,
            QcProgram::SparrowPy,
        ] {
            state.set_program(program);
            assert!(!state.has_input_file(), "{program:?}");
        }
    }

    /// Each program's input is offered under the name it is actually written
    /// as, so a saved file can be run again by hand without renaming.
    #[test]
    fn the_suggested_name_is_the_one_it_is_written_as() {
        let mut state = panel();
        for (program, name) in [
            (QcProgram::Orca, "orca_job.inp"),
            (QcProgram::Psi4, "input.dat"),
            (QcProgram::Psi4Py, "run_psi4.py"),
            (QcProgram::PySCF, "run_pyscf.py"),
        ] {
            state.set_program(program);
            assert_eq!(state.input_file_name(), name, "{program:?}");
        }
    }

    /// A Python-driven program's input is a Python script, and ORCA's is not.
    /// Getting these the wrong way round would hand a program a file it cannot
    /// read.
    #[test]
    fn the_generated_input_is_in_the_right_language() {
        let atoms: Vec<String> = ["O", "H", "H"].iter().map(|s| s.to_string()).collect();
        let water = [0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692];

        let mut state = panel();
        state.set_program(QcProgram::PySCF);
        state.method = Some(QcMethod::Dft);
        let pyscf = state.generated_input(&atoms, &water);
        assert!(pyscf.contains("gto.M("), "{pyscf}");
        assert!(!pyscf.starts_with('!'), "{pyscf}");

        state.set_program(QcProgram::Orca);
        state.method = Some(QcMethod::Dft);
        let orca = state.generated_input(&atoms, &water);
        assert!(orca.contains("%pal nprocs"), "{orca}");
        assert!(!orca.contains("gto.M("), "{orca}");
    }

    /// Switching program resets the level of theory but keeps what describes
    /// the molecule and the machine.
    #[test]
    fn switching_program_keeps_the_molecule_and_the_machine() {
        let mut state = panel();
        state.charge = -1;
        state.multiplicity = 3;
        state.config.nproc = 8;
        state.config.memory_mb = 8192;
        state.config.functional = "M062X".to_string();

        state.set_program(QcProgram::PySCF);

        assert_eq!(state.charge, -1, "the charge describes the molecule");
        assert_eq!(state.multiplicity, 3);
        assert_eq!(state.config.nproc, 8, "the thread count describes the machine");
        assert_eq!(state.config.memory_mb, 8192);
        assert_ne!(
            state.config.functional, "M062X",
            "ORCA's spelling means nothing to PySCF"
        );
    }

    /// A job the new program cannot do falls back to one every backend can.
    #[test]
    fn switching_to_a_program_without_the_job_falls_back() {
        let mut state = panel();
        state.job = JobType::Irc;
        assert!(job::supports(QcProgram::Orca, JobType::Irc));

        state.set_program(QcProgram::Xtb);
        assert_eq!(state.job, JobType::Optimize, "xTB has no IRC");
        assert!(job::supports(QcProgram::Xtb, state.job));
    }

    /// A job the new program *can* do is kept, rather than reset for no reason.
    #[test]
    fn a_job_the_new_program_supports_is_kept() {
        let mut state = panel();
        state.job = JobType::Frequencies;
        state.set_program(QcProgram::PySCF);
        assert_eq!(state.job, JobType::Frequencies);
    }

    /// A composite functional takes the basis row away: it brings its own.
    #[test]
    fn a_composite_hides_the_basis_row() {
        let mut state = panel();
        state.method = Some(QcMethod::Dft);
        state.config.functional = "r2SCAN-3c".to_string();
        assert!(!state.needs_basis());
        state.config.functional = "B3LYP".to_string();
        assert!(state.needs_basis());
    }

    /// A force field has no level of theory at all.
    #[test]
    fn dreiding_has_no_level_of_theory() {
        let mut state = panel();
        state.set_program(QcProgram::Dreiding);
        assert!(!state.has_level_of_theory());
        assert!(!state.needs_functional());
        assert!(!state.needs_basis());
    }

    /// The excited-state rows appear for exactly one job.
    #[test]
    fn the_excited_state_rows_belong_to_tddft() {
        let mut state = panel();
        for job in JobType::ALL {
            state.job = job;
            assert_eq!(
                state.needs_excited_state_options(),
                job == JobType::TdDft,
                "{job:?}"
            );
        }
    }
}
