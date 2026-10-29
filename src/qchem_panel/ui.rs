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
use crate::qchem_interfaces::job::{self, JobType};
use crate::qchem_interfaces::job::QcMethod;
use crate::qchem_interfaces::orca_run;
use crate::qchem_interfaces::program::QcProgram;
use crate::qchem_interfaces::valence::validate_electronic_state;
use crate::qchem_interfaces::xtb_optimize::{
    python_environment_row, resolve_program_executable, XtbPanelState,
};

#[path = "modern_ui.rs"]
mod modern_ui;
#[path = "dreiding_warning.rs"]
mod dreiding_warning;
pub(crate) use dreiding_warning::dreiding_warning_window;
pub(crate) use modern_ui::{docked_width, input_editor_window, is_docked, project_presentation, restore_project_presentation, project_input_state, restore_project_input_state};

/// Draws the panel and starts a run when the button is pressed.
#[allow(clippy::too_many_arguments)]
pub fn qc_panel(
    ctx: &egui::Context,
    host_ui: &mut egui::Ui,
    panel: &mut QcPanelState,
    task: &mut QcRunTask,
    opt_panel: &mut XtbPanelState,
    mol: &mut Molecule,
) -> Option<egui::Rect> {
    if !task.is_running() {
        panel.config.nproc = panel.config.nproc.clamp(1, 4);
    }
    modern_ui::show(ctx, host_ui, panel, task, opt_panel, mol)
}

fn program_path_row(
    ui: &mut egui::Ui,
    panel: &QcPanelState,
    opt_panel: &mut XtbPanelState,
    running: bool,
) {
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
                    .desired_width((ui.available_width() - 96.0).clamp(100.0, 210.0))
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

/// Both presentations launch through the same calculation path.
fn start_run(
    ctx: &egui::Context,
    panel: &mut QcPanelState,
    task: &mut QcRunTask,
    opt_panel: &XtbPanelState,
    mol: &mut Molecule,
) {
    if task.is_running() {
        return;
    }
    panel.config.nproc = panel.config.nproc.clamp(1, 4);
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
            if let Err(error) = crate::calculation_geometry::prepare_display(ctx, mol) {
                task.last_message = Some(error);
                task.last_is_error = true;
                return;
            }
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

/// The energy plot and the optimisation summary of the last run.
///
/// Drawn every frame whether or not the Quantum Chemistry window is open. The
/// run's results outlive the panel that started it: closing the panel must
/// not take the energy plot with it, and nothing else draws these windows now
/// that the Geometry Optimization tab is gone.
pub fn result_windows(
    ctx: &egui::Context,
    opt_panel: &mut XtbPanelState,
    opt_task: &crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask,
) {
    // The optimizer panel's own two windows, driven from its own flags. An
    // optimisation started here and one started there are the same
    // calculation, so they get the same summary and the same energy plot
    // rather than a second pair built to look almost alike.
    // For every program. The Geometry Optimization tab used to draw this for
    // xTB runs; that tab is gone, so this is now the only place, and the
    // log reader decides whether a program printed a summary at all.
    let report = opt_task
        .last_log
        .as_deref()
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
