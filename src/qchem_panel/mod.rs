//! The general quantum-chemistry panel: any job, any installed code.
//!
//! The other panels are each shaped around one question -- optimise this, give
//! me its frequencies, plot me a spectrum -- and are quicker for it. They stay.
//! This one is the general case: pick a program, a level of theory and a job,
//! and run it.
//!
//! What it deliberately does not do is display results. A finished run hands
//! its trajectory to the trajectory player, its Hessian to the Vibrations
//! panel and its excited states to the UV-Vis panel, because those already know
//! how to show those things and a second copy would drift from the first. What
//! stays here is the run summary: the iterations, and the energies at the end.

mod calculation_record;
pub mod calculation_export;
pub mod live_output;
mod program_selection;
mod input_generation;
pub mod scan_controls;
mod mecp_controls;
mod geometry_search_controls;
pub mod scan_results;
pub mod run;
pub mod summary;
pub mod ui;

use bevy::prelude::*;

use crate::qchem_interfaces::job::{ExcitedStateOptions, JobType, QcMethod};
use crate::qchem_interfaces::method::MethodConfig;
use crate::qchem_interfaces::program::QcProgram;

/// Everything the panel remembers between frames.
#[derive(Resource)]
#[derive(serde::Serialize, serde::Deserialize)]
#[derive(Clone)]
pub struct QcPanelState {
    pub open: bool,
    pub program: QcProgram,
    pub charge: i32,
    pub multiplicity: i32,
    /// `None` only for DREIDING, which has no level of theory to choose.
    pub method: Option<QcMethod>,
    /// The functional, basis, resources and ORCA keyword line, shared in shape
    /// with the other panels so one method block serves all of them.
    pub config: MethodConfig,
    pub job: JobType,
    pub excited: ExcitedStateOptions,
    /// Whether to also write the molecular orbitals in Molden format.
    ///
    /// Off by default: every program has to be asked for it, and for some it
    /// costs a second step. When it is on the orbitals are loaded into the
    /// Surface tool as soon as the run finishes, which is usually the reason
    /// for asking.
    pub wavefunction: bool,
    /// Whether the input is being written by hand rather than built from the
    /// controls above.
    ///
    /// The dropdowns cover what most runs need and nothing more. Every one of
    /// these programs can do things no set of controls will ever reach -- a
    /// solvent model, an embedding, a scan written in Python -- and the honest
    /// answer to that is a text box, not another row of widgets.
    pub advanced: bool,
    /// The input as it will be written, when [`QcPanelState::advanced`] is on.
    ///
    /// Filled from what the controls would have produced the first time the
    /// box is opened, so editing starts from something that runs rather than
    /// from a blank page.
    pub advanced_text: String,
    /// Whether the run summary window is open.
    pub summary_open: bool,
    /// Whether the "really close it?" prompt is showing for that window.
    ///
    /// The summary is the only record of a finished run that lives in the
    /// application rather than on disk, and its close button sits next to the
    /// Save buttons. Closing it by accident after a long calculation loses the
    /// one view of what came back, so it asks first. Nothing else in Beavyr
    /// does, and nothing else should: this window is the exception because
    /// what it holds cannot be got back without running the job again.
    pub summary_confirm_close: bool,
    /// Valence warnings are shown only once a run has been attempted, matching
    /// the optimizer panel: unsolicited warnings about a perfectly ordinary
    /// molecule read as noise.
    pub show_warnings: bool,
}

impl Default for QcPanelState {
    fn default() -> Self {
        let program = QcProgram::Orca;
        let mut config = MethodConfig::default();
        crate::qchem_interfaces::method::adopt_program_defaults(&mut config, program);
        Self {
            open: false,
            program,
            charge: 0,
            multiplicity: 1,
            method: crate::qchem_interfaces::job::default_method(program),
            config,
            // The job most runs start from.
            job: JobType::Optimize,
            excited: ExcitedStateOptions::default(),
            wavefunction: false,
            advanced: false,
            advanced_text: String::new(),
            summary_open: false,
            summary_confirm_close: false,
            show_warnings: false,
        }
    }
}

impl QcPanelState {
    /// Whether the level-of-theory rows apply at all.
    ///
    /// False for DREIDING, whose parameters follow from the perceived atom
    /// types rather than from anything the user picks.
    #[cfg(test)]
    pub fn has_level_of_theory(&self) -> bool {
        self.method.is_some()
    }

    /// Whether a functional has to be chosen for the current method.
    pub fn needs_functional(&self) -> bool {
        self.method.is_some_and(|m| m.needs_functional())
    }

    /// Whether a basis set has to be chosen.
    ///
    /// A composite functional takes the field away: it brings its own.
    pub fn needs_basis(&self) -> bool {
        let Some(method) = self.method else {
            return false;
        };
        if !method.needs_basis() {
            return false;
        }
        if method.needs_functional()
            && self.program == QcProgram::Orca
            && crate::qchem_interfaces::orca_method::is_composite(&self.config.functional)
        {
            return false;
        }
        true
    }

    /// Whether the excited-state rows apply.
    pub fn needs_excited_state_options(&self) -> bool {
        self.job.is_excited_state()
    }

    /// The name the input is written under, for the save dialog's suggestion.
    pub fn input_file_name(&self) -> &'static str {
        match self.program {
            QcProgram::Orca => "orca_job.inp",
            QcProgram::Psi4 => "input.dat",
            QcProgram::Psi4Py => "run_psi4.py",
            QcProgram::PySCF => "run_pyscf.py",
            _ => "input.txt",
        }
    }

    /// A one-line description of what will run, for the panel header.
    pub fn summary_line(&self) -> String {
        let method = match self.method {
            None => "force field".to_string(),
            Some(QcMethod::Dft) => {
                let functional = match self.program {
                    QcProgram::Orca => {
                        crate::qchem_interfaces::orca_method::functional_label(
                            &self.config.functional,
                        )
                    }
                    QcProgram::PySCF => {
                        crate::qchem_interfaces::pyscf_method::functional_label(
                            &self.config.functional,
                        )
                    }
                    _ => &self.config.functional,
                };
                functional.to_string()
            }
            Some(other) => other.label().to_string(),
        };
        if self.needs_basis() {
            format!(
                "{} {}/{} \u{2014} {}",
                self.program.label(),
                method,
                self.config.basis,
                self.job.label()
            )
        } else {
            format!("{} {method} \u{2014} {}", self.program.label(), self.job.label())
        }
    }
}
