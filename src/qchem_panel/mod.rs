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

pub mod run;
pub mod sections;
pub mod summary;
pub mod ui;
pub mod ui_classic;
pub mod ui_sections;

use bevy::prelude::*;

use crate::qchem_interfaces::job::{ExcitedStateOptions, JobType, QcMethod};
use crate::qchem_interfaces::method::MethodConfig;
use crate::qchem_interfaces::program::QcProgram;

/// Which arrangement of the runner panel is drawn.
///
/// The sectioned layout replaced a flat eleven-row form. The old one stays
/// reachable rather than being deleted, because a redesign is a claim about
/// what someone finds usable, and a claim like that should be testable by the
/// person it is made for rather than argued with.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PanelLayout {
    /// Six collapsed sections, each showing its current setting.
    #[default]
    Sections,
    /// The original flat form, every control visible at once.
    Classic,
}

impl PanelLayout {
    pub fn label(self) -> &'static str {
        match self {
            PanelLayout::Sections => "sections",
            PanelLayout::Classic => "the old layout",
        }
    }

    /// The other one, for a control that toggles.
    pub fn other(self) -> PanelLayout {
        match self {
            PanelLayout::Sections => PanelLayout::Classic,
            PanelLayout::Classic => PanelLayout::Sections,
        }
    }
}

/// Everything the panel remembers between frames.
#[derive(Resource)]
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
    /// Which arrangement is drawn. Remembered between sessions.
    pub layout: PanelLayout,
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
            layout: crate::qchem_interfaces::config::load_layout(),
            summary_open: false,
            summary_confirm_close: false,
            show_warnings: false,
        }
    }
}

impl QcPanelState {
    /// Switches program, resetting everything that does not carry across.
    ///
    /// The method, functional and basis all reset, because the catalogues do
    /// not translate between codes: ORCA's `M062X` is an input error to
    /// Behemoth and Behemoth's `m06-2x` is one to ORCA. The charge, the
    /// multiplicity and the resources stay, because they describe the molecule
    /// and the machine rather than the program.
    ///
    /// The job stays too where the new program can do it, and falls back to an
    /// optimisation where it cannot -- which every backend can do.
    pub fn set_program(&mut self, program: QcProgram) {
        if program == self.program {
            return;
        }
        self.program = program;
        self.method = crate::qchem_interfaces::job::default_method(program);
        crate::qchem_interfaces::method::adopt_program_defaults(&mut self.config, program);
        if !crate::qchem_interfaces::job::supports(program, self.job) {
            self.job = JobType::Optimize;
        }
    }

    /// Whether the level-of-theory rows apply at all.
    ///
    /// False for DREIDING, whose parameters follow from the perceived atom
    /// types rather than from anything the user picks.
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

    /// Whether this program is genuinely driven by an input file.
    ///
    /// Four are: ORCA by its `.inp`, Psi4 by an input file that is itself
    /// Python, and PySCF by a Python script -- which is simply how PySCF is
    /// used. Each of those is a file a chemist would recognise, write by hand
    /// and run outside Beavyr.
    ///
    /// Three are not. xTB and Behemoth take everything on their command lines,
    /// and DREIDING runs in this process. Sparrow belongs with them: it is
    /// driven by settings, through its own command line or through SCINE's
    /// bindings, and the Python script Beavyr generates for it is this
    /// program's scaffolding rather than Sparrow's input format. Offering that
    /// scaffolding as "the input" would invite editing something that is not a
    /// thing Sparrow has.
    pub fn has_input_file(&self) -> bool {
        matches!(
            self.program,
            QcProgram::Orca | QcProgram::Psi4 | QcProgram::Psi4Py | QcProgram::PySCF
        )
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

    /// The input the current settings would write.
    ///
    /// One function for three callers: the preview the editor starts from, the
    /// Save button, and the run itself. They cannot disagree about what would
    /// have been sent, which is the point -- an editor showing something other
    /// than what runs would be worse than no editor.
    pub fn generated_input(&self, atoms: &[String], positions_angstrom: &[f64]) -> String {
        use crate::qchem_interfaces::{orca_run, psi4_run, pyscf_run};

        let method = self.method.unwrap_or(QcMethod::Dft);
        match self.program {
            QcProgram::Orca => orca_run::panel_input_text(
                self.job,
                method,
                self.charge,
                self.multiplicity,
                &self.config,
                self.excited,
                &self.config.orca_extra,
                &orca_run::xyz_body(atoms, positions_angstrom),
            ),
            QcProgram::Psi4 => psi4_run::input_text(
                self.job,
                method,
                self.charge,
                self.multiplicity,
                &self.config,
                self.excited,
                self.wavefunction,
                atoms,
                positions_angstrom,
            ),
            QcProgram::Psi4Py => psi4_run::script_text(
                self.job,
                method,
                self.charge,
                self.multiplicity,
                &self.config,
                self.excited,
                self.wavefunction,
                atoms,
                positions_angstrom,
            ),
            QcProgram::PySCF => pyscf_run::script_text(
                self.job,
                self.charge,
                self.multiplicity,
                &self.config,
                self.excited,
                self.wavefunction,
                atoms,
                positions_angstrom,
            ),
            // Everything else is driven by settings rather than a file, so
            // there is nothing to show. See `has_input_file`.
            _ => String::new(),
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

#[cfg(test)]
mod tests {
    use super::*;

    /// The new layout is what a user gets unless they have chosen otherwise.
    #[test]
    fn the_new_layout_is_the_default() {
        assert_eq!(PanelLayout::default(), PanelLayout::Sections);
    }

    /// The choice survives a restart. Being dropped back into a layout you
    /// rejected is worse than having no switch at all, and the existing
    /// classic/windowed switch already has that wart; this one must not.
    #[test]
    fn the_layout_choice_is_remembered() {
        use crate::qchem_interfaces::config;

        let dir = std::env::temp_dir().join(format!(
            "beavyr_layout_test_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        // SAFETY: the suite runs single-threaded via --test-threads=1 and no
        // other test in this crate touches these variables concurrently.
        unsafe {
            std::env::remove_var("XDG_CONFIG_HOME");
            std::env::set_var("HOME", &dir);
        }

        // Nothing saved yet: the default.
        assert_eq!(config::load_layout(), PanelLayout::Sections);

        config::save_layout(PanelLayout::Classic);
        assert_eq!(config::load_layout(), PanelLayout::Classic);

        config::save_layout(PanelLayout::Sections);
        assert_eq!(config::load_layout(), PanelLayout::Sections);

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Review Focus 4. Switching program must update every header at once. A
    /// header still showing the previous program's basis is worse than no
    /// header, because it reads as current.
    #[test]
    fn switching_program_updates_every_header() {
        use crate::qchem_interfaces::job::QcMethod;
        use crate::qchem_interfaces::program::QcProgram;
        use sections::{header, Section};

        let mut state = QcPanelState::default();
        state.set_program(QcProgram::Orca);
        state.method = Some(QcMethod::Dft);
        state.config.functional = "M062X".to_string();
        let before = header(Section::LevelOfTheory, &state);
        assert!(before.contains("M06-2X"), "{before}");

        state.set_program(QcProgram::PySCF);
        let after = header(Section::LevelOfTheory, &state);
        assert!(
            !after.contains("M06-2X"),
            "ORCA's spelling must not survive the switch: {after}"
        );
        assert!(header(Section::Program, &state).contains("PySCF"));
    }

    /// The switch offers the other one, whichever you are on.
    #[test]
    fn the_switch_offers_the_other_layout() {
        assert_eq!(PanelLayout::Sections.other(), PanelLayout::Classic);
        assert_eq!(PanelLayout::Classic.other(), PanelLayout::Sections);
        assert!(!PanelLayout::Sections.label().is_empty());
        assert!(!PanelLayout::Classic.label().is_empty());
    }

    /// A job the chosen program cannot do must explain itself on screen.
    ///
    /// The reason currently lives in a hover tooltip. A user who cannot see
    /// why an option is unavailable is left guessing, and a tooltip is
    /// unreachable the moment this panel is used with a finger.
    #[test]
    fn an_unavailable_job_explains_itself_without_hovering() {
        use crate::qchem_interfaces::job::{self, JobType};
        use crate::qchem_interfaces::program::QcProgram;

        // The model already carries the reason; this test pins that it is
        // available to draw on screen rather than only as hover text.
        let reason =
            job::unsupported_reason(QcProgram::Xtb, JobType::Irc).expect("xTB has no IRC");
        assert!(reason.contains("ORCA") || reason.contains("Psi4"), "{reason}");

        // And that every refusal is short enough to show inline.
        for program in QcProgram::ALL {
            for a_job in JobType::ALL {
                if let Some(text) = job::unsupported_reason(program, a_job) {
                    assert!(
                        text.chars().count() < 200,
                        "{program:?}/{a_job:?} is too long to show inline: {} chars",
                        text.chars().count()
                    );
                }
            }
        }
    }
}
