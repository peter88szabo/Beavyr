//! What the runner panel is made of, and what each part currently says.
//!
//! No drawing here. Every rule the panel follows -- which sections exist, what
//! a header reads, whether a section applies to the chosen program -- is a
//! plain function with a test, so "a force field has no level of theory" is
//! something a test asserts rather than something only a click reveals.

use super::QcPanelState;
use crate::qchem_interfaces::program::QcProgram;

/// The longest a collapsed header may be.
///
/// Headers are built from strings the user supplies -- a hand-typed functional,
/// a deep path -- and those have no length limit. An unbounded header pushes
/// the expander off the row, and a section that cannot be opened is worse than
/// one with a truncated label.
pub const MAX_HEADER_CHARS: usize = 60;

/// One collapsible part of the panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Section {
    Program,
    LevelOfTheory,
    Job,
    Molecule,
    Resources,
    InputFile,
}

/// The sections, in the order they are drawn: what runs it, at what level of
/// theory, to do what, on what, with which resources, and the file that
/// results.
pub const ALL: [Section; 6] = [
    Section::Program,
    Section::LevelOfTheory,
    Section::Job,
    Section::Molecule,
    Section::Resources,
    Section::InputFile,
];

pub fn title(section: Section) -> &'static str {
    match section {
        Section::Program => "Program",
        Section::LevelOfTheory => "Level of theory",
        Section::Job => "Job",
        Section::Molecule => "Molecule",
        Section::Resources => "Resources",
        Section::InputFile => "Input file",
    }
}

/// Shortens `text` to fit a header, keeping the end.
///
/// The end, not the beginning: the distinctive part of a long path or a long
/// functional specification is at the end, and a row of leading directories
/// tells nobody anything.
fn fit(text: &str) -> String {
    let count = text.chars().count();
    if count <= MAX_HEADER_CHARS {
        return text.to_string();
    }
    let keep = MAX_HEADER_CHARS.saturating_sub(1);
    let skip = count - keep;
    format!("\u{2026}{}", text.chars().skip(skip).collect::<String>())
}

/// What this section's collapsed header reads for the current state.
pub fn header(section: Section, state: &QcPanelState) -> String {
    let text = match section {
        Section::Program => state.program.label().to_string(),
        Section::LevelOfTheory => {
            if !applies(section, state) {
                "none: a force field".to_string()
            } else if state.needs_basis() {
                format!("{}/{}", theory_name(state), state.config.basis)
            } else {
                theory_name(state)
            }
        }
        Section::Job => state.job.label().to_string(),
        Section::Molecule => format!(
            "charge {}, {}",
            state.charge,
            multiplicity_name(state.multiplicity)
        ),
        Section::Resources => format!(
            "{} thread{}, {} MB",
            state.config.nproc,
            if state.config.nproc == 1 { "" } else { "s" },
            state.config.memory_mb
        ),
        Section::InputFile => {
            if !applies(section, state) {
                "none: driven by settings".to_string()
            } else if state.advanced {
                "written by hand".to_string()
            } else {
                "built from the settings above".to_string()
            }
        }
    };
    fit(&text)
}

/// The level of theory as a name, without the basis set.
fn theory_name(state: &QcPanelState) -> String {
    use crate::qchem_interfaces::job::QcMethod;
    match state.method {
        None => "force field".to_string(),
        Some(QcMethod::Dft) => match state.program {
            QcProgram::Orca => {
                crate::qchem_interfaces::orca_method::functional_label(&state.config.functional)
                    .to_string()
            }
            QcProgram::PySCF => {
                crate::qchem_interfaces::pyscf_method::functional_label(&state.config.functional)
                    .to_string()
            }
            QcProgram::Psi4 | QcProgram::Psi4Py => {
                crate::qchem_interfaces::psi4_method::functional_label(&state.config.functional)
                    .to_string()
            }
            _ => state.config.functional.clone(),
        },
        Some(QcMethod::Semiempirical) => {
            crate::qchem_interfaces::sparrow_method::method_label(&state.config.functional)
                .to_string()
        }
        Some(other) => other.label().to_string(),
    }
}

/// A spin state named rather than numbered, because "singlet" is what a
/// chemist reads and "1" is what a program wants.
fn multiplicity_name(multiplicity: i32) -> String {
    match multiplicity {
        1 => "singlet".to_string(),
        2 => "doublet".to_string(),
        3 => "triplet".to_string(),
        4 => "quartet".to_string(),
        5 => "quintet".to_string(),
        other => format!("multiplicity {other}"),
    }
}

/// Whether this section has anything to set.
///
/// False is not "hide it". A section that does not apply stays on screen and
/// says why, because a form whose rows come and go moves under the cursor --
/// which is the complaint this redesign exists to answer.
pub fn applies(section: Section, state: &QcPanelState) -> bool {
    match section {
        Section::LevelOfTheory => state.has_level_of_theory(),
        Section::InputFile => state.has_input_file(),
        // The rest always have something to set.
        Section::Program | Section::Job | Section::Molecule | Section::Resources => true,
    }
}

/// Why a section has nothing to set, when it has nothing to set.
pub fn unavailable_note(section: Section, state: &QcPanelState) -> Option<&'static str> {
    if applies(section, state) {
        return None;
    }
    Some(match section {
        Section::LevelOfTheory => {
            "DREIDING is a generic force field. Its parameters follow from the atom types, \
             which are perceived from the structure, so there is no level of theory to set."
        }
        Section::InputFile => {
            "This program is driven by settings rather than by an input file, so there is \
             nothing to edit or save. xTB and Behemoth take everything on their command \
             lines, Sparrow is driven through SCINE, and DREIDING runs inside Beavyr."
        }
        _ => "Nothing to set here.",
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qchem_interfaces::job::{JobType, QcMethod};
    use crate::qchem_interfaces::program::QcProgram;
    use crate::qchem_panel::QcPanelState;

    fn orca() -> QcPanelState {
        let mut state = QcPanelState::default();
        state.set_program(QcProgram::Orca);
        state.method = Some(QcMethod::Dft);
        state.config.functional = "B3LYP".to_string();
        state.config.basis = "def2-SVP".to_string();
        state
    }

    /// Every section names itself, and the list is the six the spec gives.
    #[test]
    fn there_are_six_sections_and_each_is_named() {
        assert_eq!(ALL.len(), 6);
        for section in ALL {
            assert!(!title(section).is_empty(), "{section:?}");
        }
    }

    /// The header carries the current setting, which is what makes the whole
    /// panel readable without opening anything.
    #[test]
    fn the_header_carries_the_current_setting() {
        let state = orca();
        assert!(header(Section::Program, &state).contains("ORCA"));
        let theory = header(Section::LevelOfTheory, &state);
        assert!(theory.contains("B3LYP"), "{theory}");
        assert!(theory.contains("def2-SVP"), "{theory}");
        assert!(header(Section::Job, &state).contains("Optimization"));
    }

    /// Review Focus 1. A force field has no level of theory. The section must
    /// say so rather than disappear, or the form moves under the cursor --
    /// which is the complaint this redesign exists to fix.
    #[test]
    fn a_force_field_says_it_has_no_level_of_theory() {
        let mut state = QcPanelState::default();
        state.set_program(QcProgram::Dreiding);
        assert!(!applies(Section::LevelOfTheory, &state));
        let note = unavailable_note(Section::LevelOfTheory, &state).expect("a reason is given");
        assert!(note.contains("force field"), "{note}");
        // And it still has a header, because it is still on screen.
        assert!(!header(Section::LevelOfTheory, &state).is_empty());
    }

    /// Review Focus 2. Three programs are driven by settings, not by a file.
    #[test]
    fn programs_without_an_input_file_say_so() {
        let mut state = QcPanelState::default();
        for program in [QcProgram::Xtb, QcProgram::Behemoth, QcProgram::SparrowPy] {
            state.set_program(program);
            assert!(!applies(Section::InputFile, &state), "{program:?}");
            assert!(
                unavailable_note(Section::InputFile, &state).is_some(),
                "{program:?} must say why"
            );
        }
        state.set_program(QcProgram::Orca);
        assert!(applies(Section::InputFile, &state));
    }

    /// Review Focus 3. Only TD-DFT adds excited-state controls, and the Job
    /// header must read the same shape either way so the section does not
    /// change size as the job changes.
    #[test]
    fn the_job_header_reads_the_same_shape_for_every_job() {
        let mut state = orca();
        for job in [JobType::Optimize, JobType::TdDft, JobType::SinglePoint] {
            state.job = job;
            let text = header(Section::Job, &state);
            assert_eq!(text, job.label(), "{job:?}");
        }
    }

    /// Review Focus 5. Headers are built from strings the user supplies and
    /// those have no length limit. A header that grows without bound pushes
    /// the expander off the row and the section cannot be opened at all.
    #[test]
    fn a_header_is_bounded_however_long_the_setting_is() {
        let mut state = orca();
        state.config.functional = "x".repeat(400);
        state.config.basis = "y".repeat(400);
        for section in ALL {
            let text = header(section, &state);
            assert!(
                text.chars().count() <= MAX_HEADER_CHARS,
                "{section:?} header is {} chars",
                text.chars().count()
            );
        }
    }

    /// Every section that applies has a header, and every one that does not
    /// has a reason. Neither may be silent.
    #[test]
    fn nothing_is_silent() {
        let mut state = QcPanelState::default();
        for program in QcProgram::ALL {
            state.set_program(program);
            for section in ALL {
                assert!(!header(section, &state).is_empty(), "{program:?} {section:?}");
                if !applies(section, &state) {
                    assert!(
                        unavailable_note(section, &state).is_some(),
                        "{program:?} {section:?} is unavailable without saying why"
                    );
                }
            }
        }
    }
}
