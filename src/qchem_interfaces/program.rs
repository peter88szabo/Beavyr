//! Which backend a calculation is sent to.
//!
//! Beavyr links no quantum-chemistry code in. Each quantum-chemistry backend is an executable the
//! user points at, run in a scratch directory and read back from its files -- which keeps
//! Beavyr's own build to five dependencies, and means a program that crashes or is cancelled
//! cannot take the viewport with it.
//!
//! Two exceptions, for opposite reasons.
//!
//! [`QcProgram::Dreiding`], the DREIDING force field, *is* part of Beavyr and runs in this
//! process. It is offered alongside the external programs because from a user's point of view it
//! answers the same question -- optimise this, give me its frequencies -- but it needs no path,
//! no scratch directory and no output parsing.
//!
//! [`QcProgram::PySCF`] is a Python library with no command line at all, so there is no
//! executable to point at either. It still runs as a separate process, from a script Beavyr
//! generates; what it needs is a Python environment that can import it, which is the user's to
//! arrange before launching Beavyr.
//!
//! So "does this need a path field?" is [`QcProgram::needs_binary_path`], and is not the same
//! question as "does this run inside Beavyr?" ([`QcProgram::runs_in_process`]). What each backend
//! can actually do lives in [`super::capabilities`].
//!
//! Adding an external backend means adding a variant here, its capabilities, a path file where it
//! needs one, and the command construction and output parsing for it. Nothing else in the panels
//! changes.

/// An external program Beavyr can run a calculation with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[derive(serde::Serialize, serde::Deserialize)]
pub enum QcProgram {
    /// Grimme's xTB, driven through its own command line.
    Xtb,
    /// Behemoth, the user's own Rust quantum-chemistry code, driven through
    /// `behemoth --xyz … --method xtb …`.
    Behemoth,
    /// The DREIDING force field, built into Beavyr and run in this process.
    ///
    /// Not a quantum-chemistry method and not a substitute for one: it is a generic force field,
    /// good for cleaning up a structure and for a quick look at whether a geometry is a minimum,
    /// and it answers in milliseconds where the others take seconds to minutes.
    Dreiding,
    /// ORCA, driven by writing an input file and running the executable on it.
    Orca,
    /// PySCF, driven by generating a Python script and running it with `python3`
    /// from the environment Beavyr was launched in.
    PySCF,
    /// Psi4, driven through its own executable.
    ///
    /// That binary carries its own environment, so this route works whatever
    /// environment Beavyr was launched from -- which matters here, because Psi4
    /// and PySCF live in different ones on this machine and only one can be
    /// active at a time.
    Psi4,
    /// Psi4, driven as a Python module through `python3`.
    ///
    /// Needs Beavyr to have been launched from the environment holding Psi4.
    /// Same engine as [`QcProgram::Psi4`]; the difference is only how it is
    /// reached.
    Psi4Py,
    /// Sparrow's semi-empirical methods, through SCINE's Python bindings.
    ///
    /// The one backend with no optimiser and no frequency analysis of its own.
    /// It gives an energy, an analytic gradient and an analytic Hessian, and
    /// Beavyr does the rest with its own optimiser and its own projection and
    /// eigensolver.
    SparrowPy,
}

impl QcProgram {
    /// Every program, in the order the selector offers them.
    ///
    /// The three that were here first keep their positions, so an existing habit still lands on
    /// the same entry.
    pub const ALL: [QcProgram; 8] = [
        QcProgram::Xtb,
        QcProgram::Behemoth,
        QcProgram::Dreiding,
        QcProgram::Orca,
        QcProgram::PySCF,
        QcProgram::Psi4,
        QcProgram::Psi4Py,
        QcProgram::SparrowPy,
    ];

    pub fn label(self) -> &'static str {
        match self {
            QcProgram::Xtb => "xTB",
            QcProgram::Behemoth => "Behemoth",
            QcProgram::Dreiding => "DREIDING",
            QcProgram::Orca => "ORCA",
            QcProgram::PySCF => "PySCF",
            QcProgram::Psi4 => "Psi4",
            QcProgram::Psi4Py => "Psi4 (Python)",
            QcProgram::SparrowPy => "Sparrow (Python)",
        }
    }

    /// Whether this backend runs inside Beavyr rather than as a separate program.
    ///
    /// An in-process backend needs no scratch directory and no output parsing, and cannot be
    /// cancelled by killing a child process. Only DREIDING qualifies: PySCF needs no executable
    /// either, but it is still a separate process and is still killed to cancel it.
    pub fn runs_in_process(self) -> bool {
        matches!(self, QcProgram::Dreiding)
    }

    /// Whether the user has to point Beavyr at an executable for this backend.
    ///
    /// False for DREIDING, which is built in, and false for PySCF, which is reached through
    /// `python3` rather than through a program of its own. The panels use this to decide whether
    /// to draw a path field at all; a field that could never be filled usefully is a dead end.
    pub fn needs_binary_path(self) -> bool {
        super::capabilities::capabilities(self).launch.needs_binary_path()
    }

    /// The Python module this backend needs, when it is reached through Python.
    ///
    /// Used to check the environment before a run and to name what is missing when it is not
    /// there -- naming it matters, because "import failed" without the module and the interpreter
    /// looks like a broken installation rather than the wrong environment.
    pub fn python_module(self) -> Option<&'static str> {
        super::capabilities::capabilities(self).launch.python_module()
    }

    /// The executable's usual file name, used as the file dialog's filter and
    /// as the hint in an empty path field. Empty where there is no executable.
    pub fn binary_hint(self) -> &'static str {
        match self {
            QcProgram::Xtb => "xtb",
            QcProgram::Behemoth => "behemoth",
            QcProgram::Orca => "orca",
            QcProgram::Psi4 => "psi4",
            // Built in, or reached through Python: there is no executable to find.
            QcProgram::Dreiding
            | QcProgram::PySCF
            | QcProgram::Psi4Py
            | QcProgram::SparrowPy => "",
        }
    }

    /// The file its path is remembered in, under Beavyr's config directory.
    ///
    /// xTB keeps the name it has always used, so a configuration written by an
    /// earlier version is still found. Empty where there is no path to remember.
    pub fn config_file(self) -> &'static str {
        match self {
            QcProgram::Xtb => "xtb_path.txt",
            QcProgram::Behemoth => "behemoth_path.txt",
            QcProgram::Orca => "orca_path.txt",
            QcProgram::Psi4 => "psi4_path.txt",
            // Nothing to remember.
            QcProgram::Dreiding
            | QcProgram::PySCF
            | QcProgram::Psi4Py
            | QcProgram::SparrowPy => "",
        }
    }
}

impl Default for QcProgram {
    fn default() -> Self {
        QcProgram::Xtb
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_program_is_listed_once_and_described() {
        assert_eq!(QcProgram::ALL.len(), 8);
        for program in QcProgram::ALL {
            assert_eq!(
                QcProgram::ALL.iter().filter(|&&p| p == program).count(),
                1,
                "{program:?} appears twice"
            );
            assert!(!program.label().is_empty());
        }
    }

    /// Only a backend with an executable needs a hint and somewhere to remember its path. A
    /// backend offering a path field it cannot use would be a confusing dead end -- which is the
    /// case for PySCF as much as for DREIDING, though for a different reason.
    #[test]
    fn only_programs_with_an_executable_need_a_binary_and_a_path_file() {
        for program in QcProgram::ALL {
            if program.needs_binary_path() {
                assert!(!program.binary_hint().is_empty(), "{program:?}");
                assert!(!program.config_file().is_empty(), "{program:?}");
            } else {
                assert!(program.binary_hint().is_empty(), "{program:?}");
                assert!(program.config_file().is_empty(), "{program:?}");
            }
        }
    }

    /// PySCF needs no path but is not in process either: it is a separate program reached through
    /// Python. Conflating the two would either hide its cancellation or show it a dead path field.
    #[test]
    fn pyscf_needs_no_path_but_is_not_in_process() {
        assert!(!QcProgram::PySCF.needs_binary_path());
        assert!(!QcProgram::PySCF.runs_in_process());
        assert_eq!(QcProgram::PySCF.python_module(), Some("pyscf"));
    }

    /// A backend that resolves without a path must do so, or its Run button would be permanently
    /// disabled for want of something it does not need.
    #[test]
    fn a_backend_without_an_executable_resolves_with_no_path() {
        use crate::qchem_interfaces::xtb_optimize::resolve_program_executable;

        for program in [
            QcProgram::Dreiding,
            QcProgram::PySCF,
            QcProgram::Psi4Py,
            QcProgram::SparrowPy,
        ] {
            let resolved = resolve_program_executable(program, "")
                .unwrap_or_else(|e| panic!("{program:?} needs no path: {e}"));
            assert!(resolved.as_os_str().is_empty(), "{program:?}");
        }

        // One that does need an executable still must be told where it is.
        assert!(resolve_program_executable(QcProgram::Xtb, "").is_err());
        assert!(resolve_program_executable(QcProgram::Orca, "").is_err());
    }

    #[test]
    fn dreiding_is_the_only_in_process_backend() {
        assert!(QcProgram::Dreiding.runs_in_process());
        for program in QcProgram::ALL {
            if program != QcProgram::Dreiding {
                assert!(!program.runs_in_process(), "{program:?}");
            }
        }
    }

    /// Each program remembers its path in its own file, and xTB keeps the name
    /// earlier versions wrote so an existing configuration is not orphaned.
    #[test]
    fn each_program_has_its_own_config_file() {
        assert_eq!(QcProgram::Xtb.config_file(), "xtb_path.txt");
        let mut files: Vec<&str> = QcProgram::ALL
            .iter()
            .filter(|p| p.needs_binary_path())
            .map(|p| p.config_file())
            .collect();
        let count = files.len();
        files.sort_unstable();
        files.dedup();
        assert_eq!(files.len(), count, "two programs share a config file");
    }

    #[test]
    fn xtb_is_the_default() {
        assert_eq!(QcProgram::default(), QcProgram::Xtb);
    }
}
