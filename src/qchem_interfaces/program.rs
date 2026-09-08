//! Which backend a calculation is sent to.
//!
//! Beavyr links no quantum-chemistry code in. Each quantum-chemistry backend is an executable the
//! user points at, run in a scratch directory and read back from its files -- which keeps
//! Beavyr's own build to five dependencies, and means a program that crashes or is cancelled
//! cannot take the viewport with it.
//!
//! The exception is [`QcProgram::Dreiding`], the DREIDING force field, which *is* part of Beavyr
//! and runs in this process. It is offered alongside the external programs because from a user's
//! point of view it answers the same question -- optimise this, give me its frequencies -- but it
//! needs no path, no scratch directory and no output parsing, so
//! [`QcProgram::runs_in_process`] exists to tell the two apart.
//!
//! Adding an external backend means adding a variant here, a path file, and the command
//! construction and output parsing for it. Nothing else in the panels changes.

/// An external program Beavyr can run a calculation with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
}

impl QcProgram {
    /// Every program, in the order the selector offers them.
    pub const ALL: [QcProgram; 3] = [QcProgram::Xtb, QcProgram::Behemoth, QcProgram::Dreiding];

    pub fn label(self) -> &'static str {
        match self {
            QcProgram::Xtb => "xTB",
            QcProgram::Behemoth => "Behemoth",
            QcProgram::Dreiding => "DREIDING",
        }
    }

    /// Whether this backend runs inside Beavyr rather than as a separate program.
    ///
    /// An in-process backend needs no executable path, no scratch directory and no output
    /// parsing, and cannot be cancelled by killing a child process -- so the panels use this to
    /// decide whether to show a path field at all.
    pub fn runs_in_process(self) -> bool {
        matches!(self, QcProgram::Dreiding)
    }

    /// The executable's usual file name, used as the file dialog's filter and
    /// as the hint in an empty path field.
    pub fn binary_hint(self) -> &'static str {
        match self {
            QcProgram::Xtb => "xtb",
            QcProgram::Behemoth => "behemoth",
            // Built in: there is no executable to find.
            QcProgram::Dreiding => "",
        }
    }

    /// The file its path is remembered in, under Beavyr's config directory.
    ///
    /// xTB keeps the name it has always used, so a configuration written by an
    /// earlier version is still found.
    pub fn config_file(self) -> &'static str {
        match self {
            QcProgram::Xtb => "xtb_path.txt",
            QcProgram::Behemoth => "behemoth_path.txt",
            // Nothing to remember.
            QcProgram::Dreiding => "",
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
        assert_eq!(QcProgram::ALL.len(), 3);
        for program in QcProgram::ALL {
            assert_eq!(
                QcProgram::ALL.iter().filter(|&&p| p == program).count(),
                1,
                "{program:?} appears twice"
            );
            assert!(!program.label().is_empty());
        }
    }

    /// Only the external programs need an executable, and only they need somewhere to remember
    /// its path. A built-in backend offering a path field would be a confusing dead end.
    #[test]
    fn only_external_programs_need_a_binary_and_a_path_file() {
        for program in QcProgram::ALL {
            if program.runs_in_process() {
                assert!(program.binary_hint().is_empty(), "{program:?}");
                assert!(program.config_file().is_empty(), "{program:?}");
            } else {
                assert!(!program.binary_hint().is_empty(), "{program:?}");
                assert!(!program.config_file().is_empty(), "{program:?}");
            }
        }
    }

    /// An in-process backend must resolve without a path, or the run button would be permanently
    /// disabled for a backend that needs no executable.
    #[test]
    fn an_in_process_backend_resolves_with_no_path() {
        use crate::qchem_interfaces::xtb_optimize::resolve_program_executable;

        let resolved = resolve_program_executable(QcProgram::Dreiding, "")
            .expect("a built-in backend needs no path");
        assert!(resolved.as_os_str().is_empty());

        // An external one still must be told where it is.
        assert!(resolve_program_executable(QcProgram::Xtb, "").is_err());
    }

    #[test]
    fn dreiding_is_the_only_in_process_backend() {
        assert!(QcProgram::Dreiding.runs_in_process());
        assert!(!QcProgram::Xtb.runs_in_process());
        assert!(!QcProgram::Behemoth.runs_in_process());
    }

    /// Each program remembers its path in its own file, and xTB keeps the name
    /// earlier versions wrote so an existing configuration is not orphaned.
    #[test]
    fn each_program_has_its_own_config_file() {
        assert_eq!(QcProgram::Xtb.config_file(), "xtb_path.txt");
        let mut files: Vec<&str> = QcProgram::ALL
            .iter()
            .filter(|p| !p.runs_in_process())
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
