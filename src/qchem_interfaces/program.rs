//! Which external program a calculation is sent to.
//!
//! Beavyr does not link any quantum-chemistry code in. Each backend is an
//! executable the user points at, run in a scratch directory and read back
//! from its files -- which keeps Beavyr's own build to three dependencies, and
//! means a program that crashes or is cancelled cannot take the viewport with
//! it.
//!
//! Adding a backend means adding a variant here, a path file, and the command
//! construction and output parsing for it. Nothing else in the panels changes.

/// An external program Beavyr can run a calculation with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QcProgram {
    /// Grimme's xTB, driven through its own command line.
    Xtb,
    /// Behemoth, the user's own Rust quantum-chemistry code, driven through
    /// `behemoth --xyz … --method xtb …`.
    Behemoth,
}

impl QcProgram {
    /// Every program, in the order the selector offers them.
    pub const ALL: [QcProgram; 2] = [QcProgram::Xtb, QcProgram::Behemoth];

    pub fn label(self) -> &'static str {
        match self {
            QcProgram::Xtb => "xTB",
            QcProgram::Behemoth => "Behemoth",
        }
    }

    /// The executable's usual file name, used as the file dialog's filter and
    /// as the hint in an empty path field.
    pub fn binary_hint(self) -> &'static str {
        match self {
            QcProgram::Xtb => "xtb",
            QcProgram::Behemoth => "behemoth",
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
        assert_eq!(QcProgram::ALL.len(), 2);
        for program in QcProgram::ALL {
            assert_eq!(
                QcProgram::ALL.iter().filter(|&&p| p == program).count(),
                1,
                "{program:?} appears twice"
            );
            assert!(!program.label().is_empty());
            assert!(!program.binary_hint().is_empty());
        }
    }

    /// Each program remembers its path in its own file, and xTB keeps the name
    /// earlier versions wrote so an existing configuration is not orphaned.
    #[test]
    fn each_program_has_its_own_config_file() {
        assert_eq!(QcProgram::Xtb.config_file(), "xtb_path.txt");
        let mut files: Vec<&str> = QcProgram::ALL.iter().map(|p| p.config_file()).collect();
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
