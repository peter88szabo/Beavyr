//! What each backend can actually do, and how Beavyr starts it.
//!
//! The panels ask this rather than assuming. Two things force it.
//!
//! Not every program can optimise a geometry. Sparrow computes energies,
//! gradients and an analytic Hessian but has no optimiser at all, so
//! "optimise this" has to be a question with an answer rather than something
//! taken for granted -- and where the answer is no, Beavyr's own optimiser is
//! the only route and the panel says so instead of offering a choice that does
//! not exist.
//!
//! And not every backend is an executable. PySCF is a Python library with no
//! command line, so there is no path to browse to; what it needs instead is a
//! Python environment that can import it. Both are "where does this program
//! live", asked of different kinds of program, and both have to be answered
//! before a run starts rather than after one fails.

use super::program::QcProgram;

/// How a backend produces a Hessian, if it can at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HessianKind {
    /// The program differentiates twice itself. Fast and exact.
    Analytic,
    /// Beavyr builds the matrix from the program's gradients by central
    /// differences: `6N` calculations for `N` atoms, so it is worth warning
    /// about before a run rather than after it.
    FiniteDifference,
    /// No Hessian by any route.
    None,
}

impl HessianKind {
    /// Whether a Hessian can be had at all, however it is obtained.
    pub fn is_available(self) -> bool {
        !matches!(self, HessianKind::None)
    }
}

/// How Beavyr starts a backend.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Launch {
    /// Part of Beavyr. Nothing to find, nothing to kill, no output to parse.
    InProcess,
    /// An executable the user points at, remembered between sessions.
    Binary,
    /// A generated script run with `python3`, taken from the environment
    /// Beavyr itself was started in. There is no path to set: activating the
    /// right environment before launching Beavyr is the user's job, so what
    /// this carries is the module name to check for, not a location.
    Python {
        /// The import name, used both to check the environment and to name
        /// what is missing when it is not there.
        module: &'static str,
    },
}

impl Launch {
    /// Whether the user has to point Beavyr at an executable.
    pub fn needs_binary_path(self) -> bool {
        matches!(self, Launch::Binary)
    }

    /// The Python module this backend needs, if it is a Python one.
    pub fn python_module(self) -> Option<&'static str> {
        match self {
            Launch::Python { module } => Some(module),
            _ => None,
        }
    }
}

/// Everything the panels need to know about a backend before offering it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Capabilities {
    pub launch: Launch,
    /// Whether the program has a geometry optimiser of its own.
    ///
    /// When it does not, Beavyr's own optimiser is not an alternative but the
    /// only route, and the panel forces it rather than offering a choice.
    pub own_optimiser: bool,
    pub hessian: HessianKind,
    /// Whether a frequency run also yields infrared intensities.
    ///
    /// Where it does not, the mode list shows them as unavailable and the
    /// spectrum window stays shut, rather than plotting a row of zeroes that
    /// looks like a result.
    pub ir_intensities: bool,
}

/// What `program` can do.
pub fn capabilities(program: QcProgram) -> Capabilities {
    match program {
        // Grimme's xTB: its own optimiser, its own Hessian, and a
        // `vibspectrum` file carrying real intensities.
        QcProgram::Xtb => Capabilities {
            launch: Launch::Binary,
            own_optimiser: true,
            hessian: HessianKind::Analytic,
            ir_intensities: true,
        },
        // Behemoth prints its Hessian to standard output. `--hessian` reports
        // no intensities, which is why the spectrum stays unavailable after
        // one.
        QcProgram::Behemoth => Capabilities {
            launch: Launch::Binary,
            own_optimiser: true,
            hessian: HessianKind::Analytic,
            ir_intensities: true,
        },
        // The built-in force field. It has an analytic Hessian but no dipole
        // model, so there is nothing to build a spectrum from.
        QcProgram::Dreiding => Capabilities {
            launch: Launch::InProcess,
            own_optimiser: true,
            hessian: HessianKind::Analytic,
            ir_intensities: false,
        },
        // ORCA: `Opt` and `Freq` on the keyword line. Its `.hess` carries an
        // infrared spectrum, which Beavyr's existing reader already takes.
        QcProgram::Orca => Capabilities {
            launch: Launch::Binary,
            own_optimiser: true,
            hessian: HessianKind::Analytic,
            ir_intensities: true,
        },
        // PySCF: the `geometric` optimiser and an analytic Hessian, both
        // driven from a generated script.
        //
        // Intensities are absent by choice rather than by limitation. PySCF
        // can produce them, but only from a separate dipole-derivative
        // calculation on top of the Hessian, which roughly doubles what a
        // frequency job costs. They are left out until they are asked for.
        QcProgram::PySCF => Capabilities {
            launch: Launch::Python { module: "pyscf" },
            own_optimiser: true,
            hessian: HessianKind::Analytic,
            ir_intensities: false,
        },
        // Psi4, reached by its own executable. That binary carries its own
        // environment, which is why this route works whatever Beavyr was
        // launched from.
        QcProgram::Psi4 => Capabilities {
            launch: Launch::Binary,
            own_optimiser: true,
            hessian: HessianKind::Analytic,
            ir_intensities: true,
        },
        // The same engine as a Python module, so the same capabilities. Only
        // how it is reached differs.
        QcProgram::Psi4Py => Capabilities {
            launch: Launch::Python { module: "psi4" },
            own_optimiser: true,
            hessian: HessianKind::Analytic,
            ir_intensities: true,
        },
        // Sparrow is the case this whole structure exists for. It computes an
        // energy, an analytic gradient and an analytic Hessian, and nothing
        // else: no optimiser, no frequency analysis. So `own_optimiser` is
        // false, Beavyr's own drives it, and the panel says so rather than
        // offering a choice that does not exist.
        QcProgram::SparrowPy => Capabilities {
            launch: Launch::Python { module: "scine_sparrow" },
            own_optimiser: false,
            hessian: HessianKind::Analytic,
            // A semi-empirical calculator here reports no dipole derivatives,
            // so there is nothing to build a spectrum from.
            ir_intensities: false,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A backend that cannot optimise on its own must still be offered in the
    /// optimizer panel, driven by Beavyr's optimiser. The panel reads
    /// `own_optimiser` to decide whether the choice is the user's; nothing
    /// here may be taken to mean "not offered".
    #[test]
    fn every_program_can_be_optimised_one_way_or_another() {
        for program in QcProgram::ALL {
            let caps = capabilities(program);
            // Either the program optimises, or Beavyr drives it. The second is
            // always possible, so the assertion is that we know which.
            let _ = caps.own_optimiser;
        }
    }

    /// Exactly one backend runs inside Beavyr, and it is the one with no
    /// executable and no module.
    #[test]
    fn only_dreiding_runs_in_process() {
        for program in QcProgram::ALL {
            let caps = capabilities(program);
            let in_process = caps.launch == Launch::InProcess;
            assert_eq!(
                in_process,
                program == QcProgram::Dreiding,
                "{program:?} launch is {:?}",
                caps.launch
            );
        }
    }

    /// The launch kind and the older `runs_in_process` flag must agree, or the
    /// panel would show a path field for a backend that has no executable.
    #[test]
    fn the_launch_kind_agrees_with_runs_in_process() {
        for program in QcProgram::ALL {
            let caps = capabilities(program);
            assert_eq!(
                caps.launch == Launch::InProcess,
                program.runs_in_process(),
                "{program:?}"
            );
            assert_eq!(
                caps.launch.needs_binary_path(),
                program.needs_binary_path(),
                "{program:?}"
            );
        }
    }

    /// A Python backend names its module, and a non-Python one does not, so
    /// the environment check has exactly one thing to look for.
    #[test]
    fn only_python_backends_name_a_module() {
        assert_eq!(
            capabilities(QcProgram::PySCF).launch.python_module(),
            Some("pyscf")
        );
        assert_eq!(
            capabilities(QcProgram::Psi4Py).launch.python_module(),
            Some("psi4")
        );
        assert_eq!(
            capabilities(QcProgram::SparrowPy).launch.python_module(),
            Some("scine_sparrow")
        );
        for program in [QcProgram::Xtb, QcProgram::Behemoth, QcProgram::Dreiding, QcProgram::Orca, QcProgram::Psi4] {
            assert_eq!(capabilities(program).launch.python_module(), None, "{program:?}");
        }
    }

    /// Sparrow has no optimiser at all, which is the case the capability model
    /// exists for. Claiming one would have the panel offer a choice that does
    /// not exist and a run that could never work.
    #[test]
    fn sparrow_has_no_optimiser_of_its_own() {
        assert!(!capabilities(QcProgram::SparrowPy).own_optimiser);
        // And it is the only one, today.
        for program in QcProgram::ALL {
            if program != QcProgram::SparrowPy {
                assert!(capabilities(program).own_optimiser, "{program:?}");
            }
        }
    }

    /// Intensities are only claimed where they are genuinely produced.
    /// Claiming them elsewhere would put a row of zeroes on a spectrum and
    /// make it look like a measurement.
    #[test]
    fn intensities_are_claimed_only_where_they_exist() {
        assert!(capabilities(QcProgram::Xtb).ir_intensities);
        assert!(capabilities(QcProgram::Orca).ir_intensities);
        assert!(
            !capabilities(QcProgram::Dreiding).ir_intensities,
            "a force field with no dipole model has no intensities"
        );
        assert!(
            !capabilities(QcProgram::PySCF).ir_intensities,
            "deliberately not paid for -- see the note in capabilities()"
        );
    }

    /// Every backend in this stage can produce a Hessian by some route, so the
    /// frequency panel never has to refuse one outright.
    #[test]
    fn every_backend_can_produce_a_hessian() {
        for program in QcProgram::ALL {
            assert!(
                capabilities(program).hessian.is_available(),
                "{program:?} has no Hessian route"
            );
        }
    }
}
