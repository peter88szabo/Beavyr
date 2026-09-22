//! What calculation to run, and which programs can run it.
//!
//! The general panel asks two questions the single-purpose panels never had to:
//! *what kind of job* and *at what level of theory above the functional*. Both
//! answers depend on the program, and neither is uniform -- ORCA optimises a
//! transition state and follows a reaction coordinate, xTB does neither; PySCF
//! does real TD-DFT, Behemoth has only Grimme's simplified methods.
//!
//! Everything here was established by running the programs, not from their
//! manuals. Where a program cannot do something the panel says so and disables
//! it, rather than hiding it: someone choosing between codes is entitled to see
//! that ORCA does an IRC and xTB does not.

use super::program::QcProgram;

/// Which iteration history a finished run's summary lists.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SummaryKind {
    /// The geometry cycles: energy, change, step and gradient per cycle, with
    /// the convergence thresholds beside them.
    GeometryCycles,
    /// The self-consistent field iterations only, for a job that never moved
    /// the nuclei.
    ScfOnly,
}

/// Where a finished job's results belong.
///
/// A job hands its result to whichever viewer already owns that kind of thing,
/// rather than the panel growing its own copy of the trajectory player, the
/// mode list and the spectrum plot. The panel reads this to decide which window
/// to open when a run finishes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResultTarget {
    /// The trajectory player, and the energy-versus-step plot with it.
    Trajectory,
    /// The Vibrations panel: modes, thermochemistry, infrared spectrum.
    Vibrations,
    /// The UV-Vis panel: excited states and the absorption spectrum.
    UvVis,
}

/// The kind of calculation to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JobType {
    /// One energy at the structure as it stands.
    SinglePoint,
    /// One energy and the Cartesian gradient, written to a file you can keep.
    Gradient,
    /// The force-constant matrix at the structure as it stands, written to a
    /// file you can keep.
    ///
    /// Separate from [`JobType::Frequencies`] on purpose. The calculation is
    /// the same; what differs is what you get. This one hands you the matrix to
    /// save, reuse or feed to something else. That one analyses it, and opens
    /// the Vibrations panel with the modes in it.
    Hessian,
    Optimize,
    Frequencies,
    OptimizeThenFrequencies,
    /// A saddle-point search, not a minimisation.
    TransitionState,
    /// Intrinsic reaction coordinate, downhill from a transition state.
    Irc,
    /// Excited states and an absorption spectrum.
    TdDft,
}

impl JobType {
    pub const ALL: [JobType; 9] = [
        JobType::SinglePoint,
        JobType::Gradient,
        JobType::Hessian,
        JobType::Optimize,
        JobType::Frequencies,
        JobType::OptimizeThenFrequencies,
        JobType::TransitionState,
        JobType::Irc,
        JobType::TdDft,
    ];

    pub fn label(self) -> &'static str {
        match self {
            JobType::SinglePoint => "Single-point energy",
            JobType::Gradient => "Gradient",
            JobType::Hessian => "Hessian (save the matrix)",
            JobType::Optimize => "Optimization",
            JobType::Frequencies => "Frequencies",
            JobType::OptimizeThenFrequencies => "Optimization + Frequencies",
            JobType::TransitionState => "Transition-state optimization",
            JobType::Irc => "IRC",
            JobType::TdDft => "TD-DFT",
        }
    }

    /// Which viewers this job's results are handed to when it finishes.
    ///
    /// Empty for a single point, whose whole result is its energy and the
    /// convergence history behind it, both of which the run summary shows.
    pub fn result_targets(self) -> &'static [ResultTarget] {
        match self {
            // These three hand you a number or a file rather than something to
            // plot, so nothing opens; the summary offers the file to save.
            JobType::SinglePoint | JobType::Gradient | JobType::Hessian => &[],
            JobType::Optimize | JobType::TransitionState | JobType::Irc => {
                &[ResultTarget::Trajectory]
            }
            JobType::Frequencies => &[ResultTarget::Vibrations],
            // Both, and in this order: the trajectory is what happened, the
            // frequencies are the answer at the end of it.
            JobType::OptimizeThenFrequencies => {
                &[ResultTarget::Trajectory, ResultTarget::Vibrations]
            }
            JobType::TdDft => &[ResultTarget::UvVis],
        }
    }

    /// What the run summary should show for this job.
    ///
    /// Every job gets a summary. Which iterations it lists differs: a geometry
    /// job has a cycle table, a single point has only its self-consistent field
    /// history, and both end with the energetics.
    pub fn summary_kind(self) -> SummaryKind {
        match self {
            JobType::SinglePoint
            | JobType::Gradient
            | JobType::Hessian
            | JobType::TdDft => SummaryKind::ScfOnly,
            JobType::Optimize
            | JobType::OptimizeThenFrequencies
            | JobType::TransitionState
            | JobType::Irc => SummaryKind::GeometryCycles,
            JobType::Frequencies => SummaryKind::ScfOnly,
        }
    }

    /// One line under the selector, saying what the job is for.
    pub fn description(self) -> &'static str {
        match self {
            JobType::SinglePoint => {
                "One energy at the structure as it stands, with the self-consistent field \
                 history behind it."
            }
            JobType::Gradient => {
                "One energy and the forces on every atom, written to a file you can save \
                 and reuse."
            }
            JobType::Hessian => {
                "The force-constant matrix at the structure as it stands, written to a file \
                 you can save. Frequencies analyses the same matrix instead."
            }
            JobType::Optimize => "Relax the structure to the nearest minimum.",
            JobType::Frequencies => {
                "Second derivatives at the structure as it stands, for normal modes and \
                 thermochemistry."
            }
            JobType::OptimizeThenFrequencies => {
                "Relax first, then take frequencies at the relaxed structure -- which is \
                 where they mean something."
            }
            JobType::TransitionState => {
                "Search for a first-order saddle point. Start from a good guess: this \
                 converges to whatever saddle is nearest."
            }
            JobType::Irc => {
                "Follow the reaction path downhill from a transition state, to see which \
                 minima it connects."
            }
            JobType::TdDft => "Excited states and an absorption spectrum.",
        }
    }

    /// Whether the job produces a trajectory the player and energy plot show.
    pub fn produces_trajectory(self) -> bool {
        self.result_targets().contains(&ResultTarget::Trajectory)
    }

    /// Whether the job produces a Hessian the Vibrations panel shows.
    pub fn produces_frequencies(self) -> bool {
        self.result_targets().contains(&ResultTarget::Vibrations)
    }

    /// Whether the job needs the excited-state controls.
    pub fn is_excited_state(self) -> bool {
        self == JobType::TdDft
    }
}

/// The level of theory above the functional.
///
/// One enum rather than one per program, with a per-program list, so the panel
/// can be written once. A program simply does not list what it cannot do.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QcMethod {
    // Semi-empirical and parametrised.
    Gfn1,
    Gfn2,
    GfnFf,
    Tasi,
    /// One of Sparrow's semi-empirical parametrisations, named in the
    /// configuration's method field rather than here: there are eight of them
    /// and they differ only by parameter set.
    Semiempirical,
    // Wavefunction and density-functional.
    Hf,
    Dft,
    Mp2,
    RiMp2,
    Mp2F12,
    CcsdT,
    DlpnoCcsdT,
    CcsdTF12,
}

impl QcMethod {
    pub fn label(self) -> &'static str {
        match self {
            QcMethod::Gfn1 => "GFN1-xTB",
            QcMethod::Gfn2 => "GFN2-xTB",
            QcMethod::GfnFf => "GFN-FF",
            QcMethod::Tasi => "TASI",
            QcMethod::Semiempirical => "Semi-empirical",
            QcMethod::Hf => "HF",
            QcMethod::Dft => "DFT",
            QcMethod::Mp2 => "MP2",
            QcMethod::RiMp2 => "RI-MP2",
            QcMethod::Mp2F12 => "MP2-F12",
            QcMethod::CcsdT => "CCSD(T)",
            QcMethod::DlpnoCcsdT => "DLPNO-CCSD(T)",
            QcMethod::CcsdTF12 => "CCSD(T)-F12",
        }
    }

    /// Whether choosing this method means also choosing a functional.
    pub fn needs_functional(self) -> bool {
        // Sparrow's parametrisation travels in the same field a functional
        // does, so the row is shown for it too -- filled from Sparrow's own
        // catalogue rather than a functional list.
        matches!(self, QcMethod::Dft | QcMethod::Semiempirical)
    }

    /// Whether this method expands the wavefunction in a Gaussian basis, and so
    /// needs a basis set chosen for it.
    ///
    /// The semi-empirical methods carry their own element parameters and take
    /// no basis, which is why the field disappears for them rather than being
    /// shown and ignored.
    pub fn needs_basis(self) -> bool {
        !matches!(
            self,
            QcMethod::Gfn1
                | QcMethod::Gfn2
                | QcMethod::GfnFf
                | QcMethod::Tasi
                | QcMethod::Semiempirical
        )
    }

    /// Whether this method needs a correlation-fitting auxiliary basis.
    ///
    /// Everything correlated does. ORCA refuses the input without one, which
    /// is why this is supplied rather than left to the user.
    pub fn needs_auxiliary_basis(self) -> bool {
        matches!(
            self,
            QcMethod::Mp2
                | QcMethod::RiMp2
                | QcMethod::Mp2F12
                | QcMethod::CcsdT
                | QcMethod::DlpnoCcsdT
                | QcMethod::CcsdTF12
        )
    }

    /// Whether this is an explicitly correlated method, which additionally
    /// needs a complementary auxiliary basis and an F12-specific orbital basis.
    pub fn is_f12(self) -> bool {
        matches!(self, QcMethod::Mp2F12 | QcMethod::CcsdTF12)
    }

    /// Whether a correlated method's cost is worth warning about before the
    /// run rather than after it.
    pub fn is_expensive(self) -> bool {
        matches!(
            self,
            QcMethod::CcsdT | QcMethod::DlpnoCcsdT | QcMethod::CcsdTF12 | QcMethod::Mp2F12
        )
    }
}

/// The methods `program` offers, in the order the selector lists them.
///
/// Each was run against the installed program. DREIDING offers none: it is a
/// force field whose parameters follow from the perceived atom types, so the
/// whole level-of-theory block disappears for it.
pub fn methods_for(program: QcProgram) -> &'static [QcMethod] {
    match program {
        QcProgram::Xtb => &[QcMethod::Gfn2, QcMethod::Gfn1, QcMethod::GfnFf],
        QcProgram::Behemoth => &[
            QcMethod::Gfn1,
            QcMethod::Tasi,
            QcMethod::Hf,
            QcMethod::Dft,
            QcMethod::Mp2,
            QcMethod::RiMp2,
        ],
        QcProgram::Dreiding => &[],
        QcProgram::Orca => &[
            QcMethod::Hf,
            QcMethod::Dft,
            QcMethod::Mp2,
            QcMethod::RiMp2,
            QcMethod::Mp2F12,
            QcMethod::CcsdT,
            QcMethod::DlpnoCcsdT,
            QcMethod::CcsdTF12,
        ],
        QcProgram::PySCF => &[QcMethod::Hf, QcMethod::Dft, QcMethod::Mp2],
        // Both Psi4 routes are the same engine, so the same list.
        QcProgram::SparrowPy => &[QcMethod::Semiempirical],
        QcProgram::Psi4 | QcProgram::Psi4Py => &[
            QcMethod::Hf,
            QcMethod::Dft,
            QcMethod::Mp2,
            QcMethod::CcsdT,
        ],
    }
}

/// The method the selector starts on for `program`.
pub fn default_method(program: QcProgram) -> Option<QcMethod> {
    methods_for(program).first().copied()
}

/// Why a program cannot run a job, or `None` when it can.
///
/// The message is what the disabled entry says, so it has to name the program
/// and the alternative rather than just refusing.
pub fn unsupported_reason(program: QcProgram, job: JobType) -> Option<&'static str> {
    match (program, job) {
        // Every backend gives an energy, a gradient, a Hessian, and optimises.
        (_, JobType::SinglePoint)
        | (_, JobType::Gradient)
        | (_, JobType::Hessian)
        | (_, JobType::Optimize)
        | (_, JobType::Frequencies)
        | (_, JobType::OptimizeThenFrequencies) => None,

        // ORCA and Psi4 both search for saddle points and follow reaction
        // paths -- Psi4 through optking, whose `opt_type` takes `ts` and `irc`
        // and which ships the reaction-path modules to go with them.
        (QcProgram::Orca, JobType::TransitionState)
        | (QcProgram::Orca, JobType::Irc)
        | (QcProgram::Psi4, JobType::TransitionState)
        | (QcProgram::Psi4, JobType::Irc)
        | (QcProgram::Psi4Py, JobType::TransitionState)
        | (QcProgram::Psi4Py, JobType::Irc) => None,
        // Sparrow has neither, and no excited states either.
        (QcProgram::SparrowPy, JobType::TransitionState)
        | (QcProgram::SparrowPy, JobType::Irc)
        | (QcProgram::SparrowPy, JobType::TdDft) => Some(
            "Sparrow is a set of semi-empirical parametrisations: it gives an energy, a \
             gradient and a Hessian, and nothing else. Use ORCA or Psi4 for this.",
        ),
        (_, JobType::TransitionState) => Some(
            "Only ORCA and Psi4 have a transition-state optimizer. For the others, use the TS \
             Generation tool to make a guess, then optimize it with ORCA or Psi4.",
        ),
        (_, JobType::Irc) => Some(
            "Only ORCA and Psi4 can follow an intrinsic reaction coordinate.",
        ),

        // Real TD-DFT: ORCA, PySCF and Psi4.
        (QcProgram::Orca, JobType::TdDft)
        | (QcProgram::PySCF, JobType::TdDft)
        | (QcProgram::Psi4, JobType::TdDft)
        | (QcProgram::Psi4Py, JobType::TdDft) => None,
        (QcProgram::Behemoth, JobType::TdDft) => Some(
            "Behemoth has Grimme's simplified methods (sTDA and sTD-DFT), not full TD-DFT. \
             Run those from the UV-Vis panel, or choose ORCA or PySCF here for real TD-DFT.",
        ),
        (_, JobType::TdDft) => Some(
            "This program has no excited-state module. Use ORCA, Psi4 or PySCF for TD-DFT.",
        ),
    }
}

/// Whether `program` can run `job`.
pub fn supports(program: QcProgram, job: JobType) -> bool {
    unsupported_reason(program, job).is_none()
}

/// Excited-state settings, shown only for a TD-DFT job.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExcitedStateOptions {
    /// How many roots to solve for.
    pub states: u32,
    /// Whether to use the Tamm-Dancoff approximation.
    ///
    /// On by default, which is both programs' own default and the cheaper,
    /// more robust choice; turning it off gives the full linear response.
    pub tda: bool,
}

impl Default for ExcitedStateOptions {
    fn default() -> Self {
        ExcitedStateOptions { states: 10, tda: true }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The support table, asserted so a program gaining or losing a capability
    /// has exactly one place to change. Every entry was established by running
    /// the program.
    #[test]
    fn the_support_table_is_what_the_programs_actually_do() {
        // Every backend gives an energy, optimises and takes a Hessian.
        for program in QcProgram::ALL {
            for job in [
                JobType::SinglePoint,
                JobType::Gradient,
                JobType::Hessian,
                JobType::Optimize,
                JobType::Frequencies,
                JobType::OptimizeThenFrequencies,
            ] {
                assert!(supports(program, job), "{program:?} / {job:?}");
            }
        }

        // ORCA alone does saddle points and reaction paths. Verified by running
        // `! OptTS` and `! IRC` on water.
        for job in [JobType::TransitionState, JobType::Irc] {
            for program in [QcProgram::Orca, QcProgram::Psi4, QcProgram::Psi4Py] {
                assert!(supports(program, job), "{program:?} / {job:?}");
            }
            for program in [
                QcProgram::Xtb,
                QcProgram::Behemoth,
                QcProgram::Dreiding,
                QcProgram::PySCF,
                QcProgram::SparrowPy,
            ] {
                assert!(!supports(program, job), "{program:?} / {job:?}");
            }
        }

        // Real TD-DFT: everything but the semi-empirical codes.
        assert!(supports(QcProgram::Orca, JobType::TdDft));
        assert!(supports(QcProgram::PySCF, JobType::TdDft));
        assert!(supports(QcProgram::Psi4, JobType::TdDft));
        assert!(supports(QcProgram::Psi4Py, JobType::TdDft));
        for program in [
            QcProgram::Xtb,
            QcProgram::Behemoth,
            QcProgram::Dreiding,
            QcProgram::SparrowPy,
        ] {
            assert!(!supports(program, JobType::TdDft), "{program:?}");
        }
    }

    /// Behemoth's refusal must name what it does have, and must not call it
    /// TD-DFT. The simplified methods are a different model, and labelling them
    /// as the more impressive thing they resemble is exactly the mistake this
    /// message exists to prevent.
    #[test]
    fn behemoth_is_told_apart_from_real_tddft() {
        let reason = unsupported_reason(QcProgram::Behemoth, JobType::TdDft)
            .expect("Behemoth has no full TD-DFT");
        assert!(reason.contains("sTDA"), "{reason}");
        assert!(reason.contains("simplified"), "{reason}");
        assert!(reason.contains("UV-Vis"), "it should say where to find them: {reason}");
    }

    /// Every refusal says what to do instead, rather than only refusing.
    #[test]
    fn every_refusal_names_an_alternative() {
        for program in QcProgram::ALL {
            for job in JobType::ALL {
                if let Some(reason) = unsupported_reason(program, job) {
                    assert!(
                        reason.contains("ORCA")
                            || reason.contains("PySCF")
                            || reason.contains("Psi4"),
                        "{program:?} / {job:?} refuses without an alternative: {reason}"
                    );
                }
            }
        }
    }

    /// A force field has no level of theory to choose, so the whole block goes.
    #[test]
    fn dreiding_offers_no_methods() {
        assert!(methods_for(QcProgram::Dreiding).is_empty());
        assert_eq!(default_method(QcProgram::Dreiding), None);
        // And every other program offers at least one.
        for program in QcProgram::ALL {
            if program != QcProgram::Dreiding {
                assert!(default_method(program).is_some(), "{program:?}");
            }
        }
    }

    /// Only DFT takes a functional; the semi-empirical methods take no basis.
    #[test]
    fn the_method_decides_which_fields_apply() {
        assert!(QcMethod::Dft.needs_functional());
        for method in [QcMethod::Hf, QcMethod::Mp2, QcMethod::Gfn2, QcMethod::Tasi] {
            assert!(!method.needs_functional(), "{method:?}");
        }
        for method in [QcMethod::Gfn1, QcMethod::Gfn2, QcMethod::GfnFf, QcMethod::Tasi] {
            assert!(!method.needs_basis(), "{method:?} carries its own parameters");
        }
        for method in [QcMethod::Hf, QcMethod::Dft, QcMethod::Mp2, QcMethod::CcsdT] {
            assert!(method.needs_basis(), "{method:?}");
        }
    }

    /// Every correlated method needs a fitting basis, and ORCA refuses the
    /// input without one. Nothing uncorrelated does.
    #[test]
    fn the_correlated_methods_ask_for_an_auxiliary_basis() {
        for method in [
            QcMethod::Mp2,
            QcMethod::RiMp2,
            QcMethod::Mp2F12,
            QcMethod::CcsdT,
            QcMethod::DlpnoCcsdT,
            QcMethod::CcsdTF12,
        ] {
            assert!(method.needs_auxiliary_basis(), "{method:?}");
        }
        for method in [QcMethod::Hf, QcMethod::Dft, QcMethod::Gfn2, QcMethod::Tasi] {
            assert!(!method.needs_auxiliary_basis(), "{method:?}");
        }
    }

    /// The F12 methods need more than an ordinary fitting basis, so they are
    /// distinguishable.
    #[test]
    fn the_f12_methods_are_marked() {
        assert!(QcMethod::Mp2F12.is_f12());
        assert!(QcMethod::CcsdTF12.is_f12());
        assert!(!QcMethod::Mp2.is_f12());
        assert!(!QcMethod::CcsdT.is_f12());
    }

    /// A job knows which viewer its result belongs to, so the panel does not
    /// have to carry that mapping separately.
    #[test]
    fn each_job_knows_what_it_produces() {
        assert!(JobType::Optimize.produces_trajectory());
        assert!(JobType::TransitionState.produces_trajectory());
        assert!(JobType::Irc.produces_trajectory());
        assert!(!JobType::Frequencies.produces_trajectory());

        assert!(JobType::Frequencies.produces_frequencies());
        assert!(JobType::OptimizeThenFrequencies.produces_frequencies());
        assert!(!JobType::Optimize.produces_frequencies());

        // The combined job does both, which is the whole point of it.
        assert!(JobType::OptimizeThenFrequencies.produces_trajectory());

        assert!(JobType::TdDft.is_excited_state());
        assert!(!JobType::TdDft.produces_trajectory());
    }

    /// Every job's results land somewhere a user can see them. A single point
    /// is the one exception, and deliberately so: its whole answer is the
    /// energy and the iterations behind it, which the run summary shows.
    #[test]
    fn every_job_hands_its_result_to_a_viewer() {
        for job in JobType::ALL {
            let targets = job.result_targets();
            if matches!(job, JobType::SinglePoint | JobType::Gradient | JobType::Hessian) {
                assert!(targets.is_empty(), "{job:?} hands over a file, not a plot");
            } else {
                assert!(!targets.is_empty(), "{job:?} produces nothing anyone can see");
            }
        }

        assert_eq!(JobType::Optimize.result_targets(), &[ResultTarget::Trajectory]);
        assert_eq!(JobType::Frequencies.result_targets(), &[ResultTarget::Vibrations]);
        assert_eq!(JobType::TdDft.result_targets(), &[ResultTarget::UvVis]);
        // The combined job opens both, trajectory first: it is what happened,
        // and the frequencies are the answer at the end of it.
        assert_eq!(
            JobType::OptimizeThenFrequencies.result_targets(),
            &[ResultTarget::Trajectory, ResultTarget::Vibrations]
        );
        // A saddle search and a reaction path are both journeys, so both play.
        assert_eq!(JobType::TransitionState.result_targets(), &[ResultTarget::Trajectory]);
        assert_eq!(JobType::Irc.result_targets(), &[ResultTarget::Trajectory]);
    }

    /// A job that never moves the nuclei has no cycle table to show, so its
    /// summary lists the self-consistent field iterations instead.
    #[test]
    fn the_summary_lists_the_iterations_the_job_actually_had() {
        assert_eq!(JobType::SinglePoint.summary_kind(), SummaryKind::ScfOnly);
        assert_eq!(JobType::Gradient.summary_kind(), SummaryKind::ScfOnly);
        assert_eq!(JobType::Hessian.summary_kind(), SummaryKind::ScfOnly);
        assert_eq!(JobType::TdDft.summary_kind(), SummaryKind::ScfOnly);
        assert_eq!(JobType::Frequencies.summary_kind(), SummaryKind::ScfOnly);
        for job in [
            JobType::Optimize,
            JobType::OptimizeThenFrequencies,
            JobType::TransitionState,
            JobType::Irc,
        ] {
            assert_eq!(job.summary_kind(), SummaryKind::GeometryCycles, "{job:?}");
        }
    }

    /// Every backend can do a single point. It is the cheapest thing any of
    /// them does, so a refusal would be a bug.
    #[test]
    fn every_backend_does_a_single_point() {
        for program in QcProgram::ALL {
            assert!(supports(program, JobType::SinglePoint), "{program:?}");
        }
    }

    /// Ten states with Tamm-Dancoff on: both programs' own default, and the
    /// cheaper, more robust choice for a first look.
    #[test]
    fn the_excited_state_defaults_are_sensible() {
        let options = ExcitedStateOptions::default();
        assert_eq!(options.states, 10);
        assert!(options.tda);
    }

    /// Every job is described, or the selector would show a bare name for
    /// something like IRC that not everyone will recognise.
    #[test]
    fn every_job_describes_itself() {
        for job in JobType::ALL {
            assert!(!job.label().is_empty(), "{job:?}");
            assert!(job.description().len() > 20, "{job:?} is barely described");
        }
    }

    /// `supports` and `unsupported_reason` are two views of one fact, and the
    /// panel uses both -- one to enable an entry, the other to explain a
    /// disabled one. They must never disagree.
    #[test]
    fn supports_and_the_reason_are_two_views_of_one_fact() {
        for program in QcProgram::ALL {
            for job in JobType::ALL {
                assert_eq!(
                    supports(program, job),
                    unsupported_reason(program, job).is_none(),
                    "{program:?} / {job:?}"
                );
            }
        }
    }
}
