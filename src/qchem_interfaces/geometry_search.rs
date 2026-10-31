//! Settings for native and Beavyr-driven minimum and transition-state searches.
use super::{JobType, QcMethod};
use crate::optimizer::transition_state::{TransitionStateCoordinateMode, TransitionStateOptions};
use crate::qchem_interfaces::program::QcProgram;
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Driver {
    #[default]
    Native,
    Beavyr,
}
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Coordinates {
    Cartesian,
    #[default]
    Primitive,
    Pulay,
}
impl Coordinates {
    pub fn label(self) -> &'static str {
        match self {
            Self::Cartesian => "Cartesian",
            Self::Primitive => "Primitive internal",
            Self::Pulay => "Pulay redundant internal",
        }
    }
}
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum MinimumAlgorithm {
    #[default]
    Bfgs,
    SteepestDescent,
    Bb1,
    Bb2,
}
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct Minimum {
    pub coordinates: Coordinates,
    pub algorithm: MinimumAlgorithm,
    pub max_cycles: usize,
    pub energy_tol: f64,
    pub max_step: f64,
    pub rms_step: f64,
    pub max_gradient: f64,
    pub rms_gradient: f64,
    pub max_step_internal: f64,
    pub best_fit_iters: usize,
    pub best_fit_rms_tol: f64,
    pub min_line_search_scale: f64,
    pub g_cutoff: f64,
    pub redundant_penalty: f64,
}
impl Default for Minimum {
    fn default() -> Self {
        Self {
            coordinates: Coordinates::Primitive,
            algorithm: MinimumAlgorithm::Bfgs,
            max_cycles: 300,
            energy_tol: 5e-6,
            max_step: 4e-3,
            rms_step: 2e-3,
            max_gradient: 3e-4,
            rms_gradient: 1e-4,
            max_step_internal: 0.2,
            best_fit_iters: 5,
            best_fit_rms_tol: 1e-7,
            min_line_search_scale: 1e-4,
            g_cutoff: 1e-7,
            redundant_penalty: 1000.,
        }
    }
}
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Reaction {
    #[default]
    Lowest,
    Bond,
    Angle,
    Dihedral,
    Transfer,
    Mixed,
    Direction,
}
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct Options {
    pub driver: Driver,
    pub ts: TransitionStateOptions,
    pub minimum: Minimum,
    pub reaction: Reaction,
    pub reaction_atoms: String,
    pub reaction_coordinates: String,
    pub direction: String,
    pub extra_bonds: String,
    pub extra_angles: String,
    pub extra_dihedrals: String,
    pub orca_recalc_hessian: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            driver: Driver::Native,
            ts: TransitionStateOptions {
                coordinate_mode: TransitionStateCoordinateMode::Internal,
                trajectory_file: None,
                verbosity: 0,
                ..Default::default()
            },
            minimum: Default::default(),
            reaction: Default::default(),
            reaction_atoms: String::new(),
            reaction_coordinates: String::new(),
            direction: String::new(),
            extra_bonds: String::new(),
            extra_angles: String::new(),
            extra_dihedrals: String::new(),
            orca_recalc_hessian: 5,
        }
    }
}
pub fn active(job: JobType) -> bool {
    matches!(
        job,
        JobType::Optimize | JobType::TransitionState | JobType::OptimizeThenFrequencies
    )
}
pub fn local(program: QcProgram, job: JobType, o: &Options) -> bool {
    program != QcProgram::Behemoth
        && (job != JobType::TransitionState || ts_supported(program))
        && active(job)
        && (o.driver == Driver::Beavyr
            || (job == JobType::TransitionState
                && !matches!(
                    program,
                    QcProgram::Orca | QcProgram::Psi4 | QcProgram::Psi4Py
                )))
}
pub fn ts_supported(program: QcProgram) -> bool {
    matches!(
        program,
        QcProgram::Orca
            | QcProgram::Psi4
            | QcProgram::Psi4Py
            | QcProgram::PySCF
            | QcProgram::Behemoth
            | QcProgram::Xtb
    )
}
pub fn method_problem(
    program: QcProgram,
    method: Option<QcMethod>,
    job: JobType,
    o: &Options,
) -> Option<&'static str> {
    if job == JobType::TransitionState
        && program == QcProgram::Behemoth
        && method != Some(QcMethod::Gfn1)
    {
        return Some("Behemoth's native TS frontend currently supports GFN1-xTB. Choose GFN1-xTB or another program for this TS search.");
    }
    if program == QcProgram::Behemoth
        && job == JobType::TransitionState
        && o.reaction == Reaction::Direction
    {
        return Some("Behemoth's native TS frontend requires a bond, angle, dihedral, transfer or weighted reference; Cartesian direction vectors are available with Beavyr's driver.");
    }
    if local(program, job, o)
        && super::scan::method_problem(program, method, super::scan::Mode::Relaxed).is_some()
    {
        return Some("This geometry search requires a supported energy-and-gradient method.");
    }
    None
}

#[path = "geometry_search_coordinates.rs"]
mod coordinates;
pub use coordinates::{prepare_ts, remap_dummies, validate};
