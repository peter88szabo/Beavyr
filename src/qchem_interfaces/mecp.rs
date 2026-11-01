//! Settings and capability checks for two-spin-state MECP searches.
use super::{JobType, QcMethod};
use crate::optimizer::spin_crossing::SpinCrossingCoordinateMode;
use crate::qchem_interfaces::program::QcProgram;
#[derive(Debug, Clone, Copy, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Driver {
    #[default]
    Beavyr,
    Orca,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Coordinates {
    #[default]
    Cartesian,
    Primitive,
    Pulay,
}
impl Coordinates {
    pub fn mode(self) -> SpinCrossingCoordinateMode {
        match self {
            Self::Cartesian => SpinCrossingCoordinateMode::Cartesian,
            Self::Primitive => SpinCrossingCoordinateMode::SimpleInternal,
            Self::Pulay => SpinCrossingCoordinateMode::PulayRedundant,
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::Cartesian => "Cartesian",
            Self::Primitive => "Primitive internal",
            Self::Pulay => "Pulay redundant internal",
        }
    }
}
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct Options {
    pub driver: Driver,
    pub multiplicity_b: i32,
    pub max_cycles: usize,
    pub gap_tolerance: f64,
    /// SpinCross convergence units: Angstrom and Hartree/Angstrom.
    pub max_step: f64,
    pub rms_step: f64,
    pub max_gradient: f64,
    pub rms_gradient: f64,
    pub native_tolerances: Option<crate::qchem_interfaces::method::tolerances::Criteria>,
    pub coordinates: Coordinates,
}
impl Default for Options {
    fn default() -> Self {
        let optimizer = crate::optimizer::spin_crossing::SpinCrossingOptions::default();
        Self {
            driver: Default::default(),
            multiplicity_b: 3,
            max_cycles: 100,
            gap_tolerance: 5e-5,
            max_step: optimizer.max_step,
            rms_step: optimizer.rms_step,
            max_gradient: optimizer.max_gradient,
            rms_gradient: optimizer.rms_gradient,
            native_tolerances: None,
            coordinates: Default::default(),
        }
    }
}
pub fn program_problem(program: QcProgram) -> Option<&'static str> {
    if matches!(
        program,
        QcProgram::Orca
            | QcProgram::Behemoth
            | QcProgram::PySCF
            | QcProgram::Psi4
            | QcProgram::Psi4Py
    ) {
        None
    } else {
        Some("MECP requires two spin-dependent energy surfaces with gradients. Choose ORCA, Behemoth, PySCF or Psi4 with HF or DFT.")
    }
}
pub fn validate(
    program: QcProgram,
    method: Option<QcMethod>,
    options: &Options,
    mol: &crate::molecule::Molecule,
    charge: i32,
    multiplicity_a: i32,
) -> Result<(), String> {
    if let Some(reason) = program_problem(program) {
        return Err(reason.into());
    }
    if !matches!(method, Some(QcMethod::Hf | QcMethod::Dft)) {
        return Err("MECP currently supports HF and DFT energy/gradient surfaces.".into());
    }
    if multiplicity_a == options.multiplicity_b {
        return Err("State A and State B must have different spin multiplicities.".into());
    }
    for (name, multi) in [("A", multiplicity_a), ("B", options.multiplicity_b)] {
        crate::qchem_interfaces::valence::validate_electronic_state(mol, charge, multi)
            .map_err(|e| format!("State {name}: {e}"))?;
    }
    if options.max_cycles == 0 || !options.gap_tolerance.is_finite() || options.gap_tolerance <= 0.
    {
        return Err("MECP requires a positive step limit and energy-gap tolerance.".into());
    }
    if [
        options.max_step,
        options.rms_step,
        options.max_gradient,
        options.rms_gradient,
    ]
    .iter()
    .any(|v| !v.is_finite() || *v <= 0.)
    {
        return Err(
            "MECP displacement and gradient tolerances must be positive and finite.".into(),
        );
    }
    debug_assert!(super::supports(program, JobType::Gradient));
    if let Some(t) = &options.native_tolerances {
        t.validate()?;
    }
    Ok(())
}
