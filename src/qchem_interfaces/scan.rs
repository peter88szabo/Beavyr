//! Persisted one-dimensional potential-energy scan settings and validation.
use super::QcMethod;
use crate::optimizer::constrained::{ConstraintCoordinate, ConstraintTarget};
use crate::qchem_interfaces::program::QcProgram;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Driver {
    #[default]
    Beavyr,
    Orca,
    Behemoth,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Mode {
    Rigid,
    #[default]
    Relaxed,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum BondMotion {
    #[default]
    Second,
    First,
    Symmetric,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum RelaxationCoordinates {
    Primitive,
    #[default]
    Pulay,
}
impl RelaxationCoordinates {
    pub fn optimizer_mode(self) -> crate::optimizer::constrained::ConstrainedCoordinateMode {
        match self {
            Self::Primitive => {
                crate::optimizer::constrained::ConstrainedCoordinateMode::PrimitiveInternal
            }
            Self::Pulay => crate::optimizer::constrained::ConstrainedCoordinateMode::PulayRedundant,
        }
    }
    pub fn keyword(self) -> &'static str {
        match self {
            Self::Primitive => "internal",
            Self::Pulay => "pulay-redundant-internal",
        }
    }
}
impl BondMotion {
    pub fn keyword(self) -> &'static str {
        match self {
            Self::Second => "second",
            Self::First => "first",
            Self::Symmetric => "symmetric",
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Coordinate {
    #[default]
    Bond,
    Angle,
    Dihedral,
}
impl Coordinate {
    pub fn label(self) -> &'static str {
        match self {
            Self::Bond => "Bond",
            Self::Angle => "Angle",
            Self::Dihedral => "Dihedral",
        }
    }
    pub fn count(self) -> usize {
        match self {
            Self::Bond => 2,
            Self::Angle => 3,
            Self::Dihedral => 4,
        }
    }
    pub fn unit(self) -> &'static str {
        if self == Self::Bond {
            "Å"
        } else {
            "°"
        }
    }
    pub fn keyword(self) -> &'static str {
        match self {
            Self::Bond => "B",
            Self::Angle => "A",
            Self::Dihedral => "D",
        }
    }
    pub fn target(self, atoms: &[usize], value: f64) -> ConstraintTarget {
        match self {
            Self::Bond => ConstraintTarget::bond_angstrom([atoms[0], atoms[1]], value),
            Self::Angle => ConstraintTarget::angle_degree([atoms[0], atoms[1], atoms[2]], value),
            Self::Dihedral => {
                ConstraintTarget::dihedral_degree([atoms[0], atoms[1], atoms[2], atoms[3]], value)
            }
        }
    }
    pub fn validate_value(self, value: f64) -> Result<(), String> {
        if !value.is_finite()
            || (self == Self::Bond && value <= 0.)
            || (self == Self::Angle && (value <= 0. || value >= 180.))
        {
            return Err(format!(
                "Invalid {} value: {value} {}.",
                self.label().to_lowercase(),
                self.unit()
            ));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FixedCoordinate {
    pub coordinate: Coordinate,
    pub atoms: String,
    pub value: f64,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Options {
    pub driver: Driver,
    pub mode: Mode,
    pub coordinate: Coordinate,
    pub atoms: String,
    pub start: f64,
    pub end: f64,
    pub points: usize,
    pub move_fragment: bool,
    pub bond_motion: BondMotion,
    pub relaxation_coordinates: RelaxationCoordinates,
    pub constraints: Vec<FixedCoordinate>,
    pub max_cycles: usize,
    pub tolerances: crate::qchem_interfaces::method::tolerances::Criteria,
    pub stop_on_failure: bool,
    #[serde(skip)]
    pub picking: bool,
    #[serde(skip)]
    pub picked: Vec<usize>,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            driver: Driver::Beavyr,
            mode: Mode::Relaxed,
            coordinate: Coordinate::Bond,
            atoms: "1 2".into(),
            start: 1.,
            end: 2.,
            points: 11,
            move_fragment: true,
            bond_motion: Default::default(),
            relaxation_coordinates: Default::default(),
            constraints: vec![],
            max_cycles: 100,
            tolerances: crate::qchem_interfaces::method::tolerances::Criteria::defaults(QcProgram::Orca),
            stop_on_failure: true,
            picking: false,
            picked: vec![],
        }
    }
}
#[derive(Debug, Clone)]
pub struct Prepared {
    pub atoms: Vec<usize>,
    pub fixed: Vec<ConstraintTarget>,
}
pub fn atom_indices(text: &str, count: usize, symbols: &[String]) -> Result<Vec<usize>, String> {
    let mut atoms = Vec::new();
    for word in text
        .split(|c: char| c.is_whitespace() || c == ',')
        .filter(|s| !s.is_empty())
    {
        let atom = word
            .parse::<usize>()
            .ok()
            .filter(|i| *i > 0 && *i <= symbols.len())
            .ok_or_else(|| format!("Use atom numbers from 1 to {}.", symbols.len()))?
            - 1;
        if atoms.contains(&atom) {
            return Err("Choose different atoms for each coordinate.".into());
        }
        if crate::molecule::is_dummy(&symbols[atom]) {
            return Err(
                "A scan coordinate cannot refer to a dummy atom. Choose real atoms.".into(),
            );
        }
        atoms.push(atom);
    }
    if atoms.len() != count {
        return Err(format!("This coordinate requires {count} atom numbers."));
    }
    Ok(atoms)
}
impl Options {
    pub fn values(&self) -> Vec<f64> {
        (0..self.points)
            .map(|i| self.start + (self.end - self.start) * i as f64 / (self.points - 1) as f64)
            .collect()
    }
    pub fn prepare(&self, symbols: &[String]) -> Result<Prepared, String> {
        self.tolerances.validate()?;
        if !(2..=10000).contains(&self.points) {
            return Err("Choose 2–10000 scan points, including both endpoints.".into());
        }
        if self.max_cycles == 0 {
            return Err("Use at least one optimization step per point.".into());
        }
        self.coordinate.validate_value(self.start)?;
        self.coordinate.validate_value(self.end)?;
        if self.start == self.end {
            return Err("Scan start and end must differ.".into());
        }
        let atoms = atom_indices(&self.atoms, self.coordinate.count(), symbols)?;
        if self.mode == Mode::Rigid && !self.constraints.is_empty() {
            return Err("Additional fixed coordinates require a relaxed scan.".into());
        }
        let mut keys = vec![self.coordinate.target(&atoms, self.start).coordinate];
        let mut fixed = vec![];
        for c in &self.constraints {
            c.coordinate.validate_value(c.value)?;
            let indices = atom_indices(&c.atoms, c.coordinate.count(), symbols)?;
            let target = c.coordinate.target(&indices, c.value);
            if keys.iter().any(|k| same_coordinate(*k, target.coordinate)) {
                return Err(
                    "A coordinate cannot be both scanned and fixed, or constrained twice.".into(),
                );
            }
            keys.push(target.coordinate);
            fixed.push(target);
        }
        Ok(Prepared { atoms, fixed })
    }
    pub fn remap_dummies(&self, symbols: &[String]) -> Result<Self, String> {
        self.prepare(symbols)?;
        let remap = |text: &str, coordinate: Coordinate| -> Result<String, String> {
            Ok(atom_indices(text, coordinate.count(), symbols)?
                .iter()
                .map(|&i| {
                    symbols[..=i]
                        .iter()
                        .filter(|s| !crate::molecule::is_dummy(s))
                        .count()
                        .to_string()
                })
                .collect::<Vec<_>>()
                .join(" "))
        };
        let mut result = self.clone();
        result.atoms = remap(&self.atoms, self.coordinate)?;
        for fixed in &mut result.constraints {
            fixed.atoms = remap(&fixed.atoms, fixed.coordinate)?;
        }
        result.picking = false;
        result.picked.clear();
        Ok(result)
    }
}
fn same_coordinate(a: ConstraintCoordinate, b: ConstraintCoordinate) -> bool {
    match (a, b) {
        (ConstraintCoordinate::Bond(a), ConstraintCoordinate::Bond(b)) => {
            a == b || a == [b[1], b[0]]
        }
        (ConstraintCoordinate::Angle(a), ConstraintCoordinate::Angle(b)) => {
            a == b || a == [b[2], b[1], b[0]]
        }
        (ConstraintCoordinate::Dihedral(a), ConstraintCoordinate::Dihedral(b)) => {
            a == b || a == [b[3], b[2], b[1], b[0]]
        }
        _ => false,
    }
}
pub fn method_problem(
    program: QcProgram,
    method: Option<QcMethod>,
    mode: Mode,
) -> Option<&'static str> {
    if program == QcProgram::StdaXtb {
        return Some("sTDA-xTB calculates spectra only.");
    }
    if mode == Mode::Relaxed
        && method.is_some_and(|m| {
            matches!(
                m,
                QcMethod::CcsdT
                    | QcMethod::DlpnoCcsdT
                    | QcMethod::CcsdTF12
                    | QcMethod::Mp2F12
                    | QcMethod::Casscf
                    | QcMethod::Nevpt2
            )
        })
    {
        return Some("This method is currently available for rigid energy scans only in Beavyr. Choose a gradient-capable method for relaxed scans.");
    }
    None
}

/// Methods supported by Behemoth’s native scan frontend.
pub fn driver_problem(
    program: QcProgram,
    method: Option<QcMethod>,
    options: &Options,
) -> Option<&'static str> {
    if program != QcProgram::Behemoth || options.driver != Driver::Behemoth {
        return None;
    }
    if !matches!(
        method,
        Some(QcMethod::Hf | QcMethod::Dft | QcMethod::Mp2 | QcMethod::RiMp2)
    ) {
        return Some("Behemoth's built-in scan frontend supports HF, DFT, MP2 and RI-MP2. Choose the Beavyr driver for other methods.");
    }
    None
}
