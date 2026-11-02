//! Stöhr CPA (JCP 144, 151101, 2016), matching Behemoth properties/dispersion.rs.
//! h_A = sum_mu(on A) P_mu,mu; alpha_A = alpha_free h_A/Z_A.
//! This is a basis-dependent atomic descriptor, not a molecular response tensor.
use crate::{orbitals::molden::OrbitalData, trajectory::TrajectoryFrame};
use serde::{Deserialize, Serialize};
#[path = "polarizability_reference.rs"]
mod reference;
pub const BOHR3_TO_ANGSTROM3: f64 = 0.5291772109_f64 * 0.5291772109 * 0.5291772109;
#[derive(Clone, Serialize, Deserialize)]
pub struct Values {
    pub population: Vec<f64>,
    pub alpha: Vec<f64>,
    pub c6: Vec<f64>,
}
pub fn compute(data: &OrbitalData) -> Result<Values, String> {
    if data.atoms.is_empty() || data.alpha.is_empty() {
        return Err("Load a Molden wavefunction first.".into());
    }
    let mut owners = Vec::with_capacity(data.layout.n_basis);
    for shell in &data.shells {
        if shell.center >= data.atoms.len() {
            return Err("Invalid basis centre.".into());
        }
        owners.extend(std::iter::repeat_n(shell.center, shell.n_functions()));
    }
    if owners.len() != data.layout.n_basis {
        return Err("Inconsistent Molden basis layout.".into());
    }
    let mut population = vec![0.; data.atoms.len()];
    for orbital in data.alpha.iter().chain(&data.beta) {
        if !orbital.occupation.is_finite()
            || orbital.occupation < 0.
            || orbital.coefficients.len() != owners.len()
        {
            return Err("Invalid orbital occupations or coefficients.".into());
        }
        for (&atom, &c) in owners.iter().zip(&orbital.coefficients) {
            if !c.is_finite() {
                return Err("Nonfinite orbital coefficient.".into());
            }
            population[atom] += orbital.occupation * c * c;
        }
    }
    let mut alpha = vec![];
    let mut c6 = vec![];
    for (i, symbol) in data.atoms.iter().enumerate() {
        let (_, a, c, _) = reference::FREE_ATOM_TABLE
            .iter()
            .find(|v| v.0 == symbol)
            .ok_or_else(|| format!("No TS free-atom reference is available for {symbol}."))?;
        let z = crate::molecule::atomic_number(symbol).ok_or("Unknown element")? as f64;
        let ratio = population[i] / z;
        if !ratio.is_finite() || ratio < 0. {
            return Err("Invalid CPA population.".into());
        }
        alpha.push(a * ratio);
        c6.push(c * ratio * ratio);
    }
    Ok(Values {
        population,
        alpha,
        c6,
    })
}
pub fn atom_list(text: &str, n: usize) -> Result<Vec<usize>, String> {
    let mut result = std::collections::BTreeSet::new();
    for token in text
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|s| !s.is_empty())
    {
        let (start, end) = token.split_once('-').unwrap_or((token, token));
        let start: usize = start
            .parse()
            .map_err(|_| "Use atom numbers or ranges, e.g. 1,3-5.")?;
        let end: usize = end.parse().map_err(|_| "Invalid atom range.")?;
        if start == 0 || end < start || end > n {
            return Err(format!("Atom numbers must be between 1 and {n}."));
        }
        result.extend(start - 1..end);
    }
    if result.is_empty() {
        return Err("Select at least one atom.".into());
    }
    Ok(result.into_iter().collect())
}
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    Frame,
    Time,
    Distance,
    Angle,
    Dihedral,
}
impl Axis {
    pub fn label(self) -> &'static str {
        match self {
            Self::Frame => "Frame",
            Self::Time => "Time (fs)",
            Self::Distance => "Distance (Å)",
            Self::Angle => "Angle (°)",
            Self::Dihedral => "Dihedral (°)",
        }
    }
    pub fn count(self) -> usize {
        match self {
            Self::Distance => 2,
            Self::Angle => 3,
            Self::Dihedral => 4,
            _ => 0,
        }
    }
    pub fn value(
        self,
        frame: &TrajectoryFrame,
        index: usize,
        atoms: &[usize; 4],
        dt: f64,
        origin: f64,
    ) -> Option<f64> {
        if self == Self::Frame {
            return Some((index + 1) as f64);
        }
        if self == Self::Time {
            return Some(origin + index as f64 * dt);
        }
        let p = atoms
            .iter()
            .take(self.count())
            .map(|i| i.checked_sub(1).and_then(|j| frame.pos.get(j)).copied())
            .collect::<Option<Vec<_>>>()?;
        let v = match self {
            Self::Distance => p[0].distance(p[1]),
            Self::Angle => {
                let a = (p[0] - p[1]).try_normalize()?;
                let b = (p[2] - p[1]).try_normalize()?;
                a.dot(b).clamp(-1., 1.).acos().to_degrees()
            }
            Self::Dihedral => {
                let axis = (p[2] - p[1]).try_normalize()?;
                let a = p[0] - p[1];
                let b = p[3] - p[2];
                let a = (a - axis * a.dot(axis)).try_normalize()?;
                let b = (b - axis * b.dot(axis)).try_normalize()?;
                axis.cross(a).dot(b).atan2(a.dot(b)).to_degrees()
            }
            _ => return None,
        };
        v.is_finite().then_some(v as f64)
    }
}
pub fn radius(alpha_bohr3: f64, scale: f32) -> f32 {
    (3. * alpha_bohr3 * BOHR3_TO_ANGSTROM3 / (4. * std::f64::consts::PI)).cbrt() as f32 * scale
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn atom_groups_do_not_double_count() {
        assert_eq!(atom_list("1,3-5,3", 5).unwrap(), vec![0, 2, 3, 4]);
        assert!(atom_list("0", 5).is_err());
        assert!(atom_list("1-99999999", 5).is_err());
    }
    #[test]
    fn cpa_matches_behemoth_diagonal_and_spin_sum() {
        let text="[Molden Format]\n[Atoms] Angs\nH 1 1 0 0 0\n[GTO]\n1 0\ns 1 1.0\n1.0 1.0\n\n[MO]\nSym= 1\nEne= -0.5\nSpin= Alpha\nOccup= 0.5\n1 2.0\nSym= 1\nEne= -0.5\nSpin= Beta\nOccup= 0.5\n1 2.0\n";
        let data = crate::orbitals::molden::parse_molden(text).unwrap();
        let v = compute(&data).unwrap();
        assert_eq!(v.population, vec![4.]);
        assert_eq!(v.alpha, vec![18.]);
        assert_eq!(v.c6, vec![104.]);
        assert!((radius(8., 1.) / radius(1., 1.) - 2.).abs() < 1e-6);
    }
}
