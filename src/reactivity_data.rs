//! Finite-difference and frontier-orbital conceptual-DFT descriptors.
//! eta=(I-A)/2, S=1/(2 eta), mu=-(I+A)/2, omega=mu²/(2 eta).
pub use crate::orbitals::overlap;
use crate::{
    orbitals::{
        density::{self, DensityMatrices, Quantity},
        grid::{self, ScalarField},
        molden::{OrbitalData, Spin},
    },
    trajectory::TrajectoryFrame,
};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
pub const HARTREE_EV: f64 = 27.211386245988;
#[derive(Clone, Serialize, Deserialize)]
pub enum Wave {
    File(PathBuf),
    Embedded(String),
}
impl Wave {
    pub fn load(&self) -> Result<OrbitalData, String> {
        let text = match self {
            Self::File(p) => {
                std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()))?
            }
            Self::Embedded(t) => t.clone(),
        };
        crate::orbitals::molden::parse_molden(&text).map_err(|e| e.to_string())
    }
}
#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct Frontier {
    pub spin: Spin,
    pub index: usize,
}
#[derive(Clone, Serialize, Deserialize)]
pub enum Sources {
    Finite {
        reference: Wave,
        added: Wave,
        removed: Wave,
    },
    Frontier {
        reference: Wave,
        homo: Frontier,
        lumo: Frontier,
    },
}
#[derive(Clone, Serialize, Deserialize)]
pub struct Global {
    pub ionization: f64,
    pub affinity: f64,
    pub mu: f64,
    pub eta: f64,
    pub softness: Option<f64>,
    pub electrophilicity: Option<f64>,
}
impl Global {
    pub fn new(i: f64, a: f64) -> Result<Self, String> {
        if !i.is_finite() || !a.is_finite() {
            return Err("Nonfinite energies for reactivity analysis.".into());
        }
        let eta = (i - a) / 2.;
        let mu = -(i + a) / 2.;
        let valid = eta > 1e-8;
        Ok(Self {
            ionization: i,
            affinity: a,
            mu,
            eta,
            softness: valid.then(|| 1. / (2. * eta)),
            electrophilicity: valid.then(|| mu * mu / (2. * eta)),
        })
    }
}
#[derive(Clone, Serialize, Deserialize)]
pub struct Frame {
    pub geometry: TrajectoryFrame,
    pub plus: Vec<f64>,
    pub minus: Vec<f64>,
    pub global: Global,
    pub sources: Sources,
}
#[derive(Clone, Serialize, Deserialize)]
pub struct Results {
    pub version: u32,
    pub method: String,
    pub frames: Vec<Frame>,
}
#[derive(Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub enum Descriptor {
    Plus,
    Minus,
    Radical,
    Dual,
    SoftPlus,
    SoftMinus,
    PhilPlus,
    PhilMinus,
    Hardness,
    Softness,
    ChemicalPotential,
    Electrophilicity,
}
impl Descriptor {
    pub const ALL: [Self; 12] = [
        Self::Plus,
        Self::Minus,
        Self::Radical,
        Self::Dual,
        Self::SoftPlus,
        Self::SoftMinus,
        Self::PhilPlus,
        Self::PhilMinus,
        Self::Hardness,
        Self::Softness,
        Self::ChemicalPotential,
        Self::Electrophilicity,
    ];
    pub fn label(self) -> &'static str {
        match self {
            Self::Plus => "Fukui f+",
            Self::Minus => "Fukui f−",
            Self::Radical => "Fukui f0",
            Self::Dual => "Dual descriptor",
            Self::SoftPlus => "Local softness s+",
            Self::SoftMinus => "Local softness s−",
            Self::PhilPlus => "Local philicity ω+",
            Self::PhilMinus => "Local philicity ω−",
            Self::Hardness => "Global hardness",
            Self::Softness => "Global softness",
            Self::ChemicalPotential => "Chemical potential",
            Self::Electrophilicity => "Electrophilicity",
        }
    }
    pub fn unit(self) -> &'static str {
        match self {
            Self::SoftPlus | Self::SoftMinus | Self::Softness => "eV⁻¹",
            Self::PhilPlus
            | Self::PhilMinus
            | Self::Hardness
            | Self::ChemicalPotential
            | Self::Electrophilicity => "eV",
            _ => "dimensionless",
        }
    }
    pub fn local(self) -> bool {
        matches!(
            self,
            Self::Plus
                | Self::Minus
                | Self::Radical
                | Self::Dual
                | Self::SoftPlus
                | Self::SoftMinus
                | Self::PhilPlus
                | Self::PhilMinus
        )
    }
    pub fn weights(self, g: &Global) -> Option<(f64, f64)> {
        Some(match self {
            Self::Plus => (1., 0.),
            Self::Minus => (0., 1.),
            Self::Radical => (0.5, 0.5),
            Self::Dual => (1., -1.),
            Self::SoftPlus => (g.softness?, 0.),
            Self::SoftMinus => (0., g.softness?),
            Self::PhilPlus => (g.electrophilicity?, 0.),
            Self::PhilMinus => (0., g.electrophilicity?),
            _ => return None,
        })
    }
    pub fn atoms(self, frame: &Frame) -> Option<Vec<f64>> {
        let (a, b) = self.weights(&frame.global)?;
        Some(
            frame
                .plus
                .iter()
                .zip(&frame.minus)
                .map(|(p, m)| a * p + b * m)
                .collect(),
        )
    }
    pub fn value(self, frame: &Frame, group: &[usize]) -> Option<f64> {
        match self {
            Self::Hardness => Some(frame.global.eta),
            Self::Softness => frame.global.softness,
            Self::ChemicalPotential => Some(frame.global.mu),
            Self::Electrophilicity => frame.global.electrophilicity,
            _ => {
                let (a, b) = self.weights(&frame.global)?;
                group
                    .iter()
                    .map(|&i| Some(a * frame.plus.get(i)? + b * frame.minus.get(i)?))
                    .sum()
            }
        }
    }
}
fn electrons(d: &OrbitalData) -> f64 {
    d.alpha.iter().chain(&d.beta).map(|o| o.occupation).sum()
}
fn same_basis(a: &OrbitalData, b: &OrbitalData) -> bool {
    a.layout.n_basis == b.layout.n_basis
        && a.shells.len() == b.shells.len()
        && a.shells.iter().zip(&b.shells).all(|(x, y)| {
            x.center == y.center
                && x.l == y.l
                && x.spherical == y.spherical
                && x.exponents == y.exponents
                && x.raw_coefficients == y.raw_coefficients
        })
}
pub fn finite(
    geometry: TrajectoryFrame,
    sources: Sources,
    energies: [f64; 3],
) -> Result<Frame, String> {
    let Sources::Finite {
        reference,
        added,
        removed,
    } = &sources
    else {
        return Err("Expected three finite-difference states".into());
    };
    let n = reference.load()?;
    let add = added.load()?;
    let rem = removed.load()?;
    for d in [&n, &add, &rem] {
        crate::trajectory_wavefunctions::job::check_geometry(d, &geometry)?;
    }
    if !same_basis(&n, &add) || !same_basis(&n, &rem) {
        return Err("All three states must use the same Gaussian basis and atom order.".into());
    }
    if (electrons(&add) - electrons(&n) - 1.).abs() > 1e-4
        || (electrons(&n) - electrons(&rem) - 1.).abs() > 1e-4
    {
        return Err("Molden occupations do not represent N, N+1 and N−1 electrons.".into());
    }
    let overlap = overlap::Overlap::build(&n)?;
    let p = overlap.population(&n)?;
    let pa = overlap.population(&add)?;
    let pr = overlap.population(&rem)?;
    let plus: Vec<_> = pa.iter().zip(&p).map(|(a, n)| a - n).collect();
    let minus: Vec<_> = p.iter().zip(&pr).map(|(n, r)| n - r).collect();
    for f in [&plus, &minus] {
        let sum = f.iter().sum::<f64>();
        if !sum.is_finite() || (sum - 1.).abs() > 0.1 {
            return Err(format!("Condensed Fukui sum is {sum:.6}; the Molden density normalization is inconsistent with a one-electron change."));
        }
    }
    let global = Global::new(
        (energies[2] - energies[0]) * HARTREE_EV,
        (energies[0] - energies[1]) * HARTREE_EV,
    )?;
    Ok(Frame {
        geometry,
        plus,
        minus,
        global,
        sources,
    })
}
pub fn frontier(geometry: TrajectoryFrame, reference: Wave, spin: Spin) -> Result<Frame, String> {
    let d = reference.load()?;
    crate::trajectory_wavefunctions::job::check_geometry(&d, &geometry)?;
    let orbitals = d.orbitals(spin);
    if orbitals
        .iter()
        .any(|o| o.occupation > 1e-6 && (o.occupation - o.occupation.round()).abs() > 1e-5)
    {
        return Err("Fractional occupations make a single HOMO/LUMO assignment ambiguous. Use the finite-difference method or a wavefunction with integer occupations.".into());
    }
    let h = orbitals
        .iter()
        .enumerate()
        .filter(|(_, o)| o.occupation > 1e-6)
        .max_by(|(_, a), (_, b)| a.energy.total_cmp(&b.energy))
        .map(|(i, _)| i)
        .ok_or("No occupied orbital in selected spin channel")?;
    let l = orbitals
        .iter()
        .enumerate()
        .filter(|(_, o)| o.occupation <= 1e-6)
        .min_by(|(_, a), (_, b)| a.energy.total_cmp(&b.energy))
        .map(|(i, _)| i)
        .ok_or("No unoccupied orbital in selected spin channel. Export virtual orbitals too.")?;
    let overlap = overlap::Overlap::build(&d)?;
    let plus = overlap.orbital(&orbitals[l], d.atoms.len())?;
    let minus = overlap.orbital(&orbitals[h], d.atoms.len())?;
    let global = Global::new(
        -orbitals[h].energy * HARTREE_EV,
        -orbitals[l].energy * HARTREE_EV,
    )?;
    Ok(Frame {
        geometry,
        plus,
        minus,
        global,
        sources: Sources::Frontier {
            reference,
            homo: Frontier { spin, index: h },
            lumo: Frontier { spin, index: l },
        },
    })
}
pub fn field(frame: &Frame, descriptor: Descriptor, spacing: f32) -> Result<ScalarField, String> {
    let (a, b) = descriptor
        .weights(&frame.global)
        .ok_or("This descriptor has no spatial field (or the energy gap is nonpositive).")?;
    let spec = grid::plan_grid(&frame.geometry.pos, 4., spacing, 2_000_000);
    let mut values = vec![0.; spec.point_count()];
    match &frame.sources {
        Sources::Finite {
            reference,
            added,
            removed,
        } => {
            for (source, weight) in [(reference, b - a), (added, a), (removed, -b)] {
                if weight == 0. {
                    continue;
                }
                let d = source.load()?;
                crate::trajectory_wavefunctions::job::check_geometry(&d, &frame.geometry)?;
                let f = density::sample_density(
                    &d,
                    &DensityMatrices::build(&d),
                    Quantity::Density,
                    &spec,
                );
                for (v, x) in values.iter_mut().zip(f.values) {
                    *v += (weight * x as f64) as f32;
                }
            }
        }
        Sources::Frontier {
            reference,
            homo,
            lumo,
        } => {
            let d = reference.load()?;
            crate::trajectory_wavefunctions::job::check_geometry(&d, &frame.geometry)?;
            for (orbital, weight) in [(lumo, a), (homo, b)] {
                if orbital.index >= d.orbitals(orbital.spin).len() {
                    return Err("Missing saved frontier orbital".into());
                }
                if weight == 0. {
                    continue;
                }
                let f = grid::sample_orbital(&d, orbital.spin, orbital.index, &spec);
                for (v, x) in values.iter_mut().zip(f.values) {
                    *v += (weight * (x as f64).powi(2)) as f32;
                }
            }
        }
    }
    if values.iter().any(|v| !v.is_finite()) {
        return Err("Nonfinite reactivity field.".into());
    }
    Ok(ScalarField { spec, values })
}
impl Results {
    pub fn validate(&self) -> Result<(), String> {
        if self.version != 1 || self.frames.is_empty() {
            return Err("Unsupported/empty reactivity analysis".into());
        }
        for f in &self.frames {
            let n = f.geometry.atoms.len();
            if n == 0
                || n != f.geometry.pos.len()
                || n != f.plus.len()
                || n != f.minus.len()
                || f.geometry.pos.iter().any(|p| !p.is_finite())
                || f.plus.iter().chain(&f.minus).any(|v| !v.is_finite())
                || f.geometry.atoms != self.frames[0].geometry.atoms
            {
                return Err("Invalid reactivity analysis geometry or populations".into());
            }
            let expected = Global::new(f.global.ionization, f.global.affinity)?;
            if !f.global.mu.is_finite()
                || !f.global.eta.is_finite()
                || (expected.eta - f.global.eta).abs() > 1e-8
                || (expected.mu - f.global.mu).abs() > 1e-8
                || expected.softness != f.global.softness
                || expected.electrophilicity != f.global.electrophilicity
            {
                return Err("Invalid global reactivity descriptors".into());
            }
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hardness_convention_and_invalid_gap() {
        let g = Global::new(10., 2.).unwrap();
        assert_eq!(g.eta, 4.);
        assert_eq!(g.mu, -6.);
        assert_eq!(g.softness, Some(0.125));
        assert_eq!(g.electrophilicity, Some(4.5));
        assert!(Global::new(1., 2.).unwrap().softness.is_none());
    }
}
