//! Analytic Gaussian AO overlaps in the same normalization/order as the Molden viewer.
//! Mulliken populations require P*S, not the CPA diagonal used for polarizability.
use crate::orbitals::{
    basis::{cartesian_component_norm, double_factorial_odd, Shell},
    cart2sph,
    molden::{MolecularOrbital, OrbitalData},
};
fn choose(n: u8, k: u8) -> f64 {
    (0..k).fold(1., |v, i| v * (n - i) as f64 / (i + 1) as f64)
}
fn axis_overlap(i: u8, j: u8, a: f64, b: f64, alpha: f64, beta: f64) -> f64 {
    let p = alpha + beta;
    let centre = (alpha * a + beta * b) / p;
    let mut sum = 0.;
    for k in 0..=i {
        for l in 0..=j {
            let n = k + l;
            if n % 2 != 0 {
                continue;
            }
            let moment = double_factorial_odd((n / 2) as usize) / (2. * p).powi((n / 2) as i32);
            sum += choose(i, k)
                * choose(j, l)
                * (centre - a).powi((i - k) as i32)
                * (centre - b).powi((j - l) as i32)
                * moment;
        }
    }
    sum * (std::f64::consts::PI / p).sqrt() * (-alpha * beta / p * (a - b).powi(2)).exp()
}
fn terms(s: &Shell, slot: usize) -> Vec<([u8; 3], f64)> {
    let l = s.l as usize;
    let powers = cart2sph::cartesian_powers(l);
    if s.spherical {
        let row = cart2sph::molden_spherical_row(l, slot);
        let table = cart2sph::transform(l);
        powers
            .iter()
            .enumerate()
            .filter_map(|(c, &p)| {
                let v = table[row * powers.len() + c] * cartesian_component_norm(p);
                (v != 0.).then_some((p, v))
            })
            .collect()
    } else {
        let p = powers[cart2sph::molden_cartesian_column(l, slot)];
        vec![(p, cartesian_component_norm(p))]
    }
}
pub struct Overlap {
    pub matrix: Vec<f64>,
    pub owners: Vec<usize>,
    pub n: usize,
}
impl Overlap {
    pub fn build(data: &OrbitalData) -> Result<Self, String> {
        let n = data.layout.n_basis;
        if n == 0 || n > 3000 {
            return Err("Mulliken analysis supports 1–3000 basis functions.".into());
        }
        let mut basis = vec![];
        for (s, shell) in data.shells.iter().enumerate() {
            for slot in 0..shell.n_functions() {
                basis.push((s, terms(shell, slot)));
            }
        }
        if basis.len() != n {
            return Err("Invalid basis layout".into());
        }
        let mut matrix = vec![0.; n * n];
        let owners = basis.iter().map(|(s, _)| data.shells[*s].center).collect();
        for mu in 0..n {
            let (si, ta) = &basis[mu];
            let a = &data.shells[*si];
            let ca = data.positions[a.center].as_dvec3() / crate::orbitals::grid::BOHR;
            for nu in 0..=mu {
                let (sj, tb) = &basis[nu];
                let b = &data.shells[*sj];
                let cb = data.positions[b.center].as_dvec3() / crate::orbitals::grid::BOHR;
                let mut value = 0.;
                for (&alpha, &ac) in a.exponents.iter().zip(&a.coefficients) {
                    for (&beta, &bc) in b.exponents.iter().zip(&b.coefficients) {
                        for &(pa, na) in ta {
                            for &(pb, nb) in tb {
                                let integral = (0..3)
                                    .map(|k| axis_overlap(pa[k], pb[k], ca[k], cb[k], alpha, beta))
                                    .product::<f64>();
                                value += ac * bc * na * nb * integral;
                            }
                        }
                    }
                }
                if !value.is_finite() {
                    return Err("Nonfinite overlap integral".into());
                }
                matrix[mu * n + nu] = value;
                matrix[nu * n + mu] = value;
            }
        }
        Ok(Self { matrix, owners, n })
    }
    pub fn orbital(&self, o: &MolecularOrbital, natoms: usize) -> Result<Vec<f64>, String> {
        if o.coefficients.len() != self.n {
            return Err("Invalid MO size".into());
        }
        let mut p = vec![0.; natoms];
        for mu in 0..self.n {
            let sc = self.matrix[mu * self.n..(mu + 1) * self.n]
                .iter()
                .zip(&o.coefficients)
                .map(|(s, c)| s * c)
                .sum::<f64>();
            p[self.owners[mu]] += o.coefficients[mu] * sc;
        }
        let norm = p.iter().sum::<f64>();
        if !norm.is_finite() || (norm - 1.).abs() > 0.03 {
            return Err(format!("Molden MO norm is {norm:.6}, not 1. This file's coefficients/basis cannot support reliable Mulliken Fukui analysis. Check the Molden export convention or use another wavefunction export."));
        }
        Ok(p)
    }
    pub fn population(&self, data: &OrbitalData) -> Result<Vec<f64>, String> {
        let mut p = vec![0.; data.atoms.len()];
        for o in data
            .alpha
            .iter()
            .chain(&data.beta)
            .filter(|o| o.occupation > 1e-8)
        {
            let v = self.orbital(o, p.len())?;
            for (a, b) in p.iter_mut().zip(v) {
                *a += o.occupation * b;
            }
        }
        Ok(p)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn normalized_spherical_and_cartesian_diagonals() {
        for l in 0..=4 {
            for spherical in [false, true] {
                let mut s = Shell {
                    center: 0,
                    l,
                    spherical,
                    exponents: vec![0.7],
                    raw_coefficients: vec![1.],
                    coefficients: vec![],
                };
                s.normalize();
                let layout = crate::orbitals::basis::BasisLayout::build(&[s.clone()]);
                let data = OrbitalData {
                    atoms: vec!["H".into()],
                    positions: vec![bevy::prelude::Vec3::ZERO],
                    shells: vec![s],
                    layout,
                    ..Default::default()
                };
                let overlap = Overlap::build(&data).unwrap();
                for i in 0..overlap.n {
                    assert!(
                        (overlap.matrix[i * overlap.n + i] - 1.).abs() < 1e-10,
                        "l={l}, spherical={spherical}, i={i}"
                    );
                }
            }
        }
    }
    #[test]
    fn displaced_s_overlap_matches_gaussian_product() {
        let alpha = 0.7;
        let s = axis_overlap(0, 0, 0., 2., alpha, alpha)
            * axis_overlap(0, 0, 0., 0., alpha, alpha).powi(2)
            * (2. * alpha / std::f64::consts::PI).powf(1.5);
        assert!((s - (-alpha * 2.).exp()).abs() < 1e-12);
    }
}
