//! Length-gauge transition dipoles, referenced to the ground Hessian geometry.
//! The Cartesian-to-normal projection follows Behemoth's Herzberg–Teller module.
use super::{
    complex_math::C, gaussian_correlation::Moments, harmonic_model::HarmonicModel,
    input::StateModes,
};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PropertySource {
    Manual,
    File,
    PyScf,
    Psi4,
}
#[derive(Clone, Debug)]
pub struct PropertySettings {
    pub ht: bool,
    pub source: PropertySource,
    pub path: String,
    pub step: f64,
}
impl Default for PropertySettings {
    fn default() -> Self {
        Self {
            ht: false,
            source: PropertySource::Manual,
            path: String::new(),
            step: 0.001,
        }
    }
}
#[derive(Clone, Debug)]
pub struct Dipole {
    pub value: [f64; 3],
    pub derivative: Vec<[f64; 3]>,
}
impl Dipole {
    pub fn constant(value: [f64; 3], n: usize) -> Self {
        Self {
            value,
            derivative: vec![[0.; 3]; n],
        }
    }
    pub fn reversed(&self, model: &HarmonicModel) -> Self {
        let n = model.n();
        Self {
            value: std::array::from_fn(|a| {
                self.value[a]
                    + (0..n)
                        .map(|i| self.derivative[i][a] * model.displacement[i])
                        .sum::<f64>()
            }),
            derivative: (0..n)
                .map(|j| {
                    std::array::from_fn(|a| {
                        (0..n)
                            .map(|i| self.derivative[i][a] * model.j[(i, j)])
                            .sum()
                    })
                })
                .collect(),
        }
    }
    pub fn optical(&self, m: &Moments) -> C {
        let p = self.prepare(m);
        (0..3).fold(C::ZERO, |s, a| s + p.moment(&[], a, a))
    }
    /// Precontract the linear dipoles once per time step; Raman channels then
    /// require only Gaussian moments of at most four scalar observables.
    pub fn prepare(&self, m: &Moments) -> DipoleMoments {
        let n = self.derivative.len();
        let len = n + 6;
        let mut means = vec![C::ZERO; len];
        means[..n].copy_from_slice(&m.mean);
        let mut cov = vec![C::ZERO; len * len];
        for i in 0..n {
            for j in 0..n {
                cov[i * len + j] = m.xx[i * n + j];
            }
        }
        for a in 0..3 {
            let mean = C::from(self.value[a])
                + (0..n).fold(C::ZERO, |s, i| s + m.mean[i] * self.derivative[i][a]);
            means[n + a] = mean;
            means[n + 3 + a] = mean;
            for i in 0..n {
                let ca = (0..n).fold(C::ZERO, |s, j| s + m.xx[i * n + j] * self.derivative[j][a]);
                let cb = (0..n).fold(C::ZERO, |s, j| s + m.xy[i * n + j] * self.derivative[j][a]);
                cov[i * len + n + a] = ca;
                cov[(n + a) * len + i] = ca;
                cov[i * len + n + 3 + a] = cb;
                cov[(n + 3 + a) * len + i] = cb;
            }
            for b in 0..3 {
                let cross = (0..n).fold(C::ZERO, |s, i| {
                    s + (0..n).fold(C::ZERO, |t, j| {
                        t + m.xy[i * n + j] * (self.derivative[i][a] * self.derivative[j][b])
                    })
                });
                cov[(n + a) * len + n + 3 + b] = cross;
                cov[(n + 3 + b) * len + n + a] = cross;
            }
        }
        DipoleMoments { n, means, cov }
    }
}
pub struct DipoleMoments {
    n: usize,
    means: Vec<C>,
    cov: Vec<C>,
}
impl DipoleMoments {
    pub fn moment(&self, indices: &[usize], a: usize, b: usize) -> C {
        let mut indices = indices.to_vec();
        indices.extend([self.n + a, self.n + 3 + b]);
        wick(&indices, &self.means, &self.cov, self.n + 6)
    }
}

fn wick(indices: &[usize], means: &[C], cov: &[C], n: usize) -> C {
    if indices.is_empty() {
        return C::ONE;
    }
    let first = indices[0];
    let mut value = means[first] * wick(&indices[1..], means, cov, n);
    for k in 1..indices.len() {
        let remaining: Vec<_> = indices[1..]
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != k - 1)
            .map(|(_, v)| *v)
            .collect();
        value += cov[first * n + indices[k]] * wick(&remaining, means, cov, n);
    }
    value
}
pub fn parse(text: &str, reference: &StateModes, ht: bool) -> Result<Dipole, String> {
    let mut rows = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));
    if rows.next() != Some("BEAVYR_TRANSITION_V1") {
        return Err("Expected BEAVYR_TRANSITION_V1 property file; use Save template for its format and units.".into());
    }
    let number = |s: &str| {
        s.parse::<f64>()
            .ok()
            .filter(|x| x.is_finite())
            .ok_or_else(|| format!("Invalid transition property: {s}"))
    };
    let natoms = rows
        .next()
        .ok_or("Missing atom count")?
        .parse::<usize>()
        .map_err(|_| "Invalid atom count")?;
    if natoms != reference.atoms.len() {
        return Err("Transition property atom count differs from the Hessian".into());
    }
    for i in 0..natoms {
        let r = rows
            .next()
            .ok_or("Missing property geometry")?
            .split_whitespace()
            .collect::<Vec<_>>();
        if r.len() != 4 || !r[0].eq_ignore_ascii_case(&reference.atoms[i]) {
            return Err("Transition property atom order differs from the Hessian".into());
        }
        for a in 0..3 {
            if (number(r[a + 1])? - reference.coords[3 * i + a]).abs() > 1e-6 {
                return Err("Transition properties must use exactly the ground Hessian geometry and Cartesian frame (bohr).".into());
            }
        }
    }
    let vector = |line: &str| -> Result<[f64; 3], String> {
        let r = line
            .split_whitespace()
            .map(number)
            .collect::<Result<Vec<_>, _>>()?;
        r.try_into()
            .map_err(|_| "Expected three dipole components".into())
    };
    let value = vector(rows.next().ok_or("Missing transition dipole")?)?;
    let mut result = Dipole::constant(value, reference.selected.len());
    if ht {
        let cart = (0..3 * natoms)
            .map(|_| {
                rows.next()
                    .ok_or_else(|| {
                        "Missing Cartesian transition-dipole derivatives (3N rows)".to_string()
                    })
                    .and_then(vector)
            })
            .collect::<Result<Vec<_>, _>>()?;
        for (k, &mode) in reference.selected.iter().enumerate() {
            result.derivative[k] = std::array::from_fn(|a| {
                (0..3 * natoms)
                    .map(|i| reference.modes[(i, mode)] * cart[i][a])
                    .sum()
            });
        }
    }
    Ok(result)
}
pub fn template(reference: &StateModes) -> String {
    use std::fmt::Write;
    let mut s=String::from("# Length-gauge dipole in e*bohr; Cartesian derivative in e; coordinates in bohr.\n# Phase-align the same electronic state across all displacements.\nBEAVYR_TRANSITION_V1\n");
    let _ = writeln!(s, "{}", reference.atoms.len());
    for (i, a) in reference.atoms.iter().enumerate() {
        let _ = writeln!(
            s,
            "{} {:.12} {:.12} {:.12}",
            a,
            reference.coords[3 * i],
            reference.coords[3 * i + 1],
            reference.coords[3 * i + 2]
        );
    }
    s.push_str("# Replace the next row with transition dipole x y z\nNaN NaN NaN\n# Replace each row below with d(mu_x,mu_y,mu_z)/dx_i, atom-major x,y,z order\n");
    for _ in &reference.coords {
        s.push_str("NaN NaN NaN\n");
    }
    s
}
