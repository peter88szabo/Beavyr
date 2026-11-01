//! Harmonic thermal trace as two Gaussian kernels (Baiardi et al. 2013;
//! de Souza et al. 2018). Complex times are beta-it and it, in atomic units.
//! Sum/difference coordinates reduce the 2N Gaussian to two N matrices.
use super::{
    complex_math::{Lu, C},
    harmonic_model::HarmonicModel,
};
pub const KB_HARTREE: f64 = 3.166811563e-6;
pub struct GaussianTrace {
    pub model: HarmonicModel,
    transform: Vec<f64>,
    log_jacobian: f64,
    temperature: f64,
}
#[derive(Clone)]
pub struct Moments {
    pub base: C,
    pub mean: Vec<C>,
    pub xx: Vec<C>,
    pub xy: Vec<C>,
}
impl GaussianTrace {
    pub fn new(model: HarmonicModel, temperature: f64) -> Result<Self, String> {
        if !temperature.is_finite() || temperature < 0. {
            return Err("Temperature must be finite and nonnegative".into());
        }
        let (transform, log_jacobian) = model.inverse_rotation()?;
        Ok(Self {
            model,
            transform,
            log_jacobian,
            temperature,
        })
    }
    pub fn occupation(&self, i: usize) -> f64 {
        if self.temperature == 0. {
            0.
        } else {
            1. / (self.model.initial[i] / (KB_HARTREE * self.temperature)).exp_m1()
        }
    }
    fn transformed(&self, diagonal: &[C]) -> Vec<C> {
        let n = self.model.n();
        let mut out = vec![C::ZERO; n * n];
        for k in 0..n {
            for i in 0..n {
                let v = diagonal[k] * self.transform[k * n + i];
                for j in 0..n {
                    out[i * n + j] += v * self.transform[k * n + j];
                }
            }
        }
        out
    }
    pub fn sample(&self, time: f64, moments: bool) -> Result<Moments, String> {
        let n = self.model.n();
        if time == 0. {
            let mut cov = vec![C::ZERO; n * n];
            for i in 0..n {
                cov[i * n + i] = C::from((self.occupation(i) + 0.5) / self.model.initial[i]);
            }
            return Ok(Moments {
                base: C::ONE,
                mean: vec![C::ZERO; n],
                xx: cov.clone(),
                xy: cov,
            });
        }
        // A tiny imaginary-time displacement regularises kernel caustics; move it
        // from the initial kernel to the final one so the total thermal time is beta.
        let eps = 1e-8
            / self
                .model
                .initial
                .iter()
                .chain(&self.model.final_)
                .copied()
                .fold(0_f64, f64::max);
        let mut s = vec![C::ZERO; n];
        let mut v = s.clone();
        let mut log_pref = C::from(2. * self.log_jacobian);
        for i in 0..n {
            let w = self.model.final_[i];
            let q = C::new(-eps * w, -time * w).exp();
            s[i] = ((C::ONE - q) / (C::ONE + q)) * w;
            v[i] = ((C::ONE + q) / (C::ONE - q)) * w;
            log_pref += ((q * (2. * w)) / (C::ONE - q * q)).ln();
        }
        let se = self.transformed(&s);
        let mut u = se.clone();
        let mut v = self.transformed(&v);
        for i in 0..n {
            let w = self.model.initial[i];
            let thermal = if self.temperature == 0. {
                0.
            } else {
                (-w / (KB_HARTREE * self.temperature)).exp()
            };
            let phase = C::new(eps * w, time * w).exp();
            let q = phase * thermal;
            u[i * n + i] += ((C::ONE - q) / (C::ONE + q)) * w;
            v[i * n + i] += ((C::ONE + q) / (C::ONE - q)) * w;
            log_pref += ((phase * (2. * w * (1. - thermal).powi(2))) / (C::ONE - q * q)).ln();
        }
        let u = Lu::new(u, n)?;
        let v = Lu::new(v, n)?;
        let b: Vec<C> = (0..n)
            .map(|i| {
                (0..n).fold(C::ZERO, |a, j| {
                    a + se[i * n + j] * self.model.displacement[j]
                })
            })
            .collect();
        let mean = u.solve(&b);
        let exponent = (0..n).fold(C::ZERO, |a, i| {
            a + b[i] * mean[i] - b[i] * self.model.displacement[i]
        });
        let base = ((log_pref - u.logdet - v.logdet) * 0.5
            + exponent
            + C::new(0., self.model.zpe_difference() * time))
        .exp();
        if !base.finite() {
            return Err(
                "Non-finite harmonic correlation; check the Hessians and mode mapping".into(),
            );
        }
        let (xx, xy) = if moments {
            let a = u.inverse();
            let b = v.inverse();
            (
                a.iter().zip(&b).map(|(&a, &b)| (a + b) * 0.5).collect(),
                a.iter().zip(&b).map(|(&a, &b)| (a - b) * 0.5).collect(),
            )
        } else {
            (vec![], vec![])
        };
        Ok(Moments { base, mean, xx, xy })
    }
    /// Conservative sampling scale for the phase of the vibrational correlation.
    pub fn energy_scale(&self) -> f64 {
        let n = self.model.n();
        let shift: Vec<f64> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| self.transform[i * n + j] * self.model.displacement[j])
                    .sum()
            })
            .collect();
        let lambda = 0.5
            * self
                .model
                .final_
                .iter()
                .zip(shift)
                .map(|(w, d)| (w * d).powi(2))
                .sum::<f64>();
        let mismatch = self
            .model
            .initial
            .iter()
            .zip(&self.model.final_)
            .enumerate()
            .map(|(i, (a, b))| (a - b).abs() * (0.5 + self.occupation(i)))
            .sum::<f64>();
        lambda
            + mismatch
            + 4. * (2. * lambda * KB_HARTREE * self.temperature).sqrt()
            + 4. * self
                .model
                .initial
                .iter()
                .chain(&self.model.final_)
                .copied()
                .fold(0_f64, f64::max)
    }
}
/// Resolve the square-root branch continuously. The time grid bounds the phase
/// increment; the zero-time reference is exactly one.
pub fn continuous_branch(value: C, previous: C) -> C {
    if (value - previous).abs2() <= (-value - previous).abs2() {
        value
    } else {
        -value
    }
}
