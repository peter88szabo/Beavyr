//! Finite-temperature two-Hessian absorption (and the shared emission transform).
use super::{
    complex_math::{inverse_fft, C},
    gaussian_correlation::{continuous_branch, GaussianTrace},
    harmonic_model::HarmonicModel,
    CM1_TO_AU, HARTREE_TO_EV,
};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
#[derive(Clone, Debug)]
pub struct SpectrumGrid {
    pub temperature: f64,
    pub spacing_cm1: f64,
    pub lorentz_fwhm_cm1: f64,
    pub gaussian_fwhm_cm1: f64,
    pub max_points: usize,
}
impl Default for SpectrumGrid {
    fn default() -> Self {
        Self {
            temperature: 298.15,
            spacing_cm1: 10.,
            lorentz_fwhm_cm1: 100.,
            gaussian_fwhm_cm1: 0.,
            max_points: 65536,
        }
    }
}
#[derive(Clone, Debug)]
pub struct HarmonicSpectrum {
    pub profile: Vec<(f64, f64)>,
    pub correlation: Vec<(f64, C)>,
    pub zero_zero_ev: f64,
    pub spacing_ev: f64,
    pub temperature: f64,
    pub integral: f64,
    pub notes: Vec<String>,
}
#[cfg(test)]
pub fn calculate(
    model: HarmonicModel,
    grid: &SpectrumGrid,
    cancel: &AtomicBool,
    progress: &AtomicUsize,
) -> Result<HarmonicSpectrum, String> {
    calculate_with_dipole(model, grid, None, cancel, progress)
}
pub fn calculate_with_dipole(
    model: HarmonicModel,
    grid: &SpectrumGrid,
    dipole: Option<&super::transition_properties::Dipole>,
    cancel: &AtomicBool,
    progress: &AtomicUsize,
) -> Result<HarmonicSpectrum, String> {
    validate_grid(grid)?;
    let zero_zero = model.zero_zero();
    let trace = GaussianTrace::new(model, grid.temperature)?;
    let dw = grid.spacing_cm1 * CM1_TO_AU;
    let total_time = std::f64::consts::TAU / dw;
    let needed = (total_time * trace.energy_scale() * 2. / std::f64::consts::PI).ceil();
    if !needed.is_finite() || needed > grid.max_points as f64 {
        return Err("The time grid is too large. Increase spectral grid spacing or maximum time-grid points.".into());
    }
    let points = (needed as usize).max(512).next_power_of_two();
    if points > grid.max_points {
        return Err("The time grid exceeds the selected maximum".into());
    }
    let dt = total_time / points as f64;
    let gamma = 0.5 * grid.lorentz_fwhm_cm1 * CM1_TO_AU;
    let sigma = grid.gaussian_fwhm_cm1 * CM1_TO_AU / (8. * 2_f64.ln()).sqrt();
    let envelope = |t: f64| (-gamma * t - 0.5 * (sigma * t).powi(2)).exp();
    if envelope(total_time * 0.5) > 1e-6 {
        return Err("The requested grid spacing is too coarse for this linewidth. Reduce spacing or increase broadening so the correlation decays before the time window ends.".into());
    }
    let strength = if let Some(d) = dipole {
        d.optical(&trace.sample(0., true)?).re
    } else {
        1.
    };
    if !strength.is_finite() || strength <= 0. {
        return Err("The supplied transition dipole and derivatives give zero intensity".into());
    }
    let mut series = vec![C::ZERO; points];
    series[0] = C::from(strength);
    let mut previous = C::ONE;
    let mut correlation = vec![(0., C::from(strength))];
    for k in 1..=points / 2 {
        if cancel.load(Ordering::Relaxed) {
            return Err("Cancelled.".into());
        }
        let time = k as f64 * dt;
        let moments = trace.sample(time, dipole.is_some())?;
        let base = continuous_branch(moments.base, previous);
        previous = base;
        let value = base
            * if let Some(d) = dipole {
                d.optical(&moments)
            } else {
                C::ONE
            };
        correlation.push((time, value));
        series[k] = value * envelope(time);
        if k < points / 2 {
            series[points - k] = series[k].conj();
        } else {
            series[k].im = 0.;
        }
        progress.store(k * 100 / (points / 2), Ordering::Relaxed);
    }
    inverse_fft(&mut series);
    let negative: f64 = series.iter().map(|v| (-v.re).max(0.)).sum();
    if negative > 1e-4 * strength {
        return Err(format!("Spectrum has significant negative numerical weight ({negative:.3e}). Increase time-grid resolution and inspect the harmonic model."));
    }
    let integral = series.iter().map(|v| v.re).sum::<f64>();
    let mut profile: Vec<_> = series
        .iter()
        .enumerate()
        .map(|(k, v)| {
            let index = if k <= points / 2 {
                k as f64
            } else {
                k as f64 - points as f64
            };
            (
                (zero_zero + index * dw) * HARTREE_TO_EV,
                v.re.max(0.) / (dw * HARTREE_TO_EV),
            )
        })
        .collect();
    profile.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(HarmonicSpectrum {
        profile,
        correlation,
        zero_zero_ev: zero_zero * HARTREE_TO_EV,
        spacing_ev: dw * HARTREE_TO_EV,
        temperature: grid.temperature,
        integral,
        notes: vec![],
    })
}
pub fn validate_grid(grid: &SpectrumGrid) -> Result<(), String> {
    if !grid.temperature.is_finite()
        || grid.temperature < 0.
        || !grid.spacing_cm1.is_finite()
        || grid.spacing_cm1 <= 0.
        || !grid.lorentz_fwhm_cm1.is_finite()
        || grid.lorentz_fwhm_cm1 < 0.
        || !grid.gaussian_fwhm_cm1.is_finite()
        || grid.gaussian_fwhm_cm1 < 0.
        || (grid.lorentz_fwhm_cm1 == 0. && grid.gaussian_fwhm_cm1 == 0.)
        || grid.max_points < 512
        || grid.max_points > 262144
    {
        return Err("Invalid temperature, linewidth or time-grid settings".into());
    }
    Ok(())
}
