//! Resonant Kramers–Heisenberg tensor, one intermediate electronic surface.
//! Ground vibrational state at 0 K; final fundamentals, overtones and combinations.
use super::{
    complex_math::C,
    gaussian_correlation::{continuous_branch, GaussianTrace},
    harmonic_model::HarmonicModel,
    thermal_spectrum::{HarmonicSpectrum, SpectrumGrid},
    CM1_TO_AU, HARTREE_TO_EV,
};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
#[derive(Clone, Debug)]
pub struct RamanSettings {
    pub lasers: String,
    pub order: usize,
    pub line_fwhm: f64,
    pub dipole: [f64; 3],
    pub photon_factor: bool,
}
impl Default for RamanSettings {
    fn default() -> Self {
        Self {
            lasers: "532".into(),
            order: 1,
            line_fwhm: 10.,
            dipole: [0.; 3],
            photon_factor: false,
        }
    }
}
impl RamanSettings {
    pub fn wavelengths(&self) -> Result<Vec<f64>, String> {
        let values = self
            .lasers
            .split(|c: char| c == ',' || c == ';' || c.is_whitespace())
            .filter(|s| !s.is_empty())
            .map(|s| {
                s.parse::<f64>()
                    .map_err(|_| "Invalid laser wavelength".to_string())
            })
            .collect::<Result<Vec<_>, _>>()?;
        if values.is_empty()
            || values.len() > 32
            || values.iter().any(|v| !v.is_finite() || *v <= 0.)
        {
            return Err("Enter 1–32 positive laser wavelengths in nm".into());
        }
        Ok(values)
    }
}
#[derive(Clone, Debug)]
pub struct RamanLine {
    pub laser_nm: f64,
    pub modes: Vec<usize>,
    pub shift_cm1: f64,
    pub strength: f64,
}
struct Channel {
    modes: Vec<usize>,
    energy: f64,
}
#[cfg(test)]
pub fn calculate(
    model: HarmonicModel,
    grid: &SpectrumGrid,
    settings: &RamanSettings,
    cancel: &AtomicBool,
    progress: &AtomicUsize,
) -> Result<(Vec<(String, HarmonicSpectrum)>, Vec<RamanLine>), String> {
    let dipole = super::transition_properties::Dipole::constant(settings.dipole, model.n());
    calculate_with_dipole(model, grid, settings, &dipole, cancel, progress)
}
pub fn calculate_with_dipole(
    model: HarmonicModel,
    grid: &SpectrumGrid,
    settings: &RamanSettings,
    dipole: &super::transition_properties::Dipole,
    cancel: &AtomicBool,
    progress: &AtomicUsize,
) -> Result<(Vec<(String, HarmonicSpectrum)>, Vec<RamanLine>), String> {
    super::thermal_spectrum::validate_grid(grid)?;
    if grid.temperature != 0. {
        return Err("Resonance Raman currently requires 0 K (vibrational ground state), for both Beavyr and native ORCA ESD.".into());
    }
    if !(1..=2).contains(&settings.order)
        || !settings.line_fwhm.is_finite()
        || settings.line_fwhm <= 0.
        || settings.dipole.iter().any(|v| !v.is_finite())
    {
        return Err("Invalid Raman order, linewidth or transition dipole".into());
    }
    if dipole
        .value
        .iter()
        .chain(dipole.derivative.iter().flatten())
        .all(|v| *v == 0.)
    {
        return Err(
            "Supply the transition-dipole vector in atomic units for resonance Raman.".into(),
        );
    }
    let lasers = settings.wavelengths()?;
    let laser_au: Vec<_> = lasers.iter().map(|nm| 1e7 / nm * CM1_TO_AU).collect();
    let n = model.n();
    let e00 = model.zero_zero();
    let mut channels = vec![];
    for i in 0..n {
        channels.push(Channel {
            modes: vec![i],
            energy: model.initial[i],
        });
    }
    if settings.order == 2 {
        for i in 0..n {
            for j in i..n {
                channels.push(Channel {
                    modes: vec![i, j],
                    energy: model.initial[i] + model.initial[j],
                });
            }
        }
    }
    if channels.len() * lasers.len() > 200000 {
        return Err("Too many Raman channels; use fundamentals or fewer lasers.".into());
    }
    let trace = GaussianTrace::new(model, 0.)?;
    let gamma = grid.lorentz_fwhm_cm1 * CM1_TO_AU / 2.;
    let sigma = grid.gaussian_fwhm_cm1 * CM1_TO_AU / (8. * 2_f64.ln()).sqrt();
    let tmax = if gamma > 0. {
        24. / gamma
    } else {
        f64::INFINITY
    }
    .min(if sigma > 0. {
        (48_f64).sqrt() / sigma
    } else {
        f64::INFINITY
    });
    let detuning = laser_au
        .iter()
        .map(|e| (e00 - e).abs())
        .fold(0_f64, f64::max);
    // Oversample the fastest phase; the quadrature is checked on a nested coarse grid.
    let steps = ((tmax * (trace.energy_scale() + detuning) * 8. / std::f64::consts::PI).ceil()
        as usize)
        .max(512)
        .next_multiple_of(2);
    if steps > grid.max_points {
        return Err(format!(
            "Raman needs {steps} integration points; increase the maximum or lifetime broadening."
        ));
    }
    let dt = tmax / steps as f64;
    let mut amplitudes = vec![C::ZERO; channels.len() * lasers.len() * 9];
    let mut coarse = amplitudes.clone();
    let mut previous = C::ONE;
    for k in 0..=steps {
        if cancel.load(Ordering::Relaxed) {
            return Err("Cancelled.".into());
        }
        let t = k as f64 * dt;
        let m = trace.sample(t, true)?;
        let base = continuous_branch(m.base, previous);
        previous = base;
        let damping = (-gamma * t - 0.5 * (sigma * t).powi(2)).exp();
        let weighted = dipole.prepare(&m);
        for (ch, channel) in channels.iter().enumerate() {
            let i = channel.modes[0];
            let wi = trace.model.initial[i];
            for a in 0..3 {
                for b in 0..3 {
                    let factor = if channel.modes.len() == 1 {
                        weighted.moment(&[i], a, b) * (2. * wi).sqrt()
                    } else {
                        let j = channel.modes[1];
                        let moment = weighted.moment(&[i, j], a, b);
                        if i == j {
                            (moment * (2. * wi) - weighted.moment(&[], a, b)) / 2_f64.sqrt()
                        } else {
                            moment * (2. * (wi * trace.model.initial[j]).sqrt())
                        }
                    };
                    for (l, energy) in laser_au.iter().enumerate() {
                        let value =
                            C::I * base * factor * C::new(0., -(e00 - energy) * t).exp() * damping;
                        let at = (l * channels.len() + ch) * 9 + a * 3 + b;
                        amplitudes[at] +=
                            value * (dt * if k == 0 || k == steps { 0.5 } else { 1. });
                        if k % 2 == 0 {
                            coarse[at] +=
                                value * (2. * dt * if k == 0 || k == steps { 0.5 } else { 1. });
                        }
                    }
                }
            }
        }
        progress.store(k * 100 / steps, Ordering::Relaxed);
    }
    let norm = amplitudes.iter().map(|z| z.abs2()).sum::<f64>();
    let error = amplitudes
        .iter()
        .zip(&coarse)
        .map(|(a, b)| (*a - *b).abs2())
        .sum::<f64>();
    if error > 1e-5 * norm.max(1e-20) {
        return Err("Raman time integration did not converge. Increase lifetime broadening or inspect the model.".into());
    }

    let mut spectra = vec![];
    let mut lines = vec![];
    for (l, &laser) in lasers.iter().enumerate() {
        let mut sticks = vec![];
        for (ch, channel) in channels.iter().enumerate() {
            let scattered = laser_au[l] - channel.energy;
            if scattered <= 0. {
                continue;
            }
            let strength = scattering_factor(
                &amplitudes[(l * channels.len() + ch) * 9..(l * channels.len() + ch + 1) * 9],
            ) * if settings.photon_factor {
                scattered.powi(4)
            } else {
                1.
            };
            lines.push(RamanLine {
                laser_nm: laser,
                modes: channel
                    .modes
                    .iter()
                    .map(|i| trace.model.initial_labels[*i] + 1)
                    .collect(),
                shift_cm1: channel.energy / CM1_TO_AU,
                strength,
            });
            sticks.push((channel.energy / CM1_TO_AU, strength));
        }
        let end = sticks.iter().map(|x| x.0).fold(0_f64, f64::max) + 10. * settings.line_fwhm;
        let count = (end / grid.spacing_cm1).ceil() as usize + 1;
        if count > 262144 {
            return Err("Raman output grid is too large; increase spacing.".into());
        }
        let width = settings.line_fwhm / 2.;
        let profile = (0..count)
            .map(|k| {
                let x = k as f64 * grid.spacing_cm1;
                let y = sticks
                    .iter()
                    .map(|(s, v)| {
                        v * width / (std::f64::consts::PI * ((x - s).powi(2) + width * width))
                    })
                    .sum();
                (x * CM1_TO_AU * HARTREE_TO_EV, y)
            })
            .collect();
        spectra.push((format!("Resonance Raman · {laser:.3} nm"),HarmonicSpectrum{profile,correlation:vec![],zero_zero_ev:e00*HARTREE_TO_EV,spacing_ev:grid.spacing_cm1*CM1_TO_AU*HARTREE_TO_EV,temperature:0.,integral:f64::NAN,notes:vec!["0 K resonant KHD, one intermediate electronic state; orientationally averaged scattering factors (relative units).".into()]}));
    }
    Ok((spectra, lines))
}

pub fn scattering_factor(a: &[C]) -> f64 {
    let iso = (a[0] + a[4] + a[8]) / 3.;
    let gamma = 0.5 * ((a[0] - a[4]).abs2() + (a[4] - a[8]).abs2() + (a[8] - a[0]).abs2())
        + 0.75 * ((a[1] + a[3]).abs2() + (a[2] + a[6]).abs2() + (a[5] + a[7]).abs2());
    45. * iso.abs2() + 7. * gamma
}
