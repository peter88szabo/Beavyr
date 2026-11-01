//! The vibronic band shape a set of Huang-Rhys factors implies.
//!
//! The Franck-Condon analysis produces, for each mode, a Huang-Rhys factor
//! `S_k`. Those numbers already contain the whole absorption band at zero
//! temperature; this evaluates it.
//!
//! # The model
//!
//! Independent displaced harmonic oscillators with the ground state's
//! frequencies. For one mode, the intensity of the transition from `v = 0` to
//! `v = n` is a Poisson distribution in the Huang-Rhys factor,
//!
//! ```text
//! I_n = e^{-S} S^n / n!
//! ```
//!
//! and for several modes the band is the product over modes, placed at
//!
//! ```text
//! E = E_00 + sum_k n_k w_k          absorption
//! E = E_00 - sum_k n_k w_k          emission, the mirror image
//! ```
//!
//! The 0-0 energy is not the vertical excitation. In this model the excited
//! surface has the ground state's curvature displaced by `delta_k`, so the
//! vertical energy sits `lambda = sum_k S_k w_k` above the adiabatic one:
//!
//! ```text
//! E_00 = E_vertical - lambda
//! ```
//!
//! and the Stokes shift between absorption and emission maxima is `2 lambda`.
//!
//! # Evaluating the product
//!
//! Enumerating every `(n_1, n_2, ...)` is combinatorially hopeless past a
//! handful of modes. Instead the band is built on an energy grid by
//! convolution: start with a single line at `E_00`, and for each mode in turn
//! convolve with that mode's Poisson progression. That costs
//! `modes x grid x n_max` and is exact to the grid's resolution.
//!
//! # What is and is not included
//!
//! Zero temperature, so no hot bands. No Duschinsky mixing and no frequency
//! change on excitation -- the model that produced `S_k` assumed neither, and
//! including them here would be inconsistent with it. No Herzberg-Teller
//! intensity, which lives in a different term. The result is the Condon,
//! displaced-oscillator band, which for most rigid chromophores is the bulk of
//! the observed shape.

use anyhow::{ensure, Result};

use super::{FranckCondonAnalysis, CM1_TO_AU};

/// Whether the band is absorption or its mirror-image emission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BandKind {
    Absorption,
    #[allow(dead_code)] // Tested mirror profile; UI exposes absorption only.
    Emission,
}

/// A vibronic band on an energy grid.
#[derive(Debug, Clone)]
#[allow(dead_code)] // Retain Behemoth band metadata for callers and tests.
pub struct VibronicBand {
    pub kind: BandKind,
    /// The adiabatic 0-0 energy, in hartree.
    pub zero_zero_energy: f64,
    /// The vertical energy the analysis started from, in hartree.
    pub vertical_energy: f64,
    /// Total reorganisation energy `lambda`, in hartree.
    pub reorganisation: f64,
    /// Energies in hartree and the band intensity at each, normalised so the
    /// total integrated intensity is one.
    pub profile: Vec<(f64, f64)>,
    /// Spacing of the energy grid, in hartree.
    pub resolution: f64,
}

#[allow(dead_code)] // Retain Behemoth diagnostic helpers.
impl VibronicBand {
    /// The Stokes shift this model implies, `2 lambda`, in hartree.
    pub fn stokes_shift(&self) -> f64 {
        2.0 * self.reorganisation
    }

    /// The energy of the strongest point on the profile.
    pub fn maximum(&self) -> Option<f64> {
        self.profile
            .iter()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(energy, _)| *energy)
    }

    /// The intensity-weighted mean energy, which for absorption is the
    /// vertical excitation energy in this model -- a check the test uses.
    pub fn mean_energy(&self) -> f64 {
        let total: f64 = self.profile.iter().map(|(_, i)| i).sum();
        if total <= 0.0 {
            return 0.0;
        }
        self.profile.iter().map(|(e, i)| e * i).sum::<f64>() / total
    }

    /// The profile in electron-volts, for plotting.
    pub fn profile_ev(&self) -> Vec<(f64, f64)> {
        const HARTREE_TO_EV: f64 = 27.211_386_245_988;
        self.profile
            .iter()
            .map(|(e, i)| (e * HARTREE_TO_EV, *i))
            .collect()
    }
}

/// Builds the zero-temperature vibronic band implied by a Franck-Condon
/// analysis.
///
/// `vertical_energy` is the excitation energy the analysis was made at, in
/// hartree. `resolution_cm1` sets the energy grid; something below the
/// smallest active frequency keeps the progressions resolved. `max_quanta`
/// caps the progression per mode and needs to be a little beyond the largest
/// `S_k`, where the Poisson distribution peaks.
pub fn vibronic_band(
    analysis: &FranckCondonAnalysis,
    vertical_energy: f64,
    kind: BandKind,
    resolution_cm1: f64,
    max_quanta: usize,
) -> Result<VibronicBand> {
    ensure!(
        resolution_cm1.is_finite() && resolution_cm1 > 0.0,
        "the energy resolution must be positive"
    );
    ensure!(max_quanta >= 1, "at least one quantum per mode is needed");
    ensure!(
        vertical_energy.is_finite() && vertical_energy > 0.0,
        "the vertical energy must be finite and positive"
    );

    ensure!(
        max_quanta <= 1000,
        "At most 1000 quanta per mode are supported"
    );
    ensure!(
        analysis.modes.iter().all(|m| m.frequency_cm1.is_finite()
            && m.frequency_cm1 > 0.0
            && m.huang_rhys.is_finite()
            && m.huang_rhys >= 0.0),
        "Invalid mode frequency or Huang-Rhys factor"
    );
    let resolution = resolution_cm1 * CM1_TO_AU;
    let reorganisation = analysis.reorganisation_cm1 * CM1_TO_AU;
    let zero_zero_energy = vertical_energy - reorganisation;

    // The band extends from E_00 up (or down) by the sum of every mode's
    // longest progression, which is a generous bound.
    let span: f64 = analysis
        .modes
        .iter()
        .map(|mode| max_quanta as f64 * mode.frequency_cm1 * CM1_TO_AU)
        .sum::<f64>()
        + 2.0 * resolution;
    let grid_points = (span / resolution).ceil() + 1.0;
    ensure!(
        grid_points.is_finite() && grid_points <= 250_000.0,
        "Band grid is too large; increase the grid spacing or reduce the quanta limit"
    );
    ensure!(
        grid_points * analysis.modes.len() as f64 * (max_quanta + 1) as f64 <= 200_000_000.0,
        "Band calculation is too large; increase grid spacing or reduce the quanta limit"
    );
    let points = grid_points as usize;

    // The line at E_00, then one convolution per mode. Offsets are in grid
    // points, positive for absorption and negative for emission.
    let mut band = vec![0.0f64; points];
    band[0] = 1.0;
    for mode in &analysis.modes {
        let s = mode.huang_rhys;
        let quantum_points = mode.frequency_cm1 * CM1_TO_AU / resolution;
        // Poisson weights e^{-S} S^n / n!, built by recurrence to avoid
        // factorial overflow.
        let mut weights = Vec::with_capacity(max_quanta + 1);
        let mut weight = (-s).exp();
        for n in 0..=max_quanta {
            weights.push(weight);
            weight *= s / (n as f64 + 1.0);
        }
        let mut next = vec![0.0f64; points];
        for (index, &intensity) in band.iter().enumerate() {
            if intensity == 0.0 {
                continue;
            }
            for (n, &w) in weights.iter().enumerate() {
                // Nearest grid point for n quanta of this mode. Rounding rather
                // than interpolating keeps every line's full weight on the
                // grid; the resolution argument is what controls the error.
                let shift = (n as f64 * quantum_points).round() as usize;
                if index + shift < points {
                    next[index + shift] += intensity * w;
                }
            }
        }
        band = next;
    }

    let total: f64 = band.iter().sum();
    ensure!(total > 0.0, "the band has no intensity");
    let sign = match kind {
        BandKind::Absorption => 1.0,
        BandKind::Emission => -1.0,
    };
    let profile: Vec<(f64, f64)> = band
        .iter()
        .enumerate()
        .map(|(i, &intensity)| {
            (
                zero_zero_energy + sign * i as f64 * resolution,
                intensity / total,
            )
        })
        .collect();

    Ok(VibronicBand {
        kind,
        zero_zero_energy,
        vertical_energy,
        reorganisation,
        profile,
        resolution,
    })
}
