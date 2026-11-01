//! Franck-Condon activity from an excited-state gradient.
//!
//! The vertical-gradient picture, also called the independent-mode displaced
//! harmonic oscillator: the excited-state surface is taken to have the ground
//! state's curvature, displaced by `-g_k / w_k^2` along each normal
//! coordinate. Everything else follows from that one assumption.
//!
//! # Why a gradient is enough
//!
//! No excited-state optimisation is needed, and that is the point. An excited
//! minimum is often hard to converge, and near a conical intersection it may
//! not exist in any useful sense; a gradient at the Franck-Condon geometry
//! always does.
//!
//! # The quantities
//!
//! With `L` the ground-state modes normalised so `sum_i m_i L_ik^2 = 1`, and
//! `g` the excited-state Cartesian gradient at the ground-state geometry,
//!
//! ```text
//! g_k    = sum_i L_ik g_i            the gradient along mode k
//! delta_k = -g_k / w_k^{3/2}         displacement, dimensionless
//! S_k    = delta_k^2 / 2             Huang-Rhys factor
//! lambda_k = S_k w_k                 that mode's reorganisation energy
//! I_k    ~ w_k^2 delta_k^2           resonance-Raman intensity, short time
//! ```
//!
//! The mass weighting cancels in the projection. In mass-weighted coordinates
//! it is `sum_i (sqrt(m_i) L_ik)(g_i / sqrt(m_i))`, which is a plain dot
//! product of the stored Cartesian column with the raw gradient -- no masses
//! appear. They reappear only in the diagnostic below.
//!
//! `delta_k` is signed on purpose: the direction matters when two states are
//! compared, and squaring it early throws that away.
//!
//! # What the intensity is and is not
//!
//! `I_k` is a pure Franck-Condon intensity in the Condon approximation. No
//! transition dipole enters, so it gives **relative** intensities within one
//! spectrum and not absolute cross-sections, and it cannot describe
//! Herzberg-Teller (non-Condon) activity. Use
//! [`FranckCondonAnalysis::normalised_resonance_raman`] rather than the raw
//! number, which carries an arbitrary constant.

use anyhow::{bail, ensure, Result};
use ndarray::Array2;

use super::{AMU_TO_ELECTRON_MASS, CM1_TO_AU};

/// What one ground-state mode does when the molecule is placed on an excited
/// state's surface at the Franck-Condon geometry.
#[derive(Debug, Clone, PartialEq)]
pub struct ModeActivity {
    /// Index into the ground state's mode list.
    pub index: usize,
    pub frequency_cm1: f64,
    /// The excited-state gradient projected onto this mode, in atomic units.
    pub gradient_au: f64,
    /// Displacement of the excited-state minimum along this mode, in
    /// dimensionless normal coordinates: `-g_k / w_k^{3/2}`. Signed.
    pub displacement: f64,
    /// Huang-Rhys factor `S_k = delta_k^2 / 2`: the mean number of quanta of
    /// this mode excited in a vertical transition.
    pub huang_rhys: f64,
    /// This mode's share of the reorganisation energy, `S_k w_k`, in cm^-1.
    pub reorganisation_cm1: f64,
    /// Resonance-Raman intensity in the short-time approximation, before
    /// normalisation. Proportional to `w_k^2 delta_k^2`; the constant of
    /// proportionality is arbitrary and identical for every mode.
    pub resonance_raman: f64,
}

/// The whole Franck-Condon picture for one excited state.
#[derive(Debug, Clone, PartialEq)]
pub struct FranckCondonAnalysis {
    pub modes: Vec<ModeActivity>,
    /// Total reorganisation energy in cm^-1, summed over the modes analysed.
    pub reorganisation_cm1: f64,
    /// The part of the gradient no vibration accounts for, as a fraction of
    /// its norm in the mass-weighted metric.
    ///
    /// It should be small. A large value means the gradient points somewhere
    /// the mode list cannot represent: translations and rotations projected
    /// out of the Hessian, imaginary modes excluded from the analysis, or --
    /// the case worth catching -- a gradient computed at a different geometry
    /// from the Hessian.
    pub unaccounted_fraction: f64,
}

impl FranckCondonAnalysis {
    /// Modes ordered by resonance-Raman intensity, strongest first.
    #[allow(dead_code)] // Retained Behemoth analysis API.
    pub fn by_resonance_raman(&self) -> Vec<&ModeActivity> {
        let mut ordered: Vec<&ModeActivity> = self.modes.iter().collect();
        ordered.sort_by(|a, b| b.resonance_raman.total_cmp(&a.resonance_raman));
        ordered
    }

    /// `(frequency_cm1, intensity)` with the strongest intensity scaled to
    /// one, which is the only form comparable with a measured spectrum.
    pub fn normalised_resonance_raman(&self) -> Vec<(f64, f64)> {
        let peak = self
            .modes
            .iter()
            .map(|mode| mode.resonance_raman)
            .fold(0.0_f64, f64::max);
        self.modes
            .iter()
            .map(|mode| {
                let value = if peak > 0.0 {
                    mode.resonance_raman / peak
                } else {
                    0.0
                };
                (mode.frequency_cm1, value)
            })
            .collect()
    }
}

/// Projects an excited-state gradient onto the ground-state normal modes.
///
/// `modes` are the ground state's, Cartesian columns normalised so that
/// `sum_i m_i L_ik^2 = 1`. `gradient_cart` is the excited state's Cartesian
/// gradient in hartree per bohr **at the ground-state geometry**, flat, three
/// components per atom. `which` selects the genuine vibrations, since the
/// projected translations and rotations carry no meaning here.
pub fn franck_condon_activity(
    modes: &Array2<f64>,
    frequencies_cm1: &[f64],
    gradient_cart: &[f64],
    masses_amu: &[f64],
    which: &[usize],
) -> Result<FranckCondonAnalysis> {
    let (dimensions, columns) = modes.dim();
    ensure!(
        gradient_cart.len() == dimensions,
        "the gradient has {} components but the modes have {dimensions}",
        gradient_cart.len()
    );
    ensure!(
        frequencies_cm1.len() == columns,
        "{} frequencies for {columns} modes",
        frequencies_cm1.len()
    );
    // One mass per atom, not per coordinate. Without this a caller passing a
    // 3N-length list would index the first N entries and silently corrupt the
    // diagnostic below rather than failing.
    ensure!(
        masses_amu.len() * 3 == dimensions,
        "{} masses for {dimensions} Cartesian coordinates; one mass per atom is required",
        masses_amu.len()
    );
    ensure!(
        gradient_cart.iter().all(|value| value.is_finite()),
        "the excited-state gradient contains a non-finite component"
    );

    ensure!(
        masses_amu.iter().all(|m| m.is_finite() && *m > 0.0),
        "Masses must be finite and positive"
    );
    ensure!(
        modes.iter().chain(frequencies_cm1).all(|v| v.is_finite()),
        "Normal modes must be finite"
    );
    let mut activities = Vec::with_capacity(which.len());
    let mut recovered = 0.0_f64;
    for &mode in which {
        ensure!(mode < columns, "mode index {mode} is out of range");
        let frequency_cm1 = frequencies_cm1[mode];
        // An imaginary mode has no displaced minimum and no Franck-Condon
        // progression. Reporting one would be inventing a number.
        if frequency_cm1 <= 0.0 {
            continue;
        }
        let gradient_au: f64 = (0..dimensions)
            .map(|index| modes[(index, mode)] * gradient_cart[index])
            .sum();
        recovered += gradient_au * gradient_au;

        let frequency_au = frequency_cm1 * CM1_TO_AU;
        let displacement = -gradient_au / frequency_au.powf(1.5);
        let huang_rhys = 0.5 * displacement * displacement;
        activities.push(ModeActivity {
            index: mode,
            frequency_cm1,
            gradient_au,
            displacement,
            huang_rhys,
            reorganisation_cm1: huang_rhys * frequency_cm1,
            // `w^2 delta^2`, with the frequency in atomic units so the number
            // is a consistent atomic-unit quantity rather than a mixture.
            // Only ratios are meaningful either way, but a mixed-unit value
            // invites being read as something it is not.
            resonance_raman: frequency_au * frequency_au * displacement * displacement,
        });
    }

    // The comparison must be in the metric the modes are orthonormal in, which
    // is the mass-weighted one: `sum_k g_k^2` recovers `|M^-1/2 g|^2`, not
    // `|g|^2`. Comparing against the Cartesian norm would report every
    // gradient as almost entirely unaccounted for, whatever the projection
    // actually did.
    let total: f64 = (0..dimensions)
        .map(|index| {
            let mass = masses_amu[index / 3] * AMU_TO_ELECTRON_MASS;
            gradient_cart[index] * gradient_cart[index] / mass
        })
        .sum();
    if !total.is_finite() {
        bail!("the mass-weighted gradient norm is not finite");
    }
    let unaccounted_fraction = if total > 0.0 {
        (1.0 - recovered / total).max(0.0).sqrt()
    } else {
        0.0
    };

    let reorganisation_cm1 = activities.iter().map(|mode| mode.reorganisation_cm1).sum();
    Ok(FranckCondonAnalysis {
        modes: activities,
        reorganisation_cm1,
        unaccounted_fraction,
    })
}
