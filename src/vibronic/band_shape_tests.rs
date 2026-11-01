//! The vibronic band shape, against the closed forms of the Poisson model.
//!
//! The band is a product of Poisson distributions, and a Poisson distribution
//! has exactly known moments. That makes every property of the band checkable
//! analytically rather than by comparison with a stored profile.

use super::band_shape::{vibronic_band, BandKind};
use super::{FranckCondonAnalysis, ModeActivity, CM1_TO_AU};

/// Builds an analysis directly from `(frequency, Huang-Rhys)` pairs, filling
/// the other fields consistently. Only `frequency_cm1` and `huang_rhys` reach
/// the band, but the reorganisation energy has to agree with them or the 0-0
/// energy would be inconsistent with the progression.
fn analysis(modes: &[(f64, f64)]) -> FranckCondonAnalysis {
    let activities: Vec<ModeActivity> = modes
        .iter()
        .enumerate()
        .map(|(index, &(frequency_cm1, huang_rhys))| ModeActivity {
            index,
            frequency_cm1,
            gradient_au: 0.0,
            displacement: (2.0 * huang_rhys).sqrt(),
            huang_rhys,
            reorganisation_cm1: huang_rhys * frequency_cm1,
            resonance_raman: 0.0,
        })
        .collect();
    let reorganisation_cm1 = activities.iter().map(|m| m.reorganisation_cm1).sum();
    FranckCondonAnalysis {
        modes: activities,
        reorganisation_cm1,
        unaccounted_fraction: 0.0,
    }
}

/// An undisplaced mode puts everything in the 0-0 line.
///
/// `S = 0` makes the Poisson distribution a delta at `n = 0`, so the band is a
/// single line, and with no reorganisation energy that line sits exactly at
/// the vertical energy.
#[test]
fn an_undisplaced_mode_gives_only_the_zero_zero_line() {
    let vertical = 0.15;
    let band = vibronic_band(
        &analysis(&[(1500.0, 0.0)]),
        vertical,
        BandKind::Absorption,
        10.0,
        6,
    )
    .unwrap();

    assert!(
        (band.zero_zero_energy - vertical).abs() < 1e-15,
        "with no reorganisation the 0-0 line is the vertical energy"
    );
    assert!(
        (band.profile[0].1 - 1.0).abs() < 1e-12,
        "all intensity belongs in the first line, got {}",
        band.profile[0].1
    );
    let elsewhere: f64 = band.profile[1..].iter().map(|(_, i)| i).sum();
    assert!(
        elsewhere < 1e-12,
        "nothing may fall outside the 0-0 line, got {elsewhere}"
    );
}

/// One displaced mode reproduces the Poisson distribution exactly.
///
/// The `n`-th line of a single mode sits at `E_00 + n w` with intensity
/// `e^{-S} S^n / n!`. Checking the first several lines against that closed
/// form tests the recurrence, the grid placement, and the normalisation at
/// once.
#[test]
fn a_single_mode_follows_the_poisson_progression() {
    const FREQUENCY: f64 = 1000.0;
    const S: f64 = 1.4;
    const RESOLUTION: f64 = 10.0;
    let vertical = 0.2;
    let band = vibronic_band(
        &analysis(&[(FREQUENCY, S)]),
        vertical,
        BandKind::Absorption,
        RESOLUTION,
        12,
    )
    .unwrap();

    let per_quantum = (FREQUENCY / RESOLUTION).round() as usize;
    let mut poisson = (-S).exp();
    for n in 0..8 {
        let (energy, intensity) = band.profile[n * per_quantum];
        let expected_energy = band.zero_zero_energy + n as f64 * FREQUENCY * CM1_TO_AU;
        assert!(
            (energy - expected_energy).abs() < 1e-12,
            "line {n} sits at {energy}, expected {expected_energy}"
        );
        // Relative, not absolute: the progression is truncated at
        // `max_quanta`, and normalising by the truncated sum scales every
        // line up by `1 / (1 - tail)`. The tail is what it is -- a few parts
        // in 1e9 here -- so the model is recovered relatively, not exactly.
        assert!(
            (intensity - poisson).abs() / poisson < 1e-7,
            "line {n} has intensity {intensity}, Poisson says {poisson}"
        );
        poisson *= S / (n as f64 + 1.0);
    }
}

/// The band's first moment is the vertical excitation energy.
///
/// This is the model's central identity and the sharpest single check on it.
/// A Poisson distribution has mean `S`, so the mean energy of one mode's
/// progression is `E_00 + S w`; summed over modes that is
/// `E_00 + lambda = E_vertical`. It holds for any number of modes and any
/// displacements, and it ties the 0-0 energy, the reorganisation energy and
/// the progression together -- get any one of the three wrong and it fails.
#[test]
fn the_mean_energy_of_the_band_is_the_vertical_energy() {
    let vertical = 0.25;
    let band = vibronic_band(
        &analysis(&[(400.0, 2.1), (1200.0, 0.8), (3000.0, 0.15)]),
        vertical,
        BandKind::Absorption,
        5.0,
        20,
    )
    .unwrap();

    let mean = band.mean_energy();
    assert!(
        (mean - vertical).abs() < 1e-6,
        "the intensity-weighted mean is {mean}, the vertical energy {vertical}"
    );
}

/// Emission is the mirror image, so its mean sits as far below the 0-0 line as
/// absorption's sits above, and the gap between the two is the Stokes shift.
#[test]
fn emission_mirrors_absorption_about_the_zero_zero_line() {
    let modes = analysis(&[(800.0, 1.1), (1600.0, 0.6)]);
    let vertical = 0.18;
    let absorption = vibronic_band(&modes, vertical, BandKind::Absorption, 5.0, 16).unwrap();
    let emission = vibronic_band(&modes, vertical, BandKind::Emission, 5.0, 16).unwrap();

    assert!(
        (absorption.zero_zero_energy - emission.zero_zero_energy).abs() < 1e-15,
        "both bands start from the same 0-0 energy"
    );
    let above = absorption.mean_energy() - absorption.zero_zero_energy;
    let below = emission.zero_zero_energy - emission.mean_energy();
    assert!(
        (above - below).abs() < 1e-9,
        "emission must mirror absorption: {above} above vs {below} below"
    );
    let shift = absorption.mean_energy() - emission.mean_energy();
    assert!(
        (shift - absorption.stokes_shift()).abs() < 1e-9,
        "the separation of the two means is the Stokes shift: {shift} vs {}",
        absorption.stokes_shift()
    );
    assert!(
        (absorption.stokes_shift() - 2.0 * modes.reorganisation_cm1 * CM1_TO_AU).abs() < 1e-15,
        "the Stokes shift is twice the reorganisation energy"
    );
}

/// Several modes combine as a product, which shows up as a combination line.
///
/// With two modes the intensity at `E_00 + w_1 + w_2` must be the product of
/// each mode's one-quantum Poisson weight. A band built by adding
/// progressions rather than convolving them would put the wrong number there.
#[test]
fn modes_combine_multiplicatively() {
    const S1: f64 = 0.7;
    const S2: f64 = 1.3;
    const RESOLUTION: f64 = 10.0;
    let band = vibronic_band(
        &analysis(&[(500.0, S1), (1500.0, S2)]),
        0.2,
        BandKind::Absorption,
        RESOLUTION,
        10,
    )
    .unwrap();

    let combination = ((500.0 + 1500.0) / RESOLUTION).round() as usize;
    // At 2000 cm^-1 two routes arrive: one quantum of each mode, or four
    // quanta of the 500 mode. Both belong there, so the expectation is their
    // sum -- which is exactly what makes this a test of the convolution.
    let one_and_one = (S1 * (-S1).exp()) * (S2 * (-S2).exp());
    let four_of_the_first = (-S1).exp() * S1.powi(4) / 24.0 * (-S2).exp();
    let expected = one_and_one + four_of_the_first;
    let found = band.profile[combination].1;
    // Relative for the same truncation reason as the single-mode test.
    assert!(
        (found - expected).abs() / expected < 1e-6,
        "the combination line is {found}, the product model says {expected}"
    );
}

/// The profile is normalised, whatever the displacements.
#[test]
fn the_profile_integrates_to_one() {
    let band = vibronic_band(
        &analysis(&[(300.0, 3.5), (900.0, 1.0)]),
        0.1,
        BandKind::Absorption,
        5.0,
        24,
    )
    .unwrap();
    let total: f64 = band.profile.iter().map(|(_, i)| i).sum();
    assert!((total - 1.0).abs() < 1e-10, "total intensity {total}");
}

/// Nonsense arguments are refused rather than producing an empty or infinite
/// band.
#[test]
fn invalid_arguments_are_refused() {
    let modes = analysis(&[(1000.0, 1.0)]);
    assert!(vibronic_band(&modes, 0.2, BandKind::Absorption, 0.0, 6).is_err());
    assert!(vibronic_band(&modes, 0.2, BandKind::Absorption, 10.0, 0).is_err());
    assert!(vibronic_band(&modes, -0.2, BandKind::Absorption, 10.0, 6).is_err());
}
