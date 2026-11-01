use super::{
    complex_math::{inverse_fft, Lu, C},
    gaussian_correlation::{continuous_branch, GaussianTrace},
    harmonic_model::HarmonicModel,
    thermal_spectrum::{self, SpectrumGrid},
    CM1_TO_AU,
};
use ndarray::Array2;
use std::sync::atomic::{AtomicBool, AtomicUsize};
fn oscillator(w: f64, s: f64) -> HarmonicModel {
    HarmonicModel {
        initial: vec![w],
        final_: vec![w],
        j: Array2::eye(1),
        displacement: vec![(2. * s / w).sqrt()],
        electronic_gap: 0.1,
        initial_labels: vec![0],
        final_labels: vec![0],
    }
}
#[test]
fn complex_pivoted_solve_and_inverse() {
    let a = vec![
        C::new(0., 1.),
        C::new(2., 0.),
        C::new(3., 0.),
        C::new(4., -1.),
    ];
    let lu = Lu::new(a.clone(), 2).unwrap();
    let inv = lu.inverse();
    for i in 0..2 {
        for j in 0..2 {
            let p = (0..2).fold(C::ZERO, |v, k| v + a[i * 2 + k] * inv[k * 2 + j]);
            assert!((p - C::from(if i == j { 1. } else { 0. })).abs() < 1e-12);
        }
    }
}
#[test]
fn inverse_fourier_transform_places_known_frequency() {
    let n = 32;
    let mut a: Vec<_> = (0..n)
        .map(|k| C::new(0., -std::f64::consts::TAU * 3. * k as f64 / n as f64).exp())
        .collect();
    inverse_fft(&mut a);
    for (i, v) in a.iter().enumerate() {
        assert!((*v - C::from(if i == 3 { 1. } else { 0. })).abs() < 1e-12);
    }
}
#[test]
fn gaussian_trace_matches_displaced_oscillator_at_zero_temperature() {
    let w = 0.01;
    let s = 1.2;
    let trace = GaussianTrace::new(oscillator(w, s), 0.).unwrap();
    let mut previous = C::ONE;
    for k in 0..1000 {
        let t = k as f64 * 2.;
        let expected = ((C::new(0., -w * t).exp() - C::ONE) * s).exp();
        let actual = continuous_branch(trace.sample(t, false).unwrap().base, previous);
        previous = actual;
        assert!(
            (actual - expected).abs() < 1e-6,
            "t={t}, actual={actual:?}, expected={expected:?}"
        );
    }
}
#[test]
fn gaussian_trace_includes_thermal_hot_bands() {
    let w = 0.003;
    let s = 0.7;
    let trace = GaussianTrace::new(oscillator(w, s), 800.).unwrap();
    let occupation = trace.occupation(0);
    let mut previous = C::ONE;
    for k in 0..1000 {
        let t = k as f64 * 5.;
        let expected = (((C::new(0., -w * t).exp() - C::ONE) * (occupation + 1.)
            + (C::new(0., w * t).exp() - C::ONE) * occupation)
            * s)
            .exp();
        let actual = continuous_branch(trace.sample(t, false).unwrap().base, previous);
        previous = actual;
        assert!(
            (actual - expected).abs() < 1e-6,
            "t={t}: {actual:?} vs {expected:?}"
        );
    }
}
#[test]
fn equal_undisplaced_surfaces_remain_identity_through_caustics() {
    let trace = GaussianTrace::new(oscillator(0.01, 0.), 400.).unwrap();
    let mut previous = C::ONE;
    for k in 0..4000 {
        let t = k as f64 * 3.;
        let c = continuous_branch(trace.sample(t, false).unwrap().base, previous);
        previous = c;
        assert!((c - C::ONE).abs() < 1e-6, "{t}: {c:?}");
    }
}
#[test]
fn reversed_surface_restores_original_model() {
    let m = oscillator(0.01, 1.);
    let reversed = m.reversed().unwrap().reversed().unwrap();
    assert!((reversed.displacement[0] - m.displacement[0]).abs() < 1e-12);
    assert_eq!(reversed.electronic_gap, m.electronic_gap);
}
#[test]
fn broadened_band_is_positive_and_normalised() {
    let w = 1000. * CM1_TO_AU;
    let s = 0.8;
    let g = SpectrumGrid {
        temperature: 0.,
        spacing_cm1: 10.,
        lorentz_fwhm_cm1: 0.,
        gaussian_fwhm_cm1: 100.,
        ..Default::default()
    };
    let r = thermal_spectrum::calculate(
        oscillator(w, s),
        &g,
        &AtomicBool::new(false),
        &AtomicUsize::new(0),
    )
    .unwrap();
    let area = r.profile.iter().map(|p| p.1 * r.spacing_ev).sum::<f64>();
    assert!((area - 1.).abs() < 1e-6);
    assert!(r.profile.iter().all(|p| p.1 >= 0.));
    let mean = r
        .profile
        .iter()
        .map(|p| p.0 * p.1 * r.spacing_ev)
        .sum::<f64>();
    assert!((mean - (0.1 + s * w) * super::HARTREE_TO_EV).abs() < 1e-4);
}
fn harmonic_functions(x: f64, w: f64, count: usize) -> Vec<f64> {
    let mut a = vec![(w / std::f64::consts::PI).powf(0.25) * (-0.5 * w * x * x).exp()];
    if count > 1 {
        a.push((2. * w).sqrt() * x * a[0]);
    }
    for n in 1..count - 1 {
        a.push(
            (2. * w / (n + 1) as f64).sqrt() * x * a[n]
                - (n as f64 / (n + 1) as f64).sqrt() * a[n - 1],
        );
    }
    a
}
#[test]
fn mixed_unequal_modes_match_independent_coordinate_quadrature() {
    let angle = 0.6_f64;
    let j = Array2::from_shape_vec(
        (2, 2),
        vec![angle.cos(), -angle.sin(), angle.sin(), angle.cos()],
    )
    .unwrap();
    let model = HarmonicModel {
        initial: vec![0.8, 1.3],
        final_: vec![1.1, 0.7],
        j,
        displacement: vec![0.5, -0.3],
        electronic_gap: 3.,
        initial_labels: vec![0, 1],
        final_labels: vec![0, 1],
    };
    let trace = GaussianTrace::new(model.clone(), 0.).unwrap();
    let (inv, _) = model.inverse_rotation().unwrap();
    let levels = 14;
    let mut overlaps = vec![0.; levels * levels];
    let step = 0.14;
    for ix in -55..=55 {
        for iy in -55..=55 {
            let q = [ix as f64 * step, iy as f64 * step];
            let ground = harmonic_functions(q[0], 0.8, 1)[0] * harmonic_functions(q[1], 1.3, 1)[0];
            let qe: [f64; 2] = std::array::from_fn(|i| {
                (0..2)
                    .map(|k| inv[i * 2 + k] * (q[k] - model.displacement[k]))
                    .sum()
            });
            let a = harmonic_functions(qe[0], 1.1, levels);
            let b = harmonic_functions(qe[1], 0.7, levels);
            for m in 0..levels {
                for n in 0..levels {
                    overlaps[m * levels + n] += ground * a[m] * b[n] * step * step;
                }
            }
        }
    }
    assert!((overlaps.iter().map(|v| v * v).sum::<f64>() - 1.).abs() < 1e-6);
    let mut previous = C::ONE;
    for k in 0..200 {
        let t = k as f64 * 0.02;
        let expected = (0..levels)
            .flat_map(|m| (0..levels).map(move |n| (m, n)))
            .fold(C::ZERO, |a, (m, n)| {
                a + C::new(0., -(m as f64 * 1.1 + n as f64 * 0.7) * t).exp()
                    * overlaps[m * levels + n].powi(2)
            });
        let actual = continuous_branch(trace.sample(t, false).unwrap().base, previous);
        previous = actual;
        assert!(
            (actual - expected).abs() < 2e-6,
            "{t}: {actual:?} vs {expected:?}"
        );
    }
}
#[test]
fn native_orca_spectrum_parser_converts_only_the_energy_axis() {
    let r = super::orca_esd::parse_spectrum(
        "Energy TotalSpectrum FC HT\n20000 2 2 0\n20010 3 3 0\n",
        298.,
    )
    .unwrap();
    assert!((r.profile[0].0 - 20000. * CM1_TO_AU * super::HARTREE_TO_EV).abs() < 1e-12);
    assert_eq!(r.profile[0].1, 2.);
    assert!(super::orca_esd::parse_spectrum("20000 NaN", 0.).is_err());
}
#[test]
#[ignore = "requires BEAVYR_TEST_ORCA; runs a small native ORCA ESD job"]
fn native_orca_esd_absorption_smoke() {
    let binary = std::env::var("BEAVYR_TEST_ORCA").expect("Set BEAVYR_TEST_ORCA");
    let directory = super::process::directory().unwrap();
    let input="! HF STO-3G Opt Freq\n%pal nprocs 1 end\n* xyz 0 1\nO 0 0 0\nH 0.76 0.58 0\nH -0.76 0.58 0\n*\n";
    std::fs::write(directory.join("orca_job.inp"), input).unwrap();
    let cancel = AtomicBool::new(false);
    let output = super::process::run(
        crate::qchem_interfaces::orca_run::run_command(std::path::Path::new(&binary), &directory),
        &directory,
        "orca_job",
        &cancel,
    )
    .unwrap();
    crate::qchem_interfaces::orca_run::check_output(&output, true).unwrap();
    let source = super::input::HessianSource {
        path: directory.join("orca_job.hess").display().to_string(),
        xyz: String::new(),
    };
    let mut settings = super::state::Settings {
        kind: super::state::AnalysisKind::Spectrum,
        ground: source.clone(),
        excited: source,
        ..Default::default()
    };
    settings.spectrum.engine = super::spectral_job::Engine::Orca;
    settings.spectrum.orca_path = binary;
    settings.spectrum.qc_method = crate::qchem_interfaces::job::QcMethod::Hf;
    settings.spectrum.method.basis = "STO-3G".into();
    settings.spectrum.method.nproc = 1;
    settings.spectrum.roots = 1;
    settings.spectrum.grid.temperature = 0.;
    for process in [
        super::optical_spectrum::Process::Absorption,
        super::optical_spectrum::Process::Fluorescence,
        super::optical_spectrum::Process::Raman,
    ] {
        settings.spectrum.process = process;
        if process == super::optical_spectrum::Process::Raman {
            settings.spectrum.raman.lasers = "532, 500".into();
        }
        let result = super::spectral_job::run(&settings, &cancel, &AtomicUsize::new(0)).unwrap();
        assert_eq!(result.spectra.len(), if process == super::optical_spectrum::Process::Raman {2} else {1});
        if process != super::optical_spectrum::Process::Raman {
            assert!(result.spectra[0].1.profile.iter().any(|p| p.1 > 0.));
        }
        println!(
            "Native ESD {} smoke directory: {}",
            process.label(),
            result.workdir.unwrap()
        );
    }
}
#[test]
fn fluorescence_progression_uses_the_emitting_surface_and_correct_energy_sign() {
    let w = 1000. * CM1_TO_AU;
    let s = 0.8;
    let grid = SpectrumGrid {
        temperature: 0.,
        spacing_cm1: 10.,
        lorentz_fwhm_cm1: 0.,
        gaussian_fwhm_cm1: 100.,
        ..Default::default()
    };
    let cancel = AtomicBool::new(false);
    let progress = AtomicUsize::new(0);
    let abs = super::optical_spectrum::calculate(
        oscillator(w, s),
        super::optical_spectrum::Process::Absorption,
        &grid,
        false,
        &cancel,
        &progress,
    )
    .unwrap();
    let em = super::optical_spectrum::calculate(
        oscillator(w, s),
        super::optical_spectrum::Process::Fluorescence,
        &grid,
        false,
        &cancel,
        &progress,
    )
    .unwrap();
    let mean = |r: &super::thermal_spectrum::HarmonicSpectrum| {
        r.profile
            .iter()
            .map(|p| p.0 * p.1 * r.spacing_ev)
            .sum::<f64>()
    };
    assert!((mean(&abs) - (0.1 + s * w) * super::HARTREE_TO_EV).abs() < 1e-4);
    assert!((mean(&em) - (0.1 - s * w) * super::HARTREE_TO_EV).abs() < 1e-4);
    assert!((abs.zero_zero_ev - em.zero_zero_ev).abs() < 1e-12);
}
#[test]
fn unequal_frequency_emission_reversal_is_normalised() {
    let mut m = oscillator(0.008, 0.5);
    m.final_[0] = 0.012;
    let trace = GaussianTrace::new(m.reversed().unwrap(), 500.).unwrap();
    assert_eq!(trace.model.initial[0], 0.012);
    assert_eq!(trace.model.final_[0], 0.008);
    assert_eq!(trace.sample(0., true).unwrap().base, C::ONE);
}

#[test]
fn resonance_raman_matches_independent_sum_over_states() {
    let w = 1000. * CM1_TO_AU;
    let s = 0.7;
    let model = oscillator(w, s);
    let laser = 0.102;
    let nm = 1e7 * CM1_TO_AU / laser;
    let width = 300.;
    let grid = SpectrumGrid {
        temperature: 0.,
        lorentz_fwhm_cm1: width,
        gaussian_fwhm_cm1: 0.,
        max_points: 262144,
        ..Default::default()
    };
    let settings = super::raman_spectrum::RamanSettings {
        lasers: nm.to_string(),
        dipole: [0., 0., 1.],
        order: 2,
        ..Default::default()
    };
    let (_, lines) = super::raman_spectrum::calculate(
        model,
        &grid,
        &settings,
        &AtomicBool::new(false),
        &AtomicUsize::new(0),
    )
    .unwrap();
    let gamma = width * CM1_TO_AU / 2.;
    let mut expected = C::ZERO;
    let mut probability = (-s).exp();
    let mut previous = 0.;
    for n in 0..60 {
        expected += C::from(s.sqrt() * (probability - previous))
            / C::new(0.1 + n as f64 * w - laser, -gamma);
        previous = probability;
        probability *= s / (n + 1) as f64;
    }
    assert!(
        (lines[0].strength / (12. * expected.abs2()) - 1.).abs() < 0.001,
        "{} vs {}",
        lines[0].strength,
        12. * expected.abs2()
    );
    assert_eq!(lines[1].modes, vec![1, 1]);
}
#[test]
fn undisplaced_condon_model_has_no_inelastic_raman() {
    let grid = SpectrumGrid {
        temperature: 0.,
        lorentz_fwhm_cm1: 500.,
        ..Default::default()
    };
    let settings = super::raman_spectrum::RamanSettings {
        dipole: [1., 0., 0.],
        order: 2,
        ..Default::default()
    };
    let (_, lines) = super::raman_spectrum::calculate(
        oscillator(1000. * CM1_TO_AU, 0.),
        &grid,
        &settings,
        &AtomicBool::new(false),
        &AtomicUsize::new(0),
    )
    .unwrap();
    assert!(lines.iter().all(|l| l.strength < 1e-10));
}
#[test]
fn herzberg_teller_dark_oscillator_matches_thermal_exact_correlation() {
    let w = 0.003;
    let trace = GaussianTrace::new(oscillator(w, 0.), 800.).unwrap();
    let d = super::transition_properties::Dipole {
        value: [0.; 3],
        derivative: vec![[1., 0., 0.]],
    };
    let occupation = trace.occupation(0);
    let mut previous = C::ONE;
    for k in 0..300 {
        let t = k as f64 * 3.;
        let m = trace.sample(t, true).unwrap();
        let base = continuous_branch(m.base, previous);
        previous = base;
        let actual = base * d.optical(&m);
        let expected = (C::new(0., -w * t).exp() * (occupation + 1.)
            + C::new(0., w * t).exp() * occupation)
            / (2. * w);
        assert!(
            (actual - expected).abs() < 1e-4,
            "{actual:?} vs {expected:?}"
        );
    }
}
#[test]
fn herzberg_teller_fluorescence_transforms_dipole_origin_and_derivatives() {
    let model = oscillator(0.01, 1.);
    let d = super::transition_properties::Dipole {
        value: [1., 2., 3.],
        derivative: vec![[0.5, 0., 0.]],
    };
    let reverse = d.reversed(&model);
    assert!((reverse.value[0] - 1. - 0.5 * model.displacement[0]).abs() < 1e-12);
    let back = reverse.reversed(&model.reversed().unwrap());
    assert!((back.value[0] - 1.).abs() < 1e-12);
}
#[test]
fn raman_tensor_scattering_factor_is_rotation_invariant() {
    let d = [1., 2., 3.];
    let a: Vec<_> = d
        .iter()
        .flat_map(|x| d.iter().map(move |y| C::from(x * y)))
        .collect();
    assert!((super::raman_spectrum::scattering_factor(&a) - 12. * 14_f64.powi(2)).abs() < 1e-12);
}
