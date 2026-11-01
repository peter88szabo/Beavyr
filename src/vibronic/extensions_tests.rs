use super::{
    complex_math::C,
    gaussian_correlation::{continuous_branch, GaussianTrace},
    harmonic_model::HarmonicModel,
    spectrum_view::{self, Axis},
    thermal_spectrum::SpectrumGrid,
    transition_properties::{self, Dipole},
    CM1_TO_AU,
};
use ndarray::Array2;
use std::sync::atomic::{AtomicBool, AtomicUsize};
fn model(w: f64, s: f64) -> HarmonicModel {
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
fn experimental_data_units_and_invalid_rows() {
    let a =
        spectrum_view::parse_experimental("wavelength,intensity\n500,1\n600,2\n", Axis::Wavelength)
            .unwrap();
    assert!(a[0].0 < a[1].0);
    assert!((a[1].0 - 2.479683968664).abs() < 1e-10);
    assert!(spectrum_view::parse_experimental("0,1\n500,2", Axis::Wavelength).is_err());
    assert!(spectrum_view::parse_experimental("x,y\nwrong,3\n500,2", Axis::Ev).is_err());
    assert!(spectrum_view::parse_experimental("1,NaN\n2,3", Axis::Ev).is_err());
}
#[test]
fn display_shift_and_normalisation_leave_originals_unchanged() {
    let original = vec![(2., 3.), (3., 6.)];
    let displayed = spectrum_view::displayed(&original, Axis::Ev, 0.5, true);
    assert_eq!(displayed, vec![(2.5, 0.5), (3.5, 1.)]);
    assert_eq!(original, vec![(2., 3.), (3., 6.)]);
    let csv = spectrum_view::csv(&[("Computed", displayed)], Axis::Ev, 0.5, true);
    assert!(csv.contains("shift_eV: 0.5000000000"));
}
#[test]
fn svg_export_escapes_labels_and_has_both_series() {
    let curves = vec![
        ("A&B <test>", vec![(1., 0.), (2., 1.)]),
        ("Experiment", vec![(1., 0.2), (2., 0.5)]),
    ];
    let svg = spectrum_view::svg(&curves, Axis::Ev, "<spectrum>");
    assert!(svg.contains("&lt;spectrum&gt;"));
    assert!(svg.contains("A&amp;B &lt;test&gt;"));
    assert_eq!(svg.matches("<polyline").count(), 2);
    assert!(svg.ends_with("</svg>"));
}
#[test]
fn property_import_checks_geometry_units_and_projects_derivatives() {
    let a = super::input::StateModes {
        atoms: vec!["H".into()],
        masses: vec![1.],
        coords: vec![0., 0., 0.],
        modes: Array2::from_shape_vec((3, 1), vec![0.5, 0., 0.]).unwrap(),
        frequencies: vec![1000.],
        selected: vec![0],
        imaginary: 0,
    };
    let text = "BEAVYR_TRANSITION_V1\n1\nH 0 0 0\n1 2 3\n4 5 6\n0 0 0\n0 0 0\n";
    let d = transition_properties::parse(text, &a, true).unwrap();
    assert_eq!(d.derivative, vec![[2., 2.5, 3.]]);
    assert!(
        transition_properties::parse(&text.replace("H 0 0 0", "H 0.001 0 0"), &a, true).is_err()
    );
    assert!(transition_properties::parse(&transition_properties::template(&a), &a, true).is_err());
}
#[test]
fn herzberg_teller_interference_matches_independent_coordinate_integration() {
    let w = 0.8;
    let model = model(w, 0.7);
    let trace = GaussianTrace::new(model.clone(), 0.).unwrap();
    let dipole = Dipole {
        value: [0.6, 0., 0.],
        derivative: vec![[0.2, 0., 0.]],
    };
    let levels = 25;
    let dx = 0.03;
    let mut overlaps = vec![0.; levels];
    for k in -250..=250 {
        let x = k as f64 * dx;
        let g = (w / std::f64::consts::PI).powf(0.25) * (-w * x * x / 2.).exp();
        let y = x - model.displacement[0];
        let mut h = vec![(w / std::f64::consts::PI).powf(0.25) * (-w * y * y / 2.).exp()];
        h.push((2. * w).sqrt() * y * h[0]);
        for n in 1..levels - 1 {
            h.push(
                (2. * w / (n + 1) as f64).sqrt() * y * h[n]
                    - (n as f64 / (n + 1) as f64).sqrt() * h[n - 1],
            );
        }
        for n in 0..levels {
            overlaps[n] += dx * g * (0.6 + 0.2 * x) * h[n];
        }
    }
    let mut previous = C::ONE;
    for k in 0..150 {
        let t = k as f64 * 0.03;
        let moments = trace.sample(t, true).unwrap();
        let base = continuous_branch(moments.base, previous);
        previous = base;
        let expected = overlaps.iter().enumerate().fold(C::ZERO, |s, (n, v)| {
            s + C::new(0., -(n as f64) * w * t).exp() * (v * v)
        });
        assert!((base * dipole.optical(&moments) - expected).abs() < 1e-6);
    }
}
#[test]
fn herzberg_teller_raman_fundamental_matches_two_intermediate_levels() {
    let w = 1000. * CM1_TO_AU;
    let laser = 0.102;
    let gamma = 300. * CM1_TO_AU / 2.;
    let grid = SpectrumGrid {
        temperature: 0.,
        lorentz_fwhm_cm1: 300.,
        max_points: 262144,
        ..Default::default()
    };
    let settings = super::raman_spectrum::RamanSettings {
        lasers: (1e7 * CM1_TO_AU / laser).to_string(),
        ..Default::default()
    };
    let dipole = Dipole {
        value: [1., 0., 0.],
        derivative: vec![[0.02, 0., 0.]],
    };
    let (_, lines) = super::raman_spectrum::calculate_with_dipole(
        model(w, 0.),
        &grid,
        &settings,
        &dipole,
        &AtomicBool::new(false),
        &AtomicUsize::new(0),
    )
    .unwrap();
    let exact = (C::ONE / C::new(0.1 - laser, -gamma) + C::ONE / C::new(0.1 + w - laser, -gamma))
        * (0.02 / (2. * w).sqrt());
    assert!((lines[0].strength / (12. * exact.abs2()) - 1.).abs() < 0.001);
}
#[test]
fn spectrum_cancellation_returns_no_partial_result() {
    let g = SpectrumGrid::default();
    assert!(super::thermal_spectrum::calculate(
        model(1000. * CM1_TO_AU, 0.4),
        &g,
        &AtomicBool::new(true),
        &AtomicUsize::new(0)
    )
    .unwrap_err()
    .starts_with("Cancelled."));
}
#[test]
fn animation_closure_restores_equilibrium_and_leaves_other_trajectories() {
    use crate::trajectory::{TrajectoryFrame, TrajectoryState};
    use bevy::prelude::Vec3;
    let eq = TrajectoryFrame { metadata: Default::default(),
        atoms: vec!["H".into()],
        pos: vec![Vec3::ZERO],
    };
    let displaced = TrajectoryFrame { metadata: Default::default(),
        atoms: eq.atoms.clone(),
        pos: vec![Vec3::X],
    };
    let mut animation = super::mode_animation::Animation::default();
    let mut molecule = crate::molecule::Molecule::empty();
    let mut traj = TrajectoryState::default();
    animation.pending = Some((eq.clone(), vec![displaced.clone(); 60]));
    assert!(animation.apply(true, &mut molecule, &mut traj));
    assert!(traj.playing);
    assert!(animation.apply(false, &mut molecule, &mut traj));
    assert!(!traj.playing);
    assert_eq!(molecule.pos, eq.pos);
    animation.pending = Some((eq, vec![displaced; 60]));
    animation.apply(true, &mut molecule, &mut traj);
    traj.clear_for_structure();
    assert!(!animation.apply(false, &mut molecule, &mut traj));
}
