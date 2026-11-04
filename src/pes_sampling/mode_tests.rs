use super::{input, state::*, ts_mode};
use crate::qchem_interfaces::{hessian_io, xtb_freq};
use ndarray::Array2;

#[test]
fn eckart_modes_and_frequencies_match_frequency_analysis_for_linear_and_nonlinear_molecules() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let xyz = hessian_io::geometry(
        &std::fs::read_to_string(root.join("pyscf_pbe_water_hessian.xyz")).unwrap(),
    )
    .unwrap();
    let nonlinear =
        hessian_io::load(&root.join("pyscf_pbe_water_hessian.txt"), Some(&xyz)).unwrap();
    let linear = hessian_io::from_geometry(
        vec!["O".into(), "C".into(), "O".into()],
        vec![-2.2, 0., 0., 0., 0., 0., 2.2, 0., 0.],
        Array2::eye(9),
    )
    .unwrap();
    for (case, data) in [nonlinear, linear].iter().enumerate() {
        let f = input::from_hessian(
            data.atoms.clone(),
            data.masses_amu.clone(),
            data.coords_bohr.clone(),
            data.hessian.clone(),
            "test".into(),
        )
        .unwrap();
        let path = std::env::temp_dir().join(format!(
            "beavyr-eckart-parity-{}-{case}.hess",
            std::process::id()
        ));
        hessian_io::save(&path, data).unwrap();
        let mut panel = xtb_freq::XtbFreqPanelState::default();
        let mut task = xtb_freq::XtbFrequencyTask::default();
        assert!(
            xtb_freq::load_hessian_from_path(&path, &mut panel, &mut task),
            "{:?}",
            task.last_message
        );
        let result = task.result.unwrap();
        assert_eq!(f.frequencies.len(), if case == 0 { 3 } else { 4 });
        assert_eq!(f.frequencies.len(), result.positive_indices.len());
        for (j, &column) in result.positive_indices.iter().enumerate() {
            assert!(
                (f.frequencies[j] / input::CM_TO_AU - result.frequencies_cm1[column]).abs() < 1e-8
            );
            for row in 0..f.equilibrium.len() {
                assert!((f.modes[(row, j)] - result.modes[(row, column)]).abs() < 1e-12);
            }
            let mut translation = [0.; 3];
            let mut rotation = [0.; 3];
            for a in 0..f.atoms.len() {
                let r = &f.equilibrium[3 * a..3 * a + 3];
                let v = [
                    f.modes[(3 * a, j)],
                    f.modes[(3 * a + 1, j)],
                    f.modes[(3 * a + 2, j)],
                ];
                for k in 0..3 {
                    translation[k] += f.masses[a] * v[k];
                }
                let cross = [
                    r[1] * v[2] - r[2] * v[1],
                    r[2] * v[0] - r[0] * v[2],
                    r[0] * v[1] - r[1] * v[0],
                ];
                for k in 0..3 {
                    rotation[k] += f.masses[a] * cross[k];
                }
            }
            assert!(
                translation.iter().chain(&rotation).all(|v| v.abs() < 1e-7),
                "{translation:?} {rotation:?}"
            );
        }
        std::fs::remove_file(path).unwrap();
    }
}

fn saddle() -> input::Fragment {
    let mut h = Array2::zeros((6, 6));
    h[(0, 0)] = -0.3;
    h[(3, 3)] = -0.3;
    h[(0, 3)] = 0.3;
    h[(3, 0)] = 0.3;
    input::from_hessian(
        vec!["H".into(), "H".into()],
        vec![1.008; 2],
        vec![0., 0., 0., 1.4, 0., 0.],
        h,
        "saddle".into(),
    )
    .unwrap()
}

#[test]
fn ts_mode_is_retained_after_eckart_and_samples_both_sides_by_default() {
    let f = saddle();
    assert!(f.frequencies.is_empty());
    assert_eq!(f.imaginary_cm1.len(), 1);
    assert!(f.imaginary_cm1[0] < 0.);
    let settings = ts_mode::Settings::default();
    assert!(settings.enabled);
    assert_eq!(settings.side, ts_mode::Side::Both);
    for u in [0., 0.25, 0.5, 0.75, 1.] {
        let mut x = f.equilibrium.clone();
        ts_mode::apply(&mut x, &f, &settings, u).unwrap();
        let projection: f64 = x
            .iter()
            .zip(&f.equilibrium)
            .enumerate()
            .map(|(i, (a, b))| (a - b) * f.imaginary_modes[(i, 0)])
            .sum();
        if u != 0.5 {
            assert_eq!(projection.signum(), (2. * u - 1.).signum());
        }
        let maximum = x
            .chunks_exact(3)
            .zip(f.equilibrium.chunks_exact(3))
            .map(|(a, b)| (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>().sqrt() * input::BOHR)
            .fold(0., f64::max);
        assert!((maximum - settings.max_displacement_angstrom * (2. * u - 1.).abs()).abs() < 1e-12);
    }
}

#[test]
fn ts_one_sided_settings_and_disable_are_honoured_in_generated_trajectories() {
    use std::sync::atomic::{AtomicBool, AtomicUsize};
    for side in [
        ts_mode::Side::Both,
        ts_mode::Side::Positive,
        ts_mode::Side::Negative,
    ] {
        let f = saddle();
        let mut settings = Settings::default();
        settings.points = 50;
        settings.overlap.enabled = false;
        settings.partners[0].fragment = Some(f.clone());
        settings.partners[0].ts_mode.side = side;
        let out =
            super::sampling::generate(settings, 42, &AtomicBool::new(false), &AtomicUsize::new(0))
                .unwrap();
        let dots: Vec<f64> = out
            .samples
            .iter()
            .map(|s| {
                s.coordinates
                    .iter()
                    .enumerate()
                    .map(|(i, x)| (x / input::BOHR - f.equilibrium[i]) * f.imaginary_modes[(i, 0)])
                    .sum()
            })
            .collect();
        match side {
            ts_mode::Side::Both => {
                assert!(dots.iter().any(|v| *v > 0.) && dots.iter().any(|v| *v < 0.))
            }
            ts_mode::Side::Positive => assert!(dots.iter().all(|v| *v >= 0.)),
            ts_mode::Side::Negative => assert!(dots.iter().all(|v| *v <= 0.)),
        }
    }
    let f = saddle();
    let mut x = f.equilibrium.clone();
    let mut settings = ts_mode::Settings::default();
    settings.enabled = false;
    ts_mode::apply(&mut x, &f, &settings, 0.2).unwrap();
    assert_eq!(x, f.equilibrium);
    settings.enabled = true;
    settings.max_displacement_angstrom = f64::NAN;
    assert!(ts_mode::validate(&f, &settings).is_err());
}
