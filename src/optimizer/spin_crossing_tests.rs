//! Numerical reference from Peter's SpinCross_2022_12_09a.f90 (gfortran).
use super::*;

#[test]
fn mecp_gradient_and_bfgs_match_original_fortran() {
    let (parallel, perpendicular, effective) = effective_gradient(
        -1.234,
        -1.235,
        &[0.1, -0.2, 0.3],
        &[-0.2, 0.05, 0.15],
        1.,
        140.,
    )
    .unwrap();
    for (actual, expected) in [
        (
            parallel,
            vec![
                -0.1142857142857143,
                -0.021428571428571463,
                0.19285714285714287,
            ],
        ),
        (
            perpendicular,
            vec![
                0.00030000000000003364,
                -0.00025000000000002798,
                0.00015000000000001679,
            ],
        ),
        (
            effective,
            vec![
                -0.0722857142857096,
                -0.05642857142857538,
                0.21385714285714522,
            ],
        ),
    ] {
        for (a, b) in actual.iter().zip(expected) {
            assert!((a - b).abs() < 1e-13);
        }
    }
    let (h, updated) = bfgs_inverse_update(
        &(Array2::eye(3) * INITIAL_INVERSE_HESSIAN),
        &[0.02, -0.01, 0.03],
        &[0.03, -0.01, 0.06],
    );
    assert!(updated);
    for (a, b) in h.iter().zip([
        0.73008, -0.04304, -0.03888, -0.04304, 0.73552, -0.02256, -0.03888, -0.02256, 0.51568,
    ]) {
        assert!((a - b).abs() < 1e-13);
    }
}

#[test]
fn mecp_initial_inverse_hessian_is_point_seven_in_angstrom_coordinates() {
    struct Surface(f64);
    impl Objective for Surface {
        fn energy(&mut self, x: &[f64]) -> Result<f64> {
            Ok(0.5 * (x[0] - self.0).powi(2) + 0.5 * x[1] * x[1] + 0.5 * x[2] * x[2])
        }
        fn gradient(&mut self, x: &[f64], g: &mut [f64]) -> Result<()> {
            g.copy_from_slice(x);
            g[0] -= self.0;
            Ok(())
        }
    }
    assert_eq!(SpinCrossingOptions::default().initial_inverse_hessian, 0.7);
    let x = vec![0., 0.02, 0.03];
    let result = optimize_spin_crossing(
        x.clone(),
        None,
        &mut Surface(1.),
        &mut Surface(-1.),
        SpinCrossingOptions {
            coordinate_mode: SpinCrossingCoordinateMode::Cartesian,
            max_cycles: 1,
            trajectory_file: None,
            verbosity: 0,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(result.x[0], 0.);
    for i in [1, 2] {
        let expected_angstrom = x[i] * BOHR_TO_ANGSTROM - 0.7 * x[i] / BOHR_TO_ANGSTROM;
        assert!((result.x[i] * BOHR_TO_ANGSTROM - expected_angstrom).abs() < 1e-14);
    }
}

#[test]
fn mecp_convergence_checks_absolute_gap_and_original_units() {
    let o = SpinCrossingOptions::default();
    assert!(!converged(-0.1, &[0.; 3], &[0.; 3], &o));
    assert!(converged(0., &[0.; 3], &[0.; 3], &o));
    assert!(!converged(
        0.,
        &[0.005 / BOHR_TO_ANGSTROM, 0., 0.],
        &[0.; 3],
        &o
    ));
    assert!(!converged(
        0.,
        &[0.; 3],
        &[0.001 * BOHR_TO_ANGSTROM, 0., 0.],
        &o
    ));
}
