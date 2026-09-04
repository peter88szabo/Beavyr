#![allow(dead_code)]
//! Generic numerical (finite-difference) nuclear Hessian from analytic gradients.
//!
//! This mirrors the established GFN1-xTB finite-difference Hessian
//! (`hamiltonian::xtb_hamiltonian::fd_hessian`) but is not tied to any one
//! Hamiltonian: the caller supplies a closure that returns the flattened
//! Cartesian gradient (Hartree/bohr) at an arbitrary displaced geometry, so
//! any method that already has an analytic gradient (RHF/UHF/ROHF, MP2,
//! RI-MP2, RKS, UKS, ...) can get a Hessian, and therefore frequencies and
//! thermochemistry, without a dedicated analytic second-derivative
//! implementation. This is an O(N) times more expensive stopgap, not a
//! replacement for an eventual analytic Hessian.

use anyhow::{bail, Result};
use ndarray::Array2;

/// Central-difference Hessian built from repeated analytic-gradient evaluations.
///
/// `gradient_fn` must return the flattened (length `coords_bohr.len()`)
/// Cartesian gradient in Hartree/bohr for the displaced geometry it is given,
/// using the same atom ordering and flattening (`x0,y0,z0,x1,y1,z1,...`) as
/// `coords_bohr`. The returned matrix is explicitly symmetrized, since a
/// central difference of an independently converged calculation at each
/// displacement is not exactly symmetric numerically.
pub fn numerical_hessian_from_gradients<F>(
    coords_bohr: &[f64],
    step_bohr: f64,
    verbosity: u8,
    mut gradient_fn: F,
) -> Result<Array2<f64>>
where
    F: FnMut(&[f64]) -> Result<Vec<f64>>,
{
    if step_bohr <= 0.0 {
        bail!("numerical Hessian step must be positive");
    }
    let n = coords_bohr.len();
    let mut hessian = Array2::<f64>::zeros((n, n));

    for column in 0..n {
        if verbosity > 0 {
            println!("Numerical Hessian: displaced coordinate {}/{n}", column + 1);
        }

        let mut plus = coords_bohr.to_vec();
        plus[column] += step_bohr;
        let gradient_plus = gradient_fn(&plus)?;
        if gradient_plus.len() != n {
            bail!(
                "numerical Hessian gradient callback returned {} components, expected {n}",
                gradient_plus.len()
            );
        }

        let mut minus = coords_bohr.to_vec();
        minus[column] -= step_bohr;
        let gradient_minus = gradient_fn(&minus)?;
        if gradient_minus.len() != n {
            bail!(
                "numerical Hessian gradient callback returned {} components, expected {n}",
                gradient_minus.len()
            );
        }

        for row in 0..n {
            hessian[(row, column)] = (gradient_plus[row] - gradient_minus[row]) / (2.0 * step_bohr);
        }
    }

    symmetrize(&mut hessian);
    Ok(hessian)
}

fn symmetrize(hessian: &mut Array2<f64>) {
    let n = hessian.nrows();
    for i in 0..n {
        for j in 0..i {
            let value = 0.5 * (hessian[(i, j)] + hessian[(j, i)]);
            hessian[(i, j)] = value;
            hessian[(j, i)] = value;
        }
    }
}

/// Reports the relative Frobenius asymmetry of a Hessian before symmetrization,
/// a cheap diagnostic for how well the finite-difference step and the
/// underlying gradient's own convergence tolerance resolve the true Hessian.
pub fn hessian_asymmetry_ratio(hessian: &Array2<f64>) -> f64 {
    let (n, m) = hessian.dim();
    if n != m {
        return 0.0;
    }
    let mut diff_sq = 0.0;
    let mut norm_sq = 0.0;
    for i in 0..n {
        for j in 0..n {
            let hij = hessian[(i, j)];
            let hji = hessian[(j, i)];
            let d = hij - hji;
            diff_sq += d * d;
            norm_sq += hij * hij;
        }
    }
    let norm = norm_sq.sqrt();
    if norm > 0.0 {
        diff_sq.sqrt() / norm
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovers_the_hessian_of_a_quadratic_test_function() {
        // f(x) = 0.5 x^T A x with a fixed, hand-picked symmetric positive
        // definite A; the analytic gradient is A x, so a central-difference
        // Hessian of that gradient must reproduce A itself (up to the finite
        // step's O(h^2) truncation error).
        let a = Array2::from_shape_vec((3, 3), vec![4.0, 1.0, 0.2, 1.0, 3.0, 0.5, 0.2, 0.5, 2.0])
            .unwrap();
        let x0 = [0.3, -0.2, 0.1];
        let gradient_fn = |x: &[f64]| -> Result<Vec<f64>> {
            let mut g = vec![0.0; 3];
            for i in 0..3 {
                for j in 0..3 {
                    g[i] += a[(i, j)] * x[j];
                }
            }
            Ok(g)
        };
        let hessian = numerical_hessian_from_gradients(&x0, 1.0e-4, 0, gradient_fn).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (hessian[(i, j)] - a[(i, j)]).abs() < 1.0e-6,
                    "hessian[{i},{j}] = {} does not match A[{i},{j}] = {}",
                    hessian[(i, j)],
                    a[(i, j)]
                );
            }
        }
    }

    #[test]
    fn rejects_a_gradient_callback_with_the_wrong_length() {
        let x0 = [0.0, 0.0];
        let gradient_fn = |_: &[f64]| -> Result<Vec<f64>> { Ok(vec![0.0]) };
        assert!(numerical_hessian_from_gradients(&x0, 1.0e-3, 0, gradient_fn).is_err());
    }
}
