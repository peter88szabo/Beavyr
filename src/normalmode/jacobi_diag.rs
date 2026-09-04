#![allow(dead_code, unused_mut)]
pub fn jacobi(
    a: &Vec<Vec<f64>>,
    itmax: usize,
    tol: f64,
    iord: i32, // +1 ascending, -1 descending
) -> (Vec<Vec<f64>>, Vec<f64>) {
    let n = a.len();
    assert!(a.iter().all(|row| row.len() == n), "Matrix must be square");
    if n == 0 {
        return (vec![], vec![]);
    }

    // Work on a copy so input isn't destroyed
    let mut a = a.clone();

    // Eigenvectors start as I
    let mut v = vec![vec![0.0; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }

    // Diagonal as initial eigenvalue estimates
    let mut d = (0..n).map(|i| a[i][i]).collect::<Vec<_>>();

    // Helper: off-diagonal norm (Frobenius, but only lower triangle counted once)
    let mut off = |mat: &Vec<Vec<f64>>| -> f64 {
        let mut s = 0.0;
        for i in 0..n {
            for j in 0..i {
                let x = mat[i][j];
                s += x * x;
            }
        }
        s.sqrt()
    };

    // Convergence target
    let mut offn = off(&a);
    if offn <= tol {
        // NOTE (Beavyr deviation from the imported Behemoth original):
        // upstream returned here directly, which skips the `iord` sort at
        // the end of this function -- so an already-diagonal input came
        // back in its original diagonal order regardless of the requested
        // ordering, silently breaking this function's own documented
        // contract. `normal_modes` depends on ascending order to tell
        // translations/rotations from vibrations, so the sort is applied on
        // this path too.
        match iord {
            1 => sort_eigs(&mut d, &mut v, true),
            -1 => sort_eigs(&mut d, &mut v, false),
            _ => {}
        }
        return (v, d);
    }

    let mut iter = 0usize;

    while iter < itmax {
        let mut changed = false;

        // One **cyclic sweep** over all i>j
        for i in 1..n {
            for j in 0..i {
                let aij = a[i][j];
                if aij.abs() <= tol {
                    continue;
                }

                let aii = d[i];
                let ajj = d[j];

                // Stable Jacobi rotation (Golub–Van Loan)
                // tau = (ajj - aii) / (2 * aij)
                let tau = (ajj - aii) / (2.0 * aij);
                let t = {
                    // t = sign(tau) / (|tau| + sqrt(1 + tau^2))
                    // This minimizes cancellation and keeps |t| <= 1
                    let sign = if tau >= 0.0 { 1.0 } else { -1.0 };
                    sign / (tau.abs() + (1.0 + tau * tau).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;

                // Update diagonal eigenvalue estimates
                let aii_new = c * c * aii - 2.0 * s * c * aij + s * s * ajj;
                let ajj_new = s * s * aii + 2.0 * s * c * aij + c * c * ajj;
                d[i] = aii_new;
                d[j] = ajj_new;
                a[i][j] = 0.0;
                a[j][i] = 0.0;

                // Update off-diagonal rows/cols (preserve symmetry)
                for k in 0..n {
                    if k != i && k != j {
                        // rows k,i and k,j
                        let aki = a[k][i];
                        let akj = a[k][j];
                        let new_ki = c * aki - s * akj;
                        let new_kj = s * aki + c * akj;
                        a[k][i] = new_ki;
                        a[i][k] = new_ki;
                        a[k][j] = new_kj;
                        a[j][k] = new_kj;
                    }
                }

                // Update eigenvectors V = V * G(i,j)
                for k in 0..n {
                    let vki = v[k][i];
                    let vkj = v[k][j];
                    v[k][i] = c * vki - s * vkj;
                    v[k][j] = s * vki + c * vkj;
                }

                changed = true;
            }
        }

        iter += 1;

        // Check convergence every sweep
        offn = off(&a);
        if offn <= tol || !changed {
            break;
        }
    }

    if iter == itmax {
        eprintln!(
            "Warning: Jacobi reached itmax={} (offdiag norm = {:.3e})",
            itmax, offn
        );
    }

    // Sort eigenpairs if requested
    match iord {
        1 => sort_eigs(&mut d, &mut v, true),   // ascending
        -1 => sort_eigs(&mut d, &mut v, false), // descending
        _ => {}
    }

    // Orthonormality drift is usually tiny; optionally re-orthonormalize here.

    (v, d)
}

#[inline]
fn sort_eigs(vals: &mut Vec<f64>, vecs: &mut Vec<Vec<f64>>, ascending: bool) {
    let n = vals.len();
    // Build permutation of indices
    let mut idx = (0..n).collect::<Vec<_>>();
    if ascending {
        idx.sort_by(|&i, &j| vals[i].partial_cmp(&vals[j]).unwrap());
    } else {
        idx.sort_by(|&i, &j| vals[j].partial_cmp(&vals[i]).unwrap());
    }
    // Apply permutation: vals and columns of vecs
    let vals_old = vals.clone();
    let vecs_old = vecs.clone();
    for (k, &p) in idx.iter().enumerate() {
        vals[k] = vals_old[p];
        for r in 0..n {
            vecs[r][k] = vecs_old[r][p];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Regression test for the deviation noted above: an already-diagonal
    /// matrix takes the early-return path, which upstream left unsorted.
    #[test]
    fn an_already_diagonal_matrix_is_still_sorted_as_requested() {
        let a = vec![
            vec![3.0, 0.0, 0.0],
            vec![0.0, -1.0, 0.0],
            vec![0.0, 0.0, 2.0],
        ];
        let (_v, d) = jacobi(&a, 100, 1e-10, 1);
        assert_eq!(d, vec![-1.0, 2.0, 3.0], "ascending order requested");

        let (_v, d) = jacobi(&a, 100, 1e-10, -1);
        assert_eq!(d, vec![3.0, 2.0, -1.0], "descending order requested");
    }

    /// The eigenvectors must be permuted with their eigenvalues, not left
    /// behind -- otherwise sorting silently mismatches pairs.
    #[test]
    fn sorting_keeps_eigenvectors_with_their_eigenvalues() {
        let a = vec![
            vec![3.0, 0.0, 0.0],
            vec![0.0, -1.0, 0.0],
            vec![0.0, 0.0, 2.0],
        ];
        let (v, d) = jacobi(&a, 100, 1e-10, 1);
        // Eigenvalue -1.0 (now first) belonged to the second coordinate, so
        // its eigenvector must be the second unit vector.
        assert_eq!(d[0], -1.0);
        assert!((v[1][0].abs() - 1.0).abs() < 1e-12, "{v:?}");
        assert!(v[0][0].abs() < 1e-12);
        assert!(v[2][0].abs() < 1e-12);
    }
}
