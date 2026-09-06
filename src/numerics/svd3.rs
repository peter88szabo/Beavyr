//! Symmetric eigendecomposition and singular value decomposition for 3x3
//! real matrices, with no BLAS/LAPACK dependency.
//!
//! Only the fixed 3x3 case is handled, which is what makes this practical:
//! a cyclic Jacobi sweep on a 3x3 symmetric matrix converges in a handful of
//! rotations and needs no pivoting strategy or workspace queries. The SVD is
//! then built from the eigendecomposition of `A^T A` in the standard way,
//! which is accurate enough for the superposition problem this exists for
//! (Kabsch alignment) while staying entirely self-contained.

use super::mat3::{mat_vec, matmul, transpose, zeros, Mat3};

/// Eigendecomposition of a real symmetric 3x3 matrix.
///
/// Returns eigenvalues in descending order together with the matching
/// eigenvectors as the *columns* of `vectors`. The input is assumed
/// symmetric; only its lower/upper agreement matters, not which half is
/// read, since a Jacobi sweep touches both.
pub struct SymEigen3 {
    pub values: [f64; 3],
    pub vectors: Mat3,
}

/// Jacobi eigenvalue iteration for a symmetric 3x3 matrix.
pub fn symmetric_eigen3(input: &Mat3) -> SymEigen3 {
    // Work on a mutable copy; `v` accumulates the rotations.
    let mut a = *input;
    let mut v = super::mat3::IDENTITY;

    // A 3x3 symmetric matrix has only three off-diagonal pairs, so a small
    // fixed sweep budget is ample; the loop exits as soon as the
    // off-diagonal norm is negligible.
    const MAX_SWEEPS: usize = 64;
    const TOL: f64 = 1.0e-14;
    for _ in 0..MAX_SWEEPS {
        let off = a[0][1].abs() + a[0][2].abs() + a[1][2].abs();
        if off <= TOL {
            break;
        }
        for (p, q) in [(0usize, 1usize), (0, 2), (1, 2)] {
            let apq = a[p][q];
            if apq.abs() <= TOL {
                continue;
            }
            // Standard stable Jacobi rotation angle (Golub & Van Loan).
            let theta = (a[q][q] - a[p][p]) / (2.0 * apq);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;

            let mut rot = super::mat3::IDENTITY;
            rot[p][p] = c;
            rot[q][q] = c;
            rot[p][q] = s;
            rot[q][p] = -s;

            // a <- rot^T a rot, keeping it symmetric by construction.
            a = matmul(&transpose(&rot), &matmul(&a, &rot));
            v = matmul(&v, &rot);
        }
    }

    let mut values = [a[0][0], a[1][1], a[2][2]];
    let mut vectors = v;

    // Sort descending, carrying each eigenvector column with its value.
    let mut order = [0usize, 1, 2];
    order.sort_by(|&i, &j| values[j].total_cmp(&values[i]));
    let sorted_values = [values[order[0]], values[order[1]], values[order[2]]];
    let mut sorted_vectors = zeros();
    for (new_col, &old_col) in order.iter().enumerate() {
        for row in 0..3 {
            sorted_vectors[row][new_col] = vectors[row][old_col];
        }
    }
    values = sorted_values;
    vectors = sorted_vectors;

    SymEigen3 { values, vectors }
}

/// Singular value decomposition `A = U * diag(S) * V^T` for a 3x3 matrix,
/// with singular values in descending order.
pub struct Svd3 {
    pub u: Mat3,
    /// Singular values, descending. Kabsch superposition only needs `u` and
    /// `v`, but the singular values are part of what an SVD *is* -- this is
    /// a general-purpose numerics module, and dropping them would make it
    /// useless for rank/conditioning questions later.
    #[allow(dead_code)]
    pub s: [f64; 3],
    pub v: Mat3,
}

/// Computes the SVD of a 3x3 matrix via the eigendecomposition of `A^T A`.
///
/// `V` comes from the eigenvectors of `A^T A` and the singular values from
/// the square roots of its eigenvalues. Each column of `U` is then
/// `A * v_i / s_i`; for a (near-)zero singular value that division is
/// undefined, so the corresponding `U` column is instead completed to an
/// orthonormal basis via a cross product, which keeps `U` a proper rotation
/// even for rank-deficient input (a planar or linear molecule, in the
/// superposition case).
pub fn svd3(a: &Mat3) -> Svd3 {
    let ata = matmul(&transpose(a), a);
    let eig = symmetric_eigen3(&ata);

    let mut s = [0.0f64; 3];
    for i in 0..3 {
        // Clamp tiny negatives that come from round-off on a semidefinite
        // matrix; a singular value is never imaginary.
        s[i] = eig.values[i].max(0.0).sqrt();
    }

    let v = eig.vectors;
    let mut u = zeros();
    const SMALL: f64 = 1.0e-12;
    let mut valid = [false; 3];
    for col in 0..3 {
        if s[col] <= SMALL {
            continue;
        }
        let v_col = [v[0][col], v[1][col], v[2][col]];
        let av = mat_vec(a, v_col);
        for row in 0..3 {
            u[row][col] = av[row] / s[col];
        }
        valid[col] = true;
    }

    // Fill any degenerate columns with something orthonormal so `U` stays a
    // full orthogonal matrix rather than containing zero columns.
    for col in 0..3 {
        if valid[col] {
            continue;
        }
        let (a_idx, b_idx) = match col {
            0 => (1, 2),
            1 => (2, 0),
            _ => (0, 1),
        };
        if valid[a_idx] && valid[b_idx] {
            let x = [u[0][a_idx], u[1][a_idx], u[2][a_idx]];
            let y = [u[0][b_idx], u[1][b_idx], u[2][b_idx]];
            let cross = [
                x[1] * y[2] - x[2] * y[1],
                x[2] * y[0] - x[0] * y[2],
                x[0] * y[1] - x[1] * y[0],
            ];
            for row in 0..3 {
                u[row][col] = cross[row];
            }
        } else {
            // Two or more degenerate directions: fall back to a unit axis,
            // which is as good as any other choice for a rank <= 1 input.
            u[col][col] = 1.0;
        }
        valid[col] = true;
    }

    Svd3 { u, s, v }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numerics::mat3::{determinant, IDENTITY};

    fn approx_eq_mat(a: &Mat3, b: &Mat3, tol: f64) -> bool {
        (0..3).all(|i| (0..3).all(|j| (a[i][j] - b[i][j]).abs() < tol))
    }

    #[test]
    fn eigen_of_a_diagonal_matrix_is_its_diagonal_descending() {
        let d: Mat3 = [[1.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 3.0]];
        let e = symmetric_eigen3(&d);
        assert!((e.values[0] - 5.0).abs() < 1e-12, "{:?}", e.values);
        assert!((e.values[1] - 3.0).abs() < 1e-12, "{:?}", e.values);
        assert!((e.values[2] - 1.0).abs() < 1e-12, "{:?}", e.values);
    }

    #[test]
    fn eigen_reconstructs_the_original_matrix() {
        // A = V diag(w) V^T must hold for a symmetric input.
        let a: Mat3 = [[4.0, 1.0, 0.5], [1.0, 3.0, -0.2], [0.5, -0.2, 2.0]];
        let e = symmetric_eigen3(&a);
        let mut lambda = zeros();
        for i in 0..3 {
            lambda[i][i] = e.values[i];
        }
        let recon = matmul(&matmul(&e.vectors, &lambda), &transpose(&e.vectors));
        assert!(approx_eq_mat(&recon, &a, 1e-9), "{recon:?} vs {a:?}");
    }

    #[test]
    fn eigenvectors_are_orthonormal() {
        let a: Mat3 = [[4.0, 1.0, 0.5], [1.0, 3.0, -0.2], [0.5, -0.2, 2.0]];
        let e = symmetric_eigen3(&a);
        let should_be_identity = matmul(&transpose(&e.vectors), &e.vectors);
        assert!(approx_eq_mat(&should_be_identity, &IDENTITY, 1e-9));
    }

    #[test]
    fn svd_reconstructs_a_general_matrix() {
        let a: Mat3 = [[1.0, 2.0, 0.0], [0.0, 1.0, 3.0], [4.0, 0.0, 1.0]];
        let d = svd3(&a);
        let mut sigma = zeros();
        for i in 0..3 {
            sigma[i][i] = d.s[i];
        }
        let recon = matmul(&matmul(&d.u, &sigma), &transpose(&d.v));
        assert!(approx_eq_mat(&recon, &a, 1e-8), "{recon:?} vs {a:?}");
    }

    #[test]
    fn svd_singular_values_are_descending_and_nonnegative() {
        let a: Mat3 = [[3.0, 1.0, -2.0], [0.5, 2.0, 1.0], [1.0, -1.0, 4.0]];
        let d = svd3(&a);
        assert!(d.s[0] >= d.s[1] && d.s[1] >= d.s[2]);
        assert!(d.s.iter().all(|&v| v >= 0.0));
    }

    #[test]
    fn svd_factors_are_orthogonal() {
        let a: Mat3 = [[3.0, 1.0, -2.0], [0.5, 2.0, 1.0], [1.0, -1.0, 4.0]];
        let d = svd3(&a);
        assert!(approx_eq_mat(&matmul(&transpose(&d.u), &d.u), &IDENTITY, 1e-8));
        assert!(approx_eq_mat(&matmul(&transpose(&d.v), &d.v), &IDENTITY, 1e-8));
    }

    #[test]
    fn svd_of_a_rank_deficient_matrix_still_yields_an_orthogonal_u() {
        // A planar (rank-2) case: the third singular value is zero, which is
        // exactly where naive `A v / s` would divide by zero.
        let a: Mat3 = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.0]];
        let d = svd3(&a);
        assert!(d.s[2].abs() < 1e-12, "{:?}", d.s);
        assert!(approx_eq_mat(&matmul(&transpose(&d.u), &d.u), &IDENTITY, 1e-8));
        assert!(determinant(&d.u).abs() > 0.5, "U must stay non-singular");
    }

    #[test]
    fn svd_of_the_identity_is_all_unit_singular_values() {
        let d = svd3(&IDENTITY);
        for v in d.s {
            assert!((v - 1.0).abs() < 1e-12, "{:?}", d.s);
        }
    }
}
