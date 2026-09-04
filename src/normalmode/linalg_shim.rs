#![allow(dead_code)]
//! The small slice of Behemoth's `numeric::linalg` that its `normalmode`
//! code actually uses, backed by the pure-Rust Jacobi eigensolver rather
//! than BLAS/LAPACK.
//!
//! Behemoth's own `LinAlg` sits on `ndarray-linalg` with OpenBLAS or MKL --
//! a system-native dependency chain that a molecular viewer has no business
//! requiring to build. Its `eigh` already falls back to
//! `numeric::jacobi_diag::jacobi` when no backend is selected, and that
//! solver is dependency-free, so the imported analysis code runs unmodified
//! on top of this shim.
//!
//! The surface here is exactly what the imported files call -- `new`,
//! `eigh`, `matmul`, `transpose` -- and the bodies are Behemoth's own
//! `transpose_pure`/`matmul_pure`/`jacobi_eigh`, carried over verbatim.

use anyhow::{anyhow, Result};
use ndarray::{Array1, Array2};

use crate::normalmode::jacobi_diag::jacobi;

/// Which numerical backend to use. Only `Auto` exists here: unlike
/// Behemoth, there is no BLAS to choose. The variant is kept so the
/// imported `LinAlg::new(Backend::Auto)` call sites need no edit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    Auto,
}

pub struct LinAlg {
    _backend: Backend,
}

impl LinAlg {
    pub fn new(choice: Backend) -> Result<Self> {
        Ok(Self { _backend: choice })
    }

    pub fn transpose(&self, a: &Array2<f64>) -> Array2<f64> {
        let (m, n) = a.dim();
        let mut out = Array2::<f64>::zeros((n, m));
        for i in 0..m {
            for j in 0..n {
                out[(j, i)] = a[(i, j)];
            }
        }
        out
    }

    pub fn matmul(&self, a: &Array2<f64>, b: &Array2<f64>) -> Array2<f64> {
        let (m, k1) = a.dim();
        let (k2, n) = b.dim();
        assert_eq!(k1, k2, "matmul: inner dimensions must agree");

        // Pre-transpose B for cache locality, exactly as Behemoth's
        // `matmul_pure` does.
        let bt = self.transpose(b);
        let mut c = Array2::<f64>::zeros((m, n));
        for i in 0..m {
            for j in 0..n {
                let mut sum = 0.0f64;
                for p in 0..k1 {
                    sum += a[(i, p)] * bt[(j, p)];
                }
                c[(i, j)] = sum;
            }
        }
        c
    }

    /// Symmetric eigendecomposition, ascending eigenvalues.
    pub fn eigh(&self, a: &Array2<f64>) -> Result<(Array1<f64>, Array2<f64>)> {
        let (nrows, ncols) = a.dim();
        if nrows != ncols {
            return Err(anyhow!("eigh: matrix must be square"));
        }
        let n = nrows;

        let mut mat: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                mat[i][j] = a[(i, j)];
            }
        }

        let itmax = usize::max(100, 10 * n);
        let tol = 1e-10;
        let iord = 1; // +1 ascending

        let (v_mat, w_vec) = jacobi(&mat, itmax, tol, iord);

        let w = Array1::from_vec(w_vec);
        let mut v = Array2::<f64>::zeros((n, n));
        for i in 0..n {
            for j in 0..n {
                v[(i, j)] = v_mat[i][j];
            }
        }
        Ok((w, v))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eigh_recovers_a_known_spectrum() {
        // Diagonal matrix: eigenvalues are the diagonal, ascending.
        let mut a = Array2::<f64>::zeros((3, 3));
        a[(0, 0)] = 3.0;
        a[(1, 1)] = -1.0;
        a[(2, 2)] = 2.0;
        let la = LinAlg::new(Backend::Auto).unwrap();
        let (w, _v) = la.eigh(&a).unwrap();
        assert!((w[0] - -1.0).abs() < 1e-9, "{w:?}");
        assert!((w[1] - 2.0).abs() < 1e-9, "{w:?}");
        assert!((w[2] - 3.0).abs() < 1e-9, "{w:?}");
    }

    #[test]
    fn eigh_reconstructs_the_matrix_from_its_decomposition() {
        // A symmetric matrix must satisfy A ~= V diag(w) V^T.
        let mut a = Array2::<f64>::zeros((3, 3));
        let vals = [[4.0, 1.0, 0.5], [1.0, 3.0, -0.2], [0.5, -0.2, 2.0]];
        for i in 0..3 {
            for j in 0..3 {
                a[(i, j)] = vals[i][j];
            }
        }
        let la = LinAlg::new(Backend::Auto).unwrap();
        let (w, v) = la.eigh(&a).unwrap();

        let mut lambda = Array2::<f64>::zeros((3, 3));
        for i in 0..3 {
            lambda[(i, i)] = w[i];
        }
        let vt = la.transpose(&v);
        let reconstructed = la.matmul(&la.matmul(&v, &lambda), &vt);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (reconstructed[(i, j)] - a[(i, j)]).abs() < 1e-8,
                    "({i},{j}): {} vs {}",
                    reconstructed[(i, j)],
                    a[(i, j)]
                );
            }
        }
    }

    #[test]
    fn matmul_matches_a_hand_computed_product() {
        let mut a = Array2::<f64>::zeros((2, 3));
        let mut b = Array2::<f64>::zeros((3, 2));
        for i in 0..2 {
            for j in 0..3 {
                a[(i, j)] = (i * 3 + j) as f64 + 1.0; // 1..6
            }
        }
        for i in 0..3 {
            for j in 0..2 {
                b[(i, j)] = (i * 2 + j) as f64 + 1.0; // 1..6
            }
        }
        let la = LinAlg::new(Backend::Auto).unwrap();
        let c = la.matmul(&a, &b);
        // [[1,2,3],[4,5,6]] * [[1,2],[3,4],[5,6]] = [[22,28],[49,64]]
        assert_eq!(c[(0, 0)], 22.0);
        assert_eq!(c[(0, 1)], 28.0);
        assert_eq!(c[(1, 0)], 49.0);
        assert_eq!(c[(1, 1)], 64.0);
    }

    #[test]
    fn transpose_swaps_indices() {
        let mut a = Array2::<f64>::zeros((2, 3));
        a[(0, 2)] = 7.0;
        let la = LinAlg::new(Backend::Auto).unwrap();
        let t = la.transpose(&a);
        assert_eq!(t.dim(), (3, 2));
        assert_eq!(t[(2, 0)], 7.0);
    }
}
