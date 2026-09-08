//! Dense linear solve, imported from Behemoth's `numeric::linalg`.
//!
//! Behemoth's linear algebra sits on `ndarray-linalg`, which needs a system BLAS. These two
//! functions do not: they are plain Gauss-Jordan with partial pivoting, so they come across
//! verbatim and keep Beavyr free of that dependency chain. They complete the surface that the
//! imported optimizer modules need, alongside `normalmode::linalg_shim`.

use anyhow::{anyhow, Result};
use ndarray::{Array1, Array2};

pub fn solve_linear_system_array2(a: &Array2<f64>, b: &Array1<f64>) -> Result<Array1<f64>> {
    let (n, m) = a.dim();
    if n != m {
        return Err(anyhow!(
            "solve_linear_system_array2: A must be square, got {n}x{m}"
        ));
    }
    if b.len() != n {
        return Err(anyhow!(
            "solve_linear_system_array2: dimension mismatch: A is {n}x{n}, b has len {}",
            b.len()
        ));
    }
    // Reuse the existing robust Gauss-Jordan implementation.
    let mut a_vec = array2_to_vec2(a);
    let mut b_vec = b.to_vec();
    let x = gauss_jordan_solve(&mut a_vec, &mut b_vec)?;
    Ok(Array1::from(x))
}

/// Solve A x = b using Gauss-Jordan elimination with partial pivoting.
/// Mutates `a` and `b` and returns `x` on success.
pub fn gauss_jordan_solve(a: &mut [Vec<f64>], b: &mut [f64]) -> Result<Vec<f64>> {
    let n = b.len();
    if n == 0 {
        return Ok(vec![]);
    }
    if a.len() != n || a.iter().any(|r| r.len() != n) {
        return Err(anyhow!("gauss_jordan_solve: A must be nxn and match b"));
    }

    for i in 0..n {
        let mut pivot = i;
        let mut max = a[i][i].abs();
        for r in (i + 1)..n {
            let v = a[r][i].abs();
            if v > max {
                max = v;
                pivot = r;
            }
        }
        if max < 1e-14 {
            return Err(anyhow!("gauss_jordan_solve: linear system is singular"));
        }
        if pivot != i {
            a.swap(i, pivot);
            b.swap(i, pivot);
        }

        let diag = a[i][i];
        for j in i..n {
            a[i][j] /= diag;
        }
        b[i] /= diag;

        for r in 0..n {
            if r == i {
                continue;
            }
            let factor = a[r][i];
            for c in i..n {
                a[r][c] -= factor * a[i][c];
            }
            b[r] -= factor * b[i];
        }
    }
    Ok(b.to_vec())
}

fn array2_to_vec2(a: &Array2<f64>) -> Vec<Vec<f64>> {
    let (n, m) = a.dim();
    let mut out = vec![vec![0.0f64; m]; n];
    for i in 0..n {
        for j in 0..m {
            out[i][j] = a[(i, j)];
        }
    }
    out
}
