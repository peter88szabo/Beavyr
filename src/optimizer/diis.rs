//! Pulay's DIIS, imported from Behemoth's `solver::diis`.
//!
//! Behemoth uses this to accelerate SCF convergence, extrapolating a Fock matrix from a history of
//! error vectors. Nothing about it is specific to that: it takes (vector, error) pairs, solves the
//! standard Pulay equations for coefficients summing to one, and returns the extrapolated vector.
//!
//! That makes it directly reusable for **GDIIS** -- geometry DIIS, in the Császár-Pulay and
//! Farkas-Schlegel sense. Pushing `(xᵢ + eᵢ, eᵢ)`, where `eᵢ = -H⁻¹gᵢ` is the quasi-Newton
//! displacement at geometry `xᵢ`, makes [`Diis::extrapolate`] return `Σ cᵢ(xᵢ + eᵢ)`, which is
//! exactly the GDIIS geometry. See [`crate::optimizer::minimum_cartesian`].
//!
//! Changed from upstream: it takes `gauss_jordan_solve` from
//! [`crate::optimizer::linalg_solve`] rather than Behemoth's `numeric::linalg`, and carries its
//! own two-line `dot`. The algorithm is untouched.

use anyhow::{anyhow, Result};

use crate::optimizer::linalg_solve::gauss_jordan_solve;

pub struct Diis {
    max_vecs: usize,
    vecs: Vec<Vec<f64>>,
    errs: Vec<Vec<f64>>,
}

//Pulay's Direct Iterative Inversion of the Subspace (DIIS) to
//accelerate the SCF convergence
impl Diis {
    pub fn new(max_vecs: usize) -> Self {
        Self {
            max_vecs: max_vecs.max(2),
            vecs: Vec::new(),
            errs: Vec::new(),
        }
    }

    pub fn push(&mut self, vec: Vec<f64>, err: Vec<f64>) {
        if !self.vecs.is_empty() && vec.len() != self.vecs[0].len() {
            return;
        }
        if !self.errs.is_empty() && err.len() != self.errs[0].len() {
            return;
        }
        self.vecs.push(vec);
        self.errs.push(err);
        if self.vecs.len() > self.max_vecs {
            self.vecs.remove(0);
            self.errs.remove(0);
        }
    }

    /// Removes every stored vector and error so a guarded SCF driver can
    /// recover from a deteriorating or nearly singular extrapolation space.
    pub fn clear(&mut self) {
        self.vecs.clear();
        self.errs.clear();
    }

    /// Number of Fock/error pairs currently retained.
    pub fn len(&self) -> usize {
        self.vecs.len()
    }

    pub fn extrapolate(&self) -> Option<Vec<f64>> {
        self.extrapolate_with_damping(0.0)
    }

    pub fn extrapolate_with_damping(&self, damping: f64) -> Option<Vec<f64>> {
        let m = self.vecs.len();
        if m < 2 {
            return None;
        }
        let coeffs = solve_diis_coeffs(&self.errs).ok()?;
        let mut out = vec![0.0f64; self.vecs[0].len()];
        for (c, v) in coeffs.iter().zip(self.vecs.iter()) {
            for i in 0..out.len() {
                out[i] += c * v[i];
            }
        }
        let damp = damping.clamp(0.0, 1.0);
        if damp > 0.0 {
            let last = self.vecs.last().unwrap();
            for i in 0..out.len() {
                out[i] = (1.0 - damp) * out[i] + damp * last[i];
            }
        }
        Some(out)
    }
}

impl Diis {
    /// Extrapolates, shrinking the history from the oldest end until the fit is solvable.
    ///
    /// **Added here; not in Behemoth's original.** Its SCF driver never needs this: an SCF's error
    /// vectors stay well separated, so the Pulay system stays well conditioned and
    /// [`Diis::extrapolate`] simply works.
    ///
    /// Geometry is different. As an optimisation converges the quasi-Newton displacements become
    /// nearly parallel and shrink together, the `B` matrix goes singular, `gauss_jordan_solve`
    /// reports it, and `extrapolate` returns `None` -- silently, so GDIIS degrades to plain BFGS
    /// exactly when it should be helping most. That was measured: with three, four and five
    /// geometries in hand the fit failed every time.
    ///
    /// Farkas and Schlegel's remedy is this one: discard the oldest vectors and retry, since the
    /// recent ones carry the information that matters and the stale ones are what make the space
    /// degenerate. Returns the extrapolation from the largest usable window, or `None` if even two
    /// vectors cannot be fitted.
    pub fn extrapolate_shrinking(&self) -> Option<Vec<f64>> {
        let available = self.vecs.len();
        for keep in (2..=available).rev() {
            let first = available - keep;
            let errs = &self.errs[first..];
            let Ok(coeffs) = solve_diis_coeffs(errs) else {
                continue;
            };
            // A fit whose coefficients have run away is numerically worthless even when the solve
            // reported success, and is the other face of the same ill-conditioning.
            const MAX_COEFFICIENT: f64 = 1.0e3;
            if coeffs.iter().any(|c| !c.is_finite() || c.abs() > MAX_COEFFICIENT) {
                continue;
            }
            let mut out = vec![0.0f64; self.vecs[0].len()];
            for (c, v) in coeffs.iter().zip(self.vecs[first..].iter()) {
                for i in 0..out.len() {
                    out[i] += c * v[i];
                }
            }
            return Some(out);
        }
        None
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn solve_diis_coeffs(errs: &[Vec<f64>]) -> Result<Vec<f64>> {
    let m = errs.len();
    let mut b = vec![vec![0.0f64; m + 1]; m + 1];
    for i in 0..m {
        for j in 0..m {
            b[i][j] = dot(&errs[i], &errs[j]);
        }
        b[i][m] = -1.0;
        b[m][i] = -1.0;
    }
    let mut rhs = vec![0.0f64; m + 1];
    rhs[m] = -1.0;

    let sol = gauss_jordan_solve(&mut b, &mut rhs)?;
    Ok(sol[..m].to_vec())
}
