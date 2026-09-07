#![allow(dead_code)]
use anyhow::{anyhow, bail, Result};
use ndarray::Array2;

use crate::normalmode::linalg_shim::LinAlg;

pub fn build_bmat(mass_au: &[f64], q_bohr: &[f64]) -> Result<Array2<f64>> {
    let ndim = q_bohr.len();
    if ndim % 3 != 0 {
        bail!("q_bohr length must be 3N");
    }
    let natom = ndim / 3;
    if mass_au.len() != natom {
        bail!("mass length must match number of atoms");
    }

    let mut bmat = Array2::<f64>::zeros((ndim, 6));
    for i in 0..natom {
        let m = mass_au[i].sqrt();
        let ix = 3 * i;
        let iy = ix + 1;
        let iz = ix + 2;

        bmat[(ix, 0)] = m;
        bmat[(iy, 1)] = m;
        bmat[(iz, 2)] = m;

        bmat[(iy, 3)] = m * q_bohr[iz];
        bmat[(iz, 3)] = -m * q_bohr[iy];

        bmat[(ix, 4)] = -m * q_bohr[iz];
        bmat[(iz, 4)] = m * q_bohr[ix];

        bmat[(ix, 5)] = m * q_bohr[iy];
        bmat[(iy, 5)] = -m * q_bohr[ix];
    }

    Ok(bmat)
}

pub fn build_smat(mass_au: &[f64], q_bohr: &[f64], la: &LinAlg) -> Result<Array2<f64>> {
    let bmat = build_bmat(mass_au, q_bohr)?;
    let bmat_tr = la.transpose(&bmat);
    Ok(la.matmul(&bmat_tr, &bmat))
}

pub fn get_eckart_projector(mass_au: &[f64], q_bohr: &[f64], la: &LinAlg) -> Result<Array2<f64>> {
    let bmat = build_bmat(mass_au, q_bohr)?;
    let bmat_tr = la.transpose(&bmat);
    let smat = la.matmul(&bmat_tr, &bmat);
    let smat_inv = pinv_symmetric(&smat, la, 1e-12)?;
    let sinv_bt = la.matmul(&smat_inv, &bmat_tr);
    Ok(la.matmul(&bmat, &sinv_bt))
}

pub fn eckart_transform(
    mass_au: &[f64],
    q_bohr: &[f64],
    hess_mw: &Array2<f64>,
    la: &LinAlg,
) -> Result<Array2<f64>> {
    let ndim = q_bohr.len();
    let rmat = get_eckart_projector(mass_au, q_bohr, la)?;
    let mut proj = Array2::<f64>::zeros((ndim, ndim));
    for i in 0..ndim {
        proj[(i, i)] = 1.0;
    }
    proj = &proj - &rmat;
    let hess_proj = la.matmul(hess_mw, &proj);
    Ok(la.matmul(&proj, &hess_proj))
}

pub fn eckart_reactionpath_transform(
    mass_au: &[f64],
    q_bohr: &[f64],
    grad_mw: &[f64],
    hess_mw: &Array2<f64>,
    la: &LinAlg,
) -> Result<Array2<f64>> {
    let ndim = q_bohr.len();
    if grad_mw.len() != ndim {
        bail!("grad_mw must have length 3N");
    }
    let rmat = get_eckart_projector(mass_au, q_bohr, la)?;

    // Page and McIver define the path tangent as the normalised mass-weighted
    // gradient, Eq. (2), and the projector as P = R + v v^t, Eq. (D1). That
    // second form is a projector only if v is orthogonal to the
    // translation/rotation subspace R spans. The paper says it is: "the path
    // tangent v is proportional to the energy gradient and therefore has no
    // component corresponding to translation and rotation in the absence of
    // external forces".
    //
    // That holds for an exact gradient. It does not hold for a computed one.
    // ORCA's gradient for the 43-atom example in `examples/` carries a net
    // force of 8.9% of its own norm and a substantial net torque -- ordinary
    // residuals of a finite DFT integration grid -- which leaves the raw unit
    // gradient with a 3.5% component inside R. Building P from it gives
    // max|P^2 - P| = 4.2e-3: not a projector, and the seven modes it is meant
    // to annihilate come out at 0.3-1.0 cm^-1 instead of zero.
    //
    // So the gradient is projected into the internal subspace before it is
    // used as the tangent. This is the paper's own intent rather than a
    // departure from it -- v is meant to be the path direction, and a net
    // force or torque is not part of the path -- and it makes P idempotent by
    // construction for any input.
    let mut g_internal = vec![0.0_f64; ndim];
    for i in 0..ndim {
        let mut acc = grad_mw[i];
        for j in 0..ndim {
            acc -= rmat[(i, j)] * grad_mw[j];
        }
        g_internal[i] = acc;
    }
    let mut proj_grad = Array2::<f64>::zeros((ndim, ndim));
    let norm = g_internal.iter().map(|v| v * v).sum::<f64>().sqrt();
    let raw_norm = grad_mw.iter().map(|v| v * v).sum::<f64>().sqrt();
    // The test is relative, not against zero. Removing the translations and
    // rotations from a gradient that was made of little else leaves a residue
    // of rounding error -- around 1e-18 -- which is not zero, so an absolute
    // test passes it through and normalises noise into a unit vector. A
    // direction chosen at random is then projected out of the Hessian, which
    // silently destroys a real vibration.
    const MIN_INTERNAL_FRACTION: f64 = 1.0e-6;
    if raw_norm <= 0.0 || norm <= MIN_INTERNAL_FRACTION * raw_norm {
        return Err(anyhow!(
            "gradient norm is zero after removing translations and rotations \
             ({norm:.3e} of {raw_norm:.3e}); there is no reaction-path direction in it"
        ));
    }
    let mut g = vec![0.0_f64; ndim];
    for i in 0..ndim {
        g[i] = g_internal[i] / norm;
    }
    for i in 0..ndim {
        for j in 0..ndim {
            proj_grad[(i, j)] = g[i] * g[j];
        }
    }

    let mut proj = Array2::<f64>::zeros((ndim, ndim));
    for i in 0..ndim {
        proj[(i, i)] = 1.0;
    }
    proj = &proj - &rmat - &proj_grad;

    let hess_proj = la.matmul(hess_mw, &proj);
    Ok(la.matmul(&proj, &hess_proj))
}

fn pinv_symmetric(a: &Array2<f64>, la: &LinAlg, tol: f64) -> Result<Array2<f64>> {
    let (evals, evecs) = la.eigh(a)?;
    let n = evals.len();
    let mut d_inv = Array2::<f64>::zeros((n, n));
    for i in 0..n {
        let lam = evals[i];
        if lam.abs() > tol {
            d_inv[(i, i)] = 1.0 / lam;
        }
    }
    let vt = la.transpose(&evecs);
    let tmp = la.matmul(&evecs, &d_inv);
    Ok(la.matmul(&tmp, &vt))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalmode::linalg_shim::Backend;

    /// Water, bent, in bohr; masses in electron masses.
    fn water() -> (Vec<f64>, Vec<f64>) {
        const AMU: f64 = 1822.888_486_209;
        let mass_au = vec![15.994_915 * AMU, 1.007_825 * AMU, 1.007_825 * AMU];
        let q_bohr = vec![
            0.0, 0.0, 0.2217, 0.0, 1.4309, -0.8867, 0.0, -1.4309, -0.8867,
        ];
        (mass_au, q_bohr)
    }

    /// A gradient with a deliberate net force and net torque, which is what a
    /// finite-grid DFT gradient really looks like: ORCA's for the 43-atom
    /// example in `examples/` carries a net force of 8.9% of its own norm.
    fn contaminated_gradient() -> Vec<f64> {
        vec![
            0.01, -0.004, 0.02, // a genuine internal displacement ...
            -0.006, 0.011, -0.008, 0.003, -0.002, 0.005,
        ]
        .iter()
        .zip([0.004_f64; 9])
        // ... plus a uniform push on every atom, i.e. a pure translation.
        .map(|(g, shift)| g + shift)
        .collect()
    }

    fn max_abs_diff(a: &Array2<f64>, b: &Array2<f64>) -> f64 {
        let (n, m) = a.dim();
        let mut worst = 0.0_f64;
        for i in 0..n {
            for j in 0..m {
                worst = worst.max((a[(i, j)] - b[(i, j)]).abs());
            }
        }
        worst
    }

    /// The Eckart projector is a projector: `R^2 = R`.
    #[test]
    fn the_eckart_projector_is_idempotent() {
        let (mass_au, q_bohr) = water();
        let la = LinAlg::new(Backend::Auto).unwrap();
        let r = get_eckart_projector(&mass_au, &q_bohr, &la).unwrap();
        let r2 = la.matmul(&r, &r);
        assert!(max_abs_diff(&r, &r2) < 1e-12, "R is not idempotent");
    }

    /// Six directions for a bent triatomic: three translations, three
    /// rotations. The trace of a projector is the dimension it projects onto.
    #[test]
    fn the_eckart_projector_spans_six_directions() {
        let (mass_au, q_bohr) = water();
        let la = LinAlg::new(Backend::Auto).unwrap();
        let r = get_eckart_projector(&mass_au, &q_bohr, &la).unwrap();
        let trace: f64 = (0..9).map(|i| r[(i, i)]).sum();
        assert!((trace - 6.0).abs() < 1e-10, "trace was {trace}");
    }

    /// The fix, tested at its cause rather than its symptom.
    ///
    /// Page and McIver's `P = R + v v^t`, Eq. (D1), is a projector only when
    /// the tangent is orthogonal to the translation/rotation subspace. A real
    /// gradient is not: a net force or torque puts part of it inside R. The
    /// implementation removes that first, so `P` is idempotent whatever it is
    /// handed -- which is what makes it annihilate exactly seven directions.
    #[test]
    fn the_reaction_path_projector_is_idempotent_despite_a_contaminated_gradient() {
        let (mass_au, q_bohr) = water();
        let la = LinAlg::new(Backend::Auto).unwrap();
        let grad = contaminated_gradient();

        // The contamination is real: the raw tangent overlaps R.
        let r = get_eckart_projector(&mass_au, &q_bohr, &la).unwrap();
        let norm: f64 = grad.iter().map(|v| v * v).sum::<f64>().sqrt();
        let raw: Vec<f64> = grad.iter().map(|v| v / norm).collect();
        let mut leak = 0.0;
        for i in 0..9 {
            let mut acc = 0.0;
            for j in 0..9 {
                acc += r[(i, j)] * raw[j];
            }
            leak += acc * acc;
        }
        assert!(
            leak.sqrt() > 0.05,
            "the fixture must actually be contaminated; leak was {}",
            leak.sqrt()
        );

        // Recover I - P from the transform by feeding it the identity, then
        // check that P projects.
        let identity = Array2::<f64>::eye(9);
        let i_minus_p =
            eckart_reactionpath_transform(&mass_au, &q_bohr, &grad, &identity, &la).unwrap();
        let squared = la.matmul(&i_minus_p, &i_minus_p);
        assert!(
            max_abs_diff(&i_minus_p, &squared) < 1e-10,
            "I - P is not idempotent: {:e}",
            max_abs_diff(&i_minus_p, &squared)
        );

        // And it removes exactly seven directions: six plus the path.
        let trace: f64 = (0..9).map(|i| i_minus_p[(i, i)]).sum();
        assert!(
            (trace - 2.0).abs() < 1e-10,
            "I - P should have rank 3N - 7 = 2; trace was {trace}"
        );
    }

    /// A gradient that is nothing but a translation has no path direction left
    /// once the translations are removed, and must be refused rather than
    /// normalised into noise.
    #[test]
    fn a_gradient_that_is_pure_translation_is_refused() {
        let (mass_au, q_bohr) = water();
        let la = LinAlg::new(Backend::Auto).unwrap();
        // The same mass-weighted shift on every atom: a pure translation.
        let mut grad = vec![0.0; 9];
        for a in 0..3 {
            grad[3 * a] = mass_au[a].sqrt() * 0.01;
        }
        let identity = Array2::<f64>::eye(9);
        let err = eckart_reactionpath_transform(&mass_au, &q_bohr, &grad, &identity, &la)
            .expect_err("no internal direction survives");
        let message = err.to_string();
        assert!(
            message.contains("after removing translations and rotations"),
            "{message}"
        );
    }
}
