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

    let mut proj_grad = Array2::<f64>::zeros((ndim, ndim));
    let norm = grad_mw.iter().map(|v| v * v).sum::<f64>().sqrt();
    if norm <= 0.0 {
        return Err(anyhow!(
            "gradient norm is zero; cannot project reaction path"
        ));
    }
    let mut g = vec![0.0_f64; ndim];
    for i in 0..ndim {
        g[i] = grad_mw[i] / norm;
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
