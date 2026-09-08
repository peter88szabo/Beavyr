//! Pulay redundant internal coordinate machinery.
//!
//! This mirrors the minimum-search subset of Smite's `redundant_internals.py`:
//! primitive valence coordinates are kept redundant, while the Wilson G matrix
//! generalized inverse and projector remove null/redundant directions
//! algebraically.

use anyhow::{bail, Result};
use ndarray::Array2;

use crate::normalmode::linalg_shim::{Backend, LinAlg};

use super::internal_coords::{
    analyze_structure, bpg_matrix, compute_internals, diff_internals, dq_to_dx,
    internal_coords_from_bonds, ConnectivityModel, InternalCoords,
};

#[derive(Debug, Clone)]
pub struct GeneralizedInverse {
    pub matrix: Array2<f64>,
    pub rank: usize,
    pub eigenvalues: Vec<f64>,
    pub projector: Array2<f64>,
}

#[derive(Debug, Clone)]
pub struct RedundantInternalSystem {
    pub coordinates: InternalCoords,
    pub q: Vec<f64>,
    pub b: Vec<Vec<f64>>,
    pub g: Array2<f64>,
    pub ginv: Vec<Vec<f64>>,
    pub projector: Array2<f64>,
    pub rank: usize,
}

impl RedundantInternalSystem {
    pub fn nat(&self) -> usize {
        self.coordinates.nat
    }

    pub fn ncart(&self) -> usize {
        self.b
            .first()
            .map(|row| row.len())
            .unwrap_or(3 * self.nat())
    }

    pub fn nredundant(&self) -> usize {
        self.coordinates.nint().saturating_sub(self.rank)
    }
}

fn array_from_vecvec(rows: &[Vec<f64>]) -> Array2<f64> {
    let nr = rows.len();
    let nc = rows.first().map(|r| r.len()).unwrap_or(0);
    let mut out = Array2::<f64>::zeros((nr, nc));
    for i in 0..nr {
        for j in 0..nc {
            out[(i, j)] = rows[i][j];
        }
    }
    out
}

fn vecvec_from_array(a: &Array2<f64>) -> Vec<Vec<f64>> {
    let (nr, nc) = a.dim();
    let mut out = vec![vec![0.0; nc]; nr];
    for i in 0..nr {
        for j in 0..nc {
            out[i][j] = a[(i, j)];
        }
    }
    out
}

pub fn generalized_inverse_symmetric(
    matrix: &Array2<f64>,
    cutoff: f64,
) -> Result<GeneralizedInverse> {
    let (nr, nc) = matrix.dim();
    if nr != nc {
        bail!("generalized_inverse_symmetric: matrix must be square");
    }
    let mut sym = matrix.clone();
    for i in 0..nr {
        for j in 0..i {
            let v = 0.5 * (matrix[(i, j)] + matrix[(j, i)]);
            sym[(i, j)] = v;
            sym[(j, i)] = v;
        }
    }

    let la = LinAlg::new(Backend::Auto)?;
    let (evals, evecs) = la.eigh(&sym)?;
    let mut ginv = Array2::<f64>::zeros((nr, nr));
    let mut rank = 0usize;
    for k in 0..nr {
        let lam = evals[k];
        if lam <= cutoff {
            continue;
        }
        rank += 1;
        let inv = 1.0 / lam;
        for i in 0..nr {
            let vi = evecs[(i, k)];
            for j in 0..nr {
                ginv[(i, j)] += inv * vi * evecs[(j, k)];
            }
        }
    }

    let mut projector = la.matmul(&sym, &ginv);
    for i in 0..nr {
        for j in 0..i {
            let v = 0.5 * (projector[(i, j)] + projector[(j, i)]);
            projector[(i, j)] = v;
            projector[(j, i)] = v;
        }
    }

    Ok(GeneralizedInverse {
        matrix: ginv,
        rank,
        eigenvalues: evals.to_vec(),
        projector,
    })
}

pub fn build_redundant_internals(
    x: &[f64],
    model: Option<&ConnectivityModel>,
    coordinates: Option<InternalCoords>,
    bonds: Option<&[[usize; 2]]>,
    g_cutoff: f64,
) -> Result<RedundantInternalSystem> {
    if x.len() % 3 != 0 {
        bail!("redundant internals require 3N Cartesian coordinates");
    }
    let nat = x.len() / 3;
    let ic = if let Some(ic) = coordinates {
        if ic.nat != nat {
            bail!(
                "internal-coordinate atom count mismatch: got {}, expected {}",
                ic.nat,
                nat
            );
        }
        ic
    } else if let Some(bonds) = bonds {
        internal_coords_from_bonds(nat, bonds)?
    } else {
        let model =
            model.ok_or_else(|| anyhow::anyhow!("model, coordinates, or bonds required"))?;
        analyze_structure(x, model)?
    };

    let bpg = bpg_matrix(x, &ic)?;
    let b = array_from_vecvec(&bpg.b);
    let la = LinAlg::new(Backend::Auto)?;
    let g = la.matmul(&b, &la.transpose(&b));
    let gi = generalized_inverse_symmetric(&g, g_cutoff)?;
    let q = compute_internals(x, &ic);

    Ok(RedundantInternalSystem {
        coordinates: ic,
        q,
        b: bpg.b,
        g,
        ginv: vecvec_from_array(&gi.matrix),
        projector: gi.projector,
        rank: gi.rank,
    })
}

pub fn project_internal_vector(
    system: &RedundantInternalSystem,
    vector: &[f64],
) -> Result<Vec<f64>> {
    let nint = system.coordinates.nint();
    if vector.len() != nint {
        bail!(
            "internal vector length mismatch: got {}, expected {}",
            vector.len(),
            nint
        );
    }
    let mut out = vec![0.0; nint];
    for i in 0..nint {
        for j in 0..nint {
            out[i] += system.projector[(i, j)] * vector[j];
        }
    }
    Ok(out)
}

pub fn cartesian_gradient_to_redundant(
    system: &RedundantInternalSystem,
    grad_x: &[f64],
) -> Result<Vec<f64>> {
    if grad_x.len() != system.ncart() {
        bail!(
            "Cartesian gradient length mismatch: got {}, expected {}",
            grad_x.len(),
            system.ncart()
        );
    }
    let nint = system.coordinates.nint();
    let mut bg = vec![0.0; nint];
    for i in 0..nint {
        for k in 0..grad_x.len() {
            bg[i] += system.b[i][k] * grad_x[k];
        }
    }
    let mut gq = vec![0.0; nint];
    for i in 0..nint {
        for j in 0..nint {
            gq[i] += system.ginv[i][j] * bg[j];
        }
    }
    Ok(gq)
}

pub fn redundant_step_to_cartesian(
    system: &RedundantInternalSystem,
    dq: &[f64],
) -> Result<Vec<f64>> {
    let bpg = super::internal_coords::Bpg {
        b: system.b.clone(),
        inv_g: system.ginv.clone(),
    };
    dq_to_dx(&bpg, dq)
}

pub fn project_redundant_hessian(
    system: &RedundantInternalSystem,
    hessian_q: &Array2<f64>,
    redundant_penalty: f64,
) -> Result<Array2<f64>> {
    let nint = system.coordinates.nint();
    if hessian_q.dim() != (nint, nint) {
        bail!("internal Hessian shape mismatch");
    }
    let la = LinAlg::new(Backend::Auto)?;
    let p = &system.projector;
    let mut projected = la.matmul(&la.matmul(p, hessian_q), p);
    if redundant_penalty > 0.0 {
        for i in 0..nint {
            projected[(i, i)] += redundant_penalty * (1.0 - p[(i, i)]);
            for j in 0..nint {
                if i != j {
                    projected[(i, j)] -= redundant_penalty * p[(i, j)];
                }
            }
        }
    }
    for i in 0..nint {
        for j in 0..i {
            let v = 0.5 * (projected[(i, j)] + projected[(j, i)]);
            projected[(i, j)] = v;
            projected[(j, i)] = v;
        }
    }
    Ok(projected)
}

pub fn initial_redundant_hessian(
    system: &RedundantInternalSystem,
    hessian_q: &Array2<f64>,
    redundant_penalty: f64,
) -> Result<Array2<f64>> {
    project_redundant_hessian(system, hessian_q, redundant_penalty)
}

pub fn update_cartesian_from_redundant_step(
    x: &[f64],
    system: &RedundantInternalSystem,
    dq: &[f64],
    n_iter: usize,
) -> Result<Vec<f64>> {
    if x.len() != system.ncart() {
        bail!(
            "Cartesian coordinate length mismatch: got {}, expected {}",
            x.len(),
            system.ncart()
        );
    }
    let target_dq = project_internal_vector(system, dq)?;
    let mut x_new = x.to_vec();
    let dx0 = redundant_step_to_cartesian(system, &target_dq)?;
    for i in 0..x_new.len() {
        x_new[i] += dx0[i];
    }

    for _ in 0..n_iter {
        let trial = build_redundant_internals(
            &x_new,
            None,
            Some(system.coordinates.clone()),
            None,
            1.0e-7,
        )?;
        let achieved = diff_internals(&trial.q, &system.q, system.coordinates.ndiheds());
        let missing = diff_internals(&target_dq, &achieved, system.coordinates.ndiheds());
        let dx = redundant_step_to_cartesian(&trial, &missing)?;
        for i in 0..x_new.len() {
            x_new[i] += dx[i];
        }
    }
    Ok(x_new)
}
