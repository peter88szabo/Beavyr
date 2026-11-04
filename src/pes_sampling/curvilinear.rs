//! Curvilinear (internal-coordinate) normal-mode displacements.
//!
//! The sampled rectilinear displacement dx = Σ L_j Q_j is mapped to
//! first-order changes of bonds, angles and dihedrals, dq = B dx, with the
//! Wilson B matrix at equilibrium. The geometry is then found iteratively
//! from q(x) ≈ q_eq + dq (the usual generalized-inverse back-transformation).
//! Bending and torsional motion therefore follows arcs instead of straight
//! lines, so large-amplitude samples no longer stretch bonds artificially.
use super::input::Fragment;
use crate::normalmode::linalg_shim::{Backend, LinAlg};
use crate::optimizer::internal_coords::{
    analyze_structure, bpg_matrix, compute_internals, diff_internals, wilson_b, Bpg,
    InternalCoords,
};
use ndarray::Array2;

#[derive(Debug)]
pub struct Internal {
    ic: InternalCoords,
    q0: Vec<f64>,
    bpg: Bpg,
}

fn sub(x: &[f64], a: usize, b: usize) -> [f64; 3] {
    std::array::from_fn(|k| x[3 * a + k] - x[3 * b + k])
}
fn sin2(u: [f64; 3], v: [f64; 3]) -> f64 {
    let c = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let n = |w: [f64; 3]| w.iter().map(|x| x * x).sum::<f64>();
    n(c) / (n(u) * n(v)).max(1e-300)
}

/// Internal coordinates of the equilibrium structure, or an error when they
/// cannot represent every vibration (disconnected or linear units).
pub fn prepare(f: &Fragment) -> Result<Internal, String> {
    let nat = f.atoms.len();
    let x = &f.equilibrium;
    if nat < 3 || crate::normalmode::cartesian::is_linear(x) {
        return Err("Internal-coordinate displacement needs a nonlinear molecule with at least three atoms; use Cartesian displacement.".into());
    }
    // Covalent radii × 4/3 (the D3 convention): a bond up to ~1.27 × (rA + rB).
    let mut model = crate::optimizer::irc::default_connectivity_model(&f.atoms);
    model.rcov_disp.iter_mut().for_each(|r| *r *= 4. / 3.);
    let mut ic = analyze_structure(x, &model).map_err(|e| e.to_string())?;
    // Near-linear bends and dihedrals through them are ill-defined.
    let bent = |a: usize, b: usize, c: usize| sin2(sub(x, a, b), sub(x, c, b)) > 0.0025;
    ic.angles.retain(|&[a, b, c]| bent(a, b, c));
    ic.dihedrals
        .retain(|&[a, b, c, d]| bent(a, b, c) && bent(b, c, d));
    let bpg = bpg_matrix(x, &ic).map_err(|e| e.to_string())?;
    let nint = ic.nint();
    let nvib = 3 * nat - 6;
    let la = LinAlg::new(Backend::Auto).map_err(|e| e.to_string())?;
    let g = Array2::from_shape_fn((nint, nint), |(i, j)| {
        bpg.b[i].iter().zip(&bpg.b[j]).map(|(a, b)| a * b).sum()
    });
    let rank = if nint == 0 {
        0
    } else {
        let (values, _) = la.eigh(&g).map_err(|e| e.to_string())?;
        let max = values.iter().fold(0f64, |m, v| m.max(*v));
        values.iter().filter(|v| **v > 1e-8 * max).count()
    };
    if rank < nvib {
        return Err(format!("The detected bonds, angles and dihedrals span only {rank} of {nvib} vibrational directions (disconnected fragments or near-linear units). Use Cartesian displacement for this molecule."));
    }
    let q0 = compute_internals(x, &ic);
    Ok(Internal { ic, q0, bpg })
}

/// Equilibrium geometry `x0` (bohr) displaced curvilinearly by the
/// first-order internal-coordinate image of `dx`. The internal change is
/// applied in substeps (≤ 0.2 bohr or rad each), each solved by Gauss–Newton
/// with the B matrix of the current geometry, which stays robust for
/// large-amplitude bends and torsions.
pub fn displace(x0: &[f64], dx: &[f64], internal: &Internal) -> Result<Vec<f64>, String> {
    let b = &internal.bpg.b;
    let dq: Vec<f64> = b
        .iter()
        .map(|row| row.iter().zip(dx).map(|(a, v)| a * v).sum())
        .collect();
    let largest = dq.iter().fold(0f64, |m, v| m.max(v.abs()));
    if !largest.is_finite() {
        return Err("Non-finite internal-coordinate displacement.".into());
    }
    let substeps = ((largest / 0.2).ceil() as usize).clamp(1, 200);
    let nd = internal.ic.ndiheds();
    let la = LinAlg::new(Backend::Auto).map_err(|e| e.to_string())?;
    let mut x = x0.to_vec();
    for k in 1..=substeps {
        let s = k as f64 / substeps as f64;
        let target: Vec<f64> = internal.q0.iter().zip(&dq).map(|(q, d)| q + s * d).collect();
        // Intermediate substeps only need to stay on track; the last is exact.
        let tolerance = if k == substeps { 1e-9 } else { 1e-5 };
        let mut converged = false;
        for _ in 0..50 {
            let q = compute_internals(&x, &internal.ic);
            let residual = diff_internals(&target, &q, nd);
            let step = gauss_newton_step(&x, &internal.ic, &residual, &la)?;
            x.iter_mut().zip(&step).for_each(|(v, d)| *v += d);
            let size = step.iter().fold(0f64, |m, d| m.max(d.abs()));
            if !size.is_finite() {
                break;
            }
            if size < tolerance {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err("The internal-coordinate back-transformation did not converge for this sample. Reduce the excitation or use Cartesian displacement.".into());
        }
    }
    Ok(x)
}

/// Least-squares Cartesian step for internal residual r at x:
/// dx = (BᵀB)⁺ Bᵀ r, solved in the 3N-dimensional Cartesian space.
fn gauss_newton_step(x: &[f64], ic: &InternalCoords, r: &[f64], la: &LinAlg) -> Result<Vec<f64>, String> {
    let b = wilson_b(x, ic).map_err(|e| e.to_string())?;
    let n = x.len();
    let mut m = Array2::<f64>::zeros((n, n));
    let mut rhs = vec![0f64; n];
    for (row, ri) in b.iter().zip(r) {
        for i in 0..n {
            if row[i] == 0. {
                continue;
            }
            rhs[i] += row[i] * ri;
            for j in 0..n {
                m[(i, j)] += row[i] * row[j];
            }
        }
    }
    let (values, vectors) = la.eigh(&m).map_err(|e| e.to_string())?;
    let max = values.iter().fold(0f64, |a, v| a.max(*v));
    let mut step = vec![0f64; n];
    for (k, &value) in values.iter().enumerate() {
        // Translations/rotations (and any rank deficiency) are left out.
        if value <= 1e-10 * max {
            continue;
        }
        let c: f64 = (0..n).map(|i| vectors[(i, k)] * rhs[i]).sum::<f64>() / value;
        for i in 0..n {
            step[i] += c * vectors[(i, k)];
        }
    }
    Ok(step)
}
