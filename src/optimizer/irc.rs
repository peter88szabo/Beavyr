//! Intrinsic reaction coordinate following in mass-weighted Cartesian coordinates.
//!
//! This is a Rust port of Smite's `optimizer/irc.py` core algorithm: start from
//! a transition-state geometry, displace along the imaginary mode in both
//! directions, then follow steepest descent in mass-weighted coordinates with an
//! Euler or LQA predictor and a spherical corrector.

use anyhow::{bail, Result};
use ndarray::Array2;

use crate::normalmode::{project_hessian, EckartMode, HessianWeight};
use crate::normalmode::linalg_shim::{Backend, LinAlg};
use crate::optimizer::internal_coords::{
    analyze_structure, best_fit_dq_to_cart, bpg_matrix, compute_internals,
    internal_coords_from_bonds, ConnectivityModel, InternalCoords,
};
use crate::optimizer::redundant_internals::{
    build_redundant_internals, cartesian_gradient_to_redundant,
    update_cartesian_from_redundant_step,
};
use crate::optimizer::traits::Objective;
use crate::optimizer::atomic_masses::AtomicMasses;

const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;
const HARTREE_TO_KJMOL: f64 = 2625.499638;
const HARTREE_TO_KCALMOL: f64 = 627.509474;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IrcPredictor {
    Euler,
    Lqa,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IrcCoordinateMode {
    Cartesian,
    Internal,
    PulayRedundant,
}

#[derive(Debug, Clone)]
pub struct IrcOptions {
    pub masses_amu: Option<Vec<f64>>,
    pub coordinate_mode: IrcCoordinateMode,
    pub connectivity_model: Option<ConnectivityModel>,
    pub extra_bonds: Vec<[usize; 2]>,
    pub g_cutoff: f64,
    pub best_fit_iters: usize,
    pub imaginary_mode: Option<Vec<f64>>,
    pub mode_is_mass_weighted: bool,
    pub hessian: Option<Array2<f64>>,
    pub hessian_step: f64,
    pub eckart_project_hessian: bool,
    pub initial_displacement: f64,
    pub step_size: f64,
    pub step_forward: usize,
    pub step_backward: usize,
    pub step_size_min: f64,
    pub step_size_max: f64,
    pub predictor: IrcPredictor,
    pub lqa_substeps: usize,
    pub corrector_tol: f64,
    pub corrector_maxiter: usize,
    pub corrector_alpha: f64,
    pub endpoint_gradient_tol: f64,
    pub adaptive: bool,
    pub trajectory_prefix: Option<String>,
    pub trajectory_file: Option<String>,
    pub profile_file: Option<String>,
    pub verbosity: u8,
}

impl Default for IrcOptions {
    fn default() -> Self {
        Self {
            masses_amu: None,
            coordinate_mode: IrcCoordinateMode::Cartesian,
            connectivity_model: None,
            extra_bonds: Vec::new(),
            g_cutoff: 1.0e-7,
            best_fit_iters: 8,
            imaginary_mode: None,
            mode_is_mass_weighted: true,
            hessian: None,
            hessian_step: 1.0e-3,
            eckart_project_hessian: true,
            initial_displacement: 0.02,
            step_size: 0.05,
            step_forward: 50,
            step_backward: 50,
            step_size_min: 0.01,
            step_size_max: 0.20,
            predictor: IrcPredictor::Euler,
            lqa_substeps: 5,
            corrector_tol: 1.0e-4,
            corrector_maxiter: 20,
            corrector_alpha: 1.0,
            endpoint_gradient_tol: 1.0e-4,
            adaptive: true,
            trajectory_prefix: Some("irc".to_string()),
            trajectory_file: None,
            profile_file: Some("irc_energy_profile.dat".to_string()),
            verbosity: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct IrcPoint {
    pub s: f64,
    pub q: Vec<f64>,
    pub energy: f64,
    pub gradient: Vec<f64>,
    pub grad_norm_mw: f64,
    pub step_size: f64,
    pub corrector_iterations: usize,
}

#[derive(Debug, Clone)]
pub struct IrcBranch {
    pub direction: i32,
    pub points: Vec<IrcPoint>,
    pub converged: bool,
    pub message: String,
}

#[derive(Debug, Clone)]
pub struct IrcResult {
    pub atoms: Vec<String>,
    pub q_ts: Vec<f64>,
    pub ts_point: IrcPoint,
    pub imaginary_mode_mw: Vec<f64>,
    pub hessian_eigenvalues: Option<Vec<f64>>,
    pub forward: IrcBranch,
    pub backward: IrcBranch,
    pub trajectory: Vec<IrcPoint>,
    pub trajectory_file: Option<String>,
    pub profile_file: Option<String>,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

fn normalize(v: &[f64], label: &str) -> Result<Vec<f64>> {
    let n = norm(v);
    if n <= 0.0 || !n.is_finite() {
        bail!("Cannot normalize zero or invalid {label}");
    }
    Ok(v.iter().map(|x| x / n).collect())
}

fn symmetrize(a: &mut Array2<f64>) {
    let (n, m) = a.dim();
    if n != m {
        return;
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let v = 0.5 * (a[(i, j)] + a[(j, i)]);
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
}

fn finite_difference_hessian<O: Objective>(
    obj: &mut O,
    x: &[f64],
    step: f64,
) -> Result<Array2<f64>> {
    if step <= 0.0 {
        bail!("IRC finite-difference Hessian step must be positive");
    }
    let n = x.len();
    let mut hess = Array2::<f64>::zeros((n, n));
    let mut xp = x.to_vec();
    let mut xm = x.to_vec();
    let mut gp = vec![0.0; n];
    let mut gm = vec![0.0; n];
    for j in 0..n {
        xp[j] = x[j] + step;
        xm[j] = x[j] - step;
        obj.gradient(&xp, &mut gp)?;
        obj.gradient(&xm, &mut gm)?;
        for i in 0..n {
            hess[(i, j)] = (gp[i] - gm[i]) / (2.0 * step);
        }
        xp[j] = x[j];
        xm[j] = x[j];
    }
    symmetrize(&mut hess);
    Ok(hess)
}

pub(crate) fn covalent_radius_bohr(element: &str) -> f64 {
    let angstrom = match element.to_ascii_uppercase().as_str() {
        "H" => 0.31,
        "B" => 0.84,
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "F" => 0.57,
        "P" => 1.07,
        "S" => 1.05,
        "CL" => 1.02,
        "BR" => 1.20,
        "I" => 1.39,
        _ => 0.77,
    };
    angstrom / BOHR_TO_ANGSTROM
}

pub(crate) fn default_connectivity_model(atoms: &[String]) -> ConnectivityModel {
    ConnectivityModel {
        rcov_disp: atoms.iter().map(|a| covalent_radius_bohr(a)).collect(),
        kcn: 7.5,
        facmin: 0.6,
    }
}

fn internal_coords_for_q(atoms: &[String], q: &[f64], opts: &IrcOptions) -> Result<InternalCoords> {
    let model = opts
        .connectivity_model
        .clone()
        .unwrap_or_else(|| default_connectivity_model(atoms));
    let mut ic = analyze_structure(q, &model)?;
    if !opts.extra_bonds.is_empty() {
        let mut bonds = ic.bonds.clone();
        for &bond in &opts.extra_bonds {
            let [a, b] = bond;
            if a >= atoms.len() || b >= atoms.len() || a == b {
                bail!("IRC extra_bonds contains invalid bond ({a}, {b})");
            }
            let canonical = if a < b { [a, b] } else { [b, a] };
            if !bonds.iter().any(|existing| *existing == canonical) {
                bonds.push(canonical);
            }
        }
        ic = internal_coords_from_bonds(atoms.len(), &bonds)?;
    }
    if ic.nint() == 0 {
        bail!("IRC internal-coordinate mode found no internal coordinates");
    }
    Ok(ic)
}

fn internal_gradient(b: &[Vec<f64>], inv_g: &[Vec<f64>], grad_x: &[f64]) -> Vec<f64> {
    let nint = b.len();
    let mut bg = vec![0.0; nint];
    for i in 0..nint {
        bg[i] = b[i]
            .iter()
            .zip(grad_x.iter())
            .map(|(bij, gj)| bij * gj)
            .sum();
    }
    let mut gq = vec![0.0; nint];
    for i in 0..nint {
        for j in 0..nint {
            gq[i] += inv_g[i][j] * bg[j];
        }
    }
    gq
}

fn grad_mw_to_cart(grad_y: &[f64], sqrt_m: &[f64]) -> Vec<f64> {
    grad_y
        .iter()
        .zip(sqrt_m.iter())
        .map(|(g, m)| g * m)
        .collect()
}

fn normalize_cartesian_delta_to_mw_step(
    q_current: &[f64],
    q_trial: &[f64],
    sqrt_m: &[f64],
    step_size: f64,
) -> Result<Vec<f64>> {
    let mut delta_y: Vec<f64> = q_trial
        .iter()
        .zip(q_current.iter())
        .zip(sqrt_m.iter())
        .map(|((qt, qc), sm)| (qt - qc) * sm)
        .collect();
    let n = norm(&delta_y);
    if n <= 0.0 || !n.is_finite() {
        bail!("IRC internal-coordinate predictor produced a zero Cartesian displacement");
    }
    for v in &mut delta_y {
        *v *= step_size / n;
    }
    Ok(delta_y)
}

fn masses_from_atoms(atoms: &[String], opts: &IrcOptions) -> Result<Vec<f64>> {
    if let Some(masses) = &opts.masses_amu {
        if masses.len() != atoms.len() {
            bail!(
                "IRC mass count mismatch: got {}, expected {}",
                masses.len(),
                atoms.len()
            );
        }
        for &m in masses {
            if m <= 0.0 || !m.is_finite() {
                bail!("IRC requires positive finite atomic masses");
            }
        }
        return Ok(masses.clone());
    }
    let table = AtomicMasses::new();
    let mut masses = Vec::with_capacity(atoms.len());
    for atom in atoms {
        let m = table.get_mass(atom);
        if m <= 0.0 || !m.is_finite() {
            bail!("IRC requires positive finite atomic mass for {atom}");
        }
        masses.push(m);
    }
    Ok(masses)
}

fn sqrt_mass_by_coord(masses_amu: &[f64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(3 * masses_amu.len());
    for &m in masses_amu {
        let sm = m.sqrt();
        out.push(sm);
        out.push(sm);
        out.push(sm);
    }
    out
}

fn cart_to_mw(q: &[f64], sqrt_m: &[f64]) -> Vec<f64> {
    q.iter().zip(sqrt_m.iter()).map(|(x, m)| x * m).collect()
}

fn mw_to_cart(y: &[f64], sqrt_m: &[f64]) -> Vec<f64> {
    y.iter().zip(sqrt_m.iter()).map(|(x, m)| x / m).collect()
}

fn grad_cart_to_mw(grad_x: &[f64], sqrt_m: &[f64]) -> Vec<f64> {
    grad_x
        .iter()
        .zip(sqrt_m.iter())
        .map(|(g, m)| g / m)
        .collect()
}

fn hess_cart_to_mw(hess_x: &Array2<f64>, sqrt_m: &[f64]) -> Array2<f64> {
    let (n, m) = hess_x.dim();
    let mut out = Array2::<f64>::zeros((n, m));
    for i in 0..n {
        for j in 0..m {
            out[(i, j)] = hess_x[(i, j)] / sqrt_m[i] / sqrt_m[j];
        }
    }
    out
}

fn imaginary_mode_from_mass_weighted_hessian(
    hess_mw: &Array2<f64>,
) -> Result<(Vec<f64>, Vec<f64>)> {
    let mut hess_mw = hess_mw.clone();
    symmetrize(&mut hess_mw);
    let solver = LinAlg::new(Backend::Auto)?;
    let (evals, evecs) = solver.eigh(&hess_mw)?;
    let mut idx = 0usize;
    for i in 1..evals.len() {
        if evals[i] < evals[idx] {
            idx = i;
        }
    }
    if evals[idx] >= 0.0 {
        bail!("IRC initialization requires a Hessian with at least one negative eigenvalue");
    }
    let mode: Vec<f64> = (0..evals.len()).map(|i| evecs[(i, idx)]).collect();
    Ok((normalize(&mode, "imaginary mode")?, evals.to_vec()))
}

fn imaginary_mode_from_hessian(
    hess_x: &Array2<f64>,
    sqrt_m: &[f64],
) -> Result<(Vec<f64>, Vec<f64>)> {
    let hess_mw = hess_cart_to_mw(hess_x, sqrt_m);
    imaginary_mode_from_mass_weighted_hessian(&hess_mw)
}

fn prepare_imaginary_mode<O: Objective>(
    obj: &mut O,
    q_ts: &[f64],
    masses_amu: &[f64],
    sqrt_m: &[f64],
    opts: &IrcOptions,
) -> Result<(Vec<f64>, Option<Vec<f64>>)> {
    if let Some(mode) = &opts.imaginary_mode {
        if mode.len() != q_ts.len() {
            bail!(
                "IRC imaginary_mode length mismatch: got {}, expected {}",
                mode.len(),
                q_ts.len()
            );
        }
        let mode_mw = if opts.mode_is_mass_weighted {
            mode.clone()
        } else {
            cart_to_mw(mode, sqrt_m)
        };
        return Ok((normalize(&mode_mw, "imaginary mode")?, None));
    }
    let hess_x = if let Some(h) = &opts.hessian {
        h.clone()
    } else {
        finite_difference_hessian(obj, q_ts, opts.hessian_step)?
    };
    let (mode, evals) = if opts.eckart_project_hessian {
        let hess_mw = project_hessian(
            masses_amu,
            q_ts,
            &hess_x,
            EckartMode::VibRot,
            None,
            HessianWeight::MassWeighted,
        )?;
        imaginary_mode_from_mass_weighted_hessian(&hess_mw)?
    } else {
        imaginary_mode_from_hessian(&hess_x, sqrt_m)?
    };
    Ok((mode, Some(evals)))
}

fn point_from_y<O: Objective>(
    obj: &mut O,
    y: &[f64],
    sqrt_m: &[f64],
    s: f64,
    step_size: f64,
    corrector_iterations: usize,
) -> Result<IrcPoint> {
    let q = mw_to_cart(y, sqrt_m);
    let mut grad_x = vec![0.0; q.len()];
    let energy = obj.energy_gradient(&q, &mut grad_x)?;
    let grad_y = grad_cart_to_mw(&grad_x, sqrt_m);
    Ok(IrcPoint {
        s,
        q,
        energy,
        gradient: grad_x,
        grad_norm_mw: norm(&grad_y),
        step_size,
        corrector_iterations,
    })
}

fn mat_vec(a: &Array2<f64>, x: &[f64]) -> Vec<f64> {
    let (nr, nc) = a.dim();
    let mut out = vec![0.0; nr];
    for i in 0..nr {
        for j in 0..nc {
            out[i] += a[(i, j)] * x[j];
        }
    }
    out
}

fn lqa_predictor(
    grad_y: &[f64],
    hess_y: &Array2<f64>,
    step_size: f64,
    substeps: usize,
) -> Result<Vec<f64>> {
    let mut displacement = vec![0.0; grad_y.len()];
    let ds = step_size / substeps.max(1) as f64;
    for _ in 0..substeps.max(1) {
        let hd = mat_vec(hess_y, &displacement);
        let local_grad: Vec<f64> = grad_y.iter().zip(hd.iter()).map(|(g, h)| g + h).collect();
        let tangent = normalize(&local_grad, "LQA gradient")?;
        for i in 0..displacement.len() {
            displacement[i] -= ds * tangent[i];
        }
    }
    Ok(displacement)
}

fn cartesian_predictor_step<O: Objective>(
    obj: &mut O,
    y: &[f64],
    sqrt_m: &[f64],
    grad_y: &[f64],
    step_size: f64,
    opts: &IrcOptions,
) -> Result<Vec<f64>> {
    match opts.predictor {
        IrcPredictor::Euler => {
            let tangent = normalize(grad_y, "mass-weighted gradient")?;
            Ok(tangent.into_iter().map(|v| -step_size * v).collect())
        }
        IrcPredictor::Lqa => {
            let q = mw_to_cart(y, sqrt_m);
            let hess_x = finite_difference_hessian(obj, &q, opts.hessian_step)?;
            let hess_y = hess_cart_to_mw(&hess_x, sqrt_m);
            lqa_predictor(grad_y, &hess_y, step_size, opts.lqa_substeps)
        }
    }
}

fn predictor_step<O: Objective>(
    obj: &mut O,
    atoms: &[String],
    y: &[f64],
    sqrt_m: &[f64],
    grad_y: &[f64],
    step_size: f64,
    opts: &IrcOptions,
) -> Result<Vec<f64>> {
    match opts.coordinate_mode {
        IrcCoordinateMode::Cartesian => {
            cartesian_predictor_step(obj, y, sqrt_m, grad_y, step_size, opts)
        }
        IrcCoordinateMode::Internal => {
            let q = mw_to_cart(y, sqrt_m);
            let grad_x = grad_mw_to_cart(grad_y, sqrt_m);
            let ic = internal_coords_for_q(atoms, &q, opts)?;
            let bpg = bpg_matrix(&q, &ic)?;
            let grad_q = internal_gradient(&bpg.b, &bpg.inv_g, &grad_x);
            let tangent = normalize(&grad_q, "internal gradient")?;
            let dq: Vec<f64> = tangent.iter().map(|v| -step_size * v).collect();
            let mut q_trial = q.clone();
            let qs = compute_internals(&q, &ic);
            best_fit_dq_to_cart(&mut q_trial, &qs, &dq, &ic, &bpg, opts.best_fit_iters)?;
            normalize_cartesian_delta_to_mw_step(&q, &q_trial, sqrt_m, step_size)
        }
        IrcCoordinateMode::PulayRedundant => {
            let q = mw_to_cart(y, sqrt_m);
            let grad_x = grad_mw_to_cart(grad_y, sqrt_m);
            let ic = internal_coords_for_q(atoms, &q, opts)?;
            let system = build_redundant_internals(&q, None, Some(ic), None, opts.g_cutoff)?;
            let grad_q = cartesian_gradient_to_redundant(&system, &grad_x)?;
            let tangent = normalize(&grad_q, "Pulay redundant internal gradient")?;
            let dq: Vec<f64> = tangent.iter().map(|v| -step_size * v).collect();
            let q_trial =
                update_cartesian_from_redundant_step(&q, &system, &dq, opts.best_fit_iters)?;
            normalize_cartesian_delta_to_mw_step(&q, &q_trial, sqrt_m, step_size)
        }
    }
}

fn correct_irc_point<O: Objective>(
    obj: &mut O,
    y_previous: &[f64],
    y_initial: &[f64],
    sqrt_m: &[f64],
    step_size: f64,
    opts: &IrcOptions,
) -> Result<(Vec<f64>, usize, bool)> {
    let mut y = y_initial.to_vec();
    let ndim = y.len();
    for iteration in 1..=opts.corrector_maxiter {
        let mut radial: Vec<f64> = y
            .iter()
            .zip(y_previous.iter())
            .map(|(a, b)| a - b)
            .collect();
        let mut radial_norm = norm(&radial);
        if radial_norm <= 0.0 || !radial_norm.is_finite() {
            radial = vec![0.0; ndim];
            radial[0] = 1.0;
            radial_norm = 1.0;
        }
        let u: Vec<f64> = radial.iter().map(|v| v / radial_norm).collect();
        let q = mw_to_cart(&y, sqrt_m);
        let mut grad_x = vec![0.0; q.len()];
        obj.gradient(&q, &mut grad_x)?;
        let grad_y = grad_cart_to_mw(&grad_x, sqrt_m);
        let ug = dot(&u, &grad_y);
        let g_perp: Vec<f64> = grad_y
            .iter()
            .zip(u.iter())
            .map(|(g, u)| g - u * ug)
            .collect();
        let g_perp_norm = norm(&g_perp);
        if g_perp_norm < opts.corrector_tol {
            let corrected: Vec<f64> = y_previous
                .iter()
                .zip(u.iter())
                .map(|(yp, u)| yp + step_size * u)
                .collect();
            return Ok((corrected, iteration, true));
        }
        let mut alpha = opts.corrector_alpha;
        let max_move = 0.5 * step_size;
        if alpha * g_perp_norm > max_move {
            alpha = max_move / g_perp_norm;
        }
        for i in 0..ndim {
            y[i] -= alpha * g_perp[i];
        }
        radial = y
            .iter()
            .zip(y_previous.iter())
            .map(|(a, b)| a - b)
            .collect();
        radial_norm = norm(&radial);
        if radial_norm <= 0.0 || !radial_norm.is_finite() {
            return Ok((y_initial.to_vec(), iteration, false));
        }
        for i in 0..ndim {
            y[i] = y_previous[i] + step_size * radial[i] / radial_norm;
        }
    }
    Ok((y, opts.corrector_maxiter, false))
}

fn relative_kj(point: &IrcPoint, reference_energy: f64) -> f64 {
    (point.energy - reference_energy) * HARTREE_TO_KJMOL
}

fn relative_kcal(point: &IrcPoint, reference_energy: f64) -> f64 {
    (point.energy - reference_energy) * HARTREE_TO_KCALMOL
}

fn xyz_frame(atoms: &[String], point: &IrcPoint, branch: &str, reference_energy: f64) -> String {
    let mut out = String::new();
    out.push_str(&format!("{}\n", atoms.len()));
    out.push_str(&format!(
        "irc_part={} s={:.8} E_hartree={:.12} dE_from_TS_kJmol={:.8} grad_norm_mw={:.8e} h={:.8} corrector_iterations={}\n",
        branch,
        point.s,
        point.energy,
        relative_kj(point, reference_energy),
        point.grad_norm_mw,
        point.step_size,
        point.corrector_iterations
    ));
    for (i, atom) in atoms.iter().enumerate() {
        let k = 3 * i;
        out.push_str(&format!(
            "{:<2} {:16.8} {:16.8} {:16.8}\n",
            atom,
            point.q[k] * BOHR_TO_ANGSTROM,
            point.q[k + 1] * BOHR_TO_ANGSTROM,
            point.q[k + 2] * BOHR_TO_ANGSTROM
        ));
    }
    out
}

fn write_branch_xyz(
    filename: &Option<String>,
    atoms: &[String],
    points: &[IrcPoint],
    branch: &str,
    reference_energy: f64,
) -> Result<()> {
    let Some(filename) = filename else {
        return Ok(());
    };
    let mut out = String::new();
    for point in points {
        out.push_str(&xyz_frame(atoms, point, branch, reference_energy));
    }
    std::fs::write(filename, out)?;
    Ok(())
}

fn write_trajectory_xyz(
    filename: &Option<String>,
    atoms: &[String],
    trajectory: &[IrcPoint],
    reference_energy: f64,
) -> Result<()> {
    let Some(filename) = filename else {
        return Ok(());
    };
    let mut out = String::new();
    for point in trajectory {
        let branch = if point.s < 0.0 {
            "backward"
        } else if point.s > 0.0 {
            "forward"
        } else {
            "TS"
        };
        out.push_str(&xyz_frame(atoms, point, branch, reference_energy));
    }
    std::fs::write(filename, out)?;
    Ok(())
}

fn write_profile_dat(
    filename: &Option<String>,
    trajectory: &[IrcPoint],
    reference_energy: f64,
) -> Result<()> {
    let Some(filename) = filename else {
        return Ok(());
    };
    let mut out = String::new();
    out.push_str("# IRC energy profile\n");
    out.push_str("# Relative energies are measured from the TS reference energy.\n");
    out.push_str("# Columns: index branch s E_hartree dE_from_TS_kcalmol dE_from_TS_kJmol grad_norm_mw step_size corrector_iterations\n");
    for (index, point) in trajectory.iter().enumerate() {
        let branch = if point.s < 0.0 {
            "backward"
        } else if point.s > 0.0 {
            "forward"
        } else {
            "TS"
        };
        out.push_str(&format!(
            "{:6} {:>10} {:16.8} {:22.12} {:20.8} {:20.8} {:16.8e} {:12.6} {:6}\n",
            index,
            branch,
            point.s,
            point.energy,
            relative_kcal(point, reference_energy),
            relative_kj(point, reference_energy),
            point.grad_norm_mw,
            point.step_size,
            point.corrector_iterations
        ));
    }
    std::fs::write(filename, out)?;
    Ok(())
}

fn follow_branch<O: Objective>(
    obj: &mut O,
    atoms: &[String],
    y_start: &[f64],
    sqrt_m: &[f64],
    direction: i32,
    max_steps: usize,
    branch_file: Option<String>,
    reference_energy: f64,
    opts: &IrcOptions,
) -> Result<IrcBranch> {
    let branch_label = if direction > 0 { "forward" } else { "backward" };
    let mut points = Vec::new();
    let mut current_y = y_start.to_vec();
    let mut current_step = opts.step_size;
    let mut s = opts.initial_displacement;
    let mut point = point_from_y(obj, &current_y, sqrt_m, s, current_step, 0)?;
    points.push(point.clone());

    for istep in 1..=max_steps {
        if point.grad_norm_mw < opts.endpoint_gradient_tol {
            write_branch_xyz(&branch_file, atoms, &points, branch_label, reference_energy)?;
            return Ok(IrcBranch {
                direction,
                points,
                converged: true,
                message: "IRC branch reached the endpoint gradient threshold.".to_string(),
            });
        }
        let grad_y = grad_cart_to_mw(&point.gradient, sqrt_m);
        let mut accepted = false;
        let mut last_iterations = 0usize;
        let mut y_corr = current_y.clone();
        while current_step >= opts.step_size_min {
            let predicted_delta =
                predictor_step(obj, atoms, &current_y, sqrt_m, &grad_y, current_step, opts)?;
            let y_pred: Vec<f64> = current_y
                .iter()
                .zip(predicted_delta.iter())
                .map(|(a, b)| a + b)
                .collect();
            let (trial_corr, n_corr, corr_ok) =
                correct_irc_point(obj, &current_y, &y_pred, sqrt_m, current_step, opts)?;
            last_iterations = n_corr;
            if corr_ok {
                y_corr = trial_corr;
                accepted = true;
                break;
            }
            current_step *= 0.5;
        }
        if !accepted {
            write_branch_xyz(&branch_file, atoms, &points, branch_label, reference_energy)?;
            return Ok(IrcBranch {
                direction,
                points,
                converged: false,
                message: "IRC branch stopped because the corrector could not converge above the minimum step size.".to_string(),
            });
        }
        s += current_step;
        current_y = y_corr;
        point = point_from_y(obj, &current_y, sqrt_m, s, current_step, last_iterations)?;
        points.push(point.clone());
        if opts.verbosity > 0 {
            println!(
                "{:<8} {:4} s={:10.5} E={:18.10} dE_TS={:12.5} kJ/mol |g_mw|={:10.4e} h={:8.4} corr={:3}",
                branch_label,
                istep,
                s,
                point.energy,
                relative_kj(&point, reference_energy),
                point.grad_norm_mw,
                current_step,
                last_iterations
            );
        }
        if opts.adaptive {
            if last_iterations <= 3 {
                current_step = opts.step_size_max.min(1.2 * current_step);
            } else if last_iterations > 10 {
                current_step = opts.step_size_min.max(0.5 * current_step);
            }
        }
    }

    write_branch_xyz(&branch_file, atoms, &points, branch_label, reference_energy)?;
    Ok(IrcBranch {
        direction,
        points,
        converged: false,
        message: "Maximum IRC steps reached.".to_string(),
    })
}

pub fn follow_irc<O: Objective>(
    obj: &mut O,
    atoms: Vec<String>,
    q_ts: Vec<f64>,
    opts: IrcOptions,
) -> Result<IrcResult> {
    if q_ts.len() != 3 * atoms.len() {
        bail!("IRC TS coordinates must have length 3N");
    }
    if opts.initial_displacement <= 0.0 || opts.step_size <= 0.0 {
        bail!("IRC initial_displacement and step_size must be positive");
    }
    if opts.step_size_min <= 0.0 || opts.step_size_max < opts.step_size_min {
        bail!("IRC step-size bounds are invalid");
    }
    let masses = masses_from_atoms(&atoms, &opts)?;
    let sqrt_m = sqrt_mass_by_coord(&masses);
    let y_ts = cart_to_mw(&q_ts, &sqrt_m);
    let (imag_mode_mw, evals) = prepare_imaginary_mode(obj, &q_ts, &masses, &sqrt_m, &opts)?;
    let ts_point = point_from_y(obj, &y_ts, &sqrt_m, 0.0, opts.step_size, 0)?;

    if opts.verbosity > 0 {
        println!("Intrinsic reaction coordinate following");
        println!("TS energy [Eh]={:.12}", ts_point.energy);
        println!(
            "step_size={:.5} initial_displacement={:.5}",
            opts.step_size, opts.initial_displacement
        );
        println!(
            "step_forward={} step_backward={}",
            opts.step_forward, opts.step_backward
        );
        println!("adaptive step size: {}", opts.adaptive);
        println!(
            "profile file: {}",
            opts.profile_file.as_deref().unwrap_or("-")
        );
        let mode_source = if opts.imaginary_mode.is_some() {
            "user-provided imaginary_mode"
        } else {
            "lowest Hessian eigenvector"
        };
        println!("IRC mode source: {mode_source}");
        println!(
            "IRC Hessian Eckart projection: {}",
            opts.eckart_project_hessian
        );
        println!("IRC imaginary-mode norm: {:.8}", norm(&imag_mode_mw));
        if let Some(evals) = &evals {
            let negative_count = evals.iter().filter(|&&v| v < 0.0).count();
            let lowest = evals.iter().fold(f64::INFINITY, |m, &v| m.min(v));
            println!("IRC Hessian negative eigenvalues: {negative_count}");
            println!("IRC lowest Hessian eigenvalue [mass-weighted a.u.]: {lowest:.8e}");
        }
        println!("IRC branch convention:");
        println!("  forward  starts along + imaginary mode");
        println!("  backward starts along - imaginary mode");
    }

    let prefix = opts.trajectory_prefix.clone();
    let forward_file = prefix.as_ref().map(|p| format!("{p}_forward.xyz"));
    let backward_file = prefix.as_ref().map(|p| format!("{p}_backward.xyz"));
    let y_forward: Vec<f64> = y_ts
        .iter()
        .zip(imag_mode_mw.iter())
        .map(|(y, m)| y + opts.initial_displacement * m)
        .collect();
    let y_backward: Vec<f64> = y_ts
        .iter()
        .zip(imag_mode_mw.iter())
        .map(|(y, m)| y - opts.initial_displacement * m)
        .collect();

    let forward = follow_branch(
        obj,
        &atoms,
        &y_forward,
        &sqrt_m,
        1,
        opts.step_forward,
        forward_file,
        ts_point.energy,
        &opts,
    )?;
    let backward = follow_branch(
        obj,
        &atoms,
        &y_backward,
        &sqrt_m,
        -1,
        opts.step_backward,
        backward_file,
        ts_point.energy,
        &opts,
    )?;

    let mut trajectory = Vec::new();
    for point in backward.points.iter().rev() {
        let mut p = point.clone();
        p.s = -p.s.abs();
        trajectory.push(p);
    }
    trajectory.push(ts_point.clone());
    for point in &forward.points {
        let mut p = point.clone();
        p.s = p.s.abs();
        trajectory.push(p);
    }

    let trajectory_file = opts
        .trajectory_file
        .clone()
        .or_else(|| prefix.map(|p| format!("{p}_full.xyz")));
    write_trajectory_xyz(&trajectory_file, &atoms, &trajectory, ts_point.energy)?;
    write_profile_dat(&opts.profile_file, &trajectory, ts_point.energy)?;

    if opts.verbosity > 0 {
        println!("IRC output files:");
        println!(
            "  combined trajectory: {}",
            trajectory_file.as_deref().unwrap_or("-")
        );
        println!(
            "  energy profile: {}",
            opts.profile_file.as_deref().unwrap_or("-")
        );
    }

    Ok(IrcResult {
        atoms,
        q_ts,
        ts_point,
        imaginary_mode_mw: imag_mode_mw,
        hessian_eigenvalues: evals,
        forward,
        backward,
        trajectory,
        trajectory_file,
        profile_file: opts.profile_file,
    })
}
