//! RDA quasi-transition-state guess generation.
//!
//! This is the first Rust port of Smite's `rda_guess.py`: it implements the
//! standard RDA midpoint/beta-bracketing workflow with conditional Cartesian
//! or internal-coordinate relaxation, plus Smite's poor-man NEB constrained
//! image relaxation path.

use anyhow::{bail, Result};
use ndarray::Array2;

use crate::normalmode::linalg_shim::{Backend, LinAlg};
use crate::optimizer::constrained::{
    optimize_constrained_with_internals, ConstrainedCoordinateMode, ConstrainedOptimizationOptions,
    ConstraintCoordinate, ConstraintTarget,
};
use crate::optimizer::internal_coords::{
    best_fit_dq_to_cart, bpg_matrix, compute_internals, diff_internals, internal_coords_from_bonds,
    InternalCoords,
};
use crate::optimizer::minimum_cartesian::step_limit;
use crate::optimizer::traits::Objective;

const EV_TO_HARTREE: f64 = 0.036_749_322_175_654_99;
const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RdaDistanceMode {
    Cartesian,
    ReactiveInternal,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RdaConditionalCoordinates {
    Cartesian,
    Internal,
}

#[derive(Debug, Clone)]
pub struct RdaOptions {
    pub active_atoms: Option<Vec<usize>>,
    pub align: bool,
    pub align_atoms: Option<Vec<usize>>,
    pub distance_mode: RdaDistanceMode,
    pub conditional_coordinates: RdaConditionalCoordinates,
    pub first_energy_tol_ev: f64,
    pub later_energy_tol_ev: f64,
    pub distance_tol: f64,
    pub max_conditional_steps: usize,
    pub beta_start: f64,
    pub beta_step: f64,
    pub max_bracket_attempts: usize,
    pub gamma_values: Vec<f64>,
    pub cartesian_step_max_component: f64,
    pub internal_step_max_component: f64,
    pub internal_best_fit_iters: usize,
    pub fallback_to_cartesian: bool,
    pub reactive_bonds: Vec<[usize; 2]>,
    pub reactive_angles: Vec<[usize; 3]>,
    pub reactive_dihedrals: Vec<[usize; 4]>,
    pub reactive_coordinate_weights: Option<Vec<f64>>,
    pub extra_bonds: Vec<[usize; 2]>,
    pub extra_angles: Vec<[usize; 3]>,
    pub extra_dihedrals: Vec<[usize; 4]>,
    pub poormans_neb: bool,
    pub poormans_neb_images: usize,
    pub poormans_neb_maxiter: usize,
    pub poormans_neb_file: Option<String>,
    pub poormans_neb_profile_file: Option<String>,
    pub rda_trajectory_file: Option<String>,
    pub rda_chain_file: Option<String>,
    pub verbosity: u8,
}

impl Default for RdaOptions {
    fn default() -> Self {
        Self {
            active_atoms: None,
            align: true,
            align_atoms: None,
            distance_mode: RdaDistanceMode::Cartesian,
            conditional_coordinates: RdaConditionalCoordinates::Internal,
            first_energy_tol_ev: 0.01,
            later_energy_tol_ev: 0.05,
            distance_tol: 0.05,
            max_conditional_steps: 25,
            beta_start: 0.5,
            beta_step: 0.1,
            max_bracket_attempts: 6,
            gamma_values: (1..=9).map(|i| i as f64 * 0.1).collect(),
            cartesian_step_max_component: 0.10,
            internal_step_max_component: 0.10,
            internal_best_fit_iters: 8,
            fallback_to_cartesian: true,
            reactive_bonds: Vec::new(),
            reactive_angles: Vec::new(),
            reactive_dihedrals: Vec::new(),
            reactive_coordinate_weights: None,
            extra_bonds: Vec::new(),
            extra_angles: Vec::new(),
            extra_dihedrals: Vec::new(),
            poormans_neb: false,
            poormans_neb_images: 9,
            poormans_neb_maxiter: 50,
            poormans_neb_file: Some("final_reaction_path.xyz".to_string()),
            poormans_neb_profile_file: Some("poormans_neb_profile.dat".to_string()),
            rda_trajectory_file: Some(
                "detailed_optimization_steps_of_reaction_path.xyz".to_string(),
            ),
            rda_chain_file: Some("rda_quasi_ts_search_chain.xyz".to_string()),
            verbosity: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct RdaPoint {
    pub label: String,
    pub q_start: Vec<f64>,
    pub q_copt: Vec<f64>,
    pub energy: f64,
    pub nsteps: usize,
    pub direction: String,
    pub delta_d_is: f64,
    pub delta_d_fs: f64,
    pub beta: Option<f64>,
    pub gamma: Option<f64>,
    pub message: String,
}

#[derive(Debug, Clone)]
pub struct RdaGuessResult {
    pub atoms: Vec<String>,
    pub q: Vec<f64>,
    pub reaction_direction: Option<Vec<f64>>,
    pub points: Vec<RdaPoint>,
    pub direction: String,
    pub bracketed: bool,
    pub message: String,
    pub rda_trajectory_file: Option<String>,
    pub rda_chain_file: Option<String>,
    pub poormans_neb_file: Option<String>,
    pub poormans_neb_profile_file: Option<String>,
}

impl RdaGuessResult {
    pub fn xyz(&self) -> String {
        xyz_string(&self.atoms, &self.q, "RDA quasi-TS guess")
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

fn normalize(v: &[f64]) -> Option<Vec<f64>> {
    let n = norm(v);
    if n <= 0.0 || !n.is_finite() {
        return None;
    }
    Some(v.iter().map(|x| x / n).collect())
}

fn interpolate(a: &[f64], b: &[f64], alpha: f64) -> Vec<f64> {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (1.0 - alpha) * x + alpha * y)
        .collect()
}

fn xyz_string(atoms: &[String], q: &[f64], comment: &str) -> String {
    let mut out = format!("{}\n{}\n", atoms.len(), comment);
    for (i, atom) in atoms.iter().enumerate() {
        let k = 3 * i;
        out.push_str(&format!(
            "{:<2} {:16.10} {:16.10} {:16.10}\n",
            atom,
            q[k] * BOHR_TO_ANGSTROM,
            q[k + 1] * BOHR_TO_ANGSTROM,
            q[k + 2] * BOHR_TO_ANGSTROM
        ));
    }
    out
}

fn validate_inputs(atoms: &[String], q_is: &[f64], q_fs: &[f64]) -> Result<()> {
    if atoms.is_empty() {
        bail!("RDA requires at least one atom");
    }
    if q_is.len() != 3 * atoms.len() || q_fs.len() != 3 * atoms.len() {
        bail!("RDA endpoint coordinates must have length 3N");
    }
    Ok(())
}

fn active_atoms(opts: &RdaOptions, nat: usize) -> Result<Vec<usize>> {
    if let Some(active) = &opts.active_atoms {
        if active.is_empty() {
            bail!("RDA active_atoms must not be empty");
        }
        let mut seen = vec![false; nat];
        for &i in active {
            if i >= nat {
                bail!("RDA active atom index {i} out of range 0..{}", nat - 1);
            }
            if seen[i] {
                bail!("RDA active_atoms contains duplicate atom {i}");
            }
            seen[i] = true;
        }
        Ok(active.clone())
    } else {
        Ok((0..nat).collect())
    }
}

fn covalent_radius_bohr(element: &str) -> f64 {
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

fn distance_between(q: &[f64], i: usize, j: usize) -> f64 {
    let dx = q[3 * i] - q[3 * j];
    let dy = q[3 * i + 1] - q[3 * j + 1];
    let dz = q[3 * i + 2] - q[3 * j + 2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn union_endpoint_internal_coordinates(
    atoms: &[String],
    q_is: &[f64],
    q_fs: &[f64],
    extra_bonds: &[[usize; 2]],
) -> Result<InternalCoords> {
    let nat = atoms.len();
    let mut bonds = Vec::<[usize; 2]>::new();
    let mut seen = vec![false; nat * nat];
    let mut push_bond = |a: usize, b: usize, bonds: &mut Vec<[usize; 2]>| -> Result<()> {
        if a >= nat || b >= nat || a == b {
            bail!("RDA internal conditional bond index out of range or self-bond ({a}, {b})");
        }
        let (i, j) = if a < b { (a, b) } else { (b, a) };
        let idx = i * nat + j;
        if !seen[idx] {
            seen[idx] = true;
            bonds.push([i, j]);
        }
        Ok(())
    };

    for i in 0..nat.saturating_sub(1) {
        for j in (i + 1)..nat {
            let cutoff = 1.25 * (covalent_radius_bohr(&atoms[i]) + covalent_radius_bohr(&atoms[j]));
            let d_is = distance_between(q_is, i, j);
            let d_fs = distance_between(q_fs, i, j);
            if d_is <= cutoff || d_fs <= cutoff {
                push_bond(i, j, &mut bonds)?;
            }
        }
    }
    for &[i, j] in extra_bonds {
        push_bond(i, j, &mut bonds)?;
    }
    if bonds.is_empty() {
        bail!("RDA internal conditional optimization could not build a bond graph");
    }
    internal_coords_from_bonds(nat, &bonds)
}

fn validate_atom_indices(indices: &[usize], nat: usize, label: &str) -> Result<()> {
    if indices.is_empty() {
        bail!("RDA {label} must not be empty");
    }
    let mut seen = vec![false; nat];
    for &i in indices {
        if i >= nat {
            bail!("RDA {label} atom index {i} out of range 0..{}", nat - 1);
        }
        if seen[i] {
            bail!("RDA {label} contains duplicate atom {i}");
        }
        seen[i] = true;
    }
    Ok(())
}

fn centroid(q: &[f64], active: &[usize]) -> [f64; 3] {
    let mut c = [0.0; 3];
    let n = active.len() as f64;
    for &i in active {
        for k in 0..3 {
            c[k] += q[3 * i + k] / n;
        }
    }
    c
}

fn center_align(mobile: &[f64], reference: &[f64], active: &[usize]) -> Vec<f64> {
    let cm_m = centroid(mobile, active);
    let cm_r = centroid(reference, active);
    let shift = [cm_r[0] - cm_m[0], cm_r[1] - cm_m[1], cm_r[2] - cm_m[2]];
    let mut out = mobile.to_vec();
    for i in 0..(mobile.len() / 3) {
        for k in 0..3 {
            out[3 * i + k] += shift[k];
        }
    }
    out
}

fn quaternion_to_rotation(q: [f64; 4]) -> [[f64; 3]; 3] {
    let w = q[0];
    let x = q[1];
    let y = q[2];
    let z = q[3];
    [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - z * w),
            2.0 * (x * z + y * w),
        ],
        [
            2.0 * (x * y + z * w),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - x * w),
        ],
        [
            2.0 * (x * z - y * w),
            2.0 * (y * z + x * w),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ]
}

fn kabsch_align(mobile: &[f64], reference: &[f64], active: &[usize]) -> Result<Vec<f64>> {
    if active.len() < 2 {
        return Ok(center_align(mobile, reference, active));
    }
    let cm_m = centroid(mobile, active);
    let cm_r = centroid(reference, active);
    let mut s = [[0.0; 3]; 3];
    for &i in active {
        let p = [
            mobile[3 * i] - cm_m[0],
            mobile[3 * i + 1] - cm_m[1],
            mobile[3 * i + 2] - cm_m[2],
        ];
        let q = [
            reference[3 * i] - cm_r[0],
            reference[3 * i + 1] - cm_r[1],
            reference[3 * i + 2] - cm_r[2],
        ];
        for a in 0..3 {
            for b in 0..3 {
                s[a][b] += p[a] * q[b];
            }
        }
    }

    let (sxx, sxy, sxz) = (s[0][0], s[0][1], s[0][2]);
    let (syx, syy, syz) = (s[1][0], s[1][1], s[1][2]);
    let (szx, szy, szz) = (s[2][0], s[2][1], s[2][2]);
    let mut k = Array2::<f64>::zeros((4, 4));
    k[(0, 0)] = sxx + syy + szz;
    k[(0, 1)] = syz - szy;
    k[(0, 2)] = szx - sxz;
    k[(0, 3)] = sxy - syx;
    k[(1, 0)] = k[(0, 1)];
    k[(1, 1)] = sxx - syy - szz;
    k[(1, 2)] = sxy + syx;
    k[(1, 3)] = szx + sxz;
    k[(2, 0)] = k[(0, 2)];
    k[(2, 1)] = k[(1, 2)];
    k[(2, 2)] = -sxx + syy - szz;
    k[(2, 3)] = syz + szy;
    k[(3, 0)] = k[(0, 3)];
    k[(3, 1)] = k[(1, 3)];
    k[(3, 2)] = k[(2, 3)];
    k[(3, 3)] = -sxx - syy + szz;

    let solver = LinAlg::new(Backend::Auto)?;
    let (evals, evecs) = solver.eigh(&k)?;
    let mut best = 0usize;
    for i in 1..evals.len() {
        if evals[i] > evals[best] {
            best = i;
        }
    }
    let mut quat = [
        evecs[(0, best)],
        evecs[(1, best)],
        evecs[(2, best)],
        evecs[(3, best)],
    ];
    let qnorm =
        (quat[0] * quat[0] + quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]).sqrt();
    if qnorm <= 1.0e-14 || !qnorm.is_finite() {
        return Ok(center_align(mobile, reference, active));
    }
    for x in &mut quat {
        *x /= qnorm;
    }
    let rot = quaternion_to_rotation(quat);

    let mut out = mobile.to_vec();
    for i in 0..(mobile.len() / 3) {
        let p = [
            mobile[3 * i] - cm_m[0],
            mobile[3 * i + 1] - cm_m[1],
            mobile[3 * i + 2] - cm_m[2],
        ];
        for a in 0..3 {
            out[3 * i + a] = rot[a][0] * p[0] + rot[a][1] * p[1] + rot[a][2] * p[2] + cm_r[a];
        }
    }
    Ok(out)
}

#[derive(Debug, Clone)]
enum DistanceFunction {
    Cartesian {
        active_atoms: Vec<usize>,
    },
    ReactiveInternal {
        ic: InternalCoords,
        weights: Option<Vec<f64>>,
    },
}

impl DistanceFunction {
    fn label(&self) -> String {
        match self {
            Self::Cartesian { .. } => "cartesian".to_string(),
            Self::ReactiveInternal { ic, .. } => format!(
                "reactive_internal({} bonds, {} angles, {} dihedrals)",
                ic.bonds.len(),
                ic.angles.len(),
                ic.dihedrals.len()
            ),
        }
    }

    fn distance(&self, a: &[f64], b: &[f64]) -> f64 {
        match self {
            Self::Cartesian { active_atoms } => {
                let mut sum = 0.0;
                for &i in active_atoms {
                    let dx = a[3 * i] - b[3 * i];
                    let dy = a[3 * i + 1] - b[3 * i + 1];
                    let dz = a[3 * i + 2] - b[3 * i + 2];
                    sum += dx * dx + dy * dy + dz * dz;
                }
                (sum / active_atoms.len().max(1) as f64).sqrt()
            }
            Self::ReactiveInternal { ic, weights } => {
                let qa = compute_internals(a, ic);
                let qb = compute_internals(b, ic);
                let dq = diff_internals(&qa, &qb, ic.ndiheds());
                if let Some(weights) = weights {
                    let mut sum = 0.0;
                    for (d, w) in dq.iter().zip(weights.iter()) {
                        sum += w * d * d;
                    }
                    (sum / dq.len().max(1) as f64).sqrt()
                } else {
                    norm(&dq) / (dq.len().max(1) as f64).sqrt()
                }
            }
        }
    }
}

fn make_distance(
    opts: &RdaOptions,
    active_atoms: Vec<usize>,
    nat: usize,
) -> Result<DistanceFunction> {
    match opts.distance_mode {
        RdaDistanceMode::Cartesian => Ok(DistanceFunction::Cartesian { active_atoms }),
        RdaDistanceMode::ReactiveInternal => {
            if opts.reactive_bonds.is_empty()
                && opts.reactive_angles.is_empty()
                && opts.reactive_dihedrals.is_empty()
            {
                bail!("RDA reactive_internal distance requires at least one reactive coordinate");
            }
            let nint = opts.reactive_bonds.len()
                + opts.reactive_angles.len()
                + opts.reactive_dihedrals.len();
            if let Some(weights) = &opts.reactive_coordinate_weights {
                if weights.len() != nint {
                    bail!("RDA reactive_coordinate_weights length must match number of reactive coordinates");
                }
            }
            Ok(DistanceFunction::ReactiveInternal {
                ic: InternalCoords {
                    nat,
                    bonds: opts.reactive_bonds.clone(),
                    angles: opts.reactive_angles.clone(),
                    dihedrals: opts.reactive_dihedrals.clone(),
                },
                weights: opts.reactive_coordinate_weights.clone(),
            })
        }
    }
}

fn classify_direction(
    q_start: &[f64],
    q_copt: &[f64],
    q_is: &[f64],
    q_fs: &[f64],
    distance_fn: &DistanceFunction,
    distance_tol: f64,
) -> (String, f64, f64) {
    let delta_is = distance_fn.distance(q_copt, q_is) - distance_fn.distance(q_start, q_is);
    let delta_fs = distance_fn.distance(q_copt, q_fs) - distance_fn.distance(q_start, q_fs);
    if delta_is.abs() < distance_tol && delta_fs.abs() < distance_tol {
        return ("ND".to_string(), delta_is, delta_fs);
    }
    if delta_is * delta_fs > 0.0 {
        return ("ND".to_string(), delta_is, delta_fs);
    }
    if delta_is < 0.0 && delta_fs > 0.0 {
        return ("IS".to_string(), delta_is, delta_fs);
    }
    if delta_is > 0.0 && delta_fs < 0.0 {
        return ("FS".to_string(), delta_is, delta_fs);
    }
    ("ND".to_string(), delta_is, delta_fs)
}

fn identity(n: usize) -> Array2<f64> {
    let mut out = Array2::<f64>::zeros((n, n));
    for i in 0..n {
        out[(i, i)] = 1.0;
    }
    out
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

fn bfgs_inverse_update(h_inv: &mut Array2<f64>, step: &[f64], y: &[f64]) {
    let ys = dot(y, step);
    if ys <= 1.0e-14 {
        return;
    }
    let n = step.len();
    let rho = 1.0 / ys;
    let old = h_inv.clone();
    let hy = mat_vec(&old, y);
    let yhy = dot(y, &hy);
    for i in 0..n {
        for j in 0..n {
            h_inv[(i, j)] = old[(i, j)] + (1.0 + yhy * rho) * rho * step[i] * step[j]
                - rho * (step[i] * hy[j] + hy[i] * step[j]);
        }
    }
}

fn conditional_cartesian_optimize<O: Objective>(
    obj: &mut O,
    q0: &[f64],
    energy_tol: f64,
    maxstep: usize,
    step_max_component: f64,
    verbosity: u8,
) -> Result<(Vec<f64>, f64, usize, String)> {
    let mut q = q0.to_vec();
    let mut grad = vec![0.0; q.len()];
    let mut energy = obj.energy_gradient(&q, &mut grad)?;
    let mut h_inv = identity(q.len());
    let mut previous_energy = energy;

    for istep in 1..=maxstep {
        if verbosity > 1 {
            println!("    RDA Cartesian conditional step {istep}/{maxstep}");
        }
        let mut step: Vec<f64> = mat_vec(&h_inv, &grad).into_iter().map(|v| -v).collect();
        step_limit(&mut step, step_max_component);

        let mut accepted = false;
        let mut scale = 1.0;
        let mut trial_q = q.clone();
        let mut trial_grad = vec![0.0; q.len()];
        let mut trial_energy = energy;
        while scale >= 1.0e-4 {
            for i in 0..q.len() {
                trial_q[i] = q[i] + scale * step[i];
            }
            trial_energy = obj.energy_gradient(&trial_q, &mut trial_grad)?;
            if trial_energy.is_finite()
                && trial_energy <= energy + 1.0e-4 * scale * dot(&grad, &step)
            {
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if !accepted {
            for i in 0..q.len() {
                trial_q[i] = q[i] + scale * step[i];
            }
            trial_energy = obj.energy_gradient(&trial_q, &mut trial_grad)?;
        }

        let accepted_step: Vec<f64> = trial_q.iter().zip(q.iter()).map(|(a, b)| a - b).collect();
        let y: Vec<f64> = trial_grad
            .iter()
            .zip(grad.iter())
            .map(|(a, b)| a - b)
            .collect();
        bfgs_inverse_update(&mut h_inv, &accepted_step, &y);
        let d_e = trial_energy - previous_energy;
        q = trial_q;
        energy = trial_energy;
        grad = trial_grad;
        if d_e.abs() < energy_tol {
            return Ok((
                q,
                energy,
                istep,
                "RDA conditional Cartesian optimization reached the energy-change threshold."
                    .to_string(),
            ));
        }
        previous_energy = energy;
    }
    Ok((
        q,
        energy,
        maxstep,
        "RDA conditional Cartesian optimization reached maxstep.".to_string(),
    ))
}

fn internal_gradient(bpg_b: &[Vec<f64>], inv_g: &[Vec<f64>], grad_x: &[f64]) -> Vec<f64> {
    let nint = bpg_b.len();
    let mut b_grad = vec![0.0; nint];
    for i in 0..nint {
        b_grad[i] = bpg_b[i].iter().zip(grad_x.iter()).map(|(b, g)| b * g).sum();
    }
    let mut grad_q = vec![0.0; nint];
    for i in 0..nint {
        for j in 0..nint {
            grad_q[i] += inv_g[i][j] * b_grad[j];
        }
    }
    grad_q
}

fn solve_linear(a: &Array2<f64>, b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut aug = vec![vec![0.0; n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[(i, j)];
        }
        aug[i][n] = b[i];
    }
    for col in 0..n {
        let mut pivot = col;
        let mut best = aug[col][col].abs();
        for row in (col + 1)..n {
            if aug[row][col].abs() > best {
                best = aug[row][col].abs();
                pivot = row;
            }
        }
        if best < 1.0e-12 {
            return b.to_vec();
        }
        if pivot != col {
            aug.swap(pivot, col);
        }
        let div = aug[col][col];
        for j in col..=n {
            aug[col][j] /= div;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row][col];
            for j in col..=n {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }
    aug.into_iter().map(|row| row[n]).collect()
}

fn bfgs_hessian_update(hess: &Array2<f64>, step: &[f64], y: &[f64]) -> Array2<f64> {
    let ys = dot(y, step);
    if ys <= 1.0e-12 {
        return hess.clone();
    }
    let hs = mat_vec(hess, step);
    let shs = dot(step, &hs);
    if shs.abs() <= 1.0e-12 {
        return hess.clone();
    }
    let n = step.len();
    let mut out = hess.clone();
    for i in 0..n {
        for j in 0..n {
            out[(i, j)] += y[i] * y[j] / ys - hs[i] * hs[j] / shs;
        }
    }
    out
}

fn conditional_internal_optimize<O: Objective>(
    obj: &mut O,
    q0: &[f64],
    ic: &InternalCoords,
    energy_tol: f64,
    maxstep: usize,
    max_step_internal: f64,
    best_fit_iters: usize,
    verbosity: u8,
) -> Result<(Vec<f64>, f64, usize, String)> {
    if ic.nint() == 0 {
        bail!("RDA internal conditional optimization found no internal coordinates");
    }
    let mut q = q0.to_vec();
    let mut qs = compute_internals(&q, ic);
    let mut bpg = bpg_matrix(&q, ic)?;
    let mut grad_x = vec![0.0; q.len()];
    let mut energy = obj.energy_gradient(&q, &mut grad_x)?;
    let mut grad_q = internal_gradient(&bpg.b, &bpg.inv_g, &grad_x);
    let mut hess_q = identity(ic.nint());
    let mut previous_energy = energy;

    for istep in 1..=maxstep {
        if verbosity > 1 {
            println!("    RDA internal conditional step {istep}/{maxstep}");
        }
        let rhs: Vec<f64> = grad_q.iter().map(|g| -g).collect();
        let mut dq = solve_linear(&hess_q, &rhs);
        step_limit(&mut dq, max_step_internal);

        let mut accepted = false;
        let mut scale = 1.0;
        let mut trial_q = q.clone();
        let mut trial_grad_x = vec![0.0; q.len()];
        let mut trial_energy = energy;
        let mut accepted_dq = dq.clone();
        while scale >= 1.0e-4 {
            let scaled_dq: Vec<f64> = dq.iter().map(|v| scale * v).collect();
            trial_q = q.clone();
            best_fit_dq_to_cart(&mut trial_q, &qs, &scaled_dq, ic, &bpg, best_fit_iters)?;
            trial_energy = obj.energy_gradient(&trial_q, &mut trial_grad_x)?;
            if trial_energy.is_finite()
                && trial_energy <= energy + 1.0e-4 * scale * dot(&grad_q, &dq)
            {
                accepted = true;
                accepted_dq = scaled_dq;
                break;
            }
            scale *= 0.5;
        }
        if !accepted {
            accepted_dq = dq.iter().map(|v| scale * v).collect();
            trial_q = q.clone();
            best_fit_dq_to_cart(&mut trial_q, &qs, &accepted_dq, ic, &bpg, best_fit_iters)?;
            trial_energy = obj.energy_gradient(&trial_q, &mut trial_grad_x)?;
        }

        let trial_bpg = bpg_matrix(&trial_q, ic)?;
        let trial_grad_q = internal_gradient(&trial_bpg.b, &trial_bpg.inv_g, &trial_grad_x);
        let y: Vec<f64> = trial_grad_q
            .iter()
            .zip(grad_q.iter())
            .map(|(a, b)| a - b)
            .collect();
        hess_q = bfgs_hessian_update(&hess_q, &accepted_dq, &y);
        let d_e = trial_energy - previous_energy;
        q = trial_q;
        qs = compute_internals(&q, ic);
        bpg = trial_bpg;
        energy = trial_energy;
        grad_x = trial_grad_x;
        grad_q = trial_grad_q;
        if d_e.abs() < energy_tol {
            return Ok((
                q,
                energy,
                istep,
                "RDA conditional internal optimization reached the energy-change threshold."
                    .to_string(),
            ));
        }
        previous_energy = energy;
    }
    Ok((
        q,
        energy,
        maxstep,
        "RDA conditional internal optimization reached maxstep.".to_string(),
    ))
}

fn conditional_optimize<O: Objective>(
    obj: &mut O,
    q0: &[f64],
    internal_coordinates: Option<&InternalCoords>,
    opts: &RdaOptions,
    energy_tol: f64,
) -> Result<(Vec<f64>, f64, usize, String)> {
    match opts.conditional_coordinates {
        RdaConditionalCoordinates::Cartesian => conditional_cartesian_optimize(
            obj,
            q0,
            energy_tol,
            opts.max_conditional_steps,
            opts.cartesian_step_max_component,
            opts.verbosity,
        ),
        RdaConditionalCoordinates::Internal => {
            let Some(ic) = internal_coordinates else {
                if opts.fallback_to_cartesian {
                    if opts.verbosity > 0 {
                        println!("  RDA internal conditional coordinates unavailable; falling back to Cartesian.");
                    }
                    return conditional_cartesian_optimize(
                        obj,
                        q0,
                        energy_tol,
                        opts.max_conditional_steps,
                        opts.cartesian_step_max_component,
                        opts.verbosity,
                    );
                }
                bail!("RDA internal conditional coordinates unavailable");
            };
            match conditional_internal_optimize(
                obj,
                q0,
                ic,
                energy_tol,
                opts.max_conditional_steps,
                opts.internal_step_max_component,
                opts.internal_best_fit_iters,
                opts.verbosity,
            ) {
                Ok(result) => Ok(result),
                Err(err) if opts.fallback_to_cartesian => {
                    if opts.verbosity > 0 {
                        println!("  RDA internal conditional optimization failed ({err}); falling back to Cartesian.");
                    }
                    conditional_cartesian_optimize(
                        obj,
                        q0,
                        energy_tol,
                        opts.max_conditional_steps,
                        opts.cartesian_step_max_component,
                        opts.verbosity,
                    )
                }
                Err(err) => Err(err),
            }
        }
    }
}

fn make_point(
    label: &str,
    q_start: Vec<f64>,
    q_copt: Vec<f64>,
    energy: f64,
    nsteps: usize,
    q_is: &[f64],
    q_fs: &[f64],
    distance_fn: &DistanceFunction,
    distance_tol: f64,
    beta: Option<f64>,
    gamma: Option<f64>,
    message: String,
) -> RdaPoint {
    let (direction, delta_d_is, delta_d_fs) =
        classify_direction(&q_start, &q_copt, q_is, q_fs, distance_fn, distance_tol);
    RdaPoint {
        label: label.to_string(),
        q_start,
        q_copt,
        energy,
        nsteps,
        direction,
        delta_d_is,
        delta_d_fs,
        beta,
        gamma,
        message,
    }
}

fn select_bracket_candidate(
    q_a: &[f64],
    q_b: &[f64],
    q_is: &[f64],
    q_fs: &[f64],
    distance_fn: &DistanceFunction,
    gamma_values: &[f64],
) -> (Vec<f64>, f64) {
    let mut best_q = interpolate(q_a, q_b, gamma_values.first().copied().unwrap_or(0.5));
    let mut best_gamma = gamma_values.first().copied().unwrap_or(0.5);
    let mut best_score = f64::INFINITY;
    for &gamma in gamma_values {
        let q_trial = interpolate(q_a, q_b, gamma);
        let score =
            (distance_fn.distance(&q_trial, q_is) - distance_fn.distance(&q_trial, q_fs)).abs();
        if score < best_score {
            best_score = score;
            best_q = q_trial;
            best_gamma = gamma;
        }
    }
    (best_q, best_gamma)
}

fn format_opt(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{x:.3}"),
        None => "-".to_string(),
    }
}

fn print_rda_table(points: &[RdaPoint]) {
    println!("RDA summary:");
    println!(
        "{:>16} {:>8} {:>8} {:>16} {:>12} {:>12} {:>8}",
        "label", "beta", "gamma", "E[Eh]", "dIS", "dFS", "dir"
    );
    for p in points {
        println!(
            "{:>16} {:>8} {:>8} {:16.8} {:12.4e} {:12.4e} {:>8}",
            p.label,
            format_opt(p.beta),
            format_opt(p.gamma),
            p.energy,
            p.delta_d_is,
            p.delta_d_fs,
            p.direction
        );
    }
}

fn write_xyz_frame(out: &mut String, atoms: &[String], q: &[f64], comment: &str) {
    out.push_str(&xyz_string(atoms, q, comment));
}

fn write_rda_trajectory(
    filename: &Option<String>,
    atoms: &[String],
    q_is: &[f64],
    q_fs: &[f64],
    points: &[RdaPoint],
) -> Result<()> {
    let Some(filename) = filename else {
        return Ok(());
    };
    let mut out = String::new();
    write_xyz_frame(&mut out, atoms, q_is, "label=endpoint_IS");
    for p in points {
        write_xyz_frame(
            &mut out,
            atoms,
            &p.q_start,
            &format!("label={}_start direction={} beta={} gamma={} delta_d_is={:.10e} delta_d_fs={:.10e}", p.label, p.direction, format_opt(p.beta), format_opt(p.gamma), p.delta_d_is, p.delta_d_fs),
        );
        if p.nsteps > 0 {
            write_xyz_frame(
                &mut out,
                atoms,
                &p.q_copt,
                &format!(
                    "label={}_copt direction={} energy={:.12} nsteps={}",
                    p.label, p.direction, p.energy, p.nsteps
                ),
            );
        }
    }
    write_xyz_frame(&mut out, atoms, q_fs, "label=endpoint_FS");
    std::fs::write(filename, out)?;
    Ok(())
}

fn write_rda_chain(
    filename: &Option<String>,
    atoms: &[String],
    q_is: &[f64],
    q_fs: &[f64],
    points: &[RdaPoint],
    q_guess: &[f64],
) -> Result<()> {
    let Some(filename) = filename else {
        return Ok(());
    };
    let mut out = String::new();
    write_xyz_frame(&mut out, atoms, q_is, "label=endpoint_IS image=0");
    let mut image = 1usize;
    for p in points {
        if p.nsteps > 0 || p.label == "gamma_selected" {
            write_xyz_frame(
                &mut out,
                atoms,
                &p.q_copt,
                &format!("label={} image={image}", p.label),
            );
            image += 1;
        }
    }
    write_xyz_frame(
        &mut out,
        atoms,
        q_guess,
        &format!("label=rda_qts image={image}"),
    );
    write_xyz_frame(
        &mut out,
        atoms,
        q_fs,
        &format!("label=endpoint_FS image={}", image + 1),
    );
    std::fs::write(filename, out)?;
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ReactiveCoordinateMode {
    Bond,
    Angle,
    Dihedral,
}

#[derive(Debug, Clone)]
struct ReactiveCoordinateRef {
    internal_index: usize,
    mode: ReactiveCoordinateMode,
    coordinate: ConstraintCoordinate,
    label: String,
}

#[derive(Debug, Clone)]
struct PoorMansNebImage {
    index: usize,
    alpha: f64,
    q_start: Vec<f64>,
    q: Vec<f64>,
    energy: f64,
    converged: bool,
    nsteps: usize,
    message: String,
    reactive_values: Vec<f64>,
}

fn same_bond(a: [usize; 2], b: [usize; 2]) -> bool {
    a == b || a == [b[1], b[0]]
}

fn same_angle(a: [usize; 3], b: [usize; 3]) -> bool {
    a == b || a == [b[2], b[1], b[0]]
}

fn same_dihedral(a: [usize; 4], b: [usize; 4]) -> bool {
    a == b || a == [b[3], b[2], b[1], b[0]]
}

fn wrap_periodic_angle(mut x: f64) -> f64 {
    while x > std::f64::consts::PI {
        x -= 2.0 * std::f64::consts::PI;
    }
    while x <= -std::f64::consts::PI {
        x += 2.0 * std::f64::consts::PI;
    }
    x
}

fn find_reactive_coordinate_indices(
    ic: &InternalCoords,
    reactive_bonds: &[[usize; 2]],
    reactive_angles: &[[usize; 3]],
    reactive_dihedrals: &[[usize; 4]],
) -> Result<Vec<ReactiveCoordinateRef>> {
    let mut refs = Vec::new();

    for &target in reactive_bonds {
        let Some((idx, &bond)) = ic
            .bonds
            .iter()
            .enumerate()
            .find(|(_, &bond)| same_bond(bond, target))
        else {
            bail!(
                "RDA poor-man NEB could not find reactive bond {:?} in the internal coordinate set",
                target
            );
        };
        refs.push(ReactiveCoordinateRef {
            internal_index: idx,
            mode: ReactiveCoordinateMode::Bond,
            coordinate: ConstraintCoordinate::Bond(bond),
            label: format!("bond({}-{})", bond[0] + 1, bond[1] + 1),
        });
    }

    let angle_offset = ic.bonds.len();
    for &target in reactive_angles {
        let Some((idx, &angle)) = ic
            .angles
            .iter()
            .enumerate()
            .find(|(_, &angle)| same_angle(angle, target))
        else {
            bail!("RDA poor-man NEB could not find reactive angle {:?} in the internal coordinate set", target);
        };
        refs.push(ReactiveCoordinateRef {
            internal_index: angle_offset + idx,
            mode: ReactiveCoordinateMode::Angle,
            coordinate: ConstraintCoordinate::Angle(angle),
            label: format!("angle({}-{}-{})", angle[0] + 1, angle[1] + 1, angle[2] + 1),
        });
    }

    let dihedral_offset = ic.bonds.len() + ic.angles.len();
    for &target in reactive_dihedrals {
        let Some((idx, &dihedral)) = ic
            .dihedrals
            .iter()
            .enumerate()
            .find(|(_, &dihedral)| same_dihedral(dihedral, target))
        else {
            bail!("RDA poor-man NEB could not find reactive dihedral {:?} in the internal coordinate set", target);
        };
        refs.push(ReactiveCoordinateRef {
            internal_index: dihedral_offset + idx,
            mode: ReactiveCoordinateMode::Dihedral,
            coordinate: ConstraintCoordinate::Dihedral(dihedral),
            label: format!(
                "dihedral({}-{}-{}-{})",
                dihedral[0] + 1,
                dihedral[1] + 1,
                dihedral[2] + 1,
                dihedral[3] + 1
            ),
        });
    }

    Ok(refs)
}

fn interpolate_reactive_targets(
    qs_is: &[f64],
    qs_fs: &[f64],
    reactive_refs: &[ReactiveCoordinateRef],
    alpha: f64,
) -> Vec<ConstraintTarget> {
    reactive_refs
        .iter()
        .map(|r| {
            let start = qs_is[r.internal_index];
            let end = qs_fs[r.internal_index];
            let target = match r.mode {
                ReactiveCoordinateMode::Dihedral => {
                    wrap_periodic_angle(start + alpha * wrap_periodic_angle(end - start))
                }
                ReactiveCoordinateMode::Bond | ReactiveCoordinateMode::Angle => {
                    (1.0 - alpha) * start + alpha * end
                }
            };
            ConstraintTarget {
                coordinate: r.coordinate,
                target,
            }
        })
        .collect()
}

fn write_poormans_neb_chain(
    filename: &Option<String>,
    atoms: &[String],
    images: &[PoorMansNebImage],
    selected_index: usize,
) -> Result<()> {
    let Some(filename) = filename else {
        return Ok(());
    };
    let mut out = String::new();
    for image in images {
        let selected = if image.index == selected_index {
            " selected=true"
        } else {
            ""
        };
        write_xyz_frame(
            &mut out,
            atoms,
            &image.q,
            &format!(
                "label=poormans_NEB_{:03} image={} alpha={:.6} E={:.12} converged={} nsteps={}{}",
                image.index,
                image.index,
                image.alpha,
                image.energy,
                image.converged,
                image.nsteps,
                selected
            ),
        );
    }
    std::fs::write(filename, out)?;
    Ok(())
}

fn write_poormans_neb_profile(
    filename: &Option<String>,
    images: &[PoorMansNebImage],
    labels: &[String],
) -> Result<()> {
    let Some(filename) = filename else {
        return Ok(());
    };
    let e0 = images.first().map(|i| i.energy).unwrap_or(0.0);
    let mut out = String::new();
    out.push_str("# Poor-man NEB constrained-relaxed RDA profile\n");
    out.push_str("# image alpha E_Eh dE_kJ_mol converged nsteps");
    for label in labels {
        out.push(' ');
        out.push_str(label);
    }
    out.push('\n');
    for image in images {
        out.push_str(&format!(
            "{:5} {:10.6} {:18.10} {:14.6} {:>9} {:6}",
            image.index,
            image.alpha,
            image.energy,
            (image.energy - e0) * 2625.499638,
            image.converged,
            image.nsteps
        ));
        for value in &image.reactive_values {
            out.push_str(&format!(" {:14.8}", value));
        }
        out.push('\n');
    }
    std::fs::write(filename, out)?;
    Ok(())
}

fn print_poormans_neb_profile(images: &[PoorMansNebImage], labels: &[String]) {
    let e0 = images.first().map(|i| i.energy).unwrap_or(0.0);
    println!("RDA poor-man NEB profile:");
    println!(
        "{:>5} {:>8} {:>16} {:>14} {:>10} {:>8}",
        "image", "alpha", "E[Eh]", "dE[kJ/mol]", "conv", "steps"
    );
    for image in images {
        println!(
            "{:5} {:8.3} {:16.8} {:14.6} {:>10} {:8}",
            image.index,
            image.alpha,
            image.energy,
            (image.energy - e0) * 2625.499638,
            image.converged,
            image.nsteps
        );
    }
    if !labels.is_empty() {
        println!(
            "RDA poor-man NEB constrained coordinates: {}",
            labels.join(", ")
        );
    }
}

fn run_poormans_neb<O: Objective>(
    obj: &mut O,
    atoms: &[String],
    q_is: &[f64],
    q_fs: &[f64],
    distance_fn: &DistanceFunction,
    opts: &RdaOptions,
) -> Result<(Vec<PoorMansNebImage>, usize, Vec<RdaPoint>, Vec<f64>)> {
    if opts.reactive_bonds.is_empty()
        && opts.reactive_angles.is_empty()
        && opts.reactive_dihedrals.is_empty()
    {
        bail!("RDA poor-man NEB requires reactive_bonds, reactive_angles, or reactive_dihedrals");
    }
    if opts.poormans_neb_images < 3 {
        bail!("RDA poor-man NEB requires at least 3 images");
    }

    let mut extra_bonds = opts.reactive_bonds.clone();
    extra_bonds.extend(opts.extra_bonds.iter().copied());
    for angle in &opts.reactive_angles {
        extra_bonds.push([angle[0], angle[1]]);
        extra_bonds.push([angle[1], angle[2]]);
    }
    for angle in &opts.extra_angles {
        extra_bonds.push([angle[0], angle[1]]);
        extra_bonds.push([angle[1], angle[2]]);
    }
    for dihedral in &opts.reactive_dihedrals {
        extra_bonds.push([dihedral[0], dihedral[1]]);
        extra_bonds.push([dihedral[1], dihedral[2]]);
        extra_bonds.push([dihedral[2], dihedral[3]]);
    }
    for dihedral in &opts.extra_dihedrals {
        extra_bonds.push([dihedral[0], dihedral[1]]);
        extra_bonds.push([dihedral[1], dihedral[2]]);
        extra_bonds.push([dihedral[2], dihedral[3]]);
    }
    let ic = union_endpoint_internal_coordinates(atoms, q_is, q_fs, &extra_bonds)?;
    let qs_is = compute_internals(q_is, &ic);
    let qs_fs = compute_internals(q_fs, &ic);
    let reactive_refs = find_reactive_coordinate_indices(
        &ic,
        &opts.reactive_bonds,
        &opts.reactive_angles,
        &opts.reactive_dihedrals,
    )?;
    let labels: Vec<String> = reactive_refs.iter().map(|r| r.label.clone()).collect();

    if opts.verbosity > 0 {
        println!("RDA poor-man NEB: constrained relaxation of endpoint-interpolated images");
        println!("RDA poor-man NEB: images={}", opts.poormans_neb_images);
        println!("RDA poor-man NEB: constrained coordinates={:?}", labels);
    }

    let mut images = Vec::new();
    for image_index in 0..opts.poormans_neb_images {
        let alpha = image_index as f64 / (opts.poormans_neb_images - 1) as f64;
        let q_start = interpolate(q_is, q_fs, alpha);
        let (q_relaxed, energy, converged, nsteps, message) = if image_index == 0 {
            let mut grad = vec![0.0; q_is.len()];
            let energy = obj.energy_gradient(q_is, &mut grad)?;
            (
                q_is.to_vec(),
                energy,
                true,
                0,
                "Endpoint image.".to_string(),
            )
        } else if image_index == opts.poormans_neb_images - 1 {
            let mut grad = vec![0.0; q_fs.len()];
            let energy = obj.energy_gradient(q_fs, &mut grad)?;
            (
                q_fs.to_vec(),
                energy,
                true,
                0,
                "Endpoint image.".to_string(),
            )
        } else {
            if opts.verbosity > 0 {
                println!(
                    "RDA poor-man NEB: relaxing image {}/{}, alpha={:.3}",
                    image_index,
                    opts.poormans_neb_images - 1,
                    alpha
                );
            }
            let constraints = interpolate_reactive_targets(&qs_is, &qs_fs, &reactive_refs, alpha);
            let mut opt_opts = ConstrainedOptimizationOptions::default();
            opt_opts.max_cycles = opts.poormans_neb_maxiter;
            opt_opts.energy_tol = opts.later_energy_tol_ev * EV_TO_HARTREE;
            opt_opts.max_step_internal = opts.internal_step_max_component;
            opt_opts.best_fit_iters = opts.internal_best_fit_iters;
            opt_opts.coordinate_mode = ConstrainedCoordinateMode::PulayRedundant;
            opt_opts.verbosity = if opts.verbosity > 1 {
                opts.verbosity
            } else {
                0
            };
            let result = optimize_constrained_with_internals(
                obj,
                q_start.clone(),
                ic.clone(),
                constraints,
                opt_opts,
            )?;
            (
                result.x,
                result.energy,
                result.converged,
                result.cycles,
                result.message,
            )
        };
        let qs_relaxed = compute_internals(&q_relaxed, &ic);
        let reactive_values = reactive_refs
            .iter()
            .map(|r| qs_relaxed[r.internal_index])
            .collect();
        images.push(PoorMansNebImage {
            index: image_index,
            alpha,
            q_start,
            q: q_relaxed,
            energy,
            converged,
            nsteps,
            message,
            reactive_values,
        });
    }

    let selected_index = images[1..images.len() - 1]
        .iter()
        .max_by(|a, b| {
            a.energy
                .partial_cmp(&b.energy)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|image| image.index)
        .unwrap_or_else(|| {
            images
                .iter()
                .max_by(|a, b| {
                    a.energy
                        .partial_cmp(&b.energy)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|image| image.index)
                .unwrap_or(0)
        });
    let selected_q = images[selected_index].q.clone();

    let points: Vec<RdaPoint> = images
        .iter()
        .map(|image| {
            make_point(
                &format!("poormans_NEB_{:03}", image.index),
                image.q_start.clone(),
                image.q.clone(),
                image.energy,
                image.nsteps,
                q_is,
                q_fs,
                distance_fn,
                opts.distance_tol,
                None,
                None,
                image.message.clone(),
            )
        })
        .collect();

    write_poormans_neb_chain(&opts.poormans_neb_file, atoms, &images, selected_index)?;
    write_poormans_neb_profile(&opts.poormans_neb_profile_file, &images, &labels)?;
    if opts.verbosity > 0 {
        print_poormans_neb_profile(&images, &labels);
        println!(
            "RDA poor-man NEB: selected highest-energy relaxed interior image {} with E={:.10} Eh",
            selected_index, images[selected_index].energy
        );
        if let Some(file) = &opts.poormans_neb_file {
            println!("RDA poor-man NEB: wrote final relaxed reaction path to {file}");
        }
        if let Some(file) = &opts.poormans_neb_profile_file {
            println!("RDA poor-man NEB: wrote energetics profile to {file}");
        }
    }

    Ok((images, selected_index, points, selected_q))
}

fn finish_result(
    atoms: Vec<String>,
    q: Vec<f64>,
    reaction_direction: Option<Vec<f64>>,
    points: Vec<RdaPoint>,
    direction: String,
    bracketed: bool,
    message: String,
    q_is: &[f64],
    q_fs: &[f64],
    opts: &RdaOptions,
) -> Result<RdaGuessResult> {
    if opts.verbosity > 0 {
        print_rda_table(&points);
        if let Some(file) = &opts.rda_trajectory_file {
            println!("RDA: wrote detailed optimization steps to {file}");
        }
        if let Some(file) = &opts.rda_chain_file {
            println!("RDA: wrote compact quasi-TS search chain to {file}");
        }
    }
    write_rda_trajectory(&opts.rda_trajectory_file, &atoms, q_is, q_fs, &points)?;
    write_rda_chain(&opts.rda_chain_file, &atoms, q_is, q_fs, &points, &q)?;
    Ok(RdaGuessResult {
        atoms,
        q,
        reaction_direction,
        points,
        direction,
        bracketed,
        message,
        rda_trajectory_file: opts.rda_trajectory_file.clone(),
        rda_chain_file: opts.rda_chain_file.clone(),
        poormans_neb_file: opts.poormans_neb_file.clone(),
        poormans_neb_profile_file: opts.poormans_neb_profile_file.clone(),
    })
}

pub fn generate_rda_ts_guess<O: Objective>(
    obj: &mut O,
    atoms: Vec<String>,
    q_initial: Vec<f64>,
    q_final: Vec<f64>,
    opts: RdaOptions,
) -> Result<RdaGuessResult> {
    validate_inputs(&atoms, &q_initial, &q_final)?;
    let nat = atoms.len();
    let active = active_atoms(&opts, nat)?;
    let align_atoms = if let Some(indices) = &opts.align_atoms {
        validate_atom_indices(indices, nat, "align_atoms")?;
        indices.clone()
    } else {
        active.clone()
    };
    let q_is = q_initial;
    let q_fs = if opts.align {
        if opts.verbosity > 0 {
            println!("RDA: aligning final geometry to initial geometry with Kabsch rotation");
        }
        kabsch_align(&q_final, &q_is, &align_atoms)?
    } else {
        q_final
    };
    let distance_fn = make_distance(&opts, active.clone(), nat)?;
    let distance_used = distance_fn.label();
    let eps1 = opts.first_energy_tol_ev * EV_TO_HARTREE;
    let eps2 = opts.later_energy_tol_ev * EV_TO_HARTREE;

    if opts.verbosity > 0 {
        println!("RDA: starting quasi-TS guess generation");
        println!("RDA: atoms={nat} active_atoms={active:?}");
        println!("RDA: distance analysis uses {distance_used} coordinates");
        println!(
            "RDA: conditional optimization coordinates={:?}",
            opts.conditional_coordinates
        );
        let mut tmp_grad = vec![0.0; q_is.len()];
        let e_is = obj.energy_gradient(&q_is, &mut tmp_grad)?;
        let e_fs = obj.energy_gradient(&q_fs, &mut tmp_grad)?;
        println!("RDA endpoint diagnostics:");
        println!("  E(IS) [Eh] = {e_is:.10}");
        println!("  E(FS) [Eh] = {e_fs:.10}");
        println!("  dE FS-IS [Eh] = {:.10}", e_fs - e_is);
        println!(
            "  active endpoint distance [{distance_used}] = {:.6e}",
            distance_fn.distance(&q_is, &q_fs)
        );
    }

    if opts.poormans_neb {
        let (_images, _selected_index, points, q_guess) =
            run_poormans_neb(obj, &atoms, &q_is, &q_fs, &distance_fn, &opts)?;
        return finish_result(
            atoms,
            q_guess,
            normalize(
                &q_fs
                    .iter()
                    .zip(q_is.iter())
                    .map(|(a, b)| a - b)
                    .collect::<Vec<_>>(),
            ),
            points,
            "poormans_NEB".to_string(),
            false,
            "RDA poor-man NEB selected the highest-energy constrained-relaxed interior image as the quasi-TS guess.".to_string(),
            &q_is,
            &q_fs,
            &opts,
        );
    }

    let internal_coordinates =
        if opts.conditional_coordinates == RdaConditionalCoordinates::Internal {
            Some(union_endpoint_internal_coordinates(
                &atoms,
                &q_is,
                &q_fs,
                &opts.reactive_bonds,
            )?)
        } else {
            None
        };

    let mut points = Vec::new();
    let midpoint_alphas = [0.5, 0.45, 0.55, 0.40, 0.60, 0.35, 0.65, 0.30, 0.70];
    let mut midpoint_result = None;
    let mut last_midpoint_error = None;
    for alpha in midpoint_alphas {
        let q_trial = interpolate(&q_is, &q_fs, alpha);
        if opts.verbosity > 0 {
            println!("RDA: conditional optimization of midpoint alpha={alpha:.2}");
        }
        match conditional_optimize(obj, &q_trial, internal_coordinates.as_ref(), &opts, eps1) {
            Ok((q_copt, energy, nsteps, message)) => {
                midpoint_result = Some((alpha, q_trial, q_copt, energy, nsteps, message));
                break;
            }
            Err(err) => {
                last_midpoint_error = Some(format!("{err:#}"));
                if opts.verbosity > 0 {
                    println!("RDA: midpoint alpha={alpha:.2} failed: {err}");
                }
            }
        }
    }
    let Some((alpha_mid, q_mid, q_mid_copt, e_mid, n_mid, msg_mid)) = midpoint_result else {
        bail!(
            "RDA midpoint conditional optimization failed for all alpha candidates; last error: {}",
            last_midpoint_error.unwrap_or_else(|| "unknown".to_string())
        );
    };
    let mid_point = make_point(
        "midpoint",
        q_mid,
        q_mid_copt.clone(),
        e_mid,
        n_mid,
        &q_is,
        &q_fs,
        &distance_fn,
        opts.distance_tol,
        None,
        None,
        format!("{msg_mid} alpha={alpha_mid:.2} distance={distance_used}"),
    );
    if opts.verbosity > 0 {
        println!(
            "RDA: midpoint finished direction={} dIS={:.6e} dFS={:.6e} steps={}",
            mid_point.direction, mid_point.delta_d_is, mid_point.delta_d_fs, mid_point.nsteps
        );
    }
    points.push(mid_point.clone());
    if mid_point.direction == "ND" {
        if opts.verbosity > 0 {
            println!("RDA: ended, midpoint accepted as quasi-TS guess");
        }
        return finish_result(
            atoms,
            q_mid_copt,
            normalize(
                &q_fs
                    .iter()
                    .zip(q_is.iter())
                    .map(|(a, b)| a - b)
                    .collect::<Vec<_>>(),
            ),
            points,
            "ND".to_string(),
            false,
            "RDA accepted the conditionally optimized midpoint as a quasi-TS guess.".to_string(),
            &q_is,
            &q_fs,
            &opts,
        );
    }

    let reference = if mid_point.direction == "IS" {
        &q_fs
    } else {
        &q_is
    };
    let mut previous_direction = mid_point.direction.clone();
    let mut previous_copt = q_mid_copt.clone();
    let mut bracket_pair: Option<(Vec<f64>, Vec<f64>)> = None;
    let mut beta = opts.beta_start;
    for attempt in 1..=opts.max_bracket_attempts {
        let beta_clamped = beta.clamp(0.0, 1.0);
        let q_beta = interpolate(&q_mid_copt, reference, beta_clamped);
        if opts.verbosity > 0 {
            println!(
                "RDA: conditional optimization of beta candidate {attempt}, beta={beta_clamped:.3}"
            );
        }
        let (q_beta_copt, e_beta, n_beta, msg_beta) =
            match conditional_optimize(obj, &q_beta, internal_coordinates.as_ref(), &opts, eps2) {
                Ok(result) => result,
                Err(err) => {
                    if opts.verbosity > 0 {
                        println!("RDA: beta candidate {attempt} failed: {err}");
                    }
                    beta += opts.beta_step;
                    continue;
                }
            };
        let beta_point = make_point(
            &format!("beta_{attempt}"),
            q_beta,
            q_beta_copt.clone(),
            e_beta,
            n_beta,
            &q_is,
            &q_fs,
            &distance_fn,
            opts.distance_tol,
            Some(beta_clamped),
            None,
            format!("{msg_beta} distance={distance_used}"),
        );
        if opts.verbosity > 0 {
            println!(
                "RDA: beta candidate finished direction={} dIS={:.6e} dFS={:.6e} steps={}",
                beta_point.direction,
                beta_point.delta_d_is,
                beta_point.delta_d_fs,
                beta_point.nsteps
            );
        }
        if beta_point.direction == "ND" {
            points.push(beta_point);
            if opts.verbosity > 0 {
                println!("RDA: ended, nondirectional beta candidate accepted as quasi-TS guess");
            }
            return finish_result(
                atoms,
                q_beta_copt.clone(),
                normalize(
                    &q_beta_copt
                        .iter()
                        .zip(previous_copt.iter())
                        .map(|(a, b)| a - b)
                        .collect::<Vec<_>>(),
                ),
                points,
                "ND".to_string(),
                false,
                "RDA accepted a nondirectional beta candidate as a quasi-TS guess.".to_string(),
                &q_is,
                &q_fs,
                &opts,
            );
        }
        if beta_point.direction != previous_direction {
            bracket_pair = Some((previous_copt, q_beta_copt));
            points.push(beta_point);
            break;
        }
        previous_copt = q_beta_copt;
        previous_direction = beta_point.direction.clone();
        points.push(beta_point);
        beta += opts.beta_step;
    }

    let Some((q_a, q_b)) = bracket_pair else {
        if opts.verbosity > 0 {
            println!("RDA: ended without bracket; returning last candidate");
        }
        let q = points
            .last()
            .expect("RDA has midpoint point")
            .q_copt
            .clone();
        let direction = points
            .last()
            .expect("RDA has midpoint point")
            .direction
            .clone();
        return finish_result(
            atoms,
            q.clone(),
            normalize(&q.iter().zip(q_mid_copt.iter()).map(|(a, b)| a - b).collect::<Vec<_>>()),
            points,
            direction,
            false,
            "RDA did not find opposite directionality before max_bracket_attempts; returning the last candidate.".to_string(),
            &q_is,
            &q_fs,
            &opts,
        );
    };

    if opts.verbosity > 0 {
        println!("RDA: opposite directionality found; selecting gamma quasi-TS candidate");
    }
    let (q_guess, gamma) =
        select_bracket_candidate(&q_a, &q_b, &q_is, &q_fs, &distance_fn, &opts.gamma_values);
    let mut grad = vec![0.0; q_guess.len()];
    let e_guess = obj.energy_gradient(&q_guess, &mut grad)?;
    let gamma_point = make_point(
        "gamma_selected",
        q_guess.clone(),
        q_guess.clone(),
        e_guess,
        0,
        &q_is,
        &q_fs,
        &distance_fn,
        opts.distance_tol,
        None,
        Some(gamma),
        format!("Selected geometrically between opposite RDA directions; distance={distance_used}"),
    );
    points.push(gamma_point);
    if opts.verbosity > 0 {
        println!("RDA: ended with bracketed quasi-TS guess, gamma={gamma:.3}");
    }
    finish_result(
        atoms,
        q_guess,
        normalize(
            &q_b.iter()
                .zip(q_a.iter())
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>(),
        ),
        points,
        "bracketed".to_string(),
        true,
        "RDA found opposite directionality and generated a bracketed quasi-TS guess.".to_string(),
        &q_is,
        &q_fs,
        &opts,
    )
}
