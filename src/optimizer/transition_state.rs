//! Transition-state optimization with Smite-style partitioned RFO.
//!
//! This module is intentionally separate from the minimum-search optimizers.
//! It ports the practical Smite TS optimizer structure: Cartesian/internal/Pulay
//! coordinates, partitioned RFO mode following, Bofill Hessian updates, optional
//! model Hessians, trust radii, Hessian repair/recalculation controls, and the
//! same convergence table display.

use anyhow::{bail, Result};
use ndarray::{Array1, Array2};

use crate::normalmode::linalg_shim::{Backend, LinAlg};
use crate::optimizer::internal_coords::{
    analyze_structure, best_fit_dq_to_cart, bpg_matrix, compute_internals,
    internal_coords_from_bonds, ConnectivityModel, InternalCoords,
};
use crate::optimizer::redundant_internals::{
    build_redundant_internals, cartesian_gradient_to_redundant, generalized_inverse_symmetric,
    initial_redundant_hessian, update_cartesian_from_redundant_step, RedundantInternalSystem,
};
use crate::optimizer::traits::Objective;

const HARTREE_TO_KJMOL: f64 = 2625.499638;
pub const STANDARD_TS_ENERGY_TOL: f64 = 5.0e-6;
pub const STANDARD_TS_MAX_STEP: f64 = 4.0e-3;
pub const STANDARD_TS_RMS_STEP: f64 = 2.0e-3;
pub const STANDARD_TS_MAX_GRADIENT: f64 = 3.0e-4;
pub const STANDARD_TS_RMS_GRADIENT: f64 = 1.0e-4;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransitionStateCoordinateMode {
    Cartesian,
    Internal,
    PulayRedundant,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransitionStateInitialHessian {
    Exact,
    Model,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransitionStateModeTracking {
    Internal,
    Cartesian,
}

#[derive(Debug, Clone, PartialEq)]
pub enum TransitionStateReactionMode {
    Lowest,
    Direction(Vec<f64>),
    Bond([usize; 2]),
    Angle([usize; 3]),
    Dihedral([usize; 4]),
    Transfer([usize; 3]),
    Coordinates(Vec<WeightedInternalReference>),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InternalReferenceKind {
    Bond,
    Angle,
    Dihedral,
}

#[derive(Debug, Clone, PartialEq)]
pub struct WeightedInternalReference {
    pub kind: InternalReferenceKind,
    pub atoms: Vec<usize>,
    pub weight: f64,
}

impl WeightedInternalReference {
    pub fn bond(atoms: [usize; 2], weight: f64) -> Self {
        Self {
            kind: InternalReferenceKind::Bond,
            atoms: atoms.to_vec(),
            weight,
        }
    }

    pub fn angle(atoms: [usize; 3], weight: f64) -> Self {
        Self {
            kind: InternalReferenceKind::Angle,
            atoms: atoms.to_vec(),
            weight,
        }
    }

    pub fn dihedral(atoms: [usize; 4], weight: f64) -> Self {
        Self {
            kind: InternalReferenceKind::Dihedral,
            atoms: atoms.to_vec(),
            weight,
        }
    }
}

#[derive(Debug, Clone)]
pub struct TransitionStateOptions {
    pub coordinate_mode: TransitionStateCoordinateMode,
    pub max_cycles: usize,
    pub energy_tol: f64,
    pub max_step: f64,
    pub rms_step: f64,
    pub max_gradient: f64,
    pub rms_gradient: f64,
    pub max_step_internal: f64,
    pub best_fit_iters: usize,
    pub best_fit_rms_tol: f64,
    pub trajectory_file: Option<String>,
    pub hessian_step: f64,
    pub hessian_recalc_interval: usize,
    pub trust_radius: Option<f64>,
    pub trust_radius_min: f64,
    pub trust_radius_max: Option<f64>,
    pub reaction_mode: TransitionStateReactionMode,
    pub max_rejected_steps: usize,
    pub project_eckart: bool,
    pub masses_amu: Option<Vec<f64>>,
    pub final_hessian: bool,
    pub min_gradient_improvement: f64,
    pub adaptive_hessian_recalc: bool,
    pub adaptive_hessian_overlap_min: f64,
    pub adaptive_hessian_rho_min: f64,
    pub adaptive_hessian_rho_max: f64,
    pub adaptive_hessian_trust_fraction: f64,
    pub adaptive_hessian_on_negative_mode_change: bool,
    pub adaptive_hessian_on_trust_collapse: bool,
    pub adaptive_hessian_retry_recalc: bool,
    pub adaptive_hessian_recalc_cooldown: usize,
    pub skip_hessian_recalc_near_convergence: bool,
    pub hessian_recalc_near_convergence_factor: f64,
    pub repair_ts_hessian: bool,
    pub ts_hessian_eigenvalue_floor: f64,
    pub internal_hessian_model: String,
    pub ts_initial_hessian: TransitionStateInitialHessian,
    pub ts_model_negative_curvature: f64,
    pub mode_tracking_coordinates: TransitionStateModeTracking,
    pub g_cutoff: f64,
    pub redundant_penalty: f64,
    pub extra_bonds: Vec<[usize; 2]>,
    pub verbosity: u8,
}

impl Default for TransitionStateOptions {
    fn default() -> Self {
        Self {
            coordinate_mode: TransitionStateCoordinateMode::PulayRedundant,
            max_cycles: 100,
            energy_tol: STANDARD_TS_ENERGY_TOL,
            max_step: STANDARD_TS_MAX_STEP,
            rms_step: STANDARD_TS_RMS_STEP,
            max_gradient: STANDARD_TS_MAX_GRADIENT,
            rms_gradient: STANDARD_TS_RMS_GRADIENT,
            max_step_internal: 0.2,
            best_fit_iters: 5,
            best_fit_rms_tol: 1.0e-7,
            trajectory_file: Some("tsopt_traj.xyz".to_string()),
            hessian_step: 1.0e-3,
            hessian_recalc_interval: 0,
            trust_radius: None,
            trust_radius_min: 1.0e-4,
            trust_radius_max: None,
            reaction_mode: TransitionStateReactionMode::Lowest,
            max_rejected_steps: 8,
            project_eckart: true,
            masses_amu: None,
            final_hessian: true,
            min_gradient_improvement: 0.0,
            adaptive_hessian_recalc: false,
            adaptive_hessian_overlap_min: 0.55,
            adaptive_hessian_rho_min: 0.05,
            adaptive_hessian_rho_max: 5.0,
            adaptive_hessian_trust_fraction: 0.15,
            adaptive_hessian_on_negative_mode_change: true,
            adaptive_hessian_on_trust_collapse: false,
            adaptive_hessian_retry_recalc: false,
            adaptive_hessian_recalc_cooldown: 3,
            skip_hessian_recalc_near_convergence: true,
            hessian_recalc_near_convergence_factor: 3.0,
            repair_ts_hessian: true,
            ts_hessian_eigenvalue_floor: 1.0e-4,
            internal_hessian_model: "simple".to_string(),
            ts_initial_hessian: TransitionStateInitialHessian::Exact,
            ts_model_negative_curvature: 0.05,
            mode_tracking_coordinates: TransitionStateModeTracking::Internal,
            g_cutoff: 1.0e-7,
            redundant_penalty: 1000.0,
            extra_bonds: Vec::new(),
            verbosity: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct TSModeDiagnostics {
    pub index: usize,
    pub eigenvalue: f64,
    pub overlap: Option<f64>,
    pub source: String,
    pub tracking_space: String,
}

#[derive(Debug, Clone)]
pub struct TransitionStateResult {
    pub x: Vec<f64>,
    pub energy: f64,
    pub gradient: Vec<f64>,
    pub converged: bool,
    pub cycles: usize,
    pub message: String,
    pub coordinates: TransitionStateCoordinateMode,
    pub negative_modes: Option<usize>,
    pub reaction_mode: Option<Vec<f64>>,
    pub hessian: Option<Array2<f64>>,
    pub trajectory: Vec<Vec<f64>>,
    pub trajectory_energies: Vec<f64>,
}

#[derive(Debug, Clone)]
enum TSCoordinateState {
    Cartesian,
    Internal {
        ic: InternalCoords,
        q: Vec<f64>,
        b: Vec<Vec<f64>>,
        inv_g: Vec<Vec<f64>>,
    },
    Pulay(RedundantInternalSystem),
}

impl TSCoordinateState {
    fn ndim(&self, x: &[f64]) -> usize {
        match self {
            Self::Cartesian => x.len(),
            Self::Internal { q, .. } => q.len(),
            Self::Pulay(system) => system.coordinates.nint(),
        }
    }

    fn coordinate_label(&self) -> &'static str {
        match self {
            Self::Cartesian => "cartesian",
            Self::Internal { .. } => "internal coordinates",
            Self::Pulay(_) => "Pulay redundant internal coordinates",
        }
    }

    fn internal_coords(&self) -> Option<&InternalCoords> {
        match self {
            Self::Cartesian => None,
            Self::Internal { ic, .. } => Some(ic),
            Self::Pulay(system) => Some(&system.coordinates),
        }
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

fn rms(a: &[f64]) -> f64 {
    if a.is_empty() {
        0.0
    } else {
        norm(a) / (a.len() as f64).sqrt()
    }
}

fn max_abs(a: &[f64]) -> f64 {
    a.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

fn sci2(value: f64) -> String {
    let raw = format!("{value:.2e}");
    let Some((mantissa, exponent)) = raw.split_once('e') else {
        return raw;
    };
    let exponent = exponent.parse::<i32>().unwrap_or(0);
    format!("{mantissa}e{exponent:+03}")
}

fn format_threshold_guide_value(value: f64) -> String {
    sci2(value)
}

fn place_centered_text(line: &mut Vec<char>, center: f64, text: &str) {
    let mut start = (center - (text.len().saturating_sub(1) as f64) / 2.0).round() as isize;
    if start < 0 {
        start = 0;
    }
    let start = start as usize;
    let end = start + text.len();
    if end > line.len() {
        line.extend(std::iter::repeat(' ').take(end - line.len()));
    }
    for (idx, ch) in text.chars().enumerate() {
        line[start + idx] = ch;
    }
}

fn print_ts_threshold_guide(
    energy_tol: f64,
    max_step: f64,
    rms_step: f64,
    max_gradient: f64,
    rms_gradient: f64,
    header_line: &str,
) {
    let values = [
        format_threshold_guide_value(energy_tol * HARTREE_TO_KJMOL),
        format_threshold_guide_value(max_step),
        format_threshold_guide_value(rms_step),
        format_threshold_guide_value(max_gradient),
        format_threshold_guide_value(rms_gradient),
    ];
    let columns = ["dE[kJ/mol]", "max_step", "rms_step", "max_grad", "rms_grad"];
    let mut centers = Vec::with_capacity(columns.len());
    for column in columns {
        if let Some(start) = header_line.find(column) {
            centers.push(start as f64 + (column.len().saturating_sub(1) as f64) / 2.0);
        }
    }
    let mut value_line: Vec<char> = vec![' '; header_line.len()];
    let mut arrow_bar_line: Vec<char> = vec![' '; header_line.len()];
    let mut arrow_line: Vec<char> = vec![' '; header_line.len()];
    let prefix = "Threshold Values:";
    for (idx, ch) in prefix.chars().enumerate() {
        if idx < value_line.len() {
            value_line[idx] = ch;
        }
    }
    for (center, value) in centers.iter().zip(values.iter()) {
        place_centered_text(&mut value_line, *center, value);
        place_centered_text(&mut arrow_bar_line, *center, "|");
        place_centered_text(&mut arrow_line, *center, "V");
    }
    println!("{}", value_line.iter().collect::<String>().trim_end());
    println!("{}", arrow_bar_line.iter().collect::<String>().trim_end());
    println!("{}", arrow_line.iter().collect::<String>().trim_end());
}

fn print_ts_header(opts: &TransitionStateOptions) {
    println!("Transition-state convergence tolerances:");
    println!("  Energy Change:  {} Eh", sci2(opts.energy_tol));
    println!("  Max. Gradient:  {} Eh/bohr", sci2(opts.max_gradient));
    println!("  RMS Gradient:   {} Eh/bohr", sci2(opts.rms_gradient));
    println!("  Max Step:       {} bohr", sci2(opts.max_step));
    println!("  RMS Step:       {} bohr", sci2(opts.rms_step));
    let header_line = format!(
        "{:<5} {:>18} {:>18} {:>11} {:>11} {:>11} {:>11} {:>9} {:>7} {:>8} {:>9}",
        "step",
        "E[Eh]",
        "dE[kJ/mol]",
        "max_step",
        "rms_step",
        "max_grad",
        "rms_grad",
        "step_cap",
        "#imags",
        "TSmode",
        "overlap"
    );
    println!();
    print_ts_threshold_guide(
        opts.energy_tol,
        opts.max_step,
        opts.rms_step,
        opts.max_gradient,
        opts.rms_gradient,
        &header_line,
    );
    println!("{header_line}");
}

fn print_ts_step(
    istep: usize,
    energy: f64,
    energy_change: f64,
    step: &[f64],
    grad: &[f64],
    step_cap: f64,
    negative_modes: usize,
    mode_diag: Option<&TSModeDiagnostics>,
) {
    let mode_text = mode_diag
        .map(|d| d.index.to_string())
        .unwrap_or_else(|| "n/a".to_string());
    let overlap_text = mode_diag
        .and_then(|d| d.overlap)
        .map(|v| format!("{v:.2}"))
        .unwrap_or_else(|| "n/a".to_string());
    println!(
        "{:5} {:18.8} {:18.6} {:>11} {:>11} {:>11} {:>11} {:>9} {:>7} {:>8} {:>9}",
        istep,
        energy,
        energy_change * HARTREE_TO_KJMOL,
        sci2(max_abs(step)),
        sci2(rms(step)),
        sci2(max_abs(grad)),
        sci2(rms(grad)),
        sci2(step_cap),
        negative_modes,
        mode_text,
        overlap_text,
    );
}

fn center_with_fill(text: &str, width: usize, fill: char) -> String {
    if text.len() >= width {
        return text.to_string();
    }
    let pad = width - text.len();
    let left = pad / 2;
    let right = pad - left;
    format!(
        "{}{}{}",
        fill.to_string().repeat(left),
        text,
        fill.to_string().repeat(right)
    )
}

fn print_optimization_status(converged: bool, message: &str) {
    let status = if converged {
        "CONVERGED"
    } else {
        "NOT CONVERGED"
    };
    let status_text = format!(" {status} ");
    let width = (status_text.len() + 8).max(message.len());
    let bar = "=".repeat(width);
    println!();
    println!("{bar}");
    println!("{}", center_with_fill(&status_text, width, '*'));
    if !message.is_empty() {
        println!("{message}");
    }
    println!("{bar}");
    println!();
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

fn vec_from_array1(a: &Array1<f64>) -> Vec<f64> {
    a.iter().copied().collect()
}

fn symmetrize(a: &mut Array2<f64>) {
    let (nr, nc) = a.dim();
    debug_assert_eq!(nr, nc);
    for i in 0..nr {
        for j in 0..i {
            let v = 0.5 * (a[(i, j)] + a[(j, i)]);
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
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

fn outer_scaled(a: &[f64], b: &[f64], scale: f64, out: &mut Array2<f64>) {
    for i in 0..a.len() {
        for j in 0..b.len() {
            out[(i, j)] += scale * a[i] * b[j];
        }
    }
}

fn bofill_update_hessian(hess: &Array2<f64>, step: &[f64], y: &[f64]) -> (Array2<f64>, bool) {
    let h_step = mat_vec(hess, step);
    let u: Vec<f64> = y
        .iter()
        .zip(h_step.iter())
        .map(|(yi, hsi)| yi - hsi)
        .collect();
    let ss = dot(step, step);
    let us = dot(&u, step);
    let uu = dot(&u, &u);
    if ss <= 1.0e-14 || uu <= 1.0e-14 || us.abs() <= 1.0e-14 {
        return (hess.clone(), false);
    }
    let n = step.len();
    let mut updated = hess.clone();
    let phi = (us * us / (uu * ss)).clamp(0.0, 1.0);
    outer_scaled(&u, &u, phi / us, &mut updated);
    let psb_scale = (1.0 - phi) / ss;
    outer_scaled(&u, step, psb_scale, &mut updated);
    outer_scaled(step, &u, psb_scale, &mut updated);
    outer_scaled(step, step, -(1.0 - phi) * us / (ss * ss), &mut updated);
    for i in 0..n {
        for j in 0..i {
            let v = 0.5 * (updated[(i, j)] + updated[(j, i)]);
            updated[(i, j)] = v;
            updated[(j, i)] = v;
        }
    }
    let ok = updated.iter().all(|v| v.is_finite());
    if ok {
        (updated, true)
    } else {
        (hess.clone(), false)
    }
}

fn finite_difference_hessian<O: Objective>(
    obj: &mut O,
    x: &[f64],
    step: f64,
) -> Result<Array2<f64>> {
    if step <= 0.0 {
        bail!("finite-difference Hessian step must be positive");
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

fn eckart_projector(masses_amu: &[f64], coords_bohr: &[f64]) -> Result<Array2<f64>> {
    if coords_bohr.len() % 3 != 0 {
        bail!("Eckart projection requires 3N Cartesian coordinates");
    }
    let nat = coords_bohr.len() / 3;
    if masses_amu.len() != nat {
        bail!(
            "Eckart projection mass count mismatch: got {}, expected {}",
            masses_amu.len(),
            nat
        );
    }
    let mut total_mass = 0.0;
    let mut com = [0.0; 3];
    for a in 0..nat {
        let m = masses_amu[a];
        if m <= 0.0 || !m.is_finite() {
            bail!("Eckart projection requires positive finite atomic masses");
        }
        total_mass += m;
        for k in 0..3 {
            com[k] += m * coords_bohr[3 * a + k];
        }
    }
    if total_mass <= 0.0 {
        bail!("Eckart projection requires nonzero total mass");
    }
    for v in &mut com {
        *v /= total_mass;
    }

    let ndim = 3 * nat;
    let mut b = Array2::<f64>::zeros((ndim, 6));
    for a in 0..nat {
        let sm = masses_amu[a].sqrt();
        let x = coords_bohr[3 * a] - com[0];
        let y = coords_bohr[3 * a + 1] - com[1];
        let z = coords_bohr[3 * a + 2] - com[2];
        let r = 3 * a;

        b[(r, 0)] = sm;
        b[(r + 1, 1)] = sm;
        b[(r + 2, 2)] = sm;

        b[(r + 1, 3)] = sm * z;
        b[(r + 2, 3)] = -sm * y;
        b[(r, 4)] = -sm * z;
        b[(r + 2, 4)] = sm * x;
        b[(r, 5)] = sm * y;
        b[(r + 1, 5)] = -sm * x;
    }

    let mut metric = Array2::<f64>::zeros((6, 6));
    for i in 0..6 {
        for j in 0..6 {
            let mut v = 0.0;
            for r in 0..ndim {
                v += b[(r, i)] * b[(r, j)];
            }
            metric[(i, j)] = v;
        }
    }
    let metric_inv = generalized_inverse_symmetric(&metric, 1.0e-12)?.matrix;

    let mut rb = Array2::<f64>::zeros((ndim, 6));
    for r in 0..ndim {
        for j in 0..6 {
            let mut v = 0.0;
            for k in 0..6 {
                v += b[(r, k)] * metric_inv[(k, j)];
            }
            rb[(r, j)] = v;
        }
    }
    let mut projector = Array2::<f64>::zeros((ndim, ndim));
    for i in 0..ndim {
        projector[(i, i)] = 1.0;
        for j in 0..ndim {
            let mut rot_trans = 0.0;
            for k in 0..6 {
                rot_trans += rb[(i, k)] * b[(j, k)];
            }
            projector[(i, j)] -= rot_trans;
        }
    }
    symmetrize(&mut projector);
    Ok(projector)
}

fn maybe_project_ts_derivatives(
    masses_amu: Option<&[f64]>,
    coords_bohr: &[f64],
    grad: Option<&mut [f64]>,
    hess: Option<&mut Array2<f64>>,
    project_eckart: bool,
) -> Result<()> {
    if !project_eckart {
        return Ok(());
    }
    let masses = masses_amu.ok_or_else(|| {
        anyhow::anyhow!("TS Eckart projection requested but atomic masses are not available")
    })?;
    let nat = masses.len();
    if coords_bohr.len() != 3 * nat {
        bail!("Eckart projection coordinate/mass dimension mismatch");
    }
    let projector = eckart_projector(masses, coords_bohr)?;
    let ndim = coords_bohr.len();
    let sqrt_m: Vec<f64> = masses.iter().flat_map(|m| [m.sqrt(); 3]).collect();

    if let Some(g) = grad {
        if g.len() != ndim {
            bail!("Eckart projection gradient dimension mismatch");
        }
        let mut mw = vec![0.0; ndim];
        for i in 0..ndim {
            mw[i] = g[i] / sqrt_m[i];
        }
        let projected = mat_vec(&projector, &mw);
        for i in 0..ndim {
            g[i] = projected[i] * sqrt_m[i];
        }
    }

    if let Some(h) = hess {
        let (nr, nc) = h.dim();
        if nr != ndim || nc != ndim {
            bail!("Eckart projection Hessian dimension mismatch");
        }
        let mut mw = Array2::<f64>::zeros((ndim, ndim));
        for i in 0..ndim {
            for j in 0..ndim {
                mw[(i, j)] = h[(i, j)] / (sqrt_m[i] * sqrt_m[j]);
            }
        }
        let mut tmp = Array2::<f64>::zeros((ndim, ndim));
        for i in 0..ndim {
            for j in 0..ndim {
                let mut v = 0.0;
                for k in 0..ndim {
                    v += projector[(i, k)] * mw[(k, j)];
                }
                tmp[(i, j)] = v;
            }
        }
        let mut projected = Array2::<f64>::zeros((ndim, ndim));
        for i in 0..ndim {
            for j in 0..ndim {
                let mut v = 0.0;
                for k in 0..ndim {
                    v += tmp[(i, k)] * projector[(k, j)];
                }
                projected[(i, j)] = v * sqrt_m[i] * sqrt_m[j];
            }
        }
        symmetrize(&mut projected);
        *h = projected;
    }
    Ok(())
}

fn internal_gradient_from_bpg(b: &[Vec<f64>], inv_g: &[Vec<f64>], grad_x: &[f64]) -> Vec<f64> {
    let nint = b.len();
    let nx = grad_x.len();
    let mut bg = vec![0.0; nint];
    for i in 0..nint {
        for k in 0..nx {
            bg[i] += b[i][k] * grad_x[k];
        }
    }
    let mut gq = vec![0.0; nint];
    for i in 0..nint {
        for j in 0..nint {
            gq[i] += inv_g[i][j] * bg[j];
        }
    }
    gq
}

fn cartesian_hessian_to_internal(
    hess_x: &Array2<f64>,
    b: &[Vec<f64>],
    inv_g: &[Vec<f64>],
) -> Result<Array2<f64>> {
    let la = LinAlg::new(Backend::Auto)?;
    let bmat = array_from_vecvec(b);
    let invg = array_from_vecvec(inv_g);
    let tmp = la.matmul(&bmat, hess_x);
    let tmp = la.matmul(&tmp, &la.transpose(&bmat));
    let mut hq = la.matmul(&la.matmul(&invg, &tmp), &invg);
    symmetrize(&mut hq);
    Ok(hq)
}

fn initial_model_hessian(ic: &InternalCoords) -> Array2<f64> {
    let nb = ic.nbonds();
    let na = ic.nangles();
    let nd = ic.ndiheds();
    let n = nb + na + nd;
    let mut h = Array2::<f64>::zeros((n, n));
    for i in 0..nb {
        h[(i, i)] = 1.0;
    }
    for i in 0..na {
        h[(nb + i, nb + i)] = 0.25;
    }
    for i in 0..nd {
        h[(nb + na + i, nb + na + i)] = 0.125;
    }
    h
}

fn inject_model_ts_negative_curvature(
    mut hess: Array2<f64>,
    reaction_direction: Option<&[f64]>,
    magnitude: f64,
) -> Result<Array2<f64>> {
    symmetrize(&mut hess);
    let magnitude = magnitude.abs();
    if magnitude <= 0.0 {
        return Ok(hess);
    }
    if let Some(direction) = reaction_direction {
        if direction.len() == hess.nrows() {
            let nrm = norm(direction);
            if nrm > 0.0 {
                let v: Vec<f64> = direction.iter().map(|x| x / nrm).collect();
                let hv = mat_vec(&hess, &v);
                let current = dot(&v, &hv);
                outer_scaled(&v, &v, -(current + magnitude), &mut hess);
                symmetrize(&mut hess);
                return Ok(hess);
            }
        }
    }
    let la = LinAlg::new(Backend::Auto)?;
    let (mut evals, evecs) = la.eigh(&hess)?;
    if !evals.is_empty() {
        let mut idx = 0usize;
        for i in 1..evals.len() {
            if evals[i] < evals[idx] {
                idx = i;
            }
        }
        evals[idx] = -magnitude;
    }
    let mut out = Array2::<f64>::zeros(hess.dim());
    for k in 0..evals.len() {
        for i in 0..evals.len() {
            for j in 0..evals.len() {
                out[(i, j)] += evals[k] * evecs[(i, k)] * evecs[(j, k)];
            }
        }
    }
    symmetrize(&mut out);
    Ok(out)
}

fn repair_ts_minimization_eigenvalues(
    evals: &[f64],
    reaction_index: usize,
    floor: f64,
) -> Vec<f64> {
    let floor = floor.abs();
    let mut repaired = evals.to_vec();
    if floor <= 0.0 {
        return repaired;
    }
    for (idx, value) in repaired.iter_mut().enumerate() {
        if idx != reaction_index && *value < floor {
            *value = value.abs() + floor;
        }
    }
    repaired
}

fn partitioned_rfo_step_from_shifts(
    evals: &[f64],
    grad_eig: &[f64],
    reaction_index: usize,
    lambda_p: f64,
    lambda_n: f64,
) -> Vec<f64> {
    let mut step_eig = vec![0.0; grad_eig.len()];
    let denom = evals[reaction_index] - lambda_p;
    if denom.abs() > 1.0e-14 {
        step_eig[reaction_index] = -grad_eig[reaction_index] / denom;
    }
    for i in 0..evals.len() {
        if i == reaction_index {
            continue;
        }
        let denom = evals[i] - lambda_n;
        if denom.abs() > 1.0e-14 {
            step_eig[i] = -grad_eig[i] / denom;
        }
    }
    step_eig
}

fn trust_constrained_partitioned_rfo_step(
    evals: &[f64],
    grad_eig: &[f64],
    reaction_index: usize,
    lambda_p: f64,
    lambda_n: f64,
    trust_radius: f64,
) -> Vec<f64> {
    let step =
        partitioned_rfo_step_from_shifts(evals, grad_eig, reaction_index, lambda_p, lambda_n);
    if trust_radius <= 0.0 || norm(&step) <= trust_radius {
        return step;
    }
    let shifted = |gamma: f64| -> Vec<f64> {
        partitioned_rfo_step_from_shifts(
            evals,
            grad_eig,
            reaction_index,
            lambda_p + gamma,
            lambda_n - gamma,
        )
    };
    let mut hi = 1.0;
    let mut found = false;
    for _ in 0..80 {
        if norm(&shifted(hi)) <= trust_radius {
            found = true;
            break;
        }
        hi *= 2.0;
    }
    if !found {
        let nrm = norm(&step);
        return step.iter().map(|v| trust_radius * v / nrm).collect();
    }
    let mut lo = 0.0;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if norm(&shifted(mid)) > trust_radius {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    shifted(hi)
}

#[derive(Debug, Clone)]
struct PartitionedRfoResult {
    step: Vec<f64>,
    reaction_mode: Vec<f64>,
    reaction_index: usize,
    evals: Vec<f64>,
    negative_modes: usize,
    diagnostics: TSModeDiagnostics,
}

fn partitioned_rfo_step(
    hess: &Array2<f64>,
    grad: &[f64],
    previous_mode: Option<&[f64]>,
    reaction_direction: Option<&[f64]>,
    follow_lowest: bool,
    trust_radius: f64,
    repair_hessian: bool,
    hessian_eigenvalue_floor: f64,
) -> Result<PartitionedRfoResult> {
    let mut sym = hess.clone();
    symmetrize(&mut sym);
    let la = LinAlg::new(Backend::Auto)?;
    let (evals_a, evecs) = la.eigh(&sym)?;
    let evals = vec_from_array1(&evals_a);
    let n = evals.len();
    let mut mode_source = "lowest".to_string();
    let mut mode_overlap = None;
    let reaction_index = if follow_lowest {
        mode_source = "lowest".to_string();
        argmin(&evals)
    } else if let Some(ref_dir) = reaction_direction {
        if ref_dir.len() != grad.len() || norm(ref_dir) <= 0.0 {
            bail!("reaction direction must be non-zero and match optimizer coordinate length");
        }
        let ref_unit: Vec<f64> = ref_dir.iter().map(|v| v / norm(ref_dir)).collect();
        let mut best = 0usize;
        let mut best_overlap = -1.0;
        for k in 0..n {
            let ov = (0..n)
                .map(|i| evecs[(i, k)] * ref_unit[i])
                .sum::<f64>()
                .abs();
            if ov > best_overlap {
                best_overlap = ov;
                best = k;
            }
        }
        mode_source = "reference".to_string();
        mode_overlap = Some(best_overlap);
        best
    } else if let Some(prev) = previous_mode {
        if prev.len() == grad.len() && norm(prev) > 0.0 {
            let prev_unit: Vec<f64> = prev.iter().map(|v| v / norm(prev)).collect();
            let mut best = 0usize;
            let mut best_overlap = -1.0;
            for k in 0..n {
                let ov = (0..n)
                    .map(|i| evecs[(i, k)] * prev_unit[i])
                    .sum::<f64>()
                    .abs();
                if ov > best_overlap {
                    best_overlap = ov;
                    best = k;
                }
            }
            mode_source = "previous".to_string();
            mode_overlap = Some(best_overlap);
            best
        } else {
            argmin(&evals)
        }
    } else {
        argmin(&evals)
    };

    let mut reaction_mode = (0..n)
        .map(|i| evecs[(i, reaction_index)])
        .collect::<Vec<_>>();
    if let Some(prev) = previous_mode {
        if prev.len() == reaction_mode.len() && dot(&reaction_mode, prev) < 0.0 {
            for v in reaction_mode.iter_mut() {
                *v = -*v;
            }
        }
    }
    if mode_overlap.is_none() {
        if let Some(prev) = previous_mode {
            if prev.len() == reaction_mode.len() && norm(prev) > 0.0 {
                mode_overlap =
                    Some((dot(&reaction_mode, prev) / (norm(&reaction_mode) * norm(prev))).abs());
            }
        }
    }

    let mut grad_eig = vec![0.0; n];
    for k in 0..n {
        for i in 0..n {
            grad_eig[k] += evecs[(i, k)] * grad[i];
        }
    }
    let step_evals = if repair_hessian {
        repair_ts_minimization_eigenvalues(&evals, reaction_index, hessian_eigenvalue_floor)
    } else {
        evals.clone()
    };
    let lam_k = step_evals[reaction_index];
    let grad_k = grad_eig[reaction_index];
    let lambda_p = 0.5 * (lam_k + (lam_k * lam_k + 4.0 * grad_k * grad_k).sqrt());
    let mut min_indices = Vec::new();
    for i in 0..n {
        if i != reaction_index {
            min_indices.push(i);
        }
    }
    let mut lambda_n = 0.0;
    if !min_indices.is_empty() {
        let m = min_indices.len();
        let mut rfo = Array2::<f64>::zeros((m + 1, m + 1));
        for (row, &idx) in min_indices.iter().enumerate() {
            rfo[(row, row)] = step_evals[idx];
            rfo[(row, m)] = grad_eig[idx];
            rfo[(m, row)] = grad_eig[idx];
        }
        let (rfo_eval, _rfo_evec) = la.eigh(&rfo)?;
        lambda_n = rfo_eval[0];
    }
    let step_eig = trust_constrained_partitioned_rfo_step(
        &step_evals,
        &grad_eig,
        reaction_index,
        lambda_p,
        lambda_n,
        trust_radius,
    );
    let mut step = vec![0.0; n];
    for i in 0..n {
        for k in 0..n {
            step[i] += evecs[(i, k)] * step_eig[k];
        }
    }
    let negative_modes = evals.iter().filter(|v| **v < -1.0e-8).count();
    Ok(PartitionedRfoResult {
        step,
        reaction_mode,
        reaction_index,
        evals: evals.clone(),
        negative_modes,
        diagnostics: TSModeDiagnostics {
            index: reaction_index,
            eigenvalue: evals[reaction_index],
            overlap: mode_overlap,
            source: mode_source,
            tracking_space: "optimizer".to_string(),
        },
    })
}

fn argmin(values: &[f64]) -> usize {
    let mut idx = 0usize;
    for i in 1..values.len() {
        if values[i] < values[idx] {
            idx = i;
        }
    }
    idx
}

fn ts_converged(
    energy_change: f64,
    step: &[f64],
    grad: &[f64],
    opts: &TransitionStateOptions,
) -> bool {
    energy_change.abs() <= opts.energy_tol
        && max_abs(grad) <= opts.max_gradient
        && rms(grad) <= opts.rms_gradient
        && max_abs(step) <= opts.max_step
        && rms(step) <= opts.rms_step
}

fn near_ts_convergence_for_recalc(
    step: Option<&[f64]>,
    grad: &[f64],
    opts: &TransitionStateOptions,
) -> bool {
    let Some(step) = step else {
        return false;
    };
    let factor = opts.hessian_recalc_near_convergence_factor;
    if factor <= 0.0 {
        return false;
    }
    max_abs(grad) <= factor * opts.max_gradient
        && rms(grad) <= factor * opts.rms_gradient
        && max_abs(step) <= factor * opts.max_step
        && rms(step) <= factor * opts.rms_step
}

fn adaptive_hessian_reason(
    mode_diag: &TSModeDiagnostics,
    rho: f64,
    negative_modes: usize,
    previous_negative_modes: Option<usize>,
    trust_radius: f64,
    initial_trust_radius: f64,
    opts: &TransitionStateOptions,
) -> Option<String> {
    if let Some(overlap) = mode_diag.overlap {
        if overlap < opts.adaptive_hessian_overlap_min {
            return Some(format!(
                "mode overlap {overlap:.3} < {:.3}",
                opts.adaptive_hessian_overlap_min
            ));
        }
    }
    if !rho.is_finite()
        || rho < opts.adaptive_hessian_rho_min
        || rho > opts.adaptive_hessian_rho_max
    {
        return Some(format!("poor model ratio rho={rho:.3}"));
    }
    if opts.adaptive_hessian_on_negative_mode_change {
        if let Some(prev) = previous_negative_modes {
            if negative_modes != prev {
                return Some(format!(
                    "imaginary-mode count changed {prev}->{negative_modes}"
                ));
            }
        }
    }
    if opts.adaptive_hessian_on_trust_collapse {
        if trust_radius
            <= opts
                .trust_radius_min
                .max(opts.adaptive_hessian_trust_fraction * initial_trust_radius)
        {
            return Some(format!("trust radius collapsed to {}", sci2(trust_radius)));
        }
    }
    None
}

fn reaction_label(mode: &TransitionStateReactionMode) -> Option<String> {
    match mode {
        TransitionStateReactionMode::Lowest => None,
        TransitionStateReactionMode::Direction(_) => {
            Some("Following TS mode: Cartesian direction".to_string())
        }
        TransitionStateReactionMode::Bond([i, j]) => Some(format!("Following TS mode: B {i} {j}")),
        TransitionStateReactionMode::Angle([i, j, k]) => {
            Some(format!("Following TS mode: A {i} {j} {k}"))
        }
        TransitionStateReactionMode::Dihedral([i, j, k, l]) => {
            Some(format!("Following TS mode: D {i} {j} {k} {l}"))
        }
        TransitionStateReactionMode::Transfer([i, j, k]) => {
            Some(format!("Following TS mode: atom-transfer {i} {j} {k}"))
        }
        TransitionStateReactionMode::Coordinates(entries) => {
            let mut parts = Vec::new();
            for entry in entries {
                let label = match entry.kind {
                    InternalReferenceKind::Bond => "B",
                    InternalReferenceKind::Angle => "A",
                    InternalReferenceKind::Dihedral => "D",
                };
                let atoms = entry
                    .atoms
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(" ");
                parts.push(format!("{label} {atoms} * {}", entry.weight));
            }
            Some(format!("Following TS mode: combined {}", parts.join(" + ")))
        }
    }
}

fn cartesian_reaction_direction_from_bond(x: &[f64], pair: [usize; 2]) -> Result<Vec<f64>> {
    let nat = x.len() / 3;
    let [i, j] = pair;
    if i >= nat || j >= nat || i == j {
        bail!("reaction_bond atom indices are invalid");
    }
    let ri = [x[3 * i], x[3 * i + 1], x[3 * i + 2]];
    let rj = [x[3 * j], x[3 * j + 1], x[3 * j + 2]];
    let bond = [rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2]];
    let r = (bond[0] * bond[0] + bond[1] * bond[1] + bond[2] * bond[2]).sqrt();
    if r <= 0.0 {
        bail!("reaction_bond atoms occupy the same position");
    }
    let unit = [bond[0] / r, bond[1] / r, bond[2] / r];
    let mut d = vec![0.0; x.len()];
    for k in 0..3 {
        d[3 * i + k] = -unit[k];
        d[3 * j + k] = unit[k];
    }
    Ok(d)
}

fn cartesian_reaction_direction_from_transfer(x: &[f64], triple: [usize; 3]) -> Result<Vec<f64>> {
    let nat = x.len() / 3;
    let [donor, atom, acceptor] = triple;
    if donor >= nat
        || atom >= nat
        || acceptor >= nat
        || donor == atom
        || atom == acceptor
        || donor == acceptor
    {
        bail!("reaction_transfer atom indices are invalid");
    }
    let p = [x[3 * atom], x[3 * atom + 1], x[3 * atom + 2]];
    let d = [
        x[3 * donor] - p[0],
        x[3 * donor + 1] - p[1],
        x[3 * donor + 2] - p[2],
    ];
    let a = [
        x[3 * acceptor] - p[0],
        x[3 * acceptor + 1] - p[1],
        x[3 * acceptor + 2] - p[2],
    ];
    let dn = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let an = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    if dn <= 0.0 || an <= 0.0 {
        bail!("reaction_transfer atoms occupy coincident positions");
    }
    let mut out = vec![0.0; x.len()];
    for k in 0..3 {
        out[3 * atom + k] = a[k] / an - d[k] / dn;
    }
    Ok(out)
}

fn internal_coordinate_direction(
    ic: &InternalCoords,
    kind: InternalReferenceKind,
    atoms: &[usize],
) -> Result<Vec<f64>> {
    let n = ic.nint();
    let mut out = vec![0.0; n];
    match kind {
        InternalReferenceKind::Bond => {
            if atoms.len() != 2 {
                bail!("bond reference must contain two atoms");
            }
            let target = [atoms[0], atoms[1]];
            for (idx, item) in ic.bonds.iter().enumerate() {
                if *item == target || *item == [target[1], target[0]] {
                    out[idx] = 1.0;
                    return Ok(out);
                }
            }
            bail!(
                "reaction bond {:?} was not found in internal coordinates",
                target
            )
        }
        InternalReferenceKind::Angle => {
            if atoms.len() != 3 {
                bail!("angle reference must contain three atoms");
            }
            let target = [atoms[0], atoms[1], atoms[2]];
            for (i, item) in ic.angles.iter().enumerate() {
                if *item == target || *item == [target[2], target[1], target[0]] {
                    out[ic.nbonds() + i] = 1.0;
                    return Ok(out);
                }
            }
            bail!(
                "reaction angle {:?} was not found in internal coordinates",
                target
            )
        }
        InternalReferenceKind::Dihedral => {
            if atoms.len() != 4 {
                bail!("dihedral reference must contain four atoms");
            }
            let target = [atoms[0], atoms[1], atoms[2], atoms[3]];
            for (i, item) in ic.dihedrals.iter().enumerate() {
                if *item == target || *item == [target[3], target[2], target[1], target[0]] {
                    out[ic.nbonds() + ic.nangles() + i] = 1.0;
                    return Ok(out);
                }
            }
            bail!(
                "reaction dihedral {:?} was not found in internal coordinates",
                target
            )
        }
    }
}

fn internal_reaction_direction_from_mode(
    ic: &InternalCoords,
    mode: &TransitionStateReactionMode,
) -> Result<Option<Vec<f64>>> {
    match mode {
        TransitionStateReactionMode::Lowest => Ok(None),
        TransitionStateReactionMode::Direction(v) => Ok(Some(v.clone())),
        TransitionStateReactionMode::Bond(pair) => Ok(Some(internal_coordinate_direction(
            ic,
            InternalReferenceKind::Bond,
            pair,
        )?)),
        TransitionStateReactionMode::Angle(triple) => Ok(Some(internal_coordinate_direction(
            ic,
            InternalReferenceKind::Angle,
            triple,
        )?)),
        TransitionStateReactionMode::Dihedral(quad) => Ok(Some(internal_coordinate_direction(
            ic,
            InternalReferenceKind::Dihedral,
            quad,
        )?)),
        TransitionStateReactionMode::Transfer(_) => Ok(None),
        TransitionStateReactionMode::Coordinates(entries) => {
            if entries.is_empty() {
                bail!("reaction_coordinates must not be empty");
            }
            let mut combined = vec![0.0; ic.nint()];
            for entry in entries {
                let d = internal_coordinate_direction(ic, entry.kind, &entry.atoms)?;
                for i in 0..combined.len() {
                    combined[i] += entry.weight * d[i];
                }
            }
            let nrm = norm(&combined);
            if nrm <= 0.0 {
                bail!("reaction_coordinates produced a zero reaction direction");
            }
            for v in combined.iter_mut() {
                *v /= nrm;
            }
            Ok(Some(combined))
        }
    }
}

fn build_state(
    x: &[f64],
    mode: TransitionStateCoordinateMode,
    model: Option<&ConnectivityModel>,
    fixed: Option<InternalCoords>,
    extra_bonds: &[[usize; 2]],
    g_cutoff: f64,
) -> Result<TSCoordinateState> {
    match mode {
        TransitionStateCoordinateMode::Cartesian => Ok(TSCoordinateState::Cartesian),
        TransitionStateCoordinateMode::Internal => {
            let ic = if let Some(ic) = fixed {
                ic
            } else {
                let mut ic = analyze_structure(
                    x,
                    model.ok_or_else(|| {
                        anyhow::anyhow!("connectivity model required for internal TS optimization")
                    })?,
                )?;
                if !extra_bonds.is_empty() {
                    let mut bonds = ic.bonds.clone();
                    for &bond in extra_bonds {
                        let reverse = [bond[1], bond[0]];
                        if !bonds.contains(&bond) && !bonds.contains(&reverse) {
                            bonds.push(bond);
                        }
                    }
                    ic = internal_coords_from_bonds(ic.nat, &bonds)?;
                }
                ic
            };
            let q = compute_internals(x, &ic);
            let bpg = bpg_matrix(x, &ic)?;
            Ok(TSCoordinateState::Internal {
                ic,
                q,
                b: bpg.b,
                inv_g: bpg.inv_g,
            })
        }
        TransitionStateCoordinateMode::PulayRedundant => {
            let coordinates = if fixed.is_none() && !extra_bonds.is_empty() {
                let mut ic = analyze_structure(
                    x,
                    model.ok_or_else(|| {
                        anyhow::anyhow!(
                            "connectivity model required for Pulay redundant TS optimization"
                        )
                    })?,
                )?;
                let mut bonds = ic.bonds.clone();
                for &bond in extra_bonds {
                    let reverse = [bond[1], bond[0]];
                    if !bonds.contains(&bond) && !bonds.contains(&reverse) {
                        bonds.push(bond);
                    }
                }
                ic = internal_coords_from_bonds(ic.nat, &bonds)?;
                Some(ic)
            } else {
                fixed
            };
            let system = build_redundant_internals(x, model, coordinates, None, g_cutoff)?;
            Ok(TSCoordinateState::Pulay(system))
        }
    }
}

fn optimizer_gradient(state: &TSCoordinateState, grad_x: &[f64]) -> Result<Vec<f64>> {
    match state {
        TSCoordinateState::Cartesian => Ok(grad_x.to_vec()),
        TSCoordinateState::Internal { b, inv_g, .. } => {
            Ok(internal_gradient_from_bpg(b, inv_g, grad_x))
        }
        TSCoordinateState::Pulay(system) => cartesian_gradient_to_redundant(system, grad_x),
    }
}

fn optimizer_hessian_from_cartesian(
    state: &TSCoordinateState,
    hess_x: &Array2<f64>,
    redundant_penalty: f64,
) -> Result<Array2<f64>> {
    match state {
        TSCoordinateState::Cartesian => Ok(hess_x.clone()),
        TSCoordinateState::Internal { b, inv_g, .. } => {
            cartesian_hessian_to_internal(hess_x, b, inv_g)
        }
        TSCoordinateState::Pulay(system) => {
            let hq = cartesian_hessian_to_internal(hess_x, &system.b, &system.ginv)?;
            initial_redundant_hessian(system, &hq, redundant_penalty)
        }
    }
}

fn model_hessian_for_state(
    state: &TSCoordinateState,
    redundant_penalty: f64,
) -> Result<Array2<f64>> {
    match state {
        TSCoordinateState::Cartesian => {
            bail!("model initial Hessian is only available for internal-coordinate TS optimization")
        }
        TSCoordinateState::Internal { ic, .. } => Ok(initial_model_hessian(ic)),
        TSCoordinateState::Pulay(system) => {
            let hq = initial_model_hessian(&system.coordinates);
            initial_redundant_hessian(system, &hq, redundant_penalty)
        }
    }
}

fn accept_step<O: Objective>(
    obj: &mut O,
    x: &[f64],
    energy: f64,
    grad_opt: &[f64],
    hess_opt: &Array2<f64>,
    step_opt: &[f64],
    state: &TSCoordinateState,
    opts: &TransitionStateOptions,
    trust_radius: f64,
) -> Result<(Vec<f64>, f64, Vec<f64>, Vec<f64>, Vec<f64>, f64, Vec<f64>)> {
    let mut accepted_step_opt = step_opt.to_vec();
    let mut x_trial = match state {
        TSCoordinateState::Cartesian => x
            .iter()
            .zip(step_opt.iter())
            .map(|(a, b)| a + b)
            .collect::<Vec<_>>(),
        TSCoordinateState::Internal { ic, q, b, inv_g } => {
            let bpg = crate::optimizer::internal_coords::Bpg {
                b: b.clone(),
                inv_g: inv_g.clone(),
            };
            let mut xt = x.to_vec();
            best_fit_dq_to_cart(&mut xt, q, step_opt, ic, &bpg, opts.best_fit_iters)?;
            xt
        }
        TSCoordinateState::Pulay(system) => {
            update_cartesian_from_redundant_step(x, system, step_opt, opts.best_fit_iters)?
        }
    };
    if !matches!(state, TSCoordinateState::Cartesian) && trust_radius > 0.0 {
        let mut scale = 1.0;
        for _ in 0..8 {
            let cart_step: Vec<f64> = x_trial.iter().zip(x.iter()).map(|(a, b)| a - b).collect();
            let cart_norm = norm(&cart_step);
            if cart_norm <= trust_radius {
                break;
            }
            scale *= trust_radius / cart_norm;
            accepted_step_opt = step_opt.iter().map(|v| scale * v).collect();
            x_trial = match state {
                TSCoordinateState::Internal { ic, q, b, inv_g } => {
                    let bpg = crate::optimizer::internal_coords::Bpg {
                        b: b.clone(),
                        inv_g: inv_g.clone(),
                    };
                    let mut xt = x.to_vec();
                    best_fit_dq_to_cart(
                        &mut xt,
                        q,
                        &accepted_step_opt,
                        ic,
                        &bpg,
                        opts.best_fit_iters,
                    )?;
                    xt
                }
                TSCoordinateState::Pulay(system) => update_cartesian_from_redundant_step(
                    x,
                    system,
                    &accepted_step_opt,
                    opts.best_fit_iters,
                )?,
                TSCoordinateState::Cartesian => unreachable!(),
            };
        }
    }
    let mut trial_grad_x = vec![0.0; x.len()];
    let trial_energy = obj.energy_gradient(&x_trial, &mut trial_grad_x)?;
    maybe_project_ts_derivatives(
        opts.masses_amu.as_deref(),
        &x_trial,
        Some(&mut trial_grad_x),
        None,
        opts.project_eckart,
    )?;
    let cart_step: Vec<f64> = x_trial.iter().zip(x.iter()).map(|(a, b)| a - b).collect();
    let h_step = mat_vec(hess_opt, &accepted_step_opt);
    let pred = dot(grad_opt, &accepted_step_opt) + 0.5 * dot(&accepted_step_opt, &h_step);
    let actual = trial_energy - energy;
    let rho = if pred.abs() > 1.0e-14 {
        actual / pred
    } else {
        f64::NAN
    };
    Ok((
        x_trial,
        trial_energy,
        trial_grad_x,
        cart_step,
        accepted_step_opt,
        rho,
        h_step,
    ))
}

fn count_negative_modes(hess: &Array2<f64>) -> Result<usize> {
    let mut h = hess.clone();
    symmetrize(&mut h);
    let la = LinAlg::new(Backend::Auto)?;
    let (evals, _evecs) = la.eigh(&h)?;
    Ok(evals.iter().filter(|v| **v < -1.0e-8).count())
}

fn write_xyz_frame(
    handle: &mut std::fs::File,
    x: &[f64],
    step: usize,
    energy: f64,
    grad: &[f64],
) -> Result<()> {
    use std::io::Write;
    let nat = x.len() / 3;
    writeln!(handle, "{nat}")?;
    writeln!(
        handle,
        "step= {step:10} Ene = {energy:20.8} GradMax = {:20.8}",
        max_abs(grad)
    )?;
    const BOHR_TO_ANG: f64 = 0.529177210903;
    for i in 0..nat {
        writeln!(
            handle,
            "X {:15.8} {:15.8} {:15.8}",
            x[3 * i] * BOHR_TO_ANG,
            x[3 * i + 1] * BOHR_TO_ANG,
            x[3 * i + 2] * BOHR_TO_ANG,
        )?;
    }
    Ok(())
}

pub fn optimize_transition_state<O: Objective>(
    x0: Vec<f64>,
    connectivity: Option<ConnectivityModel>,
    obj: &mut O,
    mut opts: TransitionStateOptions,
) -> Result<TransitionStateResult> {
    if x0.len() % 3 != 0 {
        bail!("TS optimizer requires 3N Cartesian coordinates");
    }
    if opts.hessian_step <= 0.0 {
        bail!("hessian_step must be positive");
    }
    if opts.coordinate_mode == TransitionStateCoordinateMode::Cartesian
        && matches!(
            opts.ts_initial_hessian,
            TransitionStateInitialHessian::Model
        )
    {
        bail!("Cartesian TS optimization requires exact finite-difference initial Hessian");
    }
    if opts.coordinate_mode == TransitionStateCoordinateMode::Cartesian {
        if matches!(
            opts.reaction_mode,
            TransitionStateReactionMode::Angle(_)
                | TransitionStateReactionMode::Dihedral(_)
                | TransitionStateReactionMode::Coordinates(_)
        ) {
            bail!("angle, dihedral, and combined reaction-coordinate references require internal TS coordinates");
        }
    }
    if opts.project_eckart && opts.masses_amu.is_none() {
        bail!("TS Eckart projection requested but TransitionStateOptions::masses_amu is not set");
    }
    let mut x = x0;
    let model = connectivity.as_ref();
    let mut state = build_state(
        &x,
        opts.coordinate_mode,
        model,
        None,
        &opts.extra_bonds,
        opts.g_cutoff,
    )?;
    let fixed_ic = state.internal_coords().cloned();
    if matches!(
        opts.coordinate_mode,
        TransitionStateCoordinateMode::Internal | TransitionStateCoordinateMode::PulayRedundant
    ) && fixed_ic.is_none()
    {
        bail!("internal TS optimization found no internal coordinates");
    }
    let mut trust_radius = opts.trust_radius.unwrap_or_else(|| {
        if opts.coordinate_mode == TransitionStateCoordinateMode::Cartesian {
            0.1
        } else {
            opts.max_step_internal
        }
    });
    let trust_radius_max = opts.trust_radius_max.unwrap_or_else(|| {
        if opts.coordinate_mode == TransitionStateCoordinateMode::Cartesian {
            0.3
        } else {
            0.4
        }
    });
    trust_radius = trust_radius.clamp(opts.trust_radius_min, trust_radius_max);
    opts.trust_radius = Some(trust_radius);

    let mut grad_x = vec![0.0; x.len()];
    let mut energy = obj.energy_gradient(&x, &mut grad_x)?;
    maybe_project_ts_derivatives(
        opts.masses_amu.as_deref(),
        &x,
        Some(&mut grad_x),
        None,
        opts.project_eckart,
    )?;
    let mut grad_opt = optimizer_gradient(&state, &grad_x)?;
    let reaction_direction = match (&opts.reaction_mode, &state) {
        (TransitionStateReactionMode::Lowest, _) => None,
        (TransitionStateReactionMode::Direction(v), TSCoordinateState::Cartesian) => {
            Some(v.clone())
        }
        (TransitionStateReactionMode::Bond(pair), TSCoordinateState::Cartesian) => {
            Some(cartesian_reaction_direction_from_bond(&x, *pair)?)
        }
        (TransitionStateReactionMode::Transfer(triple), TSCoordinateState::Cartesian) => {
            Some(cartesian_reaction_direction_from_transfer(&x, *triple)?)
        }
        (_, TSCoordinateState::Cartesian) => None,
        (_, _) => internal_reaction_direction_from_mode(
            state.internal_coords().unwrap(),
            &opts.reaction_mode,
        )?,
    };

    let mut hess_opt = if opts.ts_initial_hessian == TransitionStateInitialHessian::Model {
        model_hessian_for_state(&state, opts.redundant_penalty)?
    } else {
        let mut hess_x = finite_difference_hessian(obj, &x, opts.hessian_step)?;
        maybe_project_ts_derivatives(
            opts.masses_amu.as_deref(),
            &x,
            None,
            Some(&mut hess_x),
            opts.project_eckart,
        )?;
        optimizer_hessian_from_cartesian(&state, &hess_x, opts.redundant_penalty)?
    };
    if opts.ts_initial_hessian == TransitionStateInitialHessian::Model {
        hess_opt = inject_model_ts_negative_curvature(
            hess_opt,
            reaction_direction.as_deref(),
            opts.ts_model_negative_curvature,
        )?;
    }

    let mut reaction_mode = None::<Vec<f64>>;
    let mut reaction_mode_is_lowest = false;
    let mut negative_modes = None::<usize>;
    let mut previous_negative_modes = None::<usize>;
    let mut previous_cart_step = None::<Vec<f64>>;
    let mut adaptive_recalc_next_reason = None::<String>;
    let mut last_hessian_recalc_step = 0usize;
    let initial_trust_radius = trust_radius;
    let mut trajectory = vec![x.clone()];
    let mut trajectory_energies = vec![energy];

    if opts.verbosity > 0 {
        match opts.coordinate_mode {
            TransitionStateCoordinateMode::Cartesian => {
                println!("Transition-state optimization with Cartesian partitioned RFO");
                println!("TS mode tracking coordinates: cartesian");
            }
            TransitionStateCoordinateMode::Internal
            | TransitionStateCoordinateMode::PulayRedundant => {
                println!("Transition-state optimization with internal-coordinate partitioned RFO");
                println!("Update method: Bofill");
                println!("Choice of coordinates: {}", state.coordinate_label());
                if let Some(ic) = state.internal_coords() {
                    println!(
                        "Internal coordinates: {} bonds, {} angles, {} dihedrals",
                        ic.nbonds(),
                        ic.nangles(),
                        ic.ndiheds()
                    );
                }
                println!(
                    "TS mode tracking coordinates: {:?}",
                    opts.mode_tracking_coordinates
                );
                println!("Initial TS Hessian: {:?}", opts.ts_initial_hessian);
                if opts.ts_initial_hessian == TransitionStateInitialHessian::Model {
                    println!(
                        "Initial internal Hessian model: {}",
                        opts.internal_hessian_model
                    );
                    println!(
                        "Injected model TS negative curvature: {:.6}",
                        opts.ts_model_negative_curvature.abs()
                    );
                }
            }
        }
        if let Some(label) = reaction_label(&opts.reaction_mode) {
            println!("{label}");
        }
        print_ts_header(&opts);
    }

    let mut traj_handle = if let Some(path) = &opts.trajectory_file {
        Some(std::fs::File::create(path)?)
    } else {
        None
    };
    if let Some(handle) = traj_handle.as_mut() {
        write_xyz_frame(handle, &x, 0, energy, &grad_x)?;
    }

    for cycle in 1..=opts.max_cycles {
        let mut fixed_recalc = opts.hessian_recalc_interval > 0
            && cycle > 1
            && (cycle - 1) % opts.hessian_recalc_interval == 0;
        if fixed_recalc && opts.skip_hessian_recalc_near_convergence {
            if near_ts_convergence_for_recalc(previous_cart_step.as_deref(), &grad_x, &opts) {
                if opts.verbosity > 0 {
                    println!("      TS scheduled Hessian recalculation skipped near convergence at step {}", cycle - 1);
                }
                fixed_recalc = false;
            }
        }
        if fixed_recalc || adaptive_recalc_next_reason.is_some() {
            if opts.verbosity > 0 {
                if fixed_recalc {
                    println!(
                        "      TS scheduled Hessian recalculation at step {}",
                        cycle - 1
                    );
                }
                if let Some(reason) = &adaptive_recalc_next_reason {
                    println!("      TS adaptive Hessian recalculation: {reason}");
                }
            }
            let mut hess_x = finite_difference_hessian(obj, &x, opts.hessian_step)?;
            maybe_project_ts_derivatives(
                opts.masses_amu.as_deref(),
                &x,
                None,
                Some(&mut hess_x),
                opts.project_eckart,
            )?;
            hess_opt = optimizer_hessian_from_cartesian(&state, &hess_x, opts.redundant_penalty)?;
            last_hessian_recalc_step = cycle - 1;
            adaptive_recalc_next_reason = None;
        }

        let rfo = partitioned_rfo_step(
            &hess_opt,
            &grad_opt,
            if reaction_mode_is_lowest {
                None
            } else {
                reaction_mode.as_deref()
            },
            if reaction_mode.is_none() && !reaction_mode_is_lowest {
                reaction_direction.as_deref()
            } else {
                None
            },
            reaction_mode_is_lowest,
            trust_radius,
            opts.repair_ts_hessian,
            opts.ts_hessian_eigenvalue_floor,
        )?;
        let (trial_x, trial_energy, trial_grad_x, cart_step, accepted_step_opt, rho, _h_step) =
            accept_step(
                obj,
                &x,
                energy,
                &grad_opt,
                &hess_opt,
                &rfo.step,
                &state,
                &opts,
                trust_radius,
            )?;

        let energy_change = trial_energy - energy;
        let trial_state = build_state(
            &trial_x,
            opts.coordinate_mode,
            model,
            fixed_ic.clone(),
            &opts.extra_bonds,
            opts.g_cutoff,
        )?;
        let trial_grad_opt = optimizer_gradient(&trial_state, &trial_grad_x)?;
        reaction_mode = Some(rfo.reaction_mode.clone());
        if rfo.reaction_index == argmin(&rfo.evals) {
            reaction_mode_is_lowest = true;
        }
        negative_modes = Some(rfo.negative_modes);

        if opts.verbosity > 0 {
            print_ts_step(
                cycle,
                trial_energy,
                energy_change,
                &cart_step,
                &trial_grad_x,
                trust_radius,
                rfo.negative_modes,
                Some(&rfo.diagnostics),
            );
        }

        if opts.adaptive_hessian_recalc {
            if cycle.saturating_sub(last_hessian_recalc_step)
                >= opts.adaptive_hessian_recalc_cooldown
            {
                adaptive_recalc_next_reason = adaptive_hessian_reason(
                    &rfo.diagnostics,
                    rho,
                    rfo.negative_modes,
                    previous_negative_modes,
                    trust_radius,
                    initial_trust_radius,
                    &opts,
                );
            } else {
                adaptive_recalc_next_reason = None;
            }
        }
        previous_negative_modes = Some(rfo.negative_modes);

        if let Some(handle) = traj_handle.as_mut() {
            write_xyz_frame(handle, &trial_x, cycle, trial_energy, &trial_grad_x)?;
        }
        trajectory.push(trial_x.clone());
        trajectory_energies.push(trial_energy);

        if ts_converged(energy_change, &cart_step, &trial_grad_x, &opts) {
            let mut result_hessian = None;
            let mut result_negative_modes = negative_modes;
            if opts.final_hessian {
                let mut hess_x = finite_difference_hessian(obj, &trial_x, opts.hessian_step)?;
                maybe_project_ts_derivatives(
                    opts.masses_amu.as_deref(),
                    &trial_x,
                    None,
                    Some(&mut hess_x),
                    opts.project_eckart,
                )?;
                let final_state = build_state(
                    &trial_x,
                    opts.coordinate_mode,
                    model,
                    fixed_ic.clone(),
                    &opts.extra_bonds,
                    opts.g_cutoff,
                )?;
                let hess_eval = optimizer_hessian_from_cartesian(
                    &final_state,
                    &hess_x,
                    opts.redundant_penalty,
                )?;
                result_negative_modes = Some(count_negative_modes(&hess_eval)?);
                result_hessian = Some(hess_x);
            }
            if opts.verbosity > 0 {
                print_optimization_status(true, "Transition-state convergence reached.");
            }
            return Ok(TransitionStateResult {
                x: trial_x,
                energy: trial_energy,
                gradient: trial_grad_x,
                converged: true,
                cycles: cycle,
                message: "Transition-state convergence reached.".to_string(),
                coordinates: opts.coordinate_mode,
                negative_modes: result_negative_modes,
                reaction_mode,
                hessian: result_hessian,
                trajectory,
                trajectory_energies,
            });
        }

        let y: Vec<f64> = trial_grad_opt
            .iter()
            .zip(grad_opt.iter())
            .map(|(a, b)| a - b)
            .collect();
        let (updated_hess, _updated) = bofill_update_hessian(&hess_opt, &accepted_step_opt, &y);
        hess_opt = updated_hess;

        x = trial_x;
        state = trial_state;
        energy = trial_energy;
        grad_x = trial_grad_x;
        grad_opt = trial_grad_opt;
        previous_cart_step = Some(cart_step);
    }

    if opts.verbosity > 0 {
        print_optimization_status(
            false,
            "Maximum TS optimization steps reached without convergence.",
        );
    }
    Ok(TransitionStateResult {
        x,
        energy,
        gradient: grad_x,
        converged: false,
        cycles: opts.max_cycles,
        message: "Maximum TS optimization steps reached without convergence.".to_string(),
        coordinates: opts.coordinate_mode,
        negative_modes,
        reaction_mode,
        hessian: None,
        trajectory,
        trajectory_energies,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    struct SaddleObjective;

    impl Objective for SaddleObjective {
        fn energy(&mut self, x: &[f64]) -> Result<f64> {
            Ok(-0.5 * x[0] * x[0] + 0.5 * x[1] * x[1] + 0.5 * x[2] * x[2])
        }

        fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
            grad[0] = -x[0];
            grad[1] = x[1];
            grad[2] = x[2];
            Ok(())
        }
    }

    #[test]
    fn cartesian_ts_finds_quadratic_saddle() {
        let mut obj = SaddleObjective;
        let opts = TransitionStateOptions {
            coordinate_mode: TransitionStateCoordinateMode::Cartesian,
            max_cycles: 20,
            hessian_step: 1.0e-4,
            final_hessian: false,
            project_eckart: false,
            trajectory_file: None,
            verbosity: 0,
            reaction_mode: TransitionStateReactionMode::Lowest,
            ..Default::default()
        };
        let res = optimize_transition_state(vec![0.2, 0.1, -0.1], None, &mut obj, opts).unwrap();
        assert!(res.converged);
        assert!(res.energy.abs() < 1.0e-8);
        assert!(max_abs(&res.gradient) < 1.0e-4);
    }
}
