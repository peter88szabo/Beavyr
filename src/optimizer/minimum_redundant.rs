//! Smite-style minimum search in primitive internal coordinates.
//!
//! This module is separate from the legacy `geom_opt` optimizer. It supports
//! both the existing simple redundant internals and Pulay-projected redundant
//! internals adapted from Smite.

use anyhow::{bail, Result};
use ndarray::{Array1, Array2};

use crate::optimizer::linalg_solve::solve_linear_system_array2;
use crate::optimizer::internal_coords::{
    analyze_structure, best_fit_dq_to_cart, bpg_matrix, compute_internals, diff_internals, Bpg,
    ConnectivityModel, InternalCoords,
};
use crate::optimizer::redundant_internals::{
    build_redundant_internals, cartesian_gradient_to_redundant, initial_redundant_hessian,
    update_cartesian_from_redundant_step, RedundantInternalSystem,
};
use crate::optimizer::traits::Objective;

const HARTREE_TO_KJMOL: f64 = 2625.499638;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InternalMinimumMethod {
    Bfgs,
    SteepestDescent,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InternalCoordinateMode {
    SimpleRedundant,
    PulayRedundant,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InitialInternalHessianModel {
    Simple,
}

#[derive(Debug, Clone)]
pub struct InternalMinimumOptions {
    pub method: InternalMinimumMethod,
    pub coordinate_mode: InternalCoordinateMode,
    pub hessian_model: InitialInternalHessianModel,
    pub max_cycles: usize,
    pub energy_tol: f64,
    pub max_step: f64,
    pub rms_step: f64,
    pub max_gradient: f64,
    pub rms_gradient: f64,
    pub max_step_internal: f64,
    pub best_fit_iters: usize,
    pub best_fit_rms_tol: f64,
    pub min_line_search_scale: f64,
    pub g_cutoff: f64,
    pub redundant_penalty: f64,
    pub verbosity: u8,
}

impl Default for InternalMinimumOptions {
    fn default() -> Self {
        Self {
            method: InternalMinimumMethod::Bfgs,
            coordinate_mode: InternalCoordinateMode::PulayRedundant,
            hessian_model: InitialInternalHessianModel::Simple,
            max_cycles: 100,
            energy_tol: 5.0e-6,
            max_step: 4.0e-3,
            rms_step: 2.0e-3,
            max_gradient: 3.0e-4,
            rms_gradient: 1.0e-4,
            max_step_internal: 0.2,
            best_fit_iters: 5,
            best_fit_rms_tol: 1.0e-7,
            min_line_search_scale: 1.0e-4,
            g_cutoff: 1.0e-7,
            redundant_penalty: 1000.0,
            verbosity: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct InternalMinimumResult {
    pub x: Vec<f64>,
    pub energy: f64,
    pub gradient: Vec<f64>,
    pub converged: bool,
    pub cycles: usize,
    pub message: String,
    pub coordinates: InternalCoordinateMode,
    pub trajectory: Vec<Vec<f64>>,
    pub trajectory_energies: Vec<f64>,
}

#[derive(Debug, Clone)]
struct SimpleInternalState {
    ic: InternalCoords,
    q: Vec<f64>,
    bpg: Bpg,
}

#[derive(Debug, Clone)]
enum InternalState {
    Simple(SimpleInternalState),
    Pulay(RedundantInternalSystem),
}

impl InternalState {
    fn coordinates(&self) -> &InternalCoords {
        match self {
            Self::Simple(s) => &s.ic,
            Self::Pulay(s) => &s.coordinates,
        }
    }

    fn q(&self) -> &[f64] {
        match self {
            Self::Simple(s) => &s.q,
            Self::Pulay(s) => &s.q,
        }
    }

    fn nint(&self) -> usize {
        self.coordinates().nint()
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

fn print_header(opts: &InternalMinimumOptions) {
    println!(" Energy Change:   {} Eh", sci2(opts.energy_tol));
    println!(" Max. Gradient:   {} Eh/bohr", sci2(opts.max_gradient));
    println!(" RMS Gradient:    {} Eh/bohr", sci2(opts.rms_gradient));
    println!(" Max Step:        {} bohr", sci2(opts.max_step));
    println!(" RMS Step:        {} bohr", sci2(opts.rms_step));
    println!();

    let header_line = format!(
        "{:<5} {:>18} {:>18} {:>11} {:>11} {:>11} {:>11}",
        "step", "E[Eh]", "dE[kJ/mol]", "max_step", "rms_step", "max_grad", "rms_grad"
    );
    print_threshold_guide(
        opts.energy_tol,
        opts.max_step,
        opts.rms_step,
        opts.max_gradient,
        opts.rms_gradient,
        &header_line,
    );
    println!("{header_line}");
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

fn print_threshold_guide(
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

fn max_abs(a: &[f64]) -> f64 {
    a.iter().fold(0.0, |m, v| m.max(v.abs()))
}

fn internal_gradient_simple(bpg: &Bpg, grad_x: &[f64]) -> Vec<f64> {
    let nint = bpg.b.len();
    let mut bg = vec![0.0; nint];
    for i in 0..nint {
        for k in 0..grad_x.len() {
            bg[i] += bpg.b[i][k] * grad_x[k];
        }
    }
    let mut gq = vec![0.0; nint];
    for i in 0..nint {
        for j in 0..nint {
            gq[i] += bpg.inv_g[i][j] * bg[j];
        }
    }
    gq
}

fn internal_gradient(state: &InternalState, grad_x: &[f64]) -> Result<Vec<f64>> {
    match state {
        InternalState::Simple(s) => Ok(internal_gradient_simple(&s.bpg, grad_x)),
        InternalState::Pulay(s) => cartesian_gradient_to_redundant(s, grad_x),
    }
}

fn build_state(
    x: &[f64],
    mode: InternalCoordinateMode,
    model: Option<&ConnectivityModel>,
    coordinates: Option<InternalCoords>,
    g_cutoff: f64,
) -> Result<InternalState> {
    match mode {
        InternalCoordinateMode::SimpleRedundant => {
            let ic = if let Some(ic) = coordinates {
                ic
            } else {
                analyze_structure(
                    x,
                    model.ok_or_else(|| anyhow::anyhow!("connectivity model required"))?,
                )?
            };
            let q = compute_internals(x, &ic);
            let bpg = bpg_matrix(x, &ic)?;
            Ok(InternalState::Simple(SimpleInternalState { ic, q, bpg }))
        }
        InternalCoordinateMode::PulayRedundant => Ok(InternalState::Pulay(
            build_redundant_internals(x, model, coordinates, None, g_cutoff)?,
        )),
    }
}

fn initial_hessian(ic: &InternalCoords, _model: InitialInternalHessianModel) -> Array2<f64> {
    let nb = ic.nbonds();
    let na = ic.nangles();
    let nd = ic.ndiheds();
    let nint = nb + na + nd;
    let mut h = Array2::<f64>::zeros((nint, nint));
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

fn solve_step(hess: &Array2<f64>, grad: &[f64]) -> Result<Vec<f64>> {
    let rhs = Array1::from(grad.iter().map(|g| -g).collect::<Vec<_>>());
    Ok(solve_linear_system_array2(hess, &rhs)?.to_vec())
}

fn limit_internal_step(dq: &mut [f64], max_step_internal: f64) {
    let stp = norm(dq);
    if stp > max_step_internal && stp > 0.0 {
        let scale = max_step_internal / stp;
        for v in dq.iter_mut() {
            *v *= scale;
        }
    }
}

fn converged(
    energy_change: f64,
    cart_step: &[f64],
    grad_x: &[f64],
    opts: &InternalMinimumOptions,
) -> bool {
    max_abs(grad_x) <= opts.max_gradient
        && rms(grad_x) <= opts.rms_gradient
        && max_abs(cart_step) <= opts.max_step
        && rms(cart_step) <= opts.rms_step
        && energy_change.abs() <= opts.energy_tol
}

fn apply_internal_step(
    state: &InternalState,
    x: &[f64],
    dq: &[f64],
    opts: &InternalMinimumOptions,
) -> Result<Vec<f64>> {
    match state {
        InternalState::Simple(s) => {
            let mut trial = x.to_vec();
            best_fit_dq_to_cart(&mut trial, &s.q, dq, &s.ic, &s.bpg, opts.best_fit_iters)?;
            Ok(trial)
        }
        InternalState::Pulay(s) => {
            update_cartesian_from_redundant_step(x, s, dq, opts.best_fit_iters)
        }
    }
}

fn line_search_internal<O: Objective>(
    obj: &mut O,
    x: &[f64],
    state: &InternalState,
    energy: f64,
    dq: &[f64],
    opts: &InternalMinimumOptions,
) -> Result<(Vec<f64>, f64, Vec<f64>, Vec<f64>, Vec<f64>, f64)> {
    let mut scale = 1.0;
    while scale >= opts.min_line_search_scale {
        let scaled_dq: Vec<f64> = dq.iter().map(|v| scale * v).collect();
        let trial_x = apply_internal_step(state, x, &scaled_dq, opts)?;
        let trial_energy = obj.energy(&trial_x)?;
        if trial_energy < energy {
            let mut trial_grad = vec![0.0; x.len()];
            let trial_energy = obj.energy_gradient(&trial_x, &mut trial_grad)?;
            let cart_step: Vec<f64> = trial_x.iter().zip(x.iter()).map(|(a, b)| a - b).collect();
            return Ok((
                trial_x,
                trial_energy,
                trial_grad,
                scaled_dq,
                cart_step,
                scale,
            ));
        }
        scale *= 0.5;
    }

    let scaled_dq: Vec<f64> = dq.iter().map(|v| scale * v).collect();
    let trial_x = apply_internal_step(state, x, &scaled_dq, opts)?;
    let mut trial_grad = vec![0.0; x.len()];
    let trial_energy = obj.energy_gradient(&trial_x, &mut trial_grad)?;
    let cart_step: Vec<f64> = trial_x.iter().zip(x.iter()).map(|(a, b)| a - b).collect();
    Ok((
        trial_x,
        trial_energy,
        trial_grad,
        scaled_dq,
        cart_step,
        scale,
    ))
}

fn bfgs_update_hessian(hess: &Array2<f64>, step: &[f64], y: &[f64]) -> Array2<f64> {
    let ys = dot(y, step);
    let h_step = hess.dot(&Array1::from(step.to_vec())).to_vec();
    let shs = dot(step, &h_step);
    if ys <= 1.0e-14 || shs <= 1.0e-14 {
        return hess.clone();
    }
    let mut out = hess.clone();
    let n = step.len();
    for i in 0..n {
        for j in 0..n {
            out[(i, j)] += (y[i] * y[j]) / ys - (h_step[i] * h_step[j]) / shs;
        }
    }
    for i in 0..n {
        for j in 0..i {
            let v = 0.5 * (out[(i, j)] + out[(j, i)]);
            out[(i, j)] = v;
            out[(j, i)] = v;
        }
    }
    out
}

pub fn minimize_internal<O: Objective>(
    mut x: Vec<f64>,
    conn: ConnectivityModel,
    obj: &mut O,
    opts: InternalMinimumOptions,
) -> Result<InternalMinimumResult> {
    let nat = conn.rcov_disp.len();
    if x.len() != 3 * nat {
        bail!("internal optimizer expected 3N coordinates");
    }
    if nat < 2 {
        let mut grad = vec![0.0; x.len()];
        let energy = obj.energy_gradient(&x, &mut grad)?;
        let ok = max_abs(&grad) <= opts.max_gradient;
        return Ok(InternalMinimumResult {
            x,
            energy,
            gradient: grad,
            converged: ok,
            cycles: 0,
            message: "Single-atom system has no internal coordinates.".to_string(),
            coordinates: opts.coordinate_mode,
            trajectory: Vec::new(),
            trajectory_energies: Vec::new(),
        });
    }

    let mut state = build_state(&x, opts.coordinate_mode, Some(&conn), None, opts.g_cutoff)?;
    if state.nint() == 0 {
        bail!("no internal coordinates found");
    }
    let mut hess_q = initial_hessian(state.coordinates(), opts.hessian_model);
    if let InternalState::Pulay(system) = &state {
        hess_q = initial_redundant_hessian(system, &hess_q, opts.redundant_penalty)?;
    }

    let mut grad_x = vec![0.0; x.len()];
    let mut energy = obj.energy_gradient(&x, &mut grad_x)?;
    let mut grad_q = internal_gradient(&state, &grad_x)?;
    let mut trajectory = vec![x.clone()];
    let mut trajectory_energies = vec![energy];

    if opts.verbosity > 0 {
        match opts.coordinate_mode {
            InternalCoordinateMode::SimpleRedundant => {
                println!("Choice of coordinates: simple redundant internal coordinates")
            }
            InternalCoordinateMode::PulayRedundant => {
                println!("Choice of coordinates: Pulay redundant internal coordinates")
            }
        }
        println!(
            "Internal coordinates: {} bonds, {} angles, {} dihedrals",
            state.coordinates().nbonds(),
            state.coordinates().nangles(),
            state.coordinates().ndiheds()
        );
        print_header(&opts);
    }

    for cycle in 1..=opts.max_cycles {
        let mut dq = match opts.method {
            InternalMinimumMethod::SteepestDescent => grad_q.iter().map(|g| -0.7 * g).collect(),
            InternalMinimumMethod::Bfgs => solve_step(&hess_q, &grad_q)?,
        };
        limit_internal_step(&mut dq, opts.max_step_internal);

        let old_x = x.clone();
        let old_energy = energy;
        let old_grad_q = grad_q.clone();
        let old_state = state.clone();
        let old_hess = hess_q.clone();

        let (trial_x, trial_energy, trial_grad_x, accepted_dq, cart_step, scale) =
            line_search_internal(obj, &x, &state, energy, &dq, &opts)?;
        let energy_change = trial_energy - energy;

        x = trial_x;
        energy = trial_energy;
        grad_x = trial_grad_x;
        state = build_state(
            &x,
            opts.coordinate_mode,
            None,
            Some(old_state.coordinates().clone()),
            opts.g_cutoff,
        )?;
        grad_q = internal_gradient(&state, &grad_x)?;

        if opts.verbosity > 0 {
            println!(
                "{:5} {:18.8} {:18.6} {:>11} {:>11} {:>11} {:>11}",
                cycle,
                energy,
                energy_change * HARTREE_TO_KJMOL,
                sci2(max_abs(&cart_step)),
                sci2(rms(&cart_step)),
                sci2(max_abs(&grad_x)),
                sci2(rms(&grad_x))
            );
        }

        trajectory.push(x.clone());
        trajectory_energies.push(energy);

        if converged(energy_change, &cart_step, &grad_x, &opts) {
            if opts.verbosity > 0 {
                print_optimization_status(true, "Convergence reached.");
            }
            return Ok(InternalMinimumResult {
                x,
                energy,
                gradient: grad_x,
                converged: true,
                cycles: cycle,
                message: "Convergence reached.".to_string(),
                coordinates: opts.coordinate_mode,
                trajectory,
                trajectory_energies,
            });
        }

        if opts.method == InternalMinimumMethod::Bfgs {
            let y = diff_internals(&grad_q, &old_grad_q, state.coordinates().ndiheds());
            let candidate = bfgs_update_hessian(&hess_q, &accepted_dq, &y);
            if candidate.iter().all(|v| v.is_finite()) {
                hess_q = candidate;
                if let InternalState::Pulay(system) = &state {
                    hess_q = crate::optimizer::redundant_internals::project_redundant_hessian(
                        system,
                        &hess_q,
                        opts.redundant_penalty,
                    )?;
                }
            }
        }

        if scale <= opts.min_line_search_scale && energy > old_energy {
            x = old_x;
            energy = old_energy;
            state = old_state;
            hess_q = old_hess;
            grad_x = vec![0.0; x.len()];
            energy = obj.energy_gradient(&x, &mut grad_x)?;
            grad_q = internal_gradient(&state, &grad_x)?;
        }
    }

    if opts.verbosity > 0 {
        print_optimization_status(
            false,
            "Maximum optimization steps reached without convergence.",
        );
    }

    Ok(InternalMinimumResult {
        x,
        energy,
        gradient: grad_x,
        converged: false,
        cycles: opts.max_cycles,
        message: "Maximum optimization steps reached without convergence.".to_string(),
        coordinates: opts.coordinate_mode,
        trajectory,
        trajectory_energies,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    struct DiatomicBondHarmonic {
        r0: f64,
    }

    impl Objective for DiatomicBondHarmonic {
        fn energy(&mut self, x: &[f64]) -> Result<f64> {
            let dx = x[3] - x[0];
            let dy = x[4] - x[1];
            let dz = x[5] - x[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            Ok(0.5 * (r - self.r0) * (r - self.r0))
        }

        fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
            let dx = x[3] - x[0];
            let dy = x[4] - x[1];
            let dz = x[5] - x[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            let pref = if r > 0.0 { (r - self.r0) / r } else { 0.0 };
            grad[0] = -pref * dx;
            grad[1] = -pref * dy;
            grad[2] = -pref * dz;
            grad[3] = pref * dx;
            grad[4] = pref * dy;
            grad[5] = pref * dz;
            Ok(())
        }
    }

    #[test]
    fn pulay_redundant_internals_minimize_diatomic_bond() {
        let conn = ConnectivityModel {
            rcov_disp: vec![1.0, 1.0],
            kcn: 16.0,
            facmin: 0.1,
        };
        let mut obj = DiatomicBondHarmonic { r0: 1.4 };
        let opts = InternalMinimumOptions {
            verbosity: 0,
            coordinate_mode: InternalCoordinateMode::PulayRedundant,
            max_gradient: 1.0e-8,
            rms_gradient: 1.0e-8,
            max_step: 1.0e-8,
            rms_step: 1.0e-8,
            energy_tol: 1.0e-12,
            ..Default::default()
        };
        let res =
            minimize_internal(vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0], conn, &mut obj, opts).unwrap();
        assert!(res.converged);
        assert!(res.energy < 1.0e-14);
    }
}
