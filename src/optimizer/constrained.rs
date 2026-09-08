//! Smite-style constrained internal-coordinate optimization.
//!
//! This module is separate from the legacy optimizer. It implements the core
//! constrained relaxation used by Smite coordinate scans: selected primitive
//! bond/angle/dihedral coordinates are held at target values, while all other
//! internal coordinates are relaxed with a BFGS model Hessian.

use anyhow::{bail, Result};
use ndarray::Array2;

use crate::optimizer::internal_coords::{
    analyze_structure, best_fit_dq_to_cart, bpg_matrix, compute_internals, diff_internals,
    ConnectivityModel, InternalCoords,
};
use crate::optimizer::minimum_cartesian::step_limit;
use crate::optimizer::redundant_internals::{
    build_redundant_internals, cartesian_gradient_to_redundant, project_redundant_hessian,
    update_cartesian_from_redundant_step, RedundantInternalSystem,
};
use crate::optimizer::traits::Objective;

const HARTREE_TO_KJMOL: f64 = 2625.499638;
// Matches `xyzparser::parse_xyz`'s ANGSTROM_TO_BOHR deliberately (see the
// note there) rather than the more precise CODATA 2018 value, so a
// constraint target given in Angstrom lands on the same bohr geometry the
// rest of the pipeline would parse it to.
const ANGSTROM_TO_BOHR: f64 = 1.8897259886;
const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;
const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstrainedCoordinateMode {
    PrimitiveInternal,
    PulayRedundant,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstraintCoordinate {
    Bond([usize; 2]),
    Angle([usize; 3]),
    Dihedral([usize; 4]),
}

#[derive(Debug, Clone, Copy)]
pub struct ConstraintTarget {
    pub coordinate: ConstraintCoordinate,
    pub target: f64,
}

impl ConstraintTarget {
    pub fn bond_bohr(atoms: [usize; 2], target_bohr: f64) -> Self {
        Self {
            coordinate: ConstraintCoordinate::Bond(atoms),
            target: target_bohr,
        }
    }

    pub fn bond_angstrom(atoms: [usize; 2], target_angstrom: f64) -> Self {
        Self::bond_bohr(atoms, target_angstrom * ANGSTROM_TO_BOHR)
    }

    pub fn angle_radian(atoms: [usize; 3], target_radian: f64) -> Self {
        Self {
            coordinate: ConstraintCoordinate::Angle(atoms),
            target: target_radian,
        }
    }

    pub fn angle_degree(atoms: [usize; 3], target_degree: f64) -> Self {
        Self::angle_radian(atoms, target_degree * DEG_TO_RAD)
    }

    pub fn dihedral_radian(atoms: [usize; 4], target_radian: f64) -> Self {
        Self {
            coordinate: ConstraintCoordinate::Dihedral(atoms),
            target: target_radian,
        }
    }

    pub fn dihedral_degree(atoms: [usize; 4], target_degree: f64) -> Self {
        Self::dihedral_radian(atoms, target_degree * DEG_TO_RAD)
    }
}

#[derive(Debug, Clone)]
pub struct ConstrainedOptimizationOptions {
    pub max_cycles: usize,
    pub energy_tol: f64,
    pub max_step_internal: f64,
    pub max_step: f64,
    pub rms_step: f64,
    pub max_gradient: f64,
    pub rms_gradient: f64,
    pub best_fit_iters: usize,
    pub coordinate_mode: ConstrainedCoordinateMode,
    pub g_cutoff: f64,
    pub redundant_penalty: f64,
    pub verbosity: u8,
}

impl Default for ConstrainedOptimizationOptions {
    fn default() -> Self {
        Self {
            max_cycles: 50,
            energy_tol: 5.0e-6,
            max_step_internal: 0.10,
            max_step: 4.0e-3,
            rms_step: 2.0e-3,
            max_gradient: 3.0e-4,
            rms_gradient: 1.0e-4,
            best_fit_iters: 20,
            coordinate_mode: ConstrainedCoordinateMode::PulayRedundant,
            g_cutoff: 1.0e-7,
            redundant_penalty: 1000.0,
            verbosity: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct ConstrainedOptimizationStep {
    pub cycle: usize,
    pub energy: f64,
    pub energy_change: f64,
    pub max_step: f64,
    pub rms_step: f64,
    pub max_gradient: f64,
    pub rms_gradient: f64,
}

#[derive(Debug, Clone)]
pub struct ConstrainedOptimizationResult {
    pub x: Vec<f64>,
    pub energy: f64,
    pub gradient: Vec<f64>,
    pub converged: bool,
    pub cycles: usize,
    pub message: String,
    pub constraints: Vec<ConstraintTarget>,
    pub final_values: Vec<f64>,
    pub steps: Vec<ConstrainedOptimizationStep>,
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

fn place_centered_text(line: &mut Vec<char>, center: f64, text: &str) {
    let len = line.len();
    if len == 0 || text.is_empty() {
        return;
    }
    let mut start = center.round() as isize - text.len() as isize / 2;
    if start < 0 {
        start = 0;
    }
    if start as usize + text.len() > len {
        start = len.saturating_sub(text.len()) as isize;
    }
    for (offset, ch) in text.chars().enumerate() {
        let idx = start as usize + offset;
        if idx < len {
            line[idx] = ch;
        }
    }
}

fn print_threshold_guide(opts: &ConstrainedOptimizationOptions, header_line: &str) {
    let columns = [
        ("dE[kJ/mol]", opts.energy_tol * HARTREE_TO_KJMOL),
        ("max_step", opts.max_step),
        ("rms_step", opts.rms_step),
        ("max_grad", opts.max_gradient),
        ("rms_grad", opts.rms_gradient),
    ];
    let mut value_line = vec![' '; header_line.len().max(80)];
    let mut bar_line = vec![' '; value_line.len()];
    let mut arrow_line = vec![' '; value_line.len()];
    for (name, value) in columns {
        if let Some(pos) = header_line.find(name) {
            let center = pos as f64 + name.len() as f64 / 2.0;
            place_centered_text(&mut value_line, center, &sci2(value));
            place_centered_text(&mut bar_line, center, "|");
            place_centered_text(&mut arrow_line, center, "V");
        }
    }
    println!(
        "Threshold Values:{}",
        value_line.into_iter().collect::<String>().trim_end()
    );
    println!(
        "                 {}",
        bar_line.into_iter().collect::<String>().trim_end()
    );
    println!(
        "                 {}",
        arrow_line.into_iter().collect::<String>().trim_end()
    );
}

fn print_header(opts: &ConstrainedOptimizationOptions, constraints: &[ConstraintTarget]) {
    println!("Constrained internal-coordinate optimization");
    println!("Coordinate mode: {:?}", opts.coordinate_mode);
    println!("Fixed coordinates: {}", constraints.len());
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
    print_threshold_guide(opts, &header_line);
    println!("{header_line}");
}

fn print_step(step: &ConstrainedOptimizationStep) {
    println!(
        "{:5} {:18.8} {:18.6} {:11.2e} {:11.2e} {:11.2e} {:11.2e}",
        step.cycle,
        step.energy,
        step.energy_change * HARTREE_TO_KJMOL,
        step.max_step,
        step.rms_step,
        step.max_gradient,
        step.rms_gradient
    );
}

fn wrap_angle_delta(delta: f64) -> f64 {
    let mut d = delta;
    if d > std::f64::consts::PI {
        d -= 2.0 * std::f64::consts::PI;
    }
    if d < -std::f64::consts::PI {
        d += 2.0 * std::f64::consts::PI;
    }
    d
}

fn target_delta(coordinate: ConstraintCoordinate, target: f64, current: f64) -> f64 {
    match coordinate {
        ConstraintCoordinate::Dihedral(_) => wrap_angle_delta(target - current),
        _ => target - current,
    }
}

fn find_constraint_index(ic: &InternalCoords, constraint: ConstraintCoordinate) -> Result<usize> {
    match constraint {
        ConstraintCoordinate::Bond([a, b]) => {
            for (i, &[x, y]) in ic.bonds.iter().enumerate() {
                if (a == x && b == y) || (a == y && b == x) {
                    return Ok(i);
                }
            }
            bail!("constrained bond {a}-{b} is not present in the internal coordinate set")
        }
        ConstraintCoordinate::Angle([a, b, c]) => {
            let offset = ic.nbonds();
            for (i, &[x, y, z]) in ic.angles.iter().enumerate() {
                if b == y && ((a == x && c == z) || (a == z && c == x)) {
                    return Ok(offset + i);
                }
            }
            bail!("constrained angle {a}-{b}-{c} is not present in the internal coordinate set")
        }
        ConstraintCoordinate::Dihedral([a, b, c, d]) => {
            let offset = ic.nbonds() + ic.nangles();
            for (i, &[w, x, y, z]) in ic.dihedrals.iter().enumerate() {
                if (a == w && b == x && c == y && d == z) || (a == z && b == y && c == x && d == w)
                {
                    return Ok(offset + i);
                }
            }
            bail!("constrained dihedral {a}-{b}-{c}-{d} is not present in the internal coordinate set")
        }
    }
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

fn identity(n: usize) -> Array2<f64> {
    let mut out = Array2::<f64>::zeros((n, n));
    for i in 0..n {
        out[(i, i)] = 1.0;
    }
    out
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

fn bfgs_update_hessian(hess: &Array2<f64>, step: &[f64], y: &[f64]) -> Array2<f64> {
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

#[derive(Debug, Clone)]
enum ConstraintState {
    Primitive {
        ic: InternalCoords,
        qs: Vec<f64>,
        b: Vec<Vec<f64>>,
        inv_g: Vec<Vec<f64>>,
    },
    Pulay(RedundantInternalSystem),
}

impl ConstraintState {
    fn coordinates(&self) -> &InternalCoords {
        match self {
            Self::Primitive { ic, .. } => ic,
            Self::Pulay(system) => &system.coordinates,
        }
    }

    fn q(&self) -> &[f64] {
        match self {
            Self::Primitive { qs, .. } => qs,
            Self::Pulay(system) => &system.q,
        }
    }
}

fn build_constraint_state(
    x: &[f64],
    ic: &InternalCoords,
    mode: ConstrainedCoordinateMode,
    g_cutoff: f64,
) -> Result<ConstraintState> {
    match mode {
        ConstrainedCoordinateMode::PrimitiveInternal => {
            let qs = compute_internals(x, ic);
            let bpg = bpg_matrix(x, ic)?;
            Ok(ConstraintState::Primitive {
                ic: ic.clone(),
                qs,
                b: bpg.b,
                inv_g: bpg.inv_g,
            })
        }
        ConstrainedCoordinateMode::PulayRedundant => Ok(ConstraintState::Pulay(
            build_redundant_internals(x, None, Some(ic.clone()), None, g_cutoff)?,
        )),
    }
}

fn state_gradient(state: &ConstraintState, grad_x: &[f64]) -> Result<Vec<f64>> {
    match state {
        ConstraintState::Primitive { b, inv_g, .. } => Ok(internal_gradient(b, inv_g, grad_x)),
        ConstraintState::Pulay(system) => cartesian_gradient_to_redundant(system, grad_x),
    }
}

fn apply_state_step(
    x: &[f64],
    state: &ConstraintState,
    dq: &[f64],
    best_fit_iters: usize,
) -> Result<Vec<f64>> {
    match state {
        ConstraintState::Primitive { ic, qs, b, inv_g } => {
            let mut trial = x.to_vec();
            let bpg = crate::optimizer::internal_coords::Bpg {
                b: b.clone(),
                inv_g: inv_g.clone(),
            };
            best_fit_dq_to_cart(&mut trial, qs, dq, ic, &bpg, best_fit_iters)?;
            Ok(trial)
        }
        ConstraintState::Pulay(system) => {
            update_cartesian_from_redundant_step(x, system, dq, best_fit_iters)
        }
    }
}

fn apply_constraint_targets(
    x: &mut [f64],
    ic: &InternalCoords,
    constraints: &[ConstraintTarget],
    indices: &[usize],
    best_fit_iters: usize,
) -> Result<()> {
    let qs = compute_internals(x, ic);
    let bpg = bpg_matrix(x, ic)?;
    let mut dq = vec![0.0; ic.nint()];
    for (constraint, &idx) in constraints.iter().zip(indices.iter()) {
        dq[idx] = target_delta(constraint.coordinate, constraint.target, qs[idx]);
    }
    best_fit_dq_to_cart(x, &qs, &dq, ic, &bpg, best_fit_iters)
}

pub fn optimize_constrained_internals<O: Objective>(
    obj: &mut O,
    x0: Vec<f64>,
    model: &ConnectivityModel,
    constraints: Vec<ConstraintTarget>,
    opts: ConstrainedOptimizationOptions,
) -> Result<ConstrainedOptimizationResult> {
    if constraints.is_empty() {
        bail!("constrained optimization requires at least one fixed coordinate");
    }
    let ic = analyze_structure(&x0, model)?;
    optimize_constrained_with_internals(obj, x0, ic, constraints, opts)
}

pub fn optimize_constrained_with_internals<O: Objective>(
    obj: &mut O,
    mut x: Vec<f64>,
    ic: InternalCoords,
    constraints: Vec<ConstraintTarget>,
    opts: ConstrainedOptimizationOptions,
) -> Result<ConstrainedOptimizationResult> {
    if x.len() != 3 * ic.nat {
        bail!("coordinate length does not match internal coordinate atom count");
    }
    let nint = ic.nint();
    let indices: Vec<usize> = constraints
        .iter()
        .map(|c| find_constraint_index(&ic, c.coordinate))
        .collect::<Result<Vec<_>>>()?;
    let mut constrained = vec![false; nint];
    for &idx in &indices {
        constrained[idx] = true;
    }
    let free_indices: Vec<usize> = (0..nint).filter(|&i| !constrained[i]).collect();
    if free_indices.is_empty() {
        bail!("all internal coordinates are constrained; nothing can relax");
    }

    apply_constraint_targets(&mut x, &ic, &constraints, &indices, opts.best_fit_iters)?;
    let mut grad_x = vec![0.0; x.len()];
    let mut energy = obj.energy_gradient(&x, &mut grad_x)?;
    let mut hess_q = identity(nint);
    let mut steps = Vec::new();

    if opts.verbosity > 0 {
        print_header(&opts, &constraints);
    }

    if opts.coordinate_mode == ConstrainedCoordinateMode::PulayRedundant {
        let state = build_constraint_state(&x, &ic, opts.coordinate_mode, opts.g_cutoff)?;
        hess_q = project_redundant_hessian(
            match &state {
                ConstraintState::Pulay(system) => system,
                _ => unreachable!(),
            },
            &hess_q,
            opts.redundant_penalty,
        )?;
    }

    let mut message = "Maximum constrained optimization cycles reached.".to_string();
    for cycle in 1..=opts.max_cycles {
        let state = build_constraint_state(&x, &ic, opts.coordinate_mode, opts.g_cutoff)?;
        let qs = state.q().to_vec();
        let grad_q = state_gradient(&state, &grad_x)?;
        let grad_free: Vec<f64> = free_indices.iter().map(|&i| grad_q[i]).collect();
        let mut hess_free = Array2::<f64>::zeros((free_indices.len(), free_indices.len()));
        for (ii, &i) in free_indices.iter().enumerate() {
            for (jj, &j) in free_indices.iter().enumerate() {
                hess_free[(ii, jj)] = hess_q[(i, j)];
            }
        }
        let rhs: Vec<f64> = grad_free.iter().map(|g| -g).collect();
        let mut dq_free = solve_linear(&hess_free, &rhs);
        step_limit(&mut dq_free, opts.max_step_internal);

        let mut dq = vec![0.0; nint];
        for (k, &idx) in free_indices.iter().enumerate() {
            dq[idx] = dq_free[k];
        }
        for (constraint, &idx) in constraints.iter().zip(indices.iter()) {
            dq[idx] = target_delta(constraint.coordinate, constraint.target, qs[idx]);
        }

        let trial_x = apply_state_step(&x, &state, &dq, opts.best_fit_iters)?;
        let mut trial_grad_x = vec![0.0; trial_x.len()];
        let trial_energy = obj.energy_gradient(&trial_x, &mut trial_grad_x)?;
        let trial_state =
            build_constraint_state(&trial_x, &ic, opts.coordinate_mode, opts.g_cutoff)?;
        let trial_grad_q = state_gradient(&trial_state, &trial_grad_x)?;
        let accepted_dq = diff_internals(trial_state.q(), state.q(), ic.ndiheds());
        let y: Vec<f64> = trial_grad_q
            .iter()
            .zip(grad_q.iter())
            .map(|(a, b)| a - b)
            .collect();
        hess_q = bfgs_update_hessian(&hess_q, &accepted_dq, &y);
        if let ConstraintState::Pulay(system) = &trial_state {
            hess_q = project_redundant_hessian(system, &hess_q, opts.redundant_penalty)?;
        }

        let cart_step: Vec<f64> = trial_x.iter().zip(x.iter()).map(|(a, b)| a - b).collect();
        let energy_change = trial_energy - energy;
        let grad_free_new: Vec<f64> = free_indices.iter().map(|&i| trial_grad_q[i]).collect();
        let step = ConstrainedOptimizationStep {
            cycle,
            energy: trial_energy,
            energy_change,
            max_step: max_abs(&cart_step),
            rms_step: rms(&cart_step),
            max_gradient: max_abs(&grad_free_new),
            rms_gradient: rms(&grad_free_new),
        };
        if opts.verbosity > 0 {
            print_step(&step);
        }

        x = trial_x;
        energy = trial_energy;
        grad_x = trial_grad_x;
        let converged = energy_change.abs() <= opts.energy_tol
            && step.max_gradient <= opts.max_gradient
            && step.rms_gradient <= opts.rms_gradient
            && step.max_step <= opts.max_step
            && step.rms_step <= opts.rms_step;
        steps.push(step);
        if converged {
            message = "Constrained optimization converged.".to_string();
            let final_values = final_constraint_values(&x, &ic, &constraints, &indices);
            return Ok(ConstrainedOptimizationResult {
                x,
                energy,
                gradient: grad_x,
                converged: true,
                cycles: cycle,
                message,
                constraints,
                final_values,
                steps,
            });
        }
    }

    let final_values = final_constraint_values(&x, &ic, &constraints, &indices);
    Ok(ConstrainedOptimizationResult {
        x,
        energy,
        gradient: grad_x,
        converged: false,
        cycles: opts.max_cycles,
        message,
        constraints,
        final_values,
        steps,
    })
}

fn final_constraint_values(
    x: &[f64],
    ic: &InternalCoords,
    constraints: &[ConstraintTarget],
    indices: &[usize],
) -> Vec<f64> {
    let qs = compute_internals(x, ic);
    constraints
        .iter()
        .zip(indices.iter())
        .map(|(constraint, &idx)| match constraint.coordinate {
            ConstraintCoordinate::Bond(_) => qs[idx] / ANGSTROM_TO_BOHR,
            ConstraintCoordinate::Angle(_) | ConstraintCoordinate::Dihedral(_) => {
                qs[idx] * RAD_TO_DEG
            }
        })
        .collect()
}

#[derive(Debug, Clone)]
pub struct RelaxedScanOptions {
    pub targets: Vec<f64>,
    pub coordinate_mode: ConstrainedCoordinateMode,
    pub optimizer: ConstrainedOptimizationOptions,
    pub atom_labels: Option<Vec<String>>,
    pub trajectory_file: Option<String>,
    pub verbosity: u8,
}

impl Default for RelaxedScanOptions {
    fn default() -> Self {
        Self {
            targets: Vec::new(),
            coordinate_mode: ConstrainedCoordinateMode::PulayRedundant,
            optimizer: ConstrainedOptimizationOptions::default(),
            atom_labels: None,
            trajectory_file: Some("relaxed_constraint_scan.xyz".to_string()),
            verbosity: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct RelaxedScanPoint {
    pub index: usize,
    pub target_value: f64,
    pub actual_value: f64,
    pub energy: f64,
    pub relative_energy_kcal_mol: f64,
    pub converged: bool,
    pub optimization_steps: usize,
    pub x: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct RelaxedScanResult {
    pub coordinate: ConstraintCoordinate,
    pub points: Vec<RelaxedScanPoint>,
    pub trajectory_file: Option<String>,
}

impl RelaxedScanResult {
    pub fn print_table(&self) {
        println!("{}", self.format_table());
    }

    pub fn format_table(&self) -> String {
        let mut out = String::new();
        let value_header = coordinate_value_header(self.coordinate);
        out.push_str("Relaxed constrained scan summary\n");
        out.push_str(&format!(
            "{:>5} {:>14} {:>18} {:>18} {:>8} {:>8}\n",
            "step", value_header, "E[Eh]", "dE[kcal/mol]", "opt", "conv"
        ));
        for point in &self.points {
            out.push_str(&format!(
                "{:5} {:14.6} {:18.10} {:18.6} {:8} {:>8}\n",
                point.index,
                point.actual_value,
                point.energy,
                point.relative_energy_kcal_mol,
                point.optimization_steps,
                point.converged
            ));
        }
        out
    }
}

fn coordinate_value_header(coordinate: ConstraintCoordinate) -> &'static str {
    match coordinate {
        ConstraintCoordinate::Bond(_) => "R[Ang]",
        ConstraintCoordinate::Angle(_) => "angle[deg]",
        ConstraintCoordinate::Dihedral(_) => "dihedral[deg]",
    }
}

fn coordinate_display_value(coordinate: ConstraintCoordinate, internal_value: f64) -> f64 {
    match coordinate {
        ConstraintCoordinate::Bond(_) => internal_value / ANGSTROM_TO_BOHR,
        ConstraintCoordinate::Angle(_) | ConstraintCoordinate::Dihedral(_) => {
            internal_value * RAD_TO_DEG
        }
    }
}

fn target_from_internal(coordinate: ConstraintCoordinate, target: f64) -> ConstraintTarget {
    ConstraintTarget { coordinate, target }
}

fn xyz_frame(atoms: &[String], x: &[f64], comment: &str) -> String {
    let bohr_to_angstrom = 1.0 / ANGSTROM_TO_BOHR;
    let mut out = format!("{}\n{}\n", atoms.len(), comment);
    for (i, atom) in atoms.iter().enumerate() {
        out.push_str(&format!(
            "{:<2} {:16.10} {:16.10} {:16.10}\n",
            atom,
            x[3 * i] * bohr_to_angstrom,
            x[3 * i + 1] * bohr_to_angstrom,
            x[3 * i + 2] * bohr_to_angstrom
        ));
    }
    out
}

pub fn relaxed_constraint_scan<O: Objective>(
    obj: &mut O,
    x0: Vec<f64>,
    model: &ConnectivityModel,
    coordinate: ConstraintCoordinate,
    opts: RelaxedScanOptions,
) -> Result<RelaxedScanResult> {
    if opts.targets.is_empty() {
        bail!("relaxed constrained scan requires at least one target value");
    }
    if let Some(atoms) = &opts.atom_labels {
        if atoms.len() * 3 != x0.len() {
            bail!("relaxed scan atom_labels length does not match coordinate vector");
        }
    }
    let ic = analyze_structure(&x0, model)?;
    find_constraint_index(&ic, coordinate)?;
    let mut x = x0;
    let mut points = Vec::with_capacity(opts.targets.len());
    let mut reference_energy = None;
    let mut traj = String::new();

    if opts.verbosity > 0 {
        println!(
            "Relaxed constrained scan: {}",
            coordinate_value_header(coordinate)
        );
        println!("Scan points: {}", opts.targets.len());
    }

    for (i, &target) in opts.targets.iter().enumerate() {
        let step_index = i + 1;
        let target_display = coordinate_display_value(coordinate, target);
        if opts.verbosity > 0 {
            println!(
                "scan step {:4}/{:<4} setting {}={:10.5}",
                step_index,
                opts.targets.len(),
                coordinate_value_header(coordinate),
                target_display
            );
        }
        let constraint = target_from_internal(coordinate, target);
        let mut opt_opts = opts.optimizer.clone();
        opt_opts.coordinate_mode = opts.coordinate_mode;
        opt_opts.verbosity = opts.optimizer.verbosity.min(opts.verbosity);
        let result =
            optimize_constrained_with_internals(obj, x, ic.clone(), vec![constraint], opt_opts)?;
        x = result.x.clone();
        let actual_value = result.final_values[0];
        let first = *reference_energy.get_or_insert(result.energy);
        let rel = (result.energy - first) * 627.509_474_0631;
        if let Some(atoms) = &opts.atom_labels {
            traj.push_str(&xyz_frame(
                atoms,
                &x,
                &format!(
                    "scan_step={} target={:.10} actual={:.10} E={:.12} converged={} opt_steps={}",
                    step_index,
                    target_display,
                    actual_value,
                    result.energy,
                    result.converged,
                    result.cycles
                ),
            ));
        }
        if opts.verbosity > 0 {
            println!(
                "scan     {:4} {}={:10.5} E={:18.10} dE_scan={:12.5} kcal/mol opt={:3} conv={}",
                step_index,
                coordinate_value_header(coordinate),
                actual_value,
                result.energy,
                rel,
                result.cycles,
                result.converged
            );
        }
        points.push(RelaxedScanPoint {
            index: step_index,
            target_value: target_display,
            actual_value,
            energy: result.energy,
            relative_energy_kcal_mol: rel,
            converged: result.converged,
            optimization_steps: result.cycles,
            x: x.clone(),
        });
    }

    if let Some(file) = &opts.trajectory_file {
        if opts.atom_labels.is_some() {
            std::fs::write(file, traj)?;
            if opts.verbosity > 0 {
                println!("Relaxed constrained scan trajectory written to {file}");
            }
        }
    }

    let result = RelaxedScanResult {
        coordinate,
        points,
        trajectory_file: opts.trajectory_file.clone(),
    };
    if opts.verbosity > 0 {
        println!();
        result.print_table();
    }
    Ok(result)
}
