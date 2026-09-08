//! Spin-crossing / MECP optimizer.
//!
//! This is the Rust counterpart of Smite's `optimizer.spincrossing`: two
//! independent objective surfaces are evaluated at the same geometry, an
//! effective MECP gradient is built from the two gradients and the energy gap,
//! and a quasi-Newton optimizer searches the lowest point where the two
//! surfaces cross.

use anyhow::{bail, Context, Result};
use ndarray::{Array1, Array2};
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::optimizer::internal_coords::{
    analyze_structure, best_fit_dq_to_cart, bpg_matrix, compute_internals, ConnectivityModel,
    InternalCoords,
};
use crate::optimizer::minimum_cartesian::step_limit;
use crate::optimizer::redundant_internals::{
    build_redundant_internals, cartesian_gradient_to_redundant,
    update_cartesian_from_redundant_step, RedundantInternalSystem,
};
use crate::optimizer::traits::Objective;

const HARTREE_TO_KJMOL: f64 = 2625.499638;
const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903; // CODATA 2018

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpinCrossingMethod {
    Bfgs,
    PlainBfgs,
    BfgsOld,
    Dfp,
    Sr1,
    Psb,
    Bofill,
    Broyden,
    BbGrad1,
    BbGrad2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpinCrossingCoordinateMode {
    Cartesian,
    SimpleInternal,
    PulayRedundant,
}

#[derive(Debug, Clone)]
pub struct SpinCrossingOptions {
    pub method: SpinCrossingMethod,
    pub coordinate_mode: SpinCrossingCoordinateMode,
    pub multiplicities: Option<(usize, usize)>,
    pub max_cycles: usize,
    pub energy_tol: f64,
    pub max_step: f64,
    pub rms_step: f64,
    pub max_gradient: f64,
    pub rms_gradient: f64,
    pub step_max_component: f64,
    pub initial_inverse_hessian: f64,
    pub max_step_internal: f64,
    pub best_fit_iters: usize,
    pub g_cutoff: f64,
    pub fac_parallel: f64,
    pub fac_perpendicular: f64,
    pub trajectory_file: Option<String>,
    pub atom_labels: Option<Vec<String>>,
    pub verbosity: u8,
}

impl Default for SpinCrossingOptions {
    fn default() -> Self {
        Self {
            method: SpinCrossingMethod::Bfgs,
            coordinate_mode: SpinCrossingCoordinateMode::PulayRedundant,
            multiplicities: None,
            max_cycles: 40,
            energy_tol: 5.0e-5,
            max_step: 4.0e-3,
            rms_step: 2.5e-3,
            max_gradient: 7.0e-4,
            rms_gradient: 5.0e-4,
            step_max_component: 0.1,
            initial_inverse_hessian: 0.7,
            max_step_internal: 0.2,
            best_fit_iters: 5,
            g_cutoff: 1.0e-7,
            fac_parallel: 1.0,
            fac_perpendicular: 140.0,
            trajectory_file: Some("spincross_traj.xyz".to_string()),
            atom_labels: None,
            verbosity: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SpinCrossingPoint {
    pub cycle: usize,
    pub x: Vec<f64>,
    pub energy_a: f64,
    pub energy_b: f64,
    pub gradient_a: Vec<f64>,
    pub gradient_b: Vec<f64>,
    pub effective_gradient: Vec<f64>,
    pub parallel_gradient: Vec<f64>,
    pub perpendicular_gradient: Vec<f64>,
    pub step: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct SpinCrossingResult {
    pub x: Vec<f64>,
    pub energy_a: f64,
    pub energy_b: f64,
    pub gradient_a: Vec<f64>,
    pub gradient_b: Vec<f64>,
    pub effective_gradient: Vec<f64>,
    pub converged: bool,
    pub cycles: usize,
    pub method: SpinCrossingMethod,
    pub coordinates: SpinCrossingCoordinateMode,
    pub multiplicities: Option<(usize, usize)>,
    pub trajectory_file: Option<String>,
    pub points: Vec<SpinCrossingPoint>,
    pub message: String,
}

impl SpinCrossingResult {
    pub fn energy_gap(&self) -> f64 {
        self.energy_a - self.energy_b
    }
}

#[derive(Debug, Clone)]
struct SimpleInternalState {
    ic: InternalCoords,
    q: Vec<f64>,
    bpg: crate::optimizer::internal_coords::Bpg,
}

#[derive(Debug, Clone)]
enum InternalState {
    Simple(SimpleInternalState),
    Pulay(RedundantInternalSystem),
}

impl InternalState {
    fn nint(&self) -> usize {
        match self {
            Self::Simple(s) => s.ic.nint(),
            Self::Pulay(s) => s.coordinates.nint(),
        }
    }

    fn coordinates(&self) -> &InternalCoords {
        match self {
            Self::Simple(s) => &s.ic,
            Self::Pulay(s) => &s.coordinates,
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
    a.iter().fold(0.0, |m, v| m.max(v.abs()))
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
        sci2(energy_tol * HARTREE_TO_KJMOL),
        sci2(max_step),
        sci2(rms_step),
        sci2(max_gradient),
        sci2(rms_gradient),
    ];
    let columns = ["dE[kJ/mol]", "max_step", "rms_step", "max_geff", "rms_geff"];
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

fn print_header(opts: &SpinCrossingOptions) {
    println!(" Energy Gap:      {} Eh", sci2(opts.energy_tol));
    println!(" Max. Eff. Gradient: {} Eh/bohr", sci2(opts.max_gradient));
    println!(" RMS Eff. Gradient:  {} Eh/bohr", sci2(opts.rms_gradient));
    println!(" Max Step:        {} bohr", sci2(opts.max_step));
    println!(" RMS Step:        {} bohr", sci2(opts.rms_step));
    println!();
    let header_line = format!(
        "{:<5} {:>18} {:>18} {:>18} {:>11} {:>11} {:>11} {:>11}",
        "step", "Ea[Eh]", "Eb[Eh]", "dE[kJ/mol]", "max_step", "rms_step", "max_geff", "rms_geff"
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

fn outer(a: &[f64], b: &[f64]) -> Array2<f64> {
    let mut out = Array2::<f64>::zeros((a.len(), b.len()));
    for i in 0..a.len() {
        for j in 0..b.len() {
            out[(i, j)] = a[i] * b[j];
        }
    }
    out
}

fn symmetrize(a: &mut Array2<f64>) {
    let (n, m) = a.dim();
    if n != m {
        return;
    }
    for i in 0..n {
        for j in 0..i {
            let v = 0.5 * (a[(i, j)] + a[(j, i)]);
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
}

fn bfgs_inverse_update(inv: &Array2<f64>, dx: &[f64], dg: &[f64]) -> (Array2<f64>, bool) {
    let overlap = dot(dg, dx);
    let scale = norm(dx) * norm(dg);
    if overlap <= 1.0e-14_f64.max(1.0e-8 * scale) {
        return (inv.clone(), false);
    }
    let n = dx.len();
    let rho = 1.0 / overlap;
    let eye = Array2::<f64>::eye(n);
    let xg = outer(dx, dg);
    let gx = outer(dg, dx);
    let left = &eye - &(xg * rho);
    let right = &eye - &(gx * rho);
    let mut updated = left.dot(inv).dot(&right) + outer(dx, dx) * rho;
    if updated.iter().all(|v| v.is_finite()) {
        symmetrize(&mut updated);
        (updated, true)
    } else {
        (inv.clone(), false)
    }
}

fn dfp_inverse_update(inv: &Array2<f64>, dx: &[f64], dg: &[f64]) -> (Array2<f64>, bool) {
    let h_dg = mat_vec(inv, dg);
    let overlap = dot(dx, dg);
    let ghg = dot(dg, &h_dg);
    if overlap.abs() <= 1.0e-14 || ghg.abs() <= 1.0e-14 {
        return (inv.clone(), false);
    }
    let mut updated =
        inv + &(outer(dx, dx) * (1.0 / overlap)) - &(outer(&h_dg, &h_dg) * (1.0 / ghg));
    if updated.iter().all(|v| v.is_finite()) {
        symmetrize(&mut updated);
        (updated, true)
    } else {
        (inv.clone(), false)
    }
}

fn sr1_inverse_update(inv: &Array2<f64>, dx: &[f64], dg: &[f64]) -> (Array2<f64>, bool) {
    let h_dg = mat_vec(inv, dg);
    let y: Vec<f64> = dx.iter().zip(h_dg.iter()).map(|(x, h)| x - h).collect();
    let overlap = dot(&y, dg);
    let scale = norm(&y) * norm(dg);
    if overlap.abs() <= 1.0e-14_f64.max(1.0e-8 * scale) {
        return (inv.clone(), false);
    }
    let mut updated = inv + &(outer(&y, &y) * (1.0 / overlap));
    if updated.iter().all(|v| v.is_finite()) {
        symmetrize(&mut updated);
        (updated, true)
    } else {
        (inv.clone(), false)
    }
}

fn psb_inverse_update(inv: &Array2<f64>, dx: &[f64], dg: &[f64]) -> (Array2<f64>, bool) {
    let h_dg = mat_vec(inv, dg);
    let y: Vec<f64> = dx.iter().zip(h_dg.iter()).map(|(x, h)| x - h).collect();
    let dd = dot(dg, dg);
    let yd = dot(&y, dg);
    if dd <= 1.0e-14 {
        return (inv.clone(), false);
    }
    let mut updated =
        inv + &((outer(&y, dg) + outer(dg, &y)) * (1.0 / dd)) - &(outer(dg, dg) * (yd / (dd * dd)));
    if updated.iter().all(|v| v.is_finite()) {
        symmetrize(&mut updated);
        (updated, true)
    } else {
        (inv.clone(), false)
    }
}

fn broyden_inverse_update(inv: &Array2<f64>, dx: &[f64], dg: &[f64]) -> (Array2<f64>, bool) {
    let h_dg = mat_vec(inv, dg);
    let y: Vec<f64> = dx.iter().zip(h_dg.iter()).map(|(x, h)| x - h).collect();
    let dd = dot(dg, dg);
    if dd <= 1.0e-14 {
        return (inv.clone(), false);
    }
    let updated = inv + &(outer(&y, dg) * (1.0 / dd));
    if updated.iter().all(|v| v.is_finite()) {
        (updated, true)
    } else {
        (inv.clone(), false)
    }
}

fn bofill_inverse_update(inv: &Array2<f64>, dx: &[f64], dg: &[f64]) -> (Array2<f64>, bool) {
    let (sr1, ok_sr1) = sr1_inverse_update(inv, dx, dg);
    let (psb, ok_psb) = psb_inverse_update(inv, dx, dg);
    if !ok_sr1 {
        return (psb, ok_psb);
    }
    if !ok_psb {
        return (sr1, ok_sr1);
    }
    let h_dg = mat_vec(inv, dg);
    let y: Vec<f64> = dx.iter().zip(h_dg.iter()).map(|(x, h)| x - h).collect();
    let yy = dot(&y, &y);
    let dd = dot(dg, dg);
    let yd = dot(&y, dg);
    let phi = if yy > 1.0e-24 && dd > 1.0e-24 {
        ((yd * yd) / (yy * dd)).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let mut updated = sr1 * phi + psb * (1.0 - phi);
    symmetrize(&mut updated);
    (updated, true)
}

fn bb_step(eq: usize, x1: &[f64], x2: &[f64], g1: &[f64], g2: &[f64]) -> Vec<f64> {
    let dg: Vec<f64> = g2.iter().zip(g1.iter()).map(|(a, b)| a - b).collect();
    let dx: Vec<f64> = x2.iter().zip(x1.iter()).map(|(a, b)| a - b).collect();
    let xg = dot(&dg, &dx);
    let alpha = if eq == 1 {
        let denom = dot(&dg, &dg);
        if denom.abs() > 1.0e-16 {
            xg / denom
        } else {
            0.7
        }
    } else if xg.abs() > 1.0e-16 {
        dot(&dx, &dx) / xg
    } else {
        0.7
    };
    g2.iter().map(|g| -alpha * g).collect()
}

fn quasi_newton_step(
    method: SpinCrossingMethod,
    x1: &[f64],
    x2: &[f64],
    g1: &[f64],
    g2: &[f64],
    inv_hess: &Array2<f64>,
) -> (Vec<f64>, Array2<f64>, bool) {
    let dx: Vec<f64> = x2.iter().zip(x1.iter()).map(|(a, b)| a - b).collect();
    let dg: Vec<f64> = g2.iter().zip(g1.iter()).map(|(a, b)| a - b).collect();
    let (new_inv, updated) = match method {
        SpinCrossingMethod::Bfgs | SpinCrossingMethod::PlainBfgs => {
            bfgs_inverse_update(inv_hess, &dx, &dg)
        }
        SpinCrossingMethod::BfgsOld | SpinCrossingMethod::Dfp => {
            dfp_inverse_update(inv_hess, &dx, &dg)
        }
        SpinCrossingMethod::Sr1 => sr1_inverse_update(inv_hess, &dx, &dg),
        SpinCrossingMethod::Psb => psb_inverse_update(inv_hess, &dx, &dg),
        SpinCrossingMethod::Bofill => bofill_inverse_update(inv_hess, &dx, &dg),
        SpinCrossingMethod::Broyden => broyden_inverse_update(inv_hess, &dx, &dg),
        SpinCrossingMethod::BbGrad1 => return (bb_step(1, x1, x2, g1, g2), inv_hess.clone(), true),
        SpinCrossingMethod::BbGrad2 => return (bb_step(2, x1, x2, g1, g2), inv_hess.clone(), true),
    };
    let step = mat_vec(&new_inv, g2).into_iter().map(|v| -v).collect();
    (step, new_inv, updated)
}

fn effective_gradient(
    energy_a: f64,
    energy_b: f64,
    grad_a: &[f64],
    grad_b: &[f64],
    fac_parallel: f64,
    fac_perpendicular: f64,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    if grad_a.len() != grad_b.len() {
        bail!("state gradient length mismatch");
    }
    let dgrad: Vec<f64> = grad_a
        .iter()
        .zip(grad_b.iter())
        .map(|(a, b)| a - b)
        .collect();
    let denom = dot(&dgrad, &dgrad);
    if denom <= 1.0e-24 {
        bail!("cannot build spin-crossing effective gradient because state gradients are parallel/equal");
    }
    let scale = dot(grad_a, &dgrad) / denom;
    let gap = energy_a - energy_b;
    let perpendicular: Vec<f64> = dgrad.iter().map(|g| gap * g).collect();
    let parallel: Vec<f64> = grad_a
        .iter()
        .zip(dgrad.iter())
        .map(|(g, d)| g - scale * d)
        .collect();
    let effective: Vec<f64> = parallel
        .iter()
        .zip(perpendicular.iter())
        .map(|(p, q)| fac_parallel * p + fac_perpendicular * q)
        .collect();
    Ok((parallel, perpendicular, effective))
}

fn simple_internal_gradient(
    bpg: &crate::optimizer::internal_coords::Bpg,
    grad_x: &[f64],
) -> Vec<f64> {
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
        InternalState::Simple(s) => Ok(simple_internal_gradient(&s.bpg, grad_x)),
        InternalState::Pulay(s) => cartesian_gradient_to_redundant(s, grad_x),
    }
}

fn build_state(
    x: &[f64],
    mode: SpinCrossingCoordinateMode,
    model: Option<&ConnectivityModel>,
    coordinates: Option<InternalCoords>,
    g_cutoff: f64,
) -> Result<Option<InternalState>> {
    match mode {
        SpinCrossingCoordinateMode::Cartesian => Ok(None),
        SpinCrossingCoordinateMode::SimpleInternal => {
            let ic = if let Some(ic) = coordinates {
                ic
            } else {
                analyze_structure(
                    x,
                    model.ok_or_else(|| {
                        anyhow::anyhow!("connectivity model required for internal spin-crossing")
                    })?,
                )?
            };
            let q = compute_internals(x, &ic);
            let bpg = bpg_matrix(x, &ic)?;
            Ok(Some(InternalState::Simple(SimpleInternalState {
                ic,
                q,
                bpg,
            })))
        }
        SpinCrossingCoordinateMode::PulayRedundant => Ok(Some(InternalState::Pulay(
            build_redundant_internals(x, model, coordinates, None, g_cutoff)?,
        ))),
    }
}

fn opt_coordinate(state: &Option<InternalState>, x: &[f64]) -> Vec<f64> {
    match state {
        None => x.to_vec(),
        Some(InternalState::Simple(s)) => s.q.clone(),
        Some(InternalState::Pulay(s)) => s.q.clone(),
    }
}

fn apply_internal_step(
    x: &[f64],
    state: &InternalState,
    dq: &[f64],
    opts: &SpinCrossingOptions,
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

fn limit_internal_step(dq: &mut [f64], max_step_internal: f64) {
    let n = norm(dq);
    if n > max_step_internal && n > 0.0 {
        let scale = max_step_internal / n;
        for v in dq.iter_mut() {
            *v *= scale;
        }
    }
}

fn converged(gap: f64, step: &[f64], geff: &[f64], opts: &SpinCrossingOptions) -> bool {
    gap.abs() <= opts.energy_tol
        && max_abs(step) <= opts.max_step
        && rms(step) <= opts.rms_step
        && max_abs(geff) <= opts.max_gradient
        && rms(geff) <= opts.rms_gradient
}

fn write_frame(
    writer: &mut BufWriter<File>,
    labels: &[String],
    x: &[f64],
    cycle: usize,
    energy_a: f64,
    energy_b: f64,
    geff: &[f64],
) -> Result<()> {
    writeln!(writer, "{}", labels.len())?;
    writeln!(
        writer,
        "spin-crossing step={} Ea={:.12} Eb={:.12} gap={:.12} geff_rms={:.6e}",
        cycle,
        energy_a,
        energy_b,
        energy_a - energy_b,
        rms(geff)
    )?;
    for i in 0..labels.len() {
        let k = 3 * i;
        writeln!(
            writer,
            "{:<2} {:16.10} {:16.10} {:16.10}",
            labels[i],
            x[k] * BOHR_TO_ANGSTROM,
            x[k + 1] * BOHR_TO_ANGSTROM,
            x[k + 2] * BOHR_TO_ANGSTROM
        )?;
    }
    Ok(())
}

fn evaluate_pair<A: Objective, B: Objective>(
    obj_a: &mut A,
    obj_b: &mut B,
    x: &[f64],
) -> Result<(f64, f64, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> {
    let mut ga = vec![0.0; x.len()];
    let mut gb = vec![0.0; x.len()];
    let ea = obj_a
        .energy_gradient(x, &mut ga)
        .context("spin-crossing state A energy/gradient failed")?;
    let eb = obj_b
        .energy_gradient(x, &mut gb)
        .context("spin-crossing state B energy/gradient failed")?;
    let (parallel, perpendicular, geff) = effective_gradient(ea, eb, &ga, &gb, 1.0, 140.0)?;
    Ok((ea, eb, ga, gb, parallel, perpendicular, geff))
}

pub fn optimize_spin_crossing<A: Objective, B: Objective>(
    mut x: Vec<f64>,
    connectivity: Option<ConnectivityModel>,
    obj_a: &mut A,
    obj_b: &mut B,
    opts: SpinCrossingOptions,
) -> Result<SpinCrossingResult> {
    if x.is_empty() || x.len() % 3 != 0 {
        bail!("spin-crossing optimizer requires a non-empty 3N coordinate vector");
    }
    if opts.coordinate_mode != SpinCrossingCoordinateMode::Cartesian && connectivity.is_none() {
        bail!("internal spin-crossing optimization requires a connectivity model");
    }
    let labels = opts.atom_labels.clone().unwrap_or_else(|| {
        (0..x.len() / 3)
            .map(|i| format!("X{}", i + 1))
            .collect::<Vec<_>>()
    });
    if labels.len() != x.len() / 3 {
        bail!("atom_labels length does not match coordinate vector");
    }

    let mut state = build_state(
        &x,
        opts.coordinate_mode,
        connectivity.as_ref(),
        None,
        opts.g_cutoff,
    )?;
    if let Some(s) = &state {
        if s.nint() == 0 {
            state = None;
        }
    }
    let actual_mode = if state.is_none() {
        SpinCrossingCoordinateMode::Cartesian
    } else {
        opts.coordinate_mode
    };
    let inv_dim = state.as_ref().map(|s| s.nint()).unwrap_or(x.len());
    let mut inv_hess = Array2::<f64>::eye(inv_dim) * opts.initial_inverse_hessian;
    let mut previous_coord = opt_coordinate(&state, &x);
    let mut previous_geff_opt = vec![0.0; inv_dim];
    let mut previous_x = x.clone();
    let mut previous_geff = vec![0.0; x.len()];
    let mut points = Vec::new();

    let mut traj = if let Some(path) = &opts.trajectory_file {
        Some(BufWriter::new(File::create(path)?))
    } else {
        None
    };

    if opts.verbosity > 0 {
        println!("Spin-crossing optimization");
        if let Some((m1, m2)) = opts.multiplicities {
            println!("Multiplicities: {m1}, {m2}");
        }
        println!("Update method: {:?}", opts.method);
        println!("Coordinates: {:?}", actual_mode);
        if let Some(s) = &state {
            println!(
                "Internal coordinates: {} bonds, {} angles, {} dihedrals",
                s.coordinates().nbonds(),
                s.coordinates().nangles(),
                s.coordinates().ndiheds()
            );
        }
        print_header(&opts);
    }

    for cycle in 1..=opts.max_cycles {
        let mut ga = vec![0.0; x.len()];
        let mut gb = vec![0.0; x.len()];
        let ea = obj_a
            .energy_gradient(&x, &mut ga)
            .context("spin-crossing state A energy/gradient failed")?;
        let eb = obj_b
            .energy_gradient(&x, &mut gb)
            .context("spin-crossing state B energy/gradient failed")?;
        let (parallel, perpendicular, geff) =
            effective_gradient(ea, eb, &ga, &gb, opts.fac_parallel, opts.fac_perpendicular)?;

        let current_state = if actual_mode == SpinCrossingCoordinateMode::Cartesian {
            None
        } else {
            build_state(
                &x,
                actual_mode,
                connectivity.as_ref(),
                state.as_ref().map(|s| s.coordinates().clone()),
                opts.g_cutoff,
            )?
        };
        let current_coord = opt_coordinate(&current_state, &x);
        let geff_opt = if let Some(s) = &current_state {
            internal_gradient(s, &geff)?
        } else {
            geff.clone()
        };

        let mut opt_step = if cycle == 1 {
            geff_opt
                .iter()
                .map(|g| -opts.initial_inverse_hessian * g)
                .collect::<Vec<_>>()
        } else {
            let (step, new_inv, _updated) = quasi_newton_step(
                opts.method,
                &previous_coord,
                &current_coord,
                &previous_geff_opt,
                &geff_opt,
                &inv_hess,
            );
            inv_hess = new_inv;
            step
        };

        let trial_x = if let Some(s) = &current_state {
            limit_internal_step(&mut opt_step, opts.max_step_internal);
            let mut scale = 1.0;
            let current_gap = (ea - eb).abs();
            let mut accepted = apply_internal_step(&x, s, &opt_step, &opts)?;
            while scale >= 1.0e-4 {
                let scaled: Vec<f64> = opt_step.iter().map(|v| scale * v).collect();
                let candidate = apply_internal_step(&x, s, &scaled, &opts)?;
                let ca = obj_a
                    .energy(&candidate)
                    .context("spin-crossing state A trial energy failed")?;
                let cb = obj_b
                    .energy(&candidate)
                    .context("spin-crossing state B trial energy failed")?;
                accepted = candidate;
                if (ca - cb).abs() <= current_gap {
                    break;
                }
                scale *= 0.5;
            }
            accepted
        } else {
            step_limit(&mut opt_step, opts.step_max_component);
            x.iter()
                .zip(opt_step.iter())
                .map(|(xi, si)| xi + si)
                .collect::<Vec<_>>()
        };
        let step: Vec<f64> = trial_x.iter().zip(x.iter()).map(|(a, b)| a - b).collect();

        if opts.verbosity > 0 {
            println!(
                "{:5} {:18.8} {:18.8} {:18.6} {:>11} {:>11} {:>11} {:>11}",
                cycle,
                ea,
                eb,
                (ea - eb) * HARTREE_TO_KJMOL,
                sci2(max_abs(&step)),
                sci2(rms(&step)),
                sci2(max_abs(&geff)),
                sci2(rms(&geff))
            );
        }
        if let Some(writer) = traj.as_mut() {
            write_frame(writer, &labels, &x, cycle, ea, eb, &geff)?;
        }

        points.push(SpinCrossingPoint {
            cycle,
            x: x.clone(),
            energy_a: ea,
            energy_b: eb,
            gradient_a: ga.clone(),
            gradient_b: gb.clone(),
            effective_gradient: geff.clone(),
            parallel_gradient: parallel,
            perpendicular_gradient: perpendicular,
            step: step.clone(),
        });

        if converged(ea - eb, &step, &geff, &opts) {
            let mut final_ga = vec![0.0; x.len()];
            let mut final_gb = vec![0.0; x.len()];
            let final_ea = obj_a
                .energy_gradient(&trial_x, &mut final_ga)
                .context("spin-crossing final state A energy/gradient failed")?;
            let final_eb = obj_b
                .energy_gradient(&trial_x, &mut final_gb)
                .context("spin-crossing final state B energy/gradient failed")?;
            let (_p, _q, final_geff) = effective_gradient(
                final_ea,
                final_eb,
                &final_ga,
                &final_gb,
                opts.fac_parallel,
                opts.fac_perpendicular,
            )?;
            if opts.verbosity > 0 {
                print_optimization_status(true, "Spin-crossing convergence reached.");
            }
            return Ok(SpinCrossingResult {
                x: trial_x,
                energy_a: final_ea,
                energy_b: final_eb,
                gradient_a: final_ga,
                gradient_b: final_gb,
                effective_gradient: final_geff,
                converged: true,
                cycles: cycle,
                method: opts.method,
                coordinates: actual_mode,
                multiplicities: opts.multiplicities,
                trajectory_file: opts.trajectory_file,
                points,
                message: "Spin-crossing convergence reached.".to_string(),
            });
        }

        previous_x = x;
        previous_geff = geff;
        previous_coord = current_coord;
        previous_geff_opt = geff_opt;
        x = trial_x;
        state = current_state;
    }

    let mut final_ga = vec![0.0; x.len()];
    let mut final_gb = vec![0.0; x.len()];
    let final_ea = obj_a
        .energy_gradient(&x, &mut final_ga)
        .context("spin-crossing final state A energy/gradient failed")?;
    let final_eb = obj_b
        .energy_gradient(&x, &mut final_gb)
        .context("spin-crossing final state B energy/gradient failed")?;
    let (_p, _q, final_geff) = effective_gradient(
        final_ea,
        final_eb,
        &final_ga,
        &final_gb,
        opts.fac_parallel,
        opts.fac_perpendicular,
    )?;
    if opts.verbosity > 0 {
        print_optimization_status(
            false,
            "Maximum spin-crossing optimization steps reached without convergence.",
        );
    }
    Ok(SpinCrossingResult {
        x,
        energy_a: final_ea,
        energy_b: final_eb,
        gradient_a: final_ga,
        gradient_b: final_gb,
        effective_gradient: final_geff,
        converged: false,
        cycles: opts.max_cycles,
        method: opts.method,
        coordinates: actual_mode,
        multiplicities: opts.multiplicities,
        trajectory_file: opts.trajectory_file,
        points,
        message: "Maximum spin-crossing optimization steps reached without convergence."
            .to_string(),
    })
}
