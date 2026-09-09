//! Smite-style Cartesian minimum-search optimizer.
//!
//! This is intentionally separate from the legacy `geom_opt` module so both
//! optimizers can be compared side by side.

use anyhow::{bail, Result};
use ndarray::{Array1, Array2};

use super::traits::Objective;

const HARTREE_TO_KJMOL: f64 = 2625.499638;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CartesianMinimumMethod {
    Bfgs,
    SteepestDescent,
    TwoPointGradient1,
    TwoPointGradient2,
}

#[derive(Debug, Clone)]
pub struct CartesianMinimumOptions {
    pub method: CartesianMinimumMethod,
    pub max_cycles: usize,
    pub energy_tol: f64,
    pub max_step: f64,
    pub rms_step: f64,
    pub max_gradient: f64,
    pub rms_gradient: f64,
    pub step_max_component: f64,
    pub min_line_search_scale: f64,
    pub verbosity: u8,
}

impl Default for CartesianMinimumOptions {
    fn default() -> Self {
        Self {
            method: CartesianMinimumMethod::Bfgs,
            max_cycles: 100,
            energy_tol: 5.0e-6,
            max_step: 4.0e-3,
            rms_step: 2.0e-3,
            max_gradient: 3.0e-4,
            rms_gradient: 1.0e-4,
            step_max_component: 0.1,
            min_line_search_scale: 1.0e-4,
            verbosity: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct CartesianMinimumResult {
    pub x: Vec<f64>,
    pub energy: f64,
    pub gradient: Vec<f64>,
    pub converged: bool,
    pub cycles: usize,
    pub message: String,
    pub trajectory: Vec<Vec<f64>>,
    pub trajectory_energies: Vec<f64>,
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

fn print_header(opts: &CartesianMinimumOptions) {
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

pub fn step_limit(step: &mut [f64], max_component: f64) {
    if step.is_empty() || max_component <= 0.0 {
        return;
    }
    let stpmax = (step.len() as f64) * max_component;
    let nrm = norm(step);
    if nrm > stpmax && nrm > 0.0 {
        let scale = stpmax / nrm;
        for v in step.iter_mut() {
            *v *= scale;
        }
    }
    let largest = max_abs(step);
    if largest > max_component && largest > 0.0 {
        let scale = max_component / largest;
        for v in step.iter_mut() {
            *v *= scale;
        }
    }
}

fn converged(
    energy_change: f64,
    step: &[f64],
    grad: &[f64],
    opts: &CartesianMinimumOptions,
) -> bool {
    max_abs(grad) <= opts.max_gradient
        && rms(grad) <= opts.rms_gradient
        && max_abs(step) <= opts.max_step
        && rms(step) <= opts.rms_step
        && energy_change.abs() <= opts.energy_tol
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
    let mut left = Array2::<f64>::eye(n);
    let mut right = Array2::<f64>::eye(n);
    for i in 0..n {
        for j in 0..n {
            left[(i, j)] -= rho * step[i] * y[j];
            right[(i, j)] -= rho * y[i] * step[j];
        }
    }
    let updated = left.dot(&old).dot(&right)
        + rho
            * Array1::from(step.to_vec())
                .insert_axis(ndarray::Axis(1))
                .dot(&Array1::from(step.to_vec()).insert_axis(ndarray::Axis(0)));
    if updated.iter().all(|v| v.is_finite()) {
        *h_inv = updated;
    }
}

fn line_search_cartesian<O: Objective>(
    obj: &mut O,
    x: &[f64],
    energy: f64,
    grad: &[f64],
    step: &[f64],
    min_scale: f64,
) -> Result<(Vec<f64>, f64, Vec<f64>, Vec<f64>, f64)> {
    let mut scale = 1.0;
    let directional = dot(grad, step);
    while scale >= min_scale {
        let trial_x: Vec<f64> = x
            .iter()
            .zip(step.iter())
            .map(|(xi, si)| xi + scale * si)
            .collect::<Vec<_>>();
        let mut trial_grad = vec![0.0; x.len()];
        let trial_energy = obj.energy_gradient(&trial_x, &mut trial_grad)?;
        if trial_energy <= energy + 1.0e-4 * scale * directional || trial_energy < energy {
            let accepted: Vec<f64> = step.iter().map(|v| scale * v).collect::<Vec<_>>();
            return Ok((trial_x, trial_energy, trial_grad, accepted, scale));
        }
        scale *= 0.5;
    }

    let trial_x: Vec<f64> = x
        .iter()
        .zip(step.iter())
        .map(|(xi, si)| xi + scale * si)
        .collect::<Vec<_>>();
    let mut trial_grad = vec![0.0; x.len()];
    let trial_energy = obj.energy_gradient(&trial_x, &mut trial_grad)?;
    let accepted: Vec<f64> = step.iter().map(|v| scale * v).collect::<Vec<_>>();
    Ok((trial_x, trial_energy, trial_grad, accepted, scale))
}

pub fn minimize_cartesian<O: Objective>(
    mut x: Vec<f64>,
    obj: &mut O,
    opts: CartesianMinimumOptions,
) -> Result<CartesianMinimumResult> {
    if x.is_empty() || x.len() % 3 != 0 {
        bail!("Cartesian optimizer requires a non-empty 3N coordinate vector");
    }
    let ndim = x.len();
    let mut h_inv = Array2::<f64>::eye(ndim);
    let mut grad = vec![0.0; ndim];
    let mut energy = obj.energy_gradient(&x, &mut grad)?;
    let mut previous_x: Option<Vec<f64>> = None;
    let mut previous_grad: Option<Vec<f64>> = None;
    let mut trajectory = vec![x.clone()];
    let mut trajectory_energies = vec![energy];

    if opts.verbosity > 0 {
        print_header(&opts);
    }

    for cycle in 1..=opts.max_cycles {
        let mut step: Vec<f64> = match opts.method {
            CartesianMinimumMethod::SteepestDescent => {
                grad.iter().map(|g| -0.7 * g).collect::<Vec<_>>()
            }
            CartesianMinimumMethod::TwoPointGradient1
            | CartesianMinimumMethod::TwoPointGradient2 => {
                if let (Some(px), Some(pg)) = (&previous_x, &previous_grad) {
                    let dg: Vec<f64> = grad
                        .iter()
                        .zip(pg.iter())
                        .map(|(g, p)| g - p)
                        .collect::<Vec<_>>();
                    let dx: Vec<f64> = x
                        .iter()
                        .zip(px.iter())
                        .map(|(xi, pi)| xi - pi)
                        .collect::<Vec<_>>();
                    let denom1 = dot(&dg, &dg);
                    let denom2 = dot(&dg, &dx);
                    let alpha = match opts.method {
                        CartesianMinimumMethod::TwoPointGradient1 => {
                            if denom1.abs() > 1.0e-16 {
                                denom2 / denom1
                            } else {
                                0.7
                            }
                        }
                        _ => {
                            if denom2.abs() > 1.0e-16 {
                                dot(&dx, &dx) / denom2
                            } else {
                                0.7
                            }
                        }
                    };
                    grad.iter().map(|g| -alpha * g).collect::<Vec<_>>()
                } else {
                    grad.iter().map(|g| -0.7 * g).collect::<Vec<_>>()
                }
            }
            CartesianMinimumMethod::Bfgs => mat_vec(&h_inv, &grad)
                .into_iter()
                .map(|v| -v)
                .collect::<Vec<_>>(),
        };
        step_limit(&mut step, opts.step_max_component);

        let old_x = x.clone();
        let old_grad = grad.clone();
        let old_energy = energy;
        let (trial_x, trial_energy, trial_grad, accepted_step, _scale) =
            line_search_cartesian(obj, &x, energy, &grad, &step, opts.min_line_search_scale)?;
        let energy_change = trial_energy - energy;

        if opts.verbosity > 0 {
            println!(
                "{:5} {:18.8} {:18.6} {:>11} {:>11} {:>11} {:>11}",
                cycle,
                trial_energy,
                energy_change * HARTREE_TO_KJMOL,
                sci2(max_abs(&accepted_step)),
                sci2(rms(&accepted_step)),
                sci2(max_abs(&trial_grad)),
                sci2(rms(&trial_grad))
            );
        }

        if opts.method == CartesianMinimumMethod::Bfgs {
            let y: Vec<f64> = trial_grad
                .iter()
                .zip(old_grad.iter())
                .map(|(g, og)| g - og)
                .collect::<Vec<_>>();
            bfgs_inverse_update(&mut h_inv, &accepted_step, &y);
        }

        previous_x = Some(old_x);
        previous_grad = Some(old_grad);
        x = trial_x;
        energy = trial_energy;
        grad = trial_grad;
        trajectory.push(x.clone());
        trajectory_energies.push(energy);

        if converged(energy_change, &accepted_step, &grad, &opts) {
            if opts.verbosity > 0 {
                print_optimization_status(true, "Convergence reached.");
            }
            return Ok(CartesianMinimumResult {
                x,
                energy,
                gradient: grad,
                converged: true,
                cycles: cycle,
                message: "Convergence reached.".to_string(),
                trajectory,
                trajectory_energies,
            });
        }

        if opts.min_line_search_scale > 0.0 && energy > old_energy {
            h_inv = Array2::<f64>::eye(ndim);
        }
    }

    if opts.verbosity > 0 {
        print_optimization_status(
            false,
            "Maximum optimization steps reached without convergence.",
        );
    }

    Ok(CartesianMinimumResult {
        x,
        energy,
        gradient: grad,
        converged: false,
        cycles: opts.max_cycles,
        message: "Maximum optimization steps reached without convergence.".to_string(),
        trajectory,
        trajectory_energies,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    struct HarmonicCartesian {
        target: Vec<f64>,
    }

    impl Objective for HarmonicCartesian {
        fn energy(&mut self, x: &[f64]) -> Result<f64> {
            Ok(0.5
                * x.iter()
                    .zip(self.target.iter())
                    .map(|(xi, ti)| (xi - ti) * (xi - ti))
                    .sum::<f64>())
        }

        fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
            for ((g, xi), ti) in grad.iter_mut().zip(x.iter()).zip(self.target.iter()) {
                *g = xi - ti;
            }
            Ok(())
        }
    }

    #[test]
    fn cartesian_bfgs_minimizes_harmonic_objective() {
        let mut obj = HarmonicCartesian {
            target: vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        };
        let opts = CartesianMinimumOptions {
            verbosity: 0,
            max_gradient: 1.0e-8,
            rms_gradient: 1.0e-8,
            max_step: 1.0e-8,
            rms_step: 1.0e-8,
            energy_tol: 1.0e-12,
            ..Default::default()
        };
        let res = minimize_cartesian(vec![0.3, 0.0, 0.0, 1.4, 0.0, 0.0], &mut obj, opts).unwrap();
        assert!(res.converged);
        assert!(res.energy < 1.0e-14);
    }
}
