//! Internal-coordinate BFGS geometry optimization.
//!
//! - build a fixed internal-coordinate set from the initial geometry
//! - transform Cartesian gradients into internal gradients
//! - take a bounded internal step `dq`
//! - map `dq -> dx` with a small best-fit iteration
//! - update a Hessian approximation with the BFGS Hessian update

use anyhow::{bail, Result};
use ndarray::{Array1, Array2};

use super::internal_coords::{
    analyze_structure, best_fit_dq_to_cart, bpg_matrix, compute_internals, diff_internals,
    ConnectivityModel, InternalCoords,
};
use super::traits::Objective;

use crate::optimizer::linalg_solve::solve_linear_system_array2;

const HARTREE_TO_KJMOL: f64 = 2625.499638;

#[derive(Debug, Clone)]
pub struct GeomOptOptions {
    pub max_cycles: usize,
    pub max_step_internal: f64,
    pub energy_tol: f64,
    pub max_step: f64,
    pub rms_step: f64,
    pub max_gradient: f64,
    pub g_rms_tol: f64,
    pub best_fit_iters: usize,
    pub verbosity: u8,
}

impl Default for GeomOptOptions {
    fn default() -> Self {
        Self {
            max_cycles: 300,
            max_step_internal: 0.2,
            energy_tol: 5.0e-6,
            max_step: 4.0e-3,
            rms_step: 2.0e-3,
            max_gradient: 3.0e-4,
            g_rms_tol: 1.0e-4,
            best_fit_iters: 5,
            verbosity: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GeomOptResult {
    pub x: Vec<f64>,
    pub energy: f64,
    /// Cartesian nuclear gradient at `x`, in Hartree/bohr.
    pub gradient: Vec<f64>,
    pub grad_rms: f64,
    pub cycles: usize,
    pub converged: bool,
    pub trajectory: Vec<Vec<f64>>,
    pub trajectory_energies: Vec<f64>,
}

fn vec_dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn vec_norm2(a: &[f64]) -> f64 {
    vec_dot(a, a)
}

fn vec_norm(a: &[f64]) -> f64 {
    vec_norm2(a).sqrt()
}

fn vec_rms(a: &[f64]) -> f64 {
    if a.is_empty() {
        0.0
    } else {
        vec_norm(a) / (a.len() as f64).sqrt()
    }
}

fn vec_max_abs(a: &[f64]) -> f64 {
    a.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

fn print_header(opts: &GeomOptOptions) {
    println!(" Energy Change:   {} Eh", sci2(opts.energy_tol));
    println!(" Max. Gradient:   {} Eh/bohr", sci2(opts.max_gradient));
    println!(" RMS Gradient:    {} Eh/bohr", sci2(opts.g_rms_tol));
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
        opts.g_rms_tol,
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

fn internal_gradient(
    ic: &InternalCoords,
    b: &[Vec<f64>],
    inv_g: &[Vec<f64>],
    grad_x: &[f64],
) -> Vec<f64> {
    // gq = invG * B * grad
    let nint = ic.nint();
    let nx = grad_x.len();

    let mut bg = vec![0.0f64; nint];
    for i in 0..nint {
        let mut acc = 0.0;
        for k in 0..nx {
            acc += b[i][k] * grad_x[k];
        }
        bg[i] = acc;
    }

    let mut gq = vec![0.0f64; nint];
    for i in 0..nint {
        let mut acc = 0.0;
        for j in 0..nint {
            acc += inv_g[i][j] * bg[j];
        }
        gq[i] = acc;
    }
    gq
}

fn initial_hq(ic: &InternalCoords) -> Array2<f64> {
    let nb = ic.nbonds();
    let na = ic.nangles();
    let nd = ic.ndiheds();
    let nint = nb + na + nd;
    let mut hq = Array2::<f64>::zeros((nint, nint));
    for i in 0..nb {
        hq[(i, i)] = 1.0;
    }
    for i in 0..na {
        hq[(nb + i, nb + i)] = 0.25;
    }
    for i in 0..nd {
        hq[(nb + na + i, nb + na + i)] = 0.125;
    }
    hq
}

/// Run internal-coordinate BFGS geometry optimizer.
///
/// Inputs:
/// - `x0`: initial geometry (bohr), length `3*nat`.
/// - `conn`: connectivity model (for internal coordinate discovery).
/// - `obj`: any energy/gradient provider.
pub fn geom_opt_internal_bfgs<O: Objective>(
    x: Vec<f64>,
    conn: ConnectivityModel,
    obj: &mut O,
    opts: GeomOptOptions,
) -> Result<GeomOptResult> {
    let nat = conn.rcov_disp.len();
    if x.len() != 3 * nat {
        bail!(
            "geom_opt_internal_bfgs: expected x length {}, got {}",
            3 * nat,
            x.len()
        );
    }

    let ic = analyze_structure(&x, &conn)?;
    geom_opt_internal_bfgs_with_internals(x, ic, obj, opts)
}

/// Run the internal-coordinate BFGS optimizer with a pre-defined internal coordinate set.
pub fn geom_opt_internal_bfgs_with_internals<O: Objective>(
    mut x: Vec<f64>,
    ic: InternalCoords,
    obj: &mut O,
    opts: GeomOptOptions,
) -> Result<GeomOptResult> {
    let nat = ic.nat;
    if x.len() != 3 * nat {
        bail!(
            "geom_opt_internal_bfgs_with_internals: expected x length {}, got {}",
            3 * nat,
            x.len()
        );
    }
    let nint = ic.nint();
    if nint == 0 {
        bail!("geom_opt_internal_bfgs_with_internals: no internal coordinates found");
    }

    let mut qs = compute_internals(&x, &ic);
    let mut bpg = bpg_matrix(&x, &ic)?;

    let mut grad_x = vec![0.0f64; x.len()];
    let mut energy = obj.energy_gradient(&x, &mut grad_x)?;
    let mut gq = internal_gradient(&ic, &bpg.b, &bpg.inv_g, &grad_x);

    let mut hq = initial_hq(&ic);
    let mut converged = false;
    let mut cycles_done = 0usize;
    let mut trajectory = Vec::new();
    let mut trajectory_energies = Vec::new();

    if opts.verbosity > 0 {
        print_header(&opts);
    }

    for cyc in 1..=opts.max_cycles {
        cycles_done = cyc;

        // dq = -Hq^{-1} gq  => solve(Hq, -gq)
        let rhs = Array1::from(gq.iter().map(|&v| -v).collect::<Vec<_>>());
        let mut dq = solve_linear_system_array2(&hq, &rhs)?.to_vec();

        // Step limiting ( maxstep = 0.2 on internal coords).
        let stp2 = vec_norm2(&dq);
        let stp = stp2.sqrt();
        if stp > opts.max_step_internal && stp > 0.0 {
            let s = opts.max_step_internal / stp;
            for v in &mut dq {
                *v *= s;
            }
        }

        let x_old = x.clone();
        let grad_x_old = grad_x.clone();
        let gq_old = gq.clone();
        let e_old = energy;

        best_fit_dq_to_cart(&mut x, &qs, &dq, &ic, &bpg, opts.best_fit_iters)?;
        let cart_step: Vec<f64> = x.iter().zip(x_old.iter()).map(|(a, b)| a - b).collect();

        energy = obj.energy_gradient(&x, &mut grad_x)?;

        qs = compute_internals(&x, &ic);
        bpg = bpg_matrix(&x, &ic)?;
        gq = internal_gradient(&ic, &bpg.b, &bpg.inv_g, &grad_x);

        let energy_change = energy - e_old;
        let gnorm = vec_norm(&grad_x);
        let grms = vec_rms(&grad_x);

        if opts.verbosity > 0 {
            println!(
                "{:5} {:18.8} {:18.6} {:>11} {:>11} {:>11} {:>11}",
                cyc,
                energy,
                energy_change * HARTREE_TO_KJMOL,
                sci2(vec_max_abs(&cart_step)),
                sci2(vec_rms(&cart_step)),
                sci2(vec_max_abs(&grad_x)),
                sci2(grms)
            );
        }

        trajectory.push(x.clone());
        trajectory_energies.push(energy);

        if vec_max_abs(&grad_x) <= opts.max_gradient
            && grms <= opts.g_rms_tol
            && vec_max_abs(&cart_step) <= opts.max_step
            && vec_rms(&cart_step) <= opts.rms_step
            && energy_change.abs() <= opts.energy_tol
        {
            if opts.verbosity > 0 {
                print_optimization_status(true, "Convergence reached.");
            }
            converged = true;
            break;
        }

        // BFGS Hessian update:
        // yk = gq - gqold
        // fac1 = dq·yk
        // vec = Hq*dq
        // fac2 = dq·vec
        // Hq += yk yk^T / fac1 - vec vec^T / fac2
        let yk = diff_internals(&gq, &gq_old, ic.ndiheds());
        let fac1 = vec_dot(&dq, &yk);
        let vec_h = hq.dot(&Array1::from(dq.clone())).to_vec();
        let fac2 = vec_dot(&dq, &vec_h);

        // Guard: if update would blow up, skip and keep Hq.
        if fac1.abs() > 1e-14 && fac2.abs() > 1e-14 {
            for i in 0..nint {
                for j in 0..nint {
                    hq[(i, j)] += (yk[i] * yk[j]) / fac1 - (vec_h[i] * vec_h[j]) / fac2;
                }
            }
        } else {
            // If this happens, also revert the last step (safer) and continue.
            x = x_old;
            energy = e_old;
            grad_x = grad_x_old;
            qs = compute_internals(&x, &ic);
            bpg = bpg_matrix(&x, &ic)?;
            gq = gq_old;
        }
    }

    if opts.verbosity > 0 && !converged {
        print_optimization_status(
            false,
            "Maximum optimization steps reached without convergence.",
        );
    }

    Ok(GeomOptResult {
        x,
        energy,
        gradient: grad_x.clone(),
        grad_rms: (vec_norm2(&grad_x) / (grad_x.len() as f64)).sqrt(),
        cycles: cycles_done,
        converged,
        trajectory,
        trajectory_energies,
    })
}
