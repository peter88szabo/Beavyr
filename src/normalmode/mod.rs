#![allow(dead_code, unused_imports)]
//! Normal mode analysis (mass-weighted Hessian).
//!
//! This module implements a general normal mode analysis that is not
//! tied to a specific Hamiltonian. Eckart projection can be enabled
//! to remove translations/rotations (or also the reaction-path mode).

use anyhow::{anyhow, bail, Result};
use ndarray::Array2;

use crate::normalmode::bmat_smat::{eckart_reactionpath_transform, eckart_transform};
use crate::normalmode::linalg_shim::{Backend, LinAlg};

pub mod bmat_smat;
pub mod fd_hessian;
pub mod inertia;
pub mod ir_intensity;
// Imported alongside the analysis code so it needs no BLAS: see
// `linalg_shim`'s own docs.
pub mod jacobi_diag;
pub mod linalg_shim;
mod print;
pub mod print_matrix_in_ao;
pub mod thermofuncs;

pub use print::{coord_labels, format_hessian_title, print_hessian_matrix};

const C1: f64 = 1.0 / 0.529_177_210_903; // Angstrom -> bohr (CODATA 2018 Bohr radius)
const C2: f64 = 1.0 / 627.51; // kcal/mol -> Hartree (unused but kept for parity)
const C3: f64 = 1_822.888_486_209; // g/mol -> electron mass unit (amu to au)
const C4: f64 = 27.2114; // Hartree -> eV (unused but kept for parity)
const C5: f64 = 219_474.0; // Hartree -> cm^-1 (unused but kept for parity)
const C6: f64 = 41.341_105; // fs -> time in au (unused but kept for parity)
const C9: f64 = 1.0e8 * C1; // frequency in cm^-1 -> bohr^-1
const C10: f64 = 137.035_999_084; // speed of light in atomic units (CODATA 2018 inverse fine-structure constant; was CODATA 2010's 137.035999074)

fn eigval_to_freq(lambd: f64) -> f64 {
    if lambd < 0.0 {
        -lambd.abs().sqrt()
    } else {
        lambd.sqrt()
    }
}

fn freq_to_cm1(omega: f64) -> f64 {
    omega / C10 * C9 / (std::f64::consts::PI * 2.0)
}

#[derive(Debug, Clone)]
pub struct NormalModeResult {
    pub frequencies_au: Vec<f64>,
    pub low_frequencies_au: Vec<f64>,
    pub modes: Array2<f64>, // columns: normal modes (Cartesian, not mass-weighted)
    pub negative_indices: Vec<usize>,
    pub zero_indices: Vec<usize>,
    pub positive_indices: Vec<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EckartMode {
    Off,
    VibRot,
    ReactionPath,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HessianWeight {
    Cartesian,
    MassWeighted,
}

fn build_mass_weight_matrix(mass_au: &[f64]) -> Result<Array2<f64>> {
    if mass_au.is_empty() {
        bail!("mass array is empty");
    }
    let ncoord = 3 * mass_au.len();
    let mut m = Array2::<f64>::zeros((ncoord, ncoord));
    for (i, &mi) in mass_au.iter().enumerate() {
        if mi <= 0.0 {
            bail!("mass must be positive; got {mi}");
        }
        let inv = 1.0 / mi.sqrt();
        for k in 0..3 {
            let idx = 3 * i + k;
            m[(idx, idx)] = inv;
        }
    }
    Ok(m)
}

/// Compute normal modes and frequencies from a Cartesian Hessian.
///
/// Inputs:
/// - `mass`: atomic masses (amu or any consistent unit, used as weights).
/// - `hessian`: Cartesian Hessian matrix (3N x 3N).
/// - `linear`: whether the molecule is linear (affects 5 vs 6 low modes).
pub fn normal_modes(
    mass: &[f64],
    hessian: &Array2<f64>,
    linear: bool,
    coords_bohr: &[f64],
) -> Result<NormalModeResult> {
    normal_modes_with_projection(
        mass,
        hessian,
        linear,
        Some(coords_bohr),
        EckartMode::VibRot,
        None,
    )
}

/// Compute normal modes and frequencies from a Cartesian Hessian with
/// optional Eckart projection.
///
/// Inputs:
/// - `mass`: atomic masses in amu.
/// - `hessian`: Cartesian Hessian matrix (3N x 3N) in Eh/bohr^2.
/// - `linear`: whether the molecule is linear (affects 5 vs 6 low modes).
/// - `coords_bohr`: Cartesian coordinates (bohr) for Eckart projection.
/// - `eckart`: projection mode (off, vib/rot, or reaction-path).
/// - `grad_cart`: Cartesian gradient (Eh/bohr) required for reaction-path mode.
pub fn normal_modes_with_projection(
    mass: &[f64],
    hessian: &Array2<f64>,
    linear: bool,
    coords_bohr: Option<&[f64]>,
    eckart: EckartMode,
    grad_cart: Option<&[f64]>,
) -> Result<NormalModeResult> {
    let nat = mass.len();
    let ncoord = 3 * nat;
    if ncoord == 0 {
        bail!("empty system");
    }
    let (nr, nc) = hessian.dim();
    if nr != ncoord || nc != ncoord {
        bail!("hessian must be 3N x 3N (N = number of atoms); got {nr} x {nc}");
    }

    let mass_au: Vec<f64> = mass.iter().map(|&m| m * C3).collect();

    let la = LinAlg::new(Backend::Auto)?;
    let m = build_mass_weight_matrix(&mass_au)?;
    let mut h_mw = la.matmul(&la.matmul(&m, hessian), &m);

    if eckart != EckartMode::Off {
        let coords = coords_bohr
            .ok_or_else(|| anyhow::anyhow!("coords_bohr is required for Eckart projection"))?;
        if coords.len() != ncoord {
            bail!("coords_bohr must have length 3N");
        }
        println!();
        println!("!!!!!              Eckart transformation switched on            !!!!!");
        println!("!!!!! Translations and Rotations projected out from the Hessian !!!!!");
        println!();
        match eckart {
            EckartMode::VibRot => {
                h_mw = eckart_transform(&mass_au, coords, &h_mw, &la)?;
            }
            EckartMode::ReactionPath => {
                let grad = grad_cart
                    .ok_or_else(|| anyhow!("grad_cart is required for reaction-path projection"))?;
                if grad.len() != ncoord {
                    bail!("grad_cart must have length 3N");
                }
                let mut grad_mw = vec![0.0_f64; ncoord];
                for at in 0..nat {
                    let inv = 1.0 / mass_au[at].sqrt();
                    grad_mw[3 * at] = grad[3 * at] * inv;
                    grad_mw[3 * at + 1] = grad[3 * at + 1] * inv;
                    grad_mw[3 * at + 2] = grad[3 * at + 2] * inv;
                }
                h_mw = eckart_reactionpath_transform(&mass_au, coords, &grad_mw, &h_mw, &la)?;
            }
            EckartMode::Off => {}
        }
    }

    let (evals, evecs) = la.eigh(&h_mw)?;

    let mut frequencies_au = Vec::with_capacity(ncoord);
    for i in 0..ncoord {
        frequencies_au.push(eigval_to_freq(evals[i]));
    }

    // Un-mass-weight the eigenvectors.
    let mut modes = evecs;
    for at in 0..nat {
        let inv = 1.0 / mass_au[at].sqrt();
        for k in 0..3 {
            let row = 3 * at + k;
            for col in 0..ncoord {
                modes[(row, col)] *= inv;
            }
        }
    }

    let nlow = 6 - usize::from(linear);
    let tolerance = (0.01 * C10 / C9 * (std::f64::consts::PI * 2.0)).powi(2);
    let mut negative_indices = Vec::new();
    let mut zero_indices = Vec::new();
    let mut positive_indices = Vec::new();

    for (idx, &eigenvalue) in evals.iter().enumerate() {
        if eigenvalue < 0.0 && eigenvalue.abs() > tolerance {
            negative_indices.push(idx);
        } else if eigenvalue.abs() <= tolerance {
            zero_indices.push(idx);
        } else {
            positive_indices.push(idx);
        }
    }

    // DEVIATION from the imported Behemoth code, deliberate: the absolute
    // tolerance above corresponds to about 0.01 cm^-1, which holds for small
    // molecules but not larger ones. On the 19-atom transition state in
    // examples/, the six projected modes come out at 0.02-0.04 cm^-1 --
    // ordinary double-precision noise from a 57x57 diagonalisation -- and were
    // counted as real vibrations, putting six spurious near-zero modes into
    // the mode list and, worse, into the vibrational partition function.
    //
    // When a projection was actually applied we know exactly how many modes it
    // removed, so identify them by rank rather than by an absolute threshold:
    // the projected ones are always the smallest in magnitude, by orders of
    // magnitude, so this cannot swallow a genuine low-frequency torsion.
    let nprojected = match eckart {
        EckartMode::Off => 0,
        EckartMode::VibRot => nlow,
        // The reaction-path projection also removes the gradient direction.
        EckartMode::ReactionPath => nlow + 1,
    };
    if nprojected > 0 && nprojected <= ncoord && zero_indices.len() != nprojected {
        let mut by_magnitude: Vec<usize> = (0..ncoord).collect();
        by_magnitude.sort_by(|&a, &b| evals[a].abs().total_cmp(&evals[b].abs()));
        let projected = &by_magnitude[..nprojected];
        negative_indices.clear();
        zero_indices.clear();
        positive_indices.clear();
        for (idx, &eigenvalue) in evals.iter().enumerate() {
            if projected.contains(&idx) {
                zero_indices.push(idx);
            } else if eigenvalue < 0.0 {
                negative_indices.push(idx);
            } else {
                positive_indices.push(idx);
            }
        }
    }

    let mut low_frequencies_au = Vec::new();
    if zero_indices.len() == nlow {
        for i in 0..ncoord {
            if zero_indices.contains(&i) {
                low_frequencies_au.push(frequencies_au[i]);
            }
        }
    } else {
        for i in 0..nlow.min(frequencies_au.len()) {
            low_frequencies_au.push(frequencies_au[i]);
        }
    }

    Ok(NormalModeResult {
        frequencies_au,
        low_frequencies_au,
        modes,
        negative_indices,
        zero_indices,
        positive_indices,
    })
}

pub fn project_hessian(
    mass_amu: &[f64],
    coords_bohr: &[f64],
    hessian: &Array2<f64>,
    eckart: EckartMode,
    grad_cart: Option<&[f64]>,
    weight: HessianWeight,
) -> Result<Array2<f64>> {
    let nat = mass_amu.len();
    let ncoord = 3 * nat;
    if coords_bohr.len() != ncoord {
        bail!("coords_bohr must have length 3N");
    }
    let (nr, nc) = hessian.dim();
    if nr != ncoord || nc != ncoord {
        bail!("hessian must be 3N x 3N (N = number of atoms); got {nr} x {nc}");
    }

    let mass_au: Vec<f64> = mass_amu.iter().map(|&m| m * C3).collect();
    let la = LinAlg::new(Backend::Auto)?;
    let m = build_mass_weight_matrix(&mass_au)?;
    let mut h_mw = la.matmul(&la.matmul(&m, hessian), &m);

    if eckart != EckartMode::Off {
        match eckart {
            EckartMode::VibRot => {
                h_mw = eckart_transform(&mass_au, coords_bohr, &h_mw, &la)?;
            }
            EckartMode::ReactionPath => {
                let grad = grad_cart
                    .ok_or_else(|| anyhow!("grad_cart is required for reaction-path projection"))?;
                if grad.len() != ncoord {
                    bail!("grad_cart must have length 3N");
                }
                let mut grad_mw = vec![0.0_f64; ncoord];
                for at in 0..nat {
                    let inv = 1.0 / mass_au[at].sqrt();
                    grad_mw[3 * at] = grad[3 * at] * inv;
                    grad_mw[3 * at + 1] = grad[3 * at + 1] * inv;
                    grad_mw[3 * at + 2] = grad[3 * at + 2] * inv;
                }
                h_mw = eckart_reactionpath_transform(&mass_au, coords_bohr, &grad_mw, &h_mw, &la)?;
            }
            EckartMode::Off => {}
        }
    }

    let out = match weight {
        HessianWeight::MassWeighted => h_mw,
        HessianWeight::Cartesian => {
            let mut weights = vec![0.0_f64; ncoord];
            for at in 0..nat {
                let w = mass_au[at].sqrt();
                weights[3 * at] = w;
                weights[3 * at + 1] = w;
                weights[3 * at + 2] = w;
            }
            let mut out = Array2::<f64>::zeros((ncoord, ncoord));
            for i in 0..ncoord {
                for j in 0..ncoord {
                    out[(i, j)] = h_mw[(i, j)] * weights[i] * weights[j];
                }
            }
            out
        }
    };
    Ok(out)
}

pub fn project_hessian_vibrot(
    mass_amu: &[f64],
    coords_bohr: &[f64],
    hessian: &Array2<f64>,
) -> Result<Array2<f64>> {
    let nat = mass_amu.len();
    let ncoord = 3 * nat;
    if coords_bohr.len() != ncoord {
        bail!("coords_bohr must have length 3N");
    }
    let (nr, nc) = hessian.dim();
    if nr != ncoord || nc != ncoord {
        bail!("hessian must be 3N x 3N (N = number of atoms); got {nr} x {nc}");
    }

    let mass_au: Vec<f64> = mass_amu.iter().map(|&m| m * C3).collect();
    let la = LinAlg::new(Backend::Auto)?;
    let m = build_mass_weight_matrix(&mass_au)?;
    let h_mw = la.matmul(&la.matmul(&m, hessian), &m);
    let h_mw_proj = eckart_transform(&mass_au, coords_bohr, &h_mw, &la)?;

    let mut weights = vec![0.0_f64; ncoord];
    for at in 0..nat {
        let w = mass_au[at].sqrt();
        weights[3 * at] = w;
        weights[3 * at + 1] = w;
        weights[3 * at + 2] = w;
    }

    let mut out = Array2::<f64>::zeros((ncoord, ncoord));
    for i in 0..ncoord {
        for j in 0..ncoord {
            out[(i, j)] = h_mw_proj[(i, j)] * weights[i] * weights[j];
        }
    }
    Ok(out)
}

/// Print frequencies to stdout and to `vibrational_freq_{name}.dat`.
pub fn print_frequencies(
    name: &str,
    frequencies_au: &[f64],
    linear: bool,
    reduced_masses: Option<&[f64]>,
    ir_intensities: Option<&[f64]>,
) -> Result<()> {
    let nlow = 6 - usize::from(linear);
    if frequencies_au.len() < nlow {
        bail!("not enough frequencies to print");
    }
    if let Some(mu) = reduced_masses {
        if mu.len() != frequencies_au.len() {
            bail!("reduced_masses length must match frequencies");
        }
    }
    if let Some(ir) = ir_intensities {
        if ir.len() != frequencies_au.len() {
            bail!("ir_intensities length must match frequencies");
        }
    }

    let fmt_opt = |v: Option<f64>, width: usize, prec: usize| -> String {
        match v {
            Some(val) => format!("{:>width$.prec$}", val, width = width, prec = prec),
            None => format!("{:>width$}", "n/a", width = width),
        }
    };

    println!();
    println!("---------- Low Frequencies ----------");
    println!(
        "{:6} {:6} {:12} {:12} {:12}",
        "index1", "index2", "freq[cm-1]", "mu_red[amu]", "IR[km/mol]"
    );
    for i in 0..nlow {
        let omega = frequencies_au[i];
        let freq = freq_to_cm1(omega);
        let mu = reduced_masses.and_then(|v| v.get(i)).copied();
        let mu_str = fmt_opt(mu, 12, 4);
        let ir_str = fmt_opt(None, 12, 4);
        if omega < 0.0 && omega.abs() > 0.01 * C10 / C9 * (std::f64::consts::PI * 2.0) {
            println!(
                "{:6} {:6} {:12.2} {mu_str} {ir_str} {:>10}",
                i as isize,
                i as isize - nlow as isize,
                freq,
                "<-- Imag"
            );
        } else {
            println!(
                "{:6} {:6} {:12.2} {mu_str} {ir_str}",
                i as isize,
                i as isize - nlow as isize,
                freq
            );
        }
    }
    println!("---------- High Frequencies ---------");
    for i in nlow..frequencies_au.len() {
        let omega = frequencies_au[i];
        let freq = freq_to_cm1(omega);
        let mu = reduced_masses.and_then(|v| v.get(i)).copied();
        let ir = ir_intensities.and_then(|v| v.get(i)).copied();
        let mu_str = fmt_opt(mu, 12, 4);
        let ir_str = fmt_opt(ir, 12, 4);
        if omega < 0.0 {
            println!(
                "{:6} {:6} {:12.2} {mu_str} {ir_str} {:>10}",
                i as isize,
                i as isize - nlow as isize,
                freq,
                "<-- Imag"
            );
        } else {
            println!(
                "{:6} {:6} {:12.2} {mu_str} {ir_str}",
                i as isize,
                i as isize - nlow as isize,
                freq
            );
        }
    }
    println!("-----------------------------------");
    println!();

    let filename = format!("vibrational_freq_{name}.dat");
    let mut out = String::new();
    out.push_str("---------- Low Frequencies ----------\n");
    out.push_str(&format!(
        "{:6} {:6} {:12} {:12} {:12}\n",
        "index1", "index2", "freq[cm-1]", "mu_red[amu]", "IR[km/mol]"
    ));
    for i in 0..nlow {
        let omega = frequencies_au[i];
        let freq = freq_to_cm1(omega);
        let mu = reduced_masses.and_then(|v| v.get(i)).copied();
        let mu_str = fmt_opt(mu, 12, 4);
        let ir_str = fmt_opt(None, 12, 4);
        if omega < 0.0 && omega.abs() > 0.01 * C10 / C9 * (std::f64::consts::PI * 2.0) {
            out.push_str(&format!(
                "{:6} {:6} {:12.2} {mu_str} {ir_str} {:>10}\n",
                i as isize,
                i as isize - nlow as isize,
                freq,
                "<-- Imag"
            ));
        } else {
            out.push_str(&format!(
                "{:6} {:6} {:12.2} {mu_str} {ir_str}\n",
                i as isize,
                i as isize - nlow as isize,
                freq
            ));
        }
    }
    out.push_str("---------- High Frequencies --------- \n");
    for i in nlow..frequencies_au.len() {
        let omega = frequencies_au[i];
        let freq = freq_to_cm1(omega);
        let mu = reduced_masses.and_then(|v| v.get(i)).copied();
        let ir = ir_intensities.and_then(|v| v.get(i)).copied();
        let mu_str = fmt_opt(mu, 12, 4);
        let ir_str = fmt_opt(ir, 12, 4);
        if omega < 0.0 {
            out.push_str(&format!(
                "{:6} {:6} {:12.2} {mu_str} {ir_str} {:>10}\n",
                i as isize,
                i as isize - nlow as isize,
                freq,
                "<-- Imag"
            ));
        } else {
            out.push_str(&format!(
                "{:6} {:6} {:12.2} {mu_str} {ir_str}\n",
                i as isize,
                i as isize - nlow as isize,
                freq
            ));
        }
    }
    out.push_str("-----------------------------------\n");
    std::fs::write(&filename, out)?;
    Ok(())
}
