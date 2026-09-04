#![allow(dead_code)]
use anyhow::{bail, Result};
use ndarray::Array2;

// Behemoth's `ir_intensities_fd_xtb` (and its `dipole_from_xtb_geometry`
// helper) supplied dipoles from Behemoth's *own* xTB Hamiltonian and SCF.
// Beavyr shells out to the xtb binary and has no internal Hamiltonian, so
// importing them would mean dragging in Behemoth's Hamiltonian, properties
// and SCF stacks. They are therefore the only part of this file left
// behind; the generic `ir_intensities_fd` below takes any dipole closure,
// and IR intensities in practice come from parsing xtb's own `vibspectrum`.
// See docs/superpowers/specs/2026-09-04-frequency-analysis-design.md.

/// From Behemoth's `properties::dipolemoment`, inlined so this file no
/// longer depends on that module.
const AU_TO_DEBYE: f64 = 2.541_746_473;

const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;
const AMU_TO_ELECTRON_MASS: f64 = 1_822.888_486_209;
const IR_INTENSITY_DEBYE_ANG_TO_KM_MOL: f64 = 42.255;

pub const DEFAULT_IR_FD_STEP_Q: f64 = 1.0e-2;

/// Standard Wilson-Decius-Cross normal-mode reduced mass,
/// `mu_k = 1 / sum_i L_ik^2`, where `L` is the *un-normalized* Cartesian
/// displacement eigenvector satisfying `L^T M L = I` (this is exactly what
/// `modes` already is: the mass-weighted eigenvector divided by
/// `sqrt(mass_au)`, not renormalized to a unit Cartesian norm). Renormalizing
/// each mode to unit Cartesian length first and then averaging `1/mass` over
/// its components — as an earlier version of this function did — computes a
/// different, non-standard quantity that silently collapses towards the
/// lightest atom's mass for every mode; the fix was validated against
/// PySCF's `pyscf.hessian.thermo.harmonic_analysis` (which uses the same
/// `1/sum(norm_mode**2)` formula) on water/HF/STO-3G, where the old formula
/// gave ~1.01 amu for all three modes and the correct one reproduces PySCF's
/// 1.08165/1.04647/1.08202 amu.
///
/// `modes` is built from a Hessian mass-weighted with `mass_amu * `
/// [`AMU_TO_ELECTRON_MASS`] (atomic units), so `sum_i L_ik^2` comes out in
/// electron-mass units; dividing by `AMU_TO_ELECTRON_MASS` converts the
/// result back to amu.
pub fn reduced_masses_amu(masses_amu: &[f64], modes: &Array2<f64>) -> Result<Vec<f64>> {
    let (ncoord, nmode) = modes.dim();
    if ncoord == 0 {
        bail!("normal modes are empty");
    }
    if ncoord != 3 * masses_amu.len() {
        bail!("modes length must be 3N (N = number of atoms)");
    }
    let mut reduced = vec![0.0f64; nmode];
    for mode in 0..nmode {
        let mut sum_sq = 0.0;
        for i in 0..ncoord {
            let v = modes[(i, mode)];
            sum_sq += v * v;
        }
        reduced[mode] = if sum_sq > 0.0 {
            1.0 / (sum_sq * AMU_TO_ELECTRON_MASS)
        } else {
            0.0
        };
    }
    Ok(reduced)
}


pub fn ir_intensities_fd<F>(
    coords_bohr: &[f64],
    modes: &Array2<f64>,
    fd_step_q: f64,
    mut dipole_fn: F,
) -> Result<Vec<f64>>
where
    F: FnMut(&[f64]) -> Result<[f64; 3]>,
{
    if fd_step_q <= 0.0 {
        bail!("fd_step_q must be positive");
    }
    let ncoord = coords_bohr.len();
    let (nr, nmode) = modes.dim();
    if nr != ncoord {
        bail!("modes row count must match coords length");
    }
    let mut intensities = vec![0.0f64; nmode];
    for mode in 0..nmode {
        let mut coords_p = coords_bohr.to_vec();
        let mut coords_m = coords_bohr.to_vec();
        for i in 0..ncoord {
            let disp = modes[(i, mode)] * fd_step_q;
            coords_p[i] += disp;
            coords_m[i] -= disp;
        }
        let mu_p = dipole_fn(&coords_p)?;
        let mu_m = dipole_fn(&coords_m)?;
        let dmu = [
            (mu_p[0] - mu_m[0]) / (2.0 * fd_step_q),
            (mu_p[1] - mu_m[1]) / (2.0 * fd_step_q),
            (mu_p[2] - mu_m[2]) / (2.0 * fd_step_q),
        ];
        let dmu_norm = (dmu[0] * dmu[0] + dmu[1] * dmu[1] + dmu[2] * dmu[2]).sqrt();
        let dmu_debye_per_ang_sqrt_amu =
            dmu_norm * AU_TO_DEBYE * AMU_TO_ELECTRON_MASS.sqrt() / BOHR_TO_ANGSTROM;
        intensities[mode] = IR_INTENSITY_DEBYE_ANG_TO_KM_MOL
            * dmu_debye_per_ang_sqrt_amu
            * dmu_debye_per_ang_sqrt_amu;
    }
    Ok(intensities)
}



#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn ir_intensity_converts_mass_weighted_au_derivative_to_practical_units() {
        let coords = vec![0.0, 0.0, 0.0];
        let modes = Array2::from_shape_vec((3, 1), vec![1.0, 0.0, 0.0]).unwrap();
        let intensities =
            ir_intensities_fd(&coords, &modes, 1.0e-3, |x| Ok([x[0], 0.0, 0.0])).unwrap();
        let derivative = AU_TO_DEBYE * AMU_TO_ELECTRON_MASS.sqrt() / BOHR_TO_ANGSTROM;
        let expected = IR_INTENSITY_DEBYE_ANG_TO_KM_MOL * derivative * derivative;
        assert!((intensities[0] - expected).abs() / expected < 1.0e-12);
    }

    #[test]
    fn reduced_mass_of_a_diatomic_stretch_matches_the_wilson_decius_cross_formula() {
        // Two masses m1, m2 (amu) connected by a spring along z. Center-of-
        // mass conservation alone fixes the stretch mode's Cartesian
        // direction (independent of the force constant): m1*x1 + m2*x2 = 0,
        // so the mass-weighted eigenvector points along
        // (-m2/sqrt(m1), sqrt(m2)). For this direction the standard
        // (Wilson-Decius-Cross) normal-mode reduced mass has the closed
        // form mu = m1*m2*(m1+m2)/(m1^2+m2^2), independently confirmed
        // against pyscf.hessian.thermo.harmonic_analysis on real HF/STO-3G
        // (masses 1.008/18.998403163 amu give 1.058501806766436, matching
        // PySCF's reported 1.05850181 to 8 digits).
        let m1 = 2.0_f64;
        let m2 = 3.0_f64;
        let expected_mu = m1 * m2 * (m1 + m2) / (m1 * m1 + m2 * m2);

        let a = (m1 / (m2 * (m1 + m2))).sqrt();
        let l1_amu_convention = -m2 * a / m1; // atom 0, z
        let l2_amu_convention = a; // atom 1, z

        // `modes` (as produced by `normal_modes`) is un-mass-weighted with
        // mass_au = mass_amu * AMU_TO_ELECTRON_MASS, not mass_amu directly,
        // so convert the amu-convention eigenvector accordingly.
        let scale = 1.0 / AMU_TO_ELECTRON_MASS.sqrt();
        let mut modes = Array2::<f64>::zeros((6, 1));
        modes[(2, 0)] = l1_amu_convention * scale;
        modes[(5, 0)] = l2_amu_convention * scale;

        let reduced = reduced_masses_amu(&[m1, m2], &modes).unwrap();
        assert!((reduced[0] - expected_mu).abs() / expected_mu < 1.0e-12);
    }
}
