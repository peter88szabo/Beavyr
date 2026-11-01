//! File inputs carry their own geometry; never infer it from the current viewer.
use super::{duschinsky::mass_weighted_superposition, CM1_TO_AU};
use crate::normalmode::{normal_modes_with_projection, EckartMode};
use crate::qchem_interfaces::{hessian_file, orca_engrad, orca_hess};
use ndarray::Array2;

#[derive(Clone, Default)]
pub struct HessianSource {
    pub path: String,
    pub xyz: String,
}
#[derive(Debug)]
pub struct StateModes {
    pub atoms: Vec<String>,
    pub masses: Vec<f64>,
    pub coords: Vec<f64>,
    pub modes: Array2<f64>,
    pub frequencies: Vec<f64>,
    pub selected: Vec<usize>,
    pub imaginary: usize,
}
fn read(path: &str) -> Result<String, String> {
    std::fs::read_to_string(path.trim()).map_err(|e| format!("Cannot read {path}: {e}"))
}
fn number(s: &str) -> Result<f64, String> {
    s.replace(['D', 'd'], "E")
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| format!("Invalid number: {s}"))
}
fn xyz(text: &str) -> Result<(Vec<String>, Vec<f64>), String> {
    let mut lines = text.lines();
    let n = lines
        .next()
        .ok_or("Empty XYZ")?
        .trim()
        .parse::<usize>()
        .map_err(|_| "Invalid XYZ atom count")?;
    lines.next().ok_or("Missing XYZ comment")?;
    let mut atoms = Vec::new();
    let mut coords = Vec::new();
    for _ in 0..n {
        let f: Vec<_> = lines
            .next()
            .ok_or("Incomplete XYZ")?
            .split_whitespace()
            .collect();
        if f.len() != 4 {
            return Err("XYZ needs an element and three coordinates per atom".into());
        }
        atoms.push(f[0].to_string());
        for s in &f[1..] {
            coords.push(number(s)? / 0.529177210903);
        }
    }
    if lines.any(|l| !l.trim().is_empty()) {
        return Err("Choose a single geometry, not an XYZ trajectory".into());
    }
    Ok((atoms, coords))
}
fn linear(x: &[f64]) -> bool {
    if x.len() <= 6 {
        return true;
    }
    let axis: Vec<f64> = x
        .chunks_exact(3)
        .map(|p| (0..3).map(|k| p[k] - x[k]).collect::<Vec<_>>())
        .max_by(|a, b| {
            a.iter()
                .map(|v| v * v)
                .sum::<f64>()
                .total_cmp(&b.iter().map(|v| v * v).sum::<f64>())
        })
        .unwrap();
    let norm: f64 = axis.iter().map(|v| v * v).sum();
    norm > 0.
        && x.chunks_exact(3).all(|p| {
            let a: [f64; 3] = std::array::from_fn(|k| p[k] - x[k]);
            let cross = [
                a[1] * axis[2] - a[2] * axis[1],
                a[2] * axis[0] - a[0] * axis[2],
                a[0] * axis[1] - a[1] * axis[0],
            ];
            cross.iter().map(|v| v * v).sum::<f64>() / norm < 1e-10
        })
}
pub fn load(source: &HessianSource, cutoff: f64) -> Result<StateModes, String> {
    let text = read(&source.path)?;
    let (atoms, masses, coords, hessian) = if text.lines().any(|l| l.trim() == "$atoms") {
        let h = orca_hess::parse_orca_hess(&text)?;
        (h.atoms, h.masses_amu, h.coords_bohr, h.hessian)
    } else {
        if source.xyz.trim().is_empty() {
            return Err("This Hessian has no geometry. Select its matching XYZ file (Å).".into());
        }
        let (atoms, coords) = xyz(&read(&source.xyz)?)?;
        let masses = atoms
            .iter()
            .map(|a| {
                crate::molecule::atomic_mass_amu(a).ok_or_else(|| format!("Unknown element: {a}"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let h = hessian_file::parse_turbomole_hessian(&text.replace(['D', 'd'], "E"), atoms.len())?;
        (atoms, masses, coords, h)
    };
    if atoms.len() < 2
        || masses.len() != atoms.len()
        || coords.len() != 3 * atoms.len()
        || hessian.dim() != (coords.len(), coords.len())
        || atoms.iter().any(|a| crate::molecule::is_dummy(a))
        || masses.iter().any(|v| !v.is_finite() || *v <= 0.)
        || coords.iter().chain(hessian.iter()).any(|v| !v.is_finite())
    {
        return Err(
            "Invalid Hessian, atom masses or geometry. Hessians must contain real atoms only."
                .into(),
        );
    }
    let scale = hessian.iter().fold(1_f64, |a, b| a.max(b.abs()));
    for i in 0..hessian.nrows() {
        for j in 0..i {
            if (hessian[(i, j)] - hessian[(j, i)]).abs() > 1e-7 * scale {
                return Err("Hessian is not symmetric".into());
            }
        }
    }
    let nm = normal_modes_with_projection(
        &masses,
        &hessian,
        linear(&coords),
        Some(&coords),
        EckartMode::VibRot,
        None,
    )
    .map_err(|e| e.to_string())?;
    let frequencies: Vec<_> = nm.frequencies_au.iter().map(|f| f / CM1_TO_AU).collect();
    let selected: Vec<_> = nm
        .positive_indices
        .iter()
        .copied()
        .filter(|&i| frequencies[i] > cutoff)
        .collect();
    if selected.is_empty() {
        return Err("No positive vibrational modes remain above the frequency cutoff".into());
    }
    Ok(StateModes {
        atoms,
        masses,
        coords,
        modes: nm.modes,
        frequencies,
        selected,
        imaginary: nm.negative_indices.len(),
    })
}
pub fn same_atoms(a: &[String], b: &[String]) -> Result<(), String> {
    if a.len() != b.len() || a.iter().zip(b).any(|(a, b)| !a.eq_ignore_ascii_case(b)) {
        Err("The files must contain the same atoms in the same order.".into())
    } else {
        Ok(())
    }
}
/// Parse one Turbomole gradient cycle, retaining its geometry for validation.
pub(super) fn turbomole_gradient(text: &str) -> Result<orca_engrad::OrcaGradient, String> {
    if !text.lines().any(|l| l.trim() == "$grad") {
        return Err("Expected an ORCA .engrad or a Turbomole $grad file".into());
    }
    let rows: Vec<Vec<&str>> = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('$') && !l.starts_with("cycle"))
        .map(|l| l.split_whitespace().collect())
        .collect();
    let n = rows.len() / 2;
    if n == 0 || rows.len() != 2 * n {
        return Err("Gradient must contain one geometry and one gradient cycle".into());
    }
    let mut atoms = vec![];
    let mut coords = vec![];
    let mut gradient = vec![];
    for row in &rows[..n] {
        if row.len() != 4 {
            return Err("Gradient coordinate rows need x y z element, in bohr".into());
        }
        atoms.push(row[3].to_string());
        for v in &row[..3] {
            coords.push(number(v)?);
        }
    }
    for row in &rows[n..] {
        if row.len() != 3 {
            return Err("Gradient rows need three components in Hartree/bohr".into());
        }
        for v in row {
            gradient.push(number(v)?);
        }
    }
    Ok(orca_engrad::OrcaGradient {
        atoms,
        coords_bohr: coords,
        gradient_bohr: gradient,
        energy_hartree: 0.,
    })
}
pub fn gradient(path: &str, reference: &StateModes) -> Result<Vec<f64>, String> {
    let text = read(path)?;
    let grad = if text.lines().any(|l| l.trim() == "$grad") {
        turbomole_gradient(&text)?
    } else {
        orca_engrad::parse_orca_engrad(&text)?
    };
    aligned_gradient(&grad, reference)
}
pub(super) fn aligned_gradient(
    grad: &orca_engrad::OrcaGradient,
    reference: &StateModes,
) -> Result<Vec<f64>, String> {
    same_atoms(&reference.atoms, &grad.atoms)?;
    if grad.coords_bohr.len() != reference.coords.len()
        || grad.gradient_bohr.len() != reference.coords.len()
        || grad
            .coords_bohr
            .iter()
            .chain(&grad.gradient_bohr)
            .any(|v| !v.is_finite())
    {
        return Err("Invalid gradient geometry or components".into());
    }
    let (r, mobile, fixed) =
        mass_weighted_superposition(&grad.coords_bohr, &reference.coords, &reference.masses)
            .map_err(|e| e.to_string())?;
    let max = mobile
        .chunks_exact(3)
        .zip(fixed.chunks_exact(3))
        .map(|(a, b)| {
            a.iter()
                .zip(b)
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .fold(0_f64, f64::max);
    if max > orca_engrad::GEOMETRY_TOLERANCE_BOHR {
        return Err(format!("The excited-state gradient was not calculated at the ground-state Hessian geometry (maximum aligned difference {max:.6} bohr). Use a gradient at that same geometry."));
    }
    Ok(grad
        .gradient_bohr
        .chunks_exact(3)
        .flat_map(|p| (0..3).map(move |i| (0..3).map(|j| r[i][j] * p[j]).sum::<f64>()))
        .collect())
}
