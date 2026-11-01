//! Hessians retain their own geometry and isotope masses; no display normalization.
use super::state::{Partner, Source};
use crate::normalmode::{normal_modes_with_projection, EckartMode};
use ndarray::Array2;
use std::path::Path;
pub const BOHR: f64 = 0.529177210903;
pub const AMU: f64 = 1822.888486209;
pub const CM_TO_AU: f64 = 137.035999084 * BOHR / 1e8 * std::f64::consts::TAU;
#[derive(Clone, Debug)]
pub struct Fragment {
    pub atoms: Vec<String>,
    pub masses: Vec<f64>,
    pub equilibrium: Vec<f64>,
    pub frequencies: Vec<f64>,
    pub modes: Array2<f64>,
    pub excluded: Vec<f64>,
    pub name: String,
    pub diatomic: Option<(f64, f64, Option<f64>)>,
}
pub fn masses(atoms: &[String]) -> Result<Vec<f64>, String> {
    atoms
        .iter()
        .map(|s| crate::molecule::atomic_mass_amu(s).ok_or_else(|| format!("Unknown element: {s}")))
        .collect()
}
pub fn center(x: &mut [f64], m: &[f64]) {
    let total: f64 = m.iter().sum();
    let mut c = [0.; 3];
    for (p, w) in x.chunks_exact(3).zip(m) {
        for k in 0..3 {
            c[k] += p[k] * w / total;
        }
    }
    for p in x.chunks_exact_mut(3) {
        for k in 0..3 {
            p[k] -= c[k];
        }
    }
}
pub fn geometry(text: &str) -> Result<(Vec<String>, Vec<f64>), String> {
    let mut lines = text.lines();
    let n = lines
        .next()
        .ok_or("Empty XYZ")?
        .trim()
        .parse::<usize>()
        .map_err(|_| "Invalid XYZ atom count")?;
    lines.next().ok_or("XYZ needs a comment line")?;
    let mut atoms = vec![];
    let mut x = vec![];
    for _ in 0..n {
        let line = lines.next().ok_or("Incomplete XYZ")?;
        let f: Vec<_> = line.split_whitespace().collect();
        if f.len() < 4 {
            return Err("XYZ rows need an element and three coordinates".into());
        }
        let p = f[1..4]
            .iter()
            .map(|v| v.parse::<f64>().map_err(|_| "Invalid XYZ coordinate"))
            .collect::<Result<Vec<_>, _>>()?;
        if p.iter().any(|v| !v.is_finite()) {
            return Err("XYZ coordinates must be finite".into());
        }
        if !crate::molecule::is_dummy(f[0]) {
            atoms.push(f[0].to_string());
            x.extend(p);
        }
    }
    if atoms.is_empty() {
        return Err("The structure contains no real atoms".into());
    }
    if lines.any(|s| !s.trim().is_empty()) {
        return Err("Choose a single equilibrium geometry, not a multi-frame XYZ".into());
    }
    masses(&atoms)?;
    Ok((atoms, x))
}
pub fn linear(x: &[f64]) -> bool {
    let n = x.len() / 3;
    if n <= 2 {
        return true;
    }
    let mut axis = [0.; 3];
    let mut norm = 0.;
    for p in x.chunks_exact(3) {
        let a = [p[0] - x[0], p[1] - x[1], p[2] - x[2]];
        let l = a.iter().map(|v| v * v).sum::<f64>();
        if l > norm {
            norm = l;
            axis = a;
        }
    }
    norm > 0.
        && x.chunks_exact(3).all(|p| {
            let a = [p[0] - x[0], p[1] - x[1], p[2] - x[2]];
            let cross = [
                a[1] * axis[2] - a[2] * axis[1],
                a[2] * axis[0] - a[0] * axis[2],
                a[0] * axis[1] - a[1] * axis[0],
            ];
            cross.iter().map(|v| v * v).sum::<f64>() / norm < 1e-10
        })
}
pub fn from_hessian(
    atoms: Vec<String>,
    masses: Vec<f64>,
    mut x: Vec<f64>,
    hessian: Array2<f64>,
    name: String,
) -> Result<Fragment, String> {
    if atoms.len() < 2
        || masses.len() != atoms.len()
        || x.len() != 3 * atoms.len()
        || masses.iter().any(|m| !m.is_finite() || *m <= 0.)
        || x.iter().chain(hessian.iter()).any(|v| !v.is_finite())
    {
        return Err("Invalid geometry, masses or Hessian".into());
    }
    if atoms.iter().any(|s| crate::molecule::is_dummy(s)) {
        return Err("Use a Hessian and matching geometry containing real atoms only".into());
    }
    let n = x.len();
    if hessian.dim() != (n, n) {
        return Err("Hessian dimensions must match the equilibrium geometry (3N × 3N)".into());
    }
    let max = hessian.iter().fold(0_f64, |a, b| a.max(b.abs()));
    for i in 0..n {
        for j in 0..i {
            if (hessian[(i, j)] - hessian[(j, i)]).abs() > 1e-7 * max.max(1.) {
                return Err("Hessian is not symmetric".into());
            }
        }
    }
    center(&mut x, &masses);
    let r = normal_modes_with_projection(
        &masses,
        &hessian,
        linear(&x),
        Some(&x),
        EckartMode::VibRot,
        None,
    )
    .map_err(|e| e.to_string())?;
    let keep = &r.positive_indices;
    if keep.is_empty() {
        return Err(
            "No positive vibrational modes remain after translation/rotation projection".into(),
        );
    }
    let modes = Array2::from_shape_fn((n, keep.len()), |(i, j)| r.modes[(i, keep[j])]);
    Ok(Fragment {
        atoms,
        masses,
        equilibrium: x,
        frequencies: keep.iter().map(|i| r.frequencies_au[*i]).collect(),
        modes,
        excluded: r
            .negative_indices
            .iter()
            .map(|i| r.frequencies_au[*i] / CM_TO_AU)
            .collect(),
        name,
        diatomic: None,
    })
}
pub fn load_hessian(
    path: &Path,
    geometry: Option<&(Vec<String>, Vec<f64>)>,
) -> Result<Fragment, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("Cannot read {}: {e}", path.display()))?;
    let name = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .into_owned();
    if text.lines().any(|l| l.trim() == "$atoms") {
        let h = crate::qchem_interfaces::orca_hess::parse_orca_hess(&text)?;
        from_hessian(h.atoms, h.masses_amu, h.coords_bohr, h.hessian, name)
    } else {
        let (atoms,x)=geometry.ok_or("This Hessian contains no geometry. Load its matching equilibrium XYZ first (coordinates in Å; Hessian in Eh/bohr²).")?;
        let h = crate::qchem_interfaces::hessian_file::parse_turbomole_hessian(&text, atoms.len())?;
        from_hessian(
            atoms.clone(),
            masses(atoms)?,
            x.iter().map(|v| v / BOHR).collect(),
            h,
            name,
        )
    }
}
pub fn prepare(p: &Partner) -> Result<Fragment, String> {
    match p.source {
        Source::Hessian => p
            .fragment
            .clone()
            .ok_or_else(|| "Load an equilibrium Hessian for this partner".into()),
        Source::Atom => {
            let atoms = vec![p.symbols[0].trim().to_string()];
            let masses = masses(&atoms)?;
            Ok(Fragment {
                atoms,
                masses,
                equilibrium: vec![0.; 3],
                frequencies: vec![],
                modes: Array2::zeros((3, 0)),
                excluded: vec![],
                name: "Atom".into(),
                diatomic: None,
            })
        }
        Source::Diatomic => {
            if !p.equilibrium.is_finite()
                || p.equilibrium <= 0.
                || !p.frequency.is_finite()
                || p.frequency <= 0.
                || p.morse && (!p.dissociation_ev.is_finite() || p.dissociation_ev <= 0.)
            {
                return Err("Diatomic equilibrium distance, frequency and Morse dissociation energy must be positive and finite".into());
            }
            let atoms: Vec<_> = p.symbols.iter().map(|s| s.trim().to_string()).collect();
            let masses = masses(&atoms)?;
            let re = p.equilibrium / BOHR;
            let mut x = vec![re, 0., 0., 0., 0., 0.];
            center(&mut x, &masses);
            let mu = masses[0] * masses[1] / (masses[0] + masses[1]) * AMU;
            let mut modes = Array2::zeros((6, 1));
            modes[(0, 0)] = masses[1] / (masses[0] + masses[1]) / mu.sqrt();
            modes[(3, 0)] = -masses[0] / (masses[0] + masses[1]) / mu.sqrt();
            Ok(Fragment {
                atoms,
                masses,
                equilibrium: x,
                frequencies: vec![p.frequency * CM_TO_AU],
                modes,
                excluded: vec![],
                name: "Diatomic parameters".into(),
                diatomic: Some((re, mu, p.morse.then_some(p.dissociation_ev / 27.2114))),
            })
        }
    }
}
