//! Hessians retain their own geometry and isotope masses; no display normalization.
use super::state::{Partner, Source};
use crate::normalmode::{cartesian, EckartMode};
use ndarray::Array2;
use std::path::Path;
pub const BOHR: f64 = 0.529177210903;
pub const AMU: f64 = 1822.888486209;
pub const CM_TO_AU: f64 = 137.035999084 * BOHR / 1e8 * std::f64::consts::TAU;
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Fragment {
    pub atoms: Vec<String>,
    pub masses: Vec<f64>,
    pub equilibrium: Vec<f64>,
    pub frequencies: Vec<f64>,
    pub modes: Array2<f64>,
    pub imaginary_cm1: Vec<f64>,
    pub imaginary_modes: Array2<f64>,
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
pub fn from_hessian(
    atoms: Vec<String>,
    masses: Vec<f64>,
    x: Vec<f64>,
    hessian: Array2<f64>,
    name: String,
) -> Result<Fragment, String> {
    from_hessian_projected(atoms, masses, x, hessian, name, None)
}
/// A path point: the Cartesian path tangent (bohr) is projected out together
/// with translations and rotations, leaving 3N-7 modes orthogonal to the path.
pub fn from_hessian_on_path(
    atoms: Vec<String>,
    masses: Vec<f64>,
    x: Vec<f64>,
    hessian: Array2<f64>,
    name: String,
    tangent: &[f64],
) -> Result<Fragment, String> {
    if tangent.len() != x.len() {
        return Err("Path tangent does not match the geometry".into());
    }
    // The projector mass-weights a gradient as g/sqrt(m); m·dx gives sqrt(m)·dx,
    // the mass-weighted tangent. Its scale is irrelevant.
    let gradient: Vec<f64> = tangent
        .iter()
        .enumerate()
        .map(|(i, t)| t * masses[i / 3])
        .collect();
    from_hessian_projected(atoms, masses, x, hessian, name, Some(&gradient))
}
fn from_hessian_projected(
    atoms: Vec<String>,
    masses: Vec<f64>,
    mut x: Vec<f64>,
    hessian: Array2<f64>,
    name: String,
    path_gradient: Option<&[f64]>,
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
    let mode = if path_gradient.is_some() {
        EckartMode::ReactionPath
    } else {
        EckartMode::VibRot
    };
    let r = cartesian::analyze(&masses, &x, &hessian, mode, path_gradient)?;
    center(&mut x, &masses);
    let keep = &r.positive_indices;
    if keep.is_empty() && r.negative_indices.len() != 1 {
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
        imaginary_modes: Array2::from_shape_fn((n, r.negative_indices.len()), |(i, j)| {
            r.modes[(i, r.negative_indices[j])]
        }),
        imaginary_cm1: r
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
    let h = crate::qchem_interfaces::hessian_io::load(path, geometry)?;
    let name = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .into_owned();
    from_hessian(h.atoms, h.masses_amu, h.coords_bohr, h.hessian, name)
}

/// Load a Hessian for a partner. Embedded coordinates always win. A
/// separately loaded equilibrium XYZ is used only for a coordinate-free
/// matrix, and only once, so it can never silently pair with a later file.
/// Records where the coordinates came from.
pub fn load_partner_hessian(p: &mut Partner, path: &Path) -> Result<Fragment, String> {
    let (f, source) = match load_hessian(path, None) {
        Ok(f) => (f, "coordinates from the Hessian file".to_string()),
        Err(e) => match &p.geometry {
            Some(g) => (
                load_hessian(path, Some(g))?,
                format!(
                    "coordinates from {}",
                    p.geometry_name.as_deref().unwrap_or("the loaded XYZ")
                ),
            ),
            None => return Err(e),
        },
    };
    p.geometry = None;
    p.geometry_name = None;
    p.coordinates_source = Some(source);
    Ok(f)
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
                imaginary_cm1: vec![],
                imaginary_modes: Array2::zeros((0, 0)),
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
                imaginary_cm1: vec![],
                imaginary_modes: Array2::zeros((0, 0)),
                name: "Diatomic parameters".into(),
                diatomic: Some((re, mu, p.morse.then_some(p.dissociation_ev / 27.2114))),
            })
        }
    }
}
