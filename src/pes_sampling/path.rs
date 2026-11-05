//! Reaction-path sampling: normal modes orthogonal to an IRC.
//!
//! Each path point is a self-contained Hessian (geometry and masses embedded).
//! The path tangent at point k is the difference of its neighbours'
//! geometries after mass-weighted superposition onto point k (central
//! difference inside, one-sided at the ends). It is projected out together
//! with translations and rotations (Miller-Handy-Adams reaction-path
//! Hamiltonian), so the remaining 3N-7 modes are sampled at every point.
use super::input::{self, Fragment};
use crate::normalmode::linalg_shim::{Backend, LinAlg};
use ndarray::Array2;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct PathSettings {
    pub points: Vec<Fragment>,
    pub samples_per_point: usize,
    pub excitation: super::state::Excitation,
}

impl Default for PathSettings {
    fn default() -> Self {
        Self {
            points: vec![],
            samples_per_point: 20,
            excitation: Default::default(),
        }
    }
}

/// "irc_2" sorts before "irc_10".
pub fn natural_key(name: &str) -> Vec<(bool, u128, String)> {
    let mut key = vec![];
    let mut chars = name.chars().peekable();
    while let Some(&c) = chars.peek() {
        let digit = c.is_ascii_digit();
        let mut chunk = String::new();
        while let Some(&c) = chars.peek() {
            if c.is_ascii_digit() != digit {
                break;
            }
            chunk.push(c);
            chars.next();
        }
        let number = if digit { chunk.parse().unwrap_or(u128::MAX) } else { 0 };
        key.push((digit, number, if digit { String::new() } else { chunk.to_lowercase() }));
    }
    key
}

fn centered(x: &[f64], m: &[f64]) -> Vec<f64> {
    let mut x = x.to_vec();
    input::center(&mut x, m);
    x
}

/// Rotate centered `moving` onto centered `reference` (mass-weighted Kabsch,
/// via Horn's quaternion eigenproblem).
pub fn superpose(moving: &[f64], reference: &[f64], m: &[f64]) -> Result<Vec<f64>, String> {
    let mut s = [[0f64; 3]; 3];
    for ((p, q), w) in moving.chunks_exact(3).zip(reference.chunks_exact(3)).zip(m) {
        for a in 0..3 {
            for b in 0..3 {
                s[a][b] += w * p[a] * q[b];
            }
        }
    }
    let [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]] = s;
    let n = Array2::from_shape_vec(
        (4, 4),
        vec![
            xx + yy + zz, yz - zy, zx - xz, xy - yx,
            yz - zy, xx - yy - zz, xy + yx, zx + xz,
            zx - xz, xy + yx, -xx + yy - zz, yz + zy,
            xy - yx, zx + xz, yz + zy, -xx - yy + zz,
        ],
    )
    .unwrap();
    let la = LinAlg::new(Backend::Auto).map_err(|e| e.to_string())?;
    let (values, vectors) = la.eigh(&n).map_err(|e| e.to_string())?;
    let k = (0..4).max_by(|&a, &b| values[a].total_cmp(&values[b])).unwrap();
    let (w, x, y, z) = (vectors[(0, k)], vectors[(1, k)], vectors[(2, k)], vectors[(3, k)]);
    let r = [
        [w * w + x * x - y * y - z * z, 2. * (x * y - w * z), 2. * (x * z + w * y)],
        [2. * (x * y + w * z), w * w - x * x + y * y - z * z, 2. * (y * z - w * x)],
        [2. * (x * z - w * y), 2. * (y * z + w * x), w * w - x * x - y * y + z * z],
    ];
    Ok(moving
        .chunks_exact(3)
        .flat_map(|p| (0..3).map(move |i| (0..3).map(|k| r[i][k] * p[k]).sum::<f64>()))
        .collect())
}

/// Path tangents (Cartesian, bohr) at each point, in that point's own frame.
pub fn tangents(geometries: &[Vec<f64>], m: &[f64]) -> Result<Vec<Vec<f64>>, String> {
    let n = geometries.len();
    if n < 2 {
        return Err("Reaction-path sampling needs at least two path points.".into());
    }
    let mut result = vec![];
    for k in 0..n {
        let here = centered(&geometries[k], m);
        let (a, b) = (k.saturating_sub(1), (k + 1).min(n - 1));
        let after = superpose(&centered(&geometries[b], m), &here, m)?;
        let before = superpose(&centered(&geometries[a], m), &here, m)?;
        let t: Vec<f64> = after.iter().zip(&before).map(|(p, q)| p - q).collect();
        let norm = t.iter().map(|v| v * v).sum::<f64>().sqrt();
        if !(norm > 1e-6) {
            return Err(format!("Path points {} and {} have the same geometry; remove duplicates.", a + 1, b + 1));
        }
        result.push(t);
    }
    Ok(result)
}

/// Load self-contained Hessians along a path, ordered naturally by file name.
pub fn load(files: &[PathBuf]) -> Result<Vec<Fragment>, String> {
    let mut files = files.to_vec();
    files.sort_by_key(|p| natural_key(&p.file_name().unwrap_or_default().to_string_lossy()));
    let mut data = vec![];
    for path in &files {
        let h = crate::qchem_interfaces::hessian_io::load(path, None).map_err(|e| {
            format!("{e}\nReaction-path points need self-contained Hessians (ORCA .hess or Beavyr Hessian exports).")
        })?;
        data.push((name(path), h));
    }
    build(data)
}

fn name(path: &Path) -> String {
    path.file_name().unwrap_or_default().to_string_lossy().into_owned()
}

pub fn build(data: Vec<(String, crate::qchem_interfaces::hessian_io::HessianData)>) -> Result<Vec<Fragment>, String> {
    let Some((_, first)) = data.first() else {
        return Err("Choose the Hessian files along the reaction path.".into());
    };
    let (atoms, masses) = (first.atoms.clone(), first.masses_amu.clone());
    if let Some((n, _)) = data.iter().find(|(_, h)| h.atoms != atoms || h.masses_amu != masses) {
        return Err(format!("{n}: all path points must have the same atoms, order and masses."));
    }
    let geometries: Vec<Vec<f64>> = data.iter().map(|(_, h)| h.coords_bohr.clone()).collect();
    let tangents = tangents(&geometries, &masses)?;
    data.into_iter()
        .zip(tangents)
        .enumerate()
        .map(|(k, ((n, h), t))| {
            input::from_hessian_on_path(h.atoms, h.masses_amu, h.coords_bohr, h.hessian, n, &t)
                .map_err(|e| format!("Path point {}: {e}", k + 1))
        })
        .collect()
}
