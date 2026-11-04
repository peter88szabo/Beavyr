//! Common Cartesian Hessian input and Smite-compatible export.
//! Numeric rows are unweighted Eh/bohr², atom-major x/y/z. Comment metadata
//! makes exports self-contained without changing numpy.loadtxt compatibility.
use ndarray::Array2;
use std::{fmt::Write, path::Path};

pub const BOHR: f64 = 0.529177210903;
const HEADER: &str = "# Beavyr Cartesian Hessian v1";

#[derive(Clone, Debug)]
pub struct HessianData {
    pub atoms: Vec<String>,
    pub masses_amu: Vec<f64>,
    pub coords_bohr: Vec<f64>,
    pub hessian: Array2<f64>,
    pub frequency_scale_factor: Option<f64>,
    pub ir_intensities: Vec<f64>,
}

/// Use the same symmetric matrix for PES sampling, frequency analysis and
/// export. Do not mass-weight, project or scale the saved force constants.
pub fn symmetric(matrix: &Array2<f64>) -> Result<Array2<f64>, String> {
    let n = matrix.nrows();
    if n == 0 || matrix.ncols() != n || n % 3 != 0 {
        return Err("Hessian must be a nonempty 3N × 3N Cartesian matrix.".into());
    }
    if matrix.iter().any(|v| !v.is_finite()) {
        return Err("Hessian contains a non-finite value.".into());
    }
    let scale = matrix.iter().fold(0_f64, |a, b| a.max(b.abs()));
    let mut difference = 0_f64;
    for i in 0..n {
        for j in 0..i {
            difference = difference.max((matrix[(i, j)] - matrix[(j, i)]).abs());
        }
    }
    // A generous corruption guard, not the old 1e-7 symmetry requirement.
    // Numerical derivatives (including older Behemoth output) can have
    // appreciable antisymmetric residuals. Smite averages the two triangles.
    if difference > 1e-7 + 0.25 * scale {
        return Err(format!("Hessian has a large asymmetry (max |H-Hᵀ| = {difference:.6e} Eh/bohr²; largest element = {scale:.6e}). Check the matrix ordering and calculation before sampling."));
    }
    Ok(Array2::from_shape_fn((n, n), |(i, j)| {
        0.5 * matrix[(i, j)] + 0.5 * matrix[(j, i)]
    }))
}

fn number(s: &str) -> Result<f64, String> {
    s.replace(['D', 'd'], "E")
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| format!("Invalid Hessian number: {s}"))
}

pub fn geometry(text: &str) -> Result<(Vec<String>, Vec<f64>), String> {
    let mut lines = text.lines();
    let n = lines
        .next()
        .ok_or("Empty XYZ")?
        .trim()
        .parse::<usize>()
        .map_err(|_| "Invalid XYZ atom count")?;
    lines.next().ok_or("Missing XYZ comment")?;
    let mut atoms = Vec::new();
    let mut x = Vec::new();
    for _ in 0..n {
        let f: Vec<_> = lines
            .next()
            .ok_or("Incomplete XYZ")?
            .split_whitespace()
            .collect();
        if f.len() < 4 || crate::molecule::is_dummy(f[0]) {
            return Err("Hessian geometry must contain real atoms with three coordinates.".into());
        }
        atoms.push(f[0].to_string());
        for v in &f[1..4] {
            x.push(number(v)?);
        }
    }
    if n == 0 || lines.any(|l| !l.trim().is_empty()) {
        return Err("Select one matching equilibrium XYZ, not a trajectory.".into());
    }
    Ok((atoms, x))
}

pub fn from_geometry(
    atoms: Vec<String>,
    coords_bohr: Vec<f64>,
    hessian: Array2<f64>,
) -> Result<HessianData, String> {
    let masses_amu = atoms
        .iter()
        .map(|s| crate::molecule::atomic_mass_amu(s).ok_or_else(|| format!("Unknown element: {s}")))
        .collect::<Result<Vec<_>, _>>()?;
    validate(HessianData {
        atoms,
        coords_bohr,
        masses_amu,
        hessian,
        frequency_scale_factor: None,
        ir_intensities: vec![],
    })
}

fn validate(mut data: HessianData) -> Result<HessianData, String> {
    let n = data.atoms.len();
    if n == 0
        || data.coords_bohr.len() != 3 * n
        || data.masses_amu.len() != n
        || data.hessian.dim() != (3 * n, 3 * n)
        || data.coords_bohr.iter().any(|v| !v.is_finite())
        || data.masses_amu.iter().any(|v| !v.is_finite() || *v <= 0.)
        || data
            .atoms
            .iter()
            .any(|s| crate::molecule::atomic_number(s).is_none() || crate::molecule::is_dummy(s))
    {
        return Err(
            "Hessian dimensions, real atoms, coordinates and positive isotope masses must match."
                .into(),
        );
    }
    data.hessian = symmetric(&data.hessian)?;
    Ok(data)
}

pub fn is_standard(text: &str) -> bool {
    text.lines().any(|l| l.trim() == HEADER)
}

pub fn is_self_contained(text: &str) -> bool {
    text.lines()
        .any(|l| l.trim() == "$atoms" || l.trim() == HEADER)
}

/// Explicit XYZ coordinates are in Å. Embedded coordinates and masses take
/// precedence, so a different molecule in the viewer can never override them.
pub fn parse(text: &str, xyz: Option<&(Vec<String>, Vec<f64>)>) -> Result<HessianData, String> {
    if text.lines().any(|l| l.trim() == "$atoms") {
        let h = super::orca_hess::parse_orca_hess(text)?;
        let ir_intensities = h.ir_intensities_in_eigenvalue_order();
        return validate(HessianData {
            atoms: h.atoms,
            masses_amu: h.masses_amu,
            coords_bohr: h.coords_bohr,
            hessian: h.hessian,
            frequency_scale_factor: h.frequency_scale_factor,
            ir_intensities,
        });
    }
    let embedded = is_standard(text);
    let mut data = if embedded {
        let mut atoms = vec![];
        let mut masses_amu = vec![];
        let mut coords_bohr = vec![];
        let mut scale = None;
        let mut ir = vec![];
        if !text
            .lines()
            .any(|l| l.trim() == "# units Eh/bohr^2; coordinates bohr; masses amu")
        {
            return Err("Missing or unsupported Hessian units.".into());
        }
        for line in text.lines() {
            let f: Vec<_> = line.split_whitespace().collect();
            if f.get(0..2) == Some(&["#", "atom"]) {
                if f.len() != 7 {
                    return Err("Invalid Hessian atom metadata.".into());
                }
                atoms.push(f[2].to_string());
                masses_amu.push(number(f[3])?);
                for v in &f[4..7] {
                    coords_bohr.push(number(v)?);
                }
            } else if f.get(0..2) == Some(&["#", "frequency_scale"]) {
                if f.len() != 3 {
                    return Err("Invalid frequency scale metadata.".into());
                }
                let value = number(f[2])?;
                if value <= 0. {
                    return Err("Frequency scale must be positive.".into());
                }
                scale = Some(value);
            } else if f.get(0..2) == Some(&["#", "ir_km_mol"]) {
                ir = f[2..].iter().map(|s| number(s)).collect::<Result<_, _>>()?;
            }
        }
        HessianData {
            atoms,
            masses_amu,
            coords_bohr,
            hessian: Array2::zeros((0, 0)),
            frequency_scale_factor: scale,
            ir_intensities: ir,
        }
    } else {
        let (atoms, coords) = xyz.ok_or("This matrix has no geometry. Load its matching XYZ equilibrium geometry (Å) first. Supported matrix formats are $hessian, Smite and Behemoth output. New Beavyr exports include geometry and masses.")?;
        let masses_amu = atoms
            .iter()
            .map(|s| {
                crate::molecule::atomic_mass_amu(s).ok_or_else(|| format!("Unknown element: {s}"))
            })
            .collect::<Result<_, _>>()?;
        HessianData {
            atoms: atoms.clone(),
            masses_amu,
            coords_bohr: coords.iter().map(|v| v / BOHR).collect(),
            hessian: Array2::zeros((0, 0)),
            frequency_scale_factor: None,
            ir_intensities: vec![],
        }
    };
    let n = 3 * data.atoms.len();
    let lines: Vec<_> = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect();
    let labels = super::behemoth::expected_labels(&data.atoms);
    data.hessian = if lines.first().is_some_and(|l| l.starts_with("$hessian")) {
        super::hessian_file::parse_turbomole_hessian(
            &lines.join("\n").replace(['D', 'd'], "E"),
            data.atoms.len(),
        )?
    } else if lines
        .iter()
        .any(|l| l.split_whitespace().next() == labels.first().map(String::as_str))
    {
        super::behemoth::parse_hessian_for(text, &data.atoms)?
    } else {
        // Smite's numeric matrix, or the older internal PySCF dimension + flat matrix.
        let dimensioned = lines
            .first()
            .is_some_and(|l| l.parse::<usize>().ok() == Some(n));
        let rows = &lines[usize::from(dimensioned)..];
        if !dimensioned
            && (rows.len() != n || rows.iter().any(|l| l.split_whitespace().count() != n))
        {
            return Err(format!("Expected {n} rows of {n} Hessian values (Smite format), a $hessian matrix, ORCA .hess, or Behemoth Hessian output."));
        }
        let values = rows
            .iter()
            .flat_map(|l| l.split_whitespace())
            .map(number)
            .collect::<Result<Vec<_>, _>>()?;
        if values.len() != n * n {
            return Err(format!(
                "Hessian has {} values; expected {}.",
                values.len(),
                n * n
            ));
        }
        Array2::from_shape_vec((n, n), values).map_err(|e| e.to_string())?
    };
    validate(data)
}

pub fn load(path: &Path, xyz: Option<&(Vec<String>, Vec<f64>)>) -> Result<HessianData, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("Cannot read {}: {e}", path.display()))?;
    parse(&text, xyz).map_err(|e| format!("{}: {e}", path.display()))
}

pub fn standard_text(data: &HessianData) -> Result<String, String> {
    let data = validate(data.clone())?;
    let mut out = format!("{HEADER}\n# units Eh/bohr^2; coordinates bohr; masses amu\n# unweighted Cartesian matrix; atom-major x y z; symmetric (H + H^T)/2\n");
    for ((atom, mass), xyz) in data
        .atoms
        .iter()
        .zip(&data.masses_amu)
        .zip(data.coords_bohr.chunks_exact(3))
    {
        writeln!(
            out,
            "# atom {atom} {mass:.17e} {:.17e} {:.17e} {:.17e}",
            xyz[0], xyz[1], xyz[2]
        )
        .unwrap();
    }
    if let Some(scale) = data.frequency_scale_factor {
        writeln!(out, "# frequency_scale {scale:.17e}").unwrap();
    }
    if !data.ir_intensities.is_empty() {
        out.push_str("# ir_km_mol");
        for v in &data.ir_intensities {
            write!(out, " {v:.17e}").unwrap();
        }
        out.push('\n');
    }
    for row in data.hessian.rows() {
        for (i, v) in row.iter().enumerate() {
            if i > 0 {
                out.push(' ');
            }
            write!(out, "{v:.17e}").unwrap();
        }
        out.push('\n');
    }
    Ok(out)
}

pub fn save(path: &Path, data: &HessianData) -> Result<(), String> {
    std::fs::write(path, standard_text(data)?)
        .map_err(|e| format!("Cannot save {}: {e}", path.display()))
}

#[cfg(test)]
#[path = "hessian_io_tests.rs"]
mod tests;
