//! Shared preparation for Frequency Analysis, thermochemistry and PES sampling.
use super::{EckartMode, NormalModeResult};
use ndarray::Array2;

/// Double precision, translation invariant linearity test. Use the farthest
/// atom for the axis so coincident/nearby initial atoms do not define it.
pub fn is_linear(x: &[f64]) -> bool {
    if x.len() <= 6 {
        return true;
    }
    let vectors: Vec<[f64; 3]> = x
        .chunks_exact(3)
        .map(|p| std::array::from_fn(|k| p[k] - x[k]))
        .collect();
    let norm2 = |v: &[f64; 3]| v.iter().map(|x| x * x).sum::<f64>();
    let axis = vectors
        .iter()
        .max_by(|a, b| norm2(a).total_cmp(&norm2(b)))
        .unwrap();
    let length = norm2(axis);
    if length < 1e-20 {
        return true;
    }
    vectors.iter().all(|v| {
        let cross = [
            axis[1] * v[2] - axis[2] * v[1],
            axis[2] * v[0] - axis[0] * v[2],
            axis[0] * v[1] - axis[1] * v[0],
        ];
        norm2(v) < 1e-20 || norm2(&cross) / (length * norm2(v)) < 1e-6
    })
}

pub fn analyze(
    masses: &[f64],
    coords_bohr: &[f64],
    hessian: &Array2<f64>,
    eckart: EckartMode,
    gradient: Option<&[f64]>,
) -> Result<NormalModeResult, String> {
    if coords_bohr.len() != 3 * masses.len()
        || coords_bohr.iter().any(|v| !v.is_finite())
        || masses.iter().any(|v| !v.is_finite() || *v <= 0.)
    {
        return Err("Normal modes need matching finite coordinates and positive masses.".into());
    }
    let hessian = crate::qchem_interfaces::hessian_io::symmetric(hessian)?;
    super::normal_modes_with_projection(
        masses,
        &hessian,
        is_linear(coords_bohr),
        Some(coords_bohr),
        eckart,
        gradient,
    )
    .map_err(|e| e.to_string())
}
