//! Parses the extra files xTB writes for a Hessian run (`xtb --hess` /
//! `--ohess`): the Turbomole-style `hessian` matrix and the `vibspectrum`
//! block. Both formats were confirmed by actually running xtb 6.7.1 against
//! `examples/h2o2.xyz` (see the design spec); `tests/fixtures/h2o2_hessian`
//! and `tests/fixtures/h2o2_vibspectrum` are that real, captured output.

use ndarray::Array2;

/// Parses a Turbomole-style `$hessian` file: a header line, then `(3N)^2`
/// free-format values (row-major) for an N-atom Cartesian Hessian, in
/// Eh/bohr^2.
pub fn parse_turbomole_hessian(text: &str, natoms: usize) -> Result<Array2<f64>, String> {
    if natoms == 0 {
        return Err("Cannot parse a Hessian for zero atoms.".to_string());
    }
    let ndim = 3 * natoms;
    let expected = ndim * ndim;

    let mut values = Vec::with_capacity(expected);
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Err("Hessian file is empty.".to_string());
    };
    if !header.trim_start().starts_with('$') {
        return Err(format!(
            "Expected a Turbomole-style '$hessian' header, found: {:?}",
            header.trim()
        ));
    }

    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('$') {
            continue;
        }
        for tok in trimmed.split_whitespace() {
            let v: f64 = tok
                .parse()
                .map_err(|_| format!("Could not parse Hessian value {tok:?}"))?;
            values.push(v);
        }
    }

    if values.len() != expected {
        return Err(format!(
            "Hessian has {} values; expected {expected} for {natoms} atom(s) (3N x 3N).",
            values.len()
        ));
    }

    Array2::from_shape_vec((ndim, ndim), values)
        .map_err(|e| format!("Failed to shape the Hessian matrix: {e}"))
}

/// Parses xTB's Turbomole-style `$grad` file into a flat Cartesian gradient
/// (Eh/bohr, row-major x,y,z per atom). The file lists N coordinate lines
/// (ignored here -- the geometry comes from `xtbhess.xyz` instead, since
/// that is the frame the Hessian and gradient both share) followed by N
/// gradient lines, bracketed by `$grad` / `$end`.
pub fn parse_turbomole_gradient(text: &str, natoms: usize) -> Result<Vec<f64>, String> {
    let data_lines: Vec<&str> = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('$') && !l.starts_with("cycle"))
        .collect();
    if data_lines.len() != 2 * natoms {
        return Err(format!(
            "Gradient file has {} data line(s); expected {} (coords + gradient) for {natoms} atom(s).",
            data_lines.len(),
            2 * natoms
        ));
    }

    let mut gradient = Vec::with_capacity(3 * natoms);
    for line in &data_lines[natoms..] {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() < 3 {
            return Err(format!("Malformed gradient line: {line:?}"));
        }
        for tok in &tokens[..3] {
            let v: f64 = tok
                .parse()
                .map_err(|_| format!("Could not parse gradient value {tok:?}"))?;
            gradient.push(v);
        }
    }
    Ok(gradient)
}

/// Parses the "Standard orientation" geometry block from xTB's Gaussian-98
/// style `g98.out`: element (from atomic number) and Cartesian coordinates
/// (Angstrom), in the exact frame the Hessian was computed in.
///
/// This is the reliable source for that geometry -- `xtbhess.xyz`, tried
/// first in an earlier version of this code, turned out to be written only
/// for some point groups. Reproduced directly: `xtb --hess` on
/// `examples/h2o2.xyz` (a small, symmetric molecule) writes `xtbhess.xyz`,
/// but the identical command on `examples/phenyl-OCH3.xyz` (larger, Cs
/// symmetry) does not, while both runs write `g98.out`. Parsing is
/// pattern-based (six-token rows of `center, Z, type, x, y, z`) rather than
/// counting fixed header lines, so it does not depend on the exact banner
/// text staying byte-for-byte the same across xTB versions.
pub fn parse_g98_geometry(text: &str) -> Result<(Vec<String>, Vec<f64>), String> {
    let Some(start) = text.find("Standard orientation:") else {
        return Err("g98.out has no 'Standard orientation' geometry block.".to_string());
    };

    let mut symbols = Vec::new();
    let mut coords = Vec::with_capacity(3);
    let mut in_rows = false;
    for line in text[start..].lines().skip(1) {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        let row = tokens.len() == 6
            && tokens[..3].iter().all(|t| t.parse::<i64>().is_ok())
            && tokens[3..].iter().all(|t| t.parse::<f64>().is_ok());
        if row {
            in_rows = true;
            let z: u32 = tokens[1].parse().unwrap();
            let symbol = crate::molecule::element_symbol(z)
                .ok_or_else(|| format!("g98.out uses an unknown atomic number {z}"))?;
            symbols.push(symbol.to_string());
            for tok in &tokens[3..] {
                coords.push(tok.parse::<f64>().unwrap());
            }
        } else if in_rows {
            // The closing dashed line right after the last data row.
            break;
        }
    }

    if symbols.is_empty() {
        return Err("g98.out's geometry block had no atom rows.".to_string());
    }
    Ok((symbols, coords))
}

/// One line of xTB's `$vibrational spectrum` block.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VibSpectrumLine {
    pub mode: usize,
    pub wavenumber_cm1: f64,
    pub ir_intensity_km_mol: f64,
}

/// Parses xTB's `vibspectrum` file. Each data row is either
/// `mode  wavenumber  intensity  selection` (translations/rotations, no
/// symmetry label) or `mode  symmetry  wavenumber  intensity  selection`
/// (vibrations) -- confirmed against the real file, both shapes present.
/// Read from the right rather than assuming a fixed column count: the last
/// three tokens are always selection, intensity, wavenumber regardless of
/// whether the symmetry label is present.
pub fn parse_vibspectrum(text: &str) -> Result<Vec<VibSpectrumLine>, String> {
    let mut out = Vec::new();
    let mut in_block = false;
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.starts_with("$vibrational spectrum") {
            in_block = true;
            continue;
        }
        if !in_block {
            continue;
        }
        if trimmed.starts_with("$end") {
            break;
        }
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        if tokens.len() < 4 {
            continue;
        }
        let n = tokens.len();
        let mode: usize = tokens[0]
            .parse()
            .map_err(|_| format!("Bad mode index in vibspectrum line: {trimmed:?}"))?;
        let wavenumber_cm1: f64 = tokens[n - 3]
            .parse()
            .map_err(|_| format!("Bad wavenumber in vibspectrum line: {trimmed:?}"))?;
        let ir_intensity_km_mol: f64 = tokens[n - 2]
            .parse()
            .map_err(|_| format!("Bad IR intensity in vibspectrum line: {trimmed:?}"))?;
        out.push(VibSpectrumLine {
            mode,
            wavenumber_cm1,
            ir_intensity_km_mol,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join(name);
        std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("cannot read fixture {}: {e}", path.display()))
    }

    #[test]
    fn parses_the_real_h2o2_hessian() {
        let text = fixture("h2o2_hessian");
        let h = parse_turbomole_hessian(&text, 4).unwrap();
        assert_eq!(h.dim(), (12, 12));
        // Spot-check the first row against the fixture's own first values.
        assert!((h[(0, 0)] - 0.3592143426).abs() < 1e-9);
        assert!((h[(0, 1)] - -0.2069309954).abs() < 1e-9);
    }

    #[test]
    fn the_real_hessian_is_symmetric() {
        // A physical Hessian must be (numerically near-)symmetric; this is
        // a property of the real data, not something the parser enforces.
        let text = fixture("h2o2_hessian");
        let h = parse_turbomole_hessian(&text, 4).unwrap();
        for i in 0..12 {
            for j in 0..12 {
                assert!(
                    (h[(i, j)] - h[(j, i)]).abs() < 1e-6,
                    "asymmetric at ({i},{j}): {} vs {}",
                    h[(i, j)],
                    h[(j, i)]
                );
            }
        }
    }

    #[test]
    fn rejects_a_wrong_atom_count() {
        let text = fixture("h2o2_hessian");
        assert!(parse_turbomole_hessian(&text, 3).is_err());
    }

    #[test]
    fn rejects_a_missing_header() {
        assert!(parse_turbomole_hessian("1.0 2.0 3.0 4.0", 1).is_err());
    }

    #[test]
    fn parses_the_real_h2o2_vibspectrum() {
        let text = fixture("h2o2_vibspectrum");
        let lines = parse_vibspectrum(&text).unwrap();
        assert_eq!(lines.len(), 12);
        // The six near-zero translation/rotation modes (no symmetry label).
        for l in &lines[0..6] {
            assert!(l.wavenumber_cm1.abs() < 0.01, "{l:?}");
        }
        // The six genuine vibrations, exact values from the captured file.
        let expected = [
            (259.56, 413.48255),
            (1108.86, 0.00066),
            (1158.96, 213.77306),
            (1355.45, 0.00005),
            (3531.45, 0.14239),
            (3534.68, 101.16131),
        ];
        for (line, &(w, i)) in lines[6..].iter().zip(&expected) {
            assert!((line.wavenumber_cm1 - w).abs() < 1e-6, "{line:?}");
            assert!((line.ir_intensity_km_mol - i).abs() < 1e-6, "{line:?}");
        }
    }

    #[test]
    fn vibspectrum_of_empty_text_is_empty_not_an_error() {
        assert_eq!(parse_vibspectrum("").unwrap(), vec![]);
    }

    #[test]
    fn parses_the_real_h2o2_g98_geometry() {
        let (symbols, coords) = parse_g98_geometry(&fixture("h2o2_g98.out")).unwrap();
        assert_eq!(symbols, vec!["O", "O", "H", "H"]);
        assert_eq!(coords.len(), 12);
        assert!((coords[0] - 0.0).abs() < 1e-6);
        assert!((coords[4] - 0.000002).abs() < 1e-5);
    }

    /// The exact regression: a larger, Cs-symmetry molecule for which xTB
    /// does not write `xtbhess.xyz` at all (confirmed by running
    /// `xtb --hess` on it directly), but does write `g98.out`.
    #[test]
    fn parses_the_real_phenyl_och3_g98_geometry_even_though_xtbhess_is_absent() {
        let (symbols, coords) = parse_g98_geometry(&fixture("phenyl_och3_g98.out")).unwrap();
        assert_eq!(symbols.len(), 16);
        assert_eq!(coords.len(), 48);
        assert_eq!(symbols[0], "C");
        assert_eq!(symbols[11], "O");
    }

    #[test]
    fn g98_geometry_of_text_without_the_block_is_an_error() {
        assert!(parse_g98_geometry("nothing relevant here").is_err());
    }

    #[test]
    fn parses_the_real_h2o2_gradient() {
        let text = fixture("h2o2_gradient");
        let g = parse_turbomole_gradient(&text, 4).unwrap();
        assert_eq!(g.len(), 12);
        assert!((g[0] - -7.4602070474193E-03).abs() < 1e-12);
        assert!((g[11] - -4.9994414736491E-03).abs() < 1e-12);
    }

    #[test]
    fn rejects_a_gradient_with_the_wrong_atom_count() {
        let text = fixture("h2o2_gradient");
        assert!(parse_turbomole_gradient(&text, 3).is_err());
    }

    /// The decisive test for this whole pipeline: feed the real, captured
    /// Hessian through Beavyr's own imported normal-mode analysis (parser,
    /// mass table, the pure-Rust Jacobi eigensolver shim, and the Eckart
    /// vib/rot projection all at once) and check the resulting frequencies
    /// against xTB's own independently-computed `vibspectrum` -- not merely
    /// that the code runs, but that it agrees with an outside implementation.
    #[test]
    fn normal_mode_analysis_reproduces_xtbs_own_frequencies() {
        const ANGSTROM_TO_BOHR: f64 = 1.0 / 0.529_177_210_903;

        // The Eckart projector needs the SAME geometry the Hessian was
        // computed at. `--ohess` optimizes before computing the Hessian, so
        // this is xtb's own post-optimization geometry, not the original
        // `examples/h2o2.xyz` -- using the wrong one was the first version
        // of this test's own bug, caught by exactly the cross-check this
        // test exists to do: several frequencies came out wrong (up to 12%)
        // while others matched almost exactly, which is the signature of an
        // Eckart projector built at the wrong geometry only partially
        // failing to remove rotational contamination.
        let xyz = fixture("h2o2_hessian_geometry.xyz");
        let (_n, atoms, coords_angstrom) = crate::molecule::parse_xyz_angstrom(&xyz);
        let coords_bohr: Vec<f64> = coords_angstrom
            .into_iter()
            .flatten()
            .map(|c| c * ANGSTROM_TO_BOHR)
            .collect();
        let mass: Vec<f64> = atoms
            .iter()
            .map(|sym| crate::molecule::atomic_mass_amu(sym).expect("known element"))
            .collect();

        let hessian = parse_turbomole_hessian(&fixture("h2o2_hessian"), atoms.len()).unwrap();
        let result =
            crate::normalmode::normal_modes(&mass, &hessian, false, &coords_bohr).unwrap();

        // Six near-zero (translation/rotation) modes, six genuine vibrations.
        assert_eq!(
            result.zero_indices.len(),
            6,
            "H2O2 is non-linear: exactly 6 low modes expected, got {:?}",
            result.zero_indices
        );
        assert_eq!(result.positive_indices.len(), 6, "{:?}", result.positive_indices);
        assert!(result.negative_indices.is_empty(), "a converged minimum has no imaginary modes");

        let mut computed_cm1: Vec<f64> = result
            .positive_indices
            .iter()
            .map(|&i| result.frequencies_au[i])
            .collect();
        computed_cm1 = crate::normalmode::thermofuncs::freqs_au_to_cm1(&computed_cm1);
        computed_cm1.sort_by(|a, b| a.total_cmp(b));

        let vib_lines = parse_vibspectrum(&fixture("h2o2_vibspectrum")).unwrap();
        let mut expected_cm1: Vec<f64> = vib_lines
            .iter()
            .map(|l| l.wavenumber_cm1)
            .filter(|w| w.abs() > 1.0) // drop the six ~0 translation/rotation rows
            .collect();
        expected_cm1.sort_by(|a, b| a.total_cmp(b));

        assert_eq!(computed_cm1.len(), expected_cm1.len());
        for (computed, expected) in computed_cm1.iter().zip(&expected_cm1) {
            let rel_err = (computed - expected).abs() / expected.abs();
            assert!(
                rel_err < 0.02,
                "computed {computed:.2} cm^-1 vs xTB's own {expected:.2} cm^-1 \
                 ({:.1}% off)",
                rel_err * 100.0
            );
        }
    }
}
