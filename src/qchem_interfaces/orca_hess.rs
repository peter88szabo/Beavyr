//! Reads an ORCA `.hess` file.
//!
//! The file is a list of `$name` sections ending at `$end`. We take the
//! Hessian, the geometry and masses, and ORCA's own frequencies and IR
//! intensities; the normal modes ORCA also stores are deliberately ignored,
//! because the point of loading the file is to run its Hessian through our own
//! projection and eigensolver.
//!
//! The `$hessian` block layout matches the reader in Smite's
//! `qchem_interfaces/orcarun.py`: five columns per block, each block preceded
//! by a row of column indices, each row prefixed by its own index.

use ndarray::Array2;

/// Everything read out of a `.hess`.
///
/// The whole file's contents are exposed, including fields the analysis does
/// not currently consult: they cost nothing to keep, and dropping them would
/// mean re-parsing the file the first time something wants them.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct OrcaHessian {
    pub natoms: usize,
    pub atoms: Vec<String>,
    /// As recorded in the file, so an isotope substitution is respected rather
    /// than overwritten by the standard atomic weight.
    pub masses_amu: Vec<f64>,
    /// Flattened `3N`, in bohr -- the unit `$atoms` is written in.
    pub coords_bohr: Vec<f64>,
    /// `3N x 3N`, in Eh/bohr².
    pub hessian: Array2<f64>,
    /// ORCA's own frequencies, in the order it printed them: the six projected
    /// translations/rotations first as exact zeros, then any imaginary mode,
    /// then the ascending positives. Kept for cross-checking, not for display.
    pub frequencies_cm1: Vec<f64>,
    /// `(frequency, intensity)` pairs sorted ascending by frequency, which is
    /// the order our own eigensolver returns modes in. See
    /// `ir_intensities_in_eigenvalue_order`.
    pub ir_spectrum: Vec<(f64, f64)>,
    pub multiplicity: Option<i32>,
    pub energy_hartree: Option<f64>,
    pub frequency_scale_factor: Option<f64>,
}

impl OrcaHessian {
    /// IR intensities in the order our eigensolver produces modes: ascending
    /// by eigenvalue, so an imaginary mode comes first.
    ///
    /// ORCA prints the six projected translations/rotations first and puts an
    /// imaginary mode seventh, so pairing by position -- which is valid for
    /// xTB's `vibspectrum` -- would put the imaginary mode's intensity on a
    /// translation. Sorting by frequency reproduces our order exactly, since
    /// "most negative first, then the zeros, then ascending" is what both
    /// orderings agree on once ORCA's list is sorted.
    pub fn ir_intensities_in_eigenvalue_order(&self) -> Vec<f64> {
        self.ir_spectrum.iter().map(|&(_, i)| i).collect()
    }

    /// ORCA's own frequencies in that same ascending order, for comparing
    /// against the ones we compute.
    #[allow(dead_code)]
    pub fn frequencies_in_eigenvalue_order(&self) -> Vec<f64> {
        let mut sorted = self.frequencies_cm1.clone();
        sorted.sort_by(f64::total_cmp);
        sorted
    }
}

/// Parses a `.hess` file.
///
/// `$hessian` and `$atoms` are required; the frequency and IR sections are
/// not, since the analysis recomputes the former and can proceed without the
/// latter.
pub fn parse_orca_hess(text: &str) -> Result<OrcaHessian, String> {
    let lines: Vec<&str> = text.lines().collect();

    let hessian = {
        let start = find_section(&lines, "$hessian")
            .ok_or_else(|| "no $hessian section in this file".to_string())?;
        parse_block_matrix(&lines, start)?
    };
    let n = hessian.nrows();

    let atoms_start = find_section(&lines, "$atoms")
        .ok_or_else(|| "no $atoms section in this file".to_string())?;
    let (atoms, masses_amu, coords_bohr) = parse_atoms(&lines, atoms_start)?;
    let natoms = atoms.len();

    if 3 * natoms != n {
        return Err(format!(
            "$atoms lists {natoms} atoms ({} coordinates) but $hessian is {n}x{n}",
            3 * natoms
        ));
    }

    let frequencies_cm1 = match find_section(&lines, "$vibrational_frequencies") {
        Some(start) => parse_indexed_values(&lines, start, "$vibrational_frequencies")?,
        None => Vec::new(),
    };
    let ir_spectrum = match find_section(&lines, "$ir_spectrum") {
        Some(start) => parse_ir_spectrum(&lines, start)?,
        None => Vec::new(),
    };

    Ok(OrcaHessian {
        natoms,
        atoms,
        masses_amu,
        coords_bohr,
        hessian,
        frequencies_cm1,
        ir_spectrum,
        multiplicity: find_section(&lines, "$multiplicity")
            .and_then(|i| parse_scalar(&lines, i))
            .map(|v| v as i32),
        energy_hartree: find_section(&lines, "$act_energy").and_then(|i| parse_scalar(&lines, i)),
        frequency_scale_factor: find_section(&lines, "$frequency_scale_factor")
            .and_then(|i| parse_scalar(&lines, i)),
    })
}

/// The index of the line holding `$name`, matched exactly so `$hessian` does
/// not also find a section whose name merely starts the same way.
fn find_section(lines: &[&str], name: &str) -> Option<usize> {
    lines.iter().position(|l| l.trim() == name)
}

/// The first non-blank line after `start`, and its index.
fn next_content(lines: &[&str], start: usize) -> Option<(usize, String)> {
    lines[start + 1..]
        .iter()
        .enumerate()
        .find(|(_, l)| !l.trim().is_empty())
        .map(|(offset, l)| (start + 1 + offset, l.trim().to_string()))
}

/// A section holding a single number on its own line.
fn parse_scalar(lines: &[&str], start: usize) -> Option<f64> {
    next_content(lines, start)?.1.parse().ok()
}

/// A square matrix written in blocks of five columns:
///
/// ```text
/// $hessian
/// 57
///              0        1        2        3        4
///     0    3.5E-01  7.1E-03  ...
/// ```
///
/// Each block repeats the column-index header, then all `n` rows.
fn parse_block_matrix(lines: &[&str], start: usize) -> Result<Array2<f64>, String> {
    const BLOCK: usize = 5;
    let (dim_line, dim_text) = next_content(lines, start)
        .ok_or_else(|| "$hessian ends before its dimension".to_string())?;
    let n: usize = dim_text
        .parse()
        .map_err(|_| format!("$hessian dimension is not a number: {dim_text:?}"))?;
    if n == 0 {
        return Err("$hessian has zero dimension".to_string());
    }

    let mut matrix = Array2::<f64>::zeros((n, n));
    let mut cursor = dim_line + 1;
    let mut first_col = 0;
    while first_col < n {
        let ncols = BLOCK.min(n - first_col);
        // Skip the block's column-index header.
        let (header, _) = next_content(lines, cursor - 1)
            .ok_or_else(|| format!("$hessian ends before the block at column {first_col}"))?;
        cursor = header + 1;

        for row in 0..n {
            let (index, content) = next_content(lines, cursor - 1).ok_or_else(|| {
                format!("$hessian ends inside the block at column {first_col}, row {row}")
            })?;
            cursor = index + 1;
            let values: Vec<&str> = content.split_whitespace().collect();
            // The leading token is the row's own index.
            if values.len() < ncols + 1 {
                return Err(format!(
                    "$hessian row {row} of the block at column {first_col} has \
                     {} values, expected {ncols}",
                    values.len().saturating_sub(1)
                ));
            }
            for (offset, token) in values[1..=ncols].iter().enumerate() {
                matrix[[row, first_col + offset]] = token.parse::<f64>().map_err(|_| {
                    format!("$hessian row {row} holds a non-numeric value {token:?}")
                })?;
            }
        }
        first_col += ncols;
    }
    Ok(matrix)
}

/// ```text
/// $atoms
/// 19
///  H      1.00800     -2.073404472303     1.655160222253     2.146818172856
/// ```
///
/// Coordinates are in bohr, which is what the rest of the pipeline wants.
type AtomsSection = (Vec<String>, Vec<f64>, Vec<f64>);

fn parse_atoms(lines: &[&str], start: usize) -> Result<AtomsSection, String> {
    let (count_line, count_text) =
        next_content(lines, start).ok_or_else(|| "$atoms ends before its count".to_string())?;
    let natoms: usize = count_text
        .parse()
        .map_err(|_| format!("$atoms count is not a number: {count_text:?}"))?;

    let mut atoms = Vec::with_capacity(natoms);
    let mut masses = Vec::with_capacity(natoms);
    let mut coords = Vec::with_capacity(3 * natoms);
    let mut cursor = count_line + 1;
    for i in 0..natoms {
        let (index, content) = next_content(lines, cursor - 1)
            .ok_or_else(|| format!("$atoms ends at atom {i} of {natoms}"))?;
        cursor = index + 1;
        let fields: Vec<&str> = content.split_whitespace().collect();
        if fields.len() < 5 {
            return Err(format!(
                "$atoms line for atom {i} needs symbol, mass and three coordinates: {content:?}"
            ));
        }
        atoms.push(fields[0].to_string());
        let parse = |token: &str, what: &str| -> Result<f64, String> {
            token
                .parse::<f64>()
                .map_err(|_| format!("$atoms atom {i}: {what} is not a number ({token:?})"))
        };
        masses.push(parse(fields[1], "mass")?);
        coords.push(parse(fields[2], "x")?);
        coords.push(parse(fields[3], "y")?);
        coords.push(parse(fields[4], "z")?);
    }
    Ok((atoms, masses, coords))
}

/// A section of `index  value` rows under a count, such as
/// `$vibrational_frequencies`.
fn parse_indexed_values(
    lines: &[&str],
    start: usize,
    name: &str,
) -> Result<Vec<f64>, String> {
    let (count_line, count_text) =
        next_content(lines, start).ok_or_else(|| format!("{name} ends before its count"))?;
    let count: usize = count_text
        .parse()
        .map_err(|_| format!("{name} count is not a number: {count_text:?}"))?;

    let mut values = Vec::with_capacity(count);
    let mut cursor = count_line + 1;
    for i in 0..count {
        let (index, content) =
            next_content(lines, cursor - 1).ok_or_else(|| format!("{name} ends at row {i}"))?;
        cursor = index + 1;
        let fields: Vec<&str> = content.split_whitespace().collect();
        let value = fields
            .get(1)
            .ok_or_else(|| format!("{name} row {i} has no value: {content:?}"))?;
        values.push(
            value
                .parse::<f64>()
                .map_err(|_| format!("{name} row {i} is not a number: {value:?}"))?,
        );
    }
    Ok(values)
}

/// ```text
/// $ir_spectrum
/// 57
///    106.90   0.00073996   3.73945848    -0.041449     0.013271     0.016306
/// ```
///
/// Columns are frequency, epsilon, intensity in km/mol, then the transition
/// dipole. Returned sorted ascending by frequency -- see
/// `OrcaHessian::ir_intensities_in_eigenvalue_order` for why that matters.
fn parse_ir_spectrum(lines: &[&str], start: usize) -> Result<Vec<(f64, f64)>, String> {
    /// Which column holds the intensity in km/mol, in the six-column layout
    /// ORCA 5 and 6 both write.
    const INTENSITY_COLUMN: usize = 2;

    let (count_line, count_text) =
        next_content(lines, start).ok_or_else(|| "$ir_spectrum ends before its count".to_string())?;
    let count: usize = count_text
        .parse()
        .map_err(|_| format!("$ir_spectrum count is not a number: {count_text:?}"))?;

    let mut rows = Vec::with_capacity(count);
    let mut cursor = count_line + 1;
    for i in 0..count {
        let (index, content) =
            next_content(lines, cursor - 1).ok_or_else(|| format!("$ir_spectrum ends at row {i}"))?;
        cursor = index + 1;
        let fields: Vec<&str> = content.split_whitespace().collect();
        if fields.len() <= INTENSITY_COLUMN {
            return Err(format!(
                "$ir_spectrum row {i} has {} columns, expected at least {}: {content:?}",
                fields.len(),
                INTENSITY_COLUMN + 1
            ));
        }
        let number = |col: usize| -> Result<f64, String> {
            fields[col]
                .parse::<f64>()
                .map_err(|_| format!("$ir_spectrum row {i} column {col} is not a number"))
        };
        rows.push((number(0)?, number(INTENSITY_COLUMN)?));
    }
    rows.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;

    fn example() -> String {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("TS_Gamma-3-6_ZZ-S_11_Compound_1.hess");
        std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
    }

    fn parsed() -> OrcaHessian {
        parse_orca_hess(&example()).expect("the example must parse")
    }

    #[test]
    fn the_geometry_and_masses_come_from_the_atoms_section() {
        let h = parsed();
        assert_eq!(h.natoms, 19);
        assert_eq!(h.atoms[0], "H");
        assert_eq!(h.atoms[1], "C");
        assert_eq!(h.atoms[5], "O");
        assert!((h.masses_amu[0] - 1.008).abs() < 1e-6);
        assert!((h.masses_amu[1] - 12.011).abs() < 1e-6);
        assert!((h.masses_amu[5] - 15.999).abs() < 1e-6);
        assert_eq!(h.coords_bohr.len(), 57);
        assert!((h.coords_bohr[0] + 2.073_404_472_303).abs() < 1e-12);
        assert!((h.coords_bohr[2] - 2.146_818_172_856).abs() < 1e-12);
    }

    /// Coordinates are in bohr, not angstrom -- reading them as angstrom would
    /// give a 2.05 A "C-H bond" and frequencies off by a factor of ~3.5.
    #[test]
    fn coordinates_are_in_bohr() {
        let h = parsed();
        let d: f64 = (0..3)
            .map(|k| (h.coords_bohr[k] - h.coords_bohr[3 + k]).powi(2))
            .sum::<f64>()
            .sqrt();
        let angstrom = d * BOHR_TO_ANGSTROM;
        assert!(
            (angstrom - 1.084).abs() < 0.01,
            "H-C came out {angstrom:.3} A, which is not a C-H bond"
        );
    }

    #[test]
    fn the_hessian_is_square_and_symmetric() {
        let h = parsed();
        assert_eq!(h.hessian.dim(), (57, 57));
        assert!((h.hessian[[0, 0]] - 3.5457342187E-01).abs() < 1e-12);
        assert!((h.hessian[[0, 3]] + 3.4036544410E-01).abs() < 1e-12);
        for i in 0..57 {
            for j in 0..57 {
                assert!(
                    (h.hessian[[i, j]] - h.hessian[[j, i]]).abs() < 1e-10,
                    "asymmetric at ({i},{j})"
                );
            }
        }
    }

    /// The last column block is a partial one (57 = 11 x 5 + 2), which is
    /// where an off-by-one in the block reader would show up.
    #[test]
    fn the_final_partial_column_block_is_read() {
        let h = parsed();
        assert_ne!(h.hessian[[0, 56]], 0.0);
        assert_ne!(h.hessian[[56, 56]], 0.0);
        assert!((h.hessian[[56, 0]] - h.hessian[[0, 56]]).abs() < 1e-10);
    }

    #[test]
    fn orcas_own_frequencies_are_read_in_file_order() {
        let h = parsed();
        assert_eq!(h.frequencies_cm1.len(), 57);
        // Six projected zeros, then the imaginary mode seventh.
        assert!(h.frequencies_cm1[..6].iter().all(|&f| f == 0.0));
        assert!((h.frequencies_cm1[6] + 602.605_246_697_391_6).abs() < 1e-9);
        assert!((h.frequencies_cm1[7] - 106.904_927_010_273_45).abs() < 1e-9);
        assert!((h.frequencies_cm1[56] - 3781.37).abs() < 0.01);
    }

    /// The ordering trap: ORCA puts the zeros first and the imaginary mode
    /// seventh, our eigensolver puts the imaginary mode first. Sorting is what
    /// reconciles them.
    #[test]
    fn frequencies_and_intensities_are_reordered_to_match_our_eigensolver() {
        let h = parsed();
        let freqs = h.frequencies_in_eigenvalue_order();
        assert!((freqs[0] + 602.605_246_697_391_6).abs() < 1e-9, "imaginary first");
        assert!(freqs[1..7].iter().all(|&f| f == 0.0), "then the six zeros");
        assert!(freqs[7] > 0.0, "then the positives");
        assert!(freqs.windows(2).all(|w| w[0] <= w[1]), "ascending");

        let intensities = h.ir_intensities_in_eigenvalue_order();
        assert_eq!(intensities.len(), 57);
        // A transition state's imaginary mode has no IR intensity listed, and
        // the projected modes have none either.
        assert_eq!(intensities[0], 0.0);
        assert!(intensities[1..7].iter().all(|&i| i == 0.0));
        // The highest mode is the O-H stretch at 3781.37 cm^-1.
        assert!((intensities[56] - 133.380_688_60).abs() < 1e-6);
    }

    #[test]
    fn the_ir_spectrum_pairs_each_frequency_with_its_intensity() {
        let h = parsed();
        assert_eq!(h.ir_spectrum.len(), 57);
        let (freq, intensity) = h.ir_spectrum[56];
        assert!((freq - 3781.37).abs() < 0.01);
        assert!((intensity - 133.380_688_60).abs() < 1e-6);
    }

    #[test]
    fn the_scalar_sections_are_read() {
        let h = parsed();
        assert_eq!(h.multiplicity, Some(1));
        assert_eq!(h.energy_hartree, Some(0.0));
        assert!((h.frequency_scale_factor.unwrap() - 0.974).abs() < 1e-9);
    }

    // ---- synthetic cases ----

    /// Two atoms, six coordinates: one full block of five columns plus a
    /// partial block of one.
    const MINIMAL: &str = "\
$orca_hessian_file

$hessian
6
          0        1        2        3        4
    0    1.0      0.0      0.0     -1.0      0.0
    1    0.0      2.0      0.0      0.0     -2.0
    2    0.0      0.0      3.0      0.0      0.0
    3   -1.0      0.0      0.0      1.0      0.0
    4    0.0     -2.0      0.0      0.0      2.0
    5    0.0      0.0     -3.0      0.0      0.0
          5
    0    0.0
    1    0.0
    2   -3.0
    3    0.0
    4    0.0
    5    3.0

$atoms
2
 H      1.00800      0.000000      0.000000      0.000000
 H      1.00800      0.000000      0.000000      1.400000

$end
";

    #[test]
    fn a_two_block_matrix_is_stitched_together_correctly() {
        let h = parse_orca_hess(MINIMAL).unwrap();
        assert_eq!(h.hessian.dim(), (6, 6));
        assert_eq!(h.hessian[[0, 0]], 1.0);
        assert_eq!(h.hessian[[0, 3]], -1.0);
        // Column 5 comes from the second, partial block.
        assert_eq!(h.hessian[[2, 5]], -3.0);
        assert_eq!(h.hessian[[5, 5]], 3.0);
    }

    /// Frequencies and IR are optional: the analysis recomputes the former and
    /// can run without the latter.
    #[test]
    fn a_file_without_frequencies_or_ir_still_parses() {
        let h = parse_orca_hess(MINIMAL).unwrap();
        assert!(h.frequencies_cm1.is_empty());
        assert!(h.ir_spectrum.is_empty());
        assert!(h.ir_intensities_in_eigenvalue_order().is_empty());
        assert_eq!(h.natoms, 2);
    }

    #[test]
    fn a_missing_hessian_names_the_section() {
        let err = parse_orca_hess("$orca_hessian_file\n\n$atoms\n0\n\n$end\n").unwrap_err();
        assert!(err.contains("$hessian"), "{err}");
    }

    #[test]
    fn a_missing_atoms_section_names_it() {
        let text = MINIMAL.replace("$atoms", "$other_atoms");
        let err = parse_orca_hess(&text).unwrap_err();
        assert!(err.contains("$atoms"), "{err}");
    }

    /// A Hessian that does not match the atom count is a corrupt pairing, not
    /// something to analyse anyway.
    #[test]
    fn a_hessian_that_does_not_match_the_atom_count_is_rejected() {
        let text = MINIMAL.replace(
            " H      1.00800      0.000000      0.000000      1.400000\n",
            "",
        );
        let text = text.replace("$atoms\n2\n", "$atoms\n1\n");
        let err = parse_orca_hess(&text).unwrap_err();
        assert!(err.contains("$hessian is 6x6"), "{err}");
    }

    #[test]
    fn a_truncated_block_is_reported_rather_than_silently_zero_filled() {
        let text = MINIMAL.replace("    5    0.0      0.0     -3.0      0.0      0.0\n", "");
        let err = parse_orca_hess(&text).unwrap_err();
        assert!(err.contains("$hessian"), "{err}");
    }

    #[test]
    fn a_non_numeric_hessian_entry_is_reported() {
        let text = MINIMAL.replace("    1.0      0.0      0.0     -1.0      0.0", "    x  0 0 0 0");
        let err = parse_orca_hess(&text).unwrap_err();
        assert!(err.contains("non-numeric"), "{err}");
    }
}
