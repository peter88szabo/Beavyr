//! Reads frequencies, IR intensities and normal modes from a Gaussian output.
//!
//! A Gaussian `.log` from an ordinary `Freq` job does **not** contain the
//! force-constant matrix -- that needs `iop(7/33=1)` or a formatted
//! checkpoint. What it does contain is the finished result: the frequencies,
//! their IR intensities, and the normal coordinates themselves. So this is not
//! a Hessian reader: the modes are taken as given rather than recomputed, and
//! only the thermochemistry is calculated here.
//!
//! ## Normalisation
//!
//! Gaussian prints normal coordinates normalised in plain Cartesian space,
//! `sum |d_i|^2 = 1`, and prints the reduced mass `mu` alongside. Our own
//! modes come out of a mass-weighted eigenproblem and satisfy
//! `sum m_i |d_i|^2 = 1`. Since `mu = sum m_i |d_i|^2` under Gaussian's
//! convention, dividing an imported vector by `sqrt(mu)` converts it to ours
//! -- which is what makes the animation amplitude mean the same thing whatever
//! the modes came from.

use ndarray::Array2;

use crate::molecule::element_symbol;

/// What a Gaussian frequency job reports.
#[derive(Debug, Clone)]
pub struct GaussianFrequencies {
    pub atoms: Vec<String>,
    /// Flattened `3N`, in angstrom, from the frame the modes belong to.
    pub coords_angstrom: Vec<f64>,
    /// One per printed mode: `3N - 6` of them, already excluding the
    /// translations and rotations Gaussian projected out.
    pub frequencies_cm1: Vec<f64>,
    pub ir_intensities_km_mol: Vec<f64>,
    pub reduced_masses_amu: Vec<f64>,
    /// `3N x nmodes`, columns in our own normalisation (`sum m|d|^2 = 1`).
    pub modes: Array2<f64>,
}

/// Parses a Gaussian frequency output.
pub fn parse_gaussian_log(text: &str) -> Result<GaussianFrequencies, String> {
    let lines: Vec<&str> = text.lines().collect();

    let (atoms, coords_angstrom) = parse_geometry(&lines)?;
    let natoms = atoms.len();
    let blocks = parse_mode_blocks(&lines, natoms)?;
    if blocks.frequencies.is_empty() {
        return Err(
            "no frequencies in this file: expected 'Frequencies --' lines, as a \
             Gaussian Freq job writes"
                .to_string(),
        );
    }

    let nmodes = blocks.frequencies.len();
    let mut modes = Array2::<f64>::zeros((3 * natoms, nmodes));
    for (k, column) in blocks.displacements.iter().enumerate() {
        // Gaussian's convention is sum |d|^2 = 1 with mu = sum m|d|^2, so
        // dividing by sqrt(mu) lands on our sum m|d|^2 = 1. Falling back to
        // the vector's own norm keeps a file that somehow lacks reduced masses
        // usable rather than producing a mode of arbitrary scale.
        let mu = blocks.reduced_masses.get(k).copied().unwrap_or(0.0);
        let scale = if mu > 0.0 {
            mu.sqrt()
        } else {
            let norm: f64 = column.iter().map(|d| d * d).sum::<f64>().sqrt();
            if norm > 0.0 { norm } else { 1.0 }
        };
        for (row, value) in column.iter().enumerate() {
            modes[[row, k]] = value / scale;
        }
    }

    Ok(GaussianFrequencies {
        atoms,
        coords_angstrom,
        frequencies_cm1: blocks.frequencies,
        ir_intensities_km_mol: blocks.intensities,
        reduced_masses_amu: blocks.reduced_masses,
        modes,
    })
}

/// Whether this text looks like a Gaussian output at all.
pub fn is_gaussian_log(text: &str) -> bool {
    text.contains("Entering Gaussian System")
        || text.contains("Gaussian, Inc.")
        || text.contains(" Frequencies --")
}

/// The geometry the normal modes belong to.
///
/// Gaussian reports the frequency analysis in the standard orientation, so
/// that is preferred; a `nosymm` job prints only the input orientation, which
/// is then the same frame. The last of each is used, since a job may report
/// several along the way.
fn parse_geometry(lines: &[&str]) -> Result<(Vec<String>, Vec<f64>), String> {
    let start = last_index_of(lines, "Standard orientation:")
        .or_else(|| last_index_of(lines, "Input orientation:"))
        .ok_or_else(|| {
            "no geometry in this file: expected a 'Standard orientation:' or \
             'Input orientation:' table"
                .to_string()
        })?;

    let mut atoms = Vec::new();
    let mut coords = Vec::new();
    // The table opens with a rule, two header lines and another rule, then the
    // atom rows, then a closing rule.
    let mut seen_rule = 0;
    for line in &lines[start + 1..] {
        let trimmed = line.trim();
        if trimmed.starts_with("---") {
            seen_rule += 1;
            if seen_rule >= 3 {
                break;
            }
            continue;
        }
        if seen_rule < 2 {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        // center, atomic number, atomic type, x, y, z
        if fields.len() < 6 {
            continue;
        }
        let Ok(z) = fields[1].parse::<u32>() else {
            continue;
        };
        let symbol = element_symbol(z)
            .ok_or_else(|| format!("unknown atomic number {z} in the geometry"))?;
        let mut xyz = [0.0f64; 3];
        for (i, field) in fields[3..6].iter().enumerate() {
            xyz[i] = field
                .parse()
                .map_err(|_| format!("non-numeric coordinate {field:?} in the geometry"))?;
        }
        atoms.push(symbol.to_string());
        coords.extend_from_slice(&xyz);
    }
    if atoms.is_empty() {
        return Err("the geometry table lists no atoms".to_string());
    }
    Ok((atoms, coords))
}

fn last_index_of(lines: &[&str], needle: &str) -> Option<usize> {
    lines.iter().rposition(|l| l.trim() == needle)
}

/// Accumulated across every printed block.
struct ModeBlocks {
    frequencies: Vec<f64>,
    intensities: Vec<f64>,
    reduced_masses: Vec<f64>,
    /// One `3N` vector per mode, as printed (before renormalisation).
    displacements: Vec<Vec<f64>>,
}

/// Reads the frequency blocks, each covering up to three modes side by side:
///
/// ```text
///  Frequencies --  -2468.9062                78.6012               133.6459
///  Red. masses --      1.1301                 3.8666                 5.1644
///  IR Inten    --   6843.5507                 1.4167                 0.4868
///   Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z
///      1   6    -0.00   0.00  -0.00    -0.00   0.00  -0.01    -0.02   0.01  -0.04
/// ```
///
/// A `freq=hpmodes` job prints a five-column high-precision block first, whose
/// atom rows are headed `Coord Atom Element:` instead; keying on the `Atom  AN`
/// header means those are skipped and the standard blocks read.
fn parse_mode_blocks(lines: &[&str], natoms: usize) -> Result<ModeBlocks, String> {
    let mut blocks = ModeBlocks {
        frequencies: Vec::new(),
        intensities: Vec::new(),
        reduced_masses: Vec::new(),
        displacements: Vec::new(),
    };

    for (i, line) in lines.iter().enumerate() {
        let Some(rest) = line.trim().strip_prefix("Frequencies --") else {
            continue;
        };
        let freqs = parse_row(rest)?;
        let count = freqs.len();
        if count == 0 {
            continue;
        }

        // The rows between here and the atom table, in whatever order and
        // with whatever extras the job printed.
        let mut reduced = vec![0.0; count];
        let mut intensities = vec![0.0; count];
        let mut atom_header = None;
        for (offset, follow) in lines[i + 1..].iter().enumerate().take(12) {
            let trimmed = follow.trim();
            if let Some(rest) = trimmed.strip_prefix("Red. masses --") {
                reduced = fit(parse_row(rest)?, count);
            } else if let Some(rest) = trimmed.strip_prefix("IR Inten    --") {
                intensities = fit(parse_row(rest)?, count);
            } else if trimmed.starts_with("Atom  AN") {
                atom_header = Some(i + 1 + offset);
                break;
            }
        }
        let Some(header) = atom_header else {
            return Err(format!(
                "the block at line {} has no 'Atom  AN' displacement table",
                i + 1
            ));
        };

        let mut columns = vec![vec![0.0f64; 3 * natoms]; count];
        for atom in 0..natoms {
            let row = lines.get(header + 1 + atom).ok_or_else(|| {
                format!("the block at line {} ends at atom {atom} of {natoms}", i + 1)
            })?;
            let fields: Vec<&str> = row.split_whitespace().collect();
            // center, atomic number, then three components per mode.
            let expected = 2 + 3 * count;
            if fields.len() < expected {
                return Err(format!(
                    "atom row {} of the block at line {} has {} fields, expected {expected}",
                    atom + 1,
                    i + 1,
                    fields.len()
                ));
            }
            for (k, column) in columns.iter_mut().enumerate() {
                for axis in 0..3 {
                    let field = fields[2 + 3 * k + axis];
                    column[3 * atom + axis] = field.parse().map_err(|_| {
                        format!("non-numeric displacement {field:?} at atom {}", atom + 1)
                    })?;
                }
            }
        }

        blocks.frequencies.extend(freqs);
        blocks.reduced_masses.extend(reduced);
        blocks.intensities.extend(intensities);
        blocks.displacements.extend(columns);
    }
    Ok(blocks)
}

/// The numbers on one `... --` row.
fn parse_row(rest: &str) -> Result<Vec<f64>, String> {
    rest.split_whitespace()
        .map(|token| {
            token
                .parse::<f64>()
                .map_err(|_| format!("non-numeric value {token:?} in a frequency block"))
        })
        .collect()
}

/// Pads or truncates a row to the block's mode count, so a malformed extra
/// column cannot shift every later mode.
fn fit(mut values: Vec<f64>, count: usize) -> Vec<f64> {
    values.resize(count, 0.0);
    values
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::atomic_mass_amu;

    fn example() -> String {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("c164-ts1-2-1001.log");
        std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
    }

    fn parsed() -> GaussianFrequencies {
        parse_gaussian_log(&example()).expect("the example must parse")
    }

    #[test]
    fn it_is_recognised_as_a_gaussian_log() {
        assert!(is_gaussian_log(&example()));
        assert!(!is_gaussian_log("$orca_hessian_file\n$hessian\n"));
    }

    #[test]
    fn the_geometry_comes_from_the_standard_orientation() {
        let g = parsed();
        assert_eq!(g.atoms.len(), 19);
        assert_eq!(g.coords_angstrom.len(), 57);
        assert_eq!(g.atoms[0], "C");
        assert_eq!(g.atoms[3], "H");
        assert_eq!(g.atoms[7], "O");
        // Standard orientation, not input: -1.244180 rather than -1.244301.
        assert!((g.coords_angstrom[0] + 1.244180).abs() < 1e-9, "{}", g.coords_angstrom[0]);
        assert!((g.coords_angstrom[2] + 0.011685).abs() < 1e-9);
    }

    #[test]
    fn all_three_n_minus_six_modes_are_read() {
        let g = parsed();
        assert_eq!(g.frequencies_cm1.len(), 51, "3*19 - 6");
        assert_eq!(g.ir_intensities_km_mol.len(), 51);
        assert_eq!(g.reduced_masses_amu.len(), 51);
        assert_eq!(g.modes.dim(), (57, 51));
    }

    #[test]
    fn the_frequencies_and_intensities_match_the_file() {
        let g = parsed();
        assert!((g.frequencies_cm1[0] + 2468.9062).abs() < 1e-6, "the imaginary mode");
        assert!((g.frequencies_cm1[1] - 78.6012).abs() < 1e-6);
        assert!((g.frequencies_cm1[50] - 3867.0399).abs() < 1e-6);
        assert!((g.ir_intensities_km_mol[0] - 6843.5507).abs() < 1e-4);
        assert!((g.ir_intensities_km_mol[2] - 0.4868).abs() < 1e-4);
        assert!((g.reduced_masses_amu[0] - 1.1301).abs() < 1e-4);
    }

    /// Exactly one imaginary mode, and it is the first: this is a transition
    /// state, and Gaussian lists modes in ascending order.
    #[test]
    fn there_is_exactly_one_imaginary_mode() {
        let g = parsed();
        assert_eq!(g.frequencies_cm1.iter().filter(|f| **f < 0.0).count(), 1);
        assert!(g.frequencies_cm1[0] < 0.0);
        assert!(g.frequencies_cm1[1..].iter().all(|f| *f > 0.0));
        assert!(
            g.frequencies_cm1.windows(2).all(|w| w[0] <= w[1]),
            "ascending"
        );
    }

    /// The imported modes must arrive in our normalisation, `sum m|d|^2 = 1`,
    /// or the animation amplitude would mean something different for an
    /// imported mode than for one we computed.
    #[test]
    fn modes_are_converted_to_our_mass_weighted_normalisation() {
        let g = parsed();
        let masses: Vec<f64> = g
            .atoms
            .iter()
            .map(|s| atomic_mass_amu(s).unwrap())
            .collect();
        for k in 0..g.frequencies_cm1.len() {
            let norm: f64 = (0..g.atoms.len())
                .map(|atom| {
                    let m = masses[atom];
                    (0..3)
                        .map(|axis| {
                            let d = g.modes[[3 * atom + axis, k]];
                            m * d * d
                        })
                        .sum::<f64>()
                })
                .sum();
            // The printed vectors carry only two decimals, so the recovered
            // norm is approximate; a convention error would be off by a factor
            // of the reduced mass (1.1 to 12 here), not by a few percent.
            assert!(
                (norm - 1.0).abs() < 0.05,
                "mode {k} has sum m|d|^2 = {norm:.4}, expected 1"
            );
        }
    }

    /// The mode dominated by a single hydrogen is the H-transfer coordinate of
    /// this transition state -- a good check that displacements land on the
    /// right atoms rather than being shifted by a column.
    #[test]
    fn the_imaginary_mode_is_a_hydrogen_transfer() {
        let g = parsed();
        let mover = (0..19)
            .max_by(|&a, &b| {
                let d = |i: usize| {
                    (0..3)
                        .map(|axis| g.modes[[3 * i + axis, 0]].powi(2))
                        .sum::<f64>()
                };
                d(a).total_cmp(&d(b))
            })
            .unwrap();
        assert_eq!(g.atoms[mover], "H", "atom {} moves most", mover + 1);
        assert_eq!(mover, 17, "atom 18 in the file's numbering");
    }

    // ---- failure cases ----

    #[test]
    fn a_file_with_no_frequencies_is_an_error() {
        let text = " Standard orientation:\n ---\n a\n b\n ---\n      1  6  0  0.0 0.0 0.0\n ---\n";
        let err = parse_gaussian_log(text).unwrap_err();
        assert!(err.contains("Frequencies"), "{err}");
    }

    #[test]
    fn a_file_with_no_geometry_is_an_error() {
        let err = parse_gaussian_log(" Frequencies --   100.0\n").unwrap_err();
        assert!(err.contains("orientation"), "{err}");
    }

    #[test]
    fn a_truncated_displacement_table_is_reported() {
        let mut text = example();
        // Remove the last atom row of the first block.
        let marker = "  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z";
        let at = text.find(marker).unwrap();
        let tail_start = at + marker.len();
        let mut rows: Vec<&str> = text[tail_start..].lines().collect();
        rows.remove(1);
        text = format!("{}{}", &text[..tail_start], rows.join("\n"));
        let err = parse_gaussian_log(&text).unwrap_err();
        assert!(err.contains("fields") || err.contains("atom"), "{err}");
    }
}
