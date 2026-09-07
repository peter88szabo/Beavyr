//! Driving Behemoth through its command line, and reading what it writes.
//!
//! Behemoth is run exactly as xTB is: an executable the user points at, given
//! an XYZ file in a scratch directory. Only the flags and the output formats
//! differ, and both live here so the panels stay program-agnostic.
//!
//! Verified against `behemoth --help` and real runs; the fixtures in
//! `tests/fixtures/behemoth_water_*` are captured output.

use std::path::Path;
use std::process::Command;

use ndarray::Array2;

/// The file Behemoth writes the optimized geometry to.
///
/// The spelling is Behemoth's own and has to be matched exactly -- looking for
/// the correctly spelled name would silently find nothing.
pub const OPTIMIZED_GEOMETRY_FILE: &str = "optimized_geoemtry.xyz";
/// The multi-frame XYZ it writes the optimization path to.
pub const TRAJECTORY_FILE: &str = "optimization_trajectory.xyz";

/// The electronic method to ask Behemoth for.
///
/// Only xTB is offered for now: the others need a basis set, and the DFT ones
/// a functional, which is a separate piece of UI rather than a hidden default.
pub const METHOD: &str = "xtb";

/// A geometry optimization, run in `run_dir`.
pub fn optimize_command(
    binary: &Path,
    run_dir: &Path,
    xyz_file: &str,
    charge: i32,
    multiplicity: i32,
) -> Command {
    let mut command = base_command(binary, run_dir, xyz_file, charge, multiplicity);
    command.arg("--opt");
    command
}

/// A finite-difference Hessian, printed to stdout.
///
/// `--hessian` rather than `--freq`: the matrix can be put through Beavyr's
/// own projection and eigensolver, which yields the normal modes an animation
/// needs and lets the thermochemistry follow the panel's temperature, cutoff
/// and symmetry number. `--freq` prints a finished frequency list with no
/// modes, so nothing could be animated from it.
pub fn hessian_command(
    binary: &Path,
    run_dir: &Path,
    xyz_file: &str,
    charge: i32,
    multiplicity: i32,
) -> Command {
    let mut command = base_command(binary, run_dir, xyz_file, charge, multiplicity);
    command.arg("--hessian");
    command
}

fn base_command(
    binary: &Path,
    run_dir: &Path,
    xyz_file: &str,
    charge: i32,
    multiplicity: i32,
) -> Command {
    let mut command = Command::new(binary);
    command
        .current_dir(run_dir)
        .arg("--xyz")
        .arg(xyz_file)
        .arg("--method")
        .arg(METHOD)
        .arg("--charge")
        .arg(charge.to_string())
        .arg("--multi")
        .arg(multiplicity.max(1).to_string());
    command
}

/// What Behemoth's optimization summary says.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OptimizationSummary {
    pub converged: bool,
    pub cycles: Option<usize>,
    pub final_energy_hartree: Option<f64>,
}

/// Reads the summary block Behemoth prints at the end of an optimization:
///
/// ```text
///   converged : true
///   cycles    : 3
///   E(final)  :      -5.7687749333 Eh
/// ```
pub fn parse_summary(stdout: &str) -> OptimizationSummary {
    let field = |name: &str| -> Option<String> {
        stdout.lines().find_map(|line| {
            let trimmed = line.trim();
            let rest = trimmed.strip_prefix(name)?;
            let value = rest.trim_start().strip_prefix(':')?;
            Some(value.trim().to_string())
        })
    };
    OptimizationSummary {
        converged: field("converged").as_deref() == Some("true"),
        cycles: field("cycles").and_then(|v| v.parse().ok()),
        final_energy_hartree: field("E(final)")
            .and_then(|v| v.split_whitespace().next()?.parse().ok()),
    }
}

/// The per-cycle energies of an optimization, from the trajectory's comment
/// lines: `cycle 1 E = -5.7686678546 Eh`.
///
/// Behemoth puts them there rather than in its log, which suits Beavyr -- the
/// energy plot reads them from the same file the trajectory player loads, so a
/// partial trajectory from a cancelled or failed run still plots.
pub fn parse_trajectory_energies(xyz_text: &str) -> Vec<f64> {
    xyz_text
        .lines()
        .filter_map(|line| {
            let trimmed = line.trim();
            if !trimmed.starts_with("cycle") {
                return None;
            }
            // `cycle <n> E = <value> Eh`
            let after_equals = trimmed.split('=').nth(1)?;
            after_equals.split_whitespace().next()?.parse().ok()
        })
        .collect()
}

/// Columns per printed block.
const BLOCK_COLUMNS: usize = 8;

/// The labels as Behemoth actually writes them, given the element symbols:
/// `O1x`, `H2y`. Used to check that the matrix belongs to the geometry that
/// was sent, rather than to whatever was in the directory before.
pub fn expected_labels(symbols: &[String]) -> Vec<String> {
    symbols
        .iter()
        .enumerate()
        .flat_map(|(atom, symbol)| {
            ["x", "y", "z"].into_iter().map(move |axis| {
                format!("{symbol}{}{axis}", atom + 1)
            })
        })
        .collect()
}

/// Reads the Hessian, checking its labels against the elements that were sent.
///
/// The labels carry the element symbols, so a matrix left over from another
/// molecule can be caught rather than silently analysed.
pub fn parse_hessian_for(stdout: &str, symbols: &[String]) -> Result<Array2<f64>, String> {
    let labels = expected_labels(symbols);
    let n = labels.len();
    if n == 0 {
        return Err("no atoms to read a Hessian for".to_string());
    }
    let mut matrix = Array2::<f64>::zeros((n, n));
    let mut seen = vec![false; n];
    let lines: Vec<&str> = stdout.lines().collect();
    let mut first_column = 0usize;
    let mut i = 0usize;

    while i < lines.len() && first_column < n {
        let fields: Vec<&str> = lines[i].split_whitespace().collect();
        let expected = &labels[first_column..n.min(first_column + BLOCK_COLUMNS)];
        if fields.is_empty() || fields.len() > BLOCK_COLUMNS || fields != expected {
            i += 1;
            continue;
        }
        let block_width = fields.len();
        let mut rows_read = 0usize;
        let mut j = i + 1;
        while j < lines.len() && rows_read < n {
            let row: Vec<&str> = lines[j].split_whitespace().collect();
            if row.is_empty() {
                j += 1;
                continue;
            }
            let Some(row_index) = labels.iter().position(|l| l == row[0]) else {
                break;
            };
            if row.len() != block_width + 1 {
                return Err(format!(
                    "Hessian row {} has {} values, expected {block_width}",
                    row[0],
                    row.len().saturating_sub(1)
                ));
            }
            for (offset, token) in row[1..].iter().enumerate() {
                matrix[[row_index, first_column + offset]] = token.parse::<f64>().map_err(|_| {
                    format!("Hessian row {} holds a non-numeric value {token:?}", row[0])
                })?;
            }
            seen[row_index] = true;
            rows_read += 1;
            j += 1;
        }
        if rows_read != n {
            return Err(format!(
                "the Hessian block at column {first_column} has {rows_read} rows, expected {n}"
            ));
        }
        first_column += block_width;
        i = j;
    }

    if first_column != n {
        return Err(format!(
            "Behemoth printed {first_column} of {n} expected Hessian columns -- \
             the output may belong to a different molecule"
        ));
    }
    if let Some(missing) = seen.iter().position(|s| !s) {
        return Err(format!("no Hessian row for {}", labels[missing]));
    }
    Ok(matrix)
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
            .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
    }

    fn water() -> Vec<String> {
        ["O", "H", "H"].iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn the_command_asks_for_what_it_should() {
        let command = optimize_command(
            Path::new("/bin/behemoth"),
            Path::new("/tmp/run"),
            "input.xyz",
            -1,
            2,
        );
        let args: Vec<String> = command
            .get_args()
            .map(|a| a.to_string_lossy().to_string())
            .collect();
        assert_eq!(
            args,
            vec!["--xyz", "input.xyz", "--method", "xtb", "--charge", "-1", "--multi", "2", "--opt"]
        );
        assert_eq!(command.get_current_dir(), Some(Path::new("/tmp/run")));
    }

    #[test]
    fn the_hessian_command_asks_for_the_matrix_not_the_frequency_list() {
        let command = hessian_command(
            Path::new("/bin/behemoth"),
            Path::new("/tmp/run"),
            "input.xyz",
            0,
            1,
        );
        let args: Vec<String> = command
            .get_args()
            .map(|a| a.to_string_lossy().to_string())
            .collect();
        assert!(args.contains(&"--hessian".to_string()));
        assert!(!args.contains(&"--freq".to_string()), "--freq has no modes");
    }

    /// A multiplicity of zero is not a thing; clamping keeps a bad panel value
    /// from reaching the command line.
    #[test]
    fn a_multiplicity_below_one_is_clamped() {
        let command =
            optimize_command(Path::new("b"), Path::new("."), "in.xyz", 0, 0);
        let args: Vec<String> = command
            .get_args()
            .map(|a| a.to_string_lossy().to_string())
            .collect();
        let multi = args.iter().position(|a| a == "--multi").unwrap();
        assert_eq!(args[multi + 1], "1");
    }

    // ---- summary ----

    #[test]
    fn the_summary_is_read_from_a_real_log() {
        let summary = parse_summary(&fixture("behemoth_water_opt.log"));
        assert!(summary.converged);
        assert_eq!(summary.cycles, Some(3));
        let energy = summary.final_energy_hartree.expect("a final energy");
        assert!((energy + 5.7687749333).abs() < 1e-9, "got {energy}");
    }

    #[test]
    fn a_log_without_a_summary_reports_nothing_rather_than_guessing() {
        let summary = parse_summary("Run resources\nCompute cores: 16\n");
        assert!(!summary.converged);
        assert!(summary.cycles.is_none());
        assert!(summary.final_energy_hartree.is_none());
    }

    // ---- trajectory energies ----

    #[test]
    fn the_per_cycle_energies_come_from_the_trajectory_comments() {
        let energies = parse_trajectory_energies(&fixture("behemoth_water_trajectory.xyz"));
        assert_eq!(energies.len(), 3, "three cycles: {energies:?}");
        assert!((energies[0] + 5.7686678546).abs() < 1e-9);
        assert!((energies[2] + 5.7687749333).abs() < 1e-6);
        // Monotonically downhill, as an optimization should be.
        assert!(energies.windows(2).all(|w| w[1] <= w[0] + 1e-9), "{energies:?}");
    }

    #[test]
    fn a_trajectory_with_no_energy_comments_yields_none() {
        let energies = parse_trajectory_energies("3\nwater\n O 0 0 0\n H 0 0 1\n H 0 1 0\n");
        assert!(energies.is_empty());
    }

    // ---- the Hessian ----

    #[test]
    fn the_water_hessian_is_read_whole_and_symmetric() {
        let h = parse_hessian_for(&fixture("behemoth_water_hessian.log"), &water())
            .expect("the fixture must parse");
        assert_eq!(h.dim(), (9, 9));
        // Corner values, straight from the printed matrix.
        assert!((h[[0, 0]] - 0.01353746).abs() < 1e-9);
        assert!((h[[1, 1]] - 0.59082867).abs() < 1e-9);
        assert!((h[[0, 3]] + 0.00568132).abs() < 1e-9);
        // The last column comes from the final, one-column block -- where an
        // off-by-one in the block walk would show up.
        assert!((h[[8, 8]] - 0.18980046).abs() < 1e-9);
        assert!((h[[1, 8]] + 0.22551677).abs() < 1e-9);
        for i in 0..9 {
            for j in 0..9 {
                assert!(
                    (h[[i, j]] - h[[j, i]]).abs() < 1e-8,
                    "asymmetric at ({i},{j})"
                );
            }
        }
    }

    /// The labels carry the element symbols, so a matrix belonging to another
    /// molecule is caught instead of being analysed as this one's.
    #[test]
    fn a_hessian_for_different_elements_is_rejected() {
        let wrong: Vec<String> = ["C", "H", "H"].iter().map(|s| s.to_string()).collect();
        let err = parse_hessian_for(&fixture("behemoth_water_hessian.log"), &wrong).unwrap_err();
        assert!(err.contains("different molecule") || err.contains("columns"), "{err}");
    }

    #[test]
    fn a_hessian_for_the_wrong_atom_count_is_rejected() {
        let short: Vec<String> = ["O", "H"].iter().map(|s| s.to_string()).collect();
        assert!(parse_hessian_for(&fixture("behemoth_water_hessian.log"), &short).is_err());
    }

    #[test]
    fn a_truncated_hessian_is_reported() {
        let text = fixture("behemoth_water_hessian.log");
        // Cut the final one-column block off.
        let cut = text.find("                         H3z").expect("the last block header");
        let err = parse_hessian_for(&text[..cut], &water()).unwrap_err();
        assert!(err.contains("columns"), "{err}");
    }

    #[test]
    fn no_atoms_is_an_error_not_an_empty_matrix() {
        assert!(parse_hessian_for("", &[]).is_err());
    }

    /// Runs the real binary, if it is where the user said it is, and checks
    /// the whole path: command, files written, energies and the Hessian.
    ///
    /// Ignored by default so the suite does not depend on a build of
    /// Behemoth being present. Run with
    /// `cargo test behemoth_binary -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn the_real_behemoth_binary_optimizes_and_hessians_water() {
        let binary = std::path::PathBuf::from(
            "/home/peter/Dropbox/Research_Leuven/Behemoth/target/release/behemoth",
        );
        if !binary.is_file() {
            eprintln!("no Behemoth build at {}, skipping", binary.display());
            return;
        }
        let dir = std::env::temp_dir().join(format!("beavyr_behemoth_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("input.xyz"),
            "3\n\nO 0.000000 0.000000 0.119262\n\
             H 0.000000 0.763239 -0.477047\n\
             H 0.000000 -0.763239 -0.477047\n",
        )
        .unwrap();

        // Optimization.
        let status = optimize_command(&binary, &dir, "input.xyz", 0, 1)
            .stdout(std::process::Stdio::piped())
            .output()
            .expect("Behemoth should run");
        assert!(status.status.success(), "Behemoth exited with an error");
        let log = String::from_utf8_lossy(&status.stdout).to_string();
        let summary = parse_summary(&log);
        assert!(summary.converged, "the optimization should converge");
        assert!(summary.final_energy_hartree.is_some());
        assert!(
            dir.join(OPTIMIZED_GEOMETRY_FILE).is_file(),
            "the optimized geometry file name must match Behemoth's spelling"
        );
        let trajectory = std::fs::read_to_string(dir.join(TRAJECTORY_FILE)).unwrap();
        let energies = parse_trajectory_energies(&trajectory);
        assert!(!energies.is_empty(), "the plot needs per-cycle energies");

        // Hessian.
        let hess = hessian_command(&binary, &dir, "input.xyz", 0, 1)
            .stdout(std::process::Stdio::piped())
            .output()
            .expect("Behemoth should run");
        assert!(hess.status.success());
        let atoms: Vec<String> = ["O", "H", "H"].iter().map(|s| s.to_string()).collect();
        let matrix =
            parse_hessian_for(&String::from_utf8_lossy(&hess.stdout), &atoms).expect("a Hessian");
        assert_eq!(matrix.dim(), (9, 9));
        for i in 0..9 {
            for j in 0..9 {
                assert!((matrix[[i, j]] - matrix[[j, i]]).abs() < 1e-6);
            }
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn the_labels_match_what_behemoth_writes() {
        let labels = expected_labels(&water());
        assert_eq!(labels.len(), 9);
        assert_eq!(labels[0], "O1x");
        assert_eq!(labels[4], "H2y");
        assert_eq!(labels[8], "H3z");
    }
}
