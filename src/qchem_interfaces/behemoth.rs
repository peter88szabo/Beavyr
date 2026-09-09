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

use super::method::{
    behemoth_basis_flag, behemoth_dispersion_flag, behemoth_functional_flag,
    behemoth_method_flag, MethodConfig,
};

/// The file Behemoth writes the optimized geometry to.
///
/// The spelling is Behemoth's own and has to be matched exactly -- looking for
/// the correctly spelled name would silently find nothing.
pub const OPTIMIZED_GEOMETRY_FILE: &str = "optimized_geoemtry.xyz";
/// The multi-frame XYZ it writes the optimization path to.
pub const TRAJECTORY_FILE: &str = "optimization_trajectory.xyz";


/// A geometry optimization, run in `run_dir`.
/// A single-point energy and nuclear gradient.
///
/// The command a conformer search needs: one call per gradient, with the geometry written into
/// `run_dir` beforehand. Behemoth is fast enough for this to be practical -- a TASI gradient on
/// water is nine milliseconds of compute, so the process launch dominates -- which is why driving
/// Beavyr's own search with Behemoth gradients is worth doing rather than delegating the search.
pub fn gradient_command(
    binary: &Path,
    run_dir: &Path,
    xyz_file: &str,
    charge: i32,
    multiplicity: i32,
    method: &MethodConfig,
) -> Command {
    let mut command = base_command(binary, run_dir, xyz_file, charge, multiplicity, method);
    command.arg("--gradient");
    command
}

/// Reads the energy and Cartesian gradient from a `--gradient` run's output.
///
/// The block looks like this, whichever method produced it:
///
/// ```text
/// Analytical TASI nuclear gradient
/// Units: Hartree/bohr
/// Energy:     -68.4051345926 Hartree
///      1     -0.0020104     -0.0021028     -0.0000000
///      2     -0.0026798      0.0040516      0.0000000
///      3      0.0046902     -0.0019488     -0.0000000
/// ```
///
/// Located by the "nuclear gradient" heading rather than by a method name, so a method whose
/// heading reads differently still parses. Energies are Hartree and gradients Hartree/bohr, which
/// is what the optimiser wants, so nothing is converted.
pub fn parse_gradient(stdout: &str, natoms: usize) -> Result<(f64, Vec<f64>), String> {
    let lines: Vec<&str> = stdout.lines().collect();
    let heading = lines
        .iter()
        .rposition(|line| line.to_ascii_lowercase().contains("nuclear gradient"))
        .ok_or_else(|| "Behemoth printed no nuclear gradient".to_string())?;

    let mut energy = None;
    let mut gradient = Vec::with_capacity(natoms * 3);
    for line in &lines[heading + 1..] {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix("Energy:") {
            energy = rest
                .split_whitespace()
                .next()
                .and_then(|token| token.parse::<f64>().ok());
            continue;
        }
        // An atom row: an index followed by three components. Anything else ends the block, but
        // only once rows have started -- the heading is followed by a units line first.
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() == 4 {
            if let (Ok(_index), Ok(x), Ok(y), Ok(z)) = (
                fields[0].parse::<usize>(),
                fields[1].parse::<f64>(),
                fields[2].parse::<f64>(),
                fields[3].parse::<f64>(),
            ) {
                gradient.extend_from_slice(&[x, y, z]);
                if gradient.len() == natoms * 3 {
                    break;
                }
                continue;
            }
        }
        if !gradient.is_empty() {
            break;
        }
    }

    let energy = energy.ok_or_else(|| "Behemoth printed no gradient energy".to_string())?;
    if gradient.len() != natoms * 3 {
        return Err(format!(
            "Behemoth printed {} gradient components for {natoms} atoms",
            gradient.len()
        ));
    }
    Ok((energy, gradient))
}

pub fn optimize_command(
    binary: &Path,
    run_dir: &Path,
    xyz_file: &str,
    charge: i32,
    multiplicity: i32,
    method: &MethodConfig,
) -> Command {
    let mut command = base_command(binary, run_dir, xyz_file, charge, multiplicity, method);
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
    method: &MethodConfig,
) -> Command {
    let mut command = base_command(binary, run_dir, xyz_file, charge, multiplicity, method);
    command.arg("--hessian");
    command
}

/// The flags every task shares: the geometry, the level of theory and the
/// resources.
///
/// Which of the level-of-theory flags apply is `method`'s decision, not this
/// function's -- a composite functional passes no basis set, HF-3c travels as
/// the method name rather than a functional, and an open shell is always the
/// unrestricted method. All of that lives in `super::method`, with its tests.
fn base_command(
    binary: &Path,
    run_dir: &Path,
    xyz_file: &str,
    charge: i32,
    multiplicity: i32,
    method: &MethodConfig,
) -> Command {
    let multiplicity = multiplicity.max(1);
    let mut command = Command::new(binary);
    command
        .current_dir(run_dir)
        .arg("--xyz")
        .arg(xyz_file)
        .arg("--method")
        .arg(behemoth_method_flag(method, multiplicity))
        .arg("--charge")
        .arg(charge.to_string())
        .arg("--multi")
        .arg(multiplicity.to_string());
    if let Some(functional) = behemoth_functional_flag(method) {
        command.arg("--functional").arg(functional);
    }
    if let Some(basis) = behemoth_basis_flag(method) {
        command.arg("--basis").arg(basis);
    }
    if let Some(dispersion) = behemoth_dispersion_flag(method) {
        command.arg("--dispersion").arg(dispersion);
    }
    command
        .arg("--memory")
        .arg(method.memory_mb.to_string())
        .arg("--nproc")
        .arg(method.nproc.max(1).to_string());
    command
}

/// What Behemoth's optimization summary says.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OptimizationSummary {
    pub converged: bool,
    pub cycles: Option<usize>,
    pub final_energy_hartree: Option<f64>,
    /// The RMS gradient it finished at, in Eh/bohr.
    pub rms_gradient: Option<f64>,
}

/// Reads the summary block Behemoth prints at the end of an optimization:
///
/// ```text
///   converged : true
///   cycles    : 5
///   E(final)  :      -5.7687748782 Eh
///   RMS(g)    :       0.0000601941 Eh/bohr
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
    let number = |name: &str| -> Option<f64> {
        field(name).and_then(|v| v.split_whitespace().next()?.parse().ok())
    };
    OptimizationSummary {
        converged: field("converged").as_deref() == Some("true"),
        cycles: field("cycles").and_then(|v| v.parse().ok()),
        final_energy_hartree: number("E(final)"),
        rms_gradient: number("RMS(g)"),
    }
}

/// One row of Behemoth's optimization convergence table.
///
/// The five convergence measures are exactly the ones it tests against its
/// thresholds, in the order it prints them, so a row can be compared straight
/// against `ConvergenceThresholds` column by column.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OptimizationStep {
    pub step: usize,
    pub energy_hartree: f64,
    /// Energy change from the previous step, in kJ/mol as Behemoth prints it.
    pub de_kj_mol: f64,
    pub max_step: f64,
    pub rms_step: f64,
    pub max_grad: f64,
    pub rms_grad: f64,
}

/// The convergence criteria the run was held to, in the same column order as
/// an `OptimizationStep`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConvergenceThresholds {
    pub de_kj_mol: f64,
    pub max_step: f64,
    pub rms_step: f64,
    pub max_grad: f64,
    pub rms_grad: f64,
}

impl ConvergenceThresholds {
    /// Whether each of a step's five measures is inside its threshold, in
    /// column order. `dE` is compared on magnitude, since it is signed and a
    /// descent prints it negative.
    pub fn met_by(&self, step: &OptimizationStep) -> [bool; 5] {
        [
            step.de_kj_mol.abs() <= self.de_kj_mol,
            step.max_step <= self.max_step,
            step.rms_step <= self.rms_step,
            step.max_grad <= self.max_grad,
            step.rms_grad <= self.rms_grad,
        ]
    }
}

/// The column headings, for a table that wants to label itself without
/// repeating the order in two places.
pub const STEP_COLUMNS: [&str; 5] =
    ["dE [kJ/mol]", "max step", "rms step", "max |g|", "rms |g|"];

/// Reads the per-cycle convergence table:
///
/// ```text
/// step               E[Eh]         dE[kJ/mol]    max_step    rms_step    max_grad    rms_grad
///     1        -5.76633308          -9.521052    2.70e-02    1.68e-02    4.07e-02    1.98e-02
///     2        -5.76858971          -5.924783    5.86e-02    3.30e-02    9.65e-03    5.59e-03
/// ```
///
/// Rows are recognised by shape -- seven numbers, the first a bare integer --
/// rather than by counting lines after the header, so the surrounding banner
/// text can change without silently reading the wrong lines. A partial table
/// from a cancelled run yields the rows it did print.
pub fn parse_optimization_steps(stdout: &str) -> Vec<OptimizationStep> {
    let mut steps = Vec::new();
    for line in stdout.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != 7 {
            continue;
        }
        // The step number is a bare integer; the header's "step" is not, which
        // is what keeps the header itself out of the table.
        let Ok(step) = fields[0].parse::<usize>() else {
            continue;
        };
        let mut values = [0.0f64; 6];
        let mut parsed = true;
        for (slot, text) in values.iter_mut().zip(&fields[1..]) {
            match text.parse::<f64>() {
                Ok(value) => *slot = value,
                Err(_) => {
                    parsed = false;
                    break;
                }
            }
        }
        if !parsed {
            continue;
        }
        steps.push(OptimizationStep {
            step,
            energy_hartree: values[0],
            de_kj_mol: values[1],
            max_step: values[2],
            rms_step: values[3],
            max_grad: values[4],
            rms_grad: values[5],
        });
    }
    steps
}

/// Reads the thresholds line that sits above the table:
///
/// ```text
/// Threshold Values:                 1.31e-02     4.00e-03    2.00e-03    3.00e-04    1.00e-04
/// ```
///
/// Its five numbers are in the table's own column order -- dE, max step, rms
/// step, max gradient, rms gradient -- so they line up with a row directly.
pub fn parse_convergence_thresholds(stdout: &str) -> Option<ConvergenceThresholds> {
    let line = stdout
        .lines()
        .find(|line| line.trim_start().starts_with("Threshold Values"))?;
    let after_colon = line.split(':').nth(1)?;
    let values: Vec<f64> = after_colon
        .split_whitespace()
        .filter_map(|token| token.parse().ok())
        .collect();
    if values.len() != 5 {
        return None;
    }
    Some(ConvergenceThresholds {
        de_kj_mol: values[0],
        max_step: values[1],
        rms_step: values[2],
        max_grad: values[3],
        rms_grad: values[4],
    })
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
    /// Behemoth's own output, verbatim from a TASI gradient run on water.
    ///
    /// The energy is not the last thing on its line and the atom rows carry an index, so both are
    /// easy to get wrong; this is the real text rather than a reconstruction.
    #[test]
    fn a_gradient_block_is_read_from_real_output() {
        let stdout = concat!(
            "TASI SEOEM energy\n",
            "Total energy (all terms):        -68.4051345926 Hartree\n",
            "\n",
            "Analytical TASI nuclear gradient\n",
            "Units: Hartree/bohr\n",
            "Energy:     -68.4051345926 Hartree\n",
            "     1     -0.0020104     -0.0021028     -0.0000000\n",
            "     2     -0.0026798      0.0040516      0.0000000\n",
            "     3      0.0046902     -0.0019488     -0.0000000\n",
            "\n",
            "Total run time (s): 0.009\n"
        );
        let (energy, gradient) = parse_gradient(stdout, 3).expect("the block parses");
        assert!((energy + 68.405_134_592_6).abs() < 1.0e-10);
        assert_eq!(gradient.len(), 9);
        assert!((gradient[0] + 0.002_010_4).abs() < 1.0e-9);
        assert!((gradient[4] - 0.004_051_6).abs() < 1.0e-9);
        assert!((gradient[6] - 0.004_690_2).abs() < 1.0e-9);
    }

    /// A run that printed no gradient must be an error, not a zero gradient -- which the optimiser
    /// would read as a converged structure and accept.
    #[test]
    fn output_without_a_gradient_is_an_error() {
        assert!(parse_gradient("Total energy: -1.0 Hartree\n", 3).is_err());
        // A truncated block is an error too, rather than a short gradient.
        let truncated = concat!(
            "Analytical TASI nuclear gradient\n",
            "Units: Hartree/bohr\n",
            "Energy:     -68.4051345926 Hartree\n",
            "     1     -0.0020104     -0.0021028     -0.0000000\n"
        );
        assert!(parse_gradient(truncated, 3).is_err());
    }

    use super::*;
    use super::super::method::{BehemothMethod, Dispersion};

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
            &MethodConfig::default(),
        );
        let args: Vec<String> = command
            .get_args()
            .map(|a| a.to_string_lossy().to_string())
            .collect();
        assert_eq!(
            args,
            vec![
                "--xyz",
                "input.xyz",
                "--method",
                "xtb",
                "--charge",
                "-1",
                "--multi",
                "2",
                "--memory",
                "1024",
                "--nproc",
                "1",
                "--opt"
            ]
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
            &MethodConfig::default(),
        );
        let args: Vec<String> = command
            .get_args()
            .map(|a| a.to_string_lossy().to_string())
            .collect();
        assert!(args.contains(&"--hessian".to_string()));
        assert!(!args.contains(&"--freq".to_string()), "--freq has no modes");
    }

    /// The flags for a plain DFT job: the functional, the basis and the
    /// dispersion correction all travel, and a singlet is `rks`.
    #[test]
    fn a_dft_command_carries_the_functional_basis_and_dispersion() {
        let method = MethodConfig {
            behemoth: BehemothMethod::Dft,
            functional: "pbe0".to_string(),
            basis: "def2-tzvp".to_string(),
            dispersion: Dispersion::D4,
            memory_mb: 4096,
            nproc: 8,
            ..Default::default()
        };
        let args = args_of(optimize_command(
            Path::new("b"),
            Path::new("."),
            "in.xyz",
            0,
            1,
            &method,
        ));
        assert_eq!(value_after(&args, "--method"), Some("rks"));
        assert_eq!(value_after(&args, "--functional"), Some("pbe0"));
        assert_eq!(value_after(&args, "--basis"), Some("def2-tzvp"));
        assert_eq!(value_after(&args, "--dispersion"), Some("d4"));
        assert_eq!(value_after(&args, "--memory"), Some("4096"));
        assert_eq!(value_after(&args, "--nproc"), Some("8"));
    }

    /// An open shell is always the unrestricted method, as asked: `uks` for
    /// DFT -- which Behemoth also enforces, refusing `rks` above a singlet --
    /// and `uhf` for Hartree-Fock.
    #[test]
    fn an_open_shell_command_is_unrestricted() {
        let dft = MethodConfig { behemoth: BehemothMethod::Dft, ..Default::default() };
        let args = args_of(hessian_command(Path::new("b"), Path::new("."), "in.xyz", 0, 3, &dft));
        assert_eq!(value_after(&args, "--method"), Some("uks"));

        let hf = MethodConfig { behemoth: BehemothMethod::Hf, ..Default::default() };
        let args = args_of(hessian_command(Path::new("b"), Path::new("."), "in.xyz", 0, 2, &hf));
        assert_eq!(value_after(&args, "--method"), Some("uhf"));

        let singlet =
            args_of(hessian_command(Path::new("b"), Path::new("."), "in.xyz", 0, 1, &hf));
        assert_eq!(value_after(&singlet, "--method"), Some("hf"));
    }

    /// A composite defines its own basis set and dispersion correction, so
    /// neither may be passed alongside it -- Behemoth would be given two
    /// answers to the same question.
    #[test]
    fn a_composite_command_passes_no_basis_or_dispersion() {
        let method = MethodConfig {
            behemoth: BehemothMethod::Dft,
            functional: "r2scan-3c".to_string(),
            basis: "cc-pvtz".to_string(),
            dispersion: Dispersion::D3Bj,
            ..Default::default()
        };
        let args = args_of(optimize_command(
            Path::new("b"),
            Path::new("."),
            "in.xyz",
            0,
            1,
            &method,
        ));
        assert_eq!(value_after(&args, "--functional"), Some("r2scan-3c"));
        assert!(!args.iter().any(|a| a == "--basis"), "{args:?}");
        assert!(!args.iter().any(|a| a == "--dispersion"), "{args:?}");
    }

    /// HF-3c is the one composite that is a method rather than a functional.
    #[test]
    fn hf3c_travels_as_the_method_with_no_functional() {
        let method = MethodConfig {
            behemoth: BehemothMethod::Dft,
            functional: "hf-3c".to_string(),
            ..Default::default()
        };
        let args = args_of(optimize_command(
            Path::new("b"),
            Path::new("."),
            "in.xyz",
            0,
            1,
            &method,
        ));
        assert_eq!(value_after(&args, "--method"), Some("hf-3c"));
        assert!(!args.iter().any(|a| a == "--functional"), "{args:?}");
        assert!(!args.iter().any(|a| a == "--basis"), "{args:?}");
    }

    /// A wavefunction method takes a basis but no functional or dispersion.
    #[test]
    fn a_correlated_command_carries_only_the_basis() {
        for (method, expected) in
            [(BehemothMethod::Mp2, "mp2"), (BehemothMethod::RiMp2, "ri-mp2")]
        {
            let config = MethodConfig {
                behemoth: method,
                dispersion: Dispersion::D3Bj,
                ..Default::default()
            };
            let args = args_of(optimize_command(
                Path::new("b"),
                Path::new("."),
                "in.xyz",
                0,
                1,
                &config,
            ));
            assert_eq!(value_after(&args, "--method"), Some(expected));
            assert_eq!(value_after(&args, "--basis"), Some("def2-svp"));
            assert!(!args.iter().any(|a| a == "--functional"), "{method:?}");
            assert!(!args.iter().any(|a| a == "--dispersion"), "{method:?}");
        }
    }

    /// No dispersion selected means no flag, not an empty one.
    #[test]
    fn no_dispersion_selected_passes_no_flag() {
        let method = MethodConfig { behemoth: BehemothMethod::Dft, ..Default::default() };
        let args = args_of(optimize_command(
            Path::new("b"),
            Path::new("."),
            "in.xyz",
            0,
            1,
            &method,
        ));
        assert!(!args.iter().any(|a| a == "--dispersion"), "{args:?}");
    }

    /// The convergence table, from the captured log of a real run.
    #[test]
    fn the_convergence_table_is_read_row_by_row() {
        let log = include_str!("../../tests/fixtures/behemoth_water_opt.log");
        let steps = parse_optimization_steps(log);
        assert_eq!(steps.len(), 3, "three cycles in the fixture");

        assert_eq!(steps[0].step, 1);
        assert!((steps[0].energy_hartree - -5.76866785).abs() < 1e-8);
        assert!((steps[0].de_kj_mol - -0.573303).abs() < 1e-6);
        assert!((steps[0].max_step - 1.49e-02).abs() < 1e-10);
        assert!((steps[0].rms_step - 6.88e-03).abs() < 1e-10);
        assert!((steps[0].max_grad - 8.81e-03).abs() < 1e-10);
        assert!((steps[0].rms_grad - 3.83e-03).abs() < 1e-10);

        // The last row is the converged one, and the summary's own numbers
        // have to agree with it.
        let last = steps.last().unwrap();
        assert_eq!(last.step, 3);
        assert!((last.rms_grad - 5.69e-06).abs() < 1e-10);
        let summary = parse_summary(log);
        assert_eq!(summary.cycles, Some(steps.len()));
        assert!((summary.rms_gradient.unwrap() - 0.0000056858).abs() < 1e-10);
    }

    /// The header line contains the word "step", not a number, which is what
    /// keeps it from being read as a row.
    #[test]
    fn the_table_header_is_not_read_as_a_step() {
        let steps = parse_optimization_steps(
            "step               E[Eh]         dE[kJ/mol]    max_step    rms_step    max_grad    rms_grad\n                 1        -5.76866785          -0.573303    1.49e-02    6.88e-03    8.81e-03    3.83e-03\n",
        );
        assert_eq!(steps.len(), 1);
        assert_eq!(steps[0].step, 1);
    }

    /// A run that was cancelled or died mid-table still yields the cycles it
    /// managed, so the summary can show how far it got.
    #[test]
    fn a_partial_table_yields_the_rows_it_printed() {
        let steps = parse_optimization_steps(
            "step               E[Eh]         dE[kJ/mol]    max_step    rms_step    max_grad    rms_grad\n                 1        -5.76866785          -0.573303    1.49e-02    6.88e-03    8.81e-03    3.83e-03\n                 2        -5.76877493   \n",
        );
        assert_eq!(steps.len(), 1, "the truncated row is not half-read");
    }

    #[test]
    fn a_log_with_no_table_yields_no_steps() {
        assert!(parse_optimization_steps("nothing here\n").is_empty());
        assert!(parse_optimization_steps("").is_empty());
    }

    #[test]
    fn the_thresholds_are_read_in_the_tables_column_order() {
        let log = include_str!("../../tests/fixtures/behemoth_water_opt.log");
        let t = parse_convergence_thresholds(log).expect("a thresholds line");
        assert!((t.de_kj_mol - 1.31e-02).abs() < 1e-10);
        assert!((t.max_step - 4.00e-03).abs() < 1e-10);
        assert!((t.rms_step - 2.00e-03).abs() < 1e-10);
        assert!((t.max_grad - 3.00e-04).abs() < 1e-10);
        assert!((t.rms_grad - 1.00e-04).abs() < 1e-10);
    }

    #[test]
    fn a_log_with_no_thresholds_line_reports_none() {
        assert!(parse_convergence_thresholds("no thresholds here").is_none());
        // Five numbers are expected; a truncated line is not guessed at.
        assert!(parse_convergence_thresholds("Threshold Values:  1.0  2.0").is_none());
    }

    /// The converged final step must satisfy every criterion, and an early one
    /// must not -- otherwise the per-column marks in the summary window would
    /// be meaningless.
    #[test]
    fn the_final_step_meets_every_threshold_and_the_first_does_not() {
        let log = include_str!("../../tests/fixtures/behemoth_water_opt.log");
        let steps = parse_optimization_steps(log);
        let thresholds = parse_convergence_thresholds(log).unwrap();

        assert_eq!(
            thresholds.met_by(steps.last().unwrap()),
            [true; 5],
            "the run converged, so its last cycle is inside every threshold"
        );
        assert!(
            thresholds.met_by(&steps[0]).iter().any(|met| !met),
            "the first cycle cannot already be converged"
        );
    }

    /// dE is signed -- a descent prints it negative -- so it is the magnitude
    /// that is compared, not the value.
    #[test]
    fn a_large_energy_drop_is_not_counted_as_converged() {
        let thresholds = ConvergenceThresholds {
            de_kj_mol: 1.31e-02,
            max_step: 4.00e-03,
            rms_step: 2.00e-03,
            max_grad: 3.00e-04,
            rms_grad: 1.00e-04,
        };
        let step = OptimizationStep {
            step: 1,
            energy_hartree: -5.0,
            de_kj_mol: -9.52,
            max_step: 1.0e-03,
            rms_step: 1.0e-03,
            max_grad: 1.0e-04,
            rms_grad: 1.0e-05,
        };
        assert_eq!(thresholds.met_by(&step)[0], false, "a 9.5 kJ/mol drop is not converged");
    }

    fn args_of(command: Command) -> Vec<String> {
        command
            .get_args()
            .map(|a| a.to_string_lossy().to_string())
            .collect()
    }

    /// The value a flag was given, so a test names the flag rather than
    /// hard-coding a position that shifts when a flag is added.
    fn value_after<'a>(args: &'a [String], flag: &str) -> Option<&'a str> {
        let at = args.iter().position(|a| a == flag)?;
        args.get(at + 1).map(|s| s.as_str())
    }

    /// A multiplicity of zero is not a thing; clamping keeps a bad panel value
    /// from reaching the command line.
    #[test]
    fn a_multiplicity_below_one_is_clamped() {
        let command = optimize_command(
            Path::new("b"),
            Path::new("."),
            "in.xyz",
            0,
            0,
            &MethodConfig::default(),
        );
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
        let status = optimize_command(&binary, &dir, "input.xyz", 0, 1, &MethodConfig::default())
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
        let hess = hessian_command(&binary, &dir, "input.xyz", 0, 1, &MethodConfig::default())
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
