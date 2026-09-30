//! Driving Sparrow through SCINE's Python bindings.
//!
//! Sparrow computes three things: an energy, an analytic gradient and an
//! analytic Hessian. It has no geometry optimiser and no frequency analysis, so
//! those are Beavyr's to do -- its optimiser drives the gradients, and its own
//! projection and eigensolver turn the Hessian into modes and thermochemistry.
//!
//! That division is why the generated script only ever asks for a property. It
//! never loops and never analyses, which also makes it the simplest of the
//! Python scripts here.
//!
//! The SCINE calls follow the user's own `sparrowpy.py`: the module manager
//! hands out a calculator by parametrisation name, the structure is read from a
//! file, and the properties wanted are declared before `calculate()`. Checked
//! against the installed build, where all eight parametrisations resolve and
//! both `Property.Gradients` and `Property.Hessian` exist.

use std::path::{Path, PathBuf};
use std::process::Command;

use super::job::JobType;
use super::method::MethodConfig;
use super::python_env::PYTHON;

/// The script Beavyr writes into the run directory.
pub const SCRIPT_FILE: &str = "run_sparrow.py";
/// The geometry the script reads, since SCINE reads a structure from a file.
pub const GEOMETRY_FILE: &str = "geometry.xyz";
/// The gradient, as a Turbomole `$grad` block Beavyr already reads.
pub const GRADIENT_FILE: &str = "sparrow.grad";
/// The Hessian, behind a `$hessian` header.
pub const HESSIAN_FILE: &str = "hessian.txt";
/// The molecular orbitals, when asked for.
pub const MOLDEN_FILE: &str = "orbitals.molden.input";

/// Where each file of a run in `run_dir` lives.
#[derive(Debug, Clone)]
pub struct JobFiles {
    pub script: PathBuf,
    pub geometry: PathBuf,
    pub gradient: PathBuf,
    pub hessian: PathBuf,
    pub molden: PathBuf,
}

pub fn job_files(run_dir: &Path) -> JobFiles {
    JobFiles {
        script: run_dir.join(SCRIPT_FILE),
        geometry: run_dir.join(GEOMETRY_FILE),
        gradient: run_dir.join(GRADIENT_FILE),
        hessian: run_dir.join(HESSIAN_FILE),
        molden: run_dir.join(MOLDEN_FILE),
    }
}

/// The command for a prepared script.
pub fn run_command(run_dir: &Path) -> Command {
    let mut command = Command::new(PYTHON);
    command.current_dir(run_dir).arg(SCRIPT_FILE);
    command
}

/// Which SCINE properties this job needs.
///
/// A Hessian is far the most expensive of the three, so it is asked for only
/// when it is wanted. A gradient is nearly free beside it and comes along with
/// anything that needs one.
fn properties_for(job: JobType) -> &'static str {
    match job {
        JobType::SinglePoint => "[su.Property.Energy]",
        JobType::Gradient => "[su.Property.Energy, su.Property.Gradients]",
        // Beavyr does the frequency analysis; what it needs from Sparrow is
        // the matrix.
        JobType::Hessian | JobType::Frequencies | JobType::OptimizeThenFrequencies => {
            "[su.Property.Energy, su.Property.Gradients, su.Property.Hessian]"
        }
        // An optimisation step: Beavyr's optimiser asks for one of these per
        // geometry and decides where to go next.
        JobType::Optimize => "[su.Property.Energy, su.Property.Gradients]",
        // Neither exists for a semi-empirical calculator; the panel refuses
        // these before a script is ever written.
        JobType::TransitionState | JobType::Irc | JobType::TdDft => {
            "[su.Property.Energy]"
        }
    }
}

/// The geometry file SCINE reads.
pub fn geometry_text(atoms: &[String], positions_angstrom: &[f64]) -> String {
    let mut text = format!("{}\nWritten by Beavyr for Sparrow\n", atoms.len());
    for (index, symbol) in atoms.iter().enumerate() {
        let base = index * 3;
        text.push_str(&format!(
            "{symbol:<3} {:18.10} {:18.10} {:18.10}\n",
            positions_angstrom[base],
            positions_angstrom[base + 1],
            positions_angstrom[base + 2]
        ));
    }
    text
}

/// Builds the script for one job.
pub fn script_text(
    job: JobType,
    charge: i32,
    multiplicity: i32,
    method: &MethodConfig,
    wavefunction: bool,
) -> String {
    let multiplicity = multiplicity.max(1);
    // Sparrow's parametrisation is stored in the functional field, which is the
    // one string the shared config already carries for "which method".
    let parametrisation = super::sparrow_method::method(&method.functional)
        .map_or(method.functional.as_str(), |m| m.name);
    let properties = properties_for(job);
    let wants_gradient = matches!(
        job,
        JobType::Gradient
            | JobType::Optimize
            | JobType::Hessian
            | JobType::Frequencies
            | JobType::OptimizeThenFrequencies
    );
    let wants_hessian = matches!(
        job,
        JobType::Hessian | JobType::Frequencies | JobType::OptimizeThenFrequencies
    );

    let mut text = format!(
        r#"# Written by Beavyr. One Sparrow calculation, results beside this file.
import scine_utilities as su
import scine_sparrow  # noqa: F401  (registers the Sparrow module with SCINE)

manager = su.core.ModuleManager.get_instance()
calculator = manager.get("calculator", {parametrisation:?})
calculator.structure = su.io.read({geometry:?})[0]
calculator.set_required_properties({properties})

calculator.settings["molecular_charge"] = {charge}
calculator.settings["spin_multiplicity"] = {multiplicity}
if {multiplicity} != 1:
    # An open shell has to be unrestricted or the calculation describes a
    # different electronic state from the one asked for.
    calculator.settings["spin_mode"] = "unrestricted"

results = calculator.calculate()
energy = float(results.energy)
symbols = [str(a.element) for a in calculator.structure]
coords = calculator.structure.positions  # bohr, which is what both files want
"#,
        geometry = GEOMETRY_FILE,
    );

    if wants_gradient {
        text.push_str(
            r#"
# Turbomole `$grad`, the layout Beavyr already reads.
g = results.gradients
with open("sparrow.grad", "w") as handle:
    handle.write("$grad          cartesian gradients\n")
    handle.write("  cycle =      1    SCF energy =  %20.14f   |dE/dxyz| =  0.000000\n" % energy)
    for symbol, xyz in zip(symbols, coords):
        handle.write("%22.14f %22.14f %22.14f      %s\n" % (xyz[0], xyz[1], xyz[2], symbol))
    for row in g:
        handle.write("%22.14e %22.14e %22.14e\n" % (row[0], row[1], row[2]))
    handle.write("$end\n")
"#,
        );
    }

    if wants_hessian {
        text.push_str(
            r#"
# Behind a `$hessian` header, so Beavyr's Turbomole reader takes it and its own
# projection and eigensolver do the frequency analysis. Sparrow has none.
H = results.hessian
with open("hessian.txt", "w") as handle:
    handle.write("$hessian\n")
    for row in H:
        for value in row:
            handle.write(" %22.14e" % value)
        handle.write("\n")
    handle.write("$end\n")
"#,
        );
    }

    if wavefunction {
        // SCINE's own wavefunction writer, the route the user's `sparrowpy.py`
        // already uses. A semi-empirical calculator has a minimal basis, so the
        // file is small and the step is cheap.
        text.push_str(
            r#"
wf_gen = su.core.to_wf_generator(calculator)
wf_gen.wavefunction2file("orbitals.molden.input")
"#,
        );
    }

    text.push_str("\nprint(\"energy %.12f\" % energy)\n");
    text
}

/// The energy the script printed, in Hartree.
pub fn parse_energy(text: &str) -> Option<f64> {
    text.lines()
        .filter_map(|line| line.strip_prefix("energy "))
        .next_back()
        .and_then(|rest| rest.trim().parse().ok())
}

/// The last thing said before giving up.
pub fn failure_detail(stderr: &str) -> String {
    stderr
        .lines()
        .rev()
        .find(|line| !line.trim().is_empty())
        .unwrap_or("see the run directory")
        .trim()
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config(parametrisation: &str) -> MethodConfig {
        MethodConfig {
            functional: parametrisation.to_string(),
            ..Default::default()
        }
    }

    #[test]
    fn the_parametrisation_reaches_scine() {
        let text = script_text(JobType::SinglePoint, 0, 1, &config("PM6"), false);
        assert!(text.contains(r#"manager.get("calculator", "PM6")"#), "{text}");
    }

    /// The label carries "non-SCC" for readability; SCINE is given the bare
    /// name, and sending the label would fail.
    #[test]
    fn the_bare_name_is_sent_not_the_label() {
        let text = script_text(JobType::SinglePoint, 0, 1, &config("DFTB0"), false);
        assert!(text.contains(r#""DFTB0""#), "{text}");
        assert!(!text.contains("non-SCC"), "{text}");
    }

    /// An open shell must be unrestricted.
    #[test]
    fn an_open_shell_is_unrestricted() {
        let text = script_text(JobType::SinglePoint, -1, 2, &config("PM6"), false);
        assert!(text.contains(r#"["molecular_charge"] = -1"#), "{text}");
        assert!(text.contains(r#"["spin_multiplicity"] = 2"#), "{text}");
        assert!(text.contains("unrestricted"), "{text}");
    }

    /// A Hessian is far the most expensive of the three, so it is asked for
    /// only where it is wanted. A single point must not pay for one.
    #[test]
    fn each_job_asks_only_for_what_it_needs() {
        let sp = script_text(JobType::SinglePoint, 0, 1, &config("PM6"), false);
        assert!(!sp.contains("Property.Hessian"), "{sp}");
        assert!(!sp.contains("Property.Gradients"), "{sp}");
        assert!(!sp.contains("hessian.txt"), "{sp}");

        let grad = script_text(JobType::Gradient, 0, 1, &config("PM6"), false);
        assert!(grad.contains("Property.Gradients"), "{grad}");
        assert!(!grad.contains("Property.Hessian"), "{grad}");
        assert!(grad.contains("$grad"), "{grad}");

        let hess = script_text(JobType::Hessian, 0, 1, &config("PM6"), false);
        assert!(hess.contains("Property.Hessian"), "{hess}");
        assert!(hess.contains("$hessian"), "{hess}");
    }

    /// Beavyr does the frequency analysis, so what a frequency job needs from
    /// Sparrow is the matrix and nothing more. The script never analyses.
    #[test]
    fn a_frequency_job_only_fetches_the_matrix() {
        let text = script_text(JobType::Frequencies, 0, 1, &config("PM6"), false);
        assert!(text.contains("Property.Hessian"), "{text}");
        assert!(text.contains("$hessian"), "{text}");
        // Nothing in the script diagonalises anything or writes frequencies.
        // Asserted on calls rather than on words, because the comments in the
        // script legitimately mention the frequency analysis Beavyr does.
        for call in ["eigh(", "eigvals", "numpy.linalg", "frequencies(", "normal_modes"] {
            assert!(!text.contains(call), "the script must not analyse: {call}");
        }
        // The matrix is written out and that is the end of it.
        assert!(text.contains(r#"handle.write("$hessian\n")"#), "{text}");
    }

    /// The script never optimises: Sparrow has no optimiser, and Beavyr's own
    /// drives it one gradient at a time instead.
    #[test]
    fn the_script_never_optimises() {
        let text = script_text(JobType::Optimize, 0, 1, &config("PM6"), false);
        assert!(text.contains("Property.Gradients"), "{text}");
        assert!(!text.contains("optimi"), "{text}");
        // One calculation, not a loop.
        assert_eq!(text.matches("calculator.calculate()").count(), 1, "{text}");
    }

    /// The orbitals are written only when asked for: it is a second step, and
    /// nobody wants the file unless they said so.
    #[test]
    fn the_orbitals_are_written_only_on_request() {
        let without = script_text(JobType::SinglePoint, 0, 1, &config("PM6"), false);
        assert!(!without.contains("wavefunction2file"), "{without}");

        let with = script_text(JobType::SinglePoint, 0, 1, &config("PM6"), true);
        assert!(with.contains("to_wf_generator"), "{with}");
        assert!(with.contains("orbitals.molden.input"), "{with}");
    }

    #[test]
    fn the_geometry_file_is_an_ordinary_xyz() {
        let atoms: Vec<String> = ["O", "H"].iter().map(|s| s.to_string()).collect();
        let text = geometry_text(&atoms, &[0.0, 0.0, 0.0, 0.0, 0.0, 0.96]);
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0].trim(), "2");
        assert!(lines[2].starts_with('O'), "{text}");
        assert!(lines[3].contains("0.9600000000"), "{text}");
    }

    #[test]
    fn the_energy_is_read_by_key() {
        assert!(parse_energy("chatter\nenergy -11.123456789012\n").is_some());
        assert_eq!(parse_energy("nothing\n"), None);
    }
}
