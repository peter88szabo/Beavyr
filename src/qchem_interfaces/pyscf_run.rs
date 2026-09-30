//! Driving PySCF by generating a Python script and running it.
//!
//! PySCF has no command line, so unlike every other backend there is nothing
//! to point at and no flags to build. What Beavyr writes instead is a small
//! script, run with `python3` from the environment Beavyr itself was launched
//! in. Whether that environment can import PySCF is settled before the run, by
//! [`super::python_env`].
//!
//! The script writes its results into files rather than printing them, and it
//! writes them in formats Beavyr already reads:
//!
//! * the optimisation path as a multi-frame XYZ with `cycle <n> E = <v> Eh`
//!   comments, which is what the trajectory player and energy plot take;
//! * the Hessian behind a `$hessian` header, which is Turbomole's layout and
//!   is read by [`super::hessian_file::parse_turbomole_hessian`];
//! * the final geometry as an ordinary XYZ.
//!
//! So this module generates a script and names files. It parses no numbers,
//! because nothing here needed a new format.
//!
//! Checked against PySCF 2.12.1: the `geometric` optimiser is present, its
//! callback carries `coords`, `energy` and `gradients` per step, and the
//! analytic Hessian comes back symmetric to eleven decimal places on water.

use std::path::{Path, PathBuf};
use std::process::Command;

use super::job::JobType;
use super::method::MethodConfig;
use super::pyscf_method;
use super::python_env::PYTHON;

/// The script Beavyr writes into the run directory.
pub const SCRIPT_FILE: &str = "run_pyscf.py";
/// The optimisation path, already in Beavyr's own trajectory format.
pub const TRAJECTORY_FILE: &str = "trajectory.xyz";
/// The final geometry, as an ordinary XYZ.
pub const GEOMETRY_FILE: &str = "final.xyz";
/// The Hessian, behind a `$hessian` header so the Turbomole reader takes it.
pub const HESSIAN_FILE: &str = "hessian.txt";
/// The molecular orbitals, when asked for.
pub const MOLDEN_FILE: &str = "orbitals.molden.input";

// What the script prints is captured by the caller into the run directory's
// shared `xtb.stdout`/`xtb.stderr`, the same two files every backend's output
// goes to, so there are no names for them here.

/// The gradient, as a Turbomole `$grad` block Beavyr already reads.
pub const GRADIENT_FILE: &str = "pyscf.grad";
/// The excited states, written by the script in a format of our own.
pub const SPECTRUM_FILE: &str = "spectrum.txt";

/// Where each file of a run in `run_dir` lives.
#[derive(Debug, Clone)]
pub struct JobFiles {
    pub script: PathBuf,
    pub trajectory: PathBuf,
    pub geometry: PathBuf,
    pub hessian: PathBuf,
    pub gradient: PathBuf,
    pub spectrum: PathBuf,
    pub molden: PathBuf,
}

pub fn job_files(run_dir: &Path) -> JobFiles {
    JobFiles {
        script: run_dir.join(SCRIPT_FILE),
        trajectory: run_dir.join(TRAJECTORY_FILE),
        geometry: run_dir.join(GEOMETRY_FILE),
        hessian: run_dir.join(HESSIAN_FILE),
        gradient: run_dir.join(GRADIENT_FILE),
        spectrum: run_dir.join(SPECTRUM_FILE),
        molden: run_dir.join(MOLDEN_FILE),
    }
}

/// The command that runs a prepared script.
pub fn run_command(run_dir: &Path) -> Command {
    let mut command = Command::new(PYTHON);
    command.current_dir(run_dir).arg(SCRIPT_FILE);
    command
}

/// The geometry block the script embeds, one `symbol x y z` line per atom in
/// angstrom, which is the unit `gto.M` is told to read.
fn atom_block(atoms: &[String], positions_angstrom: &[f64]) -> String {
    let mut block = String::new();
    for (index, symbol) in atoms.iter().enumerate() {
        let base = index * 3;
        block.push_str(&format!(
            "{symbol} {:.10} {:.10} {:.10}\\n",
            positions_angstrom[base],
            positions_angstrom[base + 1],
            positions_angstrom[base + 2]
        ));
    }
    block
}

/// Builds the script for one job.
///
/// Every value from Beavyr is written as a Python literal rather than
/// interpolated into an expression, so a method name with awkward characters
/// cannot become code.
#[allow(clippy::too_many_arguments)]
pub fn script_text(
    job: JobType,
    charge: i32,
    multiplicity: i32,
    method: &MethodConfig,
    excited: super::job::ExcitedStateOptions,
    wavefunction: bool,
    atoms: &[String],
    positions_angstrom: &[f64],
) -> String {
    let multiplicity = multiplicity.max(1);
    // PySCF counts unpaired electrons, not multiplicity.
    let spin = multiplicity - 1;
    let nproc = method.nproc.max(1);
    // Only the branch that applies is written. Emitting both behind a constant
    // condition would leave `mf.xc = "hf"` sitting in dead code in a
    // Hartree-Fock script, which is exactly the sort of thing someone reading
    // the run directory to work out what ran would be misled by.
    let build_body = if pyscf_method::is_hartree_fock(&method.functional) {
        "    from pyscf import scf\n\
         \x20   return scf.UHF(molecule) if molecule.spin else scf.RHF(molecule)\n"
            .to_string()
    } else {
        format!(
            "    from pyscf import dft\n\
             \x20   mf = dft.UKS(molecule) if molecule.spin else dft.RKS(molecule)\n\
             \x20   mf.xc = {functional:?}\n\
             \x20   return mf\n",
            functional = method.functional
        )
    };
    let geometry = atom_block(atoms, positions_angstrom);
    let basis = &method.basis;

    let mut script = String::new();
    script.push_str(&format!(
        r#"# Written by Beavyr. Runs one PySCF job and leaves its results beside this file.
import sys
import numpy
from pyscf import gto, lib

lib.num_threads({nproc})

ATOMS = {atoms:?}
GEOMETRY = "{geometry}"

mol = gto.M(atom=GEOMETRY, basis={basis:?}, charge={charge}, spin={spin},
            unit="Angstrom", verbose=0, output=None)

BOHR = 0.52917721092

def build(molecule):
    """The mean-field object for this molecule, restricted or not by spin."""
{build_body}

def write_xyz(path, molecule, comment):
    coords = molecule.atom_coords() * BOHR
    with open(path, "w") as handle:
        handle.write("%d\n%s\n" % (len(ATOMS), comment))
        for symbol, xyz in zip(ATOMS, coords):
            handle.write("%-3s %18.10f %18.10f %18.10f\n" % (symbol, xyz[0], xyz[1], xyz[2]))

"#
    ));

    if job.produces_trajectory() {
        script.push_str(
            r#"# The optimisation path, written straight into Beavyr's own trajectory
# format so the trajectory player and the energy plot read it unchanged.
frames = []

def record(envs):
    """One frame per optimiser step: geometry in angstrom, energy in Hartree."""
    try:
        coords = numpy.asarray(envs["coords"], dtype=float).reshape(len(ATOMS), 3) * BOHR
        energy = float(envs["energy"])
    except Exception:
        return
    frames.append((coords, energy))

from pyscf.geomopt.geometric_solver import optimize
mf = build(mol)
optimised = optimize(mf, callback=record)

with open("trajectory.xyz", "w") as handle:
    for index, (coords, energy) in enumerate(frames, start=1):
        handle.write("%d\ncycle %d E = %.10f Eh\n" % (len(ATOMS), index, energy))
        for symbol, xyz in zip(ATOMS, coords):
            handle.write("%-3s %18.10f %18.10f %18.10f\n" % (symbol, xyz[0], xyz[1], xyz[2]))

mol = optimised
write_xyz("final.xyz", mol, "PySCF optimised geometry")
"#,
        );
    }

    // The energy is always wanted, and after an optimisation it has to be the
    // optimised geometry's, so the mean field is rebuilt on whichever molecule
    // we ended up with rather than reusing the one the optimiser started from.
    script.push_str(
        r#"mf = build(mol)
energy = float(mf.kernel())
if not mf.converged:
    sys.stderr.write("PySCF SCF did not converge.\n")
    sys.exit(2)
if not numpy.isfinite(energy):
    sys.stderr.write("PySCF returned a non-finite energy: %r\n" % (energy,))
    sys.exit(2)
"#,
    );

    if !job.produces_trajectory() {
        script.push_str("write_xyz(\"final.xyz\", mol, \"PySCF geometry\")\n");
    }

    // The gradient, written where Beavyr already reads one: Turbomole's
    // `$grad`, the same layout the Psi4 route writes and the same one
    // `parse_turbomole_gradient` takes.
    if job == JobType::Gradient {
        script.push_str(
            r#"
g = mf.nuc_grad_method().kernel()
coords = mol.atom_coords()
with open("pyscf.grad", "w") as handle:
    handle.write("$grad          cartesian gradients\n")
    handle.write("  cycle =      1    SCF energy =  %20.14f   |dE/dxyz| =  0.000000\n" % energy)
    for symbol, xyz in zip(ATOMS, coords):
        handle.write("%22.14f %22.14f %22.14f      %s\n" % (xyz[0], xyz[1], xyz[2], symbol))
    for row in g:
        handle.write("%22.14e %22.14e %22.14e\n" % (row[0], row[1], row[2]))
    handle.write("$end\n")
"#,
        );
    }

    if job.produces_frequencies() || job == JobType::Hessian {
        script.push_str(
            r#"
# Behind a `$hessian` header, which is Turbomole's layout: Beavyr already has a
# reader for it, so this needs no format of its own.
h = mf.Hessian().kernel()
n = mol.natm
H = numpy.asarray(h).transpose(0, 2, 1, 3).reshape(n * 3, n * 3)
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

    if job == JobType::TdDft {
        script.push_str(&format!(
            r##"
from pyscf import tdscf

td = tdscf.TDA(mf) if {tda} else tdscf.TDDFT(mf)
td.nstates = {states}
td.kernel()
strengths = td.oscillator_strength()

# A small format of our own: the script is ours, so there is no reason to
# print a table and read it back with a pattern. Energies in eV.
HARTREE_TO_EV = 27.211386245988
with open("spectrum.txt", "w") as handle:
    handle.write("# root  energy_eV  oscillator_strength\n")
    for index, (e, f) in enumerate(zip(td.e, strengths), start=1):
        handle.write("%d %.10f %.10f\n" % (index, float(e) * HARTREE_TO_EV, float(f)))

# And in PySCF's own wording as well, so the UV-Vis reader written against
# `analyze()` takes this run exactly as it takes a log the user produced.
print("\n** Singlet excitation energies and oscillator strengths **")
for index, (e, f) in enumerate(zip(td.e, strengths), start=1):
    ev = float(e) * HARTREE_TO_EV
    nm = 1239.8419843320026 / ev if ev > 0 else 0.0
    print("Excited State %3d: %12.5f eV %9.2f nm  f=%.4f" % (index, ev, nm, float(f)))
"##,
            states = excited.states.max(1),
            tda = if excited.tda { "True" } else { "False" }
        ));
    }

    if wavefunction {
        // PySCF's own Molden writer, the route the user's `pyscfrun.py` uses.
        script.push_str(
            r#"
from pyscf.tools import molden
with open("orbitals.molden.input", "w") as handle:
    molden.header(mol, handle)
    molden.orbital_coeff(mol, handle, mf.mo_coeff, ene=mf.mo_energy, occ=mf.mo_occ)
"#,
        );
    }

    script.push_str(
        r#"
print("energy %.12f" % energy)
print("natoms %d" % mol.natm)
"#,
    );

    script
}

/// The energy the script printed, in Hartree.
///
/// Keyed by name rather than taken from a fixed line, so an extra line printed
/// by a future version cannot shift it.
pub fn parse_energy(stdout: &str) -> Option<f64> {
    stdout
        .lines()
        .find_map(|line| line.strip_prefix("energy "))
        .and_then(|rest| rest.trim().parse().ok())
}

/// The last thing Python said before giving up, for a message worth reading.
///
/// A traceback's final line is the exception and its text, which is the useful
/// part; everything above it is the call stack, which stays in the run
/// directory for when it is wanted.
pub fn failure_detail(stderr: &str) -> String {
    stderr
        .lines()
        .rev()
        .find(|line| !line.trim().is_empty())
        .unwrap_or("see xtb.stderr in the run directory")
        .trim()
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atoms() -> Vec<String> {
        ["O", "H", "H"].iter().map(|s| s.to_string()).collect()
    }

    fn water() -> Vec<f64> {
        vec![0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692]
    }

    fn method() -> MethodConfig {
        MethodConfig {
            functional: "b3lyp".to_string(),
            basis: "def2-svp".to_string(),
            nproc: 4,
            ..Default::default()
        }
    }

    #[test]
    fn the_script_carries_the_method_the_basis_and_the_threads() {
        let script = script_text(JobType::Optimize, 0, 1, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(script.contains("mf.xc = \"b3lyp\""), "{script}");
        assert!(script.contains("basis=\"def2-svp\""), "{script}");
        assert!(script.contains("lib.num_threads(4)"), "{script}");
    }

    /// PySCF counts unpaired electrons, where the panel asks for multiplicity.
    /// Passing the multiplicity straight through would make every triplet a
    /// quintet.
    #[test]
    fn multiplicity_is_converted_to_unpaired_electrons() {
        let singlet = script_text(JobType::Gradient, 0, 1, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(singlet.contains("spin=0"), "{singlet}");
        let triplet = script_text(JobType::Gradient, 0, 3, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(triplet.contains("spin=2"), "{triplet}");
        let anion = script_text(JobType::Gradient, -1, 2, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(anion.contains("charge=-1"), "{anion}");
        assert!(anion.contains("spin=1"), "{anion}");
    }

    /// Hartree-Fock builds a different PySCF object and takes no `xc`, so it
    /// cannot travel as a functional name.
    #[test]
    fn hartree_fock_builds_an_scf_object_and_sets_no_functional() {
        let mut hf = method();
        hf.functional = "hf".to_string();
        let script = script_text(JobType::Gradient, 0, 1, &hf, super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(script.contains("scf.RHF"), "{script}");
        assert!(
            !script.contains("mf.xc"),
            "Hartree-Fock takes no functional, and none should appear even unreached: {script}"
        );
        assert!(!script.contains("dft.RKS"), "{script}");

        let dft = script_text(JobType::Gradient, 0, 1, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(dft.contains("dft.RKS"), "{dft}");
        assert!(dft.contains("mf.xc = \"b3lyp\""), "{dft}");
        assert!(!dft.contains("scf.RHF"), "{dft}");
    }

    /// Only a job that optimises writes a trajectory, and only one that wants
    /// a Hessian computes one -- a Hessian is far the most expensive part, so
    /// asking for it by accident is costly.
    #[test]
    fn each_job_asks_for_exactly_what_it_needs() {
        let engrad = script_text(JobType::Gradient, 0, 1, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(!engrad.contains("geometric_solver"), "{engrad}");
        assert!(!engrad.contains("mf.Hessian()"), "{engrad}");

        let opt = script_text(JobType::Optimize, 0, 1, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(opt.contains("geometric_solver"), "{opt}");
        assert!(!opt.contains("mf.Hessian()"), "{opt}");

        let hess = script_text(JobType::Hessian, 0, 1, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(!hess.contains("geometric_solver"), "{hess}");
        assert!(hess.contains("mf.Hessian()"), "{hess}");

        let both = script_text(JobType::OptimizeThenFrequencies, 0, 1, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(both.contains("geometric_solver"), "{both}");
        assert!(both.contains("mf.Hessian()"), "{both}");
    }

    /// The Hessian is written behind a Turbomole header so Beavyr's existing
    /// reader takes it and no second format enters the codebase.
    #[test]
    fn the_hessian_is_written_in_a_format_beavyr_already_reads() {
        let script = script_text(JobType::Hessian, 0, 1, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(script.contains(r#"handle.write("$hessian\n")"#), "{script}");
    }

    /// The trajectory is written in Beavyr's own format, so a real one parses
    /// with the reader every other backend's trajectory uses.
    #[test]
    fn the_trajectory_format_is_the_one_beavyr_reads() {
        let script = script_text(JobType::Optimize, 0, 1, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(script.contains("cycle %d E = %.10f Eh"), "{script}");

        // What that produces, parsed by the shared reader.
        let produced = "3\ncycle 1 E = -75.3125218863 Eh\nO 0 0 0\nH 0 0 0\nH 0 0 0\n";
        let energies = super::super::behemoth::parse_trajectory_energies(produced);
        assert_eq!(energies.len(), 1);
        assert!((energies[0] + 75.312_521_886_3).abs() < 1.0e-9);
    }

    /// After an optimisation the energy must belong to the optimised geometry,
    /// not to the structure the optimiser started from.
    #[test]
    fn the_reported_energy_belongs_to_the_final_geometry() {
        let script = script_text(JobType::Optimize, 0, 1, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        let rebind = script.find("mol = optimised").expect("the molecule is replaced");
        let rebuild = script.find("mf = build(mol)\nenergy").expect("the energy is recomputed");
        assert!(rebind < rebuild, "the energy must be taken after the optimisation");
    }

    /// Method names are written as Python literals, so a name with a quote or
    /// a comma in it cannot become code.
    #[test]
    fn method_names_are_written_as_literals() {
        let mut awkward = method();
        awkward.functional = "0.2*HF + 0.8*B88, LYP".to_string();
        let script = script_text(JobType::Gradient, 0, 1, &awkward, super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(
            script.contains(r#"mf.xc = "0.2*HF + 0.8*B88, LYP""#),
            "{script}"
        );
    }

    /// A non-converged SCF stops the run rather than returning a number, and
    /// says so on standard error where the panel reads it.
    #[test]
    fn a_failed_scf_stops_the_script() {
        let script = script_text(JobType::Gradient, 0, 1, &method(), super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        assert!(script.contains("if not mf.converged"), "{script}");
        assert!(script.contains("sys.exit(2)"), "{script}");
        assert!(script.contains("isfinite"), "{script}");
    }

    /// The orbitals are written only when asked for.
    #[test]
    fn the_orbitals_are_written_only_on_request() {
        let without = script_text(
            JobType::SinglePoint, 0, 1, &method(),
            super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water(),
        );
        assert!(!without.contains("molden"), "{without}");

        let with = script_text(
            JobType::SinglePoint, 0, 1, &method(),
            super::super::job::ExcitedStateOptions::default(), true, &atoms(), &water(),
        );
        assert!(with.contains("molden.orbital_coeff"), "{with}");
        assert!(with.contains("orbitals.molden.input"), "{with}");
    }

    #[test]
    fn the_energy_is_read_by_key() {
        let stdout = "some chatter\nenergy -75.322774764508\nnatoms 3\n";
        let energy = parse_energy(stdout).expect("the energy is printed");
        assert!((energy + 75.322_774_764_508).abs() < 1.0e-12);
        assert_eq!(parse_energy("natoms 3\n"), None);
    }

    /// The exception is the useful line, not the top of the traceback.
    #[test]
    fn the_failure_detail_is_the_last_line_of_the_traceback() {
        let stderr = concat!(
            "Traceback (most recent call last):\n",
            "  File \"run_pyscf.py\", line 12, in <module>\n",
            "    mf.kernel()\n",
            "ValueError: basis not found for element Uuo\n"
        );
        assert_eq!(
            failure_detail(stderr),
            "ValueError: basis not found for element Uuo"
        );
        assert!(failure_detail("").contains("xtb.stderr"));
    }

    /// Generates a real script, runs it with the real PySCF, and checks what
    /// comes back. Ignored by default because it needs an environment that can
    /// import PySCF and takes a few seconds; run it with
    /// `cargo test -- --ignored pyscf_really_runs`.
    ///
    /// Worth having despite that: every other test here asserts on the text of
    /// a script, and text that looks right can still be a Python syntax error
    /// or call an API that has moved. This is the only test that would notice.
    #[test]
    #[ignore = "needs a Python that can import PySCF"]
    fn a_generated_script_really_runs() {
        use std::fs;

        if !super::super::python_env::probe("pyscf").is_ready() {
            eprintln!("no PySCF in this environment; skipping");
            return;
        }

        let dir = std::env::temp_dir().join(format!("beavyr_pyscf_e2e_{}", std::process::id()));
        fs::create_dir_all(&dir).expect("scratch directory");

        let mut small = method();
        small.basis = "sto-3g".to_string();
        small.nproc = 1;
        let script = script_text(JobType::OptimizeThenFrequencies, 0, 1, &small, super::super::job::ExcitedStateOptions::default(), false, &atoms(), &water());
        let files = job_files(&dir);
        fs::write(&files.script, script).expect("write the script");

        let output = run_command(&dir).output().expect("run python3");
        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            output.status.success(),
            "the generated script failed.\nstdout:\n{stdout}\nstderr:\n{stderr}"
        );

        // An energy, and a sane one for water in a minimal basis.
        let energy = parse_energy(&stdout).expect("an energy is printed");
        assert!(
            (-80.0..-70.0).contains(&energy),
            "water/STO-3G should be near -75 Eh, got {energy}"
        );

        // A trajectory in the shared format, read by the shared reader.
        let trajectory = fs::read_to_string(&files.trajectory).expect("a trajectory");
        let energies = super::super::behemoth::parse_trajectory_energies(&trajectory);
        assert!(!energies.is_empty(), "no optimiser steps recorded:\n{trajectory}");
        assert!(
            energies.last().unwrap() <= energies.first().unwrap(),
            "an optimisation should not end higher than it began: {energies:?}"
        );

        // A Hessian the existing Turbomole reader takes, square and symmetric.
        let hessian_text = fs::read_to_string(&files.hessian).expect("a Hessian");
        let hessian = super::super::hessian_file::parse_turbomole_hessian(&hessian_text, 3)
            .expect("the Hessian parses with Beavyr's own reader");
        assert_eq!(hessian.dim(), (9, 9));
        let asymmetry = (0..9)
            .flat_map(|i| (0..9).map(move |j| (i, j)))
            .map(|(i, j)| (hessian[[i, j]] - hessian[[j, i]]).abs())
            .fold(0.0f64, f64::max);
        assert!(asymmetry < 1.0e-6, "the Hessian is not symmetric: {asymmetry}");

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn the_command_runs_the_script_in_its_own_directory() {
        let command = run_command(Path::new("/tmp/run"));
        let args: Vec<String> = command
            .get_args()
            .map(|a| a.to_string_lossy().to_string())
            .collect();
        assert_eq!(args, vec!["run_pyscf.py"]);
        assert_eq!(command.get_current_dir(), Some(Path::new("/tmp/run")));
    }
}
