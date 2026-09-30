//! Driving Psi4, by its own executable or through Python.
//!
//! Psi4 is the one backend with two genuine routes, and they share almost
//! everything. A Psi4 input file is Python: the sugar blocks in its
//! documentation are a preprocessor convenience, and plain `psi4.` API calls
//! pass straight through. Checked directly against the installed 1.11.
//!
//! So Beavyr generates one body and wraps it twice:
//!
//! * the **executable** route writes it as `input.dat` and runs
//!   `psi4 -i input.dat -o output.dat`. That binary carries its own
//!   environment -- it reports its version with `PYTHONPATH` and
//!   `CONDA_PREFIX` cleared -- so this route works whatever environment Beavyr
//!   was launched from. It is the one to reach for when Psi4 lives in an
//!   environment other than the one holding PySCF.
//! * the **Python** route prepends `import psi4` and runs it with `python3`
//!   from Beavyr's own environment, which therefore has to be the one with
//!   Psi4 in it.
//!
//! Results are written in the formats Beavyr already reads, as the PySCF
//! script does: the trajectory in Beavyr's own `cycle <n> E = <v> Eh` layout,
//! the Hessian behind a `$hessian` header, the gradient as a Turbomole `$grad`
//! block, and the geometry as an ordinary XYZ.

use std::path::{Path, PathBuf};
use std::process::Command;

use super::job::{ExcitedStateOptions, JobType, QcMethod};
use super::method::MethodConfig;
use super::psi4_method;
use super::python_env::PYTHON;

/// The input file the executable route runs.
pub const INPUT_FILE: &str = "input.dat";
/// The script the Python route runs.
pub const SCRIPT_FILE: &str = "run_psi4.py";
/// Where Psi4's own output goes.
pub const OUTPUT_FILE: &str = "psi4.out";
/// The optimisation path, already in Beavyr's trajectory format.
pub const TRAJECTORY_FILE: &str = "trajectory.xyz";
/// The final geometry, as an ordinary XYZ.
pub const GEOMETRY_FILE: &str = "final.xyz";
/// The Hessian, behind a `$hessian` header.
pub const HESSIAN_FILE: &str = "hessian.txt";
/// The gradient, as a Turbomole `$grad` block.
pub const GRADIENT_FILE: &str = "psi4.grad";
/// The excited states, written by the script in a format of our own.
pub const SPECTRUM_FILE: &str = "spectrum.txt";
/// The molecular orbitals, when asked for.
pub const MOLDEN_FILE: &str = "orbitals.molden";

/// Which of the two routes to take.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Route {
    /// `psi4 -i input.dat`, carrying its own environment.
    Executable,
    /// `python3 run_psi4.py`, using Beavyr's environment.
    Python,
}

/// Where each file of a run in `run_dir` lives.
#[derive(Debug, Clone)]
pub struct JobFiles {
    pub input: PathBuf,
    pub script: PathBuf,
    pub output: PathBuf,
    pub trajectory: PathBuf,
    pub geometry: PathBuf,
    pub hessian: PathBuf,
    pub gradient: PathBuf,
    pub spectrum: PathBuf,
    pub molden: PathBuf,
}

pub fn job_files(run_dir: &Path) -> JobFiles {
    JobFiles {
        input: run_dir.join(INPUT_FILE),
        script: run_dir.join(SCRIPT_FILE),
        output: run_dir.join(OUTPUT_FILE),
        trajectory: run_dir.join(TRAJECTORY_FILE),
        geometry: run_dir.join(GEOMETRY_FILE),
        hessian: run_dir.join(HESSIAN_FILE),
        gradient: run_dir.join(GRADIENT_FILE),
        spectrum: run_dir.join(SPECTRUM_FILE),
        molden: run_dir.join(MOLDEN_FILE),
    }
}

/// The command for a prepared run.
pub fn run_command(route: Route, binary: &Path, run_dir: &Path, nproc: u32) -> Command {
    match route {
        Route::Executable => {
            let mut command = Command::new(binary);
            command
                .current_dir(run_dir)
                .arg("-i")
                .arg(INPUT_FILE)
                .arg("-o")
                .arg(OUTPUT_FILE)
                .arg("-n")
                .arg(nproc.max(1).to_string());
            command
        }
        Route::Python => {
            let mut command = Command::new(PYTHON);
            command.current_dir(run_dir).arg(SCRIPT_FILE);
            command
        }
    }
}

/// The method name Psi4 is given.
///
/// Hartree-Fock is `scf` here, not `hf`; a dispersion correction is part of the
/// functional's own name rather than a separate keyword, which is why there is
/// no dispersion control for this program.
pub fn method_name<'a>(method: QcMethod, functional: &'a str) -> &'a str {
    match method {
        QcMethod::Dft => functional,
        QcMethod::Hf => "scf",
        QcMethod::Mp2 => "mp2",
        QcMethod::RiMp2 => "mp2",
        QcMethod::Mp2F12 => "mp2-f12",
        QcMethod::CcsdT => "ccsd(t)",
        QcMethod::DlpnoCcsdT => "ccsd(t)",
        QcMethod::CcsdTF12 => "ccsd(t)-f12",
        QcMethod::Gfn1
        | QcMethod::Gfn2
        | QcMethod::GfnFf
        | QcMethod::Tasi
        | QcMethod::Semiempirical => "scf",
    }
}

/// The geometry block, one atom per line, with the frame pinned.
///
/// `no_com` and `no_reorient` matter and are not decoration. Psi4 would
/// otherwise move the molecule to its centre of mass and rotate it onto its
/// principal axes before every calculation, which is a geometry-dependent
/// transformation. A Hessian computed in a frame Beavyr does not know about
/// cannot be projected correctly, and the frequencies that come out look
/// plausible and are wrong -- the same trap that made xTB's reoriented
/// geometry worth chasing down. The Smite interface pins the frame for the
/// same reason.
fn geometry_block(charge: i32, multiplicity: i32, atoms: &[String], positions: &[f64]) -> String {
    let mut block = format!("{charge} {}\\n", multiplicity.max(1));
    for (index, symbol) in atoms.iter().enumerate() {
        // A position the user built against, not an atom. Psi4 would stop on
        // it, and it contributes no electrons to the charge and multiplicity
        // written above.
        if crate::molecule::is_dummy(symbol) {
            continue;
        }
        let base = index * 3;
        block.push_str(&format!(
            "{symbol} {:.10} {:.10} {:.10}\\n",
            positions[base],
            positions[base + 1],
            positions[base + 2]
        ));
    }
    block.push_str("units angstrom\\nno_com\\nno_reorient\\n");
    block
}

/// Builds the body both routes share.
#[allow(clippy::too_many_arguments)]
pub fn body_text(
    job: JobType,
    method: QcMethod,
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    excited: ExcitedStateOptions,
    wavefunction: bool,
    atoms: &[String],
    positions_angstrom: &[f64],
) -> String {
    let name = method_name(method, &config.functional);
    let geometry = geometry_block(charge, multiplicity, atoms, positions_angstrom);
    // Psi4 asks for the reference by name, and an open shell must be
    // unrestricted or the calculation describes a different state.
    let reference = if multiplicity > 1 { "uks" } else { "rks" };
    let reference = if psi4_method::is_hartree_fock(name) {
        if multiplicity > 1 { "uhf" } else { "rhf" }
    } else {
        reference
    };
    let basis = &config.basis;
    let memory = config.memory_mb.max(256);
    let nproc = config.nproc.max(1);

    let mut text = String::new();
    text.push_str(&format!(
        r#"# Written by Beavyr. One Psi4 job, results left beside this file.
psi4.set_memory("{memory} MB")
psi4.set_num_threads({nproc})

ATOMS = {atoms:?}
BOHR = 0.52917721092

mol = psi4.geometry("""
{geometry}""")

psi4.set_options({{
    "basis": {basis:?},
    "reference": {reference:?},
}})

METHOD = {name:?}

def write_xyz(path, molecule, comment):
    coords = molecule.geometry().to_array() * BOHR
    with open(path, "w") as handle:
        handle.write("%d\n%s\n" % (len(ATOMS), comment))
        for symbol, xyz in zip(ATOMS, coords):
            handle.write("%-3s %18.10f %18.10f %18.10f\n" % (symbol, xyz[0], xyz[1], xyz[2]))

"#
    ));

    // A composite brings its own basis, so the option set above is overridden.
    if psi4_method::is_composite(&config.functional) && method == QcMethod::Dft {
        text.push_str("# A composite defines its own basis; whatever was set is not used.\n");
    }

    match job {
        JobType::SinglePoint => text.push_str(
            "energy, wfn = psi4.energy(METHOD, return_wfn=True)\n             write_xyz(\"final.xyz\", mol, \"Psi4 geometry\")\n",
        ),
        JobType::Gradient => text.push_str(
            r#"grad, wfn = psi4.gradient(METHOD, return_wfn=True)
energy = wfn.energy()
g = grad.to_array()
coords = mol.geometry().to_array()
# Turbomole `$grad`, which Beavyr already reads: the coordinates first, then
# one gradient line per atom, both in bohr.
with open("psi4.grad", "w") as handle:
    handle.write("$grad          cartesian gradients\n")
    handle.write("  cycle =      1    SCF energy =  %20.14f   |dE/dxyz| =  0.000000\n" % energy)
    for symbol, xyz in zip(ATOMS, coords):
        handle.write("%22.14f %22.14f %22.14f      %s\n" % (xyz[0], xyz[1], xyz[2], symbol))
    for row in g:
        handle.write("%22.14e %22.14e %22.14e\n" % (row[0], row[1], row[2]))
    handle.write("$end\n")
write_xyz("final.xyz", mol, "Psi4 geometry")
"#,
        ),
        JobType::Hessian | JobType::Frequencies => text.push_str(
            r#"hess, wfn = psi4.hessian(METHOD, return_wfn=True)
energy = wfn.energy()
H = hess.to_array()
# Behind a `$hessian` header, which is Turbomole's layout and what Beavyr's
# existing reader takes.
with open("hessian.txt", "w") as handle:
    handle.write("$hessian\n")
    for row in H:
        for value in row:
            handle.write(" %22.14e" % value)
        handle.write("\n")
    handle.write("$end\n")
write_xyz("final.xyz", mol, "Psi4 geometry")
"#,
        ),
        JobType::Optimize
        | JobType::OptimizeThenFrequencies
        | JobType::TransitionState
        | JobType::Irc => {
            if job == JobType::TransitionState {
                // A saddle search, with a real Hessian to start from and again
                // as it goes -- without it the search wanders to whatever
                // saddle is nearest rather than the one aimed at.
                text.push_str(
                    "psi4.set_options({\"opt_type\": \"ts\", \"full_hess_every\": 5})\n",
                );
            } else if job == JobType::Irc {
                text.push_str("psi4.set_options({\"opt_type\": \"irc\"})\n");
            }
            text.push_str(
                r#"frames = []

def record(**kwargs):
    """One frame per optimiser step."""
    try:
        coords = mol.geometry().to_array() * BOHR
        frames.append((coords, float(psi4.variable("CURRENT ENERGY"))))
    except Exception:
        pass

energy = psi4.optimize(METHOD, opt_coordinates="both")
record()

with open("trajectory.xyz", "w") as handle:
    if not frames:
        frames = [(mol.geometry().to_array() * BOHR, energy)]
    for index, (coords, step_energy) in enumerate(frames, start=1):
        handle.write("%d\ncycle %d E = %.10f Eh\n" % (len(ATOMS), index, step_energy))
        for symbol, xyz in zip(ATOMS, coords):
            handle.write("%-3s %18.10f %18.10f %18.10f\n" % (symbol, xyz[0], xyz[1], xyz[2]))

write_xyz("final.xyz", mol, "Psi4 optimized geometry")
"#,
            );
            if job == JobType::OptimizeThenFrequencies {
                text.push_str(
                    r#"
hess, wfn = psi4.hessian(METHOD, return_wfn=True)
H = hess.to_array()
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
        }
        JobType::TdDft => {
            text.push_str(&format!(
                r##"from psi4.driver.procrouting.response.scf_response import tdscf_excitations

energy, wfn = psi4.energy(METHOD, return_wfn=True)
roots = tdscf_excitations(wfn, states={states}, tda={tda})

# Written in a small format of our own: the script is ours, so there is no
# reason to print a table and read it back with a pattern.
with open("spectrum.txt", "w") as handle:
    handle.write("# root  energy_eV  oscillator_strength\n")
    for index, root in enumerate(roots, start=1):
        ev = float(root["EXCITATION ENERGY"]) * 27.211386245988
        try:
            f = float(root["OSCILLATOR STRENGTH (LEN)"])
        except Exception:
            f = 0.0
        handle.write("%d %.10f %.10f\n" % (index, ev, f))
write_xyz("final.xyz", mol, "Psi4 geometry")
"##,
                states = excited.states.max(1),
                tda = if excited.tda { "True" } else { "False" }
            ));
        }
    }

    if wavefunction {
        // Psi4's own Molden writer, the route the user's `psi4run.py` uses.
        // It needs the wavefunction object, which is why every job above asks
        // for one back rather than the energy alone.
        text.push_str(
            r#"
try:
    psi4.molden(wfn, "orbitals.molden")
except NameError:
    _e, wfn = psi4.energy(METHOD, return_wfn=True)
    psi4.molden(wfn, "orbitals.molden")
"#,
        );
    }

    text.push_str(
        r#"
psi4.core.print_out("BEAVYR_ENERGY %.12f\n" % energy)
print("energy %.12f" % energy)
"#,
    );

    text
}

/// The executable route's input file: the shared body, unchanged.
///
/// Psi4's preprocessor passes plain Python through and has `psi4` already in
/// scope, so nothing has to be added.
pub fn input_text(
    job: JobType,
    method: QcMethod,
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    excited: ExcitedStateOptions,
    wavefunction: bool,
    atoms: &[String],
    positions_angstrom: &[f64],
) -> String {
    body_text(
        job, method, charge, multiplicity, config, excited, wavefunction, atoms,
        positions_angstrom,
    )
}

/// The Python route's script: the same body, with the import it needs.
pub fn script_text(
    job: JobType,
    method: QcMethod,
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    excited: ExcitedStateOptions,
    wavefunction: bool,
    atoms: &[String],
    positions_angstrom: &[f64],
) -> String {
    let mut text = String::from("import psi4\npsi4.core.set_output_file(\"psi4.out\", False)\n");
    text.push_str(&body_text(
        job, method, charge, multiplicity, config, excited, wavefunction, atoms,
        positions_angstrom,
    ));
    text
}

/// The energy the run printed, in Hartree.
pub fn parse_energy(text: &str) -> Option<f64> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix("energy ")
                .or_else(|| line.strip_prefix("BEAVYR_ENERGY "))
        })
        .next_back()
        .and_then(|rest| rest.trim().parse().ok())
}

/// The last thing said before giving up, for a message worth reading.
pub fn failure_detail(stderr: &str, output: &str) -> String {
    for text in [stderr, output] {
        if let Some(line) = text
            .lines()
            .rev()
            .find(|line| !line.trim().is_empty() && !line.trim().starts_with('*'))
        {
            if !line.trim().is_empty() {
                return line.trim().to_string();
            }
        }
    }
    "see psi4.out in the run directory".to_string()
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

    fn config() -> MethodConfig {
        MethodConfig {
            functional: "b3lyp".to_string(),
            basis: "def2-SVP".to_string(),
            nproc: 4,
            memory_mb: 4096,
            ..Default::default()
        }
    }

    fn build(job: JobType, method: QcMethod, multiplicity: i32) -> String {
        body_text(
            job,
            method,
            0,
            multiplicity,
            &config(),
            ExcitedStateOptions::default(),
            false,
            &atoms(),
            &water(),
        )
    }

    /// Psi4 would otherwise translate to the centre of mass and rotate onto
    /// its principal axes before every calculation. A Hessian built in a frame
    /// Beavyr does not know about cannot be projected correctly, and the
    /// frequencies come out plausible and wrong.
    #[test]
    fn the_laboratory_frame_is_pinned() {
        let text = build(JobType::Frequencies, QcMethod::Dft, 1);
        assert!(text.contains("no_com"), "{text}");
        assert!(text.contains("no_reorient"), "{text}");
    }

    /// Hartree-Fock is `scf` to Psi4, and its reference is not a Kohn-Sham one.
    #[test]
    fn hartree_fock_uses_psi4s_own_names() {
        let text = build(JobType::SinglePoint, QcMethod::Hf, 1);
        assert!(text.contains(r#"METHOD = "scf""#), "{text}");
        assert!(text.contains(r#""reference": "rhf""#), "{text}");
    }

    /// An open shell must be unrestricted, or the calculation describes a
    /// different electronic state from the one asked for.
    #[test]
    fn an_open_shell_is_unrestricted() {
        let dft = build(JobType::SinglePoint, QcMethod::Dft, 3);
        assert!(dft.contains(r#""reference": "uks""#), "{dft}");
        let hf = build(JobType::SinglePoint, QcMethod::Hf, 2);
        assert!(hf.contains(r#""reference": "uhf""#), "{hf}");
    }

    /// Each job calls the right Psi4 driver and nothing else. A Hessian is by
    /// far the most expensive part, so asking for one by accident is costly.
    #[test]
    fn each_job_calls_the_right_driver() {
        let sp = build(JobType::SinglePoint, QcMethod::Dft, 1);
        // The wavefunction comes back too, because Psi4's Molden writer needs
        // that object and asking for it costs nothing.
        assert!(sp.contains("psi4.energy(METHOD, return_wfn=True)"), "{sp}");
        assert!(!sp.contains("psi4.hessian"), "{sp}");
        assert!(!sp.contains("psi4.optimize"), "{sp}");

        let grad = build(JobType::Gradient, QcMethod::Dft, 1);
        assert!(grad.contains("psi4.gradient(METHOD"), "{grad}");
        assert!(grad.contains("$grad"), "written where Beavyr reads it: {grad}");

        let hess = build(JobType::Hessian, QcMethod::Dft, 1);
        assert!(hess.contains("psi4.hessian(METHOD"), "{hess}");
        assert!(hess.contains("$hessian"), "{hess}");

        let opt = build(JobType::Optimize, QcMethod::Dft, 1);
        assert!(opt.contains("psi4.optimize(METHOD"), "{opt}");
        assert!(!opt.contains("psi4.hessian"), "{opt}");

        let both = build(JobType::OptimizeThenFrequencies, QcMethod::Dft, 1);
        assert!(both.contains("psi4.optimize(METHOD"), "{both}");
        assert!(both.contains("psi4.hessian(METHOD"), "{both}");
    }

    /// A saddle search and a reaction path are ordinary optimisations with the
    /// optimiser told what to look for.
    #[test]
    fn the_saddle_and_path_searches_set_the_optimiser_mode() {
        let ts = build(JobType::TransitionState, QcMethod::Dft, 1);
        assert!(ts.contains(r#""opt_type": "ts""#), "{ts}");
        assert!(ts.contains("full_hess_every"), "a saddle search needs Hessians: {ts}");

        let irc = build(JobType::Irc, QcMethod::Dft, 1);
        assert!(irc.contains(r#""opt_type": "irc""#), "{irc}");
    }

    /// Excited states come from the entry point Psi4 actually has. The obvious
    /// guess, `psi4.tdscf_excitations`, does not exist in this build.
    #[test]
    fn tddft_uses_the_entry_point_that_exists() {
        let text = build(JobType::TdDft, QcMethod::Dft, 1);
        assert!(
            text.contains("from psi4.driver.procrouting.response.scf_response import tdscf_excitations"),
            "{text}"
        );
        assert!(text.contains("states=10"), "{text}");
        assert!(text.contains("tda=True"), "{text}");
    }

    /// The trajectory is written in the one format every backend shares, so
    /// the player and the energy plot read it unchanged.
    #[test]
    fn the_trajectory_format_is_the_shared_one() {
        let text = build(JobType::Optimize, QcMethod::Dft, 1);
        assert!(text.contains("cycle %d E = %.10f Eh"), "{text}");
    }

    /// The two routes differ by the import and nothing else, which is the whole
    /// reason one body serves both.
    #[test]
    fn the_two_routes_share_their_body() {
        let body = build(JobType::SinglePoint, QcMethod::Dft, 1);
        let script = script_text(
            JobType::SinglePoint,
            QcMethod::Dft,
            0,
            1,
            &config(),
            ExcitedStateOptions::default(),
            false,
            &atoms(),
            &water(),
        );
        let input = input_text(
            JobType::SinglePoint,
            QcMethod::Dft,
            0,
            1,
            &config(),
            ExcitedStateOptions::default(),
            false,
            &atoms(),
            &water(),
        );
        assert_eq!(input, body, "the executable route is the body unchanged");
        assert!(script.ends_with(&body), "the Python route is the body plus a header");
        assert!(script.starts_with("import psi4"), "{script}");
    }

    /// The executable is told where to put its output and how many threads to
    /// use; the Python route needs neither on its command line.
    #[test]
    fn each_route_is_launched_the_way_it_expects() {
        let exe = run_command(Route::Executable, Path::new("/x/psi4"), Path::new("/tmp/r"), 4);
        let args: Vec<String> = exe.get_args().map(|a| a.to_string_lossy().into()).collect();
        assert_eq!(args, vec!["-i", "input.dat", "-o", "psi4.out", "-n", "4"]);

        let py = run_command(Route::Python, Path::new(""), Path::new("/tmp/r"), 4);
        let args: Vec<String> = py.get_args().map(|a| a.to_string_lossy().into()).collect();
        assert_eq!(args, vec!["run_psi4.py"]);
    }

    /// Runs a real Psi4 single point by its executable, on water in a minimal
    /// basis. Ignored by default; run with
    /// `cargo test -- --ignored psi4_really_runs`.
    ///
    /// Every other test here asserts on generated text, and text that reads
    /// correctly can still be rejected by the program. This is the only one
    /// that would notice.
    #[test]
    #[ignore = "needs Psi4 installed"]
    fn a_generated_input_really_runs() {
        use std::fs;

        let binary = std::env::var("BEAVYR_PSI4_PATH").unwrap_or_else(|_| {
            "/home/peter/micromamba/envs/psi4env/bin/psi4".to_string()
        });
        if !Path::new(&binary).is_file() {
            eprintln!("no Psi4 at {binary}; skipping");
            return;
        }

        let dir = std::env::temp_dir().join(format!("beavyr_psi4_e2e_{}", std::process::id()));
        fs::create_dir_all(&dir).expect("scratch directory");

        let mut small = config();
        small.basis = "sto-3g".to_string();
        small.nproc = 1;
        small.memory_mb = 1024;

        let files = job_files(&dir);
        let text = input_text(
            JobType::Gradient,
            QcMethod::Dft,
            0,
            1,
            &small,
            ExcitedStateOptions::default(),
            false,
            &atoms(),
            &water(),
        );
        fs::write(&files.input, &text).expect("write the input");

        let status = run_command(Route::Executable, Path::new(&binary), &dir, 1)
            .status()
            .expect("run Psi4");
        let output = fs::read_to_string(&files.output).unwrap_or_default();
        assert!(status.success(), "Psi4 failed:\n{output}");

        // An energy, and a sane one for water in a minimal basis.
        let energy = parse_energy(&output)
            .unwrap_or_else(|| panic!("no energy printed:\n{output}"));
        assert!(
            (-80.0..-70.0).contains(&energy),
            "water/STO-3G should be near -75 Eh, got {energy}"
        );

        // The gradient, in the Turbomole layout Beavyr already reads.
        let gradient_text = fs::read_to_string(&files.gradient).expect("a gradient file");
        let gradient =
            super::super::hessian_file::parse_turbomole_gradient(&gradient_text, 3)
                .expect("Beavyr's own reader takes it");
        assert_eq!(gradient.len(), 9);
        assert!(gradient.iter().all(|v| v.is_finite()));

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn the_energy_is_read_by_key() {
        assert!(parse_energy("energy -76.123456789012\n").is_some());
        assert!(parse_energy("BEAVYR_ENERGY -76.123456789012\n").is_some());
        assert_eq!(parse_energy("nothing here\n"), None);
    }

    /// A composite brings its own basis, and the body says so rather than
    /// silently setting one that is then ignored.
    #[test]
    fn a_composite_says_it_brings_its_own_basis() {
        let mut composite = config();
        composite.functional = "r2scan-3c".to_string();
        let text = body_text(
            JobType::SinglePoint,
            QcMethod::Dft,
            0,
            1,
            &composite,
            ExcitedStateOptions::default(),
            false,
            &atoms(),
            &water(),
        );
        assert!(text.contains("composite defines its own basis"), "{text}");
    }

    /// Resources reach the run.
    #[test]
    fn the_resources_reach_the_run() {
        let text = build(JobType::SinglePoint, QcMethod::Dft, 1);
        assert!(text.contains(r#"psi4.set_memory("4096 MB")"#), "{text}");
        assert!(text.contains("psi4.set_num_threads(4)"), "{text}");
    }
}
