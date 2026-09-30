//! One energy and one Cartesian gradient, from whichever program was chosen.
//!
//! Three tools need exactly this and nothing else: the transition-state
//! generator, which walks a reaction path one gradient at a time; Beavyr's own
//! optimiser, which is the only route a backend without one has; and, when they
//! arrive, relaxed scans and constrained optimisations. Each of them used to
//! carry its own copy of "write the input, run it, read the numbers back", and
//! two copies had already drifted.
//!
//! So it lives here once, split into the three steps a caller needs to
//! interleave with its own cancellation:
//!
//! * [`prepare`] writes whatever input this program wants,
//! * [`command`] builds the process,
//! * [`read`] takes the energy and gradient back out.
//!
//! Spawning and waiting stay with the caller, because how a run is cancelled
//! differs between the tools and is not this module's business.
//!
//! Units throughout are Hartree and Hartree per bohr, which is what every
//! optimiser here wants, so nothing is converted on the way out.

use std::path::{Path, PathBuf};
use std::process::Command;

use super::job::{ExcitedStateOptions, JobType, QcMethod};
use super::method::MethodConfig;
use super::program::QcProgram;

const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;

/// The geometry file the command-line programs are given.
const GEOMETRY_FILE: &str = "geometry.xyz";

/// Why a program cannot be used for a gradient-driven job, or `None`.
///
/// DREIDING is the only refusal. It is a force field, and a reaction path or a
/// saddle search on one is not a result worth having; the other tools that
/// drive gradients are quantum-chemical by intent.
pub fn unsupported_reason(program: QcProgram) -> Option<&'static str> {
    match program {
        QcProgram::Dreiding => Some(
            "DREIDING is a generic force field. A reaction path on one is not worth the run; \
             choose a quantum-chemistry program.",
        ),
        _ => None,
    }
}

/// An XYZ file, from coordinates in bohr.
/// The structure as an XYZ file, for a program to read.
///
/// Dummies are left out, and the header counts what is written: they are
/// positions the user built against, and no program here accepts one.
pub fn xyz_text(atoms: &[String], coords_bohr: &[f64], comment: &str) -> String {
    let kept: Vec<(&String, &[f64])> = atoms
        .iter()
        .zip(coords_bohr.chunks_exact(3))
        .filter(|(symbol, _)| !crate::molecule::is_dummy(symbol))
        .collect();
    let mut text = format!("{}\n{comment}\n", kept.len());
    for (symbol, xyz) in kept {
        text.push_str(&format!(
            "{symbol:<3} {:16.10} {:16.10} {:16.10}\n",
            xyz[0] * BOHR_TO_ANGSTROM,
            xyz[1] * BOHR_TO_ANGSTROM,
            xyz[2] * BOHR_TO_ANGSTROM
        ));
    }
    text
}

/// Writes whatever input `program` needs for one energy and gradient.
///
/// Also clears the previous run's gradient file. That matters: several of these
/// programs write their gradient to a fixed name, and a run that fails after
/// the check but before writing would otherwise leave the last geometry's
/// gradient to be read as this one's -- a wrong answer that looks entirely
/// ordinary.
pub fn prepare(
    program: QcProgram,
    workdir: &Path,
    atoms: &[String],
    coords_bohr: &[f64],
    charge: i32,
    multiplicity: i32,
    method: &MethodConfig,
) -> Result<(), String> {
    if let Some(reason) = unsupported_reason(program) {
        return Err(reason.to_string());
    }
    std::fs::create_dir_all(workdir)
        .map_err(|e| format!("Could not create the run directory: {e}"))?;

    for stale in [gradient_path(program, workdir), energy_source(program, workdir)] {
        if stale.is_file() {
            let _ = std::fs::remove_file(&stale);
        }
    }

    let angstrom: Vec<f64> = coords_bohr.iter().map(|c| c * BOHR_TO_ANGSTROM).collect();

    match program {
        QcProgram::Xtb | QcProgram::Behemoth => {
            std::fs::write(
                workdir.join(GEOMETRY_FILE),
                xyz_text(atoms, coords_bohr, "Beavyr energy/gradient"),
            )
            .map_err(|e| format!("Could not write the geometry: {e}"))?;
        }
        QcProgram::Orca => {
            let text = super::orca_run::panel_input_text(
                JobType::Gradient,
                method_for(program, method),
                charge,
                multiplicity,
                method,
                ExcitedStateOptions::default(),
                &method.orca_extra,
                &super::orca_run::xyz_body(atoms, &angstrom),
            );
            std::fs::write(super::orca_run::job_files(workdir).input, text)
                .map_err(|e| format!("Could not write the ORCA input: {e}"))?;
        }
        QcProgram::Psi4 | QcProgram::Psi4Py => {
            let files = super::psi4_run::job_files(workdir);
            let (path, text) = if program == QcProgram::Psi4Py {
                (
                    files.script,
                    super::psi4_run::script_text(
                        JobType::Gradient,
                        method_for(program, method),
                        charge,
                        multiplicity,
                        method,
                        ExcitedStateOptions::default(),
                        false,
                        atoms,
                        &angstrom,
                    ),
                )
            } else {
                (
                    files.input,
                    super::psi4_run::input_text(
                        JobType::Gradient,
                        method_for(program, method),
                        charge,
                        multiplicity,
                        method,
                        ExcitedStateOptions::default(),
                        false,
                        atoms,
                        &angstrom,
                    ),
                )
            };
            std::fs::write(path, text)
                .map_err(|e| format!("Could not write the Psi4 job: {e}"))?;
        }
        QcProgram::PySCF => {
            let text = super::pyscf_run::script_text(
                JobType::Gradient,
                charge,
                multiplicity,
                method,
                ExcitedStateOptions::default(),
                false,
                atoms,
                &angstrom,
            );
            std::fs::write(super::pyscf_run::job_files(workdir).script, text)
                .map_err(|e| format!("Could not write the PySCF script: {e}"))?;
        }
        QcProgram::SparrowPy => {
            let files = super::sparrow_run::job_files(workdir);
            std::fs::write(
                &files.geometry,
                super::sparrow_run::geometry_text(atoms, &angstrom),
            )
            .map_err(|e| format!("Could not write the Sparrow geometry: {e}"))?;
            std::fs::write(
                &files.script,
                super::sparrow_run::script_text(
                    JobType::Gradient,
                    charge,
                    multiplicity,
                    method,
                    false,
                ),
            )
            .map_err(|e| format!("Could not write the Sparrow script: {e}"))?;
        }
        QcProgram::Dreiding => unreachable!("refused above"),
    }
    Ok(())
}

/// Which level of theory to send, for the programs that take one above the
/// functional.
///
/// The panel's own method choice is not carried into these tools, which ask
/// only for a surface to walk on, so DFT is the assumption -- the same one the
/// optimizer and frequency panels make.
fn method_for(program: QcProgram, _method: &MethodConfig) -> QcMethod {
    match program {
        QcProgram::SparrowPy => QcMethod::Semiempirical,
        _ => QcMethod::Dft,
    }
}

/// The process that computes it.
pub fn command(
    program: QcProgram,
    binary: &Path,
    workdir: &Path,
    charge: i32,
    multiplicity: i32,
    method: &MethodConfig,
) -> Result<Command, String> {
    if let Some(reason) = unsupported_reason(program) {
        return Err(reason.to_string());
    }
    if program.needs_binary_path() && binary.as_os_str().is_empty() {
        return Err(format!("Set the {} executable path first.", program.label()));
    }

    Ok(match program {
        QcProgram::Xtb => {
            let mut command = Command::new(binary);
            command
                .current_dir(workdir)
                .arg(GEOMETRY_FILE)
                .arg("--grad")
                .arg("--chrg")
                .arg(charge.to_string());
            // xTB counts unpaired electrons, not multiplicity.
            if multiplicity > 1 {
                command.arg("--uhf").arg((multiplicity - 1).to_string());
            }
            command
                .args(super::method::xtb_method_args(method))
                .arg("-P")
                .arg(method.nproc.max(1).to_string());
            command
        }
        QcProgram::Behemoth => super::behemoth::gradient_command(
            binary,
            workdir,
            GEOMETRY_FILE,
            charge,
            multiplicity,
            method,
        ),
        QcProgram::Orca => super::orca_run::run_command(binary, workdir),
        QcProgram::Psi4 | QcProgram::Psi4Py => super::psi4_run::run_command(
            if program == QcProgram::Psi4Py {
                super::psi4_run::Route::Python
            } else {
                super::psi4_run::Route::Executable
            },
            binary,
            workdir,
            method.nproc,
        ),
        QcProgram::PySCF => super::pyscf_run::run_command(workdir),
        QcProgram::SparrowPy => super::sparrow_run::run_command(workdir),
        QcProgram::Dreiding => unreachable!("refused above"),
    })
}

/// Where this program leaves its gradient.
fn gradient_path(program: QcProgram, workdir: &Path) -> PathBuf {
    match program {
        QcProgram::Orca => super::orca_run::job_files(workdir).engrad,
        QcProgram::Psi4 | QcProgram::Psi4Py => super::psi4_run::job_files(workdir).gradient,
        QcProgram::PySCF => super::pyscf_run::job_files(workdir).gradient,
        QcProgram::SparrowPy => super::sparrow_run::job_files(workdir).gradient,
        // xTB writes a Turbomole gradient beside its input; Behemoth prints to
        // standard output and writes no file.
        QcProgram::Xtb => workdir.join("gradient"),
        QcProgram::Behemoth | QcProgram::Dreiding => workdir.join("gradient"),
    }
}

/// Where this program leaves the text its energy is read from, when that is a
/// file of its own rather than what the caller captured.
fn energy_source(program: QcProgram, workdir: &Path) -> PathBuf {
    match program {
        QcProgram::Orca => super::orca_run::job_files(workdir).output,
        QcProgram::Psi4 | QcProgram::Psi4Py => super::psi4_run::job_files(workdir).output,
        _ => workdir.join("__no_such_file__"),
    }
}

/// Reads the energy and gradient back, in Hartree and Hartree per bohr.
///
/// `captured` is whatever the caller redirected the process's output to. Some
/// programs put the energy there and some in a file of their own, which is why
/// both are consulted.
pub fn read(
    program: QcProgram,
    workdir: &Path,
    natoms: usize,
    captured: &str,
) -> Result<(f64, Vec<f64>), String> {
    let gradient_file = gradient_path(program, workdir);

    let (energy, gradient) = match program {
        QcProgram::Xtb => (
            super::xtbrun::parse_xtb_energy(captured)?,
            super::xtbrun::parse_xtb_grad(&gradient_file, natoms)?,
        ),
        QcProgram::Behemoth => super::behemoth::parse_gradient(captured, natoms)?,
        QcProgram::Orca => {
            // ORCA's `.engrad` carries both, and Beavyr's reader for it also
            // checks the geometry it was computed at.
            let text = std::fs::read_to_string(&gradient_file)
                .map_err(|_| "ORCA wrote no .engrad file".to_string())?;
            let parsed = super::orca_engrad::parse_orca_engrad(&text)?;
            (parsed.energy_hartree, parsed.gradient_bohr)
        }
        QcProgram::Psi4 | QcProgram::Psi4Py => {
            let printed = std::fs::read_to_string(energy_source(program, workdir))
                .unwrap_or_else(|_| captured.to_string());
            let energy = super::psi4_run::parse_energy(&printed)
                .or_else(|| super::psi4_run::parse_energy(captured))
                .ok_or_else(|| "Psi4 printed no energy".to_string())?;
            let text = std::fs::read_to_string(&gradient_file)
                .map_err(|_| "Psi4 wrote no gradient".to_string())?;
            (
                energy,
                super::hessian_file::parse_turbomole_gradient(&text, natoms)?,
            )
        }
        QcProgram::PySCF => {
            let energy = super::pyscf_run::parse_energy(captured)
                .ok_or_else(|| "PySCF printed no energy".to_string())?;
            let text = std::fs::read_to_string(&gradient_file)
                .map_err(|_| "PySCF wrote no gradient".to_string())?;
            (
                energy,
                super::hessian_file::parse_turbomole_gradient(&text, natoms)?,
            )
        }
        QcProgram::SparrowPy => {
            let energy = super::sparrow_run::parse_energy(captured)
                .ok_or_else(|| "Sparrow printed no energy".to_string())?;
            let text = std::fs::read_to_string(&gradient_file)
                .map_err(|_| "Sparrow wrote no gradient".to_string())?;
            (
                energy,
                super::hessian_file::parse_turbomole_gradient(&text, natoms)?,
            )
        }
        QcProgram::Dreiding => return Err(unsupported_reason(program).unwrap().to_string()),
    };

    // A gradient of the wrong length, or one with a non-finite value in it,
    // would send an optimiser somewhere meaningless rather than fail. Caught
    // here so every caller gets the check.
    if !energy.is_finite() {
        return Err(format!("{} returned a non-finite energy", program.label()));
    }
    if gradient.len() != natoms * 3 {
        return Err(format!(
            "{} returned {} gradient components for {natoms} atoms",
            program.label(),
            gradient.len()
        ));
    }
    if gradient.iter().any(|v| !v.is_finite()) {
        return Err(format!(
            "{} returned a gradient with a non-finite value in it",
            program.label()
        ));
    }
    Ok((energy, gradient))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atoms() -> Vec<String> {
        ["O", "H", "H"].iter().map(|s| s.to_string()).collect()
    }

    /// Bohr in, angstrom out. Getting this backwards would scale every
    /// structure by a factor of two and still look like a molecule.
    #[test]
    fn the_geometry_is_written_in_angstrom() {
        let text = xyz_text(&atoms(), &[0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0], "x");
        let second: Vec<&str> = text.lines().nth(2).unwrap().split_whitespace().collect();
        let z: f64 = second[3].parse().unwrap();
        assert!((z - BOHR_TO_ANGSTROM).abs() < 1.0e-9, "one bohr is 0.529 Å, got {z}");
    }

    /// Every program but the force field can provide a surface to walk on.
    #[test]
    fn only_the_force_field_is_refused() {
        for program in QcProgram::ALL {
            let refused = unsupported_reason(program).is_some();
            assert_eq!(refused, program == QcProgram::Dreiding, "{program:?}");
        }
        let reason = unsupported_reason(QcProgram::Dreiding).unwrap();
        assert!(reason.contains("force field"), "{reason}");
        assert!(reason.contains("quantum-chemistry"), "the alternative: {reason}");
    }

    /// A command cannot be built without a path for a program that needs one,
    /// and must be buildable without one for a program that does not.
    #[test]
    fn a_missing_executable_is_caught_before_spawning() {
        let config = MethodConfig::default();
        for program in [QcProgram::Xtb, QcProgram::Orca, QcProgram::Psi4] {
            assert!(
                command(program, Path::new(""), Path::new("/tmp"), 0, 1, &config).is_err(),
                "{program:?} needs a path"
            );
        }
        for program in [QcProgram::PySCF, QcProgram::Psi4Py, QcProgram::SparrowPy] {
            assert!(
                command(program, Path::new(""), Path::new("/tmp"), 0, 1, &config).is_ok(),
                "{program:?} is reached through python3"
            );
        }
    }

    /// The charge and the spin state reach the command, or every ion and every
    /// radical would silently be computed as a neutral singlet.
    #[test]
    fn the_charge_and_spin_state_reach_the_command() {
        let config = MethodConfig::default();
        let cmd = command(
            QcProgram::Xtb,
            Path::new("/bin/xtb"),
            Path::new("/tmp"),
            -1,
            3,
            &config,
        )
        .expect("xTB with a path");
        let args: Vec<String> = cmd.get_args().map(|a| a.to_string_lossy().into()).collect();
        assert!(args.contains(&"-1".to_string()), "{args:?}");
        // Two unpaired electrons for a triplet, not three.
        assert!(args.contains(&"--uhf".to_string()), "{args:?}");
        assert!(args.contains(&"2".to_string()), "{args:?}");
    }

    /// Each program's gradient is looked for where that program leaves it. Two
    /// of them sharing a path would have one read the other's.
    #[test]
    fn each_program_has_its_own_gradient_file() {
        let dir = Path::new("/tmp/run");
        let mut paths: Vec<PathBuf> = [
            QcProgram::Orca,
            QcProgram::Psi4,
            QcProgram::PySCF,
            QcProgram::SparrowPy,
        ]
        .into_iter()
        .map(|p| gradient_path(p, dir))
        .collect();
        let count = paths.len();
        paths.sort();
        paths.dedup();
        assert_eq!(paths.len(), count, "two programs share a gradient file");
        // ORCA's is its own format, not a Turbomole one.
        assert!(gradient_path(QcProgram::Orca, dir).ends_with("orca_job.engrad"));
        // xTB's name is xTB's own and cannot be changed; the three Beavyr
        // generates are named apart from it and from each other, so a stale
        // file from one run can never be read as another program's.
        assert!(gradient_path(QcProgram::Xtb, dir).ends_with("gradient"));
        for program in [QcProgram::Psi4, QcProgram::PySCF, QcProgram::SparrowPy] {
            assert!(
                !gradient_path(program, dir).ends_with("gradient"),
                "{program:?} must not share xTB's fixed name"
            );
        }
    }

    /// A short or non-finite gradient is refused rather than handed to an
    /// optimiser, which would follow it somewhere meaningless and report
    /// convergence.
    #[test]
    fn a_malformed_gradient_is_refused() {
        let dir = std::env::temp_dir().join(format!("beavyr_grad_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        // Behemoth's reader is the one that takes everything from the captured
        // text, so it can be exercised without writing files.
        let short = concat!(
            "Analytical TASI nuclear gradient\n",
            "Units: Hartree/bohr\n",
            "Energy:     -68.4051345926 Hartree\n",
            "     1     -0.0020104     -0.0021028     -0.0000000\n"
        );
        assert!(read(QcProgram::Behemoth, &dir, 3, short).is_err());

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Sparrow is semi-empirical and takes no functional; everything else here
    /// is asked for as DFT, the same assumption the other panels make.
    #[test]
    fn the_level_of_theory_matches_the_program() {
        let config = MethodConfig::default();
        assert_eq!(
            method_for(QcProgram::SparrowPy, &config),
            QcMethod::Semiempirical
        );
        for program in [QcProgram::Orca, QcProgram::Psi4, QcProgram::PySCF] {
            assert_eq!(method_for(program, &config), QcMethod::Dft, "{program:?}");
        }
    }
}
