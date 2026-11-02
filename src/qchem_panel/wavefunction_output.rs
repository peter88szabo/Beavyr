//! Wavefunction output and matching-geometry handoff to Surface Tools.
use super::*;
use crate::orbitals::molden::{parse_molden, OrbitalData};
pub fn read(path: &Path, program: QcProgram) -> Result<(OrbitalData, String), String> {
    let mut text = fs::read_to_string(path)
        .map_err(|e| format!("Cannot open wavefunction {}: {e}", path.display()))?;
    if program == QcProgram::Xtb && !text.contains("Beavyr source: xTB") {
        text = text.replacen(
            "[Molden Format]",
            "[Molden Format]\n[Title]\nBeavyr source: xTB (primitive-normalized coefficients)",
            1,
        );
    }
    let data = parse_molden(&text)
        .map_err(|e| format!("Cannot read wavefunction {}: {e}", path.display()))?;
    Ok((data, text))
}
pub fn matches(data: &OrbitalData, atoms: &[String], positions: &[Vec3]) -> bool {
    data.atoms == atoms
        && data.positions.len() == positions.len()
        && data
            .positions
            .iter()
            .zip(positions)
            .all(|(a, b)| a.distance(*b) <= 1e-3)
}
pub fn request(command: &mut std::process::Command, program: QcProgram) {
    command.arg(if program == QcProgram::Xtb {
        "--molden"
    } else {
        "--wf2molden"
    });
}
pub fn produced(program: QcProgram, dir: &Path) -> Result<PathBuf, String> {
    let path = dir.join(if program == QcProgram::Xtb {
        "molden.input"
    } else {
        "canonicalMO.molden"
    });
    if path.is_file() {
        Ok(path)
    } else {
        Err(format!("{} completed, but did not produce the requested wavefunction {}. Calculation files are retained.",program.label(),path.display()))
    }
}
#[allow(clippy::too_many_arguments)]
pub fn final_xtb_family(
    program: QcProgram,
    binary: &Path,
    workdir: &Path,
    geometry: &str,
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    cancel: &AtomicBool,
    child: &Mutex<Option<Child>>,
) -> Result<PathBuf, String> {
    let (_, atoms, positions) = crate::molecule::parse_xyz_angstrom(geometry);
    let bohr: Vec<f64> = positions
        .into_iter()
        .flatten()
        .map(|x| x / crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM)
        .collect();
    let dir = workdir.join("final_orbitals");
    crate::qchem_interfaces::gradient::prepare(
        program,
        &dir,
        &atoms,
        &bohr,
        charge,
        multiplicity,
        config,
    )?;
    let mut command = crate::qchem_interfaces::gradient::command(
        program,
        binary,
        &dir,
        charge,
        multiplicity,
        config,
    )?;
    request(&mut command, program);
    command.stdout(Stdio::from(
        fs::File::create(dir.join("program.stdout")).map_err(|e| e.to_string())?,
    ));
    command.stderr(Stdio::from(
        fs::File::create(dir.join("program.stderr")).map_err(|e| e.to_string())?,
    ));
    let status = spawn_and_wait(command, program.label(), cancel, child)?;
    if !status.success() {
        return Err(format!(
            "{} final-wavefunction calculation failed ({status}). See {}",
            program.label(),
            dir.display()
        ));
    }
    let log = fs::read_to_string(dir.join("program.stdout")).map_err(|e| e.to_string())?;
    crate::qchem_interfaces::gradient::read(program, &dir, atoms.len(), &log)?;
    produced(program, &dir)
}
