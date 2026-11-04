//! MRCC execution uses the selected dmrcc directly in its own MINP directory.
use super::*;
use crate::qchem_interfaces::mrcc;
#[path = "../vibronic/process.rs"]
#[allow(dead_code)]
mod process;
#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    method: QcMethod,
    job: JobType,
    binary: &Path,
    dir: &Path,
    atoms: &[String],
    xyz: &[f64],
    charge: i32,
    mult: i32,
    config: &MethodConfig,
    excited: ExcitedStateOptions,
    orbitals: bool,
    custom: Option<&str>,
    cancel: &AtomicBool,
) -> Result<QcRunOutput, String> {
    let generated = mrcc::input(
        job, method, charge, mult, config, excited, orbitals, atoms, xyz,
    )?;
    if let Some(custom) = custom {
        mrcc::validate_custom(&generated, custom)?;
    }
    fs::create_dir_all(dir).map_err(|e| e.to_string())?;
    fs::write(dir.join("MINP"), custom.unwrap_or(&generated)).map_err(|e| e.to_string())?;
    let mut log = process::run(
        mrcc::command(binary, dir, config.nproc)?,
        dir,
        "mrcc",
        cancel,
    )?;
    let energy = mrcc::parse_energy(&log, method).map_err(|e| {
        format!(
            "{e}\nOutput: {}\n{}",
            dir.join("mrcc.out").display(),
            log.lines()
                .rev()
                .take(30)
                .collect::<Vec<_>>()
                .into_iter()
                .rev()
                .collect::<Vec<_>>()
                .join("\n")
        )
    })?;
    let mut saveable = vec![
        ("MRCC input", dir.join("MINP")),
        ("MRCC output", dir.join("mrcc.out")),
    ];
    if job == JobType::Gradient {
        let g = mrcc::parse_gradient(&log, atoms.len())?;
        let mut text = String::from("# MRCC Cartesian gradient, Hartree/bohr, input orientation\n");
        for (a, r) in atoms.iter().zip(g.chunks_exact(3)) {
            text.push_str(&format!("{a} {:.14e} {:.14e} {:.14e}\n", r[0], r[1], r[2]));
        }
        fs::write(dir.join("mrcc_gradient.txt"), text).map_err(|e| e.to_string())?;
        saveable.push(("Cartesian gradient", dir.join("mrcc_gradient.txt")));
    }
    if job == JobType::TdDft {
        let parsed = mrcc::parse_spectrum(&log, &dir.join("mrcc.out"))?;
        if parsed.states.len() != excited.states as usize {
            return Err(format!(
                "MRCC returned {} excited states; {} were requested.",
                parsed.states.len(),
                excited.states
            ));
        }
    }
    let molden = orbitals.then(|| dir.join("MOLDEN")).filter(|p| p.is_file());
    if let Some(path) = &molden {
        saveable.push(("Molden wavefunction", path.clone()));
    }
    if orbitals && molden.is_none() {
        log.push_str("\nMRCC did not produce the requested MOLDEN wavefunction.\n");
    }
    // Machine-readable energy is also used by the optimizer gradient handoff.
    log.push_str(&format!("\nBeavyr MRCC selected energy: {energy:.16}\n"));
    Ok(QcRunOutput {
        program: QcProgram::Mrcc,
        job,
        log,
        trajectory_text: None,
        hessian_path: None,
        molden_path: molden,
        kept_dir: Some(dir.to_owned()),
        saveable,
        energetics: Energetics {
            final_energy_hartree: Some(energy),
            ..Default::default()
        },
        scf: vec![],
        converged: None,
        scan: None,
    })
}

#[cfg(test)]
#[path = "mrcc_job_tests.rs"]
mod tests;
