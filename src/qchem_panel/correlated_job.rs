//! Energy-only Behemoth coupled-cluster and multireference calculations.
use super::*;
use crate::qchem_interfaces::job::correlated;

pub(super) fn input_text(
    method: QcMethod,
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    wavefunction: bool,
) -> String {
    let reference = if method == QcMethod::CcsdT {
        // Closed-shell RHF is Behemoth's plain "HF"; open shells take UHF
        // (default) or ROHF as chosen in the panel.
        match correlated::reference(config.cc_reference, multiplicity) {
            correlated::Reference::Rhf => "HF",
            correlated::Reference::Rohf => "ROHF",
            correlated::Reference::Uhf => "UHF",
        }
    } else if multiplicity == 1 {
        "HF"
    } else {
        "ROHF"
    };
    let correlation = if method == QcMethod::CcsdT {
        "CCSD(T)".to_owned()
    } else {
        format!(
            "CASSCF[{}e,{}o]{}",
            config.active_space.electrons,
            config.active_space.orbitals,
            if method == QcMethod::Nevpt2 {
                " NEVPT2"
            } else {
                ""
            }
        )
    };
    format!("nproc = {}\nmemory = {}\nmethod = {reference} {correlation}\nbasis = {}\ncharge = {charge}\nmulti = {multiplicity}\nxyzfile = input.xyz\n{}",
        config.nproc.max(1), config.memory_mb, config.basis,
        if wavefunction { "wf2molden = true\n" } else { "" })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn run_behemoth(
    method: QcMethod,
    binary: &Path,
    workdir: &Path,
    atoms: &[String],
    coordinates: &[f64],
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    wavefunction: bool,
    cancel: &AtomicBool,
    child: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    correlated::validate_active_space(method, &config.active_space, atoms, charge, multiplicity)?;
    if cancel.load(Ordering::SeqCst) {
        return Err("Cancelled.".into());
    }
    if config.basis.contains(['\n', '\r']) {
        return Err("The basis name must be a single line.".into());
    }
    let input = workdir.join("correlated.inp");
    let geometry = workdir.join("input.xyz");
    fs::write(
        &input,
        input_text(method, charge, multiplicity, config, wavefunction),
    )
    .map_err(|e| e.to_string())?;
    fs::write(
        &geometry,
        format!(
            "{}\nBeavyr {} energy\n{}",
            atoms.len(),
            method.label(),
            orca_run::xyz_body(atoms, coordinates)
        ),
    )
    .map_err(|e| e.to_string())?;
    let out = workdir.join("behemoth.out");
    let err = workdir.join("behemoth.err");
    let mut command = std::process::Command::new(binary);
    command
        .arg("correlated.inp")
        .current_dir(workdir)
        .env("OMP_NUM_THREADS", config.nproc.max(1).to_string())
        .env("OPENBLAS_NUM_THREADS", config.nproc.max(1).to_string())
        .env("MKL_NUM_THREADS", config.nproc.max(1).to_string())
        .stdout(Stdio::from(
            fs::File::create(&out).map_err(|e| e.to_string())?,
        ))
        .stderr(Stdio::from(
            fs::File::create(&err).map_err(|e| e.to_string())?,
        ));
    let status = spawn_and_wait(command, "Behemoth", cancel, child)?;
    let log = fs::read_to_string(&out).unwrap_or_default();
    let errors = fs::read_to_string(&err).unwrap_or_default();
    if !status.success() {
        return Err(format!(
            "Behemoth {} failed ({status}).\n{log}\n{errors}",
            method.label()
        ));
    }
    let energy = parse_energy(method, &log).map_err(|error| format!("{error}\n{log}\n{errors}"))?;
    if cancel.load(Ordering::SeqCst) {
        return Err("Cancelled.".into());
    }
    let mut saveable = vec![
        ("Input", input),
        ("Input geometry", geometry),
        ("Behemoth output", out),
    ];
    let molden_path = if wavefunction {
        ["correlated.casscf_canonicalMO.molden", "canonicalMO.molden"]
            .iter()
            .map(|name| workdir.join(name))
            .find(|path| path.is_file())
    } else {
        None
    };
    if let Some(path) = &molden_path {
        saveable.push(("Orbitals (Molden)", path.clone()));
    }
    Ok(QcRunOutput {
        scan: None,
        job: JobType::SinglePoint,
        program: QcProgram::Behemoth,
        log,
        trajectory_text: None,
        hessian_path: None,
        molden_path,
        saveable,
        kept_dir: Some(workdir.to_owned()),
        energetics: Energetics {
            final_energy_hartree: Some(energy),
            ..Default::default()
        },
        scf: vec![],
        converged: None,
    })
}

pub(super) fn parse_energy(method: QcMethod, log: &str) -> Result<f64, String> {
    if log.contains("NOT converged") || log.contains("did not converge") {
        return Err(format!(
            "Behemoth {} did not converge; no energy is accepted.",
            method.label()
        ));
    }
    let mut section = false;
    let mut energy = None;
    for line in log.lines().map(str::trim) {
        let value = match method {
            QcMethod::CcsdT => line
                .strip_prefix("Total CCSD(T) energy:")
                .or_else(|| line.strip_prefix("Total UCCSD(T) energy:")),
            QcMethod::Casscf => {
                if line == "CASSCF summary" {
                    section = true;
                }
                if section {
                    line.strip_prefix("total energy")
                } else {
                    None
                }
            }
            QcMethod::Nevpt2 => {
                if line.starts_with("SC-NEVPT2") && line.contains("root") {
                    section = true;
                }
                if section {
                    line.strip_prefix("E(total)")
                } else {
                    None
                }
            }
            _ => None,
        };
        if let Some(value) = value {
            energy = value
                .split_whitespace()
                .next()
                .and_then(|n| n.parse::<f64>().ok());
        }
    }
    energy.filter(|e| e.is_finite()).ok_or_else(|| {
        format!(
            "Behemoth printed no finite {} total energy.",
            method.label()
        )
    })
}

#[cfg(test)]
#[path = "correlated_job_tests.rs"]
mod tests;
