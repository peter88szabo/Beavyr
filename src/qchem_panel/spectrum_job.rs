//! Spectrum adapters for Behemoth and standalone Grimme xtb4stda + stda.
//! Standalone CLI: https://github.com/grimme-lab/xtb4stda and the sTDA manual,
//! section 7.2. Output columns follow stda.f: state, eV, nm, fL, Rv.
use super::*;
use crate::uvvis::excited_state::{self, ExcitedStateConfig, ExcitedStateMethod};
use crate::uvvis::types::{ExcitedState, TddftResult, HARTREE_TO_EV};

pub fn generator_config_path() -> Option<PathBuf> {
    crate::qchem_interfaces::config::config_path_for(QcProgram::StdaXtb)
        .map(|path| path.with_file_name("xtb4stda_path.txt"))
}
pub fn generator_path() -> String {
    generator_config_path()
        .and_then(|path| fs::read_to_string(path).ok())
        .unwrap_or_default()
        .trim()
        .to_owned()
}
pub fn save_generator_path(path: &str) -> Result<(), String> {
    let file = generator_config_path().ok_or("No configuration directory is available")?;
    fs::create_dir_all(file.parent().unwrap()).map_err(|e| e.to_string())?;
    fs::write(file, path.trim()).map_err(|e| format!("Cannot save xtb4stda path: {e}"))
}

pub fn behemoth_config(options: ExcitedStateOptions) -> ExcitedStateConfig {
    ExcitedStateConfig {
        method: options.spectrum_method,
        roots: options.states,
        emax_ev: options.energy_window_ev as f64,
        tda: options.tda,
        xtb4stda_path: String::new(),
    }
}

pub fn validate(
    program: QcProgram,
    options: ExcitedStateOptions,
    reference: &MethodConfig,
    multiplicity: i32,
    atoms: &[String],
) -> Result<(), String> {
    if options.energy_window_ev == 0 {
        return Err("The spectrum energy window must be positive.".into());
    }
    if program == QcProgram::StdaXtb {
        if multiplicity != 1 {
            return Err(
                "The standalone sTDA-xTB route requires a closed-shell singlet (multiplicity 1)."
                    .into(),
            );
        }
        return Ok(());
    }
    if options.spectrum_method == ExcitedStateMethod::StdaXtbOriginal {
        return Err("Choose GFN1, TASI or a DFT reference for Behemoth; use sTDA-xTB (Grimme) for the standalone programs.".into());
    }
    let real_atoms: Vec<_> = atoms
        .iter()
        .filter(|a| !crate::molecule::is_dummy(a))
        .cloned()
        .collect();
    let issues = excited_state::validate(
        &behemoth_config(options),
        reference,
        multiplicity,
        &real_atoms,
    );
    if let Some(issue) = issues
        .iter()
        .find(|i| i.severity == crate::qchem_interfaces::method::Severity::Block)
    {
        return Err(issue.message.clone());
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    program: QcProgram,
    job: JobType,
    binary: &Path,
    workdir: &Path,
    atoms: &[String],
    coordinates: &[f64],
    charge: i32,
    multiplicity: i32,
    reference: &MethodConfig,
    options: ExcitedStateOptions,
    cancel: &AtomicBool,
    child: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    run_with_generator(
        program,
        job,
        binary,
        workdir,
        atoms,
        coordinates,
        charge,
        multiplicity,
        reference,
        options,
        cancel,
        child,
        &generator_path(),
    )
}

#[allow(clippy::too_many_arguments)]
fn run_with_generator(
    program: QcProgram,
    job: JobType,
    binary: &Path,
    workdir: &Path,
    atoms: &[String],
    coordinates: &[f64],
    charge: i32,
    multiplicity: i32,
    reference: &MethodConfig,
    options: ExcitedStateOptions,
    cancel: &AtomicBool,
    child: &Mutex<Option<Child>>,
    generator: &str,
) -> Result<QcRunOutput, String> {
    if job != JobType::TdDft {
        return Err("sTDA-xTB calculates spectra only.".into());
    }
    validate(program, options, reference, multiplicity, atoms)?;
    let body = orca_run::xyz_body(atoms, coordinates);
    fs::write(
        workdir.join("input.xyz"),
        format!("{}\nBeavyr spectrum\n{body}", atoms.len()),
    )
    .map_err(|e| e.to_string())?;
    let log = if program == QcProgram::StdaXtb {
        let generator = crate::qchem_interfaces::xtb_optimize::resolve_program_executable(
            QcProgram::StdaXtb,
            generator,
        )
        .map_err(|e| format!("Set the xtb4stda executable path in Configure: {e}"))?;
        // xtb4stda's documented charge and spin files, in this run directory only.
        fs::write(workdir.join(".CHRG"), format!("{charge}\n")).map_err(|e| e.to_string())?;
        fs::write(workdir.join(".UHF"), "0\n").map_err(|e| e.to_string())?;
        let mut ground = std::process::Command::new(generator);
        ground.arg("input.xyz").current_dir(workdir);
        let ground_log = execute(ground, "xtb4stda", workdir, reference.nproc, cancel, child)?;
        if !workdir.join("wfn.xtb").is_file() {
            return Err(format!("xtb4stda did not produce wfn.xtb.\n{ground_log}"));
        }
        let mut spectrum = std::process::Command::new(binary);
        spectrum
            .args(["-xtb", "-e", &options.energy_window_ev.to_string()])
            .current_dir(workdir);
        let spectrum_log = execute(spectrum, "stda", workdir, reference.nproc, cancel, child)?;
        parse_stda(&spectrum_log, &workdir.join("stda.stdout"))?;
        format!("{ground_log}\n{spectrum_log}")
    } else {
        let command = excited_state::spectrum_command(
            binary,
            workdir,
            "input.xyz",
            charge,
            multiplicity,
            &behemoth_config(options),
            reference,
        );
        let log = execute(command, "behemoth", workdir, reference.nproc, cancel, child)?;
        crate::uvvis::behemoth::parse_behemoth_spectrum(
            &log,
            &workdir.join("behemoth.stdout"),
            Some(options.states as usize),
        )?;
        log
    };
    if cancel.load(Ordering::SeqCst) {
        return Err("Cancelled.".into());
    }
    let molden = workdir.join(excited_state::CANONICAL_MOLDEN_FILE);
    let molden_path = molden.is_file().then_some(molden);
    let mut saveable = vec![("Input geometry", workdir.join("input.xyz"))];
    if let Some(path) = &molden_path {
        saveable.push(("Orbitals (Molden)", path.clone()));
    }
    if workdir.join("tda.dat").is_file() {
        saveable.push(("Spectrum data", workdir.join("tda.dat")));
    }
    Ok(QcRunOutput {
        scan: None,
        job,
        program,
        log,
        trajectory_text: None,
        hessian_path: None,
        molden_path,
        saveable,
        kept_dir: Some(workdir.to_owned()),
        energetics: Energetics::default(),
        scf: vec![],
        converged: None,
    })
}

fn execute(
    mut command: std::process::Command,
    name: &str,
    workdir: &Path,
    cores: u32,
    cancel: &AtomicBool,
    child: &Mutex<Option<Child>>,
) -> Result<String, String> {
    if cancel.load(Ordering::SeqCst) {
        return Err("Cancelled.".into());
    }
    let stdout = workdir.join(format!("{name}.stdout"));
    let stderr = workdir.join(format!("{name}.stderr"));
    command
        .env("OMP_NUM_THREADS", cores.max(1).to_string())
        .env("MKL_NUM_THREADS", cores.max(1).to_string())
        .env("OPENBLAS_NUM_THREADS", cores.max(1).to_string())
        .stdout(Stdio::from(
            fs::File::create(&stdout).map_err(|e| e.to_string())?,
        ))
        .stderr(Stdio::from(
            fs::File::create(&stderr).map_err(|e| e.to_string())?,
        ));
    let status = spawn_and_wait(command, name, cancel, child)?;
    let log = fs::read_to_string(stdout).unwrap_or_default();
    let errors = fs::read_to_string(stderr).unwrap_or_default();
    if !status.success() {
        return Err(format!("{name} failed ({status}).\n{log}\n{errors}"));
    }
    Ok(format!("{log}\n{errors}"))
}

pub fn parse_stda(text: &str, source: &Path) -> Result<TddftResult, String> {
    let mut reading = false;
    let mut states = Vec::new();
    for line in text.lines() {
        let columns: Vec<_> = line.split_whitespace().collect();
        if columns.starts_with(&["state", "eV", "nm", "fL"]) {
            reading = true;
            states.clear();
            continue;
        }
        if !reading {
            continue;
        }
        if columns.is_empty() {
            if !states.is_empty() {
                break;
            } else {
                continue;
            }
        }
        let Ok(root) = columns[0].parse::<usize>() else {
            if !states.is_empty() {
                break;
            } else {
                continue;
            }
        };
        if columns.len() < 4 {
            return Err(format!("Incomplete sTDA spectrum row: {line}"));
        }
        let number = |value: &str| {
            value
                .replace(['D', 'd'], "E")
                .parse::<f64>()
                .map_err(|_| format!("Invalid number in sTDA spectrum: {line}"))
        };
        let energy_ev = number(columns[1])?;
        let wavelength_nm = number(columns[2])?;
        let strength = number(columns[3])?;
        if root == 0
            || !energy_ev.is_finite()
            || energy_ev <= 0.0
            || !wavelength_nm.is_finite()
            || wavelength_nm <= 0.0
            || !strength.is_finite()
            || strength < 0.0
        {
            return Err(format!("Invalid sTDA transition: {line}"));
        }
        states.push(ExcitedState {
            root,
            energy_au: energy_ev / HARTREE_TO_EV,
            energy_ev,
            energy_cm1: energy_ev * 8065.54393734921,
            wavelength_nm,
            s2: None,
            multiplicity: Some(1),
            oscillator_strength: Some(strength),
            d2_au2: None,
            transition_dipole_au: None,
            excitations: vec![],
        });
    }
    if states.is_empty() {
        return Err(format!(
            "sTDA produced no readable excitation table (state / eV / nm / fL).\n{text}"
        ));
    }
    Ok(TddftResult {
        source: source.to_owned(),
        program: "sTDA (Grimme)".into(),
        method: "sTDA-xTB".into(),
        unrestricted: false,
        ground_state_s2: None,
        states,
        notes: vec![],
    })
}

#[cfg(test)]
#[path = "spectrum_job_tests.rs"]
mod tests;
