//! ORCA-native surface crossing, with both spin energies retained.
use super::*;
fn input(
    method: QcMethod,
    atoms: &[String],
    coordinates: &[f64],
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
) -> String {
    // A crossing search has its own tolerances, independent of minimum jobs.
    let mut config = config.clone();
    config.optimization.set(crate::qchem_interfaces::program::QcProgram::Orca,
        config.mecp.native_tolerances.clone());
    let base = orca_run::panel_input_text(
        JobType::Optimize,
        method,
        charge,
        multiplicity,
        &config,
        Default::default(),
        &config.orca_extra,
        &orca_run::xyz_body(atoms, coordinates),
    );
    let split = base
        .find("* xyz ")
        .expect("generated ORCA input has coordinates");
    format!(
        "{}%mecp\n SurfCrossOpt true\n Mult {}\nend\n%geom\n MaxIter {}\nend\n{}",
        &base[..split],
        config.mecp.multiplicity_b,
        config.mecp.max_cycles,
        &base[split..]
    )
}
fn energies(log: &str) -> Result<Vec<(f64, f64)>, String> {
    let mut pending = vec![];
    let mut pairs = vec![];
    for line in log.lines() {
        if let Some(value) = line.strip_prefix("FINAL SINGLE POINT ENERGY") {
            let e = value
                .trim()
                .parse::<f64>()
                .map_err(|_| "Invalid ORCA MECP energy")?;
            if !e.is_finite() {
                return Err("Non-finite ORCA MECP energy.".into());
            }
            pending.push(e);
        }
        if line
            .trim_start()
            .starts_with("Energy difference between both states")
        {
            let gap = line
                .split_whitespace()
                .last()
                .and_then(|s| s.parse::<f64>().ok())
                .ok_or("Missing ORCA MECP gap")?;
            if pending.len() != 2 || !gap.is_finite() {
                return Err("ORCA did not report two unambiguous MECP state energies.".into());
            }
            let (a, b) = (pending[0], pending[1]);
            if (a - b - gap).abs() > 1e-6 {
                return Err("ORCA MECP state energies do not match its reported gap.".into());
            }
            pairs.push((a, b));
            pending.clear();
        }
    }
    if pairs.is_empty() {
        return Err("ORCA returned no MECP state energies.".into());
    }
    Ok(pairs)
}
#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    method: QcMethod,
    binary: &Path,
    dir: &Path,
    atoms: &[String],
    coordinates: &[f64],
    charge: i32,
    multiplicity_a: i32,
    config: &MethodConfig,
    cancel: &AtomicBool,
    child: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    fs::write(
        dir.join("mecp.inp"),
        input(method, atoms, coordinates, charge, multiplicity_a, config),
    )
    .map_err(|e| e.to_string())?;
    let mut command = std::process::Command::new(binary);
    command
        .arg("mecp.inp")
        .current_dir(dir)
        .stdout(Stdio::from(
            fs::File::create(dir.join("mecp.out")).map_err(|e| e.to_string())?,
        ))
        .stderr(Stdio::from(
            fs::File::create(dir.join("mecp.err")).map_err(|e| e.to_string())?,
        ));
    let status = spawn_and_wait(command, "ORCA MECP", cancel, child)?;
    if cancel.load(Ordering::SeqCst) {
        return Err("Cancelled.".into());
    }
    let raw = fs::read_to_string(dir.join("mecp.out")).map_err(|e| e.to_string())?;
    orca_run::check_output(&raw, status.success())?;
    let pairs = energies(&raw)?;
    let frames = crate::trajectory::parse_multi_xyz(
        &fs::read_to_string(dir.join("mecp_trj.xyz")).map_err(|e| e.to_string())?,
    )?;
    if frames.len() != pairs.len() {
        return Err(format!(
            "ORCA returned {} trajectory frames but {} MECP energy pairs; raw files retained.",
            frames.len(),
            pairs.len()
        ));
    }
    let mut xyz = String::new();
    let mut csv = String::from("frame,state_a_energy_Eh,state_b_energy_Eh,gap_Eh,mean_energy_Eh\n");
    for (i, (frame, &(a, b))) in frames.iter().zip(&pairs).enumerate() {
        if frame.atoms != atoms {
            return Err("ORCA changed the MECP atom order.".into());
        }
        let mean = (a + b) * 0.5;
        let coordinates: Vec<_> = frame
            .pos
            .iter()
            .flat_map(|p| p.to_array().map(f64::from))
            .collect();
        xyz.push_str(&format!(
            "{}\ncycle {} E = {mean:.12} Eh Ea={a:.12} Eb={b:.12} gap={:.12}\n{}",
            atoms.len(),
            i + 1,
            a - b,
            orca_run::xyz_body(atoms, &coordinates)
        ));
        csv.push_str(&format!(
            "{},{a:.12},{b:.12},{:.12},{mean:.12}\n",
            i + 1,
            a - b
        ));
    }
    let &(a, b) = pairs.last().unwrap();
    let converged = raw.contains("THE OPTIMIZATION HAS CONVERGED")
        && !raw.contains("The optimization did not converge")
        && (a - b).abs() <= config.mecp.gap_tolerance;
    let log = format!("ORCA built-in MECP: {}\nState A multiplicity: {}\nState B multiplicity: {}\nState A energy: {a:.12} Eh\nState B energy: {b:.12} Eh\nEnergy gap (A − B): {:.10} Eh\nGap tolerance: {:.10} Eh\nEvaluated geometries: {}",
        if converged { "Converged" } else { "Not converged to the requested crossing criteria" }, multiplicity_a, config.mecp.multiplicity_b, a-b, config.mecp.gap_tolerance, pairs.len());
    fs::write(dir.join("mecp_summary.txt"), &log).map_err(|e| e.to_string())?;
    fs::write(dir.join("mecp.csv"), csv).map_err(|e| e.to_string())?;
    fs::write(dir.join("mecp_trajectory.xyz"), &xyz).map_err(|e| e.to_string())?;
    Ok(QcRunOutput {
        job: JobType::Mecp,
        program: QcProgram::Orca,
        log,
        trajectory_text: Some(xyz),
        hessian_path: None,
        molden_path: None,
        saveable: vec![
            ("ORCA MECP input", dir.join("mecp.inp")),
            ("ORCA MECP output", dir.join("mecp.out")),
            ("MECP energies", dir.join("mecp.csv")),
            ("MECP trajectory", dir.join("mecp_trajectory.xyz")),
            ("MECP summary", dir.join("mecp_summary.txt")),
        ],
        kept_dir: Some(dir.to_owned()),
        energetics: Default::default(),
        scf: vec![],
        converged: Some(converged),
        scan: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mecp_orca_input_sets_both_states_and_parser_checks_gap() {
        let mut config = MethodConfig::default();
        config.mecp.multiplicity_b = 5;
        config.mecp.max_cycles = 27;
        let text = input(QcMethod::Hf, &["Fe".into()], &[0., 0., 0.], 0, 3, &config);
        assert!(
            text.contains("* xyz 0 3")
                && text.contains("Mult 5")
                && text.contains("MaxIter 27")
                && text.contains("SurfCrossOpt true")
        );
        let good = "FINAL SINGLE POINT ENERGY -1.1\nFINAL SINGLE POINT ENERGY -1.2\nEnergy difference between both states 0.1\n";
        assert_eq!(energies(good).unwrap(), vec![(-1.1, -1.2)]);
        assert!(energies(&good.replace("states 0.1", "states 0.9")).is_err());
    }
}
