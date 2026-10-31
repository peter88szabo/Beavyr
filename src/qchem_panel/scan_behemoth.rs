//! Behemoth's native scan CLI; read its trajectory and per-point convergence markers.
use super::*;
use crate::qchem_interfaces::method::{
    behemoth_basis_flag, behemoth_dispersion_flag, behemoth_functional_flag, behemoth_method_flag,
};

#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    binary: &Path,
    dir: &Path,
    start: &[f64],
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    report: &mut ScanReport,
    cancel: &AtomicBool,
    child: &Mutex<Option<Child>>,
) -> Result<(), String> {
    let options = &report.options;
    let prepared = options.prepare(&report.atoms)?;
    fs::write(
        dir.join("geometry.xyz"),
        format!(
            "{}\nBeavyr scan\n{}",
            report.atoms.len(),
            orca_run::xyz_body(&report.atoms, start)
        ),
    )
    .map_err(|e| e.to_string())?;
    let mut command = std::process::Command::new(binary);
    command
        .current_dir(dir)
        .arg("--xyz")
        .arg("geometry.xyz")
        .arg("--method")
        .arg(behemoth_method_flag(config, multiplicity))
        .arg("--charge")
        .arg(charge.to_string())
        .arg("--multi")
        .arg(multiplicity.to_string())
        .arg("--memory")
        .arg(config.memory_mb.to_string())
        .arg("--nproc")
        .arg(config.nproc.max(1).to_string());
    if let Some(functional) = behemoth_functional_flag(config) {
        command.arg("--functional").arg(functional);
    }
    if let Some(basis) = behemoth_basis_flag(config) {
        command.arg("--basis").arg(basis);
    }
    if let Some(dispersion) = behemoth_dispersion_flag(config) {
        command.arg("--dispersion").arg(dispersion);
    }
    command
        .arg("--scan")
        .arg("--scan-coordinate")
        .arg(options.coordinate.label().to_lowercase())
        .arg("--index-base")
        .arg("1")
        .arg("--scan-atoms")
        .arg(
            prepared
                .atoms
                .iter()
                .map(|a| (a + 1).to_string())
                .collect::<Vec<_>>()
                .join(","),
        )
        .arg("--scan-start")
        .arg(options.start.to_string())
        .arg("--scan-end")
        .arg(options.end.to_string())
        .arg("--scan-points")
        .arg(options.points.to_string())
        .arg("--scan-mode")
        .arg(if options.mode == Mode::Rigid {
            "rigid"
        } else {
            "relaxed"
        })
        .arg("--opt-coord")
        .arg(options.relaxation_coordinates.keyword());
    command
        .arg("--scan-max-cycles")
        .arg(options.max_cycles.to_string())
        .arg("--scan-move-fragment")
        .arg(options.move_fragment.to_string())
        .arg("--scan-bond-motion")
        .arg(options.bond_motion.keyword())
        .arg("--scan-stop-on-failure")
        .arg(options.stop_on_failure.to_string());
    if !options.constraints.is_empty() {
        let fixed = options
            .constraints
            .iter()
            .map(|c| {
                let atoms = scan::atom_indices(&c.atoms, c.coordinate.count(), &report.atoms)?;
                Ok(format!(
                    "{}:{}:{}",
                    c.coordinate.label().to_lowercase(),
                    atoms
                        .iter()
                        .map(|a| (a + 1).to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                    c.value
                ))
            })
            .collect::<Result<Vec<String>, String>>()?;
        command.arg("--scan-constraints").arg(fixed.join(";"));
    }
    fs::write(
        dir.join("behemoth_scan_command.txt"),
        format!("{command:?}\n"),
    )
    .map_err(|e| e.to_string())?;
    command
        .stdout(Stdio::from(
            fs::File::create(dir.join("behemoth_scan.out")).map_err(|e| e.to_string())?,
        ))
        .stderr(Stdio::from(
            fs::File::create(dir.join("behemoth_scan.err")).map_err(|e| e.to_string())?,
        ));
    fs::write(
        dir.join("scan_progress.txt"),
        format!("Behemoth built-in scan running · {} points", options.points),
    )
    .map_err(|e| e.to_string())?;
    check(cancel)?;
    let status = spawn_and_wait(command, "Behemoth scan", cancel, child)?;
    check(cancel)?;
    let log = fs::read_to_string(dir.join("behemoth_scan.out")).unwrap_or_default();
    let errors = fs::read_to_string(dir.join("behemoth_scan.err")).unwrap_or_default();
    let file = if options.mode == Mode::Rigid {
        "rigid_scan.xyz"
    } else {
        "relaxed_scan.xyz"
    };
    if let Ok(text) = fs::read_to_string(dir.join(file)) {
        read_results(&text, report)?;
        persist(report, dir)?;
    }
    if !status.success() {
        return Err(format!(
            "Behemoth scan failed ({status}).\n{errors}\nOutput: {}",
            dir.join("behemoth_scan.out").display()
        ));
    }
    if log.contains("SCF did not converge") {
        return Err("Behemoth reported an unconverged SCF calculation.".into());
    }
    if report.points.len() != report.options.points {
        return Err(format!(
            "Behemoth returned {} readable points; {} were requested. See {}",
            report.points.len(),
            report.options.points,
            dir.join("behemoth_scan.out").display()
        ));
    }
    Ok(())
}
fn field<'a>(text: &'a str, key: &str) -> Option<&'a str> {
    text.split_whitespace().find_map(|s| s.strip_prefix(key))
}
fn read_results(text: &str, report: &mut ScanReport) -> Result<(), String> {
    let frames = crate::trajectory::parse_multi_xyz(text)?;
    let lines: Vec<_> = text.lines().collect();
    let prepared = report.options.prepare(&report.atoms)?;
    if frames.len() > report.options.points {
        return Err("Behemoth returned more scan frames than requested.".into());
    }
    for (i, (frame, target)) in frames.into_iter().zip(report.options.values()).enumerate() {
        if frame.atoms != report.atoms {
            return Err("Behemoth changed the scan atom order.".into());
        }
        let comment = lines
            .get(i * (report.atoms.len() + 2) + 1)
            .ok_or("Incomplete Behemoth trajectory")?;
        let energy = field(comment, "E=")
            .or_else(|| {
                comment
                    .split("E =")
                    .nth(1)
                    .and_then(|s| s.split_whitespace().next())
            })
            .and_then(|s| s.parse::<f64>().ok())
            .filter(|e| e.is_finite())
            .ok_or("Behemoth scan point has no finite energy")?;
        let coordinates: Vec<_> = frame
            .pos
            .iter()
            .flat_map(|p| p.to_array().map(f64::from))
            .collect();
        let actual = geometry::measure(report.options.coordinate, &prepared.atoms, &coordinates)?;
        let mut valid = geometry::residual(report.options.coordinate, actual, target).abs() < 1e-3;
        for fixed in &report.options.constraints {
            let atoms = scan::atom_indices(&fixed.atoms, fixed.coordinate.count(), &report.atoms)?;
            valid &= geometry::residual(
                fixed.coordinate,
                geometry::measure(fixed.coordinate, &atoms, &coordinates)?,
                fixed.value,
            )
            .abs()
                < 1e-3;
        }
        let converged = if report.options.mode == Mode::Relaxed {
            Some(
                field(comment, "converged=")
                    .and_then(|s| s.parse::<bool>().ok())
                    .ok_or("Behemoth scan point has no convergence status")?
                    && valid,
            )
        } else {
            if !valid {
                return Err(
                    "Behemoth rigid scan coordinate differs from the requested target.".into(),
                );
            }
            None
        };
        let steps = if report.options.mode == Mode::Relaxed {
            field(comment, "opt_steps=")
                .and_then(|s| s.parse().ok())
                .ok_or("Behemoth scan point has no optimization step count")?
        } else {
            0
        };
        report.points.push(ScanPoint {
            target,
            actual,
            energy,
            converged,
            steps,
            coordinates,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn scan_behemoth_native_requires_explicit_convergence_and_matching_geometry() {
        let mut r = ScanReport {
            program: "Behemoth".into(),
            method: "HF".into(),
            driver: "Behemoth built-in".into(),
            options: scan::Options {
                start: 1.,
                end: 1.1,
                points: 2,
                ..Default::default()
            },
            atoms: vec!["H".into(), "H".into()],
            points: vec![],
            complete: false,
            message: String::new(),
        };
        let text = "2\nscan_step=1 E=-1.0 converged=false opt_steps=50\nH 0 0 0\nH 1 0 0\n";
        read_results(text, &mut r).unwrap();
        assert_eq!(r.points[0].converged, Some(false));
        r.points.clear();
        assert!(read_results(&text.replace("converged=false ", ""), &mut r).is_err());
        r.points.clear();
        read_results(
            &text
                .replace("converged=false", "converged=true")
                .replace("H 1 0 0", "H 1.4 0 0"),
            &mut r,
        )
        .unwrap();
        assert_eq!(r.points[0].converged, Some(false));
    }
}
