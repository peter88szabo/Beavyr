//! ORCA's relaxed surface scan and native multiple-XYZ energy scan routes.
use super::*;

#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    method: QcMethod,
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
    let prepared = report.options.prepare(&report.atoms)?;
    let options = report.options.clone();
    let values = options.values();
    let mut rigid_frames = vec![];
    let job = if options.mode == Mode::Relaxed {
        JobType::Optimize
    } else {
        JobType::SinglePoint
    };
    let mut config = config.clone();
    config.optimization.set(QcProgram::Orca,
        (options.mode == Mode::Relaxed && options.tolerances != crate::qchem_interfaces::method::tolerances::Criteria::defaults(QcProgram::Orca))
            .then(|| options.tolerances.clone()));
    let base = orca_run::panel_input_text(
        job,
        method,
        charge,
        multiplicity,
        &config,
        Default::default(),
        &config.orca_extra,
        &orca_run::xyz_body(&report.atoms, start),
    );
    let coordinate_start = base
        .find("* xyz ")
        .ok_or("ORCA scan input has no coordinate block")?;
    let mut input = base[..coordinate_start].to_owned();
    if options.mode == Mode::Relaxed {
        input.push_str(&format!(
            "%geom\n MaxIter {}\n Scan\n  {} {} = {:.10}, {:.10}, {}\n end\n",
            options.max_cycles,
            options.coordinate.keyword(),
            prepared
                .atoms
                .iter()
                .map(|i| i.to_string())
                .collect::<Vec<_>>()
                .join(" "),
            options.start,
            options.end,
            options.points
        ));
        if !options.constraints.is_empty() {
            input.push_str(" Constraints\n");
            for c in &options.constraints {
                let a = scan::atom_indices(&c.atoms, c.coordinate.count(), &report.atoms)?;
                input.push_str(&format!(
                    "  {{ {} {} {:.10} C }}\n",
                    c.coordinate.keyword(),
                    a.iter()
                        .map(|i| i.to_string())
                        .collect::<Vec<_>>()
                        .join(" "),
                    c.value
                ));
            }
            input.push_str(" end\n");
        }
        input.push_str("end\n");
        input.push_str(&base[coordinate_start..]);
    } else {
        let mut xyz = String::new();
        for (i, &target) in values.iter().enumerate() {
            check(cancel)?;
            let coords = geometry::rigid(&report.atoms, start, &options, &prepared, target)?;
            xyz.push_str(&format!(
                "{}\nscan point {} target={target}\n{}",
                report.atoms.len(),
                i + 1,
                orca_run::xyz_body(&report.atoms, &coords)
            ));
            rigid_frames.push(coords);
        }
        fs::write(dir.join("scan_geometries.xyz"), xyz).map_err(|e| e.to_string())?;
        input.push_str(&format!(
            "* xyzfile {charge} {multiplicity} scan_geometries.xyz\n"
        ));
    }
    fs::write(dir.join("orca_scan.inp"), input).map_err(|e| e.to_string())?;
    fs::write(
        dir.join("scan_progress.txt"),
        format!("ORCA built-in scan running · {} points", options.points),
    )
    .map_err(|e| e.to_string())?;
    check(cancel)?;
    let mut command = std::process::Command::new(binary);
    command
        .arg("orca_scan.inp")
        .current_dir(dir)
        .stdout(Stdio::from(
            fs::File::create(dir.join("orca_scan.out")).map_err(|e| e.to_string())?,
        ))
        .stderr(Stdio::from(
            fs::File::create(dir.join("orca_scan.err")).map_err(|e| e.to_string())?,
        ));
    let status = spawn_and_wait(command, "ORCA scan", cancel, child)?;
    check(cancel)?;
    let log = fs::read_to_string(dir.join("orca_scan.out")).unwrap_or_default();
    let segments: Vec<_> = log.split("RELAXED SURFACE SCAN STEP").skip(1).collect();
    let rigid_energies: Vec<_> = log
        .lines()
        .filter(|l| l.contains("FINAL SINGLE POINT ENERGY"))
        .filter_map(|l| l.split_whitespace().last()?.parse::<f64>().ok())
        .collect();
    // The native table reports the actual method energy, which can differ from SCF.
    let table_name = if options.mode == Mode::Relaxed {
        "orca_scan.relaxscanact.dat"
    } else {
        "orca_scan.xyzact.dat"
    };
    let table: Vec<Vec<f64>> = fs::read_to_string(dir.join(table_name))
        .unwrap_or_default()
        .lines()
        .filter_map(|l| {
            l.split_whitespace()
                .map(str::parse::<f64>)
                .collect::<Result<Vec<_>, _>>()
                .ok()
        })
        .filter(|row| row.len() == 2 && row.iter().all(|x| x.is_finite()))
        .collect();
    for (i, &target) in values.iter().enumerate() {
        let (coordinates, energy, converged, steps) = if options.mode == Mode::Relaxed {
            let Some(segment) = segments.get(i) else {
                break;
            };
            let Some(energy) = orca_run::parse_final_energy(segment).filter(|e| e.is_finite())
            else {
                break;
            };
            let text = match fs::read_to_string(dir.join(format!("orca_scan.{:03}.xyz", i + 1))) {
                Ok(t) => t,
                Err(_) => break,
            };
            let frames = crate::trajectory::parse_multi_xyz(&text)
                .map_err(|e| format!("Cannot read scan point {}: {e}", i + 1))?;
            let frame = frames.first().ok_or("Empty ORCA scan geometry")?;
            if frame.atoms != report.atoms {
                return Err("ORCA changed the scan atom order.".into());
            }
            let coords: Vec<f64> = frame
                .pos
                .iter()
                .flat_map(|p| p.to_array().map(f64::from))
                .collect();
            let valid = orca_run::converged(segment) && !segment.contains("SCF NOT CONVERGED");
            (
                coords,
                energy,
                Some(valid),
                segment.matches("GEOMETRY OPTIMIZATION CYCLE").count(),
            )
        } else {
            let Some(row) = table.get(i) else {
                break;
            };
            let energy = if rigid_energies.len() == table.len() {
                rigid_energies[i]
            } else {
                row[1]
            };
            (rigid_frames[i].clone(), energy, None, 0)
        };
        let actual = geometry::measure(options.coordinate, &prepared.atoms, &coordinates)?;
        let mut valid = geometry::residual(options.coordinate, actual, target).abs() < 1e-3;
        for c in &options.constraints {
            let atoms = scan::atom_indices(&c.atoms, c.coordinate.count(), &report.atoms)?;
            valid &= geometry::residual(
                c.coordinate,
                geometry::measure(c.coordinate, &atoms, &coordinates)?,
                c.value,
            )
            .abs()
                < 1e-3;
        }
        if !energy.is_finite() {
            return Err(format!("Non-finite energy at ORCA scan point {}", i + 1));
        }
        report.points.push(ScanPoint {
            target,
            actual,
            energy,
            converged: converged.map(|c| c && valid),
            steps,
            coordinates,
        });
        persist(report, dir)?;
    }
    if !status.success() {
        return Err(format!(
            "ORCA scan exited with {status}. See {}",
            dir.join("orca_scan.out").display()
        ));
    }
    if log.contains("SCF NOT CONVERGED") {
        return Err("ORCA reported an unconverged SCF calculation during the scan.".into());
    }
    orca_run::check_output(&log, status.success()).map_err(|e| {
        format!(
            "{e}\n{}",
            fs::read_to_string(dir.join("orca_scan.err")).unwrap_or_default()
        )
    })?;
    if report.points.len() != options.points {
        return Err(format!(
            "ORCA returned {} readable points; {} were requested. See {}",
            report.points.len(),
            options.points,
            dir.join("orca_scan.out").display()
        ));
    }
    Ok(())
}
