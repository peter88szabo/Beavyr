//! Cancellable one-dimensional scans, sharing the calculator's selected backends.
#[cfg(test)]
use super::calculation_objective::pyscf_script;
use super::calculation_objective::{configure_method, Evaluator};
use super::*;
use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
use crate::optimizer::constrained::optimize_constrained_with_internals;
use crate::optimizer::traits::Objective;
use crate::qchem_interfaces::job::scan::{self, Driver, Mode};
use crate::qchem_panel::scan_results::{ScanPoint, ScanReport};
#[path = "scan_behemoth.rs"]
mod behemoth_native;
#[path = "scan_geometry.rs"]
pub(crate) mod geometry;
#[path = "scan_orca.rs"]
mod native;

pub(crate) fn progress(task: &QcRunTask) -> Option<String> {
    if !task.is_running() {
        return None;
    }
    task.run_dir
        .as_ref()
        .and_then(|dir| fs::read_to_string(dir.join("scan_progress.txt")).ok())
}
fn check(cancel: &AtomicBool) -> Result<(), String> {
    if cancel.load(Ordering::SeqCst) {
        Err("Cancelled.".into())
    } else {
        Ok(())
    }
}
fn persist(report: &ScanReport, dir: &Path) -> Result<(), String> {
    fs::write(dir.join("scan.xyz"), report.xyz()).map_err(|e| e.to_string())?;
    fs::write(dir.join("scan.csv"), report.csv()).map_err(|e| e.to_string())?;
    fs::write(
        dir.join("scan_results.ron"),
        ron::ser::to_string_pretty(report, ron::ser::PrettyConfig::default())
            .map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())
}
#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    program: QcProgram,
    method: Option<QcMethod>,
    binary: &Path,
    dir: &Path,
    original_atoms: &[String],
    original_coordinates: &[f64],
    charge: i32,
    multiplicity: i32,
    original_config: &MethodConfig,
    cancel: &AtomicBool,
    child: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    check(cancel)?;
    let mut config = original_config.clone();
    config.scan = config.scan.remap_dummies(original_atoms)?;
    let (atoms, coordinates) =
        crate::calculation_geometry::real_coordinates(original_atoms, original_coordinates)?;
    let molecule = crate::molecule::Molecule::from_xyz(&format!(
        "{}\nscan\n{}",
        atoms.len(),
        orca_run::xyz_body(&atoms, &coordinates)
    ));
    crate::qchem_interfaces::valence::validate_electronic_state(&molecule, charge, multiplicity)
        .map_err(|e| e.to_string())?;
    if let Some(error) = scan::driver_problem(program, method, &config.scan) {
        return Err(error.into());
    }
    if let Some(error) = scan::method_problem(program, method, config.scan.mode) {
        return Err(error.into());
    }
    if !crate::qchem_interfaces::job::methods_for(program).is_empty()
        && !method.is_some_and(|m| crate::qchem_interfaces::job::methods_for(program).contains(&m))
    {
        return Err("Choose a supported scan method.".into());
    }
    if let Some(selected) = method {
        crate::qchem_interfaces::job::correlated::validate_active_space(
            selected,
            &config.active_space,
            &atoms,
            charge,
            multiplicity,
        )?;
    }
    if let Some(module) = program.python_module() {
        let env = crate::qchem_interfaces::python_env::cached(module);
        if !env.is_ready() {
            return Err(env.detail());
        }
    }
    if program.needs_binary_path() && binary.as_os_str().is_empty() {
        return Err(format!("Set the {} executable path.", program.label()));
    }
    // Resolve a relative file path before changing the child's working directory.
    // Keep the exact selected executable; do not search for a substitute.
    let absolute_binary = if !binary.is_absolute() && binary.components().count() > 1 {
        Some(
            std::env::current_dir()
                .map_err(|e| e.to_string())?
                .join(binary),
        )
    } else {
        None
    };
    let binary = absolute_binary.as_deref().unwrap_or(binary);
    fs::create_dir_all(dir).map_err(|e| e.to_string())?;
    configure_method(&mut config, method);
    let behemoth_native = program == QcProgram::Behemoth && config.scan.driver == Driver::Behemoth;
    let native = program == QcProgram::Orca && config.scan.driver == Driver::Orca;
    let mut report = ScanReport {
        program: program.label().into(),
        method: method
            .map(|m| {
                if m == QcMethod::Dft {
                    config.functional.clone()
                } else {
                    m.label().into()
                }
            })
            .unwrap_or_else(|| "DREIDING".into()),
        driver: if native {
            "ORCA built-in"
        } else if behemoth_native {
            "Behemoth built-in"
        } else {
            "Beavyr"
        }
        .into(),
        options: config.scan.clone(),
        atoms: atoms.clone(),
        points: vec![],
        complete: false,
        message: String::new(),
    };
    persist(&report, dir)?;
    let result = if native {
        native::run(
            method.unwrap(),
            binary,
            dir,
            &coordinates,
            charge,
            multiplicity,
            &config,
            &mut report,
            cancel,
            child,
        )
    } else if behemoth_native {
        behemoth_native::run(
            binary,
            dir,
            &coordinates,
            charge,
            multiplicity,
            &config,
            &mut report,
            cancel,
            child,
        )
    } else if program == QcProgram::Dreiding {
        use crate::forcefield::dreiding::{
            objective::{cancellation::CancellableObjective, DreidingObjective},
            DreidingTopology,
        };
        let xyz = format!(
            "{}\nscan\n{}",
            atoms.len(),
            orca_run::xyz_body(&atoms, &coordinates)
        );
        let mol = crate::qchem_interfaces::dreiding_run::molecule_for_force_field(&xyz);
        let topology = DreidingTopology::build(&mol)
            .map_err(|e| crate::forcefield::dreiding::diagnostics::message(&e))?;
        let mut objective =
            CancellableObjective::new(DreidingObjective::new(&topology), Some(cancel));
        compute(&mut objective, &coordinates, &mut report, dir, cancel)
    } else {
        let mut objective = Evaluator {
            program,
            method,
            binary,
            dir,
            atoms: &atoms,
            charge,
            multiplicity,
            config: &config,
            cancel,
            child,
            evaluations: 0,
            cached: None,
        };
        compute(&mut objective, &coordinates, &mut report, dir, cancel)
    };
    check(cancel)?;
    match result {
        Ok(()) => {
            report.complete = report.points.len() == report.options.points
                && report.points.iter().all(|p| p.converged != Some(false));
            report.message = if report.complete {
                format!("{} scan points completed.", report.points.len())
            } else {
                format!(
                    "{} of {} points calculated; at least one point did not converge.",
                    report.points.len(),
                    report.options.points
                )
            };
        }
        Err(error) => {
            report.message = format!(
                "Scan stopped after {} of {} points: {error}",
                report.points.len(),
                report.options.points
            )
        }
    }
    persist(&report, dir)?;
    let mut saveable = vec![
        ("Scan trajectory", dir.join("scan.xyz")),
        ("Scan table", dir.join("scan.csv")),
        ("Scan results", dir.join("scan_results.ron")),
    ];
    if dir.join("orca_scan.inp").is_file() {
        saveable.push(("ORCA scan input", dir.join("orca_scan.inp")));
        saveable.push(("ORCA scan output", dir.join("orca_scan.out")));
    }
    if dir.join("behemoth_scan.out").is_file() {
        saveable.push(("Behemoth scan output", dir.join("behemoth_scan.out")));
        saveable.push((
            "Behemoth scan command",
            dir.join("behemoth_scan_command.txt"),
        ));
    }
    Ok(QcRunOutput {
        job: JobType::Scan,
        program,
        log: report.message.clone(),
        trajectory_text: None,
        hessian_path: None,
        molden_path: None,
        saveable,
        kept_dir: Some(dir.into()),
        energetics: Energetics::default(),
        scf: vec![],
        converged: None,
        scan: Some(Arc::new(report)),
    })
}
fn compute<O: Objective>(
    objective: &mut O,
    start: &[f64],
    report: &mut ScanReport,
    dir: &Path,
    cancel: &AtomicBool,
) -> Result<(), String> {
    let options = report.options.clone();
    let prepared = options.prepare(&report.atoms)?;
    let internals = if options.mode == Mode::Relaxed {
        Some(geometry::internals(
            &report.atoms,
            start,
            &options,
            &prepared,
        )?)
    } else {
        None
    };
    let mut previous: Vec<f64> = start.iter().map(|v| v / BOHR_TO_ANGSTROM).collect();
    for (i, target) in options.values().into_iter().enumerate() {
        check(cancel)?;
        fs::write(
            dir.join("scan_progress.txt"),
            format!(
                "Scan point {} / {} · target {:.5} {}",
                i + 1,
                options.points,
                target,
                options.coordinate.unit()
            ),
        )
        .map_err(|e| e.to_string())?;
        let (coords, energy, converged, steps) = if options.mode == Mode::Rigid {
            let coords = geometry::rigid(&report.atoms, start, &options, &prepared, target)?;
            let bohr: Vec<f64> = coords.iter().map(|v| v / BOHR_TO_ANGSTROM).collect();
            let energy = objective
                .energy(&bohr)
                .map_err(|e| format!("Point {}: {e:#}", i + 1))?;
            (coords, energy, None, 0)
        } else {
            let mut constraints = vec![options.coordinate.target(&prepared.atoms, target)];
            constraints.extend(prepared.fixed.clone());
            let result = optimize_constrained_with_internals(
                objective,
                previous.clone(),
                internals.as_ref().unwrap().clone(),
                constraints,
                crate::optimizer::constrained::ConstrainedOptimizationOptions {
                    max_cycles: options.max_cycles,
                    coordinate_mode: options.relaxation_coordinates.optimizer_mode(),
                    verbosity: 0,
                    ..Default::default()
                },
            )
            .map_err(|e| format!("Point {}: {e:#}", i + 1))?;
            previous = result.x.clone();
            let coords: Vec<f64> = result.x.iter().map(|v| v * BOHR_TO_ANGSTROM).collect();
            let mut valid = result.converged
                && geometry::residual(
                    options.coordinate,
                    geometry::measure(options.coordinate, &prepared.atoms, &coords)?,
                    target,
                )
                .abs()
                    < 1e-3;
            for fixed in &options.constraints {
                let a = scan::atom_indices(&fixed.atoms, fixed.coordinate.count(), &report.atoms)?;
                valid &= geometry::residual(
                    fixed.coordinate,
                    geometry::measure(fixed.coordinate, &a, &coords)?,
                    fixed.value,
                )
                .abs()
                    < 1e-3;
            }
            (coords, result.energy, Some(valid), result.cycles)
        };
        check(cancel)?;
        if !energy.is_finite() || coords.iter().any(|v| !v.is_finite()) {
            return Err(format!("Non-finite result at point {}.", i + 1));
        }
        let actual = geometry::measure(options.coordinate, &prepared.atoms, &coords)?;
        report.points.push(ScanPoint {
            target,
            actual,
            energy,
            converged,
            steps,
            coordinates: coords,
        });
        persist(report, dir)?;
        if converged == Some(false) && options.stop_on_failure {
            break;
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "scan_tests.rs"]
mod tests;
