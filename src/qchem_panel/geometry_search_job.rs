//! Calculator connections for existing Behemoth-derived geometry optimizers.
use super::calculation_objective::{configure_method, Evaluator};
use super::*;
use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
use crate::optimizer::internal_coords::ConnectivityModel;
use crate::optimizer::traits::Objective;
use crate::qchem_interfaces::job::search::{self, Coordinates, MinimumAlgorithm};
#[path = "geometry_search_minimum.rs"]
mod minimum;
#[path = "geometry_search_behemoth.rs"]
mod native;
struct Outcome {
    x: Vec<f64>,
    energy: f64,
    converged: bool,
    cycles: usize,
    message: String,
    frames: Vec<Vec<f64>>,
    energies: Vec<f64>,
    hessian: Option<ndarray::Array2<f64>>,
}
#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    program: QcProgram,
    method: Option<QcMethod>,
    job: JobType,
    binary: &Path,
    dir: &Path,
    original_atoms: &[String],
    original_coordinates: &[f64],
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    cancel: &AtomicBool,
    child: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    check(cancel)?;
    let mut config = config.clone();
    if job == JobType::TransitionState {
        config.search = search::remap_dummies(&config.search, original_atoms)?;
    }
    let (atoms, coordinates) =
        crate::calculation_geometry::real_coordinates(original_atoms, original_coordinates)?;
    let mol = crate::molecule::Molecule::from_xyz(&format!(
        "{}\nsearch\n{}",
        atoms.len(),
        orca_run::xyz_body(&atoms, &coordinates)
    ));
    crate::qchem_interfaces::valence::validate_electronic_state(&mol, charge, multiplicity)
        .map_err(|e| e.to_string())?;
    if let Some(e) = search::method_problem(program, method, job, &config.search) {
        return Err(e.into());
    }
    search::validate(&config.search, job, &atoms)?;
    if let Some(module) = program.python_module() {
        let env = crate::qchem_interfaces::python_env::cached(module);
        if !env.is_ready() {
            return Err(env.detail());
        }
    }
    if program.needs_binary_path() && binary.as_os_str().is_empty() {
        return Err(format!("Set the {} executable path.", program.label()));
    }
    let binary = if !binary.is_absolute() && binary.components().count() > 1 {
        std::env::current_dir()
            .map_err(|e| e.to_string())?
            .join(binary)
    } else {
        binary.to_owned()
    };
    fs::create_dir_all(dir).map_err(|e| e.to_string())?;
    configure_method(&mut config, method);
    if program == QcProgram::Behemoth {
        return native::run(
            &binary,
            dir,
            &atoms,
            &coordinates,
            charge,
            multiplicity,
            &config,
            cancel,
            child,
        );
    }
    let model = ConnectivityModel {
        rcov_disp: atoms
            .iter()
            .map(|a| {
                crate::molecule::covalent_radius_angstrom(a) as f64 * (4. / 3.) / BOHR_TO_ANGSTROM
            })
            .collect(),
        kcn: 7.5,
        facmin: 0.6,
    };
    let x = coordinates.iter().map(|x| x / BOHR_TO_ANGSTROM).collect();
    crate::qchem_panel::calculation_record::append(dir, "Starting Beavyr geometry optimization…");
    let outcome = if program == QcProgram::Dreiding {
        use crate::forcefield::dreiding::{
            objective::{cancellation::CancellableObjective, DreidingObjective},
            DreidingTopology,
        };
        let topology = DreidingTopology::build(&mol)
            .map_err(|e| crate::forcefield::dreiding::diagnostics::message(&e))?;
        let mut objective =
            CancellableObjective::new(DreidingObjective::new(&topology), Some(cancel));
        compute(&mut objective, x, model, &atoms, &config.search, job, dir)?
    } else {
        let mut objective = Evaluator {
            program,
            method,
            binary: &binary,
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
        compute(&mut objective, x, model, &atoms, &config.search, job, dir)?
    };
    check(cancel)?;
    let final_coords: Vec<_> = outcome.x.iter().map(|v| v * BOHR_TO_ANGSTROM).collect();
    let mut output = output(program, job, dir, &atoms, outcome)?;
    if job == JobType::OptimizeThenFrequencies && output.converged == Some(true) {
        let mut frequencies = run_cancellable(
            program,
            method,
            JobType::Frequencies,
            &binary,
            &dir.join("frequencies"),
            &atoms,
            &final_coords,
            charge,
            multiplicity,
            &config,
            Default::default(),
            false,
            None,
            cancel,
            child,
        )?;
        frequencies.job = job;
        frequencies.trajectory_text = output.trajectory_text;
        frequencies.log = format!("{}\n{}", output.log, frequencies.log);
        frequencies.converged = output.converged;
        frequencies.saveable.extend(output.saveable);
        frequencies.kept_dir = Some(dir.to_owned());
        return Ok(frequencies);
    } else if job == JobType::OptimizeThenFrequencies {
        output
            .log
            .push_str("\nFrequencies were skipped because optimization did not converge.");
    }
    Ok(output)
}
fn check(cancel: &AtomicBool) -> Result<(), String> {
    if cancel.load(Ordering::SeqCst) {
        Err("Cancelled.".into())
    } else {
        Ok(())
    }
}
fn compute<O: Objective>(
    objective: &mut O,
    x: Vec<f64>,
    model: ConnectivityModel,
    atoms: &[String],
    settings: &search::Options,
    job: JobType,
    dir: &Path,
) -> Result<Outcome, String> {
    if job != JobType::TransitionState {
        return minimum::compute(objective, x, model, &settings.minimum);
    }
    let mut opts = search::prepare_ts(settings, atoms)?;
    opts.masses_amu = Some(
        atoms
            .iter()
            .map(|a| {
                crate::molecule::atomic_mass_amu(a)
                    .ok_or_else(|| format!("Unknown atomic mass for {a}"))
            })
            .collect::<Result<Vec<_>, _>>()?,
    );
    opts.trajectory_file = Some(dir.join("ts_iterations.xyz").to_string_lossy().into_owned());
    opts.verbosity = 0;
    let result = crate::optimizer::transition_state::optimize_transition_state(
        x,
        Some(model),
        objective,
        opts,
    )
    .map_err(|e| format!("TS search failed: {e:#}"))?;
    let converged =
        result.converged && result.hessian.is_some() && result.negative_modes == Some(1);
    let message = format!(
        "{}\nNegative modes: {}\nFinal Hessian: {}",
        if result.converged && !converged {
            "Stationary geometry reached, but it is not a verified first-order saddle."
        } else {
            &result.message
        },
        result
            .negative_modes
            .map(|n| n.to_string())
            .unwrap_or_else(|| "not determined".into()),
        if result.hessian.is_some() {
            "calculated"
        } else {
            "not recalculated"
        }
    );
    Ok(Outcome {
        x: result.x,
        energy: result.energy,
        converged,
        cycles: result.cycles,
        message,
        frames: result.trajectory,
        energies: result.trajectory_energies,
        hessian: result.hessian,
    })
}
fn output(
    program: QcProgram,
    job: JobType,
    dir: &Path,
    atoms: &[String],
    mut result: Outcome,
) -> Result<QcRunOutput, String> {
    if result.frames.last() != Some(&result.x) {
        result.frames.push(result.x.clone());
        result.energies.push(result.energy);
    }
    if result.frames.len() != result.energies.len() {
        return Err("Geometry and energy histories have different lengths.".into());
    }
    let mut xyz = String::new();
    let mut csv = String::from("step,energy_Eh\n");
    for (i, (x, e)) in result.frames.iter().zip(&result.energies).enumerate() {
        xyz.push_str(&format!(
            "{}\ncycle {} E = {e:.12} Eh\n{}",
            atoms.len(),
            i + 1,
            orca_run::xyz_body(
                atoms,
                &x.iter().map(|x| x * BOHR_TO_ANGSTROM).collect::<Vec<_>>()
            )
        ));
        csv.push_str(&format!("{},{e:.12}\n", i + 1));
    }
    let log = format!(
        "{}\n{}\nSteps: {}\nFinal energy: {:.12} Eh",
        job.label(),
        result.message,
        result.cycles,
        result.energy
    );
    fs::write(dir.join("search.xyz"), &xyz).map_err(|e| e.to_string())?;
    fs::write(dir.join("search.csv"), csv).map_err(|e| e.to_string())?;
    fs::write(dir.join("search_summary.txt"), &log).map_err(|e| e.to_string())?;
    let mut saveable = vec![
        ("Search trajectory", dir.join("search.xyz")),
        ("Search energies", dir.join("search.csv")),
        ("Search summary", dir.join("search_summary.txt")),
    ];
    if let Some(h) = result.hessian {
        let mut text = String::from("$hessian\n");
        for row in h.rows() {
            for value in row {
                text.push_str(&format!(" {value:.14e}"));
            }
            text.push('\n');
        }
        text.push_str("$end\n");
        let path = dir.join("ts_hessian.txt");
        fs::write(&path, text).map_err(|e| e.to_string())?;
        saveable.push(("TS Hessian", path));
    }
    Ok(QcRunOutput {
        job,
        program,
        log,
        trajectory_text: Some(xyz),
        hessian_path: None,
        molden_path: None,
        saveable,
        kept_dir: Some(dir.to_owned()),
        energetics: Energetics {
            final_energy_hartree: Some(result.energy),
            ..Default::default()
        },
        scf: vec![],
        converged: Some(result.converged),
        scan: None,
    })
}

#[cfg(test)]
#[path = "geometry_search_tests.rs"]
mod tests;

/// Honor the existing orbital-export request at the optimized geometry.
#[allow(clippy::too_many_arguments)]
pub(super) fn write_orbitals(
    mut output: QcRunOutput,
    method: Option<QcMethod>,
    binary: &Path,
    dir: &Path,
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    cancel: &AtomicBool,
    child: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    let frames = crate::trajectory::parse_multi_xyz(
        output
            .trajectory_text
            .as_deref()
            .ok_or("No final geometry for orbital export")?,
    )?;
    let frame = frames
        .last()
        .ok_or("No final geometry for orbital export")?;
    let coordinates: Vec<f64> = frame
        .pos
        .iter()
        .flat_map(|p| p.to_array().map(f64::from))
        .collect();
    let result = run_cancellable(
        output.program,
        method,
        JobType::SinglePoint,
        binary,
        &dir.join("final_orbitals"),
        &frame.atoms,
        &coordinates,
        charge,
        multiplicity,
        config,
        Default::default(),
        true,
        None,
        cancel,
        child,
    )?;
    check(cancel)?;
    if result.converged == Some(false) {
        return Err("The geometry search converged, but the final orbital calculation did not converge. Search files remain in the calculation directory.".into());
    }
    output.molden_path = result.molden_path;
    output.saveable.extend(result.saveable);
    Ok(output)
}
