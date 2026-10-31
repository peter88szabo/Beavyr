//! Two independently evaluated spin surfaces and the existing MECP optimizer.
use super::calculation_objective::{configure_method, Evaluator};
use super::*;
#[path = "mecp_orca.rs"]
mod native;
use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
use crate::optimizer::internal_coords::ConnectivityModel;
use crate::optimizer::spin_crossing::{
    optimize_spin_crossing, SpinCrossingOptions, SpinCrossingResult,
};

#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    program: QcProgram,
    method: Option<QcMethod>,
    binary: &Path,
    dir: &Path,
    original_atoms: &[String],
    original_coordinates: &[f64],
    charge: i32,
    multiplicity_a: i32,
    config: &MethodConfig,
    cancel: &AtomicBool,
    child: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    let check = || {
        if cancel.load(Ordering::SeqCst) {
            Err("Cancelled.".to_owned())
        } else {
            Ok(())
        }
    };
    check()?;
    let (atoms, coordinates) =
        crate::calculation_geometry::real_coordinates(original_atoms, original_coordinates)?;
    let mol = crate::molecule::Molecule::from_xyz(&format!(
        "{}\nMECP\n{}",
        atoms.len(),
        orca_run::xyz_body(&atoms, &coordinates)
    ));
    crate::qchem_interfaces::job::mecp::validate(
        program,
        method,
        &config.mecp,
        &mol,
        charge,
        multiplicity_a,
    )?;
    if let Some(module) = program.python_module() {
        let environment = crate::qchem_interfaces::python_env::cached(module);
        if !environment.is_ready() {
            return Err(environment.detail());
        }
    }
    if program.needs_binary_path() && binary.as_os_str().is_empty() {
        return Err(format!("Set the {} executable path.", program.label()));
    }
    let absolute_binary = if !binary.is_absolute() && binary.components().count() > 1 {
        std::env::current_dir()
            .map_err(|e| e.to_string())?
            .join(binary)
    } else {
        binary.to_path_buf()
    };
    let mut config = config.clone();
    configure_method(&mut config, method);
    fs::create_dir_all(dir).map_err(|e| e.to_string())?;
    if program == QcProgram::Orca
        && config.mecp.driver == crate::qchem_interfaces::job::mecp::Driver::Orca
    {
        return native::run(
            method.unwrap(),
            &absolute_binary,
            dir,
            &atoms,
            &coordinates,
            charge,
            multiplicity_a,
            &config,
            cancel,
            child,
        );
    }
    let a_dir = dir.join("state_a");
    let b_dir = dir.join("state_b");
    let objective = |state_dir, multiplicity| Evaluator {
        program,
        method,
        binary: &absolute_binary,
        dir: state_dir,
        atoms: &atoms,
        charge,
        multiplicity,
        config: &config,
        cancel,
        child,
        evaluations: 0,
        cached: None,
    };
    let mut a = objective(&a_dir, multiplicity_a);
    let mut b = objective(&b_dir, config.mecp.multiplicity_b);
    let options = SpinCrossingOptions {
        coordinate_mode: config.mecp.coordinates.mode(),
        multiplicities: Some((multiplicity_a as usize, config.mecp.multiplicity_b as usize)),
        max_cycles: config.mecp.max_cycles,
        energy_tol: config.mecp.gap_tolerance,
        trajectory_file: Some(
            dir.join("mecp_iterations.xyz")
                .to_string_lossy()
                .into_owned(),
        ),
        atom_labels: Some(atoms.clone()),
        verbosity: 0,
        ..Default::default()
    };
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
    let result = optimize_spin_crossing(
        coordinates.iter().map(|x| x / BOHR_TO_ANGSTROM).collect(),
        Some(model),
        &mut a,
        &mut b,
        options,
    )
    .map_err(|e| format!("MECP failed: {e:#}"))?;
    check()?;
    let (trajectory, table) = exports(&result, &atoms);
    fs::write(dir.join("mecp.xyz"), &trajectory).map_err(|e| e.to_string())?;
    fs::write(dir.join("mecp.csv"), &table).map_err(|e| e.to_string())?;
    let log = format!("{}\nState A multiplicity: {}\nState B multiplicity: {}\nState A energy: {:.12} Eh\nState B energy: {:.12} Eh\nEnergy gap (A − B): {:.10} Eh\nGap tolerance: {:.10} Eh\nSteps: {}\nCoordinates: {:?}",
        result.message, multiplicity_a, config.mecp.multiplicity_b, result.energy_a, result.energy_b,
        result.energy_a - result.energy_b, config.mecp.gap_tolerance, result.cycles, result.coordinates);
    fs::write(dir.join("mecp_summary.txt"), &log).map_err(|e| e.to_string())?;
    Ok(QcRunOutput {
        job: JobType::Mecp,
        program,
        log,
        trajectory_text: Some(trajectory),
        hessian_path: None,
        molden_path: None,
        saveable: vec![
            ("MECP trajectory", dir.join("mecp.xyz")),
            ("MECP energies", dir.join("mecp.csv")),
            ("MECP summary", dir.join("mecp_summary.txt")),
        ],
        kept_dir: Some(dir.to_path_buf()),
        energetics: Default::default(),
        scf: vec![],
        converged: Some(result.converged),
        scan: None,
    })
}
fn exports(result: &SpinCrossingResult, atoms: &[String]) -> (String, String) {
    let mut trajectory = String::new();
    let mut table =
        String::from("frame,state_a_energy_Eh,state_b_energy_Eh,gap_Eh,mean_energy_Eh\n");
    let points = result
        .points
        .iter()
        .map(|p| (p.x.as_slice(), p.energy_a, p.energy_b))
        .chain(std::iter::once((
            result.x.as_slice(),
            result.energy_a,
            result.energy_b,
        )));
    for (i, (x, ea, eb)) in points.enumerate() {
        let mean = (ea + eb) * 0.5;
        table.push_str(&format!(
            "{},{ea:.12},{eb:.12},{:.12},{mean:.12}\n",
            i + 1,
            ea - eb
        ));
        trajectory.push_str(&format!(
            "{}\ncycle {} E = {mean:.12} Eh Ea={ea:.12} Eb={eb:.12} gap={:.12}\n{}",
            atoms.len(),
            i + 1,
            ea - eb,
            orca_run::xyz_body(
                atoms,
                &x.iter().map(|x| x * BOHR_TO_ANGSTROM).collect::<Vec<_>>()
            )
        ));
    }
    (trajectory, table)
}

#[cfg(test)]
#[path = "mecp_tests.rs"]
mod tests;
