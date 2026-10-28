//! DREIDING jobs for the general calculator, including their result metadata.

use super::{last_frame_xyz, Energetics, QcRunOutput};
use crate::forcefield::dreiding::objective::cancellation::{check_cancelled, CancellableObjective};
use crate::optimizer::minimum_cartesian::CartesianMinimumOptions;
use crate::optimizer::traits::Objective;
use crate::qchem_interfaces::{dreiding_run, job::JobType, program::QcProgram};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};

/// Runs a DREIDING job, here in this process.
///
/// Cancellation is checked between force-field evaluations. Results use the
/// same files as the external backends for the trajectory, spectra and exports.
pub(super) fn run(
    job: JobType,
    workdir: &Path,
    atoms: &[String],
    positions_angstrom: &[f64],
    cancel: &AtomicBool,
) -> Result<QcRunOutput, String> {
    let result = run_with_options(
        job,
        workdir,
        atoms,
        positions_angstrom,
        dreiding_run::optimization_options(),
        cancel,
    );
    if cancel.load(Ordering::SeqCst) {
        Err("Cancelled.".into())
    } else {
        result
    }
}

fn run_with_options(
    job: JobType,
    workdir: &Path,
    atoms: &[String],
    positions_angstrom: &[f64],
    options: CartesianMinimumOptions,
    cancel: &AtomicBool,
) -> Result<QcRunOutput, String> {
    use crate::forcefield::dreiding::hessian::hessian_and_gradient_cancellable;
    use crate::forcefield::dreiding::objective::{DreidingObjective, BOHR_TO_ANGSTROM};
    use crate::forcefield::dreiding::DreidingTopology;
    check_cancelled(Some(cancel)).map_err(|e| e.to_string())?;

    // A saddle search, a reaction path and excited states are all beyond a
    // generic force field. The panel disables them, so reaching here would be
    // a bug rather than a user's mistake -- but it is said plainly either way.
    if matches!(
        job,
        JobType::TransitionState | JobType::Irc | JobType::TdDft
    ) {
        return Err(format!(
            "DREIDING is a force field and cannot do {}. Choose a quantum-chemistry program \
             for that.",
            job.label()
        ));
    }
    if atoms.is_empty()
        || positions_angstrom.len() != atoms.len() * 3
        || positions_angstrom.iter().any(|value| !value.is_finite())
    {
        return Err("DREIDING requires a nonempty structure with finite coordinates.".into());
    }

    let input_xyz = crate::qchem_interfaces::xtb_optimize::write_xyz_string(
        atoms,
        &positions_angstrom
            .chunks_exact(3)
            .map(|p| bevy::prelude::Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32))
            .collect::<Vec<_>>(),
    );

    let mut trajectory_text = None;
    let mut final_energy = None;
    let mut converged = None;
    let mut log = String::new();
    let mut saveable: Vec<(&'static str, PathBuf)> = Vec::new();

    if matches!(job, JobType::SinglePoint | JobType::Gradient) {
        let molecule = dreiding_run::molecule_for_force_field(&input_xyz);
        let topology = DreidingTopology::build(&molecule)
            .map_err(|e| crate::forcefield::dreiding::diagnostics::message(&e))?;
        let coords_bohr: Vec<f64> = molecule
            .pos
            .iter()
            .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
            .map(|value| value / BOHR_TO_ANGSTROM)
            .collect();
        let mut objective =
            CancellableObjective::new(DreidingObjective::new(&topology), Some(cancel));
        let mut gradient = vec![0.0; coords_bohr.len()];
        let energy = if job == JobType::Gradient {
            objective.energy_gradient(&coords_bohr, &mut gradient)
        } else {
            objective.energy(&coords_bohr)
        }
        .map_err(|e| format!("DREIDING calculation failed: {e}"))?;
        if !energy.is_finite() || gradient.iter().any(|value| !value.is_finite()) {
            return Err("DREIDING produced a non-finite energy or gradient.".into());
        }
        final_energy = Some(energy);
        log = format!("DREIDING {}\nEnergy: {energy:.14} Hartree\n", job.label());
        if job == JobType::Gradient {
            let path = workdir.join("dreiding.grad");
            let text = gradient_text(&molecule.atoms, &coords_bohr, energy, &gradient);
            fs::write(&path, text).map_err(|e| format!("Failed to write the gradient: {e}"))?;
            saveable.push(("Gradient (Turbomole $grad)", path));
            log.push_str("Cartesian gradient: Hartree/bohr; coordinates: bohr\n");
        }
    }

    if job.produces_trajectory() {
        let molecule = dreiding_run::molecule_for_force_field(&input_xyz);
        let result =
            dreiding_run::optimize_with_options(workdir, &molecule, options, Some(cancel))?;
        converged = Some(result.converged);
        let output = result.output;
        log = output.log;
        final_energy = output.final_energy_hartree;
        trajectory_text = Some(output.trajectory_text);
    }

    // The Hessian, at whichever geometry we ended up with.
    let hessian_path = if job.produces_frequencies() || job == JobType::Hessian {
        let geometry = trajectory_text
            .as_deref()
            .and_then(last_frame_xyz)
            .unwrap_or_else(|| input_xyz.clone());
        let molecule = crate::qchem_interfaces::dreiding_run::molecule_for_force_field(&geometry);
        let topology = DreidingTopology::build(&molecule)
            .map_err(|e| crate::forcefield::dreiding::diagnostics::message(&e))?;
        let coords_bohr: Vec<f64> = molecule
            .pos
            .iter()
            .flat_map(|p| {
                [
                    p.x as f64 / BOHR_TO_ANGSTROM,
                    p.y as f64 / BOHR_TO_ANGSTROM,
                    p.z as f64 / BOHR_TO_ANGSTROM,
                ]
            })
            .collect();
        let (hessian, _gradient) =
            hessian_and_gradient_cancellable(&topology, &coords_bohr, Some(cancel))
                .map_err(|e| e.to_string())?;

        // Behind a `$hessian` header, the format every other backend's matrix
        // is written in and the one Beavyr's reader takes.
        let mut text = String::from("$hessian\n");
        for row in hessian.rows() {
            for value in row {
                text.push_str(&format!(" {value:22.14e}"));
            }
            text.push('\n');
        }
        text.push_str("$end\n");
        let path = workdir.join("hessian.txt");
        fs::write(&path, text).map_err(|e| format!("Failed to write the Hessian: {e}"))?;
        saveable.push(("Hessian ($hessian)", path.clone()));
        Some(path)
    } else {
        None
    };

    check_cancelled(Some(cancel)).map_err(|e| e.to_string())?;
    Ok(QcRunOutput {
        // A force field has no orbitals.
        molden_path: None,
        job,
        program: QcProgram::Dreiding,
        log,
        energetics: Energetics {
            final_energy_hartree: final_energy,
            ..Default::default()
        },
        // A force field has no self-consistent field to iterate.
        scf: Vec::new(),
        converged,
        trajectory_text,
        kept_dir: (hessian_path.is_some() || !saveable.is_empty()).then(|| workdir.to_path_buf()),
        hessian_path: job
            .produces_frequencies()
            .then(|| hessian_path.clone())
            .flatten(),
        saveable,
    })
}

/// The same Turbomole layout used by the other calculator backends and reader.
fn gradient_text(atoms: &[String], coords_bohr: &[f64], energy: f64, gradient: &[f64]) -> String {
    let norm = gradient.iter().map(|g| g * g).sum::<f64>().sqrt();
    let mut text = format!(
        "$grad          cartesian gradients\n  cycle =      1    SCF energy =  {energy:20.14}   |dE/dxyz| =  {norm:.10}\n"
    );
    for (atom, xyz) in atoms.iter().zip(coords_bohr.chunks_exact(3)) {
        text.push_str(&format!(
            "{:22.14} {:22.14} {:22.14}      {atom}\n",
            xyz[0], xyz[1], xyz[2]
        ));
    }
    for xyz in gradient.chunks_exact(3) {
        text.push_str(&format!(
            "{:22.14e} {:22.14e} {:22.14e}\n",
            xyz[0], xyz[1], xyz[2]
        ));
    }
    text.push_str("$end\n");
    text
}

#[cfg(test)]
#[path = "dreiding_job_tests.rs"]
mod tests;
