//! Minimum-search dispatch to algorithms already imported from Behemoth.
use super::*;
pub(super) fn compute<O: Objective>(
    objective: &mut O,
    x: Vec<f64>,
    model: ConnectivityModel,
    o: &search::Minimum,
) -> Result<Outcome, String> {
    if o.coordinates == Coordinates::Cartesian {
        use crate::optimizer::minimum_cartesian::*;
        let r = minimize_cartesian(
            x,
            objective,
            CartesianMinimumOptions {
                method: match o.algorithm {
                    MinimumAlgorithm::Bfgs => CartesianMinimumMethod::Bfgs,
                    MinimumAlgorithm::SteepestDescent => CartesianMinimumMethod::SteepestDescent,
                    MinimumAlgorithm::Bb1 => CartesianMinimumMethod::TwoPointGradient1,
                    MinimumAlgorithm::Bb2 => CartesianMinimumMethod::TwoPointGradient2,
                },
                max_cycles: o.max_cycles,
                energy_tol: o.energy_tol,
                max_step: o.max_step,
                rms_step: o.rms_step,
                max_gradient: o.max_gradient,
                rms_gradient: o.rms_gradient,
                min_line_search_scale: o.min_line_search_scale,
                verbosity: 0,
                ..Default::default()
            },
        )
        .map_err(|e| format!("Minimum search failed: {e:#}"))?;
        return Ok(Outcome {
            x: r.x,
            energy: r.energy,
            converged: r.converged,
            cycles: r.cycles,
            message: r.message,
            frames: r.trajectory,
            energies: r.trajectory_energies,
            hessian: None,
        });
    }
    if o.coordinates == Coordinates::Primitive && o.algorithm == MinimumAlgorithm::Bfgs {
        use crate::optimizer::geom_opt::*;
        let r = geom_opt_internal_bfgs(
            x,
            model,
            objective,
            GeomOptOptions {
                max_cycles: o.max_cycles,
                energy_tol: o.energy_tol,
                max_step: o.max_step,
                rms_step: o.rms_step,
                max_gradient: o.max_gradient,
                g_rms_tol: o.rms_gradient,
                max_step_internal: o.max_step_internal,
                best_fit_iters: o.best_fit_iters,
                verbosity: 0,
            },
        )
        .map_err(|e| format!("Minimum search failed: {e:#}"))?;
        return Ok(Outcome {
            x: r.x,
            energy: r.energy,
            converged: r.converged,
            cycles: r.cycles,
            message: if r.converged {
                "Minimum-search convergence reached."
            } else {
                "Maximum steps reached without convergence."
            }
            .into(),
            frames: r.trajectory,
            energies: r.trajectory_energies,
            hessian: None,
        });
    }
    use crate::optimizer::minimum_redundant::*;
    let r = minimize_internal(
        x,
        model,
        objective,
        InternalMinimumOptions {
            method: if o.algorithm == MinimumAlgorithm::SteepestDescent {
                InternalMinimumMethod::SteepestDescent
            } else {
                InternalMinimumMethod::Bfgs
            },
            coordinate_mode: if o.coordinates == Coordinates::Primitive {
                InternalCoordinateMode::SimpleRedundant
            } else {
                InternalCoordinateMode::PulayRedundant
            },
            max_cycles: o.max_cycles,
            energy_tol: o.energy_tol,
            max_step: o.max_step,
            rms_step: o.rms_step,
            max_gradient: o.max_gradient,
            rms_gradient: o.rms_gradient,
            max_step_internal: o.max_step_internal,
            best_fit_iters: o.best_fit_iters,
            best_fit_rms_tol: o.best_fit_rms_tol,
            min_line_search_scale: o.min_line_search_scale,
            g_cutoff: o.g_cutoff,
            redundant_penalty: o.redundant_penalty,
            verbosity: 0,
            ..Default::default()
        },
    )
    .map_err(|e| format!("Minimum search failed: {e:#}"))?;
    Ok(Outcome {
        x: r.x,
        energy: r.energy,
        converged: r.converged,
        cycles: r.cycles,
        message: r.message,
        frames: r.trajectory,
        energies: r.trajectory_energies,
        hessian: None,
    })
}
