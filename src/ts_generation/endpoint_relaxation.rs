//! Optional endpoint minimization using the same objective as the TS search.
use super::{backend::RunControl, parameters::Algorithm, Endpoint, GenerationOutput};
use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
use crate::optimizer::minimum_cartesian::{minimize_cartesian, CartesianMinimumOptions};
use crate::optimizer::{rda::RdaOptions, traits::Objective};
use anyhow::{bail, Context, Result};
use bevy::prelude::Vec3;
use std::path::Path;

fn relax<O: Objective>(objective: &mut O, endpoint: &Endpoint, label: &str) -> Result<Endpoint> {
    let result = minimize_cartesian(
        endpoint.coordinates_bohr(),
        objective,
        CartesianMinimumOptions {
            verbosity: 0,
            ..Default::default()
        },
    )
    .with_context(|| format!("{label} endpoint relaxation failed"))?;
    if !result.converged {
        bail!("{label} endpoint relaxation did not converge after {} steps. TS generation was not started. {}",
            result.cycles, result.message);
    }
    if !result.energy.is_finite() || result.x.iter().any(|x| !x.is_finite()) {
        bail!("{label} endpoint relaxation returned non-finite coordinates or energy.");
    }
    let mut relaxed = endpoint.clone();
    relaxed.positions = result
        .x
        .chunks_exact(3)
        .map(|q| {
            Vec3::new(
                (q[0] * BOHR_TO_ANGSTROM) as f32,
                (q[1] * BOHR_TO_ANGSTROM) as f32,
                (q[2] * BOHR_TO_ANGSTROM) as f32,
            )
        })
        .collect();
    Ok(relaxed)
}

fn stage(control: &RunControl, label: &'static str) -> Result<()> {
    control.check()?;
    if let Ok(mut stage) = control.stage.lock() {
        *stage = Some(label);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn generate<O: Objective>(
    objective: &mut O,
    reactant: &Endpoint,
    product: &Endpoint,
    options: RdaOptions,
    algorithm: Algorithm,
    method: String,
    run_dir: &Path,
    relax_endpoints: bool,
    control: &RunControl,
) -> Result<GenerationOutput> {
    let (reactant, product) = if relax_endpoints {
        stage(control, "Relaxing reactant endpoint (A)…")?;
        let reactant = relax(objective, reactant, "Reactant")?;
        stage(control, "Relaxing product endpoint (B)…")?;
        let product = relax(objective, product, "Product")?;
        super::validate_endpoints(&reactant, &product)
            .context("Relaxed endpoints cannot define a reaction path")?;
        (reactant, product)
    } else {
        (reactant.clone(), product.clone())
    };
    stage(control, "Generating TS guess…")?;
    super::generate(
        objective, &reactant, &product, options, algorithm, method, run_dir,
    )
}
