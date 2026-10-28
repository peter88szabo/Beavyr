//! Select the requested energy model without requiring another model's atom types.

use super::*;
use std::path::PathBuf;

pub(super) fn objective<'a>(
    settings: &ConformerSettings,
    topology: Option<&'a DreidingTopology>,
    atoms: &'a [String],
    workdir: PathBuf,
) -> anyhow::Result<SearchObjective<'a>> {
    if let Some(args) = settings.engine.xtb_args() {
        let binary = settings.xtb_binary.as_ref().ok_or_else(|| {
            anyhow::anyhow!(
                "{} needs the xTB executable, which is not set",
                settings.engine.label()
            )
        })?;
        return Ok(SearchObjective::Xtb(refine::XtbObjective::new(
            workdir,
            atoms,
            binary,
            args,
            settings.charge,
            settings.multiplicity,
        )));
    }
    if settings.engine.needs_behemoth() {
        let binary = settings.behemoth_binary.as_ref().ok_or_else(|| {
            anyhow::anyhow!(
                "{} needs the Behemoth executable, which is not set",
                settings.engine.label()
            )
        })?;
        let method = crate::qchem_interfaces::method::MethodConfig {
            behemoth: crate::qchem_interfaces::method::BehemothMethod::Tasi,
            nproc: 1,
            ..Default::default()
        };
        return Ok(SearchObjective::Behemoth(refine::BehemothObjective::new(
            workdir,
            atoms,
            binary,
            method,
            settings.charge,
            settings.multiplicity,
        )));
    }
    let topology = topology.ok_or_else(|| anyhow::anyhow!("DREIDING topology is missing"))?;
    Ok(SearchObjective::Dreiding(DreidingObjective::new(topology)))
}

/// Classify rigid/ring-only inputs independently of the selected energy model.
pub(super) fn no_torsions(mol: &Molecule, topology: Option<&DreidingTopology>) -> ConformerError {
    if let Some(topology) = topology {
        return nothing_to_rotate(topology);
    }
    let mut adjacency = vec![Vec::new(); mol.atoms.len()];
    for &(a, b, _) in &mol.bonds {
        adjacency[a].push(b);
        adjacency[b].push(a);
    }
    let rings = mol
        .bonds
        .iter()
        .filter(|&&(a, b, _)| {
            let mut seen = vec![false; adjacency.len()];
            seen[a] = true;
            let mut stack = vec![a];
            while let Some(i) = stack.pop() {
                for &j in &adjacency[i] {
                    if (i == a && j == b) || (i == b && j == a) {
                        continue;
                    }
                    if j == b {
                        return true;
                    }
                    if !seen[j] {
                        seen[j] = true;
                        stack.push(j);
                    }
                }
            }
            false
        })
        .count();
    if rings > 0 {
        ConformerError::RingOnly { rings }
    } else {
        ConformerError::Rigid
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn single_external_conformer(
    settings: &ConformerSettings,
    start_bohr: &[f64],
    atoms: &[String],
    elapsed: Duration,
    workers: usize,
    rotors: usize,
    rings: usize,
) -> Result<ConformerOutcome, ConformerError> {
    use crate::optimizer::minimum_cartesian::minimize_cartesian;
    let mut objective = objective(
        settings,
        None,
        atoms,
        external_scratch().join("single_conformer"),
    )
    .map_err(|e| ConformerError::Search(e.to_string()))?;
    let result = minimize_cartesian(start_bohr.to_vec(), &mut objective, cleanup_options())
        .map_err(|e| ConformerError::Search(e.to_string()))?;
    if !result.converged {
        return Err(ConformerError::Search(format!(
            "{} single-conformer optimization did not converge",
            settings.engine.label()
        )));
    }
    Ok(ConformerOutcome {
        conformers: vec![Conformer {
            positions_angstrom: result.x.iter().map(|x| x * BOHR_TO_ANGSTROM).collect(),
            relative_energy_kcal: 0.0,
            population_fraction: 1.0,
        }],
        rotors,
        rings,
        generations: 0,
        structures_tried: 0,
        local_optimizations: 1,
        termination: "only one distinct minimum".into(),
        elapsed,
        workers,
        local_optimizer: LocalOptimizer::Cartesian,
        constrained: 0,
        best_found_in_generation: 0,
        atoms: atoms.to_vec(),
        radical_warning: None,
        refined_with_xtb: false,
    })
}

#[cfg(test)]
#[path = "search_engine_tests.rs"]
mod tests;
