//! Remove placeholders from both reaction endpoints and remap atom references.

use super::parameters::{atom_indices, Algorithm, Parameters, ReactionKind};
use super::{Endpoint, TsGeneration};
use anyhow::{anyhow, Context, Result};

fn endpoint(source: &Endpoint) -> Result<(Endpoint, Vec<Option<usize>>)> {
    let molecule = crate::molecule::Molecule {
        atoms: source.atoms.clone(),
        pos: source.positions.clone(),
        bonds: Vec::new(),
        hydrogen_bonds: Vec::new(),
        topology: None,
    };
    let (real, mapping) =
        crate::calculation_geometry::prepare(&molecule).map_err(anyhow::Error::msg)?;
    Ok((
        Endpoint {
            name: source.name.clone(),
            atoms: real.atoms,
            positions: real.pos,
        },
        mapping,
    ))
}

fn indices(text: &str, mapping: &[Option<usize>]) -> Result<String> {
    text.lines().map(|line| {
        atom_indices(line, mapping.len())?.into_iter().map(|i| {
            mapping[i].map(|j| (j+1).to_string()).ok_or_else(|| anyhow!(
                "Reaction atom {} is a dummy placeholder. Select real atoms for the reaction coordinates.", i+1
            ))
        }).collect::<Result<Vec<_>>>().map(|line| line.join(" "))
    }).collect::<Result<Vec<_>>>().map(|lines| lines.join("\n"))
}

pub(super) fn prepare(state: &TsGeneration) -> Result<(Endpoint, Endpoint, Parameters)> {
    let (reactant, mapping) = endpoint(
        state
            .reactant
            .as_ref()
            .context("Capture or load the reactant structure.")?,
    )?;
    let (product, _) = endpoint(
        state
            .product
            .as_ref()
            .context("Capture or load the product structure.")?,
    )?;
    let mut parameters = state.parameters.clone();
    if mapping.iter().any(Option::is_none) {
        match parameters.reaction {
            ReactionKind::ActiveAtoms => {}
            ReactionKind::Custom => {
                parameters.bonds = indices(&parameters.bonds, &mapping)?;
                parameters.angles = indices(&parameters.angles, &mapping)?;
                parameters.dihedrals = indices(&parameters.dihedrals, &mapping)?;
            }
            _ => parameters.indices = indices(&parameters.indices, &mapping)?,
        }
        parameters.active_atoms = indices(&parameters.active_atoms, &mapping)?;
        if parameters.options.align {
            parameters.align_atoms = indices(&parameters.align_atoms, &mapping)?;
        }
        if parameters.algorithm == Algorithm::PoorMansNeb {
            parameters.extra_bonds = indices(&parameters.extra_bonds, &mapping)?;
            parameters.extra_angles = indices(&parameters.extra_angles, &mapping)?;
            parameters.extra_dihedrals = indices(&parameters.extra_dihedrals, &mapping)?;
        }
    }
    Ok((reactant, product, parameters))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn endpoint_removal_remaps_reaction_coordinates() {
        let endpoint = Endpoint {
            name: "with drawing helper".into(),
            atoms: vec!["X".into(), "H".into(), "H".into()],
            positions: vec![
                bevy::prelude::Vec3::ZERO,
                bevy::prelude::Vec3::X,
                bevy::prelude::Vec3::Y,
            ],
        };
        let mut state = TsGeneration::default();
        state.reactant = Some(endpoint.clone());
        state.product = Some(endpoint);
        state.parameters.indices = "2 3".into();
        let (reactant, product, parameters) = prepare(&state).unwrap();
        assert_eq!(reactant.atoms, ["H", "H"]);
        assert_eq!(product.atoms, reactant.atoms);
        assert_eq!(parameters.indices, "1 2");
        state.parameters.indices = "1 3".into();
        assert!(prepare(&state)
            .err()
            .unwrap()
            .to_string()
            .contains("dummy placeholder"));
    }
}
