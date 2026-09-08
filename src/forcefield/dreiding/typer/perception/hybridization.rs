//! Assigns VSEPR-consistent hybridization states after conjugation and aromaticity analysis.
//!
//! This perception stage translates the structural annotations produced by earlier passes
//! (degree, lone pairs, conjugation, aromatic flags) into concrete `Hybridization` labels and
//! the corresponding steric numbers required by later typing decisions.

use super::model::{AnnotatedAtom, AnnotatedMolecule};
use crate::forcefield::dreiding::typer::core::error::PerceptionError;
use crate::forcefield::dreiding::typer::core::properties::{Element, GraphBondOrder, Hybridization};

/// Updates every atom with its final hybridization and steric number assignments.
///
/// The procedure evaluates the VSEPR steric number for each annotated atom, reconciles it
/// with conjugation and aromatic flags, and records the resulting `Hybridization` along with
/// the adjusted steric number in-place.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule containing degrees, lone pairs, and resonance flags.
///
/// # Errors
///
/// Returns [`PerceptionError::HybridizationInference`] when an atom presents an unsupported
/// steric environment.
pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    for atom in &mut molecule.atoms {
        atom.hybridization = initial_hybridization(atom)?;
    }

    loop {
        let mut changes = 0;
        for i in 0..molecule.atoms.len() {
            if molecule.atoms[i].hybridization == Hybridization::SP3
                && molecule.atoms[i].lone_pairs > 0
                && matches!(molecule.atoms[i].element, Element::O | Element::N)
            {
                let is_adjacent_to_pi_system = molecule.adjacency[i]
                    .iter()
                    .any(|&(neighbor_id, _)| supports_delocalization(molecule, neighbor_id));

                if is_adjacent_to_pi_system {
                    molecule.atoms[i].hybridization = Hybridization::Resonant;
                    molecule.atoms[i].is_resonant = true;
                    changes += 1;
                }
            }
        }
        if changes == 0 {
            break;
        }
    }

    for atom in &mut molecule.atoms {
        atom.steric_number = match atom.hybridization {
            Hybridization::Resonant | Hybridization::SP2 => 3,
            Hybridization::SP3 => 4,
            Hybridization::SP => 2,
            Hybridization::None | Hybridization::Unknown => atom.degree + atom.lone_pairs,
        };
    }

    Ok(())
}

/// Determines the initial hybridization for a given atom, respecting resonance flags
/// before applying pure VSEPR steric-number logic.
fn initial_hybridization(atom: &AnnotatedAtom) -> Result<Hybridization, PerceptionError> {
    if is_non_hybridized_element(atom.element) {
        return Ok(Hybridization::None);
    }

    if atom.is_resonant && !atom.is_anti_aromatic {
        return Ok(Hybridization::Resonant);
    }

    let steric_number = atom.degree + atom.lone_pairs;
    match steric_number {
        4 => Ok(Hybridization::SP3),
        3 => Ok(Hybridization::SP2),
        2 => Ok(Hybridization::SP),
        0 | 1 => Ok(Hybridization::None),
        _ => Err(PerceptionError::HybridizationInference { atom_id: atom.id }),
    }
}

/// Detects elements that should stay in the `None` hybridization state regardless of geometry.
fn is_non_hybridized_element(element: Element) -> bool {
    matches!(
        element,
        Element::Li
            | Element::Na
            | Element::K
            | Element::Rb
            | Element::Cs
            | Element::Fr
            | Element::Be
            | Element::Mg
            | Element::Ca
            | Element::Sr
            | Element::Ba
            | Element::Ra
            | Element::F
            | Element::Cl
            | Element::Br
            | Element::I
            | Element::At
            | Element::He
            | Element::Ne
            | Element::Ar
            | Element::Kr
            | Element::Xe
            | Element::Rn
            | Element::Sc
            | Element::Ti
            | Element::V
            | Element::Cr
            | Element::Mn
            | Element::Fe
            | Element::Co
            | Element::Ni
            | Element::Cu
            | Element::Zn
            | Element::Y
            | Element::Zr
            | Element::Nb
            | Element::Mo
            | Element::Tc
            | Element::Ru
            | Element::Rh
            | Element::Pd
            | Element::Ag
            | Element::Cd
            | Element::Hf
            | Element::Ta
            | Element::W
            | Element::Re
            | Element::Os
            | Element::Ir
            | Element::Pt
            | Element::Au
            | Element::Hg
    )
}

/// Checks whether a neighboring atom can participate in π-delocalization.
///
/// Atoms with SP2, SP, or Resonant hybridization can generally extend conjugation,
/// except for carbonyl-like carbons (C=O or C=S) which do not promote adjacent
/// heteroatoms to resonant hybridization.
fn supports_delocalization(molecule: &AnnotatedMolecule, neighbor_id: usize) -> bool {
    let neighbor = &molecule.atoms[neighbor_id];
    if !matches!(
        neighbor.hybridization,
        Hybridization::SP2 | Hybridization::SP | Hybridization::Resonant
    ) {
        return false;
    }

    if neighbor.element == Element::C {
        let is_carbonyl_like = molecule.adjacency[neighbor_id]
            .iter()
            .any(|&(other_id, order)| {
                order == GraphBondOrder::Double
                    && matches!(molecule.atoms[other_id].element, Element::O | Element::S)
            });
        if is_carbonyl_like {
            return false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::dreiding::typer::core::graph::MolecularGraph;
    use crate::forcefield::dreiding::typer::core::properties::{Element, GraphBondOrder};

    fn build_molecule(
        elements: &[Element],
        bonds: &[(usize, usize, GraphBondOrder)],
        setup: impl FnOnce(&mut AnnotatedMolecule),
    ) -> AnnotatedMolecule {
        let mut graph = MolecularGraph::new();
        for &element in elements {
            graph.add_atom(element);
        }
        for &(u, v, order) in bonds {
            graph.add_bond(u, v, order).unwrap();
        }
        let mut molecule = AnnotatedMolecule::new(&graph).unwrap();
        setup(&mut molecule);
        molecule
    }

    #[test]
    fn non_hybridized_elements_remain_none() {
        let mut molecule = build_molecule(
            &[Element::Na, Element::C],
            &[(0, 1, GraphBondOrder::Single)],
            |_| {},
        );
        perceive(&mut molecule).unwrap();
        assert_eq!(molecule.atoms[0].hybridization, Hybridization::None);
    }

    #[test]
    fn pre_marked_resonant_atoms_are_honored() {
        let mut molecule = build_molecule(&[Element::C], &[], |mol| {
            mol.atoms[0].is_resonant = true;
            mol.atoms[0].degree = 3;
        });
        perceive(&mut molecule).unwrap();
        assert_eq!(molecule.atoms[0].hybridization, Hybridization::Resonant);
        assert_eq!(molecule.atoms[0].steric_number, 3);
    }

    #[test]
    fn enol_ether_oxygen_is_corrected_to_resonant() {
        let mut molecule = build_molecule(
            &[Element::C, Element::C, Element::O, Element::C],
            &[
                (0, 1, GraphBondOrder::Double),
                (1, 2, GraphBondOrder::Single),
                (2, 3, GraphBondOrder::Single),
            ],
            |mol| {
                mol.atoms[0].degree = 3;
                mol.atoms[1].degree = 3;
                mol.atoms[2].degree = 2;
                mol.atoms[3].degree = 4;
                mol.atoms[2].lone_pairs = 2;
            },
        );
        perceive(&mut molecule).unwrap();

        assert_eq!(molecule.atoms[0].hybridization, Hybridization::SP2);
        assert_eq!(molecule.atoms[1].hybridization, Hybridization::SP2);
        assert_eq!(
            molecule.atoms[2].hybridization,
            Hybridization::Resonant,
            "Oxygen should be promoted to Resonant"
        );
        assert_eq!(molecule.atoms[3].hybridization, Hybridization::SP3);

        assert!(
            molecule.atoms[2].is_resonant,
            "is_resonant flag should be set for the oxygen"
        );
        assert_eq!(
            molecule.atoms[2].steric_number, 3,
            "Steric number of resonant oxygen should be 3"
        );
    }

    #[test]
    fn aniline_nitrogen_is_corrected_to_resonant() {
        let mut molecule = build_molecule(
            &[Element::C, Element::N],
            &[(0, 1, GraphBondOrder::Single)],
            |mol| {
                mol.atoms[0].is_resonant = true;
                mol.atoms[1].lone_pairs = 1;
                mol.atoms[1].degree = 3;
            },
        );
        perceive(&mut molecule).unwrap();

        assert_eq!(molecule.atoms[1].hybridization, Hybridization::Resonant);
        assert!(molecule.atoms[1].is_resonant);
        assert_eq!(molecule.atoms[1].steric_number, 3);
    }

    #[test]
    fn vsepr_rules_assign_expected_hybridizations() {
        let mut molecule = build_molecule(
            &[Element::C, Element::C, Element::C, Element::H],
            &[],
            |mol| {
                mol.atoms[0].degree = 4;
                mol.atoms[0].lone_pairs = 0;
                mol.atoms[1].degree = 3;
                mol.atoms[1].lone_pairs = 0;
                mol.atoms[2].degree = 2;
                mol.atoms[2].lone_pairs = 0;
                mol.atoms[3].degree = 1;
                mol.atoms[3].lone_pairs = 0;
            },
        );
        perceive(&mut molecule).unwrap();

        assert_eq!(molecule.atoms[0].hybridization, Hybridization::SP3);
        assert_eq!(molecule.atoms[0].steric_number, 4);

        assert_eq!(molecule.atoms[1].hybridization, Hybridization::SP2);
        assert_eq!(molecule.atoms[1].steric_number, 3);

        assert_eq!(molecule.atoms[2].hybridization, Hybridization::SP);
        assert_eq!(molecule.atoms[2].steric_number, 2);

        assert_eq!(molecule.atoms[3].hybridization, Hybridization::None);
    }

    #[test]
    fn anti_aromatic_atoms_are_not_promoted_to_resonant() {
        let mut molecule = build_molecule(&[Element::C], &[], |mol| {
            mol.atoms[0].is_resonant = true;
            mol.atoms[0].is_anti_aromatic = true;
            mol.atoms[0].degree = 3;
        });
        perceive(&mut molecule).unwrap();

        assert_eq!(molecule.atoms[0].hybridization, Hybridization::SP2);
        assert_eq!(molecule.atoms[0].steric_number, 3);
    }

    #[test]
    fn steric_numbers_above_four_raise_an_error() {
        let mut molecule = build_molecule(&[Element::S], &[], |mol| {
            mol.atoms[0].degree = 6;
        });
        let err = perceive(&mut molecule).expect_err("steric 6 should fail");

        match err {
            PerceptionError::HybridizationInference { atom_id } => assert_eq!(atom_id, 0),
            other => panic!("unexpected error returned: {other:?}"),
        }
    }

    #[test]
    fn carbonyl_carbon_does_not_propagate_resonance_to_adjacent_oxygen() {
        let mut molecule = build_molecule(
            &[Element::O, Element::C, Element::O],
            &[
                (0, 1, GraphBondOrder::Single),
                (1, 2, GraphBondOrder::Double),
            ],
            |mol| {
                mol.atoms[0].degree = 2;
                mol.atoms[0].lone_pairs = 2;
                mol.atoms[1].degree = 3;
                mol.atoms[1].lone_pairs = 0;
                mol.atoms[2].degree = 1;
                mol.atoms[2].lone_pairs = 2;
            },
        );
        perceive(&mut molecule).unwrap();

        assert_eq!(
            molecule.atoms[1].hybridization,
            Hybridization::SP2,
            "Carbonyl carbon should be SP2"
        );
        assert_eq!(
            molecule.atoms[0].hybridization,
            Hybridization::SP3,
            "Ether oxygen adjacent to C=O should NOT be promoted to Resonant"
        );
        assert!(
            !molecule.atoms[0].is_resonant,
            "Ether oxygen should not be flagged as resonant"
        );
    }

    #[test]
    fn thioester_carbon_does_not_propagate_resonance() {
        let mut molecule = build_molecule(
            &[Element::S, Element::C, Element::S],
            &[
                (0, 1, GraphBondOrder::Single),
                (1, 2, GraphBondOrder::Double),
            ],
            |mol| {
                mol.atoms[0].degree = 2;
                mol.atoms[0].lone_pairs = 2;
                mol.atoms[1].degree = 3;
                mol.atoms[1].lone_pairs = 0;
                mol.atoms[2].degree = 1;
                mol.atoms[2].lone_pairs = 2;
            },
        );
        perceive(&mut molecule).unwrap();

        assert_eq!(
            molecule.atoms[1].hybridization,
            Hybridization::SP2,
            "Thiocarbonyl carbon should be SP2"
        );
        assert_eq!(
            molecule.atoms[0].hybridization,
            Hybridization::SP3,
            "Sulfur adjacent to C=S should NOT be promoted to Resonant"
        );
    }

    #[test]
    fn sp3_neighbor_does_not_propagate_resonance() {
        let mut molecule = build_molecule(
            &[Element::C, Element::O],
            &[(0, 1, GraphBondOrder::Single)],
            |mol| {
                mol.atoms[0].degree = 4;
                mol.atoms[0].lone_pairs = 0;
                mol.atoms[1].degree = 2;
                mol.atoms[1].lone_pairs = 2;
            },
        );
        perceive(&mut molecule).unwrap();

        assert_eq!(molecule.atoms[0].hybridization, Hybridization::SP3);
        assert_eq!(
            molecule.atoms[1].hybridization,
            Hybridization::SP3,
            "Oxygen attached to SP3 carbon should remain SP3"
        );
        assert!(!molecule.atoms[1].is_resonant);
    }

    #[test]
    fn nitrogen_with_lone_pair_adjacent_to_sp2_carbon_becomes_resonant() {
        let mut molecule = build_molecule(
            &[Element::C, Element::C, Element::N],
            &[
                (0, 1, GraphBondOrder::Double),
                (1, 2, GraphBondOrder::Single),
            ],
            |mol| {
                mol.atoms[0].degree = 3;
                mol.atoms[0].lone_pairs = 0;
                mol.atoms[1].degree = 3;
                mol.atoms[1].lone_pairs = 0;
                mol.atoms[2].degree = 3;
                mol.atoms[2].lone_pairs = 1;
            },
        );
        perceive(&mut molecule).unwrap();

        assert_eq!(molecule.atoms[0].hybridization, Hybridization::SP2);
        assert_eq!(molecule.atoms[1].hybridization, Hybridization::SP2);
        assert_eq!(
            molecule.atoms[2].hybridization,
            Hybridization::Resonant,
            "Enamine nitrogen should be promoted to Resonant"
        );
        assert!(molecule.atoms[2].is_resonant);
    }

    #[test]
    fn resonance_propagates_through_chain() {
        let mut molecule = build_molecule(
            &[Element::C, Element::C, Element::O, Element::C],
            &[
                (0, 1, GraphBondOrder::Double),
                (1, 2, GraphBondOrder::Single),
                (2, 3, GraphBondOrder::Single),
            ],
            |mol| {
                mol.atoms[0].degree = 3;
                mol.atoms[1].degree = 3;
                mol.atoms[2].degree = 2;
                mol.atoms[2].lone_pairs = 2;
                mol.atoms[3].degree = 4;
            },
        );
        perceive(&mut molecule).unwrap();

        assert_eq!(molecule.atoms[0].hybridization, Hybridization::SP2);
        assert_eq!(molecule.atoms[1].hybridization, Hybridization::SP2);
        assert_eq!(
            molecule.atoms[2].hybridization,
            Hybridization::Resonant,
            "Oxygen adjacent to SP2 carbon should become Resonant"
        );
        assert_eq!(molecule.atoms[3].hybridization, Hybridization::SP3);
    }

    #[test]
    fn only_oxygen_and_nitrogen_are_promoted_to_resonant() {
        let mut molecule = build_molecule(
            &[Element::C, Element::C],
            &[(0, 1, GraphBondOrder::Single)],
            |mol| {
                mol.atoms[0].degree = 3;
                mol.atoms[0].lone_pairs = 0;
                mol.atoms[0].is_resonant = true;
                mol.atoms[1].degree = 3;
                mol.atoms[1].lone_pairs = 1;
            },
        );
        perceive(&mut molecule).unwrap();

        assert_eq!(molecule.atoms[0].hybridization, Hybridization::Resonant);
        assert_eq!(
            molecule.atoms[1].hybridization,
            Hybridization::SP3,
            "Carbon with lone pair should NOT be promoted to Resonant"
        );
        assert!(
            !molecule.atoms[1].is_resonant,
            "Carbanion carbon should not gain is_resonant flag"
        );
    }

    #[test]
    fn heteroatom_without_lone_pairs_is_not_promoted() {
        let mut molecule = build_molecule(
            &[Element::C, Element::N],
            &[(0, 1, GraphBondOrder::Single)],
            |mol| {
                mol.atoms[0].degree = 3;
                mol.atoms[0].is_resonant = true;
                mol.atoms[1].degree = 4;
                mol.atoms[1].lone_pairs = 0;
            },
        );
        perceive(&mut molecule).unwrap();

        assert_eq!(
            molecule.atoms[1].hybridization,
            Hybridization::SP3,
            "Quaternary nitrogen should remain SP3"
        );
        assert!(!molecule.atoms[1].is_resonant);
    }

    #[test]
    fn adjacent_to_resonant_atom_promotes_heteroatom() {
        let mut molecule = build_molecule(
            &[Element::C, Element::O],
            &[(0, 1, GraphBondOrder::Single)],
            |mol| {
                mol.atoms[0].degree = 3;
                mol.atoms[0].is_resonant = true;
                mol.atoms[0].hybridization = Hybridization::Resonant;
                mol.atoms[1].degree = 2;
                mol.atoms[1].lone_pairs = 2;
            },
        );
        perceive(&mut molecule).unwrap();

        assert_eq!(
            molecule.atoms[1].hybridization,
            Hybridization::Resonant,
            "Phenol oxygen should be promoted to Resonant"
        );
        assert!(molecule.atoms[1].is_resonant);
        assert_eq!(molecule.atoms[1].steric_number, 3);
    }

    #[test]
    fn adjacent_to_sp_carbon_promotes_heteroatom() {
        let mut molecule = build_molecule(
            &[Element::C, Element::C, Element::O],
            &[
                (0, 1, GraphBondOrder::Triple),
                (1, 2, GraphBondOrder::Single),
            ],
            |mol| {
                mol.atoms[0].degree = 2;
                mol.atoms[0].lone_pairs = 0;
                mol.atoms[1].degree = 2;
                mol.atoms[1].lone_pairs = 0;
                mol.atoms[2].degree = 2;
                mol.atoms[2].lone_pairs = 2;
            },
        );
        perceive(&mut molecule).unwrap();

        assert_eq!(molecule.atoms[0].hybridization, Hybridization::SP);
        assert_eq!(molecule.atoms[1].hybridization, Hybridization::SP);
        assert_eq!(
            molecule.atoms[2].hybridization,
            Hybridization::Resonant,
            "Oxygen adjacent to SP carbon should be promoted"
        );
    }
}
