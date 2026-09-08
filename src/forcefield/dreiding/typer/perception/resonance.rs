//! Detects specific, strong resonance systems via strict substructure matching.
//!
//! Unlike generalized conjugation detection, this module uses an allowlist of
//! chemically significant motifs (Carboxylate, Nitro, Guanidinium, Amide).
//! When a motif is found, its atoms are marked `is_resonant`, and the system
//! (atoms + bonds) is recorded to ensure the correct bond order in the topology.

use super::model::{AnnotatedMolecule, ResonanceSystem};
use crate::forcefield::dreiding::typer::core::error::PerceptionError;
use crate::forcefield::dreiding::typer::core::properties::{Element, GraphBondOrder};

/// Runs strict resonance perception in two phases.
///
/// First, it identifies core resonance systems (aromatics are handled upstream) like
/// carboxylates and amides, registering both their atoms and bonds. Second, it
//  propagates the `is_resonant` flag to adjacent heteroatoms that can participate
//  in conjugation.
///
/// Note: Aromatic resonance is handled in the `aromaticity` module. This module
/// focuses on non-aromatic conjugated systems.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule whose atoms will be tagged and systems recorded.
///
/// # Returns
///
/// `Ok(())` always, as this process is infallible.
pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    detect_core_functional_groups(molecule);
    propagate_resonance_to_periphery(molecule);
    Ok(())
}

/// Detects core resonance systems via substructure matching.
fn detect_core_functional_groups(molecule: &mut AnnotatedMolecule) {
    let mut processed = vec![false; molecule.atoms.len()];

    detect_carboxylate_groups(molecule, &mut processed);
    detect_nitro_groups(molecule, &mut processed);
    detect_guanidinium_groups(molecule, &mut processed);
    detect_thiourea_groups(molecule, &mut processed);
    detect_amide_groups(molecule, &mut processed);
    detect_phosphate_groups(molecule, &mut processed);
}

/// Propagates resonance flags to peripheral heteroatoms bonded to resonant systems.
fn propagate_resonance_to_periphery(molecule: &mut AnnotatedMolecule) {
    let mut newly_resonant = Vec::new();

    for i in 0..molecule.atoms.len() {
        let atom = &molecule.atoms[i];

        if atom.is_resonant
            || !matches!(atom.element, Element::O | Element::N | Element::S)
            || atom.lone_pairs == 0
        {
            continue;
        }

        let is_bonded_to_resonant_atom = molecule.adjacency[i]
            .iter()
            .any(|&(neighbor_id, _)| molecule.atoms[neighbor_id].is_resonant);

        if is_bonded_to_resonant_atom {
            newly_resonant.push(i);
        }
    }

    for atom_id in newly_resonant {
        molecule.atoms[atom_id].is_resonant = true;
    }
}

/// Helper to record a detected resonance system for later topology emission.
fn push_resonance_system(molecule: &mut AnnotatedMolecule, atoms: &[usize], bonds: &[usize]) {
    let mut sys_atoms = atoms.to_vec();
    let mut sys_bonds = bonds.to_vec();
    sys_atoms.sort_unstable();
    sys_bonds.sort_unstable();

    molecule.resonance_systems.push(ResonanceSystem {
        atom_ids: sys_atoms,
        bond_ids: sys_bonds,
    });
}

/// Finds the bond ID connecting two atoms. Panics if not found (internal consistency check).
fn find_bond_id(molecule: &AnnotatedMolecule, u: usize, v: usize) -> usize {
    molecule
        .bonds
        .iter()
        .find(|b| {
            (b.atom_ids.0 == u && b.atom_ids.1 == v) || (b.atom_ids.0 == v && b.atom_ids.1 == u)
        })
        .map(|b| b.id)
        .expect("Bond must exist between adjacent atoms")
}

/// Detects Carboxylate groups: C(=O)O-
fn detect_carboxylate_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for c_idx in 0..molecule.atoms.len() {
        if processed[c_idx] || molecule.atoms[c_idx].element != Element::C {
            continue;
        }

        let mut double_o = None;
        let mut single_o = None;

        for &(neighbor_id, order) in &molecule.adjacency[c_idx] {
            if molecule.atoms[neighbor_id].element == Element::O {
                match order {
                    GraphBondOrder::Double => double_o = Some(neighbor_id),
                    GraphBondOrder::Single if molecule.atoms[neighbor_id].degree == 1 => {
                        single_o = Some(neighbor_id);
                    }
                    _ => {}
                }
            }
        }

        if let (Some(o1), Some(o2)) = (double_o, single_o) {
            let b1 = find_bond_id(molecule, c_idx, o1);
            let b2 = find_bond_id(molecule, c_idx, o2);
            for &atom_id in &[c_idx, o1, o2] {
                molecule.atoms[atom_id].is_resonant = true;
                processed[atom_id] = true;
            }
            push_resonance_system(molecule, &[c_idx, o1, o2], &[b1, b2]);
        }
    }
}

/// Detects Nitro groups: N(=O)O-
fn detect_nitro_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for n_idx in 0..molecule.atoms.len() {
        if processed[n_idx] || molecule.atoms[n_idx].element != Element::N {
            continue;
        }

        let mut oxygen_neighbors = Vec::new();
        for &(neighbor_id, _) in &molecule.adjacency[n_idx] {
            if molecule.atoms[neighbor_id].element == Element::O {
                oxygen_neighbors.push(neighbor_id);
            }
        }

        if oxygen_neighbors.len() == 2 {
            let o1 = oxygen_neighbors[0];
            let o2 = oxygen_neighbors[1];
            let b1 = find_bond_id(molecule, n_idx, o1);
            let b2 = find_bond_id(molecule, n_idx, o2);
            for &atom_id in &[n_idx, o1, o2] {
                molecule.atoms[atom_id].is_resonant = true;
                processed[atom_id] = true;
            }
            push_resonance_system(molecule, &[n_idx, o1, o2], &[b1, b2]);
        }
    }
}

/// Detects Guanidinium groups: C(N)(N)N+
fn detect_guanidinium_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for c_idx in 0..molecule.atoms.len() {
        if processed[c_idx]
            || molecule.atoms[c_idx].element != Element::C
            || molecule.atoms[c_idx].is_in_ring
        {
            continue;
        }

        let n_neighbors: Vec<_> = molecule.adjacency[c_idx]
            .iter()
            .filter(|(id, _)| molecule.atoms[*id].element == Element::N)
            .map(|(id, _)| *id)
            .collect();

        if n_neighbors.len() == 3 {
            let bonds: Vec<_> = n_neighbors
                .iter()
                .map(|&n_id| find_bond_id(molecule, c_idx, n_id))
                .collect();
            let mut atoms = vec![c_idx];
            atoms.extend(n_neighbors);
            for &atom_id in &atoms {
                molecule.atoms[atom_id].is_resonant = true;
                processed[atom_id] = true;
            }
            push_resonance_system(molecule, &atoms, &bonds);
        }
    }
}

/// Detects Thiourea cores: C(=S)(N)(N)
fn detect_thiourea_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for c_idx in 0..molecule.atoms.len() {
        if processed[c_idx] || molecule.atoms[c_idx].element != Element::C {
            continue;
        }

        let mut sulfur_neighbor = None;
        let mut nitrogen_neighbors = Vec::new();

        for &(neighbor_id, order) in &molecule.adjacency[c_idx] {
            match (molecule.atoms[neighbor_id].element, order) {
                (Element::S, GraphBondOrder::Double) => sulfur_neighbor = Some(neighbor_id),
                (Element::N, GraphBondOrder::Single) => nitrogen_neighbors.push(neighbor_id),
                _ => {}
            }
        }

        if sulfur_neighbor.is_some() && nitrogen_neighbors.len() == 2 {
            let n1 = nitrogen_neighbors[0];
            let n2 = nitrogen_neighbors[1];
            let b1 = find_bond_id(molecule, c_idx, n1);
            let b2 = find_bond_id(molecule, c_idx, n2);

            for &atom_id in &[c_idx, n1, n2] {
                molecule.atoms[atom_id].is_resonant = true;
                processed[atom_id] = true;
            }

            push_resonance_system(molecule, &[c_idx, n1, n2], &[b1, b2]);
        }
    }
}

/// Detects Amide groups: O=C-N
fn detect_amide_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for c_idx in 0..molecule.atoms.len() {
        if processed[c_idx] || molecule.atoms[c_idx].element != Element::C {
            continue;
        }

        let mut double_o = None;
        let mut single_n = None;

        for &(neighbor_id, order) in &molecule.adjacency[c_idx] {
            match (molecule.atoms[neighbor_id].element, order) {
                (Element::O, GraphBondOrder::Double) => double_o = Some(neighbor_id),
                (Element::N, GraphBondOrder::Single) => single_n = Some(neighbor_id),
                _ => {}
            }
        }

        if let (Some(o), Some(n)) = (double_o, single_n) {
            if molecule.atoms[c_idx].is_resonant || molecule.atoms[n].is_resonant {
                continue;
            }
            let b_co = find_bond_id(molecule, c_idx, o);
            let b_cn = find_bond_id(molecule, c_idx, n);
            for &atom_id in &[c_idx, o, n] {
                molecule.atoms[atom_id].is_resonant = true;
                processed[atom_id] = true;
            }
            push_resonance_system(molecule, &[c_idx, o, n], &[b_co, b_cn]);
        }
    }
}

/// Detects Phosphate groups: P(=O)(O-)
fn detect_phosphate_groups(molecule: &mut AnnotatedMolecule, processed: &mut [bool]) {
    for p_idx in 0..molecule.atoms.len() {
        if molecule.atoms[p_idx].element != Element::P {
            continue;
        }

        let mut terminal_oxygens = Vec::new();

        for &(neighbor_id, _) in &molecule.adjacency[p_idx] {
            if molecule.atoms[neighbor_id].element == Element::O
                && molecule.atoms[neighbor_id].degree == 1
            {
                terminal_oxygens.push(neighbor_id);
            }
        }

        if terminal_oxygens.len() >= 2 {
            let o1 = terminal_oxygens[0];
            let o2 = terminal_oxygens[1];

            let b1 = find_bond_id(molecule, p_idx, o1);
            let b2 = find_bond_id(molecule, p_idx, o2);

            molecule.atoms[o1].is_resonant = true;
            molecule.atoms[o2].is_resonant = true;
            processed[o1] = true;
            processed[o2] = true;
            processed[p_idx] = true;

            push_resonance_system(molecule, &[p_idx, o1, o2], &[b1, b2]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::dreiding::typer::core::graph::MolecularGraph;
    use crate::forcefield::dreiding::typer::core::properties::Element;
    use std::collections::HashSet;

    fn build_molecule(
        elements: &[Element],
        bonds: &[(usize, usize, GraphBondOrder)],
        lone_pairs: &[(usize, u8)],
    ) -> AnnotatedMolecule {
        let mut graph = MolecularGraph::new();
        for &elem in elements {
            graph.add_atom(elem);
        }
        for &(a, b, order) in bonds {
            graph.add_bond(a, b, order).expect("valid bond");
        }
        let mut molecule = AnnotatedMolecule::new(&graph).expect("valid graph");
        for &(atom_id, lp) in lone_pairs {
            molecule.atoms[atom_id].lone_pairs = lp;
        }
        molecule
    }

    fn run_resonance_perception(mut molecule: AnnotatedMolecule) -> AnnotatedMolecule {
        perceive(&mut molecule).expect("resonance perception should succeed");
        molecule
    }

    fn resonant_atom_ids(molecule: &AnnotatedMolecule) -> HashSet<usize> {
        molecule
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(i, a)| a.is_resonant.then_some(i))
            .collect()
    }

    fn assert_resonant_atoms(molecule: &AnnotatedMolecule, expected: &[usize]) {
        let actual = resonant_atom_ids(molecule);
        let expected_set: HashSet<_> = expected.iter().copied().collect();
        assert_eq!(
            actual, expected_set,
            "resonant atoms mismatch: got {:?}, expected {:?}",
            actual, expected_set
        );
    }

    fn assert_resonance_system_count(molecule: &AnnotatedMolecule, expected: usize) {
        assert_eq!(
            molecule.resonance_systems.len(),
            expected,
            "expected {} resonance system(s), got {}",
            expected,
            molecule.resonance_systems.len()
        );
    }

    fn assert_system_contains_atoms(
        molecule: &AnnotatedMolecule,
        system_idx: usize,
        atoms: &[usize],
    ) {
        let sys = &molecule.resonance_systems[system_idx];
        let expected: HashSet<_> = atoms.iter().copied().collect();
        let actual: HashSet<_> = sys.atom_ids.iter().copied().collect();
        assert_eq!(
            actual, expected,
            "system {} atoms mismatch: got {:?}, expected {:?}",
            system_idx, actual, expected
        );
    }

    #[test]
    fn carboxylate_group_is_detected() {
        let elements = [
            Element::C,
            Element::C,
            Element::O,
            Element::O,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = [
            (0, 1, GraphBondOrder::Single),
            (1, 2, GraphBondOrder::Double),
            (1, 3, GraphBondOrder::Single),
            (0, 4, GraphBondOrder::Single),
            (0, 5, GraphBondOrder::Single),
            (0, 6, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonant_atoms(&molecule, &[1, 2, 3]);
        assert_resonance_system_count(&molecule, 1);
        assert_system_contains_atoms(&molecule, 0, &[1, 2, 3]);
    }

    #[test]
    fn ester_is_not_detected_as_carboxylate() {
        let elements = [Element::C, Element::C, Element::O, Element::O, Element::C];
        let bonds = [
            (0, 1, GraphBondOrder::Single),
            (1, 2, GraphBondOrder::Double),
            (1, 3, GraphBondOrder::Single),
            (3, 4, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        let resonant = resonant_atom_ids(&molecule);
        assert!(
            !resonant.contains(&3),
            "ester oxygen (degree 2) should not be in resonance system"
        );
    }

    #[test]
    fn nitro_group_is_detected() {
        let elements = [
            Element::C,
            Element::N,
            Element::O,
            Element::O,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = [
            (0, 1, GraphBondOrder::Single),
            (1, 2, GraphBondOrder::Double),
            (1, 3, GraphBondOrder::Single),
            (0, 4, GraphBondOrder::Single),
            (0, 5, GraphBondOrder::Single),
            (0, 6, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonant_atoms(&molecule, &[1, 2, 3]);
        assert_resonance_system_count(&molecule, 1);
        assert_system_contains_atoms(&molecule, 0, &[1, 2, 3]);
    }

    #[test]
    fn nitro_detection_is_kekule_invariant() {
        let elements = [Element::C, Element::N, Element::O, Element::O];
        let bonds = [
            (0, 1, GraphBondOrder::Single),
            (1, 2, GraphBondOrder::Single),
            (1, 3, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonant_atoms(&molecule, &[1, 2, 3]);
    }

    #[test]
    fn guanidinium_group_is_detected() {
        let elements = [
            Element::C,
            Element::N,
            Element::N,
            Element::N,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (0, 2, GraphBondOrder::Single),
            (0, 3, GraphBondOrder::Single),
            (1, 4, GraphBondOrder::Single),
            (1, 5, GraphBondOrder::Single),
            (2, 6, GraphBondOrder::Single),
            (2, 7, GraphBondOrder::Single),
            (3, 8, GraphBondOrder::Single),
            (3, 9, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonant_atoms(&molecule, &[0, 1, 2, 3]);
        assert_resonance_system_count(&molecule, 1);
        assert_system_contains_atoms(&molecule, 0, &[0, 1, 2, 3]);
    }

    #[test]
    fn ring_guanidine_is_not_detected() {
        let elements = [
            Element::C,
            Element::N,
            Element::N,
            Element::N,
            Element::C,
            Element::C,
        ];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (0, 2, GraphBondOrder::Single),
            (0, 3, GraphBondOrder::Single),
            (2, 4, GraphBondOrder::Single),
            (3, 5, GraphBondOrder::Single),
            (4, 5, GraphBondOrder::Single),
        ];
        let mut molecule = build_molecule(&elements, &bonds, &[]);
        molecule.atoms[0].is_in_ring = true;
        let molecule = run_resonance_perception(molecule);

        let resonant = resonant_atom_ids(&molecule);
        assert!(
            !resonant.contains(&0),
            "ring guanidine central carbon should not be detected"
        );
    }

    #[test]
    fn thiourea_group_is_detected() {
        let elements = [
            Element::S,
            Element::C,
            Element::N,
            Element::N,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (1, 2, GraphBondOrder::Single),
            (1, 3, GraphBondOrder::Single),
            (2, 4, GraphBondOrder::Single),
            (2, 5, GraphBondOrder::Single),
            (3, 6, GraphBondOrder::Single),
            (3, 7, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonant_atoms(&molecule, &[1, 2, 3]);
        assert_resonance_system_count(&molecule, 1);
        assert_system_contains_atoms(&molecule, 0, &[1, 2, 3]);
    }

    #[test]
    fn thioamide_with_one_nitrogen_is_not_thiourea() {
        let elements = [Element::S, Element::C, Element::N, Element::C];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (1, 2, GraphBondOrder::Single),
            (1, 3, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonance_system_count(&molecule, 0);
    }

    #[test]
    fn amide_group_is_detected() {
        let elements = [
            Element::C,
            Element::C,
            Element::O,
            Element::N,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = [
            (0, 1, GraphBondOrder::Single),
            (1, 2, GraphBondOrder::Double),
            (1, 3, GraphBondOrder::Single),
            (0, 4, GraphBondOrder::Single),
            (0, 5, GraphBondOrder::Single),
            (0, 6, GraphBondOrder::Single),
            (3, 7, GraphBondOrder::Single),
            (3, 8, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonant_atoms(&molecule, &[1, 2, 3]);
        assert_resonance_system_count(&molecule, 1);
        assert_system_contains_atoms(&molecule, 0, &[1, 2, 3]);
    }

    #[test]
    fn carboxylate_takes_priority_over_amide() {
        let elements = [Element::C, Element::O, Element::O, Element::N];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (0, 2, GraphBondOrder::Single),
            (0, 3, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonant_atoms(&molecule, &[0, 1, 2]);
        assert_resonance_system_count(&molecule, 1);
    }

    #[test]
    fn phosphate_group_is_detected() {
        let elements = [
            Element::P,
            Element::O,
            Element::O,
            Element::O,
            Element::O,
            Element::C,
            Element::C,
        ];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (0, 2, GraphBondOrder::Single),
            (0, 3, GraphBondOrder::Single),
            (0, 4, GraphBondOrder::Single),
            (3, 5, GraphBondOrder::Single),
            (4, 6, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert!(
            resonant_atom_ids(&molecule).contains(&1),
            "terminal oxygen O1 should be resonant"
        );
        assert!(
            resonant_atom_ids(&molecule).contains(&2),
            "terminal oxygen O2 should be resonant"
        );
        assert_resonance_system_count(&molecule, 1);
        assert_system_contains_atoms(&molecule, 0, &[0, 1, 2]);
    }

    #[test]
    fn phosphate_requires_two_terminal_oxygens() {
        let elements = [
            Element::P,
            Element::O,
            Element::O,
            Element::O,
            Element::C,
            Element::C,
            Element::C,
        ];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (0, 2, GraphBondOrder::Single),
            (0, 3, GraphBondOrder::Single),
            (2, 4, GraphBondOrder::Single),
            (3, 5, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonance_system_count(&molecule, 0);
    }

    #[test]
    fn peripheral_oxygen_with_lone_pairs_is_promoted() {
        let elements = [Element::C, Element::O, Element::O, Element::O];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (0, 2, GraphBondOrder::Single),
            (0, 3, GraphBondOrder::Single),
        ];
        let mut molecule = build_molecule(&elements, &bonds, &[]);
        molecule.atoms[3].lone_pairs = 2;
        let molecule = run_resonance_perception(molecule);

        assert!(
            resonant_atom_ids(&molecule).contains(&3),
            "peripheral oxygen with lone pairs should be promoted"
        );
    }

    #[test]
    fn peripheral_atom_without_lone_pairs_is_not_promoted() {
        let elements = [Element::C, Element::O, Element::O, Element::N];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (0, 2, GraphBondOrder::Single),
            (0, 3, GraphBondOrder::Single),
        ];
        let mut molecule = build_molecule(&elements, &bonds, &[]);
        molecule.atoms[3].lone_pairs = 0;
        let molecule = run_resonance_perception(molecule);

        assert!(
            !resonant_atom_ids(&molecule).contains(&3),
            "nitrogen without lone pairs should not be promoted"
        );
    }

    #[test]
    fn carbon_is_not_promoted_to_resonant() {
        let elements = [Element::C, Element::O, Element::O, Element::C];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (0, 2, GraphBondOrder::Single),
            (0, 3, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert!(
            !resonant_atom_ids(&molecule).contains(&3),
            "carbon should never be promoted via peripheral propagation"
        );
    }

    #[test]
    fn peripheral_sulfur_with_lone_pairs_is_promoted() {
        let elements = [Element::C, Element::O, Element::O, Element::S];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (0, 2, GraphBondOrder::Single),
            (0, 3, GraphBondOrder::Single),
        ];
        let mut molecule = build_molecule(&elements, &bonds, &[]);
        molecule.atoms[3].lone_pairs = 2;
        let molecule = run_resonance_perception(molecule);

        assert!(
            resonant_atom_ids(&molecule).contains(&3),
            "sulfur with lone pairs should be promoted"
        );
    }

    #[test]
    fn alkane_has_no_resonance_systems() {
        let elements = [
            Element::C,
            Element::C,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
            Element::H,
        ];
        let bonds = [
            (0, 1, GraphBondOrder::Single),
            (0, 2, GraphBondOrder::Single),
            (0, 3, GraphBondOrder::Single),
            (0, 4, GraphBondOrder::Single),
            (1, 5, GraphBondOrder::Single),
            (1, 6, GraphBondOrder::Single),
            (1, 7, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonance_system_count(&molecule, 0);
        assert!(resonant_atom_ids(&molecule).is_empty());
    }

    #[test]
    fn ketone_is_not_a_resonance_system() {
        let elements = [Element::C, Element::C, Element::C, Element::O];
        let bonds = [
            (0, 1, GraphBondOrder::Single),
            (1, 2, GraphBondOrder::Single),
            (1, 3, GraphBondOrder::Double),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonance_system_count(&molecule, 0);
    }

    #[test]
    fn multiple_functional_groups_create_separate_systems() {
        let elements = [
            Element::C,
            Element::C,
            Element::O,
            Element::O,
            Element::C,
            Element::O,
            Element::O,
        ];
        let bonds = [
            (0, 1, GraphBondOrder::Single),
            (1, 2, GraphBondOrder::Double),
            (1, 3, GraphBondOrder::Single),
            (0, 4, GraphBondOrder::Single),
            (4, 5, GraphBondOrder::Double),
            (4, 6, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonance_system_count(&molecule, 2);
        assert_resonant_atoms(&molecule, &[1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn amide_skips_already_resonant_atoms() {
        let elements = [Element::C, Element::O, Element::N];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (0, 2, GraphBondOrder::Single),
        ];
        let mut molecule = build_molecule(&elements, &bonds, &[]);
        molecule.atoms[0].is_resonant = true;
        let molecule = run_resonance_perception(molecule);

        assert_resonance_system_count(&molecule, 0);
    }

    #[test]
    fn resonance_system_ids_are_sorted() {
        let elements = [Element::C, Element::O, Element::O];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (0, 2, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        assert_resonance_system_count(&molecule, 1);
        let sys = &molecule.resonance_systems[0];

        let mut sorted_atoms = sys.atom_ids.clone();
        sorted_atoms.sort_unstable();
        assert_eq!(sys.atom_ids, sorted_atoms, "atom_ids should be sorted");

        let mut sorted_bonds = sys.bond_ids.clone();
        sorted_bonds.sort_unstable();
        assert_eq!(sys.bond_ids, sorted_bonds, "bond_ids should be sorted");
    }

    #[test]
    fn resonance_system_bond_ids_are_valid() {
        let elements = [Element::N, Element::O, Element::O];
        let bonds = [
            (0, 1, GraphBondOrder::Double),
            (0, 2, GraphBondOrder::Single),
        ];
        let molecule = run_resonance_perception(build_molecule(&elements, &bonds, &[]));

        for sys in &molecule.resonance_systems {
            for &bond_id in &sys.bond_ids {
                assert!(
                    bond_id < molecule.bonds.len(),
                    "bond_id {} exceeds bond count {}",
                    bond_id,
                    molecule.bonds.len()
                );
            }
        }
    }
}
