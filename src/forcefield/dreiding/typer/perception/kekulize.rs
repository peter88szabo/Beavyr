//! Resolves aromatic bonds into concrete single/double assignments via a Kekulé solver.
//!
//! The logic here isolates aromatic systems, validates that they sit inside rings, and runs a
//! backtracking search that respects valence limits and heteroatom allowances before updating the
//! annotated molecule in-place.

use super::model::AnnotatedMolecule;
use crate::forcefield::dreiding::typer::core::error::PerceptionError;
use crate::forcefield::dreiding::typer::core::properties::{Element, GraphBondOrder};
use std::collections::{HashMap, VecDeque, hash_map::Entry};

/// Converts aromatic bonds inside the molecule to alternating single/double assignments.
///
/// The pass validates that every aromatic bond belongs to a ring, partitions the bonds into
/// connected systems, and then runs a Kekulé solver per system.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule whose bond orders and adjacency lists are mutated in place.
///
/// # Returns
///
/// `Ok(())` when every aromatic system receives a valid assignment.
///
/// # Errors
///
/// Returns [`PerceptionError::KekulizationFailed`] when an aromatic bond lies outside a ring or no
/// valid alternating assignment exists for a system.
pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    let mut aromatic_bonds = Vec::new();
    let mut aromatic_atom_flags = vec![false; molecule.atoms.len()];
    for bond in &molecule.bonds {
        if bond.order == GraphBondOrder::Aromatic {
            aromatic_atom_flags[bond.atom_ids.0] = true;
            aromatic_atom_flags[bond.atom_ids.1] = true;
            aromatic_bonds.push(bond.id);
        }
    }

    for (atom, flag) in molecule
        .atoms
        .iter_mut()
        .zip(aromatic_atom_flags.into_iter())
    {
        if flag {
            atom.has_aromatic_edge = true;
        }
    }

    if aromatic_bonds.is_empty() {
        return Ok(());
    }

    validate_aromatic_bonds_in_rings(molecule, &aromatic_bonds)?;

    let aromatic_systems = find_aromatic_systems(molecule, &aromatic_bonds);

    let mut new_bond_orders = HashMap::new();

    for system_bonds in aromatic_systems {
        let mut solver = KekuleSolver::new(molecule, &system_bonds);
        match solver.solve() {
            Some(solution) => {
                new_bond_orders.extend(solution);
            }
            None => {
                return Err(PerceptionError::KekulizationFailed {
                    message: "could not find a valid Kekulé structure for an aromatic system"
                        .to_string(),
                });
            }
        }
    }

    for (bond_id, new_order) in new_bond_orders {
        let bond = molecule.bonds.iter_mut().find(|b| b.id == bond_id).unwrap();
        let (u, v) = bond.atom_ids;
        bond.order = new_order;

        for (neighbor_id, order) in molecule.adjacency[u].iter_mut() {
            if *neighbor_id == v {
                *order = new_order;
                break;
            }
        }
        for (neighbor_id, order) in molecule.adjacency[v].iter_mut() {
            if *neighbor_id == u {
                *order = new_order;
                break;
            }
        }
    }

    Ok(())
}

/// Backtracking assignment helper that finds valid bond orders for one aromatic system.
struct KekuleSolver<'a> {
    molecule: &'a AnnotatedMolecule,
    bond_indices: Vec<usize>,
    assignments: Vec<Option<GraphBondOrder>>,
}

impl<'a> KekuleSolver<'a> {
    /// Creates a solver scoped to the provided bond identifiers.
    ///
    /// # Arguments
    ///
    /// * `molecule` - Annotated molecule providing bond/atom metadata.
    /// * `system_bond_ids` - Aromatic bond IDs belonging to one connected system.
    fn new(molecule: &'a AnnotatedMolecule, system_bond_ids: &[usize]) -> Self {
        let bond_indices: Vec<usize> = system_bond_ids
            .iter()
            .map(|id| molecule.bonds.iter().position(|b| b.id == *id).unwrap())
            .collect();

        Self {
            molecule,
            bond_indices: bond_indices.clone(),
            assignments: vec![None; bond_indices.len()],
        }
    }

    /// Attempts to assign single/double orders to every bond in the system.
    ///
    /// # Returns
    ///
    /// Map of bond IDs to resolved orders, or `None` if no assignment satisfies the constraints.
    fn solve(&mut self) -> Option<HashMap<usize, GraphBondOrder>> {
        if self.backtrack(0) {
            let solution = self
                .bond_indices
                .iter()
                .enumerate()
                .map(|(i, &bond_idx)| {
                    let bond_id = self.molecule.bonds[bond_idx].id;
                    let order = self.assignments[i].unwrap();
                    (bond_id, order)
                })
                .collect();
            Some(solution)
        } else {
            None
        }
    }

    /// Recursively assigns orders while pruning inconsistent branches.
    ///
    /// # Arguments
    ///
    /// * `k` - Index of the bond currently being assigned.
    ///
    /// # Returns
    ///
    /// `true` when a full assignment is found downstream.
    fn backtrack(&mut self, k: usize) -> bool {
        if k == self.assignments.len() {
            return true;
        }

        for &order_choice in &[GraphBondOrder::Double, GraphBondOrder::Single] {
            self.assignments[k] = Some(order_choice);

            if self.is_consistent(k) && self.backtrack(k + 1) {
                return true;
            }
        }

        self.assignments[k] = None;
        false
    }

    /// Checks whether the partial assignment at index `k` respects valence constraints.
    fn is_consistent(&self, k: usize) -> bool {
        let bond_idx = self.bond_indices[k];
        let (u, v) = self.molecule.bonds[bond_idx].atom_ids;

        self.is_valence_ok(u) && self.is_valence_ok(v)
    }

    /// Verifies that the atom's accumulated valence does not exceed its maximum allowed value.
    ///
    /// Lone allowances let aromatic nitrogens/phosphorus carry at most one double bond while still
    /// being counted as having valence one toward aromatic contributions.
    fn is_valence_ok(&self, atom_id: usize) -> bool {
        let max_valence = get_max_valence(self.molecule.atoms[atom_id].element);
        let mut current_valence = 0;
        let mut aromatic_double_allowance = 0u8;
        let element = self.molecule.atoms[atom_id].element;

        for (neighbor_id, initial_order) in &self.molecule.adjacency[atom_id] {
            let bond = self
                .molecule
                .bonds
                .iter()
                .find(|b| {
                    (b.atom_ids.0 == atom_id && b.atom_ids.1 == *neighbor_id)
                        || (b.atom_ids.0 == *neighbor_id && b.atom_ids.1 == atom_id)
                })
                .unwrap();

            if *initial_order == GraphBondOrder::Aromatic {
                if let Some(assigned_order) = self
                    .bond_indices
                    .iter()
                    .position(|&idx| self.molecule.bonds[idx].id == bond.id)
                    .and_then(|pos| self.assignments[pos])
                {
                    let contribution = if matches!(element, Element::N | Element::P)
                        && assigned_order == GraphBondOrder::Double
                    {
                        if aromatic_double_allowance >= 1 {
                            return false;
                        }
                        aromatic_double_allowance += 1;
                        1
                    } else {
                        bond_order_to_valence(assigned_order)
                    };
                    current_valence += contribution;
                }
            } else {
                current_valence += bond_order_to_valence(*initial_order);
            }
        }

        current_valence <= max_valence
    }
}

/// Groups aromatic bonds into connected systems for independent solving.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule providing adjacency and bond metadata.
/// * `aromatic_bonds` - All bond IDs currently marked as aromatic.
///
/// # Returns
///
/// Vector of bond-ID lists, each representing one connected aromatic system.
fn find_aromatic_systems(
    molecule: &AnnotatedMolecule,
    aromatic_bonds: &[usize],
) -> Vec<Vec<usize>> {
    let mut systems = Vec::new();
    let mut visited_bonds = HashMap::new();

    for &start_bond_id in aromatic_bonds {
        if visited_bonds.contains_key(&start_bond_id) {
            continue;
        }

        let mut current_system = Vec::new();
        let mut queue = VecDeque::new();

        queue.push_back(start_bond_id);
        visited_bonds.insert(start_bond_id, true);

        while let Some(bond_id) = queue.pop_front() {
            current_system.push(bond_id);
            let bond = molecule.bonds.iter().find(|b| b.id == bond_id).unwrap();
            let (u, v) = bond.atom_ids;

            for atom_id in [u, v] {
                for (neighbor_id, order) in &molecule.adjacency[atom_id] {
                    if *order == GraphBondOrder::Aromatic {
                        let neighbor_bond = molecule
                            .bonds
                            .iter()
                            .find(|b| {
                                (b.atom_ids.0 == atom_id && b.atom_ids.1 == *neighbor_id)
                                    || (b.atom_ids.0 == *neighbor_id && b.atom_ids.1 == atom_id)
                            })
                            .unwrap();

                        if let Entry::Vacant(entry) = visited_bonds.entry(neighbor_bond.id) {
                            entry.insert(true);
                            queue.push_back(neighbor_bond.id);
                        }
                    }
                }
            }
        }
        systems.push(current_system);
    }
    systems
}

/// Ensures every aromatic bond resides entirely inside a perceived ring.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule containing ring flags.
/// * `aromatic_bonds` - Bond IDs that must be validated.
///
/// # Errors
///
/// Returns [`PerceptionError::KekulizationFailed`] when any bond touches a non-ring atom.
fn validate_aromatic_bonds_in_rings(
    molecule: &AnnotatedMolecule,
    aromatic_bonds: &[usize],
) -> Result<(), PerceptionError> {
    for &bond_id in aromatic_bonds {
        let bond = molecule.bonds.iter().find(|b| b.id == bond_id).unwrap();
        let (u, v) = bond.atom_ids;
        if !molecule.atoms[u].is_in_ring || !molecule.atoms[v].is_in_ring {
            return Err(PerceptionError::KekulizationFailed {
                message: format!(
                    "aromatic bond (ID {}) found with at least one atom not in a ring",
                    bond_id
                ),
            });
        }
    }
    Ok(())
}

/// Returns the maximum valence allowed for the provided element.
///
/// # Arguments
///
/// * `element` - Element whose typical valence limit should be enforced.
///
/// # Returns
///
/// Maximum valence counted toward aromatic assignments.
fn get_max_valence(element: Element) -> u8 {
    match element {
        Element::H | Element::F | Element::Cl | Element::Br | Element::I => 1,
        Element::O | Element::S => 2,
        Element::N | Element::P => 3,
        Element::C | Element::Si => 4,
        Element::B => 3,
        _ => 8,
    }
}

/// Converts a bond order into its valence contribution.
///
/// # Arguments
///
/// * `order` - Bond order being counted toward valence.
///
/// # Returns
///
/// Integer contribution consistent with typical valence bookkeeping.
fn bond_order_to_valence(order: GraphBondOrder) -> u8 {
    match order {
        GraphBondOrder::Single => 1,
        GraphBondOrder::Double => 2,
        GraphBondOrder::Triple => 3,
        GraphBondOrder::Aromatic => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::dreiding::typer::core::graph::MolecularGraph;
    use crate::forcefield::dreiding::typer::core::properties::Element;

    const BENZENE_ELEMENTS: [Element; 6] = [Element::C; 6];
    const BENZENE_BONDS: [(usize, usize); 6] = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
    const BENZENE_RING: [usize; 6] = [0, 1, 2, 3, 4, 5];
    const BENZENE_H_POSITIONS: [usize; 6] = [0, 1, 2, 3, 4, 5];

    const CYCLOBUTADIENE_ELEMENTS: [Element; 4] = [Element::C; 4];
    const CYCLOBUTADIENE_BONDS: [(usize, usize); 4] = [(0, 1), (1, 2), (2, 3), (3, 0)];
    const CYCLOBUTADIENE_RING: [usize; 4] = [0, 1, 2, 3];
    const CYCLOBUTADIENE_H_POSITIONS: [usize; 4] = [0, 1, 2, 3];

    const PYRIDINE_ELEMENTS: [Element; 6] = [
        Element::N,
        Element::C,
        Element::C,
        Element::C,
        Element::C,
        Element::C,
    ];
    const PYRIDINE_BONDS: [(usize, usize); 6] = BENZENE_BONDS;
    const PYRIDINE_RING: [usize; 6] = BENZENE_RING;
    const PYRIDINE_H_POSITIONS: [usize; 5] = [1, 2, 3, 4, 5];

    const NAPHTHALENE_ELEMENTS: [Element; 10] = [Element::C; 10];
    const NAPHTHALENE_BONDS: [(usize, usize); 11] = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 0),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 4),
    ];
    const NAPHTHALENE_RING_LEFT: [usize; 6] = [0, 1, 2, 3, 4, 5];
    const NAPHTHALENE_RING_RIGHT: [usize; 6] = [4, 5, 6, 7, 8, 9];
    const NAPHTHALENE_H_POSITIONS: [usize; 8] = [0, 1, 2, 3, 6, 7, 8, 9];

    const ANTHRACENE_ELEMENTS: [Element; 14] = [Element::C; 14];
    const ANTHRACENE_BONDS: [(usize, usize); 16] = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 0),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 4),
        (9, 10),
        (10, 11),
        (11, 12),
        (12, 13),
        (13, 8),
    ];
    const ANTHRACENE_RING_LEFT: [usize; 6] = [0, 1, 2, 3, 4, 5];
    const ANTHRACENE_RING_MIDDLE: [usize; 6] = [4, 5, 6, 7, 8, 9];
    const ANTHRACENE_RING_RIGHT: [usize; 6] = [8, 9, 10, 11, 12, 13];
    const ANTHRACENE_H_POSITIONS: [usize; 10] = [0, 1, 2, 3, 6, 7, 10, 11, 12, 13];

    fn aromatic_fixture(
        heavy_elements: &[Element],
        aromatic_bonds: &[(usize, usize)],
        hydrogens_on: &[usize],
        rings: &[&[usize]],
    ) -> AnnotatedMolecule {
        let mut graph = MolecularGraph::new();
        for &element in heavy_elements {
            graph.add_atom(element);
        }
        for &(u, v) in aromatic_bonds {
            graph
                .add_bond(u, v, GraphBondOrder::Aromatic)
                .expect("valid aromatic bond");
        }
        for &atom_id in hydrogens_on {
            let hydrogen_id = graph.add_atom(Element::H);
            graph
                .add_bond(atom_id, hydrogen_id, GraphBondOrder::Single)
                .expect("valid X-H bond");
        }

        let mut molecule = AnnotatedMolecule::new(&graph).expect("graph should be valid");
        annotate_ring_flags(&mut molecule, rings);
        molecule
    }

    fn annotate_ring_flags(molecule: &mut AnnotatedMolecule, rings: &[&[usize]]) {
        for ring in rings {
            let ring_size = ring.len() as u8;
            for &atom_id in *ring {
                let props = &mut molecule.atoms[atom_id];
                props.is_in_ring = true;
                match props.smallest_ring_size {
                    Some(current) if current <= ring_size => {}
                    _ => props.smallest_ring_size = Some(ring_size),
                }
            }
        }
    }

    fn assert_kekule_solution(molecule: &mut AnnotatedMolecule, rings: &[&[usize]]) {
        perceive(molecule).expect("kekulization should succeed");
        assert_no_aromatic_bonds(molecule);
        for &ring in rings {
            assert_alternating_cycle(molecule, ring);
        }
    }

    fn assert_no_aromatic_bonds(molecule: &AnnotatedMolecule) {
        assert!(
            molecule
                .bonds
                .iter()
                .all(|bond| bond.order == GraphBondOrder::Single
                    || bond.order == GraphBondOrder::Double),
            "all aromatic bonds should be resolved into concrete orders"
        );
    }

    fn assert_alternating_cycle(molecule: &AnnotatedMolecule, ring: &[usize]) {
        assert!(ring.len() >= 3, "rings must contain at least three atoms");
        let first_order = bond_order_between(molecule, ring[0], ring[1]);
        assert!(matches!(
            first_order,
            GraphBondOrder::Single | GraphBondOrder::Double
        ));

        let mut prev_order = first_order;
        for i in 1..ring.len() {
            let u = ring[i];
            let v = ring[(i + 1) % ring.len()];
            let order = bond_order_between(molecule, u, v);
            assert!(matches!(
                order,
                GraphBondOrder::Single | GraphBondOrder::Double
            ));
            assert_ne!(
                order, prev_order,
                "bond orders should alternate around the ring"
            );
            prev_order = order;
        }
    }

    fn bond_order_between(molecule: &AnnotatedMolecule, u: usize, v: usize) -> GraphBondOrder {
        molecule
            .bonds
            .iter()
            .find(|bond| {
                (bond.atom_ids.0 == u && bond.atom_ids.1 == v)
                    || (bond.atom_ids.0 == v && bond.atom_ids.1 == u)
            })
            .expect("bond should exist")
            .order
    }

    #[test]
    fn perceive_returns_error_when_aromatic_bond_is_not_in_ring() {
        let rings: [&[usize]; 0] = [];
        let mut molecule = aromatic_fixture(
            &BENZENE_ELEMENTS,
            &BENZENE_BONDS,
            &BENZENE_H_POSITIONS,
            &rings,
        );

        let err = perceive(&mut molecule).expect_err("atoms must be flagged as ring members");
        match err {
            PerceptionError::KekulizationFailed { message } => {
                assert!(message.contains("not in a ring"));
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }

    #[test]
    fn benzene_produces_classic_alternating_pattern() {
        let rings = [&BENZENE_RING[..]];
        let mut molecule = aromatic_fixture(
            &BENZENE_ELEMENTS,
            &BENZENE_BONDS,
            &BENZENE_H_POSITIONS,
            &rings,
        );
        assert_kekule_solution(&mut molecule, &rings);
    }

    #[test]
    fn cyclobutadiene_is_kekulized_even_though_antiaromatic() {
        let rings = [&CYCLOBUTADIENE_RING[..]];
        let mut molecule = aromatic_fixture(
            &CYCLOBUTADIENE_ELEMENTS,
            &CYCLOBUTADIENE_BONDS,
            &CYCLOBUTADIENE_H_POSITIONS,
            &rings,
        );
        assert_kekule_solution(&mut molecule, &rings);
    }

    #[test]
    fn pyridine_handles_heteroatom_with_correct_assignment() {
        let rings = [&PYRIDINE_RING[..]];
        let mut molecule = aromatic_fixture(
            &PYRIDINE_ELEMENTS,
            &PYRIDINE_BONDS,
            &PYRIDINE_H_POSITIONS,
            &rings,
        );
        assert_kekule_solution(&mut molecule, &rings);

        let double_bonds_to_n: usize = molecule
            .bonds
            .iter()
            .filter(|bond| {
                (bond.atom_ids.0 == 0 || bond.atom_ids.1 == 0)
                    && bond.order == GraphBondOrder::Double
            })
            .count();
        assert_eq!(double_bonds_to_n, 1);
    }

    #[test]
    fn naphthalene_kekulization_handles_fused_rings() {
        let rings = [&NAPHTHALENE_RING_LEFT[..], &NAPHTHALENE_RING_RIGHT[..]];
        let mut molecule = aromatic_fixture(
            &NAPHTHALENE_ELEMENTS,
            &NAPHTHALENE_BONDS,
            &NAPHTHALENE_H_POSITIONS,
            &rings,
        );
        assert_kekule_solution(&mut molecule, &rings);
    }

    #[test]
    fn anthracene_kekulization_covers_three_linearly_fused_rings() {
        let rings = [
            &ANTHRACENE_RING_LEFT[..],
            &ANTHRACENE_RING_MIDDLE[..],
            &ANTHRACENE_RING_RIGHT[..],
        ];
        let mut molecule = aromatic_fixture(
            &ANTHRACENE_ELEMENTS,
            &ANTHRACENE_BONDS,
            &ANTHRACENE_H_POSITIONS,
            &rings,
        );
        assert_kekule_solution(&mut molecule, &rings);
    }
}
