//! Evaluates fused ring systems to determine whether the atoms are aromatic, anti-aromatic, or neither.
//!
//! The module groups rings, builds localized models that count π-electrons under planarity
//! assumptions, and sets per-atom flags. Importantly, it also registers aromatic rings
//! as `ResonanceSystem`s so that their bonds are treated as resonant in the final topology.

use super::model::{AnnotatedAtom, AnnotatedMolecule, ResonanceSystem, Ring};
use crate::forcefield::dreiding::typer::core::error::PerceptionError;
use crate::forcefield::dreiding::typer::core::properties::GraphBondOrder;
use std::collections::{HashMap, HashSet};

/// Runs aromaticity perception over all ring systems present in the molecule.
///
/// The procedure clusters rings that share atoms, evaluates each cluster as a whole, falls back to
/// ring-by-ring evaluation when mixed behavior occurs, and annotates atoms as aromatic or
/// anti-aromatic accordingly. Confirmed aromatic systems are added to the molecule's
/// resonance systems list.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule whose atom flags should be updated.
///
/// # Returns
///
/// `Ok(())` once every ring system has been processed or when no rings exist.
pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    if molecule.rings.is_empty() {
        return Ok(());
    }

    let ring_systems_indices = find_ring_systems(&molecule.rings);

    for system_indices in ring_systems_indices {
        let system_atoms: HashSet<usize> = system_indices
            .iter()
            .flat_map(|&i| molecule.rings[i].iter())
            .copied()
            .collect();

        let model = AromaticityModel::new(molecule, &system_atoms);

        if model.is_aromatic() {
            apply_aromaticity(molecule, &system_atoms);
        } else if model.is_anti_aromatic() {
            for &atom_id in &system_atoms {
                molecule.atoms[atom_id].is_anti_aromatic = true;
            }
        } else {
            evaluate_rings_individually(molecule, &system_indices);
        }
    }

    Ok(())
}

/// Marks all atoms and bonds in the system as aromatic and resonant.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule to mutate.
/// * `system_atoms` - Atom IDs representing the aromatic system.
fn apply_aromaticity(molecule: &mut AnnotatedMolecule, system_atoms: &HashSet<usize>) {
    for &atom_id in system_atoms {
        let atom = &mut molecule.atoms[atom_id];
        atom.is_aromatic = true;
        atom.is_resonant = true;
    }

    let mut bond_ids = Vec::new();
    for bond in &molecule.bonds {
        if system_atoms.contains(&bond.atom_ids.0) && system_atoms.contains(&bond.atom_ids.1) {
            bond_ids.push(bond.id);
        }
    }

    let mut atom_ids: Vec<usize> = system_atoms.iter().copied().collect();
    atom_ids.sort_unstable();
    bond_ids.sort_unstable();

    molecule
        .resonance_systems
        .push(ResonanceSystem { atom_ids, bond_ids });
}

/// Evaluates each ring independently when a fused system lacks uniform behavior.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule to mutate.
/// * `system_indices` - Indices of rings belonging to the fused system.
fn evaluate_rings_individually(molecule: &mut AnnotatedMolecule, system_indices: &[usize]) {
    for &ring_idx in system_indices {
        let ring_atoms: HashSet<_> = molecule.rings[ring_idx].iter().copied().collect();
        let ring_model = AromaticityModel::new(molecule, &ring_atoms);

        if ring_model.is_aromatic() {
            apply_aromaticity(molecule, &ring_atoms);
        } else if ring_model.is_anti_aromatic() {
            for &atom_id in &ring_atoms {
                molecule.atoms[atom_id].is_anti_aromatic = true;
            }
        }
    }
}

/// Local model capturing the atoms and π-electron count for a ring system.
struct AromaticityModel<'a> {
    /// Annotated molecule providing adjacency information.
    molecule: &'a AnnotatedMolecule,
    /// Atom IDs forming the current system under evaluation.
    atoms: HashSet<usize>,
    /// Computed π-electron count, if evaluation succeeded.
    pi_electrons: Option<u32>,
    /// Flag describing whether the atoms satisfy the planarity heuristic.
    is_potentially_planar: bool,
}

impl<'a> AromaticityModel<'a> {
    /// Constructs the model and immediately evaluates planarity and π-electrons.
    ///
    /// # Arguments
    ///
    /// * `molecule` - Annotated molecule backing the model.
    /// * `system_atoms` - Atom IDs representing a ring system.
    fn new(molecule: &'a AnnotatedMolecule, system_atoms: &HashSet<usize>) -> Self {
        let mut model = Self {
            molecule,
            atoms: system_atoms.clone(),
            pi_electrons: None,
            is_potentially_planar: false,
        };
        model.evaluate();
        model
    }

    /// Returns `true` when the Huckel 4n+2 rule is satisfied.
    fn is_aromatic(&self) -> bool {
        if !self.is_potentially_planar {
            return false;
        }

        let all_from_aromatic_input = self
            .atoms
            .iter()
            .all(|&id| self.molecule.atoms[id].has_aromatic_edge);
        if all_from_aromatic_input {
            return true;
        }

        matches!(self.pi_electrons, Some(pi) if pi > 0 && (pi - 2) % 4 == 0)
    }

    /// Returns `true` when the Huckel 4n rule indicates anti-aromaticity.
    fn is_anti_aromatic(&self) -> bool {
        if !self.is_potentially_planar || self.has_cross_conjugation() {
            return false;
        }
        matches!(self.pi_electrons, Some(pi) if pi > 0 && pi % 4 == 0)
    }

    /// Computes planarity and π-electron counts, caching the results on the struct.
    fn evaluate(&mut self) {
        if !self
            .atoms
            .iter()
            .all(|&id| is_potentially_planar(&self.molecule.atoms[id]))
        {
            self.is_potentially_planar = false;
            return;
        }
        self.is_potentially_planar = true;

        let mut pi_count = 0;
        for &atom_id in &self.atoms {
            if let Some(contribution) = self.count_pi_contribution(atom_id) {
                pi_count += contribution;
            } else {
                self.is_potentially_planar = false;
                self.pi_electrons = None;
                return;
            }
        }
        self.pi_electrons = Some(pi_count);
    }

    /// Checks if an atom participates in a double bond within the ring system.
    fn atom_has_endocyclic_double(&self, atom_id: usize) -> bool {
        self.molecule.adjacency[atom_id]
            .iter()
            .any(|&(n_id, order)| order == GraphBondOrder::Double && self.atoms.contains(&n_id))
    }

    /// Checks if an atom carries a double bond outside the ring system.
    fn atom_has_exocyclic_double(&self, atom_id: usize) -> bool {
        self.molecule.adjacency[atom_id]
            .iter()
            .any(|&(n_id, order)| order == GraphBondOrder::Double && !self.atoms.contains(&n_id))
    }

    /// Detects cross-conjugation, which prevents anti-aromatic classification.
    fn has_cross_conjugation(&self) -> bool {
        self.atoms
            .iter()
            .any(|&atom_id| self.atom_has_exocyclic_double(atom_id))
    }

    /// Computes each atom's π contribution using bond, lone-pair, and resonance flags.
    fn count_pi_contribution(&self, atom_id: usize) -> Option<u32> {
        let atom = &self.molecule.atoms[atom_id];
        let has_endocyclic_double_bond = self.atom_has_endocyclic_double(atom_id);
        let has_exocyclic_double_bond = self.atom_has_exocyclic_double(atom_id);

        if has_endocyclic_double_bond {
            return Some(1);
        }

        if !has_exocyclic_double_bond && atom.lone_pairs > 0 {
            return Some(2);
        }

        if atom.formal_charge == -1 {
            return Some(2);
        }
        if atom.formal_charge == 1 {
            return Some(0);
        }

        if has_exocyclic_double_bond {
            return Some(1);
        }

        if atom.is_resonant && atom.is_in_ring {
            return Some(1);
        }

        if atom.has_aromatic_edge && atom.is_in_ring {
            return Some(1);
        }

        None
    }
}

/// Heuristic planarity test derived from steric number rules.
fn is_potentially_planar(atom: &AnnotatedAtom) -> bool {
    let steric_number = atom.degree + atom.lone_pairs;
    match steric_number {
        0..=3 => true,
        4 => atom.lone_pairs > 0,
        _ => false,
    }
}

/// Groups rings into fused systems via shared atoms.
fn find_ring_systems(rings: &[Ring]) -> Vec<Vec<usize>> {
    if rings.is_empty() {
        return vec![];
    }

    let ring_adj = build_ring_adjacency(rings);
    let mut systems = Vec::new();
    let mut visited = vec![false; rings.len()];

    for i in 0..rings.len() {
        if !visited[i] {
            let mut current_system_indices = Vec::new();
            let mut stack = vec![i];
            visited[i] = true;

            while let Some(ring_idx) = stack.pop() {
                current_system_indices.push(ring_idx);
                for &neighbor_idx in &ring_adj[ring_idx] {
                    if !visited[neighbor_idx] {
                        visited[neighbor_idx] = true;
                        stack.push(neighbor_idx);
                    }
                }
            }
            systems.push(current_system_indices);
        }
    }
    systems
}

/// Builds an adjacency list between rings that share at least one atom.
fn build_ring_adjacency(rings: &[Ring]) -> Vec<Vec<usize>> {
    let mut atom_to_rings: HashMap<usize, Vec<usize>> = HashMap::new();
    for (ring_idx, ring) in rings.iter().enumerate() {
        for &atom_id in ring {
            atom_to_rings.entry(atom_id).or_default().push(ring_idx);
        }
    }

    let mut adj = vec![vec![]; rings.len()];
    for ring_indices in atom_to_rings.values() {
        if ring_indices.len() > 1 {
            for i in 0..ring_indices.len() {
                for j in (i + 1)..ring_indices.len() {
                    let r1 = ring_indices[i];
                    let r2 = ring_indices[j];
                    adj[r1].push(r2);
                    adj[r2].push(r1);
                }
            }
        }
    }

    for neighbors in adj.iter_mut() {
        neighbors.sort_unstable();
        neighbors.dedup();
    }

    adj
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::dreiding::typer::core::graph::MolecularGraph;
    use crate::forcefield::dreiding::typer::core::properties::Element;

    #[derive(Clone, Copy)]
    struct AtomSpec {
        element: Element,
        formal_charge: i8,
        lone_pairs: u8,
        degree_override: Option<u8>,
    }

    impl AtomSpec {
        fn new(element: Element) -> Self {
            Self {
                element,
                formal_charge: 0,
                lone_pairs: 0,
                degree_override: None,
            }
        }

        fn with_charge(mut self, charge: i8) -> Self {
            self.formal_charge = charge;
            self
        }

        fn with_lone_pairs(mut self, lone_pairs: u8) -> Self {
            self.lone_pairs = lone_pairs;
            self
        }

        fn with_degree(mut self, degree: u8) -> Self {
            self.degree_override = Some(degree);
            self
        }
    }

    fn build_test_molecule(
        atom_specs: &[AtomSpec],
        bonds: &[(usize, usize, GraphBondOrder)],
        rings: &[&[usize]],
    ) -> AnnotatedMolecule {
        let mut graph = MolecularGraph::new();
        for spec in atom_specs {
            graph.add_atom(spec.element);
        }
        for &(a, b, order) in bonds {
            graph.add_bond(a, b, order).expect("valid bond definition");
        }

        let mut molecule = AnnotatedMolecule::new(&graph).expect("graph must be valid");
        molecule.rings = rings.iter().map(|ring| ring.to_vec()).collect();
        annotate_ring_flags(&mut molecule);
        apply_atom_specs(&mut molecule, atom_specs);
        molecule
    }

    fn annotate_ring_flags(molecule: &mut AnnotatedMolecule) {
        for ring in &molecule.rings {
            for &atom_id in ring {
                let atom = &mut molecule.atoms[atom_id];
                atom.is_in_ring = true;
            }
        }
    }

    fn apply_atom_specs(molecule: &mut AnnotatedMolecule, specs: &[AtomSpec]) {
        for (i, spec) in specs.iter().enumerate() {
            let atom = &mut molecule.atoms[i];
            atom.formal_charge = spec.formal_charge;
            atom.lone_pairs = spec.lone_pairs;
            if let Some(override_degree) = spec.degree_override {
                atom.degree = override_degree;
            }
            atom.steric_number = atom.degree + atom.lone_pairs;
        }
    }

    fn perceive_aromaticity(mut molecule: AnnotatedMolecule) -> AnnotatedMolecule {
        perceive(&mut molecule).expect("aromaticity perception should succeed");
        molecule
    }

    fn assert_flag_sets(
        molecule: &AnnotatedMolecule,
        expected_aromatic: &[usize],
        expected_anti: &[usize],
    ) {
        use std::collections::HashSet;
        let aromatic: HashSet<_> = molecule
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(idx, atom)| atom.is_aromatic.then_some(idx))
            .collect();
        let anti: HashSet<_> = molecule
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(idx, atom)| atom.is_anti_aromatic.then_some(idx))
            .collect();

        assert_eq!(
            aromatic,
            expected_aromatic.iter().copied().collect(),
            "unexpected aromatic atom assignment"
        );
        assert_eq!(
            anti,
            expected_anti.iter().copied().collect(),
            "unexpected anti-aromatic atom assignment"
        );

        for (idx, atom) in molecule.atoms.iter().enumerate() {
            assert!(
                !(atom.is_aromatic && atom.is_anti_aromatic),
                "atom {idx} cannot be aromatic and anti-aromatic simultaneously"
            );
        }
    }

    fn h() -> AtomSpec {
        AtomSpec::new(Element::H)
    }
    fn c() -> AtomSpec {
        AtomSpec::new(Element::C)
    }
    fn n_pyridine() -> AtomSpec {
        AtomSpec::new(Element::N).with_lone_pairs(1)
    }
    fn n_pyrrole() -> AtomSpec {
        AtomSpec::new(Element::N).with_lone_pairs(1)
    }
    fn o_furan() -> AtomSpec {
        AtomSpec::new(Element::O).with_lone_pairs(2)
    }
    fn o_carbonyl() -> AtomSpec {
        AtomSpec::new(Element::O).with_lone_pairs(2)
    }
    fn b_anion() -> AtomSpec {
        AtomSpec::new(Element::B).with_charge(-1)
    }

    fn benzene() -> AnnotatedMolecule {
        let atoms = vec![c(), c(), c(), c(), c(), c(), h(), h(), h(), h(), h(), h()];
        let bonds = vec![
            (0, 1, GraphBondOrder::Double),
            (1, 2, GraphBondOrder::Single),
            (2, 3, GraphBondOrder::Double),
            (3, 4, GraphBondOrder::Single),
            (4, 5, GraphBondOrder::Double),
            (5, 0, GraphBondOrder::Single),
            (0, 6, GraphBondOrder::Single),
            (1, 7, GraphBondOrder::Single),
            (2, 8, GraphBondOrder::Single),
            (3, 9, GraphBondOrder::Single),
            (4, 10, GraphBondOrder::Single),
            (5, 11, GraphBondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5]])
    }

    fn pyrrole() -> AnnotatedMolecule {
        let atoms = vec![n_pyrrole(), c(), c(), c(), c(), h(), h(), h(), h(), h()];
        let bonds = vec![
            (0, 1, GraphBondOrder::Single),
            (1, 2, GraphBondOrder::Double),
            (2, 3, GraphBondOrder::Single),
            (3, 4, GraphBondOrder::Double),
            (4, 0, GraphBondOrder::Single),
            (0, 5, GraphBondOrder::Single),
            (1, 6, GraphBondOrder::Single),
            (2, 7, GraphBondOrder::Single),
            (3, 8, GraphBondOrder::Single),
            (4, 9, GraphBondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4]])
    }

    fn borabenzene_anion() -> AnnotatedMolecule {
        let atoms = vec![b_anion(), c(), c(), c(), c(), c(), h(), h(), h(), h(), h()];
        let bonds = vec![
            (0, 1, GraphBondOrder::Double),
            (1, 2, GraphBondOrder::Single),
            (2, 3, GraphBondOrder::Double),
            (3, 4, GraphBondOrder::Single),
            (4, 5, GraphBondOrder::Double),
            (5, 0, GraphBondOrder::Single),
            (1, 6, GraphBondOrder::Single),
            (2, 7, GraphBondOrder::Single),
            (3, 8, GraphBondOrder::Single),
            (4, 9, GraphBondOrder::Single),
            (5, 10, GraphBondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5]])
    }

    fn cyclobutadiene() -> AnnotatedMolecule {
        let atoms = vec![c(), c(), c(), c(), h(), h(), h(), h()];
        let bonds = vec![
            (0, 1, GraphBondOrder::Double),
            (1, 2, GraphBondOrder::Single),
            (2, 3, GraphBondOrder::Double),
            (3, 0, GraphBondOrder::Single),
            (0, 4, GraphBondOrder::Single),
            (1, 5, GraphBondOrder::Single),
            (2, 6, GraphBondOrder::Single),
            (3, 7, GraphBondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3]])
    }

    fn cyclooctatetraene_nonplanar() -> AnnotatedMolecule {
        let atoms: Vec<_> = (0..8)
            .map(|idx| {
                if idx % 2 == 0 {
                    c().with_degree(4)
                } else {
                    c()
                }
            })
            .chain((8..16).map(|_| h()))
            .collect();
        let bonds = vec![
            (0, 1, GraphBondOrder::Double),
            (1, 2, GraphBondOrder::Single),
            (2, 3, GraphBondOrder::Double),
            (3, 4, GraphBondOrder::Single),
            (4, 5, GraphBondOrder::Double),
            (5, 6, GraphBondOrder::Single),
            (6, 7, GraphBondOrder::Double),
            (7, 0, GraphBondOrder::Single),
            (0, 8, GraphBondOrder::Single),
            (1, 9, GraphBondOrder::Single),
            (2, 10, GraphBondOrder::Single),
            (3, 11, GraphBondOrder::Single),
            (4, 12, GraphBondOrder::Single),
            (5, 13, GraphBondOrder::Single),
            (6, 14, GraphBondOrder::Single),
            (7, 15, GraphBondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5, 6, 7]])
    }

    fn cyclohexane() -> AnnotatedMolecule {
        let atoms = (0..18)
            .map(|i| if i < 6 { c() } else { h() })
            .collect::<Vec<_>>();
        let bonds = vec![
            (0, 1, GraphBondOrder::Single),
            (1, 2, GraphBondOrder::Single),
            (2, 3, GraphBondOrder::Single),
            (3, 4, GraphBondOrder::Single),
            (4, 5, GraphBondOrder::Single),
            (5, 0, GraphBondOrder::Single),
            (0, 6, GraphBondOrder::Single),
            (0, 7, GraphBondOrder::Single),
            (1, 8, GraphBondOrder::Single),
            (1, 9, GraphBondOrder::Single),
            (2, 10, GraphBondOrder::Single),
            (2, 11, GraphBondOrder::Single),
            (3, 12, GraphBondOrder::Single),
            (3, 13, GraphBondOrder::Single),
            (4, 14, GraphBondOrder::Single),
            (4, 15, GraphBondOrder::Single),
            (5, 16, GraphBondOrder::Single),
            (5, 17, GraphBondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5]])
    }

    fn naphthalene() -> AnnotatedMolecule {
        let atoms = (0..18)
            .map(|i| if i < 10 { c() } else { h() })
            .collect::<Vec<_>>();
        let bonds = vec![
            (0, 1, GraphBondOrder::Double),
            (1, 2, GraphBondOrder::Single),
            (2, 3, GraphBondOrder::Double),
            (3, 4, GraphBondOrder::Single),
            (4, 9, GraphBondOrder::Single),
            (9, 8, GraphBondOrder::Double),
            (8, 7, GraphBondOrder::Single),
            (7, 6, GraphBondOrder::Double),
            (6, 5, GraphBondOrder::Single),
            (5, 0, GraphBondOrder::Single),
            (4, 5, GraphBondOrder::Double),
            (0, 10, GraphBondOrder::Single),
            (1, 11, GraphBondOrder::Single),
            (2, 12, GraphBondOrder::Single),
            (3, 13, GraphBondOrder::Single),
            (6, 14, GraphBondOrder::Single),
            (7, 15, GraphBondOrder::Single),
            (8, 16, GraphBondOrder::Single),
            (9, 17, GraphBondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5], &[4, 5, 6, 7, 8, 9]])
    }

    fn anthracene() -> AnnotatedMolecule {
        let atoms = (0..24)
            .map(|i| if i < 14 { c() } else { h() })
            .collect::<Vec<_>>();
        let bonds = vec![
            (0, 1, GraphBondOrder::Double),
            (1, 2, GraphBondOrder::Single),
            (2, 3, GraphBondOrder::Double),
            (3, 4, GraphBondOrder::Single),
            (4, 13, GraphBondOrder::Single),
            (13, 12, GraphBondOrder::Double),
            (12, 11, GraphBondOrder::Single),
            (11, 10, GraphBondOrder::Double),
            (10, 5, GraphBondOrder::Single),
            (5, 0, GraphBondOrder::Single),
            (4, 5, GraphBondOrder::Double),
            (10, 9, GraphBondOrder::Single),
            (9, 8, GraphBondOrder::Double),
            (8, 7, GraphBondOrder::Single),
            (7, 6, GraphBondOrder::Double),
            (6, 11, GraphBondOrder::Single),
            (0, 14, GraphBondOrder::Single),
            (1, 15, GraphBondOrder::Single),
            (2, 16, GraphBondOrder::Single),
            (3, 17, GraphBondOrder::Single),
            (6, 18, GraphBondOrder::Single),
            (7, 19, GraphBondOrder::Single),
            (8, 20, GraphBondOrder::Single),
            (9, 21, GraphBondOrder::Single),
            (12, 22, GraphBondOrder::Single),
            (13, 23, GraphBondOrder::Single),
        ];
        build_test_molecule(
            &atoms,
            &bonds,
            &[
                &[0, 1, 2, 3, 4, 5],
                &[4, 5, 10, 11, 6, 7, 8, 9, 13, 12],
                &[6, 7, 8, 9, 10, 11],
            ],
        )
    }

    fn purine() -> AnnotatedMolecule {
        let atoms = vec![
            n_pyridine(),
            c(),
            n_pyridine(),
            c(),
            c(),
            n_pyrrole(),
            c(),
            n_pyridine(),
            c(),
            h(),
            h(),
            h(),
            h(),
        ];
        let bonds = vec![
            (0, 1, GraphBondOrder::Double),
            (1, 2, GraphBondOrder::Single),
            (2, 3, GraphBondOrder::Double),
            (3, 4, GraphBondOrder::Single),
            (4, 5, GraphBondOrder::Single),
            (5, 0, GraphBondOrder::Single),
            (3, 8, GraphBondOrder::Single),
            (8, 7, GraphBondOrder::Double),
            (7, 6, GraphBondOrder::Single),
            (6, 4, GraphBondOrder::Double),
            (1, 9, GraphBondOrder::Single),
            (5, 10, GraphBondOrder::Single),
            (6, 11, GraphBondOrder::Single),
            (8, 12, GraphBondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5], &[3, 4, 6, 7, 8]])
    }

    fn alpha_pyrone() -> AnnotatedMolecule {
        let atoms = vec![
            c(),
            o_furan(),
            c(),
            c(),
            c(),
            c(),
            o_carbonyl(),
            h(),
            h(),
            h(),
            h(),
        ];
        let bonds = vec![
            (0, 1, GraphBondOrder::Single),
            (1, 2, GraphBondOrder::Single),
            (2, 3, GraphBondOrder::Double),
            (3, 4, GraphBondOrder::Single),
            (4, 5, GraphBondOrder::Double),
            (5, 0, GraphBondOrder::Single),
            (0, 6, GraphBondOrder::Double),
            (2, 7, GraphBondOrder::Single),
            (3, 8, GraphBondOrder::Single),
            (4, 9, GraphBondOrder::Single),
            (5, 10, GraphBondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4, 5]])
    }

    fn pyrazole() -> AnnotatedMolecule {
        let atoms = vec![c(), n_pyrrole(), n_pyridine(), c(), c(), h(), h(), h(), h()];
        let bonds = vec![
            (0, 1, GraphBondOrder::Single),
            (1, 2, GraphBondOrder::Single),
            (2, 3, GraphBondOrder::Double),
            (3, 4, GraphBondOrder::Single),
            (4, 0, GraphBondOrder::Double),
            (0, 5, GraphBondOrder::Single),
            (1, 6, GraphBondOrder::Single),
            (3, 7, GraphBondOrder::Single),
            (4, 8, GraphBondOrder::Single),
        ];
        build_test_molecule(&atoms, &bonds, &[&[0, 1, 2, 3, 4]])
    }

    #[test]
    fn benzene_ring_is_aromatic() {
        let molecule = perceive_aromaticity(benzene());
        assert_flag_sets(&molecule, &[0, 1, 2, 3, 4, 5], &[]);

        assert_eq!(
            molecule.resonance_systems.len(),
            1,
            "Benzene should form 1 resonance system"
        );
        assert_eq!(molecule.resonance_systems[0].atom_ids.len(), 6);
        assert_eq!(molecule.resonance_systems[0].bond_ids.len(), 6);
    }

    #[test]
    fn pyrrole_lone_pair_contributes_to_aromaticity() {
        let molecule = perceive_aromaticity(pyrrole());
        assert_flag_sets(&molecule, &[0, 1, 2, 3, 4], &[]);
    }

    #[test]
    fn borabenzene_anion_is_aromatic() {
        let molecule = perceive_aromaticity(borabenzene_anion());
        assert_flag_sets(&molecule, &[0, 1, 2, 3, 4, 5], &[]);
    }

    #[test]
    fn cyclobutadiene_detected_as_antiaromatic() {
        let molecule = perceive_aromaticity(cyclobutadiene());
        assert_flag_sets(&molecule, &[], &[0, 1, 2, 3]);
    }

    #[test]
    fn cyclooctatetraene_rejected_due_to_non_planarity() {
        let molecule = perceive_aromaticity(cyclooctatetraene_nonplanar());
        assert_flag_sets(&molecule, &[], &[]);
    }

    #[test]
    fn cyclohexane_is_non_aromatic() {
        let molecule = perceive_aromaticity(cyclohexane());
        assert_flag_sets(&molecule, &[], &[]);
    }

    #[test]
    fn naphthalene_fused_rings_are_aromatic() {
        let molecule = perceive_aromaticity(naphthalene());
        assert_flag_sets(&molecule, &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], &[]);

        assert!(!molecule.resonance_systems.is_empty());
        let sys = &molecule.resonance_systems[0];
        assert_eq!(sys.atom_ids.len(), 10);
        assert_eq!(sys.bond_ids.len(), 11);
    }

    #[test]
    fn anthracene_three_ring_system_is_aromatic() {
        let molecule = perceive_aromaticity(anthracene());
        assert_flag_sets(&molecule, &(0..14).collect::<Vec<_>>(), &[]);
    }

    #[test]
    fn purine_dual_ring_system_is_aromatic() {
        let molecule = perceive_aromaticity(purine());
        assert_flag_sets(&molecule, &[0, 1, 2, 3, 4, 5, 6, 7, 8], &[]);
    }

    #[test]
    fn alpha_pyrone_is_non_aromatic() {
        let molecule = perceive_aromaticity(alpha_pyrone());
        assert_flag_sets(&molecule, &[], &[]);
    }

    #[test]
    fn pyrazole_is_aromatic() {
        let molecule = perceive_aromaticity(pyrazole());
        assert_flag_sets(&molecule, &[0, 1, 2, 3, 4], &[]);
    }

    #[test]
    fn biphenyl_registers_two_separate_resonance_systems() {
        let atoms = vec![
            c(),
            c(),
            c(),
            c(),
            c(),
            c(),
            c(),
            c(),
            c(),
            c(),
            c(),
            c(),
            h(),
            h(),
            h(),
            h(),
            h(),
            h(),
            h(),
            h(),
            h(),
            h(),
        ];
        let bonds = vec![
            (0, 1, GraphBondOrder::Double),
            (1, 2, GraphBondOrder::Single),
            (2, 3, GraphBondOrder::Double),
            (3, 4, GraphBondOrder::Single),
            (4, 5, GraphBondOrder::Double),
            (5, 0, GraphBondOrder::Single),
            (5, 6, GraphBondOrder::Single),
            (6, 7, GraphBondOrder::Double),
            (7, 8, GraphBondOrder::Single),
            (8, 9, GraphBondOrder::Double),
            (9, 10, GraphBondOrder::Single),
            (10, 11, GraphBondOrder::Double),
            (11, 6, GraphBondOrder::Single),
        ];

        let r1 = vec![0, 1, 2, 3, 4, 5];
        let r2 = vec![6, 7, 8, 9, 10, 11];

        let molecule = build_test_molecule(&atoms, &bonds, &[&r1, &r2]);
        let molecule = perceive_aromaticity(molecule);

        assert_eq!(
            molecule.resonance_systems.len(),
            2,
            "Should have 2 independent resonance systems"
        );

        let mut resonant_bond_ids = HashSet::new();
        for sys in &molecule.resonance_systems {
            for &bid in &sys.bond_ids {
                resonant_bond_ids.insert(bid);
            }
        }

        let connecting_bond_id = molecule
            .bonds
            .iter()
            .find(|b| {
                (b.atom_ids.0 == 5 && b.atom_ids.1 == 6) || (b.atom_ids.0 == 6 && b.atom_ids.1 == 5)
            })
            .map(|b| b.id)
            .expect("Connecting bond must exist");

        assert!(
            !resonant_bond_ids.contains(&connecting_bond_id),
            "Connecting bond should NOT be resonant"
        );
    }
}
