//! Detects rings and records small-set cycle representatives for subsequent perception stages.
//!
//! This module builds a minimal cycle basis from the molecular graph so aromaticity, resonance,
//! and hybridization passes can quickly determine ring membership and sizes.

use super::model::{AnnotatedMolecule, NeighborBond, Ring};
use crate::forcefield::dreiding::typer::core::error::PerceptionError;
use crate::forcefield::dreiding::typer::core::properties::GraphBondOrder;
use std::collections::{HashMap, VecDeque};

/// Computes ring information for the supplied annotated molecule.
///
/// Runs connected-component counting, enumerates simple cycle candidates, chooses a minimal cycle
/// basis, and marks atoms with ring membership metadata.
///
/// # Arguments
///
/// * `molecule` - Mutable annotated molecule that will receive ring annotations.
///
/// # Returns
///
/// `Ok(())` whether rings exist or not; the step keeps the `Result` signature to match the broader
/// perception pipeline but currently never emits [`PerceptionError`].
pub fn perceive(molecule: &mut AnnotatedMolecule) -> Result<(), PerceptionError> {
    let num_atoms = molecule.atoms.len();
    if num_atoms == 0 {
        return Ok(());
    }
    let num_bonds = molecule.bonds.len();
    let num_components = count_components(num_atoms, &molecule.adjacency);
    let cyclomatic_number = num_bonds as isize - num_atoms as isize + num_components as isize;

    if cyclomatic_number <= 0 {
        return Ok(());
    }

    let bond_id_to_index: HashMap<usize, usize> = molecule
        .bonds
        .iter()
        .map(|b| b.id)
        .enumerate()
        .map(|(i, id)| (id, i))
        .collect();

    let mut workspace = RingSearchWorkspace::new(num_atoms);
    let candidates = enumerate_cycle_candidates(molecule, &mut workspace);

    let sssr_candidates =
        select_minimal_cycle_basis(candidates, cyclomatic_number as usize, &bond_id_to_index);

    let final_rings: Vec<Ring> = sssr_candidates
        .into_iter()
        .map(|c| {
            let mut atom_ids = c.atom_ids;
            atom_ids.sort_unstable();
            atom_ids
        })
        .collect();
    molecule.rings = final_rings;

    annotate_atoms_with_ring_info(molecule);

    Ok(())
}

/// Reusable scratch buffers to avoid per-bond allocations during ring search.
struct RingSearchWorkspace {
    queue: VecDeque<usize>,
    visited: Vec<bool>,
    parent: Vec<Option<(usize, usize)>>,
}

impl RingSearchWorkspace {
    fn new(num_atoms: usize) -> Self {
        Self {
            queue: VecDeque::with_capacity(num_atoms),
            visited: vec![false; num_atoms],
            parent: vec![None; num_atoms],
        }
    }

    fn reset(&mut self) {
        self.queue.clear();
        self.visited.fill(false);
        self.parent.fill(None);
    }
}

/// Cycle descriptor storing both atom and bond identifiers.
struct RingCandidate {
    /// Ordered atom identifiers along the candidate cycle.
    atom_ids: Vec<usize>,
    /// Ordered bond identifiers along the candidate cycle.
    bond_ids: Vec<usize>,
    /// Cycle length measured in edges.
    len: usize,
}

/// Enumerates simple cycles by removing each bond and searching for alternate paths.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule whose adjacency and bonds will be analyzed.
/// * `workspace` - Reusable BFS buffers to avoid per-bond allocations.
///
/// # Returns
///
/// Collection of candidate rings containing atom and bond identifiers.
fn enumerate_cycle_candidates(
    molecule: &AnnotatedMolecule,
    workspace: &mut RingSearchWorkspace,
) -> Vec<RingCandidate> {
    let mut candidates = Vec::new();

    for bond_to_remove in &molecule.bonds {
        if let Some(path) = shortest_path_bfs(
            molecule,
            bond_to_remove.atom_ids.0,
            bond_to_remove.atom_ids.1,
            Some(bond_to_remove.id),
            workspace,
        ) {
            let mut atom_ids = path.atom_ids;
            let mut bond_ids = path.bond_ids;
            atom_ids.push(bond_to_remove.atom_ids.1);
            bond_ids.push(bond_to_remove.id);

            candidates.push(RingCandidate {
                atom_ids,
                bond_ids,
                len: path.len + 1,
            });
        }
    }
    candidates
}

/// Selects up to `cyclomatic_number` cycles forming a minimal basis using Gaussian elimination.
///
/// # Arguments
///
/// * `candidates` - Candidate cycles sorted by length.
/// * `cyclomatic_number` - Target number of independent cycles to keep.
/// * `bond_id_to_index` - Mapping from bond IDs to dense indices for bit-vector math.
///
/// # Returns
///
/// Minimal cycle basis expressed as a subset of the input candidates.
fn select_minimal_cycle_basis(
    mut candidates: Vec<RingCandidate>,
    cyclomatic_number: usize,
    bond_id_to_index: &HashMap<usize, usize>,
) -> Vec<RingCandidate> {
    candidates.sort_by_key(|c| c.len);

    let mut selected_rings = Vec::new();
    let mut basis: Vec<(BitVec, usize)> = Vec::new();

    for ring in candidates {
        let mut bitvec = BitVec::from_bond_ids(&ring.bond_ids, bond_id_to_index);

        for (basis_vec, pivot) in &basis {
            if bitvec.test(*pivot) {
                bitvec.xor(basis_vec);
            }
        }

        if let Some(pivot) = bitvec.leading_one() {
            basis.push((bitvec, pivot));
            basis.sort_by_key(|&(_, p)| p);
            selected_rings.push(ring);

            if selected_rings.len() == cyclomatic_number {
                break;
            }
        }
    }
    selected_rings
}

/// Marks atoms as ring members and records their smallest ring size.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule updated in-place.
fn annotate_atoms_with_ring_info(molecule: &mut AnnotatedMolecule) {
    for ring in &molecule.rings {
        let ring_size = ring.len() as u8;
        for &atom_id in ring {
            if let Some(props) = molecule.atoms.get_mut(atom_id) {
                props.is_in_ring = true;

                let current_smallest = props.smallest_ring_size.get_or_insert(ring_size);
                if ring_size < *current_smallest {
                    *current_smallest = ring_size;
                }
            }
        }
    }
}

/// Stores the path discovered between two atoms when a bond is removed.
struct PathData {
    /// Atom identifiers along the path (excluding the destination, which is implied).
    atom_ids: Vec<usize>,
    /// Bond identifiers traversed along the path.
    bond_ids: Vec<usize>,
    /// Path length measured in edges.
    len: usize,
}

/// BFS-based shortest path search that optionally excludes one bond from consideration.
///
/// # Arguments
///
/// * `molecule` - Annotated molecule providing adjacency and bond metadata.
/// * `start_id` - Starting atom identifier.
/// * `end_id` - Destination atom identifier.
/// * `excluded_bond_id` - Optional bond ID to ignore, simulating its removal.
/// * `workspace` - Reusable BFS buffers to avoid per-call allocations.
///
/// # Returns
///
/// A [`PathData`] instance if a route exists or `None` when the nodes disconnect.
fn shortest_path_bfs(
    molecule: &AnnotatedMolecule,
    start_id: usize,
    end_id: usize,
    excluded_bond_id: Option<usize>,
    workspace: &mut RingSearchWorkspace,
) -> Option<PathData> {
    workspace.reset();

    if start_id >= molecule.adjacency_with_bonds.len()
        || end_id >= molecule.adjacency_with_bonds.len()
    {
        return None;
    }

    workspace.visited[start_id] = true;
    workspace.queue.push_back(start_id);

    'outer: while let Some(current_id) = workspace.queue.pop_front() {
        for NeighborBond {
            neighbor_id,
            bond_id,
            ..
        } in &molecule.adjacency_with_bonds[current_id]
        {
            if Some(*bond_id) == excluded_bond_id {
                continue;
            }
            if !workspace.visited[*neighbor_id] {
                workspace.visited[*neighbor_id] = true;
                workspace.parent[*neighbor_id] = Some((current_id, *bond_id));
                workspace.queue.push_back(*neighbor_id);

                if *neighbor_id == end_id {
                    break 'outer;
                }
            }
        }
    }

    if !workspace.visited[end_id] {
        return None;
    }

    let mut atom_ids = Vec::new();
    let mut bond_ids = Vec::new();
    let mut len = 0;
    let mut cursor = end_id;

    while let Some((prev_id, bond_id)) = workspace.parent[cursor] {
        atom_ids.push(prev_id);
        bond_ids.push(bond_id);
        len += 1;
        cursor = prev_id;
        if cursor == start_id {
            break;
        }
    }
    atom_ids.reverse();
    bond_ids.reverse();

    Some(PathData {
        atom_ids,
        bond_ids,
        len,
    })
}

/// Counts the number of connected components in the molecular graph.
///
/// # Arguments
///
/// * `num_atoms` - Number of atoms in the graph.
/// * `adjacency` - Neighbor list for each atom.
///
/// # Returns
///
/// Count of disjoint components.
fn count_components(num_atoms: usize, adjacency: &[Vec<(usize, GraphBondOrder)>]) -> usize {
    let mut visited = vec![false; num_atoms];
    let mut components = 0;
    for i in 0..num_atoms {
        if !visited[i] {
            components += 1;
            let mut stack = vec![i];
            visited[i] = true;
            while let Some(current) = stack.pop() {
                for &(neighbor_id, _) in &adjacency[current] {
                    if !visited[neighbor_id] {
                        visited[neighbor_id] = true;
                        stack.push(neighbor_id);
                    }
                }
            }
        }
    }
    components
}

/// Sparse bit vector used for Gaussian elimination over GF(2).
#[derive(Clone, Debug, PartialEq, Eq)]
struct BitVec {
    /// Packed bit storage.
    data: Vec<u64>,
    /// Number of significant bits represented by the vector.
    size: usize,
}

impl BitVec {
    /// Creates a zero-initialized bit vector with the requested number of bits.
    ///
    /// # Arguments
    ///
    /// * `size` - Number of significant bits to represent.
    fn new(size: usize) -> Self {
        let words = size.div_ceil(64);
        Self {
            data: vec![0; words],
            size,
        }
    }

    /// Builds a bit vector with ones at positions mapped from the provided bond identifiers.
    ///
    /// # Arguments
    ///
    /// * `bond_ids` - Bonds participating in the cycle.
    /// * `bond_id_to_index` - Dense index lookup for bonds.
    fn from_bond_ids(bond_ids: &[usize], bond_id_to_index: &HashMap<usize, usize>) -> Self {
        let mut bitvec = Self::new(bond_id_to_index.len());
        for bond_id in bond_ids {
            if let Some(&index) = bond_id_to_index.get(bond_id) {
                let word_idx = index / 64;
                let bit_idx = index % 64;
                if word_idx < bitvec.data.len() {
                    bitvec.data[word_idx] |= 1u64 << bit_idx;
                }
            }
        }
        bitvec
    }

    /// XORs the bit vector with another vector of equal length.
    ///
    /// # Arguments
    ///
    /// * `other` - Right-hand operand.
    fn xor(&mut self, other: &Self) {
        debug_assert_eq!(self.data.len(), other.data.len());
        for (a, b) in self.data.iter_mut().zip(&other.data) {
            *a ^= *b;
        }
    }

    /// Tests whether the bit at `index` is set.
    ///
    /// # Arguments
    ///
    /// * `index` - Position to inspect.
    fn test(&self, index: usize) -> bool {
        if index >= self.size {
            return false;
        }
        let word_idx = index / 64;
        let bit_idx = index % 64;
        (self.data[word_idx] & (1u64 << bit_idx)) != 0
    }

    /// Returns the index of the most significant set bit, if any.
    ///
    /// # Returns
    ///
    /// `Some(index)` when a set bit exists or `None` when the vector is zero.
    fn leading_one(&self) -> Option<usize> {
        for (word_idx_rev, &word) in self.data.iter().enumerate().rev() {
            if word != 0 {
                let bit_pos_in_word = 63 - word.leading_zeros() as usize;
                return Some(word_idx_rev * 64 + bit_pos_in_word);
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::dreiding::typer::core::graph::MolecularGraph;
    use crate::forcefield::dreiding::typer::core::properties::{Element, GraphBondOrder};

    fn chain_graph(len: usize) -> MolecularGraph {
        let mut graph = MolecularGraph::new();
        for _ in 0..len {
            graph.add_atom(Element::C);
        }
        for i in 0..len.saturating_sub(1) {
            graph
                .add_bond(i, i + 1, GraphBondOrder::Single)
                .expect("valid bond");
        }
        graph
    }

    fn cycle_graph(len: usize) -> MolecularGraph {
        assert!(len >= 3, "cycles require at least three atoms");
        let mut graph = MolecularGraph::new();
        for _ in 0..len {
            graph.add_atom(Element::C);
        }
        for i in 0..len {
            let next = (i + 1) % len;
            graph
                .add_bond(i, next, GraphBondOrder::Single)
                .expect("valid cycle bond");
        }
        graph
    }

    fn fused_square_graph() -> MolecularGraph {
        let mut graph = MolecularGraph::new();
        for _ in 0..6 {
            graph.add_atom(Element::C);
        }

        let edges = [(0, 1), (1, 2), (2, 3), (3, 0), (2, 5), (5, 4), (4, 3)];

        for (u, v) in edges {
            graph
                .add_bond(u, v, GraphBondOrder::Single)
                .expect("valid edge");
        }

        graph
    }

    #[test]
    fn perceive_skips_acyclic_molecules() {
        let chain = chain_graph(4);
        let mut molecule = AnnotatedMolecule::new(&chain).expect("graph is valid");

        perceive(&mut molecule).expect("perception should succeed");

        assert!(
            molecule
                .atoms
                .iter()
                .all(|atom| !atom.is_in_ring && atom.smallest_ring_size.is_none())
        );
    }

    #[test]
    fn perceive_marks_ring_atoms_with_smallest_ring_size() {
        let square = cycle_graph(4);
        let mut molecule = AnnotatedMolecule::new(&square).expect("graph is valid");

        perceive(&mut molecule).expect("perception should succeed");

        for atom in &molecule.atoms {
            assert!(atom.is_in_ring, "atom {} should be in ring", atom.id);
            assert_eq!(atom.smallest_ring_size, Some(4));
        }
    }

    #[test]
    fn perceive_identifies_fused_squares_basis() {
        let graph = fused_square_graph();
        let mut molecule = AnnotatedMolecule::new(&graph).expect("graph is valid");

        perceive(&mut molecule).expect("perception should succeed");

        assert_eq!(molecule.rings.len(), 2, "expected two 4-cycles in basis");
        for ring in &molecule.rings {
            assert_eq!(ring.len(), 4, "each ring should have length 4");
        }

        for atom in &molecule.atoms {
            assert!(atom.is_in_ring, "atom {} should be in a ring", atom.id);
            assert_eq!(atom.smallest_ring_size, Some(4));
        }
    }

    #[test]
    fn shortest_path_bfs_finds_alternative_route_when_edge_removed() {
        let triangle = cycle_graph(3);
        let molecule = AnnotatedMolecule::new(&triangle).expect("graph is valid");

        let mut workspace = RingSearchWorkspace::new(molecule.atoms.len());

        let removed_bond_id = molecule
            .bonds
            .iter()
            .find(|bond| bond.atom_ids == (0, 1) || bond.atom_ids == (1, 0))
            .map(|bond| bond.id)
            .expect("triangle should contain 0-1 bond");

        let path = shortest_path_bfs(&molecule, 0, 1, Some(removed_bond_id), &mut workspace)
            .expect("path exists through third atom");

        assert_eq!(path.len, 2);
        assert_eq!(path.atom_ids, vec![0, 2]);
        assert_eq!(path.bond_ids.len(), 2);
    }

    #[test]
    fn count_components_detects_disconnected_fragments() {
        let adjacency = vec![
            vec![(1, GraphBondOrder::Single)],
            vec![(0, GraphBondOrder::Single)],
            vec![(3, GraphBondOrder::Single)],
            vec![(2, GraphBondOrder::Single)],
        ];

        assert_eq!(count_components(4, &adjacency), 2);
    }

    #[test]
    fn bitvec_supports_xor_and_leading_one() {
        let mut bond_map = HashMap::new();
        bond_map.insert(10usize, 0usize);
        bond_map.insert(20usize, 1usize);
        bond_map.insert(30usize, 2usize);

        let mut a = BitVec::from_bond_ids(&[10, 30], &bond_map);
        let b = BitVec::from_bond_ids(&[20, 30], &bond_map);

        assert!(a.test(0));
        assert!(a.test(2));
        assert!(b.test(1));
        assert!(b.test(2));

        a.xor(&b);
        assert!(a.test(0));
        assert!(a.test(1));
        assert!(!a.test(2));

        assert_eq!(a.leading_one(), Some(1));
    }
}
