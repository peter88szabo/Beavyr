//! Internal model shared across perception stages to annotate atoms, bonds, rings, and adjacency.
//!
//! The structures defined here wrap the raw `MolecularGraph` with mutable fields that each
//! perception pass enriches before the typing engine consumes them.

use crate::forcefield::dreiding::typer::core::error::GraphValidationError;
use crate::forcefield::dreiding::typer::core::graph::{BondEdge, MolecularGraph};
use crate::forcefield::dreiding::typer::core::properties::{Element, GraphBondOrder, Hybridization};

/// Neighbor descriptor bundling atom connectivity with the originating bond ID.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct NeighborBond {
    /// Adjacent atom identifier.
    pub neighbor_id: usize,
    /// Source bond identifier matching [`BondEdge::id`].
    pub bond_id: usize,
    /// Bond order for the edge.
    pub order: GraphBondOrder,
}

/// Perception-friendly atom record that stores both graph identity and inferred properties.
#[derive(Debug, Clone)]
pub struct AnnotatedAtom {
    /// Zero-based identifier matching the source [`MolecularGraph`].
    pub id: usize,
    /// Chemical element of the atom.
    pub element: Element,

    /// Current formal charge assigned by electron perception.
    pub formal_charge: i8,
    /// Number of lone pairs tracked for hybridization and resonance logic.
    pub lone_pairs: u8,
    /// Graph degree computed during adjacency building.
    pub degree: u8,

    /// Whether the atom lies on any ring identified so far.
    pub is_in_ring: bool,
    /// Size of the smallest ring containing the atom, if any.
    pub smallest_ring_size: Option<u8>,

    /// Flag set once aromaticity perception confirms Huckel criteria for this atom.
    pub is_aromatic: bool,
    /// Flag set for anti-aromatic atoms that should avoid resonance promotion.
    pub is_anti_aromatic: bool,

    /// Marks atoms that participate in a resonance system (aromatic or functional group).
    pub is_resonant: bool,
    /// Remembers whether the atom participated in any aromatic input bonds pre-Kekulé.
    pub has_aromatic_edge: bool,

    /// Steric number derived from lone pairs and neighbors for VSEPR calculations.
    pub steric_number: u8,
    /// Current hybridization assignment, defaulting to [`Hybridization::Unknown`].
    pub hybridization: Hybridization,
}

/// Convenience alias representing a ring as a list of atom identifiers.
pub type Ring = Vec<usize>;

/// Represents a subset of atoms and bonds that form a delocalized electron system.
///
/// This structure is populated by the `aromaticity` (for rings) and `resonance`
/// (for linear conjugated groups like carboxylates) passes. It is used by the
/// builder to identify which bonds should be emitted with `TopologyBondOrder::Resonant`.
#[derive(Debug, Clone)]
pub struct ResonanceSystem {
    /// IDs of atoms participating in this resonance system.
    #[allow(dead_code)]
    pub atom_ids: Vec<usize>,
    /// IDs of bonds (`BondEdge::id`) that are part of the delocalized system.
    pub bond_ids: Vec<usize>,
}

/// Annotated molecular container combining atom metadata, bonds, adjacency, and ring sets.
#[derive(Debug, Clone)]
pub struct AnnotatedMolecule {
    /// All atoms with perception-specific annotations.
    pub atoms: Vec<AnnotatedAtom>,
    /// Copy of the graph bonds to provide stable IDs and connectivity.
    pub bonds: Vec<BondEdge>,
    /// Adjacency list capturing neighbor IDs and bond orders.
    pub adjacency: Vec<Vec<(usize, GraphBondOrder)>>,
    /// Adjacency list that also records the bond ID for each neighbor edge.
    pub adjacency_with_bonds: Vec<Vec<NeighborBond>>,
    /// Collection of rings discovered during perception.
    pub rings: Vec<Ring>,
    /// Collection of all identified resonance systems.
    pub resonance_systems: Vec<ResonanceSystem>,
}

impl AnnotatedMolecule {
    /// Builds an annotated molecule from a validated [`MolecularGraph`].
    ///
    /// Initializes adjacency lists for every atom, clones the bond list, and seeds default
    /// annotations that later perception passes will populate.
    ///
    /// # Arguments
    ///
    /// * `graph` - Source molecular graph supplied by the caller.
    ///
    /// # Returns
    ///
    /// A new [`AnnotatedMolecule`] containing adjacency, cloned bonds, and zeroed annotations.
    ///
    /// # Errors
    ///
    /// Returns [`GraphValidationError::MissingAtom`] if any bond endpoint references an atom index
    /// outside the graph's atom list.
    pub fn new(graph: &MolecularGraph) -> Result<Self, GraphValidationError> {
        let mut adjacency = vec![vec![]; graph.atoms.len()];
        let mut adjacency_with_bonds = vec![vec![]; graph.atoms.len()];
        for bond in &graph.bonds {
            let (u, v) = bond.atom_ids;

            if u >= graph.atoms.len() {
                return Err(GraphValidationError::MissingAtom { atom_id: u });
            }
            if v >= graph.atoms.len() {
                return Err(GraphValidationError::MissingAtom { atom_id: v });
            }

            adjacency[u].push((v, bond.order));
            adjacency[v].push((u, bond.order));

            let neighbor_entry_u = NeighborBond {
                neighbor_id: v,
                bond_id: bond.id,
                order: bond.order,
            };
            let neighbor_entry_v = NeighborBond {
                neighbor_id: u,
                bond_id: bond.id,
                order: bond.order,
            };

            adjacency_with_bonds[u].push(neighbor_entry_u);
            adjacency_with_bonds[v].push(neighbor_entry_v);
        }

        let atoms = graph
            .atoms
            .iter()
            .map(|node| AnnotatedAtom {
                id: node.id,
                element: node.element,
                degree: adjacency[node.id].len() as u8,
                formal_charge: 0,
                lone_pairs: 0,
                is_in_ring: false,
                smallest_ring_size: None,
                is_aromatic: false,
                is_anti_aromatic: false,
                is_resonant: false,
                has_aromatic_edge: false,
                steric_number: 0,
                hybridization: Hybridization::Unknown,
            })
            .collect();

        Ok(Self {
            atoms,
            bonds: graph.bonds.clone(),
            adjacency,
            adjacency_with_bonds,
            rings: Vec::new(),
            resonance_systems: Vec::new(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::dreiding::typer::core::graph::AtomNode;
    use crate::forcefield::dreiding::typer::core::properties::GraphBondOrder;

    fn water_like_graph() -> MolecularGraph {
        let mut graph = MolecularGraph::new();
        let oxygen = graph.add_atom(Element::O);
        let hydrogen1 = graph.add_atom(Element::H);
        let hydrogen2 = graph.add_atom(Element::H);

        graph
            .add_bond(oxygen, hydrogen1, GraphBondOrder::Single)
            .expect("valid bond should be created");
        graph
            .add_bond(oxygen, hydrogen2, GraphBondOrder::Single)
            .expect("valid bond should be created");

        graph
    }

    #[test]
    fn annotated_molecule_new_populates_adjacency_and_atoms() {
        let graph = water_like_graph();
        let molecule = AnnotatedMolecule::new(&graph).expect("graph should be valid");

        let adjacency_sizes: Vec<_> = molecule
            .adjacency
            .iter()
            .map(|neighbors| neighbors.len())
            .collect();
        assert_eq!(adjacency_sizes, vec![2, 1, 1]);

        let adjacency_with_bonds_sizes: Vec<_> = molecule
            .adjacency_with_bonds
            .iter()
            .map(|neighbors| neighbors.len())
            .collect();
        assert_eq!(adjacency_with_bonds_sizes, vec![2, 1, 1]);

        let oxygen_neighbors = &molecule.adjacency[0];
        assert!(oxygen_neighbors.contains(&(1, GraphBondOrder::Single)));
        assert!(oxygen_neighbors.contains(&(2, GraphBondOrder::Single)));

        let oxygen_neighbors_with_bonds = &molecule.adjacency_with_bonds[0];
        assert_eq!(oxygen_neighbors_with_bonds.len(), 2);
        for neighbor in oxygen_neighbors_with_bonds {
            assert_eq!(neighbor.order, GraphBondOrder::Single);
            assert!(neighbor.neighbor_id == 1 || neighbor.neighbor_id == 2);
        }

        let h1_neighbors_with_bonds = &molecule.adjacency_with_bonds[1];
        assert_eq!(h1_neighbors_with_bonds.len(), 1);
        let h1_neighbor = &h1_neighbors_with_bonds[0];
        assert_eq!(h1_neighbor.neighbor_id, 0);
        assert_eq!(h1_neighbor.order, GraphBondOrder::Single);

        let h2_neighbors_with_bonds = &molecule.adjacency_with_bonds[2];
        assert_eq!(h2_neighbors_with_bonds.len(), 1);
        let h2_neighbor = &h2_neighbors_with_bonds[0];
        assert_eq!(h2_neighbor.neighbor_id, 0);
        assert_eq!(h2_neighbor.order, GraphBondOrder::Single);

        let oxygen = &molecule.atoms[0];
        assert_eq!(oxygen.element, Element::O);
        assert_eq!(oxygen.degree, 2);
        assert_eq!(oxygen.formal_charge, 0);
        assert_eq!(oxygen.lone_pairs, 0);
        assert!(!oxygen.is_in_ring);
        assert_eq!(oxygen.smallest_ring_size, None);
        assert!(!oxygen.is_aromatic);
        assert!(!oxygen.is_anti_aromatic);
        assert!(!oxygen.is_resonant);
        assert_eq!(oxygen.steric_number, 0);
        assert_eq!(oxygen.hybridization, Hybridization::Unknown);
        assert!(molecule.resonance_systems.is_empty());
    }

    #[test]
    fn annotated_molecule_new_detects_invalid_bond_endpoints() {
        let graph = MolecularGraph {
            atoms: vec![AtomNode {
                id: 0,
                element: Element::C,
            }],
            bonds: vec![BondEdge {
                id: 0,
                atom_ids: (0, 2),
                order: GraphBondOrder::Single,
            }],
        };

        let err = AnnotatedMolecule::new(&graph).expect_err("invalid bond must fail");

        match err {
            GraphValidationError::MissingAtom { atom_id } => assert_eq!(atom_id, 2),
            _ => panic!("unexpected error returned: {err:?}"),
        }
    }
}
