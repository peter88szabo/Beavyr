//! Defines the molecular graph structure used for building molecules.
//!
//! This module provides the `MolecularGraph` container, which serves as the
//! primary input interface for the library. It captures atoms and bonds with
//! `GraphBondOrder` connectivity (Single, Double, Triple, Aromatic) before
//! perception begins.

use super::error::GraphValidationError;
use super::properties::{Element, GraphBondOrder};

/// Stores the identifier and element for a single atom within a
/// [`MolecularGraph`].
#[derive(Debug, Clone)]
pub struct AtomNode {
    /// Zero-based identifier assigned when the atom is inserted into the graph.
    pub id: usize,
    /// Chemical element represented by this node.
    pub element: Element,
}

/// Captures a bond between two atoms inside a [`MolecularGraph`].
#[derive(Debug, Clone)]
pub struct BondEdge {
    /// Zero-based identifier assigned sequentially as bonds are inserted.
    pub id: usize,
    /// Tuple of atom IDs that this bond connects.
    pub atom_ids: (usize, usize),
    /// Bond multiplicity recorded for the edge.
    pub order: GraphBondOrder,
}

/// Mutable graph of atoms and bonds supplied to the perception pipeline.
#[derive(Debug, Clone, Default)]
pub struct MolecularGraph {
    /// Collection of all atoms currently present in the graph.
    pub atoms: Vec<AtomNode>,
    /// Collection of all bonds currently present in the graph.
    pub bonds: Vec<BondEdge>,
}

impl MolecularGraph {
    /// Creates an empty graph ready for atom and bond insertion.
    ///
    /// # Returns
    ///
    /// A `MolecularGraph` with zero atoms and bonds.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::MolecularGraph;
    /// let graph = MolecularGraph::new();
    /// assert!(graph.atoms.is_empty());
    /// assert!(graph.bonds.is_empty());
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a new atom with the provided [`Element`] and returns its ID.
    ///
    /// # Arguments
    ///
    /// * `element` - Chemical element to assign to the node.
    ///
    /// # Returns
    ///
    /// The zero-based identifier for the newly inserted atom.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::{Element, MolecularGraph};
    /// let mut graph = MolecularGraph::new();
    /// let carbon_id = graph.add_atom(Element::C);
    /// assert_eq!(carbon_id, 0);
    /// ```
    pub fn add_atom(&mut self, element: Element) -> usize {
        let id = self.atoms.len();
        self.atoms.push(AtomNode { id, element });
        id
    }

    /// Adds a bond between two existing atoms.
    ///
    /// # Arguments
    ///
    /// * `atom1_id` - Identifier of the first atom.
    /// * `atom2_id` - Identifier of the second atom.
    /// * `order` - Bond multiplicity to record.
    ///
    /// # Returns
    ///
    /// The zero-based identifier assigned to the new bond.
    ///
    /// # Errors
    ///
    /// Returns [`GraphValidationError::MissingAtom`] if either atom ID has not
    /// been inserted, or [`GraphValidationError::SelfBondingAtom`] if both IDs
    /// refer to the same atom.
    ///
    /// # Examples
    ///
    /// ```
    /// use dreid_typer::{GraphBondOrder, Element, MolecularGraph};
    /// let mut graph = MolecularGraph::new();
    /// let c = graph.add_atom(Element::C);
    /// let o = graph.add_atom(Element::O);
    /// let bond_id = graph.add_bond(c, o, GraphBondOrder::Double).unwrap();
    /// assert_eq!(bond_id, 0);
    /// ```
    pub fn add_bond(
        &mut self,
        atom1_id: usize,
        atom2_id: usize,
        order: GraphBondOrder,
    ) -> Result<usize, GraphValidationError> {
        if atom1_id >= self.atoms.len() {
            return Err(GraphValidationError::MissingAtom { atom_id: atom1_id });
        }
        if atom2_id >= self.atoms.len() {
            return Err(GraphValidationError::MissingAtom { atom_id: atom2_id });
        }
        if atom1_id == atom2_id {
            return Err(GraphValidationError::SelfBondingAtom { atom_id: atom1_id });
        }

        let id = self.bonds.len();
        self.bonds.push(BondEdge {
            id,
            atom_ids: (atom1_id, atom2_id),
            order,
        });
        Ok(id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::dreiding::typer::core::error::GraphValidationError;
    use crate::forcefield::dreiding::typer::core::properties::{Element, GraphBondOrder};

    fn graph_with_atoms(elements: &[Element]) -> MolecularGraph {
        let mut graph = MolecularGraph::new();
        for &element in elements {
            graph.add_atom(element);
        }
        graph
    }

    #[test]
    fn molecular_graph_add_atom_assigns_sequential_ids() {
        let mut graph = MolecularGraph::new();

        let carbon_id = graph.add_atom(Element::C);
        let oxygen_id = graph.add_atom(Element::O);
        let nitrogen_id = graph.add_atom(Element::N);

        assert_eq!(carbon_id, 0);
        assert_eq!(oxygen_id, 1);
        assert_eq!(nitrogen_id, 2);
        assert_eq!(graph.atoms.len(), 3);
        assert_eq!(graph.atoms[0].element, Element::C);
        assert_eq!(graph.atoms[1].element, Element::O);
        assert_eq!(graph.atoms[2].element, Element::N);
    }

    #[test]
    fn molecular_graph_add_bond_registers_edge() {
        let mut graph = graph_with_atoms(&[Element::C, Element::O]);

        let bond_id = graph
            .add_bond(0, 1, GraphBondOrder::Double)
            .expect("adding a valid bond should succeed");

        assert_eq!(bond_id, 0);
        assert_eq!(graph.bonds.len(), 1);
        assert_eq!(graph.bonds[0].atom_ids, (0, 1));
        assert_eq!(graph.bonds[0].order, GraphBondOrder::Double);
    }

    #[test]
    fn molecular_graph_rejects_bond_with_missing_atom() {
        let mut graph = graph_with_atoms(&[Element::C]);

        let err = graph
            .add_bond(0, 1, GraphBondOrder::Single)
            .expect_err("bonding to a missing atom should fail");

        match err {
            GraphValidationError::MissingAtom { atom_id } => assert_eq!(atom_id, 1),
            _ => panic!("unexpected error returned: {err:?}"),
        }
    }

    #[test]
    fn molecular_graph_rejects_self_bonding() {
        let mut graph = graph_with_atoms(&[Element::C]);

        let err = graph
            .add_bond(0, 0, GraphBondOrder::Single)
            .expect_err("self bonds should be rejected");

        match err {
            GraphValidationError::SelfBondingAtom { atom_id } => assert_eq!(atom_id, 0),
            _ => panic!("unexpected error returned: {err:?}"),
        }
    }
}
