//! Represents the canonical topology emitted after typing completes.
//!
//! This module defines the output structures (`MolecularTopology`, `Bond`, etc.)
//! that are ready for serialization or consumption by force field engines.
//! These types use `TopologyBondOrder`, which includes physical properties like
//! resonance, unlike the input graph.

use super::properties::{Element, Hybridization, TopologyBondOrder};

/// Canonical topology produced after the typer assigns atom types and torsions.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MolecularTopology {
    /// A list of all atoms with their final assigned properties.
    pub atoms: Vec<Atom>,
    /// A list of all bonds.
    pub bonds: Vec<Bond>,
    /// A list of all three-atom angles.
    pub angles: Vec<Angle>,
    /// A list of all proper dihedral angles (torsions).
    pub propers: Vec<ProperDihedral>,
    /// A list of all improper dihedral angles (out-of-plane bends).
    pub impropers: Vec<ImproperDihedral>,
}

/// Atom entry emitted in the final topology, combining identity and typing.
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    /// The unique identifier of the atom.
    pub id: usize,
    /// The chemical element.
    pub element: Element,
    /// The final, assigned DREIDING atom type string.
    pub atom_type: String,
    /// The perceived hybridization state.
    pub hybridization: Hybridization,
}

/// Bond entry emitted in the final topology.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Bond {
    /// The IDs of the two atoms, sorted to ensure a canonical representation.
    pub atom_ids: (usize, usize),
    /// The order of the bond.
    pub order: TopologyBondOrder,
}

impl Bond {
    /// Creates a new bond with atom IDs sorted to a canonical order.
    pub fn new(id1: usize, id2: usize, order: TopologyBondOrder) -> Self {
        let atom_ids = if id1 < id2 { (id1, id2) } else { (id2, id1) };
        Self { atom_ids, order }
    }
}

/// Angle entry emitted in the final topology.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Angle {
    /// The IDs of the three atoms (`end1`, `center`, `end2`), with end atoms sorted.
    pub atom_ids: (usize, usize, usize),
}

impl Angle {
    /// Creates a new angle with end atoms sorted to a canonical order.
    pub fn new(id1: usize, center_id: usize, id2: usize) -> Self {
        let atom_ids = if id1 < id2 {
            (id1, center_id, id2)
        } else {
            (id2, center_id, id1)
        };
        Self { atom_ids }
    }
}

/// Proper dihedral entry emitted in the final topology.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ProperDihedral {
    /// The IDs of the four atoms (`a-b-c-d`), sorted lexicographically.
    pub atom_ids: (usize, usize, usize, usize),
}

impl ProperDihedral {
    /// Creates a new proper dihedral with atom IDs sorted lexicographically.
    pub fn new(a: usize, b: usize, c: usize, d: usize) -> Self {
        let fwd = (a, b, c, d);
        let rev = (d, c, b, a);
        let atom_ids = if fwd <= rev { fwd } else { rev };
        Self { atom_ids }
    }
}

/// Improper dihedral entry emitted in the final topology.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ImproperDihedral {
    /// The IDs of the four atoms (`plane1`, `plane2`, `center`, `plane3`),
    /// with plane atoms sorted.
    pub atom_ids: (usize, usize, usize, usize),
}

impl ImproperDihedral {
    /// Creates a new improper dihedral with plane atoms sorted.
    pub fn new(p1: usize, p2: usize, center: usize, p3: usize) -> Self {
        let mut plane_ids = [p1, p2, p3];
        plane_ids.sort_unstable();
        let atom_ids = (plane_ids[0], plane_ids[1], center, plane_ids[2]);
        Self { atom_ids }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::dreiding::typer::core::properties::TopologyBondOrder;

    #[test]
    fn bond_new_sorts_atom_ids() {
        let bond = Bond::new(4, 1, TopologyBondOrder::Triple);
        assert_eq!(bond.atom_ids, (1, 4));
        assert_eq!(bond.order, TopologyBondOrder::Triple);
    }

    #[test]
    fn angle_new_orders_terminal_atoms() {
        let angle = Angle::new(7, 3, 2);
        assert_eq!(angle.atom_ids, (2, 3, 7));
    }

    #[test]
    fn proper_dihedral_new_canonicalizes_orientation() {
        let forward = ProperDihedral::new(1, 2, 3, 4);
        let reversed = ProperDihedral::new(4, 3, 2, 1);

        assert_eq!(forward.atom_ids, reversed.atom_ids);
        assert_eq!(forward.atom_ids, (1, 2, 3, 4));
    }

    #[test]
    fn improper_dihedral_new_sorts_plane_atoms() {
        let improper = ImproperDihedral::new(9, 1, 5, 4);
        assert_eq!(improper.atom_ids, (1, 4, 5, 9));
    }

    #[test]
    fn molecular_topology_default_is_empty() {
        let topology = MolecularTopology::default();

        assert!(topology.atoms.is_empty());
        assert!(topology.bonds.is_empty());
        assert!(topology.angles.is_empty());
        assert!(topology.propers.is_empty());
        assert!(topology.impropers.is_empty());
    }
}
