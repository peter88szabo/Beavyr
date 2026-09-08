//! Converts annotated molecules and assigned atom types into a full molecular topology.
//!
//! The builder stage takes the perception output and typing assignments, emitting atoms, bonds,
//! angles, proper dihedrals, and improper dihedrals expected by downstream force-field tooling.

use crate::forcefield::dreiding::typer::core::properties::{GraphBondOrder, Hybridization, TopologyBondOrder};
use crate::forcefield::dreiding::typer::core::topology::{
    Angle, Atom, Bond, ImproperDihedral, MolecularTopology, ProperDihedral,
};
use crate::forcefield::dreiding::typer::perception::{AnnotatedMolecule, ResonanceSystem};
use std::collections::HashSet;

/// Builds the `MolecularTopology` aggregate from perception results and atom-type labels.
///
/// This function effectively serializes the `AnnotatedMolecule` into the graph structures used
/// by force-field consumers by delegating to specialized helpers for each topology term.
///
/// # Arguments
///
/// * `annotated_molecule` - Molecule carrying ring, hybridization, and bonding metadata.
/// * `atom_types` - Slice of final atom-type names aligned with the molecule's atom ordering.
///
/// # Returns
///
/// A populated [`MolecularTopology`] containing atoms, bonds, angles, proper dihedrals, and
/// improper dihedrals.
pub fn build_topology(
    annotated_molecule: &AnnotatedMolecule,
    atom_types: &[String],
) -> MolecularTopology {
    let atoms = build_atoms(annotated_molecule, atom_types);
    let bonds = build_bonds(annotated_molecule);
    let angles = build_angles(annotated_molecule);
    let propers = build_propers(annotated_molecule);
    let impropers = build_impropers(annotated_molecule);

    MolecularTopology {
        atoms,
        bonds: bonds.into_iter().collect(),
        angles: angles.into_iter().collect(),
        propers: propers.into_iter().collect(),
        impropers: impropers.into_iter().collect(),
    }
}

/// Creates the atom list with element, type, and hybridization copies.
///
/// # Arguments
///
/// * `annotated_molecule` - Source molecule whose atoms provide structural metadata.
/// * `atom_types` - Slice of assigned atom-type labels.
fn build_atoms(annotated_molecule: &AnnotatedMolecule, atom_types: &[String]) -> Vec<Atom> {
    annotated_molecule
        .atoms
        .iter()
        .map(|ann_atom| Atom {
            id: ann_atom.id,
            element: ann_atom.element,
            atom_type: atom_types[ann_atom.id].clone(),
            hybridization: ann_atom.hybridization,
        })
        .collect()
}

/// Extracts unique bonds from the annotated molecule.
///
/// This function determines the final `TopologyBondOrder` by checking if a bond belongs to
/// any detected `ResonanceSystem`. If so, it is promoted to `Resonant`. Otherwise, the
/// Kekulized order (`Single`, `Double`, `Triple`) is used.
fn build_bonds(annotated_molecule: &AnnotatedMolecule) -> HashSet<Bond> {
    let resonant_bond_ids: HashSet<usize> = annotated_molecule
        .resonance_systems
        .iter()
        .flat_map(|sys: &ResonanceSystem| sys.bond_ids.iter())
        .copied()
        .collect();

    annotated_molecule
        .bonds
        .iter()
        .map(|edge| {
            let topology_order = if resonant_bond_ids.contains(&edge.id) {
                TopologyBondOrder::Resonant
            } else {
                match edge.order {
                    GraphBondOrder::Single => TopologyBondOrder::Single,
                    GraphBondOrder::Double => TopologyBondOrder::Double,
                    GraphBondOrder::Triple => TopologyBondOrder::Triple,
                    GraphBondOrder::Aromatic => TopologyBondOrder::Single, // Fallback; should not occur here
                }
            };

            Bond::new(edge.atom_ids.0, edge.atom_ids.1, topology_order)
        })
        .collect()
}

/// Generates all angle triplets by enumerating neighbor pairs around each atom.
fn build_angles(annotated_molecule: &AnnotatedMolecule) -> HashSet<Angle> {
    let mut angles = HashSet::new();
    for j in 0..annotated_molecule.atoms.len() {
        let neighbors = &annotated_molecule.adjacency[j];
        if neighbors.len() < 2 {
            continue;
        }
        for i in 0..neighbors.len() {
            for k in (i + 1)..neighbors.len() {
                let atom_i_id = neighbors[i].0;
                let atom_k_id = neighbors[k].0;
                angles.insert(Angle::new(atom_i_id, j, atom_k_id));
            }
        }
    }
    angles
}

/// Builds proper dihedrals by extending each bond to its neighboring atoms.
fn build_propers(annotated_molecule: &AnnotatedMolecule) -> HashSet<ProperDihedral> {
    let mut propers = HashSet::new();
    for bond_jk in &annotated_molecule.bonds {
        let (j, k) = bond_jk.atom_ids;

        for &(i, _) in &annotated_molecule.adjacency[j] {
            if i == k {
                continue;
            }
            for &(l, _) in &annotated_molecule.adjacency[k] {
                if l == j || l == i {
                    continue;
                }
                propers.insert(ProperDihedral::new(i, j, k, l));
            }
        }
    }
    propers
}

/// Builds improper dihedrals for planar degree-three centers with SP2-like hybridization.
fn build_impropers(annotated_molecule: &AnnotatedMolecule) -> HashSet<ImproperDihedral> {
    let mut impropers = HashSet::new();
    for atom in &annotated_molecule.atoms {
        if atom.degree == 3
            && matches!(
                atom.hybridization,
                Hybridization::SP2 | Hybridization::Resonant
            )
        {
            let neighbors = &annotated_molecule.adjacency[atom.id];
            let p1 = neighbors[0].0;
            let p2 = neighbors[1].0;
            let p3 = neighbors[2].0;
            impropers.insert(ImproperDihedral::new(p1, p2, atom.id, p3));
        }
    }
    impropers
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::dreiding::typer::core::graph::MolecularGraph;
    use crate::forcefield::dreiding::typer::core::properties::{Element, GraphBondOrder, TopologyBondOrder};
    use crate::forcefield::dreiding::typer::perception::ResonanceSystem;
    use std::collections::HashSet;

    fn planar_fragment() -> (AnnotatedMolecule, Vec<String>) {
        let mut graph = MolecularGraph::new();
        let c_left = graph.add_atom(Element::C);
        let c_center = graph.add_atom(Element::C);
        let c_right = graph.add_atom(Element::C);
        let n_cap = graph.add_atom(Element::N);
        let c_tail = graph.add_atom(Element::C);
        let h_tail = graph.add_atom(Element::H);

        graph
            .add_bond(c_left, c_center, GraphBondOrder::Single)
            .expect("valid bond");
        graph
            .add_bond(c_center, c_right, GraphBondOrder::Double)
            .expect("valid bond");
        graph
            .add_bond(c_center, n_cap, GraphBondOrder::Single)
            .expect("valid bond");
        graph
            .add_bond(c_right, c_tail, GraphBondOrder::Single)
            .expect("valid bond");
        graph
            .add_bond(c_tail, h_tail, GraphBondOrder::Single)
            .expect("valid bond");

        let mut molecule = AnnotatedMolecule::new(&graph).expect("graph should be valid");
        molecule.atoms[c_center].hybridization = Hybridization::SP2;

        molecule.resonance_systems.push(ResonanceSystem {
            atom_ids: vec![c_center, c_right, n_cap],
            bond_ids: vec![1, 2],
        });

        let atom_types = vec![
            "C_SP2_EDGE".to_string(),
            "C_R".to_string(),
            "C_SP3".to_string(),
            "N_R".to_string(),
            "C_ALK".to_string(),
            "H_".to_string(),
        ];

        (molecule, atom_types)
    }

    #[test]
    fn build_atoms_uses_atom_ids_to_assign_types() {
        let (molecule, atom_types) = planar_fragment();

        let atoms = build_atoms(&molecule, &atom_types);

        assert_eq!(atoms.len(), molecule.atoms.len());
        assert_eq!(atoms[1].atom_type, "C_R");
        assert_eq!(atoms[5].atom_type, "H_");
        assert_eq!(atoms[1].hybridization, Hybridization::SP2);
    }

    #[test]
    fn build_bonds_assigns_resonant_order_to_system_bonds() {
        let (molecule, _) = planar_fragment();

        let bonds = build_bonds(&molecule);

        assert_eq!(bonds.len(), molecule.bonds.len());

        assert!(bonds.contains(&Bond::new(0, 1, TopologyBondOrder::Single)));
        assert!(bonds.contains(&Bond::new(1, 2, TopologyBondOrder::Resonant)));
        assert!(bonds.contains(&Bond::new(1, 3, TopologyBondOrder::Resonant)));
        assert!(bonds.contains(&Bond::new(2, 4, TopologyBondOrder::Single)));
    }

    #[test]
    fn build_angles_generates_all_neighbor_pairs() {
        let (molecule, _) = planar_fragment();

        let angles = build_angles(&molecule);
        let expected: HashSet<_> = vec![
            Angle::new(0, 1, 2),
            Angle::new(0, 1, 3),
            Angle::new(2, 1, 3),
            Angle::new(1, 2, 4),
            Angle::new(2, 4, 5),
        ]
        .into_iter()
        .collect();

        assert_eq!(angles, expected);
    }

    #[test]
    fn build_propers_emits_all_valid_dihedrals() {
        let (molecule, _) = planar_fragment();

        let propers = build_propers(&molecule);
        let expected: HashSet<_> = vec![
            ProperDihedral::new(0, 1, 2, 4),
            ProperDihedral::new(3, 1, 2, 4),
            ProperDihedral::new(1, 2, 4, 5),
        ]
        .into_iter()
        .collect();

        assert_eq!(propers, expected);
    }

    #[test]
    fn build_impropers_targets_planar_degree_three_centers() {
        let (molecule, _) = planar_fragment();

        let impropers = build_impropers(&molecule);
        let expected: HashSet<_> = vec![ImproperDihedral::new(0, 2, 1, 3)]
            .into_iter()
            .collect();

        assert_eq!(impropers, expected);
    }
}
