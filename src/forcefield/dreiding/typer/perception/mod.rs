//! Coordinates the sequential perception pipeline that annotates molecules prior to typing.
//!
//! This module wires the specialized perception stages—ring detection, Kekulé expansion,
//! electron bookkeeping, aromaticity, resonance, and hybridization—into a single pass that
//! populates an [`AnnotatedMolecule`] for downstream typing.

mod aromaticity;
mod electrons;
mod hybridization;
mod kekulize;
mod model;
mod resonance;
mod rings;

pub use model::{AnnotatedAtom, AnnotatedMolecule, ResonanceSystem};

use crate::forcefield::dreiding::typer::core::error::{PerceptionError, TyperError};
use crate::forcefield::dreiding::typer::core::graph::MolecularGraph;

type PerceptionStepFn = fn(&mut AnnotatedMolecule) -> Result<(), PerceptionError>;
type PerceptionStep = (&'static str, PerceptionStepFn);

/// Runs the full perception pipeline and returns an annotated molecule.
///
/// The function constructs an [`AnnotatedMolecule`] from the input graph, executes the fixed set
/// of perception stages in order, and reports any failure with the offending step name folded
/// into the [`TyperError`].
///
/// # Arguments
///
/// * `graph` - Validated molecular graph containing atoms and bonds.
///
/// # Returns
///
/// Fully annotated molecule that records ring membership, electron bookkeeping, aromaticity,
/// resonance, and hybridization properties for every atom.
///
/// # Errors
///
/// Returns [`TyperError::InvalidInput`] when the graph contains invalid bonding, or
/// [`TyperError::PerceptionFailed`] when any perception stage emits a [`PerceptionError`].
pub fn perceive(graph: &MolecularGraph) -> Result<AnnotatedMolecule, TyperError> {
    let mut molecule = AnnotatedMolecule::new(graph).map_err(TyperError::InvalidInput)?;

    let pipeline: [PerceptionStep; 6] = [
        ("Rings", rings::perceive),
        ("Kekulization", kekulize::perceive),
        ("Electrons", electrons::perceive),
        ("Aromaticity", aromaticity::perceive),
        ("Resonance", resonance::perceive),
        ("Hybridization", hybridization::perceive),
    ];

    for (name, step_fn) in pipeline {
        step_fn(&mut molecule).map_err(|source| TyperError::PerceptionFailed {
            step: name.to_string(),
            source,
        })?;
    }

    Ok(molecule)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::dreiding::typer::core::properties::{Element, GraphBondOrder, Hybridization};

    fn benzene_graph() -> MolecularGraph {
        let mut graph = MolecularGraph::new();
        let carbons: Vec<_> = (0..6).map(|_| graph.add_atom(Element::C)).collect();
        let hydrogens: Vec<_> = (0..6).map(|_| graph.add_atom(Element::H)).collect();

        let ring_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
        for &(u, v) in &ring_edges {
            graph
                .add_bond(carbons[u], carbons[v], GraphBondOrder::Aromatic)
                .expect("valid aromatic bond in benzene ring");
        }
        for i in 0..6 {
            graph
                .add_bond(carbons[i], hydrogens[i], GraphBondOrder::Single)
                .expect("valid C-H bond");
        }

        graph
    }

    fn acridine_graph() -> MolecularGraph {
        let mut graph = MolecularGraph::new();
        let c1 = graph.add_atom(Element::C);
        let c2 = graph.add_atom(Element::C);
        let c3 = graph.add_atom(Element::C);
        let c4 = graph.add_atom(Element::C);
        let c5 = graph.add_atom(Element::C);
        let c6 = graph.add_atom(Element::C);
        let c7 = graph.add_atom(Element::C);
        let c8 = graph.add_atom(Element::C);
        let c9 = graph.add_atom(Element::C);
        let n10 = graph.add_atom(Element::N);
        let c11 = graph.add_atom(Element::C);
        let c12 = graph.add_atom(Element::C);
        let c13 = graph.add_atom(Element::C);
        let c14 = graph.add_atom(Element::C);

        let h1 = graph.add_atom(Element::H);
        let h2 = graph.add_atom(Element::H);
        let h3 = graph.add_atom(Element::H);
        let h4 = graph.add_atom(Element::H);
        let h7 = graph.add_atom(Element::H);
        let h8 = graph.add_atom(Element::H);
        let h11 = graph.add_atom(Element::H);
        let h12 = graph.add_atom(Element::H);
        let h13 = graph.add_atom(Element::H);
        let h14 = graph.add_atom(Element::H);

        let aromatic_edges = [
            (c1, c2),
            (c2, c3),
            (c3, c4),
            (c4, c5),
            (c5, c6),
            (c6, c1),
            (c5, c7),
            (c7, c8),
            (c8, c9),
            (c9, n10),
            (n10, c6),
            (n10, c11),
            (c11, c12),
            (c12, c13),
            (c13, c14),
            (c14, c9),
        ];
        for &(u, v) in &aromatic_edges {
            graph
                .add_bond(u, v, GraphBondOrder::Aromatic)
                .expect("valid aromatic bond in acridine core");
        }

        let single_edges = [
            (c1, h1),
            (c2, h2),
            (c3, h3),
            (c4, h4),
            (c7, h7),
            (c8, h8),
            (c11, h11),
            (c12, h12),
            (c13, h13),
            (c14, h14),
        ];
        for &(u, v) in &single_edges {
            graph
                .add_bond(u, v, GraphBondOrder::Single)
                .expect("valid sigma bond in acridine");
        }

        graph
    }

    fn aromatic_bond_outside_ring_graph() -> MolecularGraph {
        let mut graph = MolecularGraph::new();
        let c1 = graph.add_atom(Element::C);
        let c2 = graph.add_atom(Element::C);
        let h1 = graph.add_atom(Element::H);
        let h2 = graph.add_atom(Element::H);

        graph
            .add_bond(c1, c2, GraphBondOrder::Aromatic)
            .expect("valid aromatic bond");
        graph
            .add_bond(c1, h1, GraphBondOrder::Single)
            .expect("valid C-H bond");
        graph
            .add_bond(c2, h2, GraphBondOrder::Single)
            .expect("valid C-H bond");

        graph
    }

    fn pyrimidine_aromatic_graph() -> MolecularGraph {
        let mut graph = MolecularGraph::new();

        let n1 = graph.add_atom(Element::N);
        let c2 = graph.add_atom(Element::C);
        let n3 = graph.add_atom(Element::N);
        let c4 = graph.add_atom(Element::C);
        let c5 = graph.add_atom(Element::C);
        let c6 = graph.add_atom(Element::C);
        let h2 = graph.add_atom(Element::H);
        let h4 = graph.add_atom(Element::H);
        let h5 = graph.add_atom(Element::H);
        let h6 = graph.add_atom(Element::H);

        let aromatic_edges = [(n1, c2), (c2, n3), (n3, c4), (c4, c5), (c5, c6), (c6, n1)];
        for &(u, v) in &aromatic_edges {
            graph
                .add_bond(u, v, GraphBondOrder::Aromatic)
                .expect("valid aromatic bond in pyrimidine ring");
        }

        graph
            .add_bond(c2, h2, GraphBondOrder::Single)
            .expect("C2-H bond");
        graph
            .add_bond(c4, h4, GraphBondOrder::Single)
            .expect("C4-H bond");
        graph
            .add_bond(c5, h5, GraphBondOrder::Single)
            .expect("C5-H bond");
        graph
            .add_bond(c6, h6, GraphBondOrder::Single)
            .expect("C6-H bond");

        graph
    }

    #[test]
    fn perception_pipeline_assigns_benzene_properties() {
        let graph = benzene_graph();
        let molecule = perceive(&graph).expect("perception pipeline should succeed");

        assert_eq!(molecule.rings.len(), 1, "benzene must yield a single ring");
        for (idx, atom) in molecule.atoms.iter().enumerate() {
            match atom.element {
                Element::C => {
                    assert!(atom.is_in_ring, "carbon {idx} must be in the ring");
                    assert!(atom.is_aromatic, "carbon {idx} must be aromatic");
                    assert!(atom.is_resonant, "carbon {idx} must be resonant");
                    assert_eq!(
                        atom.hybridization,
                        Hybridization::Resonant,
                        "carbon {idx} should end Resonant"
                    );
                    assert_eq!(atom.steric_number, 3);
                }
                Element::H => {
                    assert_eq!(atom.hybridization, Hybridization::None);
                    assert_eq!(atom.steric_number, 1);
                }
                other => panic!("unexpected element in benzene fixture: {other:?}"),
            }
        }

        assert!(
            molecule
                .bonds
                .iter()
                .all(|bond| bond.order != GraphBondOrder::Aromatic),
            "all aromatic bonds should be Kekulé-expanded"
        );

        assert!(
            !molecule.resonance_systems.is_empty(),
            "Benzene should have a resonance system"
        );
    }

    #[test]
    fn perception_pipeline_marks_acridine_as_aromatic() {
        let graph = acridine_graph();
        let molecule = perceive(&graph).expect("perception pipeline should succeed");

        assert!(
            molecule.rings.len() >= 3,
            "acridine should yield multiple rings, found {}",
            molecule.rings.len()
        );

        let doubleless: Vec<_> = molecule
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, atom)| atom.element != Element::H)
            .filter(|&(idx, _)| {
                !molecule.adjacency[idx]
                    .iter()
                    .filter(|&&(neighbor_id, _)| molecule.atoms[neighbor_id].is_in_ring)
                    .any(|&(_, order)| order == GraphBondOrder::Double)
            })
            .map(|(idx, _)| idx)
            .collect();
        assert!(
            doubleless
                .iter()
                .all(|&idx| molecule.atoms[idx].is_resonant),
            "heavy ring atoms without double bonds must originate from aromatic input: {:?}",
            doubleless
        );

        let resonant_heavy: Vec<_> = molecule
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, atom)| atom.element != Element::H)
            .filter(|(_, atom)| atom.is_resonant)
            .map(|(idx, _)| idx)
            .collect();

        assert_eq!(
            resonant_heavy.len(),
            14,
            "all heavy atoms in acridine should be flagged resonant"
        );

        let non_ring_heavy: Vec<_> = molecule
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, atom)| atom.element != Element::H)
            .filter(|(_, atom)| !atom.is_in_ring)
            .map(|(idx, _)| idx)
            .collect();

        assert!(
            non_ring_heavy.is_empty(),
            "every heavy atom in acridine should belong to a ring: {:?}",
            non_ring_heavy
        );

        let aromatic_heavy: Vec<_> = molecule
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, atom)| atom.element != Element::H)
            .filter(|(_, atom)| atom.is_aromatic)
            .map(|(idx, _)| idx)
            .collect();

        assert_eq!(
            aromatic_heavy.len(),
            14,
            "all heavy atoms in acridine should be aromatic"
        );

        for idx in aromatic_heavy {
            assert!(
                molecule.atoms[idx].is_resonant,
                "aromatic atom {idx} should be marked resonant"
            );
        }
    }

    #[test]
    fn pipeline_reports_step_name_when_kekulization_fails() {
        let graph = aromatic_bond_outside_ring_graph();
        let err = perceive(&graph).expect_err("pipeline should fail before completion");

        match err {
            TyperError::PerceptionFailed { step, source } => {
                assert_eq!(step, "Kekulization");
                match source {
                    PerceptionError::KekulizationFailed { message } => {
                        assert!(
                            message.contains("not in a ring"),
                            "error message should describe the missing ring context"
                        );
                    }
                    other => panic!("unexpected perception error emitted: {other:?}"),
                }
            }
            other => panic!("expected PerceptionFailed error, got {other:?}"),
        }
    }

    #[test]
    fn pyrimidine_aromatic_input_is_detected() {
        let graph = pyrimidine_aromatic_graph();
        let molecule = perceive(&graph).expect("perception pipeline should succeed");

        let ring_atoms = [0usize, 1, 2, 3, 4, 5];

        for &atom_id in &ring_atoms {
            assert!(
                molecule.atoms[atom_id].is_resonant,
                "Ring atom {atom_id} should be resonant"
            );
            assert!(
                molecule.atoms[atom_id].is_aromatic,
                "Ring atom {atom_id} should be aromatic"
            );
            assert_eq!(
                molecule.atoms[atom_id].hybridization,
                Hybridization::Resonant,
                "Ring atom {atom_id} should be classified as resonant"
            );
        }
    }
}
