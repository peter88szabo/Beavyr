//! Ported upstream test suite for the vendored typer.
//!
//! These are `dreid-typer`'s own tests, carried across the vendoring essentially unchanged. They
//! are what makes "validated" a fact rather than an assumption here -- in particular they are the
//! proof that translating upstream's `default.rules.toml` into the Rust deck in
//! `typing::rules_data` preserved it, since a dropped or mangled rule shows up immediately as a
//! wrong atom type on one of the molecules from the DREIDING paper.
//!
//! Changes from upstream: `crate::forcefield::dreiding::typer::` paths point at this module tree instead, and the
//! separate `tests/integration_tests.rs` binary is folded in at the bottom (Beavyr is a binary
//! crate, so it has no integration-test target).

pub mod cases;

use crate::forcefield::dreiding::typer::{
    Element, GraphBondOrder, MolecularGraph, MolecularTopology, TopologyBondOrder, assign_topology,
};
use std::collections::{HashMap, HashSet};

#[derive(Debug)]
pub struct AtomBlueprint {
    pub label: &'static str,
    pub element: Element,
    pub expected_type: &'static str,
}

#[derive(Debug)]
pub struct InputBondBlueprint {
    pub atom1_label: &'static str,
    pub atom2_label: &'static str,
    pub order: GraphBondOrder,
}

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct OutputBondBlueprint {
    pub atom1_label: &'static str,
    pub atom2_label: &'static str,
    pub order: TopologyBondOrder,
}

#[derive(Debug)]
pub struct MoleculeTestCase {
    pub name: &'static str,
    pub atoms: &'static [AtomBlueprint],
    pub bonds: &'static [InputBondBlueprint],
    pub expected_bonds: &'static [OutputBondBlueprint],
}

pub struct LabeledMolecule {
    graph: MolecularGraph,
    labels: HashMap<&'static str, usize>,
}

impl LabeledMolecule {
    pub fn graph(&self) -> &MolecularGraph {
        &self.graph
    }

    pub fn id(&self, label: &'static str) -> usize {
        *self
            .labels
            .get(label)
            .unwrap_or_else(|| panic!("Unknown atom label: {}", label))
    }
}

pub fn run_molecule_test_case(case: &MoleculeTestCase) {
    let molecule = build_from_blueprint(case);

    let topology = assign_topology(molecule.graph())
        .unwrap_or_else(|err| panic!("Topology assignment failed for '{}': {:?}", case.name, err));

    verify_atom_types(&topology, &molecule, case);
    verify_bond_orders(&topology, &molecule, case);
}

fn build_from_blueprint(case: &MoleculeTestCase) -> LabeledMolecule {
    let mut graph = MolecularGraph::new();
    let mut labels = HashMap::new();

    for atom_bp in case.atoms {
        let id = graph.add_atom(atom_bp.element);
        if labels.insert(atom_bp.label, id).is_some() {
            panic!(
                "Molecule '{}': Duplicate atom label '{}'",
                case.name, atom_bp.label
            );
        }
    }

    for bond_bp in case.bonds {
        let id1 = *labels
            .get(bond_bp.atom1_label)
            .unwrap_or_else(|| panic!("Label '{}' not found", bond_bp.atom1_label));
        let id2 = *labels
            .get(bond_bp.atom2_label)
            .unwrap_or_else(|| panic!("Label '{}' not found", bond_bp.atom2_label));
        graph.add_bond(id1, id2, bond_bp.order).unwrap();
    }

    LabeledMolecule { graph, labels }
}

fn verify_atom_types(
    topology: &MolecularTopology,
    molecule: &LabeledMolecule,
    case: &MoleculeTestCase,
) {
    let actual_types: HashMap<usize, &str> = topology
        .atoms
        .iter()
        .map(|a| (a.id, a.atom_type.as_str()))
        .collect();

    let mut all_heavy_atoms_tested = true;

    for atom_bp in case.atoms {
        let atom_id = molecule.id(atom_bp.label);
        let actual_type = actual_types.get(&atom_id).unwrap_or(&"UNTYPED");

        assert_eq!(
            *actual_type, atom_bp.expected_type,
            "\n --- Test Failure ---\nMolecule: '{}'\nAtom Label: '{}'\n  Expected Type: '{}'\n  Actual Type:   '{}'\n -------------------- \n",
            case.name, atom_bp.label, atom_bp.expected_type, actual_type
        );
    }

    for atom in &topology.atoms {
        if atom.element != Element::H
            && !case.atoms.iter().any(|ab| molecule.id(ab.label) == atom.id)
        {
            all_heavy_atoms_tested = false;
            eprintln!(
                "Warning: Heavy atom with ID {} was not checked in test case '{}'",
                atom.id, case.name
            );
        }
    }
    assert!(
        all_heavy_atoms_tested,
        "One or more heavy atoms were not defined in the test case blueprint for '{}'",
        case.name
    );
}

fn verify_bond_orders(
    topology: &MolecularTopology,
    molecule: &LabeledMolecule,
    case: &MoleculeTestCase,
) {
    let expected_bonds: HashSet<_> = case
        .expected_bonds
        .iter()
        .map(|bp| {
            let mut ids = [molecule.id(bp.atom1_label), molecule.id(bp.atom2_label)];
            ids.sort_unstable();
            ((ids[0], ids[1]), bp.order)
        })
        .collect();

    let actual_bonds: HashSet<_> = topology
        .bonds
        .iter()
        .map(|bond| (bond.atom_ids, bond.order))
        .collect();

    assert_eq!(
        actual_bonds, expected_bonds,
        "\n --- Test Failure ---\nMolecule: '{}'\nBond order mismatch.\n  Expected: {:?}\n  Actual:   {:?}\n -------------------- \n",
        case.name, expected_bonds, actual_bonds
    );
}


#[cfg(test)]
mod generated {
    use super::cases;


    use cases::amino_acids::*;
    use cases::dreiding_paper::*;
    use cases::nucleic_acids::*;
    use super::run_molecule_test_case;

    macro_rules! generate_molecule_test {
        ($test_name:ident, $molecule_case:expr) => {
            #[test]
            fn $test_name() {
                run_molecule_test_case(&$molecule_case);
            }
        };
    }

    generate_molecule_test!(glycine_zwitterion_is_typed_correctly, GLYCINE_ZWITTERION);
    generate_molecule_test!(alanine_zwitterion_is_typed_correctly, ALANINE_ZWITTERION);
    generate_molecule_test!(valine_zwitterion_is_typed_correctly, VALINE_ZWITTERION);
    generate_molecule_test!(leucine_zwitterion_is_typed_correctly, LEUCINE_ZWITTERION);
    generate_molecule_test!(
        isoleucine_zwitterion_is_typed_correctly,
        ISOLEUCINE_ZWITTERION
    );
    generate_molecule_test!(proline_zwitterion_is_typed_correctly, PROLINE_ZWITTERION);
    generate_molecule_test!(serine_zwitterion_is_typed_correctly, SERINE_ZWITTERION);
    generate_molecule_test!(
        threonine_zwitterion_is_typed_correctly,
        THREONINE_ZWITTERION
    );
    generate_molecule_test!(cysteine_zwitterion_is_typed_correctly, CYSTEINE_ZWITTERION);
    generate_molecule_test!(
        methionine_zwitterion_is_typed_correctly,
        METHIONINE_ZWITTERION
    );
    generate_molecule_test!(
        aspartate_zwitterion_is_typed_correctly,
        ASPARTATE_ZWITTERION
    );
    generate_molecule_test!(
        asparagine_zwitterion_is_typed_correctly,
        ASPARAGINE_ZWITTERION
    );
    generate_molecule_test!(
        glutamate_zwitterion_is_typed_correctly,
        GLUTAMATE_ZWITTERION
    );
    generate_molecule_test!(
        glutamine_zwitterion_is_typed_correctly,
        GLUTAMINE_ZWITTERION
    );
    generate_molecule_test!(lysine_zwitterion_is_typed_correctly, LYSINE_ZWITTERION);
    generate_molecule_test!(arginine_zwitterion_is_typed_correctly, ARGININE_ZWITTERION);
    generate_molecule_test!(
        phenylalanine_zwitterion_is_typed_correctly,
        PHENYLALANINE_ZWITTERION
    );
    generate_molecule_test!(tyrosine_zwitterion_is_typed_correctly, TYROSINE_ZWITTERION);
    generate_molecule_test!(
        histidine_zwitterion_is_typed_correctly,
        HISTIDINE_ZWITTERION
    );
    generate_molecule_test!(
        tryptophan_zwitterion_is_typed_correctly,
        TRYPTOPHAN_ZWITTERION
    );
    generate_molecule_test!(diglycine_peptide_is_typed_correctly, DIGLYCINE);

    generate_molecule_test!(uracil_is_typed_correctly, URACIL);
    generate_molecule_test!(thymine_is_typed_correctly, THYMINE);
    generate_molecule_test!(cytosine_is_typed_correctly, CYTOSINE);
    generate_molecule_test!(adenine_is_typed_correctly, ADENINE);
    generate_molecule_test!(guanine_is_typed_correctly, GUANINE);
    generate_molecule_test!(deoxyadenosine_is_typed_correctly, DEOXYADENOSINE);
    generate_molecule_test!(
        dinucleotide_backbone_fragment_is_typed_correctly,
        DINUCLEOTIDE_BACKBONE
    );

    generate_molecule_test!(cyclobutane_lactone_is_typed_correctly, ACBUOL);
    generate_molecule_test!(cage_ether_is_typed_correctly, AFURPO10);
    generate_molecule_test!(adamantane_is_typed_correctly, ADAMANTANE);
    generate_molecule_test!(decalin_is_typed_correctly, DECALIN);
    generate_molecule_test!(dithiacyclohexene_is_typed_correctly, ABAXES);
    generate_molecule_test!(methyl_pyridine_is_typed_correctly, METHYL_PYRIDINE);
    generate_molecule_test!(acridine_is_typed_correctly, ACRIDINE);
    generate_molecule_test!(phosphinane_is_typed_correctly, PHOSPHINANE);
    generate_molecule_test!(trinitrobenzene_is_typed_correctly, TRINITROBENZENE);
    generate_molecule_test!(tetramethylthiourea_is_typed_correctly, TETRAMETHYLTHIOUREA);
    generate_molecule_test!(dimethyl_sulfoxide_is_typed_correctly, DIMETHYL_SULFOXIDE);
    generate_molecule_test!(methanesulfonamide_is_typed_correctly, METHANESULFONAMIDE);
    generate_molecule_test!(
        chloroacetyl_chloride_is_typed_correctly,
        CHLOROACETYL_CHLORIDE
    );
    generate_molecule_test!(adenosine_is_typed_correctly, ADENOSINE);
    generate_molecule_test!(phosphate_ester_is_typed_correctly, PHOSPHATE_ESTER);
    generate_molecule_test!(choline_cation_is_typed_correctly, CHOLINE_CATION);
    generate_molecule_test!(perchlorate_anion_is_typed_correctly, PERCHLORATE_ANION);

}
