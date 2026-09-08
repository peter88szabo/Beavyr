use super::super::{AtomBlueprint, InputBondBlueprint, MoleculeTestCase, OutputBondBlueprint};
use crate::forcefield::dreiding::typer::{Element, GraphBondOrder, TopologyBondOrder};

pub const ACBUOL: MoleculeTestCase = MoleculeTestCase {
    name: "ACBUOL - Cyclobutane Lactone",
    atoms: &[
        AtomBlueprint {
            label: "C1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C3",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C4",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C5",
            element: Element::C,
            expected_type: "C_2",
        },
        AtomBlueprint {
            label: "O1",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O2",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "C6",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "O3",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "C7",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "O4",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "H_C1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C3",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C4",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C6a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C6b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C7a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C7b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C7c",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_O4",
            element: Element::H,
            expected_type: "H_HB",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C5",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "O3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "O4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C6",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "O1",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "O2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "O2",
            atom2_label: "C6",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "O3",
            atom2_label: "C7",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H_C7a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H_C7b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H_C7c",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "O4",
            atom2_label: "H_O4",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C5",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "O3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "O4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C6",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "O1",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "O2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "O2",
            atom2_label: "C6",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "O3",
            atom2_label: "C7",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H_C7a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H_C7b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H_C7c",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "O4",
            atom2_label: "H_O4",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const AFURPO10: MoleculeTestCase = MoleculeTestCase {
    name: "AFURPO10 - Cage Ether",
    atoms: &[
        AtomBlueprint {
            label: "C1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C3",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C4",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C5",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C6",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "O1",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "O2",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "H_C1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C3",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C4",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C5",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C6a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C6b",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C6",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "O2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "O2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H_C5",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6b",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C6",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "O1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "O2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "O1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "O2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H_C5",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6b",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const ADAMANTANE: MoleculeTestCase = MoleculeTestCase {
    name: "Adamantane",
    atoms: &[
        AtomBlueprint {
            label: "C1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C3",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C4",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C5",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C6",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C7",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C8",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C9",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C10",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H_C1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C3",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C5",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C7",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C4a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C4b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C6a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C6b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C8a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C8b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C9a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C9b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C10a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C10b",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C6",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C8",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C10",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C8",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C9",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "C6",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "C9",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "C10",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H_C5",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H_C7",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H_C8a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H_C8b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C9",
            atom2_label: "H_C9a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C9",
            atom2_label: "H_C9b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C10",
            atom2_label: "H_C10a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C10",
            atom2_label: "H_C10b",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C6",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C8",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C10",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C8",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C9",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "C6",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "C9",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "C10",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H_C5",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H_C7",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H_C8a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H_C8b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C9",
            atom2_label: "H_C9a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C9",
            atom2_label: "H_C9b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C10",
            atom2_label: "H_C10a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C10",
            atom2_label: "H_C10b",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const DECALIN: MoleculeTestCase = MoleculeTestCase {
    name: "Decalin",
    atoms: &[
        AtomBlueprint {
            label: "C1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C3",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C4",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C5",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C6",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C7",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C8",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C9_fused",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C10_fused",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H_C1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C3a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C3b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C4a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C4b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C5a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C5b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C6",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C7a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C7b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C8a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C8b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C9a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C9b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C10a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C10b",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "C1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "C7",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "C8",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "C9_fused",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C9_fused",
            atom2_label: "C10_fused",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C10_fused",
            atom2_label: "C1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H_C5a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H_C5b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H_C7a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H_C7b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H_C8a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H_C8b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C9_fused",
            atom2_label: "H_C9a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C9_fused",
            atom2_label: "H_C9b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C10_fused",
            atom2_label: "H_C10a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C10_fused",
            atom2_label: "H_C10b",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "C1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "C7",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "C8",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "C9_fused",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C9_fused",
            atom2_label: "C10_fused",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C10_fused",
            atom2_label: "C1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H_C5a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H_C5b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H_C7a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H_C7b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H_C8a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H_C8b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C9_fused",
            atom2_label: "H_C9a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C9_fused",
            atom2_label: "H_C9b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C10_fused",
            atom2_label: "H_C10a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C10_fused",
            atom2_label: "H_C10b",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const ABAXES: MoleculeTestCase = MoleculeTestCase {
    name: "ABAXES - 1,4-dithiacyclohex-2-ene",
    atoms: &[
        AtomBlueprint {
            label: "S1",
            element: Element::S,
            expected_type: "S_3",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_2",
        },
        AtomBlueprint {
            label: "C3",
            element: Element::C,
            expected_type: "C_2",
        },
        AtomBlueprint {
            label: "S4",
            element: Element::S,
            expected_type: "S_3",
        },
        AtomBlueprint {
            label: "C5",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C6",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H51",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H52",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H61",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H62",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "S1",
            atom2_label: "C2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "S4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "S4",
            atom2_label: "C5",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "S1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H51",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H52",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H61",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H62",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "S1",
            atom2_label: "C2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "S4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "S4",
            atom2_label: "C5",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "S1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H51",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H52",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H61",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H62",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const METHYL_PYRIDINE: MoleculeTestCase = MoleculeTestCase {
    name: "ACPYNS - 4-Methyl-Pyridine",
    atoms: &[
        AtomBlueprint {
            label: "N1",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C3",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C4",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C5",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C6",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C7",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H5",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H6",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H71",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H72",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H73",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "N1",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C7",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H5",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H6",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H71",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H72",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H73",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C7",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H5",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H6",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H71",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H72",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H73",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "N1",
            order: TopologyBondOrder::Resonant,
        },
    ],
};

pub const ACRIDINE: MoleculeTestCase = MoleculeTestCase {
    name: "ACRAMS - Acridine",
    atoms: &[
        AtomBlueprint {
            label: "C1",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C3",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C4",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C5",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C6",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C7",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C8",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C9",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "N10",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "C11",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C12",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C13",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C14",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H4",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H7",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H8",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H11",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H12",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H13",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H14",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "C1",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C7",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "C8",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "C9",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C9",
            atom2_label: "N10",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "N10",
            atom2_label: "C6",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "N10",
            atom2_label: "C11",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C11",
            atom2_label: "C12",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C12",
            atom2_label: "C13",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C13",
            atom2_label: "C14",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C14",
            atom2_label: "C9",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H7",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H8",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C11",
            atom2_label: "H11",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C12",
            atom2_label: "H12",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C13",
            atom2_label: "H13",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C14",
            atom2_label: "H14",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "H7",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H8",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C11",
            atom2_label: "H11",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C12",
            atom2_label: "H12",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C13",
            atom2_label: "H13",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C14",
            atom2_label: "H14",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "C1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C7",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C7",
            atom2_label: "C8",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "C9",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C9",
            atom2_label: "N10",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N10",
            atom2_label: "C6",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N10",
            atom2_label: "C11",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C11",
            atom2_label: "C12",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C12",
            atom2_label: "C13",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C13",
            atom2_label: "C14",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C14",
            atom2_label: "C9",
            order: TopologyBondOrder::Resonant,
        },
    ],
};

pub const PHOSPHINANE: MoleculeTestCase = MoleculeTestCase {
    name: "AFINDS - Phosphinane",
    atoms: &[
        AtomBlueprint {
            label: "P1",
            element: Element::P,
            expected_type: "P_3",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C3",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C4",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C5",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C6",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H_P1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C3a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C3b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C4a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C4b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C5a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C5b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C6a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C6b",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "P1",
            atom2_label: "C2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "P1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "P1",
            atom2_label: "H_P1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H_C5a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H_C5b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6b",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "P1",
            atom2_label: "C2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "P1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "P1",
            atom2_label: "H_P1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "H_C3b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H_C5a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "H_C5b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6b",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const TRINITROBENZENE: MoleculeTestCase = MoleculeTestCase {
    name: "ACNPEC - 1,3,5-Trinitrobenzene",
    atoms: &[
        AtomBlueprint {
            label: "C1",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C3",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C4",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C5",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C6",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "N1",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "O1A",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O1B",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "N3",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "O3A",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O3B",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "N5",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "O5A",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O5B",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "H_C2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C4",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C6",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "C1",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "N1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "N3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "N5",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "O1A",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "O1B",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "N3",
            atom2_label: "O3A",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N3",
            atom2_label: "O3B",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "N5",
            atom2_label: "O5A",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N5",
            atom2_label: "O5B",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "C2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "C3",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "C4",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "C5",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C6",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "C1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "O1A",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "O1B",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N3",
            atom2_label: "O3A",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N3",
            atom2_label: "O3B",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N5",
            atom2_label: "O5A",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N5",
            atom2_label: "O5B",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "N1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3",
            atom2_label: "N3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "N5",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "H_C4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "H_C6",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const TETRAMETHYLTHIOUREA: MoleculeTestCase = MoleculeTestCase {
    name: "ACPRET03 - Tetramethylthiourea",
    atoms: &[
        AtomBlueprint {
            label: "CS",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "S",
            element: Element::S,
            expected_type: "S_2",
        },
        AtomBlueprint {
            label: "N1",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "N2",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "C11",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C12",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C21",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C22",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H_C11a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C11b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C11c",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C12a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C12b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C12c",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C21a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C21b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C21c",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C22a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C22b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C22c",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "CS",
            atom2_label: "S",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CS",
            atom2_label: "N1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CS",
            atom2_label: "N2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C11",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C12",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N2",
            atom2_label: "C21",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N2",
            atom2_label: "C22",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C11",
            atom2_label: "H_C11a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C11",
            atom2_label: "H_C11b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C11",
            atom2_label: "H_C11c",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C12",
            atom2_label: "H_C12a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C12",
            atom2_label: "H_C12b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C12",
            atom2_label: "H_C12c",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C21",
            atom2_label: "H_C21a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C21",
            atom2_label: "H_C21b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C21",
            atom2_label: "H_C21c",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C22",
            atom2_label: "H_C22a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C22",
            atom2_label: "H_C22b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C22",
            atom2_label: "H_C22c",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "CS",
            atom2_label: "S",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "CS",
            atom2_label: "N1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CS",
            atom2_label: "N2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C11",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C12",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N2",
            atom2_label: "C21",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N2",
            atom2_label: "C22",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C11",
            atom2_label: "H_C11a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C11",
            atom2_label: "H_C11b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C11",
            atom2_label: "H_C11c",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C12",
            atom2_label: "H_C12a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C12",
            atom2_label: "H_C12b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C12",
            atom2_label: "H_C12c",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C21",
            atom2_label: "H_C21a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C21",
            atom2_label: "H_C21b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C21",
            atom2_label: "H_C21c",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C22",
            atom2_label: "H_C22a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C22",
            atom2_label: "H_C22b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C22",
            atom2_label: "H_C22c",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const DIMETHYL_SULFOXIDE: MoleculeTestCase = MoleculeTestCase {
    name: "ACSESO10 - Dimethyl Sulfoxide",
    atoms: &[
        AtomBlueprint {
            label: "S",
            element: Element::S,
            expected_type: "S_3",
        },
        AtomBlueprint {
            label: "O",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "C1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H_C1a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C1b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C1c",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2c",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "S",
            atom2_label: "O",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "S",
            atom2_label: "C1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "S",
            atom2_label: "C2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1c",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2c",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "S",
            atom2_label: "O",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "S",
            atom2_label: "C1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "S",
            atom2_label: "C2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1c",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2c",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const METHANESULFONAMIDE: MoleculeTestCase = MoleculeTestCase {
    name: "ACISUL - Methanesulfonamide",
    atoms: &[
        AtomBlueprint {
            label: "S",
            element: Element::S,
            expected_type: "S_3",
        },
        AtomBlueprint {
            label: "O1",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O2",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H_Na",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H_Nb",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H_Ca",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_Cb",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_Cc",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "S",
            atom2_label: "O1",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "S",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "S",
            atom2_label: "N",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "S",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H_Na",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H_Nb",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "H_Ca",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "H_Cb",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "H_Cc",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "S",
            atom2_label: "O1",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "S",
            atom2_label: "O2",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "S",
            atom2_label: "N",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "S",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H_Na",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H_Nb",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "H_Ca",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "H_Cb",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "H_Cc",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const CHLOROACETYL_CHLORIDE: MoleculeTestCase = MoleculeTestCase {
    name: "ACHNAP10 - Chloroacetyl chloride",
    atoms: &[
        AtomBlueprint {
            label: "C_alpha",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "Clalpha",
            element: Element::Cl,
            expected_type: "Cl",
        },
        AtomBlueprint {
            label: "C_carbonyl",
            element: Element::C,
            expected_type: "C_2",
        },
        AtomBlueprint {
            label: "O_carbonyl",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "Clcarbonyl",
            element: Element::Cl,
            expected_type: "Cl",
        },
        AtomBlueprint {
            label: "H_alpha1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_alpha2",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "C_alpha",
            atom2_label: "Clalpha",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_alpha",
            atom2_label: "C_carbonyl",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_carbonyl",
            atom2_label: "O_carbonyl",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "C_carbonyl",
            atom2_label: "Clcarbonyl",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_alpha",
            atom2_label: "H_alpha1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_alpha",
            atom2_label: "H_alpha2",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C_alpha",
            atom2_label: "Clalpha",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_alpha",
            atom2_label: "C_carbonyl",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_carbonyl",
            atom2_label: "O_carbonyl",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "C_carbonyl",
            atom2_label: "Clcarbonyl",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_alpha",
            atom2_label: "H_alpha1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_alpha",
            atom2_label: "H_alpha2",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const ADENOSINE: MoleculeTestCase = MoleculeTestCase {
    name: "ADENOS10 - Adenosine",
    atoms: &[
        AtomBlueprint {
            label: "N1",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "N3",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "C4",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C5",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "C6",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "N7",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "C8",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "N9",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "N6",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "C1'",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C2'",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C3'",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C4'",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "O4'",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "C5'",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "O2'",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "O3'",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "O5'",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H8",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H61",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H62",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H1'",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H2'",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H3'",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H4'",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H5'1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H5'2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_O2'",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H_O3'",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H_O5'",
            element: Element::H,
            expected_type: "H_HB",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N9",
            atom2_label: "C4",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "N3",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "N3",
            atom2_label: "C2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "N1",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C6",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "C5",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C4",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "N7",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "N7",
            atom2_label: "C8",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "N9",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "N6",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H8",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N6",
            atom2_label: "H61",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N6",
            atom2_label: "H62",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N9",
            atom2_label: "C1'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1'",
            atom2_label: "C2'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2'",
            atom2_label: "C3'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3'",
            atom2_label: "C4'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4'",
            atom2_label: "O4'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "O4'",
            atom2_label: "C1'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2'",
            atom2_label: "O2'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "O2'",
            atom2_label: "H_O2'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3'",
            atom2_label: "O3'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "O3'",
            atom2_label: "H_O3'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4'",
            atom2_label: "C5'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5'",
            atom2_label: "O5'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "O5'",
            atom2_label: "H_O5'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1'",
            atom2_label: "H1'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2'",
            atom2_label: "H2'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C3'",
            atom2_label: "H3'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C4'",
            atom2_label: "H4'",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5'",
            atom2_label: "H5'1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C5'",
            atom2_label: "H5'2",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N9",
            atom2_label: "C1'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1'",
            atom2_label: "C2'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2'",
            atom2_label: "C3'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3'",
            atom2_label: "C4'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4'",
            atom2_label: "O4'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "O4'",
            atom2_label: "C1'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2'",
            atom2_label: "O2'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "O2'",
            atom2_label: "H_O2'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3'",
            atom2_label: "O3'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "O3'",
            atom2_label: "H_O3'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4'",
            atom2_label: "C5'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5'",
            atom2_label: "O5'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "O5'",
            atom2_label: "H_O5'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1'",
            atom2_label: "H1'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2'",
            atom2_label: "H2'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C3'",
            atom2_label: "H3'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C4'",
            atom2_label: "H4'",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5'",
            atom2_label: "H5'1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C5'",
            atom2_label: "H5'2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H8",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N6",
            atom2_label: "H61",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N6",
            atom2_label: "H62",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N9",
            atom2_label: "C4",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "N3",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N3",
            atom2_label: "C2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "N1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C6",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "C5",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C4",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "N7",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N7",
            atom2_label: "C8",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "N9",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "N6",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const PHOSPHATE_ESTER: MoleculeTestCase = MoleculeTestCase {
    name: "AEBDOD10 - Phosphate Ester",
    atoms: &[
        AtomBlueprint {
            label: "P",
            element: Element::P,
            expected_type: "P_3",
        },
        AtomBlueprint {
            label: "O_double",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O_single_neg",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O_bridge1",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "O_bridge2",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "C1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H_C1a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C1b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C1c",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_C2c",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O_double",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O_single_neg",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O_bridge1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O_bridge2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "O_bridge1",
            atom2_label: "C1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "O_bridge2",
            atom2_label: "C2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1c",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2c",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O_double",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O_single_neg",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O_bridge1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O_bridge2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "O_bridge1",
            atom2_label: "C1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "O_bridge2",
            atom2_label: "C2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "H_C1c",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "H_C2c",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const CHOLINE_CATION: MoleculeTestCase = MoleculeTestCase {
    name: "ACHTAR10 - Choline Cation",
    atoms: &[
        AtomBlueprint {
            label: "N_quat",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "C_Me1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C_Me2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C_Me3",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C_ethyl1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C_ethyl2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "O_hydroxyl",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "H_Me1a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_Me1b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_Me1c",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_Me2a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_Me2b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_Me2c",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_Me3a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_Me3b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_Me3c",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_CE1a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_CE1b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_CE2a",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_CE2b",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_OH",
            element: Element::H,
            expected_type: "H_HB",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N_quat",
            atom2_label: "C_Me1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N_quat",
            atom2_label: "C_Me2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N_quat",
            atom2_label: "C_Me3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N_quat",
            atom2_label: "C_ethyl1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_Me1",
            atom2_label: "H_Me1a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_Me1",
            atom2_label: "H_Me1b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_Me1",
            atom2_label: "H_Me1c",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_Me2",
            atom2_label: "H_Me2a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_Me2",
            atom2_label: "H_Me2b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_Me2",
            atom2_label: "H_Me2c",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_Me3",
            atom2_label: "H_Me3a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_Me3",
            atom2_label: "H_Me3b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_Me3",
            atom2_label: "H_Me3c",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_ethyl1",
            atom2_label: "C_ethyl2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_ethyl1",
            atom2_label: "H_CE1a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_ethyl1",
            atom2_label: "H_CE1b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_ethyl2",
            atom2_label: "O_hydroxyl",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_ethyl2",
            atom2_label: "H_CE2a",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_ethyl2",
            atom2_label: "H_CE2b",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "O_hydroxyl",
            atom2_label: "H_OH",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N_quat",
            atom2_label: "C_Me1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N_quat",
            atom2_label: "C_Me2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N_quat",
            atom2_label: "C_Me3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N_quat",
            atom2_label: "C_ethyl1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_Me1",
            atom2_label: "H_Me1a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_Me1",
            atom2_label: "H_Me1b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_Me1",
            atom2_label: "H_Me1c",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_Me2",
            atom2_label: "H_Me2a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_Me2",
            atom2_label: "H_Me2b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_Me2",
            atom2_label: "H_Me2c",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_Me3",
            atom2_label: "H_Me3a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_Me3",
            atom2_label: "H_Me3b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_Me3",
            atom2_label: "H_Me3c",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_ethyl1",
            atom2_label: "C_ethyl2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_ethyl1",
            atom2_label: "H_CE1a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_ethyl1",
            atom2_label: "H_CE1b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_ethyl2",
            atom2_label: "O_hydroxyl",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_ethyl2",
            atom2_label: "H_CE2a",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_ethyl2",
            atom2_label: "H_CE2b",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "O_hydroxyl",
            atom2_label: "H_OH",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const PERCHLORATE_ANION: MoleculeTestCase = MoleculeTestCase {
    name: "ACMBPN - Perchlorate Anion",
    atoms: &[
        AtomBlueprint {
            label: "CL",
            element: Element::Cl,
            expected_type: "Cl",
        },
        AtomBlueprint {
            label: "O1",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O2",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O3",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O4",
            element: Element::O,
            expected_type: "O_3",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "CL",
            atom2_label: "O1",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CL",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CL",
            atom2_label: "O3",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CL",
            atom2_label: "O4",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "CL",
            atom2_label: "O1",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "CL",
            atom2_label: "O2",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "CL",
            atom2_label: "O3",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "CL",
            atom2_label: "O4",
            order: TopologyBondOrder::Single,
        },
    ],
};
