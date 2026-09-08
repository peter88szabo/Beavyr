use super::super::{AtomBlueprint, InputBondBlueprint, MoleculeTestCase, OutputBondBlueprint};
use crate::forcefield::dreiding::typer::{Element, GraphBondOrder, TopologyBondOrder};

pub const URACIL: MoleculeTestCase = MoleculeTestCase {
    name: "Uracil",
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
            label: "O2",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O4",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
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
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "N3",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "N3",
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
            atom1_label: "C2",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "O4",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N3",
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
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N3",
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
            atom1_label: "C2",
            atom2_label: "O2",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "O4",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "N3",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N3",
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

pub const THYMINE: MoleculeTestCase = MoleculeTestCase {
    name: "Thymine",
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
            label: "O2",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O4",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "C7",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
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
            atom2_label: "N3",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "N3",
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
            atom1_label: "C2",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "O4",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "C5",
            atom2_label: "C7",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N3",
            atom2_label: "H3",
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
            atom1_label: "C5",
            atom2_label: "C7",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N3",
            atom2_label: "H3",
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
            atom1_label: "C2",
            atom2_label: "O2",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "O4",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "N3",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N3",
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

pub const CYTOSINE: MoleculeTestCase = MoleculeTestCase {
    name: "Cytosine",
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
            label: "O2",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "N4",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H41",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H42",
            element: Element::H,
            expected_type: "H_HB",
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
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "N3",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "N3",
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
            atom1_label: "C2",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "N4",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N4",
            atom2_label: "H41",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N4",
            atom2_label: "H42",
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
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C4",
            atom2_label: "N4",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N4",
            atom2_label: "H41",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N4",
            atom2_label: "H42",
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
            atom1_label: "C2",
            atom2_label: "O2",
            order: TopologyBondOrder::Double,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "C2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "N3",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N3",
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

pub const ADENINE: MoleculeTestCase = MoleculeTestCase {
    name: "Adenine",
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
            label: "H9",
            element: Element::H,
            expected_type: "H_HB",
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
            atom1_label: "N9",
            atom2_label: "H9",
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
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "N6",
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
            atom1_label: "N9",
            atom2_label: "H9",
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
    ],
};

pub const GUANINE: MoleculeTestCase = MoleculeTestCase {
    name: "Guanine",
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
            label: "O6",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "N2",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H8",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H9",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H21",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H22",
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
            atom2_label: "O6",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "N2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H8",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N9",
            atom2_label: "H9",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N2",
            atom2_label: "H21",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N2",
            atom2_label: "H22",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "N2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C8",
            atom2_label: "H8",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N9",
            atom2_label: "H9",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N2",
            atom2_label: "H21",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N2",
            atom2_label: "H22",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C6",
            atom2_label: "O6",
            order: TopologyBondOrder::Double,
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
    ],
};

pub const DEOXYADENOSINE: MoleculeTestCase = MoleculeTestCase {
    name: "Deoxyadenosine",
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
            label: "C5'",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "O4'",
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
            label: "H1'",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H2'1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H2'2",
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
            atom2_label: "H2'1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2'",
            atom2_label: "H2'2",
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
        OutputBondBlueprint {
            atom1_label: "N9",
            atom2_label: "C1'",
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
            atom2_label: "H2'1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2'",
            atom2_label: "H2'2",
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
    ],
};

pub const DINUCLEOTIDE_BACKBONE: MoleculeTestCase = MoleculeTestCase {
    name: "Dinucleotide Backbone Fragment",
    atoms: &[
        AtomBlueprint {
            label: "C_up",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "O3'_up",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "P",
            element: Element::P,
            expected_type: "P_3",
        },
        AtomBlueprint {
            label: "O1P",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O2P",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "O5'_down",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "C_down",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H_up1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_up2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_up3",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_down1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_down2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "H_down3",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "C_up",
            atom2_label: "O3'_up",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "O3'_up",
            atom2_label: "P",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O5'_down",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "O5'_down",
            atom2_label: "C_down",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O1P",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O2P",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_up",
            atom2_label: "H_up1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_up",
            atom2_label: "H_up2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_up",
            atom2_label: "H_up3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_down",
            atom2_label: "H_down1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_down",
            atom2_label: "H_down2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C_down",
            atom2_label: "H_down3",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "C_up",
            atom2_label: "O3'_up",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "O3'_up",
            atom2_label: "P",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O5'_down",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "O5'_down",
            atom2_label: "C_down",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_up",
            atom2_label: "H_up1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_up",
            atom2_label: "H_up2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_up",
            atom2_label: "H_up3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_down",
            atom2_label: "H_down1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_down",
            atom2_label: "H_down2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C_down",
            atom2_label: "H_down3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O1P",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "P",
            atom2_label: "O2P",
            order: TopologyBondOrder::Resonant,
        },
    ],
};
