use super::super::{AtomBlueprint, InputBondBlueprint, MoleculeTestCase, OutputBondBlueprint};
use crate::forcefield::dreiding::typer::{Element, GraphBondOrder, TopologyBondOrder};

pub const GLYCINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Glycine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HA2",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA2",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA2",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const ALANINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Alanine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB3",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB3",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB3",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const VALINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Valine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG11",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG12",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG13",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG21",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG22",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG23",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG11",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG12",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG13",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG21",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG22",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG23",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG11",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG12",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG13",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG21",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG22",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG23",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const LEUCINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Leucine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CD1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CD2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD11",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD12",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD13",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD21",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD22",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD23",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD11",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD12",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD13",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD21",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD22",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD23",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD11",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD12",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD13",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD21",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD22",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD23",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const ISOLEUCINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Isoleucine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CD1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG11",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG12",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG21",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG22",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG23",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD11",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD12",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD13",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG1",
            atom2_label: "CD1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG11",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG12",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG21",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG22",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG23",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD11",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD12",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD13",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG1",
            atom2_label: "CD1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG11",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG1",
            atom2_label: "HG12",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG21",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG22",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG23",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD11",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD12",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD13",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const PROLINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Proline Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CD",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD2",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "N",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "HD1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "HD2",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "N",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "HD1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "HD2",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const SERINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Serine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "OG",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG",
            element: Element::H,
            expected_type: "H_HB",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "OG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "OG",
            atom2_label: "HG",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "OG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "OG",
            atom2_label: "HG",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const THREONINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Threonine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "OG1",
            element: Element::O,
            expected_type: "O_3",
        },
        AtomBlueprint {
            label: "CG2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HG21",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG22",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG23",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "OG1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "OG1",
            atom2_label: "HG1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG21",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG22",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG23",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "OG1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "OG1",
            atom2_label: "HG1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG21",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG22",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG2",
            atom2_label: "HG23",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const CYSTEINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Cysteine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "SG",
            element: Element::S,
            expected_type: "S_3",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG",
            element: Element::H,
            expected_type: "H_HB",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "SG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "SG",
            atom2_label: "HG",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "SG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "SG",
            atom2_label: "HG",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const METHIONINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Methionine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "SD",
            element: Element::S,
            expected_type: "S_3",
        },
        AtomBlueprint {
            label: "CE",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HE1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HE2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HE3",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "SD",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "SD",
            atom2_label: "CE",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CE",
            atom2_label: "HE1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CE",
            atom2_label: "HE2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CE",
            atom2_label: "HE3",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "SD",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "SD",
            atom2_label: "CE",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CE",
            atom2_label: "HE1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CE",
            atom2_label: "HE2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CE",
            atom2_label: "HE3",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const ASPARTATE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Aspartate Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "OD1",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "OD2",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "OD1",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "OD2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "OD1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "OD2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const ASPARAGINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Asparagine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "OD1",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "ND2",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD21",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HD22",
            element: Element::H,
            expected_type: "H_HB",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "OD1",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "ND2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "ND2",
            atom2_label: "HD21",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "ND2",
            atom2_label: "HD22",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "OD1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "ND2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "ND2",
            atom2_label: "HD21",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "ND2",
            atom2_label: "HD22",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const GLUTAMATE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Glutamate Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CD",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "OE1",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "OE2",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG2",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "OE1",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "OE2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG2",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "OE1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "OE2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG2",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const GLUTAMINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Glutamine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CD",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "OE1",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "NE2",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HE21",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HE22",
            element: Element::H,
            expected_type: "H_HB",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "OE1",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "NE2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "NE2",
            atom2_label: "HE21",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "NE2",
            atom2_label: "HE22",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "OE1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "NE2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "NE2",
            atom2_label: "HE21",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "NE2",
            atom2_label: "HE22",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const LYSINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Lysine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CD",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CE",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "NZ",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HE1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HE2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HZ1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HZ2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HZ3",
            element: Element::H,
            expected_type: "H_HB",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "CE",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CE",
            atom2_label: "NZ",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "HD1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "HD2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CE",
            atom2_label: "HE1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CE",
            atom2_label: "HE2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "NZ",
            atom2_label: "HZ1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "NZ",
            atom2_label: "HZ2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "NZ",
            atom2_label: "HZ3",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "CE",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CE",
            atom2_label: "NZ",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "HD1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "HD2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CE",
            atom2_label: "HE1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CE",
            atom2_label: "HE2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "NZ",
            atom2_label: "HZ1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "NZ",
            atom2_label: "HZ2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "NZ",
            atom2_label: "HZ3",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const ARGININE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Arginine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CD",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "NE",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "CZ",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "NH1",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "NH2",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HG2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HE",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HH11",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HH12",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HH21",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HH22",
            element: Element::H,
            expected_type: "H_HB",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "NE",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "NE",
            atom2_label: "CZ",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CZ",
            atom2_label: "NH1",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CZ",
            atom2_label: "NH2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "HD1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "HD2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "NE",
            atom2_label: "HE",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "NH1",
            atom2_label: "HH11",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "NH1",
            atom2_label: "HH12",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "NH2",
            atom2_label: "HH21",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "NH2",
            atom2_label: "HH22",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "NE",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "NE",
            atom2_label: "CZ",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CZ",
            atom2_label: "NH1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CZ",
            atom2_label: "NH2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "HG2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "HD1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD",
            atom2_label: "HD2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "NE",
            atom2_label: "HE",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "NH1",
            atom2_label: "HH11",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "NH1",
            atom2_label: "HH12",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "NH2",
            atom2_label: "HH21",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "NH2",
            atom2_label: "HH22",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const PHENYLALANINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Phenylalanine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CD1",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CD2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CE1",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CE2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CZ",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HE1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HE2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HZ",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD1",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "CE1",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CE1",
            atom2_label: "CZ",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CZ",
            atom2_label: "CE2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CE2",
            atom2_label: "CD2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "CG",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CE1",
            atom2_label: "HE1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CE2",
            atom2_label: "HE2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CZ",
            atom2_label: "HZ",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "CE1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CE1",
            atom2_label: "CZ",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CZ",
            atom2_label: "CE2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CE2",
            atom2_label: "CD2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "CG",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CE1",
            atom2_label: "HE1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CE2",
            atom2_label: "HE2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CZ",
            atom2_label: "HZ",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const TYROSINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Tyrosine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CD1",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CD2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CE1",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CE2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CZ",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "OH",
            element: Element::O,
            expected_type: "O_R",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HE1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HE2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HH",
            element: Element::H,
            expected_type: "H_HB",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD1",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "CE1",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CE1",
            atom2_label: "CZ",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CZ",
            atom2_label: "CE2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CE2",
            atom2_label: "CD2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "CG",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CZ",
            atom2_label: "OH",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CE1",
            atom2_label: "HE1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CE2",
            atom2_label: "HE2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "OH",
            atom2_label: "HH",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "CE1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CE1",
            atom2_label: "CZ",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CZ",
            atom2_label: "CE2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CE2",
            atom2_label: "CD2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "CG",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CZ",
            atom2_label: "OH",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CE1",
            atom2_label: "HE1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CE2",
            atom2_label: "HE2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "OH",
            atom2_label: "HH",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const HISTIDINE_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Histidine Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "ND1",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "CE1",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "NE2",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "CD2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HE1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD2",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "ND1",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "ND1",
            atom2_label: "CE1",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CE1",
            atom2_label: "NE2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "NE2",
            atom2_label: "CD2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "CG",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "ND1",
            atom2_label: "HD1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CE1",
            atom2_label: "HE1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD2",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "ND1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "ND1",
            atom2_label: "CE1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CE1",
            atom2_label: "NE2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "NE2",
            atom2_label: "CD2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "CG",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "ND1",
            atom2_label: "HD1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CE1",
            atom2_label: "HE1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "HD2",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const TRYPTOPHAN_ZWITTERION: MoleculeTestCase = MoleculeTestCase {
    name: "Tryptophan Zwitterion",
    atoms: &[
        AtomBlueprint {
            label: "N",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C",
            element: Element::C,
            expected_type: "C_R",
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
            label: "CB",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "CG",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CD1",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "NE1",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "CE2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CD2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CE3",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CZ3",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CH2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "CZ2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HB2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HD1",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HE1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HE3",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HZ3",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HH2",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HZ2",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD1",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "NE1",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "NE1",
            atom2_label: "CE2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CE2",
            atom2_label: "CD2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "CG",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "CE3",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CE3",
            atom2_label: "CZ3",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CZ3",
            atom2_label: "CH2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CH2",
            atom2_label: "CZ2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "CZ2",
            atom2_label: "CE2",
            order: GraphBondOrder::Aromatic,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "NE1",
            atom2_label: "HE1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CE3",
            atom2_label: "HE3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CZ3",
            atom2_label: "HZ3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CH2",
            atom2_label: "HH2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CZ2",
            atom2_label: "HZ2",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "CA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "C",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "CB",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "CG",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CG",
            atom2_label: "CD1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "NE1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "NE1",
            atom2_label: "CE2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CE2",
            atom2_label: "CD2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "CG",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CD2",
            atom2_label: "CE3",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CE3",
            atom2_label: "CZ3",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CZ3",
            atom2_label: "CH2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CH2",
            atom2_label: "CZ2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "CZ2",
            atom2_label: "CE2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA",
            atom2_label: "HA",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CB",
            atom2_label: "HB2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CD1",
            atom2_label: "HD1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "NE1",
            atom2_label: "HE1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CE3",
            atom2_label: "HE3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CZ3",
            atom2_label: "HZ3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CH2",
            atom2_label: "HH2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CZ2",
            atom2_label: "HZ2",
            order: TopologyBondOrder::Single,
        },
    ],
};

pub const DIGLYCINE: MoleculeTestCase = MoleculeTestCase {
    name: "Diglycine (Gly-Gly)",
    atoms: &[
        AtomBlueprint {
            label: "N1",
            element: Element::N,
            expected_type: "N_3",
        },
        AtomBlueprint {
            label: "CA1",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C1",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "O1",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "N2",
            element: Element::N,
            expected_type: "N_R",
        },
        AtomBlueprint {
            label: "CA2",
            element: Element::C,
            expected_type: "C_3",
        },
        AtomBlueprint {
            label: "C2",
            element: Element::C,
            expected_type: "C_R",
        },
        AtomBlueprint {
            label: "O2",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "OXT",
            element: Element::O,
            expected_type: "O_2",
        },
        AtomBlueprint {
            label: "H1",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "H3",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA11",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HA12",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HN2",
            element: Element::H,
            expected_type: "H_HB",
        },
        AtomBlueprint {
            label: "HA21",
            element: Element::H,
            expected_type: "H_",
        },
        AtomBlueprint {
            label: "HA22",
            element: Element::H,
            expected_type: "H_",
        },
    ],
    bonds: &[
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "CA1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA1",
            atom2_label: "C1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "O1",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "N2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N2",
            atom2_label: "CA2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA2",
            atom2_label: "C2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "O2",
            order: GraphBondOrder::Double,
        },
        InputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "OXT",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H1",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H3",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA1",
            atom2_label: "HA11",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA1",
            atom2_label: "HA12",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "N2",
            atom2_label: "HN2",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA2",
            atom2_label: "HA21",
            order: GraphBondOrder::Single,
        },
        InputBondBlueprint {
            atom1_label: "CA2",
            atom2_label: "HA22",
            order: GraphBondOrder::Single,
        },
    ],
    expected_bonds: &[
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "CA1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA1",
            atom2_label: "C1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "O1",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C1",
            atom2_label: "N2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N2",
            atom2_label: "CA2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA2",
            atom2_label: "C2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "O2",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "C2",
            atom2_label: "OXT",
            order: TopologyBondOrder::Resonant,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H1",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N1",
            atom2_label: "H3",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA1",
            atom2_label: "HA11",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA1",
            atom2_label: "HA12",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "N2",
            atom2_label: "HN2",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA2",
            atom2_label: "HA21",
            order: TopologyBondOrder::Single,
        },
        OutputBondBlueprint {
            atom1_label: "CA2",
            atom2_label: "HA22",
            order: TopologyBondOrder::Single,
        },
    ],
};
