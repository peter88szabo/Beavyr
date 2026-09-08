//! DREIDING atom-typing rule deck, as native Rust data.
//!
//! GENERATED -- do not edit by hand.
//!
//! Upstream ships this deck as `resources/default.rules.toml` and parses it at
//! runtime with `serde` + `toml`. Beavyr vendors the typer and keeps its
//! dependency count at five, so the deck is translated to Rust literals here
//! instead and the TOML parser is gone. The translation was produced
//! mechanically from the upstream file; the ported upstream test suite is what
//! proves it faithful.
//!
//! Ruleset version 2.2.0, by Tony Kan and William A. Goddard III.

use std::collections::HashMap;
use std::sync::OnceLock;

use super::rules::{Conditions, Rule};
use crate::forcefield::dreiding::typer::core::properties::{Element, Hybridization};

static DEFAULT_RULES: OnceLock<Vec<Rule>> = OnceLock::new();

/// Returns the embedded DREIDING rule deck, built once and cached.
pub fn get_default_rules() -> &'static [Rule] {
    DEFAULT_RULES.get_or_init(build_default_rules)
}

#[rustfmt::skip]
fn build_default_rules() -> Vec<Rule> {
    vec![
        Rule {
            name: "H_Bridge_Diborane".to_string(),
            priority: 500,
            result_type: "H_b".to_string(),
            conditions: Conditions {
                element: Some(Element::H),
                degree: Some(2),
                neighbor_elements: HashMap::from([(Element::B, 2)]),
                ..Conditions::default()
            },
        },
        Rule {
            name: "C_Resonant_Terminal".to_string(),
            priority: 401,
            result_type: "C_2".to_string(),
            conditions: Conditions {
                element: Some(Element::C),
                degree: Some(1),
                hybridization: Some(Hybridization::Resonant),
                ..Conditions::default()
            },
        },
        Rule {
            name: "N_Resonant_Terminal".to_string(),
            priority: 401,
            result_type: "N_2".to_string(),
            conditions: Conditions {
                element: Some(Element::N),
                degree: Some(1),
                hybridization: Some(Hybridization::Resonant),
                ..Conditions::default()
            },
        },
        Rule {
            name: "O_Resonant_Terminal".to_string(),
            priority: 401,
            result_type: "O_2".to_string(),
            conditions: Conditions {
                element: Some(Element::O),
                degree: Some(1),
                hybridization: Some(Hybridization::Resonant),
                ..Conditions::default()
            },
        },
        Rule {
            name: "S_Resonant_Terminal".to_string(),
            priority: 401,
            result_type: "S_2".to_string(),
            conditions: Conditions {
                element: Some(Element::S),
                degree: Some(1),
                hybridization: Some(Hybridization::Resonant),
                ..Conditions::default()
            },
        },
        Rule {
            name: "C_Resonant".to_string(),
            priority: 400,
            result_type: "C_R".to_string(),
            conditions: Conditions {
                element: Some(Element::C),
                hybridization: Some(Hybridization::Resonant),
                ..Conditions::default()
            },
        },
        Rule {
            name: "N_Resonant".to_string(),
            priority: 400,
            result_type: "N_R".to_string(),
            conditions: Conditions {
                element: Some(Element::N),
                hybridization: Some(Hybridization::Resonant),
                ..Conditions::default()
            },
        },
        Rule {
            name: "O_Resonant".to_string(),
            priority: 400,
            result_type: "O_R".to_string(),
            conditions: Conditions {
                element: Some(Element::O),
                hybridization: Some(Hybridization::Resonant),
                ..Conditions::default()
            },
        },
        Rule {
            name: "S_Resonant".to_string(),
            priority: 400,
            result_type: "S_R".to_string(),
            conditions: Conditions {
                element: Some(Element::S),
                hybridization: Some(Hybridization::Resonant),
                ..Conditions::default()
            },
        },
        Rule {
            name: "C_Linear_SP".to_string(),
            priority: 300,
            result_type: "C_1".to_string(),
            conditions: Conditions {
                element: Some(Element::C),
                hybridization: Some(Hybridization::SP),
                ..Conditions::default()
            },
        },
        Rule {
            name: "N_Linear_SP".to_string(),
            priority: 300,
            result_type: "N_1".to_string(),
            conditions: Conditions {
                element: Some(Element::N),
                hybridization: Some(Hybridization::SP),
                ..Conditions::default()
            },
        },
        Rule {
            name: "C_Trigonal_SP2".to_string(),
            priority: 200,
            result_type: "C_2".to_string(),
            conditions: Conditions {
                element: Some(Element::C),
                hybridization: Some(Hybridization::SP2),
                ..Conditions::default()
            },
        },
        Rule {
            name: "N_Trigonal_SP2".to_string(),
            priority: 200,
            result_type: "N_2".to_string(),
            conditions: Conditions {
                element: Some(Element::N),
                hybridization: Some(Hybridization::SP2),
                ..Conditions::default()
            },
        },
        Rule {
            name: "O_Trigonal_SP2".to_string(),
            priority: 200,
            result_type: "O_2".to_string(),
            conditions: Conditions {
                element: Some(Element::O),
                hybridization: Some(Hybridization::SP2),
                ..Conditions::default()
            },
        },
        Rule {
            name: "S_Trigonal_SP2".to_string(),
            priority: 200,
            result_type: "S_2".to_string(),
            conditions: Conditions {
                element: Some(Element::S),
                hybridization: Some(Hybridization::SP2),
                ..Conditions::default()
            },
        },
        Rule {
            name: "B_Trigonal_SP2".to_string(),
            priority: 200,
            result_type: "B_2".to_string(),
            conditions: Conditions {
                element: Some(Element::B),
                hybridization: Some(Hybridization::SP2),
                ..Conditions::default()
            },
        },
        Rule {
            name: "C_Tetrahedral_SP3".to_string(),
            priority: 100,
            result_type: "C_3".to_string(),
            conditions: Conditions {
                element: Some(Element::C),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "N_Tetrahedral_SP3".to_string(),
            priority: 100,
            result_type: "N_3".to_string(),
            conditions: Conditions {
                element: Some(Element::N),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "O_Tetrahedral_SP3".to_string(),
            priority: 100,
            result_type: "O_3".to_string(),
            conditions: Conditions {
                element: Some(Element::O),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "P_Tetrahedral_SP3".to_string(),
            priority: 100,
            result_type: "P_3".to_string(),
            conditions: Conditions {
                element: Some(Element::P),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "S_Tetrahedral_SP3".to_string(),
            priority: 100,
            result_type: "S_3".to_string(),
            conditions: Conditions {
                element: Some(Element::S),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "B_Tetrahedral_SP3".to_string(),
            priority: 100,
            result_type: "B_3".to_string(),
            conditions: Conditions {
                element: Some(Element::B),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Si_Tetrahedral_SP3".to_string(),
            priority: 101,
            result_type: "Si3".to_string(),
            conditions: Conditions {
                element: Some(Element::Si),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Al_Tetrahedral_SP3".to_string(),
            priority: 101,
            result_type: "Al3".to_string(),
            conditions: Conditions {
                element: Some(Element::Al),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Ge_Tetrahedral_SP3".to_string(),
            priority: 101,
            result_type: "Ge3".to_string(),
            conditions: Conditions {
                element: Some(Element::Ge),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "As_Tetrahedral_SP3".to_string(),
            priority: 101,
            result_type: "As3".to_string(),
            conditions: Conditions {
                element: Some(Element::As),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Se_Tetrahedral_SP3".to_string(),
            priority: 101,
            result_type: "Se3".to_string(),
            conditions: Conditions {
                element: Some(Element::Se),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Sn_Tetrahedral_SP3".to_string(),
            priority: 101,
            result_type: "Sn3".to_string(),
            conditions: Conditions {
                element: Some(Element::Sn),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Sb_Tetrahedral_SP3".to_string(),
            priority: 101,
            result_type: "Sb3".to_string(),
            conditions: Conditions {
                element: Some(Element::Sb),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Te_Tetrahedral_SP3".to_string(),
            priority: 101,
            result_type: "Te3".to_string(),
            conditions: Conditions {
                element: Some(Element::Te),
                hybridization: Some(Hybridization::SP3),
                ..Conditions::default()
            },
        },
        Rule {
            name: "H_Donor_On_Oxygen".to_string(),
            priority: 80,
            result_type: "H_HB".to_string(),
            conditions: Conditions {
                element: Some(Element::H),
                neighbor_elements: HashMap::from([(Element::O, 1)]),
                ..Conditions::default()
            },
        },
        Rule {
            name: "H_Donor_On_Nitrogen".to_string(),
            priority: 79,
            result_type: "H_HB".to_string(),
            conditions: Conditions {
                element: Some(Element::H),
                neighbor_elements: HashMap::from([(Element::N, 1)]),
                ..Conditions::default()
            },
        },
        Rule {
            name: "H_Donor_On_Sulfur".to_string(),
            priority: 78,
            result_type: "H_HB".to_string(),
            conditions: Conditions {
                element: Some(Element::H),
                neighbor_elements: HashMap::from([(Element::S, 1)]),
                ..Conditions::default()
            },
        },
        Rule {
            name: "H_Standard_Default".to_string(),
            priority: 1,
            result_type: "H_".to_string(),
            conditions: Conditions {
                element: Some(Element::H),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Halogen_F".to_string(),
            priority: 50,
            result_type: "F_".to_string(),
            conditions: Conditions {
                element: Some(Element::F),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Halogen_Cl".to_string(),
            priority: 50,
            result_type: "Cl".to_string(),
            conditions: Conditions {
                element: Some(Element::Cl),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Halogen_Br".to_string(),
            priority: 50,
            result_type: "Br".to_string(),
            conditions: Conditions {
                element: Some(Element::Br),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Halogen_I".to_string(),
            priority: 50,
            result_type: "I_".to_string(),
            conditions: Conditions {
                element: Some(Element::I),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Ag".to_string(),
            priority: 20,
            result_type: "Ag".to_string(),
            conditions: Conditions {
                element: Some(Element::Ag),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Au".to_string(),
            priority: 20,
            result_type: "Au".to_string(),
            conditions: Conditions {
                element: Some(Element::Au),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Ba".to_string(),
            priority: 20,
            result_type: "Ba".to_string(),
            conditions: Conditions {
                element: Some(Element::Ba),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Ca".to_string(),
            priority: 20,
            result_type: "Ca".to_string(),
            conditions: Conditions {
                element: Some(Element::Ca),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Cd".to_string(),
            priority: 20,
            result_type: "Cd".to_string(),
            conditions: Conditions {
                element: Some(Element::Cd),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Co".to_string(),
            priority: 20,
            result_type: "Co".to_string(),
            conditions: Conditions {
                element: Some(Element::Co),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Cs".to_string(),
            priority: 20,
            result_type: "Cs".to_string(),
            conditions: Conditions {
                element: Some(Element::Cs),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Cu".to_string(),
            priority: 20,
            result_type: "Cu".to_string(),
            conditions: Conditions {
                element: Some(Element::Cu),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Fe".to_string(),
            priority: 20,
            result_type: "Fe".to_string(),
            conditions: Conditions {
                element: Some(Element::Fe),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Hg".to_string(),
            priority: 20,
            result_type: "Hg".to_string(),
            conditions: Conditions {
                element: Some(Element::Hg),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_K".to_string(),
            priority: 20,
            result_type: "K".to_string(),
            conditions: Conditions {
                element: Some(Element::K),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Li".to_string(),
            priority: 20,
            result_type: "Li".to_string(),
            conditions: Conditions {
                element: Some(Element::Li),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Mg".to_string(),
            priority: 20,
            result_type: "Mg".to_string(),
            conditions: Conditions {
                element: Some(Element::Mg),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Mn".to_string(),
            priority: 20,
            result_type: "Mn".to_string(),
            conditions: Conditions {
                element: Some(Element::Mn),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Na".to_string(),
            priority: 20,
            result_type: "Na".to_string(),
            conditions: Conditions {
                element: Some(Element::Na),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Ni".to_string(),
            priority: 20,
            result_type: "Ni".to_string(),
            conditions: Conditions {
                element: Some(Element::Ni),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Pd".to_string(),
            priority: 20,
            result_type: "Pd".to_string(),
            conditions: Conditions {
                element: Some(Element::Pd),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Pt".to_string(),
            priority: 20,
            result_type: "Pt".to_string(),
            conditions: Conditions {
                element: Some(Element::Pt),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Rb".to_string(),
            priority: 20,
            result_type: "Rb".to_string(),
            conditions: Conditions {
                element: Some(Element::Rb),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Ru".to_string(),
            priority: 20,
            result_type: "Ru".to_string(),
            conditions: Conditions {
                element: Some(Element::Ru),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Sr".to_string(),
            priority: 20,
            result_type: "Sr".to_string(),
            conditions: Conditions {
                element: Some(Element::Sr),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Tc".to_string(),
            priority: 20,
            result_type: "Tc".to_string(),
            conditions: Conditions {
                element: Some(Element::Tc),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Ti".to_string(),
            priority: 20,
            result_type: "Ti".to_string(),
            conditions: Conditions {
                element: Some(Element::Ti),
                ..Conditions::default()
            },
        },
        Rule {
            name: "Metal_Zn".to_string(),
            priority: 20,
            result_type: "Zn".to_string(),
            conditions: Conditions {
                element: Some(Element::Zn),
                ..Conditions::default()
            },
        },
    ]
}
