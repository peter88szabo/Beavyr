//! Estimating bond orders from geometry, and the valences they add up to.
//!
//! Only coordinates are available -- no Lewis structure, no formal charges --
//! so a bond's order is inferred from how short it is relative to the sum of
//! its atoms' covalent radii. That is enough to tell a double bond from a
//! single one, which is the difference between "this carbon is missing a
//! hydrogen" and "this carbon is part of a C=C".
//!
//! Counting *neighbours* instead, as a naive check does, calls every sp2
//! carbon under-coordinated and every transferring hydrogen over-coordinated.

use crate::molecule::{covalent_radius_angstrom, Molecule};

/// Covalent radii for bond-order work, from Cordero et al., *Dalton Trans.*
/// 2008, 2832 -- the single-bond radii the length ratios below are calibrated
/// against.
///
/// These are deliberately **not** `molecule::covalent_radius_angstrom`, whose
/// values are inflated so that bond *perception* errs towards finding a bond
/// (hydrogen there is 0.45 A against Cordero's 0.31). Inflated radii shrink
/// every ratio, and with them an ordinary C-H bond reads as a double and an
/// O-H bond as a triple -- so a water molecule's oxygen comes out hexavalent.
fn bond_order_radius(symbol: &str) -> f32 {
    match symbol {
        "H" => 0.31,
        "He" => 0.28,
        "Li" => 1.28,
        "Be" => 0.96,
        "B" => 0.84,
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "F" => 0.57,
        "Ne" => 0.58,
        "Na" => 1.66,
        "Mg" => 1.41,
        "Al" => 1.21,
        "Si" => 1.11,
        "P" => 1.07,
        "S" => 1.05,
        "Cl" => 1.02,
        "Ar" => 1.06,
        "K" => 2.03,
        "Ca" => 1.76,
        "Ge" => 1.20,
        "As" => 1.19,
        "Se" => 1.20,
        "Br" => 1.20,
        "Te" => 1.38,
        "I" => 1.39,
        // Anything not tabulated falls back to the perception radius. The
        // ratio will be a little small, but a missing element is worse.
        other => covalent_radius_angstrom(other),
    }
}

/// Bond-order tiers by length ratio `d / (r_a + r_b)`.
///
/// The reference ratios are the empirically calibrated ones already used in
/// `molecule_builder::attach`: triple ~0.80, double ~0.89, aromatic ~0.93,
/// single ~1.0. Boundaries sit at the midpoints between them.
const TRIPLE_MAX: f32 = 0.845;
const DOUBLE_MAX: f32 = 0.91;
const AROMATIC_MAX: f32 = 0.965;
/// Beyond this a bond is stretched well past a normal single bond: a partial
/// bond, forming or breaking. Counting it as a whole bond is what makes a
/// transition state's transferring hydrogen look divalent.
const SINGLE_MAX: f32 = 1.05;
/// What a stretched bond contributes. Half a bond is the right order of
/// magnitude at a transition state, where the transferring atom is roughly
/// half bonded to each partner -- and the two halves then sum to one whole
/// bond, which is what the atom actually has.
pub const PARTIAL_BOND_ORDER: f64 = 0.5;

/// The length a single bond between these two elements would have: the sum of
/// their covalent radii. Every bond-order judgement is a ratio against this.
pub fn single_bond_length(symbol_a: &str, symbol_b: &str) -> f32 {
    bond_order_radius(symbol_a) + bond_order_radius(symbol_b)
}

/// The estimated order of a bond of length `distance` between two elements.
pub fn estimate_bond_order(symbol_a: &str, symbol_b: &str, distance: f32) -> f64 {
    let single_bond_length = single_bond_length(symbol_a, symbol_b);
    if single_bond_length <= 0.0 {
        return 1.0;
    }
    let ratio = distance / single_bond_length;
    if ratio <= TRIPLE_MAX {
        3.0
    } else if ratio <= DOUBLE_MAX {
        2.0
    } else if ratio <= AROMATIC_MAX {
        1.5
    } else if ratio <= SINGLE_MAX {
        1.0
    } else {
        PARTIAL_BOND_ORDER
    }
}

/// The summed bond order at every atom, in atom order.
pub fn total_bond_orders(mol: &Molecule) -> Vec<f64> {
    let mut totals = vec![0.0; mol.atoms.len()];
    for &(i, j, distance) in &mol.bonds {
        if i >= totals.len() || j >= totals.len() {
            continue;
        }
        let order = estimate_bond_order(&mol.atoms[i], &mol.atoms[j], distance);
        totals[i] += order;
        totals[j] += order;
    }
    totals
}

/// The total bond orders an element is ordinarily found with, ascending.
///
/// `None` for elements with no useful rule -- transition metals above all,
/// whose coordination number is not a valence in this sense. Reporting
/// nothing is better than reporting a number that means nothing.
pub fn common_valences(symbol: &str) -> Option<&'static [f64]> {
    Some(match symbol {
        "H" => &[1.0],
        "B" | "Al" => &[3.0],
        "C" | "Si" => &[4.0],
        "N" => &[3.0],
        "O" => &[2.0],
        "F" => &[1.0],
        // Hypervalent: phosphorus as phosphine or phosphate, sulfur as
        // sulfide, sulfoxide or sulfone, the heavier halogens in their oxides.
        "P" => &[3.0, 5.0],
        "S" => &[2.0, 4.0, 6.0],
        "Cl" | "Br" | "I" => &[1.0, 3.0, 5.0, 7.0],
        _ => return None,
    })
}

/// How far a summed bond order may sit from a common valence and still count
/// as that valence.
///
/// Generous on purpose: the orders are inferred from bond lengths, so one bond
/// misread as aromatic rather than double already costs 0.5. The point is to
/// catch a whole missing hydrogen, not to police half a unit.
pub const VALENCE_TOLERANCE: f64 = 0.5;

/// What an atom's summed bond order says about it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValenceVerdict {
    /// No rule for this element, or the total matches a common valence.
    Satisfied,
    /// Below the lowest common valence: a hydrogen is probably missing, or
    /// the atom is a radical -- geometry alone cannot tell those apart.
    Short,
    /// Above the highest common valence: too many or too short bonds.
    Excess,
    /// Between two common valences, e.g. a sulfur summing to 3.
    Unusual,
}

/// Judges one atom's summed bond order against its element's common valences.
pub fn classify_valence(symbol: &str, total: f64) -> (ValenceVerdict, Option<f64>) {
    let Some(targets) = common_valences(symbol) else {
        return (ValenceVerdict::Satisfied, None);
    };
    // An isolated atom is not an under-coordinated one: a lone ion or a stray
    // fragment is a deliberate thing to have on screen.
    if total <= 0.0 {
        return (ValenceVerdict::Satisfied, None);
    }
    let nearest = targets
        .iter()
        .copied()
        .min_by(|a, b| (a - total).abs().total_cmp(&(b - total).abs()));
    if let Some(target) = nearest {
        if (total - target).abs() <= VALENCE_TOLERANCE {
            return (ValenceVerdict::Satisfied, Some(target));
        }
    }
    let lowest = targets[0];
    let highest = targets[targets.len() - 1];
    if total < lowest {
        (ValenceVerdict::Short, Some(lowest))
    } else if total > highest {
        (ValenceVerdict::Excess, Some(highest))
    } else {
        (ValenceVerdict::Unusual, nearest)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The calibrated reference lengths for a C-C bond, radius sum 1.52 A.
    #[test]
    fn the_reference_ratios_give_the_expected_orders() {
        assert_eq!(estimate_bond_order("C", "C", 1.54), 1.0, "single");
        assert_eq!(estimate_bond_order("C", "C", 1.40), 1.5, "aromatic");
        assert_eq!(estimate_bond_order("C", "C", 1.34), 2.0, "double");
        assert_eq!(estimate_bond_order("C", "C", 1.20), 3.0, "triple");
    }

    /// The bug this table exists to prevent: with bond-perception radii an
    /// ordinary C-H reads as a double bond and an O-H as a triple, making
    /// water's oxygen hexavalent.
    #[test]
    fn ordinary_single_bonds_to_hydrogen_are_single() {
        assert_eq!(estimate_bond_order("C", "H", 1.09), 1.0, "C-H");
        assert_eq!(estimate_bond_order("O", "H", 0.97), 1.0, "O-H");
        assert_eq!(estimate_bond_order("N", "H", 1.01), 1.0, "N-H");
        assert_eq!(estimate_bond_order("S", "H", 1.34), 1.0, "S-H");
    }

    /// A water oxygen is divalent, which is the end-to-end version of the
    /// same check.
    #[test]
    fn a_water_oxygen_is_divalent() {
        let total = 2.0 * estimate_bond_order("O", "H", 0.97);
        assert_eq!(total, 2.0);
        assert_eq!(classify_valence("O", total).0, ValenceVerdict::Satisfied);
    }

    /// Common multiple bonds keep their orders with these radii.
    #[test]
    fn characteristic_multiple_bonds_are_recognised() {
        assert_eq!(estimate_bond_order("C", "O", 1.23), 2.0, "carbonyl C=O");
        assert_eq!(estimate_bond_order("C", "O", 1.43), 1.0, "ether C-O");
        assert_eq!(estimate_bond_order("C", "N", 1.16), 3.0, "nitrile");
        assert_eq!(estimate_bond_order("C", "C", 1.40), 1.5, "aromatic");
    }

    /// A bond stretched past a normal single bond counts as partial. This is
    /// what stops a transferring hydrogen from looking divalent.
    #[test]
    fn a_stretched_bond_counts_as_partial() {
        // C-H radius sum is 1.07 A; a normal bond is 1.09.
        assert_eq!(estimate_bond_order("C", "H", 1.09), 1.0);
        // Halfway to the acceptor at a transfer transition state.
        assert_eq!(estimate_bond_order("C", "H", 1.30), PARTIAL_BOND_ORDER);
        assert_eq!(estimate_bond_order("O", "H", 1.30), PARTIAL_BOND_ORDER);
    }

    /// Two half bonds are one whole bond -- which is what the transferring
    /// atom actually has.
    #[test]
    fn two_partial_bonds_sum_to_one_valence() {
        let sum = estimate_bond_order("C", "H", 1.30) + estimate_bond_order("O", "H", 1.30);
        assert_eq!(sum, 1.0);
        assert_eq!(classify_valence("H", sum).0, ValenceVerdict::Satisfied);
    }

    #[test]
    fn an_unknown_element_has_no_rule_and_no_verdict() {
        assert!(common_valences("Fe").is_none());
        assert_eq!(classify_valence("Fe", 6.0).0, ValenceVerdict::Satisfied);
    }

    /// Ordinary saturated and unsaturated carbons alike are satisfied: 4 from
    /// four singles, and 4 from two singles plus a double.
    #[test]
    fn carbon_is_satisfied_by_singles_or_by_a_double() {
        assert_eq!(classify_valence("C", 4.0).0, ValenceVerdict::Satisfied);
        assert_eq!(classify_valence("C", 1.0 + 1.0 + 2.0).0, ValenceVerdict::Satisfied);
        // An aromatic ring carbon: two aromatic bonds and one single.
        assert_eq!(classify_valence("C", 1.5 + 1.5 + 1.0).0, ValenceVerdict::Satisfied);
    }

    #[test]
    fn a_carbon_missing_a_hydrogen_is_short() {
        assert_eq!(classify_valence("C", 3.0).0, ValenceVerdict::Short);
        assert_eq!(classify_valence("C", 3.0).1, Some(4.0));
        assert_eq!(classify_valence("N", 2.0).0, ValenceVerdict::Short);
        assert_eq!(classify_valence("O", 1.0).0, ValenceVerdict::Short);
    }

    #[test]
    fn a_carbon_with_five_bonds_is_in_excess() {
        assert_eq!(classify_valence("C", 5.0).0, ValenceVerdict::Excess);
        assert_eq!(classify_valence("H", 2.0).0, ValenceVerdict::Excess);
    }

    /// Hypervalent elements must not be flagged in any of their normal states.
    #[test]
    fn hypervalent_elements_are_satisfied_at_each_common_valence() {
        for total in [3.0, 5.0] {
            assert_eq!(classify_valence("P", total).0, ValenceVerdict::Satisfied, "P{total}");
        }
        for total in [2.0, 4.0, 6.0] {
            assert_eq!(classify_valence("S", total).0, ValenceVerdict::Satisfied, "S{total}");
        }
        for total in [1.0, 3.0, 5.0, 7.0] {
            assert_eq!(classify_valence("Cl", total).0, ValenceVerdict::Satisfied, "Cl{total}");
            assert_eq!(classify_valence("I", total).0, ValenceVerdict::Satisfied, "I{total}");
        }
        // ...but beyond the highest is still an error.
        assert_eq!(classify_valence("S", 8.0).0, ValenceVerdict::Excess);
        assert_eq!(classify_valence("P", 7.0).0, ValenceVerdict::Excess);
    }

    /// The elements the user named all have rules.
    #[test]
    fn every_requested_element_has_a_rule() {
        for symbol in ["C", "O", "N", "P", "S", "F", "Cl", "Br", "I", "B", "Al", "Si", "H"] {
            assert!(common_valences(symbol).is_some(), "{symbol} has no rule");
        }
    }

    /// Boron and aluminium are trivalent; a BF3 boron is complete at 3 and
    /// must not be reported as missing a fourth bond.
    #[test]
    fn trivalent_elements_are_complete_at_three() {
        assert_eq!(classify_valence("B", 3.0).0, ValenceVerdict::Satisfied);
        assert_eq!(classify_valence("Al", 3.0).0, ValenceVerdict::Satisfied);
        assert_eq!(classify_valence("B", 2.0).0, ValenceVerdict::Short);
    }

    /// A sulfur at 3 sits between its common 2 and 4.
    #[test]
    fn a_total_between_two_common_valences_is_unusual_not_short() {
        assert_eq!(classify_valence("S", 3.2).0, ValenceVerdict::Unusual);
    }

    /// A lone atom on screen is deliberate, not under-coordinated.
    #[test]
    fn an_isolated_atom_is_not_reported() {
        assert_eq!(classify_valence("C", 0.0).0, ValenceVerdict::Satisfied);
    }

    /// The tolerance must absorb one bond misread as aromatic rather than
    /// double, but not a whole missing hydrogen.
    #[test]
    fn the_tolerance_absorbs_a_misread_bond_but_not_a_missing_one() {
        assert_eq!(classify_valence("C", 3.5).0, ValenceVerdict::Satisfied, "one tier off");
        assert_eq!(classify_valence("C", 3.0).0, ValenceVerdict::Short, "a whole bond");
    }
}
