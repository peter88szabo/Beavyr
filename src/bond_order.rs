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
/// The same limit for a single bond between two second-period atoms that both
/// carry lone pairs -- N, O and F. Those lone pairs repel across the bond and
/// stretch it well past the radius sum, so the general limit calls a perfectly
/// ordinary bond half-formed. See `single_max_ratio`.
const LONE_PAIR_SINGLE_MAX: f32 = 1.16;
/// F-F is the extreme of the same effect and needs its own allowance: 1.412 A
/// against a 1.14 A radius sum is a ratio of 1.24.
const FLUORINE_SINGLE_MAX: f32 = 1.28;
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

/// How far past the radius sum a bond between these two elements can stretch
/// and still be a whole single bond.
///
/// Only the upper boundary moves: the tiers below it are calibrated on bonds
/// that are *shorter* than the radius sum, where nothing about this applies.
///
/// Lone pairs on both atoms repel across the bond and lengthen it. Hydrogen
/// peroxide's O-O is 1.458 A against a 1.32 A radius sum -- a ratio of 1.105,
/// past the general 1.05 -- and dialkyl peroxides reach 1.48. So every
/// peroxide read as two half bonds, and the diagnostics duly reported both
/// oxygens as radicals with a valence of 1.5. The same error hit hydroxylamine
/// (N-O 1.453 against 1.37) and NF3 (N-F 1.371 against 1.28).
///
/// This does not blunt the partial-bond detection that matters for reactive
/// structures: a peroxide bond in the middle of forming or breaking is 1.8 A
/// or more, a ratio above 1.36, so it still reads as partial.
fn single_max_ratio(symbol_a: &str, symbol_b: &str) -> f32 {
    let has_lone_pairs = |s: &str| matches!(s, "N" | "O" | "F");
    if !has_lone_pairs(symbol_a) || !has_lone_pairs(symbol_b) {
        return SINGLE_MAX;
    }
    // Fluorine's radius is the smallest of the three, so a bond to it is
    // under-estimated by the most: HOF's O-F is 1.442 A against a 1.23 A sum,
    // a ratio of 1.17, and F2's is 1.24.
    if symbol_a == "F" || symbol_b == "F" {
        FLUORINE_SINGLE_MAX
    } else {
        LONE_PAIR_SINGLE_MAX
    }
}

/// The estimated order of a bond of length `distance` between two elements.
pub fn estimate_bond_order(symbol_a: &str, symbol_b: &str, distance: f32) -> f64 {
    let single_bond_length = single_bond_length(symbol_a, symbol_b);
    if single_bond_length <= 0.0 {
        return 1.0;
    }
    let ratio = distance / single_bond_length;

    // Hydrogen has one electron and forms one bond, so no distance makes a bond to it a double.
    // The cap is only from above: a *stretched* H bond is still partial, which is how a hydrogen
    // in flight at a proton-transfer transition state gets its two half-bonds summing to one.
    //
    // Without the cap, a drawn O-H a few hundredths short -- 0.93 Å against 0.97 -- is promoted to
    // 1.5, and the error does not stay local: the oxygen types as a carbonyl, its neighbour picks
    // up resonance from it, and a hydroperoxide comes out `O_R, O_2` instead of `O_3, O_3`.
    if symbol_a == "H" || symbol_b == "H" {
        return if ratio <= single_max_ratio(symbol_a, symbol_b) {
            1.0
        } else {
            PARTIAL_BOND_ORDER
        };
    }

    if ratio <= TRIPLE_MAX {
        3.0
    } else if ratio <= DOUBLE_MAX {
        2.0
    } else if ratio <= AROMATIC_MAX {
        1.5
    } else if ratio <= single_max_ratio(symbol_a, symbol_b) {
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

/// The bond angle an atom of this element adopts once it has `sigma_bonds`
/// sigma bonds in total, in degrees.
///
/// VSEPR, reduced to the cases that arise when completing a partly-built
/// centre: the steric number is the sigma-bond count plus the lone pairs the
/// element carries, and the angle follows from that. Used only when geometry
/// cannot settle the direction on its own -- with two or more neighbours
/// already present, the existing arrangement determines where the next bond
/// goes far more reliably than a table does.
pub fn ideal_bond_angle_deg(symbol: &str, sigma_bonds: usize) -> f64 {
    const TRIGONAL: f64 = 120.0;
    const TETRAHEDRAL: f64 = 109.471;

    match symbol {
        // sp3 unless the geometry proves otherwise. This table is only
        // consulted when the existing arrangement is too sparse or too
        // distorted to measure a direction from, and in that situation a
        // partly-built carbon is overwhelmingly more likely to be heading for
        // a tetrahedral centre than a linear or trigonal one -- a carbon with
        // one neighbour is far more often a future CH3 than an alkyne. Real
        // sp and sp2 centres are caught before this by their bond lengths and
        // their neighbours' angles.
        "C" | "Si" => TETRAHEDRAL,
        // Boron and aluminium are the exception: with no lone pair they are
        // genuinely trigonal planar, and BH3 or AlCl3 is the common case
        // rather than a tetrahedral borate.
        "B" | "Al" => TRIGONAL,
        // One lone pair: pyramidal, slightly tighter than tetrahedral.
        "N" | "P" => match sigma_bonds {
            0 | 1 | 2 | 3 => 107.0,
            _ => TETRAHEDRAL,
        },
        // Two lone pairs: the familiar bent geometry of water and ethers.
        "O" | "S" | "Se" => match sigma_bonds {
            0 | 1 | 2 => 104.5,
            _ => TETRAHEDRAL,
        },
        // Everything else, including a hydrogen or halogen that should not be
        // gaining a second bond at all: tetrahedral is the least surprising.
        _ => TETRAHEDRAL,
    }
}

/// Whether this element has room for another bond, given what it already has.
///
/// `Some(deficit)` when it does, rounded to whole bonds; `None` when the
/// element has no rule or is already at its highest common valence.
pub fn room_for_another_bond(symbol: &str, current_total: f64) -> Option<f64> {
    let targets = common_valences(symbol)?;
    let highest = targets[targets.len() - 1];
    let deficit = highest - current_total;
    (deficit >= 1.0 - VALENCE_TOLERANCE).then_some(deficit)
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
    /// A bond to hydrogen is never *more* than single, however short it is drawn. Hydrogen has
    /// one electron; no distance makes it form two bonds.
    ///
    /// The cap is one-sided on purpose. A stretched H bond stays partial, which is what gives a
    /// hydrogen in flight at a proton-transfer transition state two half-bonds -- see
    /// `two_partial_bonds_sum_to_one_valence`.
    #[test]
    fn a_bond_to_hydrogen_is_never_more_than_single() {
        // Short of the normal length, and well short of it: still single, never promoted.
        for distance in [0.60_f32, 0.80, 0.93, 0.97] {
            assert_eq!(
                estimate_bond_order("O", "H", distance),
                1.0,
                "O-H at {distance} Å"
            );
        }
        for distance in [0.70_f32, 0.90, 1.09] {
            assert_eq!(estimate_bond_order("H", "C", distance), 1.0);
        }
        // And a stretched one is still partial, not forced to single.
        assert_eq!(estimate_bond_order("O", "H", 1.30), PARTIAL_BOND_ORDER);
        assert_eq!(estimate_bond_order("C", "H", 1.30), PARTIAL_BOND_ORDER);
    }

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

    /// The reported bug: an -OOH group's O-O read as two half bonds, so both
    /// oxygens were reported with a valence of 1.5 and flagged as possible
    /// radicals. A peroxide bond is a perfectly ordinary single bond that
    /// happens to be longer than the sum of the covalent radii.
    #[test]
    fn a_peroxide_bond_is_a_whole_single_bond() {
        // H2O2, gas phase.
        assert_eq!(estimate_bond_order("O", "O", 1.458), 1.0);
        // Dialkyl and dialkyl-substituted peroxides span roughly 1.45-1.48.
        assert_eq!(estimate_bond_order("O", "O", 1.45), 1.0);
        assert_eq!(estimate_bond_order("O", "O", 1.48), 1.0);

        // So a hydroperoxide oxygen is divalent, and says nothing about
        // radicals: one bond to carbon, one to the other oxygen.
        let valence = estimate_bond_order("C", "O", 1.43) + estimate_bond_order("O", "O", 1.458);
        assert_eq!(valence, 2.0);
        assert_eq!(classify_valence("O", valence).0, ValenceVerdict::Satisfied);
    }

    /// The same lone-pair lengthening in the other bonds it affects.
    #[test]
    fn single_bonds_between_lone_pair_atoms_are_not_partial() {
        assert_eq!(estimate_bond_order("N", "N", 1.447), 1.0, "hydrazine");
        assert_eq!(estimate_bond_order("N", "O", 1.453), 1.0, "hydroxylamine");
        assert_eq!(estimate_bond_order("N", "F", 1.371), 1.0, "NF3");
        assert_eq!(estimate_bond_order("O", "F", 1.442), 1.0, "HOF");
        assert_eq!(estimate_bond_order("F", "F", 1.412), 1.0, "F2");
    }

    /// The wider limit is granted only to the pairs that need it, so an
    /// ordinary stretched bond is still detected as partial.
    #[test]
    fn the_lone_pair_allowance_does_not_leak_to_other_elements() {
        // A carbon or hydrogen partner gets the general limit, however
        // lone-pair-rich the other atom is. A C-O at 1.55 A is only a ratio
        // of 1.09 -- inside the lone-pair allowance, outside the general one.
        assert_eq!(estimate_bond_order("C", "O", 1.55), PARTIAL_BOND_ORDER);
        assert_eq!(estimate_bond_order("O", "H", 1.30), PARTIAL_BOND_ORDER);
        assert_eq!(estimate_bond_order("C", "N", 1.70), PARTIAL_BOND_ORDER);
    }

    /// A peroxide bond being broken or formed must still read as partial:
    /// this is the case the wider limit could have cost us.
    #[test]
    fn a_forming_or_breaking_peroxide_bond_is_still_partial() {
        assert_eq!(estimate_bond_order("O", "O", 1.75), PARTIAL_BOND_ORDER);
        assert_eq!(estimate_bond_order("O", "O", 2.10), PARTIAL_BOND_ORDER);
        assert_eq!(estimate_bond_order("N", "O", 1.90), PARTIAL_BOND_ORDER);
    }

    /// The wider limit moved only the stretched boundary, so a *short* O-O is
    /// read exactly as it was before: dioxygen's 1.208 A still reads as
    /// multiply bonded.
    ///
    /// Ozone, whose terminal O-O is 1.278 A and whose true order is 1.5, sits
    /// just inside the single tier at a ratio of 0.968. That is a limitation
    /// of judging delocalised bonds by length alone -- the O-O single bond is
    /// anomalously long and the double anomalously short, so no single ratio
    /// scale fits both ends -- and it is the same answer this gave before the
    /// stretched limit was widened. Recorded rather than asserted as correct.
    #[test]
    fn widening_the_stretched_limit_left_short_bonds_alone() {
        assert_eq!(estimate_bond_order("O", "O", 1.208), 1.5, "dioxygen");
        assert_eq!(estimate_bond_order("O", "O", 1.278), 1.0, "ozone");
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

    /// Carbon defaults to sp3 at every coordination, because this table is
    /// only reached when the geometry could not be measured -- and a
    /// half-built carbon is far more often a future CH3 than an alkyne.
    #[test]
    fn carbon_always_falls_back_to_sp3() {
        for sigma in 0..=4 {
            assert!(
                (ideal_bond_angle_deg("C", sigma) - 109.471).abs() < 1e-3,
                "C with {sigma} sigma bonds gave {}",
                ideal_bond_angle_deg("C", sigma)
            );
        }
        assert!((ideal_bond_angle_deg("Si", 2) - 109.471).abs() < 1e-3);
    }

    /// The others keep the geometry their lone pairs dictate.
    #[test]
    fn the_ideal_angles_follow_vsepr() {
        assert_eq!(ideal_bond_angle_deg("O", 2), 104.5, "water");
        assert_eq!(ideal_bond_angle_deg("N", 3), 107.0, "ammonia");
        assert_eq!(ideal_bond_angle_deg("N", 2), 107.0, "a half-built amine");
        assert_eq!(ideal_bond_angle_deg("B", 3), 120.0, "borane, trigonal planar");
        assert_eq!(ideal_bond_angle_deg("Al", 2), 120.0);
    }

    /// A completed centre has no room; a partly built one does.
    #[test]
    fn room_is_reported_only_where_a_bond_can_go() {
        assert!(room_for_another_bond("C", 4.0).is_none(), "methane carbon is full");
        assert!(room_for_another_bond("C", 3.0).is_some(), "a methyl radical has room");
        assert!(room_for_another_bond("H", 1.0).is_none(), "hydrogen is monovalent");
        assert!(room_for_another_bond("H", 0.0).is_some(), "a bare H can bond once");
        assert!(room_for_another_bond("O", 2.0).is_none(), "water oxygen is full");
        // Hypervalent elements are judged against their highest state.
        assert!(room_for_another_bond("S", 2.0).is_some(), "sulfide can reach 6");
        assert!(room_for_another_bond("S", 6.0).is_none(), "a sulfone is full");
        assert!(room_for_another_bond("Fe", 3.0).is_none(), "no rule, no claim");
    }

    /// Half a bond short is within tolerance, so not "room" for a whole one.
    #[test]
    fn a_fraction_of_a_bond_is_not_room_for_another() {
        assert!(room_for_another_bond("C", 3.6).is_none());
        assert!(room_for_another_bond("C", 3.4).is_some());
    }

    /// The tolerance must absorb one bond misread as aromatic rather than
    /// double, but not a whole missing hydrogen.
    #[test]
    fn the_tolerance_absorbs_a_misread_bond_but_not_a_missing_one() {
        assert_eq!(classify_valence("C", 3.5).0, ValenceVerdict::Satisfied, "one tier off");
        assert_eq!(classify_valence("C", 3.0).0, ValenceVerdict::Short, "a whole bond");
    }
}
