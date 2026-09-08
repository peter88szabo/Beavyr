//! Defines the rule language that drives the DREIDING typing engine.
//!
//! This module owns the `Rule` schema. Upstream also parses that schema from TOML with `serde`;
//! Beavyr vendors the typer without a TOML parser, so the built-in deck lives as native Rust data
//! in [`super::rules_data`] and the parsing helpers are gone. `Rule` and `Conditions` are still
//! public, so a caller can build a custom deck in Rust and pass it to `assign_topology_with_rules`.

use crate::forcefield::dreiding::typer::core::properties::{Element, Hybridization};
use std::collections::HashMap;

/// User-facing rule that assigns a DREIDING type when its conditions match an atom.
///
/// Rules are ordered by priority and evaluated by the typing engine until the first rule whose
/// `conditions` filters match succeeds, at which point the `result_type` label is emitted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Rule {
    /// Human-readable identifier that makes debugging rule decks easier.
    pub name: String,
    /// Ordering value where larger priorities are evaluated before smaller ones.
    pub priority: i32,
    /// Resulting atom type string written into the topology output.
    /// (Named `type` in the upstream TOML schema, which is a Rust keyword.)
    pub result_type: String,
    /// Set of optional property filters that determine when the rule fires.
    pub conditions: Conditions,
}

/// Optional property filters attached to a [`Rule`].
///
/// Each field defaults to `None` or an empty map, meaning the corresponding property is not
/// considered when deciding whether a rule applies. When populated, the field describes an exact
/// value that must match the target atom.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Conditions {
    /// Required element identity.
    pub element: Option<Element>,
    /// Required formal charge assigned during perception.
    pub formal_charge: Option<i8>,

    /// Required degree (number of σ bonds) for the atom.
    pub degree: Option<u8>,
    /// Whether the atom must belong (or not belong) to a ring system.
    pub is_in_ring: Option<bool>,

    /// Required lone-pair count after electron perception.
    pub lone_pairs: Option<u8>,
    /// Required hybridization assignment (SP, SP2, SP3, etc.).
    pub hybridization: Option<Hybridization>,

    /// Whether the atom must be aromatic.
    pub is_aromatic: Option<bool>,
    /// Whether the atom must be anti-aromatic.
    pub is_anti_aromatic: Option<bool>,
    /// Whether the atom must participate in any resonance system.
    pub is_resonant: Option<bool>,

    /// Minimum counts for neighbor elements keyed by element symbol strings.
    pub neighbor_elements: HashMap<Element, u8>,
    /// Minimum counts for neighbor atom types identified by their DREIDING labels.
    pub neighbor_types: HashMap<String, u8>,
}


#[cfg(test)]
mod tests {
    use super::super::rules_data::get_default_rules;
    use super::*;

    /// The deck is generated from upstream's `default.rules.toml`. If the generator ever drops or
    /// mangles entries, the count moves -- and the ported upstream suite (which exercises the deck
    /// against the molecules in the DREIDING paper) is what proves the contents faithful.
    #[test]
    fn default_deck_has_every_upstream_rule() {
        assert_eq!(get_default_rules().len(), 62);
    }

    #[test]
    fn default_deck_is_cached_not_rebuilt() {
        assert!(std::ptr::eq(get_default_rules(), get_default_rules()));
    }

    /// Every rule must name a type and constrain something; a rule with empty conditions would
    /// match every atom and silently shadow the rest of the deck.
    #[test]
    fn every_rule_names_a_type_and_constrains_something() {
        for rule in get_default_rules() {
            assert!(!rule.result_type.is_empty(), "{} has no type", rule.name);
            assert_ne!(
                rule.conditions,
                Conditions::default(),
                "{} would match every atom",
                rule.name
            );
        }
    }

    /// The engine resolves ties by priority, so two rules sharing a priority *and* a condition set
    /// would make typing depend on deck order rather than on chemistry.
    #[test]
    fn no_two_rules_share_a_priority_and_condition_set() {
        let rules = get_default_rules();
        for (i, a) in rules.iter().enumerate() {
            for b in &rules[i + 1..] {
                assert!(
                    a.priority != b.priority || a.conditions != b.conditions,
                    "{} and {} are indistinguishable",
                    a.name,
                    b.name
                );
            }
        }
    }
}
