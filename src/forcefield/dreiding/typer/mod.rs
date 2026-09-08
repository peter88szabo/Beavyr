//! DREIDING force field atom typing and molecular topology perception.
//!
//! # Provenance
//!
//! This is a **vendored copy** of `dreid-typer` v0.4.0, MIT licensed:
//!
//! * upstream: <https://github.com/caltechmsc/dreid-typer>
//! * fork used: <https://github.com/peter88szabo/dreid-typer>
//! * `Copyright (c) 2025 California Institute of Technology, Materials and
//!   Process Simulation Center (MSC)`
//! * authors: Tony Kan and William A. Goddard III -- a co-author of the
//!   DREIDING paper this force field implements
//!
//! The licence text is in `LICENSE-MIT` beside this file. `README.md` records
//! exactly what was changed relative to upstream and why.
//!
//! It is copied in rather than declared as a dependency so that Beavyr builds
//! from its own checkout with nothing to fetch or link.
//!
//! # What it does
//!
//! A deterministic three-phase pipeline (**perceive -> type -> build**) turns a
//! plain connectivity graph into a simulation-ready topology: DREIDING atom
//! types plus the complete lists of bonds, angles, proper dihedrals and
//! inversion centres. Getting these right by hand is tedious and the failure
//! mode is silent -- a mistyped atom yields a plausible but wrong geometry --
//! which is why validated code is used here instead of a fresh implementation.
//!
//! # Example
//!
//! ```ignore
//! use crate::forcefield::dreiding::typer::{
//!     assign_topology, Element, GraphBondOrder, MolecularGraph,
//! };
//!
//! let mut graph = MolecularGraph::new();
//! let o = graph.add_atom(Element::O);
//! let h1 = graph.add_atom(Element::H);
//! let h2 = graph.add_atom(Element::H);
//! graph.add_bond(o, h1, GraphBondOrder::Single).unwrap();
//! graph.add_bond(o, h2, GraphBondOrder::Single).unwrap();
//!
//! let topology = assign_topology(&graph).unwrap();
//! assert_eq!(topology.atoms[o].atom_type, "O_3");
//! assert_eq!(topology.angles.len(), 1);
//! ```

pub mod builder;
pub mod core;
pub mod perception;
pub mod typing;

pub use self::core::error::{AssignmentError, GraphValidationError, PerceptionError, TyperError};
pub use self::core::graph::{AtomNode, BondEdge, MolecularGraph};
pub use self::core::properties::{
    Element, GraphBondOrder, Hybridization, ParseBondOrderError, ParseElementError,
    ParseHybridizationError, TopologyBondOrder,
};
pub use self::core::topology::{
    Angle, Atom, Bond, ImproperDihedral, MolecularTopology, ProperDihedral,
};

/// Rule customization utilities.
///
/// Upstream parses rules from TOML; here the built-in deck is native Rust data
/// and a custom deck is built by constructing [`Rule`](rules::Rule) values
/// directly.
pub mod rules {
    pub use super::typing::rules::{Conditions, Rule};
    pub use super::typing::rules_data::get_default_rules;
}

/// Assigns a full molecular topology using the built-in DREIDING ruleset.
///
/// This is the primary entry point: it runs chemical perception, atom typing and
/// topology construction, returning types plus every bonded term.
///
/// # Errors
///
/// Returns a [`TyperError`] if the input graph is inconsistent, if perception
/// fails, or if the rule engine cannot type every atom.
pub fn assign_topology(graph: &MolecularGraph) -> Result<MolecularTopology, TyperError> {
    assign_topology_internal(graph, rules::get_default_rules())
}

/// Assigns a full molecular topology using a custom set of typing rules.
///
/// # Errors
///
/// As [`assign_topology`].
pub fn assign_topology_with_rules(
    graph: &MolecularGraph,
    rules: &[rules::Rule],
) -> Result<MolecularTopology, TyperError> {
    assign_topology_internal(graph, rules)
}

/// Internal core function that executes the perception, typing, and building pipeline.
fn assign_topology_internal(
    graph: &MolecularGraph,
    rules: &[rules::Rule],
) -> Result<MolecularTopology, TyperError> {
    let annotated_molecule = perception::perceive(graph)?;

    let atom_types = typing::engine::assign_types(&annotated_molecule, rules)
        .map_err(TyperError::AssignmentFailed)?;

    let topology = builder::build_topology(&annotated_molecule, &atom_types);

    Ok(topology)
}

#[cfg(test)]
mod tests;
