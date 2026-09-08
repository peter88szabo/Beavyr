//! Hosts the DREIDING typing pipeline: the rule schema, the built-in rule deck, and the engine.
//!
//! This namespace exposes the rule schema (`rules`), the generated default deck (`rules_data`),
//! and the iterative assignment engine (`engine`) used by `assign_topology`.

/// Typing engine that evaluates rules over annotated molecules.
pub mod engine;
/// Rule definitions.
pub mod rules;
/// The built-in DREIDING rule deck, as native Rust data.
pub mod rules_data;
