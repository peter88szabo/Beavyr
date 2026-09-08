//! The DREIDING generic force field.
//!
//! S. L. Mayo, B. D. Olafson and W. A. Goddard III, "DREIDING: A Generic Force Field for
//! Molecular Simulations", *J. Phys. Chem.* **1990**, *94*, 8897-8909. Equation and table
//! numbers in this module refer to that paper.
//!
//! DREIDING derives its parameters from element and hybridization rather than tabulating them per
//! functional group, so a small number of rules covers the whole main group and several metals.
//! That generality is what makes it the right fit for cleaning up structures a user has just
//! sketched, where the chemistry is not known in advance.
//!
//! # Layout
//!
//! * [`typer`] -- atom typing and topology perception (vendored; see its `README.md`)
//!
//! # Design
//!
//! See `docs/superpowers/specs/2026-09-09-dreiding-force-field-design.md`. In outline: an
//! immutable `DreidingTopology`, built once and shareable across threads, plus a per-worker
//! mutable `Workspace` holding the neighbour list and minimiser history. Geometry cleanup,
//! molecular dynamics and conformational search are then all drivers over one
//! `energy_and_forces`, rather than three engines.

pub mod typer;

/// DREIDING parameters, transcribed from the paper.
pub mod params;
