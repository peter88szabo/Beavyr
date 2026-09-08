//! Error types describing the failure modes of graph validation, chemical perception, and typing.
//!
//! These enums aggregate lower-level issues so that library consumers can bubble up a single
//! `TyperError` while still inspecting fine-grained context when needed.
//!
//! Upstream derives these with `thiserror`. Beavyr vendors the typer and keeps its dependency
//! count at five, so the `Display`, `Error` and `From` impls are written out by hand here. The
//! messages are unchanged.

use std::error::Error;
use std::fmt;

/// Root error emitted by every fallible operation in the typing pipeline.
///
/// Each variant wraps a more specific error that pinpoints the subsystem that failed, allowing
/// callers to recover or log richer diagnostics without losing ergonomic `Result` signatures.
///
/// Upstream carries a `RuleParse(toml::de::Error)` variant. There is no TOML parser here -- the
/// rule deck is native Rust data in `typing::rules_data` -- so that variant does not exist.
#[derive(Debug)]
pub enum TyperError {
    /// Input validation of the `MolecularGraph` failed before perception could start.
    InvalidInput(GraphValidationError),

    /// A specific chemical perception stage reported a failure.
    PerceptionFailed {
        /// Name of the perception step (e.g., "aromaticity" or "hybridization").
        step: String,
        /// Root perception error that triggered the failure.
        source: PerceptionError,
    },

    /// The typing engine exhausted its rounds before assigning all atom types.
    AssignmentFailed(AssignmentError),
}

impl fmt::Display for TyperError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(_) => write!(f, "invalid input graph"),
            Self::PerceptionFailed { step, .. } => {
                write!(f, "chemical perception failed during '{step}' step")
            }
            Self::AssignmentFailed(_) => write!(f, "atom typing failed"),
        }
    }
}

impl Error for TyperError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::InvalidInput(e) => Some(e),
            Self::PerceptionFailed { source, .. } => Some(source),
            Self::AssignmentFailed(e) => Some(e),
        }
    }
}

impl From<GraphValidationError> for TyperError {
    fn from(e: GraphValidationError) -> Self {
        Self::InvalidInput(e)
    }
}

impl From<AssignmentError> for TyperError {
    fn from(e: AssignmentError) -> Self {
        Self::AssignmentFailed(e)
    }
}

/// Errors that describe structural or logical issues with the input `MolecularGraph`.
///
/// These failures are detected before any chemical reasoning is attempted so that malformed inputs
/// can be rejected early with precise diagnostics.
#[derive(Debug)]
pub enum GraphValidationError {
    /// A bond references an atom identifier that is missing from the graph.
    MissingAtom {
        /// Identifier of the atom that could not be found.
        atom_id: usize,
    },

    /// An atom record lists itself as one of its bonded neighbors.
    SelfBondingAtom {
        /// Identifier of the atom that incorrectly lists a self-bond.
        atom_id: usize,
    },
}

impl fmt::Display for GraphValidationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingAtom { atom_id } => {
                write!(f, "bond references a non-existent atom with ID {atom_id}")
            }
            Self::SelfBondingAtom { atom_id } => {
                write!(f, "atom with ID {atom_id} is bonded to itself")
            }
        }
    }
}

impl Error for GraphValidationError {}

/// Errors raised while running the staged chemical perception pipeline.
///
/// Each variant corresponds to a logical section of perception so that downstream callers can
/// attribute failures to the relevant chemical heuristic.
#[derive(Debug)]
pub enum PerceptionError {
    /// No valid Kekulé structure satisfied the aromatic subgraph constraints.
    KekulizationFailed {
        /// Human-readable reason supplied by the Kekulé resolver.
        message: String,
    },

    /// Hybridization inference could not determine an sp/sp2/sp3 class for an atom.
    HybridizationInference {
        /// Identifier of the atom whose steric number could not be mapped.
        atom_id: usize,
    },

    /// Catch-all variant for perception failures that do not fit the other buckets.
    Other(String),
}

impl fmt::Display for PerceptionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::KekulizationFailed { message } => write!(
                f,
                "failed to assign a valid Kekulé structure to the aromatic systems: {message}"
            ),
            Self::HybridizationInference { atom_id } => write!(
                f,
                "could not infer hybridization for atom ID {atom_id}: unhandled steric number or chemical environment"
            ),
            Self::Other(message) => {
                write!(f, "an unexpected perception error occurred: {message}")
            }
        }
    }
}

impl Error for PerceptionError {}

/// Error reported when the typing engine stalls before all atoms receive types.
///
/// This typically indicates that the ruleset lacks coverage for the perceived environments or that
/// earlier perception output was incomplete.
#[derive(Debug)]
pub struct AssignmentError {
    /// Unique identifiers of atoms that never converged to a final type.
    pub untyped_atom_ids: Vec<usize>,
    /// Total number of engine rounds completed before stalling.
    pub rounds_completed: u32,
}

impl fmt::Display for AssignmentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "engine stalled after {} rounds with {:?} still untyped",
            self.rounds_completed, self.untyped_atom_ids
        )
    }
}

impl Error for AssignmentError {}
