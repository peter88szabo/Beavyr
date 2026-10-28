//! User guidance when DREIDING cannot describe a system.

use super::BuildError;

pub(crate) const UNSUPPORTED_SYSTEM: &str = "DREIDING cannot handle this system.";

pub(crate) fn message(error: &BuildError) -> String {
    match error {
        BuildError::Typing(_) | BuildError::Unparameterized { .. } => format!(
            "{UNSUPPORTED_SYSTEM}\n\
             Open the Quantum Chemistry panel, select xTB, and run a geometry optimization. \
             Set the correct total charge and spin multiplicity for your molecule before starting.\n\
             Details: {error}"
        ),
        _ => error.to_string(),
    }
}
