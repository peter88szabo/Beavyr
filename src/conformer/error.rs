//! Conformer-search failures and actionable user messages.

use super::BuildError;

/// Why a search could not run or finish.
#[derive(Debug)]
pub enum ConformerError {
    /// DREIDING could not describe the molecule. Its message says which atom and why.
    ForceField(BuildError),
    /// Nothing to search: no rotatable bond, so the molecule has one shape.
    Rigid,
    /// The constraints cannot be used as given.
    Constraints(String),
    /// The molecule's only flexibility is in a ring of a size this search does not pucker.
    ///
    /// Five- and six-membered rings *are* sampled, through their Cremer-Pople puckering
    /// coordinates -- see [`crate::optimizer::conformer_search::rings`]. Other sizes are not:
    /// three- and four-rings are effectively rigid, and a larger ring's conformational space is
    /// not the small sphere or circle that treatment relies on. Since no ring bond can be turned
    /// on its own without prising the ring open, such a molecule has flexibility that goes
    /// unsampled, and saying "one conformer" would be wrong.
    RingOnly {
        rings: usize,
    },
    Search(String),
}

impl std::fmt::Display for ConformerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ForceField(e) => write!(
                f,
                "{}",
                crate::forcefield::dreiding::diagnostics::message(e)
            ),
            Self::Rigid => write!(
                f,
                "this molecule has no rotatable bonds, so it has only one conformer. \
                 Double bonds do not rotate, and a methyl group's rotation gives the same \
                 shape back."
            ),
            Self::RingOnly { rings } => write!(
                f,
                "this molecule's only flexibility is in its {}, and not of a size this search \
                 can pucker. Five- and six-membered rings are sampled -- chair against \
                 twist-boat, or a sugar's puckers -- but a three- or four-ring is effectively \
                 rigid, and a larger ring needs a different treatment. Turning a single ring \
                 bond is not an option: it would prise the ring open rather than give another \
                 conformer.",
                if *rings == 1 { "ring" } else { "rings" }
            ),
            Self::Constraints(message) => write!(f, "the constraints are not usable: {message}"),
            Self::Search(message) => write!(f, "the conformer search failed: {message}"),
        }
    }
}
