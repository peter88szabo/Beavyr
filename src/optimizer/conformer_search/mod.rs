//! First-principles conformer search with a torsion-space genetic algorithm.
//!
//! The implementation follows the Lamarckian search described by Supady,
//! Blum, and Baldauf (J. Chem. Inf. Model. 2015, 55, 2338): trial torsion
//! vectors are screened before a full local geometry optimization, and both
//! starting structures and optimized minima are blacklisted.

mod genetic;
pub(crate) mod local;
mod geometry;
mod options;
mod output;
mod parallel;
mod rmsd;
mod topology;

pub use genetic::{
    genetic_conformer_search, Conformer, ConformerSearchResult, ConformerSearchStatistics,
    TerminationReason,
};
pub use local::LocalOptimizer;
pub use options::{ConformerSearchOptions, SelectionMode};
pub use output::{print_conformer_search_result, write_conformer_search_files};
pub use parallel::genetic_conformer_search_parallel;
pub use topology::{discover_torsions, PreparedTorsion, TorsionKind};

#[cfg(test)]
mod tests;
