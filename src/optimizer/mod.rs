// The module is imported whole and deliberately: transition-state search, RDA with its
// poor-man's nudged elastic band, IRC following, spin-crossing optimisation and the conformer
// search's parallel and reporting entry points are all here awaiting panels of their own. Keeping
// the public surface intact is what makes re-vendoring from Behemoth a copy rather than a merge,
// so unused items are expected rather than a sign of dead code.
#![allow(dead_code, unused_imports)]

//! Geometry optimization, imported from Behemoth.
//!
//! # Provenance
//!
//! Copied from `~/Dropbox/Research_Leuven/Behemoth/src/optimizer/`. Behemoth's optimizer was
//! written to be independent of its quantum chemistry -- its own words: "any backend that can
//! provide `E(x)` and `∇E(x)` can be plugged in" -- so it transfers here essentially unchanged,
//! and Beavyr gets one tested optimiser rather than a second, weaker one.
//!
//! The import is BLAS-free. These modules need only `anyhow` and `ndarray`, both of which Beavyr
//! already has; Behemoth's `ndarray-linalg` is not required. Where they reached outside the
//! optimizer, four small redirections cover it:
//!
//! | Behemoth | here |
//! |---|---|
//! | `numeric::linalg::{Backend, LinAlg}` | `normalmode::linalg_shim`, the pure-Rust Jacobi eigensolver |
//! | `numeric::linalg::solve_linear_system_array2` | [`linalg_solve`], Behemoth's own Gauss-Jordan, which needed no BLAS |
//! | `utils::atomic_masses::AtomicMasses` | [`atomic_masses`], forwarding to Beavyr's mass table |
//! | `rand` | [`prng`], a seedable xoshiro256** in the tree |

//!
//! Two things were deliberately left behind. `xtb_gfn1` and `tasi_eht` are the quantum-chemistry
//! engine rather than optimiser logic, and Beavyr reaches xTB and Behemoth as subprocesses
//! instead. The conformer search's one-thread Rayon pool went too: it exists so that Rayon inside
//! an electronic-structure method cannot steal the other workers' threads, and Beavyr's force
//! field uses neither Rayon nor BLAS.
//!
//! # Conventions
//!
//! [`traits::Objective`] fixes the units, and they are **not** Beavyr's usual ones:
//!
//! | quantity | unit |
//! |---|---|
//! | coordinates | bohr |
//! | energy | Hartree |
//! | gradient | Hartree/bohr |
//!
//! Beavyr works in Å throughout, and the DREIDING force field in kcal/mol, so anything plugged in
//! here converts at the boundary. See `forcefield::dreiding::objective`.

/// Backend-agnostic energy and gradient interface.
pub mod traits;

/// Central-difference gradients, for objectives that supply only energies.
pub mod finite_diff;

/// Cartesian minimisation: BFGS, steepest descent, and two-point gradient methods.
pub mod minimum_cartesian;

/// Redundant internal coordinates, and the B matrix relating them to Cartesians.
pub mod internal_coords;

/// Building the redundant internal coordinate set from connectivity.
pub mod redundant_internals;

/// Minimisation in redundant internal coordinates, which converges faster on torsions.
pub mod minimum_redundant;

/// The general geometry-optimisation driver.
pub mod geom_opt;

/// Optimisation with geometric constraints held fixed.
pub mod constrained;

/// Transition-state search.
pub mod transition_state;

/// Reaction-driven adiabatic mapping, including the poor-man's nudged elastic band, used to
/// generate an initial guess for a transition-state search.
pub mod rda;

/// Intrinsic reaction coordinate following.
pub mod irc;

/// Minimum-energy crossing-point optimisation between two spin states.
pub mod spin_crossing;

/// Conformer search by genetic algorithm over the rotatable torsions.
pub mod conformer_search;

/// Dense Gauss-Jordan solve, imported alongside the optimizer.
pub mod linalg_solve;


/// Atomic masses, forwarding to Beavyr's own table.
pub mod atomic_masses;

/// Seedable pseudo-random generator, in place of the `rand` crate.
pub mod prng;

/// Rigid bond and angle scans.
pub mod scan;

/// The slice of Behemoth's `Molecule` that the scan code needs.
pub mod scan_molecule;

// Behemoth's own re-export surface, so imported call sites and future code use the same paths.
// `analytic_gradient`, `tasi_eht` and `xtb_gfn1` are absent: those bind to Behemoth's
// quantum-chemistry engine. `scan` is handled separately, being partly bound to it as well.
pub use conformer_search::{
    discover_torsions, genetic_conformer_search, genetic_conformer_search_parallel,
    print_conformer_search_result, write_conformer_search_files, Conformer, ConformerSearchOptions,
    ConformerSearchResult, ConformerSearchStatistics, PreparedTorsion, SelectionMode,
    TerminationReason, TorsionKind,
};
pub use constrained::{
    optimize_constrained_internals, optimize_constrained_with_internals, relaxed_constraint_scan,
    ConstrainedCoordinateMode, ConstrainedOptimizationOptions, ConstrainedOptimizationResult,
    ConstraintCoordinate, ConstraintTarget, RelaxedScanOptions, RelaxedScanPoint,
    RelaxedScanResult,
};
pub use finite_diff::FiniteDifference;
pub use geom_opt::{
    geom_opt_internal_bfgs, geom_opt_internal_bfgs_with_internals, GeomOptOptions, GeomOptResult,
};
pub use internal_coords::internal_coords_from_bonds;
pub use irc::{
    follow_irc, IrcBranch, IrcCoordinateMode, IrcOptions, IrcPoint, IrcPredictor, IrcResult,
};
pub use minimum_cartesian::{
    minimize_cartesian, CartesianMinimumMethod, CartesianMinimumOptions, CartesianMinimumResult,
};
pub use minimum_redundant::{
    minimize_internal, InitialInternalHessianModel, InternalCoordinateMode, InternalMinimumMethod,
    InternalMinimumOptions, InternalMinimumResult,
};
pub use rda::{
    generate_rda_ts_guess, RdaConditionalCoordinates, RdaDistanceMode, RdaGuessResult, RdaOptions,
    RdaPoint,
};
pub use redundant_internals::{build_redundant_internals, RedundantInternalSystem};
pub use scan::{
    bond_angle_degrees, bond_distance_angstrom, rigid_angle_scan, rigid_bond_scan,
    AngleScanOptions, AngleScanPoint, AngleScanResult, BondScanMoveMode, BondScanOptions,
    BondScanPoint, BondScanResult,
};
pub use scan_molecule::{ScanAtom, ScanMolecule};
pub use spin_crossing::{
    optimize_spin_crossing, SpinCrossingCoordinateMode, SpinCrossingMethod, SpinCrossingOptions,
    SpinCrossingPoint, SpinCrossingResult,
};
pub use traits::Objective;
pub use transition_state::{
    optimize_transition_state, InternalReferenceKind, TSModeDiagnostics,
    TransitionStateCoordinateMode, TransitionStateInitialHessian, TransitionStateModeTracking,
    TransitionStateOptions, TransitionStateReactionMode, TransitionStateResult,
    WeightedInternalReference,
};
