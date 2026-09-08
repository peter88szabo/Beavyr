//! Geometry optimization, imported from Behemoth.
//!
//! # Provenance
//!
//! Copied from `~/Dropbox/Research_Leuven/Behemoth/src/optimizer/`. Behemoth's optimizer was
//! written to be independent of its quantum chemistry -- its own words: "any backend that can
//! provide `E(x)` and `∇E(x)` can be plugged in" -- so it transfers here essentially unchanged,
//! and Beavyr gets one tested optimiser rather than a second, weaker one.
//!
//! The import is BLAS-free: these modules need only `anyhow` and `ndarray`, both of which Beavyr
//! already has. Behemoth's `ndarray-linalg` is not required.
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

/// Cartesian minimisation: BFGS, steepest descent, and two-point gradient methods.
pub mod minimum_cartesian;

/// Central-difference gradients, for objectives that supply only energies.
pub mod finite_diff;
