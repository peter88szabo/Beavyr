//! Optimizer interfaces.
//!
//! The intent is to keep this crate’s geometry optimization independent from
//! the quantum-chemistry code: any backend that can provide `E(x)` and `∇E(x)`
//! can be plugged in.

use anyhow::{bail, Result};

/// Objective function interface for geometry optimization.
///
/// Conventions:
/// - `x` is a flat vector of length `3*nat`, storing Cartesian coordinates in **bohr**.
/// - Energies are in **Hartree**.
/// - Gradients are in **Hartree/bohr**.
pub trait Objective {
    /// Compute the energy `E(x)`.
    fn energy(&mut self, x: &[f64]) -> Result<f64>;

    /// Compute the gradient `∇E(x)` into `grad`.
    ///
    /// Default: returns an error. Use `FiniteDifference` if you only have energies.
    fn gradient(&mut self, _x: &[f64], _grad: &mut [f64]) -> Result<()> {
        bail!("Objective::gradient not implemented")
    }

    /// Convenience: compute energy and gradient in one call.
    fn energy_gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<f64> {
        let e = self.energy(x)?;
        self.gradient(x, grad)?;
        Ok(e)
    }
}
