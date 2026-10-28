//! Cooperative cancellation at force-field evaluation boundaries.

use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{bail, Result};

use crate::optimizer::traits::Objective;

pub(crate) fn check_cancelled(cancel: Option<&AtomicBool>) -> Result<()> {
    if cancel.is_some_and(|flag| flag.load(Ordering::SeqCst)) {
        bail!("Cancelled.");
    }
    Ok(())
}

/// Keeps cancellation out of the force-field equations, while covering every
/// optimizer step, line-search evaluation and finite-difference displacement.
pub(crate) struct CancellableObjective<'a, O> {
    inner: O,
    cancel: Option<&'a AtomicBool>,
}

impl<'a, O> CancellableObjective<'a, O> {
    pub(crate) fn new(inner: O, cancel: Option<&'a AtomicBool>) -> Self {
        Self { inner, cancel }
    }
}

impl<O: Objective> Objective for CancellableObjective<'_, O> {
    fn energy(&mut self, x: &[f64]) -> Result<f64> {
        check_cancelled(self.cancel)?;
        let energy = self.inner.energy(x)?;
        check_cancelled(self.cancel)?;
        Ok(energy)
    }

    fn gradient(&mut self, x: &[f64], gradient: &mut [f64]) -> Result<()> {
        check_cancelled(self.cancel)?;
        self.inner.gradient(x, gradient)?;
        check_cancelled(self.cancel)
    }

    fn energy_gradient(&mut self, x: &[f64], gradient: &mut [f64]) -> Result<f64> {
        check_cancelled(self.cancel)?;
        let energy = self.inner.energy_gradient(x, gradient)?;
        check_cancelled(self.cancel)?;
        Ok(energy)
    }
}

#[cfg(test)]
#[path = "cancellation_tests.rs"]
mod tests;
