//! Finite-difference gradients for objectives that only provide energies.
//!
//! This is the “trick” that keeps the optimizer usable even before analytic
//! gradients are implemented in the Hamiltonians.

use anyhow::{bail, Result};

use super::traits::Objective;

/// Wrap an `Objective` and provide a numerical gradient via central differences.
#[derive(Debug)]
pub struct FiniteDifference<'a, O: Objective> {
    pub inner: &'a mut O,
    pub step: f64,
}

impl<'a, O: Objective> FiniteDifference<'a, O> {
    pub fn new(inner: &'a mut O, step: f64) -> Self {
        Self { inner, step }
    }
}

impl<'a, O: Objective> Objective for FiniteDifference<'a, O> {
    fn energy(&mut self, x: &[f64]) -> Result<f64> {
        self.inner.energy(x)
    }

    fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
        if self.step <= 0.0 {
            bail!("finite difference step must be > 0");
        }
        if grad.len() != x.len() {
            bail!(
                "gradient length mismatch: grad={} x={}",
                grad.len(),
                x.len()
            );
        }

        let h = self.step;
        let mut xp = x.to_vec();
        let mut xm = x.to_vec();

        for i in 0..x.len() {
            xp[i] = x[i] + h;
            xm[i] = x[i] - h;
            let ep = self.inner.energy(&xp)?;
            let em = self.inner.energy(&xm)?;
            grad[i] = (ep - em) / (2.0 * h);
            xp[i] = x[i];
            xm[i] = x[i];
        }
        Ok(())
    }
}
