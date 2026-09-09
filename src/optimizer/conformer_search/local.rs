//! Which optimiser relaxes each candidate conformer.
//!
//! A conformer search is thousands of small local optimisations, so this choice sets the cost of
//! the whole run -- and the right answer depends entirely on what an energy costs.
//!
//! [`LocalOptimizer::RedundantInternal`] is Behemoth's, and the right choice for quantum
//! chemistry: redundant internal coordinates follow bond stretches and torsions directly, so a
//! structure converges in few cycles. Each cycle pays for it, rebuilding the Wilson B matrix and
//! inverting the `nint × nint` G matrix. When an energy takes seconds, buying fewer cycles at
//! almost any per-cycle price is obviously right.
//!
//! [`LocalOptimizer::Cartesian`] inverts that trade, and is the right choice for a force field.
//! With DREIDING an energy and gradient cost microseconds, so the coordinate algebra dominates
//! completely: for a 32-atom molecule with ~150 internal coordinates, one cycle spends roughly
//! 6.7 Mflop on B and G against 0.01 Mflop on the force field -- about 99.8% of the run is spent
//! deciding where to step rather than working out the energy. Cartesian BFGS keeps its Hessian in
//! the 3N Cartesian space, so a cycle costs O((3N)²) with no matrix inversion, and needing more
//! cycles is a bargain at that price.

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use anyhow::Result;

use crate::optimizer::geom_opt::{geom_opt_internal_bfgs, GeomOptOptions, GeomOptResult};
use crate::optimizer::internal_coords::ConnectivityModel;
use crate::optimizer::minimum_cartesian::{
    minimize_cartesian, CartesianMinimumMethod, CartesianMinimumOptions,
};
use crate::optimizer::traits::Objective;

/// How each candidate conformer is relaxed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LocalOptimizer {
    /// BFGS in redundant internal coordinates. Fewest cycles, most work per cycle.
    ///
    /// Behemoth's own choice, and the default here so that importing the search does not silently
    /// change its behaviour.
    #[default]
    RedundantInternal,
    /// BFGS in Cartesian coordinates. More cycles, far less work per cycle.
    Cartesian,
}

impl LocalOptimizer {
    pub const ALL: [LocalOptimizer; 2] =
        [LocalOptimizer::Cartesian, LocalOptimizer::RedundantInternal];

    pub fn label(self) -> &'static str {
        match self {
            Self::Cartesian => "Cartesian BFGS",
            Self::RedundantInternal => "Internal coordinates",
        }
    }

    pub fn description(self) -> &'static str {
        match self {
            Self::Cartesian => {
                "Far faster with a force field, where an energy is nearly free and the \
                 coordinate algebra is the whole cost."
            }
            Self::RedundantInternal => {
                "Converges in fewer cycles but pays a matrix inversion for each one. The right \
                 choice when an energy is expensive."
            }
        }
    }
}

/// Total time spent relaxing candidates, in nanoseconds, and how many were relaxed.
///
/// A search is local optimisation plus bookkeeping -- generating trial torsions, rejecting
/// nonsense geometries, screening against the blacklist by RMSD -- and knowing the split is what
/// says whether making the optimiser faster is still worth anything. Two atomic adds per
/// optimisation costs nothing next to the optimisation itself.
static SPENT_NANOS: AtomicU64 = AtomicU64::new(0);
static COUNT: AtomicU64 = AtomicU64::new(0);

/// Time spent in local optimisation, and the number of optimisations, since [`reset_accounting`].
pub fn accounting() -> (std::time::Duration, u64) {
    (
        std::time::Duration::from_nanos(SPENT_NANOS.load(Ordering::Relaxed)),
        COUNT.load(Ordering::Relaxed),
    )
}

/// Zeroes the accounting, so one search can be measured without the previous one's total.
pub fn reset_accounting() {
    SPENT_NANOS.store(0, Ordering::Relaxed);
    COUNT.store(0, Ordering::Relaxed);
}

/// Relaxes one candidate, by whichever route was chosen.
pub fn local_optimize<O: Objective>(
    method: LocalOptimizer,
    coordinates: Vec<f64>,
    connectivity: ConnectivityModel,
    objective: &mut O,
    options: GeomOptOptions,
) -> Result<GeomOptResult> {
    let started = Instant::now();
    let result = dispatch(method, coordinates, connectivity, objective, options);
    SPENT_NANOS.fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
    COUNT.fetch_add(1, Ordering::Relaxed);
    result
}

fn dispatch<O: Objective>(
    method: LocalOptimizer,
    coordinates: Vec<f64>,
    connectivity: ConnectivityModel,
    objective: &mut O,
    options: GeomOptOptions,
) -> Result<GeomOptResult> {
    match method {
        LocalOptimizer::RedundantInternal => {
            geom_opt_internal_bfgs(coordinates, connectivity, objective, options)
        }
        LocalOptimizer::Cartesian => {
            cartesian(CartesianMinimumMethod::Bfgs, coordinates, objective, options)
        }
    }
}

/// The Cartesian route, presented as a [`GeomOptResult`] so the search need not care which ran.
fn cartesian<O: Objective>(
    method: CartesianMinimumMethod,
    coordinates: Vec<f64>,
    objective: &mut O,
    options: GeomOptOptions,
) -> Result<GeomOptResult> {
    let ndim = coordinates.len();
    let result = minimize_cartesian(
        coordinates,
        objective,
        CartesianMinimumOptions {
            method,
            // The convergence criteria carry across directly; they are stated on the geometry and
            // the gradient, neither of which depends on the coordinate system.
            max_cycles: options.max_cycles,
            energy_tol: options.energy_tol,
            max_step: options.max_step,
            rms_step: options.rms_step,
            max_gradient: options.max_gradient,
            rms_gradient: options.g_rms_tol,
            // `max_step_internal` limits a step in internal coordinates, where a component is an
            // angle as readily as a length. The Cartesian equivalent is a displacement, so the
            // optimiser's own default applies rather than a reinterpreted number.
            step_max_component: CartesianMinimumOptions::default().step_max_component,
            min_line_search_scale: CartesianMinimumOptions::default().min_line_search_scale,
            verbosity: options.verbosity,
        },
    )?;

    let grad_rms = if ndim == 0 {
        0.0
    } else {
        (result.gradient.iter().map(|g| g * g).sum::<f64>() / ndim as f64).sqrt()
    };

    Ok(GeomOptResult {
        x: result.x,
        energy: result.energy,
        gradient: result.gradient,
        grad_rms,
        cycles: result.cycles,
        converged: result.converged,
        trajectory: result.trajectory,
        trajectory_energies: result.trajectory_energies,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A quadratic bowl in 6 coordinates, minimum at the origin. Enough to check that both routes
    /// reach the same place and report it the same way.
    struct Bowl {
        calls: usize,
    }

    impl Objective for Bowl {
        fn energy(&mut self, x: &[f64]) -> Result<f64> {
            self.calls += 1;
            Ok(x.iter().map(|v| v * v).sum())
        }

        fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
            for (g, v) in grad.iter_mut().zip(x) {
                *g = 2.0 * v;
            }
            Ok(())
        }
    }

    fn options() -> GeomOptOptions {
        GeomOptOptions {
            verbosity: 0,
            ..GeomOptOptions::default()
        }
    }

    #[test]
    fn the_cartesian_route_finds_the_minimum() {
        let mut bowl = Bowl { calls: 0 };
        let start = vec![0.3, -0.2, 0.1, 0.25, -0.15, 0.05];
        let result = cartesian(CartesianMinimumMethod::Bfgs, start, &mut bowl, options())
            .expect("a bowl is minimisable");

        assert!(result.converged, "did not converge");
        for value in &result.x {
            assert!(value.abs() < 1.0e-3, "ended at {:?}", result.x);
        }
        assert!(result.energy < 1.0e-6);
        assert!(result.grad_rms < 1.0e-3);
    }

    /// `grad_rms` must be the root-mean-square, not the norm: the search compares it against
    /// `g_rms_tol`, so a factor of sqrt(N) would silently change when a conformer is accepted.
    #[test]
    fn grad_rms_is_a_root_mean_square() {
        struct Fixed;
        impl Objective for Fixed {
            fn energy(&mut self, _x: &[f64]) -> Result<f64> {
                Ok(0.0)
            }
            fn gradient(&mut self, _x: &[f64], grad: &mut [f64]) -> Result<()> {
                grad.iter_mut().for_each(|g| *g = 2.0);
                Ok(())
            }
        }
        // Six components of 2.0 -- two atoms, since the optimiser requires a 3N vector. The RMS
        // is 2.0 where the norm would be about 4.9.
        let result =
            cartesian(CartesianMinimumMethod::Bfgs, vec![0.0; 6], &mut Fixed, options()).unwrap();
        assert!(
            (result.grad_rms - 2.0).abs() < 1.0e-12,
            "grad_rms is {}",
            result.grad_rms
        );
    }

    #[test]
    fn both_routes_are_offered_and_described() {
        assert_eq!(LocalOptimizer::ALL.len(), 2);
        for method in LocalOptimizer::ALL {
            assert!(!method.label().is_empty());
            assert!(!method.description().is_empty());
        }
        // Behemoth's behaviour is the default, so importing the search changes nothing by itself.
        assert_eq!(LocalOptimizer::default(), LocalOptimizer::RedundantInternal);
    }
}
