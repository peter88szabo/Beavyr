use super::*;
use std::cell::Cell;

/// Requests cancellation from inside a running evaluation, without sleeps or
/// racing another thread. The optimizer and Hessian must propagate that stop.
struct Quadratic<'a> {
    cancel: &'a AtomicBool,
    evaluations: &'a Cell<usize>,
    stop_after: usize,
}

impl Objective for Quadratic<'_> {
    fn energy(&mut self, x: &[f64]) -> Result<f64> {
        let count = self.evaluations.get() + 1;
        self.evaluations.set(count);
        if count == self.stop_after {
            self.cancel.store(true, Ordering::SeqCst);
        }
        Ok(x.iter().map(|value| value * value / 2.0).sum())
    }

    fn gradient(&mut self, x: &[f64], gradient: &mut [f64]) -> Result<()> {
        self.energy_gradient(x, gradient).map(|_| ())
    }

    fn energy_gradient(&mut self, x: &[f64], gradient: &mut [f64]) -> Result<f64> {
        let energy = self.energy(x)?;
        gradient.copy_from_slice(x);
        Ok(energy)
    }
}

#[test]
fn cancellation_during_each_evaluation_discards_its_result() {
    for operation in 0..3 {
        let cancel = AtomicBool::new(false);
        let evaluations = Cell::new(0);
        let inner = Quadratic {
            cancel: &cancel,
            evaluations: &evaluations,
            stop_after: 1,
        };
        let mut objective = CancellableObjective::new(inner, Some(&cancel));
        let x = [4.0, 3.0, 2.0];
        let mut gradient = [0.0; 3];
        let result = match operation {
            0 => objective.energy(&x).map(|_| ()),
            1 => objective.gradient(&x, &mut gradient),
            _ => objective.energy_gradient(&x, &mut gradient).map(|_| ()),
        };
        assert_eq!(result.unwrap_err().to_string(), "Cancelled.");
        assert_eq!(evaluations.get(), 1);
        assert!(objective.energy(&x).is_err());
        assert_eq!(
            evaluations.get(),
            1,
            "do not start another evaluation after cancellation"
        );
    }
}

#[test]
fn cancellation_stops_an_optimizer_that_is_already_running() {
    use crate::optimizer::minimum_cartesian::{minimize_cartesian, CartesianMinimumOptions};
    let cancel = AtomicBool::new(false);
    let evaluations = Cell::new(0);
    let inner = Quadratic {
        cancel: &cancel,
        evaluations: &evaluations,
        stop_after: 4,
    };
    let mut objective = CancellableObjective::new(inner, Some(&cancel));
    let result = minimize_cartesian(
        vec![4.0, 3.0, 2.0],
        &mut objective,
        CartesianMinimumOptions {
            verbosity: 0,
            ..Default::default()
        },
    );
    assert_eq!(result.unwrap_err().to_string(), "Cancelled.");
    assert_eq!(evaluations.get(), 4);
}

#[test]
fn cancellation_stops_a_hessian_partway_through_its_displacements() {
    use crate::normalmode::fd_hessian::numerical_hessian_from_gradients;
    let cancel = AtomicBool::new(false);
    let evaluations = Cell::new(0);
    let inner = Quadratic {
        cancel: &cancel,
        evaluations: &evaluations,
        stop_after: 4,
    };
    let mut objective = CancellableObjective::new(inner, Some(&cancel));
    let result = numerical_hessian_from_gradients(&[4.0, 3.0, 2.0], 0.005, 0, |x| {
        let mut gradient = vec![0.0; x.len()];
        objective.gradient(x, &mut gradient)?;
        Ok(gradient)
    });
    assert_eq!(result.unwrap_err().to_string(), "Cancelled.");
    assert_eq!(
        evaluations.get(),
        4,
        "must not evaluate all six displacements"
    );
}

#[test]
fn dreiding_optimization_and_hessian_entry_points_honor_cancellation() {
    use crate::forcefield::dreiding::objective::{
        cleanup_options, relax_angstrom_cancellable, BOHR_TO_ANGSTROM,
    };
    use crate::forcefield::dreiding::{
        hessian::hessian_and_gradient_cancellable, DreidingTopology,
    };
    let mol = crate::qchem_interfaces::dreiding_run::molecule_for_force_field(
        "3\n\nO 0 0 0\nH 1.3 0 0\nH 0 1.3 0\n",
    );
    let topology = DreidingTopology::build(&mol).unwrap();
    let positions: Vec<f64> = mol
        .pos
        .iter()
        .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
        .collect();
    let cancel = AtomicBool::new(true);
    assert_eq!(
        relax_angstrom_cancellable(&topology, &positions, cleanup_options(), Some(&cancel))
            .unwrap_err()
            .to_string(),
        "Cancelled."
    );
    let bohr: Vec<_> = positions.iter().map(|v| v / BOHR_TO_ANGSTROM).collect();
    assert_eq!(
        hessian_and_gradient_cancellable(&topology, &bohr, Some(&cancel))
            .unwrap_err()
            .to_string(),
        "Cancelled."
    );
}
