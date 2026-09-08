//! Bridges the DREIDING force field to Beavyr's imported optimizer.
//!
//! The two speak different units. DREIDING is defined in kcal/mol and Å (the paper's units, and
//! the ones [`super::terms`] works in); [`crate::optimizer::traits::Objective`] is Behemoth's
//! interface and fixes bohr, Hartree and Hartree/bohr. All conversion happens here, at the
//! boundary, so neither side has to know about the other's convention.

use anyhow::Result;

use crate::optimizer::minimum_cartesian::{
    minimize_cartesian, CartesianMinimumMethod, CartesianMinimumOptions,
};
use crate::optimizer::traits::Objective;

use super::{DreidingTopology, EnergyBreakdown, Workspace};

/// Å per bohr.
pub const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;

/// kcal/mol per Hartree.
pub const HARTREE_TO_KCAL: f64 = 627.509_474_063_1;

/// A DREIDING topology presented as an optimisable objective.
///
/// Holds the [`Workspace`] because [`Objective`] takes `&mut self`; the topology itself stays
/// borrowed and immutable, so several objectives can share one topology across threads.
pub struct DreidingObjective<'a> {
    topology: &'a DreidingTopology,
    work: Workspace,
    /// Scratch for the kcal/mol/Å gradient, so no allocation happens per evaluation.
    scratch: Vec<f64>,
    /// Positions in Å, likewise reused.
    angstrom: Vec<f64>,
}

impl<'a> DreidingObjective<'a> {
    pub fn new(topology: &'a DreidingTopology) -> Self {
        let n = topology.natoms() * 3;
        Self {
            topology,
            work: Workspace::new(),
            scratch: vec![0.0; n],
            angstrom: vec![0.0; n],
        }
    }

    /// Energy in kcal/mol and `dE/dx` in kcal/mol/Å, from coordinates in bohr.
    fn evaluate(&mut self, x_bohr: &[f64]) -> f64 {
        for (a, b) in self.angstrom.iter_mut().zip(x_bohr) {
            *a = b * BOHR_TO_ANGSTROM;
        }
        // `energy_and_forces` returns forces; the optimizer wants the gradient.
        let energy = self
            .topology
            .energy_and_forces(&self.angstrom, &mut self.work, &mut self.scratch);
        for g in self.scratch.iter_mut() {
            *g = -*g;
        }
        energy
    }
}

impl Objective for DreidingObjective<'_> {
    fn energy(&mut self, x: &[f64]) -> Result<f64> {
        Ok(self.evaluate(x) / HARTREE_TO_KCAL)
    }

    fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
        self.evaluate(x);
        // dE[Ha]/dx[bohr] = dE[kcal]/dx[Å] × (Å/bohr) / (kcal/Ha)
        let factor = BOHR_TO_ANGSTROM / HARTREE_TO_KCAL;
        for (g, s) in grad.iter_mut().zip(&self.scratch) {
            *g = s * factor;
        }
        Ok(())
    }

    fn energy_gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<f64> {
        let energy = self.evaluate(x);
        let factor = BOHR_TO_ANGSTROM / HARTREE_TO_KCAL;
        for (g, s) in grad.iter_mut().zip(&self.scratch) {
            *g = s * factor;
        }
        Ok(energy / HARTREE_TO_KCAL)
    }
}

/// What a cleanup did, in the units a chemist reads.
#[derive(Debug, Clone)]
pub struct Relaxation {
    pub cycles: usize,
    /// kcal/mol.
    pub initial_energy: f64,
    /// kcal/mol.
    pub final_energy: f64,
    pub breakdown: EnergyBreakdown,
    /// Largest remaining force on any atom, kcal/mol/Å.
    pub max_force: f64,
    pub converged: bool,
    pub message: String,
    /// Every accepted geometry, in Å, for the optimisation plot and trajectory view.
    pub trajectory: Vec<Vec<f64>>,
    /// kcal/mol, one per trajectory frame.
    pub trajectory_energies: Vec<f64>,
}

/// Convergence suited to cleaning up a sketch rather than to publishing a geometry.
///
/// Deliberately looser than Behemoth's quantum-chemistry defaults: this is a structure to hand to
/// a real optimiser, and the last decimal costs iterations the user waits through. `max_cycles` is
/// raised instead, because a sketched structure can start much further from a minimum than a
/// quantum-chemistry job ever does.
pub fn cleanup_options() -> CartesianMinimumOptions {
    CartesianMinimumOptions {
        method: CartesianMinimumMethod::Bfgs,
        max_cycles: 400,
        // 4e-4 Ha/bohr is about 0.47 kcal/mol/Å.
        max_gradient: 4.0e-4,
        rms_gradient: 2.0e-4,
        // Silent: this runs behind a button, not a terminal.
        verbosity: 0,
        ..CartesianMinimumOptions::default()
    }
}

/// Relaxes a geometry given in Å, returning the relaxed coordinates and what happened.
///
/// This is the entry point the Fragment Editor's cleanup button uses.
pub fn relax_angstrom(
    topology: &DreidingTopology,
    pos_angstrom: &[f64],
    options: CartesianMinimumOptions,
) -> Result<(Vec<f64>, Relaxation)> {
    let mut objective = DreidingObjective::new(topology);

    let initial_energy = objective.evaluate(
        &pos_angstrom
            .iter()
            .map(|a| a / BOHR_TO_ANGSTROM)
            .collect::<Vec<_>>(),
    );

    let x_bohr: Vec<f64> = pos_angstrom.iter().map(|a| a / BOHR_TO_ANGSTROM).collect();
    let result = minimize_cartesian(x_bohr, &mut objective, options)?;

    let relaxed: Vec<f64> = result.x.iter().map(|b| b * BOHR_TO_ANGSTROM).collect();

    // Report forces the way the paper does, not in Hartree/bohr.
    let gradient_factor = HARTREE_TO_KCAL / BOHR_TO_ANGSTROM;
    let max_force = result
        .gradient
        .iter()
        .fold(0.0_f64, |acc, g| acc.max(g.abs()))
        * gradient_factor;

    let mut work = Workspace::new();
    let breakdown = topology.energy_breakdown(&relaxed, &mut work);

    let relaxation = Relaxation {
        cycles: result.cycles,
        initial_energy,
        final_energy: result.energy * HARTREE_TO_KCAL,
        breakdown,
        max_force,
        converged: result.converged,
        message: result.message,
        trajectory: result
            .trajectory
            .iter()
            .map(|frame| frame.iter().map(|b| b * BOHR_TO_ANGSTROM).collect())
            .collect(),
        trajectory_energies: result
            .trajectory_energies
            .iter()
            .map(|e| e * HARTREE_TO_KCAL)
            .collect(),
    };
    Ok((relaxed, relaxation))
}
