//! AFIR: pulling a reaction over its barrier with an artificial force.
//!
//! The method is Maeda and Morokuma's artificial force induced reaction. A constant force is added
//! along the bonds that are about to form, and the geometry is then optimised *normally*. The force
//! drags the reacting atoms together, the optimiser relaxes everything else as it goes, and the
//! resulting sequence of geometries is a reaction pathway: reactant, over the barrier, to product.
//! The highest point on that path, measured on the **true** energy rather than the forced one, is a
//! transition-state-like structure.
//!
//! # Why this and not constraints
//!
//! The obvious way to search a transition state's conformers is to hold its active bonds at fixed
//! lengths and sample the rest. It works, and it has two problems. Holding a coordinate requires
//! redundant internal coordinates, which cost roughly a hundred times a Cartesian cycle -- so a
//! search of hundreds of candidates becomes hours. And every conformer it finds has the *same*
//! active bond length, which is only right if the real transition states share one.
//!
//! AFIR has neither problem. An artificial force is a change to the potential, not a constraint,
//! so the fast Cartesian optimiser runs unchanged. And because each pathway starts from a
//! different reactant conformer at its own long distance, the transition states it reaches have
//! *different* active bond lengths -- which is the point of TStrail (Matsutani, Shomura, Tsuji &
//! Maeda, *J. Comput. Chem.* **2026**, *47*, e70480). Their enolate example found 157 transition
//! state conformers spanning C-C from 1.96 to 2.42 Å; a fixed-length search would have found the
//! ones at whichever length it was given.
//!
//! # What this is not
//!
//! A force field cannot describe a bond forming. DREIDING has fixed bond orders and a harmonic
//! stretch, so it has no barrier to cross and its "pathway" is meaningless. AFIR here is meant to
//! be driven by xTB, as TStrail drives it with GFN2 -- see [`super::refine`] for how xTB is
//! reached. The objective is generic, so anything implementing
//! [`crate::optimizer::traits::Objective`] can be pushed over a barrier by it.

use anyhow::Result;

use crate::optimizer::minimum_cartesian::{
    minimize_cartesian, CartesianMinimumMethod, CartesianMinimumOptions,
};
use crate::optimizer::traits::Objective;

/// A bond that is forming or breaking, named by the two atoms whose distance the force acts on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ActiveBond {
    pub a: usize,
    pub b: usize,
}

/// Artificial force applied along each active bond, in kcal/mol/Å.
///
/// Strong enough to pull two atoms together against the repulsion between them, gentle enough that
/// the rest of the molecule keeps up rather than being dragged through itself. AFIR is usually
/// parametrised by a model collision energy γ instead; this is the force directly, because a force
/// is the thing a chemist can reason about when a pathway comes out wrong -- too weak and it never
/// crosses, too strong and it crushes the structure.
pub const DEFAULT_FORCE_KCAL_PER_ANGSTROM: f64 = 100.0;

/// kcal/mol/Å expressed in Hartree/bohr, the units [`Objective`] works in.
fn to_atomic(force_kcal_per_angstrom: f64) -> f64 {
    use crate::forcefield::dreiding::objective::{BOHR_TO_ANGSTROM, HARTREE_TO_KCAL};
    force_kcal_per_angstrom / HARTREE_TO_KCAL * BOHR_TO_ANGSTROM
}

/// Any objective, with a constant pull added along the active bonds.
///
/// `E_AFIR(x) = E(x) + γ Σ r_ij` over the active bonds. The added term's gradient is exactly
/// `γ` along each bond, which is why this needs no approximation and no extra energy evaluations.
pub struct AfirObjective<'a, O: Objective> {
    inner: &'a mut O,
    active: Vec<ActiveBond>,
    /// Hartree/bohr.
    force: f64,
    /// The true energy at the last point evaluated, before the artificial term was added.
    ///
    /// Kept because the pathway's maximum has to be found on the real surface: the forced energy
    /// falls monotonically as the atoms are pulled together, so its maximum is the starting point
    /// and tells you nothing.
    pub last_true_energy: f64,
}

impl<'a, O: Objective> AfirObjective<'a, O> {
    pub fn new(inner: &'a mut O, active: Vec<ActiveBond>, force_kcal_per_angstrom: f64) -> Self {
        Self {
            inner,
            active,
            force: to_atomic(force_kcal_per_angstrom),
            last_true_energy: 0.0,
        }
    }

    /// The artificial term and its gradient contribution at `x`.
    fn artificial(&self, x: &[f64], grad: Option<&mut [f64]>) -> f64 {
        let mut total = 0.0;
        let mut grad = grad;
        for bond in &self.active {
            let d = [
                x[3 * bond.a] - x[3 * bond.b],
                x[3 * bond.a + 1] - x[3 * bond.b + 1],
                x[3 * bond.a + 2] - x[3 * bond.b + 2],
            ];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if r < 1.0e-12 {
                continue;
            }
            total += self.force * r;
            if let Some(g) = grad.as_deref_mut() {
                for k in 0..3 {
                    let component = self.force * d[k] / r;
                    g[3 * bond.a + k] += component;
                    g[3 * bond.b + k] -= component;
                }
            }
        }
        total
    }
}

impl<O: Objective> Objective for AfirObjective<'_, O> {
    fn energy(&mut self, x: &[f64]) -> Result<f64> {
        let true_energy = self.inner.energy(x)?;
        self.last_true_energy = true_energy;
        Ok(true_energy + self.artificial(x, None))
    }

    fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
        self.inner.gradient(x, grad)?;
        self.artificial(x, Some(grad));
        Ok(())
    }

    fn energy_gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<f64> {
        let true_energy = self.inner.energy_gradient(x, grad)?;
        self.last_true_energy = true_energy;
        Ok(true_energy + self.artificial(x, Some(grad)))
    }
}

/// One point on an AFIR pathway.
#[derive(Debug, Clone)]
pub struct PathPoint {
    /// bohr.
    pub coordinates: Vec<f64>,
    /// The **true** energy in Hartree, without the artificial term.
    pub energy: f64,
    /// The active bond lengths at this point, in Å, in the order they were given.
    pub active_lengths: Vec<f64>,
}

/// A pathway pushed out by an artificial force.
#[derive(Debug, Clone)]
pub struct AfirPath {
    pub points: Vec<PathPoint>,
    /// Index of the highest-energy point: the transition-state-like structure.
    pub maximum: usize,
    /// Whether the path actually rose and fell rather than climbing to its last point.
    ///
    /// A maximum at the final point means the force was still dragging the structure uphill when
    /// the run ended -- so the barrier was not crossed, and the "maximum" is just where it stopped.
    /// Worth reporting rather than passing off as a transition state.
    pub crossed: bool,
}

impl AfirPath {
    /// The transition-state-like structure, in bohr.
    pub fn transition_state_like(&self) -> &PathPoint {
        &self.points[self.maximum]
    }
}

/// Runs one AFIR pathway from `start`, returning the path and its highest point.
///
/// `cycles` caps the forced optimisation. TStrail stops at 30 rather than converging, because the
/// path only has to pass over the barrier -- what happens after that is the product, which a later
/// step handles.
pub fn afir_path<O: Objective>(
    objective: &mut O,
    start: Vec<f64>,
    active: &[ActiveBond],
    force_kcal_per_angstrom: f64,
    cycles: usize,
) -> Result<AfirPath> {
    use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;

    let mut forced = AfirObjective::new(objective, active.to_vec(), force_kcal_per_angstrom);
    let result = minimize_cartesian(
        start,
        &mut forced,
        CartesianMinimumOptions {
            method: CartesianMinimumMethod::Bfgs,
            max_cycles: cycles,
            // The forced surface has no minimum to reach -- the force keeps pulling -- so the
            // convergence tests are left wide and the cycle count is what stops it.
            max_gradient: 0.0,
            rms_gradient: 0.0,
            energy_tol: 0.0,
            verbosity: 0,
            ..CartesianMinimumOptions::default()
        },
    )?;

    // The true energy along the path. The forced energy falls monotonically, so its maximum is
    // the start; the real surface is where the barrier shows.
    let mut points = Vec::with_capacity(result.trajectory.len());
    for frame in &result.trajectory {
        let energy = forced.inner.energy(frame)?;
        let active_lengths = active
            .iter()
            .map(|bond| {
                let d = [
                    frame[3 * bond.a] - frame[3 * bond.b],
                    frame[3 * bond.a + 1] - frame[3 * bond.b + 1],
                    frame[3 * bond.a + 2] - frame[3 * bond.b + 2],
                ];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() * BOHR_TO_ANGSTROM
            })
            .collect();
        points.push(PathPoint {
            coordinates: frame.clone(),
            energy,
            active_lengths,
        });
    }

    if points.is_empty() {
        return Err(anyhow::anyhow!("the AFIR run produced no path points"));
    }
    let maximum = points
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.energy.total_cmp(&b.energy))
        .map(|(index, _)| index)
        .unwrap_or(0);
    // A barrier was crossed only if the path came back down afterwards.
    let crossed = maximum + 1 < points.len() && maximum > 0;

    Ok(AfirPath {
        points,
        maximum,
        crossed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two atoms on a double well along their separation: minima at 2.5 and 5.5 bohr with a
    /// barrier between them at 4.0.
    ///
    /// A Lennard-Jones pair will not do. It has a well and a repulsive wall but no *barrier*, so
    /// there is no transition state to find and no crossing to detect -- which is exactly what
    /// AFIR is for. This is the simplest thing shaped like a reaction.
    struct DoubleWell;

    impl DoubleWell {
        /// Minima at `MIDDLE ± HALF_WIDTH`, barrier at `MIDDLE`.
        const MIDDLE: f64 = 4.0;
        const HALF_WIDTH: f64 = 1.5;
        const STIFFNESS: f64 = 0.002;

        fn separation(x: &[f64]) -> f64 {
            ((x[0] - x[3]).powi(2) + (x[1] - x[4]).powi(2) + (x[2] - x[5]).powi(2)).sqrt()
        }

        fn energy_of(r: f64) -> f64 {
            let offset = (r - Self::MIDDLE).powi(2) - Self::HALF_WIDTH.powi(2);
            Self::STIFFNESS * offset * offset
        }

        fn slope_of(r: f64) -> f64 {
            let offset = (r - Self::MIDDLE).powi(2) - Self::HALF_WIDTH.powi(2);
            4.0 * Self::STIFFNESS * offset * (r - Self::MIDDLE)
        }
    }

    impl Objective for DoubleWell {
        fn energy(&mut self, x: &[f64]) -> Result<f64> {
            Ok(Self::energy_of(Self::separation(x).max(0.2)))
        }

        fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
            let r = Self::separation(x).max(0.2);
            let de_dr = Self::slope_of(r);
            for k in 0..3 {
                let unit = (x[k] - x[3 + k]) / r;
                grad[k] = de_dr * unit;
                grad[3 + k] = -de_dr * unit;
            }
            Ok(())
        }
    }

    /// The forced gradient must match central differences of the forced energy, or the optimiser
    /// walks the wrong way and the whole pathway is meaningless.
    #[test]
    fn the_forced_gradient_matches_finite_differences() {
        let start = vec![0.0, 0.0, 0.0, 5.3, 0.4, -0.2];
        let mut inner = DoubleWell;
        let mut forced = AfirObjective::new(
            &mut inner,
            vec![ActiveBond { a: 0, b: 1 }],
            DEFAULT_FORCE_KCAL_PER_ANGSTROM,
        );

        let mut analytic = vec![0.0; start.len()];
        forced.gradient(&start, &mut analytic).unwrap();

        const H: f64 = 1.0e-6;
        for c in 0..start.len() {
            let mut plus = start.clone();
            let mut minus = start.clone();
            plus[c] += H;
            minus[c] -= H;
            let numeric =
                (forced.energy(&plus).unwrap() - forced.energy(&minus).unwrap()) / (2.0 * H);
            let scale = analytic[c].abs().max(numeric.abs()).max(1.0);
            assert!(
                (analytic[c] - numeric).abs() / scale < 1.0e-5,
                "coordinate {c}: analytic {} vs numeric {numeric}",
                analytic[c]
            );
        }
    }

    /// The artificial term must sum to zero force on the pair as a whole, or the molecule drifts.
    #[test]
    fn the_artificial_force_is_internal() {
        let start = vec![0.0, 0.0, 0.0, 5.3, 0.4, -0.2];
        let mut inner = DoubleWell;
        let mut forced =
            AfirObjective::new(&mut inner, vec![ActiveBond { a: 0, b: 1 }], 100.0);
        let mut grad = vec![0.0; 6];
        forced.gradient(&start, &mut grad).unwrap();
        for axis in 0..3 {
            let net = grad[axis] + grad[3 + axis];
            assert!(net.abs() < 1.0e-12, "net force on axis {axis} is {net}");
        }
    }

    /// The force must actually pull the pair together -- that is its whole job.
    #[test]
    fn the_force_drives_the_pair_together() {
        let start = vec![0.0, 0.0, 0.0, 5.5, 0.0, 0.0];
        let before = DoubleWell::separation(&start);
        let mut inner = DoubleWell;
        let path = afir_path(
            &mut inner,
            start,
            &[ActiveBond { a: 0, b: 1 }],
            DEFAULT_FORCE_KCAL_PER_ANGSTROM,
            30,
        )
        .expect("a path is produced");

        assert!(!path.points.is_empty());
        let after = DoubleWell::separation(&path.points.last().unwrap().coordinates);
        assert!(
            after < before - 0.5,
            "the pair went from {before:.2} to {after:.2} bohr; the force did not pull"
        );
        // Active bond lengths are reported in Angstrom, so they must be smaller than the bohr
        // separation they came from.
        let reported = path.points.last().unwrap().active_lengths[0];
        assert!((reported - after * 0.529_177_210_903).abs() < 1.0e-9);
    }

    /// The maximum must be found on the true energy. The forced energy falls the whole way, so
    /// taking its maximum would always return the starting structure.
    #[test]
    fn the_maximum_is_on_the_true_energy() {
        // Start at the outer minimum and pull inward: over the barrier at 4.0 bohr and down into
        // the inner minimum at 2.5 -- a reaction path with a real transition state on it.
        let start = vec![0.0, 0.0, 0.0, 5.5, 0.0, 0.0];
        let mut inner = DoubleWell;
        let path = afir_path(&mut inner, start, &[ActiveBond { a: 0, b: 1 }], 20.0, 60)
            .expect("a path is produced");

        // The reported maximum really is the largest true energy in the path.
        let largest = path
            .points
            .iter()
            .map(|p| p.energy)
            .fold(f64::NEG_INFINITY, f64::max);
        assert_eq!(path.transition_state_like().energy, largest);
        // And it is not simply the first point, which is what a forced-energy maximum would give.
        assert!(
            path.maximum > 0,
            "the maximum landed on the starting structure"
        );
        // The barrier of this surface sits at 4.0 bohr, so that is where the maximum belongs.
        let at_maximum = DoubleWell::separation(&path.transition_state_like().coordinates);
        assert!(
            (at_maximum - DoubleWell::MIDDLE).abs() < 0.6,
            "the maximum is at {at_maximum:.2} bohr; the barrier is at {:.1}",
            DoubleWell::MIDDLE
        );
        // Having gone over and come down, the path must report a crossing.
        assert!(path.crossed, "the path crossed the barrier but did not say so");
    }

    /// A run that never came back down must say so rather than offer its last structure as a
    /// transition state.
    #[test]
    fn a_path_that_only_climbs_is_not_marked_as_crossed() {
        // A short run under a strong force: it is still being pushed uphill when it stops.
        let start = vec![0.0, 0.0, 0.0, 5.0, 0.0, 0.0];
        let mut inner = DoubleWell;
        let path = afir_path(&mut inner, start, &[ActiveBond { a: 0, b: 1 }], 20.0, 2)
            .expect("a path is produced");
        if path.maximum + 1 == path.points.len() {
            assert!(!path.crossed, "a still-climbing path must not claim a crossing");
        }
    }
}
