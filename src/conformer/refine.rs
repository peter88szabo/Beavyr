//! Re-ranking conformers with xTB's GFN2.
//!
//! The search itself has to run on DREIDING: it needs tens of thousands of energies, and an
//! external program costs a process launch for each. Re-ranking is the opposite shape of problem.
//! Only the conformers that survived are re-optimised -- a few dozen, a few hundred xTB calls,
//! minutes -- and that is worth doing, because relative conformer energies are precisely where a
//! generic force field is weakest. DREIDING is reliable about which shapes *exist* and much less
//! so about their order.
//!
//! GFN2-xTB does the re-ranking. It is a semiempirical tight-binding method with a self-consistent
//! charge treatment, so unlike any force field it accounts for polarisation and for dispersion
//! from the electronic structure rather than from fitted pair terms. That is exactly what decides
//! conformer energies: an intramolecular hydrogen bond, a stacked pair of rings, a dipole that
//! prefers one torsion over another.
//!
//! Each conformer is optimised in its own scratch directory, so they can run at the same time.
//! That is what [`crate::qchem_interfaces::xtbrun::call_xtb_in`] exists for.

use std::path::{Path, PathBuf};

use anyhow::{anyhow, Result};

use crate::optimizer::minimum_cartesian::{
    minimize_cartesian, CartesianMinimumMethod, CartesianMinimumOptions,
};
use crate::optimizer::traits::Objective;
use crate::qchem_interfaces::xtbrun::{call_xtb_in, QcInput};

use super::{Conformer, ConformerOutcome};

/// How many of the search's conformers to re-optimise.
///
/// A cap, because the point is to sharpen the ordering of the shapes worth caring about, not to
/// re-run the whole search at a hundred times the cost. Anything past this is far enough up in
/// energy that its exact place in the list does not matter.
pub const MAX_REFINED: usize = 24;

/// An xTB calculation presented as an optimisable objective.
///
/// xTB already speaks bohr and Hartree, so unlike DREIDING there is nothing to convert.
struct XtbObjective<'a> {
    workdir: PathBuf,
    atoms: &'a [String],
    input: QcInput,
    /// Counted so a caller can report how much external work a refinement cost.
    calls: usize,
}

impl Objective for XtbObjective<'_> {
    fn energy(&mut self, x: &[f64]) -> Result<f64> {
        self.calls += 1;
        // `--dipole` is the cheapest single-point argument that still prints the total energy.
        call_xtb_in(&self.workdir, x, self.atoms, &self.input, "--dipole")
            .map(|(energy, _)| energy)
            .map_err(|e| anyhow!("{e}"))
    }

    fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
        self.calls += 1;
        let (_, gradient) = call_xtb_in(&self.workdir, x, self.atoms, &self.input, "--grad")
            .map_err(|e| anyhow!("{e}"))?;
        if gradient.len() != grad.len() {
            return Err(anyhow!(
                "xTB returned {} gradient components for {} coordinates",
                gradient.len(),
                grad.len()
            ));
        }
        grad.copy_from_slice(&gradient);
        Ok(())
    }

    fn energy_gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<f64> {
        self.calls += 1;
        // One call for both, which halves the number of process launches.
        let (energy, gradient) = call_xtb_in(&self.workdir, x, self.atoms, &self.input, "--grad")
            .map_err(|e| anyhow!("{e}"))?;
        if gradient.len() != grad.len() {
            return Err(anyhow!(
                "xTB returned {} gradient components for {} coordinates",
                gradient.len(),
                grad.len()
            ));
        }
        grad.copy_from_slice(&gradient);
        Ok(energy)
    }
}

/// Convergence for a refinement.
///
/// Looser than a publication geometry and deliberately so: the conformers arrive already relaxed
/// on DREIDING, so this only has to move them to the nearest GFN2 minimum and report its energy.
/// Every extra cycle is another process launch.
fn refine_options() -> CartesianMinimumOptions {
    CartesianMinimumOptions {
        method: CartesianMinimumMethod::Bfgs,
        max_cycles: 60,
        max_gradient: 5.0e-4,
        rms_gradient: 2.5e-4,
        verbosity: 0,
        ..CartesianMinimumOptions::default()
    }
}

/// What a refinement did.
#[derive(Debug, Clone)]
pub struct Refinement {
    /// How many conformers were re-optimised.
    pub refined: usize,
    /// How many xTB calls that took.
    pub xtb_calls: usize,
    /// How many conformers changed places in the ranking.
    pub reordered: usize,
}

/// Re-optimises the leading conformers with GFN2-xTB and re-ranks the set.
///
/// `binary` is the xTB executable. Conformers past [`MAX_REFINED`] are left where they are, after
/// the refined ones, and keep their DREIDING energies -- which is why the outcome records that the
/// energies are mixed rather than pretending otherwise.
pub fn refine_with_xtb(
    outcome: &mut ConformerOutcome,
    binary: &Path,
    scratch: &Path,
    charge: i32,
    multiplicity: i32,
) -> Result<Refinement> {
    if outcome.conformers.is_empty() {
        return Ok(Refinement {
            refined: 0,
            xtb_calls: 0,
            reordered: 0,
        });
    }

    let input = QcInput {
        path: binary.to_path_buf(),
        functional: String::new(),
        charge,
        multiplicity,
        // GFN2-xTB, one thread each.
        //
        // GFN2 rather than GFN-FF: it is a semiempirical tight-binding method with an SCF, not a
        // force field, so it describes polarisation and dispersion properly -- which is the whole
        // reason to re-rank at all. It costs more per structure than GFN-FF, but the count is
        // dozens of conformers rather than the tens of thousands the search itself needs, so the
        // cost lands in minutes either way.
        //
        // The thread limit is so that running several conformers at once does not oversubscribe
        // the machine.
        additional: "--gfn 2 -P 1".to_string(),
        wfu: false,
    };

    let take = outcome.conformers.len().min(MAX_REFINED);
    let before: Vec<Vec<f64>> = outcome.conformers[..take]
        .iter()
        .map(|c| c.positions_angstrom.clone())
        .collect();

    let mut energies = Vec::with_capacity(take);
    let mut calls = 0usize;
    for (index, positions) in before.iter().enumerate() {
        let mut objective = XtbObjective {
            workdir: scratch.join(format!("conformer_{index:03}")),
            atoms: &outcome.atoms,
            input: input.clone(),
            calls: 0,
        };
        let x_bohr: Vec<f64> = positions
            .iter()
            .map(|a| a / crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM)
            .collect();
        let result = minimize_cartesian(x_bohr, &mut objective, refine_options())?;
        calls += objective.calls;
        energies.push((index, result.energy, result.x));
    }

    // Re-rank on the GFN2 energies, lowest first.
    let mut order: Vec<usize> = (0..take).collect();
    order.sort_by(|a, b| energies[*a].1.total_cmp(&energies[*b].1));
    let reordered = order.iter().enumerate().filter(|(at, was)| *at != **was).count();

    let lowest = energies
        .iter()
        .map(|(_, energy, _)| *energy)
        .fold(f64::INFINITY, f64::min);
    let refined: Vec<Conformer> = order
        .iter()
        .map(|&slot| {
            let (_, energy, ref x) = energies[slot];
            Conformer {
                positions_angstrom: x
                    .iter()
                    .map(|b| b * crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM)
                    .collect(),
                relative_energy_kcal: (energy - lowest)
                    * crate::forcefield::dreiding::objective::HARTREE_TO_KCAL,
                // Filled in below, once the whole set is known.
                population_fraction: 0.0,
            }
        })
        .collect();

    // Everything past the cap keeps its place and its DREIDING energy, shifted so it still sits
    // above the refined set rather than appearing to compete with it.
    let tail_offset = refined
        .last()
        .map(|c| c.relative_energy_kcal)
        .unwrap_or(0.0);
    let tail: Vec<Conformer> = outcome.conformers[take..]
        .iter()
        .map(|c| Conformer {
            positions_angstrom: c.positions_angstrom.clone(),
            relative_energy_kcal: c.relative_energy_kcal + tail_offset,
            population_fraction: 0.0,
        })
        .collect();

    outcome.conformers = refined.into_iter().chain(tail).collect();
    let relative: Vec<f64> = outcome
        .conformers
        .iter()
        .map(|c| c.relative_energy_kcal)
        .collect();
    for (conformer, fraction) in outcome
        .conformers
        .iter_mut()
        .zip(super::populations(&relative))
    {
        conformer.population_fraction = fraction;
    }
    outcome.refined_with_xtb = true;

    Ok(Refinement {
        refined: take,
        xtb_calls: calls,
        reordered,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Runs the whole refinement against the real xTB the user has configured. Ignored, since it
    /// needs an external program and takes seconds:
    ///
    /// ```text
    /// cargo test --release --bins refinement_against_real_xtb -- --ignored --nocapture
    /// ```
    #[test]
    #[ignore = "needs the configured xTB executable"]
    fn refinement_against_real_xtb() {
        use crate::conformer::{search, ConformerSettings, Thoroughness};
        use crate::molecule::Molecule;
        use crate::qchem_interfaces::program::QcProgram;

        let configured = crate::qchem_interfaces::config::load_path(QcProgram::Xtb);
        let Ok(binary) =
            crate::qchem_interfaces::xtb_optimize::resolve_program_executable(QcProgram::Xtb, &configured)
        else {
            println!("  no xTB configured; nothing to test against");
            return;
        };
        println!("  xTB: {}", binary.display());

        // n-Butane: two conformers, anti and gauche, about 0.9 kcal/mol apart experimentally.
        let mut mol = Molecule::from_xyz(
            "14\n\n\
             C -1.9200  0.2600  0.0000\n\
             C -0.6400 -0.5700  0.0000\n\
             C  0.6400  0.2600  0.0000\n\
             C  1.9200 -0.5700  0.0000\n\
             H -2.8200 -0.3600  0.0000\n\
             H -1.9600  0.9000  0.8900\n\
             H -1.9600  0.9000 -0.8900\n\
             H -0.6200 -1.2200  0.8900\n\
             H -0.6200 -1.2200 -0.8900\n\
             H  0.6200  0.9100  0.8900\n\
             H  0.6200  0.9100 -0.8900\n\
             H  2.8200  0.0500  0.0000\n\
             H  1.9600 -1.2100  0.8900\n\
             H  1.9600 -1.2100 -0.8900\n",
        );
        mol.recompute_bonds(1.2, 2.5);

        let mut outcome = search(
            &mol,
            ConformerSettings {
                thoroughness: Thoroughness::Quick,
                ncore: 2,
                ..ConformerSettings::default()
            },
        )
        .expect("butane is searchable");
        println!("  DREIDING:");
        for (index, c) in outcome.conformers.iter().enumerate() {
            println!("    {}: {:>6.2} kcal/mol", index + 1, c.relative_energy_kcal);
        }

        let scratch = std::env::temp_dir().join(format!("beavyr_refine_{}", std::process::id()));
        let report = refine_with_xtb(&mut outcome, &binary, &scratch, 0, 1)
            .expect("the refinement runs");
        let _ = std::fs::remove_dir_all(&scratch);

        println!(
            "  GFN2 ({} conformers, {} xTB calls, {} moved):",
            report.refined, report.xtb_calls, report.reordered
        );
        for (index, c) in outcome.conformers.iter().enumerate() {
            println!("    {}: {:>6.2} kcal/mol", index + 1, c.relative_energy_kcal);
        }

        assert!(outcome.refined_with_xtb);
        assert!(report.xtb_calls > 0, "no xTB call was made");
        // The lowest is still the reference, and populations still sum to one.
        assert_eq!(outcome.conformers[0].relative_energy_kcal, 0.0);
        let total: f64 = outcome.conformers.iter().map(|c| c.population_fraction).sum();
        assert!((total - 1.0).abs() < 1.0e-9, "populations sum to {total}");
        // And the list is ordered, or the table would read wrongly.
        for pair in outcome.conformers.windows(2) {
            assert!(pair[0].relative_energy_kcal <= pair[1].relative_energy_kcal + 1.0e-9);
        }
    }

    #[test]
    fn the_refinement_cap_is_a_sensible_size() {
        // Enough to cover the thermally accessible set, small enough that the external cost stays
        // in minutes rather than hours.
        assert!(MAX_REFINED >= 10 && MAX_REFINED <= 50);
    }

    /// An empty set must be a no-op rather than an error: a search that found nothing is not a
    /// refinement failure.
    #[test]
    fn refining_nothing_is_not_an_error() {
        let mut outcome = ConformerOutcome {
            conformers: Vec::new(),
            rotors: 0,
            rings: 0,
            generations: 0,
            structures_tried: 0,
            local_optimizations: 0,
            termination: String::new(),
            elapsed: std::time::Duration::ZERO,
            workers: 1,
            local_optimizer: crate::optimizer::conformer_search::LocalOptimizer::Cartesian,
            best_found_in_generation: 0,
            atoms: Vec::new(),
            radical_warning: None,
            refined_with_xtb: false,
        };
        let report = refine_with_xtb(
            &mut outcome,
            Path::new("/nonexistent/xtb"),
            Path::new("/tmp"),
            0,
            1,
        )
        .expect("an empty set is a no-op");
        assert_eq!(report.refined, 0);
        assert!(!outcome.refined_with_xtb);
    }
}
