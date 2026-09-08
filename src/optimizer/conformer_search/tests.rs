use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Duration;

use anyhow::Result;

use crate::optimizer::internal_coords::ConnectivityModel;
use crate::optimizer::Objective;

use super::geometry::{coordinates_from_genome, dihedral_angle, wrap_angle};
use super::rmsd::test_aligned_rmsd;
use super::topology::{PreparedTorsion, TorsionKind};
use super::{genetic_conformer_search, genetic_conformer_search_parallel, ConformerSearchOptions};

fn chain_geometry() -> Vec<f64> {
    vec![
        0.0, 0.0, 0.0, // C1
        2.5, 0.0, 0.0, // C2
        3.5, 2.0, 0.0, // C3
        4.5, 2.0, 2.0, // C4
    ]
}

fn chain_connectivity() -> ConnectivityModel {
    ConnectivityModel {
        rcov_disp: vec![1.3; 4],
        kcn: 7.5,
        facmin: 0.45,
    }
}

#[test]
fn torsion_builder_sets_the_requested_angle() {
    let template = chain_geometry();
    let torsion = PreparedTorsion {
        atoms: [0, 1, 2, 3],
        central_bond: [1, 2],
        moving_atoms: vec![2, 3],
        kind: TorsionKind::Rotatable,
        user_selected: true,
    };
    let target = -1.2;
    let coordinates = coordinates_from_genome(&template, &[torsion.clone()], &[target]).unwrap();
    assert!(wrap_angle(dihedral_angle(&coordinates, torsion.atoms) - target).abs() < 1.0e-8);
}

#[test]
fn heavy_atom_rmsd_removes_translation_and_rotation() {
    let reference = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0];
    let transformed = vec![3.0, -2.0, 1.0, 3.0, -1.0, 1.0, 1.0, -2.0, 1.0];
    assert!(test_aligned_rmsd(&reference, &transformed, &[0, 1, 2]) < 1.0e-10);
}

struct FlatObjective;

impl Objective for FlatObjective {
    fn energy(&mut self, _coordinates: &[f64]) -> Result<f64> {
        Ok(-1.0)
    }

    fn gradient(&mut self, _coordinates: &[f64], gradient: &mut [f64]) -> Result<()> {
        gradient.fill(0.0);
        Ok(())
    }

    fn energy_gradient(&mut self, _coordinates: &[f64], gradient: &mut [f64]) -> Result<f64> {
        gradient.fill(0.0);
        Ok(-1.0)
    }
}

#[test]
fn genetic_search_builds_a_reproducible_population() {
    let mut options = ConformerSearchOptions::default();
    options.population_size = 2;
    options.max_generations = 1;
    options.stagnation_generations = 2;
    options.rotatable_dihedrals = vec![[0, 1, 2, 3]];
    options.minimum_nonbonded_distance_angstrom = 0.0;
    options.maximum_bonded_distance_angstrom = 10.0;
    options.rmsd_cutoff_angstrom = 0.001;
    options.torsion_mutation_probability = 1.0;
    options.random_seed = Some(7);
    options.verbosity = 0;
    let elements = vec!["C".to_string(); 4];
    let result = genetic_conformer_search(
        chain_geometry(),
        &elements,
        chain_connectivity(),
        &mut FlatObjective,
        options,
    )
    .unwrap();
    assert_eq!(result.final_population.len(), 2);
    assert!(result.conformers.len() >= 2);
    assert!(result
        .conformers
        .iter()
        .all(|conformer| conformer.energy_hartree == -1.0));
}

struct CountingObjective {
    active: Arc<AtomicUsize>,
    maximum: Arc<AtomicUsize>,
}

impl CountingObjective {
    fn evaluate(&self, gradient: Option<&mut [f64]>) -> f64 {
        let active = self.active.fetch_add(1, Ordering::SeqCst) + 1;
        self.maximum.fetch_max(active, Ordering::SeqCst);
        // Keep the deliberately trivial evaluation alive long enough for the
        // test to observe overlapping initial-population jobs.
        std::thread::sleep(Duration::from_millis(20));
        if let Some(gradient) = gradient {
            gradient.fill(0.0);
        }
        self.active.fetch_sub(1, Ordering::SeqCst);
        -1.0
    }
}

impl Objective for CountingObjective {
    fn energy(&mut self, _coordinates: &[f64]) -> Result<f64> {
        Ok(self.evaluate(None))
    }

    fn gradient(&mut self, _coordinates: &[f64], gradient: &mut [f64]) -> Result<()> {
        self.evaluate(Some(gradient));
        Ok(())
    }

    fn energy_gradient(&mut self, _coordinates: &[f64], gradient: &mut [f64]) -> Result<f64> {
        Ok(self.evaluate(Some(gradient)))
    }
}

#[test]
fn parallel_search_overlaps_initial_local_optimizations() {
    let mut options = ConformerSearchOptions::default();
    options.parallel_workers = Some(4);
    options.population_size = 4;
    options.max_generations = 1;
    options.stagnation_generations = 2;
    options.rotatable_dihedrals = vec![[0, 1, 2, 3]];
    options.minimum_nonbonded_distance_angstrom = 0.0;
    options.maximum_bonded_distance_angstrom = 10.0;
    options.rmsd_cutoff_angstrom = 0.001;
    options.torsion_mutation_probability = 1.0;
    options.random_seed = Some(17);
    options.verbosity = 0;

    let active = Arc::new(AtomicUsize::new(0));
    let maximum = Arc::new(AtomicUsize::new(0));
    let factory_active = Arc::clone(&active);
    let factory_maximum = Arc::clone(&maximum);
    let result = genetic_conformer_search_parallel(
        chain_geometry(),
        &vec!["C".to_string(); 4],
        chain_connectivity(),
        || {
            Ok(CountingObjective {
                active: Arc::clone(&factory_active),
                maximum: Arc::clone(&factory_maximum),
            })
        },
        options,
    )
    .unwrap();

    assert_eq!(result.final_population.len(), 4);
    let expected = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).min(4);
    assert_eq!(maximum.load(Ordering::SeqCst), expected);
}
