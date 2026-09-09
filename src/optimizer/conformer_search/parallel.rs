//! Parallel local-optimization driver for the genetic conformer search.
//!
//! Trial generation, blacklist updates, and population selection remain
//! serial and deterministic. Only the expensive, independent local geometry
//! optimizations are evaluated concurrently, each with a fresh objective.

use anyhow::{bail, ensure, Result};
use crate::optimizer::prng::StdRng;

use crate::optimizer::internal_coords::ConnectivityModel;
use super::local::local_optimize;
use crate::optimizer::{GeomOptResult, Objective};

use super::genetic::{
    crossover, energy_order, mutated_candidate, random_genome, select_parents, validate_options,
    Conformer, ConformerSearchResult, ConformerSearchStatistics, TerminationReason,
};
use super::geometry::{coordinates_from_genome, genome_from_coordinates, sensible_geometry};
use super::options::ConformerSearchOptions;
use super::output::{print_generation, print_search_header};
use super::rmsd::{heavy_atom_indices, unique_against_blacklist};
use super::topology::{discover_torsions, PreparedTorsion};

struct LocalOutcome {
    result: GeomOptResult,
    generation: usize,
}

fn optimize_batch<O, F>(
    starts: Vec<(Vec<f64>, usize)>,
    connectivity: &ConnectivityModel,
    options: &ConformerSearchOptions,
    objective_factory: &F,
) -> Result<Vec<LocalOutcome>>
where
    O: Objective + Send,
    F: Fn() -> Result<O> + Sync,
{
    // Behemoth relies on its callers to size each batch. Beavyr hands it a whole generation at
    // once, so the cap lives here instead: the batch is worked through `workers` at a time, in
    // dispatch order, which is what keeps the result independent of how many cores were used.
    let workers = options.effective_parallel_workers();
    let mut outcomes = Vec::with_capacity(starts.len());
    let mut remaining = starts;
    while !remaining.is_empty() {
        let take = workers.min(remaining.len());
        let chunk: Vec<(Vec<f64>, usize)> = remaining.drain(..take).collect();
        let mut finished = optimize_chunk::<O, F>(chunk, connectivity, options, objective_factory)?;
        outcomes.append(&mut finished);
    }
    Ok(outcomes)
}

/// One batch of at most `workers` structures, each on its own thread.
fn optimize_chunk<O, F>(
    starts: Vec<(Vec<f64>, usize)>,
    connectivity: &ConnectivityModel,
    options: &ConformerSearchOptions,
    objective_factory: &F,
) -> Result<Vec<LocalOutcome>>
where
    O: Objective + Send,
    F: Fn() -> Result<O> + Sync,
{
    std::thread::scope(|scope| {
        let handles = starts
            .into_iter()
            .map(|(coordinates, generation)| {
                let connectivity = connectivity.clone();
                let local_options = options.local_optimization.clone();
                let local_method = options.local_optimizer;
                let local_constraints = options.constraints.clone();
                // Behemoth confines each worker to a one-thread Rayon pool here, so that Rayon
                // inside its electronic-structure method cannot borrow the other conformer
                // workers' threads. Beavyr's force field uses neither Rayon nor BLAS, so one
                // worker is already one thread and the pool is unnecessary.
                scope.spawn(move || {
                    let mut objective = objective_factory()?;
                    let result = local_optimize(
                        local_method,
                        coordinates,
                        connectivity,
                        &mut objective,
                        local_options,
                        &local_constraints,
                    )?;
                    Ok::<LocalOutcome, anyhow::Error>(LocalOutcome { result, generation })
                })
            })
            .collect::<Vec<_>>();

        let mut outcomes = Vec::with_capacity(handles.len());
        for handle in handles {
            let outcome = handle
                .join()
                .map_err(|_| anyhow::anyhow!("parallel conformer worker panicked"))??;
            outcomes.push(outcome);
        }
        Ok(outcomes)
    })
}

fn accept_outcomes(
    outcomes: Vec<LocalOutcome>,
    torsions: &[PreparedTorsion],
    optimized_blacklist: &mut Vec<Vec<f64>>,
    blacklist: &mut Vec<Vec<f64>>,
    rmsd_atoms: &[usize],
    options: &ConformerSearchOptions,
    statistics: &mut ConformerSearchStatistics,
) -> Vec<Conformer> {
    let mut accepted = Vec::new();
    for outcome in outcomes {
        let result = outcome.result;
        blacklist.push(result.x.clone());
        if !result.converged {
            statistics.nonconverged_optimizations += 1;
            continue;
        }
        if !unique_against_blacklist(
            &result.x,
            optimized_blacklist,
            rmsd_atoms,
            options.rmsd_cutoff_angstrom,
            options.chiral,
        ) {
            statistics.duplicate_minima += 1;
            continue;
        }
        optimized_blacklist.push(result.x.clone());
        accepted.push(Conformer {
            torsion_angles_radians: genome_from_coordinates(&result.x, torsions),
            coordinates_bohr: result.x,
            energy_hartree: result.energy,
            local_optimization_cycles: result.cycles,
            generation: outcome.generation,
        });
    }
    accepted
}

/// Run the genetic conformer search with independent local optimizations in
/// parallel. The factory must return a fresh electronic-structure objective,
/// including an independent molecule and SCF state, for every call.
pub fn genetic_conformer_search_parallel<O, F>(
    template_coordinates_bohr: Vec<f64>,
    elements: &[String],
    connectivity: ConnectivityModel,
    objective_factory: F,
    options: ConformerSearchOptions,
) -> Result<ConformerSearchResult>
where
    O: Objective + Send,
    F: Fn() -> Result<O> + Sync,
{
    validate_options(&options)?;
    ensure!(
        template_coordinates_bohr.len() == 3 * elements.len(),
        "conformer-search coordinate/element count mismatch"
    );
    let workers = options.effective_parallel_workers();
    // Each job receives one compute thread; keep process-global BLAS at one
    // thread as well to prevent nested oversubscription.

    let (bonds, torsions) = discover_torsions(
        &template_coordinates_bohr,
        elements,
        &connectivity,
        &options.rotatable_dihedrals,
        &options.cis_trans_dihedrals,
        options.exclude_methyl_rotors,
        &options.constraints,
    )?;
    let rmsd_atoms = heavy_atom_indices(elements);
    let mut rng = match options.random_seed {
        Some(seed) => StdRng::seed_from_u64(seed),
        None => StdRng::from_entropy(),
    };
    if options.verbosity > 0 {
        print_search_header(&options, &torsions, workers);
    }

    let mut statistics = ConformerSearchStatistics::default();
    let mut blacklist = Vec::new();
    let mut optimized_blacklist = Vec::new();
    let mut population = Vec::new();
    let mut all_conformers = Vec::new();
    let mut initial_attempts = 0usize;

    while population.len() < options.population_size
        && initial_attempts < options.max_initial_attempts
    {
        let batch_capacity = workers.min(options.population_size - population.len());
        let mut starts = Vec::with_capacity(batch_capacity);
        while starts.len() < batch_capacity && initial_attempts < options.max_initial_attempts {
            initial_attempts += 1;
            statistics.attempted_structures += 1;
            let genome = random_genome(&mut rng, &torsions);
            let Ok(coordinates) =
                coordinates_from_genome(&template_coordinates_bohr, &torsions, &genome)
            else {
                statistics.rejected_insensible += 1;
                continue;
            };
            if !sensible_geometry(
                &coordinates,
                elements.len(),
                &bonds,
                options.minimum_nonbonded_distance_angstrom,
                options.maximum_bonded_distance_angstrom,
            ) {
                statistics.rejected_insensible += 1;
                continue;
            }
            if options.unique_initial_population
                && !unique_against_blacklist(
                    &coordinates,
                    &blacklist,
                    &rmsd_atoms,
                    options.rmsd_cutoff_angstrom,
                    options.chiral,
                )
            {
                statistics.rejected_blacklisted += 1;
                continue;
            }
            // Reserve the start before generating the rest of the batch, so
            // two workers never optimize equivalent starting structures.
            blacklist.push(coordinates.clone());
            starts.push((coordinates, 0));
        }
        if starts.is_empty() {
            break;
        }
        statistics.local_optimizations += starts.len();
        let outcomes = optimize_batch::<O, F>(starts, &connectivity, &options, &objective_factory)?;
        for conformer in accept_outcomes(
            outcomes,
            &torsions,
            &mut optimized_blacklist,
            &mut blacklist,
            &rmsd_atoms,
            &options,
            &mut statistics,
        ) {
            population.push(conformer.clone());
            all_conformers.push(conformer);
        }
    }
    if population.len() < options.population_size {
        bail!(
            "conformer-search initialization produced only {} distinct converged minima for a population of {}; increase `max_initial_attempts`, relax screening/RMSD cutoffs, or reduce `population_size`",
            population.len(), options.population_size
        );
    }
    population.sort_by(energy_order);
    if options.verbosity > 0 {
        print_generation(0, &population, all_conformers.len(), &statistics);
    }

    let mut generations = 0;
    let mut stagnant = 0;
    let mut previous_best = population[0].energy_hartree;
    let mut termination = TerminationReason::MaximumGenerations;

    for generation in 1..=options.max_generations {
        generations = generation;
        if let Some(target) = options.target_energy_hartree {
            if population[0].energy_hartree <= target {
                termination = TerminationReason::TargetEnergyReached {
                    target_hartree: target,
                };
                break;
            }
        }
        // Behemoth breeds one crossover pair per generation, which leaves at most two structures
        // to optimise however many cores are available. `offspring_per_generation` widens that so
        // a generation can actually fill a machine. It is a property of the search, deliberately
        // not of the core count: children are bred from the RNG in sequence before any is
        // dispatched, so the same seed breeds the same generation whatever `workers` is, and only
        // the concurrency of the optimisations changes.
        let wanted = options.offspring_per_generation.max(1);
        let mut starts = Vec::with_capacity(wanted);
        while starts.len() < wanted {
            let (first, second) = select_parents(&mut rng, &population, &options);
            let (base_a, base_b) = crossover(
                &mut rng,
                &template_coordinates_bohr,
                &torsions,
                &population[first].torsion_angles_radians,
                &population[second].torsion_angles_radians,
                &bonds,
                &options,
            );
            let mut bred_any = false;
            for base in [&base_a, &base_b] {
                if starts.len() >= wanted {
                    break;
                }
                let Some((_genome, coordinates)) = mutated_candidate(
                    &mut rng,
                    &template_coordinates_bohr,
                    &torsions,
                    base,
                    &bonds,
                    &blacklist,
                    &rmsd_atoms,
                    &options,
                    &mut statistics,
                ) else {
                    continue;
                };
                blacklist.push(coordinates.clone());
                starts.push((coordinates, generation));
                bred_any = true;
            }
            // A crossover that yields nothing usable means the neighbourhood is exhausted; take
            // whatever the generation managed rather than spinning.
            if !bred_any {
                break;
            }
        }
        if starts.is_empty() {
            termination = TerminationReason::MutationAttemptsExhausted;
            break;
        }
        statistics.local_optimizations += starts.len();
        let outcomes = optimize_batch::<O, F>(starts, &connectivity, &options, &objective_factory)?;
        for conformer in accept_outcomes(
            outcomes,
            &torsions,
            &mut optimized_blacklist,
            &mut blacklist,
            &rmsd_atoms,
            &options,
            &mut statistics,
        ) {
            population.push(conformer.clone());
            all_conformers.push(conformer);
        }

        population.sort_by(energy_order);
        population.truncate(options.population_size);
        let best = population[0].energy_hartree;
        if (best - previous_best).abs() <= options.energy_convergence_hartree {
            stagnant += 1;
        } else {
            stagnant = 0;
        }
        previous_best = best;
        if options.verbosity > 0 {
            print_generation(generation, &population, all_conformers.len(), &statistics);
        }
        if stagnant >= options.stagnation_generations {
            termination = TerminationReason::EnergyStagnation {
                generations: stagnant,
            };
            break;
        }
    }

    all_conformers.sort_by(energy_order);
    population.sort_by(energy_order);
    Ok(ConformerSearchResult {
        conformers: all_conformers,
        final_population: population,
        torsions,
        generations,
        termination,
        statistics,
        output_energy_window_hartree: options.output_energy_window_hartree,
    })
}
