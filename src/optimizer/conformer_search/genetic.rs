use std::cmp::Ordering;
use std::fmt;

use anyhow::{bail, ensure, Result};
use crate::optimizer::prng::StdRng;


use crate::optimizer::internal_coords::ConnectivityModel;
use super::local::local_optimize;
use crate::optimizer::Objective;

use super::geometry::{
    coordinates_from_genome, genome_from_coordinates, sensible_geometry, wrap_angle,
};
use super::options::{ConformerSearchOptions, SelectionMode};
use super::output::{print_generation, print_search_header};
use super::rmsd::{heavy_atom_indices, unique_against_blacklist};
use super::topology::{discover_torsions, PreparedTorsion, TorsionKind};

#[derive(Debug, Clone)]
pub struct Conformer {
    pub coordinates_bohr: Vec<f64>,
    pub energy_hartree: f64,
    pub torsion_angles_radians: Vec<f64>,
    pub local_optimization_cycles: usize,
    pub generation: usize,
}

#[derive(Debug, Clone, Default)]
pub struct ConformerSearchStatistics {
    pub attempted_structures: usize,
    pub local_optimizations: usize,
    pub nonconverged_optimizations: usize,
    pub rejected_insensible: usize,
    pub rejected_blacklisted: usize,
    pub duplicate_minima: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum TerminationReason {
    MaximumGenerations,
    EnergyStagnation { generations: usize },
    TargetEnergyReached { target_hartree: f64 },
    MutationAttemptsExhausted,
}

impl fmt::Display for TerminationReason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MaximumGenerations => formatter.write_str("maximum generation count reached"),
            Self::EnergyStagnation { generations } => write!(
                formatter,
                "lowest energy unchanged within tolerance for {generations} generations"
            ),
            Self::TargetEnergyReached { target_hartree } => {
                write!(formatter, "target energy {target_hartree:.12} Eh reached")
            }
            Self::MutationAttemptsExhausted => {
                formatter.write_str("no new sensible, unique mutation could be generated")
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct ConformerSearchResult {
    pub conformers: Vec<Conformer>,
    pub final_population: Vec<Conformer>,
    pub torsions: Vec<PreparedTorsion>,
    pub generations: usize,
    pub termination: TerminationReason,
    pub statistics: ConformerSearchStatistics,
    pub output_energy_window_hartree: f64,
}

pub(super) fn energy_order(left: &Conformer, right: &Conformer) -> Ordering {
    left.energy_hartree
        .partial_cmp(&right.energy_hartree)
        .unwrap_or(Ordering::Equal)
}

pub(super) fn random_genome(rng: &mut StdRng, torsions: &[PreparedTorsion]) -> Vec<f64> {
    torsions
        .iter()
        .map(|torsion| match torsion.kind {
            TorsionKind::Rotatable => rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI),
            TorsionKind::CisTrans => {
                if rng.gen_bool(0.5) {
                    0.0
                } else {
                    std::f64::consts::PI
                }
            }
        })
        .collect()
}

fn roulette_index(rng: &mut StdRng, weights: &[f64]) -> usize {
    let sum = weights.iter().copied().sum::<f64>();
    if sum <= 0.0 || !sum.is_finite() {
        return rng.gen_range(0..weights.len());
    }
    let mut draw = rng.gen_range(0.0..sum);
    for (index, &weight) in weights.iter().enumerate() {
        draw -= weight.max(0.0);
        if draw <= 0.0 {
            return index;
        }
    }
    weights.len() - 1
}

pub(super) fn select_parents(
    rng: &mut StdRng,
    population: &[Conformer],
    options: &ConformerSearchOptions,
) -> (usize, usize) {
    if options.selection == SelectionMode::Random {
        let first = rng.gen_range(0..population.len());
        let mut second = rng.gen_range(0..population.len() - 1);
        if second >= first {
            second += 1;
        }
        return (first, second);
    }

    let minimum = population.first().unwrap().energy_hartree;
    let maximum = population.last().unwrap().energy_hartree;
    let variance = maximum - minimum;
    let mut fitness = if variance < options.energy_variance_hartree {
        vec![1.0; population.len()]
    } else {
        population
            .iter()
            .map(|individual| (maximum - individual.energy_hartree) / variance)
            .collect::<Vec<_>>()
    };
    if options.selection == SelectionMode::ReverseRoulette {
        fitness.reverse();
    }
    if fitness.iter().sum::<f64>() < options.fitness_sum_limit {
        let first = 0;
        let second = rng.gen_range(1..population.len());
        return (first, second);
    }
    let first = roulette_index(rng, &fitness);
    let mut second = roulette_index(rng, &fitness);
    for _ in 0..32 {
        if second != first {
            return (first, second);
        }
        second = roulette_index(rng, &fitness);
    }
    second = (first + 1 + rng.gen_range(0..population.len() - 1)) % population.len();
    (first, second)
}

pub(super) fn crossover(
    rng: &mut StdRng,
    template: &[f64],
    torsions: &[PreparedTorsion],
    first: &[f64],
    second: &[f64],
    bonds: &[[usize; 2]],
    options: &ConformerSearchOptions,
) -> (Vec<f64>, Vec<f64>) {
    if torsions.len() < 2 || !rng.gen_bool(options.crossover_probability) {
        return (first.to_vec(), second.to_vec());
    }
    for _ in 0..options.crossover_attempts {
        let cut = rng.gen_range(1..torsions.len());
        let mut child_a = first[..cut].to_vec();
        child_a.extend_from_slice(&second[cut..]);
        let mut child_b = second[..cut].to_vec();
        child_b.extend_from_slice(&first[cut..]);
        let sensible = |genome: &[f64]| {
            coordinates_from_genome(template, torsions, genome)
                .ok()
                .is_some_and(|coordinates| {
                    sensible_geometry(
                        &coordinates,
                        template.len() / 3,
                        bonds,
                        options.minimum_nonbonded_distance_angstrom,
                        options.maximum_bonded_distance_angstrom,
                    )
                })
        };
        if sensible(&child_a) && sensible(&child_b) {
            return (child_a, child_b);
        }
    }
    (first.to_vec(), second.to_vec())
}

fn mutate_positions(
    rng: &mut StdRng,
    genome: &mut [f64],
    torsions: &[PreparedTorsion],
    kind: TorsionKind,
    probability: f64,
    maximum: usize,
) -> bool {
    let positions = torsions
        .iter()
        .enumerate()
        .filter_map(|(index, torsion)| (torsion.kind == kind).then_some(index))
        .collect::<Vec<_>>();
    if positions.is_empty() || maximum == 0 || !rng.gen_bool(probability) {
        return false;
    }
    let count = rng.gen_range(1..=maximum.min(positions.len()));
    let mut available = positions;
    for _ in 0..count {
        let selected = rng.gen_range(0..available.len());
        let index = available.swap_remove(selected);
        genome[index] = match kind {
            TorsionKind::Rotatable => {
                let degrees: i32 = rng.gen_range(-179..=180);
                (degrees as f64).to_radians()
            }
            TorsionKind::CisTrans => {
                if wrap_angle(genome[index]).abs() <= std::f64::consts::FRAC_PI_2 {
                    std::f64::consts::PI
                } else {
                    0.0
                }
            }
        };
    }
    true
}

pub(super) fn mutated_candidate(
    rng: &mut StdRng,
    template: &[f64],
    torsions: &[PreparedTorsion],
    base: &[f64],
    bonds: &[[usize; 2]],
    blacklist: &[Vec<f64>],
    rmsd_atoms: &[usize],
    options: &ConformerSearchOptions,
    statistics: &mut ConformerSearchStatistics,
) -> Option<(Vec<f64>, Vec<f64>)> {
    for _ in 0..options.mutation_attempts {
        statistics.attempted_structures += 1;
        let mut genome = base.to_vec();
        let rot_changed = mutate_positions(
            rng,
            &mut genome,
            torsions,
            TorsionKind::Rotatable,
            options.torsion_mutation_probability,
            options.max_torsion_mutations,
        );
        let ct_changed = mutate_positions(
            rng,
            &mut genome,
            torsions,
            TorsionKind::CisTrans,
            options.cis_trans_mutation_probability,
            options.max_cis_trans_mutations,
        );
        let _mutation_applied = rot_changed || ct_changed;
        let Ok(coordinates) = coordinates_from_genome(template, torsions, &genome) else {
            statistics.rejected_insensible += 1;
            continue;
        };
        if !sensible_geometry(
            &coordinates,
            template.len() / 3,
            bonds,
            options.minimum_nonbonded_distance_angstrom,
            options.maximum_bonded_distance_angstrom,
        ) {
            statistics.rejected_insensible += 1;
            continue;
        }
        if !unique_against_blacklist(
            &coordinates,
            blacklist,
            rmsd_atoms,
            options.rmsd_cutoff_angstrom,
            options.chiral,
        ) {
            statistics.rejected_blacklisted += 1;
            continue;
        }
        return Some((genome, coordinates));
    }
    None
}

fn locally_optimize<O: Objective>(
    objective: &mut O,
    connectivity: &ConnectivityModel,
    torsions: &[PreparedTorsion],
    starting_coordinates: Vec<f64>,
    generation: usize,
    blacklist: &mut Vec<Vec<f64>>,
    optimized_blacklist: &mut Vec<Vec<f64>>,
    rmsd_atoms: &[usize],
    options: &ConformerSearchOptions,
    statistics: &mut ConformerSearchStatistics,
) -> Result<Option<Conformer>> {
    blacklist.push(starting_coordinates.clone());
    statistics.local_optimizations += 1;
    let result = local_optimize(
        options.local_optimizer,
        starting_coordinates,
        connectivity.clone(),
        objective,
        options.local_optimization.clone(),
    )?;
    let is_distinct_minimum = unique_against_blacklist(
        &result.x,
        optimized_blacklist,
        rmsd_atoms,
        options.rmsd_cutoff_angstrom,
        options.chiral,
    );
    blacklist.push(result.x.clone());
    if !result.converged {
        statistics.nonconverged_optimizations += 1;
        return Ok(None);
    }
    if !is_distinct_minimum {
        statistics.duplicate_minima += 1;
        return Ok(None);
    }
    optimized_blacklist.push(result.x.clone());
    Ok(Some(Conformer {
        torsion_angles_radians: genome_from_coordinates(&result.x, torsions),
        coordinates_bohr: result.x,
        energy_hartree: result.energy,
        local_optimization_cycles: result.cycles,
        generation,
    }))
}

pub(super) fn validate_options(options: &ConformerSearchOptions) -> Result<()> {
    ensure!(
        options.population_size >= 2,
        "conformer-search population size must be at least 2"
    );
    ensure!(
        options.parallel_workers.is_none_or(|workers| workers > 0),
        "conformer-search parallel worker count must be at least 1"
    );
    ensure!(
        options.max_generations > 0,
        "conformer-search maximum generations must be positive"
    );
    ensure!(
        options.stagnation_generations > 0,
        "conformer-search stagnation generations must be positive"
    );
    ensure!(
        options.mutation_attempts > 0,
        "conformer-search mutation attempts must be positive"
    );
    ensure!(
        options.max_initial_attempts >= options.population_size,
        "conformer-search initial-attempt limit must be at least the population size"
    );
    for (name, probability) in [
        ("crossover", options.crossover_probability),
        ("torsion mutation", options.torsion_mutation_probability),
        ("cis/trans mutation", options.cis_trans_mutation_probability),
    ] {
        ensure!(
            (0.0..=1.0).contains(&probability),
            "{name} probability must be between 0 and 1"
        );
    }
    ensure!(
        options.rmsd_cutoff_angstrom > 0.0,
        "conformer-search RMSD cutoff must be positive"
    );
    ensure!(
        options.minimum_nonbonded_distance_angstrom >= 0.0,
        "minimum nonbonded distance cannot be negative"
    );
    ensure!(
        options.maximum_bonded_distance_angstrom > 0.0,
        "maximum bonded distance must be positive"
    );
    Ok(())
}

/// Run the first-principles torsional genetic conformer search.
pub fn genetic_conformer_search<O: Objective>(
    template_coordinates_bohr: Vec<f64>,
    elements: &[String],
    connectivity: ConnectivityModel,
    objective: &mut O,
    options: ConformerSearchOptions,
) -> Result<ConformerSearchResult> {
    validate_options(&options)?;
    ensure!(
        template_coordinates_bohr.len() == 3 * elements.len(),
        "conformer-search coordinate/element count mismatch"
    );
    let (bonds, torsions) = discover_torsions(
        &template_coordinates_bohr,
        elements,
        &connectivity,
        &options.rotatable_dihedrals,
        &options.cis_trans_dihedrals,
        options.exclude_methyl_rotors,
    )?;
    let rmsd_atoms = heavy_atom_indices(elements);
    let mut rng = match options.random_seed {
        Some(seed) => StdRng::seed_from_u64(seed),
        None => StdRng::from_entropy(),
    };
    if options.verbosity > 0 {
        print_search_header(&options, &torsions, 1);
    }

    let mut statistics = ConformerSearchStatistics::default();
    let mut blacklist = Vec::new();
    let mut optimized_blacklist = Vec::new();
    let mut population = Vec::new();
    let mut all_conformers = Vec::new();

    for _ in 0..options.max_initial_attempts {
        if population.len() == options.population_size {
            break;
        }
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
        if let Some(conformer) = locally_optimize(
            objective,
            &connectivity,
            &torsions,
            coordinates,
            0,
            &mut blacklist,
            &mut optimized_blacklist,
            &rmsd_atoms,
            &options,
            &mut statistics,
        )? {
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
        let mut generated_any = false;
        for base in [&base_a, &base_b] {
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
            generated_any = true;
            if let Some(conformer) = locally_optimize(
                objective,
                &connectivity,
                &torsions,
                coordinates,
                generation,
                &mut blacklist,
                &mut optimized_blacklist,
                &rmsd_atoms,
                &options,
                &mut statistics,
            )? {
                population.push(conformer.clone());
                all_conformers.push(conformer);
            }
        }
        if !generated_any {
            termination = TerminationReason::MutationAttemptsExhausted;
            break;
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
        if generation == options.max_generations {
            termination = TerminationReason::MaximumGenerations;
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
