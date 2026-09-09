use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{ensure, Context, Result};

use super::genetic::{Conformer, ConformerSearchResult, ConformerSearchStatistics};
use super::options::ConformerSearchOptions;
use super::topology::{PreparedTorsion, TorsionKind};

const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;
const HARTREE_TO_KCAL_MOL: f64 = 627.509_474_063_1;

pub(crate) fn print_search_header(
    options: &ConformerSearchOptions,
    torsions: &[PreparedTorsion],
    parallel_workers: usize,
) {
    println!("Genetic conformer search");
    println!("Parallel local jobs:    {parallel_workers} (one objective per job)");
    if parallel_workers > 2 {
        println!(
            "Generation parallelism: at most 2 jobs (the published algorithm creates two children)"
        );
    }
    println!("Population size:       {}", options.population_size);
    println!("Maximum generations:   {}", options.max_generations);
    println!("Selection:             {}", options.selection.label());
    println!(
        "Crossover probability: {:.3}",
        options.crossover_probability
    );
    println!(
        "Torsion mutation:      {:.3} (up to {} torsions)",
        options.torsion_mutation_probability, options.max_torsion_mutations
    );
    println!(
        "Uniqueness threshold:  {:.3} Angstrom heavy-atom RMSD",
        options.rmsd_cutoff_angstrom
    );
    println!(
        "Mirror equivalence:    {}",
        if options.chiral {
            "disabled"
        } else {
            "enabled"
        }
    );
    println!("Torsional degrees of freedom: {}", torsions.len());
    for (index, torsion) in torsions.iter().enumerate() {
        let kind = match torsion.kind {
            TorsionKind::Rotatable => "rotatable",
            TorsionKind::CisTrans => "cis/trans",
            TorsionKind::RingPucker => "ring pucker",
        };
        let selection = if torsion.user_selected {
            "user-selected"
        } else {
            "automatically selected"
        };
        println!(
            "  {:>3}: [{}, {}, {}, {}]  {kind}, {selection} (1-based atom numbering)",
            index + 1,
            torsion.atoms[0] + 1,
            torsion.atoms[1] + 1,
            torsion.atoms[2] + 1,
            torsion.atoms[3] + 1
        );
    }
    if torsions.iter().any(|torsion| !torsion.user_selected) {
        println!(
            "Automatic torsions are inferred from XYZ connectivity; inspect this list and use the explicit `torsions` setting for chemically curated bonds."
        );
    }
    println!();
    println!(
        "{:>5} {:>18} {:>15} {:>10} {:>12} {:>12}",
        "gen", "best E (Eh)", "spread (kcal/mol)", "population", "unique", "local opts"
    );
}

pub(crate) fn print_generation(
    generation: usize,
    population: &[Conformer],
    unique_count: usize,
    statistics: &ConformerSearchStatistics,
) {
    let spread = (population.last().unwrap().energy_hartree
        - population.first().unwrap().energy_hartree)
        * HARTREE_TO_KCAL_MOL;
    println!(
        "{:>5} {:>18.10} {:>15.5} {:>10} {:>12} {:>12}",
        generation,
        population[0].energy_hartree,
        spread,
        population.len(),
        unique_count,
        statistics.local_optimizations
    );
}

pub fn print_conformer_search_result(result: &ConformerSearchResult) {
    let best = result.conformers.first();
    println!();
    println!("Conformer-search result");
    println!("Termination:            {}", result.termination);
    println!("Generations completed:  {}", result.generations);
    println!("Unique minima found:    {}", result.conformers.len());
    println!(
        "Local optimizations:    {} ({} nonconverged)",
        result.statistics.local_optimizations, result.statistics.nonconverged_optimizations
    );
    println!(
        "Rejected trials:        {} steric/connectivity, {} blacklisted",
        result.statistics.rejected_insensible, result.statistics.rejected_blacklisted
    );
    println!(
        "Duplicate minima:       {}",
        result.statistics.duplicate_minima
    );
    if let Some(best) = best {
        println!("------------------------------------------------------");
        println!("Lowest conformer energy: {:>20.12} Eh", best.energy_hartree);
        println!("======================================================");
    }
}

fn write_xyz_frame<W: Write>(
    writer: &mut W,
    elements: &[String],
    conformer: &Conformer,
    rank: usize,
    reference_energy: f64,
) -> Result<()> {
    ensure!(
        conformer.coordinates_bohr.len() == 3 * elements.len(),
        "conformer coordinate count does not match element count"
    );
    let relative = (conformer.energy_hartree - reference_energy) * HARTREE_TO_KCAL_MOL;
    writeln!(writer, "{}", elements.len())?;
    writeln!(
        writer,
        "conformer={} energy={:.12} Eh relative={:.6} kcal/mol generation={} local_cycles={}",
        rank,
        conformer.energy_hartree,
        relative,
        conformer.generation,
        conformer.local_optimization_cycles
    )?;
    for (element, coordinate) in elements
        .iter()
        .zip(conformer.coordinates_bohr.chunks_exact(3))
    {
        writeln!(
            writer,
            "{:<3} {:>18.10} {:>18.10} {:>18.10}",
            element,
            coordinate[0] * BOHR_TO_ANGSTROM,
            coordinate[1] * BOHR_TO_ANGSTROM,
            coordinate[2] * BOHR_TO_ANGSTROM
        )?;
    }
    Ok(())
}

/// Write the low-energy ensemble, the lowest conformer, and a CSV energy list.
pub fn write_conformer_search_files(
    result: &ConformerSearchResult,
    elements: &[String],
    ensemble_path: impl AsRef<Path>,
    lowest_path: impl AsRef<Path>,
    energy_table_path: impl AsRef<Path>,
) -> Result<usize> {
    let best = result
        .conformers
        .first()
        .context("conformer search returned no optimized minimum")?;
    let selected = result
        .conformers
        .iter()
        .take_while(|conformer| {
            conformer.energy_hartree - best.energy_hartree <= result.output_energy_window_hartree
        })
        .collect::<Vec<_>>();

    let mut ensemble = BufWriter::new(File::create(ensemble_path.as_ref()).with_context(|| {
        format!(
            "failed to create conformer ensemble `{}`",
            ensemble_path.as_ref().display()
        )
    })?);
    for (index, conformer) in selected.iter().enumerate() {
        write_xyz_frame(
            &mut ensemble,
            elements,
            conformer,
            index + 1,
            best.energy_hartree,
        )?;
    }

    let mut lowest = BufWriter::new(File::create(lowest_path.as_ref()).with_context(|| {
        format!(
            "failed to create lowest-conformer file `{}`",
            lowest_path.as_ref().display()
        )
    })?);
    write_xyz_frame(&mut lowest, elements, best, 1, best.energy_hartree)?;

    let mut table =
        BufWriter::new(File::create(energy_table_path.as_ref()).with_context(|| {
            format!(
                "failed to create conformer energy table `{}`",
                energy_table_path.as_ref().display()
            )
        })?);
    writeln!(
        table,
        "rank,energy_hartree,relative_kcal_mol,generation,local_optimization_cycles"
    )?;
    for (index, conformer) in result.conformers.iter().enumerate() {
        writeln!(
            table,
            "{},{:.12},{:.8},{},{}",
            index + 1,
            conformer.energy_hartree,
            (conformer.energy_hartree - best.energy_hartree) * HARTREE_TO_KCAL_MOL,
            conformer.generation,
            conformer.local_optimization_cycles
        )?;
    }
    Ok(selected.len())
}
