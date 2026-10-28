//! Execute a conformer search on real atoms only.

use super::*;

#[path = "search_engine.rs"]
mod engine;

/// Runs a conformer search. Blocking: call it on a background task.
pub fn search(
    mol: &Molecule,
    mut settings: ConformerSettings,
) -> Result<ConformerOutcome, ConformerError> {
    let (prepared, mapping) =
        crate::calculation_geometry::prepare(mol).map_err(ConformerError::Search)?;
    settings.constraints = super::preparation::remap_constraints(&settings.constraints, &mapping)
        .map_err(ConformerError::Constraints)?;
    let mol = &prepared;
    let started = Instant::now();
    let topology = if settings.engine == SearchEngine::Dreiding {
        Some(DreidingTopology::build(mol).map_err(ConformerError::ForceField)?)
    } else {
        None
    };

    let atoms = mol.atoms.clone();
    let coordinates_bohr: Vec<f64> = mol
        .pos
        .iter()
        .flat_map(|p| {
            [
                p.x as f64 / BOHR_TO_ANGSTROM,
                p.y as f64 / BOHR_TO_ANGSTROM,
                p.z as f64 / BOHR_TO_ANGSTROM,
            ]
        })
        .collect();

    let connectivity = connectivity_model(&atoms);
    let constraint_targets: Vec<_> = settings.constraints.iter().map(|c| c.target()).collect();
    if let Err(problem) = constraints::check(&settings.constraints, &atoms) {
        return Err(ConformerError::Constraints(problem));
    }

    // The rotor count must come from the search's own torsion perception, not from the force
    // field's `rotatable_bonds`. The two disagree on purpose: DREIDING counts every acyclic single
    // bond, while the search excludes methyl rotors, whose rotation returns the same shape. For
    // butane that is three bonds against one -- and budgeting from the larger number asks for a
    // population bigger than the molecule has distinct conformers.
    let degrees_of_freedom = match discover_torsions(
        &coordinates_bohr,
        &atoms,
        &connectivity,
        &[],
        &[],
        true,
        &constraint_targets,
    ) {
        Ok((_bonds, torsions)) => torsions,
        // Perception only refuses when there is nothing to sample at all. Whether that is
        // because the molecule is genuinely rigid or because its flexibility is in a ring too
        // large for the pucker treatment changes what the user should be told.
        Err(_) => return Err(engine::no_torsions(mol, topology.as_ref())),
    };
    let rings = degrees_of_freedom
        .iter()
        .filter(|dof| dof.kind == TorsionKind::RingPucker)
        .count();
    let rotors = degrees_of_freedom.len() - rings;
    if degrees_of_freedom.is_empty() {
        return Err(engine::no_torsions(mol, topology.as_ref()));
    }

    let (requested_population, max_generations) =
        settings.thoroughness.budget(degrees_of_freedom.len());
    let space = search_space(&degrees_of_freedom);
    let population_size = requested_population.min(space).max(2);

    // A small molecule can still have fewer distinct minima than even that estimate -- symmetry
    // makes conformers coincide. Rather than report the search's "could not fill the population"
    // as a failure, come down to what the molecule actually has. Halving converges in a few
    // steps and each attempt on a small molecule is milliseconds.
    // Resolved once, so the report cannot claim a different number from the one used.
    let workers = resolve_workers(settings.ncore);

    let mut population = population_size;
    let mut last_error;
    loop {
        let options = ConformerSearchOptions {
            population_size: population,
            max_generations,
            offspring_per_generation: OFFSPRING_PER_GENERATION,
            local_optimizer: settings.local_optimizer,
            constraints: constraint_targets.clone(),
            random_seed: settings.reproducible.then_some(FIXED_SEED),
            parallel_workers: Some(workers),
            output_energy_window_hartree: settings.energy_window_kcal / HARTREE_TO_KCAL,
            // Methyl rotation does not make a distinct conformer, only a redundant copy.
            exclude_methyl_rotors: true,
            verbosity: 0,
            ..ConformerSearchOptions::default()
        };
        // One objective per worker. For DREIDING each borrows the one shared topology, which is
        // exactly what the topology/workspace split was for. For an external engine each gets its
        // own scratch directory, since xTB is handed fixed file names and two calculations
        // sharing a directory would overwrite each other.
        let worker = std::sync::atomic::AtomicUsize::new(0);
        match genetic_conformer_search_parallel(
            coordinates_bohr.clone(),
            &atoms,
            connectivity.clone(),
            || {
                let slot = worker.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                engine::objective(
                    &settings,
                    topology.as_ref(),
                    &atoms,
                    external_scratch().join(format!("worker_{slot:03}")),
                )
            },
            options,
        ) {
            Ok(result) => {
                let mut outcome = assemble(
                    result,
                    atoms,
                    started.elapsed(),
                    workers,
                    settings.local_optimizer,
                    rings,
                );
                outcome.radical_warning = topology.as_ref().and_then(|t| t.radical_warning());
                outcome.constrained = settings.constraints.len();
                return Ok(outcome);
            }
            Err(error) => {
                last_error = error.to_string();
                if population <= 2 {
                    break;
                }
                population = (population / 2).max(2);
            }
        }
    }

    // A genetic search needs two distinct minima to breed from. Failing at a population of two
    // means the molecule has only one -- true of unsubstituted cyclohexane, whose chair is the
    // only minimum this force field keeps -- and that is an answer, not a failure. Relax the given
    // structure and report it as the single conformer.
    // A constrained search must never fall back to relaxing the structure freely: the geometry is
    // a saddle point, and an unconstrained relaxation would slide it down to reactant or product
    // and report the result as a conformer.
    if last_error.contains("distinct converged minima") && settings.constraints.is_empty() {
        if let Some(topology) = topology.as_ref() {
            return single_conformer(
                topology,
                mol,
                &atoms,
                started.elapsed(),
                workers,
                rotors,
                rings,
            );
        }
        return engine::single_external_conformer(
            &settings,
            &coordinates_bohr,
            &atoms,
            started.elapsed(),
            workers,
            rotors,
            rings,
        );
    }
    Err(ConformerError::Search(last_error))
}
