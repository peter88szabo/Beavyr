//! Conformer search: find the low-energy shapes a flexible molecule can take.
//!
//! The search itself is Behemoth's torsion-space genetic algorithm
//! ([`crate::optimizer::conformer_search`], following Supady, Blum and Baldauf, *J. Chem. Inf.
//! Model.* **2015**, *55*, 2338; local copy at
//! `~/Dropbox/Research_Leuven/Behemoth/Papers/ConfoSearch/Genetic/ci5b00243.pdf`). What this
//! module adds is everything needed to drive it from a structure on screen: the energy source,
//! automatic settings, and a background task.
//!
//! # Why the force field has to be in-process
//!
//! A genetic search runs thousands of local optimisations, each needing tens of gradient
//! evaluations -- so 10⁴ to 10⁵ energy calls for a molecule with a handful of rotors. The DREIDING
//! force field answers each in microseconds, in this process. An external program answers each in
//! a process launch: xTB's GFN2 is fast once running, but at ~100 ms of startup per call the
//! same search would take hours to days rather than seconds.
//!
//! That is why the search runs on DREIDING, and why xTB is offered as a *refinement* of the
//! conformers that survive instead. Re-optimising the best few dozen is a few hundred xTB calls --
//! a minute or two -- and it buys a markedly better energy ordering than a generic force field
//! gives, because relative conformer energies are where a generic force field is weakest.

pub mod afir;
pub mod constraints;
pub mod refine;
pub mod tstrail;
pub mod ui;

use std::time::{Duration, Instant};

use bevy::prelude::*;
use bevy::tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task};

use crate::forcefield::dreiding::objective::{
    cleanup_options, relax_angstrom, DreidingObjective, BOHR_TO_ANGSTROM, HARTREE_TO_KCAL,
};
use crate::forcefield::dreiding::{BuildError, DreidingTopology};
use crate::molecule::{covalent_radius_angstrom, Molecule};
use crate::optimizer::conformer_search::{
    discover_torsions, genetic_conformer_search_parallel, ConformerSearchOptions,
    ConformerSearchResult, LocalOptimizer, PreparedTorsion, TorsionKind,
};
use crate::optimizer::internal_coords::ConnectivityModel;

/// How hard to search. The only knob offered, because everything else can be derived.
///
/// The counts scale with the number of rotatable bonds the search finds, since a molecule with two
/// rotors needs a fraction of the sampling one with eight does. Behemoth's own defaults (a
/// population of 5 over 10 generations) are tuned for quantum-chemistry energies where every
/// evaluation is expensive; against DREIDING they are needlessly timid, so these are far more
/// generous.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Thoroughness {
    /// A quick look, for checking whether a molecule is flexible at all.
    Quick,
    #[default]
    Normal,
    /// For a molecule with many rotors, or when the answer matters.
    Thorough,
}

impl Thoroughness {
    pub fn label(self) -> &'static str {
        match self {
            Self::Quick => "Quick",
            Self::Normal => "Normal",
            Self::Thorough => "Thorough",
        }
    }

    pub fn description(self) -> &'static str {
        match self {
            Self::Quick => "A small population over few generations. Seconds.",
            Self::Normal => "The usual choice.",
            Self::Thorough => "Large population, many generations. Best for many rotors.",
        }
    }

    /// Population size and generation count, given how many rotors were found.
    ///
    /// Each generation breeds [`OFFSPRING_PER_GENERATION`] children, so the total number of local
    /// optimisations is roughly `population + generations × offspring`.
    fn budget(self, rotors: usize) -> (usize, usize) {
        let rotors = rotors.max(1);
        match self {
            Self::Quick => ((6 * rotors).clamp(12, 40), 12),
            Self::Normal => ((10 * rotors).clamp(20, 80), 30),
            Self::Thorough => ((16 * rotors).clamp(32, 150), 60),
        }
    }
}

/// How many children each generation breeds, and so how wide the search can run.
///
/// Behemoth breeds one crossover pair -- two children -- per generation, which caps the
/// generation loop at two concurrent optimisations no matter how many cores are available. Eight
/// keeps a machine of that size busy while leaving enough generations for selection to actually
/// work: a genetic algorithm that is all width and no depth is just random sampling.
///
/// Deliberately a constant rather than the core count. If the generation width followed `ncore`,
/// the same molecule and seed would breed differently on different machines, and a conformer
/// someone found could not be reproduced. This way `ncore` changes only how fast the generation
/// is worked through.
pub const OFFSPRING_PER_GENERATION: usize = 8;

/// What the user chooses. Deliberately short: rotors, population and generations are all worked
/// out from the structure rather than asked for.
#[derive(Debug, Clone, PartialEq)]
pub struct ConformerSettings {
    pub thoroughness: Thoroughness,
    /// Fix the random seed, so a search can be repeated exactly.
    pub reproducible: bool,
    /// Which optimiser relaxes each candidate conformer.
    ///
    /// Cartesian by default. The search's cost is almost entirely the coordinate algebra rather
    /// than the force field -- for a 32-atom molecule an internal-coordinate cycle spends roughly
    /// 6.7 Mflop building and inverting matrices against 0.01 Mflop on the energy -- so the
    /// optimiser that does least per cycle wins by a wide margin, even needing more cycles.
    pub local_optimizer: LocalOptimizer,
    /// How many structures to optimise at once.
    ///
    /// A conformer search is embarrassingly parallel: every local optimisation is independent, so
    /// this scales almost linearly until it runs out of cores. It does **not** change the answer
    /// -- the search reserves each starting structure before handing it out and collects results
    /// in the order it dispatched them, so a run with a fixed seed gives the same conformers on
    /// any number of cores.
    pub ncore: usize,
    /// Keep conformers within this window of the lowest, in kcal/mol.
    pub energy_window_kcal: f64,
    /// Coordinates held fixed, for searching the conformers of a transition state.
    ///
    /// Empty for an ordinary search. When set, every candidate is relaxed by the constrained
    /// optimiser and no degree of freedom the constraints already hold is sampled -- see
    /// [`constraints`].
    pub constraints: Vec<constraints::Constraint>,
}

impl Default for ConformerSettings {
    fn default() -> Self {
        Self {
            thoroughness: Thoroughness::default(),
            reproducible: true,
            // Not Behemoth's default. Its internal-coordinate optimiser is right when an energy
            // costs seconds; here an energy costs microseconds and the trade inverts.
            local_optimizer: LocalOptimizer::Cartesian,
            ncore: available_cores(),
            // Wide enough to keep everything thermally accessible, and then some: at room
            // temperature 10 kcal/mol is already a population of ~1e-7.
            energy_window_kcal: 10.0,
            constraints: Vec::new(),
        }
    }
}

/// One conformer, in the units a chemist reads.
#[derive(Debug, Clone)]
pub struct Conformer {
    /// Flat `[x0, y0, z0, ...]` in Å.
    pub positions_angstrom: Vec<f64>,
    /// kcal/mol above the lowest conformer found.
    pub relative_energy_kcal: f64,
    /// Fraction of molecules expected in this conformer at 298.15 K, from the relative energies
    /// alone. Ignores vibrational entropy and any symmetry factor, so it is indicative only.
    pub population_fraction: f64,
}

/// A finished search.
#[derive(Debug, Clone)]
pub struct ConformerOutcome {
    pub conformers: Vec<Conformer>,
    /// Rotatable bonds the search detected and sampled.
    pub rotors: usize,
    /// Rings whose pucker the search sampled, five- and six-membered.
    pub rings: usize,
    pub generations: usize,
    pub structures_tried: usize,
    pub local_optimizations: usize,
    pub termination: String,
    pub elapsed: Duration,
    /// How many structures were optimised at once.
    pub workers: usize,
    /// Which optimiser relaxed each candidate.
    pub local_optimizer: LocalOptimizer,
    /// The generation that produced the lowest-energy conformer.
    ///
    /// Worth showing: if the best structure only turned up in the last generation, the search was
    /// still improving when it stopped and more effort would probably find something better.
    pub best_found_in_generation: usize,
    /// Atom symbols, so a conformer can be written out or loaded back.
    pub atoms: Vec<String>,
    /// Set when the structure looks open-shell, saying what that means for the result.
    pub radical_warning: Option<String>,
    /// How many coordinates were held fixed. Non-zero means these are transition-state
    /// conformers, whose energies are saddle-point energies rather than minima.
    pub constrained: usize,
    /// Set once the leading conformers have been re-optimised and re-ranked with GFN2.
    ///
    /// When set, the energies of the refined conformers are GFN2 and any beyond the refinement
    /// cap are still DREIDING -- so the list is honestly mixed, and says so.
    pub refined_with_xtb: bool,
}

/// How many cores the machine has, as the default worker count.
///
/// A search is the one thing in Beavyr worth taking the whole machine for -- it is pure compute
/// with no interface to keep responsive, and it finishes sooner the more it is given. The field in
/// the panel is there to give some back when something else needs it.
pub fn available_cores() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .max(1)
}

/// Roughly how many distinct shapes the sampled degrees of freedom can make between them.
///
/// Used only to cap the population: the search needs that many *distinct* starting minima, and
/// asking for more than exist cannot be satisfied.
///
/// Each kind of degree of freedom contributes what it actually offers. A rotatable bond has about
/// three staggered minima and a cis/trans bond two, but a ring has as many canonical conformers as
/// its size allows -- 38 for a six-ring, 20 for a five-ring. Treating a ring as though it had
/// three, which an earlier version did, capped methylcyclohexane's population at three and found
/// one chair where there are two.
fn search_space(degrees_of_freedom: &[PreparedTorsion]) -> usize {
    // Capped so a floppy molecule does not overflow or ask for an absurd population; the
    // thoroughness budget is the real limit on size.
    const CEILING: usize = 100_000;
    degrees_of_freedom
        .iter()
        .fold(1usize, |total, dof| {
            let multiplicity = match dof.kind {
                TorsionKind::Rotatable => 3,
                TorsionKind::CisTrans => 2,
                TorsionKind::RingPucker => dof
                    .ring
                    .as_ref()
                    .map(|ring| ring.conformers.len())
                    .unwrap_or(1),
            };
            total.saturating_mul(multiplicity.max(1))
        })
        .min(CEILING)
}

/// The result for a molecule that turns out to have exactly one minimum.
///
/// Reached when the search cannot assemble two distinct starting minima to breed from. The
/// molecule is flexible on paper -- it has rotors or a ring -- but every trial geometry relaxes to
/// the same place, so there is one conformer and this reports it rather than an error.
fn single_conformer(
    topology: &DreidingTopology,
    mol: &Molecule,
    atoms: &[String],
    elapsed: Duration,
    workers: usize,
    rotors: usize,
    rings: usize,
) -> Result<ConformerOutcome, ConformerError> {
    let start: Vec<f64> = mol
        .pos
        .iter()
        .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
        .collect();
    let (relaxed, _) = relax_angstrom(topology, &start, cleanup_options())
        .map_err(|e| ConformerError::Search(e.to_string()))?;

    Ok(ConformerOutcome {
        conformers: vec![Conformer {
            positions_angstrom: relaxed,
            relative_energy_kcal: 0.0,
            population_fraction: 1.0,
        }],
        rotors,
        rings,
        generations: 0,
        structures_tried: 0,
        local_optimizations: 1,
        termination: "only one distinct minimum".into(),
        elapsed,
        workers,
        local_optimizer: LocalOptimizer::Cartesian,
        constrained: 0,
        best_found_in_generation: 0,
        atoms: atoms.to_vec(),
        radical_warning: topology.radical_warning(),
        refined_with_xtb: false,
    })
}

/// Why a molecule offered no torsion to sample: genuinely rigid, or flexible only in a ring.
///
/// The distinction matters to the user. "One conformer" is right for water and wrong for
/// cyclohexane, and a chemist told the latter would rightly not believe the tool.
fn nothing_to_rotate(topology: &DreidingTopology) -> ConformerError {
    // Reaching here with ring bonds present means none of the rings was five- or six-membered,
    // since those would have become degrees of freedom in their own right.
    let rings = topology.ring_bond_count();
    if rings > 0 {
        ConformerError::RingOnly { rings }
    } else {
        ConformerError::Rigid
    }
}

/// The worker count a request resolves to: never more than the machine has, never fewer than one.
pub fn resolve_workers(requested: usize) -> usize {
    requested.clamp(1, available_cores())
}

/// Room temperature, for the population estimate.
const ROOM_TEMPERATURE_K: f64 = 298.15;

/// kcal/mol per K, the gas constant.
const R_KCAL: f64 = 0.001_987_204_1;

/// Boltzmann populations from relative energies at room temperature.
pub(crate) fn populations(relative_kcal: &[f64]) -> Vec<f64> {
    let weights: Vec<f64> = relative_kcal
        .iter()
        .map(|e| (-e / (R_KCAL * ROOM_TEMPERATURE_K)).exp())
        .collect();
    let total: f64 = weights.iter().sum();
    if total <= 0.0 {
        return vec![0.0; relative_kcal.len()];
    }
    weights.iter().map(|w| w / total).collect()
}

/// Runs a conformer search. Blocking: call it on a background task.
pub fn search(
    mol: &Molecule,
    settings: ConformerSettings,
) -> Result<ConformerOutcome, ConformerError> {
    let started = Instant::now();
    let topology = DreidingTopology::build(mol).map_err(ConformerError::ForceField)?;

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
    let constraint_targets: Vec<_> = settings
        .constraints
        .iter()
        .map(|c| c.target())
        .collect();
    if let Err(problem) = constraints::check(&settings.constraints, &atoms) {
        return Err(ConformerError::Constraints(problem));
    }

    // The rotor count must come from the search's own torsion perception, not from the force
    // field's `rotatable_bonds`. The two disagree on purpose: DREIDING counts every acyclic single
    // bond, while the search excludes methyl rotors, whose rotation returns the same shape. For
    // butane that is three bonds against one -- and budgeting from the larger number asks for a
    // population bigger than the molecule has distinct conformers.
    let degrees_of_freedom =
        match discover_torsions(
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
            Err(_) => return Err(nothing_to_rotate(&topology)),
        };
    let rings = degrees_of_freedom
        .iter()
        .filter(|dof| dof.kind == TorsionKind::RingPucker)
        .count();
    let rotors = degrees_of_freedom.len() - rings;
    if degrees_of_freedom.is_empty() {
        return Err(nothing_to_rotate(&topology));
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
        // One objective per worker, each borrowing the one shared topology. This is exactly what
        // the topology/workspace split was for: the parameters and term lists are immutable and
        // shared, while the neighbour list and scratch buffers live per worker, so N threads
        // optimise N structures without contending on anything.
        match genetic_conformer_search_parallel(
            coordinates_bohr.clone(),
            &atoms,
            connectivity.clone(),
            || Ok(DreidingObjective::new(&topology)),
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
                outcome.radical_warning = topology.radical_warning();
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
        return single_conformer(
            &topology,
            mol,
            &atoms,
            started.elapsed(),
            workers,
            rotors,
            rings,
        );
    }
    Err(ConformerError::Search(last_error))
}

/// Coordination-number connectivity for the search's own bond perception.
///
/// The search finds its torsions from this rather than from Beavyr's displayed bond list, using
/// the smooth counting function of `analyze_structure`: a pair counts as bonded when
/// `1/(1 + exp(-k(R⁰/R − 1)))` exceeds `facmin`.
///
/// **The radii must carry the D3 scaling.** This is the detail that matters: `ConnectivityModel`
/// documents its `rcov_disp` as already scaled, and with `k = 7.5` and `facmin = 0.6` the
/// criterion is calibrated for that. Feeding it plain covalent radii puts an ordinary C-C single
/// bond at `fac = 0.494`, just under the threshold, so *no* carbon skeleton is perceived and the
/// search reports that a flexible molecule has no rotatable bond. Scaled, the same bond sits at
/// 0.922 and a 1-3 contact at 0.195, which is the separation the threshold expects.
///
/// (Behemoth's `optimizer::irc::default_connectivity_model` builds this same struct from unscaled
/// radii. It is only a fallback there -- IRC normally receives a properly scaled model from the
/// xTB Hamiltonian -- but it is worth knowing about.)
fn connectivity_model(atoms: &[String]) -> ConnectivityModel {
    ConnectivityModel {
        rcov_disp: atoms
            .iter()
            .map(|a| covalent_radius_angstrom(a) as f64 * D3_RADIUS_SCALE / BOHR_TO_ANGSTROM)
            .collect(),
        kcn: 7.5,
        facmin: 0.6,
    }
}

/// The scaling D3 applies to covalent radii before using them in a coordination number.
const D3_RADIUS_SCALE: f64 = 4.0 / 3.0;

/// The seed used when a search is asked to be reproducible.
///
/// Any fixed value does; this one is arbitrary. What matters is that the same structure and
/// settings give the same conformers, so that a result can be reproduced when someone reports
/// that a search missed something.
const FIXED_SEED: u64 = 0x_C0FF_EE_15_600D;

fn assemble(
    result: ConformerSearchResult,
    atoms: Vec<String>,
    elapsed: Duration,
    workers: usize,
    local_optimizer: LocalOptimizer,
    rings: usize,
) -> ConformerOutcome {
    let lowest = result
        .conformers
        .iter()
        .map(|c| c.energy_hartree)
        .fold(f64::INFINITY, f64::min);

    let relative: Vec<f64> = result
        .conformers
        .iter()
        .map(|c| (c.energy_hartree - lowest) * HARTREE_TO_KCAL)
        .collect();
    let fractions = populations(&relative);

    let conformers = result
        .conformers
        .iter()
        .zip(relative)
        .zip(fractions)
        .map(|((c, relative_energy_kcal), population_fraction)| Conformer {
            positions_angstrom: c
                .coordinates_bohr
                .iter()
                .map(|b| b * BOHR_TO_ANGSTROM)
                .collect(),
            relative_energy_kcal,
            population_fraction,
        })
        .collect();

    ConformerOutcome {
        conformers,
        rotors: result.torsions.len().saturating_sub(rings),
        rings,
        generations: result.generations,
        structures_tried: result.statistics.attempted_structures,
        local_optimizations: result.statistics.local_optimizations,
        termination: format!("{:?}", result.termination),
        elapsed,
        workers,
        local_optimizer,
        constrained: 0,
        radical_warning: None,
        refined_with_xtb: false,
        best_found_in_generation: result
            .conformers
            .iter()
            .min_by(|a, b| a.energy_hartree.total_cmp(&b.energy_hartree))
            .map(|c| c.generation)
            .unwrap_or(0),
        atoms,
    }
}

/// Why a search could not run or finish.
#[derive(Debug)]
pub enum ConformerError {
    /// DREIDING could not describe the molecule. Its message says which atom and why.
    ForceField(BuildError),
    /// Nothing to search: no rotatable bond, so the molecule has one shape.
    Rigid,
    /// The constraints cannot be used as given.
    Constraints(String),
    /// The molecule's only flexibility is in a ring of a size this search does not pucker.
    ///
    /// Five- and six-membered rings *are* sampled, through their Cremer-Pople puckering
    /// coordinates -- see [`crate::optimizer::conformer_search::rings`]. Other sizes are not:
    /// three- and four-rings are effectively rigid, and a larger ring's conformational space is
    /// not the small sphere or circle that treatment relies on. Since no ring bond can be turned
    /// on its own without prising the ring open, such a molecule has flexibility that goes
    /// unsampled, and saying "one conformer" would be wrong.
    RingOnly { rings: usize },
    Search(String),
}

impl std::fmt::Display for ConformerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ForceField(e) => write!(f, "{e}"),
            Self::Rigid => write!(
                f,
                "this molecule has no rotatable bonds, so it has only one conformer. \
                 Double bonds do not rotate, and a methyl group's rotation gives the same \
                 shape back."
            ),
            Self::RingOnly { rings } => write!(
                f,
                "this molecule's only flexibility is in its {}, and not of a size this search \
                 can pucker. Five- and six-membered rings are sampled -- chair against \
                 twist-boat, or a sugar's puckers -- but a three- or four-ring is effectively \
                 rigid, and a larger ring needs a different treatment. Turning a single ring \
                 bond is not an option: it would prise the ring open rather than give another \
                 conformer.",
                if *rings == 1 { "ring" } else { "rings" }
            ),
            Self::Constraints(message) => write!(f, "the constraints are not usable: {message}"),
            Self::Search(message) => write!(f, "the conformer search failed: {message}"),
        }
    }
}

/// The running search, and the last result. A Bevy resource so the panel can poll it.
#[derive(Resource, Default)]
pub struct ConformerRun {
    task: Option<Task<Result<ConformerOutcome, ConformerError>>>,
    /// The GFN2 re-ranking, which runs after a search rather than inside it.
    ///
    /// Separate because it is a different kind of wait: a search is seconds and a refinement is
    /// minutes of external process launches, so the user starts it deliberately and watches the
    /// results they already have while it runs.
    refine_task: Option<Task<Result<(ConformerOutcome, refine::Refinement), String>>>,
    /// What the last refinement did, or why it could not run.
    pub refinement: Option<Result<refine::Refinement, String>>,
    pub settings: ConformerSettings,
    pub outcome: Option<ConformerOutcome>,
    pub error: Option<String>,
    /// Which conformer is currently shown in the viewer, if any.
    pub previewing: Option<usize>,
    started: Option<Instant>,
}

impl ConformerRun {
    pub fn is_running(&self) -> bool {
        self.task.is_some()
    }

    pub fn is_refining(&self) -> bool {
        self.refine_task.is_some()
    }

    /// Starts re-ranking the current conformers with xTB's GFN2.
    ///
    /// `binary` comes from the path the Geometry Optimization panel already keeps for xTB, so
    /// there is no second place to configure it.
    pub fn start_refinement(&mut self, binary: std::path::PathBuf) {
        if self.is_running() || self.is_refining() {
            return;
        }
        let Some(outcome) = self.outcome.clone() else {
            return;
        };
        self.refinement = None;

        let scratch = crate::qchem_interfaces::xtb_optimize::xtb_scratch_dir()
            .join(format!("conformers_{}", std::process::id()));
        self.refine_task = Some(AsyncComputeTaskPool::get().spawn(async move {
            let mut outcome = outcome;
            match refine::refine_with_xtb(&mut outcome, &binary, &scratch, 0, 1) {
                Ok(report) => {
                    // The scratch tree holds one directory per conformer and is of no further
                    // use once the energies are read back.
                    let _ = std::fs::remove_dir_all(&scratch);
                    Ok((outcome, report))
                }
                Err(e) => {
                    let _ = std::fs::remove_dir_all(&scratch);
                    Err(e.to_string())
                }
            }
        }));
    }

    pub fn elapsed(&self) -> Option<Duration> {
        self.started.map(|t| t.elapsed())
    }

    /// Starts a search on a background task, so a long one cannot stall the interface.
    pub fn start(&mut self, mol: &Molecule) {
        if self.is_running() {
            return;
        }
        self.error = None;
        self.outcome = None;
        self.previewing = None;
        self.started = Some(Instant::now());

        // The molecule is copied rather than borrowed: the search runs for seconds, and the user
        // is free to keep editing meanwhile.
        let snapshot = mol.clone();
        let settings = self.settings.clone();
        self.task = Some(
            AsyncComputeTaskPool::get().spawn(async move { search(&snapshot, settings) }),
        );
    }
}

/// Collects a finished search. Registered by [`crate::conformer::ui`].
pub fn poll_conformer_search(mut run: ResMut<ConformerRun>) {
    if let Some(task) = run.task.as_mut() {
        if let Some(result) = block_on(future::poll_once(task)) {
            run.task = None;
            match result {
                Ok(outcome) => run.outcome = Some(outcome),
                Err(error) => run.error = Some(error.to_string()),
            }
        }
    }
    if let Some(task) = run.refine_task.as_mut() {
        if let Some(result) = block_on(future::poll_once(task)) {
            run.refine_task = None;
            match result {
                Ok((outcome, report)) => {
                    run.outcome = Some(outcome);
                    run.refinement = Some(Ok(report));
                    // The ranking changed, so whichever conformer was on screen may no longer be
                    // the row the user is looking at.
                    run.previewing = None;
                }
                Err(error) => run.refinement = Some(Err(error)),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn molecule(xyz: &str) -> Molecule {
        let mut mol = Molecule::from_xyz(xyz);
        mol.recompute_bonds(1.2, 2.5);
        mol
    }

    /// Cores a test search may use.
    ///
    /// `ConformerSettings::default()` takes the whole machine, which is right for the application
    /// and wrong for a test suite: several searches running at once would each try to fill every
    /// core. Tests take this instead, so the suite stays polite on a machine doing other work.
    const TEST_CORES: usize = 2;

    /// Settings for a test search: quick, repeatable, and modest about cores.
    fn test_settings() -> ConformerSettings {
        ConformerSettings {
            thoroughness: Thoroughness::Quick,
            reproducible: true,
            ncore: TEST_CORES,
            ..ConformerSettings::default()
        }
    }

    /// n-Butane: one central C-C rotor once the two methyls are excluded, and the textbook case
    /// of anti versus gauche.
    const BUTANE: &str = "14\n\n\
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
        H  1.9600 -1.2100 -0.8900\n";

    #[test]
    fn a_rigid_molecule_is_reported_not_searched() {
        // Water has no rotatable bond at all.
        let water = molecule("3\n\nO 0.000 0.000 0.000\nH 0.960 0.000 0.000\nH -0.240 0.929 0.000\n");
        let error = search(&water, test_settings())
            .expect_err("water has one shape");
        assert!(matches!(error, ConformerError::Rigid), "{error:?}");
        // The message must explain why rather than just refuse.
        assert!(error.to_string().contains("rotatable"));
    }

    /// A molecule DREIDING cannot type must fail with the force field's own message, naming the
    /// atom, rather than a generic search failure.
    #[test]
    fn an_unparameterised_molecule_reports_the_force_field_error() {
        let lif = molecule("2\n\nLi 0.000 0.000 0.000\nF 0.000 0.000 1.564\n");
        match search(&lif, test_settings()) {
            Err(ConformerError::ForceField(_)) => {}
            other => panic!("expected a force-field error, got {other:?}"),
        }
    }

    /// The connectivity radii must carry the D3 scaling, or the smooth counting function puts an
    /// ordinary C-C single bond just below its threshold and no carbon skeleton is perceived at
    /// all -- after which a flexible molecule is reported as rigid.
    #[test]
    fn connectivity_perceives_a_normal_carbon_skeleton() {
        use crate::optimizer::internal_coords::analyze_structure;

        let butane = molecule(BUTANE);
        let atoms = butane.atoms.clone();
        let coordinates_bohr: Vec<f64> = butane
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

        let internals = analyze_structure(&coordinates_bohr, &connectivity_model(&atoms))
            .expect("butane is a valid structure");
        // Butane has 13 bonds: three C-C and ten C-H.
        assert_eq!(internals.bonds.len(), 13, "bonds: {:?}", internals.bonds);

        let carbon_carbon = internals
            .bonds
            .iter()
            .filter(|[i, j]| atoms[*i] == "C" && atoms[*j] == "C")
            .count();
        assert_eq!(carbon_carbon, 3, "the carbon skeleton was not perceived");
    }

    /// The rotor count must match the search's own torsion perception, not the force field's
    /// broader `rotatable_bonds`. Butane has three acyclic C-C single bonds but only one rotor
    /// that changes its shape, the other two being methyl rotations.
    #[test]
    fn methyl_rotors_do_not_count_towards_the_budget() {
        let butane = molecule(BUTANE);
        let topology = DreidingTopology::build(&butane).unwrap();
        assert_eq!(
            topology.rotatable_bonds().len(),
            3,
            "the force field counts every acyclic single bond"
        );

        let outcome = search(
            &butane,
            ConformerSettings {
                thoroughness: Thoroughness::Quick,
                ..test_settings()
            },
        )
        .expect("butane is searchable");
        assert_eq!(
            outcome.rotors, 1,
            "the search should sample one rotor, the methyls excluded"
        );
    }

    #[test]
    fn butane_search_finds_more_than_one_conformer() {
        let butane = molecule(BUTANE);
        let outcome = search(&butane, ConformerSettings {
            thoroughness: Thoroughness::Quick,
            ..test_settings()
        })
        .expect("butane is searchable");

        assert!(outcome.rotors >= 1, "no rotor found in butane");
        assert!(!outcome.conformers.is_empty());
        // The lowest is the reference, so it must sit at zero.
        assert!(outcome.conformers[0].relative_energy_kcal.abs() < 1.0e-9);
        // Energies must come back ordered, since the table shows them as ranked.
        for pair in outcome.conformers.windows(2) {
            assert!(pair[0].relative_energy_kcal <= pair[1].relative_energy_kcal + 1.0e-9);
        }
        // Every conformer must have a full set of coordinates.
        for c in &outcome.conformers {
            assert_eq!(c.positions_angstrom.len(), butane.atoms.len() * 3);
        }
    }

    /// A reproducible search must give the same answer twice, which is the point of fixing the
    /// seed -- otherwise a report that a search missed a conformer cannot be chased.
    #[test]
    fn a_reproducible_search_repeats_exactly() {
        let butane = molecule(BUTANE);
        let settings = ConformerSettings {
            thoroughness: Thoroughness::Quick,
            reproducible: true,
            ..test_settings()
        };
        let first = search(&butane, settings.clone()).expect("searchable");
        let second = search(&butane, settings).expect("searchable");

        assert_eq!(first.conformers.len(), second.conformers.len());
        for (a, b) in first.conformers.iter().zip(&second.conformers) {
            assert!((a.relative_energy_kcal - b.relative_energy_kcal).abs() < 1.0e-9);
        }
    }

    /// The core count must not change the answer. The panel says so, and it is the whole reason
    /// parallelising a search is safe: the search reserves each starting structure before handing
    /// it to a worker, and collects results in the order it dispatched them, so nothing depends on
    /// which thread finished first.
    ///
    /// Kept to a couple of workers rather than the whole machine so the test itself stays polite.
    #[test]
    fn the_core_count_does_not_change_the_result() {
        let butane = molecule(BUTANE);
        let run = |ncore: usize| {
            search(
                &butane,
                ConformerSettings {
                    thoroughness: Thoroughness::Quick,
                    reproducible: true,
                    ncore,
                    ..test_settings()
                },
            )
            .expect("butane is searchable")
        };

        let serial = run(1);
        let parallel = run(3);

        assert_eq!(serial.workers, 1);
        assert_eq!(parallel.workers, 3.min(available_cores()));
        assert_eq!(
            serial.conformers.len(),
            parallel.conformers.len(),
            "different number of conformers on a different core count"
        );
        for (a, b) in serial.conformers.iter().zip(&parallel.conformers) {
            assert!(
                (a.relative_energy_kcal - b.relative_energy_kcal).abs() < 1.0e-9,
                "energies differ: {} vs {}",
                a.relative_energy_kcal,
                b.relative_energy_kcal
            );
        }
    }

    /// A request for more cores than the machine has must be clamped rather than oversubscribing,
    /// and zero must not mean "no workers".
    ///
    /// Tested on the resolver directly rather than by running a search: proving the clamp by
    /// actually taking every core is the one thing the clamp exists to prevent, and a test suite
    /// has no business doing it on a machine that is busy with something else.
    #[test]
    fn the_worker_count_is_clamped_to_the_machine() {
        let cores = available_cores();
        assert_eq!(resolve_workers(4096), cores, "a huge request must be clamped");
        assert_eq!(resolve_workers(0), 1, "zero cores must still run");
        assert_eq!(resolve_workers(1), 1);
        if cores > 1 {
            assert_eq!(resolve_workers(cores - 1), cores - 1, "a modest request stands");
        }
    }

    /// n-Decane, all-anti: 32 atoms and seven rotatable bonds once the methyls are excluded, so
    /// a genuinely flexible molecule rather than a toy.
    const DECANE: &str = "32\n\
        n-decane, all-anti\n\
        C 0.0000 0.0000 0.0000\n\
        C 1.2684 0.8556 0.0000\n\
        C 2.5369 0.0000 0.0000\n\
        C 3.8053 0.8556 0.0000\n\
        C 5.0737 0.0000 0.0000\n\
        C 6.3421 0.8556 0.0000\n\
        C 7.6106 0.0000 0.0000\n\
        C 8.8790 0.8556 0.0000\n\
        C 10.1474 0.0000 0.0000\n\
        C 11.4158 0.8556 0.0000\n\
        H -0.8756 0.6491 0.0000\n\
        H -0.0135 -0.6290 -0.8901\n\
        H -0.0135 -0.6290 0.8901\n\
        H 1.2684 1.4006 0.9439\n\
        H 1.2684 1.4006 -0.9439\n\
        H 2.5369 -0.5450 0.9439\n\
        H 2.5369 -0.5450 -0.9439\n\
        H 3.8053 1.4006 0.9439\n\
        H 3.8053 1.4006 -0.9439\n\
        H 5.0737 -0.5450 0.9439\n\
        H 5.0737 -0.5450 -0.9439\n\
        H 6.3421 1.4006 0.9439\n\
        H 6.3421 1.4006 -0.9439\n\
        H 7.6106 -0.5450 0.9439\n\
        H 7.6106 -0.5450 -0.9439\n\
        H 8.8790 1.4006 0.9439\n\
        H 8.8790 1.4006 -0.9439\n\
        H 10.1474 -0.5450 0.9439\n\
        H 10.1474 -0.5450 -0.9439\n\
        H 12.2915 0.2064 0.0000\n\
        H 11.4294 1.4846 -0.8901\n\
        H 11.4294 1.4846 0.8901\n";

    /// Measures what the two speed choices actually buy. Ignored by default -- it is a timing
    /// measurement, which has no business failing a test suite on a loaded machine -- but kept in
    /// the tree so the numbers can be re-checked, and **must be run in release**: a debug build
    /// makes the matrix algebra look far worse than it is.
    ///
    /// ```text
    /// cargo test --release --bins conformer_search_speed -- --ignored --nocapture
    /// ```
    #[test]
    #[ignore = "timing measurement, not a correctness check"]
    fn conformer_search_speed() {
        let decane = molecule(DECANE);
        let settings = |method: LocalOptimizer, ncore: usize| ConformerSettings {
            thoroughness: Thoroughness::Normal,
            reproducible: true,
            local_optimizer: method,
            ncore,
            ..test_settings()
        };

        let mut rows = Vec::new();
        // Only the Cartesian route is timed end to end. A whole search on internal coordinates
        // takes minutes on this molecule, and `one_local_optimisation_each_way` already measures
        // that difference in isolation, which is the cleaner comparison anyway.
        for (method, ncore) in [
            (LocalOptimizer::Cartesian, 1),
            (LocalOptimizer::Cartesian, 4),
        ] {
            crate::optimizer::conformer_search::local::reset_accounting();
            let outcome = search(&decane, settings(method, ncore)).expect("decane is searchable");
            let (in_optimiser, optimisations) =
                crate::optimizer::conformer_search::local::accounting();
            println!(
                "  {:<22} {} core{}: {:>7.2} s   {} tried, {} optimised, {} conformers",
                method.label(),
                ncore,
                if ncore == 1 { " " } else { "s" },
                outcome.elapsed.as_secs_f64(),
                outcome.structures_tried,
                outcome.local_optimizations,
                outcome.conformers.len()
            );
            // With several workers the optimiser time is a sum over threads, so divide by the
            // worker count to compare it against wall-clock.
            let optimiser_share = in_optimiser.as_secs_f64() / ncore as f64;
            println!(
                "  {:>28} {:.2} s in the optimiser ({:.0}% of wall-clock), {:.1} ms each over {} calls",
                "",
                optimiser_share,
                100.0 * optimiser_share / outcome.elapsed.as_secs_f64().max(1e-9),
                1000.0 * in_optimiser.as_secs_f64() / optimisations.max(1) as f64,
                optimisations
            );
            rows.push((method, ncore, outcome.elapsed.as_secs_f64()));
        }

        println!();
        for (_, ncore, elapsed) in &rows[1..] {
            println!(
                "  {ncore} cores vs 1: {:.2}x",
                rows[0].2 / elapsed.max(1e-9)
            );
        }
    }

    /// Times a single local optimisation each way, which is the claim that matters: the whole
    /// search is thousands of these, so their ratio is the search's ratio. Ignored, and must be
    /// run in release.
    ///
    /// ```text
    /// cargo test --release --bins one_local_optimisation -- --ignored --nocapture
    /// ```
    #[test]
    #[ignore = "timing measurement, not a correctness check"]
    fn one_local_optimisation_each_way() {
        use crate::forcefield::dreiding::objective::DreidingObjective;
        use crate::optimizer::conformer_search::local::local_optimize;
        use crate::optimizer::geom_opt::GeomOptOptions;

        let decane = molecule(DECANE);
        let topology = DreidingTopology::build(&decane).unwrap();
        let coords: Vec<f64> = decane
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
        let connectivity = connectivity_model(&decane.atoms);
        let options = GeomOptOptions {
            verbosity: 0,
            ..GeomOptOptions::default()
        };

        println!("  n-decane, {} atoms", decane.atoms.len());
        let mut timings = Vec::new();
        for method in LocalOptimizer::ALL {
            let mut objective = DreidingObjective::new(&topology);
            let started = Instant::now();
            let result = local_optimize(
                method,
                coords.clone(),
                connectivity.clone(),
                &mut objective,
                options.clone(),
                &[],
            )
            .expect("decane relaxes");
            let elapsed = started.elapsed().as_secs_f64();
            println!(
                "  {:<22} {:>8.3} s  {:>4} cycles  converged {}  E = {:.6} Eh",
                method.label(),
                elapsed,
                result.cycles,
                result.converged,
                result.energy
            );
            timings.push(elapsed);
        }
        println!();
        let slowest = timings.iter().copied().fold(0.0_f64, f64::max);
        for (method, elapsed) in LocalOptimizer::ALL.iter().zip(&timings) {
            println!(
                "  {:<22} {:>6.0}x faster than the slowest",
                method.label(),
                slowest / elapsed.max(1e-9)
            );
        }
    }

    /// Cyclohexane, chair. 18 atoms, six C-C bonds -- every one of them in the ring.
    const CYCLOHEXANE: &str = "18\n\
        cyclohexane chair\n\
        C 1.4600 0.0000 0.2500\n\
        C 0.7300 1.2644 -0.2500\n\
        C -0.7300 1.2644 0.2500\n\
        C -1.4600 0.0000 -0.2500\n\
        C -0.7300 -1.2644 0.2500\n\
        C 0.7300 -1.2644 -0.2500\n\
        H 1.4600 0.0000 1.3400\n\
        H 2.4000 0.0000 -0.3000\n\
        H 0.7300 1.2644 -1.3400\n\
        H 1.2000 2.0785 0.3000\n\
        H -0.7300 1.2644 1.3400\n\
        H -1.2000 2.0785 -0.3000\n\
        H -1.4600 0.0000 -1.3400\n\
        H -2.4000 0.0000 0.3000\n\
        H -0.7300 -1.2644 1.3400\n\
        H -1.2000 -2.0785 -0.3000\n\
        H 0.7300 -1.2644 -1.3400\n\
        H 1.2000 -2.0785 0.3000\n";

    /// Cyclohexane's ring is a searchable degree of freedom, even though not one of its bonds
    /// rotates. This is the whole point of the pucker coordinates: a torsion-space search alone
    /// would call the molecule rigid.
    #[test]
    fn a_ring_is_a_searchable_degree_of_freedom() {
        let cyclohexane = molecule(CYCLOHEXANE);
        let topology = DreidingTopology::build(&cyclohexane).expect("cyclohexane types fine");
        // No rotatable bond: every C-C is in the ring.
        assert!(topology.rotatable_bonds().is_empty());

        let outcome = search(&cyclohexane, test_settings()).expect("the ring is searchable");
        // Counted as a ring, not a rotatable bond: none of its bonds turns.
        assert_eq!(outcome.rotors, 0);
        assert_eq!(outcome.rings, 1, "the ring itself is the degree of freedom");
        assert!(!outcome.conformers.is_empty());
        for conformer in &outcome.conformers {
            assert_eq!(
                conformer.positions_angstrom.len(),
                cyclohexane.atoms.len() * 3
            );
        }
    }

    /// A constrained relaxation must hold what it was told to hold. This is the whole point for
    /// transition states: the reacting part stays put while the rest of the molecule moves.
    ///
    /// Tested on one relaxation rather than a whole search, and deliberately. Constraints can only
    /// be expressed in redundant internal coordinates, so a constrained candidate costs what an
    /// internal-coordinate cycle costs -- roughly a hundred times a Cartesian one -- and a search
    /// made of hundreds of them has no business inside a test suite. What is being tested is that
    /// the constraint survives a relaxation, and one relaxation shows that.
    #[test]
    fn a_constrained_relaxation_holds_its_coordinates() {
        use crate::conformer::constraints::{measure_angle, measure_bond, ActiveSiteLevel};
        use crate::forcefield::dreiding::objective::DreidingObjective;
        use crate::optimizer::conformer_search::local::local_optimize;
        use crate::optimizer::geom_opt::GeomOptOptions;

        let butane = molecule(BUTANE);
        let topology = DreidingTopology::build(&butane).unwrap();
        let coords_bohr: Vec<f64> = butane
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
        let adjacency = {
            let mut adjacency = vec![Vec::new(); butane.atoms.len()];
            for &(i, j, _) in &butane.bonds {
                adjacency[i].push(j);
                adjacency[j].push(i);
            }
            adjacency
        };

        // Treat C1-C2-C3 as an active site: hold both lengths and the angle between them, well
        // away from the values the force field would relax them to.
        let mut stretched = coords_bohr.clone();
        // Pull C1 out along x so the held bond starts at a length DREIDING does not want.
        stretched[0] -= 0.6;
        let held = constraints::active_site(
            &stretched,
            0,
            1,
            2,
            ActiveSiteLevel::BondsAndAngle,
            &adjacency,
        );
        assert_eq!(held.len(), 3);
        let wanted_bond = measure_bond(&stretched, 0, 1);
        let wanted_angle = measure_angle(&stretched, 0, 1, 2);
        // The constraint has to be doing work, or the test proves nothing.
        assert!(
            (wanted_bond - measure_bond(&coords_bohr, 0, 1)).abs() > 0.3,
            "the test geometry was not actually distorted"
        );

        let targets: Vec<_> = held.iter().map(|c| c.target()).collect();
        let mut objective = DreidingObjective::new(&topology);
        let relaxed = local_optimize(
            crate::optimizer::conformer_search::LocalOptimizer::RedundantInternal,
            stretched,
            connectivity_model(&butane.atoms),
            &mut objective,
            GeomOptOptions {
                verbosity: 0,
                // Few cycles on purpose. The constrained optimiser satisfies its constraints
                // early and then relaxes the rest, so holding them is visible within a handful of
                // steps -- and internal coordinates are expensive enough that letting this run to
                // full convergence would put half a minute into the test suite.
                max_cycles: 20,
                ..GeomOptOptions::default()
            },
            &targets,
        )
        .expect("a constrained relaxation runs");

        let bond = measure_bond(&relaxed.x, 0, 1);
        let angle = measure_angle(&relaxed.x, 0, 1, 2);
        assert!(
            (bond - wanted_bond).abs() < 0.02,
            "the held bond drifted from {wanted_bond:.4} to {bond:.4} bohr"
        );
        assert!(
            (angle - wanted_angle).to_degrees().abs() < 2.0,
            "the held angle drifted from {:.1}° to {:.1}°",
            wanted_angle.to_degrees(),
            angle.to_degrees()
        );
    }

    /// The same claim through a whole search. Ignored: a constrained search is hundreds of
    /// internal-coordinate relaxations and takes minutes.
    ///
    /// ```text
    /// cargo test --release --bins a_constrained_search_holds -- --ignored --nocapture
    /// ```
    #[test]
    #[ignore = "a constrained search takes minutes"]
    fn a_constrained_search_holds_its_coordinates() {
        use crate::conformer::constraints::{measure_bond, ActiveSiteLevel};

        let butane = molecule(BUTANE);
        let coords_bohr: Vec<f64> = butane
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
        let adjacency = {
            let mut adjacency = vec![Vec::new(); butane.atoms.len()];
            for &(i, j, _) in &butane.bonds {
                adjacency[i].push(j);
                adjacency[j].push(i);
            }
            adjacency
        };
        let held =
            constraints::active_site(&coords_bohr, 0, 1, 2, ActiveSiteLevel::Bonds, &adjacency);
        let wanted = measure_bond(&coords_bohr, 0, 1);

        let outcome = search(
            &butane,
            ConformerSettings {
                constraints: held,
                ..test_settings()
            },
        )
        .expect("a constrained butane is searchable");
        assert_eq!(outcome.constrained, 2);
        println!(
            "  {} conformers in {:.1} s",
            outcome.conformers.len(),
            outcome.elapsed.as_secs_f64()
        );
        for conformer in &outcome.conformers {
            let bohr: Vec<f64> = conformer
                .positions_angstrom
                .iter()
                .map(|a| a / BOHR_TO_ANGSTROM)
                .collect();
            assert!((measure_bond(&bohr, 0, 1) - wanted).abs() < 0.02);
        }
    }

    /// A constraint that names the same coordinate twice must be refused with a message, not sent
    /// to a linear solve that reports a singular system.
    #[test]
    fn a_duplicated_constraint_is_refused_before_searching() {
        let butane = molecule(BUTANE);
        let one = constraints::Constraint {
            coordinate: crate::optimizer::constrained::ConstraintCoordinate::Bond([0, 1]),
            value: 2.9,
        };
        let error = search(
            &butane,
            ConformerSettings {
                constraints: vec![one.clone(), one],
                ..test_settings()
            },
        )
        .expect_err("a duplicate must be refused");
        match &error {
            ConformerError::Constraints(message) => {
                assert!(message.contains("fixed twice"), "{message}")
            }
            other => panic!("expected a constraint error, got {other:?}"),
        }
    }

    /// Cyclohexane's conformers: the chair lowest, twist-boats well above it. Getting the
    /// *ordering* right is the chemistry that matters -- a search that put a boat below the chair
    /// would be worthless.
    #[test]
    fn cyclohexane_puts_the_chair_below_the_boats() {
        use crate::optimizer::conformer_search::rings;

        let outcome = search(&molecule(CYCLOHEXANE), test_settings()).expect("searchable");
        assert!(!outcome.conformers.is_empty());

        let ring = [0usize, 1, 2, 3, 4, 5];
        let pucker_of = |conformer: &Conformer| {
            let bohr: Vec<f64> = conformer
                .positions_angstrom
                .iter()
                .map(|a| a / BOHR_TO_ANGSTROM)
                .collect();
            rings::measure(&bohr, &ring).expect("six-ring")
        };

        // The global minimum is a chair, at a pole of the puckering sphere.
        let lowest = pucker_of(&outcome.conformers[0]);
        assert!(
            lowest.theta_deg < 25.0 || lowest.theta_deg > 155.0,
            "the lowest conformer is not a chair: θ = {:.1}°",
            lowest.theta_deg
        );
        // And properly puckered rather than flattened: planar cyclohexane is a maximum.
        assert!(
            lowest.amplitude > 0.7,
            "the ring came back nearly flat, amplitude {:.2} bohr",
            lowest.amplitude
        );

        // Anything that is not a chair must cost real energy. The experimental chair-to-twist-boat
        // gap is about 5.5 kcal/mol; a generic force field overshoots, but the sign is not in
        // doubt.
        for (conformer, pucker) in outcome.conformers.iter().zip(outcome.conformers.iter().map(pucker_of)) {
            let is_chair = pucker.theta_deg < 25.0 || pucker.theta_deg > 155.0;
            if !is_chair {
                assert!(
                    conformer.relative_energy_kcal > 3.0,
                    "a boat at only {:.2} kcal/mol above the chair",
                    conformer.relative_energy_kcal
                );
            }
        }
    }

    /// Methylcyclohexane, equatorial chair. C7H14, 21 atoms.
    const METHYLCYCLOHEXANE: &str = "21\n\
        methylcyclohexane, equatorial chair\n\
        C 1.4600 0.0000 0.2500\n\
        C 0.7300 1.2644 -0.2500\n\
        C -0.7300 1.2644 0.2500\n\
        C -1.4600 0.0000 -0.2500\n\
        C -0.7300 -1.2644 0.2500\n\
        C 0.7300 -1.2644 -0.2500\n\
        H 1.4600 0.0000 1.3400\n\
        H 0.7300 1.2644 -1.3400\n\
        H 1.2000 2.0785 0.3000\n\
        H -0.7300 1.2644 1.3400\n\
        H -1.2000 2.0785 -0.3000\n\
        H -1.4600 0.0000 -1.3400\n\
        H -2.4000 0.0000 0.3000\n\
        H -0.7300 -1.2644 1.3400\n\
        H -1.2000 -2.0785 -0.3000\n\
        H 0.7300 -1.2644 -1.3400\n\
        H 1.2000 -2.0785 0.3000\n\
        C 2.7806 0.0000 -0.5227\n\
        H 3.0938 -1.0278 -0.7060\n\
        H 3.5433 0.5139 0.0623\n\
        H 2.6443 0.5139 -1.4742\n";

    /// The textbook ring case. Methylcyclohexane has two chairs -- methyl equatorial and methyl
    /// axial, about 1.8 kcal/mol apart experimentally -- and, higher up, twist-boats. Finding them
    /// is the proof that ring puckering works on real chemistry: no bond in the ring rotates, so a
    /// torsion-space search would find only the methyl spinning and call the molecule solved.
    #[test]
    fn methylcyclohexane_finds_both_chairs() {
        let mecy = molecule(METHYLCYCLOHEXANE);
        let outcome = search(&mecy, test_settings()).expect("methylcyclohexane is searchable");
        // The ring puckers; the methyl's own rotation is excluded as redundant.
        assert_eq!(outcome.rings, 1, "the ring should be sampled");

        let ring = [0usize, 1, 2, 3, 4, 5];
        let pucker_of = |conformer: &Conformer| {
            let bohr: Vec<f64> = conformer
                .positions_angstrom
                .iter()
                .map(|a| a / BOHR_TO_ANGSTROM)
                .collect();
            crate::optimizer::conformer_search::rings::measure(&bohr, &ring).expect("six-ring")
        };
        let is_chair = |p: &crate::optimizer::conformer_search::rings::Pucker| {
            p.theta_deg < 30.0 || p.theta_deg > 150.0
        };

        let puckers: Vec<_> = outcome.conformers.iter().map(pucker_of).collect();
        for (index, p) in puckers.iter().enumerate() {
            println!(
                "  conformer {}: dE {:>6.2} kcal/mol, θ {:>6.1}°, φ {:>6.1}°, Q {:.2} Å  {}",
                index + 1,
                outcome.conformers[index].relative_energy_kcal,
                p.theta_deg,
                p.phi_deg,
                p.amplitude,
                if is_chair(p) { "chair" } else { "boat/twist" }
            );
        }

        // The global minimum must be a chair. Getting that wrong would be a real error.
        assert!(
            is_chair(&puckers[0]),
            "the lowest conformer is not a chair: θ = {:.1}°",
            puckers[0].theta_deg
        );

        // And both chairs must be there -- that is the axial/equatorial pair.
        let chairs: Vec<usize> = puckers
            .iter()
            .enumerate()
            .filter(|(_, p)| is_chair(p))
            .map(|(index, _)| index)
            .collect();
        assert!(
            chairs.len() >= 2,
            "found {} chair(s) among {} conformers; axial and equatorial should both be there",
            chairs.len(),
            puckers.len()
        );

        // The two chairs differ by where the methyl sits, so a couple of kcal/mol, not tens.
        let gap = outcome.conformers[chairs[1]].relative_energy_kcal;
        assert!(
            gap > 0.0 && gap < 8.0,
            "the two chairs are {gap:.2} kcal/mol apart"
        );
    }

    /// Cyclopentane. C5H10, 15 atoms.
    const CYCLOPENTANE: &str = "15\n\
        cyclopentane\n\
        C 1.2500 0.0000 0.2000\n\
        C 0.3863 1.1888 0.0000\n\
        C -1.0113 0.7347 0.0000\n\
        C -1.0113 -0.7347 0.0000\n\
        C 0.3863 -1.1888 0.0000\n\
        H 1.8500 0.0000 1.1100\n\
        H 1.8500 0.0000 -0.7100\n\
        H 0.5717 1.7595 0.9100\n\
        H 0.5717 1.7595 -0.9100\n\
        H -1.4967 1.0874 0.9100\n\
        H -1.4967 1.0874 -0.9100\n\
        H -1.4967 -1.0874 0.9100\n\
        H -1.4967 -1.0874 -0.9100\n\
        H 0.5717 -1.7595 0.9100\n\
        H 0.5717 -1.7595 -0.9100\n";

    /// Five-rings pucker on a circle rather than a sphere -- the pseudorotation itinerary of
    /// envelopes and twists -- which is the coordinate sugars are described in. Checked separately
    /// from the six-ring because it is a different expansion, with one puckering coordinate
    /// instead of two.
    #[test]
    fn a_five_ring_is_searched_on_its_pseudorotation_itinerary() {
        use crate::optimizer::conformer_search::rings;

        let cyclopentane = molecule(CYCLOPENTANE);
        let outcome = search(&cyclopentane, test_settings()).expect("the ring is searchable");
        assert_eq!(outcome.rotors, 0);
        assert_eq!(outcome.rings, 1, "the ring is the degree of freedom");
        assert!(!outcome.conformers.is_empty());

        // Every conformer must be a genuinely puckered five-ring, not a flattened one: a planar
        // cyclopentane is the maximum, not a minimum.
        for conformer in &outcome.conformers {
            let bohr: Vec<f64> = conformer
                .positions_angstrom
                .iter()
                .map(|a| a / BOHR_TO_ANGSTROM)
                .collect();
            let pucker = rings::measure(&bohr, &[0, 1, 2, 3, 4]).expect("five-ring");
            assert!(
                pucker.amplitude > 0.3,
                "ring came back nearly planar, amplitude {:.2} bohr",
                pucker.amplitude
            );
        }
    }

    /// The species atmospheric and carbohydrate chemistry actually brings: radicals, peroxides,
    /// and rings with an ether or ester oxygen in them.
    ///
    /// Typing these correctly is what everything downstream rests on -- a peroxy oxygen typed as a
    /// carbonyl would be given a bond about 0.2 Å too short -- and the failures are silent, so the
    /// expectations are pinned here rather than left to a probe.
    #[test]
    fn the_species_of_atmospheric_and_sugar_chemistry_type_correctly() {
        struct Case {
            name: &'static str,
            xyz: &'static str,
            /// Types expected at particular atoms.
            expect: &'static [(usize, &'static str)],
            radical: bool,
        }

        let cases = [
            Case {
                name: "CH3OO* peroxy radical",
                xyz: "6\n\nC 0.000 0.000 0.000\nH 0.630 0.890 0.100\nH 0.630 -0.890 0.100\nH -0.500 0.000 -0.970\nO -1.030 0.000 1.010\nO -0.430 0.000 2.230\n",
                // Both oxygens sp3: a peroxy O typed as a carbonyl would want an O-O near 1.11 Å
                // against the real 1.32.
                expect: &[(0, "C_3"), (4, "O_3"), (5, "O_3")],
                radical: true,
            },
            Case {
                name: "CH3OOH hydroperoxide",
                xyz: "7\n\nC 0.000 0.000 0.000\nH 0.630 0.890 0.100\nH 0.630 -0.890 0.100\nH -0.500 0.000 -0.970\nO -1.030 0.000 1.010\nO -0.430 0.000 2.230\nH 0.070 0.760 2.430\n",
                expect: &[(4, "O_3"), (5, "O_3"), (6, "H_HB")],
                radical: false,
            },
            Case {
                name: "OH* hydroxyl radical",
                xyz: "2\n\nO 0.000 0.000 0.000\nH 0.970 0.000 0.000\n",
                expect: &[(0, "O_3")],
                radical: true,
            },
            Case {
                name: "(CH3)3CO* alkoxy radical",
                xyz: "14\n\
            tert-butoxy radical\n\
            C 0.0000 0.0000 0.0000\n\
            O 0.0000 0.0000 1.3800\n\
            C 1.4428 0.0000 -0.5095\n\
            H 1.4428 0.0000 -1.5995\n\
            H 1.9562 -0.8901 -0.1458\n\
            H 1.9562 0.8901 -0.1458\n\
            C -0.7206 1.2485 -0.5095\n\
            H -0.7205 1.2482 -1.5995\n\
            H -0.2061 2.1380 -0.1461\n\
            H -1.7483 1.2479 -0.1461\n\
            C -0.7206 -1.2485 -0.5095\n\
            H -0.7205 -1.2482 -1.5995\n\
            H -1.7483 -1.2479 -0.1461\n\
            H -0.2061 -2.1380 -0.1461\n",
                expect: &[(0, "C_3"), (1, "O_3")],
                radical: true,
            },
            Case {
                name: "1,2-dioxane, ring -O-O-",
                xyz: "14\n\
            1,2-dioxane (cyclic peroxide)\n\
            O 1.4600 0.0000 0.2500\n\
            O 0.7300 1.2644 -0.2500\n\
            C -0.7300 1.2644 0.2500\n\
            C -1.4600 0.0000 -0.2500\n\
            C -0.7300 -1.2644 0.2500\n\
            C 0.7300 -1.2644 -0.2500\n\
            H -0.7300 1.2644 1.3400\n\
            H -1.2000 2.0785 -0.3000\n\
            H -1.4600 0.0000 -1.3400\n\
            H -2.4000 0.0000 0.3000\n\
            H -0.7300 -1.2644 1.3400\n\
            H -1.2000 -2.0785 -0.3000\n\
            H 0.7300 -1.2644 -1.3400\n\
            H 1.2000 -2.0785 0.3000\n",
                // A peroxide O-O is a single bond even inside a ring: never aromatic.
                expect: &[(0, "O_3"), (1, "O_3"), (2, "C_3")],
                radical: false,
            },
            Case {
                name: "tetrahydropyran, the pyranose ring",
                xyz: "16\n\
            thp\n\
            O 1.4400 0.0000 0.2400\n\
            C 0.7200 1.2471 -0.2400\n\
            C -0.7200 1.2471 0.2400\n\
            C -1.4400 0.0000 -0.2400\n\
            C -0.7200 -1.2471 0.2400\n\
            C 0.7200 -1.2471 -0.2400\n\
            H 0.7200 1.2471 -1.3300\n\
            H 1.1900 2.0611 0.3100\n\
            H -0.7200 1.2471 1.3300\n\
            H -1.1900 2.0611 -0.3100\n\
            H -1.4400 0.0000 -1.3300\n\
            H -2.3800 0.0000 0.3100\n\
            H -0.7200 -1.2471 1.3300\n\
            H -1.1900 -2.0611 -0.3100\n\
            H 0.7200 -1.2471 -1.3300\n\
            H 1.1900 -2.0611 0.3100\n",
                expect: &[(0, "O_3"), (1, "C_3")],
                radical: false,
            },
            Case {
                name: "tetrahydrofuran, the furanose ring",
                xyz: "13\n\
            thf\n\
            O 1.2300 0.0000 0.2000\n\
            C 0.3801 1.1698 -0.2000\n\
            C -0.9951 0.7230 0.2000\n\
            C -0.9951 -0.7230 -0.2000\n\
            C 0.3801 -1.1698 0.2000\n\
            H 0.3801 1.1698 -1.2900\n\
            H 0.6706 2.0638 0.3500\n\
            H -0.9951 0.7230 1.2900\n\
            H -1.7556 1.2755 -0.3500\n\
            H -0.9951 -0.7230 -1.2900\n\
            H -1.7556 -1.2755 0.3500\n\
            H 0.3801 -1.1698 1.2900\n\
            H 0.6706 -2.0638 -0.3500\n",
                expect: &[(0, "O_3"), (1, "C_3")],
                radical: false,
            },
        ];

        for case in cases {
            let mol = molecule(case.xyz);
            let topology = match DreidingTopology::build(&mol) {
                Ok(t) => t,
                Err(e) => panic!("{}: {e}", case.name),
            };
            let types = topology.atom_types();
            for (index, wanted) in case.expect {
                assert_eq!(
                    &types[*index], wanted,
                    "{}: atom {} typed {} but should be {wanted}. All: {types:?}",
                    case.name,
                    index + 1,
                    types[*index]
                );
            }
            assert_eq!(
                topology.open_shell().odd_electron_count,
                case.radical,
                "{}: electron parity",
                case.name
            );
            // A radical must produce a warning; a closed-shell molecule must not be told it has
            // an odd electron count.
            if case.radical {
                let warning = topology
                    .radical_warning()
                    .unwrap_or_else(|| panic!("{}: no radical warning", case.name));
                assert!(warning.contains("radical"), "{}: {warning}", case.name);
            }
        }
    }

    /// A lactone's ring carbonyl must come out as a carbonyl -- sp2 carbon and sp2 oxygen -- while
    /// the ester oxygen in the ring stays sp3. Getting this wrong would put the ring in the wrong
    /// shape entirely.
    #[test]
    fn a_ring_ester_types_as_a_carbonyl_and_an_ether() {
        let mol = molecule("15\n\
            delta-valerolactone\n\
            O 1.4400 0.0000 0.1000\n\
            C 0.7200 1.2471 -0.1000\n\
            C -0.7200 1.2471 0.1000\n\
            C -1.4400 0.0000 -0.1000\n\
            C -0.7200 -1.2471 0.1000\n\
            C 0.7200 -1.2471 -0.1000\n\
            O 1.3250 2.2950 -0.1000\n\
            H -0.7200 1.2471 1.1900\n\
            H -1.1900 2.0611 -0.4500\n\
            H -1.4400 0.0000 -1.1900\n\
            H -2.3800 0.0000 0.4500\n\
            H -0.7200 -1.2471 1.1900\n\
            H -1.1900 -2.0611 -0.4500\n\
            H 0.7200 -1.2471 -1.1900\n\
            H 1.1900 -2.0611 0.4500\n");
        let topology = DreidingTopology::build(&mol).expect("a lactone types");
        let types = topology.atom_types();
        // Atom 1 is the ester carbon, atom 7 its carbonyl oxygen, atom 0 the ring ether oxygen.
        assert!(
            types[1] == "C_2" || types[1] == "C_R",
            "the ester carbon typed {} -- expected sp2 or resonant. All: {types:?}",
            types[1]
        );
        assert!(
            types[6] == "O_2" || types[6] == "O_R",
            "the carbonyl oxygen typed {} -- expected sp2 or resonant. All: {types:?}",
            types[6]
        );
    }

    /// Water really is rigid, so it must still get the plain message rather than the ring one.
    #[test]
    fn a_genuinely_rigid_molecule_still_says_so() {
        let water = molecule("3\n\nO 0.000 0.000 0.000\nH 0.960 0.000 0.000\nH -0.240 0.929 0.000\n");
        let error = search(&water, test_settings()).expect_err("water has one shape");
        assert!(matches!(error, ConformerError::Rigid), "{error:?}");
        assert!(error.to_string().contains("only one conformer"));
    }

    #[test]
    fn populations_sum_to_one_and_favour_the_lowest() {
        let fractions = populations(&[0.0, 1.0, 2.0]);
        let total: f64 = fractions.iter().sum();
        assert!((total - 1.0).abs() < 1.0e-12);
        assert!(fractions[0] > fractions[1] && fractions[1] > fractions[2]);
    }

    /// Degenerate conformers must share the population equally, and 1.36 kcal/mol is about a
    /// factor of ten at room temperature -- a useful sanity check on the constant.
    #[test]
    fn populations_match_the_boltzmann_factor() {
        let equal = populations(&[0.0, 0.0]);
        assert!((equal[0] - 0.5).abs() < 1.0e-12);

        let decade = populations(&[0.0, 1.364]);
        assert!(
            (decade[0] / decade[1] - 10.0).abs() < 0.1,
            "ratio is {}",
            decade[0] / decade[1]
        );
    }

    /// A generation must be wide enough to keep several cores busy, or `ncore` buys nothing in
    /// the main loop -- which is exactly the limitation this replaced.
    #[test]
    fn a_generation_is_wide_enough_to_parallelise() {
        assert!(
            OFFSPRING_PER_GENERATION >= 4,
            "a generation of {OFFSPRING_PER_GENERATION} cannot fill a machine"
        );
        // And there must still be enough generations for selection to mean anything.
        for level in [Thoroughness::Quick, Thoroughness::Normal, Thoroughness::Thorough] {
            let (_, generations) = level.budget(3);
            assert!(
                generations >= 10,
                "{} has only {generations} generations",
                level.label()
            );
        }
    }

    /// The budget must grow with the rotor count and with thoroughness, and stay bounded so a
    /// large molecule cannot ask for an unbounded search.
    #[test]
    fn the_search_budget_scales_with_rotors_and_thoroughness() {
        let (quick_pop, quick_gen) = Thoroughness::Quick.budget(3);
        let (normal_pop, normal_gen) = Thoroughness::Normal.budget(3);
        assert!(quick_pop < normal_pop && quick_gen < normal_gen);

        // More rotors, more sampling.
        assert!(Thoroughness::Normal.budget(8).0 > Thoroughness::Normal.budget(2).0);
        // But bounded, and never zero even for a nominally rotor-free call.
        assert!(Thoroughness::Thorough.budget(1000).0 <= 150);
        assert!(Thoroughness::Quick.budget(0).0 >= 12);
    }
}
