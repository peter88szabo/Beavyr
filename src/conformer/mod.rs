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
//! a process launch: xTB's GFN-FF is fast once running, but at ~100 ms of startup per call the
//! same search would take hours to days rather than seconds.
//!
//! That is why the search runs on DREIDING, and why xTB is offered as a *refinement* of the
//! conformers that survive instead. Re-optimising the best few dozen is a few hundred xTB calls --
//! a minute or two -- and it buys a markedly better energy ordering than a generic force field
//! gives, because relative conformer energies are where a generic force field is weakest.

pub mod ui;

use std::time::{Duration, Instant};

use bevy::prelude::*;
use bevy::tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task};

use crate::forcefield::dreiding::objective::{
    DreidingObjective, BOHR_TO_ANGSTROM, HARTREE_TO_KCAL,
};
use crate::forcefield::dreiding::{BuildError, DreidingTopology};
use crate::molecule::{covalent_radius_angstrom, Molecule};
use crate::optimizer::conformer_search::{
    discover_torsions, genetic_conformer_search, ConformerSearchOptions, ConformerSearchResult,
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
    fn budget(self, rotors: usize) -> (usize, usize) {
        let rotors = rotors.max(1);
        match self {
            Self::Quick => ((6 * rotors).clamp(12, 40), 30),
            Self::Normal => ((10 * rotors).clamp(20, 80), 80),
            Self::Thorough => ((16 * rotors).clamp(32, 150), 200),
        }
    }
}

/// What the user chooses. Deliberately short: rotors, population and generations are all worked
/// out from the structure rather than asked for.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConformerSettings {
    pub thoroughness: Thoroughness,
    /// Re-optimise the surviving conformers with xTB's GFN-FF and re-rank them.
    ///
    /// **Not yet implemented**, and disabled in the panel accordingly. The blocker is not the
    /// cost -- refining a few dozen survivors is a few hundred xTB calls, minutes, and well worth
    /// it since relative conformer energies are exactly where a generic force field is weakest.
    /// It is that `qchem_interfaces::xtbrun::call_xtb` writes `xtb_geom_file_for_abinitioMD.xyz`
    /// and reads `gradient` from the *current working directory*, so calling it in a loop would
    /// litter the user's directory and could not be made concurrent. Driving it needs a
    /// scratch-directory-aware variant first; `xtb_scratch_dir()` already exists for that.
    pub refine_with_gfnff: bool,
    /// Fix the random seed, so a search can be repeated exactly.
    pub reproducible: bool,
    /// Keep conformers within this window of the lowest, in kcal/mol.
    pub energy_window_kcal: f64,
}

impl Default for ConformerSettings {
    fn default() -> Self {
        Self {
            thoroughness: Thoroughness::default(),
            refine_with_gfnff: false,
            reproducible: true,
            // Wide enough to keep everything thermally accessible, and then some: at room
            // temperature 10 kcal/mol is already a population of ~1e-7.
            energy_window_kcal: 10.0,
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
    pub generations: usize,
    pub structures_tried: usize,
    pub local_optimizations: usize,
    pub termination: String,
    pub elapsed: Duration,
    /// The generation that produced the lowest-energy conformer.
    ///
    /// Worth showing: if the best structure only turned up in the last generation, the search was
    /// still improving when it stopped and more effort would probably find something better.
    pub best_found_in_generation: usize,
    /// Atom symbols, so a conformer can be written out or loaded back.
    pub atoms: Vec<String>,
}

/// Room temperature, for the population estimate.
const ROOM_TEMPERATURE_K: f64 = 298.15;

/// kcal/mol per K, the gas constant.
const R_KCAL: f64 = 0.001_987_204_1;

/// Boltzmann populations from relative energies at room temperature.
fn populations(relative_kcal: &[f64]) -> Vec<f64> {
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

    // The rotor count must come from the search's own torsion perception, not from the force
    // field's `rotatable_bonds`. The two disagree on purpose: DREIDING counts every acyclic single
    // bond, while the search excludes methyl rotors, whose rotation returns the same shape. For
    // butane that is three bonds against one -- and budgeting from the larger number asks for a
    // population bigger than the molecule has distinct conformers.
    let rotors = match discover_torsions(&coordinates_bohr, &atoms, &connectivity, &[], &[], true) {
        Ok((_bonds, torsions)) => torsions.len(),
        // The only reason perception refuses is that nothing can rotate.
        Err(_) => return Err(ConformerError::Rigid),
    };
    if rotors == 0 {
        return Err(ConformerError::Rigid);
    }

    let (requested_population, max_generations) = settings.thoroughness.budget(rotors);
    // Each rotor contributes roughly three staggered minima, so the whole torsion space is about
    // 3^rotors. Asking for a population larger than that cannot be satisfied: the search needs
    // distinct starting minima and there are only so many to find. Butane has one rotor and, with
    // its two gauche forms treated as equivalent, just two conformers.
    let space = 3usize.saturating_pow(rotors.min(10) as u32);
    let population_size = requested_population.min(space).max(2);

    // A small molecule can still have fewer distinct minima than even that estimate -- symmetry
    // makes conformers coincide. Rather than report the search's "could not fill the population"
    // as a failure, come down to what the molecule actually has. Halving converges in a few
    // steps and each attempt on a small molecule is milliseconds.
    let mut population = population_size;
    let mut last_error;
    loop {
        let options = ConformerSearchOptions {
            population_size: population,
            max_generations,
            random_seed: settings.reproducible.then_some(FIXED_SEED),
            output_energy_window_hartree: settings.energy_window_kcal / HARTREE_TO_KCAL,
            // Methyl rotation does not make a distinct conformer, only a redundant copy.
            exclude_methyl_rotors: true,
            verbosity: 0,
            ..ConformerSearchOptions::default()
        };
        let mut objective = DreidingObjective::new(&topology);
        match genetic_conformer_search(
            coordinates_bohr.clone(),
            &atoms,
            connectivity.clone(),
            &mut objective,
            options,
        ) {
            Ok(result) => return Ok(assemble(result, atoms, started.elapsed())),
            Err(error) => {
                last_error = error.to_string();
                if population <= 2 {
                    break;
                }
                population = (population / 2).max(2);
            }
        }
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
        rotors: result.torsions.len(),
        generations: result.generations,
        structures_tried: result.statistics.attempted_structures,
        local_optimizations: result.statistics.local_optimizations,
        termination: format!("{:?}", result.termination),
        elapsed,
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
    Search(String),
}

impl std::fmt::Display for ConformerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ForceField(e) => write!(f, "{e}"),
            Self::Rigid => write!(
                f,
                "this molecule has no rotatable bonds, so it has only one conformer. \
                 Rings and double bonds do not rotate, and a methyl group's rotation \
                 gives the same shape back."
            ),
            Self::Search(message) => write!(f, "the conformer search failed: {message}"),
        }
    }
}

/// The running search, and the last result. A Bevy resource so the panel can poll it.
#[derive(Resource, Default)]
pub struct ConformerRun {
    task: Option<Task<Result<ConformerOutcome, ConformerError>>>,
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
        let settings = self.settings;
        self.task = Some(
            AsyncComputeTaskPool::get().spawn(async move { search(&snapshot, settings) }),
        );
    }
}

/// Collects a finished search. Registered by [`crate::conformer::ui`].
pub fn poll_conformer_search(mut run: ResMut<ConformerRun>) {
    let Some(task) = run.task.as_mut() else {
        return;
    };
    let Some(result) = block_on(future::poll_once(task)) else {
        return;
    };
    run.task = None;
    match result {
        Ok(outcome) => run.outcome = Some(outcome),
        Err(error) => run.error = Some(error.to_string()),
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
        let error = search(&water, ConformerSettings::default())
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
        match search(&lif, ConformerSettings::default()) {
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
                ..ConformerSettings::default()
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
            ..ConformerSettings::default()
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
            ..ConformerSettings::default()
        };
        let first = search(&butane, settings).expect("searchable");
        let second = search(&butane, settings).expect("searchable");

        assert_eq!(first.conformers.len(), second.conformers.len());
        for (a, b) in first.conformers.iter().zip(&second.conformers) {
            assert!((a.relative_energy_kcal - b.relative_energy_kcal).abs() < 1.0e-9);
        }
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
