use crate::optimizer::GeomOptOptions;

const EV_TO_HARTREE: f64 = 1.0 / 27.211_386_245_988;
const KCAL_MOL_TO_HARTREE: f64 = 1.0 / 627.509_474_063_1;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SelectionMode {
    #[default]
    Roulette,
    ReverseRoulette,
    Random,
}

impl SelectionMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::Roulette => "roulette",
            Self::ReverseRoulette => "reverse roulette",
            Self::Random => "random",
        }
    }
}

#[derive(Debug, Clone)]
pub struct ConformerSearchOptions {
    /// Maximum number of independent local optimizations evaluated at once.
    /// `None` uses the configured Rayon worker count (normally `nproc`).
    pub parallel_workers: Option<usize>,
    pub population_size: usize,
    pub max_generations: usize,
    pub stagnation_generations: usize,
    pub energy_convergence_hartree: f64,
    pub target_energy_hartree: Option<f64>,
    pub energy_variance_hartree: f64,
    pub selection: SelectionMode,
    pub fitness_sum_limit: f64,
    pub crossover_probability: f64,
    pub crossover_attempts: usize,
    pub torsion_mutation_probability: f64,
    pub cis_trans_mutation_probability: f64,
    pub max_torsion_mutations: usize,
    pub max_cis_trans_mutations: usize,
    pub mutation_attempts: usize,
    pub max_initial_attempts: usize,
    pub minimum_nonbonded_distance_angstrom: f64,
    pub maximum_bonded_distance_angstrom: f64,
    pub rmsd_cutoff_angstrom: f64,
    pub chiral: bool,
    pub unique_initial_population: bool,
    pub exclude_methyl_rotors: bool,
    /// Explicit rotatable dihedrals. Library indices are zero-based; the CLI
    /// translates its documented one-based input before calling the library.
    pub rotatable_dihedrals: Vec<[usize; 4]>,
    /// Explicit two-state cis/trans dihedrals, using the same index convention.
    pub cis_trans_dihedrals: Vec<[usize; 4]>,
    pub random_seed: Option<u64>,
    pub output_energy_window_hartree: f64,
    pub local_optimization: GeomOptOptions,
    pub verbosity: u8,
}

impl Default for ConformerSearchOptions {
    fn default() -> Self {
        Self {
            parallel_workers: None,
            population_size: 5,
            max_generations: 10,
            stagnation_generations: 10,
            energy_convergence_hartree: 0.001 * EV_TO_HARTREE,
            target_energy_hartree: None,
            energy_variance_hartree: 0.001 * EV_TO_HARTREE,
            selection: SelectionMode::Roulette,
            fitness_sum_limit: 1.2,
            crossover_probability: 0.95,
            crossover_attempts: 20,
            torsion_mutation_probability: 0.5,
            cis_trans_mutation_probability: 0.5,
            max_torsion_mutations: 2,
            max_cis_trans_mutations: 1,
            mutation_attempts: 100,
            max_initial_attempts: 1_000,
            minimum_nonbonded_distance_angstrom: 1.3,
            maximum_bonded_distance_angstrom: 2.15,
            rmsd_cutoff_angstrom: 0.2,
            chiral: false,
            unique_initial_population: true,
            exclude_methyl_rotors: true,
            rotatable_dihedrals: Vec::new(),
            cis_trans_dihedrals: Vec::new(),
            random_seed: None,
            output_energy_window_hartree: 20.0 * KCAL_MOL_TO_HARTREE,
            local_optimization: GeomOptOptions {
                verbosity: 0,
                ..GeomOptOptions::default()
            },
            verbosity: 1,
        }
    }
}

impl ConformerSearchOptions {
    /// Effective conformer-level worker count after applying the optional cap.
    pub fn effective_parallel_workers(&self) -> usize {
        let available = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1).max(1);
        self.parallel_workers
            .unwrap_or(available)
            .min(available)
            .max(1)
    }
}
