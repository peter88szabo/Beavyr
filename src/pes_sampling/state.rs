use super::{input::Fragment, output::SamplingOutput};
use bevy::prelude::Resource;
use bevy::tasks::Task;
use std::sync::{
    atomic::{AtomicBool, AtomicUsize},
    Arc,
};
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Workflow {
    Unimolecular,
    Bimolecular,
    ReactionPath,
}
/// Source of the uniform numbers behind every random choice.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum RandomNumbers {
    /// Scrambled Halton low-discrepancy sequence: more even coverage.
    #[default]
    QuasiRandom,
    PseudoRandom,
}
impl RandomNumbers {
    pub fn label(self) -> &'static str {
        match self {
            Self::QuasiRandom => "Quasi-random (Halton)",
            Self::PseudoRandom => "Pseudo-random",
        }
    }
    pub fn comment(self) -> &'static str {
        match self {
            Self::QuasiRandom => "random=scrambled_halton",
            Self::PseudoRandom => "random=xoshiro256",
        }
    }
}
/// Space in which normal-mode displacements are applied.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Coordinates {
    /// Straight-line displacement along Cartesian normal modes.
    #[default]
    Cartesian,
    /// The same displacement mapped to bonds/angles/dihedrals (curvilinear).
    Internal,
}
impl Coordinates {
    pub fn label(self) -> &'static str {
        match self {
            Self::Cartesian => "Cartesian (rectilinear)",
            Self::Internal => "Internal (curvilinear)",
        }
    }
    pub fn comment(self) -> &'static str {
        match self {
            Self::Cartesian => "coordinates=cartesian",
            Self::Internal => "coordinates=internal_curvilinear",
        }
    }
}
/// Display unit for energies; values are stored in kcal/mol.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum EnergyUnit {
    #[default]
    KcalMol,
    KjMol,
    Ev,
    Wavenumber,
}
impl EnergyUnit {
    pub const ALL: [Self; 4] = [Self::KcalMol, Self::KjMol, Self::Ev, Self::Wavenumber];
    pub fn label(self) -> &'static str {
        match self {
            Self::KcalMol => "kcal/mol",
            Self::KjMol => "kJ/mol",
            Self::Ev => "eV",
            Self::Wavenumber => "cm⁻¹",
        }
    }
    /// Units per kcal/mol.
    pub fn per_kcal(self) -> f64 {
        match self {
            Self::KcalMol => 1.,
            Self::KjMol => 4.184,
            // 1 eV = 23.060548 kcal/mol; 1 kcal/mol = 349.7551 cm⁻¹.
            Self::Ev => 1. / 23.060548,
            Self::Wavenumber => 349.7551,
        }
    }
    /// Hartree → this unit (the sampler converts kcal/mol with 627.51).
    pub fn from_hartree(self, e: f64) -> f64 {
        e * 627.51 * self.per_kcal()
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Source {
    Hessian,
    Diatomic,
    Atom,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Population {
    Zpe,
    Thermal,
    Quantum,
    Energy,
}
impl Population {
    pub fn label(self) -> &'static str {
        match self {
            Self::Zpe => "ZPE",
            Self::Thermal => "Maximum temperature",
            Self::Quantum => "Maximum quantum number",
            Self::Energy => "Maximum energy per mode",
        }
    }
}
/// How a mode's displacement is drawn within its classical turning points ±A.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Displacement {
    /// Q = A·(2u−1): even coverage of the well, including its bottom.
    #[default]
    Uniform,
    /// Q = A·cos(2πu): time-weighted (arcsine) distribution, as in QCT/Smite.
    ClassicalPhase,
}
impl Displacement {
    pub fn label(self) -> &'static str {
        match self {
            Self::Uniform => "Uniform in amplitude",
            Self::ClassicalPhase => "Classical phase (cos)",
        }
    }
    pub fn comment(self) -> &'static str {
        match self {
            Self::Uniform => "displacement=uniform_within_turning_points",
            Self::ClassicalPhase => "phases=uniform_0_2pi",
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum TemperatureRule {
    BoltzmannTail,
    ZpePlusKt,
    SmiteThermal,
}
impl TemperatureRule {
    pub fn label(self) -> &'static str {
        match self {
            Self::BoltzmannTail => "Boltzmann-tail cutoff",
            Self::ZpePlusKt => "ZPE + kB × Tmax",
            Self::SmiteThermal => "Smite thermal (no hard cutoff)",
        }
    }
}
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct Excitation {
    pub population: Population,
    pub temperature: f64,
    pub quantum: u32,
    pub energy_kcal: f64,
    pub temperature_rule: TemperatureRule,
    pub retained_population: f64,
}
impl Default for Excitation {
    fn default() -> Self {
        Self {
            population: Population::Zpe,
            temperature: 300.,
            quantum: 0,
            energy_kcal: 1.,
            temperature_rule: TemperatureRule::BoltzmannTail,
            retained_population: 99.9,
        }
    }
}
#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct Partner {
    pub source: Source,
    pub fragment: Option<Fragment>,
    pub geometry: Option<(Vec<String>, Vec<f64>)>,
    pub symbols: [String; 2],
    pub equilibrium: f64,
    pub frequency: f64,
    pub morse: bool,
    pub dissociation_ev: f64,
    pub excitation: Excitation,
    pub rigid: bool,
    pub overrides: Vec<Excitation>,
    pub use_overrides: bool,
    pub ts_mode: super::ts_mode::Settings,
    pub rotors: Vec<super::rotors::Rotor>,
    /// File name of the equilibrium XYZ waiting for a coordinate-free Hessian.
    pub geometry_name: Option<String>,
    /// Where the loaded fragment's coordinates came from (shown to the user).
    pub coordinates_source: Option<String>,
}
impl Default for Partner {
    fn default() -> Self {
        Self {
            source: Source::Hessian,
            fragment: None,
            geometry: None,
            symbols: ["H".into(), "H".into()],
            equilibrium: 0.74,
            frequency: 4401.,
            morse: false,
            dissociation_ev: 4.75,
            excitation: Default::default(),
            rigid: false,
            overrides: vec![],
            use_overrides: false,
            ts_mode: Default::default(),
            rotors: vec![],
            geometry_name: None,
            coordinates_source: None,
        }
    }
}
#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct Settings {
    pub workflow: Workflow,
    pub partners: [Partner; 2],
    pub points: usize,
    pub distance_start: f64,
    pub distance_end: f64,
    pub distances: usize,
    pub samples_per_distance: usize,
    pub repeatable: bool,
    pub seed: u64,
    pub overlap: super::overlap::OverlapFilter,
    pub displacement: Displacement,
    pub random: RandomNumbers,
    pub coordinates: Coordinates,
    pub diversity: super::diverse::Diversity,
    pub path: super::path::PathSettings,
    pub energy_unit: EnergyUnit,
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            workflow: Workflow::Unimolecular,
            partners: [Partner::default(), Partner::default()],
            points: 100,
            distance_start: 3.,
            distance_end: 6.,
            distances: 10,
            samples_per_distance: 100,
            repeatable: false,
            seed: 42,
            overlap: Default::default(),
            displacement: Default::default(),
            random: Default::default(),
            coordinates: Default::default(),
            diversity: Default::default(),
            path: Default::default(),
            energy_unit: Default::default(),
        }
    }
}
/// One-click combinations of the sampling options.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Preset {
    MlTrainingSet,
    QctInitialConditions,
}
impl Preset {
    pub const ALL: [Self; 2] = [Self::MlTrainingSet, Self::QctInitialConditions];
    pub fn label(self) -> &'static str {
        match self {
            Self::MlTrainingSet => "ML training set",
            Self::QctInitialConditions => "QCT initial conditions",
        }
    }
    pub fn description(self) -> &'static str {
        match self {
            Self::MlTrainingSet => "Even coverage of the wells for fitting a PES: uniform displacement, quasi-random numbers, diverse subset (3×), Cartesian displacement, overlap filter on; every molecule (and the path) thermal with the Boltzmann-tail cutoff at its current temperature. Hessians, counts and temperatures are kept.",
            Self::QctInitialConditions => "Smite/Pebble QCT convention: classical phases (cos), pseudo-random numbers, no subset selection, Cartesian displacement, ZPE in every mode, per-mode limits off. Hessians and counts are kept.",
        }
    }
}
impl Settings {
    pub fn apply_preset(&mut self, preset: Preset) {
        let ml = preset == Preset::MlTrainingSet;
        self.displacement = if ml { Displacement::Uniform } else { Displacement::ClassicalPhase };
        self.random = if ml { RandomNumbers::QuasiRandom } else { RandomNumbers::PseudoRandom };
        self.coordinates = Coordinates::Cartesian;
        self.diversity.enabled = ml;
        if ml {
            self.diversity.oversampling = 3.;
            self.overlap.enabled = true;
        }
        let set = |e: &mut Excitation| {
            if ml {
                e.population = Population::Thermal;
                e.temperature_rule = TemperatureRule::BoltzmannTail;
            } else {
                e.population = Population::Zpe;
            }
        };
        for p in &mut self.partners {
            set(&mut p.excitation);
            p.use_overrides = false;
        }
        set(&mut self.path.excitation);
    }
}
#[derive(Resource)]
pub struct PesSampling {
    pub open: bool,
    pub settings: Settings,
    pub task: Option<Task<Result<SamplingOutput, String>>>,
    pub output: Option<SamplingOutput>,
    pub cancel: Arc<AtomicBool>,
    pub progress: Arc<AtomicUsize>,
    pub total: usize,
    pub message: Option<String>,
    pub is_error: bool,
    pub load_requested: bool,
    /// Mode animation in the viewer, started from the mode table.
    pub animation: crate::vibronic::mode_animation::Animation,
    pub histogram: super::analysis::HistogramView,
}
impl Default for PesSampling {
    fn default() -> Self {
        Self {
            open: false,
            settings: Default::default(),
            task: None,
            output: None,
            cancel: Arc::new(AtomicBool::new(false)),
            progress: Arc::new(AtomicUsize::new(0)),
            total: 0,
            message: None,
            is_error: false,
            load_requested: false,
            animation: Default::default(),
            histogram: Default::default(),
        }
    }
}
