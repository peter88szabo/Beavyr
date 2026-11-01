use super::{input::Fragment, output::SamplingOutput};
use bevy::prelude::Resource;
use bevy::tasks::Task;
use std::sync::{
    atomic::{AtomicBool, AtomicUsize},
    Arc,
};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Workflow {
    Unimolecular,
    Bimolecular,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
    Hessian,
    Diatomic,
    Atom,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
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
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
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
#[derive(Clone, Debug)]
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
#[derive(Clone)]
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
        }
    }
}
#[derive(Clone)]
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
        }
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
        }
    }
}
