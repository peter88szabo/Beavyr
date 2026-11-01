use super::{analysis::AnalysisResult, input::HessianSource};
use bevy::{prelude::Resource, tasks::Task};
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum AnalysisKind {
    FranckCondon,
    Duschinsky,
}
#[derive(Clone)]
pub struct Settings {
    pub kind: AnalysisKind,
    pub ground: HessianSource,
    pub excited: HessianSource,
    pub gradient: String,
    pub cutoff: f64,
    pub make_band: bool,
    pub energy_ev: f64,
    pub resolution: f64,
    pub quanta: usize,
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            kind: AnalysisKind::FranckCondon,
            ground: Default::default(),
            excited: Default::default(),
            gradient: String::new(),
            cutoff: 0.,
            make_band: true,
            energy_ev: 3.,
            resolution: 20.,
            quanta: 12,
        }
    }
}
#[derive(Resource, Default)]
pub struct VibronicState {
    pub open: bool,
    pub settings: Settings,
    pub task: Option<Task<Result<AnalysisResult, String>>>,
    pub result: Option<AnalysisResult>,
    pub error: Option<String>,
    pub plot: usize,
    pub heatmap: Option<bevy_egui::egui::TextureHandle>,
}
