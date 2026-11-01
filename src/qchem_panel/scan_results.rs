//! Scan data, exports, and an independent results window.
use super::run::QcRunTask;
use crate::qchem_interfaces::job::scan::{Coordinate, Options};
use crate::ui_style::ResponseExt;
use bevy_egui::egui;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ScanPoint {
    pub target: f64,
    pub actual: f64,
    pub energy: f64,
    pub converged: Option<bool>,
    pub steps: usize,
    pub coordinates: Vec<f64>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ScanReport {
    pub program: String,
    pub method: String,
    pub driver: String,
    pub options: Options,
    pub atoms: Vec<String>,
    pub points: Vec<ScanPoint>,
    pub complete: bool,
    pub message: String,
}
#[derive(Clone, Copy, Default, PartialEq, Eq)]
enum EnergyUnit {
    #[default]
    Kcal,
    Kj,
}
impl EnergyUnit {
    fn label(self) -> &'static str {
        match self {
            Self::Kcal => "kcal/mol",
            Self::Kj => "kJ/mol",
        }
    }
    fn factor(self) -> f64 {
        match self {
            Self::Kcal => 627.5094740631,
            Self::Kj => 2625.499639479163,
        }
    }
}
#[derive(Clone, Copy, Default, PartialEq, Eq)]
enum Reference {
    #[default]
    First,
    Last,
}
#[derive(Clone, Copy, Default)]
struct ScanView {
    unit: EnergyUnit,
    reference: Reference,
    reversed: bool,
}
impl ScanView {
    fn indices(self, report: &ScanReport) -> Vec<usize> {
        let mut indices: Vec<_> = (0..report.points.len()).collect();
        if self.reversed {
            indices.reverse();
        }
        indices
    }
    fn reference_index(self, report: &ScanReport) -> usize {
        if self.reversed ^ (self.reference == Reference::Last) {
            report.points.len().saturating_sub(1)
        } else {
            0
        }
    }
    fn energy(self, report: &ScanReport, i: usize) -> f64 {
        (report.points[i].energy - report.points[self.reference_index(report)].energy)
            * self.unit.factor()
    }
    fn controls(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.selectable_value(&mut self.unit, EnergyUnit::Kcal, "kcal/mol");
            ui.selectable_value(&mut self.unit, EnergyUnit::Kj, "kJ/mol");
            ui.add_space(12.);
            ui.label("Reference:");
            ui.selectable_value(&mut self.reference, Reference::First, "First point");
            ui.selectable_value(&mut self.reference, Reference::Last, "Last point");
            if ui
                .selectable_label(self.reversed, "Reverse order")
                .clicked_once()
            {
                self.reversed = !self.reversed;
            }
        });
    }
}
impl ScanReport {
    pub fn xyz(&self) -> String {
        self.xyz_ordered(ScanView::default())
    }
    fn xyz_ordered(&self, view: ScanView) -> String {
        let mut text = String::new();
        for i in view.indices(self) {
            let p = &self.points[i];
            text.push_str(&format!(
                "{}\ncycle {} E = {:.12} Eh target={:.8} actual={:.8} converged={:?}\n",
                self.atoms.len(),
                i + 1,
                p.energy,
                p.target,
                p.actual,
                p.converged
            ));
            text.push_str(&crate::qchem_interfaces::orca_run::xyz_body(
                &self.atoms,
                &p.coordinates,
            ));
        }
        text
    }
    pub fn csv(&self) -> String {
        self.csv_displayed(ScanView::default())
    }
    fn csv_displayed(&self, view: ScanView) -> String {
        let mut text=format!("point,target_{0},actual_{0},energy_Eh,relative_energy_{1},reference_point,converged,optimization_steps\n",if self.options.coordinate==Coordinate::Bond{"angstrom"}else{"degree"},view.unit.label().replace('/',"_"));
        for i in view.indices(self) {
            let p = &self.points[i];
            text.push_str(&format!(
                "{},{:.8},{:.8},{:.12},{:.8},{},{},{}\n",
                i + 1,
                p.target,
                p.actual,
                p.energy,
                view.energy(self, i),
                view.reference_index(self) + 1,
                p.converged
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "n/a".into()),
                p.steps
            ));
        }
        text
    }
}
#[derive(Clone, Default)]
struct WindowState {
    open: bool,
    plot_open: bool,
    report: Option<Arc<ScanReport>>,
    calculation_dir: Option<std::path::PathBuf>,
    selected: usize,
    view: ScanView,
}
fn state(ctx: &egui::Context) -> WindowState {
    ctx.data_mut(|d| d.get_temp(egui::Id::new("qc_scan_results")))
        .unwrap_or_default()
}
fn store(ctx: &egui::Context, s: WindowState) {
    ctx.data_mut(|d| d.insert_temp(egui::Id::new("qc_scan_results"), s));
}
impl WindowState {
    fn sync(&mut self, task: &QcRunTask) {
        if let Some(report) = task.last.as_ref().and_then(|o| o.scan.as_ref()) {
            if !self
                .report
                .as_ref()
                .is_some_and(|old| Arc::ptr_eq(old, report))
            {
                self.report = Some(report.clone());
                self.calculation_dir = task.last.as_ref().and_then(|o| o.kept_dir.clone());
                self.open = true;
                self.plot_open = !report.points.is_empty();
                self.selected = 0;
                self.view.reversed = false;
            }
        }
    }
}
pub(crate) fn is_open(ctx: &egui::Context, task: &QcRunTask) -> bool {
    let mut s = state(ctx);
    s.sync(task);
    s.open
}
pub(super) fn reopen_button(ui: &mut egui::Ui, task: &QcRunTask) {
    let mut s = state(ui.ctx());
    s.sync(task);
    if s.report.is_some() && ui.button("Open scan results").clicked_once() {
        s.open = true;
        store(ui.ctx(), s);
    }
}
#[derive(Default)]
pub(crate) struct Actions {
    pub loaded: bool,
    pub trajectory: bool,
}
fn load(
    report: &ScanReport,
    indices: impl Iterator<Item = usize>,
    mol: &mut crate::molecule::Molecule,
    traj: &mut crate::trajectory::TrajectoryState,
) -> bool {
    let frames = indices
        .map(|i| crate::trajectory::TrajectoryFrame { metadata: crate::trajectory_data::Metadata::energy(report.points[i].energy, false, format!("Scan point {} E = {:.14} Eh", i+1, report.points[i].energy)),
            atoms: report.atoms.clone(),
            pos: report.points[i]
                .coordinates
                .chunks_exact(3)
                .map(|p| bevy::prelude::Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32))
                .collect(),
        })
        .collect();
    crate::ts_generation::show_frames(frames, mol, traj)
}
fn export(ui: &mut egui::Ui, label: &str, name: &str, text: &str) {
    if ui.button(label).clicked_once() {
        if let Some(path) = crate::recent_dir::save_file(crate::recent_dir::save(name)) {
            if let Err(e) = std::fs::write(path, text) {
                crate::runtime_warning::observe(
                    ui.ctx(),
                    "scan_export",
                    "Cannot save scan results",
                    Some(&e.to_string()),
                );
            }
        }
    }
}

pub(crate) fn window(
    ctx: &egui::Context,
    task: &QcRunTask,
    mol: &mut crate::molecule::Molecule,
    traj: &mut crate::trajectory::TrajectoryState,
) -> (Vec<egui::Rect>, Actions) {
    let mut s = state(ctx);
    s.sync(task);
    let mut actions = Actions::default();
    let mut rects = vec![];
    if let Some(report) = s.report.clone() {
        if s.open {
            let mut open = s.open;
            if let Some(window) = egui::Window::new("Scan Results")
                .id(egui::Id::new("qc_scan_results_window"))
                .frame(crate::ui_style::window_frame(ctx))
                .open(&mut open)
                .resizable(true)
                .default_size([540., 460.])
                .default_pos(ctx.content_rect().left_top() + egui::vec2(20., 60.))
                .vscroll(true)
                .show(ctx, |ui| {
                    crate::ui_style::modern(ui);
                    ui.label(format!(
                        "{} · {} · {}",
                        report.program, report.method, report.driver
                    ));
                    ui.label(&report.message);
                    super::calculation_export::button(ui, s.calculation_dir.as_deref(), true);
                    if !report.points.is_empty() {
                        s.view.controls(ui);
                        ui.horizontal_wrapped(|ui| {
                            if ui.button("Show Plot").clicked_once() {
                                s.plot_open = true;
                            }
                            if ui.button("Copy data").clicked_once() {
                                ui.ctx().copy_text(report.csv_displayed(s.view));
                            }
                            export(ui, "Save CSV…", "scan.csv", &report.csv_displayed(s.view));
                            export(ui, "Save XYZ…", "scan.xyz", &report.xyz_ordered(s.view));
                        });
                        ui.horizontal_wrapped(|ui| {
                            if ui.button("Show selected").clicked_once() {
                                actions.loaded =
                                    load(&report, std::iter::once(s.selected), mol, traj);
                            }
                            if ui.button("Load as Trajectory").clicked_once() {
                                actions.loaded =
                                    load(&report, s.view.indices(&report).into_iter(), mol, traj);
                                actions.trajectory = actions.loaded;
                            }
                        });
                        egui::ScrollArea::both()
                            .max_height(ui.available_height().max(150.))
                            .show(ui, |ui| {
                                egui::Grid::new("scan_points").striped(true).show(ui, |ui| {
                                    for name in [
                                        "Point".into(),
                                        format!("Target ({})", report.options.coordinate.unit()),
                                        format!("Actual ({})", report.options.coordinate.unit()),
                                        format!("ΔE ({})", s.view.unit.label()),
                                        "Status".into(),
                                    ] {
                                        ui.strong(name);
                                    }
                                    ui.end_row();
                                    for i in s.view.indices(&report) {
                                        let p = &report.points[i];
                                        if ui
                                            .selectable_label(s.selected == i, (i + 1).to_string())
                                            .clicked_once()
                                        {
                                            s.selected = i;
                                            actions.loaded =
                                                load(&report, std::iter::once(i), mol, traj);
                                        }
                                        ui.monospace(format!("{:.5}", p.target));
                                        ui.monospace(format!("{:.5}", p.actual));
                                        ui.monospace(format!("{:.6}", s.view.energy(&report, i)));
                                        ui.label(match p.converged {
                                            Some(true) => "Converged",
                                            Some(false) => "Not converged",
                                            None => "Energy",
                                        });
                                        ui.end_row();
                                    }
                                });
                            });
                    }
                })
            {
                rects.push(window.response.rect);
            }
            s.open = open;
        }
        if s.plot_open && !report.points.is_empty() {
            let mut open = s.plot_open;
            if let Some(window) = egui::Window::new("Scan Plot")
                .id(egui::Id::new("qc_scan_plot_window"))
                .frame(crate::ui_style::window_frame(ctx))
                .open(&mut open)
                .resizable(true)
                .default_size([600., 420.])
                .default_pos(egui::pos2(
                    (ctx.content_rect().right() - 640.).max(20.),
                    ctx.content_rect().top() + 60.,
                ))
                .show(ctx, |ui| {
                    crate::ui_style::modern(ui);
                    s.view.controls(ui);
                    if let Some(index) = plot::draw(ui, &report, s.selected, s.view) {
                        s.selected = index;
                    }
                })
            {
                rects.push(window.response.rect);
            }
            s.plot_open = open;
        }
    }
    store(ctx, s);
    (rects, actions)
}
#[path = "scan_plot.rs"]
mod plot;
#[cfg(test)]
#[path = "scan_results_tests.rs"]
mod tests;
