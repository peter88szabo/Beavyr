//! Dataset inspection changes a selection, never the source sampling output.
use crate::{
    dataset_filter::{self, Dataset, Options, Row},
    trajectory::TrajectoryState,
    ui_style::ResponseExt,
};
use bevy::{
    prelude::*,
    tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task},
};
use bevy_egui::{egui, EguiContexts, EguiPrimaryContextPass};
use std::sync::{
    atomic::{AtomicBool, AtomicUsize, Ordering},
    Arc,
};
const OPEN: &str = "open_dataset_inspector";
pub struct Plugin;
impl bevy::prelude::Plugin for Plugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<State>()
            .add_systems(EguiPrimaryContextPass, ui.after(crate::ui::ui_panel));
    }
}
enum ResultData {
    Rows(Vec<Row>),
    Saved(std::path::PathBuf),
}
#[derive(Resource)]
struct State {
    open: bool,
    data: Option<Dataset>,
    rows: Vec<Row>,
    undo: Option<Vec<Row>>,
    options: Options,
    target: usize,
    task: Option<Task<Result<ResultData, String>>>,
    cancel: Arc<AtomicBool>,
    progress: Arc<AtomicUsize>,
    total: usize,
    saving: bool,
    message: String,
}
impl Default for State {
    fn default() -> Self {
        Self {
            open: false,
            data: None,
            rows: vec![],
            undo: None,
            options: Default::default(),
            target: 100,
            task: None,
            cancel: Default::default(),
            progress: Default::default(),
            total: 0,
            saving: false,
            message: String::new(),
        }
    }
}
pub fn button(ui: &mut egui::Ui, pes: bool) {
    if ui.button("Inspect dataset…").clicked_once() {
        ui.ctx()
            .data_mut(|d| d.insert_temp(egui::Id::new(OPEN), if pes { 2u8 } else { 1u8 }));
    }
}
fn analyze(state: &mut State, diversity: bool) {
    let Some(data) = state.data.clone() else {
        return;
    };
    let rows = state.rows.clone();
    let options = state.options.clone();
    let target = state.target;
    state.cancel = Default::default();
    state.progress = Default::default();
    let cancel = state.cancel.clone();
    let progress = state.progress.clone();
    state.saving = false;
    state.total = if diversity {
        target.min(rows.iter().filter(|r| r.keep).count())
    } else {
        data.frames.len()
    };
    state.task = Some(AsyncComputeTaskPool::get().spawn(async move {
        if diversity {
            dataset_filter::diverse(&data, &rows, target, &cancel, &progress)
        } else {
            dataset_filter::analyze(&data, &options, &cancel, &progress)
        }
        .map(ResultData::Rows)
    }));
}
fn ui(
    mut contexts: EguiContexts,
    mut state: ResMut<State>,
    pes: Res<crate::pes_sampling::PesSampling>,
    mut traj: ResMut<TrajectoryState>,
    mut mol: ResMut<crate::molecule::Molecule>,
    settings: Res<crate::settings::MolSettings>,
    mut events: MessageWriter<crate::events::MoleculeChanged>,
    mut layout: ResMut<crate::ui_layout::UiLayout>,
    mut history: ResMut<crate::structure_history::StructureHistory>,
    mut regions: ResMut<crate::ui::UiPanelRegions>,
) {
    let Ok(ctx) = contexts.ctx_mut() else {
        return;
    };
    if let Some(task) = &mut state.task {
        if let Some(result) = block_on(future::poll_once(task)) {
            state.task = None;
            match result {
                Ok(ResultData::Rows(rows)) => {
                    state.undo = Some(std::mem::replace(&mut state.rows, rows));
                    state.message = "Analysis complete. Selection can be edited or undone.".into();
                }
                Ok(ResultData::Saved(path)) => state.message = format!("Saved {}", path.display()),
                Err(e) => crate::runtime_warning::push(ctx, "Dataset inspector", &e),
            }
        } else {
            ctx.request_repaint();
        }
    }
    if let Some(source) = ctx.data_mut(|d| d.remove_temp::<u8>(egui::Id::new(OPEN))) {
        state.open = true;
        if state.task.is_none() {
            let data = if source == 2 {
                children_pes(&pes)
            } else {
                children_pes(&pes)
                    .ok()
                    .filter(|data| {
                        crate::trajectory_data::fingerprint(&data.frames)
                            == crate::trajectory_data::fingerprint(&traj.frames)
                    })
                    .map(Ok)
                    .unwrap_or_else(|| Dataset::from_trajectory(&traj))
            };
            match data {
                Ok(data) => {
                    state.rows = vec![Row::default(); data.frames.len()];
                    state.data = Some(data);
                    state.undo = None;
                    state.message.clear();
                    analyze(&mut state, false);
                }
                Err(e) => crate::runtime_warning::push(ctx, "Dataset inspector", &e),
            }
        }
    }
    if !state.open {
        return;
    }
    let mut open = true;
    let mut view = None;
    let mut load = false;
    let window = egui::Window::new("PES dataset inspector")
        .open(&mut open)
        .resizable(true)
        .default_size([760., 560.])
        .show(ctx, |ui| {
            crate::ui_style::modern(ui);
            let busy = state.task.is_some();
            ui.add_enabled_ui(!busy, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.checkbox(&mut state.options.overlap.enabled, "Overlap filter");
                    for (label, value) in [("H–H", &mut state.options.overlap.hydrogen_hydrogen)]
                    {
                        ui.label(label);
                        ui.add(egui::DragValue::new(value).range(0.0..=5.).speed(0.01));
                    }
                    ui.label("H–other");
                    ui.add(
                        egui::DragValue::new(&mut state.options.overlap.hydrogen_other)
                            .range(0.0..=5.)
                            .speed(0.01),
                    );
                    ui.label("Other–other (Å)");
                    ui.add(
                        egui::DragValue::new(&mut state.options.overlap.other_other)
                            .range(0.0..=5.)
                            .speed(0.01),
                    );
                });
                ui.horizontal_wrapped(|ui| {
                    ui.checkbox(&mut state.options.duplicates, "Remove RMSD duplicates");
                    ui.add(
                        egui::DragValue::new(&mut state.options.rmsd)
                            .range(0.0..=10.)
                            .speed(0.005)
                            .suffix(" Å"),
                    );
                    if ui.button("Analyse / reset filters").clicked_once() {
                        analyze(&mut state, false);
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Diverse subset size");
                    ui.add(egui::DragValue::new(&mut state.target).range(1..=1_000_000));
                    if ui.button("Select diverse subset").clicked_once() {
                        analyze(&mut state, true);
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    if ui.button("Keep all").clicked_once() {
                        state.undo = Some(state.rows.clone());
                        for r in &mut state.rows {
                            r.keep = true;
                        }
                    }
                    if ui.button("Keep none").clicked_once() {
                        state.undo = Some(state.rows.clone());
                        for r in &mut state.rows {
                            r.keep = false;
                        }
                    }
                    if ui
                        .add_enabled(state.undo.is_some(), egui::Button::new("Undo selection"))
                        .clicked_once()
                    {
                        if let Some(rows) = state.undo.take() {
                            let current = std::mem::replace(&mut state.rows, rows);
                            state.undo = Some(current);
                        }
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    if ui.button("Load selected as Trajectory").clicked_once() {
                        load = true;
                    }
                    if ui.button("Save selected dataset…").clicked_once() {
                        if let Some(data) = state.data.clone() {
                            if let Some(parent) = rfd::FileDialog::new()
                                .set_title("Choose a folder for a new dataset export")
                                .pick_folder()
                            {
                                let rows = state.rows.clone();
                                state.saving = true;
                                state.task = Some(AsyncComputeTaskPool::get().spawn(async move {
                                    dataset_filter::save(&data, &rows, &parent)
                                        .map(ResultData::Saved)
                                }));
                            }
                        }
                    }
                    if ui.button("Copy selection CSV").clicked_once() {
                        if let Some(data) = &state.data {
                            ctx.copy_text(data.manifest(&state.rows));
                        }
                    }
                });
            });
            if busy {
                ui.horizontal(|ui| {
                    ui.spinner();
                    if state.saving {
                        ui.label("Saving dataset…");
                    } else {
                        ui.label(format!(
                            "{} / {}",
                            state.progress.load(Ordering::Relaxed),
                            state.total
                        ));
                        if ui.button("Cancel").clicked_once() {
                            state.cancel.store(true, Ordering::Relaxed);
                        }
                    }
                });
            }
            if !state.message.is_empty() {
                ui.add(egui::Label::new(&state.message).selectable(true));
            }
            ui.label(format!(
                "{} / {} frames selected · RMSD uses matching atom order",
                state.rows.iter().filter(|r| r.keep).count(),
                state.rows.len()
            ));
            if let Some(data) = state.data.clone() {
                ui.small(&data.source);
                let mut changed = false;
                egui::ScrollArea::both()
                    .auto_shrink([false, false])
                    .show_rows(ui, 24., state.rows.len(), |ui, range| {
                        for i in range {
                            ui.horizontal(|ui| {
                                let mut keep = state.rows[i].keep;
                                if ui
                                    .add_enabled(!busy, egui::Checkbox::without_text(&mut keep))
                                    .changed()
                                {
                                    if !changed {
                                        state.undo = Some(state.rows.clone());
                                        changed = true;
                                    }
                                    state.rows[i].keep = keep;
                                }
                                if ui
                                    .selectable_label(false, format!("Frame {}", i + 1))
                                    .clicked_once()
                                {
                                    view = Some(i);
                                }
                                ui.label(&state.rows[i].reason);
                                if let Some(j) = state.rows[i].duplicate_of {
                                    ui.label(format!("of {}", j + 1));
                                }
                                if let Some(r) = state.rows[i].nearest_rmsd {
                                    ui.monospace(format!("{r:.5} Å"));
                                }
                                ui.add(
                                    egui::Label::new(&data.frames[i].metadata.comment).truncate(),
                                );
                            });
                        }
                    });
            }
        });
    state.open = open;
    if let Some(w) = window {
        let r = w.response.rect;
        regions
            .windows
            .push(Rect::new(r.min.x, r.min.y, r.max.x, r.max.y));
    }
    if load || view.is_some() {
        if let Some(data) = &state.data {
            let indices: Vec<_> = if let Some(i) = view {
                vec![i]
            } else {
                state
                    .rows
                    .iter()
                    .enumerate()
                    .filter(|(_, r)| r.keep)
                    .map(|(i, _)| i)
                    .collect()
            };
            if indices.is_empty() {
                crate::runtime_warning::push(ctx, "Dataset inspector", "No frames are selected.");
                return;
            }
            let previous = history.before(&mol, &traj);
            let frames = indices
                .into_iter()
                .map(|i| data.frames[i].clone())
                .collect();
            if crate::ts_generation::show_frames(frames, &mut mol, &mut traj) {
                history.replaced(previous, &mol, &traj, "Dataset selection");
                mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                events.write(crate::events::MoleculeChanged::parse_xyz(false));
                layout.open[crate::ui_layout::Tab::Trajectory.index()] = true;
            }
        }
    }
}
fn children_pes(pes: &crate::pes_sampling::PesSampling) -> Result<Dataset, String> {
    let out = pes.output.as_ref().ok_or("Generate PES samples first.")?;
    let mut bytes = vec![];
    out.write_xyz(&mut bytes)?;
    Dataset::from_xyz(
        format!("PES sampling seed {}", out.seed),
        std::str::from_utf8(&bytes).map_err(|e| e.to_string())?,
    )
}
