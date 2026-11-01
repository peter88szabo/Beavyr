//! Completed RDA and Poor Man's NEB searches have an independent results window.
#[cfg(test)]
#[path = "results_window_tests.rs"]
mod tests;
use super::{backend, backend::RunControl, parameters::Algorithm, show_frames, TsGeneration};
use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
use crate::molecule::Molecule;
use crate::trajectory::{TrajectoryFrame, TrajectoryState};
use crate::ui_style::ResponseExt;
use bevy::prelude::Vec3;
use bevy_egui::egui;
use std::sync::{Arc, Weak};

#[derive(Clone, Default)]
struct WindowState {
    open: bool,
    completed_run: Weak<RunControl>,
}
impl WindowState {
    fn sync(&mut self, state: &TsGeneration) {
        if state.output.is_none() {
            *self = Self::default();
        } else {
            // Each start creates a new control. Keeping a weak reference identifies
            // the completed run without retaining its worker or reopening on each frame.
            let run = Arc::downgrade(&state.control);
            if !self.completed_run.ptr_eq(&run) {
                self.completed_run = run;
                self.open = true;
            }
        }
    }
}
fn state(ctx: &egui::Context) -> WindowState {
    ctx.data_mut(|data| data.get_temp(egui::Id::new("ts_results_window_state")))
        .unwrap_or_default()
}
fn store(ctx: &egui::Context, value: WindowState) {
    ctx.data_mut(|data| data.insert_temp(egui::Id::new("ts_results_window_state"), value));
}
pub(crate) fn is_open(ctx: &egui::Context, generation: &TsGeneration) -> bool {
    let mut value = state(ctx);
    value.sync(generation);
    value.open
}
pub(super) fn reopen_button(ui: &mut egui::Ui, generation: &TsGeneration) {
    if generation.output.is_some() && ui.button("Open TS results").clicked_once() {
        let mut value = state(ui.ctx());
        value.sync(generation);
        value.open = true;
        store(ui.ctx(), value);
    }
}

/// Called independently of the setup panel; returns whether an explicit result
/// action changed the displayed structure, plus the area excluded from picking.
pub(crate) fn window(
    ctx: &egui::Context,
    generation: &mut TsGeneration,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
) -> (Option<egui::Rect>, bool) {
    let mut value = state(ctx);
    value.sync(generation);
    let mut displayed = false;
    let rect = if value.open {
        let available = ctx.content_rect();
        let size = egui::vec2(
            560.0_f32.min(available.width()),
            600.0_f32.min(available.height()),
        );
        egui::Window::new("TS Generation Results")
            .id(egui::Id::new("ts_generation_results"))
            .frame(crate::ui_style::window_frame(ctx))
            .open(&mut value.open)
            .resizable(true)
            .default_size(size)
            .default_pos(available.center() - size / 2.0)
            .vscroll(true)
            .show(ctx, |ui| {
                crate::ui_style::modern(ui);
                displayed = results_body(ui, generation, mol, traj);
            })
            .map(|response| response.response.rect)
    } else {
        None
    };
    store(ctx, value);
    (rect, displayed)
}

fn results_body(
    ui: &mut egui::Ui,
    state: &mut TsGeneration,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
) -> bool {
    let mut displayed = false;
    if let Some(output) = &state.output {
        ui.strong(format!(
            "{} result · {}",
            output.algorithm.label(),
            output.method
        ));
        ui.label(format!("TS guess energy: {:.9} Eh", output.energy_hartree));
        ui.label(&output.guess.message);
        ui.small("This is a TS guess. Optimize the saddle point and check its imaginary mode before treating it as a confirmed transition state.");
        ui.horizontal_wrapped(|ui| {
            if ui.button("Show TS guess").clicked_once() {
                let positions = output
                    .guess
                    .q
                    .as_chunks::<3>()
                    .0
                    .iter()
                    .map(|q| {
                        Vec3::new(
                            (q[0] * BOHR_TO_ANGSTROM) as f32,
                            (q[1] * BOHR_TO_ANGSTROM) as f32,
                            (q[2] * BOHR_TO_ANGSTROM) as f32,
                        )
                    })
                    .collect();
                displayed = show_frames(
                    vec![TrajectoryFrame { metadata: Default::default(),
                        atoms: output.guess.atoms.clone(),
                        pos: positions,
                    }],
                    mol,
                    traj,
                );
            }
            save(
                ui,
                "Save TS guess XYZ…",
                &format!("{}_ts_guess.xyz", output.algorithm.file_stem()),
                &output.guess_xyz,
                "xyz",
                &mut state.message,
                &mut state.is_error,
            );
        });
        ui.horizontal_wrapped(|ui| {
            if ui.button("Show path in Trajectory").clicked_once() {
                displayed = show_frames(output.frames.clone(), mol, traj);
            }
            save(
                ui,
                if output.algorithm == Algorithm::PoorMansNeb {
                    "Save NEB trajectory XYZ…"
                } else {
                    "Save reaction path XYZ…"
                },
                &format!("{}_path.xyz", output.algorithm.file_stem()),
                &output.path_xyz,
                "xyz",
                &mut state.message,
                &mut state.is_error,
            );
        });
        ui.label(format!(
            "{} XYZ frames, including both endpoints.",
            output.frames.len()
        ));
        ui.horizontal_wrapped(|ui| {
            if let Some(profile) = &output.profile_dat {
                if ui.button("Plot NEB energy").clicked_once() {
                    state.energy_plot_open = true;
                }
                save(
                    ui,
                    "Save NEB energy profile…",
                    "poormans_neb_profile.dat",
                    profile,
                    "dat",
                    &mut state.message,
                    &mut state.is_error,
                );
            } else {
                save(
                    ui,
                    "Save RDA search XYZ…",
                    "rda_search.xyz",
                    &output.search_xyz,
                    "xyz",
                    &mut state.message,
                    &mut state.is_error,
                );
            }
        });
        egui::CollapsingHeader::new("Search points / relaxed images").show(ui, |ui| {
            egui::ScrollArea::vertical()
                .max_height(230.0)
                .show(ui, |ui| {
                    egui::Grid::new("ts_result_points")
                        .striped(true)
                        .show(ui, |ui| {
                            ui.strong("Point");
                            ui.strong("Energy (Eh)");
                            ui.strong("Steps");
                            ui.end_row();
                            for point in &output.guess.points {
                                ui.label(&point.label).on_hover_text(&point.message);
                                ui.monospace(format!("{:.8}", point.energy));
                                ui.label(point.nsteps.to_string());
                                ui.end_row();
                            }
                        });
                });
        });
    }
    displayed
}

fn save(
    ui: &mut egui::Ui,
    label: &str,
    name: &str,
    text: &str,
    extension: &str,
    message: &mut Option<String>,
    is_error: &mut bool,
) {
    if !ui.button(label).clicked_once() {
        return;
    }
    let dialog = crate::recent_dir::save(name).add_filter(extension, &[extension]);
    if let Some(mut path) = crate::recent_dir::save_file(dialog) {
        if !path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case(extension))
        {
            path.set_extension(extension);
        }
        match backend::save_text(&path, text) {
            Ok(()) => {
                *message = Some(format!("Saved {}", path.display()));
                *is_error = false;
            }
            Err(error) => {
                *message = Some(error.to_string());
                *is_error = true;
            }
        }
    }
}
