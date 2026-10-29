//! Completed conformer searches have their own window, independent of the setup panel.

use super::*;
use std::time::Instant;

#[derive(Default)]
pub(crate) struct ResultActions {
    pub loaded: bool,
    pub open_trajectory: bool,
}

#[derive(Clone, Default)]
struct WindowState {
    open: bool,
    seen_result: bool,
    completed_run: Option<Instant>,
}

impl WindowState {
    fn sync(&mut self, run: &ConformerRun) {
        if run.outcome.is_none() {
            self.open = false;
            self.seen_result = false;
        } else if !self.seen_result || self.completed_run != run.started {
            self.open = true;
            self.seen_result = true;
            self.completed_run = run.started;
        }
    }
}

fn state(ctx: &egui::Context) -> WindowState {
    ctx.data_mut(|data| {
        data.get_temp::<WindowState>(egui::Id::new("conformer_results_window_state"))
    })
    .unwrap_or_default()
}

fn store(ctx: &egui::Context, value: WindowState) {
    ctx.data_mut(|data| data.insert_temp(egui::Id::new("conformer_results_window_state"), value));
}

pub(crate) fn is_open(ctx: &egui::Context, run: &ConformerRun) -> bool {
    let mut value = state(ctx);
    value.sync(run);
    value.open
}

pub(super) fn reopen_button(ui: &mut egui::Ui, run: &ConformerRun) {
    if run.outcome.is_some() && ui.button("Open conformer results").clicked_once() {
        let mut value = state(ui.ctx());
        value.sync(run);
        value.open = true;
        store(ui.ctx(), value);
    }
}

pub(crate) fn window(
    ctx: &egui::Context,
    run: &mut ConformerRun,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
    settings: &MolSettings,
) -> (Option<egui::Rect>, ResultActions) {
    let mut value = state(ctx);
    value.sync(run);
    let mut actions = ResultActions::default();
    let rect = if value.open {
        let outcome = run
            .outcome
            .clone()
            .expect("an open result window has an outcome");
        let available = ctx.content_rect();
        let size = egui::vec2(
            520.0_f32.min((available.width() - 32.0).max(300.0)),
            680.0_f32.min((available.height() - 32.0).max(240.0)),
        );
        egui::Window::new("Conformer Search Results")
            .id(egui::Id::new("conformer_results_window_v2"))
            .open(&mut value.open)
            .resizable(true)
            .default_size(size)
            .default_pos(available.center() - size / 2.0)
            .vscroll(true)
            .show(ctx, |ui| {
                crate::ui_style::modern(ui);
                actions = results_section(ui, run, &outcome, mol, traj, settings);
            })
            .map(|response| response.response.rect)
    } else {
        None
    };
    store(ctx, value);
    (rect, actions)
}

fn results_section(
    ui: &mut egui::Ui,
    run: &mut ConformerRun,
    outcome: &ConformerOutcome,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
    settings: &MolSettings,
) -> ResultActions {
    let mut actions = ResultActions::default();

    // Rings and rotatable bonds are both degrees of freedom but they are not the same thing, and
    // calling a ring a "rotatable bond" would be wrong in a way a chemist would notice.
    let freedom = match (outcome.rotors, outcome.rings) {
        (0, 1) => "1 ring".to_string(),
        (0, rings) => format!("{rings} rings"),
        (rotors, 0) => format!(
            "{rotors} rotatable bond{}",
            if rotors == 1 { "" } else { "s" }
        ),
        (rotors, rings) => format!(
            "{rotors} rotatable bond{} and {rings} ring{}",
            if rotors == 1 { "" } else { "s" },
            if rings == 1 { "" } else { "s" }
        ),
    };
    ui.label(format!(
        "{} conformer{} from {freedom}, in {:.1} s on {} core{}.",
        outcome.conformers.len(),
        if outcome.conformers.len() == 1 {
            ""
        } else {
            "s"
        },
        outcome.elapsed.as_secs_f64(),
        outcome.workers,
        if outcome.workers == 1 { "" } else { "s" }
    ));
    ui.label(
        egui::RichText::new(format!(
            "{} generations, {} structures tried, {} local optimisations. Finished: {}.",
            outcome.generations,
            outcome.structures_tried,
            outcome.local_optimizations,
            outcome.termination
        ))
        .small()
        .weak(),
    );
    // If the winner only appeared at the very end, the search had not settled and more effort is
    // likely to find something better. Saying so is more use than the raw generation number.
    if outcome.generations > 1 && outcome.best_found_in_generation >= outcome.generations {
        ui.label(
            egui::RichText::new(
                "The best conformer was found in the last generation, so the search was still \
                 improving when it stopped -- try a higher search effort.",
            )
            .small()
            .color(egui::Color32::from_rgb(200, 160, 70)),
        );
    }

    if let Some(warning) = &outcome.radical_warning {
        ui.add_space(4.0);
        ui.label(
            egui::RichText::new(warning)
                .small()
                .color(egui::Color32::from_rgb(200, 160, 70)),
        );
    }

    refinement_row(ui, run, outcome);

    ui.add_space(8.0);
    ui.horizontal(|ui| {
        if ui
            .button("Load as trajectory")
            .on_hover_text(
                "Puts every conformer in the Trajectory tool, so they can be stepped through \
                 or played as an animation.",
            )
            .clicked_once()
        {
            load_as_trajectory(outcome, mol, traj, settings);
            actions.loaded = true;
            actions.open_trajectory = true;
        }
        if ui
            .button("Save as XYZ trajectory")
            .on_hover_text("One XYZ frame per conformer, ordered by energy.")
            .clicked_once()
        {
            save_trajectory(ui.ctx(), outcome);
        }
    });

    ui.add_space(8.0);
    egui::ScrollArea::vertical()
        .max_height(440.0)
        .show(ui, |ui| {
            egui::Grid::new("conformer_table")
                .num_columns(4)
                .spacing([18.0, 8.0])
                .striped(true)
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("#").strong());
                    ui.label(egui::RichText::new("ΔE (kcal/mol)").strong());
                    ui.label(egui::RichText::new("Population").strong());
                    ui.label("");
                    ui.end_row();

                    for (index, conformer) in outcome.conformers.iter().enumerate() {
                        let showing = run.previewing == Some(index);
                        ui.label(format!("{}", index + 1));
                        ui.label(format!("{:.2}", conformer.relative_energy_kcal));
                        // Below a tenth of a percent the number stops being informative.
                        if conformer.population_fraction >= 0.001 {
                            ui.label(format!("{:.1} %", conformer.population_fraction * 100.0));
                        } else {
                            ui.label(egui::RichText::new("< 0.1 %").weak());
                        }
                        if ui
                            .add(
                                egui::Button::new(if showing { "Viewing" } else { "View" })
                                    .selected(showing),
                            )
                            .clicked_once()
                        {
                            // Results retain the searched structure even if the user edited
                            // or replaced the molecule while the search was running.
                            mol.topology = None;
                            mol.atoms = outcome.atoms.clone();
                            traj.clear_for_structure();
                            mol.set_pos(
                                conformer
                                    .positions_angstrom
                                    .chunks_exact(3)
                                    .map(|c| Vec3::new(c[0] as f32, c[1] as f32, c[2] as f32))
                                    .collect(),
                            );
                            mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                            run.previewing = Some(index);
                            actions.loaded = true;
                        }
                        ui.end_row();
                    }
                });
        });

    ui.add_space(4.0);
    ui.label(
        egui::RichText::new(
            "Populations are Boltzmann factors at 298 K, without vibrational entropy or symmetry.",
        )
        .small()
        .weak(),
    );

    ui.add_space(8.0);
    ui.separator();
    ui.add_space(6.0);
    ui.collapsing("References", |ui| citation(ui));

    actions
}

#[cfg(test)]
#[path = "results_window_tests.rs"]
mod tests;
