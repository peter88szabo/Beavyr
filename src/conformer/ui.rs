//! The Conformer Search panel.
//!
//! Deliberately short on controls. The number of rotatable bonds, which torsions to sample, the
//! population size and the generation count are all worked out from the structure, because they
//! are not things a chemist should have to guess at -- and a wrong guess quietly limits what the
//! search can find. What is left is how hard to look, and whether to re-rank the answers with
//! xTB.

use bevy::prelude::*;
use bevy_egui::egui;

use crate::events::MoleculeChanged;
use crate::molecule::Molecule;
use crate::settings::MolSettings;
use crate::trajectory::{TrajectoryFrame, TrajectoryState};

use super::{available_cores, ConformerOutcome, ConformerRun, Thoroughness};

/// Draws the panel. Returns `true` if the conformers were just loaded into the viewer, so the
/// caller can announce the structure change.
pub fn conformer_panel(
    ui: &mut egui::Ui,
    run: &mut ConformerRun,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
    settings: &MolSettings,
) -> bool {
    let mut loaded = false;

    ui.label(
        egui::RichText::new(
            "Finds the low-energy shapes a flexible molecule can take, by searching over its \
             rotatable bonds.",
        )
        .small(),
    );
    ui.add_space(6.0);

    // --- Method -------------------------------------------------------------------------------
    ui.label(egui::RichText::new("Energies from").strong());
    ui.horizontal(|ui| {
        ui.label("DREIDING force field");
        ui.label(egui::RichText::new("(built in)").weak().small());
    });
    ui.label(
        egui::RichText::new(
            "A search runs tens of thousands of energy evaluations. DREIDING answers each in \
             microseconds because it runs inside Beavyr; an external program costs a process \
             launch per evaluation, which would turn seconds into hours.",
        )
        .small()
        .weak(),
    );

    ui.add_space(6.0);
    // Shown but not yet available, rather than hidden: the plan is visible, and the checkbox
    // cannot mislead by appearing to do something. See `conformer::mod` for what it needs.
    ui.add_enabled_ui(false, |ui| {
        ui.checkbox(
            &mut run.settings.refine_with_gfnff,
            "Re-rank the results with xTB GFN-FF",
        );
    });
    ui.label(
        egui::RichText::new(
            "Not yet available. Re-optimising only the surviving conformers would be a few \
             hundred xTB calls rather than tens of thousands, and relative conformer energies \
             are where a generic force field is weakest -- so it is worth doing, but xTB is \
             currently driven through a path that writes fixed filenames into the working \
             directory and cannot safely be called in a loop.",
        )
        .small()
        .weak(),
    );

    ui.add_space(8.0);
    ui.separator();
    ui.add_space(8.0);

    // --- How hard to look ---------------------------------------------------------------------
    ui.label(egui::RichText::new("Search effort").strong());
    ui.horizontal(|ui| {
        for level in [
            Thoroughness::Quick,
            Thoroughness::Normal,
            Thoroughness::Thorough,
        ] {
            if ui
                .selectable_label(run.settings.thoroughness == level, level.label())
                .on_hover_text(level.description())
                .clicked()
            {
                run.settings.thoroughness = level;
            }
        }
    });
    ui.label(
        egui::RichText::new(run.settings.thoroughness.description())
            .small()
            .weak(),
    );

    ui.add_space(4.0);
    ui.checkbox(&mut run.settings.reproducible, "Repeatable")
        .on_hover_text(
            "Fixes the random seed, so running the search again on the same structure gives \
             the same conformers.",
        );

    ui.add_space(4.0);
    let cores = available_cores();
    ui.horizontal(|ui| {
        ui.label("Cores");
        ui.add_enabled(
            !run.is_running(),
            egui::DragValue::new(&mut run.settings.ncore)
                .speed(1.0)
                .range(1..=cores),
        )
        .on_hover_text(
            "How many structures to optimise at once. Every local optimisation is independent, \
             so this scales almost linearly. It does not change the answer: with a fixed seed the \
             same conformers come back on any number of cores.",
        );
        ui.label(egui::RichText::new(format!("of {cores} available")).weak().small());
    });

    ui.add_space(4.0);
    ui.horizontal(|ui| {
        ui.label("Keep within");
        ui.add(
            egui::DragValue::new(&mut run.settings.energy_window_kcal)
                .speed(0.5)
                .range(1.0..=50.0)
                .suffix(" kcal/mol"),
        );
        ui.label("of the best");
    });

    ui.add_space(8.0);
    ui.separator();
    ui.add_space(8.0);

    // --- Run ----------------------------------------------------------------------------------
    ui.horizontal(|ui| {
        let can_run = !run.is_running() && !mol.atoms.is_empty();
        if ui
            .add_enabled(can_run, egui::Button::new("Search for conformers"))
            .clicked()
        {
            run.start(mol);
        }
        if run.is_running() {
            ui.spinner();
            if let Some(elapsed) = run.elapsed() {
                ui.label(format!("{:.1} s", elapsed.as_secs_f64()));
            }
        }
    });

    if mol.atoms.is_empty() {
        ui.label(
            egui::RichText::new("Load or build a structure first.")
                .small()
                .weak(),
        );
    }

    if let Some(error) = &run.error {
        ui.add_space(6.0);
        ui.colored_label(egui::Color32::from_rgb(220, 120, 90), error);
    }

    // --- Results ------------------------------------------------------------------------------
    if let Some(outcome) = run.outcome.clone() {
        ui.add_space(8.0);
        ui.separator();
        ui.add_space(8.0);
        loaded |= results_section(ui, run, &outcome, mol, traj, settings);
    }

    loaded
}

fn results_section(
    ui: &mut egui::Ui,
    run: &mut ConformerRun,
    outcome: &ConformerOutcome,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
    settings: &MolSettings,
) -> bool {
    let mut loaded = false;

    ui.label(egui::RichText::new("Results").strong());
    ui.label(format!(
        "{} conformer{} from {} rotatable bond{}, in {:.1} s on {} core{}.",
        outcome.conformers.len(),
        if outcome.conformers.len() == 1 { "" } else { "s" },
        outcome.rotors,
        if outcome.rotors == 1 { "" } else { "s" },
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

    ui.add_space(8.0);
    ui.horizontal(|ui| {
        if ui
            .button("Load as trajectory")
            .on_hover_text(
                "Puts every conformer in the Trajectory tool, so they can be stepped through \
                 or played as an animation.",
            )
            .clicked()
        {
            load_as_trajectory(outcome, mol, traj, settings);
            loaded = true;
        }
        if ui
            .button("Save as XYZ trajectory")
            .on_hover_text("One XYZ frame per conformer, ordered by energy.")
            .clicked()
        {
            save_trajectory(outcome);
        }
    });

    ui.add_space(8.0);
    egui::ScrollArea::vertical()
        .max_height(260.0)
        .show(ui, |ui| {
            egui::Grid::new("conformer_table")
                .num_columns(4)
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
                        if ui.selectable_label(showing, "Show").clicked() {
                            mol.set_pos(
                                conformer
                                    .positions_angstrom
                                    .chunks_exact(3)
                                    .map(|c| Vec3::new(c[0] as f32, c[1] as f32, c[2] as f32))
                                    .collect(),
                            );
                            mol.recompute_bonds(
                                settings.bond_thresh_scale,
                                settings.hbond_cutoff,
                            );
                            run.previewing = Some(index);
                            loaded = true;
                        }
                        ui.end_row();
                    }
                });
        });

    ui.add_space(4.0);
    ui.label(
        egui::RichText::new(
            "Populations are Boltzmann factors at 298 K from these energies alone. They leave \
             out vibrational entropy and symmetry, so read them as indicative.",
        )
        .small()
        .weak(),
    );

    ui.add_space(8.0);
    ui.separator();
    ui.add_space(6.0);
    citation(ui);

    loaded
}

/// The reference for the search algorithm, shown once a search has produced results.
///
/// The genetic algorithm here is someone else's published work, not ours, so anyone using these
/// conformers in a paper should cite it. Putting this next to the results rather than in a
/// documentation file elsewhere is the only way it reliably gets seen.
fn citation(ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Please cite").strong());
    ui.label(
        egui::RichText::new(
            "The search algorithm is not ours. If you publish these conformers, cite the paper \
             it comes from:",
        )
        .small(),
    );
    ui.add_space(2.0);
    ui.label(
        egui::RichText::new(
            "A. Supady, V. Blum, C. Baldauf, \"First-Principles Molecular Structure Search with \
             a Genetic Algorithm\", J. Chem. Inf. Model. 2015, 55, 2338-2348.",
        )
        .small()
        .italics(),
    );
    ui.add_space(2.0);
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new("doi:10.1021/acs.jcim.5b00243").small().weak());
        if ui
            .small_button("Copy")
            .on_hover_text("Copies the full reference to the clipboard.")
            .clicked()
        {
            ui.ctx().copy_text(CITATION.to_string());
        }
    });
    ui.add_space(4.0);
    ui.label(
        egui::RichText::new(
            "The DREIDING force field used for the energies has its own reference: \
             S. L. Mayo, B. D. Olafson, W. A. Goddard III, J. Phys. Chem. 1990, 94, 8897-8909.",
        )
        .small()
        .weak(),
    );
}

/// The reference in full, for the clipboard.
const CITATION: &str = "A. Supady, V. Blum, C. Baldauf, \"First-Principles Molecular Structure \
Search with a Genetic Algorithm\", J. Chem. Inf. Model. 2015, 55, 2338-2348. \
doi:10.1021/acs.jcim.5b00243\n\
S. L. Mayo, B. D. Olafson, W. A. Goddard III, \"DREIDING: A Generic Force Field for Molecular \
Simulations\", J. Phys. Chem. 1990, 94, 8897-8909. doi:10.1021/j100389a010";

/// Puts the conformers into the Trajectory tool, lowest energy first.
fn load_as_trajectory(
    outcome: &ConformerOutcome,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
    settings: &MolSettings,
) {
    traj.frames = outcome
        .conformers
        .iter()
        .map(|conformer| TrajectoryFrame {
            atoms: outcome.atoms.clone(),
            pos: conformer
                .positions_angstrom
                .chunks_exact(3)
                .map(|c| Vec3::new(c[0] as f32, c[1] as f32, c[2] as f32))
                .collect(),
        })
        .collect();
    traj.current_frame = 0;
    traj.last_applied = None;
    traj.playing = false;
    traj.overlay_dirty = true;
    traj.current_file = None;
    // Conformers differ only in torsion angles, so the bond graph is the same in every frame.
    // This is the case `fixed_bonds` exists for, and it avoids re-perceiving bonds per frame.
    traj.fixed_bonds = true;

    if let Some(first) = traj.frames.first() {
        mol.atoms = first.atoms.clone();
        mol.pos = first.pos.clone();
        mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
    }
}

/// Writes one XYZ frame per conformer, with the relative energy on each comment line.
fn save_trajectory(outcome: &ConformerOutcome) {
    let dialog = rfd::FileDialog::new()
        .set_file_name("conformers.xyz")
        .add_filter("XYZ trajectory", &["xyz"]);
    let Some(path) = crate::recent_dir::save_file(dialog) else {
        return;
    };

    let mut text = String::new();
    for (index, conformer) in outcome.conformers.iter().enumerate() {
        text += &format!("{}\n", outcome.atoms.len());
        // The comment line carries the energy, so the file is still meaningful on its own.
        text += &format!(
            "conformer {} of {}   dE = {:.3} kcal/mol   population = {:.2} %\n",
            index + 1,
            outcome.conformers.len(),
            conformer.relative_energy_kcal,
            conformer.population_fraction * 100.0
        );
        for (symbol, xyz) in outcome
            .atoms
            .iter()
            .zip(conformer.positions_angstrom.chunks_exact(3))
        {
            text += &format!(
                "{:<3} {:>14.8} {:>14.8} {:>14.8}\n",
                symbol, xyz[0], xyz[1], xyz[2]
            );
        }
    }

    if let Err(error) = std::fs::write(&path, text) {
        warn!("could not write {}: {error}", path.display());
    }
}

/// Announces a structure change after the panel replaced the geometry.
pub fn announce(loaded: bool, ev_changed: &mut MessageWriter<MoleculeChanged>) {
    if loaded {
        ev_changed.write(MoleculeChanged::parse_xyz(true));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformer::Conformer;
    use std::time::Duration;

    fn outcome() -> ConformerOutcome {
        ConformerOutcome {
            conformers: vec![
                Conformer {
                    positions_angstrom: vec![0.0, 0.0, 0.0, 1.1, 0.0, 0.0],
                    relative_energy_kcal: 0.0,
                    population_fraction: 0.8,
                },
                Conformer {
                    positions_angstrom: vec![0.0, 0.0, 0.0, 0.0, 1.1, 0.0],
                    relative_energy_kcal: 0.9,
                    population_fraction: 0.2,
                },
            ],
            rotors: 1,
            generations: 4,
            structures_tried: 40,
            local_optimizations: 30,
            termination: "Converged".into(),
            elapsed: Duration::from_millis(1500),
            workers: 4,
            best_found_in_generation: 3,
            atoms: vec!["C".into(), "H".into()],
        }
    }

    #[test]
    fn loading_as_a_trajectory_keeps_every_conformer_in_order() {
        let mut mol = Molecule::empty();
        let mut traj = TrajectoryState::default();
        let settings = MolSettings::default();
        let outcome = outcome();

        load_as_trajectory(&outcome, &mut mol, &mut traj, &settings);

        assert_eq!(traj.frames.len(), 2);
        assert_eq!(traj.current_frame, 0);
        assert!(!traj.playing, "a fresh load should not start animating");
        // The first frame is the lowest-energy conformer, and it is what the viewer shows.
        assert_eq!(mol.atoms, outcome.atoms);
        assert!((mol.pos[1].x - 1.1).abs() < 1.0e-5);
        // Conformers share a bond graph, so it need not be re-perceived per frame.
        assert!(traj.fixed_bonds);
        // A stale frame index from a previous trajectory would show the wrong structure.
        assert_eq!(traj.last_applied, None);
    }

    #[test]
    fn every_frame_has_one_position_per_atom() {
        let mut mol = Molecule::empty();
        let mut traj = TrajectoryState::default();
        load_as_trajectory(&outcome(), &mut mol, &mut traj, &MolSettings::default());
        for frame in &traj.frames {
            assert_eq!(frame.pos.len(), frame.atoms.len());
        }
    }
}
