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

use crate::optimizer::conformer_search::LocalOptimizer;

use super::{available_cores, ConformerOutcome, ConformerRun, SearchEngine, Thoroughness};

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

    // --- Engine -------------------------------------------------------------------------------
    ui.horizontal(|ui| {
        ui.label("Energies");
        egui::ComboBox::from_id_salt("conformer_engine")
            .selected_text(run.settings.engine.label())
            .show_ui(ui, |ui| {
                for engine in SearchEngine::ALL {
                    ui.selectable_value(&mut run.settings.engine, engine, engine.label());
                }
            });
        ui.label(
            egui::RichText::new(run.settings.engine.speed_note())
                .small()
                .weak(),
        );
    });

    if run.settings.engine.needs_xtb() {
        xtb_path_row(ui, run);
        ui.horizontal(|ui| {
            ui.label("Charge");
            ui.add_enabled(
                !run.is_running(),
                egui::DragValue::new(&mut run.settings.charge).range(-10..=10),
            );
            ui.label("Multiplicity");
            ui.add_enabled(
                !run.is_running(),
                egui::DragValue::new(&mut run.settings.multiplicity).range(1..=10),
            );
        });
    }

    ui.add_space(6.0);

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
    ui.add_space(4.0);
    ui.checkbox(&mut run.settings.reproducible, "Repeatable")
        .on_hover_text(
            "Fixes the random seed, so running the search again on the same structure gives \
             the same conformers.",
        );

    ui.add_space(4.0);
    ui.horizontal(|ui| {
        ui.label("Optimiser");
        egui::ComboBox::from_id_salt("conformer_local_optimiser")
            .selected_text(run.settings.local_optimizer.label())
            .show_ui(ui, |ui| {
                for method in LocalOptimizer::ALL {
                    ui.selectable_value(
                        &mut run.settings.local_optimizer,
                        method,
                        method.label(),
                    )
                    .on_hover_text(method.description());
                }
            });
    });

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
    //
    // An external engine needs the xTB executable. Resolved here rather than inside the search, so
    // that module never has to know how paths are configured, and so the button can simply be
    // disabled with the reason shown when it is missing.
    let xtb = run
        .settings
        .engine
        .needs_xtb()
        .then(|| resolve_xtb(&run.xtb_path));
    run.settings.xtb_binary = match &xtb {
        Some(Ok(path)) => Some(path.clone()),
        _ => None,
    };

    ui.horizontal(|ui| {
        let have_engine = !run.settings.engine.needs_xtb() || run.settings.xtb_binary.is_some();
        let can_run = !run.is_running() && !mol.atoms.is_empty() && have_engine;
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
    if let Some(Err(problem)) = &xtb {
        ui.label(egui::RichText::new(problem).small().weak());
    }

    if let Some(error) = &run.error {
        ui.add_space(6.0);
        ui.colored_label(egui::Color32::from_rgb(220, 120, 90), error);
    }

    trail_section(ui, run, mol);

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
        if outcome.conformers.len() == 1 { "" } else { "s" },
        outcome.elapsed.as_secs_f64(),
        outcome.workers,
        if outcome.workers == 1 { "" } else { "s" }
    ));
    ui.label(
        egui::RichText::new(format!(
            "Local optimiser: {}. Energies from {}.",
            outcome.local_optimizer.label(),
            if outcome.refined_with_xtb {
                "xTB GFN2"
            } else {
                "the DREIDING force field"
            }
        ))
        .small()
        .weak(),
    );
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
            "Populations are Boltzmann factors at 298 K, without vibrational entropy or symmetry.",
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

/// The xTB executable path, editable here as well as in the optimizer panel.
///
/// Both read and write the same configuration file, so a path set in either -- or in an earlier
/// session -- is found by the other. Offering it here matters because an external engine simply
/// does not run without it, and being told to go and set it somewhere else is a dead end.
fn xtb_path_row(ui: &mut egui::Ui, run: &mut ConformerRun) {
    use crate::qchem_interfaces::program::QcProgram;

    ui.horizontal(|ui| {
        ui.label("xTB path");
        let response = ui.add_enabled(
            !run.is_running(),
            egui::TextEdit::singleline(&mut run.xtb_path)
                .desired_width(200.0)
                .hint_text("xtb"),
        );
        if response.lost_focus() {
            crate::qchem_interfaces::config::save_path(QcProgram::Xtb, &run.xtb_path);
        }
        if ui
            .add_enabled(!run.is_running(), egui::Button::new("Browse…"))
            .clicked()
        {
            let mut dialog = rfd::FileDialog::new();
            if let Some(directory) = std::path::Path::new(&run.xtb_path)
                .parent()
                .filter(|p| p.is_dir())
            {
                dialog = dialog.set_directory(directory);
            }
            if let Some(path) = crate::recent_dir::pick_file(dialog) {
                run.xtb_path = path.display().to_string();
                crate::qchem_interfaces::config::save_path(QcProgram::Xtb, &run.xtb_path);
            }
        }
        // A tick as soon as it resolves, so there is no need to press Search to find out.
        if resolve_xtb(&run.xtb_path).is_ok() {
            ui.label(egui::RichText::new("found").small().weak());
        }
    });
}

/// Resolves an xTB path the way the optimizer panel does, so both agree on what counts as set.
fn resolve_xtb(path: &str) -> Result<std::path::PathBuf, String> {
    crate::qchem_interfaces::xtb_optimize::resolve_program_executable(
        crate::qchem_interfaces::program::QcProgram::Xtb,
        path,
    )
}

/// Transition-state conformers, by TStrail.
///
/// Its own section rather than part of the ordinary search, because it answers a different
/// question and costs differently: minutes of xTB rather than seconds of force field, and it needs
/// to be told which bonds are reacting.
fn trail_section(ui: &mut egui::Ui, run: &mut ConformerRun, mol: &Molecule) {
    ui.add_space(8.0);
    egui::CollapsingHeader::new("Transition-state conformers")
        .id_salt("conformer_tstrail")
        .default_open(false)
        .show(ui, |ui| {
            ui.label(
                egui::RichText::new(
                    "A transition state cannot simply be relaxed -- an optimiser slides it to \
                     reactant or product. This samples the conformers of the reactant instead, \
                     where the reacting atoms are still apart, then pushes each one over the \
                     barrier with an artificial force. Each pathway reaches its own transition \
                     state, so the results have a spread of active bond lengths rather than all \
                     sharing one.",
                )
                .small(),
            );

            ui.add_space(6.0);
            active_bond_editor(ui, run, mol);

            ui.add_space(6.0);
            ui.horizontal(|ui| {
                ui.label("Force");
                ui.add_enabled(
                    !run.is_trailing(),
                    egui::DragValue::new(&mut run.trail.force_kcal_per_angstrom)
                        .speed(5.0)
                        .range(10.0..=500.0)
                        .suffix(" kcal/mol/Å"),
                )
                .on_hover_text(
                    "How hard the reacting atoms are pulled together. Too weak and the barrier is \
                     never crossed; too strong and the structure is crushed on the way.",
                );
            });
            ui.horizontal(|ui| {
                ui.label("Steps per pathway");
                ui.add_enabled(
                    !run.is_trailing(),
                    egui::DragValue::new(&mut run.trail.cycles).range(5..=200),
                );
                ui.label("Pathways");
                ui.add_enabled(
                    !run.is_trailing(),
                    egui::DragValue::new(&mut run.trail.max_pathways).range(1..=100),
                )
                .on_hover_text("One per reactant conformer, so this caps the external cost.");
            });
            ui.horizontal(|ui| {
                ui.label("Charge");
                ui.add_enabled(
                    !run.is_trailing(),
                    egui::DragValue::new(&mut run.trail.charge).range(-10..=10),
                );
                ui.label("Multiplicity");
                ui.add_enabled(
                    !run.is_trailing(),
                    egui::DragValue::new(&mut run.trail.multiplicity).range(1..=10),
                );
            });

            xtb_path_row(ui, run);
            let resolved = resolve_xtb(&run.xtb_path);

            ui.add_space(6.0);
            ui.horizontal(|ui| {
                let ready = resolved.is_ok()
                    && !run.trail.active.is_empty()
                    && !mol.atoms.is_empty()
                    && !run.is_trailing()
                    && !run.is_running();
                if ui
                    .add_enabled(ready, egui::Button::new("Search transition-state conformers"))
                    .clicked()
                {
                    if let Ok(binary) = &resolved {
                        run.start_trail(mol, binary.clone());
                    }
                }
                if run.is_trailing() {
                    ui.spinner();
                    if let Some(elapsed) = run.trail_elapsed() {
                        ui.label(format!("{:.0} s", elapsed.as_secs_f64()));
                    }
                }
            });
            if run.trail.active.is_empty() {
                ui.label(
                    egui::RichText::new("Add at least one active bond first.")
                        .small()
                        .weak(),
                );
            }
            if let Err(problem) = &resolved {
                ui.label(
                    egui::RichText::new(format!(
                        "{problem} This needs xTB: a force field cannot break a bond."
                    ))
                    .small()
                    .weak(),
                );
            }
            if let Some(problem) = &run.trail_error {
                ui.colored_label(
                    egui::Color32::from_rgb(220, 120, 90),
                    egui::RichText::new(problem).small(),
                );
            }

            if let Some(outcome) = run.trail_outcome.clone() {
                trail_results(ui, &outcome);
            }
        });
}

/// Entering and listing the bonds that are forming or breaking.
fn active_bond_editor(ui: &mut egui::Ui, run: &mut ConformerRun, mol: &Molecule) {
    ui.label(egui::RichText::new("Active bonds").strong());
    ui.label(
        egui::RichText::new(
            "The atoms whose bond is forming or breaking, by the numbers shown on screen.",
        )
        .small()
        .weak(),
    );

    let count = mol.atoms.len();
    ui.horizontal(|ui| {
        ui.label("Atom");
        ui.add_enabled(
            !run.is_trailing(),
            egui::DragValue::new(&mut run.pending_active.0).range(1..=count.max(1)),
        );
        ui.label("to");
        ui.add_enabled(
            !run.is_trailing(),
            egui::DragValue::new(&mut run.pending_active.1).range(1..=count.max(1)),
        );

        // Show what those numbers mean, so a mistyped index is obvious before it is added.
        let name = |one_based: usize| match mol.atoms.get(one_based.wrapping_sub(1)) {
            Some(symbol) => format!("{symbol}{one_based}"),
            None => format!("#{one_based}"),
        };
        let (a, b) = run.pending_active;
        let sensible = a != b && a >= 1 && b >= 1 && a <= count && b <= count;
        if sensible {
            let distance = (mol.pos[a - 1] - mol.pos[b - 1]).length();
            ui.label(
                egui::RichText::new(format!(
                    "{} - {} ({:.2} Å)",
                    name(a),
                    name(b),
                    distance
                ))
                .small(),
            );
        }
        if ui
            .add_enabled(sensible && !run.is_trailing(), egui::Button::new("Add"))
            .clicked()
        {
            let bond = crate::conformer::afir::ActiveBond { a: a - 1, b: b - 1 };
            // Adding the same pair twice would double the force on it.
            let already = run
                .trail
                .active
                .iter()
                .any(|existing| {
                    (existing.a, existing.b) == (bond.a, bond.b)
                        || (existing.a, existing.b) == (bond.b, bond.a)
                });
            if !already {
                run.trail.active.push(bond);
            }
        }
    });

    let mut remove = None;
    for (index, bond) in run.trail.active.iter().enumerate() {
        ui.horizontal(|ui| {
            let name = |atom: usize| match mol.atoms.get(atom) {
                Some(symbol) => format!("{symbol}{}", atom + 1),
                None => format!("#{}", atom + 1),
            };
            let distance = if bond.a < mol.pos.len() && bond.b < mol.pos.len() {
                format!(" ({:.2} Å)", (mol.pos[bond.a] - mol.pos[bond.b]).length())
            } else {
                String::new()
            };
            ui.label(
                egui::RichText::new(format!("{} - {}{distance}", name(bond.a), name(bond.b)))
                    .small(),
            );
            if ui.small_button("Remove").clicked() {
                remove = Some(index);
            }
        });
    }
    if let Some(index) = remove {
        run.trail.active.remove(index);
    }
}

/// What a pathway search found.
fn trail_results(ui: &mut egui::Ui, outcome: &crate::conformer::tstrail::TrailOutcome) {
    ui.add_space(8.0);
    ui.separator();
    ui.add_space(6.0);
    ui.label(egui::RichText::new("Transition-state guesses").strong());
    ui.label(format!(
        "{} guess{} from {} pathway{} ({} reactant conformers, {} discarded as broken), {} xTB          calls.",
        outcome.guesses.len(),
        if outcome.guesses.len() == 1 { "" } else { "es" },
        outcome.pathways,
        if outcome.pathways == 1 { "" } else { "s" },
        outcome.reactant_conformers,
        outcome.discarded_broken,
        outcome.xtb_calls
    ));
    if outcome.failed > 0 {
        ui.label(
            egui::RichText::new(format!(
                "{} pathway{} failed{}",
                outcome.failed,
                if outcome.failed == 1 { "" } else { "s" },
                match &outcome.first_failure {
                    Some(problem) => format!(": {problem}"),
                    None => ".".into(),
                }
            ))
            .small()
            .weak(),
        );
    }

    // The spread is the point of the method, so it is stated rather than left to be read off the
    // table: a single value means a fixed-length search would have found the same thing.
    for (which, (lo, hi)) in outcome.active_length_range().iter().enumerate() {
        ui.label(
            egui::RichText::new(format!(
                "Active bond {} spans {lo:.2} to {hi:.2} Å across the guesses.",
                which + 1
            ))
            .small(),
        );
    }

    ui.add_space(6.0);
    egui::ScrollArea::vertical()
        .id_salt("tstrail_table")
        .max_height(220.0)
        .show(ui, |ui| {
            egui::Grid::new("tstrail_grid")
                .num_columns(4)
                .striped(true)
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("#").strong());
                    ui.label(egui::RichText::new("ΔE (kcal/mol)").strong());
                    ui.label(egui::RichText::new("Active bonds (Å)").strong());
                    ui.label(egui::RichText::new("Crossed").strong());
                    ui.end_row();
                    for (index, guess) in outcome.guesses.iter().enumerate() {
                        ui.label(format!("{}", index + 1));
                        ui.label(format!("{:.2}", guess.relative_energy_kcal));
                        ui.label(
                            guess
                                .active_lengths
                                .iter()
                                .map(|l| format!("{l:.2}"))
                                .collect::<Vec<_>>()
                                .join(", "),
                        );
                        // A pathway still climbing when it stopped is a weaker starting point,
                        // and saying so is more use than hiding it.
                        if guess.crossed {
                            ui.label("yes");
                        } else {
                            ui.label(egui::RichText::new("no").weak());
                        }
                        ui.end_row();
                    }
                });
        });

    ui.add_space(4.0);
    ui.label(
        egui::RichText::new(
            "These are transition-state guesses, not transition states: the pathway maximum is \
             not refined and no saddle-point optimisation or IRC has been run. Take them to a \
             transition-state optimisation.",
        )
        .small()
        .weak(),
    );

    ui.add_space(6.0);
    ui.label(egui::RichText::new("Please cite").strong());
    ui.label(
        egui::RichText::new(
            "This method is not ours. If you publish transition states found with it, cite:",
        )
        .small(),
    );
    ui.label(
        egui::RichText::new(crate::conformer::tstrail::CITATION)
            .small()
            .italics(),
    );
    if ui
        .small_button("Copy")
        .on_hover_text("Copies the reference to the clipboard.")
        .clicked()
    {
        ui.ctx()
            .copy_text(crate::conformer::tstrail::CITATION.to_string());
    }
}

/// The GFN2 re-ranking: the button, and what the last one did.
///
/// Offered on the results rather than as a search setting, because it is a different kind of
/// wait. The search is seconds; this is minutes of process launches, so it is started
/// deliberately and the results already on screen stay readable while it runs.
fn refinement_row(ui: &mut egui::Ui, run: &mut ConformerRun, outcome: &ConformerOutcome) {
    ui.add_space(8.0);
    ui.separator();
    ui.add_space(6.0);

    let resolved = resolve_xtb(&run.xtb_path);
    let available = outcome.conformers.len();
    ui.horizontal(|ui| {
        ui.label("Re-rank the lowest");
        ui.add_enabled(
            !run.is_refining() && !outcome.refined_with_xtb,
            egui::DragValue::new(&mut run.refine_top).range(1..=available.max(1)),
        )
        .on_hover_text(
            "How many conformers to re-optimise. The cost is roughly linear in this, so it is a \
             judgement about how much the ordering matters -- refine all of them if it matters a \
             lot.",
        );
        ui.label(format!("of {available}"));
        if ui
            .add_enabled(
                !run.is_refining() && !outcome.refined_with_xtb,
                egui::Button::new("All"),
            )
            .clicked()
        {
            run.refine_top = available.max(1);
        }
    });

    ui.horizontal(|ui| {
        let count = available.min(run.refine_top.max(1));
        let can_refine = resolved.is_ok() && !run.is_refining() && !outcome.refined_with_xtb;
        if ui
            .add_enabled(
                can_refine,
                egui::Button::new(format!("Re-rank top {count} with xTB GFN2")),
            )
            .on_hover_text(
                "Re-optimises the leading conformers with GFN2 and re-orders them. Relative \
                 conformer energies are where a generic force field is weakest, so this is \
                 usually worth the wait.",
            )
            .clicked()
        {
            if let Ok(binary) = &resolved {
                run.start_refinement(binary.clone());
            }
        }
        if run.is_refining() {
            ui.spinner();
            ui.label("running xTB…");
        }
    });

    if let Err(problem) = &resolved {
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new(problem).small().weak());
        });
        xtb_path_row(ui, run);
    }

    match &run.refinement {
        Some(Ok(report)) => {
            ui.label(
                egui::RichText::new(format!(
                    "Re-ranked {} conformer{} with {} xTB calls; {} changed place{}.",
                    report.refined,
                    if report.refined == 1 { "" } else { "s" },
                    report.xtb_calls,
                    report.reordered,
                    if report.reordered == 1 { "" } else { "s" }
                ))
                .small(),
            );
        }
        Some(Err(problem)) => {
            ui.colored_label(
                egui::Color32::from_rgb(220, 120, 90),
                egui::RichText::new(format!("xTB re-ranking failed: {problem}")).small(),
            );
        }
        None => {}
    }

    if let Some(Ok(report)) = &run.refinement {
        if outcome.refined_with_xtb && report.refined < outcome.conformers.len() {
            // The list is honestly mixed, so say so rather than let the numbers look comparable.
            ui.label(
                egui::RichText::new(format!(
                    "The first {} energies are GFN2; the remaining {} are still DREIDING and are \
                     only placed after them.",
                    report.refined,
                    outcome.conformers.len() - report.refined
                ))
                .small()
                .weak(),
            );
        }
    }
}

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
            rings: 0,
            generations: 4,
            structures_tried: 40,
            local_optimizations: 30,
            termination: "Converged".into(),
            elapsed: Duration::from_millis(1500),
            workers: 4,
            local_optimizer: LocalOptimizer::Cartesian,
            constrained: 0,
            best_found_in_generation: 3,
            radical_warning: None,
            refined_with_xtb: false,
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
