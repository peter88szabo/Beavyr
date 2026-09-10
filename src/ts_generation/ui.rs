//! The TS Generation tool, shared by the floating and classic layouts.

use std::sync::atomic::Ordering;

use bevy::prelude::Vec3;
use bevy_egui::egui;

use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
use crate::molecule::Molecule;
use crate::optimizer::rda::{RdaConditionalCoordinates, RdaDistanceMode};
use crate::qchem_interfaces::{method_ui, program::QcProgram};
use crate::trajectory::{TrajectoryFrame, TrajectoryState};

use super::parameters::{Algorithm, Parameters, ReactionKind};
use super::{backend, show_frames, Endpoint, TsGeneration};

fn endpoint_row(
    ui: &mut egui::Ui,
    slot: &str,
    endpoint: &mut Option<Endpoint>,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
    message: &mut Option<String>,
    is_error: &mut bool,
) -> bool {
    let mut displayed = false;
    ui.group(|ui| {
        ui.horizontal(|ui| {
            ui.strong(if slot == "A" { "A — Reactant" } else { "B — Product" });
            if endpoint.is_some() {
                ui.colored_label(
                    egui::Color32::from_rgb(75, 190, 120),
                    format!("● {slot} loaded"),
                );
            } else {
                ui.weak(format!("○ {slot} not loaded"));
            }
        });
        if let Some(endpoint) = endpoint.as_ref() {
            ui.label(format!("{} — {} atoms", endpoint.name, endpoint.atoms.len()));
        }
        ui.horizontal_wrapped(|ui| {
            if ui
                .add_enabled(
                    !mol.atoms.is_empty() && !traj.playing,
                    egui::Button::new(format!("Use current as {slot}")),
                )
                .on_hover_text("Store a copy of the displayed structure, including any drawing or optimization changes. Pause trajectory playback first.")
                .clicked()
            {
                match Endpoint::capture(mol, &format!("{slot} from viewer")) {
                    Ok(captured) => {
                        *endpoint = Some(captured);
                        *message = None;
                        *is_error = false;
                    }
                    Err(error) => {
                        *message = Some(error.to_string());
                        *is_error = true;
                    }
                }
            }
            if ui.button(format!("Load {slot} XYZ…")).clicked() {
                let dialog = crate::recent_dir::open().add_filter("XYZ structure", &["xyz"]);
                if let Some(path) = crate::recent_dir::pick_file(dialog) {
                    match Endpoint::from_path(&path) {
                        Ok(loaded) => {
                            *endpoint = Some(loaded);
                            *message = None;
                            *is_error = false;
                        }
                        Err(error) => {
                            *message = Some(error.to_string());
                            *is_error = true;
                        }
                    }
                }
            }
            if ui
                .add_enabled(endpoint.is_some(), egui::Button::new(format!("Show {slot}")))
                .on_hover_text("Display the stored endpoint in the viewer for inspection, editing or optimization.")
                .clicked()
            {
                if let Some(endpoint) = endpoint.as_ref() {
                    displayed = show_frames(
                        vec![TrajectoryFrame {
                            atoms: endpoint.atoms.clone(),
                            pos: endpoint.positions.clone(),
                        }],
                        mol,
                        traj,
                    );
                }
            }
        });
    });
    displayed
}

fn text_row(ui: &mut egui::Ui, label: &str, text: &mut String, hint: &str) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(
            egui::TextEdit::singleline(text)
                .desired_width(260.0)
                .hint_text(hint),
        );
    });
}

fn reaction_controls(ui: &mut egui::Ui, parameters: &mut Parameters) {
    ui.separator();
    ui.strong("Reaction definition");
    egui::ComboBox::from_id_salt("ts_reaction_kind")
        .selected_text(parameters.reaction.label())
        .show_ui(ui, |ui| {
            for kind in ReactionKind::ALL {
                if ui
                    .selectable_value(&mut parameters.reaction, kind, kind.label())
                    .changed()
                {
                    parameters.indices = match kind {
                        ReactionKind::Bond => "1 2",
                        ReactionKind::Transfer | ReactionKind::Angle => "1 2 3",
                        ReactionKind::Dihedral => "1 2 3 4",
                        _ => "",
                    }
                    .into();
                    parameters.weights.clear();
                    parameters.options.distance_mode = if kind == ReactionKind::ActiveAtoms {
                        RdaDistanceMode::Cartesian
                    } else {
                        RdaDistanceMode::ReactiveInternal
                    };
                }
            }
        });
    ui.label(parameters.reaction.hint());
    ui.weak("Atom numbers start at 1 and must match in both endpoint structures.");
    match parameters.reaction {
        ReactionKind::ActiveAtoms => {}
        ReactionKind::Custom => {
            text_row(ui, "Bonds", &mut parameters.bonds, "1 2; 2 3");
            text_row(ui, "Angles", &mut parameters.angles, "1 2 3");
            text_row(ui, "Dihedrals", &mut parameters.dihedrals, "1 2 3 4");
        }
        _ => text_row(ui, "Atom numbers", &mut parameters.indices, "1 2"),
    }
    text_row(
        ui,
        "Active atoms",
        &mut parameters.active_atoms,
        "blank = all atoms",
    );
    ui.small("For Cartesian distance, active atoms select which displacements are measured. They also define the default alignment subset. All atoms remain free to relax.");
    if parameters.reaction != ReactionKind::ActiveAtoms {
        egui::ComboBox::from_id_salt("ts_distance_metric")
            .selected_text(match parameters.options.distance_mode {
                RdaDistanceMode::Cartesian => "Cartesian distance on active atoms",
                RdaDistanceMode::ReactiveInternal => "Distance in reactive coordinates",
            })
            .show_ui(ui, |ui| {
                ui.selectable_value(
                    &mut parameters.options.distance_mode,
                    RdaDistanceMode::Cartesian,
                    "Cartesian distance on active atoms",
                );
                ui.selectable_value(
                    &mut parameters.options.distance_mode,
                    RdaDistanceMode::ReactiveInternal,
                    "Distance in reactive coordinates",
                );
            });
        if parameters.options.distance_mode == RdaDistanceMode::ReactiveInternal {
            text_row(
                ui,
                "Coordinate weights",
                &mut parameters.weights,
                "blank = equal; bonds, angles, dihedrals",
            );
        }
    }
    ui.checkbox(
        &mut parameters.options.align,
        "Align product to reactant before generating the path",
    );
    if parameters.options.align {
        text_row(
            ui,
            "Alignment atoms",
            &mut parameters.align_atoms,
            "blank = active atoms (or all)",
        );
    }
}

fn scalar_row(
    ui: &mut egui::Ui,
    label: &str,
    value: &mut f64,
    speed: f64,
    range: std::ops::RangeInclusive<f64>,
) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(value).speed(speed).range(range));
    });
}

fn count_row(ui: &mut egui::Ui, label: &str, value: &mut usize, min: usize) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(value).range(min..=10000));
    });
}

fn algorithm_controls(ui: &mut egui::Ui, parameters: &mut Parameters) {
    let neb = parameters.algorithm == Algorithm::PoorMansNeb;
    if neb {
        ui.separator();
        ui.strong("Poor Man's NEB images");
        count_row(
            ui,
            "Images (including both endpoints)",
            &mut parameters.options.poormans_neb_images,
            3,
        );
        count_row(
            ui,
            "Maximum relaxation steps per image",
            &mut parameters.options.poormans_neb_maxiter,
            1,
        );
        ui.small("Reactive bond lengths, angles and dihedrals are interpolated between the endpoints and constrained while each interior image relaxes.");
    }
    egui::CollapsingHeader::new(if neb { "Relaxation controls" } else { "RDA controls" })
        .id_salt("ts_algorithm_controls").default_open(false).show(ui, |ui| {
            let opts = &mut parameters.options;
            if !neb {
                egui::ComboBox::from_id_salt("ts_conditional_coordinates")
                    .selected_text(match opts.conditional_coordinates {
                        RdaConditionalCoordinates::Internal => "Conditional relaxation: internal coordinates",
                        RdaConditionalCoordinates::Cartesian => "Conditional relaxation: Cartesian coordinates",
                    }).show_ui(ui, |ui| {
                        ui.selectable_value(&mut opts.conditional_coordinates, RdaConditionalCoordinates::Internal, "Internal coordinates");
                        ui.selectable_value(&mut opts.conditional_coordinates, RdaConditionalCoordinates::Cartesian, "Cartesian coordinates");
                    });
                count_row(ui, "Maximum conditional steps", &mut opts.max_conditional_steps, 1);
                scalar_row(ui, "First energy tolerance (eV)", &mut opts.first_energy_tol_ev, 0.001, 1.0e-8..=10.0);
                scalar_row(ui, "Starting beta", &mut opts.beta_start, 0.01, 0.0..=1.0);
                scalar_row(ui, "Beta increment", &mut opts.beta_step, 0.01, 0.001..=1.0);
                count_row(ui, "Maximum bracket attempts", &mut opts.max_bracket_attempts, 1);
                text_row(ui, "Gamma values", &mut parameters.gamma_values, "0.1 0.2 … 0.9");
                scalar_row(ui, "Distance tolerance (bohr / rad)", &mut opts.distance_tol, 0.005, 1.0e-8..=10.0);
                ui.small("The distance metric uses bohr for Cartesian coordinates and bonds, radians for angles and dihedrals.");
                ui.checkbox(&mut opts.fallback_to_cartesian, "Allow Cartesian relaxation if internal coordinates fail");
                scalar_row(ui, "Cartesian step limit (bohr)", &mut opts.cartesian_step_max_component, 0.01, 1.0e-6..=10.0);
            }
            scalar_row(ui, "Later energy tolerance (eV)", &mut opts.later_energy_tol_ev, 0.001, 1.0e-8..=10.0);
            scalar_row(ui, "Internal step limit (bohr / rad)", &mut opts.internal_step_max_component, 0.01, 1.0e-6..=10.0);
            count_row(ui, "Internal-coordinate fit iterations", &mut opts.internal_best_fit_iters, 1);
            if neb {
                ui.small("Additional connectivity for the NEB internal-coordinate set:");
                text_row(ui, "Extra bonds", &mut parameters.extra_bonds, "1 2; 3 4");
                text_row(ui, "Extra angles", &mut parameters.extra_angles, "1 2 3");
                text_row(ui, "Extra dihedrals", &mut parameters.extra_dihedrals, "1 2 3 4");
            }
        });
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
    if !ui.button(label).clicked() {
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

/// Returns true only when the user explicitly changes the displayed structure.
pub fn panel(
    ui: &mut egui::Ui,
    state: &mut TsGeneration,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
) -> bool {
    let mut displayed = false;
    let running = state.is_running();
    ui.label("Generate a TS guess between reactant A and product B using RDA or Poor Man's NEB.");
    if ui
        .add(
            egui::Button::new(
                egui::RichText::new("How To Use")
                    .strong()
                    .color(egui::Color32::WHITE),
            )
            .fill(egui::Color32::from_rgb(35, 105, 170))
            .stroke(egui::Stroke::new(
                1.0,
                egui::Color32::from_rgb(95, 178, 240),
            ))
            .min_size(egui::vec2(140.0, 32.0)),
        )
        .on_hover_text("Guide to endpoints, Reaction definition, and RDA / NEB controls.")
        .clicked()
    {
        state.help_open = true;
    }
    ui.small("Draw or optimize A, capture it, then do the same for B. Captures stay stored while you edit the viewer. Use Show A/B to display a stored endpoint; capture again to keep subsequent edits.");
    ui.add_enabled_ui(!running, |ui| {
        ui.horizontal(|ui| {
            ui.selectable_value(&mut state.parameters.algorithm, Algorithm::Rda, "RDA");
            ui.selectable_value(
                &mut state.parameters.algorithm,
                Algorithm::PoorMansNeb,
                "Poor Man's NEB",
            );
        });
        displayed |= endpoint_row(
            ui,
            "A",
            &mut state.reactant,
            mol,
            traj,
            &mut state.message,
            &mut state.is_error,
        );
        displayed |= endpoint_row(
            ui,
            "B",
            &mut state.product,
            mol,
            traj,
            &mut state.message,
            &mut state.is_error,
        );
        if ui
            .add_enabled(
                state.reactant.is_some() && state.product.is_some(),
                egui::Button::new("Swap reactant and product"),
            )
            .clicked()
        {
            std::mem::swap(&mut state.reactant, &mut state.product);
        }
        reaction_controls(ui, &mut state.parameters);
        algorithm_controls(ui, &mut state.parameters);
        ui.separator();
        ui.strong("Energy and gradient calculation");
        ui.horizontal(|ui| {
            egui::ComboBox::from_id_salt("ts_program")
                .selected_text(state.program.label())
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut state.program, QcProgram::Xtb, "xTB");
                    ui.selectable_value(&mut state.program, QcProgram::Behemoth, "Behemoth");
                });
            ui.label("Charge");
            ui.add(egui::DragValue::new(&mut state.charge).range(-100..=100));
            ui.label("Multiplicity");
            ui.add(egui::DragValue::new(&mut state.multiplicity).range(1..=100));
        });
        let executable = if state.program == QcProgram::Xtb {
            &mut state.xtb_path
        } else {
            &mut state.behemoth_path
        };
        ui.horizontal(|ui| {
            ui.label("Executable");
            ui.add(egui::TextEdit::singleline(executable).desired_width(270.0));
            if ui.button("Browse…").clicked() {
                if let Some(path) = crate::recent_dir::pick_file(crate::recent_dir::open()) {
                    *executable = path.to_string_lossy().into_owned();
                }
            }
        });
        let atoms = state
            .reactant
            .as_ref()
            .map(|e| e.atoms.as_slice())
            .unwrap_or(&[]);
        method_ui::method_config_ui(
            ui,
            state.program,
            &mut state.method,
            state.multiplicity,
            atoms,
            running,
            "ts_generation",
        );
    });

    ui.separator();
    if running {
        ui.horizontal(|ui| {
            ui.spinner();
            let evaluations = state.control.evaluations.load(Ordering::Relaxed);
            let elapsed = state
                .elapsed()
                .map(crate::qchem_interfaces::xtb_optimize::format_elapsed)
                .unwrap_or_default();
            ui.label(format!("{evaluations} energy/gradient calls · {elapsed}"));
            let cancelling = state.control.cancelled.load(Ordering::Relaxed);
            if ui
                .add_enabled(
                    !cancelling,
                    egui::Button::new(if cancelling {
                        "Cancelling…"
                    } else {
                        "Cancel"
                    }),
                )
                .clicked()
            {
                state.control.cancel();
            }
        });
        if let Ok(energy) = state.control.latest_energy.lock() {
            if let Some(energy) = *energy {
                ui.label(format!("Last energy: {energy:.9} Eh"));
            }
        }
        ui.ctx()
            .request_repaint_after(std::time::Duration::from_millis(100));
    } else {
        let problem = state.input_options().err().map(|e| e.to_string());
        if let Some(problem) = &problem {
            ui.label(problem);
        }
        let label = if state.parameters.algorithm == Algorithm::PoorMansNeb {
            "Generate NEB path and TS guess"
        } else {
            "Generate RDA TS guess"
        };
        if ui
            .add_enabled(problem.is_none(), egui::Button::new(label))
            .clicked()
        {
            if let Err(error) = state.start() {
                state.message = Some(error.to_string());
                state.is_error = true;
            }
        }
    }
    if let Some(message) = &state.message {
        if state.is_error {
            ui.colored_label(egui::Color32::from_rgb(225, 110, 90), message);
        } else {
            ui.label(message);
        }
    }
    if let Some(output) = &state.output {
        ui.separator();
        ui.strong(format!(
            "{} result · {}",
            output.algorithm.label(),
            output.method
        ));
        ui.label(format!("TS guess energy: {:.9} Eh", output.energy_hartree));
        ui.label(&output.guess.message);
        ui.small("This is a TS guess. Optimize the saddle point and check its imaginary mode before treating it as a confirmed transition state.");
        ui.horizontal_wrapped(|ui| {
            if ui.button("Show TS guess").clicked() {
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
                    vec![TrajectoryFrame {
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
            if ui.button("Show path in Trajectory").clicked() {
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
                if ui.button("Plot NEB energy").clicked() {
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
