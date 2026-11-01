//! Existing Behemoth optimizer controls, shared with the calculator.
use super::{run::QcRunTask, QcPanelState};
use crate::optimizer::transition_state::*;
use crate::qchem_interfaces::{
    job::{search::*, JobType},
    program::QcProgram,
};
use bevy_egui::egui;
fn number(ui: &mut egui::Ui, label: &str, value: &mut f64) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(value).speed(0.0001).max_decimals(10));
    });
}
fn count(ui: &mut egui::Ui, label: &str, value: &mut usize) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(value));
    });
}
fn optional(ui: &mut egui::Ui, label: &str, value: &mut Option<f64>, default: f64) {
    ui.horizontal(|ui| {
        let mut set = value.is_some();
        if ui.checkbox(&mut set, label).changed() {
            *value = if set { Some(default) } else { None };
        }
        if let Some(v) = value {
            ui.add(egui::DragValue::new(v).speed(0.001));
        }
    });
}
pub(super) fn settings(ui: &mut egui::Ui, panel: &mut QcPanelState, task: &QcRunTask) {
    let ts = panel.job == JobType::TransitionState;
    crate::ui_style::group(ui, if ts { "TS search" } else { "Minimum search" }, |ui| {
        ui.add_enabled_ui(!task.is_running(), |ui| {
            let program = panel.program;
            let o = &mut panel.config.search;
            if matches!(
                program,
                QcProgram::Orca | QcProgram::Psi4 | QcProgram::Psi4Py
            ) || (program == QcProgram::PySCF && !ts)
                || (program == QcProgram::Xtb && !ts)
            {
                ui.horizontal_wrapped(|ui| {
                    let native = match program {
                        QcProgram::Orca => "ORCA built-in",
                        QcProgram::Psi4 | QcProgram::Psi4Py => "Psi4 built-in",
                        QcProgram::PySCF => "PySCF / geomeTRIC",
                        _ => "xTB built-in",
                    };
                    ui.selectable_value(&mut o.driver, Driver::Native, native);
                    ui.selectable_value(&mut o.driver, Driver::Beavyr, "Beavyr");
                });
            }
            let is_local = local(program, panel.job, o);
            let behemoth = program == QcProgram::Behemoth;
            if ts && program == QcProgram::Orca && !is_local {
                count(
                    ui,
                    "Recalculate Hessian every N steps",
                    &mut o.orca_recalc_hessian,
                );
                return;
            }
            if !is_local && !(behemoth && ts) {
                return;
            }
            if !ts {
                minimum(ui, &mut o.minimum);
                return;
            }
            let t = &mut o.ts;
            egui::ComboBox::from_id_salt("ts_coordinates")
                .selected_text(match t.coordinate_mode {
                    TransitionStateCoordinateMode::Cartesian => "Cartesian",
                    TransitionStateCoordinateMode::Internal => "Primitive internal",
                    TransitionStateCoordinateMode::PulayRedundant => "Pulay redundant internal",
                })
                .show_ui(ui, |ui| {
                    for (v, l) in [
                        (TransitionStateCoordinateMode::Cartesian, "Cartesian"),
                        (
                            TransitionStateCoordinateMode::Internal,
                            "Primitive internal",
                        ),
                        (
                            TransitionStateCoordinateMode::PulayRedundant,
                            "Pulay redundant internal",
                        ),
                    ] {
                        ui.selectable_value(&mut t.coordinate_mode, v, l);
                    }
                });
            count(ui, "Max steps", &mut t.max_cycles);
            count(
                ui,
                "Recalculate Hessian every N steps",
                &mut t.hessian_recalc_interval,
            );
            ui.horizontal(|ui| {
                ui.label("Initial Hessian");
                ui.selectable_value(
                    &mut t.ts_initial_hessian,
                    TransitionStateInitialHessian::Exact,
                    "Calculated",
                );
                ui.selectable_value(
                    &mut t.ts_initial_hessian,
                    TransitionStateInitialHessian::Model,
                    "Model",
                );
            });
            ui.horizontal(|ui| {
                ui.label("Mode tracking");
                ui.selectable_value(
                    &mut t.mode_tracking_coordinates,
                    TransitionStateModeTracking::Internal,
                    "Internal",
                );
                ui.selectable_value(
                    &mut t.mode_tracking_coordinates,
                    TransitionStateModeTracking::Cartesian,
                    "Cartesian",
                );
            });
            egui::ComboBox::from_id_salt("ts_reference")
                .selected_text(format!("{:?}", o.reaction))
                .show_ui(ui, |ui| {
                    for (r, l) in [
                        (Reaction::Lowest, "Lowest mode"),
                        (Reaction::Bond, "Bond"),
                        (Reaction::Angle, "Angle"),
                        (Reaction::Dihedral, "Dihedral"),
                        (Reaction::Transfer, "Atom transfer"),
                        (Reaction::Mixed, "Weighted coordinates"),
                        (Reaction::Direction, "Cartesian direction"),
                    ] {
                        if !behemoth || r != Reaction::Direction {
                            ui.selectable_value(&mut o.reaction, r, l);
                        }
                    }
                });
            match o.reaction {
                Reaction::Lowest => {}
                Reaction::Mixed => {
                    ui.add(
                        egui::TextEdit::singleline(&mut o.reaction_coordinates)
                            .hint_text("bond:1,2:1;bond:2,3:-1"),
                    );
                }
                Reaction::Direction => {
                    ui.add(
                        egui::TextEdit::multiline(&mut o.direction)
                            .desired_rows(2)
                            .hint_text("Three Cartesian components per atom"),
                    );
                }
                _ => {
                    ui.add(
                        egui::TextEdit::singleline(&mut o.reaction_atoms)
                            .hint_text("Atom numbers (starting at 1)"),
                    );
                }
            }
            egui::CollapsingHeader::new("Additional coordinates")
                .id_salt("ts_extra_coordinates")
                .show(ui, |ui| {
                    for (label, text) in [
                        ("Bonds", &mut o.extra_bonds),
                        ("Angles", &mut o.extra_angles),
                        ("Dihedrals", &mut o.extra_dihedrals),
                    ] {
                        ui.horizontal(|ui| {
                            ui.label(label);
                            ui.add(egui::TextEdit::singleline(text).hint_text("1,2;3,4"));
                        });
                    }
                });
            egui::CollapsingHeader::new("Advanced TS settings")
                .id_salt("ts_advanced")
                .show(ui, |ui| {
                    let t = &mut o.ts;
                    optional(ui, "Initial trust radius", &mut t.trust_radius, 0.2);
                    optional(ui, "Maximum trust radius", &mut t.trust_radius_max, 0.5);
                    number(ui, "Minimum trust radius", &mut t.trust_radius_min);
                    number(
                        ui,
                        "Model negative curvature",
                        &mut t.ts_model_negative_curvature,
                    );
                    count(ui, "Maximum rejected steps", &mut t.max_rejected_steps);
                    ui.checkbox(&mut t.project_eckart, "Project rotations and translations");
                    if behemoth {
                        return;
                    }
                    ui.checkbox(&mut t.final_hessian, "Calculate final Hessian");
                    number(ui, "Hessian difference step (bohr)", &mut t.hessian_step);
                    for (l, v) in [
                        ("Maximum internal step", &mut t.max_step_internal),
                        (
                            "Minimum gradient improvement",
                            &mut t.min_gradient_improvement,
                        ),
                        ("Metric cutoff", &mut t.g_cutoff),
                        ("Redundant-coordinate penalty", &mut t.redundant_penalty),
                    ] {
                        number(ui, l, v);
                    }
                    count(ui, "Coordinate fitting iterations", &mut t.best_fit_iters);
                    ui.checkbox(&mut t.repair_ts_hessian, "Repair TS model Hessian");
                    number(
                        ui,
                        "Hessian eigenvalue floor",
                        &mut t.ts_hessian_eigenvalue_floor,
                    );
                    ui.checkbox(
                        &mut t.adaptive_hessian_recalc,
                        "Adaptive Hessian recalculation",
                    );
                    for (l, v) in [
                        ("Minimum mode overlap", &mut t.adaptive_hessian_overlap_min),
                        ("Minimum step agreement", &mut t.adaptive_hessian_rho_min),
                        ("Maximum step agreement", &mut t.adaptive_hessian_rho_max),
                        (
                            "Trust-radius fraction",
                            &mut t.adaptive_hessian_trust_fraction,
                        ),
                        (
                            "Near-convergence factor",
                            &mut t.hessian_recalc_near_convergence_factor,
                        ),
                    ] {
                        number(ui, l, v);
                    }
                    count(
                        ui,
                        "Adaptive recalculation cooldown",
                        &mut t.adaptive_hessian_recalc_cooldown,
                    );
                    ui.checkbox(
                        &mut t.adaptive_hessian_on_negative_mode_change,
                        "Recalculate when negative-mode count changes",
                    );
                    ui.checkbox(
                        &mut t.adaptive_hessian_on_trust_collapse,
                        "Recalculate on trust-radius collapse",
                    );
                    ui.checkbox(
                        &mut t.adaptive_hessian_retry_recalc,
                        "Recalculate after rejected steps",
                    );
                    ui.checkbox(
                        &mut t.skip_hessian_recalc_near_convergence,
                        "Skip scheduled Hessian near convergence",
                    );
                });
        });
    });
}
fn minimum(ui: &mut egui::Ui, m: &mut Minimum) {
    egui::ComboBox::from_id_salt("minimum_coordinates")
        .selected_text(m.coordinates.label())
        .show_ui(ui, |ui| {
            for c in [
                Coordinates::Cartesian,
                Coordinates::Primitive,
                Coordinates::Pulay,
            ] {
                ui.selectable_value(&mut m.coordinates, c, c.label());
            }
        });
    egui::ComboBox::from_id_salt("minimum_algorithm")
        .selected_text(format!("{:?}", m.algorithm))
        .show_ui(ui, |ui| {
            for (a, l) in [
                (MinimumAlgorithm::Bfgs, "BFGS"),
                (MinimumAlgorithm::SteepestDescent, "Steepest descent"),
                (MinimumAlgorithm::Bb1, "Two-point gradient 1"),
                (MinimumAlgorithm::Bb2, "Two-point gradient 2"),
            ] {
                if m.coordinates == Coordinates::Cartesian
                    || matches!(
                        a,
                        MinimumAlgorithm::Bfgs | MinimumAlgorithm::SteepestDescent
                    )
                {
                    ui.selectable_value(&mut m.algorithm, a, l);
                }
            }
        });
    count(ui, "Max steps", &mut m.max_cycles);
    egui::CollapsingHeader::new("Advanced minimum settings").show(ui, |ui| {
        if m.coordinates != Coordinates::Cartesian {
            number(ui, "Maximum internal step", &mut m.max_step_internal);
            count(ui, "Coordinate fitting iterations", &mut m.best_fit_iters);
        }
        if !(m.coordinates == Coordinates::Primitive && m.algorithm == MinimumAlgorithm::Bfgs) {
            number(
                ui,
                "Minimum line-search scale",
                &mut m.min_line_search_scale,
            );
            if m.coordinates == Coordinates::Pulay {
                number(ui, "Metric cutoff", &mut m.g_cutoff);
                number(ui, "Redundant-coordinate penalty", &mut m.redundant_penalty);
            }
        }
    });
}
