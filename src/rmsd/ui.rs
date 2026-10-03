//! The Structure Comparison panel: pin a reference frame, compare the
//! current frame against it, and optionally superpose.

use crate::ui_style::ResponseExt;
use bevy::prelude::*;
use bevy_egui::egui;

use crate::molecule::Molecule;
use crate::trajectory::TrajectoryState;

use super::kabsch::{compare, same_atom_ordering, superpose};
use super::RmsdState;

/// Draws the panel. Called from `ui.rs`, matching how the other feature
/// panels are free functions the main panel invokes.
pub fn rmsd_panel(
    ui: &mut egui::Ui,
    state: &mut RmsdState,
    traj: &mut TrajectoryState,
    mol: &mut Molecule,
) {
    ui.scope(|ui| {
        crate::ui_style::modern(ui);
        crate::ui_style::heading(
            ui,
            "Compare molecular structures",
            "Pin a reference, measure differences and align frames.",
        );
        panel_body(ui, state, traj, mol);
    });
}

fn panel_body(
    ui: &mut egui::Ui,
    state: &mut RmsdState,
    traj: &mut TrajectoryState,
    mol: &mut Molecule,
) {
    if traj.frames.is_empty() {
        ui.weak("Load or generate a trajectory to compare frames.");
        return;
    }

    let total = traj.frames.len();
    let current = traj.current_frame.min(total - 1);
    let reference_index = state.effective_reference();

    crate::ui_style::group(ui, "Reference structure", |ui| {
        ui.label(format!(
            "Reference: frame {} of {}{}",
            reference_index + 1,
            total,
            if state.reference_frame.is_none() {
                "  (default: first frame)"
            } else {
                "  (pinned)"
            }
        ));

        ui.horizontal_wrapped(|ui| {
            if ui.button("Use current frame as reference").clicked_once() {
                pin_reference(state, traj, current);
            }
            if state.reference_frame.is_some() && ui.button("Reset to first frame").clicked_once() {
                state.reference_frame = None;
                state.reference_atoms.clear();
                state.reference_pos.clear();
                state.last_report = None;
                state.last_message = Some("Reference reset to the first frame.".to_string());
                state.last_error = None;
            }
        });
    });

    // Whatever is pinned, or the first frame if nothing is.
    let (reference_atoms, reference_pos) = if state.reference_pos.is_empty() {
        let frame = &traj.frames[reference_index.min(total - 1)];
        (frame.atoms.clone(), frame.pos.clone())
    } else {
        (state.reference_atoms.clone(), state.reference_pos.clone())
    };

    if ui
        .add(
            crate::ui_style::primary("Compare current frame")
                .min_size(egui::vec2(ui.available_width(), 36.0)),
        )
        .clicked_once()
    {
        let frame = &traj.frames[current];
        run_comparison(
            state,
            &frame.atoms,
            &frame.pos,
            &reference_atoms,
            &reference_pos,
        );
    }
    crate::ui_style::group(ui, "Alignment", |ui| {
        ui.small("Superpose the displayed frame or the whole trajectory onto the reference.");
        ui.horizontal_wrapped(|ui| {
            if ui.button("Align current frame").clicked_once() {
                align_frames(
                    state,
                    traj,
                    mol,
                    &reference_atoms,
                    &reference_pos,
                    AlignScope::Current(current),
                );
            }
            if ui.button("Align all frames").clicked_once() {
                align_frames(
                    state,
                    traj,
                    mol,
                    &reference_atoms,
                    &reference_pos,
                    AlignScope::All,
                );
            }
        });
    });

    if let Some(err) = &state.last_error {
        ui.colored_label(egui::Color32::from_rgb(220, 80, 80), err);
    }
    if let Some(msg) = &state.last_message {
        ui.colored_label(egui::Color32::from_rgb(80, 180, 100), msg);
    }

    if let Some(report) = &state.last_report {
        crate::ui_style::group(ui, "Comparison result", |ui| {
            ui.label(
                egui::RichText::new(format!("{:.4} Å", report.aligned_rmsd))
                    .size(24.0)
                    .strong()
                    .color(crate::ui_style::ACCENT),
            );
            ui.small("RMSD after best-fit alignment");
            egui::Grid::new("rmsd_summary")
                .striped(true)
                .show(ui, |ui| {
                    ui.label("Atoms compared");
                    ui.label(format!("{}", report.atom_count));
                    ui.end_row();
                    ui.label("RMSD (as placed)");
                    ui.label(format!("{:.4} Å", report.raw_rmsd));
                    ui.end_row();
                    ui.label("RMSD (best fit)");
                    ui.strong(format!("{:.4} Å", report.aligned_rmsd));
                    ui.end_row();
                    ui.label("Largest atom shift");
                    ui.label(format!("{:.4} Å", report.max_atom_deviation));
                    ui.end_row();
                });
            ui.weak(
                "\"As placed\" includes any overall drift or rotation; \"best fit\" \
             is after removing it, so it reflects real structural change only.",
            );
        });
    }
}

fn pin_reference(state: &mut RmsdState, traj: &TrajectoryState, frame_index: usize) {
    let frame = &traj.frames[frame_index.min(traj.frames.len() - 1)];
    state.reference_frame = Some(frame_index);
    state.reference_atoms = frame.atoms.clone();
    // Captured by value, not by index: realigning the trajectory in place
    // would otherwise silently move the reference out from under the user.
    state.reference_pos = frame.pos.clone();
    state.last_report = None;
    state.last_error = None;
    state.last_message = Some(format!(
        "Frame {} pinned as the reference.",
        frame_index + 1
    ));
}

fn run_comparison(
    state: &mut RmsdState,
    atoms: &[String],
    pos: &[Vec3],
    reference_atoms: &[String],
    reference_pos: &[Vec3],
) {
    if !same_atom_ordering(atoms, reference_atoms) {
        state.last_report = None;
        state.last_message = None;
        state.last_error = Some(
            "The two structures do not list the same elements in the same order, \
             so an atom-by-atom RMSD would be meaningless."
                .to_string(),
        );
        return;
    }
    match compare(pos, reference_pos) {
        Ok(report) => {
            state.last_report = Some(report);
            state.last_error = None;
            state.last_message = None;
        }
        Err(err) => {
            state.last_report = None;
            state.last_message = None;
            state.last_error = Some(err);
        }
    }
}

enum AlignScope {
    Current(usize),
    All,
}

fn align_frames(
    state: &mut RmsdState,
    traj: &mut TrajectoryState,
    mol: &mut Molecule,
    reference_atoms: &[String],
    reference_pos: &[Vec3],
    scope: AlignScope,
) {
    let indices: Vec<usize> = match scope {
        AlignScope::Current(i) => vec![i],
        AlignScope::All => (0..traj.frames.len()).collect(),
    };

    for &i in &indices {
        let frame = &traj.frames[i];
        if !same_atom_ordering(&frame.atoms, reference_atoms) {
            state.last_error = Some(format!(
                "Frame {} does not match the reference's atom ordering; nothing \
                 was aligned.",
                i + 1
            ));
            state.last_message = None;
            return;
        }
        let superposition = match superpose(&frame.pos, reference_pos) {
            Ok(s) => s,
            Err(err) => {
                state.last_error = Some(err);
                state.last_message = None;
                return;
            }
        };
        let aligned = superposition.apply_all(&frame.pos);
        traj.frames[i].pos = aligned;
    }

    // The displayed molecule is a copy of the current frame, so it has to be
    // refreshed too or the viewport keeps showing the pre-alignment pose
    // until the next scrub.
    let current = traj.current_frame.min(traj.frames.len() - 1);
    mol.pos = traj.frames[current].pos.clone();
    mol.recompute_bonds(2.0, 3.0);
    traj.last_applied = None;
    traj.overlay_dirty = true;

    // Report where things stand now that the fit has been applied.
    let frame = &traj.frames[current];
    run_comparison(
        state,
        &frame.atoms,
        &frame.pos,
        reference_atoms,
        reference_pos,
    );
    state.last_error = None;
    state.last_message = Some(match scope {
        AlignScope::Current(i) => format!("Frame {} aligned to the reference.", i + 1),
        AlignScope::All => format!("All {} frames aligned to the reference.", indices.len()),
    });
}
