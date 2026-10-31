//! Scan setup and ordered viewport atom selection.
use super::{run::QcRunTask, QcPanelState};
use crate::qchem_interfaces::job::{scan::*, JobType};
use crate::qchem_interfaces::program::QcProgram;
use crate::ui_style::ResponseExt;
use bevy::prelude::*;
use bevy_egui::egui;

pub fn collect_picks(
    mut clicks: MessageReader<crate::events::ViewportClicked>,
    mut panel: ResMut<QcPanelState>,
    task: Res<QcRunTask>,
    mol: Res<crate::molecule::Molecule>,
) {
    for click in clicks.read() {
        if !panel.open
            || panel.job != JobType::Scan
            || task.is_running()
            || !panel.config.scan.picking
        {
            continue;
        }
        let scan = &mut panel.config.scan;
        if let Some(i) = click
            .hit
            .filter(|i| *i < mol.atoms.len() && !crate::molecule::is_dummy(&mol.atoms[*i]))
        {
            if !scan.picked.contains(&i) {
                scan.picked.push(i);
            }
            scan.atoms = scan
                .picked
                .iter()
                .map(|i| (i + 1).to_string())
                .collect::<Vec<_>>()
                .join(" ");
            if scan.picked.len() == scan.coordinate.count() {
                scan.picking = false;
            }
        }
    }
}
fn coordinate_picker(
    ui: &mut egui::Ui,
    id: impl std::hash::Hash + std::fmt::Debug,
    value: &mut Coordinate,
) {
    egui::ComboBox::from_id_salt(id)
        .selected_text(value.label())
        .show_ui(ui, |ui| {
            for c in [Coordinate::Bond, Coordinate::Angle, Coordinate::Dihedral] {
                ui.selectable_value(value, c, c.label());
            }
        });
}
pub(super) fn settings(
    ui: &mut egui::Ui,
    panel: &mut QcPanelState,
    mol: &crate::molecule::Molecule,
    task: &QcRunTask,
) {
    crate::ui_style::group(ui, "Scan settings", |ui| {
        ui.add_enabled_ui(!task.is_running(), |ui| {
            let scan = &mut panel.config.scan;
            if panel.program == QcProgram::Orca {
                ui.horizontal(|ui| { ui.label("Driver"); ui.selectable_value(&mut scan.driver, Driver::Orca, "ORCA built-in"); ui.selectable_value(&mut scan.driver, Driver::Beavyr, "Beavyr"); });
            }
            if panel.program == QcProgram::Behemoth {
                ui.horizontal(|ui| { ui.label("Driver"); ui.selectable_value(&mut scan.driver, Driver::Behemoth, "Behemoth built-in"); ui.selectable_value(&mut scan.driver, Driver::Beavyr, "Beavyr"); });
            }
            ui.horizontal(|ui| {
                let rigid = ui.add_enabled(scan.constraints.is_empty(), egui::Button::new("Rigid").selected(scan.mode == Mode::Rigid));
                if rigid.clicked_once() { scan.mode = Mode::Rigid; }
                rigid.on_disabled_hover_text("Remove fixed coordinates before choosing a rigid scan.");
                ui.selectable_value(&mut scan.mode, Mode::Relaxed, "Relaxed");
                let previous = scan.coordinate;
                coordinate_picker(ui, "scan_coordinate", &mut scan.coordinate);
                if previous != scan.coordinate {
                    scan.atoms = (1..=scan.coordinate.count()).map(|i| i.to_string()).collect::<Vec<_>>().join(" ");
                    (scan.start, scan.end) = match scan.coordinate { Coordinate::Bond => (1.,2.), Coordinate::Angle => (90.,120.), Coordinate::Dihedral => (0.,180.) };
                    scan.picking = false; scan.picked.clear();
                }
            });
            ui.horizontal(|ui| {
                ui.label("Atoms"); ui.add(egui::TextEdit::singleline(&mut scan.atoms).desired_width(110.));
                if ui.button(if scan.picking { "Stop picking" } else { "Pick atoms" }).clicked_once() {
                    scan.picking = !scan.picking;
                    if scan.picking { scan.picked.clear(); scan.atoms.clear(); }
                }
            });
            if scan.picking { ui.label(format!("Pick {} atoms in order: {} selected", scan.coordinate.count(), scan.picked.len())); }
            ui.horizontal(|ui| {
                ui.label(format!("Start ({})", scan.coordinate.unit())); ui.add(egui::DragValue::new(&mut scan.start).speed(0.01));
                ui.label("End"); ui.add(egui::DragValue::new(&mut scan.end).speed(0.01));
            });
            ui.horizontal(|ui| { ui.label("Points"); ui.add(egui::DragValue::new(&mut scan.points).range(2..=10000)); });
            if ui.button("Use current value as start").clicked_once() {
                if let Ok(indices) = atom_indices(&scan.atoms, scan.coordinate.count(), &mol.atoms) {
                    let coords: Vec<f64> = mol.pos.iter().flat_map(|p| p.to_array().map(f64::from)).collect();
                    if let Ok(value) = super::run::scan_job::geometry::measure(scan.coordinate, &indices, &coords) { scan.start = value; }
                }
            }
            if scan.mode == Mode::Rigid {
                if scan.coordinate == Coordinate::Bond {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Move:");
                        ui.selectable_value(&mut scan.bond_motion, BondMotion::Second, "Second");
                        ui.selectable_value(&mut scan.bond_motion, BondMotion::First, "First");
                        ui.selectable_value(&mut scan.bond_motion, BondMotion::Symmetric, "Both symmetrically");
                    });
                }
                ui.checkbox(&mut scan.move_fragment, "Move connected fragment")
                    .on_hover_text("Move the attached group together. Disable to move only the final selected atom. Rings cannot be split into a moving fragment.");
            } else {
                if panel.program != QcProgram::Orca || scan.driver == Driver::Beavyr {
                    egui::ComboBox::from_id_salt("scan_relaxation_coordinates")
                        .selected_text(match scan.relaxation_coordinates { RelaxationCoordinates::Primitive => "Primitive internal", RelaxationCoordinates::Pulay => "Pulay redundant internal" })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut scan.relaxation_coordinates, RelaxationCoordinates::Primitive, "Primitive internal");
                            ui.selectable_value(&mut scan.relaxation_coordinates, RelaxationCoordinates::Pulay, "Pulay redundant internal");
                        });
                }
                ui.horizontal(|ui| { ui.label("Max steps / point"); ui.add(egui::DragValue::new(&mut scan.max_cycles).range(1..=10000)); });
                // ORCA controls its own continuation; its per-point statuses are still reported.
                if panel.program != QcProgram::Orca || scan.driver == Driver::Beavyr {
                    ui.checkbox(&mut scan.stop_on_failure, "Stop at unconverged point");
                }
                egui::CollapsingHeader::new("Fixed coordinates").id_salt("scan_fixed_coordinates").show(ui, |ui| {
                    let mut remove = None;
                    for (i, fixed) in scan.constraints.iter_mut().enumerate() {
                        ui.push_id(i, |ui| {
                            ui.horizontal(|ui| {
                                coordinate_picker(ui, "coordinate", &mut fixed.coordinate);
                                ui.add(egui::TextEdit::singleline(&mut fixed.atoms).desired_width(80.));
                                ui.add(egui::DragValue::new(&mut fixed.value).speed(0.01).suffix(fixed.coordinate.unit()));
                                if ui.small_button("×").clicked_once() { remove = Some(i); }
                            });
                        });
                    }
                    if let Some(i) = remove { scan.constraints.remove(i); }
                    if ui.button("Add fixed coordinate").clicked_once() {
                        scan.constraints.push(FixedCoordinate { coordinate: Coordinate::Bond, atoms: String::new(), value: 1. });
                    }
                });
            }
        });
        if let Some(progress) = super::run::scan_job::progress(task) {
            ui.label(progress);
        }
    });
}

#[cfg(test)]
#[path = "scan_controls_tests.rs"]
mod tests;
