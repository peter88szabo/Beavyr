//! Two-state controls for a minimum-energy crossing point calculation.
use super::{run::QcRunTask, QcPanelState};
use crate::qchem_interfaces::job::mecp::{Coordinates, Driver};
use crate::qchem_interfaces::program::QcProgram;
use bevy_egui::egui;
pub(super) fn settings(ui: &mut egui::Ui, panel: &mut QcPanelState, task: &QcRunTask) {
    crate::ui_style::group(ui, "MECP / Spin crossing", |ui| {
        ui.add_enabled_ui(!task.is_running(), |ui| {
            let o = &mut panel.config.mecp;
            if panel.program == QcProgram::Orca {
                ui.horizontal(|ui| {
                    ui.label("Driver");
                    ui.selectable_value(&mut o.driver, Driver::Orca, "ORCA built-in");
                    ui.selectable_value(&mut o.driver, Driver::Beavyr, "Beavyr MECP");
                });
            }
            ui.horizontal(|ui| {
                ui.label("State B multiplicity");
                ui.add(egui::DragValue::new(&mut o.multiplicity_b).range(1..=99));
            });
            ui.horizontal(|ui| {
                ui.label("Max steps");
                ui.add(egui::DragValue::new(&mut o.max_cycles).range(1..=10000));
            });
            ui.horizontal(|ui| {
                ui.label(if panel.program == QcProgram::Orca && o.driver == Driver::Orca { "Final gap limit (Eh)" } else { "Gap tolerance (Eh)" }).on_hover_text("ORCA built-in uses its own optimization convergence criteria; Beavyr also checks the final gap against this limit.");
                ui.add(
                    egui::DragValue::new(&mut o.gap_tolerance)
                        .speed(1e-6)
                        .range(1e-10..=1.)
                        .max_decimals(10),
                );
            });
            if panel.program != QcProgram::Orca || o.driver == Driver::Beavyr {
            egui::ComboBox::from_id_salt("mecp_coordinates")
                .selected_text(o.coordinates.label())
                .show_ui(ui, |ui| {
                    for c in [
                        Coordinates::Cartesian,
                        Coordinates::Primitive,
                        Coordinates::Pulay,
                    ] {
                        ui.selectable_value(&mut o.coordinates, c, c.label());
                    }
                });
            }
        });
    });
}
