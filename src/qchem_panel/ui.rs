//! Drawing the general quantum-chemistry panel.
//!
//! This file owns the window and the choice of layout. The two layouts
//! themselves live beside it: [`super::ui_sections`] is the current one and
//! [`super::ui_classic`] is the flat form it replaced, kept reachable.

use bevy_egui::egui;

use super::run::QcRunTask;
use super::PanelLayout;
use super::QcPanelState;
use crate::molecule::Molecule;
use crate::qchem_interfaces::program::QcProgram;
use crate::qchem_interfaces::xtb_optimize::XtbPanelState;

/// The control that swaps layouts, at the foot of the panel.
///
/// Worded as what it does rather than as a setting, and placed last so it is
/// out of the way of ordinary use without being hidden.
fn layout_switch(ui: &mut egui::Ui, panel: &mut QcPanelState) {
    ui.separator();
    ui.horizontal(|ui| {
        if ui
            .small_button(format!("use {}", panel.layout.other().label()))
            .clicked()
        {
            panel.layout = panel.layout.other();
            crate::qchem_interfaces::config::save_layout(panel.layout);
        }
    });
}

/// Draws the panel and starts a run when the button is pressed.
#[allow(clippy::too_many_arguments)]
pub fn qc_panel(
    ctx: &egui::Context,
    panel: &mut QcPanelState,
    task: &mut QcRunTask,
    opt_panel: &mut XtbPanelState,
    opt_task: &crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask,
    mol: &Molecule,
) {
    let mut open = panel.open;
    egui::Window::new("Quantum Chemistry")
        .open(&mut open)
        .default_size([460.0, 640.0])
        .vscroll(true)
        .show(ctx, |ui| {
            match panel.layout {
                // Two layouts, one set of controls: the sections call into
                // ui_classic's rows rather than copying them, so a capability
                // added once appears in both.
                PanelLayout::Sections => super::ui_sections::body(ui, panel, task, opt_panel, mol),
                PanelLayout::Classic => super::ui_classic::body(ui, panel, task, opt_panel, mol),
            }
            layout_switch(ui, panel);
        });
    panel.open = open;

    // The optimizer panel's own two windows, driven from its own flags. An
    // optimisation started here and one started there are the same
    // calculation, so they get the same summary and the same energy plot
    // rather than a second pair built to look almost alike.
    let report = opt_task
        .last_log
        .as_deref()
        .filter(|_| {
            opt_task.last_program == Some(QcProgram::Behemoth)
                || opt_task.last_program == Some(QcProgram::Dreiding)
        })
        .map(|log| {
            crate::qchem_interfaces::xtb_optimize::OptimizationReport::from_log(
                log,
                opt_task.last_duration,
            )
        });
    crate::qchem_interfaces::xtb_optimize::energy_history_window(
        ctx,
        &mut opt_panel.energy_plot_open,
        &opt_task.last_energy_history,
        opt_task.last_is_error,
    );
    if let Some(report) = &report {
        crate::qchem_interfaces::xtb_optimize::optimization_summary_window(
            ctx,
            &mut opt_panel.summary_open,
            report,
        );
    }
}
