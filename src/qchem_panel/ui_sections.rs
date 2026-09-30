//! The runner panel as six collapsed sections.
//!
//! Every rule this draws lives in [`super::sections`], which has no drawing
//! dependency and carries the tests. This file only draws.
//!
//! Two properties are the point, and both are easy to lose in a later edit:
//!
//! * the whole setup reads without opening anything, because each header
//!   carries its own current setting;
//! * nothing appears or disappears, because a section that does not apply
//!   stays on screen and says why. A form whose rows come and go moves under
//!   the cursor, which is the complaint this layout answers.
//!
//! The controls inside the sections are the same ones the old layout drew,
//! called from [`super::ui_classic`] rather than copied. Two panels drawing
//! the same control two ways is how they drift.

use bevy_egui::egui;

use super::run::QcRunTask;
use super::sections::{self, Section};
use super::QcPanelState;
use crate::molecule::Molecule;
use crate::qchem_interfaces::valence::validate_electronic_state;
use crate::qchem_interfaces::xtb_optimize::{
    python_environment_ready, python_environment_row, XtbPanelState,
};

pub(crate) fn body(
    ui: &mut egui::Ui,
    panel: &mut QcPanelState,
    task: &mut QcRunTask,
    opt_panel: &mut XtbPanelState,
    mol: &Molecule,
) {
    let running = task.is_running();

    ui.label(egui::RichText::new(panel.summary_line()).strong());
    ui.separator();

    // Both gates on the Run button are settled HERE, outside the sections,
    // and nowhere else.
    //
    // This is the one rule this file cannot break. egui does not run a
    // collapsed CollapsingHeader's body, so a value assigned inside a section
    // is not assigned at all while that section is shut -- it keeps whatever
    // default it was given, and the default for a gate is "allowed". Set
    // python_ready inside the Program section and collapsing that section
    // re-enables Run for a backend this Python cannot import; the user waits
    // for the job and gets an import error instead.
    let python_ready = python_environment_ready(panel.program);
    // level_of_theory draws its verdict but does not gate on it: every path
    // through it returns true today, and a missing method is refused by
    // run.rs with its own message. If it ever gains a false path, that verdict
    // belongs in the model beside python_environment_ready, not in a closure.
    let level_ok = true;

    for section in sections::ALL {
        let heading = format!(
            "{}  \u{2014}  {}",
            sections::title(section),
            sections::header(section, panel)
        );

        egui::CollapsingHeader::new(heading)
            .id_salt(format!("qc_section_{section:?}"))
            // Only Program starts open, so a freshly opened panel shows the
            // whole setup at a glance rather than one expanded wall of
            // controls with the rest below the fold.
            .default_open(section == Section::Program)
            .show(ui, |ui| {
                // A section with nothing to set says why, in place. It is not
                // removed: that is what made the old form jump.
                if let Some(note) = sections::unavailable_note(section, panel) {
                    ui.label(egui::RichText::new(note).small().weak());
                    return;
                }
                match section {
                    Section::Program => {
                        super::ui_classic::program_row(ui, panel, opt_panel, running);
                        // Drawn here, decided above.
                        let _ = python_environment_row(ui, panel.program, running);
                    }
                    Section::LevelOfTheory => {
                        let _ = super::ui_classic::level_of_theory(ui, panel, running);
                    }
                    Section::Job => {
                        super::ui_classic::job_row(ui, panel, running);
                        if panel.needs_excited_state_options() {
                            super::ui_classic::excited_state_row(ui, panel, running);
                        }
                    }
                    Section::Molecule => {
                        super::ui_classic::charge_and_multiplicity(ui, panel, running);
                        electronic_state(ui, panel, mol);
                    }
                    Section::Resources => {
                        super::ui_classic::resources_row(ui, panel, running);
                    }
                    Section::InputFile => {
                        super::ui_classic::wavefunction_row(ui, panel, running);
                        super::ui_classic::input_row(ui, panel, mol, running);
                    }
                }
            });
    }

    // Run, its status line and the summary button stay below the sections and
    // always visible: they are what the panel is for, and hiding the Run
    // button inside a collapsible would be its own small joke.
    let electronic = validate_electronic_state(mol, panel.charge, panel.multiplicity);

    ui.separator();
    super::ui_classic::run_row(
        ui,
        panel,
        task,
        opt_panel,
        mol,
        running,
        python_ready,
        level_ok,
        electronic.is_ok(),
    );
    super::ui_classic::results_row(ui, panel, task);
}

/// The charge and spin state's own verdict, inside the Molecule section.
///
/// On screen rather than in a tooltip: a run refused for an impossible
/// electronic state is something the user has to act on, and a reason they
/// cannot reach is no reason at all.
fn electronic_state(ui: &mut egui::Ui, panel: &QcPanelState, mol: &Molecule) {
    match validate_electronic_state(mol, panel.charge, panel.multiplicity) {
        Err(error) => {
            ui.colored_label(
                egui::Color32::from_rgb(230, 120, 90),
                format!("\u{26a0} {error}"),
            );
        }
        // Warnings wait until a run has been attempted, matching the other
        // panels: unsolicited valence warnings about a perfectly ordinary
        // molecule read as noise.
        Ok((_, warnings)) if panel.show_warnings => {
            for warning in warnings {
                ui.label(egui::RichText::new(format!("\u{26a0} {}", warning.message)).weak());
            }
        }
        Ok(_) => {}
    }
}
