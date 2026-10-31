//! Spectrum-specific controls in the general calculator.
use super::*;
use crate::qchem_panel::run::spectrum_job;
use crate::uvvis::excited_state::ExcitedStateMethod;

pub(super) fn active(panel: &QcPanelState) -> bool {
    panel.job == JobType::TdDft && matches!(panel.program, QcProgram::Behemoth | QcProgram::StdaXtb)
}

pub(super) fn settings(ui: &mut egui::Ui, panel: &mut QcPanelState, running: bool) {
    ui.add_enabled_ui(!running, |ui| {
        if panel.program == QcProgram::Behemoth {
            egui::ComboBox::from_id_salt("qc_spectrum_method")
                .selected_text(panel.excited.spectrum_method.label())
                .show_ui(ui, |ui| {
                    for method in [
                        ExcitedStateMethod::StdaXtbGfn1,
                        ExcitedStateMethod::StdaTasi,
                        ExcitedStateMethod::StdaDft,
                    ] {
                        ui.selectable_value(
                            &mut panel.excited.spectrum_method,
                            method,
                            method.label(),
                        );
                    }
                });
            ui.horizontal(|ui| {
                ui.label("States");
                ui.add(egui::DragValue::new(&mut panel.excited.states).range(1..=200));
            });
            if panel.excited.spectrum_method == ExcitedStateMethod::StdaDft {
                egui::ComboBox::from_id_salt("qc_spectrum_functional")
                    .selected_text(&panel.config.functional)
                    .show_ui(ui, |ui| {
                        for functional in crate::uvvis::excited_state::spectrum_functionals() {
                            ui.selectable_value(
                                &mut panel.config.functional,
                                functional.cli.to_owned(),
                                functional.label,
                            );
                        }
                    });
                ui.horizontal(|ui| {
                    ui.label("Basis");
                    egui::ComboBox::from_id_salt("qc_spectrum_basis")
                        .selected_text(&panel.config.basis)
                        .show_ui(ui, |ui| {
                            for group in method::BASIS_GROUPS {
                                ui.weak(group.name);
                                for (label, keyword) in group.sets {
                                    ui.selectable_value(
                                        &mut panel.config.basis,
                                        (*keyword).to_owned(),
                                        *label,
                                    );
                                }
                            }
                        });
                });
                ui.checkbox(&mut panel.excited.tda, "Tamm-Dancoff (sTDA)");
            }
        }
        ui.horizontal(|ui| {
            ui.label("Energy window (eV)");
            ui.add(egui::DragValue::new(&mut panel.excited.energy_window_ev).range(1..=100));
        });
    });
}

pub(super) fn generator_row(ui: &mut egui::Ui, state: &mut Presentation, running: bool) {
    ui.horizontal(|ui| {
        ui.label("xtb4stda path");
        let changed = ui
            .add_enabled(
                !running,
                egui::TextEdit::singleline(&mut state.generator_path)
                    .desired_width(180.0)
                    .hint_text("xtb4stda"),
            )
            .changed();
        let mut selected = false;
        if ui
            .add_enabled(!running, egui::Button::new("Browse…"))
            .clicked_once()
        {
            if let Some(path) = crate::recent_dir::pick_file(crate::recent_dir::open()) {
                state.generator_path = path.display().to_string();
                selected = true;
            }
        }
        if changed || selected {
            state.generator_error = spectrum_job::save_generator_path(&state.generator_path).err();
        }
    });
    crate::runtime_warning::observe(
        ui.ctx(),
        "xtb4stda_path_save",
        "Cannot save xtb4stda path",
        state.generator_error.as_deref(),
    );
}
