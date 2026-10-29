//! Configuration shared by the calculator, optimizer, and frequency panels.

use super::{cached, launch, notify_missing};
use crate::ui_style::ResponseExt;
use bevy_egui::egui;

#[derive(Clone)]
struct Draft {
    settings: launch::Settings,
    error: Option<String>,
}

fn path_row(ui: &mut egui::Ui, label: &str, value: &mut String, hint: &str) {
    ui.label(label);
    ui.horizontal(|ui| {
        ui.add(
            egui::TextEdit::singleline(value)
                .desired_width((ui.available_width() - 86.0).max(100.0))
                .hint_text(hint),
        );
        if ui.button("Browse…").clicked_once() {
            if let Some(path) = rfd::FileDialog::new().pick_file() {
                *value = path.display().to_string();
            }
        }
    });
}

pub(crate) fn environment_row(ui: &mut egui::Ui, module: &'static str, running: bool) -> bool {
    let id = egui::Id::new(("python_environment_draft", module));
    let saved = launch::settings(module);
    let mut draft = ui
        .ctx()
        .data_mut(|data| data.get_temp::<Draft>(id))
        .unwrap_or_else(|| Draft {
            settings: saved.as_ref().cloned().unwrap_or_default(),
            error: saved.as_ref().err().cloned(),
        });
    let mut check = false;
    ui.push_id(module, |ui| {
        ui.add_enabled_ui(!running, |ui| {
            ui.label("Environment activation");
            let hint = if cfg!(windows) {
                "& 'C:\\path\\to\\venv\\Scripts\\Activate.ps1'"
            } else {
                "source /path/to/venv/bin/activate"
            };
            ui.add(
                egui::TextEdit::singleline(&mut draft.settings.activation)
                    .desired_width(ui.available_width())
                    .hint_text(hint),
            );
            ui.small("Or: conda activate myenv / micromamba activate myenv");
            ui.collapsing("Python interpreter and activation shell", |ui| {
                path_row(
                    ui,
                    "Python interpreter",
                    &mut draft.settings.interpreter,
                    super::PYTHON,
                );
                path_row(
                    ui,
                    "Activation shell",
                    &mut draft.settings.shell,
                    if cfg!(windows) {
                        "powershell.exe"
                    } else if cfg!(target_os = "macos") {
                        "/bin/zsh"
                    } else {
                        "/bin/bash"
                    },
                );
            });
            ui.horizontal(|ui| {
                check = ui.button("Save & Check").clicked_once();
                if saved.as_ref().is_ok_and(|value| value != &draft.settings) {
                    ui.small("Unsaved changes");
                }
            });
        });
    });
    if check {
        match launch::save(module, draft.settings.clone()) {
            Ok(()) => {
                draft.error = None;
                notify_missing(ui.ctx(), &cached(module));
            }
            Err(error) => {
                crate::runtime_warning::push(ui.ctx(), "Python settings error", &error);
                draft.error = Some(error);
            }
        }
    }
    if let Some(error) = &draft.error {
        ui.colored_label(egui::Color32::from_rgb(230, 120, 90), error);
    }
    let environment = cached(module);
    let ready = environment.is_ready();
    if ready {
        ui.label(egui::RichText::new(environment.short()).small().weak())
            .on_hover_text(environment.detail());
    } else {
        ui.colored_label(egui::Color32::from_rgb(230, 120, 90), environment.short())
            .on_hover_text(environment.detail());
    }
    ui.ctx().data_mut(|data| data.insert_temp(id, draft));
    ready
}
