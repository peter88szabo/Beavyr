//! A centered, acknowledged warning for unsupported DREIDING calculations.

use crate::qchem_panel::run::QcRunTask;
use bevy_egui::egui;

pub(crate) fn dreiding_warning_window(
    ctx: &egui::Context,
    task: &mut QcRunTask,
) -> Option<egui::Rect> {
    let message = task.dreiding_warning.as_deref()?;
    let response = egui::Modal::new(egui::Id::new("unsupported_dreiding_system")).show(ctx, |ui| {
        crate::runtime_warning::message_body(ui, "DREIDING unavailable for this system", &message)
    });
    if response.inner || response.should_close() {
        task.dreiding_warning = None;
    }
    Some(response.response.rect)
}

#[cfg(test)]
#[path = "dreiding_warning_tests.rs"]
mod tests;
