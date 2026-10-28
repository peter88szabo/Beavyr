//! A centered, acknowledged warning for unsupported DREIDING calculations.

use crate::qchem_panel::run::QcRunTask;
use bevy_egui::egui;

pub(crate) fn dreiding_warning_window(
    ctx: &egui::Context,
    task: &mut QcRunTask,
) -> Option<egui::Rect> {
    let message = task.dreiding_warning.as_deref()?;
    let response = egui::Modal::new(egui::Id::new("unsupported_dreiding_system")).show(ctx, |ui| {
        ui.set_width(480.0_f32.min((ctx.content_rect().width() - 48.0).max(160.0)));
        ui.heading("DREIDING unavailable for this system");
        ui.add_space(8.0);
        egui::ScrollArea::vertical()
            .max_height((ctx.content_rect().height() - 160.0).max(100.0))
            .show(ui, |ui| {
                ui.label(message);
            });
        ui.add_space(12.0);
        ui.button("OK").clicked()
    });
    if response.inner || response.should_close() {
        task.dreiding_warning = None;
    }
    Some(response.response.rect)
}

#[cfg(test)]
#[path = "dreiding_warning_tests.rs"]
mod tests;
