//! Shared compact rows for surface and atomic-property controls.
use bevy_egui::egui;
pub fn rows(ui: &mut egui::Ui, id: &str, contents: impl FnOnce(&mut egui::Ui)) {
    ui.scope(|ui| {
        ui.spacing_mut().slider_width = 100.;
        ui.spacing_mut().interact_size.y = 23.;
        ui.style_mut().drag_value_text_style = egui::TextStyle::Small;
        egui::Grid::new(id)
            .num_columns(2)
            .spacing([14., 5.])
            .show(ui, contents);
    });
}
pub fn slider(ui: &mut egui::Ui, label: &str, slider: egui::Slider<'_>) -> bool {
    ui.label(label);
    let changed = ui.add(slider).changed();
    ui.end_row();
    changed
}
