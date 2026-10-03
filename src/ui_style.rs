//! Shared presentation for the floating tools. Styling stays local to each UI.

use bevy_egui::egui;

pub const ACCENT: egui::Color32 = egui::Color32::from_rgb(104, 174, 255);

pub trait ResponseExt {
    /// Process an action once even when egui retries layout in the same frame.
    fn clicked_once(&self) -> bool;
}

impl ResponseExt for egui::Response {
    fn clicked_once(&self) -> bool {
        if !self.clicked() {
            return false;
        }
        let frame = self.ctx.cumulative_frame_nr();
        let id = self.id.with("beavyr_action_frame");
        self.ctx.data_mut(|data| {
            if data.get_temp::<u64>(id) == Some(frame) {
                false
            } else {
                data.insert_temp(id, frame);
                true
            }
        })
    }
}

pub fn modern(ui: &mut egui::Ui) {
    let style = ui.style_mut();
    style.spacing.item_spacing = egui::vec2(8.0, 8.0);
    style.spacing.button_padding = egui::vec2(10.0, 6.0);
    style.spacing.interact_size.y = 30.0;
    style.visuals.selection.bg_fill = egui::Color32::from_rgb(35, 75, 117);
    style.visuals.selection.stroke = egui::Stroke::new(1.0, ACCENT);
    for widget in [
        &mut style.visuals.widgets.inactive,
        &mut style.visuals.widgets.hovered,
        &mut style.visuals.widgets.active,
    ] {
        widget.corner_radius = egui::CornerRadius::same(6);
    }
}

pub fn heading(ui: &mut egui::Ui, title: &str, description: &str) {
    ui.label(egui::RichText::new(title).size(18.0).strong());
    if !description.is_empty() {
        ui.label(egui::RichText::new(description).small().weak());
    }
}

pub fn card(ui: &egui::Ui) -> egui::Frame {
    egui::Frame::new()
        .fill(ui.visuals().faint_bg_color)
        .stroke(ui.visuals().widgets.noninteractive.bg_stroke)
        .corner_radius(8)
        .inner_margin(10)
}

pub fn window_frame(ctx: &egui::Context) -> egui::Frame {
    egui::Frame::window(&ctx.global_style())
        .corner_radius(10)
        .inner_margin(12)
}

pub fn collapsible_group<R>(
    ui: &mut egui::Ui,
    title: &str,
    default_open: bool,
    contents: impl FnOnce(&mut egui::Ui) -> R,
) {
    card(ui).show(ui, |ui| {
        ui.set_min_width(ui.available_width().max(0.0));
        egui::CollapsingHeader::new(egui::RichText::new(title).strong())
            .default_open(default_open)
            .show(ui, contents);
    });
}

pub fn group<R>(ui: &mut egui::Ui, title: &str, contents: impl FnOnce(&mut egui::Ui) -> R) -> R {
    card(ui)
        .show(ui, |ui| {
            ui.set_min_width(ui.available_width().max(0.0));
            if !title.is_empty() {
                ui.strong(title);
            }
            contents(ui)
        })
        .inner
}

pub fn primary(text: impl Into<egui::WidgetText>) -> egui::Button<'static> {
    egui::Button::new(text.into().color(egui::Color32::WHITE))
        .fill(egui::Color32::from_rgb(38, 96, 159))
}
