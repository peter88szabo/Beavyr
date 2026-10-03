//! Image-export controls. Capture and writing stay in the existing exporter.

use crate::ui_style::ResponseExt;

use crate::export_image::{ExportFormat, ExportSelectArea, ExportSettings};
use crate::ui_style;
use bevy_egui::egui;

fn select_area(area: &mut ExportSelectArea) {
    area.active = true;
    area.dragging = false;
    area.start = None;
    area.end = None;
    area.has_area = false;
}

pub fn panel(
    ui: &mut egui::Ui,
    settings: &mut ExportSettings,
    area: &mut ExportSelectArea,
) -> bool {
    ui.scope(|ui| {
        ui_style::modern(ui);
        ui_style::heading(
            ui,
            "Export the current view",
            "Choose a region and the output size.",
        );
        ui_style::group(ui, "File format", |ui| {
            let width = (ui.available_width() - 8.0) / 2.0;
            ui.horizontal(|ui| {
                for (format, label, hint) in [
                    (
                        ExportFormat::Png,
                        "PNG",
                        "Lossless image with sharp lines and text",
                    ),
                    (
                        ExportFormat::Jpg,
                        "JPG",
                        "Compressed image with a smaller file size",
                    ),
                ] {
                    if ui
                        .add_sized(
                            [width, 34.0],
                            egui::Button::new(label).selected(settings.format == format),
                        )
                        .on_hover_text(hint)
                        .clicked_once()
                    {
                        settings.format = format;
                    }
                }
            });
        });
        ui_style::group(ui, "Resolution", |ui| {
            ui.horizontal(|ui| {
                ui.label("Size (px)");
                ui.add(egui::DragValue::new(&mut settings.canvas_width).range(256..=8192));
                ui.label("×");
                ui.add(egui::DragValue::new(&mut settings.canvas_height).range(256..=8192));
            });
            ui.horizontal(|ui| {
                ui.label("Print resolution");
                egui::ComboBox::from_id_salt("export_dpi")
                    .selected_text(format!("{} DPI", settings.dpi))
                    .show_ui(ui, |ui| {
                        for dpi in [150_u32, 200, 300, 600, 1200] {
                            ui.selectable_value(&mut settings.dpi, dpi, format!("{dpi} DPI"));
                        }
                    });
            });
        });
        ui_style::group(ui, "Framing", |ui| {
            ui.horizontal_wrapped(|ui| {
                if ui
                    .button(if area.active {
                        "Selecting…"
                    } else if area.has_area {
                        "Redraw area"
                    } else {
                        "Select area"
                    })
                    .clicked_once()
                {
                    select_area(area);
                }
                if area.has_area && ui.button("Clear area").clicked_once() {
                    select_area(area);
                    area.active = false;
                }
            });
            if area.has_area {
                ui.label(
                    egui::RichText::new("Area selected · ready to export")
                        .small()
                        .color(ui_style::ACCENT),
                );
            } else {
                ui.small("Drag a rectangle in the 3D view to select the region.");
            }
        });
        ui.add_enabled(
            area.has_area,
            ui_style::primary(egui::RichText::new("Export image").color(egui::Color32::WHITE))
                .min_size([ui.available_width(), 36.0].into()),
        )
        .on_disabled_hover_text("Select an area in the 3D view first")
        .clicked_once()
    })
    .inner
}
