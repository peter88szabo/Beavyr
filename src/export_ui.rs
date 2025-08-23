use bevy::prelude::*;
use bevy_egui::egui;

use crate::export::{ExportParams, ExportRequestQueue, ExportUiState};

const MM_PER_INCH: f32 = 25.4;

/// Draws the export UI inline inside an existing egui panel/section.
pub fn export_section(
    ui: &mut egui::Ui,
    ui_state: &mut ExportUiState,
    queue: &mut ExportRequestQueue,
) {
    // Save button
    if ui.button("Choose folder & save").clicked() {
        let mut dlg = rfd::FileDialog::new();
        if let Some(dir) = &ui_state.last_dir {
            dlg = dlg.set_directory(dir);
        }

        if let Some(folder) = dlg.pick_folder() {
            ui_state.last_dir = Some(folder.clone());

            // Compute pixel size from mm + dpi
            let w_px = ((ui_state.width_mm / MM_PER_INCH) * ui_state.dpi as f32)
                .round()
                .max(1.0) as u32;
            let h_px = ((ui_state.height_mm / MM_PER_INCH) * ui_state.dpi as f32)
                .round()
                .max(1.0) as u32;

            queue.queue.push(ExportParams {
                output_dir: folder,
                extension: "png".into(),
                width: w_px.max(16),
                height: h_px.max(16),
            });
        }
    }

    // Small hint below button
    ui.label(egui::RichText::new("Saves numbered PNG (e.g. 00000.png)").small());

    ui.add_space(6.0);
    ui.separator();
    ui.add_space(4.0);

    // DPI block (row)
    ui.horizontal(|ui| {
        ui.label("DPI:");
        egui::ComboBox::from_id_salt("dpi_combo_inline")
            .selected_text(format!("{}", ui_state.dpi))
            .show_ui(ui, |ui| {
                for &d in &[300_u32, 600_u32, 1200_u32] {
                    ui.selectable_value(&mut ui_state.dpi, d, format!("{}", d));
                }
            });
    });

    ui.add_space(6.0);

    // Size block (row)
    ui.horizontal(|ui| {
        ui.label("Size (mm):");
        ui.add(
            egui::DragValue::new(&mut ui_state.width_mm)
                .range(5.0..=2000.0)
                .speed(1.0),
        );
        ui.label("×");
        ui.add(
            egui::DragValue::new(&mut ui_state.height_mm)
                .range(5.0..=2000.0)
                .speed(1.0),
        );
    });

    // Computed px preview
    let w_px = ((ui_state.width_mm / MM_PER_INCH) * ui_state.dpi as f32)
        .round()
        .max(1.0) as u32;
    let h_px = ((ui_state.height_mm / MM_PER_INCH) * ui_state.dpi as f32)
        .round()
        .max(1.0) as u32;
    ui.add_space(4.0);
    ui.label(egui::RichText::new(format!("→ {} × {} px", w_px, h_px)).small());
}

