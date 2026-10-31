use super::*;
use crate::molecule::composition;
use crate::ui_style::ResponseExt;
use bevy_egui::egui;
use std::time::Duration;

impl StructureHistory {
    pub(super) fn delete_entry(&mut self, entry: &Entry) {
        match fs::remove_file(&entry.path) {
            Ok(()) => {}
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => {
                self.message = Some(format!(
                    "Could not remove {} from history: {error}",
                    entry.name
                ));
                self.is_error = true;
                return;
            }
        }
        self.entries.retain(|saved| saved.path != entry.path);
        self.message = Some(format!("Removed {} from history.", entry.name));
        self.is_error = false;
    }
}

fn status(ui: &mut egui::Ui, history: &StructureHistory) {
    if let Some(message) = &history.message {
        if history.is_error {
            ui.colored_label(egui::Color32::from_rgb(225, 110, 90), message);
        } else {
            ui.weak(message);
        }
    }
}
pub fn panel(
    ui: &mut egui::Ui,
    history: &mut StructureHistory,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
) -> bool {
    if history.refreshed.elapsed() >= Duration::from_secs(2) {
        history.refresh();
    }
    let mut selected = None;
    let mut removed = None;
    ui.scope(|ui| {
        crate::ui_style::modern(ui);
        ui.horizontal(|ui| {
            ui.small(format!("{} of 20 saved", history.entries.len()));
            ui.menu_button("How history works", |ui| {
                ui.set_max_width(320.0);
                ui.label("Loading another structure or clearing the display preserves the previous one.");
                ui.label("Use Save in Molecule Editor to keep an edited version. Edits and trajectory playback do not add entries automatically.");
            });
        });
        status(ui, history);
        for (index, entry) in history.entries.iter().enumerate() {
            let comp = composition(&entry.snapshot.atoms);
            ui.push_id(&entry.path, |ui| {
                crate::ui_style::group(ui, "", |ui| {
                    ui.horizontal(|ui| {
                        ui.label(egui::RichText::new(format!("{:02}", index + 1)).small().color(crate::ui_style::ACCENT));
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            if ui
                                .add(crate::ui_style::primary(egui::RichText::new("Load").small()))
                                .clicked_once()
                            {
                                selected = Some(entry.clone());
                            }
                            if ui.small_button("Delete").on_hover_text("Remove this saved snapshot from Recent Structures").clicked_once() {
                                removed = Some(entry.clone());
                            }
                            ui.add(egui::Label::new(egui::RichText::new(&entry.name).small().strong()).truncate())
                                .on_hover_text(&entry.name);
                        });
                    });
                    ui.horizontal_wrapped(|ui| {
                        ui.label(egui::RichText::new(&comp.formula).color(crate::ui_style::ACCENT));
                        ui.small(format!("{} atoms", comp.natoms));
                        let age = SystemTime::now()
                            .duration_since(UNIX_EPOCH).unwrap_or_default().as_nanos()
                            .saturating_sub(entry.saved) / 1_000_000_000;
                        ui.weak(if age < 60 { "Just now".into() }
                            else if age < 3600 { format!("{} min ago", age / 60) }
                            else if age < 86400 { format!("{} h ago", age / 3600) }
                            else { format!("{} days ago", age / 86400) });
                    });
                });
            });
        }
        if history.entries.is_empty() {
            crate::ui_style::group(ui, "No saved structures yet", |ui| {
                ui.small("Save a version from Molecule Editor to keep it here.");
            });
        }
    });
    ui.ctx().request_repaint_after(Duration::from_secs(2));
    if let Some(entry) = removed {
        history.delete_entry(&entry);
    }
    if let Some(entry) = selected {
        history.restore(entry, mol, traj);
        return true;
    }
    false
}
