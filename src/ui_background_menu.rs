//! Compact background context menu; the caller applies geometry actions.
//! Every entry leaves the menu open; it closes with its Close button.

use crate::ui_style::ResponseExt;

use crate::orientation::{Axis, Plane};
use crate::settings::{MolSettings, RepresentationMode};
use bevy::prelude::Color;
use bevy_egui::egui;

#[derive(Clone, Copy)]
pub enum Action {
    Center,
    Orient(Plane),
    Mirror(Plane),
    Flip(Axis),
}

pub fn contents(
    ui: &mut egui::Ui,
    settings: &mut MolSettings,
    have_molecule: bool,
) -> Option<Action> {
    crate::ui_style::modern(ui);
    // As wide as the widest entry, with every button stretched to that width.
    ui.with_layout(egui::Layout::top_down_justified(egui::Align::LEFT), |ui| {
        entries(ui, settings, have_molecule)
    })
    .inner
}

/// The persisted atom-label switches the Structure panel also shows.
const LABEL_TOGGLES: [(&str, &str); 2] = [
    ("show_atom_type", "Element symbols"),
    ("show_atom_index", "Atom numbers"),
];

fn entries(ui: &mut egui::Ui, settings: &mut MolSettings, have_molecule: bool) -> Option<Action> {
    ui.label(egui::RichText::new("Molecule tools").strong().size(16.0));
    ui.separator();
    let mut action = None;
    if ui
        .add_enabled(have_molecule, egui::Button::new("Center molecule"))
        .clicked_once()
    {
        action = Some(Action::Center);
    }
    ui.menu_button("Representation", |ui| {
        for (mode, label) in [
            (RepresentationMode::BallAndStick, "Ball and stick"),
            (RepresentationMode::SpaceFilling, "Space filling (CPK)"),
            (RepresentationMode::SticksRounded, "Sticks"),
            (
                RepresentationMode::LowResBallsAndLines,
                "Low-res balls and lines",
            ),
            (RepresentationMode::LinesOnly, "Lines"),
            (RepresentationMode::BackboneTrace, "Backbone / trace"),
        ] {
            if ui
                .selectable_value(&mut settings.representation, mode, label)
                .clicked_once()
            {
                settings.geometry_dirty = true;
                ui.close();
            }
        }
    });
    ui.add_enabled_ui(have_molecule, |ui| {
        ui.menu_button("Orientation", |ui| {
            ui.menu_button("Orient to plane", |ui| {
                for plane in Plane::ALL {
                    if ui.button(plane.label()).clicked_once() {
                        action = Some(Action::Orient(plane));
                        ui.close();
                    }
                }
            });
            ui.menu_button("Mirror across plane", |ui| {
                for plane in Plane::ALL {
                    if ui.button(plane.label()).clicked_once() {
                        action = Some(Action::Mirror(plane));
                        ui.close();
                    }
                }
            })
            .response
            .on_hover_text("Reflection; produces the enantiomer of a chiral molecule");
            ui.menu_button("Flip around axis", |ui| {
                for axis in Axis::ALL {
                    if ui.button(axis.label()).clicked_once() {
                        action = Some(Action::Flip(axis));
                        ui.close();
                    }
                }
            })
            .response
            .on_hover_text("Rotate by 180° while preserving chirality");
        });
    });
    ui.separator();
    for (key, label) in LABEL_TOGGLES {
        let id = egui::Id::new(key);
        let mut on = ui.ctx().data_mut(|d| d.get_persisted::<bool>(id).unwrap_or(false));
        if ui.checkbox(&mut on, label).changed() {
            ui.ctx().data_mut(|d| d.insert_persisted(id, on));
        }
    }
    ui.checkbox(&mut settings.show_hbonds, "Hydrogen bonds");
    ui.menu_button("Background", |ui| {
        for (label, color) in [("Black", Color::BLACK), ("White", Color::WHITE)] {
            if ui.button(label).clicked_once() {
                settings.bg_color = color;
                settings.lighting_dirty = true;
                ui.close();
            }
        }
    });
    action
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn menu_is_as_wide_as_its_widest_entry_and_has_label_toggles() {
        let ctx = egui::Context::default();
        let mut settings = MolSettings::default();
        let mut output = None;
        for _ in 0..3 {
            let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
                egui::Area::new(egui::Id::new("menu_width_test"))
                    .fixed_pos(egui::pos2(10.0, 10.0))
                    .show(ui.ctx(), |ui| {
                        contents(ui, &mut settings, true);
                    });
            });
            out.textures_delta.clear();
            output = Some(out);
        }
        let output = output.unwrap();
        let mut texts = vec![];
        fn collect(shape: &egui::Shape, out: &mut Vec<(String, egui::Rect)>) {
            match shape {
                egui::Shape::Text(t) => out.push((
                    t.galley.text().to_string(),
                    egui::Rect::from_min_size(t.pos, t.galley.size()),
                )),
                egui::Shape::Vec(v) => v.iter().for_each(|s| collect(s, out)),
                _ => {}
            }
        }
        output.shapes.iter().for_each(|s| collect(&s.shape, &mut texts));
        for label in ["Element symbols", "Atom numbers", "Hydrogen bonds"] {
            assert!(texts.iter().any(|(t, _)| t == label), "missing {label}");
        }
        let widest = texts.iter().map(|(_, r)| r.width()).fold(0.0, f32::max);
        let area = ctx.memory(|m| m.area_rect(egui::Id::new("menu_width_test"))).unwrap();
        // Text plus button padding and frame, but no fixed 240 px minimum.
        assert!(area.width() < widest + 80.0, "menu {area:?} widest text {widest}");
    }
}
