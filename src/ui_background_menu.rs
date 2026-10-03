//! Compact background context menu; the caller applies geometry actions.

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
    Close,
}

pub fn contents(
    ui: &mut egui::Ui,
    settings: &mut MolSettings,
    have_molecule: bool,
) -> Option<Action> {
    crate::ui_style::modern(ui);
    ui.set_min_width(240.0);
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
                action = Some(Action::Close);
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
    if ui
        .checkbox(&mut settings.show_hbonds, "Hydrogen bonds")
        .changed()
    {
        action = Some(Action::Close);
    }
    ui.menu_button("Background", |ui| {
        for (label, color) in [("Black", Color::BLACK), ("White", Color::WHITE)] {
            if ui.button(label).clicked_once() {
                settings.bg_color = color;
                settings.lighting_dirty = true;
                action = Some(Action::Close);
                ui.close();
            }
        }
    });
    action
}
