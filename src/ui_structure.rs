//! Compact metadata and a coordinate area that fills the Structure window.

use std::path::Path;

use bevy_egui::egui;

use crate::molecule::Composition;
use crate::ui_style::{self, ResponseExt};

pub fn summary(ui: &mut egui::Ui, composition: &Composition, file: Option<&Path>) {
    let text = if composition.natoms == 0 {
        "No structure loaded".to_string()
    } else {
        format!("{} atoms · {}", composition.natoms, composition.formula)
    };
    let mut details = match composition.mass_amu {
        Some(mass) => format!("Molecular weight: {mass:.3} g/mol"),
        None if composition.natoms > 0 => format!(
            "Molecular weight unavailable — unknown elements: {}",
            composition.unknown.join(", ")
        ),
        None => "Open an XYZ file or paste coordinates to begin.".to_string(),
    };
    if let Some(path) = file {
        details.push_str(&format!("\nFile: {}", path.display()));
    }
    ui.add(egui::Label::new(egui::RichText::new(text).small().weak()).truncate())
        .on_hover_text(details);
}

#[derive(Debug, PartialEq, Eq)]
pub enum CoordinateAction {
    None,
    Apply,
    Clear,
    Reset,
}

/// The window supplies a finite height; the classic scrolling side panel
/// uses a fixed viewport instead. Only coordinates scroll in either case.
pub fn coordinates(
    ui: &mut egui::Ui,
    windowed: bool,
    edit_mode: bool,
    text: &mut String,
    live_lines: &[String],
) -> CoordinateAction {
    ui.small(if edit_mode {
        "XYZ coordinates · Å"
    } else {
        "Coordinates · Å · read only"
    });
    let footer_height = if edit_mode {
        ui.spacing().interact_size.y.max(30.0) + ui.spacing().item_spacing.y
    } else {
        0.0
    };
    let available_height = if windowed {
        ui.available_height()
    } else {
        360.0
    };
    let height = (available_height - footer_height - 16.0).max(40.0);
    egui::Frame::new()
        .fill(ui.visuals().extreme_bg_color)
        .stroke(ui.visuals().widgets.noninteractive.bg_stroke)
        .corner_radius(6)
        .inner_margin(8)
        .show(ui, |ui| {
            let width = ui.available_width();
            egui::ScrollArea::both()
                .id_salt("structure_coordinates")
                .max_height(height)
                .auto_shrink([false; 2])
                .show(ui, |ui| {
                    if edit_mode {
                        let line_height = ui.text_style_height(&egui::TextStyle::Monospace)
                            + ui.spacing().extra_text_line_spacing;
                        ui.add(
                            egui::TextEdit::multiline(text)
                                .desired_width(width)
                                .desired_rows((height / line_height).ceil().max(1.0) as usize)
                                .code_editor()
                                .frame(egui::Frame::NONE)
                                .margin(0)
                                .lock_focus(true),
                        );
                    } else {
                        ui.spacing_mut().item_spacing.y = 2.0;
                        for line in live_lines {
                            ui.monospace(line);
                        }
                    }
                });
        });

    let mut action = CoordinateAction::None;
    if edit_mode {
        ui.horizontal(|ui| {
            if ui
                .add(ui_style::primary("Apply coordinates"))
                .clicked_once()
            {
                action = CoordinateAction::Apply;
            }
            if ui
                .button("Reset")
                .on_hover_text("Replace the draft with the current structure.")
                .clicked_once()
            {
                action = CoordinateAction::Reset;
            }
            if ui
                .button("Clear")
                .on_hover_text("Clear the draft without changing the molecule.")
                .clicked_once()
            {
                action = CoordinateAction::Clear;
            }
        });
    }
    action
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame(
        ctx: &egui::Context,
        size: [f32; 2],
        events: Vec<egui::Event>,
        draft: &mut String,
        action: &mut CoordinateAction,
    ) -> egui::FullOutput {
        let mut output = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1280.0, 720.0),
                )),
                events,
                ..Default::default()
            },
            |ui| {
                crate::ui_layout::editor_section(
                    ui,
                    true,
                    &mut true,
                    &mut Vec::new(),
                    size,
                    "Structure layout check",
                    |ui| {
                        ui.horizontal_wrapped(|ui| {
                            for label in ["Open XYZ…", "View mode", "Copy XYZ", "Save XYZ"] {
                                let _ = ui.button(label);
                            }
                        });
                        ui.horizontal(|ui| {
                            summary(
                                ui,
                                &crate::molecule::composition(&vec!["C".into(); 3]),
                                None,
                            );
                            ui.menu_button("View", |_| {});
                        });
                        let next = coordinates(ui, true, true, draft, &[]);
                        if next != CoordinateAction::None {
                            *action = next;
                        }
                    },
                );
            },
        );
        output.textures_delta.clear();
        output
    }

    fn context() -> egui::Context {
        let ctx = egui::Context::default();
        ctx.global_style_mut(|style| {
            for (text_style, size) in [
                (egui::TextStyle::Body, 16.0),
                (egui::TextStyle::Monospace, 15.0),
                (egui::TextStyle::Button, 16.0),
                (egui::TextStyle::Small, 13.0),
            ] {
                let font = if text_style == egui::TextStyle::Monospace {
                    egui::FontId::monospace(size)
                } else {
                    egui::FontId::proportional(size)
                };
                style.text_styles.insert(text_style, font);
            }
        });
        ctx
    }

    fn label_rect(output: &egui::FullOutput, label: &str) -> egui::Rect {
        output
            .shapes
            .iter()
            .find_map(|shape| match &shape.shape {
                egui::Shape::Text(text) if text.galley.text() == label => {
                    Some(egui::Rect::from_min_size(text.pos, text.galley.size()))
                }
                _ => None,
            })
            .unwrap_or_else(|| panic!("Missing control: {label}"))
    }

    #[test]
    fn editor_uses_remaining_height_and_keeps_actions_visible() {
        for size in [[500.0, 620.0], [400.0, 620.0], [400.0, 470.0]] {
            let ctx = context();
            let mut draft = "1\nDraft\nC 0 0 0\n".to_string();
            let original = draft.clone();
            let mut action = CoordinateAction::None;
            for _ in 0..3 {
                frame(&ctx, size, vec![], &mut draft, &mut action);
            }
            let output = frame(&ctx, size, vec![], &mut draft, &mut action);
            let viewport = output
                .shapes
                .iter()
                .find_map(|shape| match &shape.shape {
                    egui::Shape::Rect(rect)
                        if rect.fill != egui::Color32::TRANSPARENT
                            && rect.rect.width() > 250.0
                            && rect.rect.height() > 200.0 =>
                    {
                        Some(rect.rect)
                    }
                    _ => None,
                })
                .unwrap_or_else(|| {
                    panic!(
                        "Coordinate viewport must be visible: {:?}",
                        output
                            .shapes
                            .iter()
                            .filter_map(|shape| match &shape.shape {
                                egui::Shape::Rect(rect) => Some((rect.rect, rect.fill)),
                                _ => None,
                            })
                            .collect::<Vec<_>>()
                    )
                });
            assert!(
                viewport.height() > size[1] * if size[1] >= 600.0 { 0.60 } else { 0.50 },
                "Coordinates got only {viewport:?}"
            );
            for label in ["Apply coordinates", "Reset", "Clear"] {
                let rect = label_rect(&output, label);
                assert!(
                    rect.top() >= viewport.bottom(),
                    "{label} overlaps the editor"
                );
                assert!(
                    rect.bottom() <= size[1] + 80.0 && rect.bottom() < 720.0,
                    "{label} is outside the window: {rect:?}"
                );
                assert!(
                    rect.right() < size[0] + 80.0,
                    "{label} does not fit horizontally"
                );
            }
            assert_eq!(
                draft, original,
                "Drawing the editor must preserve its draft"
            );
            assert_eq!(action, CoordinateAction::None);

            // The blank area below a small draft must also accept clicks and
            // typing, rather than merely looking like a larger text editor.
            let pos = egui::pos2(viewport.center().x, viewport.bottom() - 60.0);
            for pressed in [true, false] {
                frame(
                    &ctx,
                    size,
                    vec![
                        egui::Event::PointerMoved(pos),
                        egui::Event::PointerButton {
                            pos,
                            button: egui::PointerButton::Primary,
                            pressed,
                            modifiers: egui::Modifiers::NONE,
                        },
                    ],
                    &mut draft,
                    &mut action,
                );
            }
            frame(
                &ctx,
                size,
                vec![egui::Event::Text("#".into())],
                &mut draft,
                &mut action,
            );
            assert!(
                draft.contains('#'),
                "The expanded editor did not accept typing"
            );
            assert_eq!(action, CoordinateAction::None);
        }
    }

    #[test]
    fn apply_remains_reachable_below_a_long_scrolling_draft() {
        let ctx = context();
        let mut draft = format!("300\nLarge draft\n{}", "C 0 0 0\n".repeat(300));
        let original = draft.clone();
        let mut action = CoordinateAction::None;
        for _ in 0..3 {
            frame(&ctx, [500.0, 620.0], vec![], &mut draft, &mut action);
        }
        let output = frame(&ctx, [500.0, 620.0], vec![], &mut draft, &mut action);
        let pos = label_rect(&output, "Apply coordinates").center();
        for pressed in [true, false] {
            frame(
                &ctx,
                [500.0, 620.0],
                vec![
                    egui::Event::PointerMoved(pos),
                    egui::Event::PointerButton {
                        pos,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::NONE,
                    },
                ],
                &mut draft,
                &mut action,
            );
        }
        assert_eq!(action, CoordinateAction::Apply);
        assert_eq!(draft, original);
    }
}
