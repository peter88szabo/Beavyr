//! Presentation controls and diagnostics, independent of the renderer's work.

use bevy_egui::egui;

use crate::diagnostics::{DiagnosticSeverity, DiagnosticsReport};
use crate::settings::{MolSettings, RepresentationMode};
use crate::ui_style::{self, ResponseExt};

pub fn representation_panel(ui: &mut egui::Ui, settings: &mut MolSettings, atom_count: usize) {
    ui_style::heading(
        ui,
        "Choose how atoms are drawn",
        "Change the display without changing the structure.",
    );
    let mut changed = false;
    let mut mesh_changed = false;
    ui_style::group(ui, "Display style", |ui| {
        let width = ((ui.available_width() - 8.0) / 2.0).max(100.0);
        let choices = [
            (
                RepresentationMode::BallAndStick,
                "Ball and stick",
                "Best for small molecules and clusters",
            ),
            (
                RepresentationMode::SpaceFilling,
                "Space filling",
                "CPK spheres showing atomic van der Waals sizes",
            ),
            (
                RepresentationMode::SticksRounded,
                "Sticks",
                "Bonds with rounded joints",
            ),
            (
                RepresentationMode::LowResBallsAndLines,
                "Balls and lines",
                "Low-resolution spheres for large structures",
            ),
            (
                RepresentationMode::LinesOnly,
                "Lines",
                "Lightweight bond lines",
            ),
            (
                RepresentationMode::BackboneTrace,
                "Backbone / trace",
                "Trace through non-hydrogen atoms",
            ),
        ];
        for row in choices.chunks(2) {
            ui.horizontal(|ui| {
                for &(mode, label, hint) in row {
                    if ui
                        .add_sized(
                            [width, 36.0],
                            egui::Button::new(label).selected(settings.representation == mode),
                        )
                        .on_hover_text(hint)
                        .clicked_once()
                    {
                        settings.representation = mode;
                        changed = true;
                    }
                }
            });
        }
        if atom_count > 1000 {
            ui.small("Large structures automatically use a lighter representation when needed.");
        }
    });
    ui_style::group(ui, "Size and detail", |ui| {
        match settings.representation {
            RepresentationMode::BallAndStick => {
                ui.small("Atom and bond sizes are available in Appearance.");
            }
            RepresentationMode::SpaceFilling => {
                mesh_changed |= ui
                    .add(
                        egui::Slider::new(&mut settings.cpk_atom_resolution, 0..=10)
                            .text("Sphere detail"),
                    )
                    .changed();
                mesh_changed |= ui
                    .add(
                        egui::Slider::new(&mut settings.cpk_atom_scale, 0.6..=2.6)
                            .text("Sphere size"),
                    )
                    .changed();
            }
            RepresentationMode::SticksRounded => {
                mesh_changed |= ui
                    .add(
                        egui::Slider::new(&mut settings.stick_radius, 0.02..=0.40)
                            .text("Stick thickness"),
                    )
                    .changed();
            }
            RepresentationMode::LowResBallsAndLines => {
                mesh_changed |= ui
                    .add(
                        egui::Slider::new(&mut settings.low_res_atom_resolution, 0..=6)
                            .text("Sphere detail"),
                    )
                    .changed();
                mesh_changed |= ui
                    .add(
                        egui::Slider::new(&mut settings.low_res_atom_scale, 0.2..=1.5)
                            .text("Sphere size"),
                    )
                    .changed();
            }
            RepresentationMode::BackboneTrace => {
                changed |= ui
                    .checkbox(&mut settings.trace_show_atoms, "Show trace atoms")
                    .changed();
                if settings.trace_show_atoms {
                    mesh_changed |= ui
                        .add(
                            egui::Slider::new(&mut settings.low_res_atom_resolution, 0..=6)
                                .text("Sphere detail"),
                        )
                        .changed();
                    mesh_changed |= ui
                        .add(
                            egui::Slider::new(&mut settings.trace_atom_scale, 0.2..=1.5)
                                .text("Sphere size"),
                        )
                        .changed();
                }
            }
            RepresentationMode::LinesOnly => {}
        }
        if matches!(
            settings.representation,
            RepresentationMode::LowResBallsAndLines
                | RepresentationMode::LinesOnly
                | RepresentationMode::BackboneTrace
        ) {
            ui.add(
                egui::Slider::new(&mut settings.line_bond_thickness, 1.0..=50.0)
                    .text("Line thickness"),
            );
        }
    });
    settings.geometry_dirty |= changed;
    settings.meshes_dirty |= mesh_changed;
}

pub fn diagnostics_panel(ui: &mut egui::Ui, report: &DiagnosticsReport, have_molecule: bool) {
    ui_style::heading(
        ui,
        "Check the structure",
        "Valence, charge, isolated atoms and close contacts.",
    );
    if !have_molecule {
        ui_style::group(ui, "Load a structure to inspect", |ui| {
            ui.small("The report updates when the displayed molecule changes.");
        });
        return;
    }
    if report.items.is_empty() {
        ui_style::group(ui, "No issues detected", |ui| {
            ui.colored_label(
                egui::Color32::from_rgb(100, 195, 145),
                "The current structure passes these checks.",
            );
        });
        return;
    }
    ui.horizontal_wrapped(|ui| {
        ui.colored_label(
            egui::Color32::from_rgb(255, 170, 80),
            format!("{} warnings", report.warning_count()),
        );
        ui.weak(format!("{} notes", report.info_count()));
    });
    for (index, item) in report.items.iter().enumerate() {
        ui.push_id(index, |ui| {
            ui_style::group(ui, "", |ui| {
                let (label, color) = match item.severity {
                    DiagnosticSeverity::Warning => {
                        ("Warning", egui::Color32::from_rgb(255, 170, 80))
                    }
                    DiagnosticSeverity::Info => ("Note", ui_style::ACCENT),
                };
                ui.horizontal(|ui| {
                    ui.colored_label(color, egui::RichText::new(label).small().strong());
                    if let Some(atom) = item.atom {
                        ui.weak(format!("Atom {}", atom + 1));
                    }
                });
                ui.add(egui::Label::new(&item.message).wrap());
            });
        });
    }
    ui.small("Diagnostics use perceived single-bond connectivity. Multiple bonds and formal charges may need chemical review.");
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::Molecule;
    use crate::trajectory::{TrajectoryFrame, TrajectoryState};

    fn context() -> egui::Context {
        let ctx = egui::Context::default();
        ctx.style_mut_of(ctx.theme(), |style| {
            style.text_styles = [
                (egui::TextStyle::Heading, egui::FontId::proportional(22.0)),
                (egui::TextStyle::Body, egui::FontId::proportional(16.0)),
                (egui::TextStyle::Monospace, egui::FontId::monospace(15.0)),
                (egui::TextStyle::Button, egui::FontId::proportional(16.0)),
                (egui::TextStyle::Small, egui::FontId::proportional(13.0)),
            ]
            .into();
        });
        ctx
    }
    use bevy::prelude::Vec3;

    fn frame(
        ctx: &egui::Context,
        width: f32,
        events: Vec<egui::Event>,
        mut draw: impl FnMut(&mut egui::Ui),
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
                egui::Window::new("Tool layout check")
                    .fixed_pos([16.0, 16.0])
                    .fixed_size([width, 620.0])
                    .frame(ui_style::window_frame(ui.ctx()))
                    .vscroll(true)
                    .show(ui.ctx(), |ui| {
                        ui_style::modern(ui);
                        draw(ui);
                    });
            },
        );
        output.textures_delta.clear();
        output
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

    fn click(
        ctx: &egui::Context,
        width: f32,
        pos: egui::Pos2,
        mut draw: impl FnMut(&mut egui::Ui),
    ) {
        for pressed in [true, false] {
            frame(
                ctx,
                width,
                vec![
                    egui::Event::PointerMoved(pos),
                    egui::Event::PointerButton {
                        pos,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::NONE,
                    },
                ],
                &mut draw,
            );
        }
    }

    #[test]
    fn representation_choices_fit_and_keep_renderer_updates() {
        for width in [400.0, 500.0] {
            let ctx = context();
            let mut settings = MolSettings::default();
            settings.geometry_dirty = false;
            settings.meshes_dirty = false;
            for _ in 0..2 {
                frame(&ctx, width, vec![], |ui| {
                    representation_panel(ui, &mut settings, 8)
                });
            }
            let output = frame(&ctx, width, vec![], |ui| {
                representation_panel(ui, &mut settings, 8)
            });
            for label in [
                "Ball and stick",
                "Space filling",
                "Sticks",
                "Balls and lines",
                "Lines",
                "Backbone / trace",
            ] {
                let rect = label_rect(&output, label);
                assert!(
                    rect.right() <= width + 40.0 && rect.bottom() < 650.0,
                    "{label} is outside the tool"
                );
            }
            assert!(
                !settings.geometry_dirty && !settings.meshes_dirty,
                "Displaying choices must not rebuild geometry"
            );
            click(&ctx, width, label_rect(&output, "Sticks").center(), |ui| {
                representation_panel(ui, &mut settings, 8)
            });
            assert_eq!(settings.representation, RepresentationMode::SticksRounded);
            assert!(settings.geometry_dirty);
            assert!(!settings.meshes_dirty);
        }
    }

    #[test]
    fn measurement_actions_are_visible_and_picking_stays_exclusive() {
        let ctx = context();
        let mol = Molecule::from_xyz("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n");
        let settings = MolSettings::default();
        let mut measurements = crate::measurements::Measurements::default();
        measurements.is_active = true;
        measurements.pending = Some(0);
        for _ in 0..2 {
            frame(&ctx, 400.0, vec![], |ui| {
                crate::ui_measurements::measurements_panel(ui, &mut measurements, &mol, &settings)
            });
        }
        let output = frame(&ctx, 400.0, vec![], |ui| {
            crate::ui_measurements::measurements_panel(ui, &mut measurements, &mol, &settings)
        });
        for label in ["New distance", "New angle", "New dihedral"] {
            let rect = label_rect(&output, label);
            assert!(
                rect.right() < 440.0 && rect.bottom() < 700.0,
                "{label} should be immediately available"
            );
        }
        click(
            &ctx,
            400.0,
            label_rect(&output, "New angle").center(),
            |ui| crate::ui_measurements::measurements_panel(ui, &mut measurements, &mol, &settings),
        );
        assert!(measurements.angle_active);
        assert!(!measurements.is_active && !measurements.dihedral_active);
        assert!(measurements.pending.is_none());
    }

    #[test]
    fn comparison_reports_differences_without_moving_the_structure() {
        let ctx = context();
        let mut mol = Molecule::from_xyz("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n");
        let reference = mol.pos.clone();
        for pos in &mut mol.pos {
            *pos += Vec3::new(3.0, 1.0, 0.0);
        }
        let original = mol.pos.clone();
        let mut traj = TrajectoryState::default();
        traj.frames = vec![
            TrajectoryFrame {
                atoms: mol.atoms.clone(),
                pos: reference,
            },
            TrajectoryFrame {
                atoms: mol.atoms.clone(),
                pos: original.clone(),
            },
        ];
        traj.current_frame = 1;
        let mut state = crate::rmsd::RmsdState::default();
        for _ in 0..2 {
            frame(&ctx, 400.0, vec![], |ui| {
                crate::rmsd::ui::rmsd_panel(ui, &mut state, &mut traj, &mut mol)
            });
        }
        let output = frame(&ctx, 400.0, vec![], |ui| {
            crate::rmsd::ui::rmsd_panel(ui, &mut state, &mut traj, &mut mol)
        });
        click(
            &ctx,
            400.0,
            label_rect(&output, "Compare current frame").center(),
            |ui| crate::rmsd::ui::rmsd_panel(ui, &mut state, &mut traj, &mut mol),
        );
        let report = state
            .last_report
            .as_ref()
            .expect("comparison button should remain connected");
        assert!(report.raw_rmsd > 3.0 && report.aligned_rmsd < 1.0e-5);
        assert_eq!(mol.pos, original);
        assert_eq!(traj.frames[1].pos, original);
    }
}
