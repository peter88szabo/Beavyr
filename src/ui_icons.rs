//! Toolbar icons drawn as vectors, independent of font and emoji support.

use bevy_egui::egui;
use std::f32::consts::{PI, TAU};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Icon {
    Structure,
    History,
    Representation,
    Optimize,
    TransitionState,
    Conformers,
    Vibrations,
    Spectrum,
    Trajectory,
    Compare,
    Measure,
    Surface,
    Diagnostics,
    Appearance,
    Export,
    QuantumChemistry,
    MoleculeEditor,
}

impl Icon {
    pub fn label(self) -> &'static str {
        match self {
            Self::Structure => "Structure",
            Self::History => "Recent structures",
            Self::Representation => "Representation",
            Self::Optimize => "Geometry optimization",
            Self::TransitionState => "Transition-state generation",
            Self::Conformers => "Conformer search",
            Self::Vibrations => "Vibrations",
            Self::Spectrum => "UV-Vis spectrum",
            Self::Trajectory => "Trajectory",
            Self::Compare => "Compare structures",
            Self::Measure => "Measurements",
            Self::Surface => "Orbital and density surfaces",
            Self::Diagnostics => "Diagnostics",
            Self::Appearance => "Appearance",
            Self::Export => "Export image",
            Self::QuantumChemistry => "Quantum chemistry",
            Self::MoleculeEditor => "Molecule editor",
        }
    }

    fn paint(self, painter: &egui::Painter, rect: egui::Rect, color: egui::Color32) {
        let point = |x: f32, y: f32| rect.min + egui::vec2(x, y) * (rect.width() / 24.0);
        let scale = rect.width() / 24.0;
        let stroke = egui::Stroke::new(1.7 * scale.max(0.8), color);
        let line = |a: [f32; 2], b: [f32; 2]| {
            painter.line_segment([point(a[0], a[1]), point(b[0], b[1])], stroke)
        };
        let path = |points: &[[f32; 2]]| {
            painter.add(egui::Shape::line(
                points.iter().map(|p| point(p[0], p[1])).collect(),
                stroke,
            ))
        };
        let circle = |x: f32, y: f32, radius: f32| {
            painter.circle_stroke(point(x, y), radius * scale, stroke)
        };
        let ellipse = |x: f32, y: f32, rx: f32, ry: f32, angle: f32| {
            let points = (0..=40)
                .map(|index| {
                    let theta = TAU * index as f32 / 40.0;
                    let dx = rx * theta.cos();
                    let dy = ry * theta.sin();
                    point(
                        x + dx * angle.cos() - dy * angle.sin(),
                        y + dx * angle.sin() + dy * angle.cos(),
                    )
                })
                .collect();
            painter.add(egui::Shape::line(points, stroke));
        };
        match self {
            Self::Structure => {
                path(&[
                    [6., 3.],
                    [14., 3.],
                    [19., 8.],
                    [19., 21.],
                    [6., 21.],
                    [6., 3.],
                ]);
                path(&[[14., 3.], [14., 8.], [19., 8.]]);
                for y in [12., 16.] {
                    line([9., y], [16., y]);
                }
            }
            Self::History => {
                let points = (0..=32)
                    .map(|i| {
                        let angle = -PI * 0.75 + i as f32 * PI * 1.75 / 32.0;
                        point(12. + 8. * angle.cos(), 12. + 8. * angle.sin())
                    })
                    .collect();
                painter.add(egui::Shape::line(points, stroke));
                path(&[[3., 4.], [3., 9.], [8., 9.]]);
                path(&[[12., 7.], [12., 12.], [16., 14.]]);
            }
            Self::Representation => {
                line([6., 14.], [11., 7.]);
                line([14., 7.], [18., 12.]);
                circle(5., 17., 3.);
                circle(12., 4., 3.);
                circle(20., 15., 3.);
            }
            Self::Optimize => {
                path(&[
                    [15., 3.],
                    [13., 7.],
                    [17., 11.],
                    [21., 9.],
                    [20., 14.],
                    [16., 15.],
                    [7., 22.],
                    [3., 18.],
                    [11., 10.],
                    [11., 6.],
                    [15., 3.],
                ]);
            }
            Self::TransitionState => {
                path(&[
                    [2., 20.],
                    [6., 18.],
                    [9., 10.],
                    [12., 5.],
                    [15., 10.],
                    [18., 18.],
                    [22., 20.],
                ]);
                path(&[[2., 6.], [7., 6.], [5., 4.]]);
                line([7., 6.], [5., 8.]);
                path(&[[22., 6.], [17., 6.], [19., 4.]]);
                line([17., 6.], [19., 8.]);
            }
            Self::Conformers => {
                for (x, y) in [(8., 9.), (16., 15.)] {
                    let points = (0..=6)
                        .map(|i| {
                            let angle = PI / 3.0 * i as f32;
                            point(x + 5.5 * angle.cos(), y + 5.5 * angle.sin())
                        })
                        .collect();
                    painter.add(egui::Shape::line(points, stroke));
                }
            }
            Self::Vibrations => {
                let points = (0..=40)
                    .map(|i| {
                        let t = i as f32 / 40.0;
                        point(2. + 20. * t, 12. - 6. * (t * TAU).sin())
                    })
                    .collect();
                painter.add(egui::Shape::line(points, stroke));
                painter.circle_filled(point(2., 12.), 1.4 * scale, color);
                painter.circle_filled(point(22., 12.), 1.4 * scale, color);
            }
            Self::Spectrum => {
                path(&[[3., 3.], [3., 21.], [22., 21.]]);
                path(&[
                    [5., 18.],
                    [7., 16.],
                    [9., 10.],
                    [11., 18.],
                    [14., 16.],
                    [17., 4.],
                    [20., 18.],
                    [22., 18.],
                ]);
            }
            Self::Trajectory => {
                painter.add(egui::Shape::convex_polygon(
                    vec![point(7., 4.), point(21., 12.), point(7., 20.)],
                    color.gamma_multiply(0.18),
                    stroke,
                ));
            }
            Self::Compare => {
                for y in [5., 12., 19.] {
                    circle(5., y, 2.);
                    circle(19., y, 2.);
                    line([9., y], [11., y]);
                    line([13., y], [15., y]);
                }
                line([5., 7.], [5., 10.]);
                line([5., 14.], [5., 17.]);
                line([19., 7.], [19., 10.]);
                line([19., 14.], [19., 17.]);
            }
            Self::Measure => {
                path(&[[3., 17.], [17., 3.], [21., 7.], [7., 21.], [3., 17.]]);
                for (x, y, length) in [(6., 14., 2.), (9., 11., 3.), (12., 8., 2.), (15., 5., 3.)] {
                    line([x, y], [x + length, y + length]);
                }
            }
            Self::Surface => {
                ellipse(7.5, 7.5, 6., 3.5, PI / 4.0);
                ellipse(16.5, 16.5, 6., 3.5, PI / 4.0);
                line([5.5, 7.5], [9.5, 7.5]);
                line([7.5, 5.5], [7.5, 9.5]);
                line([14.5, 16.5], [18.5, 16.5]);
            }
            Self::Diagnostics => {
                path(&[[12., 3.], [22., 21.], [2., 21.], [12., 3.]]);
                line([12., 9.], [12., 14.]);
                painter.circle_filled(point(12., 17.5), scale, color);
            }
            Self::Appearance => {
                for (x, y) in [(8., 5.), (16., 12.), (10., 19.)] {
                    line([3., y], [x - 2.5, y]);
                    line([x + 2.5, y], [21., y]);
                    circle(x, y, 2.5);
                }
            }
            Self::Export => {
                path(&[[15., 6.], [3., 6.], [3., 21.], [17., 21.], [17., 17.]]);
                path(&[[4., 18.], [8., 13.], [12., 17.], [15., 15.]]);
                circle(8., 10., 1.);
                line([20., 15.], [20., 3.]);
                path(&[[17., 6.], [20., 3.], [23., 6.]]);
            }
            Self::QuantumChemistry => {
                for angle in [0., PI / 3., PI * 2. / 3.] {
                    ellipse(12., 12., 10., 4., angle);
                }
                painter.circle_filled(point(12., 12.), 1.8 * scale, color);
            }
            Self::MoleculeEditor => {
                line([5., 5.], [10., 4.]);
                line([5., 5.], [6., 10.]);
                for (x, y) in [(3., 5.), (12., 4.), (6., 12.)] {
                    circle(x, y, 2.);
                }
                path(&[
                    [9., 17.],
                    [18., 8.],
                    [22., 12.],
                    [13., 21.],
                    [9., 22.],
                    [9., 17.],
                ]);
                line([16., 10.], [20., 14.]);
            }
        }
    }
}

pub struct IconButton {
    icon: Icon,
    label: &'static str,
    selected: bool,
    size: f32,
}

impl IconButton {
    pub fn new(icon: Icon) -> Self {
        Self {
            icon,
            label: icon.label(),
            selected: false,
            size: 34.0,
        }
    }
    pub fn selected(mut self, selected: bool) -> Self {
        self.selected = selected;
        self
    }
    pub fn label(mut self, label: &'static str) -> Self {
        self.label = label;
        self
    }
    pub fn size(mut self, size: f32) -> Self {
        self.size = size;
        self
    }
}

impl egui::Widget for IconButton {
    fn ui(self, ui: &mut egui::Ui) -> egui::Response {
        let (rect, response) =
            ui.allocate_exact_size(egui::Vec2::splat(self.size), egui::Sense::click());
        response.widget_info(|| {
            egui::WidgetInfo::labeled(egui::WidgetType::Button, ui.is_enabled(), self.label)
        });
        if ui.is_rect_visible(rect) {
            let visuals = ui.style().interact(&response);
            let accent = egui::Color32::from_rgb(104, 174, 255);
            let fill = if self.selected {
                egui::Color32::from_rgb(35, 75, 117)
            } else if response.hovered() || response.has_focus() {
                visuals.bg_fill
            } else {
                ui.visuals().faint_bg_color
            };
            ui.painter().rect_filled(rect, 7.0, fill);
            if self.selected || response.has_focus() {
                ui.painter().rect_stroke(
                    rect,
                    7.0,
                    egui::Stroke::new(1.0, accent),
                    egui::StrokeKind::Inside,
                );
            }
            let color = if self.selected {
                accent
            } else {
                visuals.fg_stroke.color
            };
            self.icon
                .paint(ui.painter(), rect.shrink((self.size - 22.0) / 2.0), color);
        }
        response
    }
}
