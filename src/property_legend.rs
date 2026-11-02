//! Fixed, trajectory-wide colour ranges. Display units do not alter physical bounds.
use bevy::prelude::Color;
use bevy_egui::egui;
#[derive(Clone)]
pub struct Scale {
    pub show: bool,
    pub manual: bool,
    pub bounds: [f64; 2],
    auto: [f64; 2],
    key: Option<(usize, usize)>,
    pub signed: bool,
}
impl Default for Scale {
    fn default() -> Self {
        Self {
            show: true,
            manual: false,
            bounds: [0., 1.],
            auto: [0., 1.],
            key: None,
            signed: false,
        }
    }
}
impl Scale {
    pub fn update(&mut self, key: (usize, usize), signed: bool, values: impl FnOnce() -> Vec<f64>) {
        if self.key == Some(key) {
            return;
        }
        self.key = Some(key);
        self.signed = signed;
        let max = values()
            .into_iter()
            .filter(|v| v.is_finite())
            .fold(1e-12_f64, |m, v| m.max(v.abs()));
        self.auto = [if signed { -max } else { 0. }, max];
        if !self.manual {
            self.bounds = self.auto;
        }
    }
    pub fn controls(&mut self, ui: &mut egui::Ui, factor: f64, unit: &str) {
        ui.horizontal_wrapped(|ui| {
            ui.checkbox(&mut self.show, "Colour legend");
            if ui.checkbox(&mut self.manual, "Custom range").changed() && !self.manual {
                self.bounds = self.auto;
            }
        });
        if self.manual {
            ui.horizontal(|ui| {
                let mut lo = self.bounds[0] * factor;
                let mut hi = self.bounds[1] * factor;
                ui.label("Min");
                ui.add(egui::DragValue::new(&mut lo).speed(0.01));
                ui.label("Max");
                ui.add(egui::DragValue::new(&mut hi).speed(0.01));
                ui.label(unit);
                if lo.is_finite()
                    && hi.is_finite()
                    && lo < hi
                    && (!self.signed || (lo < 0. && hi > 0.))
                {
                    self.bounds = [lo / factor, hi / factor];
                }
            });
        } else {
            ui.small("Colour range: all frames");
        }
    }
    pub fn color(&self, value: f64, positive: Color, negative: Color, opacity: f32) -> Color {
        let c = if self.signed && value < 0. {
            negative
        } else {
            positive
        }
        .to_srgba();
        let t = if self.signed {
            if value < 0. {
                value / self.bounds[0].min(-1e-12)
            } else {
                value / self.bounds[1].max(1e-12)
            }
        } else {
            (value - self.bounds[0]) / (self.bounds[1] - self.bounds[0]).max(1e-12)
        }
        .clamp(0., 1.) as f32;
        Color::srgba(
            0.88 + (c.red - 0.88) * t,
            0.88 + (c.green - 0.88) * t,
            0.88 + (c.blue - 0.88) * t,
            opacity,
        )
    }
    pub fn overlay(
        &self,
        ctx: &egui::Context,
        id: &str,
        title: &str,
        factor: f64,
        unit: &str,
        colors: [Color; 2],
        offset: f32,
    ) {
        if !self.show {
            return;
        }
        egui::Area::new(egui::Id::new(id))
            .anchor(egui::Align2::RIGHT_BOTTOM, [-24., -100. - offset])
            .interactable(false)
            .show(ctx, |ui| {
                egui::Frame::popup(ui.style()).show(ui, |ui| {
                    ui.label(format!("{title} ({unit})"));
                    let (rect, _) =
                        ui.allocate_exact_size(egui::vec2(200., 12.), egui::Sense::hover());
                    for i in 0..100 {
                        let t = i as f64 / 99.;
                        let value = self.bounds[0] + t * (self.bounds[1] - self.bounds[0]);
                        let c = self.color(value, colors[0], colors[1], 1.).to_srgba();
                        ui.painter().rect_filled(
                            egui::Rect::from_min_max(
                                egui::pos2(rect.left() + i as f32 * 2., rect.top()),
                                egui::pos2(rect.left() + (i + 1) as f32 * 2., rect.bottom()),
                            ),
                            0.,
                            egui::Color32::from_rgb(
                                (c.red * 255.) as u8,
                                (c.green * 255.) as u8,
                                (c.blue * 255.) as u8,
                            ),
                        );
                    }
                    ui.horizontal(|ui| {
                        ui.label(format!("{:.4}", self.bounds[0] * factor));
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            ui.label(format!("{:.4}", self.bounds[1] * factor));
                        });
                    });
                });
            });
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fixed_range_and_signed_mapping() {
        let mut s = Scale::default();
        s.update((1, 0), true, || vec![-2., 1.]);
        assert_eq!(s.bounds, [-2., 2.]);
        s.update((1, 0), true, || panic!("frame change must not rescale"));
        s.manual = true;
        s.bounds = [-1., 3.];
        s.update((2, 0), true, || vec![20.]);
        assert_eq!(s.bounds, [-1., 3.]);
        let c = s
            .color(0., Color::srgb(0., 0., 1.), Color::srgb(1., 0., 0.), 0.5)
            .to_srgba();
        assert_eq!(c.red, c.blue);
        assert_eq!(c.alpha, 0.5);
    }
}
