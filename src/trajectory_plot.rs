//! Simple selectable plots without allocating one widget per trajectory point.
use bevy_egui::egui;
fn limits(values: impl Iterator<Item = f64>) -> (f64, f64) {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for v in values {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    if !lo.is_finite() {
        return (0., 1.);
    }
    if hi - lo < 1e-12 {
        (lo - 0.5, hi + 0.5)
    } else {
        let p = (hi - lo) * 0.06;
        (lo - p, hi + p)
    }
}
fn canvas(ui: &mut egui::Ui) -> (egui::Response, egui::Rect) {
    let (response, _) = ui.allocate_painter(
        egui::vec2(
            ui.available_width().max(200.),
            ui.available_height().max(240.) - 28.,
        ),
        egui::Sense::click(),
    );
    let r = response.rect.shrink2(egui::vec2(58., 25.));
    (response, r)
}
pub fn series(
    ui: &mut egui::Ui,
    values: &[Option<f64>],
    selected: usize,
    unit: &str,
) -> Option<usize> {
    let points:Vec<_>=values.iter().enumerate().map(|(i,v)|v.map(|y|(i as f64+1.,y))).collect();
    crate::analysis_plot::draw(ui,"trajectory_series",crate::analysis_plot::Plot{points:&points,x:"Frame",y:unit,note:"",reverse_x:false,point_ids:None},Some(selected))
}

pub fn histogram(ui: &mut egui::Ui, values: &[Option<f64>], bins: usize, unit: &str) {
    let data: Vec<_> = values.iter().flatten().copied().collect();
    let (lo, hi) = limits(data.iter().copied());
    let mut counts = vec![0usize; bins];
    for v in data {
        let i = (((v - lo) / (hi - lo) * bins as f64) as usize).min(bins - 1);
        counts[i] += 1;
    }
    let (response, r) = canvas(ui);
    let maximum = *counts.iter().max().unwrap_or(&1);
    for (i, n) in counts.iter().enumerate() {
        let left = r.left() + r.width() * i as f32 / bins as f32;
        let width = r.width() / bins as f32;
        let rect = egui::Rect::from_min_max(
            egui::pos2(
                left,
                r.bottom() - r.height() * *n as f32 / maximum.max(1) as f32,
            ),
            egui::pos2(left + width - 1., r.bottom()),
        );
        ui.painter()
            .rect_filled(rect, 0., egui::Color32::from_rgb(85, 180, 235));
        if response.hover_pos().is_some_and(|p| rect.contains(p)) {
            response.clone().on_hover_text(format!(
                "{:.6}–{:.6} {unit}\n{n} frames",
                lo + (hi - lo) * i as f64 / bins as f64,
                lo + (hi - lo) * (i + 1) as f64 / bins as f64
            ));
        }
    }
    ui.painter().text(
        egui::pos2(r.left(), r.bottom() + 5.),
        egui::Align2::LEFT_TOP,
        format!("{lo:.4}"),
        egui::FontId::monospace(12.),
        ui.visuals().text_color(),
    );
    ui.painter().text(
        egui::pos2(r.right(), r.bottom() + 5.),
        egui::Align2::RIGHT_TOP,
        format!("{hi:.4} {unit} · count"),
        egui::FontId::monospace(12.),
        ui.visuals().text_color(),
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn clicking_the_third_point_returns_frame_three_despite_an_energy_gap() {
        let ctx = egui::Context::default();
        let mut selected = None;
        let mut frame = |events| {
            let mut out = ctx.run_ui(
                egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(640., 420.),
                    )),
                    events,
                    ..Default::default()
                },
                |ui| {
                    selected = series(ui, &[Some(1.), None, Some(3.)], 0, "kcal/mol");
                },
            );
            out.textures_delta.clear();
            out
        };
        let out = frame(vec![]);
        let point = out
            .shapes
            .iter()
            .filter_map(|shape| match &shape.shape {
                egui::Shape::Circle(c) if (c.radius - 2.5).abs() < 0.01 => Some(c.center),
                _ => None,
            })
            .last()
            .expect("energy points");
        for pressed in [true, false] {
            frame(vec![
                egui::Event::PointerMoved(point),
                egui::Event::PointerButton {
                    pos: point,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ]);
        }
        assert_eq!(selected, Some(2));
    }
}
