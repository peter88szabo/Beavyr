//! XY plot preserving frame order, including repeated/nonmonotonic coordinates.
use bevy_egui::egui;
pub fn draw(
    ui: &mut egui::Ui,
    points: &[Option<(f64, f64)>],
    selected: Option<usize>,
    x_label: &str,
    unit: &str,
) -> Option<usize> {
    draw_named(
        ui,
        points,
        selected,
        x_label,
        unit,
        "Sum of atomic polarizabilities",
    )
}
pub fn draw_named(
    ui: &mut egui::Ui,
    points: &[Option<(f64, f64)>],
    selected: Option<usize>,
    x_label: &str,
    unit: &str,
    title: &str,
) -> Option<usize> {
    let (response, painter) = ui.allocate_painter(
        egui::vec2(
            ui.available_width().max(240.),
            ui.available_height().max(260.),
        ),
        egui::Sense::click(),
    );
    let r = response.rect.shrink2(egui::vec2(65., 32.));
    let limits = |axis: usize| {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for &(x, y) in points.iter().flatten() {
            let v = if axis == 0 { x } else { y };
            lo = lo.min(v);
            hi = hi.max(v);
        }
        if !lo.is_finite() {
            return (0., 1.);
        }
        let pad = ((hi - lo) * 0.05).max(0.0001);
        (lo - pad, hi + pad)
    };
    let (xl, xh) = limits(0);
    let (yl, yh) = limits(1);
    let color = ui.visuals().text_color();
    let accent = egui::Color32::from_rgb(80, 190, 240);
    for i in 0..=4 {
        let f = i as f32 / 4.;
        let x = r.left() + r.width() * f;
        let y = r.bottom() - r.height() * f;
        painter.line_segment(
            [egui::pos2(r.left(), y), egui::pos2(r.right(), y)],
            egui::Stroke::new(0.5, color.gamma_multiply(0.2)),
        );
        painter.text(
            egui::pos2(r.left() - 5., y),
            egui::Align2::RIGHT_CENTER,
            format!("{:.3}", yl + (yh - yl) * f as f64),
            egui::FontId::monospace(11.),
            color,
        );
        painter.text(
            egui::pos2(x, r.bottom() + 4.),
            egui::Align2::CENTER_TOP,
            format!("{:.3}", xl + (xh - xl) * f as f64),
            egui::FontId::monospace(11.),
            color,
        );
    }
    painter.text(
        egui::pos2(r.center().x, response.rect.bottom()),
        egui::Align2::CENTER_BOTTOM,
        x_label,
        egui::FontId::proportional(13.),
        color,
    );
    painter.text(
        response.rect.left_top(),
        egui::Align2::LEFT_TOP,
        format!("{title} ({unit})"),
        egui::FontId::proportional(13.),
        color,
    );
    let mut previous = None;
    let mut nearest = None;
    let mut distance = 225.;
    for (i, p) in points.iter().enumerate() {
        if let Some((x, y)) = p {
            let p = egui::pos2(
                r.left() + r.width() * ((x - xl) / (xh - xl)) as f32,
                r.bottom() - r.height() * ((y - yl) / (yh - yl)) as f32,
            );
            if let Some(old) = previous {
                painter.line_segment([old, p], egui::Stroke::new(1.5, accent));
            }
            previous = Some(p);
            painter.circle_filled(p, 2.5, accent);
            if selected == Some(i) {
                painter.circle_stroke(p, 5., egui::Stroke::new(2., egui::Color32::YELLOW));
            }
            if let Some(mouse) = response.hover_pos() {
                let d = mouse.distance_sq(p);
                if d < distance {
                    distance = d;
                    nearest = Some((i, *x, *y));
                }
            }
        } else {
            previous = None;
        }
    }
    if let Some((i, x, y)) = nearest {
        response.clone().on_hover_text(format!(
            "Frame {}\n{x_label}: {x:.8}\n{y:.8} {unit}\nClick to display",
            i + 1
        ));
        if response.clicked() {
            return Some(i);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reversed_coordinate_and_missing_point_keep_original_frame_indices() {
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
                    selected = draw_named(
                        ui,
                        &[Some((9., 2.)), None, Some((1., -1.))],
                        None,
                        "Distance (Å)",
                        "dimensionless",
                        "Fukui",
                    );
                },
            );
            out.textures_delta.clear();
            out
        };
        let out = frame(vec![]);
        let points: Vec<_> = out
            .shapes
            .iter()
            .filter_map(|s| match &s.shape {
                egui::Shape::Circle(c) if (c.radius - 2.5).abs() < 0.01 => Some(c.center),
                _ => None,
            })
            .collect();
        assert_eq!(points.len(), 2);
        assert!(points[1].x < points[0].x);
        for pressed in [true, false] {
            frame(vec![
                egui::Event::PointerMoved(points[1]),
                egui::Event::PointerButton {
                    pos: points[1],
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ]);
        }
        assert_eq!(selected, Some(2));
    }
}
