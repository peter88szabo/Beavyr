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
    let y=format!("{title} ({unit})");
    crate::analysis_plot::draw(ui,"atomic_property_plot",crate::analysis_plot::Plot{points,x:x_label,y:&y,note:"",reverse_x:false,point_ids:None},selected)
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
