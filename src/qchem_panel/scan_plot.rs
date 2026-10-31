//! Interactive scan plot using the same order, units and reference as the table.
use super::*;
fn coordinate(report: &ScanReport, index: usize) -> f64 {
    let p = &report.points[index];
    if report.options.coordinate == Coordinate::Dihedral {
        p.target + (p.actual - p.target + 180.).rem_euclid(360.) - 180.
    } else {
        p.actual
    }
}
fn points(report: &ScanReport, view: ScanView) -> Vec<(usize, [f64; 2])> {
    view.indices(report)
        .into_iter()
        .map(|i| (i, [coordinate(report, i), view.energy(report, i)]))
        .collect()
}
fn hover_text(report: &ScanReport, index: usize, view: ScanView) -> String {
    format!(
        "Point {}\n{}: {:.6} {}\nΔE: {:.6} {}",
        index + 1,
        report.options.coordinate.label(),
        coordinate(report, index),
        report.options.coordinate.unit(),
        view.energy(report, index),
        view.unit.label()
    )
}
pub(super) fn draw(
    ui: &mut egui::Ui,
    report: &ScanReport,
    selected: usize,
    view: ScanView,
) -> Option<usize> {
    let (rect, response) = ui.allocate_exact_size(
        egui::vec2(
            ui.available_width().max(300.),
            ui.available_height().max(220.),
        ),
        egui::Sense::click(),
    );
    let plot = rect.shrink2(egui::vec2(60., 32.));
    let painter = ui.painter_at(rect);
    let text = ui.visuals().text_color();
    let values = points(report, view);
    let xmin = values
        .iter()
        .map(|(_, p)| p[0])
        .fold(f64::INFINITY, f64::min);
    let xmax = values
        .iter()
        .map(|(_, p)| p[0])
        .fold(f64::NEG_INFINITY, f64::max);
    let mut ymin = values
        .iter()
        .map(|(_, p)| p[1])
        .fold(f64::INFINITY, f64::min);
    let mut ymax = values
        .iter()
        .map(|(_, p)| p[1])
        .fold(f64::NEG_INFINITY, f64::max);
    let pad = ((ymax - ymin) * 0.05).max(0.01);
    ymin -= pad;
    ymax += pad;
    let descending = values.first().unwrap().1[0] > values.last().unwrap().1[0];
    let screen = |p: [f64; 2]| {
        let mut t = (p[0] - xmin) / (xmax - xmin).max(1e-9);
        if descending {
            t = 1. - t;
        }
        egui::pos2(
            plot.left() + t as f32 * plot.width(),
            plot.bottom() - ((p[1] - ymin) / (ymax - ymin)) as f32 * plot.height(),
        )
    };
    painter.text(
        rect.left_top(),
        egui::Align2::LEFT_TOP,
        format!("ΔE ({})", view.unit.label()),
        egui::FontId::proportional(13.),
        text,
    );
    painter.text(
        rect.center_bottom(),
        egui::Align2::CENTER_BOTTOM,
        format!(
            "{} ({})",
            report.options.coordinate.label(),
            report.options.coordinate.unit()
        ),
        egui::FontId::proportional(13.),
        text,
    );
    for i in 0..=4 {
        let t = i as f64 / 4.;
        let x = xmin + (xmax - xmin) * t;
        let y = ymin + (ymax - ymin) * t;
        let py = screen([xmin, y]).y;
        painter.line_segment(
            [egui::pos2(plot.left(), py), egui::pos2(plot.right(), py)],
            egui::Stroke::new(1., text.gamma_multiply(0.2)),
        );
        painter.text(
            egui::pos2(plot.left() - 5., py),
            egui::Align2::RIGHT_CENTER,
            format!("{y:.2}"),
            egui::FontId::monospace(11.),
            text,
        );
        painter.text(
            egui::pos2(screen([x, ymin]).x, plot.bottom() + 4.),
            egui::Align2::CENTER_TOP,
            format!("{x:.2}"),
            egui::FontId::monospace(11.),
            text,
        );
    }
    let positions: Vec<_> = values.iter().map(|(i, p)| (*i, screen(*p))).collect();
    for pair in positions.windows(2) {
        painter.line_segment(
            [pair[0].1, pair[1].1],
            egui::Stroke::new(1.5, egui::Color32::LIGHT_BLUE),
        );
    }
    for &(i, p) in &positions {
        painter.circle_filled(
            p,
            if i == selected { 5. } else { 3.5 },
            if report.points[i].converged == Some(false) {
                egui::Color32::LIGHT_RED
            } else {
                egui::Color32::LIGHT_BLUE
            },
        );
    }
    let nearest = response.hover_pos().and_then(|p| {
        positions
            .iter()
            .min_by(|(_, a), (_, b)| a.distance_sq(p).total_cmp(&b.distance_sq(p)))
            .filter(|(_, q)| q.distance_sq(p) <= 144.)
            .map(|(i, _)| *i)
    });
    let clicked = response.clicked_once();
    if let Some(i) = nearest {
        response.on_hover_text(hover_text(report, i, view));
    }
    if clicked {
        nearest
    } else {
        None
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn scan_plot_values_hover_and_order_match_display_settings() {
        let report = super::super::tests::output().scan.unwrap();
        let mut view = ScanView::default();
        assert_eq!(points(&report, view)[0].0, 0);
        assert!(hover_text(&report, 1, view).contains("kcal/mol"));
        view.unit = EnergyUnit::Kj;
        view.reference = Reference::Last;
        view.reversed = true;
        let p = points(&report, view);
        assert_eq!(p[0].0, 2);
        assert_eq!(p.last().unwrap().0, 0);
        assert_eq!(p.last().unwrap().1[1], 0.);
        assert!(hover_text(&report, 2, view).contains("kJ/mol"));
    }
}
