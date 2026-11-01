//! Small interactive plots with numerical hover readouts.
use super::Duschinsky;
use bevy_egui::egui;
pub fn curve(ui: &mut egui::Ui, points: &[(f64, f64)], x_label: &str, y_label: &str, sticks: bool) {
    if points.is_empty() {
        return;
    }
    let (whole, response) = ui.allocate_exact_size(
        egui::vec2(ui.available_width().max(260.), 240.),
        egui::Sense::hover(),
    );
    let rect = egui::Rect::from_min_max(
        whole.min + egui::vec2(55., 22.),
        whole.max - egui::vec2(15., 40.),
    );
    let painter = ui.painter();
    let color = ui.visuals().text_color();
    let font = egui::FontId::proportional(11.);
    let xmin = points.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
    let xmax = points
        .iter()
        .map(|p| p.0)
        .fold(f64::NEG_INFINITY, f64::max)
        .max(xmin + 1e-9);
    let ymax = points.iter().map(|p| p.1).fold(0_f64, f64::max).max(1e-12);
    let pos = |x: f64, y: f64| {
        egui::pos2(
            rect.left() + ((x - xmin) / (xmax - xmin)) as f32 * rect.width(),
            rect.bottom() - (y / ymax) as f32 * rect.height(),
        )
    };
    painter.line_segment([rect.left_top(), rect.left_bottom()], (1., color));
    painter.line_segment([rect.left_bottom(), rect.right_bottom()], (1., color));
    for k in 0..=4 {
        let x = xmin + (xmax - xmin) * k as f64 / 4.;
        let p = pos(x, 0.);
        painter.text(
            p + egui::vec2(0., 5.),
            egui::Align2::CENTER_TOP,
            format!("{x:.2}"),
            font.clone(),
            color,
        );
    }
    painter.text(
        whole.center_bottom(),
        egui::Align2::CENTER_BOTTOM,
        x_label,
        font.clone(),
        color,
    );
    painter.text(
        whole.left_top(),
        egui::Align2::LEFT_TOP,
        y_label,
        font.clone(),
        color,
    );
    painter.text(
        rect.left_top() - egui::vec2(4., 0.),
        egui::Align2::RIGHT_TOP,
        format!("{ymax:.3}"),
        font.clone(),
        color,
    );
    let stroke = egui::Stroke::new(1.4, egui::Color32::from_rgb(105, 185, 250));
    if sticks {
        for &(x, y) in points {
            painter.line_segment([pos(x, 0.), pos(x, y)], stroke);
        }
    } else {
        // Preserve extrema per screen column, avoiding hundreds of thousands
        // of vertices while retaining narrow spectral peaks.
        let columns = rect.width().ceil().max(1.) as usize;
        let mut bins = vec![None::<(f64, f64, f64)>; columns + 1];
        for &(x, y) in points {
            let i = (((x - xmin) / (xmax - xmin)) * columns as f64).round() as usize;
            let b = &mut bins[i.min(columns)];
            *b = Some(match *b {
                Some((x, lo, hi)) => (x, lo.min(y), hi.max(y)),
                None => (x, y, y),
            });
        }
        let vertices = bins
            .into_iter()
            .flatten()
            .flat_map(|(x, lo, hi)| [pos(x, lo), pos(x, hi)])
            .collect::<Vec<_>>();
        painter.add(egui::Shape::line(vertices, stroke));
    }
    if let Some(p) = response.hover_pos().filter(|p| rect.contains(*p)) {
        let x = xmin + (p.x - rect.left()) as f64 / rect.width() as f64 * (xmax - xmin);
        if let Some(&(px, py)) = points
            .iter()
            .min_by(|a, b| (a.0 - x).abs().total_cmp(&(b.0 - x).abs()))
        {
            painter.circle_filled(pos(px, py), 3., color);
            response.on_hover_text(format!("{x_label}: {px:.6}\n{y_label}: {py:.8}"));
        }
    }
}
pub fn texture(ctx: &egui::Context, d: &Duschinsky) -> egui::TextureHandle {
    let (rows, cols) = d.rotation.dim();
    let pixels = d
        .rotation
        .iter()
        .map(|v| {
            let t = (v * v).clamp(0., 1.) as f32;
            egui::Color32::from_rgb(
                (15. + 240. * t) as u8,
                (25. + 165. * t) as u8,
                (55. + 15. * t) as u8,
            )
        })
        .collect();
    ctx.load_texture(
        "duschinsky_overlap",
        egui::ColorImage::new([cols, rows], pixels),
        egui::TextureOptions::NEAREST,
    )
}
pub fn matrix(ui: &mut egui::Ui, d: &Duschinsky, texture: &egui::TextureHandle) {
    ui.label("J² · rows: ground modes · columns: excited modes (dark 0 → bright 1)");
    let side = ui.available_width().clamp(240., 480.);
    let response = ui.add(
        egui::Image::new(texture)
            .fit_to_exact_size(egui::vec2(side, side))
            .sense(egui::Sense::hover()),
    );
    if let Some(p) = response.hover_pos() {
        let row = (((p.y - response.rect.top()) / response.rect.height() * d.rows.len() as f32)
            as usize)
            .min(d.rows.len() - 1);
        let col = (((p.x - response.rect.left()) / response.rect.width() * d.columns.len() as f32)
            as usize)
            .min(d.columns.len() - 1);
        let j = d.rotation[(row, col)];
        response.on_hover_text(format!(
            "Ground mode {}: {:.2} cm⁻¹\nExcited mode {}: {:.2} cm⁻¹\nJ = {j:.8}\nJ² = {:.8}",
            d.rows[row] + 1,
            d.row_frequencies_cm1[row],
            d.columns[col] + 1,
            d.column_frequencies_cm1[col],
            j * j
        ));
    }
}
