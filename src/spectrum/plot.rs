//! The shared spectrum painter: a white plot with gridlines, a broadened
//! curve or a stick spectrum, peak markers and a hover readout.
//!
//! The IR and UV-Vis windows differ only in units, axis range and the wording
//! of the hover box, all of which are parameters here.

use bevy_egui::egui;

use super::broadening::{intensity_at, sample_spectrum, BroadeningKind};

const AXIS_LABEL_SIZE: f32 = 13.0;
const AXIS_TITLE_SIZE: f32 = 16.0;
/// How close the pointer must come to a peak marker before its readout shows.
const HOVER_RADIUS: f32 = 14.0;
/// Points sampled across the plot width for the continuous curve. Enough that
/// a narrow line shape still renders smoothly at typical window sizes.
const CURVE_POINTS: usize = 400;

/// The horizontal axis of a spectrum plot.
pub struct SpectrumAxis<'a> {
    pub x_min: f64,
    pub x_max: f64,
    /// `None` picks a round step automatically. The IR spectrum pins this to
    /// 500 cm⁻¹ because that is the convention its readers expect.
    pub tick_step: Option<f64>,
    pub tick_decimals: usize,
    pub title: &'a str,
}

/// A round tick step -- 1, 2 or 5 times a power of ten -- giving roughly six
/// to ten intervals across `span`. Round numbers matter more than an exact
/// tick count: readers interpolate between labelled gridlines, which only
/// works if the labels are values they can do arithmetic with.
pub fn nice_tick_step(span: f64) -> f64 {
    const TARGET_INTERVALS: f64 = 8.0;
    if !span.is_finite() || span <= 0.0 {
        return 1.0;
    }
    let raw = span / TARGET_INTERVALS;
    let magnitude = 10f64.powf(raw.log10().floor());
    let normalized = raw / magnitude;
    let step = if normalized <= 1.5 {
        1.0
    } else if normalized <= 3.5 {
        2.0
    } else if normalized <= 7.5 {
        5.0
    } else {
        10.0
    };
    step * magnitude
}

/// The tick positions across an axis: every multiple of `step` strictly
/// inside the range. A tick exactly on `x_min` is skipped -- it would draw a
/// gridline on top of the plot's own border.
pub fn tick_positions(x_min: f64, x_max: f64, step: f64) -> Vec<f64> {
    if step <= 0.0 || !(x_max > x_min) {
        return Vec::new();
    }
    let epsilon = step * 1.0e-6;
    let mut ticks = Vec::new();
    let mut k = (x_min / step).ceil();
    loop {
        let tick = k * step;
        if tick > x_max + epsilon {
            break;
        }
        if tick > x_min + epsilon {
            ticks.push(tick);
        }
        k += 1.0;
        // A degenerate step relative to the range could otherwise spin
        // forever painting invisible gridlines.
        if ticks.len() > 1000 {
            break;
        }
    }
    ticks
}

/// Draws the spectrum, filling the remaining space in `ui`.
///
/// `peaks` are `(position, intensity)` pairs already expressed in the axis's
/// own units. `hover_label` is asked for the readout text of the peak the
/// pointer is nearest, by index into `peaks`.
pub fn draw_spectrum(
    ui: &mut egui::Ui,
    peaks: &[(f64, f64)],
    kind: BroadeningKind,
    width: f64,
    axis: &SpectrumAxis,
    hover_label: &dyn Fn(usize) -> String,
) {
    let x_min = axis.x_min;
    let x_max = axis.x_max.max(axis.x_min + 1.0e-6);
    let is_sticks = kind.is_sticks();

    // Sticks have no continuous curve to sample; the plot's y-scale then comes
    // straight from the raw intensities instead of a broadened curve's peak.
    let curve = if is_sticks {
        Vec::new()
    } else {
        sample_spectrum(peaks, kind, width, x_min, x_max, CURVE_POINTS)
    };
    let max_y = if is_sticks {
        peaks.iter().map(|&(_, y)| y).fold(0.0f64, f64::max)
    } else {
        curve.iter().map(|&(_, y)| y).fold(0.0f64, f64::max)
    }
    .max(1.0e-9);

    let desired_size = egui::vec2(
        ui.available_width().max(200.0),
        ui.available_height().max(120.0),
    );
    let (rect, response) = ui.allocate_exact_size(desired_size, egui::Sense::hover());
    let painter = ui.painter_at(rect);
    // Fixed white background regardless of the app's own theme, with dark
    // text/gridlines to stay legible against it -- a consistent look
    // independent of light/dark mode.
    painter.rect_filled(rect, 0.0, egui::Color32::WHITE);
    let plot_color = egui::Color32::from_rgb(30, 30, 30);
    let grid_color = plot_color.gamma_multiply(0.35);
    let line_color = egui::Color32::from_rgb(40, 100, 170);

    let plot_rect = egui::Rect::from_min_max(
        rect.min + egui::vec2(50.0, 10.0),
        rect.max - egui::vec2(10.0, 46.0),
    );

    let to_screen = |x: f64, y: f64| -> egui::Pos2 {
        let sx = plot_rect.left() + ((x - x_min) / (x_max - x_min)) as f32 * plot_rect.width();
        let sy = plot_rect.bottom() - (y / max_y) as f32 * plot_rect.height();
        egui::pos2(sx, sy)
    };

    let step = axis.tick_step.unwrap_or_else(|| nice_tick_step(x_max - x_min));
    for tick in tick_positions(x_min, x_max, step) {
        let sx = plot_rect.left() + ((tick - x_min) / (x_max - x_min)) as f32 * plot_rect.width();
        painter.line_segment(
            [egui::pos2(sx, plot_rect.top()), egui::pos2(sx, plot_rect.bottom())],
            egui::Stroke::new(1.0, grid_color),
        );
        painter.text(
            egui::pos2(sx, plot_rect.bottom() + 6.0),
            egui::Align2::CENTER_TOP,
            format!("{tick:.*}", axis.tick_decimals),
            egui::FontId::monospace(AXIS_LABEL_SIZE),
            plot_color,
        );
    }
    painter.rect_stroke(
        plot_rect,
        0.0,
        egui::Stroke::new(1.5, plot_color.gamma_multiply(0.6)),
        egui::StrokeKind::Outside,
    );

    let peak_points: Vec<egui::Pos2>;
    if is_sticks {
        // One vertical bar per predicted peak, from the baseline up to its own
        // raw intensity -- the theoretical stick spectrum, not a simulated
        // continuous one.
        peak_points = peaks.iter().map(|&(x, y)| to_screen(x, y)).collect();
        let baseline_y = plot_rect.bottom();
        for &p in &peak_points {
            painter.line_segment(
                [egui::pos2(p.x, baseline_y), p],
                egui::Stroke::new(2.0, line_color),
            );
        }
    } else {
        let points: Vec<egui::Pos2> = curve.iter().map(|&(x, y)| to_screen(x, y)).collect();
        if points.len() >= 2 {
            painter.add(egui::Shape::line(points, egui::Stroke::new(2.0, line_color)));
        }
        // Peak markers sit at the *curve's* height at that position, not the
        // peak's raw intensity, so a hover near a marker lands on what is
        // actually visible.
        peak_points = peaks
            .iter()
            .map(|&(x, _)| to_screen(x, intensity_at(peaks, kind, width, x)))
            .collect();
        for &p in &peak_points {
            painter.circle_filled(p, 2.5, line_color);
        }
    }

    if let Some(hover_pos) = response.hover_pos() {
        if let Some((idx, dist)) = peak_points
            .iter()
            .enumerate()
            .map(|(i, p)| (i, p.distance(hover_pos)))
            .min_by(|a, b| a.1.total_cmp(&b.1))
        {
            if dist <= HOVER_RADIUS {
                // Fixed dark colour, matching the plot's own fixed white
                // background rather than the app's theme -- the theme's
                // "strong text" is often light, which would vanish here.
                painter.circle_stroke(peak_points[idx], 5.0, egui::Stroke::new(2.0, plot_color));
                let galley = painter.layout_no_wrap(
                    hover_label(idx),
                    egui::FontId::proportional(AXIS_LABEL_SIZE),
                    plot_color,
                );
                let padding = egui::vec2(6.0, 4.0);
                let box_size = galley.size() + padding * 2.0;
                let mut box_pos = peak_points[idx] + egui::vec2(12.0, -box_size.y - 12.0);
                box_pos.x = box_pos.x.clamp(rect.left(), (rect.right() - box_size.x).max(rect.left()));
                box_pos.y = box_pos.y.clamp(rect.top(), (rect.bottom() - box_size.y).max(rect.top()));
                let box_rect = egui::Rect::from_min_size(box_pos, box_size);
                painter.rect_filled(box_rect, 4.0, egui::Color32::from_rgb(240, 240, 240));
                painter.rect_stroke(
                    box_rect,
                    4.0,
                    egui::Stroke::new(1.0, plot_color.gamma_multiply(0.5)),
                    egui::StrokeKind::Outside,
                );
                painter.galley(box_rect.min + padding, galley, plot_color);
            }
        }
    }

    painter.text(
        egui::pos2(plot_rect.center().x, rect.bottom() - 4.0),
        egui::Align2::CENTER_BOTTOM,
        axis.title,
        egui::FontId::proportional(AXIS_TITLE_SIZE),
        plot_color,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The point of a "nice" step is round numbers, so check the mantissa is
    /// one of 1/2/5 and that the tick count stays readable.
    #[test]
    fn nice_steps_are_round_numbers_at_a_readable_density() {
        for span in [1.0, 3.7, 12.0, 55.0, 500.0, 4321.0, 0.85, 0.02] {
            let step = nice_tick_step(span);
            let intervals = span / step;
            assert!(
                (3.0..=12.0).contains(&intervals),
                "span {span} gave step {step} -> {intervals} intervals"
            );
            let mantissa = step / 10f64.powf(step.log10().floor());
            assert!(
                (mantissa - 1.0).abs() < 1e-9
                    || (mantissa - 2.0).abs() < 1e-9
                    || (mantissa - 5.0).abs() < 1e-9,
                "step {step} is not 1/2/5 x 10^n (mantissa {mantissa})"
            );
        }
    }

    /// A typical UV-Vis range in each unit should come out with the steps a
    /// chemist would draw by hand.
    #[test]
    fn typical_uv_vis_ranges_get_conventional_steps() {
        assert_eq!(nice_tick_step(700.0 - 200.0), 50.0, "200-700 nm");
        assert_eq!(nice_tick_step(6.0 - 1.5), 0.5, "1.5-6 eV");
    }

    #[test]
    fn a_degenerate_span_does_not_produce_a_zero_or_infinite_step() {
        for span in [0.0, -5.0, f64::NAN, f64::INFINITY] {
            let step = nice_tick_step(span);
            assert!(step.is_finite() && step > 0.0, "span {span} gave step {step}");
        }
    }

    #[test]
    fn ticks_are_multiples_of_the_step_inside_the_range() {
        let ticks = tick_positions(200.0, 700.0, 50.0);
        assert_eq!(ticks, vec![250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0]);
    }

    /// The IR spectrum's axis starts at zero; a gridline drawn there would
    /// land exactly on the plot border and read as a doubled frame.
    #[test]
    fn a_tick_exactly_on_the_left_edge_is_skipped() {
        let ticks = tick_positions(0.0, 2000.0, 500.0);
        assert_eq!(ticks, vec![500.0, 1000.0, 1500.0, 2000.0]);
    }

    #[test]
    fn an_inverted_or_degenerate_range_yields_no_ticks() {
        assert!(tick_positions(700.0, 200.0, 50.0).is_empty());
        assert!(tick_positions(200.0, 200.0, 50.0).is_empty());
        assert!(tick_positions(200.0, 700.0, 0.0).is_empty());
    }

    /// A step far smaller than the range must not spin the tick loop forever.
    #[test]
    fn an_absurdly_fine_step_is_bounded() {
        let ticks = tick_positions(0.0, 1.0e6, 1.0e-3);
        assert!(ticks.len() <= 1001);
    }
}
