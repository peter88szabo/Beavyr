//! The relaxed Poor Man's NEB path, plotted without rerunning the backend.

use bevy_egui::egui;

use super::parameters::Algorithm;
use super::{GenerationOutput, TsGeneration};

// Matches optimizer::rda's saved NEB energy profile.
const HARTREE_TO_KJ_MOL: f64 = 2625.499638;

struct Profile {
    /// Image number (starting at 1), energy relative to A in kJ/mol.
    points: Vec<[f64; 2]>,
    selected: usize,
    lower: f64,
    upper: f64,
}

impl Profile {
    fn from_energies(energies: &[f64]) -> Option<Self> {
        if energies.len() < 3 || energies.iter().any(|e| !e.is_finite()) {
            return None;
        }
        let points: Vec<_> = energies
            .iter()
            .enumerate()
            .map(|(i, e)| [i as f64 + 1.0, (e - energies[0]) * HARTREE_TO_KJ_MOL])
            .collect();
        // The generator chooses the highest interior image, even when an
        // endpoint is higher. Match its tie-breaking (last equal maximum).
        let selected = (1..energies.len() - 1)
            .max_by(|&a, &b| energies[a].partial_cmp(&energies[b]).unwrap())?;
        let min = points.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min);
        let max = points
            .iter()
            .map(|p| p[1])
            .fold(f64::NEG_INFINITY, f64::max);
        let padding = ((max - min) * 0.08).max(0.1);
        Some(Self {
            points,
            selected,
            lower: min - padding,
            upper: max + padding,
        })
    }

    fn from_output(output: &GenerationOutput) -> Option<Self> {
        if output.algorithm != Algorithm::PoorMansNeb {
            return None;
        }
        Self::from_energies(
            &output
                .guess
                .points
                .iter()
                .map(|p| p.energy)
                .collect::<Vec<_>>(),
        )
    }
}

/// Return the window rectangle so the viewport's picking/camera logic knows
/// this is UI. Called independently of the TS panel, so the plot stays open
/// when the generation panel closes.
pub fn window(ctx: &egui::Context, state: &mut TsGeneration) -> Option<egui::Rect> {
    if !state.energy_plot_open {
        return None;
    }
    let Some(output) = state
        .output
        .as_ref()
        .filter(|o| o.algorithm == Algorithm::PoorMansNeb)
    else {
        state.energy_plot_open = false;
        return None;
    };
    egui::Window::new("Poor Man's NEB Energy")
        .id(egui::Id::new("poormans_neb_energy_window"))
        .open(&mut state.energy_plot_open)
        .resizable(true)
        .default_size([620.0, 400.0])
        .min_size([440.0, 300.0])
        .show(ctx, |ui| {
            ui.label(&output.method);
            let Some(profile) = Profile::from_output(output) else {
                ui.label("No finite NEB energy profile is available.");
                return;
            };
            ui.label(format!(
                "Energy relative to A · selected TS guess: image {} ({:+.2} kJ/mol)",
                profile.selected + 1,
                profile.points[profile.selected][1]
            ));
            ui.weak("A and B are the endpoints. Hover over an image for its energy.");
            draw(ui, output, &profile);
        })
        .map(|window| window.response.rect)
}

fn draw(ui: &mut egui::Ui, output: &GenerationOutput, profile: &Profile) {
    let size = egui::vec2(
        ui.available_width().max(380.0),
        ui.available_height().max(200.0),
    );
    let (rect, response) = ui.allocate_exact_size(size, egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let plot = egui::Rect::from_min_max(
        rect.min + egui::vec2(74.0, 25.0),
        rect.max - egui::vec2(20.0, 42.0),
    );
    let text = ui.visuals().text_color();
    let grid = text.gamma_multiply(0.18);
    let blue = egui::Color32::from_rgb(90, 170, 240);
    let endpoint = egui::Color32::from_rgb(80, 195, 130);
    let selected = egui::Color32::from_rgb(240, 145, 70);
    let count = profile.points.len();
    let to_screen = |point: [f64; 2]| {
        egui::pos2(
            plot.left() + ((point[0] - 1.0) / (count - 1) as f64) as f32 * plot.width(),
            plot.bottom()
                - ((point[1] - profile.lower) / (profile.upper - profile.lower)) as f32
                    * plot.height(),
        )
    };
    painter.text(
        egui::pos2(plot.left(), rect.top()),
        egui::Align2::LEFT_TOP,
        "ΔE from A (kJ/mol)",
        egui::FontId::proportional(14.0),
        text,
    );
    for tick in 0..=5 {
        let energy = profile.lower + (profile.upper - profile.lower) * tick as f64 / 5.0;
        let y = to_screen([1.0, energy]).y;
        painter.line_segment(
            [egui::pos2(plot.left(), y), egui::pos2(plot.right(), y)],
            egui::Stroke::new(1.0, grid),
        );
        painter.text(
            egui::pos2(plot.left() - 6.0, y),
            egui::Align2::RIGHT_CENTER,
            format!("{energy:.2}"),
            egui::FontId::monospace(12.0),
            text,
        );
    }
    let tick_count = (count - 1).min(5);
    for tick in 0..=tick_count {
        let index = tick * (count - 1) / tick_count;
        let x = to_screen(profile.points[index]).x;
        let label = if index == 0 {
            "1 (A)".into()
        } else if index == count - 1 {
            format!("{count} (B)")
        } else {
            (index + 1).to_string()
        };
        painter.line_segment(
            [egui::pos2(x, plot.top()), egui::pos2(x, plot.bottom())],
            egui::Stroke::new(1.0, grid),
        );
        painter.text(
            egui::pos2(x, plot.bottom() + 5.0),
            egui::Align2::CENTER_TOP,
            label,
            egui::FontId::monospace(12.0),
            text,
        );
    }
    painter.rect_stroke(
        plot,
        0.0,
        egui::Stroke::new(1.0, text.gamma_multiply(0.6)),
        egui::StrokeKind::Outside,
    );
    let screen: Vec<_> = profile.points.iter().copied().map(to_screen).collect();
    painter.add(egui::Shape::line(
        screen.clone(),
        egui::Stroke::new(2.0, blue),
    ));
    for (i, p) in screen.iter().enumerate() {
        let color = if i == profile.selected {
            selected
        } else if i == 0 || i == count - 1 {
            endpoint
        } else {
            blue
        };
        painter.circle_filled(*p, if i == profile.selected { 5.0 } else { 3.5 }, color);
    }
    painter.text(
        screen[profile.selected] - egui::vec2(0.0, 9.0),
        egui::Align2::CENTER_BOTTOM,
        "TS",
        egui::FontId::proportional(13.0),
        selected,
    );
    painter.text(
        egui::pos2(plot.center().x, rect.bottom()),
        egui::Align2::CENTER_BOTTOM,
        "Image along path (A → B)",
        egui::FontId::proportional(14.0),
        text,
    );
    if let Some(cursor) = response.hover_pos() {
        if let Some((index, _)) = screen
            .iter()
            .enumerate()
            .filter(|(_, p)| p.distance(cursor) <= 16.0)
            .min_by(|(_, a), (_, b)| a.distance(cursor).total_cmp(&b.distance(cursor)))
        {
            painter.circle_stroke(screen[index], 7.0, egui::Stroke::new(1.5, text));
            response.on_hover_text(format!(
                "Image {}{}\nEnergy: {:.10} Eh\nRelative to A: {:+.4} kJ/mol\n{}",
                index + 1,
                if index == profile.selected {
                    " · selected TS"
                } else if index == 0 {
                    " · A"
                } else if index == count - 1 {
                    " · B"
                } else {
                    ""
                },
                output.guess.points[index].energy,
                profile.points[index][1],
                output.guess.points[index].message,
            ));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relative_energies_use_a_and_keep_both_endpoints_and_negative_values() {
        let profile = Profile::from_energies(&[-10.0, -9.99, -10.02]).unwrap();
        assert_eq!(profile.points[0], [1.0, 0.0]);
        assert_eq!(profile.points[2][0], 3.0);
        assert!((profile.points[1][1] - 26.25499638).abs() < 1.0e-8);
        assert!((profile.points[2][1] + 52.50999276).abs() < 1.0e-8);
    }

    #[test]
    fn selected_ts_is_the_highest_interior_image_even_if_b_is_higher() {
        let profile = Profile::from_energies(&[0.0, 1.0, 2.0, 2.0, 4.0]).unwrap();
        assert_eq!(profile.selected, 3);
        assert!(profile.upper > profile.points[4][1]);
    }

    #[test]
    fn flat_profile_has_finite_nonzero_plot_range() {
        let profile = Profile::from_energies(&[-42.0; 5]).unwrap();
        assert!(profile.lower < 0.0 && profile.upper > 0.0);
        assert!(profile.points.iter().all(|p| p[1] == 0.0));
    }

    #[test]
    fn incomplete_and_nonfinite_profiles_are_not_plotted() {
        assert!(Profile::from_energies(&[]).is_none());
        assert!(Profile::from_energies(&[0.0, 1.0]).is_none());
        assert!(Profile::from_energies(&[0.0, f64::NAN, 1.0]).is_none());
    }
}
