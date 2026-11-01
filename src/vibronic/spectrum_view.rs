//! Display transformations leave computed spectra unchanged and exportable.
use super::{thermal_spectrum::HarmonicSpectrum, CM1_TO_AU, HARTREE_TO_EV};
use bevy_egui::egui;
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Axis {
    #[default]
    Ev,
    Wavenumber,
    Wavelength,
}
impl Axis {
    pub fn label(self) -> &'static str {
        match self {
            Self::Ev => "eV",
            Self::Wavenumber => "cm⁻¹",
            Self::Wavelength => "nm",
        }
    }
    pub fn from_ev(self, x: f64) -> f64 {
        match self {
            Self::Ev => x,
            Self::Wavenumber => x / (CM1_TO_AU * HARTREE_TO_EV),
            Self::Wavelength => 1239.841984332 / x,
        }
    }
    pub fn to_ev(self, x: f64) -> f64 {
        match self {
            Self::Ev => x,
            Self::Wavenumber => x * CM1_TO_AU * HARTREE_TO_EV,
            Self::Wavelength => 1239.841984332 / x,
        }
    }
}
#[derive(Default)]
pub struct View {
    pub axis: Axis,
    pub experimental_axis: Axis,
    pub shift_ev: f64,
    pub experimental: Vec<(f64, f64)>,
    pub experimental_name: String,
    pub normalise: bool,
    pub selected: usize,
    pub is_raman: bool,
    pub initialised: bool,
}
impl View {
    pub fn initialise(&mut self, raman: bool) {
        if self.initialised && self.is_raman != raman {
            self.experimental.clear();
        }
        self.is_raman = raman;
        self.experimental_axis = if raman { Axis::Wavenumber } else { Axis::Ev };
        self.axis = if raman { Axis::Wavenumber } else { Axis::Ev };
        self.shift_ev = 0.;
        self.selected = 0;
        self.normalise = true;
        self.initialised = true;
    }
}
pub fn parse_experimental(text: &str, axis: Axis) -> Result<Vec<(f64, f64)>, String> {
    let mut rows = vec![];
    let mut header = false;
    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line
            .split(|c: char| c == ',' || c == ';' || c.is_whitespace())
            .filter(|s| !s.is_empty())
            .collect();
        if fields.len() != 2 {
            return Err(format!(
                "Experimental spectrum line {} needs two columns",
                i + 1
            ));
        }
        let (x, y) = match (fields[0].parse::<f64>(), fields[1].parse::<f64>()) {
            (Ok(x), Ok(y)) => (x, y),
            _ if rows.is_empty() && !header => {
                header = true;
                continue;
            }
            _ => return Err(format!("Invalid experimental spectrum line {}", i + 1)),
        };
        if !x.is_finite() || !y.is_finite() || x < 0. || (axis == Axis::Wavelength && x == 0.) {
            return Err(format!(
                "Invalid experimental spectrum values on line {}",
                i + 1
            ));
        }
        rows.push((axis.to_ev(x), y));
    }
    if rows.len() < 2 {
        return Err("Experimental spectrum needs at least two numeric rows".into());
    }
    rows.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(rows)
}
pub fn displayed(
    points: &[(f64, f64)],
    axis: Axis,
    shift: f64,
    normalise: bool,
) -> Vec<(f64, f64)> {
    let scale = if normalise {
        points
            .iter()
            .map(|p| p.1.abs())
            .fold(0_f64, f64::max)
            .max(1e-300)
    } else {
        1.
    };
    let mut out: Vec<_> = points
        .iter()
        .filter_map(|&(x, y)| {
            let x = x + shift;
            if axis == Axis::Wavelength && x <= 0. {
                None
            } else {
                Some((axis.from_ev(x), y / scale))
            }
        })
        .collect();
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out
}
fn visible(points: &[(f64, f64)]) -> &[(f64, f64)] {
    let peak = points.iter().map(|p| p.1.abs()).fold(0_f64, f64::max);
    let first = points
        .iter()
        .position(|p| p.1.abs() > peak * 1e-5)
        .unwrap_or(0);
    let last = points
        .iter()
        .rposition(|p| p.1.abs() > peak * 1e-5)
        .map(|n| n + 1)
        .unwrap_or(points.len());
    &points[first.saturating_sub(5)..(last + 5).min(points.len())]
}
pub fn csv(curves: &[(&str, Vec<(f64, f64)>)], axis: Axis, shift: f64, normalised: bool) -> String {
    use std::fmt::Write;
    let mut out=format!("# Display axis: {}\n# Calculated energy shift_eV: {shift:.10}\n# Peak-normalised: {normalised}\n# Intensities are sampled values; axis conversion does not apply a density Jacobian.\n",axis.label());
    for (label, points) in curves {
        let _ = writeln!(out, "# {label}\nx,intensity");
        for (x, y) in points {
            let _ = writeln!(out, "{x:.10},{y:.12e}");
        }
    }
    out
}
fn bounds(curves: &[(&str, Vec<(f64, f64)>)]) -> (f64, f64, f64, f64) {
    let mut b = (f64::INFINITY, f64::NEG_INFINITY, 0_f64, 0_f64);
    for (_, points) in curves {
        for (x, y) in points {
            b.0 = b.0.min(*x);
            b.1 = b.1.max(*x);
            b.2 = b.2.min(*y);
            b.3 = b.3.max(*y);
        }
    }
    b.1 = b.1.max(b.0 + 1e-9);
    b.3 = b.3.max(b.2 + 1e-12);
    b
}
fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
pub fn svg(curves: &[(&str, Vec<(f64, f64)>)], axis: Axis, title: &str) -> String {
    use std::fmt::Write;
    let (lo, hi, bottom, top) = bounds(curves);
    let mut s=format!("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"640\" viewBox=\"0 0 1000 640\"><rect width=\"1000\" height=\"640\" fill=\"white\"/><g font-family=\"sans-serif\" font-size=\"16\" fill=\"black\"><text x=\"80\" y=\"30\">{}</text><path d=\"M80 65V565H960\" fill=\"none\" stroke=\"black\"/>",escape(title));
    for k in 0..=5 {
        let x = lo + (hi - lo) * k as f64 / 5.;
        let y = bottom + (top - bottom) * k as f64 / 5.;
        let _=write!(s,"<text x=\"{}\" y=\"590\" text-anchor=\"middle\">{x:.3}</text><text x=\"72\" y=\"{}\" text-anchor=\"end\">{y:.3}</text>",80+176*k,570-100*k);
    }
    for (k, (label, points)) in curves.iter().enumerate() {
        let color = if k == 0 { "#1261a0" } else { "#c33d22" };
        let data = points
            .iter()
            .map(|(x, y)| {
                format!(
                    "{:.3},{:.3}",
                    80. + (x - lo) / (hi - lo) * 880.,
                    565. - (y - bottom) / (top - bottom) * 500.
                )
            })
            .collect::<Vec<_>>()
            .join(" ");
        let _=write!(s,"<polyline fill=\"none\" stroke=\"{color}\" stroke-width=\"1.5\" points=\"{data}\"/><text x=\"{}\" y=\"54\" fill=\"{color}\">{}</text>",80+k*440,escape(label));
    }
    let _ = write!(
        s,
        "<text x=\"520\" y=\"625\" text-anchor=\"middle\">{}</text></g></svg>",
        axis.label()
    );
    s
}
pub fn plot(ui: &mut egui::Ui, curves: &[(&str, Vec<(f64, f64)>)], axis: Axis) {
    let (lo, hi, bottom, top) = bounds(curves);
    if !lo.is_finite() {
        return;
    }
    let (whole, response) = ui.allocate_exact_size(
        egui::vec2(ui.available_width().max(260.), 270.),
        egui::Sense::hover(),
    );
    let rect = egui::Rect::from_min_max(
        whole.min + egui::vec2(55., 25.),
        whole.max - egui::vec2(15., 40.),
    );
    let painter = ui.painter();
    let text = ui.visuals().text_color();
    let font = egui::FontId::proportional(11.);
    let position = |x: f64, y: f64| {
        egui::pos2(
            rect.left() + ((x - lo) / (hi - lo)) as f32 * rect.width(),
            rect.bottom() - ((y - bottom) / (top - bottom)) as f32 * rect.height(),
        )
    };
    painter.line_segment([rect.left_top(), rect.left_bottom()], (1., text));
    painter.line_segment([rect.left_bottom(), rect.right_bottom()], (1., text));
    for k in 0..=4 {
        let x = lo + (hi - lo) * k as f64 / 4.;
        painter.text(
            position(x, bottom) + egui::vec2(0., 5.),
            egui::Align2::CENTER_TOP,
            format!("{x:.3}"),
            font.clone(),
            text,
        );
    }
    painter.text(
        whole.center_bottom(),
        egui::Align2::CENTER_BOTTOM,
        axis.label(),
        font.clone(),
        text,
    );
    let mut hover = String::new();
    for (k, (label, points)) in curves.iter().enumerate() {
        let color = if k == 0 {
            egui::Color32::from_rgb(90, 180, 250)
        } else {
            egui::Color32::from_rgb(240, 125, 80)
        };
        painter.text(
            whole.left_top() + egui::vec2(k as f32 * 260., 0.),
            egui::Align2::LEFT_TOP,
            label,
            font.clone(),
            color,
        );
        // Retain each pixel column's extrema to preserve narrow peaks.
        let cols = rect.width().ceil() as usize;
        let mut bins = vec![None::<(f64, f64, f64)>; cols + 1];
        for &(x, y) in points {
            let at = (((x - lo) / (hi - lo)) * cols as f64).round() as usize;
            let b = &mut bins[at.min(cols)];
            *b = Some(match *b {
                Some((x, a, b)) => (x, a.min(y), b.max(y)),
                None => (x, y, y),
            });
        }
        let vertices = bins
            .into_iter()
            .flatten()
            .flat_map(|(x, a, b)| [position(x, a), position(x, b)])
            .collect();
        painter.add(egui::Shape::line(vertices, egui::Stroke::new(1.5, color)));
        if let Some(pos) = response.hover_pos().filter(|p| rect.contains(*p)) {
            let x = lo + (pos.x - rect.left()) as f64 / rect.width() as f64 * (hi - lo);
            if let Some(&(x, y)) = points
                .iter()
                .min_by(|a, b| (a.0 - x).abs().total_cmp(&(b.0 - x).abs()))
            {
                painter.circle_filled(position(x, y), 3., color);
                hover.push_str(&format!("{label}: {x:.6} {} · {y:.8}\n", axis.label()));
            }
        }
    }
    if !hover.is_empty() {
        response.on_hover_text(hover);
    }
}
pub fn results(ui: &mut egui::Ui, view: &mut View, spectra: &[(String, HarmonicSpectrum)]) {
    use crate::ui_style::ResponseExt;
    if spectra.is_empty() {
        return;
    }
    view.selected = view.selected.min(spectra.len() - 1);
    ui.horizontal_wrapped(|ui| {
        egui::ComboBox::from_id_salt("vibronic_plot_series")
            .selected_text(&spectra[view.selected].0)
            .show_ui(ui, |ui| {
                for (i, (name, _)) in spectra.iter().enumerate() {
                    ui.selectable_value(&mut view.selected, i, name);
                }
            });
        for axis in [Axis::Ev, Axis::Wavenumber, Axis::Wavelength] {
            if view.is_raman && axis == Axis::Wavelength {
                continue;
            }
            ui.selectable_value(&mut view.axis, axis, axis.label());
        }
        ui.checkbox(&mut view.normalise, "Normalise peaks");
    });
    ui.horizontal(|ui| {
        ui.label("Calculated shift");
        if view.is_raman {
            let mut shift = Axis::Wavenumber.from_ev(view.shift_ev);
            if ui
                .add(egui::DragValue::new(&mut shift).speed(1.).suffix(" cm⁻¹"))
                .changed()
            {
                view.shift_ev = Axis::Wavenumber.to_ev(shift);
            }
        } else {
            ui.add(
                egui::DragValue::new(&mut view.shift_ev)
                    .speed(0.001)
                    .suffix(" eV"),
            );
        }
        if ui.button("Reset shift").clicked_once() {
            view.shift_ev = 0.;
        }
    });
    ui.horizontal_wrapped(|ui| {
        ui.label("Experimental input");
        for axis in [Axis::Ev, Axis::Wavenumber, Axis::Wavelength] {
            if view.is_raman && axis == Axis::Wavelength {
                continue;
            }
            ui.selectable_value(&mut view.experimental_axis, axis, axis.label());
        }
        if ui.button("Load spectrum…").clicked_once() {
            if let Some(path) = rfd::FileDialog::new().pick_file() {
                match std::fs::read_to_string(&path)
                    .map_err(|e| e.to_string())
                    .and_then(|text| parse_experimental(&text, view.experimental_axis))
                {
                    Ok(data) => {
                        view.experimental = data;
                        view.experimental_name = path
                            .file_name()
                            .unwrap_or_default()
                            .to_string_lossy()
                            .into_owned();
                    }
                    Err(e) => crate::runtime_warning::push(ui.ctx(), "Experimental spectrum", &e),
                }
            }
        }
        if !view.experimental.is_empty() && ui.button("Remove overlay").clicked_once() {
            view.experimental.clear();
        }
    });
    let (label, spectrum) = &spectra[view.selected];
    let calculated = displayed(&spectrum.profile, view.axis, view.shift_ev, view.normalise);
    let mut curves = vec![(label.as_str(), visible(&calculated).to_vec())];
    if !view.experimental.is_empty() {
        curves.push((
            view.experimental_name.as_str(),
            displayed(&view.experimental, view.axis, 0., view.normalise),
        ));
    }
    plot(ui, &curves, view.axis);
    ui.horizontal(|ui| {
        if ui.button("Save figure…").clicked_once() {
            if let Some(path) = rfd::FileDialog::new()
                .set_file_name("vibronic_spectrum.svg")
                .add_filter("SVG", &["svg"])
                .save_file()
            {
                if let Err(e) = std::fs::write(path, svg(&curves, view.axis, label)) {
                    crate::runtime_warning::push(ui.ctx(), "Save figure", &e.to_string());
                }
            }
        }
        if ui.button("Copy displayed data").clicked_once() {
            ui.ctx()
                .copy_text(csv(&curves, view.axis, view.shift_ev, view.normalise));
        }
        if ui.button("Save displayed data…").clicked_once() {
            if let Some(path) = rfd::FileDialog::new()
                .set_file_name("vibronic_display.csv")
                .save_file()
            {
                if let Err(e) =
                    std::fs::write(path, csv(&curves, view.axis, view.shift_ev, view.normalise))
                {
                    crate::runtime_warning::push(ui.ctx(), "Save displayed data", &e.to_string());
                }
            }
        }
    });
    ui.collapsing("Calculation details",|ui|{ui.label(format!("0–0: {:.6} eV · {:.2} K",spectrum.zero_zero_ev,spectrum.temperature));ui.label(format!("Time samples: {} · Correlation weight: {:.6}",spectrum.correlation.len(),spectrum.integral));for note in &spectrum.notes{ui.label(note);}ui.label("Axis conversion preserves sampled intensities, not spectral densities. Raw computed data remain unchanged.");});
}
