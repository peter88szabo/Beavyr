//! The UV-Vis tool's panel and its absorption-spectrum window.

use bevy_egui::egui;

use super::types::{ExcitedState, TddftResult};
use super::{SpectrumUnit, UvVisState};
use crate::spectrum::broadening::{format_spectrum_dat, format_spectrum_sticks, BroadeningKind};
use crate::spectrum::plot::{draw_spectrum, SpectrumAxis};

/// Draws the "UV-Vis (TD-DFT)" section: load a TD-DFT output, list its roots
/// with the orbital excitations behind each, and open the spectrum plot.
pub fn uvvis_panel(ui: &mut egui::Ui, state: &mut UvVisState) {
    ui.horizontal(|ui| {
        if ui.button("Load ORCA output\u{2026}").clicked() {
            load_from_dialog(state);
        }
        if state.result.is_some() && ui.button("Clear").clicked() {
            state.result = None;
            state.load_error = None;
            state.expanded_root = None;
            state.spectrum_window_open = false;
        }
    });
    ui.weak("Reads an existing TD-DFT output. ORCA only for now.");

    if let Some(err) = &state.load_error {
        ui.add_space(4.0);
        ui.colored_label(egui::Color32::from_rgb(200, 80, 80), err);
    }

    // egui closures capture `state` mutably, so the read-only view of the
    // parsed result and the fields the controls write are split apart by
    // field here rather than borrowed through `state` as a whole.
    let UvVisState {
        result,
        expanded_root,
        spectrum_window_open,
        ..
    } = &mut *state;
    let Some(result) = result.as_ref() else {
        return;
    };

    ui.separator();
    if let Some(name) = result.source.file_name() {
        ui.label(
            egui::RichText::new(name.to_string_lossy().to_string())
                .strong()
                .monospace(),
        )
        .on_hover_text(result.source.display().to_string());
    }
    ui.label(format!(
        "{} \u{2022} {} \u{2022} {} root{}",
        result.program,
        result.method,
        result.states.len(),
        if result.states.len() == 1 { "" } else { "s" }
    ));
    if result.unrestricted {
        match result.ground_state_s2 {
            Some(s2) => ui.label(format!(
                "Unrestricted reference, ground-state \u{27e8}S\u{b2}\u{27e9} = {s2:.4}"
            )),
            None => ui.label("Unrestricted reference"),
        };
    }
    for note in &result.notes {
        ui.colored_label(egui::Color32::from_rgb(190, 150, 60), note);
    }

    let has_intensities = result.has_intensities();
    ui.add_space(6.0);
    ui.horizontal(|ui| {
        let label = if *spectrum_window_open {
            "Hide Absorption Spectrum"
        } else {
            "Show Absorption Spectrum"
        };
        ui.add_enabled_ui(has_intensities, |ui| {
            if ui.button(label).clicked() {
                *spectrum_window_open = !*spectrum_window_open;
            }
        });
        if !has_intensities {
            ui.weak("(no oscillator strengths in this file)");
        }
    });

    ui.add_space(6.0);
    ui.label(
        egui::RichText::new("Click a root to see the orbital excitations behind it.").weak(),
    );
    ui.add_space(2.0);
    state_table(ui, result, expanded_root);

    // The spectrum lives in its own window, opened from this panel -- the
    // same arrangement the IR spectrum uses.
    let ctx = ui.ctx().clone();
    uvvis_spectrum_window(&ctx, state);
}

fn load_from_dialog(state: &mut UvVisState) {
    let Some(path) = rfd::FileDialog::new()
        .add_filter("ORCA output", &["out", "log", "txt"])
        .add_filter("All files", &["*"])
        .pick_file()
    else {
        return;
    };
    match std::fs::read_to_string(&path) {
        Ok(text) => match super::orca::parse_orca_tddft(&text, &path) {
            Ok(result) => {
                state.result = Some(result);
                state.load_error = None;
                state.expanded_root = None;
            }
            Err(err) => {
                state.load_error = Some(format!("{}: {err}", path.display()));
            }
        },
        Err(err) => {
            state.load_error = Some(format!("cannot read {}: {err}", path.display()));
        }
    }
}

/// The root list.
///
/// Laid out as fixed-width monospace rows rather than an `egui::Grid`: a grid
/// sizes each column to its own contents and then *wraps* anything wider,
/// which split every "orbital pair + percentage" line in two. One preformatted
/// string per row cannot wrap, so a contribution and its percentage always
/// stay on the same line, and the columns line up by construction.
///
/// `<S²>` and `Mult` appear only for an unrestricted reference -- an empty
/// column for a closed-shell run is just noise.
fn state_table(ui: &mut egui::Ui, result: &TddftResult, expanded: &mut Option<usize>) {
    let show_spin = result.unrestricted;
    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.add(
                egui::Label::new(egui::RichText::new(header_text(show_spin)).monospace().strong())
                    .wrap_mode(egui::TextWrapMode::Extend),
            );
            ui.separator();
            for st in &result.states {
                state_row(ui, st, show_spin, expanded);
            }
        });
}

/// Column headings, spaced to match `state_row_text`.
fn header_text(show_spin: bool) -> String {
    let mut text = format!("  {:>3}  {:>8}  {:>8}  {:>10}", "#", "E (eV)", "nm", "f");
    if show_spin {
        text.push_str(&format!("  {:>8}  {:>4}", "<S^2>", "Mult"));
    }
    text
}

/// One row of the list, preformatted so every column lands at the same
/// character offset on every row.
fn state_row_text(st: &ExcitedState, show_spin: bool, is_open: bool) -> String {
    let marker = if is_open { "\u{25be}" } else { "\u{25b8}" };
    let fosc = match st.oscillator_strength {
        Some(f) => format!("{f:>10.6}"),
        None => format!("{:>10}", "--"),
    };
    let mut text = format!(
        "{marker} {:>3}  {:>8.3}  {:>8.1}  {fosc}",
        st.root, st.energy_ev, st.wavelength_nm
    );
    if show_spin {
        match st.s2 {
            Some(s2) => text.push_str(&format!("  {s2:>8.4}")),
            None => text.push_str(&format!("  {:>8}", "--")),
        }
        match st.multiplicity {
            Some(m) => text.push_str(&format!("  {m:>4}")),
            None => text.push_str(&format!("  {:>4}", "--")),
        }
    }
    text
}

fn state_row(
    ui: &mut egui::Ui,
    st: &ExcitedState,
    show_spin: bool,
    expanded: &mut Option<usize>,
) {
    let is_open = *expanded == Some(st.root);
    let row = ui.add(
        egui::Button::new(egui::RichText::new(state_row_text(st, show_spin, is_open)).monospace())
            .frame(false)
            .selected(is_open)
            .wrap_mode(egui::TextWrapMode::Extend),
    );
    if row.clicked() {
        *expanded = if is_open { None } else { Some(st.root) };
    }
    // The transition dipole rides along as a tooltip rather than a line in the
    // expanded block: it is worth keeping, but it is not what the expansion is
    // for and it crowded out the orbital list.
    if let Some(d) = st.transition_dipole_au {
        row.on_hover_text(format!(
            "\u{3bc} = ({:.5}, {:.5}, {:.5}) au{}\n{:.6} au   {:.1} cm\u{207b}\u{b9}",
            d[0],
            d[1],
            d[2],
            st.d2_au2
                .map(|d2| format!("\nD\u{b2} = {d2:.5} au\u{b2}"))
                .unwrap_or_default(),
            st.energy_au,
            st.energy_cm1
        ));
    }

    if !is_open {
        return;
    }
    if st.excitations.is_empty() {
        ui.weak("      No excitations printed above the program's threshold.");
        return;
    }
    for line in excitation_lines(st) {
        ui.add(
            egui::Label::new(egui::RichText::new(line).monospace())
                .wrap_mode(egui::TextWrapMode::Extend),
        );
    }
}

/// One line per contribution: the orbital pair and its percentage, nothing
/// else, indented to sit under the root it belongs to. Padded to a common
/// width so the arrows and the percentages each form a column.
fn excitation_lines(st: &ExcitedState) -> Vec<String> {
    let ranked = st.excitations_by_weight();
    let widest = ranked.iter().map(|e| e.label().len()).max().unwrap_or(0);
    ranked
        .iter()
        .map(|e| format!("      {:<widest$}   {:>5.1} %", e.label(), e.percent()))
        .collect()
}

/// The absorption-spectrum window: the same controls as the IR spectrum, plus
/// a choice of wavelength or energy axis.
pub fn uvvis_spectrum_window(ctx: &egui::Context, state: &mut UvVisState) {
    if !state.spectrum_window_open {
        return;
    }
    if state.result.is_none() {
        state.spectrum_window_open = false;
        return;
    }
    let UvVisState {
        result,
        spectrum_window_open,
        unit,
        broadening,
        width_nm,
        width_ev,
        ..
    } = &mut *state;
    let result = result.as_ref().expect("checked just above");
    let mut open = *spectrum_window_open;
    egui::Window::new("UV-Vis Absorption Spectrum")
        .open(&mut open)
        .resizable(true)
        .default_size([560.0, 360.0])
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Line shape");
                for kind in BroadeningKind::ALL {
                    ui.selectable_value(broadening, kind, kind.label());
                }
            });
            ui.horizontal(|ui| {
                ui.label("Axis");
                ui.selectable_value(unit, SpectrumUnit::Nanometres, "nm");
                ui.selectable_value(unit, SpectrumUnit::ElectronVolts, "eV");
                ui.separator();
                // A stick spectrum has no width -- it is drawn (and exported)
                // at the exact predicted positions.
                let sticks = broadening.is_sticks();
                let chosen = *unit;
                // Each axis keeps its own width: 20 nm and 0.3 eV are both
                // sensible, and carrying one over as the other would either
                // flatten the spectrum or reduce it to spikes.
                let (range, speed) = match chosen {
                    SpectrumUnit::Nanometres => (0.5..=200.0, 0.5),
                    SpectrumUnit::ElectronVolts => (0.005..=2.0, 0.01),
                };
                let width = chosen.width_of_mut(width_nm, width_ev);
                ui.add_enabled_ui(!sticks, |ui| {
                    ui.label(format!("Width ({})", chosen.label()));
                    ui.add(egui::DragValue::new(width).range(range).speed(speed));
                });
            });

            let chosen = *unit;
            let current_width = chosen.width_of(*width_nm, *width_ev);
            let peaks = match chosen {
                SpectrumUnit::Nanometres => result.peaks_nm(),
                SpectrumUnit::ElectronVolts => result.peaks_ev(),
            };
            if peaks.is_empty() {
                ui.weak("This file has no oscillator strengths to plot.");
                return;
            }
            let states = result.states_with_intensity();
            let (x_min, x_max) = axis_range(&peaks, chosen);

            if ui.button("Export .dat\u{2026}").clicked() {
                export_dialog(*broadening, chosen, current_width, &peaks, x_min, x_max);
            }

            let axis = SpectrumAxis {
                x_min,
                x_max,
                tick_step: None,
                tick_decimals: chosen.tick_decimals(),
                title: chosen.axis_title(),
            };
            draw_spectrum(ui, &peaks, *broadening, current_width, &axis, &|i| {
                hover_label(states[i])
            });
        });
    *spectrum_window_open = open;
}

/// What the plot shows when the pointer is on a peak. Both units are given
/// whichever axis is displayed, since a chemist reading one usually wants the
/// other, plus the dominant orbital pair -- the state's character in one line.
fn hover_label(st: &ExcitedState) -> String {
    let mut text = format!(
        "S{}  {:.1} nm  {:.3} eV\nf = {:.6}",
        st.root,
        st.wavelength_nm,
        st.energy_ev,
        st.oscillator_strength.unwrap_or(0.0)
    );
    if let Some(e) = st.dominant_excitation() {
        text.push_str(&format!("\n{}  ({:.0}%)", e.label(), e.weight * 100.0));
    }
    if let Some(s2) = st.s2 {
        text.push_str(&format!("\n\u{27e8}S\u{b2}\u{27e9} = {s2:.3}"));
    }
    text
}

/// The plotted range: the peaks plus a margin, never crossing zero. Unlike an
/// IR spectrum there is no reason to force the axis down to zero -- a UV-Vis
/// plot starting at 0 nm would waste almost all of its width.
pub fn axis_range(peaks: &[(f64, f64)], unit: SpectrumUnit) -> (f64, f64) {
    let margin = unit.axis_margin();
    let lo = peaks.iter().map(|&(x, _)| x).fold(f64::INFINITY, f64::min);
    let hi = peaks.iter().map(|&(x, _)| x).fold(f64::NEG_INFINITY, f64::max);
    if !lo.is_finite() || !hi.is_finite() {
        return (0.0, margin);
    }
    ((lo - margin).max(0.0), hi + margin)
}

fn export_dialog(
    broadening: BroadeningKind,
    unit: SpectrumUnit,
    width: f64,
    peaks: &[(f64, f64)],
    x_min: f64,
    x_max: f64,
) {
    // Sticks export the raw predicted pairs as-is -- there is no continuous
    // curve to resample, so the fixed grid the broadened shapes use does not
    // apply here.
    let text = if broadening.is_sticks() {
        format_spectrum_sticks(peaks)
    } else {
        format_spectrum_dat(
            peaks,
            broadening,
            width,
            x_min,
            x_max,
            unit.export_resolution(),
        )
    };
    let name = format!("uv_vis_spectrum_{}.dat", unit.label());
    if let Some(path) = rfd::FileDialog::new()
        .set_file_name(name)
        .add_filter("Data file", &["dat", "txt"])
        .save_file()
    {
        if let Err(err) = std::fs::write(&path, text) {
            eprintln!("Failed to write {}: {err}", path.display());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::uvvis::types::{ev_to_nm, TddftResult};
    use std::path::Path;

    fn parsed() -> TddftResult {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("AbsorptionSpectrum")
            .join("Exc_TPSSh-D4-def2-TZVP_SCMwater_16.out");
        let text = std::fs::read_to_string(&path).unwrap();
        crate::uvvis::orca::parse_orca_tddft(&text, &path).unwrap()
    }

    #[test]
    fn the_nanometre_axis_brackets_every_peak() {
        let result = parsed();
        let peaks = result.peaks_nm();
        let (lo, hi) = axis_range(&peaks, SpectrumUnit::Nanometres);
        // The example spans 242.4 to 640.8 nm.
        assert!((lo - 192.4).abs() < 1e-6, "got {lo}");
        assert!((hi - 690.8).abs() < 1e-6, "got {hi}");
        assert!(peaks.iter().all(|&(x, _)| x > lo && x < hi));
    }

    #[test]
    fn the_energy_axis_brackets_every_peak() {
        let result = parsed();
        let peaks = result.peaks_ev();
        let (lo, hi) = axis_range(&peaks, SpectrumUnit::ElectronVolts);
        assert!((lo - (1.934849 - 0.5)).abs() < 1e-3, "got {lo}");
        assert!((hi - (5.115694 + 0.5)).abs() < 1e-3, "got {hi}");
    }

    /// A peak closer to zero than the margin must not push the axis negative,
    /// which would put a nonsensical wavelength on the plot.
    #[test]
    fn the_axis_never_crosses_zero() {
        let (lo, _) = axis_range(&[(10.0, 1.0)], SpectrumUnit::Nanometres);
        assert_eq!(lo, 0.0);
    }

    #[test]
    fn an_empty_peak_list_still_gives_a_usable_range() {
        let (lo, hi) = axis_range(&[], SpectrumUnit::Nanometres);
        assert!(hi > lo, "a degenerate axis would divide by zero when scaling");
    }

    /// The hover text is what identifies a band, so it must carry the root,
    /// both units, the intensity and the state's character.
    #[test]
    fn the_hover_readout_names_the_root_and_its_dominant_excitation() {
        let result = parsed();
        let label = hover_label(&result.states[0]);
        assert!(label.starts_with("S1 "), "{label}");
        assert!(label.contains("640.8 nm"), "{label}");
        assert!(label.contains("1.935 eV"), "{label}");
        assert!(label.contains("0.002236"), "{label}");
        assert!(label.contains("b89 --> b90"), "{label}");
        assert!(label.contains("99%"), "{label}");
        assert!(label.contains("\u{27e8}S\u{b2}\u{27e9} = 8.871"), "{label}");
    }

    /// Every column must land at the same character offset on every row --
    /// that is the whole reason the rows are preformatted strings.
    #[test]
    fn every_row_is_the_same_width() {
        let result = parsed();
        let widths: Vec<usize> = result
            .states
            .iter()
            .map(|st| state_row_text(st, true, false).chars().count())
            .collect();
        let first = widths[0];
        assert!(
            widths.iter().all(|&w| w == first),
            "ragged rows would misalign every column: {widths:?}"
        );
        // ...and the heading lines up with them.
        assert_eq!(header_text(true).chars().count(), first);
    }

    #[test]
    fn a_row_carries_its_root_energy_wavelength_and_intensity() {
        let result = parsed();
        let row = state_row_text(&result.states[0], true, false);
        assert!(row.contains("1.935"), "{row}");
        assert!(row.contains("640.8"), "{row}");
        assert!(row.contains("0.002236"), "{row}");
        assert!(row.contains("8.8709"), "{row}");
        assert!(row.trim_end().ends_with('6'), "multiplicity 6: {row}");
    }

    /// A closed-shell run drops the two spin columns from both the heading and
    /// the rows, so they stay the same width as each other.
    #[test]
    fn a_restricted_run_drops_the_spin_columns_from_row_and_heading() {
        let result = parsed();
        let row = state_row_text(&result.states[0], false, false);
        assert!(!row.contains("8.8709"), "{row}");
        assert_eq!(header_text(false).chars().count(), row.chars().count());
        assert!(!header_text(false).contains("S^2"));
    }

    #[test]
    fn the_expansion_marker_does_not_change_the_row_width() {
        let result = parsed();
        let closed = state_row_text(&result.states[0], true, false);
        let open = state_row_text(&result.states[0], true, true);
        assert_eq!(closed.chars().count(), open.chars().count());
        assert_ne!(closed, open, "an expanded row must look expanded");
    }

    /// The bug this layout replaced: a contribution and its percentage were
    /// wrapped onto separate lines by the grid cell they sat in. One string
    /// per contribution cannot be split.
    #[test]
    fn each_contribution_and_its_percentage_share_one_line() {
        let result = parsed();
        for st in &result.states {
            for line in excitation_lines(st) {
                assert!(!line.contains('\n'), "{line:?}");
                assert!(line.contains("-->"), "{line:?}");
                assert!(line.trim_end().ends_with('%'), "{line:?}");
            }
        }
    }

    /// The expanded block is meant to be exactly the orbital pairs and their
    /// percentages -- no coefficients, energies or dipoles crowding them out.
    #[test]
    fn the_expanded_block_lists_only_orbital_pairs_and_percentages() {
        let result = parsed();
        let lines = excitation_lines(&result.states[0]);
        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0].trim(), "b89 --> b90    98.8 %");
    }

    /// Ranked by weight, and padded so the percentages form a column even when
    /// the orbital indices differ in width.
    #[test]
    fn the_expanded_block_is_ranked_and_column_aligned() {
        let result = parsed();
        // Root 10 mixes five contributions with three- and two-digit indices.
        let lines = excitation_lines(&result.states[9]);
        assert_eq!(lines.len(), 5);
        assert!(lines[0].trim_start().starts_with("b87 --> b92"), "{:?}", lines[0]);
        assert!(lines[0].contains("53.8 %"), "{:?}", lines[0]);
        let percent_column = |l: &String| l.rfind('%').unwrap();
        let first = percent_column(&lines[0]);
        assert!(
            lines.iter().all(|l| percent_column(l) == first),
            "percentages must line up: {lines:?}"
        );
        // Descending, so the dominant configuration reads first.
        let percents: Vec<f64> = result.states[9]
            .excitations_by_weight()
            .iter()
            .map(|e| e.percent())
            .collect();
        assert!(percents.windows(2).all(|w| w[0] >= w[1]), "{percents:?}");
    }

    #[test]
    fn a_restricted_state_shows_bare_orbital_indices() {
        let text = "TD-DFT/TDA EXCITED STATES\n\nSTATE  1:  E= 0.1 au 2.721 eV\n    45 ->  47  :  0.9\n\nend\n";
        let r = crate::uvvis::orca::parse_orca_tddft(text, Path::new("t.out")).unwrap();
        assert_eq!(excitation_lines(&r.states[0])[0].trim(), "45 --> 47    90.0 %");
    }

    /// Peak index must address the same state the plot drew, or every hover
    /// would name the wrong root.
    #[test]
    fn peak_order_matches_the_state_order_behind_it() {
        let result = parsed();
        let peaks = result.peaks_nm();
        let states = result.states_with_intensity();
        assert_eq!(peaks.len(), states.len());
        for (i, &(x, y)) in peaks.iter().enumerate() {
            assert_eq!(x, states[i].wavelength_nm);
            assert_eq!(Some(y), states[i].oscillator_strength);
        }
    }

    /// Switching axes reorders nothing: state 1 is the lowest energy but the
    /// *longest* wavelength, so the nm peak list runs the opposite way.
    #[test]
    fn the_energy_and_wavelength_peak_lists_stay_state_aligned() {
        let result = parsed();
        let nm = result.peaks_nm();
        let ev = result.peaks_ev();
        let states = result.states_with_intensity();
        for i in 0..states.len() {
            // The nm positions come from the absorption table, printed to a
            // tenth of a nanometre; the eV ones are derived from hartree. They
            // must agree to within that printed precision.
            assert!((nm[i].0 - ev_to_nm(ev[i].0)).abs() < 0.06, "state {i}");
        }
        assert!(nm[0].0 > nm[79].0, "wavelength descends with root number");
        assert!(ev[0].0 < ev[79].0, "energy ascends with root number");
    }
}
