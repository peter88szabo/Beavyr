//! UV-Vis opens this independent, horizontally and vertically resizable window.
use super::{
    analysis::{self, AnalysisData},
    input::HessianSource,
    plot,
    state::*,
    HARTREE_TO_EV,
};
use crate::ui_style::ResponseExt;
use bevy::tasks::{block_on, futures_lite::future, AsyncComputeTaskPool};
use bevy_egui::egui;
const OPEN_ID: &str = "open_vibronic_analysis";
pub fn open(ctx: &egui::Context) {
    ctx.data_mut(|d| d.insert_temp(egui::Id::new(OPEN_ID), true));
}
fn file(ui: &mut egui::Ui, label: &str, path: &mut String) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::TextEdit::singleline(path).desired_width(260.));
        if ui.button("Browse…").clicked_once() {
            if let Some(p) = rfd::FileDialog::new().set_title(label).pick_file() {
                *path = p.display().to_string();
            }
        }
    });
}
fn source(ui: &mut egui::Ui, label: &str, s: &mut HessianSource) {
    ui.push_id(label, |ui| {
        file(ui, label, &mut s.path);
        file(ui, "XYZ (only if Hessian has no geometry)", &mut s.xyz);
    });
}
pub fn poll(ctx: &egui::Context, state: &mut VibronicState) {
    if let Some(task) = &mut state.task {
        if let Some(result) = block_on(future::poll_once(task)) {
            state.task = None;
            ctx.request_repaint();
            state.open = true;
            match result {
                Ok(_) if state.cancel.load(std::sync::atomic::Ordering::Relaxed) => {}
                Ok(result) => {
                    if !result.notices.is_empty() {
                        crate::runtime_warning::push(
                            ctx,
                            "Vibronic analysis — warnings",
                            &result.notices.join("\n\n"),
                        );
                    }
                    state.heatmap = if let AnalysisData::Duschinsky(d) = &result.data {
                        Some(plot::texture(ctx, d))
                    } else {
                        None
                    };
                    state.plot = 0;
                    state.view.initialise(
                        result.settings.spectrum.process == super::optical_spectrum::Process::Raman,
                    );
                    state.result = Some(result);
                }
                Err(e) => {
                    if !e.starts_with("Cancelled.") {
                        crate::runtime_warning::push(ctx, "Vibronic analysis", &e);
                    }
                    state.error = Some(e);
                }
            }
        } else {
            ctx.request_repaint_after(std::time::Duration::from_millis(100));
        }
    }
}
pub fn window(
    ctx: &egui::Context,
    state: &mut VibronicState,
    uvvis: &crate::uvvis::UvVisState,
) -> Option<egui::Rect> {
    poll(ctx, state);
    if ctx
        .data_mut(|d| d.remove_temp::<bool>(egui::Id::new(OPEN_ID)))
        .unwrap_or(false)
    {
        state.open = true;
    }
    if !state.open {
        return None;
    }
    let mut open = true;
    let response = egui::Window::new("Vibronic analysis")
        .id(egui::Id::new("vibronic_analysis_window"))
        .open(&mut open)
        .resizable(true)
        .default_size([720., 680.])
        .min_size([360., 260.])
        .show(ctx, |ui| {
            egui::ScrollArea::both()
                .auto_shrink([false, false])
                .show(ui, |ui| body(ui, state, uvvis));
        });
    state.open = open;
    response.map(|r| r.response.rect)
}

fn body(ui: &mut egui::Ui, state: &mut VibronicState, uvvis: &crate::uvvis::UvVisState) {
    let ctx = ui.ctx().clone();

    controls(ui, state, uvvis);
    if state.task.is_some() {
        ui.horizontal(|ui| {
            ui.spinner();
            ui.label(format!(
                "Analysing… {}%",
                state.progress.load(std::sync::atomic::Ordering::Relaxed)
            ));
            if ui.button("Cancel").clicked_once() {
                state
                    .cancel
                    .store(true, std::sync::atomic::Ordering::Relaxed);
            }
        });
    }
    if let Some(error) = &state.error {
        if ui.button("Show error").clicked_once() {
            crate::runtime_warning::push(&ctx, "Vibronic analysis", error);
        }
    }
    results(ui, state);
}

fn controls(ui: &mut egui::Ui, state: &mut VibronicState, uvvis: &crate::uvvis::UvVisState) {
    let ctx = ui.ctx().clone();
    ui.add_enabled_ui(state.task.is_none(),|ui| {
                    let s=&mut state.settings;
                    ui.horizontal_wrapped(|ui| {
                        ui.selectable_value(&mut s.kind,AnalysisKind::FranckCondon,"Franck–Condon / resonance Raman");
                        ui.selectable_value(&mut s.kind,AnalysisKind::Duschinsky,"Duschinsky rotation");
                        ui.selectable_value(&mut s.kind,AnalysisKind::Spectrum,"Vibronic spectra");
                    });
                    source(ui,"Ground-state Hessian",&mut s.ground);
                    if matches!(s.kind,AnalysisKind::Duschinsky|AnalysisKind::Spectrum) {source(ui,"Excited-state Hessian",&mut s.excited);}
                    else {
                        file(ui,"Excited-state gradient",&mut s.gradient);
                        ui.small("Gradient at the ground-state geometry: ORCA .engrad or Turbomole $grad, in Hartree/bohr.");
                    }
                    if s.kind==AnalysisKind::Spectrum { super::spectral_controls::controls(ui,&mut s.spectrum,&s.ground); }
                    if s.kind!=AnalysisKind::Spectrum {ui.horizontal(|ui| {ui.label("Minimum frequency");ui.add(egui::DragValue::new(&mut s.cutoff).range(0.0..=10000.).suffix(" cm⁻¹"));});}
                    if s.kind==AnalysisKind::FranckCondon {
                        ui.checkbox(&mut s.make_band,"Absorption band");
                        if s.make_band {
                            ui.horizontal_wrapped(|ui| {
                                ui.label("Vertical energy");ui.add(egui::DragValue::new(&mut s.energy_ev).speed(0.01).range(0.0001..=1000.).suffix(" eV"));
                                if let Some(results)=&uvvis.result {
                                    egui::ComboBox::from_id_salt("vibronic_uvvis_root").selected_text("Use UV-Vis state…").show_ui(ui,|ui| {
                                        for root in &results.states {if ui.selectable_label(false,format!("State {} · {:.4} eV",root.root,root.energy_ev)).clicked_once() {s.energy_ev=root.energy_ev;}}
                                    });
                                }
                            });
                            ui.horizontal_wrapped(|ui| {
                                ui.label("Grid spacing");ui.add(egui::DragValue::new(&mut s.resolution).range(0.1..=1000.).suffix(" cm⁻¹"));
                                ui.add_space(18.);ui.label("Max. quanta per mode");ui.add(egui::DragValue::new(&mut s.quanta).range(1..=1000));
                            });
                        }
                    }
                    if ui.button("Analyse").clicked_once() {
                        let settings=state.settings.clone();state.error=None;state.animation.stop=true;
                        state.cancel=Default::default();state.progress=Default::default();
                        let cancel=state.cancel.clone();let progress=state.progress.clone();
                        state.task=Some(AsyncComputeTaskPool::get().spawn(async move {analysis::run_with_control(settings,&cancel,&progress)}));
                        ctx.request_repaint();
                    }
                });
}

fn results(ui: &mut egui::Ui, state: &mut VibronicState) {
    let ctx = ui.ctx().clone();
    if let Some(result) = &state.result {
        ui.separator();
        ui.horizontal(|ui| {
            if ui.button("Copy data").clicked_once() {
                ctx.copy_text(result.csv());
            }
            if ui.button("Save data…").clicked_once() {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("CSV", &["csv"])
                    .set_file_name("vibronic_analysis.csv")
                    .save_file()
                {
                    if let Err(e) = std::fs::write(&path, result.csv()) {
                        crate::runtime_warning::push(&ctx, "Save vibronic data", &e.to_string());
                    }
                }
            }
        });
        ui.small(format!("Calculated from: {}", result.settings.ground.path));
        if state.error.is_some() || state.task.is_some() {
            ui.label("Showing the previous completed analysis.");
        }
        for note in &result.notices {
            ui.colored_label(ui.visuals().warn_fg_color, note);
        }
        match &result.data {
            AnalysisData::Spectrum(output) => {
                super::spectrum_view::results(ui, &mut state.view, &output.spectra);
                super::mode_animation::controls(
                    ui,
                    &mut state.animation,
                    &output.ground,
                    &output.excited,
                );
                super::spectral_controls::results(ui, output);
            }
            AnalysisData::FranckCondon { activity, band } => {
                ui.label(format!(
                    "Reorganisation: {:.2} cm⁻¹ · unaccounted gradient: {:.2}%",
                    activity.reorganisation_cm1,
                    100. * activity.unaccounted_fraction
                ));
                ui.small("Displaced harmonic approximation; Raman intensities are relative (short-time limit).");
                ui.horizontal(|ui| {
                    ui.selectable_value(&mut state.plot, 0, "Resonance Raman");
                    ui.selectable_value(&mut state.plot, 1, "Huang–Rhys factors");
                    if band.is_some() {
                        ui.selectable_value(&mut state.plot, 2, "Absorption band");
                    }
                });
                match state.plot {
                    2 if band.is_some() => {
                        let b = band.as_ref().unwrap();
                        ui.small(format!(
                            "0 K · no Duschinsky mixing · 0–0: {:.4} eV",
                            b.zero_zero_energy * HARTREE_TO_EV
                        ));
                        plot::curve(ui, &b.profile_ev(), "Energy (eV)", "Bin probability", false);
                    }
                    1 => plot::curve(
                        ui,
                        &activity
                            .modes
                            .iter()
                            .map(|m| (m.frequency_cm1, m.huang_rhys))
                            .collect::<Vec<_>>(),
                        "Frequency (cm⁻¹)",
                        "Huang–Rhys factor",
                        true,
                    ),
                    _ => plot::curve(
                        ui,
                        &activity.normalised_resonance_raman(),
                        "Frequency (cm⁻¹)",
                        "Relative intensity",
                        true,
                    ),
                }
                egui::Grid::new("fc_modes").striped(true).show(ui, |ui| {
                    for h in [
                        "Mode",
                        "cm⁻¹",
                        "Displacement",
                        "S",
                        "λ (cm⁻¹)",
                        "Raman (rel.)",
                    ] {
                        ui.strong(h);
                    }
                    ui.end_row();
                    for (m, (_, rr)) in activity
                        .modes
                        .iter()
                        .zip(activity.normalised_resonance_raman())
                    {
                        for v in [
                            format!("{}", m.index + 1),
                            format!("{:.2}", m.frequency_cm1),
                            format!("{:.5}", m.displacement),
                            format!("{:.5}", m.huang_rhys),
                            format!("{:.3}", m.reorganisation_cm1),
                            format!("{rr:.6}"),
                        ] {
                            ui.label(v);
                        }
                        ui.end_row();
                    }
                });
            }
            AnalysisData::Duschinsky(d) => {
                ui.label(format!(
                    "Mean participation: {:.3} · Σ K²/2: {:.5}",
                    d.mean_participation(),
                    d.total_huang_rhys()
                ));
                if let Some(texture) = &state.heatmap {
                    plot::matrix(ui, d, texture);
                }
                egui::Grid::new("duschinsky_correspondence")
                    .striped(true)
                    .show(ui, |ui| {
                        for h in [
                            "Excited",
                            "cm⁻¹",
                            "Ground",
                            "cm⁻¹",
                            "Overlap",
                            "Participation",
                        ] {
                            ui.strong(h);
                        }
                        ui.end_row();
                        for c in d.correspondence() {
                            for v in [
                                format!("{}", c.column_index + 1),
                                format!("{:.2}", c.column_frequency_cm1),
                                format!("{}", c.row_index + 1),
                                format!("{:.2}", c.row_frequency_cm1),
                                format!("{:.5}", c.overlap),
                                format!("{:.3}", c.participation),
                            ] {
                                ui.label(v);
                            }
                            ui.end_row();
                        }
                    });
                ui.separator();
                egui::Grid::new("duschinsky_displacements")
                    .striped(true)
                    .show(ui, |ui| {
                        for h in ["Ground mode", "cm⁻¹", "K", "K²/2"] {
                            ui.strong(h);
                        }
                        ui.end_row();
                        for (i, mode) in d.rows.iter().enumerate() {
                            let k = d.displacement[i];
                            for v in [
                                format!("{}", mode + 1),
                                format!("{:.2}", d.row_frequencies_cm1[i]),
                                format!("{k:.6}"),
                                format!("{:.6}", 0.5 * k * k),
                            ] {
                                ui.label(v);
                            }
                            ui.end_row();
                        }
                    });
            }
        }
    }
}
