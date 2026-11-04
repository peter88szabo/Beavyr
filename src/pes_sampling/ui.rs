//! Independent PES-sampling setup; completion always loads the sampled trajectory.
use super::{input, sampling, state::*};
use crate::{molecule::Molecule, trajectory::TrajectoryState, ui_style::ResponseExt};
use bevy::tasks::{block_on, futures_lite::future, AsyncComputeTaskPool};
use bevy_egui::egui;
use std::sync::{
    atomic::{AtomicBool, AtomicUsize, Ordering},
    Arc,
};
fn error(ctx: &egui::Context, state: &mut PesSampling, message: String) {
    state.message = Some(message.clone());
    state.is_error = true;
    crate::runtime_warning::push(ctx, "PES sampling", &message);
}
fn number(ui: &mut egui::Ui, label: &str, value: &mut f64, speed: f64) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(value).speed(speed).max_decimals(6));
    });
}
fn count(ui: &mut egui::Ui, label: &str, value: &mut usize) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(value).range(1..=1_000_000));
    });
}
fn excitation(ui: &mut egui::Ui, e: &mut Excitation) {
    ui.horizontal_wrapped(|ui| {
        egui::ComboBox::from_id_salt("population")
            .selected_text(e.population.label())
            .show_ui(ui, |ui| {
                for p in [
                    Population::Zpe,
                    Population::Thermal,
                    Population::Quantum,
                    Population::Energy,
                ] {
                    ui.selectable_value(&mut e.population, p, p.label());
                }
            });
        match e.population {
            Population::Zpe => {}
            Population::Thermal => {
                ui.add(
                    egui::DragValue::new(&mut e.temperature)
                        .speed(1.)
                        .suffix(" K")
                        .range(0.0..=1e8),
                );
            }
            Population::Quantum => {
                ui.add(egui::DragValue::new(&mut e.quantum).prefix("n = "));
            }
            Population::Energy => {
                ui.add(
                    egui::DragValue::new(&mut e.energy_kcal)
                        .speed(0.1)
                        .suffix(" kcal/mol")
                        .range(0.0..=1e8),
                );
            }
        }
    });
    if e.population == Population::Thermal {
        egui::ComboBox::from_id_salt("temperature_rule")
            .selected_text(e.temperature_rule.label())
            .show_ui(ui, |ui| {
                for rule in [
                    TemperatureRule::BoltzmannTail,
                    TemperatureRule::ZpePlusKt,
                    TemperatureRule::SmiteThermal,
                ] {
                    ui.selectable_value(&mut e.temperature_rule, rule, rule.label());
                }
            });
        if e.temperature_rule == TemperatureRule::BoltzmannTail {
            ui.horizontal(|ui| {
                ui.label("Retained population");
                ui.add(
                    egui::DragValue::new(&mut e.retained_population)
                        .suffix(" %")
                        .speed(0.01)
                        .range(0.001..=99.999999),
                );
            });
        }
    }
}
fn partner(ui: &mut egui::Ui, p: &mut Partner, mol: &Molecule) -> Result<(), String> {
    let previous_source = p.source;
    ui.horizontal_wrapped(|ui| {
        ui.selectable_value(&mut p.source, Source::Hessian, "Hessian");
        ui.selectable_value(&mut p.source, Source::Diatomic, "Diatomic");
        ui.selectable_value(&mut p.source, Source::Atom, "Atom");
    });
    if p.source != previous_source {
        p.rotors.clear();
        // Per-mode controls belong to the loaded Hessian, not a new source.
        p.use_overrides = false;
    }
    match p.source {
        Source::Hessian => {
            let mut request = None;
            ui.horizontal_wrapped(|ui| {
                if ui.button("Load Hessian…").clicked_once() {
                    request = Some(0);
                }
                if ui
                    .button("Equilibrium XYZ…")
                    .on_hover_text("Needed only for a Hessian without its own atom coordinates")
                    .clicked_once()
                {
                    request = Some(1);
                }
                if ui.button("Use current geometry").clicked_once() {
                    request = Some(2);
                }
            });
            match request {
                Some(0) => {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("Hessian", &["hess", "hessian", "txt", "dat", "log", "out"])
                        .add_filter("All files", &["*"])
                        .pick_file()
                    {
                        let f = input::load_hessian(&path, p.geometry.as_ref())?;
                        p.overrides = vec![p.excitation.clone(); f.frequencies.len()];
                        p.use_overrides = false;
                        if f.imaginary_cm1.len() > 1 {
                            crate::runtime_warning::push(ui.ctx(),"PES sampling: excluded modes",&format!("Multiple imaginary modes are excluded from PES sampling; this is not a first-order saddle. Frequencies (cm^-1): {:?}",f.imaginary_cm1));
                        }
                        p.rotors.clear();
                        p.fragment = Some(f);
                    }
                }
                Some(1) => {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("XYZ", &["xyz"])
                        .pick_file()
                    {
                        let text = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
                        p.geometry = Some(input::geometry(&text)?);
                        p.fragment = None;
                    }
                }
                Some(2) => {
                    let mut atoms = vec![];
                    let mut x = vec![];
                    for (a, pos) in mol.atoms.iter().zip(&mol.pos) {
                        if !crate::molecule::is_dummy(a) {
                            atoms.push(a.clone());
                            x.extend(pos.to_array().map(f64::from));
                        }
                    }
                    if atoms.is_empty() {
                        return Err("No real atoms are displayed".into());
                    }
                    p.geometry = Some((atoms, x));
                    p.fragment = None;
                }
                _ => {}
            }
            if let Some(f) = &p.fragment {
                ui.label(format!(
                    "{} · {} atoms · {} modes",
                    f.name,
                    f.atoms.len(),
                    f.frequencies.len()
                ));
                if f.imaginary_cm1.len() > 1 {
                    ui.colored_label(
                        egui::Color32::YELLOW,
                        format!("{} imaginary mode(s) excluded", f.imaginary_cm1.len()),
                    );
                }
            } else if let Some((atoms, _)) = &p.geometry {
                ui.label(format!("Equilibrium geometry: {} atoms", atoms.len()));
            }
        }
        Source::Diatomic => {
            ui.horizontal(|ui| {
                ui.label("Atoms");
                for text in &mut p.symbols {
                    ui.add(egui::TextEdit::singleline(text).desired_width(40.));
                }
            });
            number(ui, "Equilibrium distance (Å)", &mut p.equilibrium, 0.001);
            number(ui, "Frequency (cm⁻¹)", &mut p.frequency, 1.);
            ui.horizontal(|ui| {
                ui.selectable_value(&mut p.morse, false, "Harmonic");
                ui.selectable_value(&mut p.morse, true, "Morse");
            });
            if p.morse {
                number(
                    ui,
                    "Dissociation energy De (eV)",
                    &mut p.dissociation_ev,
                    0.01,
                );
            }
        }
        Source::Atom => {
            ui.horizontal(|ui| {
                ui.label("Element");
                ui.add(egui::TextEdit::singleline(&mut p.symbols[0]).desired_width(50.));
            });
        }
    }
    if p.source != Source::Atom {
        ui.checkbox(&mut p.rigid, "Rigid geometry");
        ui.add_enabled_ui(!p.rigid, |ui| {
            excitation(ui, &mut p.excitation);
            if p.source == Source::Hessian {
                if let Some(f) = &p.fragment {
                    egui::CollapsingHeader::new("Mode settings").show(ui, |ui| {
                        if ui
                            .checkbox(&mut p.use_overrides, "Set limits per mode")
                            .changed()
                            && p.use_overrides
                        {
                            p.overrides = vec![p.excitation.clone(); f.frequencies.len()];
                        }
                        if p.overrides.len() != f.frequencies.len() {
                            p.overrides = vec![p.excitation.clone(); f.frequencies.len()];
                        }
                        egui::ScrollArea::vertical()
                            .max_height(200.)
                            .show(ui, |ui| {
                                for (i, (&w, e)) in
                                    f.frequencies.iter().zip(&mut p.overrides).enumerate()
                                {
                                    ui.push_id(i, |ui| {
                                        ui.vertical(|ui| {
                                            ui.label(format!(
                                                "Mode {}  {:.1} cm⁻¹",
                                                i + 1,
                                                w / input::CM_TO_AU
                                            ));
                                            ui.add_enabled_ui(p.use_overrides, |ui| {
                                                excitation(ui, e)
                                            });
                                        });
                                    });
                                }
                            });
                    });
                }
            }
        });
    }
    if p.source == Source::Hessian { super::rotor_ui::controls(ui, p); }
    Ok(())
}
fn start(ctx: &egui::Context, state: &mut PesSampling) {
    match sampling::validate(&state.settings) {
        Ok(_) => {}
        Err(e) => {
            error(ctx, state, e);
            return;
        }
    }
    state.total = sampling::total(&state.settings).unwrap();
    state.cancel = Arc::new(AtomicBool::new(false));
    state.progress = Arc::new(AtomicUsize::new(0));
    let seed = if state.settings.repeatable {
        state.settings.seed
    } else {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos() as u64)
            .unwrap_or(42)
    };
    let settings = state.settings.clone();
    let cancel = Arc::clone(&state.cancel);
    let progress = Arc::clone(&state.progress);
    state.task = Some(
        AsyncComputeTaskPool::get()
            .spawn(async move { sampling::generate(settings, seed, &cancel, &progress) }),
    );
    state.message = None;
    state.is_error = false;
}
pub(crate) fn poll(ctx: &egui::Context, state: &mut PesSampling) {
    if let Some(task) = &mut state.task {
        if let Some(result) = block_on(future::poll_once(task)) {
            state.task = None;
            match result {
                Ok(out) if !state.cancel.load(Ordering::SeqCst) => {
                    state.message = Some(format!(
                        "{} samples · seed {} · {} rejected",
                        out.samples.len(),
                        out.seed,
                        out.rejected
                    ));
                    state.is_error = false;
                    state.output = Some(out);
                    state.load_requested = true;
                }
                Ok(_) => {
                    state.message = Some("Cancelled. Previous results were retained.".into());
                }
                Err(e) => {
                    if e.starts_with("Cancelled") {
                        state.message = Some(e);
                        state.is_error = false;
                    } else {
                        error(ctx, state, e);
                    }
                }
            }
        } else {
            ctx.request_repaint_after(std::time::Duration::from_millis(50));
        }
    }
}
pub(crate) fn load_if_requested(
    state: &mut PesSampling,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
) -> bool {
    if state.load_requested {
        state.load_requested = false;
        if let Some(out) = &state.output {
            return out.load(mol, traj);
        }
    }
    false
}
pub(crate) fn window(
    ctx: &egui::Context,
    state: &mut PesSampling,
    mol: &Molecule,
) -> Option<egui::Rect> {
    if !state.open {
        return None;
    }
    let mut open = state.open;
    let response = egui::Window::new("PES Sampling")
        .id(egui::Id::new("pes_sampling_window"))
        .open(&mut open)
        .default_pos([85., 65.])
        .default_size([455., 660.])
        .min_size([300., 220.])
        .resizable(true)
        .show(ctx, |ui| {
            crate::ui_style::modern(ui);
            egui::ScrollArea::both()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    let running = state.task.is_some();
                    let mut load_error = None;
                    ui.add_enabled_ui(!running, |ui| {
                        ui.horizontal(|ui| {
                            ui.selectable_value(
                                &mut state.settings.workflow,
                                Workflow::Unimolecular,
                                "Unimolecular",
                            );
                            ui.selectable_value(
                                &mut state.settings.workflow,
                                Workflow::Bimolecular,
                                "Bimolecular",
                            );
                        });
                        let n = if state.settings.workflow == Workflow::Bimolecular {
                            2
                        } else {
                            1
                        };
                        for i in 0..n {
                            ui.push_id(i, |ui| {
                                crate::ui_style::group(
                                    ui,
                                    if n == 1 {
                                        "Molecule"
                                    } else if i == 0 {
                                        "Partner A"
                                    } else {
                                        "Partner B"
                                    },
                                    |ui| {
                                        if let Err(e) =
                                            partner(ui, &mut state.settings.partners[i], mol)
                                        {
                                            load_error = Some(e);
                                        }
                                    },
                                );
                            });
                        }
                        if n == 1 {
                            count(ui, "Number of samples", &mut state.settings.points);
                        } else {
                            ui.horizontal(|ui| {
                                ui.label("COM-distance (Å)");
                                ui.add(
                                    egui::DragValue::new(&mut state.settings.distance_start)
                                        .speed(0.1)
                                        .max_decimals(6),
                                );
                                ui.label("–");
                                ui.add(
                                    egui::DragValue::new(&mut state.settings.distance_end)
                                        .speed(0.1)
                                        .max_decimals(6),
                                );
                            });
                            count(ui, "Distance points", &mut state.settings.distances);
                            count(
                                ui,
                                "Samples per distance",
                                &mut state.settings.samples_per_distance,
                            );
                            if let Ok(n) = sampling::total(&state.settings) {
                                ui.label(format!("{n} samples total"));
                            }
                        }
                        ui.horizontal_wrapped(|ui| {
                            ui.label("Mode displacement");
                            ui.selectable_value(
                                &mut state.settings.displacement,
                                Displacement::Uniform,
                                Displacement::Uniform.label(),
                            )
                            .on_hover_text(
                                "Q = A·(2u−1): even coverage between the classical turning \
                                 points, including the bottom of the well.",
                            );
                            ui.selectable_value(
                                &mut state.settings.displacement,
                                Displacement::ClassicalPhase,
                                Displacement::ClassicalPhase.label(),
                            )
                            .on_hover_text(
                                "Q = A·cos(2πu): time-weighted (QCT) distribution; samples \
                                 crowd near the turning points.",
                            );
                        });
                        ui.horizontal(|ui| {
                            ui.checkbox(&mut state.settings.repeatable, "Repeatable");
                            if state.settings.repeatable {
                                ui.label("Seed");
                                ui.add(egui::DragValue::new(&mut state.settings.seed));
                            }
                        });
                        egui::CollapsingHeader::new("Overlap filter").show(ui, |ui| {
                            let filter = &mut state.settings.overlap;
                            ui.checkbox(&mut filter.enabled, "Enabled");
                            ui.add_enabled_ui(filter.enabled, |ui| {
                                number(ui, "H–H (Å)", &mut filter.hydrogen_hydrogen, 0.01);
                                number(ui, "H–other (Å)", &mut filter.hydrogen_other, 0.01);
                                number(ui, "Other–other (Å)", &mut filter.other_other, 0.01);
                            });
                        });
                    });
                    if let Some(e) = load_error {
                        error(ctx, state, e);
                    }
                    if running {
                        ui.add(
                            egui::ProgressBar::new(
                                state.progress.load(Ordering::Relaxed) as f32
                                    / state.total.max(1) as f32,
                            )
                            .show_percentage(),
                        );
                        if ui.button("Cancel sampling").clicked_once() {
                            state.cancel.store(true, Ordering::SeqCst);
                        }
                    } else if ui.button("Generate samples").clicked_once() {
                        start(ctx, state);
                    }
                    if let Some(message) = &state.message {
                        ui.add(egui::Label::new(message).selectable(true).wrap());
                    }
                    let mut save_error = None;
                    if let Some(out) = &state.output {
                        ui.horizontal(|ui| {
                            crate::dataset_inspector::button(ui, true);
                            if ui.button("Save XYZ…").clicked_once() {
                                if let Some(path) = rfd::FileDialog::new()
                                    .add_filter("XYZ trajectory", &["xyz"])
                                    .set_file_name("pes_samples.xyz")
                                    .save_file()
                                {
                                    if let Err(e) = out.save(&path) {
                                        save_error = Some(e);
                                    }
                                }
                            }
                            if ui.button("Load as Trajectory").clicked_once() {
                                state.load_requested = true;
                            }
                        });
                        for w in &out.warnings {
                            ui.label(w);
                        }
                    }
                    if let Some(e) = save_error {
                        error(ctx, state, e);
                    }
                });
        });
    state.open = open;
    response.map(|r| r.response.rect)
}
