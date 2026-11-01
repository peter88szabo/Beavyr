use super::spectral_job::{Controls, Engine, Output};
use crate::{
    qchem_interfaces::{job::QcMethod, program::QcProgram},
    ui_style::ResponseExt,
};
use bevy_egui::egui;
pub fn controls(ui: &mut egui::Ui, c: &mut Controls, ground: &super::input::HessianSource) {
    ui.horizontal(|ui| {
        ui.selectable_value(
            &mut c.process,
            super::optical_spectrum::Process::Absorption,
            "Absorption",
        );
        ui.selectable_value(
            &mut c.process,
            super::optical_spectrum::Process::Fluorescence,
            "Fluorescence",
        );
        if ui
            .selectable_value(
                &mut c.process,
                super::optical_spectrum::Process::Raman,
                "Resonance Raman",
            )
            .clicked()
        {
            c.grid.temperature = 0.;
        }
    });
    ui.horizontal(|ui| {
        ui.selectable_value(&mut c.engine, Engine::Internal, "Beavyr");
        ui.selectable_value(&mut c.engine, Engine::Orca, "ORCA ESD");
    });
    ui.horizontal(|ui| {
        ui.label("Electronic gap (ES minimum − GS minimum)");
        ui.add(
            egui::DragValue::new(&mut c.gap_ev)
                .speed(0.01)
                .range(0.00001..=1000.)
                .suffix(" eV"),
        );
    });
    ui.small("Electronic energies only; the Hessians supply the zero-point correction.");
    ui.horizontal_wrapped(|ui| {
        ui.label("Temperature");
        ui.add(
            egui::DragValue::new(&mut c.grid.temperature)
                .range(0.0..=10000.)
                .suffix(" K"),
        );
        ui.add_space(16.);
        ui.label("Grid spacing");
        ui.add(
            egui::DragValue::new(&mut c.grid.spacing_cm1)
                .range(0.01..=1000.)
                .suffix(" cm⁻¹"),
        );
    });
    ui.horizontal_wrapped(|ui| {
        ui.label("Lorentzian FWHM");
        ui.add(
            egui::DragValue::new(&mut c.grid.lorentz_fwhm_cm1)
                .range(0.0..=10000.)
                .suffix(" cm⁻¹"),
        );
        ui.add_space(16.);
        ui.label("Gaussian FWHM");
        ui.add(
            egui::DragValue::new(&mut c.grid.gaussian_fwhm_cm1)
                .range(0.0..=10000.)
                .suffix(" cm⁻¹"),
        );
    });
    if c.process == super::optical_spectrum::Process::Raman {
        ui.horizontal(|ui| {
            ui.label("Laser wavelengths (nm)");
            ui.text_edit_singleline(&mut c.raman.lasers);
        });
        ui.horizontal(|ui| {
            ui.selectable_value(&mut c.raman.order, 1, "Fundamentals");
            ui.selectable_value(&mut c.raman.order, 2, "+ Overtones / combinations");
            ui.label("Raman FWHM");
            ui.add(
                egui::DragValue::new(&mut c.raman.line_fwhm)
                    .range(0.01..=1000.)
                    .suffix(" cm⁻¹"),
            );
        });
        if c.engine == Engine::Internal {
            if c.properties.source == super::transition_properties::PropertySource::Manual {
                ui.horizontal(|ui| {
                    ui.label("Transition dipole (a.u.), x / y / z");
                    for v in &mut c.raman.dipole {
                        ui.add(egui::DragValue::new(v).speed(0.01));
                    }
                });
            }
            ui.checkbox(
                &mut c.raman.photon_factor,
                "Scattered-photon frequency⁴ factor",
            );
        }
    }
    property_controls(ui, c, ground);
    if c.engine == Engine::Internal {
        if c.process != super::optical_spectrum::Process::Raman {
            ui.checkbox(&mut c.frequency_factor, "Optical photon-energy prefactor");
        }
        ui.horizontal(|ui| {
            ui.label("Maximum time-grid points");
            egui::ComboBox::from_id_salt("harmonic_time_points")
                .selected_text(c.grid.max_points.to_string())
                .show_ui(ui, |ui| {
                    for p in [4096, 8192, 16384, 32768, 65536, 131072, 262144] {
                        ui.selectable_value(&mut c.grid.max_points, p, p.to_string());
                    }
                });
        });
    }
    if c.engine == Engine::Orca
        || matches!(
            c.properties.source,
            super::transition_properties::PropertySource::PyScf
                | super::transition_properties::PropertySource::Psi4
        )
    {
        if c.engine == Engine::Orca {
            ui.horizontal(|ui| {
                ui.label("ORCA executable");
                ui.text_edit_singleline(&mut c.orca_path);
                if ui.button("Browse…").clicked_once() {
                    if let Some(p) = rfd::FileDialog::new().pick_file() {
                        c.orca_path = p.display().to_string();
                        crate::qchem_interfaces::config::save_path(QcProgram::Orca, &c.orca_path);
                    }
                }
            });
        }
        ui.horizontal_wrapped(|ui| {
            ui.label("Charge");
            ui.add(egui::DragValue::new(&mut c.charge));
            ui.label("Multiplicity");
            ui.add(egui::DragValue::new(&mut c.multiplicity).range(1..=99));
            ui.label("Root");
            ui.add(egui::DragValue::new(&mut c.root).range(1..=1000));
            ui.label("Roots");
            ui.add(egui::DragValue::new(&mut c.roots).range(1..=1000));
            ui.checkbox(&mut c.tda, "TDA");
        });
        ui.horizontal(|ui| {
            ui.selectable_value(&mut c.qc_method, QcMethod::Hf, "HF");
            ui.selectable_value(&mut c.qc_method, QcMethod::Dft, "DFT");
            if c.qc_method == QcMethod::Dft {
                ui.label("Functional");
                ui.text_edit_singleline(&mut c.method.functional);
            }
            ui.label("Basis");
            ui.text_edit_singleline(&mut c.method.basis);
        });
        ui.horizontal(|ui| {
            ui.label("Cores");
            ui.add(egui::DragValue::new(&mut c.method.nproc).range(1..=4));
            ui.label("Memory");
            ui.add(
                egui::DragValue::new(&mut c.method.memory_mb)
                    .range(256..=1000000)
                    .suffix(" MB"),
            );
        });
        if c.engine == Engine::Orca {
            ui.horizontal(|ui| {
                ui.label("Extra keywords");
                ui.text_edit_singleline(&mut c.method.orca_extra);
            });
        }
    }
}
pub fn results(ui: &mut egui::Ui, output: &Output) {
    for note in &output.notes {
        ui.label(note);
    }
    if let Some(dir) = &output.workdir {
        ui.horizontal(|ui| {
            ui.label("Calculation files retained");
            if ui.button("Copy directory").clicked_once() {
                ui.ctx().copy_text(dir.clone());
            }
        });
    }
}

fn property_controls(ui: &mut egui::Ui, c: &mut Controls, ground: &super::input::HessianSource) {
    use super::transition_properties::PropertySource;
    ui.checkbox(
        &mut c.properties.ht,
        "Herzberg–Teller (transition-dipole derivatives)",
    );
    if c.engine == Engine::Orca {
        return;
    }
    if c.properties.ht || c.process == super::optical_spectrum::Process::Raman {
        ui.horizontal_wrapped(|ui| {
            ui.label("Transition properties");
            ui.selectable_value(
                &mut c.properties.source,
                PropertySource::Manual,
                "Manual dipole",
            );
            ui.selectable_value(&mut c.properties.source, PropertySource::File, "File");
            ui.selectable_value(&mut c.properties.source, PropertySource::PyScf, "PySCF");
            ui.selectable_value(&mut c.properties.source, PropertySource::Psi4, "Psi4");
        });
        if c.properties.source == PropertySource::File {
            ui.horizontal(|ui| {
                ui.text_edit_singleline(&mut c.properties.path);
                if ui.button("Browse…").clicked_once() {
                    if let Some(p) = rfd::FileDialog::new().pick_file() {
                        c.properties.path = p.display().to_string();
                    }
                }
                if ui.button("Save template…").clicked_once() {
                    let result = super::input::load(ground, 0.);
                    match result {
                        Ok(a) => {
                            if let Some(path) = rfd::FileDialog::new()
                                .set_file_name("transition.properties")
                                .save_file()
                            {
                                if let Err(e) =
                                    std::fs::write(path, super::transition_properties::template(&a))
                                {
                                    crate::runtime_warning::push(
                                        ui.ctx(),
                                        "Transition properties",
                                        &e.to_string(),
                                    );
                                }
                            }
                        }
                        Err(e) => {
                            crate::runtime_warning::push(ui.ctx(), "Transition properties", &e)
                        }
                    }
                }
            });
        }
        if matches!(
            c.properties.source,
            PropertySource::PyScf | PropertySource::Psi4
        ) {
            let module = if c.properties.source == PropertySource::PyScf {
                "pyscf"
            } else {
                "psi4"
            };
            crate::qchem_interfaces::python_env::environment_row(ui, module, false);
            if c.properties.ht {
                ui.horizontal(|ui| {
                    ui.label("Central-difference step");
                    ui.add(
                        egui::DragValue::new(&mut c.properties.step)
                            .range(0.00001..=0.1)
                            .speed(0.0001)
                            .suffix(" bohr"),
                    );
                });
            }
            ui.small("State overlap and phase are checked at each displacement; ambiguous matches stop the job.");
        }
    }
}
