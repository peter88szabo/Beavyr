//! The method block both job panels roll down under their program selector.
//!
//! One implementation, called from the optimizer and from the frequency
//! analysis, so the two cannot drift apart. Everything it decides -- which
//! fields apply, what blocks a run -- comes from `super::method`, which has no
//! GUI dependency and carries the tests; this file is only the drawing.

use bevy_egui::egui;

use super::method::{
    basis_label, functional, is_blocked, validate, BehemothMethod, Dispersion, MethodConfig,
    Severity, TwoElectron, XtbMethod, BASIS_GROUPS, FUNCTIONALS, RIJK_UNAVAILABLE,
};
use super::program::QcProgram;

/// The method line shown on the collapsed header, so the level of theory is
/// legible without opening the section.
pub fn method_summary(program: QcProgram, config: &MethodConfig) -> String {
    match program {
        QcProgram::Xtb => config.xtb.label().to_string(),
        QcProgram::Dreiding => "DREIDING force field".to_string(),
        QcProgram::Behemoth => match config.behemoth {
            BehemothMethod::Dft => {
                let name = functional(&config.functional)
                    .map(|f| f.label)
                    .unwrap_or(&config.functional);
                let fields = super::method::visible_fields(program, config);
                if fields.basis {
                    format!("{name}/{}", basis_label(&config.basis))
                } else {
                    // A composite names its own basis, so appending one would
                    // be wrong as well as redundant.
                    name.to_string()
                }
            }
            other => {
                let fields = super::method::visible_fields(program, config);
                if fields.basis {
                    format!("{}/{}", other.label(), basis_label(&config.basis))
                } else {
                    other.label().to_string()
                }
            }
        },
    }
}

/// Draws the method block and returns whether the configuration is blocked, so
/// the caller can disable its Run button.
///
/// `running` disables every control while a calculation is in flight, matching
/// the rest of both panels: changing the level of theory under a running job
/// would describe something other than what is running.
pub fn method_config_ui(
    ui: &mut egui::Ui,
    program: QcProgram,
    config: &mut MethodConfig,
    charge_multiplicity: i32,
    atoms: &[String],
    running: bool,
    id_prefix: &str,
) -> bool {
    let issues = validate(program, config, charge_multiplicity, atoms);
    let blocked = is_blocked(&issues);

    egui::CollapsingHeader::new(format!("Method \u{2014} {}", method_summary(program, config)))
        .id_salt(format!("{id_prefix}_method_block"))
        // Open by default: the level of theory is the thing most worth seeing
        // without having to go looking for it.
        .default_open(true)
        .show(ui, |ui| {
            match program {
                QcProgram::Dreiding => {
                    ui.label(
                        "A generic force field with no level of theory to set: its parameters \
                         follow from the atom types, which are perceived from the structure.",
                    );
                    ui.label(
                        egui::RichText::new(
                            "Mayo, Olafson & Goddard, J. Phys. Chem. 1990, 94, 8897.",
                        )
                        .small()
                        .weak(),
                    );
                }
                QcProgram::Xtb => {
                    ui.horizontal(|ui| {
                        ui.label("Method");
                        ui.add_enabled_ui(!running, |ui| {
                            egui::ComboBox::from_id_salt(format!("{id_prefix}_xtb_method"))
                                .selected_text(config.xtb.label())
                                .show_ui(ui, |ui| {
                                    for method in XtbMethod::ALL {
                                        ui.selectable_value(
                                            &mut config.xtb,
                                            method,
                                            method.label(),
                                        );
                                    }
                                });
                        });
                    });
                }
                QcProgram::Behemoth => {
                    ui.horizontal(|ui| {
                        ui.label("Method");
                        ui.add_enabled_ui(!running, |ui| {
                            egui::ComboBox::from_id_salt(format!("{id_prefix}_beh_method"))
                                .selected_text(config.behemoth.label())
                                .show_ui(ui, |ui| {
                                    for method in BehemothMethod::ALL {
                                        ui.selectable_value(
                                            &mut config.behemoth,
                                            method,
                                            method.label(),
                                        );
                                    }
                                });
                        });
                    });
                }
            }

            // Recomputed after the method dropdown, so a method changed this
            // frame shows the right fields immediately rather than a frame
            // late.
            let fields = super::method::visible_fields(program, config);

            if fields.functional {
                ui.horizontal(|ui| {
                    ui.label("Functional");
                    ui.add_enabled_ui(!running, |ui| {
                        let selected = functional(&config.functional)
                            .map(|f| f.label.to_string())
                            .unwrap_or_else(|| config.functional.clone());
                        egui::ComboBox::from_id_salt(format!("{id_prefix}_functional"))
                            .selected_text(selected)
                            .show_ui(ui, |ui| {
                                let mut separated = false;
                                for entry in &FUNCTIONALS {
                                    // The composites come after the common
                                    // functionals; mark where they start.
                                    if entry.composite.is_some() && !separated {
                                        ui.separator();
                                        ui.weak("composite (\u{201C}3c\u{201D})");
                                        separated = true;
                                    }
                                    let mut chosen = config.functional.clone();
                                    if ui
                                        .selectable_value(
                                            &mut chosen,
                                            entry.cli.to_string(),
                                            entry.label,
                                        )
                                        .clicked()
                                    {
                                        config.functional = entry.cli.to_string();
                                    }
                                }
                            });
                    });
                });
            }

            if fields.dispersion {
                ui.horizontal(|ui| {
                    ui.label("Dispersion");
                    ui.add_enabled_ui(!running, |ui| {
                        egui::ComboBox::from_id_salt(format!("{id_prefix}_dispersion"))
                            .selected_text(config.dispersion.label())
                            .show_ui(ui, |ui| {
                                for d in Dispersion::ALL {
                                    ui.selectable_value(&mut config.dispersion, d, d.label());
                                }
                            });
                    });
                });
            }

            if fields.basis {
                ui.horizontal(|ui| {
                    ui.label("Basis set");
                    ui.add_enabled_ui(!running, |ui| {
                        egui::ComboBox::from_id_salt(format!("{id_prefix}_basis"))
                            .selected_text(basis_label(&config.basis))
                            .show_ui(ui, |ui| {
                                for (i, group) in BASIS_GROUPS.iter().enumerate() {
                                    if i > 0 {
                                        ui.separator();
                                    }
                                    ui.weak(group.name);
                                    for (label, cli) in group.sets {
                                        let mut chosen = config.basis.clone();
                                        if ui
                                            .selectable_value(
                                                &mut chosen,
                                                cli.to_string(),
                                                *label,
                                            )
                                            .clicked()
                                        {
                                            config.basis = cli.to_string();
                                        }
                                    }
                                }
                            });
                    });
                });
            }

            if fields.two_electron {
                ui.horizontal(|ui| {
                    ui.label("Two-electron");
                    // Shown so the option is discoverable, permanently
                    // disabled because it cannot work here: RIJK has no
                    // analytic gradient, and both of these panels
                    // differentiate. The reason sits next to it rather than
                    // waiting for a failed run to explain itself.
                    ui.add_enabled_ui(false, |ui| {
                        egui::ComboBox::from_id_salt(format!("{id_prefix}_two_electron"))
                            .selected_text(config.two_electron.label())
                            .show_ui(ui, |ui| {
                                for t in TwoElectron::ALL {
                                    ui.selectable_value(&mut config.two_electron, t, t.label());
                                }
                            });
                    });
                    ui.weak("exact only \u{2014} RIJK has no gradient");
                });
            }

            if fields.memory || fields.nproc {
                ui.horizontal(|ui| {
                    if fields.memory {
                        ui.label("Memory (MB)");
                        ui.add_enabled(
                            !running,
                            egui::DragValue::new(&mut config.memory_mb)
                                .speed(64.0)
                                .range(64..=1_048_576),
                        );
                    }
                    if fields.nproc {
                        ui.label("nproc");
                        ui.add_enabled(
                            !running,
                            egui::DragValue::new(&mut config.nproc).speed(1.0).range(1..=1024),
                        );
                    }
                });
            }

            for issue in &issues {
                match issue.severity {
                    Severity::Block => {
                        ui.colored_label(
                            egui::Color32::from_rgb(230, 120, 90),
                            format!("\u{26a0} {}", issue.message),
                        );
                    }
                    Severity::Note => {
                        ui.weak(&issue.message);
                    }
                }
            }
            if fields.two_electron {
                ui.weak(RIJK_UNAVAILABLE);
            }
        });

    blocked
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::method::BehemothMethod;

    fn config(method: BehemothMethod) -> MethodConfig {
        MethodConfig { behemoth: method, ..Default::default() }
    }

    /// The collapsed header has to say what will run, since that is the whole
    /// point of showing a summary rather than just "Method".
    #[test]
    fn the_header_names_the_method_and_basis() {
        let mut dft = config(BehemothMethod::Dft);
        dft.functional = "pbe0".to_string();
        dft.basis = "def2-tzvp".to_string();
        assert_eq!(method_summary(QcProgram::Behemoth, &dft), "PBE0/def2-TZVP");

        let hf = config(BehemothMethod::Hf);
        assert_eq!(method_summary(QcProgram::Behemoth, &hf), "HF/def2-SVP");
    }

    /// A composite names its own basis set, so the summary must not append
    /// one -- it would be both redundant and wrong.
    #[test]
    fn a_composite_header_carries_no_basis() {
        let mut dft = config(BehemothMethod::Dft);
        dft.functional = "r2scan-3c".to_string();
        assert_eq!(method_summary(QcProgram::Behemoth, &dft), "r2SCAN-3c");
    }

    #[test]
    fn a_semi_empirical_header_is_just_the_method() {
        assert_eq!(
            method_summary(QcProgram::Behemoth, &config(BehemothMethod::Gfn1Xtb)),
            "GFN1-xTB"
        );
        assert_eq!(
            method_summary(QcProgram::Behemoth, &config(BehemothMethod::Tasi)),
            "TASI"
        );
    }

    #[test]
    fn the_xtb_header_is_its_parametrisation() {
        let mut c = MethodConfig::default();
        assert_eq!(method_summary(QcProgram::Xtb, &c), "GFN2-xTB");
        c.xtb = XtbMethod::GfnFf;
        assert_eq!(method_summary(QcProgram::Xtb, &c), "GFN-FF");
    }

    /// A hand-typed functional or basis still shows, rather than falling back
    /// to something misleading.
    #[test]
    fn an_unlisted_functional_shows_as_itself() {
        let mut dft = config(BehemothMethod::Dft);
        dft.functional = "gga_x_wc+gga_c_pbe".to_string();
        assert!(method_summary(QcProgram::Behemoth, &dft).starts_with("gga_x_wc+gga_c_pbe"));
    }
}
