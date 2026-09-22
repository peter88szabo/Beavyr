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
        // For ORCA and PySCF the functional *is* the method: there is no
        // separate method choice above it, so the summary is the level of
        // theory and nothing else.
        QcProgram::Orca => {
            let name = super::orca_method::functional_label(&config.functional);
            if super::orca_method::is_composite(&config.functional) {
                // A composite names its own basis, so appending one would be
                // wrong as well as redundant.
                name.to_string()
            } else {
                format!("{name}/{}", super::orca_method::basis_label(&config.basis))
            }
        }
        QcProgram::PySCF => format!(
            "{}/{}",
            super::pyscf_method::functional_label(&config.functional),
            super::pyscf_method::basis_label(&config.basis)
        ),
        QcProgram::SparrowPy => {
            super::sparrow_method::method_label(&config.functional).to_string()
        }
        QcProgram::Psi4 | QcProgram::Psi4Py => {
            let name = super::psi4_method::functional_label(&config.functional);
            if super::psi4_method::is_composite(&config.functional) {
                name.to_string()
            } else {
                format!("{name}/{}", super::psi4_method::basis_label(&config.basis))
            }
        }
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
                // Neither has a method choice above the functional: choosing
                // B3LYP *is* choosing the method. A "Method" dropdown with one
                // entry would be noise, so the functional and basis rows below
                // are the whole block.
                QcProgram::Orca => {
                    ui.label(
                        egui::RichText::new(
                            "The functional and basis set below are written onto ORCA's keyword \
                             line. Anything the lists do not cover \u{2014} a dispersion \
                             correction, tighter convergence \u{2014} goes in Extra keywords.",
                        )
                        .small()
                        .weak(),
                    );
                }
                QcProgram::PySCF => {
                    ui.label(
                        egui::RichText::new(
                            "Runs through the python3 of the environment Beavyr was started in. \
                             The functional accepts any LibXC specification, not only the \
                             listed names.",
                        )
                        .small()
                        .weak(),
                    );
                }
                QcProgram::Psi4 => {
                    ui.label(
                        egui::RichText::new(
                            "Psi4's own executable, which carries its own environment, so this \
                             route works whichever environment Beavyr was started from. A \
                             dispersion correction is part of the functional's name here, not a \
                             separate setting.",
                        )
                        .small()
                        .weak(),
                    );
                }
                QcProgram::Psi4Py => {
                    ui.label(
                        egui::RichText::new(
                            "The same engine as Psi4, reached through the python3 of the \
                             environment Beavyr was started in. Use the plain Psi4 entry instead \
                             if that environment does not have it.",
                        )
                        .small()
                        .weak(),
                    );
                }
                QcProgram::SparrowPy => {
                    ui.label(
                        egui::RichText::new(
                            "Semi-empirical parametrisations, through SCINE's Python bindings. \
                             Sparrow gives an energy, a gradient and a Hessian; Beavyr's own \
                             optimiser and frequency analysis do the rest, because Sparrow has \
                             neither.",
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

            // ORCA and PySCF each read their own catalogue, because the same
            // functional is spelled differently by different codes -- ORCA
            // rejects `M06-2X` and wants `M062X`, Behemoth is the other way
            // round. Sharing one list would mean translating between them,
            // which is exactly what per-program catalogues avoid.
            if fields.functional && program == QcProgram::Orca {
                orca_functional_row(ui, config, running, id_prefix);
            }
            if fields.functional && program == QcProgram::PySCF {
                pyscf_functional_row(ui, config, running, id_prefix);
            }
            if fields.functional && matches!(program, QcProgram::Psi4 | QcProgram::Psi4Py) {
                psi4_functional_row(ui, config, running, id_prefix);
            }
            if fields.functional && program == QcProgram::SparrowPy {
                sparrow_method_row(ui, config, running, id_prefix);
            }
            if fields.functional
                && !matches!(
                    program,
                    QcProgram::Orca
                        | QcProgram::PySCF
                        | QcProgram::Psi4
                        | QcProgram::Psi4Py
                        | QcProgram::SparrowPy
                )
            {
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

            if fields.basis && program == QcProgram::Orca {
                orca_basis_row(ui, config, running, id_prefix);
            }
            if fields.basis && program == QcProgram::PySCF {
                pyscf_basis_row(ui, config, running, id_prefix);
            }
            if fields.basis && matches!(program, QcProgram::Psi4 | QcProgram::Psi4Py) {
                psi4_basis_row(ui, config, running, id_prefix);
            }
            if fields.basis
                && !matches!(
                    program,
                    QcProgram::Orca
                        | QcProgram::PySCF
                        | QcProgram::Psi4
                        | QcProgram::Psi4Py
                        | QcProgram::SparrowPy
                )
            {
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

/// ORCA's functional dropdown, reading ORCA's own catalogue.
///
/// The composites are separated from the ordinary functionals, as in the
/// Behemoth block, because choosing one takes the basis field away and that is
/// worth signalling before the click rather than after it.
pub(crate) fn orca_functional_row(
    ui: &mut egui::Ui,
    config: &mut MethodConfig,
    running: bool,
    id_prefix: &str,
) {
    use super::orca_method;

    ui.horizontal(|ui| {
        ui.label("Functional");
        ui.add_enabled_ui(!running, |ui| {
            egui::ComboBox::from_id_salt(format!("{id_prefix}_orca_functional"))
                .selected_text(orca_method::functional_label(&config.functional))
                .show_ui(ui, |ui| {
                    // Grouped by family: forty-six entries in one flat list
                    // would be unreadable, and the family is what a chemist is
                    // choosing between before the individual name.
                    for (index, family) in orca_method::Family::ALL.into_iter().enumerate() {
                        let mut any = false;
                        for entry in orca_method::functionals_in(family) {
                            if !any {
                                if index > 0 {
                                    ui.separator();
                                }
                                ui.weak(family.heading());
                                any = true;
                            }
                            let mut chosen = config.functional.clone();
                            if ui
                                .selectable_value(
                                    &mut chosen,
                                    entry.keyword.to_string(),
                                    entry.label,
                                )
                                .clicked()
                            {
                                config.functional = entry.keyword.to_string();
                            }
                        }
                    }
                });
        });
    });

    // A double hybrid silently gains a correlation-fitting basis, which is the
    // only way ORCA will run one. Said out loud, because it changes the input
    // that gets written and a user reading the run directory should not have to
    // work out where the extra keyword came from.
    if let Some(aux) = orca_method::auxiliary_basis(&config.functional, &config.basis) {
        ui.label(
            egui::RichText::new(format!(
                "A double hybrid needs a correlation-fitting basis; {aux} is added automatically."
            ))
            .small()
            .weak(),
        );
    }
    if let Some(brings) = orca_method::composite_note(&config.functional) {
        ui.label(egui::RichText::new(format!("Brings: {brings}.")).small().weak());
    }

    orca_dispersion_row(ui, config, running, id_prefix);

    // Everything ORCA can be asked for that no dropdown covers.
    ui.horizontal(|ui| {
        ui.label("Extra keywords");
        ui.add_enabled(
            !running,
            egui::TextEdit::singleline(&mut config.orca_extra)
                .desired_width(220.0)
                .hint_text("e.g. D4 TightSCF"),
        );
    });
}

/// ORCA's dispersion dropdown.
///
/// Disabled, with the reason beside it, where the functional already carries a
/// correction: a composite brings its own, and a name ending in `-D3`, `-D4` or
/// `-V` is corrected by definition. Adding one on top would double-count, and
/// the resulting energy would look perfectly ordinary.
pub(crate) fn orca_dispersion_row(
    ui: &mut egui::Ui,
    config: &mut MethodConfig,
    running: bool,
    id_prefix: &str,
) {
    use super::orca_method::{self, Dispersion};

    let built_in = orca_method::is_composite(&config.functional)
        || orca_method::carries_own_dispersion(&config.functional);

    ui.horizontal(|ui| {
        ui.label("Dispersion");
        ui.add_enabled_ui(!running && !built_in, |ui| {
            egui::ComboBox::from_id_salt(format!("{id_prefix}_orca_dispersion"))
                .selected_text(if built_in {
                    "built in"
                } else {
                    config.orca_dispersion.label()
                })
                .show_ui(ui, |ui| {
                    for entry in Dispersion::ALL {
                        ui.selectable_value(&mut config.orca_dispersion, entry, entry.label());
                    }
                });
        });
        if built_in {
            ui.label(
                egui::RichText::new("this functional already includes one")
                    .small()
                    .weak(),
            );
        }
    });
}

/// ORCA's basis dropdown, reading ORCA's own catalogue.
pub(crate) fn orca_basis_row(ui: &mut egui::Ui, config: &mut MethodConfig, running: bool, id_prefix: &str) {
    use super::orca_method;

    ui.horizontal(|ui| {
        ui.label("Basis set");
        ui.add_enabled_ui(!running, |ui| {
            egui::ComboBox::from_id_salt(format!("{id_prefix}_orca_basis"))
                .selected_text(orca_method::basis_label(&config.basis))
                .show_ui(ui, |ui| {
                    for (index, group) in orca_method::BASIS_GROUPS.iter().enumerate() {
                        if index > 0 {
                            ui.separator();
                        }
                        ui.weak(group.name);
                        for (label, keyword) in group.sets {
                            let mut chosen = config.basis.clone();
                            if ui
                                .selectable_value(&mut chosen, keyword.to_string(), *label)
                                .clicked()
                            {
                                config.basis = keyword.to_string();
                            }
                        }
                    }
                });
        });
    });
}

/// Sparrow's parametrisation dropdown.
///
/// Labelled "Parametrisation" rather than "Functional": there is no functional
/// here and no basis set either. Sparrow is a set of parameter fits, and
/// calling the choice something it is not would be its own small lie.
pub(crate) fn sparrow_method_row(
    ui: &mut egui::Ui,
    config: &mut MethodConfig,
    running: bool,
    id_prefix: &str,
) {
    use super::sparrow_method;

    ui.horizontal(|ui| {
        ui.label("Parametrisation");
        ui.add_enabled_ui(!running, |ui| {
            egui::ComboBox::from_id_salt(format!("{id_prefix}_sparrow_method"))
                .selected_text(sparrow_method::method_label(&config.functional))
                .show_ui(ui, |ui| {
                    for (index, family) in sparrow_method::Family::ALL.into_iter().enumerate() {
                        if index > 0 {
                            ui.separator();
                        }
                        ui.weak(family.heading());
                        for entry in sparrow_method::methods_in(family) {
                            let mut chosen = config.functional.clone();
                            if ui
                                .selectable_value(&mut chosen, entry.name.to_string(), entry.label)
                                .clicked()
                            {
                                config.functional = entry.name.to_string();
                            }
                        }
                    }
                });
        });
    });
}

/// Psi4's functional dropdown, reading Psi4's own catalogue.
///
/// No dispersion row goes with it: Psi4 asks for a correction inside the
/// functional's own name (`b97-d3bj`, `wb97x-d3bj`), so a separate control
/// would be a second way to say the same thing, and the two could disagree.
pub(crate) fn psi4_functional_row(
    ui: &mut egui::Ui,
    config: &mut MethodConfig,
    running: bool,
    id_prefix: &str,
) {
    use super::psi4_method;

    ui.horizontal(|ui| {
        ui.label("Functional");
        ui.add_enabled_ui(!running, |ui| {
            egui::ComboBox::from_id_salt(format!("{id_prefix}_psi4_functional"))
                .selected_text(psi4_method::functional_label(&config.functional))
                .show_ui(ui, |ui| {
                    for (index, family) in psi4_method::Family::ALL.into_iter().enumerate() {
                        let mut any = false;
                        for entry in psi4_method::functionals_in(family) {
                            if !any {
                                if index > 0 {
                                    ui.separator();
                                }
                                ui.weak(family.heading());
                                any = true;
                            }
                            let mut chosen = config.functional.clone();
                            if ui
                                .selectable_value(&mut chosen, entry.name.to_string(), entry.label)
                                .clicked()
                            {
                                config.functional = entry.name.to_string();
                            }
                        }
                    }
                });
        });
    });

    if let Some(brings) = psi4_method::composite_note(&config.functional) {
        ui.label(egui::RichText::new(format!("Brings: {brings}.")).small().weak());
    }
}

/// Psi4's basis dropdown, reading Psi4's own catalogue.
pub(crate) fn psi4_basis_row(
    ui: &mut egui::Ui,
    config: &mut MethodConfig,
    running: bool,
    id_prefix: &str,
) {
    use super::psi4_method;

    ui.horizontal(|ui| {
        ui.label("Basis set");
        ui.add_enabled_ui(!running, |ui| {
            egui::ComboBox::from_id_salt(format!("{id_prefix}_psi4_basis"))
                .selected_text(psi4_method::basis_label(&config.basis))
                .show_ui(ui, |ui| {
                    for (index, group) in psi4_method::BASIS_GROUPS.iter().enumerate() {
                        if index > 0 {
                            ui.separator();
                        }
                        ui.weak(group.name);
                        for (label, spelling) in group.sets {
                            let mut chosen = config.basis.clone();
                            if ui
                                .selectable_value(&mut chosen, spelling.to_string(), *label)
                                .clicked()
                            {
                                config.basis = spelling.to_string();
                            }
                        }
                    }
                });
        });
    });
}

/// PySCF's functional dropdown, reading PySCF's own catalogue.
pub(crate) fn pyscf_functional_row(
    ui: &mut egui::Ui,
    config: &mut MethodConfig,
    running: bool,
    id_prefix: &str,
) {
    use super::pyscf_method;

    ui.horizontal(|ui| {
        ui.label("Functional");
        ui.add_enabled_ui(!running, |ui| {
            egui::ComboBox::from_id_salt(format!("{id_prefix}_pyscf_functional"))
                .selected_text(pyscf_method::functional_label(&config.functional))
                .show_ui(ui, |ui| {
                    // Grouped by family, as ORCA's is and for the same reason.
                    for (index, family) in pyscf_method::Family::ALL.into_iter().enumerate() {
                        let mut any = false;
                        for entry in pyscf_method::functionals_in(family) {
                            if !any {
                                if index > 0 {
                                    ui.separator();
                                }
                                ui.weak(family.heading());
                                any = true;
                            }
                            let mut chosen = config.functional.clone();
                            if ui
                                .selectable_value(&mut chosen, entry.cli.to_string(), entry.label)
                                .clicked()
                            {
                                config.functional = entry.cli.to_string();
                            }
                        }
                    }
                });
        });
    });
}

/// PySCF's basis dropdown, reading PySCF's own catalogue.
pub(crate) fn pyscf_basis_row(ui: &mut egui::Ui, config: &mut MethodConfig, running: bool, id_prefix: &str) {
    use super::pyscf_method;

    ui.horizontal(|ui| {
        ui.label("Basis set");
        ui.add_enabled_ui(!running, |ui| {
            egui::ComboBox::from_id_salt(format!("{id_prefix}_pyscf_basis"))
                .selected_text(pyscf_method::basis_label(&config.basis))
                .show_ui(ui, |ui| {
                    for (index, group) in pyscf_method::BASIS_GROUPS.iter().enumerate() {
                        if index > 0 {
                            ui.separator();
                        }
                        ui.weak(group.name);
                        for (label, cli) in group.sets {
                            let mut chosen = config.basis.clone();
                            if ui
                                .selectable_value(&mut chosen, cli.to_string(), *label)
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
