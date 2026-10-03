//! Modern calculation setup. Presentation state lives in egui; chemistry and
//! running jobs remain shared with the Classic form in `ui`.

use crate::ui_style::ResponseExt;

use super::*;
use crate::molecule::dummy_count;
use crate::qchem_interfaces::{method, orca_method, psi4_method, pyscf_method, sparrow_method};
use std::path::PathBuf;

#[path = "input_editor.rs"]
mod input_editor;
pub(crate) use input_editor::show as input_editor_window;

const ACCENT: egui::Color32 = egui::Color32::from_rgb(104, 174, 255);
const WARNING: egui::Color32 = egui::Color32::from_rgb(230, 170, 90);
const ERROR: egui::Color32 = egui::Color32::from_rgb(230, 120, 90);
const WIDTH: f32 = 500.0;
const PREFS_ID: &str = "qc_panel_presentation_v1";

#[derive(Clone)]
struct Presentation {
    modern: bool,
    docked: bool,
    configure: bool,
    method_search: String,
    functional_search: String,
    basis_search: String,
    binary_check: Option<(QcProgram, String, Result<PathBuf, String>)>,
}

impl Default for Presentation {
    fn default() -> Self {
        Self {
            modern: std::env::var_os("BEAVYR_CLASSIC_QC").is_none(),
            docked: false,
            configure: false,
            method_search: String::new(),
            functional_search: String::new(),
            basis_search: String::new(),
            binary_check: None,
        }
    }
}

fn presentation(ctx: &egui::Context) -> Presentation {
    ctx.data_mut(|data| {
        data.get_temp_mut_or_default::<Presentation>(egui::Id::new(PREFS_ID))
            .clone()
    })
}

fn store(ctx: &egui::Context, state: Presentation) {
    ctx.data_mut(|data| data.insert_temp(egui::Id::new(PREFS_ID), state));
}

pub(super) fn is_modern(ctx: &egui::Context) -> bool {
    presentation(ctx).modern
}

pub(crate) fn is_docked(ctx: &egui::Context) -> bool {
    let state = presentation(ctx);
    state.modern && state.docked
}

pub(crate) fn docked_width(ctx: &egui::Context) -> f32 {
    ctx.data_mut(|data| {
        data.get_temp::<f32>(egui::Id::new("qc_panel_docked_width"))
            .unwrap_or(WIDTH)
    })
}

pub(super) fn layout_switch(ui: &mut egui::Ui) {
    let mut state = presentation(ui.ctx());
    ui.horizontal(|ui| {
        ui.small("Layout");
        ui.selectable_value(&mut state.modern, true, "Modern");
        ui.selectable_value(&mut state.modern, false, "Classic");
    });
    store(ui.ctx(), state);
    ui.separator();
}

pub(super) fn show(
    ctx: &egui::Context,
    host_ui: &mut egui::Ui,
    panel: &mut QcPanelState,
    task: &mut QcRunTask,
    opt_panel: &mut XtbPanelState,
    mol: &Molecule,
) -> Option<egui::Rect> {
    let mut open = panel.open;
    let mut close = false;
    let rect = if is_docked(ctx) {
        let result = egui::Panel::right("modern_quantum_chemistry")
            .default_size(WIDTH)
            .min_size(440.0)
            .max_size(720.0)
            .resizable(true)
            .show(host_ui, |ui| {
                close = contents(ui, panel, task, opt_panel, mol)
            });
        ctx.data_mut(|data| {
            data.insert_temp(
                egui::Id::new("qc_panel_docked_width"),
                result.response.rect.width(),
            )
        });
        Some(result.response.rect)
    } else {
        let viewport = ctx.viewport_rect();
        egui::Window::new("Quantum Chemistry")
            .id(egui::Id::new("modern_quantum_chemistry_window_v1"))
            .open(&mut open)
            .resizable(true)
            .default_size([WIDTH, 760.0_f32.min(viewport.height() - 32.0)])
            .min_width(440.0)
            .min_height(540.0)
            .default_pos([viewport.left() + 64.0, viewport.top() + 16.0])
            .show(ctx, |ui| close = contents(ui, panel, task, opt_panel, mol))
            .map(|window| window.response.rect)
    };
    panel.open = open && !close;
    rect
}

fn style(ui: &mut egui::Ui) {
    let style = ui.style_mut();
    style.spacing.item_spacing = egui::vec2(8.0, 8.0);
    style.spacing.button_padding = egui::vec2(10.0, 6.0);
    style.spacing.interact_size.y = 30.0;
    style.visuals.selection.bg_fill = egui::Color32::from_rgb(35, 75, 117);
    style.visuals.selection.stroke = egui::Stroke::new(1.0, ACCENT);
    for widget in [
        &mut style.visuals.widgets.inactive,
        &mut style.visuals.widgets.hovered,
        &mut style.visuals.widgets.active,
    ] {
        widget.corner_radius = egui::CornerRadius::same(6);
    }
}

fn card(ui: &egui::Ui) -> egui::Frame {
    egui::Frame::new()
        .fill(ui.visuals().faint_bg_color)
        .stroke(ui.visuals().widgets.noninteractive.bg_stroke)
        .corner_radius(8)
        .inner_margin(10)
}

struct CodeStatus {
    ready: bool,
    label: &'static str,
    detail: String,
}

fn code_status(
    state: &mut Presentation,
    panel: &QcPanelState,
    opt_panel: &XtbPanelState,
) -> CodeStatus {
    if let Some(module) = panel.program.python_module() {
        let environment = crate::qchem_interfaces::python_env::cached(module);
        return CodeStatus {
            ready: environment.is_ready(),
            label: if environment.is_ready() {
                "Python ready"
            } else {
                "Setup needed"
            },
            detail: environment.detail(),
        };
    }
    if !panel.program.needs_binary_path() {
        return CodeStatus {
            ready: true,
            label: "Built in",
            detail: "DREIDING runs inside Beavyr.".into(),
        };
    }
    let path = opt_panel
        .paths
        .iter()
        .find(|(program, _)| *program == panel.program)
        .map(|(_, path)| path.as_str())
        .unwrap_or("");
    if !state
        .binary_check
        .as_ref()
        .is_some_and(|(program, checked, _)| *program == panel.program && checked == path)
    {
        state.binary_check = Some((
            panel.program,
            path.to_owned(),
            resolve_program_executable(panel.program, path),
        ));
    }
    match &state.binary_check.as_ref().unwrap().2 {
        Ok(path) => CodeStatus {
            ready: true,
            label: "Executable found",
            detail: path.display().to_string(),
        },
        Err(error) => CodeStatus {
            ready: false,
            label: "Setup needed",
            detail: error.clone(),
        },
    }
}

fn header(
    ui: &mut egui::Ui,
    state: &mut Presentation,
    panel: &mut QcPanelState,
    running: bool,
) -> bool {
    let mut close = false;
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new("Calculation setup").size(18.0).strong());
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            close = ui
                .small_button("×")
                .on_hover_text("Close setup; a running calculation continues")
                .clicked_once();
            if ui
                .small_button(if state.docked { "Float" } else { "Dock" })
                .clicked_once()
            {
                state.docked = !state.docked;
            }
            ui.selectable_value(&mut state.modern, false, "Classic")
                .on_hover_text("Use the previous form with the same settings and running job");
            ui.selectable_value(&mut state.modern, true, "Modern");
        });
    });

    // The code chooser is deliberately outside the scrolling settings.
    ui.strong("Quantum chemistry code");
    ui.horizontal(|ui| {
        let mut chosen = panel.program;
        ui.add_enabled_ui(!running, |ui| {
            egui::ComboBox::from_id_salt("modern_qc_program")
                .width((ui.available_width() - 116.0).max(180.0))
                .selected_text(panel.program.label())
                .show_ui(ui, |ui| {
                    for program in QcProgram::ALL {
                        ui.selectable_value(&mut chosen, program, program.label());
                    }
                });
        });
        if chosen != panel.program {
            panel.set_program(chosen);
            state.method_search.clear();
            state.functional_search.clear();
            state.basis_search.clear();
            state.binary_check = None;
        }
        if ui.button("Configure…").clicked_once() {
            state.configure = !state.configure;
        }
    });
    close
}

fn electronic_state_row(
    ui: &mut egui::Ui,
    panel: &mut QcPanelState,
    mol: &Molecule,
    running: bool,
) {
    charge_and_multiplicity(ui, panel, running);
    let dummies = dummy_count(&mol.atoms);
    if dummies > 0 {
        ui.small(format!(
            "{dummies} dummy atom(s) excluded from the calculation"
        ));
    }
}

fn job_choices(ui: &mut egui::Ui, panel: &mut QcPanelState, running: bool) {
    card(ui).show(ui, |ui| {
        ui.set_min_width((ui.available_width() - 1.0).max(0.0));
        ui.horizontal(|ui| {
            ui.strong("What to calculate");
            ui.add_enabled_ui(!running, |ui| {
                ui.menu_button("More…", |ui| {
                    for job in JobType::ALL {
                        let reason = job::unsupported_reason(panel.program, job);
                        let response = ui.add_enabled(
                            reason.is_none(),
                            egui::Button::new(job.label()).selected(panel.job == job),
                        );
                        if let Some(reason) = reason {
                            response.clone().on_disabled_hover_text(reason);
                        }
                        if response.clicked_once() {
                            panel.job = job;
                            ui.close();
                        }
                    }
                });
            });
        });
        let width = (ui.available_width() - 8.0) / 2.0;
        egui::Grid::new("modern_qc_job_tiles")
            .spacing([8.0, 6.0])
            .show(ui, |ui| {
                for (index, (job, label)) in [
                    (JobType::SinglePoint, "Energy"),
                    (JobType::Optimize, "Optimize"),
                    (JobType::Frequencies, "Frequencies"),
                    (JobType::OptimizeThenFrequencies, "Optimize + frequencies"),
                ]
                .into_iter()
                .enumerate()
                {
                    let reason = job::unsupported_reason(panel.program, job);
                    let response = ui.add_enabled(
                        !running && reason.is_none(),
                        egui::Button::new(label)
                            .selected(panel.job == job)
                            .min_size([width, 34.0].into()),
                    );
                    if let Some(reason) = reason {
                        response.clone().on_disabled_hover_text(reason);
                    }
                    if response.clicked_once() {
                        panel.job = job;
                    }
                    if index % 2 == 1 {
                        ui.end_row();
                    }
                }
            });
        ui.add(
            egui::Label::new(egui::RichText::new(panel.job.description()).small().weak())
                .truncate(),
        )
        .on_hover_text(panel.job.description());
        if panel.needs_excited_state_options() {
            excited_state_row(ui, panel, running);
        }
    });
}

#[derive(Clone)]
struct Choice {
    label: &'static str,
    keyword: &'static str,
    group: &'static str,
}

impl Choice {
    fn matches(&self, query: &str) -> bool {
        let query = query.trim().to_lowercase();
        self.label.to_lowercase().contains(&query)
            || self.keyword.to_lowercase().contains(&query)
            || self.group.to_lowercase().contains(&query)
    }
}

fn picker(
    ui: &mut egui::Ui,
    label: &str,
    id: &str,
    value: &mut String,
    choices: &[Choice],
    search: &mut String,
    allow_custom: bool,
    enabled: bool,
) {
    ui.horizontal(|ui| {
        ui.add_sized([86.0, 30.0], egui::Label::new(label));
        let selected = choices
            .iter()
            .find(|choice| choice.keyword == value)
            .map(|choice| choice.label)
            .unwrap_or(value.as_str())
            .to_owned();
        ui.add_enabled_ui(enabled, |ui| {
            egui::ComboBox::from_id_salt(id)
                .width(ui.available_width())
                .selected_text(selected)
                .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
                .show_ui(ui, |ui| {
                    ui.add(
                        egui::TextEdit::singleline(search)
                            .hint_text("Search by name or keyword")
                            .desired_width(f32::INFINITY),
                    );
                    egui::ScrollArea::vertical()
                        .max_height(280.0)
                        .show(ui, |ui| {
                            let mut group = None;
                            let mut found = false;
                            for (index, choice) in choices
                                .iter()
                                .enumerate()
                                .filter(|(_, choice)| choice.matches(search))
                            {
                                found = true;
                                if group != Some(choice.group) && !choice.group.is_empty() {
                                    ui.small(choice.group);
                                    group = Some(choice.group);
                                }
                                ui.push_id(index, |ui| {
                                    if ui
                                        .selectable_label(*value == choice.keyword, choice.label)
                                        .clicked_once()
                                    {
                                        *value = choice.keyword.to_owned();
                                        ui.close();
                                    }
                                });
                            }
                            if !found {
                                ui.small("No matching entries");
                            }
                            if allow_custom
                                && !search.trim().is_empty()
                                && !choices.iter().any(|choice| {
                                    choice.keyword.eq_ignore_ascii_case(search.trim())
                                })
                            {
                                if ui
                                    .button(format!("Use custom: {}", search.trim()))
                                    .clicked_once()
                                {
                                    *value = search.trim().to_owned();
                                    ui.close();
                                }
                            }
                        });
                });
        });
    });
}

fn functional_choices(program: QcProgram) -> Vec<Choice> {
    match program {
        QcProgram::Orca => orca_method::FUNCTIONALS
            .iter()
            .map(|entry| Choice {
                label: entry.label,
                keyword: entry.keyword,
                group: entry.family.heading(),
            })
            .collect(),
        QcProgram::PySCF => pyscf_method::FUNCTIONALS
            .iter()
            .map(|entry| Choice {
                label: entry.label,
                keyword: entry.cli,
                group: entry.family.heading(),
            })
            .collect(),
        QcProgram::Psi4 | QcProgram::Psi4Py => psi4_method::FUNCTIONALS
            .iter()
            .map(|entry| Choice {
                label: entry.label,
                keyword: entry.name,
                group: entry.family.heading(),
            })
            .collect(),
        QcProgram::SparrowPy => sparrow_method::METHODS
            .iter()
            .map(|entry| Choice {
                label: entry.label,
                keyword: entry.name,
                group: entry.family.heading(),
            })
            .collect(),
        _ => method::FUNCTIONALS
            .iter()
            .map(|entry| Choice {
                label: entry.label,
                keyword: entry.cli,
                group: if entry.composite.is_some() {
                    "Composite"
                } else {
                    "Functionals"
                },
            })
            .collect(),
    }
}

fn basis_choices(program: QcProgram, f12: bool) -> Vec<Choice> {
    if f12 {
        return orca_run::F12_BASIS_SETS
            .iter()
            .map(|(label, keyword)| Choice {
                label,
                keyword,
                group: "F12 basis sets",
            })
            .collect();
    }
    match program {
        QcProgram::Orca => orca_method::BASIS_GROUPS
            .iter()
            .flat_map(|group| {
                group.sets.iter().map(move |(label, keyword)| Choice {
                    label,
                    keyword,
                    group: group.name,
                })
            })
            .collect(),
        QcProgram::PySCF => pyscf_method::BASIS_GROUPS
            .iter()
            .flat_map(|group| {
                group.sets.iter().map(move |(label, keyword)| Choice {
                    label,
                    keyword,
                    group: group.name,
                })
            })
            .collect(),
        QcProgram::Psi4 | QcProgram::Psi4Py => psi4_method::BASIS_GROUPS
            .iter()
            .flat_map(|group| {
                group.sets.iter().map(move |(label, keyword)| Choice {
                    label,
                    keyword,
                    group: group.name,
                })
            })
            .collect(),
        _ => method::BASIS_GROUPS
            .iter()
            .flat_map(|group| {
                group.sets.iter().map(move |(label, keyword)| Choice {
                    label,
                    keyword,
                    group: group.name,
                })
            })
            .collect(),
    }
}

fn theory_card(
    ui: &mut egui::Ui,
    state: &mut Presentation,
    panel: &mut QcPanelState,
    running: bool,
) {
    card(ui).show(ui, |ui| {
        ui.set_min_width((ui.available_width() - 1.0).max(0.0));
        ui.strong("Method and basis");
        let methods = job::methods_for(panel.program);
        if methods.is_empty() {
            ui.small("DREIDING force field · parameters follow the atom types");
            return;
        }
        let choices: Vec<_> = methods
            .iter()
            .map(|method| Choice {
                label: method.label(),
                keyword: method.label(),
                group: "",
            })
            .collect();
        let mut selected = panel
            .method
            .map(|method| method.label())
            .unwrap_or("Choose a method")
            .to_owned();
        picker(
            ui,
            "Method",
            "modern_qc_method",
            &mut selected,
            &choices,
            &mut state.method_search,
            false,
            !running,
        );
        if let Some(&method) = methods.iter().find(|method| method.label() == selected) {
            if panel.method != Some(method) {
                panel.method = Some(method);
                adopt_method_defaults(panel, method);
                state.basis_search.clear();
            }
        }

        if panel.needs_functional() {
            picker(
                ui,
                if panel.program == QcProgram::SparrowPy {
                    "Parameters"
                } else {
                    "Functional"
                },
                "modern_qc_functional",
                &mut panel.config.functional,
                &functional_choices(panel.program),
                &mut state.functional_search,
                true,
                !running,
            );
            if panel.program == QcProgram::Orca {
                if let Some(note) = orca_method::composite_note(&panel.config.functional) {
                    ui.add(
                        egui::Label::new(
                            egui::RichText::new(format!("Includes: {note}"))
                                .small()
                                .weak(),
                        )
                        .truncate(),
                    )
                    .on_hover_text(note);
                }
                let built_in = orca_method::is_composite(&panel.config.functional)
                    || orca_method::carries_own_dispersion(&panel.config.functional);
                ui.horizontal(|ui| {
                    ui.add_sized([86.0, 30.0], egui::Label::new("Dispersion"));
                    ui.add_enabled_ui(!running && !built_in, |ui| {
                        egui::ComboBox::from_id_salt("modern_qc_dispersion")
                            .width(ui.available_width())
                            .selected_text(if built_in {
                                "Included in functional"
                            } else {
                                panel.config.orca_dispersion.label()
                            })
                            .show_ui(ui, |ui| {
                                for value in orca_method::Dispersion::ALL {
                                    ui.selectable_value(
                                        &mut panel.config.orca_dispersion,
                                        value,
                                        value.label(),
                                    );
                                }
                            });
                    })
                    .response
                    .on_disabled_hover_text(
                        "This functional carries its own dispersion correction.",
                    );
                });
            }
        }
        if panel.needs_basis() {
            let f12 = panel.method.is_some_and(|method| method.is_f12());
            picker(
                ui,
                "Basis set",
                "modern_qc_basis",
                &mut panel.config.basis,
                &basis_choices(panel.program, f12),
                &mut state.basis_search,
                !f12,
                !running,
            );
            if panel.program == QcProgram::Orca {
                let auxiliary = panel
                    .method
                    .map(|method| {
                        orca_run::basis_keywords(
                            method,
                            &panel.config.functional,
                            &panel.config.basis,
                        )
                    })
                    .unwrap_or_default();
                if auxiliary.len() > 1 {
                    let text = format!(
                        "Fitting basis added automatically: {}",
                        auxiliary[1..].join(", ")
                    );
                    ui.add(egui::Label::new(egui::RichText::new(&text).small().weak()).truncate())
                        .on_hover_text(text);
                }
            }
        }
        if panel.method.is_some_and(|method| method.is_expensive()) {
            ui.label(
                egui::RichText::new("Computationally demanding method")
                    .small()
                    .color(WARNING),
            );
        }
    });
}

fn configuration(
    ui: &mut egui::Ui,
    state: &mut Presentation,
    panel: &mut QcPanelState,
    opt_panel: &mut XtbPanelState,
    running: bool,
) {
    let response = egui::CollapsingHeader::new(format!("Configure {}", panel.program.label()))
        .id_salt("modern_qc_configuration")
        .open(Some(state.configure))
        .show(ui, |ui| {
            program_path_row(ui, panel, opt_panel, running);
            if panel.program.python_module().is_some() {
                python_environment_row(ui, panel.program, running);
            } else if !panel.program.needs_binary_path() {
                ui.small("Built into Beavyr; no external executable is needed.");
            }
            if panel.program.needs_binary_path()
                && ui
                    .add_enabled(!running, egui::Button::new("Check again"))
                    .clicked_once()
            {
                state.binary_check = None;
            }
        });
    if response.header_response.clicked_once() {
        state.configure = !state.configure;
    }
}

fn resources_row(ui: &mut egui::Ui, panel: &mut QcPanelState, running: bool) {
    ui.add_enabled_ui(!running, |ui| {
        ui.horizontal(|ui| {
            ui.label("CPU cores");
            ui.add(
                egui::DragValue::new(&mut panel.config.nproc)
                    .range(1..=4)
                    .speed(1.0),
            );
            ui.label("Memory (MB)");
            let memory = ui.add(
                egui::DragValue::new(&mut panel.config.memory_mb)
                    .range(256..=1_048_576)
                    .speed(64.0),
            );
            if panel.program == QcProgram::Orca {
                memory.on_hover_text(format!(
                    "ORCA receives {} MB per core",
                    (panel.config.memory_mb / panel.config.nproc).max(256)
                ));
            }
        });
    });
}

fn optional_settings(ui: &mut egui::Ui, panel: &mut QcPanelState, running: bool) {
    egui::CollapsingHeader::new("Advanced options")
        .id_salt("modern_qc_advanced")
        .show(ui, |ui| {
            wavefunction_row(ui, panel, running);
            if panel.program == QcProgram::Orca {
                ui.label("Extra keywords");
                ui.add_enabled(
                    !running,
                    egui::TextEdit::singleline(&mut panel.config.orca_extra)
                        .desired_width(f32::INFINITY)
                        .hint_text("e.g. TightSCF CPCM(water)"),
                );
            }
        });
}

fn start_label(job: JobType) -> &'static str {
    match job {
        JobType::SinglePoint => "Calculate energy",
        JobType::Gradient => "Calculate gradient",
        JobType::Hessian => "Calculate Hessian",
        JobType::Optimize => "Start optimization",
        JobType::Frequencies => "Calculate frequencies",
        JobType::OptimizeThenFrequencies => "Optimize + calculate frequencies",
        JobType::TransitionState => "Start transition-state search",
        JobType::Irc => "Start IRC calculation",
        JobType::TdDft => "Calculate excited states",
    }
}

fn blocking_reason(panel: &QcPanelState, mol: &Molecule, code: &CodeStatus) -> Option<String> {
    if mol.atoms.is_empty() {
        return Some("Load or build a molecule first".into());
    }
    if let Err(error) = validate_electronic_state(mol, panel.charge, panel.multiplicity) {
        return Some(error.to_string());
    }
    if !code.ready {
        return Some(code.detail.clone());
    }
    if let Some(reason) = job::unsupported_reason(panel.program, panel.job) {
        return Some(reason.into());
    }
    if !job::methods_for(panel.program).is_empty()
        && !panel
            .method
            .is_some_and(|method| job::methods_for(panel.program).contains(&method))
    {
        return Some("Choose a method supported by this code".into());
    }
    if panel.has_input_file() && panel.advanced && panel.advanced_text.trim().is_empty() {
        return Some("The custom input is empty".into());
    }
    None
}

fn footer(
    ui: &mut egui::Ui,
    panel: &mut QcPanelState,
    task: &mut QcRunTask,
    opt_panel: &XtbPanelState,
    mol: &Molecule,
    code: &CodeStatus,
) {
    ui.separator();
    input_editor::entry(ui, panel, mol, task.is_running());
    let summary = if panel.has_input_file() && panel.advanced {
        "Custom input is used exactly as written".to_owned()
    } else {
        panel.summary_line()
    };
    ui.add(
        egui::Label::new(egui::RichText::new(&summary).small().color(
            if panel.has_input_file() && panel.advanced {
                WARNING
            } else {
                ACCENT
            },
        ))
        .truncate(),
    )
    .on_hover_text(summary);
    let running = task.is_running();
    let reason = input_editor::program_mismatch(ui.ctx(), panel)
        .or_else(|| blocking_reason(panel, mol, code));
    let status = if running {
        format!(
            "Running · {}",
            task.elapsed()
                .map(crate::qchem_interfaces::xtb_optimize::format_elapsed)
                .unwrap_or_default()
        )
    } else if let Some(reason) = &reason {
        reason.clone()
    } else if task.last_is_error {
        format!(
            "Last run failed: {}",
            task.last_message
                .as_deref()
                .unwrap_or("See the calculation output")
        )
    } else if let Some(message) = &task.last_message {
        message.clone()
    } else {
        "Ready to run".to_owned()
    };
    ui.add(
        egui::Label::new(egui::RichText::new(&status).small().color(
            if reason.is_some() || task.last_is_error {
                ERROR
            } else {
                ui.visuals().text_color()
            },
        ))
        .truncate(),
    )
    .on_hover_text(&status);
    if running {
        ui.horizontal(|ui| {
            ui.spinner();
            ui.label("Calculation running");
            if ui.button("Cancel calculation").clicked_once() {
                task.cancel();
            }
        });
    } else {
        let response = ui.add_enabled(
            reason.is_none(),
            egui::Button::new(
                egui::RichText::new(start_label(panel.job)).color(egui::Color32::WHITE),
            )
            .fill(egui::Color32::from_rgb(38, 96, 159))
            .min_size([ui.available_width(), 36.0].into()),
        );
        if let Some(reason) = &reason {
            response.clone().on_disabled_hover_text(reason);
        }
        if response.clicked_once() {
            start_run(panel, task, opt_panel, mol);
        }
    }
    if task.last.is_some() {
        if ui
            .button(if panel.summary_open {
                "Hide results"
            } else {
                "View results ↗"
            })
            .clicked_once()
        {
            panel.summary_open = !panel.summary_open;
        }
    }
}

fn contents(
    ui: &mut egui::Ui,
    panel: &mut QcPanelState,
    task: &mut QcRunTask,
    opt_panel: &mut XtbPanelState,
    mol: &Molecule,
) -> bool {
    let mut state = presentation(ui.ctx());
    let mut close = false;
    ui.scope(|ui| {
        style(ui);
        close = header(ui, &mut state, panel, task.is_running());
        let code = code_status(&mut state, panel, opt_panel);
        let response = ui
            .add(
                egui::Button::new(
                    egui::RichText::new(code.label)
                        .small()
                        .color(if code.ready { ACCENT } else { WARNING }),
                )
                .small(),
            )
            .on_hover_text(&code.detail);
        if response.clicked_once() {
            state.configure = true;
        }
        electronic_state_row(ui, panel, mol, task.is_running());
        resources_row(ui, panel, task.is_running());
        if panel.has_input_file() && panel.advanced {
            ui.label(
                egui::RichText::new(
                    "Custom input active · settings and geometry do not update its text",
                )
                .small()
                .color(WARNING),
            );
        }
        let remaining = ui.available_rect_before_wrap();
        let footer_height = 112.0
            + if panel.has_input_file() { 38.0 } else { 0.0 }
            + if task.last.is_some() { 38.0 } else { 0.0 };
        let body_bottom = (remaining.bottom() - footer_height - 8.0).max(remaining.top() + 70.0);
        let body_rect =
            egui::Rect::from_min_max(remaining.min, egui::pos2(remaining.right(), body_bottom));
        let footer_rect = egui::Rect::from_min_max(
            egui::pos2(remaining.left(), body_bottom + 8.0),
            egui::pos2(
                remaining.right(),
                (body_bottom + footer_height + 8.0).max(remaining.bottom()),
            ),
        );
        let mut body = ui.new_child(
            egui::UiBuilder::new()
                .id_salt("modern_qc_body")
                .max_rect(body_rect),
        );
        egui::ScrollArea::vertical()
            .id_salt("modern_qc_scroll")
            .auto_shrink([false, false])
            .max_height(body_rect.height())
            .show(&mut body, |ui| {
                if state.configure {
                    configuration(ui, &mut state, panel, opt_panel, task.is_running());
                }
                job_choices(ui, panel, task.is_running());
                theory_card(ui, &mut state, panel, task.is_running());
                optional_settings(ui, panel, task.is_running());
                if panel.show_warnings {
                    if let Ok((_, warnings)) =
                        validate_electronic_state(mol, panel.charge, panel.multiplicity)
                    {
                        if !warnings.is_empty() {
                            ui.small(crate::qchem_interfaces::valence::BONDING_WARNING);
                        }
                    }
                }
            });
        // A path or Python check may have changed inside Configure this frame.
        let code = code_status(&mut state, panel, opt_panel);
        let mut footer_ui = ui.new_child(
            egui::UiBuilder::new()
                .id_salt("modern_qc_footer")
                .max_rect(footer_rect),
        );
        footer(&mut footer_ui, panel, task, opt_panel, mol, &code);
        ui.allocate_rect(remaining, egui::Sense::hover());
    });
    store(ui.ctx(), state);
    close
}

#[cfg(test)]
#[path = "modern_ui_tests.rs"]
mod tests;
