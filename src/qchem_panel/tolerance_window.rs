//! Transactional editing of optimizer stopping criteria.
use super::{job, QcPanelState};
use crate::qchem_interfaces::{
    job::JobType,
    method::{tolerances::Criteria, MethodConfig},
    program::QcProgram,
};
use crate::ui_style::ResponseExt;
use bevy_egui::egui;

#[derive(Clone, Copy, Debug, PartialEq)]
enum Target {
    Minimum,
    Ts,
    Mecp,
    OrcaMecp,
    Native,
    Xtb,
    Scan,
}
#[derive(Clone)]
struct Editor {
    target: Target,
    program: QcProgram,
    job: JobType,
    values: Criteria,
    gap: f64,
    native_defaults: bool,
    xtb_level: String,
}
fn id() -> egui::Id {
    egui::Id::new("optimization_tolerance_editor")
}
fn target(panel: &QcPanelState) -> Option<Target> {
    if panel.job == JobType::Mecp {
        return Some(
            if panel.program == QcProgram::Orca
                && panel.config.mecp.driver == job::mecp::Driver::Orca
            {
                Target::OrcaMecp
            } else {
                Target::Mecp
            },
        );
    }
    if panel.job == JobType::Scan && panel.config.scan.mode == job::scan::Mode::Relaxed {
        return Some(Target::Scan);
    }
    if !job::search::active(panel.job) {
        return None;
    }
    if job::search::local(panel.program, panel.job, &panel.config.search)
        || (panel.program == QcProgram::Behemoth && panel.job == JobType::TransitionState)
    {
        return Some(if panel.job == JobType::TransitionState {
            Target::Ts
        } else {
            Target::Minimum
        });
    }
    Some(if panel.program == QcProgram::Xtb {
        Target::Xtb
    } else {
        Target::Native
    })
}
impl Editor {
    fn new(panel: &QcPanelState, target: Target, defaults: bool) -> Self {
        let default_config = MethodConfig::default();
        let c = if defaults {
            &default_config
        } else {
            &panel.config
        };
        let values = match target {
            Target::Minimum => {
                let t = &c.search.minimum;
                Criteria {
                    energy: t.energy_tol,
                    max_step: t.max_step,
                    rms_step: t.rms_step,
                    max_gradient: t.max_gradient,
                    rms_gradient: t.rms_gradient,
                }
            }
            Target::Ts => {
                let t = &c.search.ts;
                Criteria {
                    energy: t.energy_tol,
                    max_step: t.max_step,
                    rms_step: t.rms_step,
                    max_gradient: t.max_gradient,
                    rms_gradient: t.rms_gradient,
                }
            }
            Target::Mecp => {
                let t = &c.mecp;
                Criteria {
                    energy: t.gap_tolerance,
                    max_step: t.max_step,
                    rms_step: t.rms_step,
                    max_gradient: t.max_gradient,
                    rms_gradient: t.rms_gradient,
                }
            }
            Target::Scan => c.scan.tolerances.clone(),
            Target::OrcaMecp => c
                .mecp
                .native_tolerances
                .clone()
                .unwrap_or_else(|| Criteria::defaults(QcProgram::Orca)),
            _ => c
                .optimization
                .get(panel.program)
                .cloned()
                .unwrap_or_else(|| Criteria::defaults(panel.program)),
        };
        Self {
            target,
            program: panel.program,
            job: panel.job,
            values,
            gap: c.mecp.gap_tolerance,
            native_defaults: if target == Target::OrcaMecp {
                c.mecp.native_tolerances.is_none()
            } else {
                c.optimization.get(panel.program).is_none()
            },
            xtb_level: c
                .optimization
                .xtb_level
                .clone()
                .unwrap_or_else(|| "normal".into()),
        }
    }
    fn save(&self, c: &mut MethodConfig) -> Result<(), String> {
        self.values.validate()?;
        if !self.gap.is_finite() || self.gap <= 0. {
            return Err("The energy-gap tolerance must be positive and finite.".into());
        }
        let v = &self.values;
        match self.target {
            Target::Minimum => {
                let t = &mut c.search.minimum;
                t.energy_tol = v.energy;
                t.max_step = v.max_step;
                t.rms_step = v.rms_step;
                t.max_gradient = v.max_gradient;
                t.rms_gradient = v.rms_gradient;
            }
            Target::Ts => {
                let t = &mut c.search.ts;
                t.energy_tol = v.energy;
                t.max_step = v.max_step;
                t.rms_step = v.rms_step;
                t.max_gradient = v.max_gradient;
                t.rms_gradient = v.rms_gradient;
            }
            Target::Mecp => {
                let t = &mut c.mecp;
                t.gap_tolerance = v.energy;
                t.max_step = v.max_step;
                t.rms_step = v.rms_step;
                t.max_gradient = v.max_gradient;
                t.rms_gradient = v.rms_gradient;
            }
            Target::Scan => c.scan.tolerances = v.clone(),
            Target::OrcaMecp => {
                c.mecp.gap_tolerance = self.gap;
                c.mecp.native_tolerances = (!self.native_defaults).then(|| v.clone());
            }
            Target::Native => c
                .optimization
                .set(self.program, (!self.native_defaults).then(|| v.clone())),
            Target::Xtb => {
                c.optimization.xtb_level =
                    (self.xtb_level != "normal").then(|| self.xtb_level.clone())
            }
        }
        Ok(())
    }
}
pub(super) fn button(ui: &mut egui::Ui, panel: &QcPanelState, running: bool) {
    if let Some(target) = target(panel) {
        if ui
            .add_enabled(!running, egui::Button::new("Tolerances"))
            .clicked_once()
        {
            ui.ctx()
                .data_mut(|d| d.insert_temp(id(), Editor::new(panel, target, false)));
        }
    }
}
fn number(ui: &mut egui::Ui, label: &str, value: &mut f64) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(
            egui::DragValue::new(value)
                .speed(1e-6)
                .max_decimals(12)
                .range(1e-14..=1.),
        );
    });
}
pub(crate) fn window(
    ctx: &egui::Context,
    panel: &mut QcPanelState,
    running: bool,
) -> Option<egui::Rect> {
    let Some(mut edit) = ctx.data_mut(|d| d.get_temp::<Editor>(id())) else {
        return None;
    };
    // A draft cannot be applied to a different optimizer after a selection change.
    if edit.program != panel.program || edit.job != panel.job || Some(edit.target) != target(panel)
    {
        ctx.data_mut(|d| d.remove::<Editor>(id()));
        return None;
    }
    let mut open = true;
    let mut close = false;
    let result = egui::Window::new("Optimization tolerances")
        .id(id().with("window"))
        .open(&mut open)
        .resizable(true)
        .default_size([460., 330.])
        .min_size([300., 200.])
        .show(ctx, |ui| {
            egui::ScrollArea::both()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.add_enabled_ui(!running, |ui| {
                        if edit.target == Target::Xtb {
                            egui::ComboBox::from_id_salt("xtb_tolerance_level")
                                .selected_text(&edit.xtb_level)
                                .show_ui(ui, |ui| {
                                    for level in [
                                        "crude",
                                        "sloppy",
                                        "loose",
                                        "normal",
                                        "tight",
                                        "verytight",
                                        "extreme",
                                    ] {
                                        ui.selectable_value(
                                            &mut edit.xtb_level,
                                            level.into(),
                                            level,
                                        );
                                    }
                                });
                        } else {
                            let native = matches!(edit.target, Target::Native | Target::OrcaMecp);
                            if native {
                                ui.checkbox(&mut edit.native_defaults, "Use backend defaults");
                            }
                            ui.add_enabled_ui(!native || !edit.native_defaults, |ui| {
                                let crossing = edit.target == Target::Mecp;
                                let psi4_internal = edit.target == Target::Native
                                    && matches!(edit.program, QcProgram::Psi4 | QcProgram::Psi4Py);
                                let v = &mut edit.values;
                                number(
                                    ui,
                                    if crossing {
                                        "Energy gap (Eh)"
                                    } else {
                                        "Energy change (Eh)"
                                    },
                                    &mut v.energy,
                                );
                                number(
                                    ui,
                                    if crossing {
                                        "Maximum displacement (Å)"
                                    } else if psi4_internal {
                                        "Maximum displacement (a.u.)"
                                    } else {
                                        "Maximum displacement (bohr)"
                                    },
                                    &mut v.max_step,
                                );
                                number(
                                    ui,
                                    if crossing {
                                        "RMS displacement (Å)"
                                    } else if psi4_internal {
                                        "RMS displacement (a.u.)"
                                    } else {
                                        "RMS displacement (bohr)"
                                    },
                                    &mut v.rms_step,
                                );
                                number(
                                    ui,
                                    if crossing {
                                        "Maximum gradient (Eh/Å)"
                                    } else if psi4_internal {
                                        "Maximum gradient (a.u.)"
                                    } else {
                                        "Maximum gradient (Eh/bohr)"
                                    },
                                    &mut v.max_gradient,
                                );
                                number(
                                    ui,
                                    if crossing {
                                        "RMS gradient (Eh/Å)"
                                    } else if psi4_internal {
                                        "RMS gradient (a.u.)"
                                    } else {
                                        "RMS gradient (Eh/bohr)"
                                    },
                                    &mut v.rms_gradient,
                                );
                            });
                            if edit.target == Target::OrcaMecp {
                                number(ui, "Final gap limit (Eh)", &mut edit.gap);
                            }
                        }
                        ui.separator();
                        ui.horizontal(|ui| {
                            if ui.button("Save").clicked_once() {
                                match edit.save(&mut panel.config) {
                                    Ok(()) => close = true,
                                    Err(e) => crate::runtime_warning::push(
                                        ctx,
                                        "Invalid optimization tolerances",
                                        &e,
                                    ),
                                }
                            }
                            if ui.button("Reset defaults").clicked_once() {
                                edit = Editor::new(panel, edit.target, true);
                            }
                            if ui.button("Cancel").clicked_once() {
                                close = true;
                            }
                        });
                    });
                });
        });
    ctx.data_mut(|d| {
        if open && !close {
            d.insert_temp(id(), edit);
        } else {
            d.remove::<Editor>(id());
        }
    });
    result.map(|r| r.response.rect)
}

#[cfg(test)]
#[path = "optimization_tolerance_tests.rs"]
mod tests;
