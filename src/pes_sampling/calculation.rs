//! Open the sampler using the Hessian and geometry of a completed calculation.
use super::{input, state::*};
use crate::qchem_interfaces::job::JobType;
#[cfg(test)]
use crate::qchem_interfaces::program::QcProgram;
use crate::qchem_panel::run::QcRunOutput;
use crate::ui_style::ResponseExt;
use bevy_egui::egui;
#[cfg(test)]
use std::{fs, path::Path};

fn supports_sampling(job: JobType) -> bool {
    matches!(
        job,
        JobType::Hessian | JobType::Frequencies | JobType::OptimizeThenFrequencies
    )
}

pub fn results_button(ui: &mut egui::Ui, output: &QcRunOutput, state: &mut PesSampling) {
    if !supports_sampling(output.job) {
        return;
    }
    if ui
        .add_enabled(state.task.is_none(), egui::Button::new("Sample PES"))
        .on_disabled_hover_text("Wait for the current PES sampling run to finish.")
        .clicked_once()
    {
        state.open = true;
        match load_calculation(output) {
            Ok(fragment) => {
                if fragment.imaginary_cm1.len() > 1 {
                    crate::runtime_warning::push(ui.ctx(), "PES sampling: excluded modes", &format!(
                        "Multiple imaginary modes are excluded from PES sampling; this is not a first-order saddle. Frequencies (cm^-1): {:?}", fragment.imaginary_cm1));
                }
                adopt(state, fragment);
            }
            Err(message) => {
                state.message = Some(message.clone());
                state.is_error = true;
                crate::runtime_warning::push(
                    ui.ctx(),
                    "Cannot load calculation for PES sampling",
                    &message,
                );
            }
        }
    }
}

fn adopt(state: &mut PesSampling, fragment: input::Fragment) {
    state.settings.workflow = Workflow::Unimolecular;
    let partner = &mut state.settings.partners[0];
    partner.source = Source::Hessian;
    partner.rotors.clear();
    partner.rigid = false;
    partner.geometry = Some((
        fragment.atoms.clone(),
        fragment
            .equilibrium
            .iter()
            .map(|x| x * input::BOHR)
            .collect(),
    ));
    partner.overrides = vec![partner.excitation.clone(); fragment.frequencies.len()];
    partner.use_overrides = false;
    partner.fragment = Some(fragment);
    state.message = None;
    state.is_error = false;
}

fn load_calculation(output: &QcRunOutput) -> Result<input::Fragment, String> {
    let path = output
        .hessian_path
        .as_deref()
        .or_else(|| {
            // A matrix-only Hessian job offers its file for saving without opening
            // Frequency Analysis, so hessian_path is intentionally empty for it.
            output
                .saveable
                .iter()
                .find(|(label, _)| label.starts_with("Hessian"))
                .map(|(_, path)| path.as_path())
        })
        .ok_or("This calculation did not produce a Hessian file.")?;
    let h = crate::qchem_panel::run::hessian_output::load(output, path)?;
    input::from_hessian(
        h.atoms,
        h.masses_amu,
        h.coords_bohr,
        h.hessian,
        path.file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned(),
    )
}

#[cfg(test)]
#[path = "calculation_tests.rs"]
mod tests;
