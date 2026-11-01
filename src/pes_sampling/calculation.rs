//! Open the sampler using the Hessian and geometry of a completed calculation.
use super::{input, state::*};
use crate::qchem_interfaces::{job::JobType, program::QcProgram};
use crate::qchem_panel::run::QcRunOutput;
use crate::ui_style::ResponseExt;
use bevy_egui::egui;
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
                if !fragment.excluded.is_empty() {
                    crate::runtime_warning::push(ui.ctx(), "PES sampling: excluded modes", &format!(
                        "As in Smite, imaginary modes are excluded from normal-mode sampling. Excluded frequencies (cm^-1): {:?}", fragment.excluded));
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

fn read(path: &Path) -> Result<String, String> {
    fs::read_to_string(path).map_err(|e| format!("Cannot read {}: {e}", path.display()))
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
    if read(path)?.lines().any(|line| line.trim() == "$atoms") {
        return input::load_hessian(path, None);
    }
    let dir = path
        .parent()
        .ok_or("The Hessian has no calculation directory.")?;
    // Use backend files from this run, never the current viewer or frequency
    // panel: either may already contain a different structure. In particular,
    // xTB's Hessian uses its reoriented g98 geometry, not the input orientation.
    let geometry = match output.program {
        QcProgram::Xtb => crate::qchem_interfaces::hessian_file::parse_g98_geometry(
            &read(&dir.join("hess/g98.out"))?,
        )?,
        QcProgram::Behemoth => input::geometry(&read(&dir.join("hess/input.xyz"))?)?,
        QcProgram::Dreiding => input::geometry(&read(&dir.join("hessian_geometry.xyz"))?)?,
        QcProgram::PySCF | QcProgram::Psi4 | QcProgram::Psi4Py =>
            input::geometry(&read(&dir.join("final.xyz"))?)?,
        QcProgram::SparrowPy => input::geometry(&read(&dir.join("geometry.xyz"))?)?,
        _ => return Err("The Hessian does not contain its matching geometry. Load the Hessian and its equilibrium XYZ in PES Sampling.".into()),
    };
    input::load_hessian(path, Some(&geometry))
}

#[cfg(test)]
#[path = "calculation_tests.rs"]
mod tests;
