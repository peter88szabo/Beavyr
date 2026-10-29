//! Presentation units for the transition table, independent of the spectrum axes.

use super::super::types::ExcitedState;
use crate::ui_style::ResponseExt;
use bevy_egui::egui;

#[derive(Clone, Copy, Default, PartialEq, Eq)]
pub(super) enum EnergyUnit {
    #[default]
    ElectronVolts,
    Wavenumbers,
}

impl EnergyUnit {
    pub fn header(self) -> &'static str {
        match self {
            Self::ElectronVolts => "E (eV)",
            Self::Wavenumbers => "E(cm-1)",
        }
    }

    pub fn value(self, state: &ExcitedState) -> String {
        match self {
            Self::ElectronVolts => format!("{:>8.3}", state.energy_ev),
            Self::Wavenumbers => format!("{:>8.1}", state.energy_cm1),
        }
    }
}

pub(super) fn selected(ctx: &egui::Context) -> EnergyUnit {
    ctx.data_mut(|data| data.get_temp::<EnergyUnit>(egui::Id::new("uvvis_table_energy_unit")))
        .unwrap_or_default()
}

pub(super) fn switcher(ui: &mut egui::Ui) {
    let mut unit = selected(ui.ctx());
    ui.horizontal(|ui| {
        for (label, value) in [
            ("eV", EnergyUnit::ElectronVolts),
            ("cm-1", EnergyUnit::Wavenumbers),
        ] {
            if ui
                .add(egui::Button::new(label).selected(unit == value))
                .clicked_once()
            {
                unit = value;
            }
        }
    });
    ui.ctx()
        .data_mut(|data| data.insert_temp(egui::Id::new("uvvis_table_energy_unit"), unit));
}
