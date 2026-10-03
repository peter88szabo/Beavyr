//! Explicit pairing prevents a spectrum from opening an unrelated wavefunction.
use std::sync::{Arc, Weak};

use super::types::{Excitation, Spin as ExcitationSpin, TddftResult};
use crate::orbitals::{
    browser,
    density::Quantity,
    molden::{OrbitalData, Spin},
    OrbitalState,
};
use crate::ui_style::ResponseExt;
use bevy_egui::egui;

#[derive(Clone, Default)]
pub struct OrbitalLink {
    data: Weak<OrbitalData>,
    pub pending_molden: Option<std::path::PathBuf>,
    pub reveal: bool,
}

impl OrbitalLink {
    pub fn bind(&mut self, orbitals: &OrbitalState) {
        self.data = orbitals
            .data
            .as_ref()
            .map(Arc::downgrade)
            .unwrap_or_default();
    }

    pub fn matches(&self, orbitals: &OrbitalState) -> bool {
        self.data
            .upgrade()
            .zip(orbitals.data.as_ref())
            .is_some_and(|(linked, current)| Arc::ptr_eq(&linked, current))
    }

    pub fn controls(&mut self, ui: &mut egui::Ui, orbitals: &OrbitalState) {
        ui.horizontal_wrapped(|ui| {
            if ui
                .button("Matching Molden…")
                .on_hover_text("Load the orbital file from this spectrum's calculation")
                .clicked_once()
            {
                self.pending_molden = crate::recent_dir::pick_file(
                    crate::recent_dir::open().add_filter("Molden", &["molden", "input", "inp"]),
                );
            }
            if self.matches(orbitals) {
                ui.weak("Orbitals linked");
                if ui.small_button("Unlink").clicked_once() {
                    self.data = Weak::new();
                }
            } else if ui
                .add_enabled(
                    orbitals.data.is_some(),
                    egui::Button::new("Link current orbitals"),
                )
                .on_hover_text(
                    "Use only when the Surface tool holds the orbitals from this calculation",
                )
                .clicked_once()
            {
                self.bind(orbitals);
            }
        });
    }

    pub fn contribution(
        &mut self,
        ui: &mut egui::Ui,
        result: &TddftResult,
        e: &Excitation,
        orbitals: &mut OrbitalState,
    ) {
        ui.horizontal(|ui| {
            ui.add_space(24.0);
            for target in [false, true] {
                if target {
                    ui.monospace("→");
                }
                let (number, spin, label) = if target {
                    (e.to_orbital, e.to_spin, e.to_label.as_deref())
                } else {
                    (e.from_orbital, e.from_spin, e.from_label.as_deref())
                };
                let mapped = if !self.matches(orbitals) {
                    Err("Link the matching Molden file first".to_string())
                } else {
                    map_orbital(result, orbitals.data.as_ref().unwrap(), number, spin, label)
                };
                let tooltip = match &mapped {
                    Ok((spin, index)) => format!(
                        "Show {:?} orbital {} in Surface\n{}",
                        spin,
                        index + 1,
                        e.described()
                    ),
                    Err(error) => error.clone(),
                };
                let clicked = ui
                    .add_enabled(
                        mapped.is_ok(),
                        egui::Button::new(
                            egui::RichText::new(format!("{}{}", spin.suffix(), number)).monospace(),
                        )
                        .frame(false),
                    )
                    .on_hover_text(&tooltip)
                    .on_disabled_hover_text(&tooltip)
                    .clicked_once();
                if clicked {
                    if let Ok((spin, index)) = mapped {
                        if orbitals.quantity != Quantity::Orbital {
                            orbitals.quantity = Quantity::Orbital;
                            orbitals.isovalue = Quantity::Orbital.default_isovalue();
                        }
                        orbitals.spin = spin;
                        orbitals.show_surface = true;
                        orbitals.browser.reset_filters();
                        browser::select(orbitals, index);
                        self.reveal = true;
                    }
                }
            }
            ui.monospace(format!("{:5.1} %", e.percent()));
        });
    }
}

fn map_orbital(
    result: &TddftResult,
    data: &OrbitalData,
    number: u32,
    spin: ExcitationSpin,
    label: Option<&str>,
) -> Result<(Spin, usize), String> {
    let spin = match spin {
        ExcitationSpin::Beta => Spin::Beta,
        ExcitationSpin::Alpha => Spin::Alpha,
        ExcitationSpin::Unspecified if !result.unrestricted => Spin::Alpha,
        _ => return Err("This output does not identify the orbital's spin".into()),
    };
    let set = data.orbitals(spin);
    let index = if result.program.starts_with("ORCA") {
        // ORCA's output numbers MOs from zero; the Surface list starts at one.
        number as usize
    } else if result.program.starts_with("Psi4") {
        let n = number.checked_sub(1).ok_or("Invalid orbital number")? as usize;
        if let Some(irrep) = label.and_then(|s| s.strip_prefix("irrep ")) {
            set.iter()
                .enumerate()
                .filter(|(_, mo)| mo.symmetry.eq_ignore_ascii_case(irrep))
                .nth(n)
                .map(|(i, _)| i)
                .ok_or("Molden symmetry labels do not match this excitation")?
        } else {
            n
        }
    } else if result.program.starts_with("PySCF") || result.program.starts_with("Behemoth") {
        number.checked_sub(1).ok_or("Invalid orbital number")? as usize
    } else {
        return Err("Orbital numbering is not supported for this program".into());
    };
    if index >= set.len() {
        return Err("This orbital is absent from the linked Molden file".into());
    }
    Ok((spin, index))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::orbitals::molden::MolecularOrbital;
    #[test]
    fn numbering_spin_symmetry_and_stale_links() {
        let mut result = super::super::detect_and_parse("ORCA\nTD-DFT EXCITED STATES\nSTATE 1: E= 0.1 au 2.721 eV 21947 cm**-1\n 0a -> 1a : 0.9 (c=0.95)\n", std::path::Path::new("test.out")).unwrap();
        let mo = |sym: &str| MolecularOrbital {
            symmetry: sym.into(),
            energy: 0.0,
            occupation: 1.0,
            coefficients: vec![],
        };
        let data = OrbitalData {
            alpha: vec![mo("A1"), mo("B2"), mo("A1")],
            ..Default::default()
        };
        assert_eq!(
            map_orbital(&result, &data, 1, ExcitationSpin::Alpha, None)
                .unwrap()
                .1,
            1
        );
        result.program = "Behemoth".into();
        assert_eq!(
            map_orbital(&result, &data, 1, ExcitationSpin::Alpha, None)
                .unwrap()
                .1,
            0
        );
        assert!(map_orbital(&result, &data, 1, ExcitationSpin::Beta, None).is_err());
        result.program = "Psi4".into();
        assert_eq!(
            map_orbital(&result, &data, 2, ExcitationSpin::Alpha, Some("irrep A1"))
                .unwrap()
                .1,
            2
        );
        assert!(map_orbital(&result, &data, 2, ExcitationSpin::Alpha, Some("irrep B2")).is_err());
        let mut orbitals = OrbitalState::default();
        orbitals.adopt(data.clone(), None, None, None);
        let mut link = OrbitalLink::default();
        link.bind(&orbitals);
        assert!(link.matches(&orbitals));
        orbitals.adopt(data, None, None, None);
        assert!(!link.matches(&orbitals));
    }
}
