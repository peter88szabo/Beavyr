//! Photon-energy axis and optical frequency factors for absorption/emission.
//! Spectra are relative; no fluorescence lifetime is inferred from an arbitrarily
//! broadened and finite spectral window.
use super::{
    harmonic_model::HarmonicModel,
    thermal_spectrum::{self, HarmonicSpectrum, SpectrumGrid},
};
use std::sync::atomic::{AtomicBool, AtomicUsize};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Process {
    Absorption,
    Fluorescence,
    Raman,
}
impl Process {
    pub fn label(self) -> &'static str {
        match self {
            Self::Absorption => "Absorption",
            Self::Fluorescence => "Fluorescence",
            Self::Raman => "Resonance Raman",
        }
    }
    pub fn orca(self) -> &'static str {
        match self {
            Self::Absorption => "ABS",
            Self::Fluorescence => "FLUOR",
            Self::Raman => "RR",
        }
    }
}
pub fn calculate(
    model: HarmonicModel,
    process: Process,
    grid: &SpectrumGrid,
    frequency_factor: bool,
    cancel: &AtomicBool,
    progress: &AtomicUsize,
) -> Result<HarmonicSpectrum, String> {
    calculate_with_dipole(
        model,
        process,
        grid,
        frequency_factor,
        None,
        cancel,
        progress,
    )
}
pub fn calculate_with_dipole(
    model: HarmonicModel,
    process: Process,
    grid: &SpectrumGrid,
    frequency_factor: bool,
    dipole: Option<&super::transition_properties::Dipole>,
    cancel: &AtomicBool,
    progress: &AtomicUsize,
) -> Result<HarmonicSpectrum, String> {
    let reversed_dipole = if process == Process::Fluorescence {
        dipole.map(|d| d.reversed(&model))
    } else {
        None
    };
    let dipole = reversed_dipole.as_ref().or(dipole);
    let model = if process == Process::Fluorescence {
        model.reversed()?
    } else {
        model
    };
    let mut spectrum =
        thermal_spectrum::calculate_with_dipole(model, grid, dipole, cancel, progress)?;
    if process == Process::Fluorescence {
        for p in &mut spectrum.profile {
            p.0 = -p.0;
        }
        spectrum.zero_zero_ev = -spectrum.zero_zero_ev;
        spectrum.profile.reverse();
    }
    spectrum.profile.retain(|p| p.0 > 0.);
    if frequency_factor {
        let power = if process == Process::Fluorescence {
            3
        } else {
            1
        };
        for p in &mut spectrum.profile {
            p.1 *= p.0.powi(power);
        }
        spectrum.notes.push(format!("Relative {} spectrum with the photon-energy factor E^{power}; temperature refers to the initial electronic state.",process.label().to_lowercase()));
    } else {
        spectrum
            .notes
            .push("Franck–Condon density without the optical photon-energy prefactor.".into());
    }
    let area = spectrum
        .profile
        .iter()
        .map(|p| p.1 * spectrum.spacing_ev)
        .sum::<f64>();
    if !area.is_finite() || area <= 0. {
        return Err("No finite positive-energy optical intensity remains".into());
    }
    for p in &mut spectrum.profile {
        p.1 /= area;
    }
    Ok(spectrum)
}
