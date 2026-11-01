//! Vibronic analysis adapted from Behemoth, with a separate UV-Vis analysis window.
mod analysis;
mod band_shape;
mod constants;
mod duschinsky;
mod franck_condon;
mod input;
mod plot;
mod state;
pub mod ui;
pub use band_shape::{vibronic_band, BandKind, VibronicBand};
pub(crate) use constants::*;
pub use duschinsky::{duschinsky, Duschinsky};
#[cfg(test)]
use franck_condon::ModeActivity;
pub use franck_condon::{franck_condon_activity, FranckCondonAnalysis};
pub use state::VibronicState;
#[cfg(test)]
mod band_shape_tests;
#[cfg(test)]
mod legacy;
#[cfg(test)]
mod tests;
#[cfg(test)]
mod workflow_tests;

mod complex_math;
mod gaussian_correlation;
mod harmonic_model;
#[cfg(test)]
mod harmonic_tests;
mod thermal_spectrum;

mod orca_esd;
mod process;
mod spectral_controls;
mod spectral_job;

mod optical_spectrum;

mod raman_spectrum;

mod transition_properties;

mod property_job;

mod mode_animation;
mod spectrum_view;

#[cfg(test)]
mod extensions_tests;
