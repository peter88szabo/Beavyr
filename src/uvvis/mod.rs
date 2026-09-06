//! UV-Vis absorption spectra from TD-DFT output files.
//!
//! The tool reads an existing quantum-chemistry output rather than running the
//! program itself -- unlike the xTB optimizer and Hessian tools, a TD-DFT job
//! is something the user has already queued and run elsewhere. ORCA is the
//! only reader for now; `types.rs` is deliberately program-independent so a
//! second one slots in beside `orca.rs` without the UI noticing.

use bevy::prelude::*;

use crate::spectrum::broadening::BroadeningKind;
use types::TddftResult;

pub mod orca;
pub mod types;
pub mod ui;

/// Which axis the absorption spectrum is drawn against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpectrumUnit {
    Nanometres,
    ElectronVolts,
}

impl SpectrumUnit {
    pub fn label(self) -> &'static str {
        match self {
            SpectrumUnit::Nanometres => "nm",
            SpectrumUnit::ElectronVolts => "eV",
        }
    }

    pub fn axis_title(self) -> &'static str {
        match self {
            SpectrumUnit::Nanometres => "Wavelength (nm)",
            SpectrumUnit::ElectronVolts => "Energy (eV)",
        }
    }

    /// How much blank axis to leave beyond the outermost peak, in this unit.
    pub fn axis_margin(self) -> f64 {
        match self {
            SpectrumUnit::Nanometres => 50.0,
            SpectrumUnit::ElectronVolts => 0.5,
        }
    }

    /// The grid an exported `.dat` is sampled on: 0.1 nm as specified, and a
    /// comparably fine 0.001 eV on the energy axis.
    pub fn export_resolution(self) -> f64 {
        match self {
            SpectrumUnit::Nanometres => 0.1,
            SpectrumUnit::ElectronVolts => 0.001,
        }
    }

    /// Decimal places for the axis tick labels.
    pub fn tick_decimals(self) -> usize {
        match self {
            SpectrumUnit::Nanometres => 0,
            SpectrumUnit::ElectronVolts => 1,
        }
    }

    /// Picks this unit's line width out of the pair the state keeps one of per
    /// axis. Kept here so the two call sites -- reading the width to broaden
    /// with, and handing it to the drag control -- cannot disagree about which
    /// field belongs to which axis.
    pub fn width_of(self, width_nm: f64, width_ev: f64) -> f64 {
        match self {
            SpectrumUnit::Nanometres => width_nm,
            SpectrumUnit::ElectronVolts => width_ev,
        }
    }

    pub fn width_of_mut<'a>(self, width_nm: &'a mut f64, width_ev: &'a mut f64) -> &'a mut f64 {
        match self {
            SpectrumUnit::Nanometres => width_nm,
            SpectrumUnit::ElectronVolts => width_ev,
        }
    }
}

/// Everything the UV-Vis tool remembers between frames.
#[derive(Resource)]
pub struct UvVisState {
    /// The parsed file, once one has been loaded.
    pub result: Option<TddftResult>,
    /// Why the last load attempt failed, shown until the next one succeeds.
    pub load_error: Option<String>,
    /// The root whose contributing excitations are expanded in the table.
    pub expanded_root: Option<usize>,
    pub spectrum_window_open: bool,
    pub unit: SpectrumUnit,
    pub broadening: BroadeningKind,
    /// Line width in nm, used while the nm axis is showing.
    pub width_nm: f64,
    /// Line width in eV, kept separately so switching axes back and forth
    /// does not mangle a width the user set in the other unit.
    pub width_ev: f64,
}

impl Default for UvVisState {
    fn default() -> Self {
        Self {
            result: None,
            load_error: None,
            expanded_root: None,
            spectrum_window_open: false,
            unit: SpectrumUnit::Nanometres,
            broadening: BroadeningKind::Gaussian,
            // Conventional starting points for a simulated UV-Vis band: broad
            // enough that a dense root list reads as bands rather than grass.
            width_nm: 20.0,
            width_ev: 0.3,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_unit_describes_its_own_axis() {
        assert_eq!(SpectrumUnit::Nanometres.label(), "nm");
        assert!(SpectrumUnit::Nanometres.axis_title().contains("nm"));
        assert!(SpectrumUnit::ElectronVolts.axis_title().contains("eV"));
    }

    /// The user asked for 0.1 nm exports specifically.
    #[test]
    fn the_nanometre_export_grid_is_a_tenth_of_a_nanometre() {
        assert_eq!(SpectrumUnit::Nanometres.export_resolution(), 0.1);
    }

    /// Switching axes must not carry a width of 20 nm over as 20 eV, which
    /// would flatten the whole spectrum into one featureless hump.
    #[test]
    fn each_axis_keeps_its_own_width() {
        let mut state = UvVisState::default();
        let width = |s: &UvVisState| s.unit.width_of(s.width_nm, s.width_ev);
        assert_eq!(width(&state), 20.0);
        state.unit = SpectrumUnit::ElectronVolts;
        assert_eq!(width(&state), 0.3);
        *state
            .unit
            .width_of_mut(&mut state.width_nm, &mut state.width_ev) = 0.8;
        state.unit = SpectrumUnit::Nanometres;
        assert_eq!(width(&state), 20.0, "the nm width is untouched");
        state.unit = SpectrumUnit::ElectronVolts;
        assert_eq!(width(&state), 0.8);
    }

    #[test]
    fn nothing_is_loaded_or_open_on_first_launch() {
        let state = UvVisState::default();
        assert!(state.result.is_none());
        assert!(state.load_error.is_none());
        assert!(!state.spectrum_window_open);
    }
}
