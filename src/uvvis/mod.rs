//! UV-Vis absorption spectra, read from a quantum-chemistry output or produced
//! by a run Beavyr started.
//!
//! Four readers fill one set of structures: ORCA, Psi4 and PySCF for real
//! TD-DFT, and Behemoth for its simplified sTDA and sTD-DFT. `types.rs` is
//! program-independent precisely so the panel never has to know which produced
//! what, and [`detect_and_parse`] is the one place that decides.

use bevy::prelude::*;

use crate::spectrum::broadening::BroadeningKind;
use types::TddftResult;

pub mod behemoth;
pub mod excited_state;
pub mod orca;
pub mod psi4;
pub mod pyscf;
pub mod run;
pub mod types;
pub mod ui;

/// Reads whichever program's output this is.
///
/// The file is identified by what it says about itself rather than by its
/// extension: every one of these programs writes `.out`, `.log` or `.dat`
/// depending only on what the user redirected to.
///
/// Each reader is tried only when its own program has been recognised, so a
/// file that is plainly ORCA's never reaches PySCF's reader and come back with
/// a confusing message about a heading it does not have. When nothing is
/// recognised, the error names every program that was considered -- which is
/// more use than "could not parse".
pub fn detect_and_parse(
    text: &str,
    source: &std::path::Path,
) -> Result<types::TddftResult, String> {
    if psi4::is_psi4_output(text) {
        return psi4::parse_psi4_tddft(text, source);
    }
    if pyscf::is_pyscf_output(text) {
        return pyscf::parse_pyscf_tddft(text, source);
    }
    // ORCA's reader is tried before Behemoth's because an ORCA file is
    // recognisable from its excited-state heading alone, where Behemoth's root
    // table needs more of the file read to be sure.
    if text.contains("EXCITED STATES") || text.contains("ORCA") {
        return orca::parse_orca_tddft(text, source);
    }
    if text.contains("sTDA") || text.contains("sTD-DFT") {
        return behemoth::parse_behemoth_spectrum(text, source, None);
    }
    // Nothing identified it. Try ORCA anyway, since that reader gives the most
    // specific message about what it expected, and fall back to naming the
    // four programs if it too finds nothing.
    orca::parse_orca_tddft(text, source).map_err(|_| {
        "This file carries no excited states Beavyr recognises. It reads TD-DFT output from \
         ORCA, Psi4 and PySCF, and sTDA/sTD-DFT output from Behemoth."
            .to_string()
    })
}

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
    /// What excited-state calculation to run, and the reference it runs on.
    /// The reference is a `MethodConfig` because that is what already holds a
    /// functional, a basis set, a memory budget and a thread count -- there is
    /// no second copy of those here.
    pub run_config: excited_state::ExcitedStateConfig,
    pub reference: crate::qchem_interfaces::method::MethodConfig,
    pub charge: i32,
    pub multiplicity: i32,
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
            run_config: excited_state::ExcitedStateConfig::default(),
            // The spectrum engine takes only global hybrids, so the reference
            // starts on one rather than on the shared default, which is not.
            reference: crate::qchem_interfaces::method::MethodConfig {
                behemoth: crate::qchem_interfaces::method::BehemothMethod::Dft,
                functional: excited_state::DEFAULT_SPECTRUM_FUNCTIONAL.to_string(),
                ..Default::default()
            },
            charge: 0,
            multiplicity: 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    /// Each program's output goes to its own reader. Getting this wrong is not
    /// a parse failure but a wrong answer: Psi4 prints its excitation energy in
    /// atomic units where PySCF prints electronvolts, so reading one with the
    /// other's reader is wrong by a factor of twenty-seven and plots a spectrum
    /// that looks entirely reasonable.
    #[test]
    fn each_program_reaches_its_own_reader() {
        let psi4 = concat!(
            "Psi4 1.11\n  Psi4 started on host\n",
            "Excited State    1 (1 B2):   0.27971 au   162.97 nm f =  0.0123\n",
        );
        let result = detect_and_parse(psi4, Path::new("a.out")).expect("Psi4 output");
        assert!(result.program.starts_with("Psi4"), "{}", result.program);
        // 0.27971 Hartree is 7.6 eV, not 0.28.
        assert!(result.states[0].energy_ev > 7.0, "{:?}", result.states[0]);

        let pyscf = concat!(
            "** Singlet excitation energies and oscillator strengths **\n",
            "Excited State   1:      7.61100 eV    162.97 nm  f=0.0123\n",
        );
        let result = detect_and_parse(pyscf, Path::new("b.log")).expect("PySCF output");
        assert_eq!(result.program, "PySCF");
        // Already electronvolts; it must not be multiplied again.
        assert!((result.states[0].energy_ev - 7.611).abs() < 1.0e-9);
    }

    /// A Psi4 run started from inside Beavyr reaches the panel by the same
    /// route as a file the user produced.
    ///
    /// Its `tdscf_excitations` prints through Psi4's own output channel, so
    /// what lands in `psi4.out` is the layout this reader was written against.
    /// Asserted here because the panel hands that log straight to
    /// `detect_and_parse`, and nothing else checks the two agree.
    #[test]
    fn a_psi4_run_started_here_reaches_the_panel() {
        // Shaped as the generated script's own output: Psi4's banner, then
        // what `tdscf_excitations` prints.
        let log = concat!(
            "    Psi4 1.11 release\n",
            "  Psi4 started on somehost\n",
            "                        Tamm-Dancoff approximation\n",
            "\n",
            "Excited State    1 (1 B2):   0.27971 au   162.97 nm f =  0.0123\n",
            "    4 (B2) -> 5 (A1)     0.68901\n",
            "Excited State    2 (1 A1):   0.34959 au   130.39 nm f =  0.0000\n",
            "BEAVYR_ENERGY -76.321074727877\n",
        );
        let result = detect_and_parse(log, Path::new("Psi4 run")).expect("a Psi4 spectrum");
        assert!(result.program.starts_with("Psi4"), "{}", result.program);
        assert_eq!(result.states.len(), 2);
        assert_eq!(result.method, "TD-DFT/TDA");
        // Atomic units in, electronvolts out.
        assert!(result.states[0].energy_ev > 7.0, "{:?}", result.states[0]);
        assert_eq!(result.states[0].oscillator_strength, Some(0.0123));
    }

    /// A file none of the readers recognises says so, and names what Beavyr
    /// does read -- more use than refusing without explanation.
    #[test]
    fn an_unrecognised_file_names_what_is_supported() {
        let err = detect_and_parse("just some text\n", Path::new("x.txt"))
            .expect_err("nothing to read here");
        for program in ["ORCA", "Psi4", "PySCF", "Behemoth"] {
            assert!(err.contains(program), "{program} is not mentioned: {err}");
        }
    }

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
