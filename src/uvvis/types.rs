//! The program-independent description of a TD-DFT excited-state calculation.
//!
//! Nothing here knows about ORCA. A parser for another program fills the same
//! structures, and the UI only ever sees these.

use std::path::PathBuf;

/// Planck's constant times the speed of light, in eV·nm -- the conversion
/// between a transition energy and its wavelength.
pub const EV_NM: f64 = 1239.841_984_3;

/// CODATA 2018 hartree-to-electronvolt conversion.
pub const HARTREE_TO_EV: f64 = 27.211_386_245_988;

/// Which spin manifold an orbital belongs to. Restricted references print no
/// spin label at all, which is `Unspecified` rather than a guess at alpha.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Spin {
    Alpha,
    Beta,
    Unspecified,
}

impl Spin {
    /// The letter the output file uses, or nothing for a restricted run.
    pub fn suffix(self) -> &'static str {
        match self {
            Spin::Alpha => "a",
            Spin::Beta => "b",
            Spin::Unspecified => "",
        }
    }

    pub fn from_letter(letter: char) -> Option<Spin> {
        match letter {
            'a' | 'A' => Some(Spin::Alpha),
            'b' | 'B' => Some(Spin::Beta),
            _ => None,
        }
    }
}

/// One orbital-to-orbital contribution to an excited state.
#[derive(Debug, Clone, PartialEq)]
pub struct Excitation {
    pub from_orbital: u32,
    pub from_spin: Spin,
    pub to_orbital: u32,
    pub to_spin: Spin,
    /// The weight |c|², as printed. Only contributions above the program's own
    /// print threshold appear, so these do not sum to one.
    pub weight: f64,
    /// The signed CI coefficient, when the output prints one.
    pub coefficient: Option<f64>,
}

impl Excitation {
    /// The orbital pair as it reads in the output, e.g. `89b → 90b`.
    pub fn label(&self) -> String {
        format!(
            "{}{} \u{2192} {}{}",
            self.from_orbital,
            self.from_spin.suffix(),
            self.to_orbital,
            self.to_spin.suffix()
        )
    }
}

/// One root of the TD-DFT problem.
#[derive(Debug, Clone, PartialEq)]
pub struct ExcitedState {
    /// The root number as the program labels it, counting from 1.
    pub root: usize,
    pub energy_au: f64,
    pub energy_ev: f64,
    pub energy_cm1: f64,
    /// Taken from the absorption table when it lists one, else derived from
    /// the transition energy.
    pub wavelength_nm: f64,
    /// The spin-squared expectation value; only unrestricted references have
    /// one to report.
    pub s2: Option<f64>,
    /// The multiplicity the program estimated from the rounded `<S²>`.
    pub multiplicity: Option<u32>,
    /// `None` when the run produced no absorption table -- distinct from a
    /// genuinely dark state, whose oscillator strength is `Some(0.0)`.
    pub oscillator_strength: Option<f64>,
    pub d2_au2: Option<f64>,
    pub transition_dipole_au: Option<[f64; 3]>,
    /// In the order printed, which is by ascending source orbital.
    pub excitations: Vec<Excitation>,
}

impl ExcitedState {
    /// Contributions from largest weight down, so the dominant configuration
    /// reads first -- the order a spectroscopist wants, as opposed to the
    /// orbital-index order the file happens to print.
    pub fn excitations_by_weight(&self) -> Vec<&Excitation> {
        let mut sorted: Vec<&Excitation> = self.excitations.iter().collect();
        sorted.sort_by(|a, b| b.weight.total_cmp(&a.weight));
        sorted
    }

    /// The single largest contribution, the one usually quoted as "the"
    /// character of the state.
    pub fn dominant_excitation(&self) -> Option<&Excitation> {
        self.excitations
            .iter()
            .max_by(|a, b| a.weight.total_cmp(&b.weight))
    }
}

/// Everything one TD-DFT output file had to say.
#[derive(Debug, Clone, PartialEq)]
pub struct TddftResult {
    pub source: PathBuf,
    /// e.g. `ORCA 6.1.0`.
    pub program: String,
    /// e.g. `TD-DFT/TDA`.
    pub method: String,
    pub unrestricted: bool,
    /// The reference determinant's own `<S²>`, for judging spin contamination
    /// of the ground state against the excited ones.
    pub ground_state_s2: Option<f64>,
    pub states: Vec<ExcitedState>,
    /// Anything the parser wants the user to know about what it could not
    /// read -- a missing absorption table, extra state blocks it skipped.
    pub notes: Vec<String>,
}

impl TddftResult {
    /// Whether an absorption table was found. Without one there are still
    /// states to list, but nothing to plot.
    pub fn has_intensities(&self) -> bool {
        self.states
            .iter()
            .any(|s| s.oscillator_strength.is_some())
    }

    /// `(wavelength_nm, oscillator_strength)` for every state that has an
    /// intensity. Dark states stay in the list -- they contribute a zero-height
    /// peak, which is correct, and keeps peak indices lined up with the table.
    pub fn peaks_nm(&self) -> Vec<(f64, f64)> {
        self.states
            .iter()
            .filter_map(|s| s.oscillator_strength.map(|f| (s.wavelength_nm, f)))
            .collect()
    }

    /// The same peaks against an energy axis.
    pub fn peaks_ev(&self) -> Vec<(f64, f64)> {
        self.states
            .iter()
            .filter_map(|s| s.oscillator_strength.map(|f| (s.energy_ev, f)))
            .collect()
    }

    /// The states behind `peaks_nm` / `peaks_ev`, in the same order, so a plot
    /// hover can name the root it is pointing at.
    pub fn states_with_intensity(&self) -> Vec<&ExcitedState> {
        self.states
            .iter()
            .filter(|s| s.oscillator_strength.is_some())
            .collect()
    }
}

/// A transition energy in eV as a wavelength in nm. Zero and negative energies
/// have no wavelength; they are reported as zero rather than an infinity that
/// would blow up an axis range.
pub fn ev_to_nm(energy_ev: f64) -> f64 {
    if energy_ev > 0.0 {
        EV_NM / energy_ev
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn excitation(from: u32, to: u32, weight: f64) -> Excitation {
        Excitation {
            from_orbital: from,
            from_spin: Spin::Beta,
            to_orbital: to,
            to_spin: Spin::Beta,
            weight,
            coefficient: None,
        }
    }

    #[test]
    fn an_excitation_reads_the_way_the_output_prints_it() {
        assert_eq!(excitation(89, 90, 0.99).label(), "89b \u{2192} 90b");
    }

    #[test]
    fn a_restricted_excitation_carries_no_spin_letter() {
        let mut e = excitation(45, 47, 0.5);
        e.from_spin = Spin::Unspecified;
        e.to_spin = Spin::Unspecified;
        assert_eq!(e.label(), "45 \u{2192} 47");
    }

    fn state(root: usize, ev: f64, fosc: Option<f64>) -> ExcitedState {
        ExcitedState {
            root,
            energy_au: ev / 27.211,
            energy_ev: ev,
            energy_cm1: ev * 8065.54,
            wavelength_nm: ev_to_nm(ev),
            s2: None,
            multiplicity: None,
            oscillator_strength: fosc,
            d2_au2: None,
            transition_dipole_au: None,
            excitations: vec![
                excitation(87, 92, 0.358),
                excitation(88, 91, 0.550),
                excitation(89, 91, 0.059),
            ],
        }
    }

    #[test]
    fn contributions_are_ranked_by_weight_not_orbital_index() {
        let s = state(8, 2.29, Some(0.045));
        let ranked: Vec<u32> = s.excitations_by_weight().iter().map(|e| e.from_orbital).collect();
        assert_eq!(ranked, vec![88, 87, 89]);
        assert_eq!(s.dominant_excitation().unwrap().from_orbital, 88);
    }

    #[test]
    fn a_state_with_no_listed_contributions_has_no_dominant_one() {
        let mut s = state(1, 1.0, None);
        s.excitations.clear();
        assert!(s.dominant_excitation().is_none());
        assert!(s.excitations_by_weight().is_empty());
    }

    fn result(states: Vec<ExcitedState>) -> TddftResult {
        TddftResult {
            source: PathBuf::from("t.out"),
            program: "ORCA 6.1.0".into(),
            method: "TD-DFT/TDA".into(),
            unrestricted: true,
            ground_state_s2: Some(8.75),
            states,
            notes: Vec::new(),
        }
    }

    #[test]
    fn a_run_without_an_absorption_table_has_no_intensities_to_plot() {
        let r = result(vec![state(1, 1.9, None), state(2, 2.1, None)]);
        assert!(!r.has_intensities());
        assert!(r.peaks_nm().is_empty());
    }

    /// A dark state is not a missing state: it belongs on the axis at zero
    /// height, and dropping it would misalign hover indices with the table.
    #[test]
    fn dark_states_are_kept_as_zero_height_peaks() {
        let r = result(vec![state(1, 1.9, Some(0.0)), state(2, 2.1, Some(0.5))]);
        assert!(r.has_intensities());
        assert_eq!(r.peaks_nm().len(), 2);
        assert_eq!(r.peaks_nm()[0].1, 0.0);
        assert_eq!(r.states_with_intensity().len(), 2);
    }

    #[test]
    fn nm_and_ev_peaks_describe_the_same_states_in_the_same_order() {
        let r = result(vec![state(1, 1.9, Some(0.1)), state(2, 3.1, Some(0.2))]);
        let nm = r.peaks_nm();
        let ev = r.peaks_ev();
        assert_eq!(nm.len(), ev.len());
        for (i, ((x_nm, y_nm), (x_ev, y_ev))) in nm.iter().zip(&ev).enumerate() {
            assert_eq!(y_nm, y_ev, "intensity differs at {i}");
            assert!((x_nm - ev_to_nm(*x_ev)).abs() < 1.0e-9);
        }
    }

    /// A 400 nm photon is 3.0996 eV; the round trip must not drift.
    #[test]
    fn the_wavelength_conversion_round_trips() {
        let nm = ev_to_nm(3.099_604_96);
        assert!((nm - 400.0).abs() < 1.0e-4, "got {nm}");
        assert!((ev_to_nm(ev_to_nm(2.5)) - 2.5).abs() < 1.0e-9);
    }

    #[test]
    fn a_nonpositive_energy_has_no_wavelength_rather_than_an_infinite_one() {
        assert_eq!(ev_to_nm(0.0), 0.0);
        assert_eq!(ev_to_nm(-1.0), 0.0);
    }
}
