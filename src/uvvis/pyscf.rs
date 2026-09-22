//! Reading PySCF's excited states into the same structures the ORCA reader
//! fills.
//!
//! The layout comes from PySCF's own `analyze()` in `pyscf/tdscf/rhf.py`, not
//! from a sample output. Its two format strings are
//!
//! ```text
//! 'Excited State %3d: %12.5f eV %9.2f nm  f=%.4f'
//! 'Excited State %3d: %4s %12.5f eV %9.2f nm  f=%.4f'   (with symmetry)
//! ```
//!
//! printing as
//!
//! ```text
//! Excited State   1:      7.61100 eV    162.97 nm  f=0.0123
//! Excited State   1:   B2      7.61100 eV    162.97 nm  f=0.0123
//! ```
//!
//! with contributions beneath as `'    %4d -> %-4d %12.5f'`.
//!
//! Note the difference from Psi4, which prints the same quantity in atomic
//! units: here it is already electronvolts. Assuming one layout for both would
//! be wrong by a factor of twenty-seven, which is exactly the kind of mistake
//! that looks plausible on a plot.

use std::path::{Path, PathBuf};

use super::types::{ExcitedState, Excitation, Spin, TddftResult, HARTREE_TO_EV};

/// Reciprocal centimetres per electronvolt.
const EV_TO_CM1: f64 = 8065.543_937_0;

/// The line every root starts with.
const STATE_PREFIX: &str = "Excited State";

/// Whether this text looks like PySCF output carrying a spectrum.
///
/// PySCF prints no banner of its own by default, so this keys on the heading
/// its `analyze()` writes above the roots.
pub fn is_pyscf_output(text: &str) -> bool {
    text.contains("excitation energies and oscillator strengths")
}

/// Reads a PySCF output's excited-state section.
pub fn parse_pyscf_tddft(text: &str, source: &Path) -> Result<TddftResult, String> {
    let lines: Vec<&str> = text.lines().collect();
    let mut states: Vec<ExcitedState> = Vec::new();
    let mut notes: Vec<String> = Vec::new();

    let mut index = 0usize;
    while index < lines.len() {
        let trimmed = lines[index].trim();
        if !trimmed.starts_with(STATE_PREFIX) {
            index += 1;
            continue;
        }
        match parse_state_line(trimmed) {
            Some(mut state) => {
                index += 1;
                while index < lines.len() {
                    let line = lines[index].trim();
                    if line.is_empty() {
                        index += 1;
                        continue;
                    }
                    if line.starts_with(STATE_PREFIX) {
                        break;
                    }
                    match parse_excitation_line(line) {
                        Some(excitation) => {
                            state.excitations.push(excitation);
                            index += 1;
                        }
                        None => break,
                    }
                }
                states.push(state);
            }
            None => index += 1,
        }
    }

    if states.is_empty() {
        return Err(
            "no PySCF excited states found: expected lines beginning \"Excited State\", as a \
             TDDFT analyze() writes"
                .to_string(),
        );
    }

    // PySCF names the manifold in its heading, and solves one at a time.
    let triplet = text.contains("Triplet excitation energies");
    // TDA and the full response are different classes there, and the class
    // name is what the log carries.
    let tda = text.contains("TDA") || text.contains("CIS");
    let method = match (triplet, tda) {
        (false, true) => "TD-DFT/TDA",
        (false, false) => "TD-DFT",
        (true, true) => "TD-DFT/TDA (triplets)",
        (true, false) => "TD-DFT (triplets)",
    };

    if states.iter().all(|s| s.oscillator_strength.is_none()) {
        notes.push(
            "This output lists no oscillator strengths, so the states are here but there is \
             nothing to plot."
                .to_string(),
        );
    }
    if triplet {
        notes.push(
            "These are triplet excitations. They carry no oscillator strength from a singlet \
             ground state, so the spectrum they plot is not an absorption spectrum."
                .to_string(),
        );
    }

    Ok(TddftResult {
        source: PathBuf::from(source),
        program: "PySCF".to_string(),
        method: method.to_string(),
        unrestricted: text.contains("UKS") || text.contains("UHF"),
        ground_state_s2: None,
        states,
        notes,
    })
}

/// One `Excited State` line.
///
/// The energy here is already in electronvolts, unlike Psi4's.
fn parse_state_line(line: &str) -> Option<ExcitedState> {
    let rest = line.strip_prefix(STATE_PREFIX)?;
    let (root_text, tail) = rest.split_once(':')?;
    let root: usize = root_text.trim().parse().ok()?;

    let energy_ev = number_before(tail, "eV")?;
    let wavelength_nm = number_before(tail, "nm").unwrap_or(0.0);
    let oscillator_strength = tail
        .split_once("f=")
        .or_else(|| tail.split_once("f ="))
        .and_then(|(_, f)| f.split_whitespace().next())
        .and_then(|token| token.parse::<f64>().ok());

    Some(ExcitedState {
        root,
        energy_au: energy_ev / HARTREE_TO_EV,
        energy_ev,
        energy_cm1: energy_ev * EV_TO_CM1,
        wavelength_nm: if wavelength_nm > 0.0 {
            wavelength_nm
        } else {
            super::types::ev_to_nm(energy_ev)
        },
        s2: None,
        multiplicity: None,
        oscillator_strength,
        d2_au2: None,
        transition_dipole_au: None,
        excitations: Vec::new(),
    })
}

/// The number immediately before `unit`.
fn number_before(text: &str, unit: &str) -> Option<f64> {
    let at = text.find(unit)?;
    text[..at]
        .split_whitespace()
        .next_back()
        .and_then(|token| token.parse().ok())
}

/// One contribution, as `    4 -> 5        0.68901`.
fn parse_excitation_line(line: &str) -> Option<Excitation> {
    let (left, right) = line.split_once("->")?;
    let from_orbital: u32 = left.split_whitespace().next()?.parse().ok()?;
    let mut tokens = right.split_whitespace();
    let to_orbital: u32 = tokens.next()?.parse().ok()?;
    let coefficient: f64 = tokens.next()?.parse().ok()?;

    Some(Excitation {
        from_orbital,
        from_spin: Spin::Unspecified,
        to_orbital,
        to_spin: Spin::Unspecified,
        // PySCF prints the coefficient; squaring gives the weight the panel
        // ranks by, the same quantity every other reader stores.
        weight: coefficient * coefficient,
        coefficient: Some(coefficient),
        from_label: None,
        to_label: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Laid out as PySCF's own `analyze()` prints it.
    const OUTPUT: &str = concat!(
        "\n** Singlet excitation energies and oscillator strengths **\n",
        "Excited State   1:      7.61100 eV    162.97 nm  f=0.0123\n",
        "       4 -> 5         0.68901\n",
        "       3 -> 6        -0.12000\n",
        "Excited State   2:      9.51300 eV    130.39 nm  f=0.0000\n",
        "       4 -> 6         0.70100\n",
    );

    #[test]
    fn the_roots_are_read_with_their_energies_and_strengths() {
        let result =
            parse_pyscf_tddft(OUTPUT, Path::new("water.log")).expect("two roots are here");
        assert_eq!(result.states.len(), 2);

        let first = &result.states[0];
        assert_eq!(first.root, 1);
        // Already electronvolts here, unlike Psi4.
        assert!((first.energy_ev - 7.611).abs() < 1.0e-9, "{first:?}");
        assert!((first.energy_au - 7.611 / HARTREE_TO_EV).abs() < 1.0e-9);
        assert!((first.wavelength_nm - 162.97).abs() < 1.0e-9);
        assert_eq!(first.oscillator_strength, Some(0.0123));
        assert_eq!(result.states[1].oscillator_strength, Some(0.0));
    }

    /// The unit trap. PySCF prints electronvolts where Psi4 prints atomic
    /// units; reading one as the other is wrong by a factor of twenty-seven
    /// and would look entirely plausible on a plot.
    #[test]
    fn the_energy_is_taken_as_electronvolts_not_hartree() {
        let result = parse_pyscf_tddft(OUTPUT, Path::new("x.log")).unwrap();
        let first = &result.states[0];
        assert!(
            first.energy_ev > 1.0,
            "7.6 eV, not 7.6 Eh: {}",
            first.energy_ev
        );
        assert!(first.energy_au < 1.0, "{}", first.energy_au);
    }

    /// A symmetry label sits between the colon and the energy, and must not
    /// shift which number is read.
    #[test]
    fn a_symmetry_label_does_not_shift_the_energy() {
        let symmetric = concat!(
            "** Singlet excitation energies and oscillator strengths **\n",
            "Excited State   1:   B2      7.61100 eV    162.97 nm  f=0.0123\n",
        );
        let result = parse_pyscf_tddft(symmetric, Path::new("x.log")).unwrap();
        assert!((result.states[0].energy_ev - 7.611).abs() < 1.0e-9);
        assert!((result.states[0].wavelength_nm - 162.97).abs() < 1.0e-9);
    }

    #[test]
    fn contributions_are_stored_as_weights() {
        let result = parse_pyscf_tddft(OUTPUT, Path::new("x.log")).unwrap();
        let first = &result.states[0];
        assert_eq!(first.excitations.len(), 2);
        assert_eq!(first.excitations[0].from_orbital, 4);
        assert_eq!(first.excitations[0].to_orbital, 5);
        assert!((first.excitations[0].weight - 0.68901f64.powi(2)).abs() < 1.0e-9);
        assert!(first.excitations[1].coefficient.unwrap() < 0.0);
    }

    /// A triplet run is labelled as one and carries a warning, because those
    /// roots have no oscillator strength from a singlet ground state and the
    /// curve they draw is not an absorption spectrum.
    #[test]
    fn a_triplet_run_is_labelled_and_flagged() {
        let triplets = OUTPUT.replace("Singlet excitation", "Triplet excitation");
        let result = parse_pyscf_tddft(&triplets, Path::new("x.log")).unwrap();
        assert!(result.method.contains("triplets"), "{}", result.method);
        assert!(
            result.notes.iter().any(|n| n.contains("not an absorption spectrum")),
            "{:?}",
            result.notes
        );
    }

    #[test]
    fn output_without_excited_states_is_an_error() {
        let err = parse_pyscf_tddft("converged SCF energy = -76.0\n", Path::new("x.log"))
            .expect_err("no roots here");
        assert!(err.contains("Excited State"), "{err}");
    }

    #[test]
    fn a_pyscf_spectrum_is_recognised() {
        assert!(is_pyscf_output(OUTPUT));
        assert!(!is_pyscf_output("ORCA TERMINATED NORMALLY\n"));
    }
}
