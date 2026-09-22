//! Reading Psi4's excited states into the same structures the ORCA reader
//! fills.
//!
//! `types.rs` was written to be program-independent for exactly this, so the
//! panel does not change: a spectrum from Psi4 and one from ORCA are the same
//! thing by the time they reach it.
//!
//! The layout is taken from Psi4's own printing code, in
//! `driver/procrouting/response/scf_response.py`, rather than from a sample
//! output. Its format string is
//!
//! ```text
//! f"\nExcited State {i+1:4d} ({mult} {irrep}):"
//! f"{E_ex_au:> 10.5f} au   {E_ex_nm: >.2f} nm f = {f_length: >.4f}\n"
//! ```
//!
//! which prints as
//!
//! ```text
//! Excited State    1 (1 B2):   0.27971 au   163.00 nm f =  0.0123
//! ```
//!
//! and the contributing excitations beneath it as
//!
//! ```text
//!    12 (B2) -> 13 (A1)     0.68901
//! ```
//!
//! Psi4 reports no per-root `<S²>`, so that stays `None`; a restricted ORCA run
//! has none either and the panel already copes.

use std::path::{Path, PathBuf};

use super::types::{ExcitedState, Excitation, Spin, TddftResult, HARTREE_TO_EV};

/// Reciprocal centimetres per electronvolt.
const EV_TO_CM1: f64 = 8065.543_937_0;

/// The line every root starts with.
const STATE_PREFIX: &str = "Excited State";

/// Whether this text looks like a Psi4 output at all.
///
/// Keyed on Psi4's own banner rather than on the excited-state line, so a Psi4
/// file with no excited states is recognised as Psi4 and refused with a useful
/// message rather than being mistaken for another program's.
pub fn is_psi4_output(text: &str) -> bool {
    text.contains("Psi4: An Open-Source Ab Initio Electronic Structure Package")
        || text.contains("Psi4 started on")
        || (text.contains("Psi4") && text.contains("Excited State"))
}

/// Reads a Psi4 output's excited-state section.
///
/// Fails only when there are no roots at all. A run cut short, or one that
/// found fewer roots than were asked for, is a normal result and comes back
/// with whatever it has.
pub fn parse_psi4_tddft(text: &str, source: &Path) -> Result<TddftResult, String> {
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
                // The contributing excitations follow, one per line, until
                // something that is not one.
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
            "no Psi4 excited states found: expected lines beginning \"Excited State\", as a \
             tdscf_excitations run writes"
                .to_string(),
        );
    }

    // Tamm-Dancoff or the full response, as Psi4 states in its own header.
    let tda = text.contains("Tamm-Dancoff") || text.contains("TDA");
    let method = if tda { "TD-DFT/TDA" } else { "TD-DFT" };

    if states.iter().all(|s| s.oscillator_strength.is_none()) {
        notes.push(
            "This output lists no oscillator strengths, so the states are here but there is \
             nothing to plot."
                .to_string(),
        );
    }

    Ok(TddftResult {
        source: PathBuf::from(source),
        program: program_name(text),
        method: method.to_string(),
        // Psi4's excited-state print does not state the reference's spin, and
        // guessing it from the roots would be a guess.
        unrestricted: text.contains("UKS") || text.contains("UHF"),
        ground_state_s2: None,
        states,
        notes,
    })
}

/// `Psi4 1.11`, or just `Psi4` when the version is not printed.
fn program_name(text: &str) -> String {
    for line in text.lines().take(80) {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix("Psi4 ") {
            if let Some(version) = rest.split_whitespace().next() {
                if version.chars().next().is_some_and(|c| c.is_ascii_digit()) {
                    return format!("Psi4 {version}");
                }
            }
        }
    }
    "Psi4".to_string()
}

/// One `Excited State` line.
///
/// The energy is in atomic units here, unlike PySCF which prints
/// electronvolts, so it is converted rather than assumed.
fn parse_state_line(line: &str) -> Option<ExcitedState> {
    let rest = line.strip_prefix(STATE_PREFIX)?;
    // `   1 (1 B2):   0.27971 au   163.00 nm f =  0.0123`
    let root: usize = rest.split_whitespace().next()?.parse().ok()?;

    // Everything after the closing bracket of the symmetry label, or after the
    // root number when there is no bracket.
    let tail = match rest.find("):") {
        Some(at) => &rest[at + 2..],
        None => rest.split_once(':').map(|(_, t)| t)?,
    };

    let energy_au = number_before(tail, "au")?;
    let wavelength_nm = number_before(tail, "nm").unwrap_or(0.0);
    let oscillator_strength = tail
        .split_once("f =")
        .or_else(|| tail.split_once("f="))
        .and_then(|(_, f)| f.split_whitespace().next())
        .and_then(|token| token.parse::<f64>().ok());

    let energy_ev = energy_au * HARTREE_TO_EV;
    Some(ExcitedState {
        root,
        energy_au,
        energy_ev,
        energy_cm1: energy_ev * EV_TO_CM1,
        wavelength_nm: if wavelength_nm > 0.0 {
            wavelength_nm
        } else {
            super::types::ev_to_nm(energy_ev)
        },
        s2: None,
        multiplicity: multiplicity_from(rest),
        oscillator_strength,
        d2_au2: None,
        transition_dipole_au: None,
        excitations: Vec::new(),
    })
}

/// The multiplicity Psi4 prints inside the bracket, as in `(1 B2)`.
fn multiplicity_from(rest: &str) -> Option<u32> {
    let open = rest.find('(')?;
    let close = rest.find(')')?;
    if close <= open {
        return None;
    }
    rest[open + 1..close].split_whitespace().next()?.parse().ok()
}

/// The number immediately before `unit` on this line.
///
/// Taken by looking back from the unit rather than by counting tokens, so a
/// symmetry label or an extra column cannot shift it.
fn number_before(text: &str, unit: &str) -> Option<f64> {
    let at = text.find(unit)?;
    text[..at]
        .split_whitespace()
        .next_back()
        .and_then(|token| token.parse().ok())
}

/// One contributing excitation, as `12 (B2) -> 13 (A1)   0.68901`.
///
/// The symmetry labels are optional; a run without symmetry prints
/// `12 -> 13   0.68901`.
fn parse_excitation_line(line: &str) -> Option<Excitation> {
    let (left, right) = line.split_once("->")?;
    let from_orbital: u32 = left.split_whitespace().next()?.parse().ok()?;

    let mut right_tokens = right.split_whitespace();
    let to_orbital: u32 = right_tokens.next()?.parse().ok()?;
    // The coefficient is the last number on the line.
    let coefficient: f64 = right
        .split_whitespace()
        .filter_map(|token| token.parse::<f64>().ok())
        .next_back()?;

    // Psi4 prints the coefficient, not the weight. Squaring it gives the same
    // quantity the other readers store, so the panel ranks contributions the
    // same way whichever program produced them.
    Some(Excitation {
        from_orbital,
        from_spin: Spin::Unspecified,
        to_orbital,
        to_spin: Spin::Unspecified,
        weight: coefficient * coefficient,
        coefficient: Some(coefficient),
        from_label: None,
        to_label: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Text laid out exactly as Psi4's own print statement produces it.
    const OUTPUT: &str = concat!(
        "    Psi4 1.11 release\n",
        "  Psi4 started on somehost\n",
        "                        Tamm-Dancoff approximation\n",
        "\n",
        "Contributing excitations\n",
        "Only contributions with abs(coeff) > 1.00e-01 will be printed:\n",
        "\n",
        "Excited State    1 (1 B2):   0.27971 au   162.97 nm f =  0.0123\n",
        "    4 (B2) -> 5 (A1)     0.68901\n",
        "    3 (A1) -> 6 (B2)    -0.12000\n",
        "\n",
        "Excited State    2 (1 A1):   0.34959 au   130.39 nm f =  0.0000\n",
        "    4 (B2) -> 6 (B2)     0.70100\n",
    );

    #[test]
    fn the_roots_are_read_with_their_energies_and_strengths() {
        let result =
            parse_psi4_tddft(OUTPUT, Path::new("water.out")).expect("two roots are here");
        assert_eq!(result.states.len(), 2);

        let first = &result.states[0];
        assert_eq!(first.root, 1);
        // Psi4 prints atomic units, so the electronvolt figure is derived.
        assert!((first.energy_au - 0.27971).abs() < 1.0e-9, "{first:?}");
        assert!((first.energy_ev - 0.27971 * HARTREE_TO_EV).abs() < 1.0e-6);
        assert!((first.wavelength_nm - 162.97).abs() < 1.0e-6);
        assert_eq!(first.oscillator_strength, Some(0.0123));
        assert_eq!(first.multiplicity, Some(1));

        // A dark state is a strength of zero, not a missing one.
        assert_eq!(result.states[1].oscillator_strength, Some(0.0));
    }

    /// Psi4 prints the coefficient; the panel ranks by weight. Squaring is what
    /// makes a Psi4 state rank the same way an ORCA one does.
    #[test]
    fn contributions_are_stored_as_weights() {
        let result = parse_psi4_tddft(OUTPUT, Path::new("water.out")).unwrap();
        let first = &result.states[0];
        assert_eq!(first.excitations.len(), 2);
        assert_eq!(first.excitations[0].from_orbital, 4);
        assert_eq!(first.excitations[0].to_orbital, 5);
        assert!((first.excitations[0].coefficient.unwrap() - 0.68901).abs() < 1.0e-9);
        assert!((first.excitations[0].weight - 0.68901f64.powi(2)).abs() < 1.0e-9);

        // The dominant one is the larger coefficient, regardless of sign.
        let dominant = first.dominant_excitation().unwrap();
        assert_eq!(dominant.from_orbital, 4);
        assert_eq!(dominant.to_orbital, 5);
    }

    /// A negative coefficient is a real contribution, not a parse failure.
    #[test]
    fn a_negative_coefficient_is_kept() {
        let result = parse_psi4_tddft(OUTPUT, Path::new("water.out")).unwrap();
        let second = &result.states[0].excitations[1];
        assert!(second.coefficient.unwrap() < 0.0);
        assert!(second.weight > 0.0, "a weight is never negative");
    }

    /// Without symmetry Psi4 prints no bracketed labels, and the same lines
    /// still have to parse.
    #[test]
    fn output_without_symmetry_labels_parses() {
        let plain = concat!(
            "Psi4 1.11\n",
            "Excited State    1 (1 A):   0.27971 au   162.97 nm f =  0.0123\n",
            "    4 -> 5     0.68901\n",
        );
        let result = parse_psi4_tddft(plain, Path::new("x.out")).unwrap();
        assert_eq!(result.states.len(), 1);
        assert_eq!(result.states[0].excitations.len(), 1);
        assert_eq!(result.states[0].excitations[0].to_orbital, 5);
    }

    /// The Tamm-Dancoff switch is reported, because it is a different model
    /// and the panel labels what actually ran.
    #[test]
    fn the_response_model_is_named() {
        let tda = parse_psi4_tddft(OUTPUT, Path::new("x.out")).unwrap();
        assert_eq!(tda.method, "TD-DFT/TDA");

        let full = OUTPUT.replace("Tamm-Dancoff approximation", "full linear response");
        let full = full.replace("TDA", "");
        let parsed = parse_psi4_tddft(&full, Path::new("x.out")).unwrap();
        assert_eq!(parsed.method, "TD-DFT");
    }

    /// The version is carried when it is printed.
    #[test]
    fn the_program_version_is_read_when_present() {
        let result = parse_psi4_tddft(OUTPUT, Path::new("x.out")).unwrap();
        assert_eq!(result.program, "Psi4 1.11");
    }

    /// A file with no excited states is an error naming what was expected,
    /// rather than an empty result that looks like a dark spectrum.
    #[test]
    fn output_without_excited_states_is_an_error() {
        let err = parse_psi4_tddft("Psi4 1.11\nTotal Energy = -76.0\n", Path::new("x.out"))
            .expect_err("no roots here");
        assert!(err.contains("Excited State"), "{err}");
    }

    /// Recognising the program is separate from finding roots, so a Psi4 file
    /// without a spectrum is still known to be Psi4.
    #[test]
    fn a_psi4_file_is_recognised_even_without_a_spectrum() {
        assert!(is_psi4_output("  Psi4 started on host\n  Total Energy = -76.0\n"));
        assert!(!is_psi4_output("ORCA TERMINATED NORMALLY\n"));
    }
}
