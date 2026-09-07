//! Reading Behemoth's sTDA / sTD-DFT spectrum into the same structures the
//! ORCA reader fills.
//!
//! `types.rs` was written to be program-independent for exactly this, so
//! nothing in the UI changes: a spectrum Beavyr ran and a spectrum it loaded
//! from a file are the same thing by the time they reach the panel.
//!
//! Two blocks are read. The root table is the spectrum:
//!
//! ```text
//! root       energy(eV)     lambda(nm)            f
//! 1            4.115545         301.26     0.000000
//! 2            8.338299         148.69     0.173763
//! ```
//!
//! and the decomposition section underneath gives each root its leading
//! configurations, with the engine's own name in the header -- `sTDA` for
//! `--stda`, `sTD-DFT` for `--stddft`:
//!
//! ```text
//! Root 1   sTDA    =     4.1155 eV      301.26 nm  f = 0.00000
//!        8 HOMO       [O lonepair] ==>    9 LUMO       [C-O antibond]  99.89%
//! ```
//!
//! Behemoth reports no per-root `<S^2>` and no transition dipole, so those
//! stay `None`. The panel already copes -- a restricted ORCA run has no
//! `<S^2>` either.

use std::path::{Path, PathBuf};

use super::types::{ExcitedState, Excitation, Spin, TddftResult, HARTREE_TO_EV};

/// Reciprocal centimetres per electronvolt, for the `energy_cm1` field.
const EV_TO_CM1: f64 = 8065.543_937_0;

/// Reads a Behemoth log containing an sTDA or sTD-DFT spectrum.
///
/// Fails only when there is no root table at all: a run that was cut short by
/// its energy window, or that produced fewer roots than were asked for, is a
/// normal result and is returned with whatever roots it has.
pub fn parse_behemoth_spectrum(
    stdout: &str,
    source: &Path,
    requested_roots: Option<usize>,
) -> Result<TddftResult, String> {
    let mut notes = Vec::new();
    let states_from_table = parse_root_table(stdout);
    if states_from_table.is_empty() {
        // Behemoth's own error is far more use than anything invented here,
        // so it is passed through when there is one.
        if let Some(error) = stdout
            .lines()
            .find(|line| line.trim_start().starts_with("Error:"))
        {
            return Err(error.trim().to_string());
        }
        return Err("No excited-state table found in the output.".to_string());
    }

    let decompositions = parse_decompositions(stdout);
    // Behemoth solves one spin channel per run and names it in its parameter
    // line, so every root shares that multiplicity. ORCA reports the same
    // thing per root from a rounded <S^2>; this is the same information from a
    // more direct source.
    let multiplicity = spin_multiplicity(stdout);
    let mut states = states_from_table;
    for state in &mut states {
        if let Some(excitations) = decompositions.get(&state.root) {
            state.excitations = excitations.clone();
        }
        state.multiplicity = multiplicity;
    }
    let missing = states.iter().filter(|s| s.excitations.is_empty()).count();
    if missing > 0 {
        notes.push(format!(
            "{missing} of {} roots have no listed configurations; the spectrum is \
             unaffected but those states cannot be assigned.",
            states.len()
        ));
    }

    // The caller knows what it asked for; Behemoth does not echo the root
    // count, so reading it back out of the log would never find one. The log
    // is still consulted, for a log that came from somewhere else.
    if let Some(requested) = requested_roots.or_else(|| echoed_roots(stdout)) {
        if states.len() < requested {
            // The commonest surprise with this engine, and not an error: the
            // energy window silently bounds the root count.
            let window = energy_window_ev(stdout)
                .map(|ev| format!(" (Emax = {ev} eV)"))
                .unwrap_or_default();
            notes.push(format!(
                "{requested} roots were requested but {} are within the energy window{window}. \
                 Raise the window to see more.",
                states.len()
            ));
        }
    }

    // Said once, so the absence of the dipole line in the state list is
    // explained rather than looking like a parse failure. Everything else the
    // ORCA reader shows -- root, energy in eV, wavelength, oscillator
    // strength, and the orbital pairs with their percentages -- is present.
    notes.push(
        "Behemoth reports no transition-dipole components or <S^2> per root, so those \
         rows are absent; everything else is as an ORCA output gives it."
            .to_string(),
    );

    Ok(TddftResult {
        source: PathBuf::from(source),
        program: "Behemoth".to_string(),
        method: method_label(stdout),
        // Behemoth prints one spin channel for a restricted reference and
        // labels an unrestricted one explicitly.
        unrestricted: stdout.contains("unrestricted sTDA") || stdout.contains("Spin = alpha"),
        ground_state_s2: None,
        states,
        notes,
    })
}

/// The engine and the orbitals it ran on, for the panel's method line.
///
/// Behemoth names the orbital source in a header line of its own:
/// `Orbital source: Behemoth RKS PBE0 (this response model ...)`. The
/// parenthetical caveat is dropped -- it belongs in the documentation, not in
/// a one-line label.
fn method_label(stdout: &str) -> String {
    let engine = if stdout.contains("sTD-DFT") {
        "sTD-DFT"
    } else {
        "sTDA"
    };
    let source = stdout
        .lines()
        .find_map(|line| line.trim().strip_prefix("Orbital source:"))
        .map(|rest| rest.split('(').next().unwrap_or(rest).trim().to_string())
        .filter(|s| !s.is_empty());
    match source {
        Some(source) => format!("{engine} / {source}"),
        None => engine.to_string(),
    }
}

/// `--stddft-roots` as the run echoed it, for a log that came from somewhere
/// other than this panel and so has no requested count to compare against.
fn echoed_roots(stdout: &str) -> Option<usize> {
    stdout.lines().find_map(|line| {
        let rest = line.trim().strip_prefix("Requested roots")?;
        rest.trim_start()
            .trim_start_matches([':', '='])
            .split_whitespace()
            .next()?
            .parse()
            .ok()
    })
}

/// The multiplicity of the channel that was solved, from the parameter line
/// `Spin = singlet, Emax = 10.000 eV, ...`.
fn spin_multiplicity(stdout: &str) -> Option<u32> {
    for line in stdout.lines() {
        for field in line.split(',') {
            let Some(rest) = field.trim().strip_prefix("Spin") else {
                continue;
            };
            let name = rest.trim_start().trim_start_matches('=').trim();
            return match name {
                "singlet" => Some(1),
                "doublet" => Some(2),
                "triplet" => Some(3),
                _ => None,
            };
        }
    }
    None
}

/// The `Emax` from the parameter line: `Spin = singlet, Emax = 10.000 eV, ...`
fn energy_window_ev(stdout: &str) -> Option<f64> {
    for line in stdout.lines() {
        for field in line.split(',') {
            let Some(rest) = field.trim().strip_prefix("Emax") else {
                continue;
            };
            let value = rest.trim_start().strip_prefix('=')?;
            return value.split_whitespace().next()?.parse().ok();
        }
    }
    None
}

/// The root table: four numbers a line, the first a bare integer.
///
/// Recognised by shape rather than by counting lines after the header, so the
/// banners around it can change without this silently reading the wrong rows.
/// The header itself begins with the word `root`, which is not an integer, and
/// is skipped by the same rule.
fn parse_root_table(stdout: &str) -> Vec<ExcitedState> {
    let mut states = Vec::new();
    for line in stdout.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() != 4 {
            continue;
        }
        let Ok(root) = fields[0].parse::<usize>() else {
            continue;
        };
        let (Ok(energy_ev), Ok(wavelength_nm), Ok(oscillator)) = (
            fields[1].parse::<f64>(),
            fields[2].parse::<f64>(),
            fields[3].parse::<f64>(),
        ) else {
            continue;
        };
        // A root's energy is positive and its wavelength finite; anything else
        // is a four-number line that happens not to be a root.
        if energy_ev <= 0.0 || wavelength_nm <= 0.0 {
            continue;
        }
        // A repeated root number means a second table -- the decomposition
        // section reprints the same header and rows -- so the first wins.
        if states.iter().any(|s: &ExcitedState| s.root == root) {
            continue;
        }
        states.push(ExcitedState {
            root,
            energy_au: energy_ev / HARTREE_TO_EV,
            energy_ev,
            energy_cm1: energy_ev * EV_TO_CM1,
            // Behemoth's own value, rather than a re-derivation, so the table
            // matches the output it came from.
            wavelength_nm,
            s2: None,
            multiplicity: None, // filled below, from the run's spin channel

            oscillator_strength: Some(oscillator),
            d2_au2: None,
            transition_dipole_au: None,
            excitations: Vec::new(),
        });
    }
    states.sort_by_key(|s| s.root);
    states
}

/// The leading configurations per root, from the decomposition section.
fn parse_decompositions(stdout: &str) -> std::collections::HashMap<usize, Vec<Excitation>> {
    let mut per_root: std::collections::HashMap<usize, Vec<Excitation>> =
        std::collections::HashMap::new();
    let mut current: Option<usize> = None;
    for line in stdout.lines() {
        let trimmed = line.trim();
        if let Some(root) = parse_root_header(trimmed) {
            current = Some(root);
            per_root.entry(root).or_default();
            continue;
        }
        if !trimmed.contains("==>") {
            // A blank line or a banner ends the block; anything else between
            // roots is prose.
            if trimmed.is_empty() {
                current = None;
            }
            continue;
        }
        let Some(root) = current else {
            continue;
        };
        if let Some(excitation) = parse_contribution(trimmed) {
            per_root.entry(root).or_default().push(excitation);
        }
    }
    per_root
}

/// `Root 1   sTDA    =     4.1155 eV      301.26 nm  f = 0.00000` -> `1`.
///
/// The engine name between the root number and the energy differs by mode, so
/// it is not matched: only the `Root <n>` prefix is, which both modes share.
fn parse_root_header(line: &str) -> Option<usize> {
    let rest = line.strip_prefix("Root ")?;
    let number = rest.split_whitespace().next()?;
    let root = number.parse().ok()?;
    // The header carries an energy; a bare "Root 3" in prose does not.
    line.contains("eV").then_some(root)
}

/// One contribution line:
///
/// ```text
/// 8 HOMO       [O lonepair] ==>    9 LUMO       [C-O antibond]  99.89%
/// ```
///
/// The bracketed character labels contain spaces, so the two sides are split
/// on the bracket rather than on whitespace.
fn parse_contribution(line: &str) -> Option<Excitation> {
    let (left, right) = line.split_once("==>")?;
    let (from_orbital, from_label) = parse_orbital(left)?;
    // The percentage is whatever follows the right-hand side's closing
    // bracket, or its last field when there is no bracket.
    let (to_orbital, to_label) = parse_orbital(right)?;
    let percent_field = match right.rsplit_once(']') {
        Some((_, tail)) => tail.split_whitespace().next(),
        None => right.split_whitespace().nth(2),
    }?;
    let percent: f64 = percent_field.trim_end_matches('%').parse().ok()?;

    Some(Excitation {
        from_orbital,
        from_spin: Spin::Unspecified,
        to_orbital,
        to_spin: Spin::Unspecified,
        // Stored as a fraction, as everywhere else in `types`, so `percent()`
        // gives back what the output printed.
        weight: percent / 100.0,
        coefficient: None,
        from_label,
        to_label,
    })
}

/// `8 HOMO       [O lonepair]` -> `(8, Some("HOMO [O lonepair]"))`.
///
/// The name is normalised to single spaces: the output pads it into columns,
/// and that padding would otherwise end up inside the label.
fn parse_orbital(side: &str) -> Option<(u32, Option<String>)> {
    let side = side.trim();
    let mut fields = side.split_whitespace();
    let orbital: u32 = fields.next()?.parse().ok()?;

    let name = match side.find('[') {
        Some(open) => {
            let close = side[open..].find(']')? + open;
            // Everything between the index and the bracket, plus the bracket.
            let relative: String = side[..open]
                .split_whitespace()
                .skip(1)
                .collect::<Vec<&str>>()
                .join(" ");
            let character = side[open..=close].to_string();
            let joined = format!("{relative} {character}").trim().to_string();
            (!joined.is_empty()).then_some(joined)
        }
        None => {
            let relative: String = fields.collect::<Vec<&str>>().join(" ");
            // A trailing percentage is not part of the name.
            let relative = relative
                .split_whitespace()
                .filter(|f| !f.ends_with('%'))
                .collect::<Vec<&str>>()
                .join(" ");
            (!relative.is_empty()).then_some(relative)
        }
    };
    Some((orbital, name))
}

/// A sanity check used by the tests: the table prints both the energy and the
/// wavelength, and the two must agree or one column is being read wrong.
#[cfg(test)]
fn table_is_self_consistent(state: &ExcitedState) -> bool {
    let derived = super::types::ev_to_nm(state.energy_ev);
    // The table rounds the wavelength to two decimals.
    (derived - state.wavelength_nm).abs() < 0.02
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::types::ev_to_nm;

    const STDA: &str = include_str!("../../tests/fixtures/behemoth_ch2o_stda.log");
    const STDDFT: &str = include_str!("../../tests/fixtures/behemoth_ch2o_stddft.log");
    const STDA_TASI: &str = include_str!("../../tests/fixtures/behemoth_ch2o_stda_tasi.log");

    fn parse(log: &str) -> TddftResult {
        parse_behemoth_spectrum(log, Path::new("run.log"), None).expect("a spectrum")
    }

    fn parse_asking_for(log: &str, roots: usize) -> TddftResult {
        parse_behemoth_spectrum(log, Path::new("run.log"), Some(roots)).expect("a spectrum")
    }

    /// The TDA run on formaldehyde, with the energies, wavelengths and
    /// oscillator strengths from the table.
    #[test]
    fn the_tda_root_table_is_read() {
        let result = parse(STDA);
        assert_eq!(result.program, "Behemoth");
        assert_eq!(result.states.len(), 6);

        let first = &result.states[0];
        assert_eq!(first.root, 1);
        assert!((first.energy_ev - 4.113967).abs() < 1e-6);
        assert!((first.wavelength_nm - 301.37).abs() < 1e-6);
        // n -> pi* in formaldehyde is symmetry forbidden: a dark state, which
        // is not the same as a state with no intensity reported.
        assert_eq!(first.oscillator_strength, Some(0.0));

        let bright = &result.states[1];
        assert!((bright.energy_ev - 8.333784).abs() < 1e-6);
        assert!((bright.oscillator_strength.unwrap() - 0.171667).abs() < 1e-6);
    }

    /// The coupled run reports six roots, and labels its engine differently,
    /// which the parser has to accept in the decomposition header.
    #[test]
    fn the_coupled_run_is_read_and_named_differently() {
        let result = parse(STDDFT);
        assert_eq!(result.states.len(), 6);
        assert!(result.method.starts_with("sTD-DFT"), "{}", result.method);
        assert!(result.method.contains("PBE0"), "{}", result.method);

        // Every root got its configurations despite the different label.
        for state in &result.states {
            assert!(
                !state.excitations.is_empty(),
                "root {} has no configurations",
                state.root
            );
        }
    }

    #[test]
    fn the_tda_run_is_named_for_its_engine_and_orbitals() {
        let result = parse(STDA);
        assert!(result.method.starts_with("sTDA"), "{}", result.method);
        assert!(result.method.contains("PBE0"), "{}", result.method);
        // The parenthetical caveat is not part of a one-line label.
        assert!(!result.method.contains("response model"), "{}", result.method);
    }

    #[test]
    fn a_tasi_run_names_tasi_as_its_orbital_source() {
        let result = parse(STDA_TASI);
        assert!(result.method.contains("TASI"), "{}", result.method);
        assert!(!result.states.is_empty());
    }

    /// The configurations, with the orbital names and character labels that
    /// are the most useful part of the line.
    #[test]
    fn a_configuration_carries_its_orbital_names_and_characters() {
        let result = parse(STDA);
        let first = result.states[0].excitations.first().expect("a contribution");
        assert_eq!(first.from_orbital, 8);
        assert_eq!(first.to_orbital, 9);
        assert!((first.percent() - 99.89).abs() < 0.01, "{}", first.percent());
        assert_eq!(first.from_label.as_deref(), Some("HOMO [O lonepair]"));
        assert_eq!(first.to_label.as_deref(), Some("LUMO [C-O antibond]"));
        assert_eq!(
            first.described(),
            "8 HOMO [O lonepair] --> 9 LUMO [C-O antibond]"
        );
    }

    /// A root with more than one leading configuration keeps all of them, in
    /// the order printed.
    #[test]
    fn a_root_with_several_configurations_keeps_them_all() {
        let result = parse(STDDFT);
        let root4 = result.states.iter().find(|s| s.root == 4).expect("root 4");
        assert_eq!(root4.excitations.len(), 2, "{:?}", root4.excitations);
        assert!((root4.excitations[0].percent() - 64.46).abs() < 0.01);
        assert!((root4.excitations[1].percent() - 33.83).abs() < 0.01);
        // Ranked by weight, the larger one leads.
        assert_eq!(root4.dominant_excitation().unwrap().from_orbital, 8);
    }

    /// The column padding must not leak into the stored name.
    #[test]
    fn column_padding_is_not_part_of_an_orbital_name() {
        let (orbital, name) = parse_orbital("       7 HOMO-1     [O-C bond]  ").unwrap();
        assert_eq!(orbital, 7);
        assert_eq!(name.as_deref(), Some("HOMO-1 [O-C bond]"));
    }

    /// A character label contains spaces, so the split cannot be on
    /// whitespace -- this is the line that proves it.
    #[test]
    fn a_multi_word_character_label_survives() {
        let e = parse_contribution(
            "       8 HOMO       [O lonepair] ==>    9 LUMO       [C-O antibond]  99.89%",
        )
        .expect("a contribution");
        assert_eq!(e.from_label.as_deref(), Some("HOMO [O lonepair]"));
        assert_eq!(e.to_label.as_deref(), Some("LUMO [C-O antibond]"));
        assert!((e.weight - 0.9989).abs() < 1e-9);
    }

    /// The energies and wavelengths Behemoth prints have to agree with each
    /// other, or one of the two columns is being read from the wrong place.
    #[test]
    fn every_roots_energy_and_wavelength_agree() {
        for log in [STDA, STDDFT, STDA_TASI] {
            for state in &parse(log).states {
                assert!(
                    table_is_self_consistent(state),
                    "root {}: {} eV is {} nm, table says {}",
                    state.root,
                    state.energy_ev,
                    ev_to_nm(state.energy_ev),
                    state.wavelength_nm
                );
            }
        }
    }

    /// The commonest surprise with this engine: `--stddft-roots` is bounded
    /// by `--stddft-emax`, so asking for more roots than fit in the window
    /// silently returns fewer. That is a note, not an error, and the note has
    /// to be able to quote the window.
    ///
    /// The fixtures were captured with a wide window (25 eV) precisely so
    /// they return the full six; the truncation itself is exercised on a
    /// synthetic log below.
    #[test]
    fn the_energy_window_is_available_to_explain_a_short_root_count() {
        assert_eq!(parse(STDA).states.len(), 6, "the wide window fits all six");
        assert_eq!(
            energy_window_ev(STDA),
            Some(25.0),
            "the window has to be quotable in the note"
        );
    }

    /// Fewer roots than requested is reported as a note against the run,
    /// rather than as a failure or in silence.
    #[test]
    fn fewer_roots_than_requested_is_a_note_not_an_error() {
        let log = "Spin = singlet, Emax = 10.000 eV, ax = 0.250\n\
                   root       energy(eV)     lambda(nm)            f\n\
                   1            4.115545         301.26     0.000000\n\
                   2            8.338299         148.69     0.173763\n";
        // The panel passes the count it asked for, which is the only place it
        // is known.
        let result = parse_asking_for(log, 6);
        assert_eq!(result.states.len(), 2);
        let note = result
            .notes
            .iter()
            .find(|n| n.contains("requested"))
            .unwrap_or_else(|| panic!("no note about the root count: {:?}", result.notes));
        assert!(note.contains('6'), "{note}");
        assert!(note.contains("10"), "the note quotes the window: {note}");
    }

    /// Getting every root asked for must not produce a note about it.
    #[test]
    fn asking_for_exactly_what_came_back_draws_no_note() {
        let result = parse_asking_for(STDA, 6);
        assert_eq!(result.states.len(), 6);
        assert!(
            !result.notes.iter().any(|n| n.contains("requested")),
            "{:?}",
            result.notes
        );
    }

    /// The same information ORCA gives per root from its <S^2>, taken from
    /// the channel Behemoth says it solved.
    #[test]
    fn every_root_carries_the_multiplicity_of_the_channel_that_was_solved() {
        for log in [STDA, STDDFT, STDA_TASI] {
            let result = parse(log);
            for state in &result.states {
                assert_eq!(
                    state.multiplicity,
                    Some(1),
                    "root {} of a singlet run",
                    state.root
                );
            }
        }
        assert_eq!(spin_multiplicity("Spin = triplet, Emax = 10.000 eV"), Some(3));
        assert_eq!(spin_multiplicity("no spin line"), None);
    }

    /// What Behemoth does not print is stated, so the missing dipole row in
    /// the state list reads as an absence rather than a parse failure.
    #[test]
    fn the_fields_behemoth_does_not_report_are_named() {
        let result = parse(STDA);
        assert!(
            result.notes.iter().any(|n| n.contains("transition-dipole")),
            "{:?}",
            result.notes
        );
        assert!(result.states.iter().all(|s| s.transition_dipole_au.is_none()));
        assert!(result.states.iter().all(|s| s.d2_au2.is_none()));
    }

    #[test]
    fn the_energy_window_is_read_from_the_parameter_line() {
        let ev = energy_window_ev(
            "Spin = singlet, Emax = 10.000 eV, ax = 0.250, Gamma^J exponent (yJ) = 0.657",
        );
        assert_eq!(ev, Some(10.0));
        assert!(energy_window_ev("no parameters here").is_none());
    }

    /// The spectrum the panel plots comes straight out of the parsed result.
    #[test]
    fn the_result_feeds_the_existing_plot() {
        let result = parse(STDDFT);
        assert!(result.has_intensities());
        let nm = result.peaks_nm();
        let ev = result.peaks_ev();
        assert_eq!(nm.len(), 6);
        assert_eq!(ev.len(), 6);
        assert_eq!(result.states_with_intensity().len(), 6);
        // Dark roots are kept as zero-height peaks rather than dropped.
        assert!(nm.iter().any(|&(_, f)| f == 0.0));
    }

    /// Behemoth's own error message is more use than a generic one.
    #[test]
    fn behemoths_own_error_is_passed_through() {
        let err = parse_behemoth_spectrum(
            "ERI algorithm: whatever\n\
             Error: sTDA/sTD-DFT with normal DFT requires a global-hybrid functional \
             with nonzero exact exchange\n",
            Path::new("run.log"),
            None,
        )
        .expect_err("this run failed");
        assert!(err.contains("global-hybrid"), "{err}");
    }

    #[test]
    fn a_log_with_no_table_and_no_error_says_so() {
        let err = parse_behemoth_spectrum("nothing useful\n", Path::new("run.log"), None)
            .expect_err("no spectrum");
        assert!(err.contains("No excited-state table"), "{err}");
    }

    /// The decomposition section reprints the root table's header and rows.
    /// Reading them twice would double every state.
    #[test]
    fn a_reprinted_root_table_does_not_duplicate_states() {
        let result = parse(STDA);
        let mut roots: Vec<usize> = result.states.iter().map(|s| s.root).collect();
        let count = roots.len();
        roots.sort_unstable();
        roots.dedup();
        assert_eq!(roots.len(), count, "a root appears twice");
        assert_eq!(roots, vec![1, 2, 3, 4, 5, 6]);
    }

    /// A four-number line that is not a root -- a Mulliken charge table, say --
    /// must not become a state.
    #[test]
    fn a_four_number_line_that_is_not_a_root_is_skipped() {
        // Behemoth's own charge table has this shape with a symbol in it.
        assert!(parse_root_table("     2 O       5.492402      16.138449\n").is_empty());
        // And a negative or zero energy is not an excitation.
        assert!(parse_root_table("1  -1.5  300.0  0.1\n").is_empty());
        assert!(parse_root_table("1   0.0    0.0  0.0\n").is_empty());
    }
}
