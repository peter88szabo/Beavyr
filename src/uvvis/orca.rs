//! Reads TD-DFT excited states and the absorption spectrum out of an ORCA
//! output file.
//!
//! Three independent passes, each tolerant of the others being absent: the
//! excited-state block, the electric-dipole absorption table, and the ground
//! state's spin contamination. A run that produced roots but no absorption
//! table still yields a usable list, because the roots and their orbital
//! character are worth reading on their own.

use std::path::{Path, PathBuf};

use super::types::{ev_to_nm, Excitation, ExcitedState, Spin, TddftResult, HARTREE_TO_EV};

/// The heading above the absorption table we want. Matched exactly, because
/// the velocity-gauge and spin-orbit-corrected tables have headings that
/// merely contain this one's wording.
const ABSORPTION_HEADING: &str = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS";

/// Parses an ORCA output file's TD-DFT section.
///
/// Fails only when there is no excited-state block at all; everything else
/// degrades to a `note` on the result.
pub fn parse_orca_tddft(text: &str, source: &Path) -> Result<TddftResult, String> {
    let lines: Vec<&str> = text.lines().collect();
    let mut notes = Vec::new();

    let blocks = find_state_blocks(&lines);
    let Some(&(header_line, ref method)) = blocks.first() else {
        return Err(
            "no TD-DFT excited states found: expected a line ending in \
             'EXCITED STATES' (as ORCA writes above its root list)"
                .to_string(),
        );
    };
    if blocks.len() > 1 {
        // Root numbering restarts in each block, so a second block's roots
        // could not be told apart from the first's in the absorption table.
        notes.push(format!(
            "{} further excited-state block(s) in this file were skipped; only the first \
             ('{method}') is shown.",
            blocks.len() - 1
        ));
    }

    let mut states = parse_state_block(&lines, header_line + 1)?;
    if states.is_empty() {
        return Err(format!(
            "the '{method} EXCITED STATES' block at line {} lists no states",
            header_line + 1
        ));
    }

    match parse_absorption_table(&lines) {
        Ok(rows) => {
            if rows.is_empty() {
                notes.push(
                    "no absorption table in this file: the states have no oscillator \
                     strengths, so there is nothing to plot."
                        .to_string(),
                );
            }
            for row in rows {
                if let Some(state) = states.iter_mut().find(|s| s.root == row.root) {
                    state.oscillator_strength = Some(row.oscillator_strength);
                    state.d2_au2 = row.d2;
                    state.transition_dipole_au = row.dipole;
                    if let Some(nm) = row.wavelength_nm {
                        state.wavelength_nm = nm;
                    }
                }
            }
        }
        Err(err) => notes.push(format!("absorption table skipped: {err}")),
    }

    let unrestricted = states.iter().any(|s| s.s2.is_some())
        || lines.iter().any(|l| l.trim() == "UHF SPIN CONTAMINATION");

    Ok(TddftResult {
        source: PathBuf::from(source),
        program: parse_program_version(&lines),
        method: method.clone(),
        unrestricted,
        ground_state_s2: parse_ground_state_s2(&lines),
        states,
        notes,
    })
}

/// `ORCA 6.1.0`, or just `ORCA` if the banner is missing.
fn parse_program_version(lines: &[&str]) -> String {
    for line in lines {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix("Program Version ") {
            if let Some(version) = rest.split_whitespace().next() {
                return format!("ORCA {version}");
            }
        }
    }
    "ORCA".to_string()
}

/// `Expectation value of <S**2>     :     8.756797`
fn parse_ground_state_s2(lines: &[&str]) -> Option<f64> {
    lines.iter().find_map(|line| {
        let trimmed = line.trim();
        if !trimmed.starts_with("Expectation value of <S**2>") {
            return None;
        }
        trimmed.rsplit(':').next()?.trim().parse().ok()
    })
}

/// Every `... EXCITED STATES` heading, as `(line index, method name)`. The
/// method is whatever precedes the words, e.g. `TD-DFT/TDA`; a trailing
/// `(SINGLETS)` / `(TRIPLETS)` qualifier is kept as part of it.
fn find_state_blocks<'a>(lines: &[&'a str]) -> Vec<(usize, String)> {
    lines
        .iter()
        .enumerate()
        .filter_map(|(i, line)| {
            let trimmed = line.trim();
            let (before, after) = trimmed.split_once("EXCITED STATES")?;
            let before = before.trim();
            let after = after.trim();
            // The heading is the whole line: `TD-DFT/TDA EXCITED STATES`,
            // optionally with a spin qualifier. Prose that merely mentions
            // excited states is not a heading.
            if before.is_empty() || before.contains(' ') {
                return None;
            }
            if !after.is_empty() && !(after.starts_with('(') && after.ends_with(')')) {
                return None;
            }
            let method = if after.is_empty() {
                before.to_string()
            } else {
                format!("{before} {after}")
            };
            Some((i, method))
        })
        .collect()
}

/// Reads `STATE` records and their excitation lines, starting just after the
/// block heading. Preamble prose before the first `STATE` is skipped; the
/// block ends at the first line after that which is neither.
fn parse_state_block(lines: &[&str], start: usize) -> Result<Vec<ExcitedState>, String> {
    let mut states: Vec<ExcitedState> = Vec::new();
    for (offset, line) in lines[start.min(lines.len())..].iter().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with("STATE") {
            states.push(parse_state_line(trimmed).map_err(|e| {
                format!("line {}: {e}", start + offset + 1)
            })?);
            continue;
        }
        if let Some(excitation) = parse_excitation_line(trimmed) {
            match states.last_mut() {
                Some(state) => state.excitations.push(excitation),
                // An orbital pair before any STATE record is not part of this
                // block; treat it as the block never having started.
                None => continue,
            }
            continue;
        }
        if states.is_empty() {
            // Still in the explanatory preamble ORCA prints under the heading.
            continue;
        }
        break;
    }
    Ok(states)
}

/// `STATE  1:  E=   0.071104 au      1.935 eV    15605.6 cm**-1 <S**2> =   8.870929 Mult 6`
///
/// The `<S**2>` and `Mult` tail is present only for unrestricted references,
/// and the cm⁻¹ column is missing in some versions, so both are optional.
fn parse_state_line(line: &str) -> Result<ExcitedState, String> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    let root: usize = tokens
        .get(1)
        .and_then(|t| t.trim_end_matches(':').parse().ok())
        .ok_or_else(|| format!("cannot read a root number from {line:?}"))?;

    // Each energy is the number immediately before its unit.
    let value_before = |unit: &str| -> Option<f64> {
        let i = tokens.iter().position(|t| *t == unit)?;
        tokens.get(i.checked_sub(1)?)?.parse().ok()
    };
    let energy_au = value_before("au").unwrap_or(0.0);
    let printed_ev = value_before("eV")
        .ok_or_else(|| format!("no excitation energy in eV in {line:?}"))?;
    // ORCA prints the eV column to three decimals but the hartree column to
    // six, so the hartree value is ~40x more precise. On a 600 nm band those
    // three decimals quantize the peak to a third of a nanometre, which is
    // visible once the spectrum is plotted against wavelength.
    let energy_ev = if energy_au > 0.0 {
        energy_au * HARTREE_TO_EV
    } else {
        printed_ev
    };
    let energy_cm1 = value_before("cm**-1").unwrap_or(energy_ev * 8065.543_937);

    let s2 = tokens
        .iter()
        .position(|t| *t == "<S**2>")
        .and_then(|i| tokens[i + 1..].iter().find_map(|t| t.parse::<f64>().ok()));
    let multiplicity = tokens
        .iter()
        .position(|t| *t == "Mult")
        .and_then(|i| tokens.get(i + 1)?.parse::<u32>().ok());

    Ok(ExcitedState {
        root,
        energy_au,
        energy_ev,
        energy_cm1,
        wavelength_nm: ev_to_nm(energy_ev),
        s2,
        multiplicity,
        oscillator_strength: None,
        d2_au2: None,
        transition_dipole_au: None,
        excitations: Vec::new(),
    })
}

/// `89b ->  90b  :     0.988112 (c= -0.99403830)`, or without the spin letters
/// and the coefficient for a restricted reference.
fn parse_excitation_line(line: &str) -> Option<Excitation> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 5 || tokens[1] != "->" {
        return None;
    }
    let colon = tokens.iter().position(|t| *t == ":")?;
    let (from_orbital, from_spin) = parse_orbital(tokens[0])?;
    let (to_orbital, to_spin) = parse_orbital(tokens[2])?;
    let weight: f64 = tokens.get(colon + 1)?.parse().ok()?;

    // `(c= -0.994)` -- the sign may be glued to the marker or stand alone.
    let coefficient = tokens.iter().position(|t| t.starts_with("(c=")).and_then(|i| {
        let inline = tokens[i].trim_start_matches("(c=").trim_end_matches(')');
        if !inline.is_empty() {
            inline.parse().ok()
        } else {
            tokens.get(i + 1)?.trim_end_matches(')').parse().ok()
        }
    });

    Some(Excitation {
        from_orbital,
        from_spin,
        to_orbital,
        to_spin,
        weight,
        coefficient,
    })
}

/// `89b` -> (89, Beta); `47` -> (47, Unspecified).
fn parse_orbital(token: &str) -> Option<(u32, Spin)> {
    let last = token.chars().last()?;
    match Spin::from_letter(last) {
        Some(spin) => Some((token[..token.len() - last.len_utf8()].parse().ok()?, spin)),
        None => Some((token.parse().ok()?, Spin::Unspecified)),
    }
}

/// One row of the electric-dipole absorption table.
#[derive(Debug, Clone, PartialEq)]
struct AbsorptionRow {
    root: usize,
    oscillator_strength: f64,
    wavelength_nm: Option<f64>,
    d2: Option<f64>,
    dipole: Option<[f64; 3]>,
}

/// Which numeric column of the table holds what. Columns are located by their
/// printed names, because ORCA 5 and ORCA 6 order and name them differently.
#[derive(Debug, Default, Clone, Copy, PartialEq)]
struct AbsorptionColumns {
    fosc: Option<usize>,
    wavelength: Option<usize>,
    d2: Option<usize>,
    dx: Option<usize>,
    dy: Option<usize>,
    dz: Option<usize>,
}

/// Reads the names row of the table. The first name labels the transition
/// column, which is text, so the numeric columns start after it.
fn map_absorption_columns(names_row: &str) -> AbsorptionColumns {
    let mut cols = AbsorptionColumns::default();
    for (i, name) in names_row.split_whitespace().skip(1).enumerate() {
        match name {
            n if n.starts_with("fosc") => cols.fosc = Some(i),
            "Wavelength" => cols.wavelength = Some(i),
            // D/T/P: electric dipole, transition dipole and velocity gauges
            // all print the same shape of table.
            "D2" | "T2" | "P2" => cols.d2 = Some(i),
            "DX" | "TX" | "PX" => cols.dx = Some(i),
            "DY" | "TY" | "PY" => cols.dy = Some(i),
            "DZ" | "TZ" | "PZ" => cols.dz = Some(i),
            _ => {}
        }
    }
    cols
}

/// Splits a data row into its root number and its numeric columns.
///
/// ORCA 6 labels the row `0-6A  ->  1-6A`; ORCA 5 just numbers it. Either way
/// the root comes from the label, never from the row's position, so a table
/// that omits a root does not shift every later assignment.
fn split_absorption_row(line: &str) -> Option<(usize, Vec<f64>)> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    let (root_token, rest) = match tokens.iter().position(|t| *t == "->") {
        Some(arrow) => (*tokens.get(arrow + 1)?, &tokens[arrow + 2..]),
        None => (*tokens.first()?, &tokens[1..]),
    };
    // `1-6A` is root 1 of multiplicity 6, symmetry A.
    let root: usize = root_token.split('-').next()?.parse().ok()?;
    let numbers: Vec<f64> = rest.iter().filter_map(|t| t.parse().ok()).collect();
    if numbers.len() != rest.len() {
        return None;
    }
    Some((root, numbers))
}

/// Finds and reads the electric-dipole absorption table. An absent table is
/// `Ok(vec![])`; a table whose columns cannot be understood is an `Err`, since
/// silently guessing them would put wrong intensities on the plot.
fn parse_absorption_table(lines: &[&str]) -> Result<Vec<AbsorptionRow>, String> {
    let Some(heading) = lines.iter().position(|l| l.trim() == ABSORPTION_HEADING) else {
        return Ok(Vec::new());
    };
    // The names row is the first line under the heading that names a column
    // we recognise; between them sits a rule of dashes.
    let names_row = lines[heading + 1..]
        .iter()
        .take(4)
        .find(|l| l.contains("Wavelength") || l.contains("fosc"))
        .ok_or_else(|| format!("no column names under the heading at line {}", heading + 1))?;
    let cols = map_absorption_columns(names_row);
    let fosc_col = cols.fosc.ok_or_else(|| {
        format!("no oscillator-strength column in {:?}", names_row.trim())
    })?;

    let mut rows = Vec::new();
    let mut started = false;
    for line in &lines[heading + 1..] {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            // Blank lines separate the table from what follows, but also sit
            // between the heading and the header rows in some versions.
            if started {
                break;
            }
            continue;
        }
        let Some((root, numbers)) = split_absorption_row(trimmed) else {
            if started {
                break;
            }
            continue;
        };
        let at = |col: Option<usize>| col.and_then(|c| numbers.get(c).copied());
        let Some(oscillator_strength) = numbers.get(fosc_col).copied() else {
            continue;
        };
        started = true;
        let dipole = match (at(cols.dx), at(cols.dy), at(cols.dz)) {
            (Some(x), Some(y), Some(z)) => Some([x, y, z]),
            _ => None,
        };
        rows.push(AbsorptionRow {
            root,
            oscillator_strength,
            wavelength_nm: at(cols.wavelength),
            d2: at(cols.d2),
            dipole,
        });
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn example() -> String {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("AbsorptionSpectrum")
            .join("Exc_TPSSh-D4-def2-TZVP_SCMwater_16.out");
        std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
    }

    fn parsed() -> TddftResult {
        parse_orca_tddft(&example(), Path::new("example.out")).expect("example must parse")
    }

    // ---- against the real captured ORCA 6.1 output ----

    #[test]
    fn all_eighty_roots_are_read() {
        let r = parsed();
        assert_eq!(r.states.len(), 80);
        assert_eq!(r.states[0].root, 1);
        assert_eq!(r.states[79].root, 80);
        assert_eq!(r.program, "ORCA 6.1.0");
        assert_eq!(r.method, "TD-DFT/TDA");
    }

    #[test]
    fn the_first_and_last_state_energies_match_the_file() {
        let r = parsed();
        let first = &r.states[0];
        assert!((first.energy_au - 0.071104).abs() < 1e-9);
        assert!((first.energy_cm1 - 15605.6).abs() < 1e-6);
        let last = &r.states[79];
        assert!((last.energy_cm1 - 41260.9).abs() < 1e-6);
        assert!((last.energy_ev - 5.115694).abs() < 1e-4, "{}", last.energy_ev);
    }

    /// The eV column is rounded to three decimals; the hartree column on the
    /// same line is not. Taking the energy from hartree agrees with the six
    /// decimals the absorption table prints independently, which the rounded
    /// column cannot.
    #[test]
    fn transition_energies_keep_the_precision_of_the_hartree_column() {
        let r = parsed();
        // Absorption table, root 1: 1.934849 eV. The STATE line says "1.935".
        let ev = r.states[0].energy_ev;
        assert!((ev - 1.934849).abs() < 1e-4, "got {ev}");
        assert!((ev - 1.935).abs() > 1e-6, "the rounded column was used: {ev}");
    }

    #[test]
    fn spin_contamination_is_read_for_this_unrestricted_run() {
        let r = parsed();
        assert!(r.unrestricted);
        assert!((r.ground_state_s2.unwrap() - 8.756797).abs() < 1e-9);
        assert!((r.states[0].s2.unwrap() - 8.870929).abs() < 1e-9);
        assert_eq!(r.states[0].multiplicity, Some(6));
        // Root 19 is the strongly contaminated one.
        assert!((r.states[18].s2.unwrap() - 10.425643).abs() < 1e-9);
    }

    #[test]
    fn a_single_configuration_state_lists_exactly_its_one_excitation() {
        let r = parsed();
        let s = &r.states[0];
        assert_eq!(s.excitations.len(), 1);
        let e = &s.excitations[0];
        assert_eq!(e.from_orbital, 89);
        assert_eq!(e.to_orbital, 90);
        assert_eq!(e.from_spin, Spin::Beta);
        assert_eq!(e.to_spin, Spin::Beta);
        assert!((e.weight - 0.988112).abs() < 1e-9);
        assert!((e.coefficient.unwrap() + 0.99403830).abs() < 1e-9);
        assert_eq!(e.label(), "89b \u{2192} 90b");
    }

    #[test]
    fn a_multiconfigurational_state_keeps_every_contribution_in_file_order() {
        let r = parsed();
        // STATE 10 lists five contributions.
        let s = &r.states[9];
        assert_eq!(s.root, 10);
        let pairs: Vec<(u32, u32)> = s
            .excitations
            .iter()
            .map(|e| (e.from_orbital, e.to_orbital))
            .collect();
        assert_eq!(pairs, vec![(84, 93), (85, 94), (87, 90), (87, 92), (88, 91)]);
        // ...and ranks them by weight for display.
        assert_eq!(s.dominant_excitation().unwrap().from_orbital, 87);
        assert!((s.dominant_excitation().unwrap().weight - 0.538456).abs() < 1e-9);
    }

    /// The last state's contributions mix alpha and beta, which is only
    /// possible for an unrestricted reference -- a good check that the spin
    /// letter is read per orbital and not assumed per state.
    #[test]
    fn alpha_and_beta_contributions_are_told_apart() {
        let r = parsed();
        let s = &r.states[79];
        assert_eq!(s.excitations[0].from_spin, Spin::Alpha);
        assert_eq!(s.excitations[0].from_orbital, 87);
        assert_eq!(s.excitations[2].from_spin, Spin::Beta);
        assert_eq!(s.excitations[2].from_orbital, 78);
    }

    #[test]
    fn the_absorption_table_supplies_intensities_and_dipoles() {
        let r = parsed();
        assert!(r.has_intensities());
        let s = &r.states[0];
        assert!((s.oscillator_strength.unwrap() - 0.002236399).abs() < 1e-12);
        assert!((s.wavelength_nm - 640.8).abs() < 1e-9);
        assert!((s.d2_au2.unwrap() - 0.04718).abs() < 1e-9);
        let d = s.transition_dipole_au.unwrap();
        assert!((d[0] - 0.00093).abs() < 1e-9);
        assert!((d[1] + 0.11230).abs() < 1e-9);
        assert!((d[2] - 0.18592).abs() < 1e-9);
    }

    /// Rows are matched to states by root number, so the last row must land on
    /// the last state and nowhere else.
    #[test]
    fn every_root_gets_its_own_row() {
        let r = parsed();
        assert!(r.states.iter().all(|s| s.oscillator_strength.is_some()));
        let last = &r.states[79];
        assert!((last.oscillator_strength.unwrap() - 0.000603496).abs() < 1e-12);
        assert!((last.wavelength_nm - 242.4).abs() < 1e-9);
    }

    /// The velocity-dipole table follows the electric-dipole one and has the
    /// same shape; reading it by mistake would silently replace every
    /// intensity with a different quantity.
    #[test]
    fn the_velocity_gauge_table_is_not_mistaken_for_the_electric_dipole_one() {
        let r = parsed();
        // Velocity-gauge fosc for root 1 is 0.001221997, electric is 0.002236399.
        assert!((r.states[0].oscillator_strength.unwrap() - 0.002236399).abs() < 1e-12);
        // ...and the velocity table's P2 for root 1 is 0.00013, not 0.04718.
        assert!((r.states[0].d2_au2.unwrap() - 0.04718).abs() < 1e-9);
    }

    #[test]
    fn the_example_parses_without_complaint() {
        assert!(parsed().notes.is_empty(), "{:?}", parsed().notes);
    }

    // ---- synthetic cases ----

    const MINIMAL: &str = "\
-------------------------
TD-DFT/TDA EXCITED STATES
-------------------------

the weight of the individual excitations are printed if larger than 1.0e-02

STATE  1:  E=   0.100000 au      2.000 eV    16000.0 cm**-1
    45 ->  47  :     0.900000

STATE  2:  E=   0.200000 au      4.000 eV    32000.0 cm**-1
    46 ->  47  :     0.800000

Storing amplitudes in GBW file ...
";

    #[test]
    fn a_restricted_run_parses_without_spin_letters_or_coefficients() {
        let r = parse_orca_tddft(MINIMAL, Path::new("t.out")).unwrap();
        assert_eq!(r.states.len(), 2);
        assert!(!r.unrestricted);
        assert!(r.ground_state_s2.is_none());
        assert!(r.states[0].s2.is_none());
        let e = &r.states[0].excitations[0];
        assert_eq!(e.from_spin, Spin::Unspecified);
        assert!(e.coefficient.is_none());
        assert_eq!(e.label(), "45 \u{2192} 47");
    }

    /// Roots with no absorption table are still worth listing -- that is the
    /// whole reason the table is parsed separately.
    #[test]
    fn states_without_an_absorption_table_still_parse_and_say_so() {
        let r = parse_orca_tddft(MINIMAL, Path::new("t.out")).unwrap();
        assert!(!r.has_intensities());
        assert!(r.states.iter().all(|s| s.oscillator_strength.is_none()));
        assert!(
            r.notes.iter().any(|n| n.contains("no absorption table")),
            "{:?}",
            r.notes
        );
        // Wavelength falls back to the transition energy, taken from the
        // hartree column: 0.1 Eh = 2.72114 eV = 455.63 nm.
        assert!((r.states[0].energy_ev - 2.721_138_6).abs() < 1e-6);
        assert!((r.states[0].wavelength_nm - 455.633).abs() < 1e-3, "{}", r.states[0].wavelength_nm);
    }

    #[test]
    fn the_block_ends_at_the_first_unrelated_line() {
        // "Storing amplitudes..." must not be read as a third state.
        let r = parse_orca_tddft(MINIMAL, Path::new("t.out")).unwrap();
        assert_eq!(r.states.len(), 2);
    }

    #[test]
    fn a_file_with_no_tddft_section_is_an_error_naming_what_was_missing() {
        let err = parse_orca_tddft("SCF CONVERGED\nTOTAL ENERGY -76.4\n", Path::new("t.out"))
            .unwrap_err();
        assert!(err.contains("EXCITED STATES"), "{err}");
    }

    #[test]
    fn an_empty_state_block_is_an_error_rather_than_an_empty_list() {
        let text = "TD-DFT/TDA EXCITED STATES\n\nStoring amplitudes\n";
        assert!(parse_orca_tddft(text, Path::new("t.out")).is_err());
    }

    #[test]
    fn a_second_state_block_is_skipped_with_a_note() {
        let text = format!(
            "{MINIMAL}\nTD-DFT/TDA EXCITED STATES (TRIPLETS)\n\n\
             STATE  1:  E=   0.050000 au      1.000 eV    8000.0 cm**-1\n    45 ->  47  :  0.9\n"
        );
        let r = parse_orca_tddft(&text, Path::new("t.out")).unwrap();
        assert_eq!(r.states.len(), 2, "only the first block's roots");
        assert!(r.notes.iter().any(|n| n.contains("skipped")), "{:?}", r.notes);
    }

    #[test]
    fn a_heading_mentioned_in_prose_is_not_taken_for_a_block() {
        let text = "we now compute the EXCITED STATES of the system\nSTATE  1:  E= 0.1 au 2.0 eV\n";
        assert!(parse_orca_tddft(text, Path::new("t.out")).is_err());
    }

    // ---- column mapping ----

    #[test]
    fn orca_6_columns_are_located_by_name() {
        let cols = map_absorption_columns(
            "     Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ",
        );
        assert_eq!(cols.fosc, Some(3));
        assert_eq!(cols.wavelength, Some(2));
        assert_eq!(cols.d2, Some(4));
        assert_eq!(cols.dx, Some(5));
        assert_eq!(cols.dz, Some(7));
    }

    /// ORCA 5 names and orders the same columns differently; mapping by name
    /// is what lets one parser read both.
    #[test]
    fn orca_5_columns_are_located_by_name_too() {
        let cols = map_absorption_columns(
            "       State   Energy  Wavelength   fosc         T2         TX    TY    TZ",
        );
        assert_eq!(cols.fosc, Some(2));
        assert_eq!(cols.wavelength, Some(1));
        assert_eq!(cols.d2, Some(3));
        assert_eq!(cols.dx, Some(4));
        assert_eq!(cols.dz, Some(6));
    }

    #[test]
    fn an_orca_6_row_is_split_into_its_root_and_numbers() {
        let (root, numbers) = split_absorption_row(
            "  0-6A  ->  1-6A    1.934849   15605.6   640.8   0.002236399   0.04718   0.00093  -0.11230   0.18592",
        )
        .unwrap();
        assert_eq!(root, 1);
        assert_eq!(numbers.len(), 8);
        assert!((numbers[2] - 640.8).abs() < 1e-9);
    }

    #[test]
    fn an_orca_5_row_is_split_into_its_root_and_numbers() {
        let (root, numbers) =
            split_absorption_row("   7   15605.6    640.8   0.002236399   0.04718   0.00093  -0.11230   0.18592")
                .unwrap();
        assert_eq!(root, 7);
        assert_eq!(numbers.len(), 7);
    }

    #[test]
    fn a_row_with_non_numeric_columns_is_rejected_rather_than_partly_read() {
        assert!(split_absorption_row("  0-6A  ->  1-6A    1.93   n/a   640.8").is_none());
        assert!(split_absorption_row("").is_none());
    }

    /// A table whose columns cannot be understood must not silently produce
    /// plausible-looking intensities.
    #[test]
    fn an_unreadable_table_is_reported_not_guessed() {
        let text = format!(
            "{MINIMAL}\n{ABSORPTION_HEADING}\n\
             ------\n     Transition   Mystery   Other\n------\n  0-1A -> 1-1A   1.0   2.0\n\n"
        );
        let r = parse_orca_tddft(&text, Path::new("t.out")).unwrap();
        assert!(!r.has_intensities());
        assert!(
            r.notes.iter().any(|n| n.contains("absorption table skipped")),
            "{:?}",
            r.notes
        );
    }

    #[test]
    fn a_table_missing_a_root_does_not_shift_the_others() {
        // Root 1 has no row; root 2's row must still land on root 2.
        let text = format!(
            "{MINIMAL}\n{ABSORPTION_HEADING}\n\
             -----\n     Transition      Energy  Wavelength fosc(D2)      D2        DX        DY        DZ\n\
                                  (eV)    (nm)                 (au**2)    (au)      (au)      (au)\n\
             -----\n  0-1A  ->  2-1A    4.000000   310.0   0.500000   1.0   0.1   0.2   0.3\n\n"
        );
        let r = parse_orca_tddft(&text, Path::new("t.out")).unwrap();
        assert!(r.states[0].oscillator_strength.is_none(), "root 1 had no row");
        assert_eq!(r.states[1].oscillator_strength, Some(0.5));
        assert!((r.states[1].wavelength_nm - 310.0).abs() < 1e-9);
    }
}
