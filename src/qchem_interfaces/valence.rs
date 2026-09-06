//! Electronic-state validation for xTB runs.
//!
//! xTB performs none of its own: an impossible spin state (verified by
//! hand, e.g. `--uhf 1` on hydrogen peroxide, an 18-electron closed-shell
//! molecule) still exits 0 with "normal termination", silently ignoring the
//! request rather than erroring. This module is what actually catches that.

use crate::bond_order::estimate_bond_order;
use crate::molecule::{atomic_number, Molecule};

/// A physically impossible charge/multiplicity combination. Blocks the run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ElectronicStateError(pub String);

impl std::fmt::Display for ElectronicStateError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// A structurally unusual atom -- an apparent radical centre or an
/// over-valent atom -- that is chemically *possible* but worth a second
/// look. Never blocks the run.
#[derive(Debug, Clone, PartialEq)]
pub struct ValenceWarning {
    pub atom_index: usize,
    pub message: String,
}

/// Typical valence (bond count) for common organic elements. An element
/// outside this table contributes no warning: silence is safer than a false
/// positive for chemistry the heuristic doesn't understand.
fn typical_valence(symbol: &str) -> Option<u32> {
    Some(match symbol {
        "H" | "F" | "Cl" | "Br" | "I" => 1,
        "O" => 2,
        "N" => 3,
        "C" => 4,
        "P" => 3,
        "S" => 2,
        _ => return None,
    })
}

/// Validate `charge`/`multiplicity` against the electrons implied by `mol`'s
/// own drawn structure (element identity plus its own bond graph -- Beavyr
/// carries no separate per-atom formal-charge field to consult).
///
/// On success, returns the xTB `--uhf` value (the number of unpaired
/// electrons) together with any non-blocking valence warnings.
pub fn validate_electronic_state(
    mol: &Molecule,
    charge: i32,
    multiplicity: i32,
) -> Result<(i32, Vec<ValenceWarning>), ElectronicStateError> {
    if multiplicity < 1 {
        return Err(ElectronicStateError(
            "Spin multiplicity must be at least 1.".to_string(),
        ));
    }
    if mol.atoms.is_empty() {
        return Err(ElectronicStateError(
            "There is no structure to optimize.".to_string(),
        ));
    }

    let mut electrons: i64 = 0;
    for symbol in &mol.atoms {
        let Some(z) = atomic_number(symbol) else {
            return Err(ElectronicStateError(format!(
                "xTB cannot use the non-element atom label {symbol:?}."
            )));
        };
        electrons += z as i64;
    }
    electrons -= charge as i64;

    if electrons <= 0 {
        return Err(ElectronicStateError(
            "Charge leaves no electrons in the selected structure.".to_string(),
        ));
    }

    let unpaired = multiplicity as i64 - 1;
    if unpaired > electrons {
        return Err(ElectronicStateError(
            "Multiplicity requires more unpaired electrons than are available.".to_string(),
        ));
    }
    if (electrons + multiplicity as i64) % 2 != 1 {
        return Err(ElectronicStateError(format!(
            "Multiplicity {multiplicity} is incompatible with {electrons} electrons; \
             choose a multiplicity with the opposite parity."
        )));
    }

    // Warning pass: how much of a "valence deficit" (typical - actual bond
    // ORDER, floored at 0) is present overall, versus how much the requested
    // multiplicity's unpaired-electron budget can explain. A deficit beyond
    // that budget is reported per atom, so the user knows which one to look
    // at, rather than as one vague aggregate message.
    //
    // Bond ORDER, not raw connection count, is essential here: `mol.bonds`
    // is a plain connectivity list with no notion of double/triple/aromatic
    // character, so an aromatic ring carbon (2 ring neighbours + 1
    // substituent = 3 connections) or an alkene carbon (=CH2: 1 double-bond
    // neighbour + 2 H = 3 connections) would otherwise look one bond short
    // of the typical 4 on every single ring or double bond in the
    // structure -- exactly the false positive this heuristic existed to
    // avoid making up in the first place.
    let mut bond_count = vec![0u32; mol.atoms.len()];
    let mut bond_order_sum = vec![0.0f64; mol.atoms.len()];
    for &(i, j, distance) in &mol.bonds {
        if i < bond_count.len() {
            bond_count[i] += 1;
        }
        if j < bond_count.len() {
            bond_count[j] += 1;
        }
        if i < bond_order_sum.len() && j < bond_order_sum.len() {
            let order = estimate_bond_order(&mol.atoms[i], &mol.atoms[j], distance);
            bond_order_sum[i] += order;
            bond_order_sum[j] += order;
        }
    }

    let mut warnings = Vec::new();
    let mut explained = unpaired.max(0) as f64;
    const EPSILON: f64 = 1.0e-6;
    for (idx, symbol) in mol.atoms.iter().enumerate() {
        let Some(typical) = typical_valence(symbol) else {
            continue;
        };
        let deficit = (typical as f64 - bond_order_sum[idx]).max(0.0);
        if deficit < EPSILON {
            continue;
        }
        if deficit <= explained + EPSILON {
            // A real, expected radical/lone-pair site given the multiplicity
            // the user asked for -- not worth flagging.
            explained -= deficit;
            continue;
        }
        warnings.push(ValenceWarning {
            atom_index: idx,
            message: format!(
                "{symbol}{} has {} bond(s) (estimated valence {:.1}); {typical} is \
                 typical for {symbol}. If this is a radical centre, increase the \
                 multiplicity; otherwise check the structure.",
                idx + 1,
                bond_count[idx],
                bond_order_sum[idx],
            ),
        });
    }

    Ok((unpaired as i32, warnings))
}


#[cfg(test)]
mod tests {
    use super::*;

    /// Every bond gets a realistic *single*-bond distance (the sum of its two
    /// atoms' own covalent radii, i.e. exactly ratio 1.0) by construction --
    /// deliberately, so an ordinary fixture reads as ordinary single bonds
    /// under `estimate_bond_order` regardless of which elements it names.
    /// `mol_with_orders` below is for a test that needs a double/triple/
    /// aromatic bond on purpose.
    fn mol(atoms: &[&str], bonds: &[(usize, usize)]) -> Molecule {
        let dated_bonds: Vec<(usize, usize, f32)> = bonds
            .iter()
            .map(|&(i, j)| {
                // A real single bond, by the same definition the estimator
                // uses -- not the inflated bond-perception radii, which would
                // make this fixture's "single" bond read as a partial one.
                let d = crate::bond_order::single_bond_length(atoms[i], atoms[j]);
                (i, j, d)
            })
            .collect();
        Molecule {
            atoms: atoms.iter().map(|s| s.to_string()).collect(),
            pos: vec![Default::default(); atoms.len()],
            bonds: dated_bonds,
            hydrogen_bonds: vec![],
        }
    }

    /// Like `mol`, but each bond's distance is given explicitly, as a
    /// fraction of the single-bond length -- for a test that specifically
    /// wants a double/triple/aromatic bond.
    fn mol_with_ratios(atoms: &[&str], bonds: &[(usize, usize, f32)]) -> Molecule {
        let dated_bonds: Vec<(usize, usize, f32)> = bonds
            .iter()
            .map(|&(i, j, ratio)| {
                let single = crate::bond_order::single_bond_length(atoms[i], atoms[j]);
                (i, j, single * ratio)
            })
            .collect();
        Molecule {
            atoms: atoms.iter().map(|s| s.to_string()).collect(),
            pos: vec![Default::default(); atoms.len()],
            bonds: dated_bonds,
            hydrogen_bonds: vec![],
        }
    }

    #[test]
    fn water_singlet_neutral_is_valid_with_no_warnings() {
        let m = mol(&["O", "H", "H"], &[(0, 1), (0, 2)]);
        let (uhf, warnings) = validate_electronic_state(&m, 0, 1).unwrap();
        assert_eq!(uhf, 0);
        assert!(warnings.is_empty());
    }

    #[test]
    fn multiplicity_below_one_is_rejected() {
        let m = mol(&["O", "H", "H"], &[(0, 1), (0, 2)]);
        assert!(validate_electronic_state(&m, 0, 0).is_err());
    }

    #[test]
    fn empty_structure_is_rejected() {
        let m = mol(&[], &[]);
        assert!(validate_electronic_state(&m, 0, 1).is_err());
    }

    #[test]
    fn unknown_element_is_rejected() {
        let m = mol(&["Xx"], &[]);
        let err = validate_electronic_state(&m, 0, 1).unwrap_err();
        assert!(err.0.contains("Xx"));
    }

    /// The exact regression case from the design spec: hydrogen peroxide is
    /// an 18-electron closed-shell molecule, so multiplicity 2 (one unpaired
    /// electron) is parity-impossible -- and xTB itself accepts it silently.
    #[test]
    fn h2o2_rejects_an_impossible_doublet() {
        let m = mol(
            &["O", "O", "H", "H"],
            &[(0, 1), (0, 2), (1, 3)],
        );
        assert!(validate_electronic_state(&m, 0, 1).is_ok());
        let err = validate_electronic_state(&m, 0, 2).unwrap_err();
        assert!(err.0.contains("parity") || err.0.contains("incompatible"));
    }

    #[test]
    fn charge_leaving_no_electrons_is_rejected() {
        let m = mol(&["H"], &[]);
        // Hydrogen has 1 electron; a +1 charge leaves none.
        assert!(validate_electronic_state(&m, 1, 1).is_err());
    }

    #[test]
    fn too_high_a_multiplicity_for_available_electrons_is_rejected() {
        let m = mol(&["H"], &[]);
        // 1 electron total, so at most 1 unpaired electron (multiplicity 2).
        assert!(validate_electronic_state(&m, 0, 4).is_err());
    }

    /// A radical with an ODD electron count (like an isolated methyl
    /// radical, 9 electrons) can never legally be a singlet at all -- that
    /// case is caught outright by the hard parity check, before the warning
    /// heuristic gets a say. To test the warning path itself, this uses a
    /// synthetic even-electron "singlet carbene"-like fragment (a 2-bonded
    /// carbon) where multiplicity 1 is parity-legal but chemically unusual.
    #[test]
    fn a_valence_deficit_is_flagged_only_when_the_multiplicity_does_not_explain_it() {
        // C with only 2 bonds (typical 4): 8 electrons total, so multiplicity
        // 1 is parity-legal, but a deficit of 2 is not explained by 0
        // unpaired electrons.
        let m = mol(&["C", "H", "H"], &[(0, 1), (0, 2)]);
        let (_uhf, warnings) = validate_electronic_state(&m, 0, 1).unwrap();
        assert_eq!(warnings.len(), 1);
        assert_eq!(warnings[0].atom_index, 0);

        // At multiplicity 3 (2 unpaired electrons, still parity-legal for 8
        // electrons), the budget exactly explains the deficit of 2.
        let (_uhf, warnings) = validate_electronic_state(&m, 0, 3).unwrap();
        assert!(warnings.is_empty(), "{warnings:?}");
    }

    /// An isolated methyl radical is a genuine, odd-electron radical: it
    /// cannot legally be a singlet at all (the hard parity check rejects it
    /// outright), and as the correct doublet it needs no warning because its
    /// one dangling bond is exactly what one unpaired electron explains.
    #[test]
    fn a_genuine_radical_is_hard_blocked_as_a_singlet_not_merely_warned() {
        let m = mol(&["C", "H", "H", "H"], &[(0, 1), (0, 2), (0, 3)]);
        let err = validate_electronic_state(&m, 0, 1).unwrap_err();
        assert!(err.0.contains("parity") || err.0.contains("incompatible"));

        let (_uhf, warnings) = validate_electronic_state(&m, 0, 2).unwrap();
        assert!(warnings.is_empty(), "{warnings:?}");
    }

    #[test]
    fn estimate_bond_order_matches_the_reference_ratios() {
        // Real ratios for a C-C bond (radius sum 1.50 A), taken from the
        // values documented in `molecule_builder::attach`: single 1.54,
        // aromatic 1.40, double 1.34, triple 1.20.
        assert_eq!(estimate_bond_order("C", "C", 1.54), 1.0);
        assert_eq!(estimate_bond_order("C", "C", 1.40), 1.5);
        assert_eq!(estimate_bond_order("C", "C", 1.34), 2.0);
        assert_eq!(estimate_bond_order("C", "C", 1.20), 3.0);
    }

    /// The exact bug the user reported: opening `examples/phenyl-OCH3.xyz`
    /// (a plain, closed-shell aromatic ring) showed a valence warning on
    /// every ring carbon, because each one has only 3 raw connections (2
    /// ring neighbours + 1 substituent) and the old check compared that
    /// count directly against a typical valence of 4 with no notion that
    /// two of those three connections are aromatic (order 1.5 each, not 1).
    #[test]
    fn an_ordinary_aromatic_ring_produces_no_warnings() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("phenyl-OCH3.xyz");
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("cannot read example {}: {e}", path.display()));
        let (_n, atoms, coords) = crate::molecule::parse_xyz_angstrom(&text);
        let pos: Vec<bevy::prelude::Vec3> = coords
            .into_iter()
            .map(|v| bevy::prelude::Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let mut m = Molecule {
            atoms,
            pos,
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        m.recompute_bonds(2.0, 3.0);

        let (_uhf, warnings) = validate_electronic_state(&m, 0, 1).unwrap();
        assert!(
            warnings.is_empty(),
            "an ordinary aromatic ring must not be flagged as under-valent: {warnings:?}"
        );
    }

    /// A plain alkene (ethylene-like: H2C=CH2) has the same shape of false
    /// positive as the aromatic case, on a plain double bond rather than a
    /// ring -- both carbons have only 3 raw connections, and both are
    /// perfectly ordinary.
    #[test]
    fn an_ordinary_double_bond_produces_no_warnings() {
        let m = mol_with_ratios(
            &["C", "C", "H", "H", "H", "H"],
            &[
                (0, 1, 0.893), // C=C, double
                (0, 2, 1.0),
                (0, 3, 1.0),
                (1, 4, 1.0),
                (1, 5, 1.0),
            ],
        );
        let (_uhf, warnings) = validate_electronic_state(&m, 0, 1).unwrap();
        assert!(warnings.is_empty(), "{warnings:?}");
    }

    /// A genuine phenoxy-methoxy radical example already in the repo (7 C /
    /// 7 H / 1 O, 57 electrons -- odd). As an odd-electron species, running
    /// it as a singlet is not merely unusual, it is impossible: the hard
    /// parity check must reject it outright. As the correct doublet, its one
    /// radical oxygen (one bond where two is typical) is exactly what the
    /// one unpaired electron explains, so no warning should fire either.
    #[test]
    fn the_repo_radical_example_is_hard_blocked_as_a_singlet_and_clean_as_a_doublet() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("phenyl-OCH3_radical.xyz");
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("cannot read example {}: {e}", path.display()));
        let (_n, atoms, coords) = crate::molecule::parse_xyz_angstrom(&text);
        let pos: Vec<bevy::prelude::Vec3> = coords
            .into_iter()
            .map(|v| bevy::prelude::Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let mut m = Molecule {
            atoms,
            pos,
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        m.recompute_bonds(2.0, 3.0);

        let err = validate_electronic_state(&m, 0, 1).unwrap_err();
        assert!(
            err.0.contains("parity") || err.0.contains("incompatible"),
            "an odd-electron radical must be hard-blocked as a singlet: {err}"
        );

        let (_uhf, warnings) = validate_electronic_state(&m, 0, 2).unwrap();
        assert!(
            warnings.is_empty(),
            "the same structure as a doublet should explain its own radical: {warnings:?}"
        );
    }
}
