//! Rotational symmetry numbers by point group.
//!
//! Sigma is the order of the group's **proper-rotation subgroup**: how many
//! indistinguishable orientations a rotation can bring the molecule into.
//! Reflections, inversions and improper rotations do not count, which is the
//! usual source of error -- C2v has four symmetry operations but sigma = 2,
//! and S4 has four but sigma = 2 as well.
//!
//! Sigma divides the rotational partition function, so an error of a factor f
//! shifts the rotational entropy by `R ln f` and every free energy with it:
//! 1.4 cal/mol/K for a factor of 2, 4.9 for a factor of 12.
//!
//! This table is a reference the user reads before typing a number; nothing
//! detects a point group automatically.

/// Every group, as `(symbol, sigma)`, grouped into the families the reference
/// table prints as rows.
struct Family {
    name: &'static str,
    members: &'static [(&'static str, f64)],
}

const FAMILIES: &[Family] = &[
    Family {
        name: "No rotation",
        members: &[
            ("C1", 1.0),
            ("Cs", 1.0),
            ("Ci", 1.0),
            ("C\u{221e}v", 1.0),
            ("Kh", 1.0),
        ],
    },
    Family {
        name: "Linear",
        members: &[("D\u{221e}h", 2.0)],
    },
    Family {
        name: "Cn",
        members: &[
            ("C2", 2.0),
            ("C3", 3.0),
            ("C4", 4.0),
            ("C5", 5.0),
            ("C6", 6.0),
            ("C7", 7.0),
            ("C8", 8.0),
        ],
    },
    Family {
        name: "Cnv",
        members: &[
            ("C2v", 2.0),
            ("C3v", 3.0),
            ("C4v", 4.0),
            ("C5v", 5.0),
            ("C6v", 6.0),
            ("C7v", 7.0),
            ("C8v", 8.0),
        ],
    },
    Family {
        name: "Cnh",
        members: &[
            ("C2h", 2.0),
            ("C3h", 3.0),
            ("C4h", 4.0),
            ("C5h", 5.0),
            ("C6h", 6.0),
            ("C7h", 7.0),
            ("C8h", 8.0),
        ],
    },
    Family {
        name: "Dn",
        members: &[
            ("D2", 4.0),
            ("D3", 6.0),
            ("D4", 8.0),
            ("D5", 10.0),
            ("D6", 12.0),
            ("D7", 14.0),
            ("D8", 16.0),
        ],
    },
    Family {
        name: "Dnh",
        members: &[
            ("D2h", 4.0),
            ("D3h", 6.0),
            ("D4h", 8.0),
            ("D5h", 10.0),
            ("D6h", 12.0),
            ("D7h", 14.0),
            ("D8h", 16.0),
        ],
    },
    Family {
        name: "Dnd",
        members: &[
            ("D2d", 4.0),
            ("D3d", 6.0),
            ("D4d", 8.0),
            ("D5d", 10.0),
            ("D6d", 12.0),
            ("D7d", 14.0),
            ("D8d", 16.0),
        ],
    },
    Family {
        // Sn's rotation subgroup is C(n/2), so sigma is n/2, not n.
        name: "Sn",
        members: &[("S4", 2.0), ("S6", 3.0), ("S8", 4.0)],
    },
    Family {
        name: "Cubic",
        members: &[
            ("T", 12.0),
            ("Th", 12.0),
            ("Td", 12.0),
            ("O", 24.0),
            ("Oh", 24.0),
            ("I", 60.0),
            ("Ih", 60.0),
        ],
    },
];

/// Sigma for a Schoenflies symbol, or `None` if the symbol is not recognised.
/// Case-insensitive, and `Cinfv` / `Dinfh` are accepted for the linear groups.
pub fn symmetry_number(symbol: &str) -> Option<f64> {
    let wanted = normalize(symbol);
    FAMILIES
        .iter()
        .flat_map(|f| f.members)
        .find(|(name, _)| normalize(name) == wanted)
        .map(|(_, sigma)| *sigma)
}

/// Lowercases and spells the infinity sign out, so `Cinfv`, `C∞v` and `CINFV`
/// all match.
fn normalize(symbol: &str) -> String {
    symbol
        .trim()
        .to_ascii_lowercase()
        .replace('\u{221e}', "inf")
}

/// The table as `(family name, members)` rows, for a UI that lays it out
/// itself rather than taking the preformatted text.
pub fn families() -> Vec<(&'static str, Vec<(&'static str, f64)>)> {
    FAMILIES
        .iter()
        .map(|f| (f.name, f.members.to_vec()))
        .collect()
}

/// The reference table as printable text, laid out one family per row.
///
/// Generated from the same data `symmetry_number` looks up, so the table a
/// user reads and the value the code would use cannot disagree.
pub fn reference_table() -> String {
    use std::fmt::Write as _;
    const LABEL_W: usize = 12;
    const CELL_W: usize = 9;

    let mut out = String::new();
    out.push_str("Rotational symmetry number by point group\n");
    out.push_str("sigma = order of the proper-rotation subgroup\n");
    out.push_str("(reflections and improper rotations do not count)\n\n");
    for family in FAMILIES {
        let _ = write!(out, "{:<LABEL_W$}", family.name);
        for (name, sigma) in family.members {
            // The infinity sign is one character but several bytes, so pad by
            // character count rather than letting {:<w$} count bytes.
            let cell = format!("{name} {sigma:.0}");
            let pad = CELL_W.saturating_sub(cell.chars().count());
            let _ = write!(out, "{cell}{:pad$}", "");
        }
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_familiar_molecules_come_out_right() {
        for (symbol, sigma) in [
            ("C1", 1.0),   // most molecules
            ("Cs", 1.0),   // HOCl
            ("C2v", 2.0),  // water
            ("D\u{221e}h", 2.0), // CO2
            ("C3v", 3.0),  // ammonia
            ("D2h", 4.0),  // ethylene
            ("D3h", 6.0),  // BF3
            ("D3d", 6.0),  // staggered ethane
            ("D4h", 8.0),  // XeF4
            ("Td", 12.0),  // methane
            ("D6h", 12.0), // benzene
            ("Oh", 24.0),  // SF6
            ("Ih", 60.0),  // C60
        ] {
            assert_eq!(symmetry_number(symbol), Some(sigma), "{symbol}");
        }
    }

    /// The five entries where the reference Fortran this table came from is
    /// wrong. Sigma counts proper rotations only: Sn's subgroup is C(n/2), and
    /// T/O/I have 12/24/60 rotations, not half that.
    #[test]
    fn the_corrected_entries_are_the_proper_rotation_counts() {
        assert_eq!(symmetry_number("S4"), Some(2.0), "not 4: rotations are E, C2");
        assert_eq!(symmetry_number("S6"), Some(3.0), "not 6");
        assert_eq!(symmetry_number("S8"), Some(4.0), "not 8");
        assert_eq!(symmetry_number("T"), Some(12.0), "not 6");
        assert_eq!(symmetry_number("O"), Some(24.0), "not 12");
        assert_eq!(symmetry_number("I"), Some(60.0), "not 30");
    }

    /// Dn, Dnh and Dnd share a rotation subgroup of order 2n.
    #[test]
    fn the_d_families_agree_at_each_order() {
        for n in 2..=8 {
            let plain = symmetry_number(&format!("D{n}")).unwrap();
            assert_eq!(plain, 2.0 * n as f64, "D{n}");
            assert_eq!(symmetry_number(&format!("D{n}h")), Some(plain));
            assert_eq!(symmetry_number(&format!("D{n}d")), Some(plain));
        }
    }

    /// Cn, Cnv and Cnh all have exactly the n rotations of Cn.
    #[test]
    fn the_c_families_agree_at_each_order() {
        for n in 2..=8 {
            assert_eq!(symmetry_number(&format!("C{n}")), Some(n as f64));
            assert_eq!(symmetry_number(&format!("C{n}v")), Some(n as f64));
            assert_eq!(symmetry_number(&format!("C{n}h")), Some(n as f64));
        }
    }

    /// C8 exists here. The Fortran this came from had a duplicated `case("C7")`
    /// where C8 was meant, so C8 fell through to its error branch.
    #[test]
    fn c8_is_present() {
        assert_eq!(symmetry_number("C8"), Some(8.0));
    }

    #[test]
    fn linear_groups_are_accepted_spelled_either_way() {
        assert_eq!(symmetry_number("Cinfv"), Some(1.0));
        assert_eq!(symmetry_number("C\u{221e}v"), Some(1.0));
        assert_eq!(symmetry_number("DINFH"), Some(2.0));
        assert_eq!(symmetry_number(" d\u{221e}h "), Some(2.0));
    }

    #[test]
    fn an_unknown_symbol_is_none_rather_than_a_guess() {
        assert_eq!(symmetry_number("C9v"), None);
        assert_eq!(symmetry_number(""), None);
        assert_eq!(symmetry_number("banana"), None);
    }

    #[test]
    fn no_symbol_is_listed_twice() {
        let mut names: Vec<String> = FAMILIES
            .iter()
            .flat_map(|f| f.members)
            .map(|(n, _)| normalize(n))
            .collect();
        let count = names.len();
        names.sort();
        names.dedup();
        assert_eq!(names.len(), count, "a duplicated symbol shadows another");
    }

    /// The printed table must contain every entry the lookup knows, or a user
    /// reading it would be missing a group the code accepts.
    #[test]
    fn the_printed_table_lists_every_group_with_its_value() {
        let text = reference_table();
        for (name, sigma) in FAMILIES.iter().flat_map(|f| f.members) {
            let cell = format!("{name} {sigma:.0}");
            assert!(text.contains(&cell), "{cell:?} missing from the table");
        }
        for family in FAMILIES {
            assert!(text.contains(family.name), "{} missing", family.name);
        }
    }

    /// Every row lines up: the infinity sign is three bytes but one column.
    #[test]
    fn the_columns_align_despite_the_infinity_sign() {
        let text = reference_table();
        let rows: Vec<&str> = text
            .lines()
            .filter(|l| FAMILIES.iter().any(|f| l.starts_with(f.name)))
            .collect();
        assert_eq!(rows.len(), FAMILIES.len());
        for row in rows {
            // Each cell starts on a 9-character boundary after the 12-wide
            // label, so a row's character count is always 12 + 9k.
            let width = row.chars().count();
            assert_eq!((width - 12) % 9, 0, "ragged row: {row:?} ({width} chars)");
        }
    }
}
