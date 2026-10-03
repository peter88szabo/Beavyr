//! The level of theory an ORCA calculation runs at.
//!
//! A catalogue of its own, not a translation of anyone else's. There is a
//! concrete reason beyond tidiness: ORCA spells the Minnesota hybrid `M062X`
//! and rejects `M06-2X` with an input error, while Behemoth requires `m06-2x`
//! and does not know `m062x`. Two codes, the same functional, incompatible
//! spellings.
//!
//! Every keyword below was run against the installed ORCA (6.1.1) on H2 and
//! kept only if the job reached `ORCA TERMINATED NORMALLY`. That test dropped
//! `M06-2X`, and it also dropped `SCAN`, `PBEsol`, `HCTH407`, `M05-2X`, `MN15`
//! and `BMK`, none of which this build knows under any spelling tried. They
//! are absent rather than offered-and-broken.
//!
//! The double hybrids are here because the same test showed they work -- but
//! only when an auxiliary basis is supplied, which is why
//! [`auxiliary_basis`] exists.

/// Which family a functional belongs to, so a list of this size can be read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Family {
    /// Hartree-Fock, which is not a functional at all.
    HartreeFock,
    LdaGga,
    MetaGga,
    Hybrid,
    RangeSeparated,
    /// Carries an MP2-like correlation term, so it needs an auxiliary basis.
    DoubleHybrid,
    /// Brings its own basis set and dispersion correction.
    Composite,
}

impl Family {
    /// The heading this family sits under in the dropdown.
    pub fn heading(self) -> &'static str {
        match self {
            Family::HartreeFock => "wavefunction",
            Family::LdaGga => "LDA / GGA",
            Family::MetaGga => "meta-GGA",
            Family::Hybrid => "hybrid",
            Family::RangeSeparated => "range-separated hybrid",
            Family::DoubleHybrid => "double hybrid",
            Family::Composite => "composite (\u{201C}3c\u{201D})",
        }
    }

    /// Every family, in the order the dropdown lists them.
    pub const ALL: [Family; 7] = [
        Family::HartreeFock,
        Family::LdaGga,
        Family::MetaGga,
        Family::Hybrid,
        Family::RangeSeparated,
        Family::DoubleHybrid,
        Family::Composite,
    ];
}

/// A functional as the panel offers it.
pub struct Functional {
    /// What the dropdown shows.
    pub label: &'static str,
    /// The keyword written onto ORCA's `!` line. ORCA is case-insensitive, but
    /// the catalogue stores the spelling its own manual uses.
    pub keyword: &'static str,
    pub family: Family,
    /// What a composite brings with it, or `None` for everything else.
    pub composite: Option<&'static str>,
}

/// Every functional offered, grouped by family.
///
/// All of these were run against ORCA 6.1.1 and terminated normally.
pub const FUNCTIONALS: [Functional; 48] = [
    Functional { label: "HF", keyword: "HF", family: Family::HartreeFock, composite: None },

    // LDA and GGA.
    Functional { label: "LDA", keyword: "LDA", family: Family::LdaGga, composite: None },
    Functional { label: "BLYP", keyword: "BLYP", family: Family::LdaGga, composite: None },
    Functional { label: "BP86", keyword: "BP86", family: Family::LdaGga, composite: None },
    Functional { label: "BPW91", keyword: "BPW91", family: Family::LdaGga, composite: None },
    Functional { label: "OLYP", keyword: "OLYP", family: Family::LdaGga, composite: None },
    Functional { label: "PBE", keyword: "PBE", family: Family::LdaGga, composite: None },
    Functional { label: "PW91", keyword: "PW91", family: Family::LdaGga, composite: None },
    Functional { label: "revPBE", keyword: "revPBE", family: Family::LdaGga, composite: None },
    Functional { label: "RPBE", keyword: "RPBE", family: Family::LdaGga, composite: None },
    Functional { label: "B97-D3", keyword: "B97-D3", family: Family::LdaGga, composite: None },

    // meta-GGA. `SCAN` is deliberately absent: this build rejects it.
    Functional { label: "M06-L", keyword: "M06L", family: Family::MetaGga, composite: None },
    Functional { label: "TPSS", keyword: "TPSS", family: Family::MetaGga, composite: None },
    Functional { label: "revTPSS", keyword: "revTPSS", family: Family::MetaGga, composite: None },
    Functional { label: "r2SCAN", keyword: "r2SCAN", family: Family::MetaGga, composite: None },
    Functional { label: "B97M-V", keyword: "B97M-V", family: Family::MetaGga, composite: None },

    // Global hybrids.
    Functional { label: "B3LYP", keyword: "B3LYP", family: Family::Hybrid, composite: None },
    // ORCA's VWN5 variant, which is what most other codes call plain B3LYP.
    Functional { label: "B3LYP(VWN5)", keyword: "B3LYP/G", family: Family::Hybrid, composite: None },
    Functional { label: "B3P86", keyword: "B3P86", family: Family::Hybrid, composite: None },
    Functional { label: "B3PW91", keyword: "B3PW91", family: Family::Hybrid, composite: None },
    Functional { label: "BHandHLYP", keyword: "BHANDHLYP", family: Family::Hybrid, composite: None },
    // ORCA's spelling. `M06-2X` is an input error here; see the module note.
    Functional { label: "M06-2X", keyword: "M062X", family: Family::Hybrid, composite: None },
    Functional { label: "M06", keyword: "M06", family: Family::Hybrid, composite: None },
    Functional { label: "mPW1PW", keyword: "MPW1PW", family: Family::Hybrid, composite: None },
    Functional { label: "O3LYP", keyword: "O3LYP", family: Family::Hybrid, composite: None },
    Functional { label: "PBE0", keyword: "PBE0", family: Family::Hybrid, composite: None },
    Functional { label: "PW6B95", keyword: "PW6B95", family: Family::Hybrid, composite: None },
    Functional { label: "TPSS0", keyword: "TPSS0", family: Family::Hybrid, composite: None },
    Functional { label: "TPSSh", keyword: "TPSSh", family: Family::Hybrid, composite: None },
    Functional { label: "X3LYP", keyword: "X3LYP", family: Family::Hybrid, composite: None },

    // Range-separated hybrids.
    Functional { label: "CAM-B3LYP", keyword: "CAM-B3LYP", family: Family::RangeSeparated, composite: None },
    Functional { label: "LC-BLYP", keyword: "LC-BLYP", family: Family::RangeSeparated, composite: None },
    Functional { label: "wB97", keyword: "wB97", family: Family::RangeSeparated, composite: None },
    Functional { label: "wB97M-V", keyword: "wB97M-V", family: Family::RangeSeparated, composite: None },
    Functional { label: "wB97X", keyword: "wB97X", family: Family::RangeSeparated, composite: None },
    Functional { label: "wB97X-D3", keyword: "wB97X-D3", family: Family::RangeSeparated, composite: None },
    Functional { label: "wB97X-V", keyword: "wB97X-V", family: Family::RangeSeparated, composite: None },

    // Double hybrids. Each needs an auxiliary basis, added automatically.
    Functional { label: "B2GP-PLYP", keyword: "B2GP-PLYP", family: Family::DoubleHybrid, composite: None },
    Functional { label: "B2PLYP", keyword: "B2PLYP", family: Family::DoubleHybrid, composite: None },
    Functional { label: "DSD-BLYP", keyword: "DSD-BLYP", family: Family::DoubleHybrid, composite: None },
    Functional { label: "DSD-PBEP86", keyword: "DSD-PBEP86", family: Family::DoubleHybrid, composite: None },
    Functional { label: "mPW2PLYP", keyword: "mPW2PLYP", family: Family::DoubleHybrid, composite: None },
    Functional { label: "PWPB95", keyword: "PWPB95", family: Family::DoubleHybrid, composite: None },

    // Composites, each bringing its own basis set and corrections.
    Functional {
        label: "B97-3c",
        keyword: "B97-3c",
        family: Family::Composite,
        composite: Some("refitted B97 GGA + def2-mTZVP + D3(BJ) + SRB"),
    },
    Functional {
        label: "HF-3c",
        keyword: "HF-3c",
        family: Family::Composite,
        composite: Some("plain HF + MINIX + D3(BJ) + gCP/SRB"),
    },
    Functional {
        label: "PBEh-3c",
        keyword: "PBEh-3c",
        family: Family::Composite,
        composite: Some("PBE hybrid (42% HF) + def2-mSVP + D3(BJ) + gCP"),
    },
    Functional {
        label: "r2SCAN-3c",
        keyword: "r2SCAN-3c",
        family: Family::Composite,
        composite: Some("r2SCAN + def2-mTZVPP + D4 + gCP"),
    },
    Functional {
        label: "wB97X-3c",
        keyword: "wB97X-3c",
        family: Family::Composite,
        composite: Some("wB97X + vDZP + built-in ECPs"),
    },
];

/// The default functional.
pub const DEFAULT_FUNCTIONAL: &str = "B3LYP";

/// A dispersion correction, as ORCA's keyword line takes it.
///
/// Its own enum rather than Behemoth's, for the usual reason: Behemoth spells
/// these `d3bj`, `d3zero`, `d4` and offers several damping variants ORCA does
/// not, so one list would have to translate.
///
/// These four come from ORCA's documentation rather than from running them.
/// A keyword this build did not know would fail as an input error within
/// seconds of starting, which is cheap to discover, and the free-text keyword
/// line is there for anything else.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[derive(serde::Serialize, serde::Deserialize)]
pub enum Dispersion {
    #[default]
    None,
    /// Grimme's D3 with zero damping.
    D3Zero,
    /// D3 with Becke-Johnson damping, the usual choice.
    D3Bj,
    /// Grimme's D4.
    D4,
}

impl Dispersion {
    pub const ALL: [Dispersion; 4] =
        [Dispersion::None, Dispersion::D3Zero, Dispersion::D3Bj, Dispersion::D4];

    pub fn label(self) -> &'static str {
        match self {
            Dispersion::None => "None",
            Dispersion::D3Zero => "D3 (zero damping)",
            Dispersion::D3Bj => "D3(BJ)",
            Dispersion::D4 => "D4",
        }
    }

    /// The keyword, or `None` when nothing should be written.
    pub fn keyword(self) -> Option<&'static str> {
        match self {
            Dispersion::None => None,
            Dispersion::D3Zero => Some("D3ZERO"),
            Dispersion::D3Bj => Some("D3BJ"),
            Dispersion::D4 => Some("D4"),
        }
    }
}

/// The dispersion keyword to write for this configuration, if any.
///
/// A composite brings its own correction, so nothing is added alongside one
/// even when the field was set before the composite was chosen -- the same rule
/// the basis follows, and for the same reason.
///
/// A functional that already carries a correction in its own name gets none
/// either: `B97-D3` is D3-corrected by definition, and `wB97X-V` and `B97M-V`
/// carry the non-local VV10 term. Adding D3 on top would double-count.
pub fn dispersion_keyword(functional: &str, dispersion: Dispersion) -> Option<&'static str> {
    if is_composite(functional) || carries_own_dispersion(functional) {
        return None;
    }
    dispersion.keyword()
}

/// Whether a functional's own definition already includes a dispersion or
/// non-local correlation term.
pub fn carries_own_dispersion(functional: &str) -> bool {
    let name = functional.to_ascii_uppercase();
    name.ends_with("-D3") || name.ends_with("-D4") || name.ends_with("-V")
}

/// A group of basis sets, as the picker shows them.
pub struct BasisGroup {
    pub name: &'static str,
    /// Display name and ORCA keyword. Identical here, since ORCA's own
    /// spellings are already readable.
    pub sets: &'static [(&'static str, &'static str)],
}

/// Every basis set offered, grouped by family.
///
/// All 43 were run against ORCA 6.1.1 and terminated normally. The families
/// match Behemoth's catalogue, so the same basis is findable in the same place
/// whichever program is selected.
pub const BASIS_GROUPS: [BasisGroup; 6] = [
    BasisGroup {
        name: "def2",
        sets: &[
            ("def2-SV(P)", "def2-SV(P)"),
            ("def2-SVP", "def2-SVP"),
            ("def2-SVPD", "def2-SVPD"),
            ("def2-TZVP", "def2-TZVP"),
            ("def2-TZVPD", "def2-TZVPD"),
            ("def2-TZVPP", "def2-TZVPP"),
            ("def2-TZVPPD", "def2-TZVPPD"),
            ("def2-QZVP", "def2-QZVP"),
            ("def2-QZVPD", "def2-QZVPD"),
            ("def2-QZVPP", "def2-QZVPP"),
            ("def2-QZVPPD", "def2-QZVPPD"),
        ],
    },
    BasisGroup {
        name: "def2, minimally augmented",
        sets: &[("ma-def2-SVP", "ma-def2-SVP"), ("ma-def2-TZVP", "ma-def2-TZVP")],
    },
    BasisGroup {
        name: "pc",
        sets: &[
            ("pc-0", "pc-0"),
            ("pc-1", "pc-1"),
            ("pc-2", "pc-2"),
            ("pc-3", "pc-3"),
            ("pc-4", "pc-4"),
            ("aug-pc-0", "aug-pc-0"),
            ("aug-pc-1", "aug-pc-1"),
            ("aug-pc-2", "aug-pc-2"),
            ("aug-pc-3", "aug-pc-3"),
            ("aug-pc-4", "aug-pc-4"),
        ],
    },
    BasisGroup {
        name: "cc",
        sets: &[
            ("cc-pVDZ", "cc-pVDZ"),
            ("cc-pVTZ", "cc-pVTZ"),
            ("cc-pVQZ", "cc-pVQZ"),
            ("cc-pV5Z", "cc-pV5Z"),
            ("aug-cc-pVDZ", "aug-cc-pVDZ"),
            ("aug-cc-pVTZ", "aug-cc-pVTZ"),
            ("aug-cc-pVQZ", "aug-cc-pVQZ"),
            ("aug-cc-pV5Z", "aug-cc-pV5Z"),
        ],
    },
    BasisGroup {
        name: "Pople",
        sets: &[
            ("3-21G", "3-21G"),
            ("6-31G", "6-31G"),
            ("6-31G*", "6-31G*"),
            ("6-31G**", "6-31G**"),
            ("6-31+G*", "6-31+G*"),
            ("6-31++G**", "6-31++G**"),
            ("6-311G", "6-311G"),
            ("6-311G*", "6-311G*"),
            ("6-311G**", "6-311G**"),
            ("6-311+G*", "6-311+G*"),
            ("6-311++G**", "6-311++G**"),
        ],
    },
    BasisGroup {
        name: "STO-nG",
        // ORCA 6.1.1 ships STO-3G only; STO-4G/5G/6G are not keywords here.
        sets: &[("STO-3G", "STO-3G")],
    },
];

/// The default basis set.
pub const DEFAULT_BASIS: &str = "def2-SVP";

/// The catalogue entry for a keyword, if it is one we offer.
pub fn functional(keyword: &str) -> Option<&'static Functional> {
    FUNCTIONALS
        .iter()
        .find(|f| f.keyword.eq_ignore_ascii_case(keyword))
}

/// Whether a keyword is a composite, which brings its own basis set.
pub fn is_composite(keyword: &str) -> bool {
    functional(keyword).is_some_and(|f| f.family == Family::Composite)
}

/// What a composite brings, for the note beside it.
pub fn composite_note(keyword: &str) -> Option<&'static str> {
    functional(keyword).and_then(|f| f.composite)
}

/// Whether a keyword names a double hybrid, which carries an MP2-like
/// correlation term.
pub fn is_double_hybrid(keyword: &str) -> bool {
    functional(keyword).is_some_and(|f| f.family == Family::DoubleHybrid)
}

/// The auxiliary basis a double hybrid needs, or `None` when none is wanted.
///
/// A double hybrid computes an MP2-like correlation energy, which ORCA will not
/// do without a correlation-fitting basis. Tested directly: `! B2PLYP def2-SVP`
/// is an input error, `! B2PLYP def2-SVP def2-SVP/C` runs. Supplying it here
/// rather than leaving it to the user is the difference between a double hybrid
/// that works and one that fails with a message about a missing auxiliary set.
///
/// A composite is excluded because it defines its own basis, so there is no
/// basis name here to build a `/C` set from.
pub fn auxiliary_basis(keyword: &str, basis: &str) -> Option<String> {
    if !is_double_hybrid(keyword) || is_composite(keyword) || basis.trim().is_empty() {
        return None;
    }
    Some(format!("{}/C", basis.trim()))
}

/// The display name for a keyword, falling back to the keyword itself so a
/// hand-typed one still shows something.
pub fn functional_label(keyword: &str) -> &str {
    functional(keyword).map_or(keyword, |f| f.label)
}

/// The display name for a basis keyword.
pub fn basis_label(keyword: &str) -> &str {
    for group in &BASIS_GROUPS {
        for (label, name) in group.sets {
            if name.eq_ignore_ascii_case(keyword) {
                return label;
            }
        }
    }
    keyword
}

/// The functionals in one family, for a dropdown that groups them.
pub fn functionals_in(family: Family) -> impl Iterator<Item = &'static Functional> {
    FUNCTIONALS.iter().filter(move |f| f.family == family)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_defaults_are_listed_in_the_catalogue() {
        assert!(functional(DEFAULT_FUNCTIONAL).is_some());
        assert_eq!(basis_label(DEFAULT_BASIS), "def2-SVP");
    }

    /// The catalogue is the size it is because every entry was run. A silent
    /// shrink would mean someone dropped entries without testing.
    #[test]
    fn the_catalogue_is_the_full_tested_list() {
        assert_eq!(FUNCTIONALS.len(), 48);
        let basis_count: usize = BASIS_GROUPS.iter().map(|g| g.sets.len()).sum();
        assert_eq!(basis_count, 43);
    }

    /// The spelling ORCA actually accepts. Tested against ORCA 6.1.1: the
    /// hyphenated form is an input error, so the catalogue must not carry it.
    #[test]
    fn the_minnesota_hybrid_uses_orcas_spelling() {
        let entry = functional("M062X").expect("M062X is ORCA's spelling");
        assert_eq!(entry.keyword, "M062X");
        assert_eq!(entry.label, "M06-2X", "the label stays readable");
        assert!(
            !FUNCTIONALS.iter().any(|f| f.keyword == "M06-2X"),
            "M06-2X is an ORCA input error and must not be sent"
        );
    }

    /// Keywords this build rejects must not be offered. Each was tried, under
    /// the spellings a user would expect, and none ran.
    #[test]
    fn keywords_this_build_rejects_are_absent() {
        for absent in ["SCAN", "PBEsol", "HCTH407", "M05-2X", "M052X", "MN15", "BMK"] {
            assert!(
                functional(absent).is_none(),
                "{absent} does not run in ORCA 6.1.1 and must not be offered"
            );
        }
        // The one that does work and is easily confused with SCAN.
        assert!(functional("r2SCAN").is_some());
    }

    /// A double hybrid needs a correlation-fitting basis or ORCA refuses the
    /// input. Tested directly against the binary.
    #[test]
    fn a_double_hybrid_gets_an_auxiliary_basis() {
        assert!(is_double_hybrid("B2PLYP"));
        assert_eq!(
            auxiliary_basis("B2PLYP", "def2-SVP").as_deref(),
            Some("def2-SVP/C")
        );
        assert_eq!(
            auxiliary_basis("DSD-PBEP86", "def2-TZVP").as_deref(),
            Some("def2-TZVP/C")
        );
    }

    /// Nothing else gets one: an auxiliary basis on an ordinary hybrid is
    /// wasted setup at best.
    #[test]
    fn nothing_else_gets_an_auxiliary_basis() {
        for keyword in ["B3LYP", "PBE0", "HF", "r2SCAN", "wB97X"] {
            assert_eq!(auxiliary_basis(keyword, "def2-SVP"), None, "{keyword}");
        }
        // A composite defines its own basis, so there is nothing to build from.
        assert_eq!(auxiliary_basis("r2SCAN-3c", "def2-SVP"), None);
        // And no basis means no auxiliary basis.
        assert_eq!(auxiliary_basis("B2PLYP", "  "), None);
    }

    /// A composite defines its own basis, so the panel must be able to tell
    /// which entries take the basis field away.
    #[test]
    fn the_composites_are_marked_and_describe_themselves() {
        for keyword in ["r2SCAN-3c", "B97-3c", "PBEh-3c", "HF-3c", "wB97X-3c"] {
            assert!(is_composite(keyword), "{keyword}");
            let note = composite_note(keyword).unwrap_or("");
            assert!(!note.is_empty(), "{keyword} says nothing about itself");
        }
        for keyword in ["B3LYP", "PBE0", "M062X", "HF", "B2PLYP"] {
            assert!(!is_composite(keyword), "{keyword}");
        }
    }

    /// Every entry belongs to a family the dropdown lists, or it would be
    /// invisible.
    #[test]
    fn every_functional_is_reachable_through_its_family() {
        let mut seen = 0usize;
        for family in Family::ALL {
            let count = functionals_in(family).count();
            assert!(!family.heading().is_empty());
            seen += count;
        }
        assert_eq!(seen, FUNCTIONALS.len(), "a functional is in no listed family");
    }

    /// No keyword appears twice, in either catalogue.
    #[test]
    fn no_keyword_is_listed_twice() {
        let mut names: Vec<String> = FUNCTIONALS
            .iter()
            .map(|f| f.keyword.to_ascii_lowercase())
            .collect();
        let count = names.len();
        names.sort();
        names.dedup();
        assert_eq!(names.len(), count, "a functional is listed twice");

        let mut bases: Vec<String> = BASIS_GROUPS
            .iter()
            .flat_map(|g| g.sets.iter().map(|(_, k)| k.to_ascii_lowercase()))
            .collect();
        let count = bases.len();
        bases.sort();
        bases.dedup();
        assert_eq!(bases.len(), count, "a basis set is listed twice");
    }

    /// The families Behemoth groups its basis sets by are all present, so the
    /// same set is found in the same place whichever program is selected.
    #[test]
    fn the_basis_families_match_the_ones_behemoth_offers() {
        let names: Vec<&str> = BASIS_GROUPS.iter().map(|g| g.name).collect();
        for family in ["def2", "pc", "cc", "Pople", "STO-nG"] {
            assert!(
                names.iter().any(|n| n.starts_with(family)),
                "no {family} group: {names:?}"
            );
        }
    }

    /// An ordinary functional takes whichever correction was chosen.
    #[test]
    fn an_ordinary_functional_takes_the_chosen_correction() {
        assert_eq!(dispersion_keyword("B3LYP", Dispersion::D3Bj), Some("D3BJ"));
        assert_eq!(dispersion_keyword("PBE0", Dispersion::D4), Some("D4"));
        assert_eq!(dispersion_keyword("TPSS", Dispersion::D3Zero), Some("D3ZERO"));
        assert_eq!(dispersion_keyword("B3LYP", Dispersion::None), None);
    }

    /// A composite brings its own correction, so none is added alongside it --
    /// even when the field was set before the composite was chosen. The same
    /// rule the basis follows.
    #[test]
    fn a_composite_takes_no_extra_correction() {
        for keyword in ["r2SCAN-3c", "B97-3c", "PBEh-3c", "HF-3c", "wB97X-3c"] {
            assert_eq!(
                dispersion_keyword(keyword, Dispersion::D4),
                None,
                "{keyword} already has one"
            );
        }
    }

    /// A functional whose own name carries a correction gets none on top.
    /// `B97-D3` is D3-corrected by definition, and the `-V` functionals carry
    /// the non-local VV10 term; adding D3 as well would double-count, and the
    /// energy that came back would look entirely ordinary.
    #[test]
    fn a_functional_that_carries_its_own_correction_gets_no_second_one() {
        for keyword in ["B97-D3", "wB97X-D3", "wB97X-V", "wB97M-V", "B97M-V"] {
            assert!(carries_own_dispersion(keyword), "{keyword}");
            assert_eq!(
                dispersion_keyword(keyword, Dispersion::D3Bj),
                None,
                "{keyword} would be double-counted"
            );
        }
        // And the ordinary ones are not caught by that rule.
        for keyword in ["B3LYP", "PBE0", "M062X", "TPSSh", "r2SCAN"] {
            assert!(!carries_own_dispersion(keyword), "{keyword}");
        }
    }

    /// Every correction has a keyword except the one that means "none".
    #[test]
    fn every_correction_but_none_has_a_keyword() {
        for entry in Dispersion::ALL {
            assert!(!entry.label().is_empty());
            assert_eq!(
                entry.keyword().is_none(),
                entry == Dispersion::None,
                "{entry:?}"
            );
        }
        assert_eq!(Dispersion::default(), Dispersion::None);
    }

    /// An unlisted keyword still shows as itself.
    #[test]
    fn an_unlisted_keyword_shows_as_itself() {
        assert_eq!(functional_label("SCAN"), "SCAN");
        assert_eq!(basis_label("ZORA-def2-TZVP"), "ZORA-def2-TZVP");
    }
}
