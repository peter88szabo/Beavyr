//! The level of theory a Psi4 calculation runs at.
//!
//! Its own catalogue, like every other program's. Every name below was read out
//! of the installed Psi4 (1.11) rather than taken from documentation: the
//! functionals from `psi4.driver.procrouting.dft.functionals`, which holds 1339
//! entries, and the basis sets from the `.gbs` files in Psi4's own data
//! directory, of which there are 538. Neither list is offered whole; a dropdown
//! of a thousand functionals would be unusable.
//!
//! Psi4's basis files are named with the punctuation flattened -- `6-31gss.gbs`
//! for `6-31G**` -- but an input file takes the readable spelling, so that is
//! what the catalogue stores.

/// Which family a functional belongs to, so a list of this size can be read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Family {
    HartreeFock,
    LdaGga,
    MetaGga,
    Hybrid,
    RangeSeparated,
    DoubleHybrid,
    /// Brings its own basis set and corrections.
    Composite,
}

impl Family {
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
    /// What the generated input passes to `energy()` and its relatives.
    pub name: &'static str,
    pub family: Family,
    /// What a composite brings with it.
    pub composite: Option<&'static str>,
}

impl Functional {
    /// Whether this is Hartree-Fock rather than a density functional.
    ///
    /// Psi4 takes `scf` as a method name like any other, so unlike PySCF this
    /// needs no different object built. It still matters: Hartree-Fock takes no
    /// dispersion correction and is not a functional.
    pub fn is_hartree_fock(&self) -> bool {
        self.family == Family::HartreeFock
    }
}

/// The functionals offered, grouped by family.
///
/// Every name was confirmed present in this build's functional dictionary.
/// Psi4's own `-d3bj` suffixed variants are how it asks for a dispersion
/// correction, which is why the dispersion choice is a suffix here rather than
/// a separate keyword as it is for ORCA.
pub const FUNCTIONALS: [Functional; 32] = [
    Functional { label: "HF", name: "scf", family: Family::HartreeFock, composite: None },

    Functional { label: "BLYP", name: "blyp", family: Family::LdaGga, composite: None },
    Functional { label: "BP86", name: "bp86", family: Family::LdaGga, composite: None },
    Functional { label: "B97-D3(BJ)", name: "b97-d3bj", family: Family::LdaGga, composite: None },
    Functional { label: "HCTH407", name: "hcth407", family: Family::LdaGga, composite: None },
    Functional { label: "PBE", name: "pbe", family: Family::LdaGga, composite: None },
    Functional { label: "revPBE", name: "revpbe", family: Family::LdaGga, composite: None },
    Functional { label: "SOGGA11", name: "sogga11", family: Family::LdaGga, composite: None },

    Functional { label: "M06-L", name: "m06-l", family: Family::MetaGga, composite: None },
    Functional { label: "r2SCAN", name: "r2scan", family: Family::MetaGga, composite: None },
    Functional { label: "SCAN", name: "scan", family: Family::MetaGga, composite: None },
    Functional { label: "TPSS", name: "tpss", family: Family::MetaGga, composite: None },

    Functional { label: "B3LYP", name: "b3lyp", family: Family::Hybrid, composite: None },
    Functional { label: "B3LYP(VWN5)", name: "b3lyp5", family: Family::Hybrid, composite: None },
    Functional { label: "B97-0", name: "b97-0", family: Family::Hybrid, composite: None },
    Functional { label: "BHandHLYP", name: "bhhlyp", family: Family::Hybrid, composite: None },
    Functional { label: "M06", name: "m06", family: Family::Hybrid, composite: None },
    Functional { label: "M06-2X", name: "m06-2x", family: Family::Hybrid, composite: None },
    Functional { label: "MN15", name: "mn15", family: Family::Hybrid, composite: None },
    Functional { label: "PBE0", name: "pbe0", family: Family::Hybrid, composite: None },
    Functional { label: "PW6B95", name: "pw6b95", family: Family::Hybrid, composite: None },
    Functional { label: "TPSSh", name: "tpssh", family: Family::Hybrid, composite: None },

    Functional {
        label: "CAM-B3LYP",
        name: "cam-b3lyp",
        family: Family::RangeSeparated,
        composite: None,
    },
    Functional { label: "wB97X", name: "wb97x", family: Family::RangeSeparated, composite: None },
    Functional { label: "wB97X-D", name: "wb97x-d", family: Family::RangeSeparated, composite: None },
    Functional {
        label: "wB97X-D3(BJ)",
        name: "wb97x-d3bj",
        family: Family::RangeSeparated,
        composite: None,
    },
    Functional { label: "wB97X-V", name: "wb97x-v", family: Family::RangeSeparated, composite: None },
    Functional { label: "wB97M-V", name: "wb97m-v", family: Family::RangeSeparated, composite: None },

    Functional { label: "B2PLYP", name: "b2plyp", family: Family::DoubleHybrid, composite: None },
    Functional {
        label: "B2PLYP-D3(BJ)",
        name: "b2plyp-d3bj",
        family: Family::DoubleHybrid,
        composite: None,
    },
    Functional {
        label: "DSD-BLYP-D3(BJ)",
        name: "dsd-blyp-d3bj",
        family: Family::DoubleHybrid,
        composite: None,
    },

    Functional {
        label: "r2SCAN-3c",
        name: "r2scan-3c",
        family: Family::Composite,
        composite: Some("r2SCAN + def2-mTZVPP + D4 + gCP"),
    },
];

/// The default functional.
pub const DEFAULT_FUNCTIONAL: &str = "b3lyp";

/// A group of basis sets, as the picker shows them.
pub struct BasisGroup {
    pub name: &'static str,
    /// Display name and the spelling written into the input. Identical here:
    /// Psi4 takes the readable form even though its files are named otherwise.
    pub sets: &'static [(&'static str, &'static str)],
}

/// The basis sets offered, grouped by family. All confirmed present in this
/// build's basis directory.
pub const BASIS_GROUPS: [BasisGroup; 5] = [
    BasisGroup {
        name: "def2",
        sets: &[
            ("def2-SVP", "def2-SVP"),
            ("def2-SVPD", "def2-SVPD"),
            ("def2-TZVP", "def2-TZVP"),
            ("def2-TZVPP", "def2-TZVPP"),
            ("def2-TZVPPD", "def2-TZVPPD"),
            ("def2-QZVP", "def2-QZVP"),
            ("def2-QZVPP", "def2-QZVPP"),
            ("def2-QZVPPD", "def2-QZVPPD"),
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
        name: "pcseg",
        sets: &[
            ("pcseg-0", "pcseg-0"),
            ("pcseg-1", "pcseg-1"),
            ("pcseg-2", "pcseg-2"),
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
        sets: &[("STO-3G", "STO-3G"), ("STO-6G", "STO-6G")],
    },
];

/// The default basis set.
pub const DEFAULT_BASIS: &str = "def2-SVP";

/// The catalogue entry for a name, if it is one we offer.
pub fn functional(name: &str) -> Option<&'static Functional> {
    FUNCTIONALS.iter().find(|f| f.name.eq_ignore_ascii_case(name))
}

/// Whether this name selects Hartree-Fock.
pub fn is_hartree_fock(name: &str) -> bool {
    functional(name).is_some_and(|f| f.is_hartree_fock())
}

/// Whether this name is a composite, which brings its own basis set.
pub fn is_composite(name: &str) -> bool {
    functional(name).is_some_and(|f| f.family == Family::Composite)
}

/// What a composite brings, for the note beside it.
pub fn composite_note(name: &str) -> Option<&'static str> {
    functional(name).and_then(|f| f.composite)
}

/// The display name for a functional, falling back to the name itself.
pub fn functional_label(name: &str) -> &str {
    functional(name).map_or(name, |f| f.label)
}

/// The display name for a basis set.
pub fn basis_label(name: &str) -> &str {
    for group in &BASIS_GROUPS {
        for (label, spelling) in group.sets {
            if spelling.eq_ignore_ascii_case(name) {
                return label;
            }
        }
    }
    name
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

    /// Psi4 takes Hartree-Fock as the method name `scf`, not `hf`. Sending the
    /// obvious spelling would fail, and the label stays readable either way.
    #[test]
    fn hartree_fock_travels_as_psi4_calls_it() {
        let entry = functional("scf").expect("Psi4 spells it scf");
        assert_eq!(entry.label, "HF");
        assert!(is_hartree_fock("scf"));
        assert!(functional("hf").is_none(), "hf is not a Psi4 method name");
    }

    /// A dispersion correction is part of the functional's name in Psi4, not a
    /// separate keyword as it is for ORCA. That is why there is no dispersion
    /// dropdown for this program.
    #[test]
    fn dispersion_is_carried_in_the_functional_name() {
        for name in ["b97-d3bj", "wb97x-d3bj", "b2plyp-d3bj", "dsd-blyp-d3bj"] {
            assert!(functional(name).is_some(), "{name}");
        }
    }

    /// Every entry belongs to a family the dropdown lists, or it is invisible.
    #[test]
    fn every_functional_is_reachable_through_its_family() {
        let seen: usize = Family::ALL.iter().map(|&f| functionals_in(f).count()).sum();
        assert_eq!(seen, FUNCTIONALS.len());
        for family in Family::ALL {
            assert!(!family.heading().is_empty());
        }
    }

    /// A composite brings its own basis, so the panel can take that field away.
    #[test]
    fn the_composite_is_marked_and_describes_itself() {
        assert!(is_composite("r2scan-3c"));
        assert!(composite_note("r2scan-3c").is_some_and(|n| !n.is_empty()));
        assert!(!is_composite("b3lyp"));
    }

    /// No name appears twice, in either catalogue.
    #[test]
    fn no_name_is_listed_twice() {
        let mut names: Vec<String> =
            FUNCTIONALS.iter().map(|f| f.name.to_ascii_lowercase()).collect();
        let count = names.len();
        names.sort();
        names.dedup();
        assert_eq!(names.len(), count, "a functional is listed twice");

        let mut bases: Vec<String> = BASIS_GROUPS
            .iter()
            .flat_map(|g| g.sets.iter().map(|(_, s)| s.to_ascii_lowercase()))
            .collect();
        let count = bases.len();
        bases.sort();
        bases.dedup();
        assert_eq!(bases.len(), count, "a basis set is listed twice");
    }

    /// An unlisted name still shows as itself.
    #[test]
    fn an_unlisted_name_shows_as_itself() {
        assert_eq!(functional_label("wb97x-d3"), "wb97x-d3");
        assert_eq!(basis_label("ano-rcc"), "ano-rcc");
    }
}
