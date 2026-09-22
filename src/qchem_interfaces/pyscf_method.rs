//! The level of theory a PySCF calculation runs at.
//!
//! A catalogue of its own rather than a translation of anyone else's:
//! switching program resets the method to this program's default, and nothing
//! tries to map a name from one code onto another.
//!
//! Every name below was checked against the installed PySCF (2.12.1) rather
//! than taken from documentation: each basis through `gto.basis.load` on
//! carbon, each functional through `dft.libxc.parse_xc`. A name that did not
//! load is not here. That test dropped `STO-4G` and `STO-5G`, which this build
//! does not carry, and it is also why `def2-SV(P)` appears under PySCF's own
//! escaped spelling.
//!
//! The double hybrids are absent on purpose. `parse_xc` rejects `b2plyp` and
//! its relatives because PySCF does not treat them as a single functional: the
//! MP2-like term has to be added by hand afterwards. Offering them here would
//! mean generating a script that silently dropped the correlation term, giving
//! a number that is not the functional it claims.
//!
//! PySCF takes no memory flag of the kind Behemoth has -- it allocates as it
//! goes -- so only the thread count is offered.

/// Which family a functional belongs to, so a list of this size can be read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Family {
    /// Hartree-Fock, which is not a functional at all.
    HartreeFock,
    LdaGga,
    MetaGga,
    Hybrid,
    RangeSeparated,
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
        }
    }

    /// Every family, in the order the dropdown lists them.
    pub const ALL: [Family; 5] = [
        Family::HartreeFock,
        Family::LdaGga,
        Family::MetaGga,
        Family::Hybrid,
        Family::RangeSeparated,
    ];
}

/// A functional as the panel offers it.
pub struct Functional {
    /// What the dropdown shows.
    pub label: &'static str,
    /// What the generated script assigns to `mf.xc`. PySCF is case-insensitive
    /// here; the catalogue stores the lowercase form its own documentation uses.
    pub cli: &'static str,
    pub family: Family,
}

impl Functional {
    /// Whether this is Hartree-Fock rather than a density functional.
    ///
    /// It changes which PySCF class the script builds, `scf.RHF`/`scf.UHF`
    /// instead of `dft.RKS`/`dft.UKS`, so it cannot simply be passed as an
    /// `xc` string.
    pub fn is_hartree_fock(&self) -> bool {
        self.family == Family::HartreeFock
    }
}

/// Every functional offered, grouped by family.
///
/// All 50 were accepted by `dft.libxc.parse_xc` in PySCF 2.12.1.
pub const FUNCTIONALS: [Functional; 50] = [
    Functional { label: "HF", cli: "hf", family: Family::HartreeFock },

    // LDA and GGA.
    Functional { label: "SVWN (LDA)", cli: "svwn", family: Family::LdaGga },
    Functional { label: "B97-D", cli: "b97-d", family: Family::LdaGga },
    Functional { label: "B97-D3", cli: "b97-d3", family: Family::LdaGga },
    Functional { label: "BLYP", cli: "blyp", family: Family::LdaGga },
    Functional { label: "BP86", cli: "bp86", family: Family::LdaGga },
    Functional { label: "GAM", cli: "gam", family: Family::LdaGga },
    Functional { label: "HCTH407", cli: "hcth407", family: Family::LdaGga },
    Functional { label: "N12", cli: "n12", family: Family::LdaGga },
    Functional { label: "OLYP", cli: "olyp", family: Family::LdaGga },
    Functional { label: "PBE", cli: "pbe", family: Family::LdaGga },
    Functional { label: "PBEsol", cli: "pbesol", family: Family::LdaGga },
    Functional { label: "PW91", cli: "pw91", family: Family::LdaGga },
    Functional { label: "revPBE", cli: "revpbe", family: Family::LdaGga },
    Functional { label: "RPBE", cli: "rpbe", family: Family::LdaGga },
    Functional { label: "SOGGA11", cli: "sogga11", family: Family::LdaGga },

    // meta-GGA.
    Functional { label: "B97M-V", cli: "b97m-v", family: Family::MetaGga },
    Functional { label: "M06-L", cli: "m06-l", family: Family::MetaGga },
    Functional { label: "M11-L", cli: "m11-l", family: Family::MetaGga },
    Functional { label: "MN12-L", cli: "mn12-l", family: Family::MetaGga },
    Functional { label: "MN15-L", cli: "mn15-l", family: Family::MetaGga },
    Functional { label: "r2SCAN", cli: "r2scan", family: Family::MetaGga },
    Functional { label: "SCAN", cli: "scan", family: Family::MetaGga },
    Functional { label: "TPSS", cli: "tpss", family: Family::MetaGga },
    Functional { label: "revTPSS", cli: "revtpss", family: Family::MetaGga },

    // Global hybrids.
    Functional { label: "B1LYP", cli: "b1lyp", family: Family::Hybrid },
    Functional { label: "B3LYP", cli: "b3lyp", family: Family::Hybrid },
    // PySCF's VWN5 variant, which is what most other codes call plain B3LYP.
    Functional { label: "B3LYP(VWN5)", cli: "b3lyp5", family: Family::Hybrid },
    Functional { label: "B3P86", cli: "b3p86", family: Family::Hybrid },
    Functional { label: "B3PW91", cli: "b3pw91", family: Family::Hybrid },
    Functional { label: "B97", cli: "b97", family: Family::Hybrid },
    Functional { label: "BHandHLYP", cli: "bhandhlyp", family: Family::Hybrid },
    Functional { label: "M05", cli: "m05", family: Family::Hybrid },
    Functional { label: "M05-2X", cli: "m05-2x", family: Family::Hybrid },
    Functional { label: "M06", cli: "m06", family: Family::Hybrid },
    Functional { label: "M06-2X", cli: "m06-2x", family: Family::Hybrid },
    Functional { label: "MN15", cli: "mn15", family: Family::Hybrid },
    Functional { label: "O3LYP", cli: "o3lyp", family: Family::Hybrid },
    Functional { label: "PBE0", cli: "pbe0", family: Family::Hybrid },
    Functional { label: "PW6B95", cli: "pw6b95", family: Family::Hybrid },
    Functional { label: "TPSS0", cli: "tpss0", family: Family::Hybrid },
    Functional { label: "TPSSh", cli: "tpssh", family: Family::Hybrid },
    Functional { label: "X3LYP", cli: "x3lyp", family: Family::Hybrid },

    // Range-separated hybrids.
    Functional { label: "CAM-B3LYP", cli: "camb3lyp", family: Family::RangeSeparated },
    Functional { label: "HSE03", cli: "hse03", family: Family::RangeSeparated },
    Functional { label: "HSE06", cli: "hse06", family: Family::RangeSeparated },
    Functional { label: "wB97", cli: "wb97", family: Family::RangeSeparated },
    Functional { label: "wB97M-V", cli: "wb97m-v", family: Family::RangeSeparated },
    Functional { label: "wB97X", cli: "wb97x", family: Family::RangeSeparated },
    Functional { label: "wB97X-D", cli: "wb97x-d", family: Family::RangeSeparated },
];

/// The default functional.
pub const DEFAULT_FUNCTIONAL: &str = "b3lyp";

/// A group of basis sets, as the picker shows them.
pub struct BasisGroup {
    pub name: &'static str,
    /// Display name and the name PySCF is given. They differ only where
    /// PySCF's own spelling is unreadable.
    pub sets: &'static [(&'static str, &'static str)],
}

/// Every basis set offered, grouped by family.
///
/// All 42 were loaded from PySCF 2.12.1's own library. The families match
/// Behemoth's catalogue, so the same basis is findable in the same place
/// whichever program is selected.
pub const BASIS_GROUPS: [BasisGroup; 5] = [
    BasisGroup {
        name: "def2",
        sets: &[
            // PySCF escapes the parentheses in its file names, and
            // `def2-sv(p)` does not load. This is the spelling that does.
            ("def2-SV(P)", "def2-sv_p_"),
            ("def2-SVP", "def2-svp"),
            ("def2-SVPD", "def2-svpd"),
            ("def2-TZVP", "def2-tzvp"),
            ("def2-TZVPD", "def2-tzvpd"),
            ("def2-TZVPP", "def2-tzvpp"),
            ("def2-TZVPPD", "def2-tzvppd"),
            ("def2-QZVP", "def2-qzvp"),
            ("def2-QZVPD", "def2-qzvpd"),
            ("def2-QZVPP", "def2-qzvpp"),
            ("def2-QZVPPD", "def2-qzvppd"),
        ],
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
            ("cc-pVDZ", "cc-pvdz"),
            ("cc-pVTZ", "cc-pvtz"),
            ("cc-pVQZ", "cc-pvqz"),
            ("cc-pV5Z", "cc-pv5z"),
            ("aug-cc-pVDZ", "aug-cc-pvdz"),
            ("aug-cc-pVTZ", "aug-cc-pvtz"),
            ("aug-cc-pVQZ", "aug-cc-pvqz"),
            ("aug-cc-pV5Z", "aug-cc-pv5z"),
        ],
    },
    BasisGroup {
        name: "Pople",
        sets: &[
            ("3-21G", "3-21g"),
            ("6-31G", "6-31g"),
            ("6-31G*", "6-31g*"),
            ("6-31G**", "6-31g**"),
            ("6-31+G*", "6-31+g*"),
            ("6-31++G**", "6-31++g**"),
            ("6-311G", "6-311g"),
            ("6-311G*", "6-311g*"),
            ("6-311G**", "6-311g**"),
            ("6-311+G*", "6-311+g*"),
            ("6-311++G**", "6-311++g**"),
        ],
    },
    BasisGroup {
        name: "STO-nG",
        // STO-4G and STO-5G are not in this build's library; they were tried.
        sets: &[("STO-3G", "sto-3g"), ("STO-6G", "sto-6g")],
    },
];

/// The default basis set. Small, so a first run finishes quickly.
pub const DEFAULT_BASIS: &str = "def2-svp";

/// The catalogue entry for a functional name, if it is one we offer.
pub fn functional(cli: &str) -> Option<&'static Functional> {
    FUNCTIONALS.iter().find(|f| f.cli.eq_ignore_ascii_case(cli))
}

/// Whether this name selects Hartree-Fock rather than a density functional.
///
/// An unlisted name is taken to be a functional, because a LibXC specification
/// typed into the free-text field is one, and Hartree-Fock has exactly one
/// spelling.
pub fn is_hartree_fock(cli: &str) -> bool {
    functional(cli).is_some_and(|f| f.is_hartree_fock())
}

/// The display name for a basis name, falling back to the name itself so a
/// hand-typed basis still shows something.
pub fn basis_label(cli: &str) -> &str {
    for group in &BASIS_GROUPS {
        for (label, name) in group.sets {
            if name.eq_ignore_ascii_case(cli) {
                return label;
            }
        }
    }
    cli
}

/// The display name for a functional, falling back to the name itself.
pub fn functional_label(cli: &str) -> &str {
    functional(cli).map_or(cli, |f| f.label)
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

    /// The catalogue is the size it is because every entry was loaded from the
    /// real PySCF. A silent shrink would mean someone dropped entries without
    /// testing.
    #[test]
    fn the_catalogue_is_the_full_tested_list() {
        assert_eq!(FUNCTIONALS.len(), 50);
        let basis_count: usize = BASIS_GROUPS.iter().map(|g| g.sets.len()).sum();
        assert_eq!(basis_count, 42);
    }

    /// Hartree-Fock is not a functional and has to build a different PySCF
    /// object, so it must be distinguishable from every density functional.
    #[test]
    fn hartree_fock_is_marked_and_nothing_else_is() {
        assert!(is_hartree_fock("hf"));
        assert!(is_hartree_fock("HF"), "matching is case-insensitive");
        for entry in &FUNCTIONALS {
            if entry.cli != "hf" {
                assert!(!is_hartree_fock(entry.cli), "{}", entry.cli);
            }
        }
    }

    /// A hand-typed LibXC specification is a functional, not Hartree-Fock.
    /// Getting this wrong would build an `RHF` object and silently ignore the
    /// functional the user asked for.
    #[test]
    fn an_unlisted_name_is_treated_as_a_functional() {
        assert!(!is_hartree_fock("0.2*HF + 0.8*B88, LYP"));
        assert!(!is_hartree_fock("wb97x-d3"));
    }

    /// Names this build does not carry must not be offered. The double hybrids
    /// in particular: PySCF needs the MP2 term added separately, so passing one
    /// as `mf.xc` would give a number that is not the functional it claims.
    #[test]
    fn names_this_build_rejects_are_absent() {
        for absent in ["b2plyp", "mpw2plyp", "b2gpplyp", "dsdblyp", "bpw91", "m11", "bmk"] {
            assert!(
                functional(absent).is_none(),
                "{absent} is not usable through mf.xc and must not be offered"
            );
        }
        for absent in ["sto-4g", "sto-5g", "def2-sv(p)"] {
            let listed = BASIS_GROUPS
                .iter()
                .any(|g| g.sets.iter().any(|(_, cli)| cli.eq_ignore_ascii_case(absent)));
            assert!(!listed, "{absent} does not load in PySCF 2.12.1");
        }
    }

    /// PySCF escapes the parentheses, and the obvious spelling does not load.
    /// The label stays readable while the name sent is the one that works.
    #[test]
    fn def2_sv_p_uses_pyscfs_escaped_spelling() {
        assert_eq!(basis_label("def2-sv_p_"), "def2-SV(P)");
    }

    /// Every entry belongs to a family the dropdown lists, or it would be
    /// invisible.
    #[test]
    fn every_functional_is_reachable_through_its_family() {
        let seen: usize = Family::ALL.iter().map(|&f| functionals_in(f).count()).sum();
        assert_eq!(seen, FUNCTIONALS.len(), "a functional is in no listed family");
        for family in Family::ALL {
            assert!(!family.heading().is_empty());
        }
    }

    /// No name appears twice, in either catalogue: a duplicate would give the
    /// dropdown two identical entries that behave differently.
    #[test]
    fn no_name_is_listed_twice() {
        let mut names: Vec<&str> = FUNCTIONALS.iter().map(|f| f.cli).collect();
        let count = names.len();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), count, "a functional is listed twice");

        let mut bases: Vec<&str> = BASIS_GROUPS
            .iter()
            .flat_map(|g| g.sets.iter().map(|(_, cli)| *cli))
            .collect();
        let count = bases.len();
        bases.sort_unstable();
        bases.dedup();
        assert_eq!(bases.len(), count, "a basis set is listed twice");
    }

    /// The families Behemoth groups its basis sets by are all present, so the
    /// same set is found in the same place whichever program is selected.
    #[test]
    fn the_basis_families_match_the_ones_behemoth_offers() {
        let names: Vec<&str> = BASIS_GROUPS.iter().map(|g| g.name).collect();
        for family in ["def2", "pc", "cc", "Pople", "STO-nG"] {
            assert!(names.contains(&family), "no {family} group: {names:?}");
        }
    }

    /// An unlisted name still shows as itself rather than falling back to
    /// something misleading.
    #[test]
    fn an_unlisted_name_shows_as_itself() {
        assert_eq!(functional_label("wb97x-d3"), "wb97x-d3");
        assert_eq!(basis_label("ano-rcc"), "ano-rcc");
    }
}
