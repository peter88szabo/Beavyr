//! The semi-empirical methods Sparrow offers.
//!
//! Sparrow is not a quantum-chemistry program in the sense the others are: it
//! is a collection of semi-empirical parametrisations. There is no functional
//! and no basis set to choose, only which parametrisation to use, so its method
//! row is one dropdown and the rest of the level-of-theory block disappears.
//!
//! Every name below was resolved through SCINE's own module manager against the
//! installed build, by asking it for a calculator and keeping the names it
//! returned one for.
//!
//! What Sparrow does **not** have is an optimiser or a frequency analysis. It
//! computes an energy, an analytic gradient and an analytic Hessian, and Beavyr
//! does the rest with its own machinery: its optimiser drives the gradients,
//! and its projection and eigensolver turn the Hessian into modes and
//! thermochemistry. That is not a workaround; it is the only route Sparrow has,
//! which is why `own_optimiser` is false for it.

/// Which family a parametrisation belongs to, so the dropdown groups them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Family {
    /// The neglect-of-differential-overlap lineage.
    Nddo,
    /// Density-functional tight binding.
    Dftb,
}

impl Family {
    pub fn heading(self) -> &'static str {
        match self {
            Family::Nddo => "NDDO",
            Family::Dftb => "DFTB",
        }
    }

    pub const ALL: [Family; 2] = [Family::Nddo, Family::Dftb];
}

/// One parametrisation.
pub struct Method {
    /// What the dropdown shows.
    pub label: &'static str,
    /// The name SCINE's module manager is given.
    pub name: &'static str,
    pub family: Family,
}

/// Every parametrisation offered.
///
/// All eight were confirmed by asking SCINE's module manager for a calculator
/// by each name against the installed build.
pub const METHODS: [Method; 8] = [
    Method { label: "MNDO", name: "MNDO", family: Family::Nddo },
    Method { label: "AM1", name: "AM1", family: Family::Nddo },
    Method { label: "RM1", name: "RM1", family: Family::Nddo },
    Method { label: "PM3", name: "PM3", family: Family::Nddo },
    Method { label: "PM6", name: "PM6", family: Family::Nddo },
    // Sparrow's own name for the non-self-consistent-charge variant.
    Method { label: "DFTB0 (non-SCC)", name: "DFTB0", family: Family::Dftb },
    Method { label: "DFTB2", name: "DFTB2", family: Family::Dftb },
    Method { label: "DFTB3", name: "DFTB3", family: Family::Dftb },
];

/// The default parametrisation. PM6 is the most broadly parameterised of the
/// NDDO set and the usual first choice.
pub const DEFAULT_METHOD: &str = "PM6";

/// The catalogue entry for a name, if it is one we offer.
pub fn method(name: &str) -> Option<&'static Method> {
    METHODS.iter().find(|m| m.name.eq_ignore_ascii_case(name))
}

/// The display name, falling back to the name itself.
pub fn method_label(name: &str) -> &str {
    method(name).map_or(name, |m| m.label)
}

/// The parametrisations in one family, for a dropdown that groups them.
pub fn methods_in(family: Family) -> impl Iterator<Item = &'static Method> {
    METHODS.iter().filter(move |m| m.family == family)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_is_listed() {
        assert!(method(DEFAULT_METHOD).is_some());
        assert_eq!(method_label("PM6"), "PM6");
    }

    /// The eight that SCINE actually resolves. A silent shrink would mean
    /// someone dropped one without checking.
    #[test]
    fn the_catalogue_is_the_full_tested_list() {
        assert_eq!(METHODS.len(), 8);
        for name in ["MNDO", "AM1", "RM1", "PM3", "PM6", "DFTB0", "DFTB2", "DFTB3"] {
            assert!(method(name).is_some(), "{name} was confirmed to resolve");
        }
    }

    /// The label says "non-SCC" where the name does not, because that is the
    /// distinction a chemist is choosing on and `DFTB0` alone does not carry it.
    #[test]
    fn the_non_self_consistent_variant_says_so() {
        let entry = method("DFTB0").expect("DFTB0 resolves");
        assert_eq!(entry.name, "DFTB0", "SCINE is given the bare name");
        assert!(entry.label.contains("non-SCC"), "{}", entry.label);
    }

    /// Every entry belongs to a family the dropdown lists, or it is invisible.
    #[test]
    fn every_method_is_reachable_through_its_family() {
        let seen: usize = Family::ALL.iter().map(|&f| methods_in(f).count()).sum();
        assert_eq!(seen, METHODS.len());
        assert_eq!(methods_in(Family::Dftb).count(), 3);
        assert_eq!(methods_in(Family::Nddo).count(), 5);
    }

    #[test]
    fn no_name_is_listed_twice() {
        let mut names: Vec<String> = METHODS.iter().map(|m| m.name.to_ascii_lowercase()).collect();
        let count = names.len();
        names.sort();
        names.dedup();
        assert_eq!(names.len(), count);
    }

    #[test]
    fn an_unlisted_name_shows_as_itself() {
        assert_eq!(method_label("PM7"), "PM7");
    }
}
