//! DREIDING parameters, transcribed from the paper's tables.
//!
//! S. L. Mayo, B. D. Olafson and W. A. Goddard III, *J. Phys. Chem.* **1990**, *94*, 8897-8909.
//!
//! | table | contents | here |
//! |---|---|---|
//! | I | bond radius `R⁰` and bond angle `θ⁰` per atom type | [`GEOMETRY`] |
//! | II | van der Waals `R⁰`, `D⁰` per *element* | [`VDW`] |
//! | III | bond, angle and inversion force constants | [`BOND_K_PER_ORDER`], [`ANGLE_K`], [`INVERSION_K`] |
//! | IV | torsion barriers | [`torsion_for_pair`] |
//! | V | hydrogen bond, no-charges convention | [`D_HB`], [`R_HB`] |
//!
//! Everything is standard DREIDING. The paper also defines **DREIDING/A** (Tables VI and VII), a
//! *different* parameter set -- its O_3 bond angle is 109.471° where Table I gives 104.51° -- and
//! the two are not mixed here. The one exception is the paper's own: Table II takes the van der
//! Waals parameters for Na⁺, Ca²⁺, Fe²⁺ and Zn²⁺ from DREIDING/A, and is followed in that.

/// Bond radius `R⁰` (Å) and equilibrium bond angle `θ⁰` (degrees) per atom type. Table I.
///
/// Table I also lists `Ga3`, `In3` and `O_1`, which the typer's ruleset never emits; they are
/// included so the table matches the paper rather than the typer.
#[rustfmt::skip]
pub const GEOMETRY: &[(&str, f64, f64)] = &[
    // type    R⁰, Å    θ⁰, deg
    ("H_",     0.330,   180.0),
    ("H_HB",   0.330,   180.0),
    ("H_b",    0.510,    90.0),
    ("B_3",    0.880,   109.471),
    ("B_2",    0.790,   120.0),
    ("C_3",    0.770,   109.471),
    ("C_R",    0.700,   120.0),
    ("C_2",    0.670,   120.0),
    ("C_1",    0.602,   180.0),
    ("N_3",    0.702,   106.7),
    ("N_R",    0.650,   120.0),
    ("N_2",    0.615,   120.0),
    ("N_1",    0.556,   180.0),
    ("O_3",    0.660,   104.51),
    ("O_R",    0.660,   120.0),
    ("O_2",    0.560,   120.0),
    ("O_1",    0.528,   180.0),
    ("F_",     0.611,   180.0),
    ("Al3",    1.047,   109.471),
    ("Si3",    0.937,   109.471),
    ("P_3",    0.890,    93.3),
    ("S_3",    1.040,    92.1),
    ("Cl",     0.997,   180.0),
    ("Ga3",    1.210,   109.471),
    ("Ge3",    1.210,   109.471),
    ("As3",    1.210,    92.1),
    ("Se3",    1.210,    90.6),
    ("Br",     1.167,   180.0),
    ("In3",    1.390,   109.471),
    ("Sn3",    1.373,   109.471),
    ("Sb3",    1.432,    91.6),
    ("Te3",    1.280,    90.3),
    ("I_",     1.360,   180.0),
    ("Na",     1.860,    90.0),
    ("Ca",     1.940,    90.0),
    ("Fe",     1.285,    90.0),
    ("Zn",     1.330,   109.471),
];

/// van der Waals `R⁰` (Å) and `D⁰` (kcal/mol), Table II.
///
/// Indexed by **element symbol**, not atom type: DREIDING's nonbonded parameters do not depend on
/// hybridisation. The two exceptions are hydrogens, which are looked up by type -- `H_HB` is given
/// an almost-zero `D⁰` because a hydrogen-bonding hydrogen's attraction is carried by the explicit
/// hydrogen-bond term instead, and double-counting it would over-bind.
///
/// The paper's `ζ` column (for the exponential-6 form) is omitted: the Lennard-Jones 12-6 form is
/// used here, and `ζ` plays no part in it.
#[rustfmt::skip]
pub const VDW: &[(&str, f64, f64)] = &[
    // element  R⁰, Å    D⁰, kcal/mol
    ("H",       3.195,   0.0152),
    ("H_b",     3.195,   0.0152),   // by type: diborane bridging hydrogen
    ("H_HB",    3.195,   0.0001),   // by type: attraction lives in the H-bond term
    ("B",       4.02,    0.095),
    ("C",       3.8983,  0.0951),
    ("N",       3.6621,  0.0774),
    ("O",       3.4046,  0.0957),
    ("F",       3.4720,  0.0725),
    ("Al",      4.39,    0.31),
    ("Si",      4.27,    0.31),
    ("P",       4.1500,  0.3200),
    ("S",       4.0300,  0.3440),
    ("Cl",      3.9503,  0.2833),
    ("Ga",      4.39,    0.40),
    ("Ge",      4.27,    0.40),
    ("As",      4.15,    0.41),
    ("Se",      4.03,    0.43),
    ("Br",      3.95,    0.37),
    ("In",      4.59,    0.55),
    ("Sn",      4.47,    0.55),
    ("Sb",      4.35,    0.55),
    ("Te",      4.23,    0.57),
    ("I",       4.15,    0.51),
    ("Na",      3.144,   0.5),      // Na⁺, from DREIDING/A per Table II
    ("Ca",      3.472,   0.05),     // Ca²⁺, ditto
    ("Fe",      4.54,    0.055),    // Fe²⁺, ditto
    ("Zn",      4.54,    0.055),    // Zn²⁺, ditto
];

/// Bond force constant `K` (kcal/mol/Å²) for bond order `n`, Table III: `K(n) = n × 700`.
pub const BOND_K_PER_ORDER: f64 = 700.0;

/// Bond dissociation energy `D` (kcal/mol) for bond order `n`, Table III: `D(n) = n × 70`.
///
/// Used only by the Morse form (`DREIDING/M`, Eq 5a), which is not implemented -- the paper itself
/// recommends the harmonic form as the default because a sketched starting geometry may be far
/// from equilibrium, which is exactly the case here. Kept so the table is complete.
pub const BOND_D_PER_ORDER: f64 = 70.0;

/// Angle bend force constant (kcal/mol/rad²), Table III and Eq 12. Independent of I, J, K.
pub const ANGLE_K: f64 = 100.0;

/// Inversion force constant (kcal/mol/rad²), Table III and Eq 29.
pub const INVERSION_K: f64 = 40.0;

/// Correction subtracted from the sum of bond radii, Eq 6: `R⁰_IJ = R⁰_I + R⁰_J − 0.01 Å`.
pub const BOND_LENGTH_CORRECTION: f64 = 0.01;

/// Hydrogen bond well depth (kcal/mol), Table V.
///
/// Table V makes this **conditional on the charge convention**: 4.0 with experimental charges,
/// 7.0 with Gasteiger charges, 9.5 with none. Beavyr has no charge model, so the no-charges
/// convention is used throughout -- which the paper supports directly, reproducing the water dimer
/// at −6.03 kcal/mol and 2.92 Å against 6.13 and 2.94 experimentally.
///
/// **If partial charges are ever added, this must change to 7.0** (Gasteiger) to stay consistent.
pub const D_HB: f64 = 9.5;

/// Hydrogen bond equilibrium donor-acceptor distance (Å), Table V.
pub const R_HB: f64 = 2.75;

/// Converts force in kcal/mol/Å divided by mass in amu to acceleration in Å/fs². For dynamics.
pub const FORCE_TO_ACCEL: f64 = 4.184e-4;

/// Boltzmann constant in kcal/mol/K. For dynamics.
pub const BOLTZMANN: f64 = 0.0019_872_041;

/// Looks up Table I geometry for an atom type.
pub fn geometry(atom_type: &str) -> Option<(f64, f64)> {
    GEOMETRY
        .iter()
        .find(|(t, _, _)| *t == atom_type)
        .map(|(_, r, a)| (*r, *a))
}

/// Looks up Table II van der Waals parameters for an atom type.
///
/// Hydrogens are matched by type first (`H_b`, `H_HB` differ), everything else by element, which
/// is the leading alphabetic run of the type string: `C_3` and `C_R` both give `C`, `Si3` gives
/// `Si`.
pub fn vdw(atom_type: &str) -> Option<(f64, f64)> {
    if let Some(hit) = VDW.iter().find(|(t, _, _)| *t == atom_type) {
        return Some((hit.1, hit.2));
    }
    let element: String = atom_type
        .chars()
        .take_while(|c| c.is_ascii_alphabetic())
        .collect();
    VDW.iter()
        .find(|(t, _, _)| *t == element)
        .map(|(_, r, d)| (*r, *d))
}

/// Whether an atom type has a complete standard-DREIDING parameter set.
///
/// The typer's ruleset can emit 56 types; the paper parameterises 34 of them. A type without
/// parameters must make the caller **refuse**, naming the atom, rather than fall back to a guess:
/// an invented radius produces a confident, wrong geometry, which is worse than no geometry.
///
/// The uncovered types are the transition and alkali/alkaline-earth metals outside
/// {Na, Ca, Fe, Zn}, plus `S_2` and `S_R`. The sulfur pair is the one that bites in practice --
/// **thiophene types as `S_R`** -- because Table I gives a bond radius and angle for sp³ `S_3`
/// only, even though sulfur's van der Waals parameters are per element and so are available.
pub fn is_parameterized(atom_type: &str) -> bool {
    geometry(atom_type).is_some() && vdw(atom_type).is_some()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every type in Table I must also resolve van der Waals parameters, or a molecule could type
    /// successfully and then fail to get an energy.
    #[test]
    fn every_geometry_type_has_vdw_parameters() {
        for (atom_type, _, _) in GEOMETRY {
            assert!(
                vdw(atom_type).is_some(),
                "{atom_type} has Table I geometry but no Table II van der Waals parameters"
            );
        }
    }

    #[test]
    fn vdw_falls_back_from_type_to_element() {
        // C_3 and C_R share carbon's parameters; DREIDING's nonbonded terms ignore hybridisation.
        assert_eq!(vdw("C_3"), vdw("C_R"));
        assert_eq!(vdw("C_3"), Some((3.8983, 0.0951)));
        assert_eq!(vdw("Si3"), Some((4.27, 0.31)));
    }

    /// H_HB must *not* fall back to plain hydrogen: its attraction is carried by the explicit
    /// hydrogen-bond term, and taking H's D⁰ as well would double-count it.
    #[test]
    fn hydrogen_variants_are_distinguished_by_type() {
        assert_eq!(vdw("H_"), Some((3.195, 0.0152)));
        assert_eq!(vdw("H_b"), Some((3.195, 0.0152)));
        assert_eq!(vdw("H_HB"), Some((3.195, 0.0001)));
    }

    #[test]
    fn tables_have_no_duplicate_entries() {
        for table_name in ["GEOMETRY", "VDW"] {
            let types: Vec<&str> = if table_name == "GEOMETRY" {
                GEOMETRY.iter().map(|(t, _, _)| *t).collect()
            } else {
                VDW.iter().map(|(t, _, _)| *t).collect()
            };
            let mut sorted = types.clone();
            sorted.sort_unstable();
            let before = sorted.len();
            sorted.dedup();
            assert_eq!(before, sorted.len(), "{table_name} has a duplicate entry");
        }
    }

    /// Guards the coverage gap so it cannot be closed by accident. Widening this needs a
    /// deliberate decision about where the new parameters came from.
    #[test]
    fn coverage_gap_is_what_the_paper_leaves_uncovered() {
        // Every type the vendored typer's ruleset can emit.
        #[rustfmt::skip]
        const EMITTABLE: &[&str] = &[
            "Ag", "Al3", "As3", "Au", "B_2", "B_3", "Ba", "Br", "C_1", "C_2", "C_3", "Ca", "Cd",
            "Cl", "Co", "C_R", "Cs", "Cu", "F_", "Fe", "Ge3", "H_", "H_b", "Hg", "H_HB", "I_",
            "K", "Li", "Mg", "Mn", "N_1", "N_2", "N_3", "Na", "Ni", "N_R", "O_2", "O_3", "O_R",
            "P_3", "Pd", "Pt", "Rb", "Ru", "S_2", "S_3", "Sb3", "Se3", "Si3", "Sn3", "Sr", "S_R",
            "Tc", "Te3", "Ti", "Zn",
        ];
        #[rustfmt::skip]
        const UNPARAMETERIZED: &[&str] = &[
            "Ag", "Au", "Ba", "Cd", "Co", "Cs", "Cu", "Hg", "K", "Li", "Mg", "Mn", "Ni", "Pd",
            "Pt", "Rb", "Ru", "S_2", "S_R", "Sr", "Tc", "Ti",
        ];

        let missing: Vec<&str> = EMITTABLE
            .iter()
            .copied()
            .filter(|t| !is_parameterized(t))
            .collect();
        assert_eq!(missing, UNPARAMETERIZED);
        assert_eq!(EMITTABLE.len(), 56);
        assert_eq!(EMITTABLE.len() - UNPARAMETERIZED.len(), 34);
    }
}
