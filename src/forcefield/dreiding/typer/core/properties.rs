//! Defines the core chemical properties shared across the typer, including elements,
//! bond orders, and hybridization states recognized by the perception pipeline.
//!
//! This module centralizes all publicly exposed enums and helpers that other
//! subsystems use to reason about atomic identity, bonding multiplicity, and VSEPR
//! classifications. Keeping these definitions in one place ensures consistent
//! parsing and documentation across the crate.
//!
//! Upstream derives the parse errors with `thiserror` and provides `serde`
//! `Deserialize` impls so the rule deck can name elements and hybridizations as
//! strings in TOML. Beavyr vendors the typer with no TOML parser -- the deck is
//! native Rust data in `typing::rules_data` -- so the `Deserialize` impls are
//! gone and the error impls are written out by hand. `FromStr`, which the
//! `Deserialize` impls merely delegated to, is unchanged and still public.

use std::error::Error;
use std::fmt;
use std::str::FromStr;

/// Enumerates every element the typer understands along with its atomic number.
///
/// The variants are grouped by periodic trends (non-metals, alkali metals, etc.)
/// so that code consuming this API can rely on exhaustive matches while still
/// understanding the chemical context of each atom.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum Element {
    // --- Non-metals ---
    /// Hydrogen (H, Z = 1).
    H = 1,
    /// Helium (He, Z = 2).
    He = 2,
    /// Boron (B, Z = 5).
    B = 5,
    /// Carbon (C, Z = 6).
    C,
    /// Nitrogen (N, Z = 7).
    N,
    /// Oxygen (O, Z = 8).
    O,
    /// Fluorine (F, Z = 9).
    F,
    /// Neon (Ne, Z = 10).
    Ne,
    /// Silicon (Si, Z = 14).
    Si = 14,
    /// Phosphorus (P, Z = 15).
    P,
    /// Sulfur (S, Z = 16).
    S,
    /// Chlorine (Cl, Z = 17).
    Cl,
    /// Argon (Ar, Z = 18).
    Ar,
    /// Germanium (Ge, Z = 32).
    Ge = 32,
    /// Arsenic (As, Z = 33).
    As,
    /// Selenium (Se, Z = 34).
    Se,
    /// Bromine (Br, Z = 35).
    Br,
    /// Krypton (Kr, Z = 36).
    Kr,
    /// Antimony (Sb, Z = 51).
    Sb = 51,
    /// Tellurium (Te, Z = 52).
    Te,
    /// Iodine (I, Z = 53).
    I,
    /// Xenon (Xe, Z = 54).
    Xe,
    /// Astatine (At, Z = 85).
    At = 85,
    /// Radon (Rn, Z = 86).
    Rn,

    // --- Alkali Metals ---
    /// Lithium (Li, Z = 3).
    Li = 3,
    /// Sodium (Na, Z = 11).
    Na = 11,
    /// Potassium (K, Z = 19).
    K = 19,
    /// Rubidium (Rb, Z = 37).
    Rb = 37,
    /// Cesium (Cs, Z = 55).
    Cs = 55,
    /// Francium (Fr, Z = 87).
    Fr = 87,

    // --- Alkaline Earth Metals ---
    /// Beryllium (Be, Z = 4).
    Be = 4,
    /// Magnesium (Mg, Z = 12).
    Mg = 12,
    /// Calcium (Ca, Z = 20).
    Ca = 20,
    /// Strontium (Sr, Z = 38).
    Sr = 38,
    /// Barium (Ba, Z = 56).
    Ba = 56,
    /// Radium (Ra, Z = 88).
    Ra = 88,

    // --- Post-transition Metals ---
    /// Aluminum (Al, Z = 13).
    Al = 13,
    /// Gallium (Ga, Z = 31).
    Ga = 31,
    /// Indium (In, Z = 49).
    In = 49,
    /// Tin (Sn, Z = 50).
    Sn = 50,
    /// Thallium (Tl, Z = 81).
    Tl = 81,
    /// Lead (Pb, Z = 82).
    Pb,
    /// Bismuth (Bi, Z = 83).
    Bi,
    /// Polonium (Po, Z = 84).
    Po,

    // --- Transition Metals ---
    /// Scandium (Sc, Z = 21).
    Sc = 21,
    /// Titanium (Ti, Z = 22).
    Ti,
    /// Vanadium (V, Z = 23).
    V,
    /// Chromium (Cr, Z = 24).
    Cr,
    /// Manganese (Mn, Z = 25).
    Mn,
    /// Iron (Fe, Z = 26).
    Fe,
    /// Cobalt (Co, Z = 27).
    Co,
    /// Nickel (Ni, Z = 28).
    Ni,
    /// Copper (Cu, Z = 29).
    Cu,
    /// Zinc (Zn, Z = 30).
    Zn,
    /// Yttrium (Y, Z = 39).
    Y = 39,
    /// Zirconium (Zr, Z = 40).
    Zr,
    /// Niobium (Nb, Z = 41).
    Nb,
    /// Molybdenum (Mo, Z = 42).
    Mo,
    /// Technetium (Tc, Z = 43).
    Tc,
    /// Ruthenium (Ru, Z = 44).
    Ru,
    /// Rhodium (Rh, Z = 45).
    Rh,
    /// Palladium (Pd, Z = 46).
    Pd,
    /// Silver (Ag, Z = 47).
    Ag,
    /// Cadmium (Cd, Z = 48).
    Cd,
    /// Hafnium (Hf, Z = 72).
    Hf = 72,
    /// Tantalum (Ta, Z = 73).
    Ta,
    /// Tungsten (W, Z = 74).
    W,
    /// Rhenium (Re, Z = 75).
    Re,
    /// Osmium (Os, Z = 76).
    Os,
    /// Iridium (Ir, Z = 77).
    Ir,
    /// Platinum (Pt, Z = 78).
    Pt,
    /// Gold (Au, Z = 79).
    Au,
    /// Mercury (Hg, Z = 80).
    Hg,

    // --- Lanthanides ---
    /// Lanthanum (La, Z = 57).
    La = 57,
    /// Cerium (Ce, Z = 58).
    Ce,
    /// Praseodymium (Pr, Z = 59).
    Pr,
    /// Neodymium (Nd, Z = 60).
    Nd,
    /// Promethium (Pm, Z = 61).
    Pm,
    /// Samarium (Sm, Z = 62).
    Sm,
    /// Europium (Eu, Z = 63).
    Eu,
    /// Gadolinium (Gd, Z = 64).
    Gd,
    /// Terbium (Tb, Z = 65).
    Tb,
    /// Dysprosium (Dy, Z = 66).
    Dy,
    /// Holmium (Ho, Z = 67).
    Ho,
    /// Erbium (Er, Z = 68).
    Er,
    /// Thulium (Tm, Z = 69).
    Tm,
    /// Ytterbium (Yb, Z = 70).
    Yb,
    /// Lutetium (Lu, Z = 71).
    Lu,

    // --- Actinides ---
    /// Actinium (Ac, Z = 89).
    Ac = 89,
    /// Thorium (Th, Z = 90).
    Th,
    /// Protactinium (Pa, Z = 91).
    Pa,
    /// Uranium (U, Z = 92).
    U,
    /// Neptunium (Np, Z = 93).
    Np,
    /// Plutonium (Pu, Z = 94).
    Pu,
    /// Americium (Am, Z = 95).
    Am,
    /// Curium (Cm, Z = 96).
    Cm,
    /// Berkelium (Bk, Z = 97).
    Bk,
    /// Californium (Cf, Z = 98).
    Cf,
    /// Einsteinium (Es, Z = 99).
    Es,
    /// Fermium (Fm, Z = 100).
    Fm,
    /// Mendelevium (Md, Z = 101).
    Md,
    /// Nobelium (No, Z = 102).
    No,
    /// Lawrencium (Lr, Z = 103).
    Lr,

    // --- Superheavy Elements ---
    /// Rutherfordium (Rf, Z = 104).
    Rf = 104,
    /// Dubnium (Db, Z = 105).
    Db,
    /// Seaborgium (Sg, Z = 106).
    Sg,
    /// Bohrium (Bh, Z = 107).
    Bh,
    /// Hassium (Hs, Z = 108).
    Hs,
    /// Meitnerium (Mt, Z = 109).
    Mt,
    /// Darmstadtium (Ds, Z = 110).
    Ds,
    /// Roentgenium (Rg, Z = 111).
    Rg,
    /// Copernicium (Cn, Z = 112).
    Cn,
    /// Nihonium (Nh, Z = 113).
    Nh,
    /// Flerovium (Fl, Z = 114).
    Fl,
    /// Moscovium (Mc, Z = 115).
    Mc,
    /// Livermorium (Lv, Z = 116).
    Lv,
    /// Tennessine (Ts, Z = 117).
    Ts,
    /// Oganesson (Og, Z = 118).
    Og,
}

/// Error returned when parsing an unknown or misspelled element symbol.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseElementError(String);

impl fmt::Display for ParseElementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid element symbol: '{}'", self.0)
    }
}

impl Error for ParseElementError {}

impl FromStr for Element {
    type Err = ParseElementError;

    /// Parses an atomic symbol into an [`Element`] variant.
    ///
    /// The parser accepts standard IUPAC symbols (e.g., `"C"`, `"Mg"`) and
    /// maps them to the corresponding enum variant without case folding.
    ///
    /// # Errors
    ///
    /// Returns [`ParseElementError`] whenever the provided symbol is not listed in
    /// the periodic table definition above.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "H" => Ok(Self::H),
            "He" => Ok(Self::He),
            "Li" => Ok(Self::Li),
            "Be" => Ok(Self::Be),
            "B" => Ok(Self::B),
            "C" => Ok(Self::C),
            "N" => Ok(Self::N),
            "O" => Ok(Self::O),
            "F" => Ok(Self::F),
            "Ne" => Ok(Self::Ne),
            "Na" => Ok(Self::Na),
            "Mg" => Ok(Self::Mg),
            "Al" => Ok(Self::Al),
            "Si" => Ok(Self::Si),
            "P" => Ok(Self::P),
            "S" => Ok(Self::S),
            "Cl" => Ok(Self::Cl),
            "Ar" => Ok(Self::Ar),
            "K" => Ok(Self::K),
            "Ca" => Ok(Self::Ca),
            "Sc" => Ok(Self::Sc),
            "Ti" => Ok(Self::Ti),
            "V" => Ok(Self::V),
            "Cr" => Ok(Self::Cr),
            "Mn" => Ok(Self::Mn),
            "Fe" => Ok(Self::Fe),
            "Co" => Ok(Self::Co),
            "Ni" => Ok(Self::Ni),
            "Cu" => Ok(Self::Cu),
            "Zn" => Ok(Self::Zn),
            "Ga" => Ok(Self::Ga),
            "Ge" => Ok(Self::Ge),
            "As" => Ok(Self::As),
            "Se" => Ok(Self::Se),
            "Br" => Ok(Self::Br),
            "Kr" => Ok(Self::Kr),
            "Rb" => Ok(Self::Rb),
            "Sr" => Ok(Self::Sr),
            "Y" => Ok(Self::Y),
            "Zr" => Ok(Self::Zr),
            "Nb" => Ok(Self::Nb),
            "Mo" => Ok(Self::Mo),
            "Tc" => Ok(Self::Tc),
            "Ru" => Ok(Self::Ru),
            "Rh" => Ok(Self::Rh),
            "Pd" => Ok(Self::Pd),
            "Ag" => Ok(Self::Ag),
            "Cd" => Ok(Self::Cd),
            "In" => Ok(Self::In),
            "Sn" => Ok(Self::Sn),
            "Sb" => Ok(Self::Sb),
            "Te" => Ok(Self::Te),
            "I" => Ok(Self::I),
            "Xe" => Ok(Self::Xe),
            "Cs" => Ok(Self::Cs),
            "Ba" => Ok(Self::Ba),
            "La" => Ok(Self::La),
            "Ce" => Ok(Self::Ce),
            "Pr" => Ok(Self::Pr),
            "Nd" => Ok(Self::Nd),
            "Pm" => Ok(Self::Pm),
            "Sm" => Ok(Self::Sm),
            "Eu" => Ok(Self::Eu),
            "Gd" => Ok(Self::Gd),
            "Tb" => Ok(Self::Tb),
            "Dy" => Ok(Self::Dy),
            "Ho" => Ok(Self::Ho),
            "Er" => Ok(Self::Er),
            "Tm" => Ok(Self::Tm),
            "Yb" => Ok(Self::Yb),
            "Lu" => Ok(Self::Lu),
            "Hf" => Ok(Self::Hf),
            "Ta" => Ok(Self::Ta),
            "W" => Ok(Self::W),
            "Re" => Ok(Self::Re),
            "Os" => Ok(Self::Os),
            "Ir" => Ok(Self::Ir),
            "Pt" => Ok(Self::Pt),
            "Au" => Ok(Self::Au),
            "Hg" => Ok(Self::Hg),
            "Tl" => Ok(Self::Tl),
            "Pb" => Ok(Self::Pb),
            "Bi" => Ok(Self::Bi),
            "Po" => Ok(Self::Po),
            "At" => Ok(Self::At),
            "Rn" => Ok(Self::Rn),
            "Fr" => Ok(Self::Fr),
            "Ra" => Ok(Self::Ra),
            "Ac" => Ok(Self::Ac),
            "Th" => Ok(Self::Th),
            "Pa" => Ok(Self::Pa),
            "U" => Ok(Self::U),
            "Np" => Ok(Self::Np),
            "Pu" => Ok(Self::Pu),
            "Am" => Ok(Self::Am),
            "Cm" => Ok(Self::Cm),
            "Bk" => Ok(Self::Bk),
            "Cf" => Ok(Self::Cf),
            "Es" => Ok(Self::Es),
            "Fm" => Ok(Self::Fm),
            "Md" => Ok(Self::Md),
            "No" => Ok(Self::No),
            "Lr" => Ok(Self::Lr),
            "Rf" => Ok(Self::Rf),
            "Db" => Ok(Self::Db),
            "Sg" => Ok(Self::Sg),
            "Bh" => Ok(Self::Bh),
            "Hs" => Ok(Self::Hs),
            "Mt" => Ok(Self::Mt),
            "Ds" => Ok(Self::Ds),
            "Rg" => Ok(Self::Rg),
            "Cn" => Ok(Self::Cn),
            "Nh" => Ok(Self::Nh),
            "Fl" => Ok(Self::Fl),
            "Mc" => Ok(Self::Mc),
            "Lv" => Ok(Self::Lv),
            "Ts" => Ok(Self::Ts),
            "Og" => Ok(Self::Og),
            _ => Err(ParseElementError(s.to_string())),
        }
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Element {
    /// Returns the number of valence electrons for main-group elements.
    ///
    /// This helper covers Groups 1–18 where valence counts follow periodic trends
    /// and returns `None` for transition metals, lanthanides, and actinides whose
    /// electron configurations require specialized treatment.
    ///
    /// # Returns
    ///
    /// `Some(count)` with the valence electron total when it can be determined or
    /// `None` for elements outside the supported groups.
    pub fn valence_electrons(&self) -> Option<u8> {
        use Element::*;
        match self {
            // Group 1
            H | Li | Na | K | Rb | Cs | Fr => Some(1),
            // Group 2
            Be | Mg | Ca | Sr | Ba | Ra => Some(2),
            // Group 13
            B | Al | Ga | In | Tl => Some(3),
            // Group 14
            C | Si | Ge | Sn | Pb => Some(4),
            // Group 15
            N | P | As | Sb | Bi => Some(5),
            // Group 16
            O | S | Se | Te | Po => Some(6),
            // Group 17
            F | Cl | Br | I | At => Some(7),
            // Group 18 (Noble Gases)
            He | Ne | Ar | Kr | Xe | Rn => Some(8),
            // For transition metals, lanthanides, actinides, valence is complex.
            _ => None,
        }
    }
}


/// Describes the discrete bond multiplicities supported in input graphs.
///
/// The values intentionally match their valence contribution so that helpers like
/// electron counting can treat `Single` as 1, `Double` as 2, etc. This enum
/// includes aromatic bonds for input parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum GraphBondOrder {
    /// A single bond contributing one unit of valence.
    Single = 1,
    /// A double bond contributing two units of valence.
    Double = 2,
    /// A triple bond contributing three units of valence.
    Triple = 3,
    /// A Kekulé-aromatic bond that gets expanded during perception.
    Aromatic = 4,
}

/// Describes the discrete bond multiplicities used in output topologies.
///
/// The values intentionally match their valence contribution so that helpers like
/// electron counting can treat `Single` as 1, `Double` as 2, etc. This enum
/// uses resonant bonds instead of aromatic for physical determination.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum TopologyBondOrder {
    /// A single bond contributing one unit of valence.
    Single = 1,
    /// A double bond contributing two units of valence.
    Double = 2,
    /// A triple bond contributing three units of valence.
    Triple = 3,
    /// A resonant bond indicating participation in a delocalized π-system.
    Resonant = 4,
}

/// Error returned when parsing a bond order string that does not match the enum.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseBondOrderError(String);

impl fmt::Display for ParseBondOrderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid bond order: '{}'", self.0)
    }
}

impl Error for ParseBondOrderError {}

impl FromStr for GraphBondOrder {
    type Err = ParseBondOrderError;

    /// Parses a textual bond order into the [`GraphBondOrder`] enum.
    ///
    /// Accepts the canonical variant names (`"Single"`, `"Double"`, etc.).
    ///
    /// # Errors
    ///
    /// Returns [`ParseBondOrderError`] if the string is not one of the supported
    /// keywords.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Single" => Ok(Self::Single),
            "Double" => Ok(Self::Double),
            "Triple" => Ok(Self::Triple),
            "Aromatic" => Ok(Self::Aromatic),
            _ => Err(ParseBondOrderError(s.to_string())),
        }
    }
}

impl FromStr for TopologyBondOrder {
    type Err = ParseBondOrderError;

    /// Parses a textual bond order into the [`TopologyBondOrder`] enum.
    ///
    /// Accepts the canonical variant names (`"Single"`, `"Double"`, etc.).
    ///
    /// # Errors
    ///
    /// Returns [`ParseBondOrderError`] if the string is not one of the supported
    /// keywords.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Single" => Ok(Self::Single),
            "Double" => Ok(Self::Double),
            "Triple" => Ok(Self::Triple),
            "Resonant" => Ok(Self::Resonant),
            _ => Err(ParseBondOrderError(s.to_string())),
        }
    }
}

impl fmt::Display for GraphBondOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl fmt::Display for TopologyBondOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// Captures the VSEPR-derived hybridization states recognized by the typer.
///
/// These variants are used both as perception outputs and as serialized atom-type
/// hints when exporting topologies.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Hybridization {
    /// sp hybridization (linear geometry).
    SP,
    /// sp² hybridization (trigonal planar geometry).
    SP2,
    /// sp³ hybridization (tetrahedral geometry).
    SP3,
    /// Resonant hybridization, indicating participation in a delocalized π-system.
    Resonant,
    /// Used for atoms where hybridization is not typically considered (e.g., ions, halogens).
    None,
    /// An initial or error state before perception is complete.
    Unknown,
}

/// Error returned when parsing an unrecognized hybridization label.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseHybridizationError(String);

impl fmt::Display for ParseHybridizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid hybridization string: '{}'", self.0)
    }
}

impl Error for ParseHybridizationError {}

impl FromStr for Hybridization {
    type Err = ParseHybridizationError;

    /// Parses a serialized hybridization label into the [`Hybridization`] enum.
    ///
    /// The parser expects the exact variant names (`"SP"`, `"SP2"`, etc.) and
    /// surfaces descriptive errors if anything else is encountered.
    ///
    /// # Errors
    ///
    /// Returns [`ParseHybridizationError`] if the string does not map to a known
    /// hybridization state.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "SP" => Ok(Self::SP),
            "SP2" => Ok(Self::SP2),
            "SP3" => Ok(Self::SP3),
            "Resonant" => Ok(Self::Resonant),
            "None" => Ok(Self::None),
            "Unknown" => Ok(Self::Unknown),
            _ => Err(ParseHybridizationError(s.to_string())),
        }
    }
}

impl fmt::Display for Hybridization {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

