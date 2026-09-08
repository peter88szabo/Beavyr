//! The level of theory a calculation runs at, and whether it can run at all.
//!
//! Both job panels -- the optimizer and the frequency analysis -- ask for the
//! same thing: a method, whatever that method needs (a functional, a basis
//! set, a dispersion correction), and the resources to run it with. All of it
//! lives here, with no dependency on egui, so every rule below is a plain
//! function with a test rather than something only reachable by clicking.
//!
//! The catalogs are not guesses. Every name was checked against the programs
//! themselves -- `xtb --help`, `behemoth --list func|basis|disp|abinitio` --
//! and the design note in
//! `docs/superpowers/specs/2026-09-07-method-selection-design.md` records what
//! that turned up, including the several places the obvious spelling was wrong.

use crate::molecule::atomic_number;

use super::program::QcProgram;

/// The lightest element that needs an effective core potential in the def2
/// family: rubidium. Everything up to krypton, transition metals included, is
/// all-electron.
///
/// Behemoth implements no ECPs -- nothing in its help, its basis catalog or
/// its ab-initio list mentions them -- so a Gaussian-basis calculation on
/// anything this heavy cannot be run, rather than merely being inadvisable.
pub const FIRST_ECP_ELEMENT: u32 = 37;

/// The elements TASI is parameterised for, from the `atoms` block of
/// `parameters/tasi_chonclsf_parameters.json`. The file name spells them out.
pub const TASI_ELEMENTS: [&str; 7] = ["C", "Cl", "F", "H", "N", "O", "S"];

/// Which xTB parametrisation to run.
///
/// GFN0 is deliberately absent: `--gfn 0` is accepted by the command line but
/// dies with "Parameter file param_gfn0-xtb.txt not found!", because that file
/// is not shipped with the build in use. Adding it back is one variant here
/// and one line in `ALL`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XtbMethod {
    Gfn1,
    Gfn2,
    GfnFf,
}

impl XtbMethod {
    pub const ALL: [XtbMethod; 3] = [XtbMethod::Gfn1, XtbMethod::Gfn2, XtbMethod::GfnFf];

    pub fn label(self) -> &'static str {
        match self {
            XtbMethod::Gfn1 => "GFN1-xTB",
            XtbMethod::Gfn2 => "GFN2-xTB",
            // The force field's name is GFN-FF, and its flag is `--gfnff`.
            XtbMethod::GfnFf => "GFN-FF",
        }
    }
}

impl Default for XtbMethod {
    /// xTB's own default is `--gfn 2`.
    fn default() -> Self {
        XtbMethod::Gfn2
    }
}

/// Which Behemoth method to run.
///
/// `Dft` is one entry rather than two because the choice between Behemoth's
/// `rks` and `uks` is not the user's to make: `--method rks` on an open shell
/// is refused outright. The multiplicity decides, in `behemoth_method_flag`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BehemothMethod {
    Tasi,
    Gfn1Xtb,
    Hf,
    Dft,
    Mp2,
    RiMp2,
}

impl BehemothMethod {
    pub const ALL: [BehemothMethod; 6] = [
        BehemothMethod::Tasi,
        BehemothMethod::Gfn1Xtb,
        BehemothMethod::Hf,
        BehemothMethod::Dft,
        BehemothMethod::Mp2,
        BehemothMethod::RiMp2,
    ];

    pub fn label(self) -> &'static str {
        match self {
            BehemothMethod::Tasi => "TASI",
            BehemothMethod::Gfn1Xtb => "GFN1-xTB",
            BehemothMethod::Hf => "HF",
            BehemothMethod::Dft => "DFT",
            BehemothMethod::Mp2 => "MP2",
            BehemothMethod::RiMp2 => "RI-MP2",
        }
    }

    /// Whether this method expands the wavefunction in a Gaussian basis, and
    /// so needs core functions for every element it is given.
    ///
    /// The semi-empirical and parametrised methods do not: they carry their
    /// own element parameters, which is why the ECP rule does not apply to
    /// them and why TASI has an element rule of its own instead.
    pub fn uses_gaussian_basis(self) -> bool {
        match self {
            BehemothMethod::Tasi | BehemothMethod::Gfn1Xtb => false,
            BehemothMethod::Hf
            | BehemothMethod::Dft
            | BehemothMethod::Mp2
            | BehemothMethod::RiMp2 => true,
        }
    }

    /// Whether Behemoth requires a closed-shell reference for it.
    pub fn requires_closed_shell(self) -> bool {
        matches!(self, BehemothMethod::Mp2 | BehemothMethod::RiMp2)
    }
}

impl Default for BehemothMethod {
    /// Behemoth's own default, and what the panels already sent before there
    /// was anything to choose -- so an existing habit does not silently turn
    /// into an expensive DFT job.
    fn default() -> Self {
        BehemothMethod::Gfn1Xtb
    }
}

/// A dispersion correction applied after the SCF, or none at all.
///
/// From `behemoth --list disp`. Tkatchenko-Scheffler is listed there as "not
/// available yet" and so is not offered.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dispersion {
    None,
    D3Bj,
    D3Zero,
    D3Bjm,
    D3Zerom,
    D3Op,
    D3Cso,
    D3Z,
    D4,
}

impl Dispersion {
    pub const ALL: [Dispersion; 9] = [
        Dispersion::None,
        Dispersion::D3Bj,
        Dispersion::D3Zero,
        Dispersion::D3Bjm,
        Dispersion::D3Zerom,
        Dispersion::D3Op,
        Dispersion::D3Cso,
        Dispersion::D3Z,
        Dispersion::D4,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Dispersion::None => "None",
            Dispersion::D3Bj => "D3(BJ)",
            Dispersion::D3Zero => "D3(0)",
            Dispersion::D3Bjm => "D3(BJ,M)",
            Dispersion::D3Zerom => "D3(0,M)",
            Dispersion::D3Op => "D3(op)",
            Dispersion::D3Cso => "D3(CSO)",
            Dispersion::D3Z => "D3(Z)",
            Dispersion::D4 => "D4",
        }
    }

    /// The `--dispersion` keyword, or `None` when no flag should be passed.
    pub fn flag(self) -> Option<&'static str> {
        Some(match self {
            Dispersion::None => return None,
            Dispersion::D3Bj => "d3bj",
            Dispersion::D3Zero => "d3zero",
            Dispersion::D3Bjm => "d3bjm",
            Dispersion::D3Zerom => "d3zerom",
            Dispersion::D3Op => "d3op",
            Dispersion::D3Cso => "d3cso",
            Dispersion::D3Z => "d3z",
            Dispersion::D4 => "d4",
        })
    }
}

impl Default for Dispersion {
    fn default() -> Self {
        Dispersion::None
    }
}

/// How the Coulomb and exchange integrals are handled.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TwoElectron {
    Exact,
    Rijk,
}

impl TwoElectron {
    pub const ALL: [TwoElectron; 2] = [TwoElectron::Exact, TwoElectron::Rijk];

    pub fn label(self) -> &'static str {
        match self {
            TwoElectron::Exact => "exact",
            TwoElectron::Rijk => "RIJK",
        }
    }
}

impl Default for TwoElectron {
    fn default() -> Self {
        TwoElectron::Exact
    }
}

/// A DFT functional as the panel offers it.
pub struct Functional {
    /// What the dropdown shows.
    pub label: &'static str,
    /// What `--functional` is given. Behemoth is case-insensitive here, but
    /// the catalog stores the lowercase form its own listing uses.
    pub cli: &'static str,
    /// A composite carries its own basis set and dispersion correction, so
    /// choosing one takes both of those fields away.
    pub composite: Option<&'static str>,
    /// Whether it is a global hybrid, i.e. carries exact exchange.
    ///
    /// From the exact-exchange fractions `behemoth --list func` prints beside
    /// each hybrid. It matters because the sTDA/sTD-DFT engine refuses a pure
    /// functional outright: "sTDA/sTD-DFT with normal DFT requires a
    /// global-hybrid functional with nonzero exact exchange". The spectrum
    /// panel offers only the hybrids for that reason.
    pub hybrid: bool,
}

/// The nine common functionals in alphabetical order, then the five
/// composites.
///
/// Behemoth's real catalog is far larger and `--functional` also accepts raw
/// LibXC specifications; a dropdown of all of it would be unusable. Every name
/// here was checked to exist in `--list func`, which is where `m06-2x` came
/// from -- `m062x` is not a name Behemoth knows.
pub const FUNCTIONALS: [Functional; 14] = [
    Functional { label: "B3LYP", cli: "b3lyp", composite: None, hybrid: true },
    Functional { label: "B97", cli: "b97", composite: None, hybrid: true },
    Functional { label: "M06-2X", cli: "m06-2x", composite: None, hybrid: true },
    Functional { label: "MN15", cli: "mn15", composite: None, hybrid: true },
    Functional { label: "PBE", cli: "pbe", composite: None, hybrid: false },
    Functional { label: "PBE0", cli: "pbe0", composite: None, hybrid: true },
    Functional { label: "R2SCAN", cli: "r2scan", composite: None, hybrid: false },
    Functional { label: "revB3LYP", cli: "revb3lyp", composite: None, hybrid: true },
    Functional { label: "TPSSh", cli: "tpssh", composite: None, hybrid: true },
    Functional {
        label: "B3LYP-3c",
        cli: "b3lyp-3c",
        composite: Some("B3LYP (VWN5) + def2-mSVP + D3(BJ)+ATM + gCP"),
        hybrid: true,
    },
    Functional {
        label: "B97-3c",
        cli: "b97-3c",
        composite: Some("refitted mB97 GGA + def2-mTZVP + D3(BJ)+ATM + SRB"),
        // A GGA refit, with no exact exchange.
        hybrid: false,
    },
    Functional {
        label: "HF-3c",
        cli: "hf-3c",
        composite: Some("plain HF + MINIX + D3(BJ) + combined gCP/SRB"),
        // Hartree-Fock is all exact exchange.
        hybrid: true,
    },
    Functional {
        label: "PBEh-3c",
        cli: "pbeh-3c",
        composite: Some("PBE hybrid (42% HF) + def2-mSVP + D3(BJ)+ATM + damped gCP"),
        hybrid: true,
    },
    Functional {
        label: "r2SCAN-3c",
        cli: "r2scan-3c",
        composite: Some("r2SCAN + def2-mTZVPP + D4(BJ)+ATM + damped gCP"),
        // Unmodified r2SCAN, which is pure.
        hybrid: false,
    },
];

/// The default functional.
pub const DEFAULT_FUNCTIONAL: &str = "b3lyp";

/// HF-3c is the one composite that is not a functional at all. Behemoth's help
/// is explicit: "HF-3c is plain Hartree-Fock, not a DFT functional, so it is
/// requested as `--method hf-3c` rather than `--functional hf-3c`". It is
/// still listed with the other composites in the dropdown, because that is
/// where a user looks for it, and translated here.
pub const HF3C: &str = "hf-3c";

/// The catalog entry for a `--functional` name, if it is one we offer.
pub fn functional(cli: &str) -> Option<&'static Functional> {
    FUNCTIONALS.iter().find(|f| f.cli.eq_ignore_ascii_case(cli))
}

/// Whether a functional name is one of the composites, which brings its own
/// basis set and dispersion correction.
pub fn is_composite(cli: &str) -> bool {
    functional(cli).is_some_and(|f| f.composite.is_some())
}

/// Whether a functional carries exact exchange. Unknown names -- a raw LibXC
/// specification, say -- are taken to be hybrids, because refusing to run
/// something on a guess is worse than letting Behemoth give its own verdict.
pub fn is_hybrid(cli: &str) -> bool {
    functional(cli).map_or(true, |f| f.hybrid)
}

/// The elements in `atoms` that need an effective core potential, which
/// nothing implements. Shared with the excited-state panel, which has the same
/// rule for the same reason.
pub fn elements_needing_ecp(atoms: &[String]) -> Vec<String> {
    offending_elements(atoms, |symbol| {
        crate::molecule::atomic_number(symbol).is_some_and(|z| z >= FIRST_ECP_ELEMENT)
    })
}

/// The elements in `atoms` that TASI has no parameters for.
pub fn elements_outside_tasi(atoms: &[String]) -> Vec<String> {
    offending_elements(atoms, |symbol| {
        !TASI_ELEMENTS.iter().any(|e| e.eq_ignore_ascii_case(symbol))
    })
}

/// A group of basis sets, as the picker shows them.
pub struct BasisGroup {
    pub name: &'static str,
    /// Display name and command-line name. They differ where the command-line
    /// spelling is unreadable: `6-31G*` is passed as `6-31g_st`.
    pub sets: &'static [(&'static str, &'static str)],
}

/// The basis sets offered, grouped by family, from the 45 in `--list basis`.
///
/// Not offered, and why:
///
/// * `lanl2dz`, `lanl2dzdp`, `lanl2tz`, `lanl2tz+`, `lanl2tzf`, `mod-lanl2dz`
///   -- ECP basis sets, and there is no ECP implementation behind them.
/// * `minix`, `def2-msvp-pbeh3c`, `def2-mtzvp-b973c`, `def2-mtzvppd-r2scan3c`
///   -- each belongs to a composite and is selected by it automatically.
/// * every `-rifit`, `-jfit` and `-jkfit` set -- auxiliary bases that Behemoth
///   chooses on its own.
pub const BASIS_GROUPS: [BasisGroup; 5] = [
    BasisGroup {
        name: "def2",
        sets: &[
            ("def2-SV(P)", "def2-sv(p)"),
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
            ("aug-cc-pVDZ", "aug-cc-pvdz"),
            ("aug-cc-pVTZ", "aug-cc-pvtz"),
            ("aug-cc-pVQZ", "aug-cc-pvqz"),
        ],
    },
    BasisGroup {
        name: "STO-nG",
        sets: &[
            ("STO-3G", "sto-3g"),
            ("STO-4G", "sto-4g"),
            ("STO-5G", "sto-5g"),
            ("STO-6G", "sto-6g"),
        ],
    },
    BasisGroup {
        name: "Pople",
        sets: &[
            ("3-21G", "3-21g"),
            ("6-31G", "6-31g"),
            ("6-31G*", "6-31g_st"),
            ("6-31G**", "6-31g_st_st"),
        ],
    },
];

/// The default basis set.
pub const DEFAULT_BASIS: &str = "def2-svp";

/// The display name for a `--basis` name, falling back to the name itself so a
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

/// Everything a calculation needs beyond the geometry, the charge and the
/// multiplicity.
///
/// One per panel, not shared: running a Hessian at a different level than the
/// optimization is a legitimate thing to want, and is the same reasoning by
/// which each panel already keeps its own program selection.
#[derive(Debug, Clone, PartialEq)]
pub struct MethodConfig {
    pub xtb: XtbMethod,
    pub behemoth: BehemothMethod,
    pub functional: String,
    pub dispersion: Dispersion,
    pub basis: String,
    pub two_electron: TwoElectron,
    /// Behemoth's `--memory`, in MiB. It has no xTB equivalent.
    pub memory_mb: u32,
    /// Threads. Behemoth's `--nproc`, xTB's `-P`.
    pub nproc: u32,
}

impl Default for MethodConfig {
    fn default() -> Self {
        MethodConfig {
            xtb: XtbMethod::default(),
            behemoth: BehemothMethod::default(),
            functional: DEFAULT_FUNCTIONAL.to_string(),
            dispersion: Dispersion::default(),
            basis: DEFAULT_BASIS.to_string(),
            two_electron: TwoElectron::default(),
            // Behemoth's own default.
            memory_mb: 1024,
            // Explicit and reproducible. Behemoth with no `--nproc` takes
            // every thread the environment offers, which competes with
            // Beavyr's own rendering; choosing more is the user's call.
            nproc: 1,
        }
    }
}

/// Which fields apply to the current program and method.
///
/// Read by the panel and by its tests, so "is the basis field showing?" has
/// exactly one answer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VisibleFields {
    pub functional: bool,
    pub basis: bool,
    pub dispersion: bool,
    /// Shown for the methods it applies to, but never enabled in these
    /// panels: see `RIJK_UNAVAILABLE`.
    pub two_electron: bool,
    pub memory: bool,
    pub nproc: bool,
}

/// Why the RIJK control is shown but never enabled here.
///
/// Behemoth's help: RIJK "density-fits both Coulomb and exchange from one
/// shared auxiliary basis (energy only, no analytic gradient yet, so
/// `--gradient`/`--opt`/`--hessian`/`--freq`/`--scan`/`--md` are refused with
/// it)". Every job either panel runs is one of those, so the option is left
/// visible and disabled rather than offered and guaranteed to fail. It becomes
/// live for free if a single-point panel is ever added.
pub const RIJK_UNAVAILABLE: &str =
    "RIJK has no analytic gradient yet, so it cannot be used for optimization or frequencies.";

pub fn visible_fields(program: QcProgram, config: &MethodConfig) -> VisibleFields {
    let hidden = VisibleFields {
        functional: false,
        basis: false,
        dispersion: false,
        two_electron: false,
        memory: false,
        nproc: false,
    };
    match program {
        // xTB takes a parametrisation and a thread count. It has no memory
        // option at all, so the field is hidden rather than shown and
        // silently ignored.
        QcProgram::Xtb => VisibleFields { nproc: true, ..hidden },
        // DREIDING has no level of theory to choose. Its parameters come from the atom types,
        // which are perceived rather than selected, and it is single-threaded because a force
        // field evaluation is already far below the cost of spawning threads for it.
        QcProgram::Dreiding => hidden,
        QcProgram::Behemoth => {
            let base = VisibleFields { memory: true, nproc: true, ..hidden };
            match config.behemoth {
                BehemothMethod::Tasi | BehemothMethod::Gfn1Xtb => base,
                BehemothMethod::Hf | BehemothMethod::Mp2 | BehemothMethod::RiMp2 => {
                    VisibleFields { basis: true, ..base }
                }
                BehemothMethod::Dft => {
                    if is_composite(&config.functional) {
                        // A composite defines its own basis set and its own
                        // dispersion correction, so neither is the user's to
                        // set. What it brings is stated instead.
                        VisibleFields { functional: true, ..base }
                    } else {
                        VisibleFields {
                            functional: true,
                            basis: true,
                            dispersion: true,
                            two_electron: true,
                            ..base
                        }
                    }
                }
            }
        }
    }
}

/// How seriously to take a configuration issue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    /// Shown as a warning, and the Run button is disabled: this configuration
    /// cannot produce a result, so it is stopped before a process is spawned
    /// rather than after the program fails.
    Block,
    /// Shown as ordinary text. Something the panel decided on the user's
    /// behalf, with nothing to fix.
    Note,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MethodIssue {
    pub severity: Severity,
    pub message: String,
}

impl MethodIssue {
    pub fn block(message: impl Into<String>) -> Self {
        MethodIssue { severity: Severity::Block, message: message.into() }
    }

    pub fn note(message: impl Into<String>) -> Self {
        MethodIssue { severity: Severity::Note, message: message.into() }
    }
}

/// Every element in `atoms`, deduplicated and in the order first seen, that
/// `keep` selects. Used to name offenders in a message rather than saying
/// "some elements".
fn offending_elements(atoms: &[String], keep: impl Fn(&str) -> bool) -> Vec<String> {
    let mut found: Vec<String> = Vec::new();
    for atom in atoms {
        if keep(atom) && !found.iter().any(|s| s == atom) {
            found.push(atom.clone());
        }
    }
    found
}

/// Everything wrong with, or worth saying about, this configuration.
pub fn validate(
    program: QcProgram,
    config: &MethodConfig,
    multiplicity: i32,
    atoms: &[String],
) -> Vec<MethodIssue> {
    let mut issues = Vec::new();
    if program != QcProgram::Behemoth {
        // xTB takes no basis set and carries parameters for elements up to
        // radon, so none of the rules below apply to it.
        return issues;
    }

    if config.behemoth == BehemothMethod::Tasi {
        let outside = offending_elements(atoms, |symbol| {
            !TASI_ELEMENTS.iter().any(|e| e.eq_ignore_ascii_case(symbol))
        });
        if !outside.is_empty() {
            issues.push(MethodIssue::block(format!(
                "TASI has no parameters for {}. It is parameterised for {} only, \
                 so this calculation cannot be run.",
                outside.join(", "),
                TASI_ELEMENTS.join(", ")
            )));
        }
    }

    if config.behemoth.uses_gaussian_basis() {
        let heavy = offending_elements(atoms, |symbol| {
            atomic_number(symbol).is_some_and(|z| z >= FIRST_ECP_ELEMENT)
        });
        if !heavy.is_empty() {
            issues.push(MethodIssue::block(format!(
                "{} need an effective core potential -- every basis set beyond krypton \
                 is ECP-based -- and no ECP is implemented, so this calculation must not \
                 be run. Use a semi-empirical method, or a lighter model system.",
                heavy.join(", ")
            )));
        }
    }

    if config.behemoth.requires_closed_shell() && multiplicity > 1 {
        issues.push(MethodIssue::block(format!(
            "{} requires a closed-shell reference; multiplicity is {multiplicity}.",
            config.behemoth.label()
        )));
    }

    if config.behemoth == BehemothMethod::RiMp2 {
        issues.push(MethodIssue::note(
            "The correlation-fitting basis is chosen automatically -- there is nothing to set.",
        ));
    }

    if config.behemoth == BehemothMethod::Dft {
        if let Some(brings) = functional(&config.functional).and_then(|f| f.composite) {
            issues.push(MethodIssue::note(format!(
                "A composite method brings its own basis set and dispersion correction: {brings}.",
            )));
        }
    }

    // Unreachable through the panel, which never enables the control, but the
    // rule belongs where every other rule lives rather than resting on the UI
    // being the only guard.
    if visible_fields(program, config).two_electron && config.two_electron == TwoElectron::Rijk {
        issues.push(MethodIssue::block(RIJK_UNAVAILABLE));
    }

    issues
}

/// Whether anything blocks this configuration from running.
pub fn is_blocked(issues: &[MethodIssue]) -> bool {
    issues.iter().any(|i| i.severity == Severity::Block)
}

/// The `--method` value for Behemoth, given the multiplicity.
///
/// An open shell is always the unrestricted method. For DFT that is forced:
/// `--method rks` on a doublet is refused with "requires `--multi 1`; use
/// `--method uks`". For Hartree-Fock it is a deliberate redundancy -- `hf`
/// already chooses UHF on its own, verified on the water cation, where `hf`
/// and `uhf` give the same energy to ten decimals -- so that the command line
/// records which reference was used instead of implying it.
///
/// HF-3c is the exception, and takes no companion method name: it selects UHF
/// for an open shell itself.
pub fn behemoth_method_flag(config: &MethodConfig, multiplicity: i32) -> &'static str {
    let open_shell = multiplicity > 1;
    match config.behemoth {
        BehemothMethod::Tasi => "tasi",
        BehemothMethod::Gfn1Xtb => "xtb",
        BehemothMethod::Hf => {
            if open_shell {
                "uhf"
            } else {
                "hf"
            }
        }
        BehemothMethod::Mp2 => "mp2",
        BehemothMethod::RiMp2 => "ri-mp2",
        BehemothMethod::Dft => {
            if is_hf3c(&config.functional) {
                HF3C
            } else if open_shell {
                "uks"
            } else {
                "rks"
            }
        }
    }
}

fn is_hf3c(functional: &str) -> bool {
    functional.eq_ignore_ascii_case(HF3C)
}

/// The `--functional` value, or `None` when the method takes none.
///
/// HF-3c returns `None` because it travels as the method name instead.
pub fn behemoth_functional_flag(config: &MethodConfig) -> Option<&str> {
    if config.behemoth != BehemothMethod::Dft || is_hf3c(&config.functional) {
        return None;
    }
    Some(&config.functional)
}

/// The `--basis` value, or `None` when the method supplies its own.
pub fn behemoth_basis_flag(config: &MethodConfig) -> Option<&str> {
    match config.behemoth {
        BehemothMethod::Tasi | BehemothMethod::Gfn1Xtb => None,
        BehemothMethod::Hf | BehemothMethod::Mp2 | BehemothMethod::RiMp2 => Some(&config.basis),
        // Every composite, HF-3c included, pulls in its own published basis.
        BehemothMethod::Dft => (!is_composite(&config.functional)).then_some(&*config.basis),
    }
}

/// The `--dispersion` value, or `None` for no flag.
///
/// A composite carries its own, so nothing is passed alongside one even if the
/// field was set before the composite was chosen.
pub fn behemoth_dispersion_flag(config: &MethodConfig) -> Option<&'static str> {
    if config.behemoth != BehemothMethod::Dft || is_composite(&config.functional) {
        return None;
    }
    config.dispersion.flag()
}

/// The xTB flags for the chosen parametrisation: `--gfn 1`, `--gfn 2`, or the
/// standalone `--gfnff`.
pub fn xtb_method_args(config: &MethodConfig) -> Vec<String> {
    match config.xtb {
        XtbMethod::Gfn1 => vec!["--gfn".to_string(), "1".to_string()],
        XtbMethod::Gfn2 => vec!["--gfn".to_string(), "2".to_string()],
        XtbMethod::GfnFf => vec!["--gfnff".to_string()],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atoms(symbols: &[&str]) -> Vec<String> {
        symbols.iter().map(|s| s.to_string()).collect()
    }

    fn behemoth(method: BehemothMethod) -> MethodConfig {
        MethodConfig { behemoth: method, ..Default::default() }
    }

    #[test]
    fn the_defaults_are_the_ones_asked_for() {
        let config = MethodConfig::default();
        assert_eq!(config.xtb, XtbMethod::Gfn2, "GFN2 is xTB's default");
        assert_eq!(config.functional, "b3lyp");
        assert_eq!(config.basis, "def2-svp");
        assert_eq!(config.dispersion, Dispersion::None);
        assert_eq!(config.memory_mb, 1024, "1 GB");
        assert_eq!(config.nproc, 1);
        assert_eq!(config.two_electron, TwoElectron::Exact);
    }

    /// GFN0 is not offered: the parameter file it needs is not shipped, and
    /// `--gfn 0` dies rather than falling back.
    #[test]
    fn the_xtb_list_is_the_three_that_run() {
        let labels: Vec<&str> = XtbMethod::ALL.iter().map(|m| m.label()).collect();
        assert_eq!(labels, vec!["GFN1-xTB", "GFN2-xTB", "GFN-FF"]);
    }

    #[test]
    fn xtb_takes_a_parametrisation_and_threads_and_nothing_else() {
        let fields = visible_fields(QcProgram::Xtb, &MethodConfig::default());
        assert!(fields.nproc);
        assert!(!fields.memory, "xTB has no memory flag");
        assert!(!fields.basis);
        assert!(!fields.functional);
        assert!(!fields.dispersion);
        assert!(!fields.two_electron);
    }

    #[test]
    fn a_semi_empirical_behemoth_method_needs_no_basis() {
        for method in [BehemothMethod::Tasi, BehemothMethod::Gfn1Xtb] {
            let fields = visible_fields(QcProgram::Behemoth, &behemoth(method));
            assert!(!fields.basis, "{method:?}");
            assert!(!fields.functional, "{method:?}");
            assert!(fields.memory && fields.nproc, "{method:?}");
        }
    }

    #[test]
    fn a_wavefunction_method_shows_only_the_basis() {
        for method in [BehemothMethod::Hf, BehemothMethod::Mp2, BehemothMethod::RiMp2] {
            let fields = visible_fields(QcProgram::Behemoth, &behemoth(method));
            assert!(fields.basis, "{method:?}");
            assert!(!fields.functional, "{method:?}");
            assert!(!fields.dispersion, "{method:?}");
        }
    }

    #[test]
    fn plain_dft_shows_the_functional_basis_dispersion_and_two_electron() {
        let fields = visible_fields(QcProgram::Behemoth, &behemoth(BehemothMethod::Dft));
        assert!(fields.functional);
        assert!(fields.basis);
        assert!(fields.dispersion);
        assert!(fields.two_electron);
    }

    /// The requested behaviour: a composite defines its own basis set and
    /// dispersion correction, so both fields go away.
    #[test]
    fn a_composite_hides_the_basis_and_dispersion_fields() {
        for name in ["hf-3c", "r2scan-3c", "pbeh-3c", "b97-3c", "b3lyp-3c"] {
            let config = MethodConfig {
                behemoth: BehemothMethod::Dft,
                functional: name.to_string(),
                ..Default::default()
            };
            let fields = visible_fields(QcProgram::Behemoth, &config);
            assert!(fields.functional, "{name}");
            assert!(!fields.basis, "{name} brings its own basis");
            assert!(!fields.dispersion, "{name} brings its own dispersion");
        }
    }

    #[test]
    fn the_functional_list_is_alphabetical_then_composite() {
        let plain: Vec<&str> = FUNCTIONALS
            .iter()
            .filter(|f| f.composite.is_none())
            .map(|f| f.label)
            .collect();
        let mut sorted = plain.clone();
        sorted.sort_by_key(|s| s.to_lowercase());
        assert_eq!(plain, sorted, "the common functionals are in abc order");

        let first_composite = FUNCTIONALS.iter().position(|f| f.composite.is_some()).unwrap();
        assert!(
            FUNCTIONALS[first_composite..].iter().all(|f| f.composite.is_some()),
            "every composite comes after every plain functional"
        );
        assert_eq!(plain.len(), 9);
        assert_eq!(FUNCTIONALS.len() - plain.len(), 5, "five composites");
    }

    /// A display name that does not match the command-line name is exactly the
    /// bug this catches: `M06-2X` is `m06-2x`, and `m062x` is not a name
    /// Behemoth knows.
    #[test]
    fn every_functional_name_is_lowercase_as_behemoth_lists_it() {
        for f in &FUNCTIONALS {
            assert_eq!(f.cli, f.cli.to_lowercase(), "{}", f.label);
            assert!(!f.cli.is_empty());
        }
        assert_eq!(functional("M06-2X").unwrap().cli, "m06-2x");
        assert!(functional("m062x").is_none(), "m062x is not a Behemoth name");
    }

    #[test]
    fn the_basis_groups_are_the_families_asked_for() {
        let names: Vec<&str> = BASIS_GROUPS.iter().map(|g| g.name).collect();
        assert_eq!(names, vec!["def2", "pc", "cc", "STO-nG", "Pople"]);
        // The default is in the catalog and is reachable from the picker.
        assert_eq!(basis_label(DEFAULT_BASIS), "def2-SVP");
    }

    /// The star spellings are the ones that would silently fail if the display
    /// name reached the command line.
    #[test]
    fn a_basis_display_name_is_not_its_command_line_name() {
        assert_eq!(basis_label("6-31g_st"), "6-31G*");
        assert_eq!(basis_label("6-31g_st_st"), "6-31G**");
        // Anything hand-typed still shows as itself.
        assert_eq!(basis_label("def2-mtzvp-b973c"), "def2-mtzvp-b973c");
    }

    /// No ECP basis set is offered, because no ECP is implemented.
    #[test]
    fn the_basis_catalog_offers_nothing_that_needs_an_ecp_or_is_auxiliary() {
        for group in &BASIS_GROUPS {
            for (label, cli) in group.sets {
                assert!(!cli.contains("lanl"), "{label} is an ECP basis");
                assert!(!cli.contains("fit"), "{label} is an auxiliary basis");
                assert!(!cli.contains("3c"), "{label} belongs to a composite");
                assert_ne!(*cli, "minix", "MINIX belongs to HF-3c");
            }
        }
    }

    #[test]
    fn an_ordinary_organic_molecule_under_dft_is_fine() {
        let issues = validate(
            QcProgram::Behemoth,
            &behemoth(BehemothMethod::Dft),
            1,
            &atoms(&["C", "H", "H", "H", "O", "H"]),
        );
        assert!(!is_blocked(&issues), "{issues:?}");
    }

    /// The requested rule: a heavy element needs an ECP, none is implemented,
    /// so the calculation must not run.
    #[test]
    fn a_heavy_element_blocks_a_gaussian_basis_calculation() {
        for method in [
            BehemothMethod::Hf,
            BehemothMethod::Dft,
            BehemothMethod::Mp2,
            BehemothMethod::RiMp2,
        ] {
            let issues = validate(
                QcProgram::Behemoth,
                &behemoth(method),
                1,
                &atoms(&["Rb", "Cl"]),
            );
            assert!(is_blocked(&issues), "{method:?} should be blocked");
            assert!(
                issues.iter().any(|i| i.message.contains("Rb")),
                "{method:?}: the message names the element"
            );
        }
    }

    /// The threshold is where def2 switches to ECPs, so transition metals --
    /// which are all-electron there -- must not be caught.
    #[test]
    fn a_transition_metal_is_all_electron_and_not_blocked() {
        for symbol in ["Fe", "Cu", "Zn", "Kr", "Br"] {
            let issues = validate(
                QcProgram::Behemoth,
                &behemoth(BehemothMethod::Dft),
                1,
                &atoms(&[symbol, "C"]),
            );
            assert!(!is_blocked(&issues), "{symbol} is all-electron in def2: {issues:?}");
        }
        // Rubidium is the first that is not.
        assert_eq!(atomic_number("Rb"), Some(FIRST_ECP_ELEMENT));
    }

    /// The parametrised methods carry their own elements and are exempt.
    #[test]
    fn a_heavy_element_does_not_block_a_parametrised_method() {
        for method in [BehemothMethod::Gfn1Xtb] {
            let issues =
                validate(QcProgram::Behemoth, &behemoth(method), 1, &atoms(&["Pb", "H"]));
            assert!(!is_blocked(&issues), "{method:?}: {issues:?}");
        }
        // And xTB the program never gets these rules at all.
        let issues = validate(QcProgram::Xtb, &MethodConfig::default(), 1, &atoms(&["Pb", "H"]));
        assert!(issues.is_empty());
    }

    /// The user's own addition: TASI is parameterised for seven elements, and
    /// anything else must be blocked rather than run.
    #[test]
    fn tasi_blocks_an_element_it_has_no_parameters_for() {
        let issues = validate(
            QcProgram::Behemoth,
            &behemoth(BehemothMethod::Tasi),
            1,
            &atoms(&["C", "H", "P", "Br"]),
        );
        assert!(is_blocked(&issues));
        let message = &issues[0].message;
        assert!(message.contains('P') && message.contains("Br"), "{message}");
    }

    #[test]
    fn tasi_runs_on_the_elements_it_does_have() {
        let issues = validate(
            QcProgram::Behemoth,
            &behemoth(BehemothMethod::Tasi),
            1,
            &atoms(&["C", "H", "N", "O", "F", "S", "Cl"]),
        );
        assert!(!is_blocked(&issues), "{issues:?}");
    }

    #[test]
    fn mp2_is_blocked_on_an_open_shell() {
        for method in [BehemothMethod::Mp2, BehemothMethod::RiMp2] {
            let triplet = validate(QcProgram::Behemoth, &behemoth(method), 3, &atoms(&["C"]));
            assert!(is_blocked(&triplet), "{method:?} needs a closed shell");
            let singlet = validate(QcProgram::Behemoth, &behemoth(method), 1, &atoms(&["C"]));
            assert!(!is_blocked(&singlet), "{method:?} at multiplicity 1: {singlet:?}");
        }
    }

    #[test]
    fn ri_mp2_says_the_fitting_basis_is_automatic() {
        let issues =
            validate(QcProgram::Behemoth, &behemoth(BehemothMethod::RiMp2), 1, &atoms(&["C"]));
        assert!(issues
            .iter()
            .any(|i| i.severity == Severity::Note && i.message.contains("automatically")));
    }

    #[test]
    fn a_composite_says_what_it_brings() {
        let config = MethodConfig {
            behemoth: BehemothMethod::Dft,
            functional: "r2scan-3c".to_string(),
            ..Default::default()
        };
        let issues = validate(QcProgram::Behemoth, &config, 1, &atoms(&["C"]));
        assert!(!is_blocked(&issues));
        assert!(issues.iter().any(|i| i.message.contains("def2-mTZVPP")), "{issues:?}");
    }

    /// An open shell is always unrestricted, as asked.
    #[test]
    fn an_open_shell_is_always_unrestricted() {
        let hf = behemoth(BehemothMethod::Hf);
        assert_eq!(behemoth_method_flag(&hf, 1), "hf");
        assert_eq!(behemoth_method_flag(&hf, 2), "uhf");
        assert_eq!(behemoth_method_flag(&hf, 3), "uhf");

        let dft = behemoth(BehemothMethod::Dft);
        assert_eq!(behemoth_method_flag(&dft, 1), "rks");
        assert_eq!(behemoth_method_flag(&dft, 2), "uks");
        assert_eq!(behemoth_method_flag(&dft, 3), "uks");
    }

    /// HF-3c travels as the method, not the functional -- Behemoth refuses it
    /// as a `--functional` -- and needs no companion method name for an open
    /// shell, choosing UHF itself.
    #[test]
    fn hf3c_is_a_method_and_carries_no_functional_or_basis() {
        let config = MethodConfig {
            behemoth: BehemothMethod::Dft,
            functional: "hf-3c".to_string(),
            ..Default::default()
        };
        assert_eq!(behemoth_method_flag(&config, 1), "hf-3c");
        assert_eq!(behemoth_method_flag(&config, 3), "hf-3c");
        assert_eq!(behemoth_functional_flag(&config), None);
        assert_eq!(behemoth_basis_flag(&config), None);
        assert_eq!(behemoth_dispersion_flag(&config), None);
    }

    #[test]
    fn a_plain_functional_carries_its_functional_basis_and_dispersion() {
        let config = MethodConfig {
            behemoth: BehemothMethod::Dft,
            functional: "pbe0".to_string(),
            dispersion: Dispersion::D4,
            basis: "def2-tzvp".to_string(),
            ..Default::default()
        };
        assert_eq!(behemoth_method_flag(&config, 1), "rks");
        assert_eq!(behemoth_functional_flag(&config), Some("pbe0"));
        assert_eq!(behemoth_basis_flag(&config), Some("def2-tzvp"));
        assert_eq!(behemoth_dispersion_flag(&config), Some("d4"));
    }

    /// A dispersion correction left over from a previous choice must not be
    /// smuggled in alongside a composite, which has its own.
    #[test]
    fn a_composite_passes_no_basis_or_dispersion_even_if_they_were_set() {
        let config = MethodConfig {
            behemoth: BehemothMethod::Dft,
            functional: "b97-3c".to_string(),
            dispersion: Dispersion::D3Bj,
            basis: "cc-pvtz".to_string(),
            ..Default::default()
        };
        assert_eq!(behemoth_functional_flag(&config), Some("b97-3c"));
        assert_eq!(behemoth_basis_flag(&config), None);
        assert_eq!(behemoth_dispersion_flag(&config), None);
    }

    #[test]
    fn a_semi_empirical_method_passes_no_basis_functional_or_dispersion() {
        for method in [BehemothMethod::Tasi, BehemothMethod::Gfn1Xtb] {
            let config = behemoth(method);
            assert_eq!(behemoth_basis_flag(&config), None, "{method:?}");
            assert_eq!(behemoth_functional_flag(&config), None, "{method:?}");
            assert_eq!(behemoth_dispersion_flag(&config), None, "{method:?}");
        }
        assert_eq!(behemoth_method_flag(&behemoth(BehemothMethod::Tasi), 1), "tasi");
        assert_eq!(behemoth_method_flag(&behemoth(BehemothMethod::Gfn1Xtb), 1), "xtb");
    }

    /// Dispersion belongs to DFT; nothing else takes the flag.
    #[test]
    fn only_dft_passes_a_dispersion_correction() {
        for method in [BehemothMethod::Hf, BehemothMethod::Mp2, BehemothMethod::RiMp2] {
            let config = MethodConfig {
                behemoth: method,
                dispersion: Dispersion::D3Bj,
                ..Default::default()
            };
            assert_eq!(behemoth_dispersion_flag(&config), None, "{method:?}");
        }
    }

    #[test]
    fn the_xtb_flags_are_gfn_levels_and_the_force_field() {
        let gfn = |m| xtb_method_args(&MethodConfig { xtb: m, ..Default::default() });
        assert_eq!(gfn(XtbMethod::Gfn1), vec!["--gfn", "1"]);
        assert_eq!(gfn(XtbMethod::Gfn2), vec!["--gfn", "2"]);
        // The force field is a standalone flag, not a `--gfn` level.
        assert_eq!(gfn(XtbMethod::GfnFf), vec!["--gfnff"]);
    }

    #[test]
    fn no_dispersion_means_no_flag() {
        assert_eq!(Dispersion::None.flag(), None);
        for d in Dispersion::ALL.iter().filter(|d| **d != Dispersion::None) {
            let flag = d.flag().expect("a keyword");
            assert_eq!(flag, flag.to_lowercase(), "{}", d.label());
        }
    }
}
