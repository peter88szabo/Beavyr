//! What excited-state calculation to run, and whether it can run.
//!
//! Behemoth only. xTB has no excited-state module, and ORCA is not driven by
//! Beavyr -- its output is read, never produced.
//!
//! Shaped like `qchem_interfaces::method`, deliberately: same split between a
//! config struct, a `visible_fields`, a `validate` and the command
//! construction, so there is one pattern for this in the codebase rather than
//! two. The functional, basis set, memory and thread count are not duplicated
//! here -- they come from the existing `MethodConfig`, which already holds
//! exactly those.
//!
//! Everything below was checked against the binary and against Behemoth's own
//! `src/utils/cli/` and `src/hamiltonian/stddft/`; the design note in
//! `docs/superpowers/specs/2026-09-07-excited-state-runs-design.md` records
//! what that turned up.

use std::path::Path;
use std::process::Command;

use crate::qchem_interfaces::method::{
    behemoth_basis_flag, elements_needing_ecp, elements_outside_tasi, is_hybrid, Functional,
    MethodConfig, MethodIssue, FUNCTIONALS, TASI_ELEMENTS,
};

/// Which route into Behemoth's sTDA/sTD-DFT engine to take.
///
/// `--stda` is the engine's Tamm-Dancoff mode and `--stddft` its coupled one;
/// which orbitals it builds on follows from `--method`, which Behemoth's own
/// CLI tests show it infers when `--stda-orbitals` is not given.
///
/// Behemoth's other excited-state engine -- `--rpa` / `--tda`, the RPA-xTB
/// response -- is not here. It is a different response model with its own
/// parameterisation and an active-space energy window in place of a root
/// count, so it needs its own controls rather than being bent into these.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExcitedStateMethod {
    /// Grimme's original sTDA-xTB, `--method xtb --stda`. Needs an external
    /// `xtb4stda` binary to generate its orbitals.
    StdaXtbOriginal,
    /// sTDA on Behemoth's own GFN1 orbitals,
    /// `--method xtb --stda --stda-orbitals gfn1`. No external binary.
    StdaXtbGfn1,
    /// sTDA on TASI orbitals, `--method tasi --stda`. No external binary, but
    /// TASI's seven elements only.
    StdaTasi,
    /// sTDA or sTD-DFT on a hybrid DFT reference.
    Tddft,
}

impl ExcitedStateMethod {
    pub const ALL: [ExcitedStateMethod; 4] = [
        ExcitedStateMethod::StdaXtbOriginal,
        ExcitedStateMethod::StdaXtbGfn1,
        ExcitedStateMethod::StdaTasi,
        ExcitedStateMethod::Tddft,
    ];

    pub fn label(self) -> &'static str {
        match self {
            ExcitedStateMethod::StdaXtbOriginal => "sTDA-xTB (Grimme's original)",
            ExcitedStateMethod::StdaXtbGfn1 => "sTDA-xTB (GFN1 orbitals)",
            ExcitedStateMethod::StdaTasi => "sTDA (TASI)",
            ExcitedStateMethod::Tddft => "TD-DFT",
        }
    }

    /// A one-line explanation for under the dropdown, since the difference
    /// between these is not obvious from their names.
    pub fn description(self) -> &'static str {
        match self {
            ExcitedStateMethod::StdaXtbOriginal => {
                "The original sTDA-xTB. Needs the external xtb4stda orbital generator."
            }
            ExcitedStateMethod::StdaXtbGfn1 => {
                "sTDA on Behemoth's own GFN1 orbitals. No external program needed."
            }
            ExcitedStateMethod::StdaTasi => {
                "sTDA on TASI orbitals. No external program needed; C, H, N, O, F, S, Cl only."
            }
            ExcitedStateMethod::Tddft => {
                "Simplified TD-DFT on a hybrid reference. Choose the functional and basis set."
            }
        }
    }

    /// The `--method` this route needs.
    fn behemoth_method(self, multiplicity: i32) -> &'static str {
        match self {
            ExcitedStateMethod::StdaXtbOriginal | ExcitedStateMethod::StdaXtbGfn1 => "xtb",
            ExcitedStateMethod::StdaTasi => "tasi",
            // An open shell is always the unrestricted method: `rks` is
            // refused above a singlet outright.
            ExcitedStateMethod::Tddft => {
                if multiplicity > 1 {
                    "uks"
                } else {
                    "rks"
                }
            }
        }
    }

    /// Whether the Tamm-Dancoff switch means anything for this route.
    ///
    /// Only TD-DFT can turn it off. `--stddft` "currently requires a
    /// global-hybrid `--method rks` or `--method uks` reference", so every
    /// other route is Tamm-Dancoff by construction and the switch is hidden
    /// for them rather than shown and ignored.
    pub fn supports_full_response(self) -> bool {
        self == ExcitedStateMethod::Tddft
    }
}

impl Default for ExcitedStateMethod {
    /// The route that needs no external binary and no basis set, and covers
    /// most organic chromophores.
    fn default() -> Self {
        ExcitedStateMethod::StdaTasi
    }
}

/// The functionals the spectrum panel offers.
///
/// Only global hybrids, and no composites. The engine refuses a pure
/// functional -- "sTDA/sTD-DFT with normal DFT requires a global-hybrid
/// functional with nonzero exact exchange" -- which rules out PBE and R2SCAN.
/// The composites are left out because two of the five are not hybrids either
/// and all of them carry their own basis set, which contradicts a panel that
/// asks for one.
pub fn spectrum_functionals() -> Vec<&'static Functional> {
    FUNCTIONALS
        .iter()
        .filter(|f| f.hybrid && f.composite.is_none())
        .collect()
}

/// The default functional for a spectrum. A hybrid, necessarily.
pub const DEFAULT_SPECTRUM_FUNCTIONAL: &str = "b3lyp";

/// The canonical molecular orbitals `--wf2molden` writes.
///
/// Canonical rather than natural: the excitation assignments are given as
/// canonical MO indices ("8 HOMO ==> 9 LUMO"), so these are the orbitals those
/// numbers refer to. `naturalMO.molden` is written in the same run and is the
/// right file for a correlated density, which is not what this panel shows.
pub const CANONICAL_MOLDEN_FILE: &str = "canonicalMO.molden";

/// Everything the excited-state calculation needs beyond the reference's own
/// level of theory.
#[derive(Debug, Clone, PartialEq)]
pub struct ExcitedStateConfig {
    pub method: ExcitedStateMethod,
    /// `--stddft-roots`. Behemoth's default is 5.
    pub roots: u32,
    /// `--stddft-emax`, in eV. Behemoth's default is 10.
    ///
    /// Not asked for, and included anyway: roots outside this window are not
    /// reported at all, so without it the root count silently does not mean
    /// what it says. Asking for six roots on formaldehyde at the default
    /// window returns three.
    pub emax_ev: f64,
    /// Tamm-Dancoff, i.e. `--stda` rather than `--stddft`. On by default.
    pub tda: bool,
    /// Where the external `xtb4stda` generator lives, for the original
    /// sTDA-xTB route. Empty until the user points at one.
    pub xtb4stda_path: String,
}

impl Default for ExcitedStateConfig {
    fn default() -> Self {
        ExcitedStateConfig {
            method: ExcitedStateMethod::default(),
            roots: 5,
            emax_ev: 10.0,
            tda: true,
            xtb4stda_path: String::new(),
        }
    }
}

/// Which fields the current method needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VisibleFields {
    pub functional: bool,
    pub basis: bool,
    pub tda: bool,
    pub xtb4stda_path: bool,
}

pub fn visible_fields(config: &ExcitedStateConfig) -> VisibleFields {
    let tddft = config.method == ExcitedStateMethod::Tddft;
    VisibleFields {
        functional: tddft,
        basis: tddft,
        tda: config.method.supports_full_response(),
        xtb4stda_path: config.method == ExcitedStateMethod::StdaXtbOriginal,
    }
}

/// What the panel says about the GFN1 route, which ran but was numerically
/// unstable on every molecule tried.
pub const GFN1_INSTABILITY: &str =
    "GFN1 orbitals were unstable on both molecules tested here (\"sTDA A' matrix is unstable\"). \
     It needs no external program, so it is worth a try, but expect that error.";

/// Everything wrong with, or worth saying about, this configuration.
pub fn validate(
    config: &ExcitedStateConfig,
    reference: &MethodConfig,
    multiplicity: i32,
    atoms: &[String],
) -> Vec<MethodIssue> {
    let mut issues = Vec::new();

    match config.method {
        ExcitedStateMethod::StdaXtbOriginal => {
            if config.xtb4stda_path.trim().is_empty() {
                issues.push(MethodIssue::block(
                    "The original sTDA-xTB needs the external xtb4stda orbital generator. \
                     Give its path, or pick a route that does not need one.",
                ));
            }
        }
        ExcitedStateMethod::StdaXtbGfn1 => issues.push(MethodIssue::note(GFN1_INSTABILITY)),
        ExcitedStateMethod::StdaTasi => {
            let outside = elements_outside_tasi(atoms);
            if !outside.is_empty() {
                issues.push(MethodIssue::block(format!(
                    "TASI has no parameters for {}. It is parameterised for {} only, \
                     so this calculation cannot be run.",
                    outside.join(", "),
                    TASI_ELEMENTS.join(", ")
                )));
            }
        }
        ExcitedStateMethod::Tddft => {
            let heavy = elements_needing_ecp(atoms);
            if !heavy.is_empty() {
                issues.push(MethodIssue::block(format!(
                    "{} need an effective core potential -- every basis set beyond krypton \
                     is ECP-based -- and no ECP is implemented, so this calculation must not \
                     be run.",
                    heavy.join(", ")
                )));
            }
            // Unreachable through the dropdown, which offers only hybrids, but
            // the rule belongs where the others live rather than resting on
            // the UI being the only guard.
            if !is_hybrid(&reference.functional) {
                issues.push(MethodIssue::block(format!(
                    "{} has no exact exchange. The sTDA/sTD-DFT engine requires a \
                     global-hybrid functional.",
                    reference.functional
                )));
            }
        }
    }

    if !config.tda && !config.method.supports_full_response() {
        issues.push(MethodIssue::note(
            "The coupled response needs a hybrid DFT reference; this route is Tamm-Dancoff only.",
        ));
    }
    if config.roots == 0 {
        issues.push(MethodIssue::block("Ask for at least one root."));
    }
    // The root count is bounded by the window, which is worth saying before a
    // run rather than after it comes back short.
    issues.push(MethodIssue::note(format!(
        "Roots above {} eV are not reported, however many are requested.",
        format_ev(config.emax_ev)
    )));

    if multiplicity > 1 && config.method == ExcitedStateMethod::Tddft {
        issues.push(MethodIssue::note(
            "An open-shell reference runs as UKS.",
        ));
    }

    issues
}

/// An energy without trailing zeros, so a note reads "10 eV" not "10.00 eV".
fn format_ev(ev: f64) -> String {
    if (ev - ev.round()).abs() < 1e-9 {
        format!("{}", ev.round() as i64)
    } else {
        format!("{ev:.2}")
    }
}

/// The command that runs the spectrum.
pub fn spectrum_command(
    binary: &Path,
    run_dir: &Path,
    xyz_file: &str,
    charge: i32,
    multiplicity: i32,
    config: &ExcitedStateConfig,
    reference: &MethodConfig,
) -> Command {
    let multiplicity = multiplicity.max(1);
    let mut command = Command::new(binary);
    command
        .current_dir(run_dir)
        .arg("--xyz")
        .arg(xyz_file)
        .arg("--method")
        .arg(config.method.behemoth_method(multiplicity))
        .arg("--charge")
        .arg(charge.to_string())
        .arg("--multi")
        .arg(multiplicity.to_string());

    if config.method == ExcitedStateMethod::Tddft {
        command.arg("--functional").arg(&reference.functional);
        // The ground-state orbitals the excitations are expressed in, written
        // as Molden in the same run -- verified to work alongside `--stda`, so
        // no second calculation is needed.
        //
        // TD-DFT only, and not merely because `--wf2molden` is restricted to
        // the Gaussian-basis references. A semi-empirical orbital set is
        // parameterised to reproduce *transitions*, not the ground-state
        // electronic structure, so its orbitals are not a description of the
        // ground state and showing them as one would be wrong. That holds even
        // where such a set could be written, which is why this is a decision
        // rather than a limitation.
        command.arg("--wf2molden");
        // Reuses the reference's own rule, so a basis is passed exactly when
        // the method takes one.
        if let Some(basis) = behemoth_basis_flag(&MethodConfig {
            behemoth: crate::qchem_interfaces::method::BehemothMethod::Dft,
            ..reference.clone()
        }) {
            command.arg("--basis").arg(basis);
        }
    }

    // GFN1 orbitals have to be asked for; the original generator and TASI are
    // what `--method` already implies, and naming them again is redundant.
    if config.method == ExcitedStateMethod::StdaXtbGfn1 {
        command.arg("--stda-orbitals").arg("gfn1");
    }
    if config.method == ExcitedStateMethod::StdaXtbOriginal
        && !config.xtb4stda_path.trim().is_empty()
    {
        command
            .arg("--xtb4stda-executable")
            .arg(config.xtb4stda_path.trim());
    }

    // The engine's mode: Tamm-Dancoff, or the coupled response.
    command.arg(if config.tda && config.method.supports_full_response() {
        "--stda"
    } else if config.method.supports_full_response() {
        "--stddft"
    } else {
        // Every non-DFT route is Tamm-Dancoff only.
        "--stda"
    });

    command
        .arg("--stddft-roots")
        .arg(config.roots.max(1).to_string())
        .arg("--stddft-emax")
        .arg(format!("{:.3}", config.emax_ev))
        .arg("--memory")
        .arg(reference.memory_mb.to_string())
        .arg("--nproc")
        .arg(reference.nproc.max(1).to_string());
    command
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qchem_interfaces::method::{is_blocked, Severity};

    fn atoms(symbols: &[&str]) -> Vec<String> {
        symbols.iter().map(|s| s.to_string()).collect()
    }

    fn config(method: ExcitedStateMethod) -> ExcitedStateConfig {
        ExcitedStateConfig { method, ..Default::default() }
    }

    fn args(command: Command) -> Vec<String> {
        command
            .get_args()
            .map(|a| a.to_string_lossy().to_string())
            .collect()
    }

    fn value_after<'a>(args: &'a [String], flag: &str) -> Option<&'a str> {
        let at = args.iter().position(|a| a == flag)?;
        args.get(at + 1).map(|s| s.as_str())
    }

    fn build(config: &ExcitedStateConfig, multiplicity: i32) -> Vec<String> {
        args(spectrum_command(
            Path::new("behemoth"),
            Path::new("."),
            "in.xyz",
            0,
            multiplicity,
            config,
            &MethodConfig::default(),
        ))
    }

    #[test]
    fn the_defaults_are_five_roots_a_ten_ev_window_and_tda_on() {
        let c = ExcitedStateConfig::default();
        assert_eq!(c.roots, 5);
        assert_eq!(c.emax_ev, 10.0);
        assert!(c.tda, "TDA is the default, as asked");
    }

    /// Only TD-DFT takes a functional and a basis set.
    #[test]
    fn a_functional_and_basis_are_shown_only_for_tddft() {
        for method in ExcitedStateMethod::ALL {
            let fields = visible_fields(&config(method));
            let tddft = method == ExcitedStateMethod::Tddft;
            assert_eq!(fields.functional, tddft, "{method:?}");
            assert_eq!(fields.basis, tddft, "{method:?}");
            assert_eq!(fields.tda, tddft, "{method:?} TDA switch");
        }
    }

    /// The path field belongs to the one route that needs an external binary.
    #[test]
    fn the_xtb4stda_path_is_shown_only_where_it_is_needed() {
        for method in ExcitedStateMethod::ALL {
            assert_eq!(
                visible_fields(&config(method)).xtb4stda_path,
                method == ExcitedStateMethod::StdaXtbOriginal,
                "{method:?}"
            );
        }
    }

    /// Only the hybrids are offered, because the engine refuses the rest.
    #[test]
    fn the_functional_list_is_hybrids_only() {
        let names: Vec<&str> = spectrum_functionals().iter().map(|f| f.label).collect();
        assert_eq!(
            names,
            vec!["B3LYP", "B97", "M06-2X", "MN15", "PBE0", "revB3LYP", "TPSSh"]
        );
        // The two the engine would reject.
        assert!(!names.contains(&"PBE"));
        assert!(!names.contains(&"R2SCAN"));
        // And no composites, which carry their own basis set.
        assert!(spectrum_functionals().iter().all(|f| f.composite.is_none()));
        assert!(is_hybrid(DEFAULT_SPECTRUM_FUNCTIONAL));
    }

    /// The requested method: the original sTDA-xTB, blocked until its
    /// external generator is pointed at.
    #[test]
    fn the_original_stda_xtb_is_blocked_without_its_generator() {
        let mut c = config(ExcitedStateMethod::StdaXtbOriginal);
        let issues = validate(&c, &MethodConfig::default(), 1, &atoms(&["C", "O"]));
        assert!(is_blocked(&issues));
        assert!(
            issues.iter().any(|i| i.message.contains("xtb4stda")),
            "{issues:?}"
        );

        c.xtb4stda_path = "/opt/xtb4stda/bin/xtb4stda".to_string();
        let issues = validate(&c, &MethodConfig::default(), 1, &atoms(&["C", "O"]));
        assert!(!is_blocked(&issues), "{issues:?}");
    }

    /// The same element rule the other panels use, for the same reason.
    #[test]
    fn tasi_is_blocked_on_an_element_it_has_no_parameters_for() {
        let issues = validate(
            &config(ExcitedStateMethod::StdaTasi),
            &MethodConfig::default(),
            1,
            &atoms(&["C", "H", "Br"]),
        );
        assert!(is_blocked(&issues));
        assert!(issues.iter().any(|i| i.message.contains("Br")), "{issues:?}");
    }

    #[test]
    fn tasi_runs_on_its_own_elements() {
        let issues = validate(
            &config(ExcitedStateMethod::StdaTasi),
            &MethodConfig::default(),
            1,
            &atoms(&["C", "H", "O", "N", "S", "F", "Cl"]),
        );
        assert!(!is_blocked(&issues), "{issues:?}");
    }

    #[test]
    fn tddft_is_blocked_on_an_element_needing_an_ecp() {
        let issues = validate(
            &config(ExcitedStateMethod::Tddft),
            &MethodConfig::default(),
            1,
            &atoms(&["Rb", "C"]),
        );
        assert!(is_blocked(&issues));
        assert!(issues.iter().any(|i| i.message.contains("Rb")), "{issues:?}");
    }

    /// A pure functional cannot reach the engine, and the rule says so.
    #[test]
    fn tddft_is_blocked_on_a_functional_with_no_exact_exchange() {
        let reference = MethodConfig { functional: "pbe".to_string(), ..Default::default() };
        let issues =
            validate(&config(ExcitedStateMethod::Tddft), &reference, 1, &atoms(&["C"]));
        assert!(is_blocked(&issues));
        assert!(
            issues.iter().any(|i| i.message.contains("global-hybrid")),
            "{issues:?}"
        );
    }

    #[test]
    fn the_gfn1_route_warns_about_the_instability_that_was_observed() {
        let issues = validate(
            &config(ExcitedStateMethod::StdaXtbGfn1),
            &MethodConfig::default(),
            1,
            &atoms(&["C", "O"]),
        );
        // Worth trying, so a note rather than a block.
        assert!(!is_blocked(&issues));
        assert!(
            issues
                .iter()
                .any(|i| i.severity == Severity::Note && i.message.contains("unstable")),
            "{issues:?}"
        );
    }

    #[test]
    fn zero_roots_is_blocked() {
        let c = ExcitedStateConfig { roots: 0, ..config(ExcitedStateMethod::StdaTasi) };
        assert!(is_blocked(&validate(&c, &MethodConfig::default(), 1, &atoms(&["C"]))));
    }

    /// The window bounds the root count, and that is said before the run.
    #[test]
    fn the_energy_window_is_explained_up_front() {
        let issues = validate(
            &config(ExcitedStateMethod::StdaTasi),
            &MethodConfig::default(),
            1,
            &atoms(&["C"]),
        );
        let note = issues
            .iter()
            .find(|i| i.message.contains("not reported"))
            .unwrap_or_else(|| panic!("{issues:?}"));
        assert!(note.message.contains("10 eV"), "{}", note.message);
    }

    // ---- command construction ----

    /// TDA on gives `--stda`, TDA off gives `--stddft`: the requested switch.
    #[test]
    fn the_tda_switch_chooses_the_engine_mode() {
        let mut c = config(ExcitedStateMethod::Tddft);
        let tda = build(&c, 1);
        assert!(tda.contains(&"--stda".to_string()), "{tda:?}");
        assert!(!tda.contains(&"--stddft".to_string()));

        c.tda = false;
        let coupled = build(&c, 1);
        assert!(coupled.contains(&"--stddft".to_string()), "{coupled:?}");
        assert!(!coupled.contains(&"--stda".to_string()));
    }

    /// Turning TDA off on a route that cannot do the coupled response must
    /// still produce a command that runs, not one Behemoth rejects.
    #[test]
    fn a_tda_only_route_stays_tamm_dancoff_even_with_the_switch_off() {
        for method in [
            ExcitedStateMethod::StdaTasi,
            ExcitedStateMethod::StdaXtbGfn1,
            ExcitedStateMethod::StdaXtbOriginal,
        ] {
            let c = ExcitedStateConfig { tda: false, ..config(method) };
            let built = build(&c, 1);
            assert!(built.contains(&"--stda".to_string()), "{method:?}: {built:?}");
            assert!(!built.contains(&"--stddft".to_string()), "{method:?}");
        }
    }

    #[test]
    fn each_route_asks_for_the_method_it_needs() {
        assert_eq!(
            value_after(&build(&config(ExcitedStateMethod::StdaTasi), 1), "--method"),
            Some("tasi")
        );
        assert_eq!(
            value_after(&build(&config(ExcitedStateMethod::StdaXtbGfn1), 1), "--method"),
            Some("xtb")
        );
        assert_eq!(
            value_after(
                &build(&config(ExcitedStateMethod::StdaXtbOriginal), 1),
                "--method"
            ),
            Some("xtb")
        );
        assert_eq!(
            value_after(&build(&config(ExcitedStateMethod::Tddft), 1), "--method"),
            Some("rks")
        );
    }

    /// An open shell is UKS, as asked -- and `rks` would be refused anyway.
    #[test]
    fn an_open_shell_tddft_reference_is_uks() {
        assert_eq!(
            value_after(&build(&config(ExcitedStateMethod::Tddft), 3), "--method"),
            Some("uks")
        );
    }

    /// The requested behaviour: a TD-DFT run also writes the ground-state
    /// orbitals, in the same run, so the Surface tool can show the orbitals a
    /// root involves. Verified against the binary: `--wf2molden` alongside
    /// `--stda` produces both the spectrum and canonicalMO.molden.
    #[test]
    fn a_tddft_run_also_writes_the_ground_state_orbitals() {
        let built = build(&config(ExcitedStateMethod::Tddft), 1);
        assert!(built.contains(&"--wf2molden".to_string()), "{built:?}");
    }

    /// The sTDA routes must not ask for orbitals, and the reason is chemical
    /// rather than mechanical: a semi-empirical set is fitted to reproduce
    /// transitions, not the ground-state electronic structure, so it is not a
    /// description of the ground state and must not be offered as one. The
    /// restriction of `--wf2molden` to the Gaussian-basis references happens
    /// to agree, but this would hold without it.
    #[test]
    fn an_stda_route_does_not_ask_for_ground_state_orbitals() {
        for method in [
            ExcitedStateMethod::StdaTasi,
            ExcitedStateMethod::StdaXtbGfn1,
            ExcitedStateMethod::StdaXtbOriginal,
        ] {
            let built = build(&config(method), 1);
            assert!(
                !built.contains(&"--wf2molden".to_string()),
                "{method:?}: {built:?}"
            );
        }
    }

    /// The canonical set is the one the assignments are written in.
    #[test]
    fn the_orbital_file_read_back_is_the_canonical_one() {
        assert_eq!(CANONICAL_MOLDEN_FILE, "canonicalMO.molden");
        assert!(!CANONICAL_MOLDEN_FILE.contains("natural"));
    }

    /// Only TD-DFT carries a functional and a basis set; the semi-empirical
    /// routes would be given two answers to the same question.
    #[test]
    fn only_tddft_passes_a_functional_and_a_basis() {
        let tddft = build(&config(ExcitedStateMethod::Tddft), 1);
        assert_eq!(value_after(&tddft, "--functional"), Some("b3lyp"));
        assert_eq!(value_after(&tddft, "--basis"), Some("def2-svp"));

        for method in [ExcitedStateMethod::StdaTasi, ExcitedStateMethod::StdaXtbGfn1] {
            let built = build(&config(method), 1);
            assert!(!built.iter().any(|a| a == "--functional"), "{method:?}");
            assert!(!built.iter().any(|a| a == "--basis"), "{method:?}");
        }
    }

    /// GFN1 orbitals must be asked for explicitly; the other two routes are
    /// what `--method` already implies.
    #[test]
    fn only_the_gfn1_route_names_its_orbital_source() {
        let gfn1 = build(&config(ExcitedStateMethod::StdaXtbGfn1), 1);
        assert_eq!(value_after(&gfn1, "--stda-orbitals"), Some("gfn1"));

        for method in [
            ExcitedStateMethod::StdaTasi,
            ExcitedStateMethod::StdaXtbOriginal,
            ExcitedStateMethod::Tddft,
        ] {
            let built = build(&config(method), 1);
            assert!(
                !built.iter().any(|a| a == "--stda-orbitals"),
                "{method:?}: {built:?}"
            );
        }
    }

    #[test]
    fn the_original_route_passes_the_generator_path_it_was_given() {
        let c = ExcitedStateConfig {
            xtb4stda_path: "  /opt/xtb4stda  ".to_string(),
            ..config(ExcitedStateMethod::StdaXtbOriginal)
        };
        let built = build(&c, 1);
        // Trimmed: a path pasted with surrounding whitespace still works.
        assert_eq!(value_after(&built, "--xtb4stda-executable"), Some("/opt/xtb4stda"));
    }

    #[test]
    fn the_roots_and_window_reach_the_command_line() {
        let c = ExcitedStateConfig {
            roots: 20,
            emax_ev: 12.5,
            ..config(ExcitedStateMethod::StdaTasi)
        };
        let built = build(&c, 1);
        assert_eq!(value_after(&built, "--stddft-roots"), Some("20"));
        assert_eq!(value_after(&built, "--stddft-emax"), Some("12.500"));
        // The reference's resources come along too.
        assert_eq!(value_after(&built, "--memory"), Some("1024"));
        assert_eq!(value_after(&built, "--nproc"), Some("1"));
    }

    /// A root count of zero would be a command Behemoth rejects; the run is
    /// blocked before that, and the clamp is the second line of defence.
    #[test]
    fn a_zero_root_count_is_clamped_rather_than_passed_on() {
        let c = ExcitedStateConfig { roots: 0, ..config(ExcitedStateMethod::StdaTasi) };
        assert_eq!(value_after(&build(&c, 1), "--stddft-roots"), Some("1"));
    }

    /// Every route offered must be describable, or the dropdown has an entry
    /// the user cannot make sense of.
    #[test]
    fn every_route_is_labelled_and_described() {
        for method in ExcitedStateMethod::ALL {
            assert!(!method.label().is_empty(), "{method:?}");
            assert!(method.description().len() > 20, "{method:?}");
        }
        assert_eq!(ExcitedStateMethod::ALL.len(), 4);
    }
}
