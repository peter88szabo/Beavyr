# Method Selection for the Optimizer and Frequency Panels — Design

Date: 2026-09-07
Status: approved for implementation

## Goal

Both job panels can now choose which program runs a calculation, but not *what
it should calculate*. xTB always runs its default parametrisation and Behemoth
always runs GFN1-xTB, because that is all the command construction ever asks
for.

This adds a method block that rolls down under the program selector in both
panels: pick a method, and the fields that method needs appear underneath it —
functional, dispersion correction, basis set, two-electron treatment, memory
and thread count. Fields that do not apply are not shown, and a configuration
that cannot run is blocked before the program is launched rather than after it
fails.

The same block appears in the optimizer and in the frequency panel, from one
implementation, so the two cannot drift apart.

## What the programs actually accept

Verified against the binaries, not assumed. Several of the requested options
turned out not to exist under the names expected, and one cannot be used in
these panels at all.

### xTB

From `xtb --help` on the build at
`~/Programs/orca_6_1_1_linux_x86-64_shared_openmpi418/xtb`:

| | |
|---|---|
| `--gfn 1` / `--gfn 2` | GFN1-xTB and GFN2-xTB; `--gfn` defaults to 2 |
| `--gfnff` (alias `--gff`) | GFN-FF, the force field. The name is **GFN-FF** |
| `-P N`, `--parallel N` | thread count |

Two findings:

* **GFN0-xTB is not available.** `--gfn 0` is parsed but dies with
  `Parameter file param_gfn0-xtb.txt not found!`. It is omitted from the list.
  Should that file be installed later, adding it back is one list entry.
* **xTB has no memory option.** The memory field is therefore Behemoth-only,
  and is hidden when xTB is selected rather than shown and silently ignored.

Confirmed working, on water: GFN1 −5.768649 Eh, GFN2 −5.070376 Eh, GFN-FF
−0.327216 Eh. The GFN1 figure agrees with Behemoth's own GFN1 energy for the
same geometry to seven decimals, which is a useful cross-check on both.

### Behemoth

From `behemoth --help` and `behemoth --list func|basis|disp|abinitio`:

| | |
|---|---|
| `--method xtb\|tasi\|hf\|uhf\|rohf\|mp2\|ri-mp2\|rks\|uks` | `xtb` is GFN1; `hf` picks restricted or unrestricted from `--multi` |
| `--functional NAME` | required with `rks`/`uks`, ignored otherwise |
| `--dispersion d3bj\|d3zero\|d3bjm\|d3zerom\|d3op\|d3cso\|d3z\|d4` | added after SCF; damping parameters looked up per functional |
| `--basis NAME` | required for `hf`, `uhf`, `rohf`, `mp2`, `ri-mp2` and the DFT methods |
| `--two-electron exact\|rijk` | `exact` is the default |
| `--jk-basis NAME` | RIJK fitting basis, default `def2-universal-jkfit` |
| `--memory MB` | default 1024 |
| `--nproc N` | default is the environment's own thread configuration |

Four findings that shape the design:

* **There is no `dft` method.** DFT is `rks` for a closed shell and `uks`
  otherwise. The UI keeps the single entry "DFT" that was asked for and picks
  the variant from the multiplicity already set in the panel: `--multi 1` gives
  `rks`, anything higher gives `uks`. This matches Behemoth's own rule, which
  refuses `rks` for an open shell.
* **An open shell is always unrestricted.** `--method rks` on a doublet is
  refused outright — "requires `--multi 1`; use `--method uks`" — so the
  multiplicity has to pick `rks` or `uks` whatever else happens. For Hartree-Fock
  it turns out not to matter: `--method hf` on the water cation doublet runs a
  UHF SCF and gives exactly the energy `--method uhf` gives, −75.5621508290 Eh,
  so `hf`'s automatic choice is already unrestricted and `--spinpol`'s
  documented ROHF default does not reach it. `uhf` is nonetheless named
  explicitly above a singlet, so the command line states which reference was
  used rather than depending on that internal choice staying as it is. HF-3c
  needs nothing: it picks UHF for an open shell by itself, verified the same
  way.
* **HF-3c is a method, not a functional.** The help is explicit: "HF-3c is
  plain Hartree-Fock, not a DFT functional, so it is requested as
  `--method hf-3c` rather than `--functional hf-3c`". It is still listed with
  the other composites in the functional dropdown, as requested, and the
  translation happens where the command is built.
* **RIJK cannot be used here at all.** It is energy-only: "no analytic
  gradient yet, so `--gradient`/`--opt`/`--hessian`/`--freq`/`--scan`/`--md`
  are refused with it". Every job either panel runs is one of those. The
  control is therefore shown for DFT and HF but **disabled**, with the reason
  next to it, so the option is discoverable and cannot be set to something that
  is guaranteed to fail. It becomes live for free if a single-point panel is
  ever added.
* **No ECP support.** Nothing in the help, the basis catalog or the ab-initio
  list mentions effective core potentials or pseudopotentials. The `lanl2dz`
  family of files is present in the catalog, but those basis sets are
  valence-only and meaningless without the ECP integrals to go with them, so
  they are not offered.

The requested functional names were checked one by one against
`--list func`. All nine exist, with one spelling correction: **M06-2X is
`m06-2x`**, not `m062x`. The rest are `b3lyp`, `pbe`, `pbe0`, `b97`, `tpssh`,
`revb3lyp`, `r2scan`, `mn15`.

TASI's parameter file, `parameters/tasi_chonclsf_parameters.json`, carries
atom parameters for exactly seven elements: **C, Cl, F, H, N, O, S**.

## Data model

`src/qchem_interfaces/method.rs`, with no dependency on egui so that every rule
below is unit-testable:

```rust
pub struct MethodConfig {
    pub xtb: XtbMethod,               // Gfn1 | Gfn2 (default) | GfnFf
    pub behemoth: BehemothMethod,     // Tasi | Gfn1Xtb (default) | Hf | Dft | Mp2 | RiMp2
    pub functional: String,           // "b3lyp"
    pub dispersion: Dispersion,       // None (default) | D3Bj | .. | D4
    pub basis: String,                // "def2-svp"
    pub two_electron: TwoElectron,    // Exact (default) | Rijk
    pub memory_mb: u32,               // 1024
    pub nproc: u32,                   // 1
}
```

Behemoth's default method stays GFN1-xTB: it is Behemoth's own default, and it
is what the panels already send today, so an existing habit does not silently
turn into an expensive DFT job.

`nproc` defaults to 1. Behemoth with no `--nproc` uses every thread the
environment offers, which competes with Beavyr's own rendering; 1 is explicit
and reproducible, and choosing more is the user's call, as they noted.

Each panel owns a `MethodConfig`. They are not shared, for the same reason the
program selection is not: running a Hessian at a different level than the
optimization is a legitimate thing to want.

## Field visibility

One function, read by both the UI and its tests, so that "is the basis field
showing?" has exactly one answer:

```rust
pub struct VisibleFields {
    pub functional: bool,
    pub basis: bool,
    pub dispersion: bool,
    pub two_electron: bool,   // always disabled in these panels; see above
    pub memory: bool,
    pub nproc: bool,
}
pub fn visible_fields(program: QcProgram, config: &MethodConfig) -> VisibleFields
```

| program and method | fields shown |
|---|---|
| xTB — GFN1, GFN2, GFN-FF | `nproc` |
| Behemoth — TASI, GFN1-xTB | `memory`, `nproc` |
| Behemoth — HF, MP2, RI-MP2 | `basis`, `memory`, `nproc` |
| Behemoth — DFT, plain functional | `functional`, `dispersion`, `basis`, `two_electron`, `memory`, `nproc` |
| Behemoth — DFT, composite functional | `functional`, `memory`, `nproc` |

A composite suppresses both the basis set and the dispersion correction,
because it defines its own. Instead of those two fields the panel states what
the composite brings: r2SCAN-3c, for instance, is "r2SCAN + def2-mTZVPP +
D4(BJ)+ATM + damped gCP", which is worth showing rather than leaving the user
to wonder where the basis went.

## Catalogs

### Functionals

The nine common ones in alphabetical order, then a separator, then the five
composites:

```
B3LYP  B97  M06-2X  MN15  PBE  PBE0  R2SCAN  revB3LYP  TPSSh
---
B3LYP-3c  B97-3c  HF-3c  PBEh-3c  r2SCAN-3c
```

Behemoth's full catalog is far larger, and `--functional` also accepts raw
LibXC specifications. Offering all of it in a dropdown would be unusable; this
is the short list that was asked for, and the field remains free text so
anything in `--list func` can still be typed.

### Basis sets

Grouped as requested, from the 45 in `--list basis`:

| group | entries |
|---|---|
| def2 | def2-SV(P), **def2-SVP**, def2-SVPD, def2-TZVP, def2-TZVPD, def2-TZVPP, def2-TZVPPD, def2-QZVP, def2-QZVPD, def2-QZVPP, def2-QZVPPD |
| pc | pc-0 … pc-4, aug-pc-0 … aug-pc-4 |
| cc | cc-pVDZ, cc-pVTZ, cc-pVQZ, aug-cc-pVDZ, aug-cc-pVTZ, aug-cc-pVQZ |
| STO-nG | STO-3G, STO-4G, STO-5G, STO-6G |
| Pople | 3-21G, 6-31G, 6-31G\*, 6-31G\*\* |

def2-SVP is the default. Display names differ from the names passed on the
command line — `6-31G*` is `6-31g_st` and `6-31G**` is `6-31g_st_st` — so the
catalog stores both and the dropdown shows the readable one.

Deliberately not offered:

* `lanl2dz`, `lanl2dzdp`, `lanl2tz`, `lanl2tz+`, `lanl2tzf`, `mod-lanl2dz` —
  ECP basis sets with no ECP implementation behind them.
* `minix`, `def2-msvp-pbeh3c`, `def2-mtzvp-b973c`, `def2-mtzvpp-r2scan3c` —
  each belongs to a composite and is selected by it automatically.
* every `-rifit`, `-jfit` and `-jkfit` set — these are auxiliary bases chosen
  automatically, which is the next section.

### Dispersion

`None` (the default, meaning no `--dispersion` flag at all), then the seven D3
dampings `d3bj`, `d3zero`, `d3bjm`, `d3zerom`, `d3op`, `d3cso`, `d3z`, then
`d4`. Tkatchenko–Scheffler is listed by Behemoth as "not available yet" and is
omitted.

## Validation

`validate(program, config, multiplicity, atoms) -> Vec<MethodIssue>`, where an
issue carries a severity and a message. There are two severities:

* `Block` — shown as a warning in the panel *and* disables the Run button, so a
  configuration that cannot work is stopped before a process is spawned rather
  than after the program fails.
* `Note` — shown as ordinary text, explaining something the panel has decided
  on the user's behalf. Nothing to fix.

| severity | condition | message |
|---|---|---|
| Block | any element with Z ≥ 37, and a Gaussian-basis method (HF, DFT, MP2, RI-MP2) | names the offending elements; explains that these need an effective core potential beyond krypton and that no ECP is implemented, so the calculation must not be run |
| Block | TASI with any element outside C, H, N, O, F, S, Cl | names the offending elements and lists the seven TASI is parameterised for |
| Block | MP2 or RI-MP2 with multiplicity > 1 | Behemoth requires a closed-shell reference for both |
| Block | RIJK with any job in these panels | unreachable through the UI, since the control is disabled, but the rule is stated where every other rule lives rather than relying on the UI to be the only guard |
| Note | RI-MP2 selected | the correlation-fitting basis is chosen automatically; there is nothing to set |
| Note | composite functional selected | what it brings: its own basis set and dispersion correction |
| Note | RIJK control | no analytic gradient yet, so it cannot be used for optimization or frequencies |

Z ≥ 37 is the threshold because that is where the def2 family switches to
ECPs — rubidium onwards. Elements up to krypton, transition metals included,
are all-electron and fine.

The ECP check is keyed on the *method*, not the basis set: every Gaussian-basis
method needs core functions for a heavy element regardless of which of the
offered sets is chosen. xTB, GFN-FF and TASI are parametrised rather than
basis-set methods and are exempt — subject, for TASI, to its own element rule.

## Command construction

Existing signatures gain a `&MethodConfig`. The mapping, in one place:

**xTB** — `--gfn 1` or `--gfn 2` or `--gfnff`, plus `-P <nproc>`. No memory
flag exists, so none is passed.

**Behemoth** — `--nproc <nproc>` and `--memory <memory_mb>` always, plus:

| method | flags |
|---|---|
| TASI | `--method tasi` |
| GFN1-xTB | `--method xtb` |
| HF, multiplicity 1 | `--method hf --basis <basis>` |
| HF, multiplicity > 1 | `--method uhf --basis <basis>` |
| MP2 / RI-MP2 | `--method mp2` / `--method ri-mp2`, `--basis <basis>` |
| DFT, plain functional | `--method rks` (multiplicity 1) or `uks` (> 1), `--functional <f>`, `--basis <b>`, and `--dispersion <d>` unless None |
| DFT, composite other than HF-3c | `--method rks` or `uks` by the same rule, `--functional <f>` — no basis, no dispersion |
| DFT, HF-3c | `--method hf-3c` — no functional, no basis, no dispersion |

The spin treatment follows the panel's multiplicity and is always the
unrestricted one above a singlet: `uhf` rather than `hf`, and `uks` rather than
`rks`. For DFT that is forced — `rks` refuses an open shell. For Hartree-Fock
it is a deliberate redundancy, since `hf` already chooses UHF on its own, so
that the command line records the reference rather than implying it. HF-3c
takes no method name alongside it and needs none: it selects UHF for an open
shell itself.

## UI

`method_ui.rs` exposes one function that both panels call directly under their
program row:

```rust
pub fn method_config_ui(ui: &mut egui::Ui, program: QcProgram,
                        config: &mut MethodConfig, multiplicity: u32,
                        atoms: &[String]) -> bool   // any Block issue?
```

It draws a collapsible "Method" section, open by default so the level of theory
is visible without hunting for it, and returns whether the configuration is
blocked so the caller can disable its Run button. The header line shows the
current method when collapsed.

Both panels already own the charge and multiplicity fields the validation
needs, so nothing new has to be threaded through.

## Testing

The catalogs and rules are pure functions, so they are tested directly:

* every catalog entry's command-line name is one Behemoth accepts — checked
  against the strings captured from `--list basis` and `--list func`, so a
  typo in a display name cannot reach the command line
* field visibility for each program-and-method combination in the table above,
  including that a composite hides the basis and dispersion fields and a plain
  functional shows them
* each Block rule fires on a structure that violates it and stays quiet on one
  that does not: a molecule containing rubidium under DFT, a TASI job on a
  molecule containing phosphorus, MP2 as a triplet
* the ECP rule does not fire for xTB or GFN-FF, and does not fire for a
  transition metal, which is all-electron in def2
* command construction for each row of the mapping table, including that HF-3c
  produces `--method hf-3c` with no `--functional`, that a composite passes no
  `--basis`, and that DFT follows the multiplicity into `rks` or `uks`
* xTB gets no memory flag

## Out of scope

Behemoth's other tasks and options — `--gcp`, `--d3-atm`, the SCF controls,
`--opt-coord`, solvation, and everything under the TS/IRC/NEB and sTDA
headings. Each wants its own UI, and this change is about the level of theory
the two existing jobs run at.

Also out of scope: a single-point energy panel, which is the place the RIJK
control would become live, and free-text entry of a raw LibXC functional
specification, which the field permits but the dropdown does not advertise.
