# Running Excited-State Calculations from the Absorption Spectrum Tool — Design

Date: 2026-09-07
Status: approved for implementation

## Goal

The UV-Vis tool reads a finished TD-DFT output and plots it. This adds the
other half: running the calculation from the panel, the way the optimizer and
the frequency analysis already do.

Behemoth only. xTB has no excited-state module, and ORCA is not driven by
Beavyr at all -- its output is read, not produced.

The controls asked for: a method (sTDA-xTB or TD-DFT), a functional and basis
set when it is TD-DFT, a root count in every case, and a TDA switch for TD-DFT
that defaults to on.

## What Behemoth actually offers

Verified by running the binary and by reading `src/utils/cli/` and
`src/hamiltonian/stddft/`, because the flag names overlap in a way the help
alone does not settle.

There are **two unrelated excited-state engines**:

| engine | flags | notes |
|---|---|---|
| sTDA / sTD-DFT | `--stda`, `--stddft` | `--stddft-roots`, `--stddft-emax`, `--stda-orbitals` |
| RPA / TDA-xTB | `--rpa`, `--tda` | its own `--rpa-*` options; root count set by an energy window, not a root count |

This panel drives the **sTDA/sTD-DFT engine only**. It is the one that takes a
root count, it is the one that covers both of the requested methods, and its
four orbital sources share one set of controls. The RPA/TDA-xTB engine has a
different parameterisation and would need its own controls; it is out of scope
below.

### The four ways into the sTDA engine

`--stda` is the Tamm-Dancoff mode of the engine and `--stddft` the coupled
one. Which orbitals it builds on is set by `--stda-orbitals`, which
Behemoth's own CLI tests show is inferred from `--method` when not given:

| requested method | command | works today? |
|---|---|---|
| sTDA-xTB (published) | `--method xtb --stda` | **needs an external `xtb4stda` binary** |
| sTDA-xTB (GFN1 orbitals) | `--method xtb --stda --stda-orbitals gfn1` | runs, but was numerically unstable in testing |
| sTDA (TASI) | `--method tasi --stda` | yes |
| TD-DFT | `--method rks\|uks --functional H --basis B --stda` or `--stddft` | yes, `H` a global hybrid |

Three findings that shape the UI:

* **The published sTDA-xTB needs a binary the machine does not have.**
  `--method xtb --stda` fails with "failed to execute original sTDA-xTB
  generator `xtb4stda`; install xtb4stda, set BEHEMOTH_XTB4STDA, or pass
  `--xtb4stda-executable`". `~/Programs/xtb4stda` is a source tree with no
  built binary. So the method is offered with a path field of its own, in the
  same shape as the xTB and Behemoth paths, and is blocked with that
  explanation until one is given.

* **GFN1 orbitals were unstable on both test molecules.** Formaldehyde gives
  "sTDA A' matrix is unstable: lowest eigenvalue = -0.0763 Hartree", benzene
  -0.0164. It is still offered, because it needs no external binary and may
  well be fine on other systems, but the panel says what was observed rather
  than letting the user discover it as an opaque failure.

* **`--stddft` requires a global-hybrid RKS/UKS reference.** With a pure
  functional: "sTDA/sTD-DFT with normal DFT requires a global-hybrid
  functional with nonzero exact exchange". With TASI: "`--stddft` currently
  requires a global-hybrid `--method rks` or `--method uks` reference". So the
  TDA switch is meaningful for TD-DFT alone; every other route is
  Tamm-Dancoff by construction, and the switch is hidden for them rather than
  shown and ignored.

### Hybrids only

`behemoth --list func` groups functionals by rung and prints each hybrid's
exact-exchange fraction. Of the nine the panels offer, seven are global
hybrids and two are not:

| | |
|---|---|
| offered for TD-DFT | B3LYP 20%, B97 19.4%, M06-2X 54%, MN15 44%, PBE0 25%, revB3LYP 20%, TPSSh 10% |
| omitted | **PBE** (GGA), **R2SCAN** (tau meta-GGA) -- no exact exchange, so the engine refuses them |

The composites are omitted too: two of the five are not hybrids either
(b97-3c, r2scan-3c), HF-3c is a method rather than a functional, and all of
them carry their own basis set, which conflicts with a panel that asks for
one. The two hybrid composites are no loss -- `pbeh-3c` and `b3lyp-3c` are
their functionals at a fixed small basis, which this panel can express
directly.

### Root count and energy window

`--stddft-roots N` (default 5) is the root count. It is not the whole story:
`--stddft-emax EV` (default 10 eV) is an energy window, and roots outside it
are not reported at all. Asking for 6 roots on formaldehyde returned 3,
because the rest lay above 10 eV.

A root count that silently does not happen is worse than no root count, so the
window is exposed beside it. It is the only field here that was not asked for,
and it is included precisely so that the field that *was* asked for means what
it says.

### Output to parse

One table gives the spectrum:

```
root       energy(eV)     lambda(nm)            f
1            4.115545         301.26     0.000000
2            8.338299         148.69     0.173763
```

and a second section gives the assignments, with the engine's name in the
header -- `sTDA` for `--stda`, `sTD-DFT` for `--stddft`, so the parser accepts
both:

```
Root 1   sTDA    =     4.1155 eV      301.26 nm  f = 0.00000
       8 HOMO       [O lonepair] ==>    9 LUMO       [C-O antibond]  99.89%
```

Each contribution carries the MO index, its HOMO/LUMO-relative name, a
bracketed Pipek-Mezey character label, and a percentage. `Excitation` already
holds the indices and the weight; the orbital names and the character labels
are extra, and are kept because they are the most useful part of the line for
reading a spectrum.

The header block also names the orbital source, which becomes the method
string shown in the panel:

```
Orbital source: Behemoth RKS PBE0 (this response model does not replace the orbital Hamiltonian)
```

### The ground-state orbitals, in the same run

A spectrum is far more use with the orbitals its assignments name. Behemoth
can write them without a second calculation: `--wf2molden` "writes
canonicalMO.molden and naturalMO.molden for Gaussian-basis HF, MP2/RI-MP2,
RKS, and UKS references", and it works alongside `--stda` -- verified on
formaldehyde, where one run produced the six-root spectrum and both Molden
files.

So the TD-DFT route always passes `--wf2molden`, the run reads
`canonicalMO.molden` back out of its scratch directory before that directory
is discarded, and the Surface tool adopts it automatically.

**Canonical, not natural.** The assignments are given as canonical MO indices
("8 HOMO ==> 9 LUMO"), so those are the orbitals the numbers refer to.
`naturalMO.molden` is the right file for a correlated density, which is not
what this shows.

**Not for the sTDA routes.** `--wf2molden` is documented for the
Gaussian-basis references only, and the xTB and TASI routes have no Gaussian
basis to write, so they do not ask for it and nothing in the Surface tool is
disturbed when one of them runs.

The file parses with the existing `orbitals::molden` reader unchanged, and the
HOMO it reports -- MO 8 -- agrees with the assignment section's own `8 HOMO`,
which is a useful cross-check on both.

`OrbitalState` gains `provenance` and `molden_text`. A computed set has no file
behind it, so the panel names the run it came from and offers a
"Save .molden..." button; the text is held in memory because the scratch
directory is gone, which is exactly why the button is needed. A loaded file
keeps showing its name and neither field. Both routes go through one
`OrbitalState::adopt`, so a set that was computed and a set that was loaded
cannot end up in different states.

## Code

### `src/uvvis/behemoth.rs` (new)

`parse_behemoth_stda(stdout, source) -> Result<TddftResult, String>`, filling
the same `TddftResult` that `orca.rs` fills. `types.rs` was written to be
program-independent for exactly this, so the UI needs no change to display
the result.

`Excitation` gains two optional fields, `from_label` and `to_label`, for the
`HOMO` / `LUMO+1` names and the bracketed characters. ORCA prints neither, so
they stay `None` there and the existing display is unaffected.

Behemoth reports no `<S²>` per root and no transition dipole, so `s2`,
`d2_au2` and `transition_dipole_au` stay `None` -- which the panel already
handles, since a restricted ORCA run has no `<S²>` either.

### `src/uvvis/excited_state.rs` (new)

```rust
pub enum ExcitedStateMethod { StdaXtbOriginal, StdaXtbGfn1, StdaTasi, Tddft }

pub struct ExcitedStateConfig {
    pub method: ExcitedStateMethod,
    pub roots: u32,          // 5
    pub emax_ev: f64,        // 10.0
    pub tda: bool,           // true
    pub xtb4stda_path: String,
}
```

with `visible_fields`, `validate` and the command construction, mirroring
`qchem_interfaces::method` in shape so there is one pattern for this in the
codebase rather than two. The functional, basis, memory and thread count come
from the existing `MethodConfig`, which already holds exactly those and needs
no duplicate.

Validation, at `Block` severity:

* the published sTDA-xTB with no `xtb4stda` path given
* TASI with any element outside its seven, reusing `method::TASI_ELEMENTS`
* TD-DFT with a non-hybrid functional -- unreachable through the dropdown,
  which offers only hybrids, but stated where the other rules live
* any element from rubidium up under TD-DFT, reusing the ECP rule

and at `Note` severity: what the GFN1 instability was, that the root count is
bounded by the energy window, and that the non-TDA switch needs a hybrid.

### `src/uvvis/run.rs` (new)

The background task: `AsyncComputeTaskPool`, `block_on(poll_once)`, a child
process that cancellation kills, a scratch run directory. Structurally the
same as `xtb_freq.rs` rather than a second pattern -- the differences are the
command line and that the result is parsed from stdout.

## UI

The Absorption Spectrum panel gains a collapsible "Run Calculation" block
above the existing load button, closed by default, in the same shape as the
frequency panel's "Run Freq Calc":

```
▼ Run Calculation
    Program  [Behemoth]        (only program with an excited-state module)
    Path     [...] [Browse…]
    Method   [TD-DFT ▾]
    Functional [B3LYP ▾]  Basis set [def2-SVP ▾]
    Roots [5]   Energy window (eV) [10]   [x] Tamm-Dancoff (TDA)
    Charge [0]  Multiplicity [1]
    Memory (MB) [1024]  nproc [1]
    [Run Spectrum]
```

The program row is a fixed label rather than a dropdown: there is one program
that can do this, and a dropdown of one is a lie about the choice available.

A finished run fills the same `UvVisState` a loaded file does, so the spectrum
plot, the state table and the `.dat` export all work on it unchanged.

## Testing

* the root table and the assignment section parsed from captured output of
  real runs, one `--stda` and one `--stddft`, as fixtures
* both engine labels (`sTDA`, `sTD-DFT`) accepted in the assignment header
* the root count in the table agreeing with the assignment section
* a run whose roots were cut short by the energy window parsing to the rows
  that are there, with a note rather than an error
* field visibility per method: no functional or basis except for TD-DFT, no
  TDA switch except for TD-DFT
* each Block rule firing and staying quiet appropriately, including the
  xtb4stda path and the TASI element rule
* command construction for all four methods, including that TD-DFT with TDA
  off uses `--stddft` and with TDA on uses `--stda`, and that only TD-DFT
  passes a functional and basis
* the parsed result feeding `peaks_nm`/`peaks_ev` so the existing plot has
  something to draw

## Out of scope

The RPA/TDA-xTB engine (`--rpa`, `--tda`), which is a different response model
with its own parameterisation and an energy window in place of a root count.
It is worth having, and it is the natural fallback for an xTB spectrum on a
molecule TASI does not cover, but it needs its own controls rather than being
bent into these.

Also out of scope: `--stddft-damping` and the custom `--stddft-yj`/`-yk`
exponents, `--stddft` for unrestricted references beyond what falls out of the
multiplicity, and building `xtb4stda` on the user's behalf.
