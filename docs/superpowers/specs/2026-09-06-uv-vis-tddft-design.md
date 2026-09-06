# UV-Vis / TD-DFT Spectrum Tool — Design

Date: 2026-09-06
Status: approved for implementation

## Goal

A new tool, reached from its own icon on the rail, that reads a quantum-chemistry
output file, extracts the TD-DFT excited states (roots), lists them with their
energies, oscillator strengths, transition dipoles, spin contamination and the
orbital-to-orbital excitations that make them up, and plots the simulated
absorption spectrum with the same line-shape / width / stick / export controls
the IR spectrum already has.

ORCA is the only supported program for now. The parser lives behind a boundary
that a second program can be added to later without touching the UI.

## Scope

In scope:

* Load an **existing** ORCA output file (`.out`, `.log`, any text file). This tool
  does not launch ORCA — unlike the xTB optimizer and Hessian tools, which run
  the binary themselves.
* Parse the `TD-DFT/TDA EXCITED STATES` (and plain `TD-DFT EXCITED STATES`) block:
  per-root energy in au / eV / cm⁻¹, `<S**2>`, estimated multiplicity, and the
  weighted list of single-orbital excitations.
* Parse the `ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS` table:
  wavelength, oscillator strength, D², and the DX/DY/DZ transition-dipole
  components.
* Show a per-root table, sortable in file order, with an expandable list of the
  contributing excitations.
* Plot the absorption spectrum in **nm or eV** (user's choice) with Gaussian /
  Lorentzian / Voigt / Sticks and an adjustable width, hover readout, and `.dat`
  export.

Out of scope (deliberately): running ORCA, CD/ECD spectra, velocity-gauge
spectra, natural transition orbitals, SOC-corrected states, and rendering the
involved orbitals in 3D. The parser keeps enough structure that CD and velocity
tables can be added later as extra tables rather than a rewrite.

## Module layout

Two new module trees. Per the standing rule, no `mod.rs` contains an algorithm —
each holds declarations (and, where needed, the Bevy resource struct) only.

### `src/spectrum/` — extracted, program-agnostic spectrum machinery

The IR spectrum's line-shape maths and plot painter currently live inside
`qchem_interfaces/xtb_freq.rs` and are hard-wired to cm⁻¹. They move out
unchanged in behaviour and become unit-agnostic, so IR and UV-Vis share one copy.

* `broadening.rs`
  * `BroadeningKind { Gaussian, Lorentzian, Voigt, Sticks }` (moved verbatim)
  * `lineshape(kind, x, x0, height, width) -> f64`
  * `sample_spectrum(peaks, kind, width, x_min, x_max, n_points) -> Vec<(f64, f64)>`
  * `format_spectrum_dat(peaks, kind, width, x_min, x_max, resolution) -> String`
  * `format_spectrum_sticks(peaks) -> String`

  The `frequencies` / `intensities` parameter pair becomes a single
  `peaks: &[(f64, f64)]`, which is what makes these callable for any x-axis.

* `plot.rs`
  * `SpectrumAxis { x_min, x_max, tick_step: Option<f64>, tick_decimals, title }`
  * `nice_tick_step(span) -> f64` — a 1/2/5×10ⁿ step giving roughly 6–10 ticks.
    IR keeps its explicit 500 cm⁻¹ step by passing `Some(500.0)`; UV-Vis passes
    `None` and gets a sensible step for whichever unit and range it is showing.
  * `draw_spectrum(ui, peaks, kind, width, axis, hover_label) `— the white-background
    painter, grid, curve-or-sticks, peak markers and hover box, lifted from
    `ir_spectrum_window`. `hover_label: &dyn Fn(usize) -> String` lets each caller
    word its own readout ("1234.5 cm⁻¹ / 12.3 km/mol" vs "S5 · 412.3 nm · f=0.012").

**Broadening convention.** The line shape is applied in the *currently displayed*
unit: in nm mode a Gaussian is Gaussian in nm, in eV mode it is Gaussian in eV.
These are not the same curve, because λ = 1239.841 984 / E is nonlinear. This is
what simple spectrum plotters do and it is what makes the width control mean what
the label says, so it is the chosen behaviour; the panel states the unit next to
the width field.

### `src/uvvis/` — the tool itself

* `mod.rs` — module declarations plus the `UvVisState` Bevy resource.
* `types.rs` — `Spin`, `Excitation`, `ExcitedState`, `TddftResult`.
* `orca.rs` — `parse_orca_tddft(text: &str) -> Result<TddftResult, String>`.
* `ui.rs` — `uvvis_panel(ui, state)` and `uvvis_spectrum_window(ctx, state)`.

## Data model

```rust
pub enum Spin { Alpha, Beta, Unspecified }

pub struct Excitation {
    pub from_orbital: u32,
    pub from_spin: Spin,
    pub to_orbital: u32,
    pub to_spin: Spin,
    /// |c|^2, the weight ORCA prints; only those above its own 1e-2 cutoff appear.
    pub weight: f64,
    /// The signed CI coefficient, when the output prints one.
    pub coefficient: Option<f64>,
}

pub struct ExcitedState {
    pub root: usize,
    pub energy_au: f64,
    pub energy_ev: f64,
    pub energy_cm1: f64,
    /// From the absorption table when present, else derived from the energy.
    pub wavelength_nm: f64,
    /// Only for unrestricted references.
    pub s2: Option<f64>,
    pub multiplicity: Option<u32>,
    /// From the electric-dipole absorption table; `None` if that table is absent.
    pub oscillator_strength: Option<f64>,
    pub d2_au2: Option<f64>,
    pub transition_dipole_au: Option<[f64; 3]>,
    pub excitations: Vec<Excitation>,
}

pub struct TddftResult {
    pub source: PathBuf,
    pub program: String,      // e.g. "ORCA 6.1.0"
    pub method: String,       // e.g. "TD-DFT/TDA"
    pub unrestricted: bool,
    pub ground_state_s2: Option<f64>,
    pub states: Vec<ExcitedState>,
}
```

`unrestricted` is true when any state carries an `<S**2>` value or the file
contains a `UHF SPIN CONTAMINATION` block. It drives whether the table shows the
`<S²>` and `Mult` columns at all — showing an empty column for a closed-shell
run is noise.

## Parsing ORCA

Three independent passes over the text, each tolerant of the others being absent,
so a run that produced states but no absorption table still lists its roots.

**1. Excited-state block.** Locate a line matching `TD-DFT/TDA EXCITED STATES`,
`TD-DFT EXCITED STATES`, or either with a ` (SINGLETS)` / ` (TRIPLETS)` suffix.
Then read `STATE` records:

```
STATE  1:  E=   0.071104 au      1.935 eV    15605.6 cm**-1 <S**2> =   8.870929 Mult 6
    89b ->  90b  :     0.988112 (c= -0.99403830)
```

The `<S**2> = … Mult …` tail is optional (restricted runs omit it in older
versions). Excitation lines are `<int><spin?> -> <int><spin?> : <weight> [(c= <coeff>)]`,
where the spin letter `a`/`b` is present for unrestricted references and absent
for restricted ones. The block ends at the first line that is neither a `STATE`
record, an excitation, nor blank.

**2. Absorption table.** Locate `ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC
DIPOLE MOMENTS`, then the column-name row, and map columns **by name** rather
than by fixed position, because ORCA 5 and ORCA 6 order and name them
differently:

```
     Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ
                      (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)
  0-6A  ->  1-6A    1.934849   15605.6   640.8   0.002236399   0.04718   0.00093  -0.11230   0.18592
```

A data row is split into a leading transition label (everything up to and
including the `->` and the token after it) and a numeric tail. The name row gives
the tail's meaning: the token starting `fosc` is the oscillator strength,
`Wavelength` is nm, `DX`/`DY`/`DZ` are the dipole components, and of the two
`Energy` columns the first is eV and the second cm⁻¹ (their order in the units
row confirms it). If the name row cannot be understood, the whole table is
skipped with a warning rather than being guessed at.

Rows are matched to states by the root number parsed out of the *destination*
label (`1-6A` → root 1), not by row order, so a table that omits a root does not
shift every later assignment.

**3. Ground-state spin.** `Expectation value of <S**2>     :     8.756797`.

### Failure behaviour

`parse_orca_tddft` returns `Err` only when the file contains no recognisable
excited-state block at all — the message names what was looked for. Partial
information (states but no absorption table, or a table whose columns could not
be mapped) yields `Ok` with `oscillator_strength: None` on the affected states
and a note the panel shows, because a root list without intensities is still
useful.

## UI

### Rail

A twelfth tab, `Tab::UvVis`, icon `🌈`, title `UV-Vis (TD-DFT)`, inserted directly
after `Vibrations` so the spectroscopy tools sit together.

### Panel

* **Load ORCA output…** — `rfd` file dialog, filtered to `out`/`log`/`txt`/all.
  The loaded path and the parsed program/method line are shown underneath.
* **Summary** — number of roots, whether the reference was unrestricted, and the
  ground-state `<S²>` when there is one.
* **State table**, one row per root, in a scroll area:
  `#`, `E (eV)`, `λ (nm)`, `f`, and — only when unrestricted — `<S²>` and `Mult`.
  Clicking a row expands the contributing excitations beneath it, each shown as
  `89b → 90b   0.988 (c = −0.994)`, sorted by descending weight so the dominant
  configuration reads first. A row also shows its transition-dipole components on
  hover.
* **Show Absorption Spectrum** — opens the plot window.

Roots with zero oscillator strength stay in the list (they are real states) but
contribute nothing to the plotted curve, which is physically correct.

### Spectrum window

Same controls as the IR window, plus a unit selector:

* Line shape: Gaussian / Lorentzian / Voigt / Sticks
* Width: in nm (default 20 nm) or eV (default 0.3 eV), disabled for Sticks
* Unit: **nm** (default) or **eV**
* Export .dat…

Axis range: in nm, from `max(0, λ_min − 50)` to `λ_max + 50`; in eV, from
`max(0, E_min − 0.5)` to `E_max + 0.5`. Ticks are the auto "nice" step. Unlike
the IR spectrum there is no reason to force the axis down to zero — a UV-Vis
plot that starts at 0 nm wastes almost all of its width.

Export resolution: **0.1 nm** in nm mode, as specified. In eV mode the equivalent
is 0.001 eV, chosen so the exported grid is comparably fine across the same
range. Sticks export the raw `(position, f)` pairs exactly as listed, no grid —
the same rule the IR sticks export already follows.

## Testing

Against the real captured output in `examples/AbsorptionSpectrum/`:

* 80 roots parsed, first and last energies and `<S**2>` values exact.
* State 1 has exactly one excitation, `89b → 90b`, weight 0.988112, c = −0.994038.
* State 10 has five excitations, in the order printed.
* Absorption row for root 1: 640.8 nm, f = 0.002236399, DZ = 0.18592.
* Root numbers from the table line up with the state block (root 80 = 5.115694 eV).
* `unrestricted` is true and `ground_state_s2` ≈ 8.756797.
* A truncated file (excited states only, table cut off) still parses, with
  `oscillator_strength: None`.
* A file with no TD-DFT block at all returns `Err`.

For the shared spectrum module: the moved line-shape and export functions keep
their existing tests, re-pointed at the new module, plus new tests that
`nice_tick_step` produces 6–10 ticks across representative nm and eV spans, and
that a nm-mode export lands on a 0.1 nm grid.
