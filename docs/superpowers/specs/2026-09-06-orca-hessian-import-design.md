# Import an ORCA `.hess` File — Design

Date: 2026-09-06
Status: approved for implementation

## Goal

A "Load Hessian…" button beside "Run Frequencies" in the Vibrations panel.
It reads a Hessian someone already computed elsewhere and runs it through the
*same* analysis the xTB path uses: Eckart projection, our own normal modes,
frequencies, IR spectrum, mode animation and thermochemistry.

ORCA's `.hess` is the only format for now. Reference file:
`examples/TS_Gamma-3-6_ZZ-S_11_Compound_1.hess` — 19 atoms, 57 coordinates, a
transition state with one imaginary mode at −602.6 cm⁻¹.

## Why this is small

`analyze()` already does all the work, from a `RawFrequencyOutput`:

```rust
struct RawFrequencyOutput {
    atoms: Vec<String>,
    coords_bohr: Vec<f64>,
    hessian: Array2<f64>,
    vib_lines: Vec<VibSpectrumLine>,
    gradient_bohr: Option<Vec<f64>>,
}
```

So the whole feature is: a parser that fills that struct from a `.hess`, plus a
button. Nothing downstream changes — the mode list, animation, thermochemistry
window and IR spectrum all work on a `FrequencyResult` and cannot tell where it
came from.

## The `.hess` format

Sections are `$name` lines; the file ends at `$end`. Only these are read:

| Section | Content |
|---|---|
| `$hessian` | `n`, then column blocks of five: an index header row, then `n` rows of `row_index  v0 v1 v2 v3 v4`. Units Eh/bohr². |
| `$atoms` | `n`, then `symbol  mass_amu  x  y  z` per atom. **Coordinates are in bohr.** |
| `$vibrational_frequencies` | `n`, then `index  wavenumber_cm1`. ORCA's own, already projected: the six translations/rotations are exactly zero. |
| `$ir_spectrum` | `n`, then `freq  eps  intensity_km_mol  TX  TY  TZ`. |
| `$multiplicity`, `$act_energy`, `$frequency_scale_factor` | scalars, read for display only. |

`$normal_modes` and `$dipole_derivatives` are parsed past but not used: we
compute our own modes, which is the point of running the file through our
pipeline rather than just displaying ORCA's numbers.

The block layout of `$hessian` is confirmed against the reader in
`Smite/src/qchem_interfaces/orcarun.py::readOrcaHessian`, which is the only
ORCA-specific section Smite reads — it gets geometry and masses from its caller
rather than from the file, so `$atoms` and `$ir_spectrum` are new here.

### Mode ordering — the one real trap

ORCA lists modes with the six projected translations/rotations **first** (as
exact zeros) and the imaginary mode **seventh**:

```
0..5   0.00
6      -602.61
7      106.90   ...
```

Our eigensolver returns modes ascending by eigenvalue, so the imaginary one
comes **first**, then the near-zeros, then the positives. The xTB path pairs
intensities to modes by position, which is valid there because xTB's
`vibspectrum` uses the same ascending order — but doing that here would put the
imaginary mode's intensity on a translation.

Fix: sort ORCA's `(frequency, intensity)` pairs ascending by frequency before
handing them over. Sorting −602.61 to the front and leaving the zeros in the
middle reproduces our eigensolver's order exactly. A test asserts the
intensities land on the right modes by comparing against ORCA's own frequency
list.

### Masses

`$atoms` carries a mass per atom. `analyze()` currently derives masses from the
element symbol, which is right for xTB but would silently ignore an isotope
substitution recorded in the file. `RawFrequencyOutput` gains
`masses_amu: Option<Vec<f64>>`; `analyze()` uses it when present and falls back
to the symbol table otherwise, so the xTB path is unchanged.

## Code

* **`src/qchem_interfaces/orca_hess.rs`** (new)
  * `pub struct OrcaHessian { natoms, atoms, masses_amu, coords_bohr, hessian, frequencies_cm1, ir_intensities_km_mol, multiplicity, energy_hartree, frequency_scale_factor }`
  * `pub fn parse_orca_hess(text: &str) -> Result<OrcaHessian, String>`
  * private: `find_section`, `parse_block_matrix`, `parse_atoms`,
    `parse_indexed_values`, `parse_ir_spectrum`, `parse_scalar`

  Errors name the section and what was expected. A missing `$hessian` or
  `$atoms` is fatal; a missing `$ir_spectrum` or `$vibrational_frequencies`
  is not — intensities become zero and the analysis still runs.

* **`src/qchem_interfaces/xtb_freq.rs`**
  * `RawFrequencyOutput` gains `masses_amu: Option<Vec<f64>>`.
  * `pub fn load_hessian_file(path, multiplicity, eckart, temp_k, cutoff) -> Result<FrequencyResult, String>`
    reads the file, dispatches on content, builds the raw struct, calls `analyze`.
  * The panel gets a **Load Hessian…** button next to Run Frequencies, plus the
    loaded file's name in the status line. It also loads the file's geometry
    into the viewport, since the program now starts empty and a Hessian with
    nothing on screen would look broken.

  Loading is synchronous: parsing a 150 KB file and diagonalising a 57×57
  matrix is milliseconds, so it needs none of the background-task machinery the
  xTB runs use.

* **`src/ui.rs`** — passes `&mut Molecule` and the `MoleculeChanged` writer to
  the panel so the loaded geometry can reach the viewport.

## Reaction-path projection

A `.hess` carries no gradient, so `EckartMode::ReactionPath` has nothing to
project against. `analyze()` already handles that: `gradient_bohr: None` falls
back to the ordinary Eckart projection with a note explaining why. Correct as
is — a `.hess` is written at a stationary point anyway.

## Testing

Against the real example file:

* 19 atoms, 57×57 Hessian, symmetric to 1e-10.
* Element symbols and masses match `$atoms`; the first atom's coordinates are
  the bohr values in the file, and the H–C distance comes out ≈1.08 Å.
* Our computed frequencies reproduce ORCA's own `$vibrational_frequencies` for
  every non-zero mode, including the −602.6 cm⁻¹ imaginary one, to within a few
  tenths of a cm⁻¹ — the real cross-validation, mirroring the xTB test that
  checks our pipeline against `vibspectrum`.
* Exactly one imaginary mode, six near-zero modes, fifty positive.
* IR intensities land on the right modes after the ascending sort: the
  imaginary mode's intensity is 0.0 and the 3781.37 cm⁻¹ mode's is 133.38 km/mol.
* Thermochemistry runs and gives a positive zero-point energy.
* A truncated file (no `$ir_spectrum`) still analyses, with zero intensities.
* A file with no `$hessian` is an error naming the section.

## Later (not in this change)

Reading an ORCA **frequency-calculation output** (`.out`) directly, taking its
frequencies, IR intensities and thermochemistry as printed rather than
recomputing anything. The thermochemistry panel would then show only what such
a file actually reports — total energies, enthalpies, entropies and the ZPE —
not the electronic / vibrational / rotational / translational breakdown we
produce ourselves, which is only available when we do the analysis.

That means `FrequencyResult` will need to distinguish "computed here" from
"read from a report", so the thermochemistry window can show the fuller table
only when it has one. Requested 2026-09-06; deferred deliberately.

**Animation is possible from a `.out`**, contrary to first impression. The
output carries a `NORMAL MODES` block with the full 3N x 3N matrix of
displacement vectors (zeros in the six projected columns, real values from
column 6 on) as well as `CARTESIAN COORDINATES (ANGSTROEM)`. So a `.out` has
frequencies, intensities, thermochemistry *and* modes to animate.

The care point is the convention: ORCA says those vectors are "the Cartesian
displacements weighted by the diagonal matrix M(i,i)=1/sqrt(m[i])... normalized
but *not* orthogonal", which is not necessarily the normalisation our
`animate_mode` assumes. Getting it wrong gives an animation that looks
plausible but has the wrong relative amplitudes -- heavy atoms moving too far
or not far enough. It is directly checkable: the `.hess` and the `.out` in
examples/ describe the same modes, so animating a mode from each must produce
the same displacements.

### Loading a `.hess` and a `.out` together

The two files are authoritative for different things, so loading both is not
redundant — it is the best case:

| | from the `.hess` | from the `.out` |
|---|---|---|
| Normal modes / animation | our own, recomputed | ORCA's, as printed |
| Frequencies | ours, at any scale factor | ORCA's, fixed |
| Thermochemistry | ours, at **any temperature and cutoff** | ORCA's, at the one temperature it ran |
| IR intensities | only if `$ir_spectrum` is present | always |

So when both are loaded: the Hessian drives the modes, the animation and our
own thermochemistry — because only it can be re-analysed at a different
temperature, cutoff or scaling factor — while the output supplies IR
intensities when the Hessian has none, and its own reported thermochemistry
shown alongside ours as a labelled reference rather than merged into it.
Two numbers that disagree should be visibly two numbers.

### Output alone: no animation

**Decided requirement:** a `.out` loaded on its own provides exactly three
things — the **frequency list**, the **IR intensities**, and the
**thermochemistry** (only the totals it prints: energies, enthalpies,
entropies, ZPE). Nothing else. There are no normal modes and therefore no
animation: the Animate buttons are disabled and the panel says "Animation
needs the Hessian; load the matching .hess file."

This is a deliberate policy choice, not a limitation discovered in the format:
the example `.out` does carry a populated `NORMAL MODES` block. But ORCA does
not always print it (reduced print levels omit it), the vectors' mass-weighting
convention differs from ours, and treating the Hessian as the one source of
animatable modes keeps a single well-tested path instead of two that can
disagree. If we later want animation from an output alone, the block is there
to read and this warning is where it would be lifted.

**The pairing must be verified, not assumed.** Nothing stops a user picking a
`.hess` and a `.out` from different calculations, and the result would be a
mode list silently labelled with another molecule's intensities. Before
combining, check that the atom count and element sequence match and that the
frequencies agree to within a cm⁻¹ or so; if they do not, keep them separate
and say why. The example pair in `examples/` is a good regression case in both
directions — they match, and a deliberately mismatched pair must be rejected.

Note the example `.out` contains **two** frequency calculations (imaginary
modes at −604.57 and −602.61 cm⁻¹); the `.hess` corresponds to the second. A
reader must take the last one, or match on frequencies, rather than the first
it encounters.
