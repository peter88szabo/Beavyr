# Frequency analysis (Hessian, normal modes, thermochemistry) — design

Status: pending written-spec review.
Author: assistant, in collaboration with Peter Szabo (peter.szabo@kuleuven.be).

## 1. Purpose

Add vibrational frequency analysis to Beavyr, on top of the xTB geometry
optimizer added earlier (`docs/superpowers/specs/2026-09-04-xtb-geometry-optimizer-design.md`):

1. Run an xTB Hessian calculation on the structure currently on screen
   (normally straight after an optimization).
2. Read the resulting Hessian and run **Behemoth's own** normal-mode
   analysis on it, imported into Beavyr rather than reimplemented — including
   the Eckart vib/rot projection, the reaction-path-out projection, and the
   thermochemistry toolset.
3. List the resulting modes, let the user pick one, and animate it in the
   existing trajectory player (the approach `Smite` uses for the same job).
4. Offer a button — mirroring the optimizer's energy-plot button — that shows
   the thermochemistry summary exactly as Behemoth prints it.

## 2. Background research (facts gathered before designing)

### 2.1 What xTB gives us (verified by running it)

Ran `xtb input.xyz --ohess --chrg 0` (xtb 6.7.1, the binary at
`/home/peter/Programs/orca_6_1_1_linux_x86-64_shared_openmpi418/xtb`) on
`examples/h2o2.xyz`; exit 0, and it wrote:

- **`hessian`** — the Cartesian Hessian. Turbomole-style: a `$hessian`
  header line, then free-format values across lines. Verified count: exactly
  144 numbers for N=4, i.e. a full 3N x 3N matrix, row-major.
- **`vibspectrum`** — a `$vibrational spectrum` block listing, per mode, the
  wavenumber in cm^-1 and the IR intensity in km/mol, ending with `$end`.
  Verified for H2O2: six ~0 modes then 259.56, 1108.86, 1158.96, 1355.45,
  3531.45, 3534.68 cm^-1. This is where IR intensities come from (§4.3).
- **`g98.out`** — normal modes in Gaussian-98 format, useful as an
  independent cross-check of our own eigenvectors.
- xTB also prints its own thermochemistry to stdout (`zero point energy`,
  `G(RRHO)`, a temperature table). Beavyr will **not** use it: the user
  asked for Behemoth's thermochemistry output specifically.

`--hess` (Hessian at the current geometry) and `--ohess` (optimize, then
Hessian) both exist; §4.1 covers which Beavyr uses.

### 2.2 What Behemoth already has (to be imported)

`/home/peter/Dropbox/Research_Leuven/Behemoth/src/normalmode/`:

| file | lines | contents |
|---|---|---|
| `mod.rs` | 484 | `normal_modes`, `normal_modes_with_projection`, `project_hessian`, `project_hessian_vibrot`, `print_frequencies`, `NormalModeResult`, `EckartMode {Off, VibRot, ReactionPath}`, `HessianWeight` |
| `bmat_smat.rs` | 125 | `build_bmat`, `build_smat`, `get_eckart_projector`, `eckart_transform`, `eckart_reactionpath_transform` |
| `inertia.rs` | 62 | `get_brot` (rotational constants from the inertia tensor) |
| `thermofuncs.rs` | 539 | `ThermoResults`, `eval_thermo`, `print_thermo`, `freqs_au_to_cm1`, `brot_from_coords` |
| `print.rs` | 39 | `coord_labels`, `print_hessian_matrix`, `format_hessian_title` |
| `fd_hessian.rs` | 152 | `numerical_hessian_from_gradients` (generic over a gradient closure), `hessian_asymmetry_ratio` |
| `ir_intensity.rs` | 227 | `reduced_masses_amu`, `ir_intensities_fd` (generic over a dipole closure), `ir_intensities_fd_xtb` |

**The dependency problem, and its solution.** Behemoth's `normalmode` calls
`crate::numeric::linalg::{Backend, LinAlg}`, and Behemoth's `Cargo.toml`
pulls `ndarray-linalg` with OpenBLAS/MKL — a system-native BLAS/LAPACK
chain. Beavyr currently depends only on `bevy`, `bevy_egui` and `rfd`, with
no native math libraries, and requiring a system OpenBLAS to build a
molecular *viewer* would be a serious regression in build simplicity.

Resolved by a fact found while reading Behemoth: `LinAlg::eigh` already
falls back to `crate::numeric::jacobi_diag::jacobi` when no BLAS backend is
selected, and **`jacobi_diag.rs` (153 lines) is pure Rust with zero
imports**, operating on `Vec<Vec<f64>>`. So Beavyr imports that solver and
provides a small local `LinAlg` shim — leaving the scientific code
essentially byte-identical while dropping the BLAS requirement entirely.

The shim's required surface is small and was measured, not guessed: grepping
every `LinAlg` call in `src/normalmode/*.rs` shows the imported code uses
exactly four things — `LinAlg::new`, `la.eigh`, `la.matmul`, `la.transpose`.
`eigh` is backed by the imported `jacobi_diag`, and `matmul` can reuse
Behemoth's own pure-Rust `LinAlg::matmul_pure` (`numeric/linalg.rs:271`),
so even the shim is largely imported rather than written.

Hessians here are 3N x 3N (a 200-atom system is 600 x 600), well within
Jacobi's practical range.

**New dependencies:** `ndarray` and `anyhow` only. Both are pure Rust with
no native/system libraries, and keeping them means the imported files need
almost no edits. (`ndarray-linalg`, the part that needs BLAS, is *not*
added.)

### 2.3 What `ir_intensity.rs` cannot bring with it

The user asked to import everything including `ir_intensity.rs`. Most of it
comes across, but one function honestly cannot:

- `reduced_masses_amu` — pure math, imports cleanly.
- `ir_intensities_fd<F>` — generic over a `dipole_fn` closure, imports
  cleanly; usable as soon as anything can supply dipoles.
- `ir_intensities_fd_xtb` — a four-line wrapper that supplies that closure
  from `crate::hamiltonian::xtb_hamiltonian` (Behemoth's *own* xTB
  implementation: `XtbParams`, `SpinPol`, `dipole_from_xtb_geometry`, plus
  `crate::properties::polarizability_finitefield::dipole_from_scf` and
  `crate::molecule::Molecule`). Beavyr shells out to the xtb *binary* and
  has no internal xTB Hamiltonian or SCF, so importing this one wrapper
  would mean dragging in Behemoth's Hamiltonian, properties and SCF stacks
  (`gfn1_baked_params.rs` alone is 23,570 lines). **It is therefore dropped,
  and only it.** IR intensities come from xTB's own `vibspectrum` instead
  (§4.3), which is the same number by a cheaper route.

`fd_hessian.rs` imports cleanly in full (it is generic over a gradient
closure). It stays unused until something wires a gradient provider to it —
`xtb --grad` could, later — and is imported now as requested.

### 2.4 How Smite animates a mode (the UX reference)

`Smite/src/normalmode/nmodeprint.py::print_single_mode` writes a multi-frame
XYZ where each frame is `q_cart = q_eq + L · (A cos θ)` for θ over one full
period (600 steps by default) — i.e. a sinusoidal displacement along the
selected mode's eigenvector. Beavyr already has a trajectory player with
frame-by-frame playback, so mode animation is generated in memory as frames
and handed to `TrajectoryState`; no new rendering code is needed.

## 3. Non-goals

- No internal xTB/SCF implementation (see §2.3); Beavyr keeps shelling out
  to the xtb binary.
- No IRC/transition-state search. The reaction-path projection is imported
  and wired, but Beavyr does not itself walk a reaction path.
- No `ndarray-linalg`/BLAS dependency (see §2.2).
- Behemoth's `print_frequencies`/`print_thermo` keep writing to stdout as
  they do today; the GUI gets the same text through a `format_*` sibling
  (§4.4) rather than by rewriting them.

## 4. Design

### 4.1 Running the Hessian

A new "Frequencies" section in the same right-hand panel as the optimizer,
reusing its machinery wholesale: the same background-task shape
(`AsyncComputeTaskPool` + per-frame `poll_once`), the same cancel handling,
the same live elapsed-time/progress readout, the same run-directory rules
(discard on success, keep on failure for `xtb.stdout`/`xtb.stderr`), and the
same validated charge/multiplicity from the optimizer's panel state (one
electronic state for the molecule, not two).

- Button: **Run Frequencies** — invokes `xtb input.xyz --hess --chrg <c>
  [--uhf <u>]` on the geometry currently on screen. `--hess`, not `--ohess`:
  the panel already has its own Optimize button, and silently re-optimizing
  under a "frequencies" button would move the user's structure without being
  asked. A frequency analysis is only meaningful at a stationary point, so
  the panel says so when the structure has not been optimized in this
  session.
- When the projection mode is **Reaction path**, the run additionally
  requests a gradient (`--grad`), since `eckart_reactionpath_transform`
  requires one. If the resulting gradient norm is essentially zero — the
  normal case at a converged minimum, where Behemoth's own code would
  `bail!("gradient norm is zero")` — the panel states that plainly and falls
  back to Vib-Rot rather than surfacing an error, per the user's choice.

### 4.2 Reading the Hessian

New `src/qchem_interfaces/hessian_file.rs`:

```rust
/// Parses xTB's Turbomole-style `hessian` file: a `$hessian` header, then
/// 3N*3N free-format values, row-major.
pub fn parse_turbomole_hessian(text: &str, natoms: usize) -> Result<Array2<f64>, String>;

/// Parses xTB's `vibspectrum`: (wavenumber_cm1, ir_intensity_km_per_mol) per mode.
pub fn parse_vibspectrum(text: &str) -> Result<Vec<(f64, f64)>, String>;
```

Both are pure functions over text, unit-tested against the **real captured
output** from the §2.1 run (saved as a fixture) rather than hand-written
stubs — including the exact H2O2 frequencies listed there.

### 4.3 Analysis

`src/normalmode/` in Beavyr, holding the imported Behemoth files plus the
shim:

```
src/normalmode/
  mod.rs           // imported from Behemoth, re-exports the public API
  bmat_smat.rs     // imported: Eckart + reaction-path projections
  inertia.rs       // imported
  thermofuncs.rs   // imported, plus `format_thermo` (§4.4)
  print.rs         // imported (+ the 70-line print_matrix_in_ao utility it needs)
  fd_hessian.rs    // imported, currently unused (§2.3)
  ir_intensity.rs  // imported minus `ir_intensities_fd_xtb` (§2.3)
  jacobi_diag.rs   // imported from Behemoth's numeric/, pure Rust eigensolver
  linalg_shim.rs   // NEW: `LinAlg` with the four methods the imported code
                   //      calls (new/eigh/matmul/transpose), eigh backed by
                   //      jacobi_diag and matmul by Behemoth's matmul_pure
```

Masses come from a new atomic-mass table keyed by element symbol (Beavyr has
`atomic_number` and `covalent_radius_angstrom` but no masses); `linear` is
determined from the inertia tensor via the imported `inertia.rs`. IR
intensities are attached to modes from `vibspectrum` (§4.2), matched by mode
order, and shown as a column beside each frequency.

### 4.4 Thermochemistry summary

Behemoth's `print_thermo(&ThermoResults, temp, freq_cutoff, elec_energy)`
writes a formatted block to stdout. To show *exactly that text* in a window,
`print_thermo`'s body moves into a new `format_thermo(...) -> String` and
`print_thermo` becomes `println!("{}", format_thermo(...))` — so the printed
and displayed forms cannot drift apart. A **Show Thermochemistry** button,
mirroring the optimizer's energy-plot button, opens a window rendering that
string in a monospace font. Temperature (default 298.15 K) and the
low-frequency cutoff are exposed as panel inputs.

### 4.5 Mode list and animation

- The mode list shows index, frequency (cm^-1), IR intensity, and a marker
  for imaginary modes (negative eigenvalues — `NormalModeResult` already
  separates `negative_indices`/`zero_indices`/`positive_indices`).
- Selecting a mode generates animation frames in memory following Smite's
  formula (§2.4), `q_eq + L·A cos θ`, with a user-adjustable amplitude and a
  fixed frame count, and loads them into `TrajectoryState` with
  `PlaybackMode::Loop` — a vibration genuinely repeats, unlike an
  optimization path, which is why the optimizer forces `Once` and this
  forces `Loop`.
- Per the user's choice this **overwrites** any loaded trajectory (e.g. the
  optimization path), and the panel says so, rather than prompting each time.

## 5. Testing plan

- `parse_turbomole_hessian` / `parse_vibspectrum`: against the real captured
  `hessian` and `vibspectrum` fixtures from §2.1 — asserting 144 values for
  H2O2, symmetry of the parsed matrix, and the exact frequency list.
- **Numerical agreement with xTB**: the decisive test. Run Behemoth's
  imported `normal_modes` on the parsed H2O2 Hessian and assert the
  resulting frequencies match `vibspectrum`'s own values to within a small
  tolerance. This validates the whole chain — parser, mass table, Jacobi
  eigensolver shim, and Eckart projection — against an independent
  implementation, rather than merely testing that the code runs.
- `jacobi_diag`: eigenvalues/vectors of small known matrices, and
  reconstruction (`A ≈ V Λ Vᵀ`).
- `thermofuncs`: `eval_thermo` on H2O2's frequencies; assert the zero-point
  energy is close to xTB's own reported `zero point energy` value from the
  same run (0.024943563139 Eh, §2.1) as an external cross-check.
- Eckart projection: assert that projecting a minimum's Hessian yields six
  near-zero modes (five for a linear molecule), and that reaction-path mode
  refuses a zero gradient the way Behemoth's code does.
- Mode animation: frame count, that frame 0 equals the equilibrium geometry,
  and that the maximum displacement respects the requested amplitude.
- The xtb Hessian run itself: an integration test gated on a resolvable xtb
  binary (`BEAVYR_XTB_PATH`/PATH), skipped with a printed reason otherwise —
  matching the optimizer's existing pattern.

## 6. Risks

- Jacobi is O(n³) with a larger constant than LAPACK; a very large system
  (>300 atoms, 900+ coordinates) will be noticeably slower to diagonalize.
  It runs on the background task pool, so the UI stays responsive, and the
  panel reports progress. If it ever becomes a real limit, the shim is the
  single place a faster backend would slot in.
- The imported files carry Behemoth's `anyhow`-based error style, which
  differs from Beavyr's `Result<_, String>` convention. Deliberately kept
  as-is inside `src/normalmode/` (the point is to import, not rewrite), with
  conversion at the boundary where Beavyr's UI calls in.
- Mass table and unit conventions (amu vs. atomic units, bohr vs. Angstrom)
  are the most likely source of a silent factor-of-something error, which is
  exactly why §5 cross-checks the final frequencies against xTB's own.
