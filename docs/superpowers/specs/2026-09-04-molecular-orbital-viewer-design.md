# 3D Molecular Orbital Viewer

**Date:** 2026-09-04
**Status:** Approved for implementation
**Component:** `src/orbitals/` (new subsystem in Beavyr)

## 1. Purpose

Add a molecular orbital viewer to Beavyr. The user loads a `.molden` file,
picks an orbital from a list showing its energy and occupation, and sees the
orbital rendered as a two-lobe isosurface in the existing 3D scene. An isovalue
slider controls how much of the orbital is enclosed by the surface.

Closed-shell files present a single orbital set. Open-shell files present
separate alpha and beta sets, selected independently.

## 2. Scope

### In scope

- Parsing the Molden format: `[Atoms]`, `[GTO]`, `[5D]/[7F]/[9G]`, `[MO]`
- Gaussian basis evaluation on a 3D grid, angular momentum s through g
- Closed-shell (single set) and open-shell (alpha + beta) orbitals
- Orbital list with index, symmetry label, energy, occupation, HOMO/LUMO marker
- Marching-cubes isosurface for the +psi and -psi lobes
- User-controlled isovalue
- Replacing the displayed molecule with the geometry from the Molden file

### Out of scope for this version

- Total electron density and spin density (the machinery generalises to them,
  but they are not built now)
- Slater-type orbitals (`[STO]`), natural orbitals, transition densities
- Cube-file import or export
- Multiple simultaneous orbitals

## 3. Constraints and observations from the reference material

### 3.1 The example file

`examples/Exc_TPSSh-D4-def2-TZVP_SCMwater_16.molden`, 15 MB, written by
`orca_2mkl`:

| Property | Value | Consequence |
|---|---|---|
| `[Atoms] AU` | Bohr | must convert; Beavyr is Angstrom throughout |
| `[5D] [7F] [9G]` | pure spherical | Cartesian-to-spherical transform required |
| 546 alpha + 546 beta | open shell | the alpha/beta requirement is exercised immediately |
| 43 atoms including Fe | up to g functions | evaluator must cover l = 0..4 |
| MO coefficients | dense `index value` pairs | straightforward parse |

15 MB is small enough to parse eagerly into memory. No streaming.

### 3.2 Closed versus open shell

Follow Molden's own rule from `rdmolf.f`: a file is open shell if and only if
any `Spin= Beta` block appears (Molden sets `iuhf = 1`). Molden distinguishes
the two by token length, `nstr == 5` for `Alpha` and `nstr == 4` for `Beta`.

Represent this as `beta: Vec<MolecularOrbital>` being empty for closed shell,
rather than a separate enum. `beta.is_empty()` then drives whether the UI shows
a spin selector at all.

### 3.3 Ordering conventions

Three orderings meet in this subsystem and disagree:

| | Cartesian order | Spherical order |
|---|---|---|
| Schlegel generator (from Behemoth) | `lx` descending, then `ly` descending: `xx, xy, xz, yy, yz, zz` | `m = -l ..= +l` |
| Molden file format | `xx, yy, zz, xy, xz, yz` | `0, +1, -1, +2, -2` |

Both require an explicit permutation. This is the single most likely source of
silent, plausible-looking wrong output, so the mapping lives in exactly two
functions (`molden_spherical_row`, `molden_cartesian_column`) and is covered by
its own tests.

### 3.4 Normalisation

Three conventions must be reconciled: primitive normalisation `N(alpha, l)`,
contraction renormalisation, and the per-Cartesian-component factor. Molden
applies the latter via its `pnorm` table in `gaussian.f` (sqrt(3), sqrt(5),
sqrt(15), sqrt(7), ...). Behemoth's `cartesian_to_spherical.rs` explicitly
documents that it uses a *different* Cartesian normalisation and that the
common reference table "cannot be copied directly".

Contraction coefficients in the example are pre-scaled by the producer
(`300784.84637  2.2063712728`), so they cannot be assumed raw.

This is the main technical risk. It is handled by measurement rather than
assumption; see section 7.

## 4. Architecture

### 4.1 Module layout

```
src/orbitals/
  mod.rs          plugin wiring, resources, systems, public types
  cart2sph.rs     Cartesian-to-spherical transform (adapted from Behemoth)
  basis.rs        shell/primitive model, normalisation, AO layout
  molden.rs       .molden parser -> OrbitalData
  evaluate.rs     contracted Gaussian AO values; psi = sum c_i phi_i
  grid.rs         bounding box, sampling, ScalarField
  isosurface.rs   marching cubes -> positions, normals, indices
  ui.rs           egui panel
```

Seven focused files rather than one large one. `builder_ui.rs` at 2,746 lines
is the cautionary example already in this repository.

### 4.2 Reuse from Behemoth

A targeted lift, not a dependency. Beavyr stays standalone.

Taken: the Schlegel and Frisch (1995) coefficient generator from
`src/integrals/cartesian_to_spherical.rs`, with its integral-block machinery
removed and its `nfac` / `n_over_k` / `n_cart` helpers reimplemented locally
(a few lines each).

Not taken: the remaining ~20,000 lines of `src/integrals`, which compute
integrals rather than point values and are irrelevant here.

Not available to take: there is no existing "evaluate a contracted Gaussian at
a point" routine in Behemoth. That is new code, though it is the easy part.

### 4.3 Data model

```rust
pub struct MoldenShell {
    center: usize,           // atom index
    l: u8,                   // 0=s .. 4=g
    spherical: bool,         // resolved per-l from [5D]/[7F]/[9G]
    exponents: Vec<f64>,
    coefficients: Vec<f64>,
}

pub struct MolecularOrbital {
    symmetry: String,        // "1a"
    energy: f64,             // Hartree
    occupation: f64,
    coefficients: Vec<f64>,  // length n_basis
}

pub struct OrbitalData {
    atoms: Vec<String>,
    positions: Vec<Vec3>,    // Angstrom, converted if the header says AU
    shells: Vec<MoldenShell>,
    n_basis: usize,
    alpha: Vec<MolecularOrbital>,
    beta:  Vec<MolecularOrbital>,   // empty => closed shell
}
```

### 4.4 Evaluation strategy

Loop over shells, not over grid points. Each shell gets a bounding sub-box from
its most diffuse exponent (`exp(-alpha r^2) < 1e-8` implies
`r > sqrt(18.4 / alpha)`) and accumulates only inside that box. Combined with
skipping AOs whose MO coefficient is negligible -- for a localised orbital most
of the 546 coefficients are around 1e-7 -- this is what makes the example
tractable.

Parallelised over z-slabs on `ComputeTaskPool`.

### 4.5 Data flow

```
pick orbital ---> AsyncComputeTaskPool::spawn(sample_field)
                       |  viewport stays at full framerate
                       v
                  ScalarField cached (current orbital only, ~2 MB)
                       |
    isovalue --------> marching cubes (~20 ms, synchronous)
                       v
               two Mesh3d entities: +psi blue, -psi red
```

Retaining the field for the currently displayed orbital costs nothing -- it
must exist to be meshed at all -- and makes the isovalue slider live rather
than deferred to drag release.

Surface normals come from the field gradient by central differences, not from
face normals. Cheap, and considerably better looking.

### 4.6 Grid sizing

Bounding box of the atoms plus 3.0 Angstrom padding, default step 0.20
Angstrom, total point count capped (2 million) to bound worst-case time.
Auto-sized; not exposed in the UI.

### 4.7 Integration with existing Beavyr

- Loading a `.molden` replaces the `Molecule` and emits
  `MoleculeChanged::parse_xyz(true)`. Measurement clearing and camera auto-fit
  are then handled by existing systems with no new code.
- Orbital lobes are ordinary `Mesh3d` on `LAYER_MAIN`, so they depth-interleave
  correctly with atoms and bonds and export through the existing screenshot
  path unchanged.
- Surface state gets its own dirty tier (`field_dirty` / `surface_dirty`),
  following the ladder established in `MolSettings`.

## 5. User interface

A new collapsing "Orbitals" section in the existing side panel.

```
Orbitals
  [ Load .molden file ]
  File: Exc_TPSSh-D4-def2-TZVP_SCMwater_16.molden
  546 basis functions, open shell

  Spin: (o) Alpha  ( ) Beta

    269   a    -0.3145 Eh   occ 1.00
    270   a    -0.2876 Eh   occ 1.00
    271   a    -0.2314 Eh   occ 1.00   <- HOMO
    272   a     0.0187 Eh   occ 0.00   <- LUMO
    273   a     0.0412 Eh   occ 0.00

  Isovalue  0.0500  [====|--------]
  [x] Show surface     Opacity [========|--]
```

- The spin selector is hidden entirely for closed-shell files, where
  occupations read `2.00`. That is the visual cue for a single set.
- HOMO and LUMO are derived from occupation, computed **per spin set**: in the
  open-shell example the alpha and beta HOMOs are different orbitals.
- A computation in progress shows a spinner; the viewport stays interactive.

## 6. Error handling

`Result<OrbitalData, MoldenError>` carrying line numbers, surfaced in the panel
the way `xyz_buf.warning` already is. Explicit failures, no silent fallbacks:

- unsupported angular momentum (l > 4)
- basis size versus MO coefficient length mismatch
- malformed or missing sections
- no `[MO]` block

## 7. Validation strategy for normalisation

The risk in section 3.4 is addressed by two checks that convert guesswork into
measurement.

**Structural.** `n_basis` derived from the shell list must equal the MO
coefficient vector length. Getting 5D versus 6D wrong makes this mismatch
immediately, so the failure is loud at parse time rather than a convincing but
wrong picture.

**Numerical.** The integral of `|psi|^2` over a fine grid must be approximately
1. If it returns 3.0 or 1/15, the ratio names the missing factor directly.
Orthogonality, integral of `psi_i psi_j` approximately 0 for i != j, catches
ordering and transform errors that self-overlap alone would miss.

Both are tests. The norm is additionally shown in the panel as a health
readout, so a convention problem with an unusual basis is visible to the user
rather than silent.

## 8. Testing

Test-driven, using small hand-written fixtures rather than the 15 MB file.

| # | Test | Catches |
|---|---|---|
| 1 | parse closed-shell H2O/STO-3G fixture | basic format |
| 2 | parse open-shell alpha+beta fixture | the spin split |
| 3 | `n_basis` versus coefficient length | 5D/6D confusion |
| 4 | s/p/d AO values against analytic expressions | evaluator |
| 5 | integral of psi^2 approximately 1 | **normalisation** |
| 6 | integral of psi_i psi_j approximately 0 | ordering, transform |
| 7 | marching cubes on an analytic sphere, area approximately 4 pi r^2 | mesher |
| 8 | transform rows orthonormal (`T Tt = I`) | coefficient generator |
| 9 | Molden slot permutations are bijective | ordering tables |
| 10 | the real example loads, 546+546, HOMO labelled | end to end |

Test 5 is the one that matters. Tests 7, 8 and 9 are self-validating and need
no external reference data.

## 9. Estimated size

Roughly 1,800 to 2,200 lines including tests, across seven files. The parser
and the mesher are mechanical. The normalisation reconciliation is the genuine
unknown and should be expected to need one or two empirical iterations against
test 5.

## 10. Decisions taken

| Decision | Choice | Rejected alternative |
|---|---|---|
| Isosurface rendering | CPU marching cubes into a Bevy `Mesh` | GPU raymarching a 3D texture; gizmo point cloud |
| Basis code | targeted lift of Behemoth's spherical transform | path dependency on Behemoth; writing everything fresh |
| Isovalue control | single slider, auto-sized grid | exposing box padding and grid step; enclosed-density percentage |
| Scope | molecular orbitals only | orbitals plus density plus spin density |
| Responsiveness | background compute, field cached for the current orbital | blocking; multi-orbital cache |
| Frame pacing | unchanged, stays `Continuous` | reactive redraw (rejected earlier for playback smoothness) |
