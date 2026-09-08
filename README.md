<p align="center">
  <img src="packaging/beavyr.svg" alt="" width="132" height="132">
</p>

<h1 align="center">Beavyr</h1>

A molecular viewer and quantum-chemistry front-end, written in Rust on
[Bevy](https://bevyengine.org/). It displays structures and trajectories at
interactive frame rates, drives external quantum-chemistry programs, and reads
back what they produce — geometries, Hessians, orbitals, excited states — into
the analyses a chemist actually wants: normal modes, thermochemistry, IR and
UV-Vis spectra, isosurfaces, RMSD.

Beavyr links no quantum-chemistry code of its own. Every backend is an
executable you point it at, run in a scratch directory and parsed back from its
output files. That keeps the dependency list to five crates — `bevy`,
`bevy_egui`, `rfd`, `ndarray`, `anyhow` — with no BLAS and no Python.

## Build and run

```sh
cargo run                                       # dev build, dynamic linking, fast rebuilds
cargo build --release --no-default-features     # standalone binary, no dynamic linking
cargo test                                      # 646 passing, 11 ignored
```

## Opening files

Beavyr takes filenames on the command line and loads each into the tool that
owns it:

```sh
beavyr ZZAllyl.xyz            # structure on screen
beavyr trajectory.xyz         # several frames: trajectory, ready to play
beavyr water.hess             # frequencies already analysed
beavyr mo.molden              # orbitals loaded
beavyr excited.out            # excited states, if the file has any
beavyr mol.xyz mo.molden      # several at once, each into its own tool
beavyr --help                 # the list of formats
```

`.xyz`, `.hess`, `.molden` and `.engrad` are recognised by name. `.out` and
`.log` are written by every program in the field, so those are decided by
looking inside: a Gaussian banner means frequencies, an absorption block means
excited states. A file that is neither says so rather than opening nothing —
an ORCA `Freq` output, in particular, is told that its force constants are in
the companion `.hess`.

Anything that fails is reported both in the terminal and at the top of the
Structure panel. Beavyr still opens: a bad argument is not a reason to refuse
to start.

The `dev` feature (on by default) enables Bevy's dynamic linking. Turn it off
for anything you intend to ship. Dependencies are built at `opt-level = 3` even
in a dev profile, because an unoptimised wgpu makes the viewer unusable.

## What it can do

The interface is an icon rail down the side; each icon opens one window. Any
number can be open at once, and each keeps its own state.

### Structure

The atom count, molecular formula in Hill notation (`C12 : H24 : O2`) and
molecular weight sit at the top of the panel — the cheapest way to notice that
a structure is one hydrogen short before spending CPU on it. An unknown element
symbol suppresses the weight and is named, rather than quietly producing a
weight that is too small.

Load or paste XYZ. Edit the coordinates as text and re-apply, or edit the
molecule in 3D. The builder adds and removes atoms, attaches fragments with
chemistry-aware placement (VSEPR geometry and valency, not just a bond vector),
sets distances and angles numerically, rotates about a bond, and undoes each of
those. Sixteen fragments ship with it: `-CH3 -OH -NH2 -Phenyl -Pyrrole -OCH3
-COOH -NO2 -OOH -CCH -CycloPentane -CycloHexane -Cl -Br -I -SH`.

A Z-matrix editor works alongside the Cartesian view — build in internal
coordinates, convert either way.

### Representation and Appearance

Six representations: ball-and-stick, rounded sticks, low-resolution balls and
lines, lines only, backbone trace, space-filling. Per-element visibility
filters, atom and bond scaling, bond radius as a fraction of the smaller
covalent radius, uniform or atom-split bond colours.

Colour schemes: CPK, Jmol, VMD, Molden, Molden0, or a custom scheme you name,
save, and optionally load at startup. PBR material controls (metallic,
roughness, reflectance, emissive), an ambient light and a key/fill/rim rig,
background colour.

Hydrogen bonds are detected on a distance cutoff and drawn as dashed lines with
their own colour, thickness and dash geometry.

### Geometry Optimization

Runs `xtb --opt` or Behemoth in the background — cancellable, non-blocking, the
viewer stays at 60 fps — and hands back the full optimization trajectory to
play through. Methods: GFN1-xTB, GFN2-xTB, GFN-FF for xTB; HF, DFT, MP2,
RI-MP2, TASI, GFN1-xTB for Behemoth, with functional, basis set, dispersion
correction (D3 in seven flavours, D4) and RIJK or exact two-electron
treatment. Charge and multiplicity are validated before the job is sent: xTB
will happily accept an impossible spin state and produce numbers, so Beavyr
checks the electron count itself.

### Vibrations

Either run a frequency calculation (`xtb --hess`, or Behemoth) or load a
Hessian someone else computed. From there:

* **Normal modes** — mass-weighted Hessian, Jacobi diagonalisation, Eckart
  projection to remove translation and rotation. Animate any mode at a chosen
  amplitude and frame rate.
* **Reaction-path projection** — supply a gradient (`.engrad`) and the mode
  along the reaction path is projected out too, which is what you want at a
  transition state.
* **IR spectrum** — intensities from the program's own dipole derivatives,
  broadened as Gaussian, Lorentzian, Voigt or plain sticks, with adjustable
  width and a frequency scaling factor. Exportable as `.dat`.
* **Thermochemistry** — ZPE, thermal corrections, entropy and free energy at a
  chosen temperature, using Grimme's quasi-RRHO treatment of low-lying modes
  with a settable cutoff. The rotational symmetry number must be stated rather
  than guessed; a point-group reference table is built into the panel.

### UV-Vis (TD-DFT)

Reads excited states — energies and oscillator strengths — from an ORCA or
Behemoth output, or launches the calculation itself through Behemoth
(functional, basis, TDA, unrestricted reference, memory). Plots absorption
against wavelength or energy with the same four line shapes, over a settable
energy window, and exports the curve.

### Surface Tools

Loads a `.molden` file, evaluates the Gaussian basis on a grid (Cartesian to
spherical harmonics handled), and meshes an isosurface by marching cubes. Three
quantities: a single molecular orbital, the total electron density, or the spin
density. Adjustable isovalue, grid resolution, opacity, per-lobe colours, and a
wireframe overlay. The sampled field is retained, so dragging the isovalue
slider only re-meshes rather than re-evaluating the basis.

### Trajectory

Plays multi-frame XYZ — MD, optimization paths, IRC — forwards or backwards,
looping, at a chosen frame rate, with a frame slider.

Frames can also be overlaid on the current structure rather than played: a
solid overlay and a transparent "ghost" overlay, each with its own stride, so a
100-frame optimization can be shown as every tenth frame ghosted behind the
final geometry.

### Structure Comparison

RMSD between frames or loaded structures, as placed and after optimal Kabsch
superposition. Pin any frame as the reference, align one frame or all of them,
and see which atom moved furthest.

### Measurements

Distances, angles and dihedrals by clicking atoms. Left-drag orbits the camera,
so a click is distinguished from a drag by how far the pointer travelled.
Per-measurement and global styling, labels drawn next to the line.

### Diagnostics

Flags structures that are probably wrong before you spend CPU on them: close
contacts, isolated atoms, and bond orders that sit between an element's common
valences (a likely missing hydrogen).

### Export Image

PNG or JPEG of the viewport, or of a selected area of it.

## File formats

| Format | Read | Write |
|---|---|---|
| XYZ (single and multi-frame) | ✅ | ✅ |
| Molden | ✅ orbitals, basis, density | ✅ |
| ORCA `.hess` | ✅ Hessian, geometry, masses, frequencies, IR | — |
| ORCA `.engrad` | ✅ Cartesian gradient | — |
| ORCA output (`.out`) | ✅ excited states | — |
| Gaussian `.log` | ✅ frequencies, IR intensities, normal modes¹ | — |
| Turbomole `hessian` / `vibspectrum` (xTB) | ✅ | — |
| Behemoth output | ✅ optimization, Hessian, sTDA/sTD-DFT | — |
| Spectra, thermochemistry | — | ✅ `.dat` |
| Viewport | — | ✅ PNG, JPEG |

¹ A plain Gaussian `Freq` job does not write the force-constant matrix; that
needs `iop(7/33=1)` or a formatted checkpoint. Beavyr reads the frequencies and
modes Gaussian printed.

## External programs

Beavyr *drives*: **xTB** and **Behemoth**. Point it at the executable once and
the path is remembered.

Beavyr *reads output from*: **ORCA**, **Gaussian**. These are jobs you queue and
run yourself — a TD-DFT or DFT Hessian job is not something to launch from a
viewer's button.

## Repository layout

```
src/
  main.rs               Bevy wiring: resources, messages, system sets
  scene.rs camera.rs    rendering, orbit camera, dirty-flag refresh tiers
  picking/              click-vs-drag atom picking
  molecule.rs           topology; bond_order.rs, hbonds.rs perception
  molecule_builder/     fragments, attachment, Z-matrix, bond rotation
  normalmode/           Hessian → modes, Eckart projection, thermochemistry
  vibronic/             Franck-Condon activity, Duschinsky rotation
  orbitals/             Molden → basis evaluation → grid → isosurface
  uvvis/                excited states from TD-DFT output
  spectrum/             shared peak list, broadening and plotting
  qchem_interfaces/     external programs: run them, parse them, validate input
  rmsd/ numerics/       Kabsch superposition, SVD3, 3x3 matrix helpers
  ui*.rs                egui panels and the icon rail
docs/superpowers/specs/ design documents, one per feature
examples/               real quantum-chemistry output used as test material
tests/fixtures/         parser fixtures
```

Two conventions worth knowing before editing:

**Dirty-flag refresh tiers.** `MolSettings` carries `coords_dirty`,
`materials_dirty`, `meshes_dirty`, `geometry_dirty`, `bond_topology_dirty`.
Each names the *smallest* work that will make the view correct again. Setting a
heavier flag than necessary is what once made every slider drag rebuild the
whole scene.

**No BLAS.** `normalmode/linalg_shim.rs` exists so analysis code imported from
Behemoth runs on plain `ndarray` without `ndarray-linalg` and a system BLAS.
Keep it that way.

## Implemented but not yet exposed

`src/vibronic/` is complete and tested but has no panel attached to it yet:

* **Franck-Condon activity** — projects an excited-state gradient taken at the
  ground-state geometry onto the ground-state normal modes, giving each mode's
  displacement, Huang-Rhys factor and resonance-Raman intensity in the
  short-time approximation. It needs no excited-state optimisation, which
  matters when an excited-state minimum is hard to converge.
* **Duschinsky rotation** — how the modes of one electronic state are built
  from the modes of another, for asking which ground-state mode an
  excited-state band descends from.

## Roadmap

`RENDERING_ROADMAP.md` covers the graphics side in detail. The short version:
the camera carries no environment map, so the `metallic` and `reflectance`
sliders cannot physically do anything yet, and the light rig is fixed in world
space so orbiting swings the molecule into shadow. Ambient occlusion and
contact shadows share a depth prepass and are the largest legibility win
available.

Chemistry side:

- Conformer search
- Simulated NMR spectra
- A UI for the vibronic analyses above
- Optimizer with a simple built-in force field
- File associations, so double-clicking an `.xyz` opens Beavyr

## Installing

`packaging/build-deb.sh` builds a Debian package that installs Beavyr into the
applications menu and puts `beavyr` on the PATH:

```sh
./packaging/build-deb.sh
sudo apt install ./target/deb/beavyr_0.2.0_amd64.deb
sudo apt remove beavyr                 # and back out again
```

See `packaging/README.md`.

## License

GPL-3.0. See `LICENSE`.
