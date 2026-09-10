<p align="center">
  <img src="packaging/beavyr.svg" alt="" width="132" height="132">
</p>

<h1 align="center">Beavyr</h1>

Beavyr is a molecular viewer for chemists. Load a structure or build one from
scratch, watch a trajectory or a reaction play out, and see what a
quantum-chemistry calculation actually looked like — normal modes, IR and
UV-Vis spectra, molecular orbitals, conformers — all in one fast, interactive
window.

## Download

Get the current release from the
[Releases page](https://github.com/peter88szabo/Beavyr/releases):

| Platform | File | Install |
|---|---|---|
| Windows | `beavyr-setup.exe` | Run it |
| macOS | `beavyr-macos.dmg` | Open it, drag Beavyr to Applications |
| Debian, Ubuntu, Mint | `beavyr_amd64.deb` | `sudo apt install ./beavyr_amd64.deb` |
| Arch, Manjaro | build with `packaging/PKGBUILD` | `makepkg -si` |

No Python, no separate runtime, nothing else to install. Beavyr starts and
shows the viewer with no structure loaded — it does not need xTB or any other
program to open.

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

## What you can do

The interface is an icon rail down the side; each icon opens one window. Any
number can be open at once, and each keeps its own state.

### Molecule Editor

The atom count, molecular formula in Hill notation (`C12 : H24 : O2`) and
molecular weight sit at the top of the panel — the cheapest way to notice that
a structure is one hydrogen short before doing anything else with it. An
unknown element symbol suppresses the weight and is named, rather than
quietly producing a weight that is too small.

Load or paste XYZ. Edit the coordinates as text and re-apply, or edit the
molecule directly in 3D. The builder adds and removes atoms, attaches
fragments with chemistry-aware placement (VSEPR geometry and valency, not just
a bond vector), sets distances and angles numerically, rotates about a bond,
and undoes each of those. Sixteen fragments ship with it: `-CH3 -OH -NH2
-Phenyl -Pyrrole -OCH3 -COOH -NO2 -OOH -CCH -CycloPentane -CycloHexane -Cl -Br
-I -SH`.

A Z-matrix editor works alongside the Cartesian view — build in internal
coordinates, convert either way.

**Clean up geometry** relaxes bond lengths and angles a fragment left
strained, using a force field built into Beavyr — no external program, no
waiting. It is a starting structure for a real optimisation, not a substitute
for one, and it says so.

Right-click the background for the whole-molecule tools: recentre the view,
toggle hydrogen bonds, and **Orient / Mirror / Flip** the structure into a
chosen coordinate plane or about a chosen axis. Orient and Flip keep a chiral
molecule as it was; Mirror deliberately gives you its enantiomer.

### Trajectory

Plays multi-frame XYZ — MD, optimisation paths, IRC, a set of conformers — 
forwards or backwards, looping, at a chosen frame rate, with a frame slider.

Frames can also be overlaid on the current structure rather than played: a
solid overlay and a transparent "ghost" overlay, each with its own stride, so
a 100-frame optimisation can be shown as every tenth frame ghosted behind the
final geometry.

### Vibrations

Animate any normal mode at a chosen amplitude and frame rate — the fastest way
to see whether a calculated frequency is the stretch, bend, or torsion you
expected. Run a frequency calculation yourself (through xTB, or the built-in
force field for a quick check) or load a Hessian someone else computed.

* **Normal modes** — mass-weighted Hessian, Eckart projection to remove
  translation and rotation.
* **Reaction-path projection** — supply a gradient and the mode along the
  reaction path is projected out too, which is what you want at a transition
  state.
* **IR spectrum** — intensities from the program's own dipole derivatives,
  broadened as Gaussian, Lorentzian, Voigt or plain sticks, with adjustable
  width and a frequency scaling factor. Exportable as `.dat`.
* **Thermochemistry** — ZPE, thermal corrections, entropy and free energy at a
  chosen temperature, using Grimme's quasi-RRHO treatment of low-lying modes
  with a settable cutoff. The rotational symmetry number must be stated rather
  than guessed; a point-group reference table is built into the panel.

### Orbital Viewer

Loads a `.molden` file, evaluates the Gaussian basis on a grid, and meshes an
isosurface by marching cubes. Three quantities: a single molecular orbital,
the total electron density, or the spin density. Adjustable isovalue, grid
resolution, opacity, per-lobe colours, and a wireframe overlay. The sampled
field is retained, so dragging the isovalue slider only re-meshes rather than
re-evaluating the basis from scratch.

### UV-Vis

Reads excited states — energies and oscillator strengths — from an ORCA
TD-DFT output, and plots the absorption spectrum against wavelength or energy
with the same four line shapes UV-Vis needs, over a settable energy window.
Export the curve for a figure.

### Conformational Analysis

Finds the low-energy shapes a flexible molecule can take, searching over its
rotatable bonds and, for a five- or six-membered ring, its pucker — a chair
against a twist-boat, or a sugar's envelopes and twists, not only what one
bond does. Runs on the built-in force field by default (seconds, no external
program), or on xTB — GFN-FF, GFN1, or GFN2 — for a slower, more careful
ranking of the same conformers.

Results load straight into the Trajectory tool to step through, or export as
a multi-frame XYZ. A "re-rank with xTB" button lets you sharpen the ordering
of just the best few conformers after a quick built-in search finds them.

For a transition state, the same tool can search over reactant conformers and
push each one over its barrier with an artificial force, collecting a set of
transition-state guesses that can have genuinely different active-bond
lengths — not the single length a conventional constrained search would be
stuck with.

### Measurements

Distances, angles and dihedrals by clicking atoms. Left-drag orbits the
camera, so a click is distinguished from a drag by how far the pointer
travelled. Per-measurement and global styling, labels drawn next to the line.

### Representation and Appearance

Six representations: ball-and-stick, rounded sticks, low-resolution balls and
lines, lines only, backbone trace, space-filling. Per-element visibility
filters, atom and bond scaling, bond radius as a fraction of the smaller
covalent radius, uniform or atom-split bond colours.

Colour schemes: CPK, Jmol, VMD, Molden, Molden0, or a custom scheme you name,
save, and optionally load at startup. Material controls (metallic, roughness,
reflectance, emissive), an ambient light and a key/fill/rim rig, background
colour.

Hydrogen bonds are detected on a distance cutoff and drawn as dashed lines
with their own colour, thickness and dash geometry.

### Also included

* **Geometry Optimization** — runs `xtb --opt`, or the built-in force field,
  in the background without blocking the viewer, and hands back the whole
  optimisation trajectory to play through.
* **Structure Comparison** — RMSD between frames or loaded structures, as
  placed and after optimal superposition. Pin any frame as the reference and
  see which atom moved furthest.
* **Diagnostics** — flags structures that are probably wrong before you spend
  time on them: close contacts, isolated atoms, and bond orders that sit
  between an element's common valences (a likely missing hydrogen), including
  a check for an odd number of electrons that means the molecule is a
  radical.
* **Export Image** — PNG or JPEG of the viewport, or of a selected area of
  it.

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
| Spectra, thermochemistry | — | ✅ `.dat` |
| Viewport | — | ✅ PNG, JPEG |

¹ A plain Gaussian `Freq` job does not write the force-constant matrix; that
needs `iop(7/33=1)` or a formatted checkpoint. Beavyr reads the frequencies and
modes Gaussian printed.

## xTB

Beavyr's own force field handles geometry cleanup and conformer searching with
nothing to install. For quantum-chemistry results — optimisation, frequencies,
a better conformer ranking — point Beavyr at an
[xTB](https://xtb-docs.readthedocs.io/) executable once and the path is
remembered. Beavyr also reads output that **ORCA** or **Gaussian** already
produced; those are jobs you run yourself, not something a viewer launches.

## License

GPL-3.0. See `LICENSE`.

## Building from source

```sh
cargo run                                       # dev build, dynamic linking, fast rebuilds
cargo build --release --no-default-features     # standalone binary, no dynamic linking
cargo test                                       # 951 passing, 17 ignored
```

The `dev` feature (on by default) enables Bevy's dynamic linking. Turn it off
for anything you intend to ship. Dependencies are built at `opt-level = 3`
even in a dev profile, because an unoptimised wgpu makes the viewer unusable.

Beavyr links no quantum-chemistry code of its own. Every external backend is
an executable you point it at, run in a scratch directory and parsed back
from its output files. That keeps the dependency list to five crates —
`bevy`, `bevy_egui`, `rfd`, `ndarray`, `anyhow` — with no BLAS and no Python.


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

- Simulated NMR spectra
- A UI for the vibronic analyses above
- File associations, so double-clicking an `.xyz` opens Beavyr
