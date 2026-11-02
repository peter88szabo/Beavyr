<p align="center">
  <img src="packaging/beavyr.svg" alt="" width="132" height="132">
</p>

<h1 align="center">Beavyr</h1>

Beavyr is a molecular viewer for chemists. Load a structure or build one from
scratch, watch a trajectory or a reaction play out, and see what a
quantum-chemistry calculation actually looked like — normal modes, IR and
UV-Vis spectra, molecular orbitals, conformers — all in one fast, interactive
window.

See the tools in action: [build a molecule](#molecule-editor),
[optimize its geometry](#geometry-optimization), [animate vibrations](#vibrations),
[view orbitals](#orbital-viewer), [inspect excited states](#uv-vis),
[find conformers](#conformational-analysis), or [measure a structure](#measurements).

**Author:** Peter Szabo
**Affiliations:** KU Leuven & BIRA-IASB, Belgium
**Join Discord:** https://discord.gg/bA32eVnCQ 

## Download

Get the current release from the
[Releases page](https://github.com/peter88szabo/Beavyr/releases):

| Platform | File | Install |
|---|---|---|
| Windows | `beavyr-setup.exe` | Run it |
| macOS | `beavyr-macos.dmg` | Open it, drag Beavyr to Applications |
| Debian, Ubuntu, Mint | `beavyr_amd64.deb` | `sudo apt install ./beavyr_amd64.deb` |
| Arch, Manjaro | build with `packaging/PKGBUILD` | `makepkg -si` |

No Python, no separate runtime, nothing else to install.

## What you can do

The interface is an icon rail down the side; each icon opens one window. Any
number can be open at once, and each keeps its own state.

### Molecule Editor

Load or paste XYZ coordinates through **Structure (XYZ)**, then open
**Molecule Editor** to edit the structure or build one by adding atoms and
fragments. Select an atom in the viewer to work on it; the Z-matrix table lets
you edit bond lengths, angles and dihedrals.

The editor opens as a floating, resizable window. Common elements
and fragments are one click away; search the fragment library by name or formula.
The selection card identifies the atom that **Add** or **Replace** will act on,
and Undo/Redo and the geometry tools stay visible while the library scrolls.
Drag it wherever you prefer, or choose **Dock** to attach it to the right edge.
Use **Float** to detach it again.

After attaching or replacing a fragment, the existing **Newly Placed Fragment**
controls are ready to rotate, stretch and bend it. **Fragment Editor ↗** opens
the existing tool with that fragment's atoms and attachment axis already set.
Its controls and Apply/Cancel behavior are unchanged.

**Peptide Builder…** opens a separate floating window with 22 amino-acid tiles,
sequence input, N- or C-terminal growth, protonation choices and a rotatable
preview. Build a separate peptide or extend a selected free amino/carboxyl group.
Insertion is undoable and prepares the existing Fragment Editor for the new
atoms. Its **Edit residues** tab adds previewed residue replacement, side-chain
fitting, ACE/NME caps, positioned disulfide links and PDB export. Residue chemistry
survives project save/load and Undo/Redo; **Add missing H** restores named peptide
hydrogens and now supports neutral amine/amide nitrogen in ordinary structures.
Peptides hide hydrogen bonds, and chains longer than 20 residues default to a
backbone trace. **Structure → Open structure…** and command-line file opening
accept PDB/mmCIF proteins, retaining residue identifiers and recorded covalent
links. The Peptide Builder's **Edit residues** tab also offers φ/ψ/ω edits,
range presets, a clickable Ramachandran plot, and clash-ranked Dunbrack rotamers
with nearby side-chain repacking. Every conformational edit previews before an
undoable Apply. See the [Peptide Builder guide](data/amino_acids/PEPTIDE_BUILDER.md).

Protein import uses `pdbtbx` 0.12; rotamers use `dunbrack` 0.1 and its embedded
2010 tables. Both are Rust dependencies; there is no Python, runtime download or
new runtime configuration. The tables increase the binary's data footprint and
make the first compilation slower. Parser parallelism is disabled; builds and
tests remain limited to four cores.

![Modern Molecular Editor with a selected atom, fragment palette, Add and Replace actions, and fixed geometry tools](Images/EditorModern.png)

The atom count and molecular formula in Hill notation (`C12H24O2`) sit at the
top of the modern panel; hover the formula for the molecular weight. This is
the cheapest way to notice that a structure is one hydrogen short before doing
anything else with it. An unknown element symbol suppresses the weight and is
named, rather than quietly producing a weight that is too small.

Load or paste XYZ. Edit the coordinates as text and re-apply, or edit the
molecule directly in 3D. The builder adds and removes atoms, attaches
fragments with chemistry-aware placement (VSEPR geometry and valency, not just
a bond vector), sets distances and angles numerically, rotates about a bond,
and undoes each of those. Sixteen fragments ship with it: `-CH3 -OH -NH2
-Phenyl -Pyrrole -OCH3 -COOH -NO2 -OOH -CCH -CycloPentane -CycloHexane -Cl -Br
-I -SH`.

A Z-matrix editor works alongside the Cartesian view — build in internal
coordinates, convert either way.

**Add missing H** fills saturated sp3 carbon to four bonds and neutral oxygen
to two, using the same chemically aware placement as **Add Atom (auto)**.
It keeps existing coordinates fixed and skips centres with detected multiple
or partial bonds. **Undo H addition** restores the structure before the addition;
use **Save to history** to keep an edited version between sessions.

The Molecular Editor's **Select** menu selects by element, connected fragment,
bonded neighbours, or the last placed fragment. Click replaces the selection;
**Shift-click** toggles atoms. The menu also offers **All visible**, **Invert**,
and **Clear selection**. The original dashed magenta rings show selected atoms. The active atom remains
the attachment target for the existing placement controls.

**Freeze** fixes selected atoms for quick cleanup; amber boxes mark them.
The **Frozen** menu unfreezes selected atoms or all atoms. **Isolate** hides
unselected atoms, bonds and labels; **Show all** restores them. Isolation affects
the view and picking only: calculations still use the complete molecule.
Selection, isolation and freezing are session tools and are reset when a new
structure is loaded. Appending fragments preserves freezes on existing atoms.
These tools do not change the Fragment Editor's ordered picks or transformations.

**Clean up geometry** opens a floating DREIDING cleanup window. Choose **Whole
structure** or **Selected atoms**, then **Preview cleanup**. Frozen atoms always
stay fixed; in selected mode all unselected atoms stay fixed too. The built-in
force field runs in the background and can be cancelled. Cyan outlines show the
proposed geometry and displacement without changing the current coordinates.
**Accept cleanup** records one Undo/Redo step; **Discard** or closing the window
leaves the molecule unchanged. Changes to the source geometry or movable atoms
invalidate the preview. Dummy atoms are retained, fixed, outside the force field.
Cleanup reports energy change and convergence, and rejects unsupported atom
types. It prepares a starting structure for optimisation; metal coordination,
unusual bonding and electronic states require a suitable quantum method.

The atom context menu also offers **Delete atom** beside **Use as rotation
center**. Deletion keeps the remaining Cartesian coordinates fixed and recreates
the Z-matrix, including when later rows referenced the deleted atom.
Clearing, replacing, or editing the geometry also removes orbital/density
surfaces that belong to the previous structure.

The **Structure (XYZ)** panel keeps file actions in a compact toolbar and gives
the remaining window space to coordinates. The atom count and formula share
one line; hover over it for molecular weight and the file path. **View** contains
atom labels, hydrogen bonds and centering. Editing keeps **Apply coordinates**,
**Reset**, and **Clear** visible below the scrollable draft.

![Compact Structure panel with a large coordinate area](Images/StructureModern.png)

Right-click the background for the whole-molecule tools: recentre the view,
choose **Representation**, toggle hydrogen bonds, or open **Orientation** for
**Orient to plane**, **Mirror across plane**, and **Flip around axis**. Orient
and Flip keep a chiral molecule as it was; Mirror gives you its enantiomer.

### Projects and custom fragments

Use **Structure (XYZ) → Project → Save project…** to write a `.beavyr` file.
It embeds the structure, trajectory, calculation settings and manual input,
Molden orbitals, frequency results and their Hessian or imported modes, UV–Vis
results, spectrum settings, camera, appearance and window layout. Original data
files are not needed to reopen it. **Open project…** validates the saved data
before replacing the session; **Restore previous session** reverses that open.
Calculations must be stopped before opening, and restored playback starts paused.
Executable paths stay local to the installation. Calculation log files, scratch
artifacts, conformer search results and measurement annotations are not part of
this project format; save those separately using their tools.

In Molecular Editor, open **My fragments → Save molecule as fragment…**.
A **name is required**: it becomes the button used to select the fragment later.
Choose **Whole molecule** or **Selected atoms**, then the **Attachment atom**.
Atom selections accept ranges such as `1-6, 9`; **Add picked atom** uses the
viewport selection, and **Use edited fragment** uses the current Fragment
Editor selection. Select a saved button, then **Place**, **Add** or **Replace**
through the existing placement controls. Attached fragments use the existing
Fragment Editor, including translation, rotation and dihedral adjustment, and
placement supports Undo/Redo. Right-click a saved button to remove it; **Undo
removal** restores the most recently removed entry.

The library is saved per user under `~/.config/beavyr/fragments` on Linux
(`$XDG_CONFIG_HOME/beavyr/fragments` when set), `%APPDATA%/Beavyr/fragments` on
Windows, or `~/Library/Application Support/Beavyr/fragments` on macOS. Names
must be unique within the library. Projects and fragments use versioned RON
serialization through `serde`/`ron`; egui persistence stores the window layout.

### Quantum Chemistry

**Tolerances…** opens a separate, resizable editor for minimum optimization,
Opt+Freq, transition-state optimization, MECP, and relaxed scans. **Save**
applies the draft values to the calculation settings; **Cancel** leaves the
previous values intact. **Reset defaults** resets the draft, then **Save**
applies it. Settings are retained in saved projects. Native optimizers keep
their own defaults until custom criteria are enabled; native xTB offers its
convergence-level presets. Beavyr MECP uses Å and Eh/Å for displacement and
gradient criteria, as in SpinCross. Other units are shown beside each field.

Behemoth's native minimum, TS, and relaxed-scan optimizers accept the optional
CLI flags `--opt-energy-tol`, `--opt-max-step`, `--opt-rms-step`,
`--opt-max-gradient`, and `--opt-rms-gradient`. Corresponding underscore names
are accepted in its `@opt` input section. These are convergence criteria,
not step-size limits; unchanged settings preserve Behemoth's existing defaults.

Open the atom-shaped launcher beside the Molecule Editor launcher. The modern
calculator first opens floating on the left; **Dock** attaches it to the right edge and **Float**
detaches it again. The **Quantum chemistry code** selector stays at the top,
with **Configure…** beside it. Its status tells
you whether the executable or Python environment is available.

For **sTDA-xTB (Grimme)**, enter both the standalone **stda** executable path
and the **xtb4stda** executable path in **Configure…**. Both paths are remembered;
Behemoth is not used for this backend. It offers only **UV-Vis / excited states**,
with a closed-shell singlet reference and an energy window in eV. The installed
programs and their parameter files must be configured for your environment.
Beavyr follows the [xtb4stda workflow](https://github.com/grimme-lab/xtb4stda)
and [sTDA manual, section 7.2](https://www.chemie.uni-bonn.de/grimme/de/software/stda/stda_manual.pdf):
run xtb4stda on the geometry, then run stda on the resulting `wfn.xtb`.

For **Behemoth**, choose **More… → UV-Vis / excited states** and select **sTDA-xTB (GFN1 orbitals)**,
**sTDA (TASI)**, or **sTDA / sTD-DFT (DFT reference)**. The DFT route offers
hybrid functionals and a basis; Tamm-Dancoff selects sTDA, and switching it off
selects sTD-DFT. Completed spectra open in the UV-Vis panel. Calculation failures
use the existing error popup with its Copy button.

**CCSD(T)** is available for ORCA, Psi4, PySCF and Behemoth. **CASSCF** and
**NEVPT2** are available for ORCA, PySCF and Behemoth, with explicit **Active
electrons** and **Orbitals** settings. The newly added methods support **Energy**
calculations and rigid scans: choosing one selects Energy unless Scan is active,
and disables geometry and spectrum jobs. Existing ORCA and Psi4 CCSD(T) job choices are unchanged. NEVPT2 reports
the corrected total energy; each backend retains its own frozen-core and
orbital-selection defaults.

Choose **More… → Scan** for a one-dimensional bond, angle, or dihedral scan.
**Scan settings** contains rigid/relaxed selection, ordered atom numbers (starting
at 1), viewport atom picking, start/end values, and the number of points including
both endpoints. Bond distances use Å; angles and dihedrals use degrees. Rigid scans
can move the connected fragment or only the last selected atom. A ring cannot be
split into a moving fragment; choose atom-only movement or a relaxed scan instead.
Relaxed scans allow additional fixed bonds, angles and dihedrals, a maximum step
limit, and a choice of primitive or Pulay redundant internal coordinates with
the Beavyr and Behemoth drivers. Rigid bond scans can move the first atom,
second atom, or both symmetrically; connected fragments follow that selection
when fragment movement is enabled. The same choices apply to rigid ORCA scans.

The Beavyr and Behemoth drivers can stop at an unconverged point.

ORCA offers **ORCA built-in** and **Beavyr** drivers. The built-in route uses
ORCA's relaxed `%geom Scan` or a multiple-XYZ energy scan for rigid geometries.
Behemoth also offers a **Behemoth built-in** driver, invoking its existing
`--scan` CLI with the same coordinate, atoms, start, end, points and mode settings.
Its native frontend supports HF/DFT/MP2/RI-MP2, rigid and relaxed bond/angle/dihedral
scans, fragment movement, extra fixed coordinates, maximum optimization steps,
and stopping at unconverged points. Rebuild Behemoth with the companion scan
CLI/input extension before using these new controls. They map to
`--scan-max-cycles`, `--scan-move-fragment`, `--scan-bond-motion`,
`--opt-coord`, `--scan-stop-on-failure`, and `--scan-constraints`; equivalent
`@scan` keys are `max_cycles`, `move_fragment`, `bond_motion`, `coordinates`,
`stop_on_failure`, and `constraints`. A constraint list uses
`bond:1,2:0.96;angle:2,1,3:104.5`.
Other backends use Beavyr's scan driver and the selected program for energies
and gradients. DREIDING retains its element/topology limitations; sTDA-xTB is
spectrum-only. Methods currently restricted to energies can be used in rigid
scans, but not relaxed scans. Scan inputs are generated from these settings;
the ordinary single-job input editor is not used. Dummy atoms are removed and
real atom numbers are remapped; a scan cannot refer to a dummy atom.

Completed scans automatically open separate **Scan Results** and **Scan Plot**
windows. **Show Plot** reopens a closed plot; **Open scan results** reopens the
results table. Energies are relative to the **First point** by default. Select
**Last point** to change the reference, and toggle **kcal/mol** (default) or
**kJ/mol**. Hover over a plotted point to see its coordinate and relative energy.
**Reverse order** reverses the numerical table, plot direction, clipboard/export
data, and loaded trajectory together; first/last reference follows that order.
The original calculated energies and geometries are retained.

**Copy data** copies the table as CSV, with absolute Hartree energies as well as
the selected relative energies, units, and reference point. **Save CSV…** and
**Save XYZ…** export the displayed order. Select a point to display its geometry
without resetting the camera, or choose **Load as Trajectory** to load the
displayed order and open the Trajectory panel. Convergence is
reported per point. Files `scan.xyz`, `scan.csv`, and `scan_results.ron` are also
saved in the calculation directory in the original calculation order.
Cancellation preserves files already written and does not apply a result to the
viewer. Only one scan coordinate is swept at a time.

Choose **More… → Transition-state optimization** for a saddle search.
ORCA offers its native OptTS driver or **Beavyr**; PySCF and xTB use Beavyr's
existing Behemoth-derived optimizer. Psi4 retains its native driver and also
offers Beavyr. Behemoth invokes its own internal TS optimizer; its current
CLI supports GFN1-xTB for this job.

Beavyr exposes Cartesian, primitive internal and Pulay redundant internal
coordinates, reaction-mode references, additional coordinates and advanced
trust-radius/Hessian controls. Defaults follow Behemoth's CLI: primitive
internal coordinates, calculated initial Hessian, internal mode tracking,
100 steps and no periodic Hessian rebuild. Both native ORCA and Beavyr offer
**Recalculate Hessian every N steps**; zero retains the initial Hessian without
periodic rebuilds. ORCA starts at its previous interval of five steps.
A local TS is reported successful only after convergence and a final
calculated Hessian with exactly one negative mode. Disabling the final
Hessian leaves the stationary structure unverified.

For **Optimize** and **Opt+Freq**, choose **Beavyr** to use the existing minimum
optimizers with PySCF, ORCA, Psi4 or xTB. The PySCF route needs PySCF energy
and gradient support, but does not require geomeTRIC. Cartesian BFGS,
steepest descent and both two-point gradient methods are available; internal
coordinates offer BFGS and steepest descent. The default is primitive
internal BFGS. Native optimizers remain selectable; Behemoth keeps its own.
Failed optimization skips follow-up frequencies and is reported explicitly.
Search trajectories, energies and summaries are saveable; requested orbital
export runs at the final converged geometry.

Choose **More… → MECP / Spin crossing** for a minimum-energy crossing search
between two spin multiplicities. Set **State A multiplicity**, **State B
multiplicity**, and the common charge, HF/DFT method, and basis. Both states
are evaluated independently at the same geometry. This workflow supports
ORCA, Behemoth, PySCF, and Psi4; spin-independent force fields are excluded.

For ORCA, select **ORCA built-in** or **Beavyr MECP**. The native driver uses
`%mecp SurfCrossOpt true` and `Mult` for State B; the coordinate block specifies
State A. It uses ORCA's optimization criteria and checks the final energy gap
against **Final gap limit (Eh)**. The Beavyr driver exposes Cartesian,
primitive internal, and Pulay redundant internal coordinates, a step limit,
and an energy-gap tolerance. Both drivers report nonconvergence explicitly.
The default Beavyr search uses Cartesian BFGS in the original SpinCross
Angstrom convention: the initial inverse Hessian is diagonal with **0.7**
on every diagonal and zero off-diagonal entries. The initial previous
effective gradient is zero. Backend geometries and gradients are converted
at the optimizer boundary; the experimental initialization is unchanged.

MECP completion opens the calculation summary and loads the trajectory into
the Trajectory panel. The summary shows both state energies, their difference,
and the two multiplicities, with a **Copy MECP results** button. CSV, XYZ,
and summary files can be saved; trajectory energy plots use the mean of the
two state energies. A converged crossing is not a substitute for checking
curvature along the crossing seam; an adapted MECP frequency calculation is
not included in this job.

For PySCF, Psi4 (Python), and Sparrow (Python), open **Configure** and enter an
**Environment activation** command, then click **Save & Check**. Settings are
remembered separately for each backend and apply without restarting Beavyr.
On Linux/macOS, examples are `source /any/path/to/venv/bin/activate`,
`conda activate myenv`, or `micromamba activate psi4env`. Environment paths may
be anywhere; quote paths containing spaces. Conda and micromamba shell hooks
must already be initialized in the selected shell, as for normal terminal use.

Under **Python interpreter and activation shell**, you can instead browse to
the environment's Python executable, leaving activation empty. Beavyr detects
the OS and defaults to Bash on Linux, zsh on macOS, and PowerShell on Windows;
an existing Bash/zsh login-shell choice is respected. Windows venv activation
uses `& 'C:\path\to\venv\Scripts\Activate.ps1'`, or select `python.exe`
directly. Activation runs in the calculation's child process and does not
change Beavyr's own environment. Failed activation prevents the calculation.

Choose **Energy**, **Optimize**, **Frequencies**, or **Optimize + frequencies**;
**More…** offers the other supported calculations. Search the program's method,
functional and basis lists by name or keyword. Charge, multiplicity, CPU cores
and memory stay visible above the scrolling settings. Additional options
expand when needed; manual input opens in a separate editor. The action and status stay at the
bottom while the settings scroll, and resources are limited to four cores.

The **Classic** option restores the previous calculator form. Both layouts use
the same code selection, settings, custom input and running calculation. For
Classic from startup, launch with `BEAVYR_CLASSIC_QC=1`.

**Edit input…** opens a large, resizable editor with **Save input…**, **Run input**,
**Sync from settings**, and **Close editor**. Closing keeps the draft active.
Sync replaces it with input generated from the current controls and molecule;
**Use settings** returns to generated input while retaining the manual draft.

![Modern Quantum Chemistry calculator with a prominent code selector and fixed calculation action](Images/QuantumChemistryModern.png)

The main application title includes the version declared in `Cargo.toml`,
for example **Beavyr 0.2.1**.

The toolbars use consistent vector icons with hover labels and blue highlights
for open tools. They stay sharp at different display scales and do not depend
on emoji fonts.

### Geometry Optimization

Relax the displayed structure with xTB or the built-in DREIDING
force field. Optimization runs in the background so the viewer stays responsive.

1. Load or draw a molecule, open **Quantum Chemistry** from the atom-shaped
   launcher, and choose **Optimize**.
2. Choose the code and method. For an external program, use **Configure…** and
   **Browse…** to select its executable. Check the structure's **Charge** and
   **Multiplicity** before starting, and set **CPU cores** and **Memory (MB)**
   directly in the setup header.
3. Click **Start optimization**. After completion, open **Trajectory** to inspect the
   optimization steps, or **Show Energy Plot** to see the energy history when
   available. Use **Save to history** in Molecule Editor to retain the result.

![Geometry Optimization with xTB settings, a completed optimization, and energy plotted against iteration](Images/Optimize.png)

### Recent Structures

The **Recent Structures** icon keeps the last 20 distinct structures, newest
first. Loading another structure or pressing **Clear Display** preserves the
structure you are leaving, including unsaved drawings and edits. The newly
loaded standalone structure is saved too. Click **Load** beside an entry to
put its coordinates back on screen.

Search by name or formula. Each saved version has a compact card showing its
name, composition and age, with **Load** alongside it.

Editing and trajectory playback do not create history entries. Use **Save to
history** at the top of **Molecule Editor** to keep a particular edited or
optimized version; an optional name helps distinguish versions. A paused
trajectory frame can also be saved explicitly with this button.

Snapshots are written immediately as XYZ coordinates in the configuration
directory, so they survive restarts and remain available in another Beavyr
window. Windows merge their saved entries and retain the newest 20; reopening
a structure moves it to the top without adding a duplicate. On Linux the files
are in `$XDG_CONFIG_HOME/beavyr/structure_history` or
`~/.config/beavyr/structure_history`; on Windows they use
`%APPDATA%/Beavyr/structure_history`, and on macOS
`~/Library/Application Support/Beavyr/structure_history`.

### Trajectory

Plays multi-frame XYZ — MD, optimisation paths, IRC, a set of conformers — 
forwards or backwards, looping, at a chosen frame rate, with a frame slider.

The panel separates **Trajectory file**, **Playback**, and **Frame overlays**.
Use **First**, **Previous**, **Next**, and **Last** to navigate, or **Play** and
**Play backward** for motion. Playback modes and speed remain beside the controls.

Frames can also be overlaid on the current structure rather than played: a
solid overlay and a transparent "ghost" overlay, each with its own stride, so
a 100-frame optimisation can be shown as every tenth frame ghosted behind the
final geometry.

### Animation export and trajectory analysis

**Export animation…**, in the Structure and Trajectory panels, saves the loaded
trajectory/vibration or a turntable rotation of the current molecule. Choose
GIF or MP4, image dimensions, frame rate and frame range/stride. GIF export is
built in; its LZW encoder uses the small `weezl` dependency (also already used
indirectly by the image stack). MP4 requires an FFmpeg executable with the
`libx264` encoder: use FFmpeg from PATH or select its executable in the export
window. No shell activation or FFmpeg installation is performed by Beavyr.
MP4 dimensions must be even. Export can be cancelled and incomplete output is
discarded. The main molecule and camera are not changed during capture.

Exports render molecular atoms/bonds and the chosen background. They do not
include UI panels, orbital surfaces, measurement overlays or custom labels.
Line bonds are approximated with thin cylinders. GIF is limited to 256 colours
per frame and centisecond timing; MP4 supports smoother colour gradients.

**Trajectory analysis…** in Trajectory plots energy, distance, angle, dihedral or aligned
RMSD, and can show distributions or copy/save CSV. Atom indices are 1-based.
Click a plotted point to display that frame without resetting the camera.
Missing energies leave gaps instead of shifting frame indices. Energies are
read from explicit XYZ comment fields/units or retained calculation results;
unlabelled numbers are not guessed to be energies. Relative energies can use
the first, minimum or last available value as reference. An existing QC energy
plot also selects a frame when its energy sequence matches the loaded path.

### Dataset inspection, saved views and comparison

**Inspect dataset…**, in Trajectory and PES Sampling results, opens a separate resizable
window. It detects short contacts using the Pebble H–H/H–other/other–other
cutoffs and aligned-RMSD duplicates, with editable thresholds. Atom identity
and order must match for RMSD comparisons; atom permutations are not searched.
Selection is editable/undoable, and a diverse subset can be selected by greedy
farthest-point sampling. Load a selection as a trajectory or save it to a new
folder containing `selected.xyz` and `selection.csv`. Original XYZ coordinate
text is retained when available; the manifest records source frames and
exclusion reasons. Filtering does not modify the source file.

**Saved views / notes…** in Structure stores named camera orientations, notes and atom
labels. Save a `.beavyr` project to keep these between sessions. Labels follow
the ordered atoms of matching frames. **Side by side…** in the RMSD panel opens
two independent molecular views: capture the current molecule, load a file or
use the current trajectory in either pane. Rotate, zoom, pan, fit or align B to
A, optionally synchronizing rotation/zoom. Alignment can be undone and never
modifies the molecule in the main editor.

### Vibrations

Animate any normal mode at a chosen amplitude and frame rate — the fastest way
to see whether a calculated frequency is the stretch, bend, or torsion you
expected. Run a frequency calculation yourself (through xTB, or the built-in
force field for a quick check) or load a Hessian someone else computed.

* **Normal modes** — mass-weighted Hessian, Eckart projection to remove
  translation and rotation.
* **Reaction-path projection** — supply a gradient and the mode along the
  reaction path is projected out too, which is what you want along a MEP 
  (reaction path). This can be handy if you do variational TST.
* **IR spectrum** — intensities from the program's own dipole derivatives,
  broadened as Gaussian, Lorentzian, Voigt or plain sticks, with adjustable
  width and a frequency scaling factor. Exportable as `.dat`.
* **Thermochemistry** — ZPE, thermal corrections, entropy and free energy at a
  chosen temperature, using Grimme's quasi-RRHO treatment of low-lying modes
  with a settable cutoff. The rotational symmetry number must be stated rather
  than guessed; a point-group reference table is built into the panel.

Use **Load Hessian…** to open a supported frequency file, then click
**Animate** beside a mode to see the corresponding molecular motion.
**Analysis settings** starts collapsed; expand it to adjust analysis and
thermochemistry. The mode table highlights the playing mode, and
**Calculate frequencies** opens the calculation settings.

![Frequency Analysis showing normal-mode frequencies and animation of an imaginary mode](Images/Vibrations.png)

### Orbital Viewer

Loads a `.molden` file, evaluates the Gaussian basis on a grid, and meshes an
isosurface by marching cubes. Three quantities: a single molecular orbital,
the total electron density, or the spin density. Adjustable isovalue, grid
resolution, opacity, per-lobe colours, and a wireframe overlay. The sampled
field is retained, so dragging the isovalue slider only re-meshes rather than
re-evaluating the basis from scratch.

Open **Surface Tools**, click **Load Molden…**, and choose **Orbitals**,
**Electron density**, or **Spin density**. For an orbital, select a row in the
list and, for an open-shell calculation, choose **Alpha** or **Beta** spin.
Type an orbital number in **#** and press **Enter** to reveal it. Use **↑/↓**
(or the previous/next buttons) to move through the visible list; **Page Up/Down**
moves ten entries. **Filter** offers occupied/virtual orbitals and an inclusive
energy range in Hartree. Filtering preserves the original orbital numbers;
a direct index jump clears filters so the requested orbital is visible.
Adjust **Isovalue** and **Surface resolution**, then choose **Solid** or **Mesh**.
The orbital list fills the remaining window height, with **HOMO** and **LUMO**
marked on their rows. The default orbital isovalue is **0.07**; isovalue,
visibility and opacity stay above the list. **Surface detail and colors** opens
the resolution, solid/mesh style and lobe color controls.

![Modern Surface Tools panel with isovalue, visibility and opacity above the orbital browser](Images/SurfaceToolsModern.png)

*Solid view with the selected orbital and separate colours for its two phases.*

![Surface Tools displaying a selected beta orbital as a solid blue and yellow isosurface](Images/OrbitalTools.png)

*Mesh view; use the mesh grid and line-width controls to adjust its appearance.*

![Surface Tools displaying an orbital as a blue and yellow mesh](Images/OrbitalTools2.png)

### UV-Vis

Reads excited states — energies and oscillator strengths — from an ORCA
TD-DFT output, and plots the absorption spectrum against wavelength or energy
with the same four line shapes UV-Vis needs, over a settable energy window.
The wavelength plot starts with a **5.0 nm** width. Export the curve for a figure.

Click **Load output…** to load excited states, expand a state to inspect
its orbital contributions, or click **Show Absorption Spectrum** to plot them.
Use **Matching Molden…** to pair the spectrum with its orbital file, or
**Link current orbitals** when Surface Tools already holds that calculation's
orbitals. Click either orbital number in an expanded contribution to open it
in Surface Tools. The display retains the output's numbering and maps it to
the appropriate file index, spin and (for Psi4) symmetry. Loading another
wavefunction invalidates the link. Runs that supply both a spectrum and
Molden orbitals link them automatically. Use the original canonical orbitals
from the same calculation.
**Calculate spectrum** opens the calculation controls; it stays open during a
run so the cancellation control remains available.

![UV-Vis tool showing excited-state energies and expanded orbital contributions alongside an orbital surface](Images/ExcitedTool.png)

### Transition-State Initial Guess Generation

The **TS Generation** button on the left icon rail opens RDA and Poor Man's
NEB. Both use the Rust optimizers with energy and
gradient calculations supplied by an installed xTB executable.

![TS initial guess generation](Images/TSgenerator.png)

**How To Use** opens a guide with atom-number examples for every Reaction
definition, active atoms, alignment, coordinate weights, and the complete
workflow. Detailed RDA option explanations are in a section that starts collapsed.

1. Use **Load A/B XYZ** to load endpoint files, or draw/optimize the molecule
   in the viewer and press **Use current as A** (reactant). Then draw/optimize
   the second structure and press **Use current as B** (product). Each stored
   endpoint has a green loaded marker and a **Show A/B** button to display it
   again. Captures stay unchanged while you edit or optimize the viewer; capture
   again to keep those changes. Both endpoints must have the same atoms in the
   same order.
2. Define the reaction using one-based atom numbers: a bond change (two atoms),
   atom transfer (donor, transferred atom, acceptor), an angle, a dihedral, or
   custom reactive bonds/angles/dihedrals. RDA also supports Cartesian motion
   on selected active atoms. Active atoms select the distance/alignment subset;
   they are not frozen-atom constraints.
3. Select RDA or Poor Man's NEB and configure the backend, method, charge and
   multiplicity. RDA exposes beta bracketing, gamma search, distance and energy
   tolerances, coordinate mode and step limits in a section that starts collapsed.
   Poor Man's NEB exposes image
   count, constrained relaxation limits and additional connectivity.
   Enable **Relax End-point Geoms** to optimize both endpoints with the same
   backend and level of theory before the search. This is off by default;
   unconverged endpoint optimizations stop the run with an error.
4. Generate, then use **Show TS guess** or **Show path in Trajectory** in the
   separate **TS Generation Results** window, which opens automatically when
   either method finishes. **Open TS results** in the setup panel reopens it.
   Calculations run in the background and can be cancelled, including during
   endpoint relaxation.

**Save TS guess XYZ** writes one structure in Å. **Save NEB trajectory XYZ**
writes the complete multi-frame XYZ path, including both endpoints and all
constrained-relaxed interior images; each frame's comment includes its energy,
convergence status and whether it is the selected TS guess. The highest-energy
interior image is selected. **Plot NEB energy** opens a resizable energy plot
of all images, including A and B, relative to A in kJ/mol. The selected TS is
marked; hovering over an image shows its absolute and relative energies.
The NEB energy profile and RDA search structures
also have separate save buttons. These outputs are TS guesses; saddle-point
optimization and a frequency check are needed to establish a transition state.

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

Open **Conformer Search**, choose the energy method and search effort, then
click **Search for conformers**. The energy model and search effort have
separate sections, with three effort presets and visible core controls.
Use **View** beside a result to inspect it,
**Load as trajectory** to browse the set, or **Save as XYZ trajectory** to export it.

![Conformer Search showing search settings, relative energies, populations, and controls to display or export conformers](Images/Confsearch.png)

### PES Sampling

The **PES Sampling** button on the left toolbar generates a sampled XYZ
trajectory. It follows Pebble's unimolecular/bimolecular workflows and Smite's
normal-mode sampling conventions. Completed samples load automatically and
open the **Trajectory** panel, stopped on the first frame. **Save XYZ…** exports
all samples; **Load as Trajectory** reopens the retained result.

For each molecule or collision partner, load an ORCA `.hess` containing its
geometry and masses, or load the matching equilibrium XYZ followed by a bare
`$hessian` matrix in Eh/bohr². Diatomics can instead use their two elements,
equilibrium distance in Å and frequency in cm⁻¹. An optional Morse oscillator
at J = 0 additionally needs De in eV. Atomic partners and rigid partners are
also supported. Imaginary modes are excluded, with an explicit notice, as in
Smite; translations and rotations are projected out. Hessian isotope masses
are retained.

Every mode uses an independent phase uniformly distributed over [0, 2π).
With a specified energy ceiling, `Q = sqrt(2 Emax) / omega * cos(phase)`;
Cartesian displacements use mass-normalized mode vectors. Choose ZPE,
maximum quantum number, maximum energy **per mode**, or temperature. The
common setting applies to every mode; **Mode settings → Set limits per mode**
overrides individual modes. Temperature choices are:

- **Boltzmann-tail cutoff:** the smallest maximum quantum number retaining
  the selected canonical population (99.9% by default).
- **ZPE + kB × Tmax:** a continuous energy ceiling above the zero-point energy.
- **Smite thermal (no hard cutoff):** draw quantum numbers from the canonical
  distribution at the specified temperature for each sample. Morse sampling
  uses its physically bound states.

For bimolecular sampling, set first/last **COM distances**, distance-point count
and samples per distance. The distance grid includes both endpoints (one point
uses the first distance). Both partners receive independent phases and
isotropic rotations. Re-centering preserves the exact chosen COM distance,
and the combined center of mass stays at the origin. XYZ comments record the
distance, distance index, sample index, partner boundary and random seed.

**Repeatable** fixes the seed for identical inputs/settings. The overlap filter
is enabled by default for every sample, using Pebble's minimum distances:
H–H 0.4 Å, H–other 0.5 Å, and other–other 0.8 Å. It checks all atom pairs,
including pairs between collision partners. The cutoffs are editable, and the
filter can be disabled. Rejected samples are redrawn, changing the accepted
distribution. Generation runs in the background and
can be cancelled. Previous results survive cancellation or failure. Preview
is limited to five million atom-frames to bound memory use. This module
exports geometry samples only; it does not assign quantum-chemistry PES
energies or propagate molecular dynamics. Save the XYZ explicitly; sampling
setup/results are not stored as additional project fields.

### Measurements

Distances, angles and dihedrals by clicking atoms. Left-drag orbits the
camera, so a click is distinguished from a drag by how far the pointer
travelled. Per-measurement and global styling, labels drawn next to the line.

Open **Measurements** and choose **New distance**, **New angle**, or
**New dihedral**, then click two, three, or four atoms respectively. Use
**Hide** or **Delete** beside a measurement to manage the labels on screen.
Each measurement type has its own section, with its picking action visible.

![Measurements panel with atom-picking controls and labelled distances on a molecular structure](Images/Measurment.png)

### Representation and Appearance

Six representations: ball-and-stick, rounded sticks, low-resolution balls and
lines, lines only, backbone trace, space-filling. Per-element visibility
filters, atom and bond scaling, bond radius as a fraction of the smaller
covalent radius, uniform or atom-split bond colours.

The Representation panel presents the six styles as a two-column chooser,
followed by the size and detail controls for the selected style. Appearance
opens **Colors and materials** first; **Atom and bond size**, **Lighting**, and
**Element visibility** have separate expandable sections.

Colour schemes: CPK, Jmol, VMD, Molden, Molden0, or a custom scheme you name,
save, and optionally load at startup. Material controls (metallic, roughness,
reflectance, emissive), an ambient light and a key/fill/rim rig, background
colour.

Hydrogen bonds are detected on a distance cutoff and drawn as dashed lines
with their own colour, thickness and dash geometry.

### Also included

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

The Export Image panel groups **File format**, **Resolution**, and **Framing**.
Choose PNG or JPG, set the pixel dimensions and print DPI, then **Select area**
and drag a rectangle in the 3D view. **Export image** becomes available once
the region is selected; **Redraw area** and **Clear area** adjust the framing.

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

### Ubuntu dependencies

Bevy 0.19 needs rustc 1.95 or newer, which is newer than Ubuntu's own `rustc`
package. Install Rust through rustup and remove the apt version, which would
otherwise sit ahead of rustup's in `PATH`:

```sh
sudo snap install rustup --classic
rustup default stable
sudo apt remove rustc cargo        # only if installed
```

Then install the system libraries Bevy builds against:

```sh
sudo apt install build-essential pkg-config libudev-dev \
                 libx11-dev libxkbcommon-dev libwayland-dev
```

`UBUNTU_BUILD_DEPENDENCIES.txt` explains what each package is for. Only
building from source needs these packages; the `.deb` does not.

### Building

```sh
cargo run                                       # dev build, dynamic linking, fast rebuilds
cargo build --release --no-default-features     # standalone binary, no dynamic linking
cargo test -- --test-threads=1                   # run the tests serially
```

Cargo is limited to four build jobs by `.cargo/config.toml`. For a hard limit
of four CPU cores on Linux, prefix build, test and run commands with
`taskset -c 0-3` (or four CPU IDs allowed on your machine). Launched calculations
inherit that CPU limit.

The `dev` feature (on by default) enables Bevy's dynamic linking. Turn it off
for anything you intend to ship. Dependencies are built at `opt-level = 3`
even in a dev profile, because an unoptimised wgpu makes the viewer unusable.

Bevy is built with its 2D, 3D and UI features, with audio disabled at compile
time. Linux builds and packages therefore do not require ALSA.

Beavyr links no quantum-chemistry code of its own. Every external backend is
an executable you point it at, run in a scratch directory and parsed back
from its output files. That keeps the dependency list to five crates —
`bevy`, `bevy_egui`, `rfd`, `ndarray`, `anyhow` — with no BLAS and no Python.


## Vibronic analysis

Open **UV-Vis → Vibronic analysis → Vibronic spectra**. The window resizes
horizontally and vertically. Select matching ground- and excited-state Hessians
(with their geometries and identical atom order) and enter the electronic
minimum-to-minimum gap in eV, without zero-point energy.

- Absorption and fluorescence use both Hessians, Duschinsky mixing, frequency
  changes and the initial state's thermal population. Internal spectra are
  relative; they do not report absolute fluorescence lifetimes.
- Resonance Raman accepts comma-separated laser wavelengths in nm and includes
  fundamentals, or fundamentals plus overtones/combinations. It currently uses
  one intermediate electronic state at 0 K. The intermediate-state linewidths
  and the displayed Raman linewidth are separate settings.
- Enable Herzberg–Teller for linear transition-dipole contributions, including
  interference with Franck–Condon amplitudes. For internal calculations, import
  a property file using **Save template**, or collect properties with PySCF or
  Psi4. Both use the existing saved Python interpreter/activation settings.
  The driver matches states by wavefunction overlap and aligns their phase
  before central differences; ambiguous matches stop with a copyable error.
- **ORCA ESD** invokes the user-configured ORCA executable directly. ORCA performs
  its native absorption, fluorescence, Raman and optional Herzberg–Teller work.
  Native ESD requires ORCA-format Hessians. Calculation directories are retained.
- Results support ground/excited mode animation, two-column experimental data
  overlays, explicit shifts, selectable axes, CSV copy/save and SVG figures.
  Axis changes preserve sampled intensity values, without a density Jacobian.
  Display transformations leave the original calculated data unchanged.

The original gradient-based Franck–Condon/short-time Raman activity analysis and
Duschinsky rotation diagnostics remain available as separate choices. The full
spectral calculation requires stable harmonic minima; imaginary modes are rejected.
Other backends can supply Hessians and property files. Automatic transition-dipole
collection is currently available for HF/DFT references in PySCF and Psi4; native
ORCA handles its own properties. Behemoth has not been modified for this feature.

Validation and method details: `reports/vibronic_extensions_2026-11-01.txt`.

## Recovery, comparison and trajectory properties

**Project → Autosave / recovery…** controls periodic recovery snapshots (60 seconds
by default, three generations per session). Recovery files are separate from manual
projects. Inactive sessions can be restored without restarting interrupted jobs;
Restore previous session undoes the restore. On Linux snapshots live in
`$XDG_STATE_HOME/beavyr/recovery` or `~/.local/state/beavyr/recovery`.

The side-by-side comparison window supports linked trajectory frames, manual or
geometric atom mapping, and displacement colours after rigid alignment. Automatic
mapping preserves elements, but equivalent atoms still require scientific review.
It supports up to 1,000 atoms / 500 atoms of one element; manual mapping remains
available for larger structures. Neither comparison camera changes the main view.

Open **Surface Tools → Trajectory wavefunctions…** to compute one wavefunction per
geometry. GFN2-xTB is the default; GFN1-xTB and the supported HF/DFT Quantum Chemistry
setup routes are also available. Set charge/multiplicity, select the executable,
choose a frame range and output folder. Each frame and an incremental manifest are
stored on disk. Load `wavefunctions.beavyr-wfn` to reopen a series. Playback can wait
for the matching surface. Orbital index playback is the default; optional character
tracking follows the selected orbital using adjacent-frame cross-AO overlaps and
aligns its phase. Weak or ambiguous matches stop tracked playback. The xTB version tested exports no separate beta orbitals, so its Molden
files cannot provide spin-density surfaces; a spin-resolved backend/file is required.

**Surface Tools → Atomic polarizabilities…** provides the Stöhr CPA/TS analysis used
in Behemoth: hA = sum of the AO density diagonal on atom A, alphaA = alphaFree hA/ZA.
It is an atomic descriptor, not a molecular dipole-response polarizability. It uses
full atomic Z consistently with Behemoth; valence-only/ECP calculations therefore
have the same basis/core limitations. Halo volume is proportional to alphaA. Values
can be shown in Å³ or bohr³. Atom/group sums can be plotted against frame, time,
distance, angle or dihedral; time requires the interval per sampled frame. Plot points
select frames without changing zoom. CSV copy/save includes the selected atoms.

**Surface Tools → Reactivity / Fukui…** offers fixed-geometry N/N+1/N−1 calculations
and a separately labelled HOMO/LUMO approximation. The finite-difference route uses
the Quantum Chemistry setup (GFN1/GFN2-xTB or HF/DFT with ORCA, PySCF or Psi4), with
reference charge and all three multiplicities specified explicitly. Atom-condensed
values use Mulliken populations with analytic Gaussian AO overlaps. Available
quantities are f+, f−, f0, the dual descriptor, local softness, local philicity, and
global hardness, softness, chemical potential and electrophilicity. The convention
is eta=(I−A)/2, S=1/(2 eta), mu=−(I+A)/2, omega=mu²/(2 eta); energies use eV.
Nonpositive gaps disable reciprocal-gap descriptors. Local hardness is not assigned
an arbitrary formulation. Fractional occupations are rejected by the single-HOMO/
LUMO approximation; finite differences remain available. Degeneracies, basis/charge
partition dependence and unbound electron-added states require interpretation.

Reactivity displays support signed surfaces, signed halos, numerical labels,
atom/group plots, global-descriptor plots, CSV copy/save, and reopenable analysis
files. Analysis files retain paths to native wavefunctions: keep their calculation
folders. Projects/recovery also retain completed polarizability and reactivity data,
measurements and the last QC result; external native files remain at their paths.
All these analysis windows resize in both directions.

Animation export can include per-frame orbital/density surfaces, polarizability and
reactivity displays, measurements recalculated at every frame, custom labels and
frame/energy captions. A changing geometry requires a matching wavefunction for
every surface frame. Exports use a private render scene and preserve the editor's
geometry/camera. Solid surfaces retain opacity/colours; exported mesh surfaces use
one-pixel lines. GIF export is native; MP4 requires a user-supplied FFmpeg executable.
The `ab_glyph` dependency rasterizes egui's bundled font into exported frames without
requiring a system font installation.

## Roadmap

- Simulated NMR spectra
- File associations, so double-clicking an `.xyz` opens Beavyr

## Surface controls, legends and trajectory workflow

Surface Tools and its atomic-property windows use aligned compact control rows.
Optional method/tracking details start collapsed. Atomic halos have an optional
colour legend with a fixed range derived from all calculated frames; a custom
range can be set. Polarizabilities use a sequential palette and signed reactivity
uses a diverging palette. Radius/volume semantics are unchanged. Display-unit
switches convert legend labels without changing the physical range. Animation
exports use the same halo colours; the viewer legend itself is not overlaid on GIFs.

**Trajectory wavefunctions → Resume calculation…** opens an existing calculation
folder. New runs save `wavefunction-plan.ron`, `plan.checksum`, settings and an
incremental manifest. Resume verifies the input plan, completed geometries and
wavefunction checksums, skips verified frames, and retries unfinished frames in
new folders. The saved method/charge/spin/range are used; current panel settings
do not replace them. A file lock prevents two workers resuming the same folder.
Older runs without a saved plan cannot be resumed. Progress shows elapsed time,
completed/remaining frames, the active source frame and an estimated time remaining.

In **Trajectory wavefunctions**, select an orbital in Surface Tools and choose
**Track selected orbital**. Tracking works forward and backward from that anchor
using normalized real-orbital overlaps `C(A)^T S(A,B) C(B)` within the same spin set.
The default minimum absolute overlap is 0.5, with a 0.1 separation from the next
candidate. Both thresholds can be adjusted. Phase signs are aligned for playback
and orbital animation export; native files are not modified. Near-degenerate
mixtures have no unique individual orbital identity, so ambiguous matches remain
unassigned. Tracking uses stored orientations (no rigid-rotation alignment) and
supports at most 1,000 basis functions per frame. Selecting another orbital turns
tracking off so a new anchor can be chosen. Tracks are session-only and can be rebuilt.

Polarizability, reactivity, trajectory-series and scan-energy plots offer a closed
**Plot appearance / export** section: font size, colour, solid/dashed/no line,
circle/square/no markers, line width, marker size and image dimensions. SVG and PNG
exports share the plotted geometry, axis units and scan reference caption. Missing
values remain gaps, and reversed scans preserve source-point selection. Spectrum
plots and trajectory histograms retain their existing separate rendering controls.
