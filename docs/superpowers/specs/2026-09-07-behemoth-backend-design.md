# Behemoth as a Second Calculation Backend — Design

Date: 2026-09-07
Status: approved for implementation

## Goal

Beavyr currently runs one external program, xTB, and asks for one binary path.
The user now has a second — Behemoth, their own Rust quantum-chemistry code —
and wants to choose which program a calculation goes to, giving a binary path
per program the same way the xTB path is given today.

Behemoth is used through its command line, exactly as xTB is. No source-level
linking: Behemoth's build needs CMake, LibXC and a BLAS, none of which Beavyr
has or wants, and a panic inside a linked-in SCF would take the GUI with it.

## What Behemoth's CLI offers

Verified against `Behemoth/target/release/behemoth --help` and real runs on
water:

```
behemoth --xyz FILE --method xtb --charge N --multi N --opt
behemoth --xyz FILE --method xtb --charge N --multi N --hessian
behemoth --xyz FILE --method xtb --charge N --multi N --freq
```

`--method` defaults to xTB-GFN1 and also accepts `tasi|hf|uhf|rohf|mp2|ri-mp2|
rks|uks`; only `xtb` is wired up initially, since the others need `--basis`
and (for DFT) `--functional`, which is a separate piece of UI.

### `--opt`

Writes two files into the working directory and a summary to stdout:

| | |
|---|---|
| `optimized_geoemtry.xyz` | final geometry (the filename's spelling is Behemoth's, and must be matched exactly) |
| `optimization_trajectory.xyz` | multi-frame XYZ, each frame's comment line reading `cycle 1 E = -5.7686678546 Eh` |
| stdout | `converged : true`, `cycles : 3`, `E(final) : -5.7687749333 Eh`, `RMS(g) : …` |

The per-cycle energies live in the trajectory's comment lines rather than in
the log, which suits Beavyr: the energy-versus-iteration plot reads them from
the same file the trajectory player loads.

### `--hessian`

Prints the full `3N x 3N` matrix to stdout in blocks of eight columns, labelled
by atom and axis:

```
                     O1x           O1y           O1z           H2x  …
       O1x    0.01353746   -0.00000000    0.00000000   -0.00568132  …
       O1y   -0.00000000    0.59082867   -0.00000000    0.00000000  …
```

This is the route to use for frequencies. Feeding the matrix into Beavyr's own
`analyze` gives the Eckart projection, the normal modes needed for animation,
and thermochemistry at whatever temperature, cutoff and symmetry number the
panel is set to — none of which is available from a printed frequency list.

### `--freq`

Prints frequencies, reduced masses and IR intensities in km/mol, then
thermochemistry in the same layout Beavyr uses (Beavyr's `format_thermo` was
imported from Behemoth). Not used: it reports no normal modes, so it cannot
animate. Kept in mind as the reference for validating the Hessian route.

## Code

### `src/qchem_interfaces/program.rs` (new)

```rust
pub enum QcProgram { Xtb, Behemoth }
```

with `ALL`, `label()`, `binary_hint()` (the file name to look for) and
`config_file()` (where its path is remembered). Adding a third program means
adding a variant and its command construction, nothing else.

### `src/qchem_interfaces/config.rs`

Generalised from one path to one per program: `load_path(program)` /
`save_path(program)`. xTB keeps reading and writing `xtb_path.txt`, so an
existing configuration survives; Behemoth gets `behemoth_path.txt`.

### `src/qchem_interfaces/behemoth.rs` (new)

* `optimize_command(binary, xyz, charge, multiplicity) -> Command`
* `hessian_command(binary, xyz, charge, multiplicity) -> Command`
* `parse_trajectory_energies(xyz_text) -> Vec<f64>` — from the `cycle N E = …`
  comments, so a partial trajectory from a failed run still plots.
* `parse_hessian(stdout, natoms) -> Result<Array2<f64>, String>` — the labelled
  block matrix. Columns are located by counting, rows by their label, and a row
  whose label does not match the expected atom and axis is an error rather than
  a guess.
* `parse_summary(stdout) -> BehemothSummary { converged, cycles, final_energy }`

### Existing modules

`xtb_optimize` and `xtb_freq` gain a program parameter. Their run directories,
cancellation, polling, progress display and result rendering are unchanged --
only the command line and the output parsing branch on the program.

## UI

Both panels get the same two controls in place of the single "xTB path" field:

```
Program  [xTB ▾]   Path [/usr/local/bin/xtb        ] [Browse…]
```

The path field shows and edits the path of the *selected* program, so switching
programs switches paths. The frequency panel's "Run Freq Calc (with xTB)"
becomes "Run Freq Calc" plus the program name.

The optimizer and the frequency panel keep their own program selection: running
a Hessian in a different program than the optimization is a legitimate thing to
want, and coupling them would prevent it.

## Validation

Behemoth's own `--freq` on water (GFN1-xTB) reports 1506.27, 3551.05 and
3646.50 cm⁻¹ with IR intensities 46.07, 6.10 and 18.26 km/mol, and a ZPE of
0.019829 Eh. Beavyr's analysis of the Hessian from the same geometry must
reproduce the frequencies to within about a cm⁻¹ and the ZPE to a fraction of
a kcal/mol, which is the same standard the ORCA and Gaussian readers are held
to. A captured fixture of the Hessian output goes into `tests/fixtures/` so the
test does not need the binary.

Parser tests cover: a full water Hessian round trip; symmetry of the result; a
matrix whose row labels disagree with the atom order; a truncated block; and
trajectory comments both well-formed and missing.

## Out of scope

Behemoth's other methods (`hf`, `mp2`, `rks`, `uks`, …) and its many task
flags (`--opt-ts`, `--rpa`, `--stddft`, `--moldyn`, scans, IRC). Each needs its
own UI for basis sets, functionals and task parameters. The point of this change
is the plumbing: a program list, a path per program, and one non-xTB backend
proving the seam works.
