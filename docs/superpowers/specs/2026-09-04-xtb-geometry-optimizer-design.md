# xTB geometry optimizer — design

Status: approved by user (chat), pending written-spec review.
Author: assistant, in collaboration with Peter Szabo (peter.szabo@kuleuven.be).

## 1. Purpose

Add a geometry-optimization tool to Beavyr, backed by Grimme's xtb binary,
modeled on the existing feature in the sibling project `Tautomer`
(`/home/peter/Dropbox/Research_Leuven/Tautomer/src/qchem_interface/xtbrun.rs`
and its driving UI in `src/ui.rs`). The user gives the path to an xtb
executable and the charge/multiplicity of the molecule currently on screen;
Beavyr runs the optimization in the background, and on success loads the
optimization trajectory into the existing trajectory player and plays it once
automatically.

Unlike Tautomer, which only reads back the final geometry (`xtbopt.xyz`),
Beavyr already has a trajectory player (`src/trajectory.rs`), so this feature
also surfaces the *path* of the optimization (`xtbopt.log`), not just its
endpoint.

Before running, Beavyr independently checks that the given charge and
multiplicity are physically possible for the drawn structure (electron count
parity, multiplicity vs. available electrons) and flags — without blocking —
valences that look chemically unusual. This check exists because xtb performs
none of its own: empirically verified below (§3), it silently "succeeds" on
an impossible spin state.

## 2. Non-goals

- No single-point energy/gradient/Hessian calculations. Only `--opt`.
- No replacement of the dead `src/qchem_interfaces/` module (it is not wired
  into `main.rs`'s module list at all, so it currently does not even compile
  into the binary; it targets ORCA/PySCF/Psi4 energy calls, not xtb
  optimization). It is left untouched.
- No support for other engines (ORCA, Psi4, PySCF) in this pass — xtb only,
  matching the user's request.
- No per-atom formal-charge UI. Beavyr's `Molecule` (`src/molecule.rs`) has no
  per-atom charge field (unlike Tautomer's `atom.charge`); electron counting
  works purely from element identity and the structure's own bond graph, per
  the user's own framing ("number of atoms + electron based on the drawn
  structure and its valency of atoms").

## 3. Background research (facts gathered before designing)

- **Tautomer's xtb runner** (`xtbrun.rs`, 174 lines) is a clean reference:
  resolves the executable (absolute path or `$PATH` lookup), creates a
  timestamped scratch directory, writes `input.xyz`, spawns
  `xtb input.xyz --opt --chrg <n> [--uhf <n>]` with stdout/stderr redirected to
  files, polls the child non-blockingly with a cancel flag checked every
  50 ms, and on success reads back `xtbopt.xyz`. This shape (background
  `Task`, `Arc<AtomicBool>` cancel, `Arc<Mutex<Option<Child>>>` for kill) is
  reused directly; see §5.
- **Tautomer's validation** (`threed_structure/mod.rs::validate_xtb_electronic_state`)
  computes total electrons from atomic numbers plus each heavy atom's
  attached-hydrogen count, subtracts the charge, and checks:
  multiplicity ≥ 1; electrons > 0 after charge; unpaired electrons
  (`multiplicity - 1`) ≤ available electrons; and electron/multiplicity
  parity (`(electrons + multiplicity) % 2 == 1`). This is the hard-block
  logic Beavyr adopts, adapted to work off Beavyr's own bond-graph
  representation instead of Tautomer's explicit per-atom charge field.
- **I ran xtb myself** (`/home/peter/Programs/orca_6_1_1_linux_x86-64_shared_openmpi418/xtb`,
  version 6.7.1) against `examples/h2o2.xyz`:
  - `xtb input.xyz --opt --chrg 0` succeeds, exit 0, and writes: `xtbopt.log`
    (the trajectory — 17 frames, standard multi-frame XYZ, each frame's
    comment line reading `energy: <Eh> gnorm: <norm> xtb: 6.7.1 (edcfbbe)`),
    `xtbopt.xyz` (final frame only), `xtbrestart`, `xtbtopo.mol`, `charges`,
    `wbo`, plus whatever is redirected to `xtb.stdout`/`xtb.stderr`.
  - `xtbopt.log`'s frame format is plain multi-frame XYZ and needs no new
    parser: `trajectory::parse_multi_xyz` (`src/trajectory.rs:96`) already
    accepts arbitrary comment-line content, verified by inspection of its
    `looks_like_coord` / count-line logic.
  - **Critically**, `xtb input.xyz --opt --chrg 0 --uhf 1` on H2O2 (18
    electrons, an even-electron closed-shell molecule for which multiplicity
    2 / one unpaired electron is impossible) exits 0, prints "normal
    termination of xtb", and its own log explicitly states
    `unpaired electrons   1` while producing the *same total energy* as the
    correct closed-shell run — i.e. it silently ignored the impossible
    request rather than erroring. This is the concrete justification for
    Beavyr doing its own hard-block validation rather than trusting xtb's
    exit code.
- **Beavyr has no on-disk persistence today** (confirmed: no `serde`/`ron`
  save-to-disk call anywhere in `src/`, `settings.rs` is in-memory only). The
  xtb executable path is the first thing worth remembering across runs, so
  this feature also introduces Beavyr's first (deliberately tiny) config
  file.
- **Beavyr already has the right background-task idiom** in
  `src/orbitals/mod.rs`: `AsyncComputeTaskPool::get().spawn(async move {...})`
  polled once per frame via `block_on(future::poll_once(&mut task))`. The xtb
  task follows this exact pattern instead of Tautomer's raw
  `bevy::tasks::Task` field (functionally the same primitive; matching
  Beavyr's own existing style rather than importing Tautomer's).
- `src/molecule.rs::covalent_radius_angstrom` already has a 103-element match
  arm in periodic-table order; a parallel `atomic_number(symbol) -> Option<u32>`
  table is added alongside it for electron counting (no existing table to
  reuse — grepped, none found).

## 4. User-facing behavior (from clarifying questions)

1. **On success**: the trajectory is loaded into `TrajectoryState`, the app
   switches to the Trajectory tool/panel, a message states the optimization
   succeeded (with final energy), and the trajectory auto-plays once from
   frame 1. The user can replay or scrub freely afterward, exactly as with
   any other loaded trajectory. The on-screen `Molecule` becomes the
   optimized (final-frame) geometry, matching normal trajectory-load
   behavior already in `src/ui.rs` (`apply_current_frame`).
2. **Valence findings** (over-valent atoms, radical-implying odd bond counts
   inconsistent with the given multiplicity, etc.) are warnings only — shown
   in the panel, never blocking the Optimize button. Only physically
   *impossible* electron/multiplicity combinations (parity, insufficient
   electrons for the requested multiplicity) are hard errors that block the
   run.
3. **xtb path** is stored in a small RON config file under the OS config
   directory (via the `dirs` crate — not currently a dependency; added as a
   minimal, no-proc-macro, single-purpose addition), read on startup and
   written whenever the user changes/confirms the path in the panel. The UI
   field is still always editable and pre-fills from
   `$BEAVYR_XTB_PATH` → saved config → plain `"xtb"` (PATH lookup), in that
   order, mirroring Tautomer's `default_xtb_executable`.

## 5. Architecture

New module `src/qchem/` (the existing `src/qchem_interfaces/` is dead code
for a different purpose and is left alone, per §2):

```
src/qchem/
  mod.rs        // declares submodules, re-exports the small public surface
  xtb_run.rs    // process spawning, background task resource, cancellation
  valence.rs    // pure validation: electron counting, hard errors, warnings
  config.rs     // XtbConfig: load/save the RON file with the executable path
```

### 5.1 `valence.rs` — pure, synchronous, fully unit-testable

```rust
pub struct ElectronicStateError(pub String);        // hard block
pub struct ValenceWarning {                          // never blocks
    pub atom_index: usize,
    pub message: String,                             // e.g. "carbon C7 has 5 bonds; expected 4"
}

/// Returns the xtb `--uhf` value (unpaired electrons) on success.
/// Electron count comes from atomic number per atom, summed, minus `charge`.
/// Bonds come from `mol.bonds` (already computed by `Molecule::recompute_bonds`),
/// not from any drawn/declared valence -- Beavyr has no separate valence model.
pub fn validate_electronic_state(
    mol: &Molecule,
    charge: i32,
    multiplicity: i32,
) -> Result<(i32, Vec<ValenceWarning>), ElectronicStateError>
```

Hard-block checks (ported from Tautomer's `validate_xtb_electronic_state`,
adapted to Beavyr's data model — no per-atom charge field, no
`is_pseudoatom`/custom-label concept, so those two Tautomer-specific checks
are dropped and replaced with "every symbol must resolve via the new
`atomic_number` table, else it's an unsupported element for xtb"):
- multiplicity < 1 → error
- any atom symbol without a known atomic number → error (xtb cannot run on it)
- electrons (after charge) ≤ 0 → error
- unpaired electrons (`multiplicity - 1`) > available electrons → error
- `(electrons + multiplicity) % 2 != 1` → error (parity)

Warning-only checks (new — Tautomer has no equivalent, since it works from
explicit user-declared charges/valence rather than inferring from bonds):
- an atom's bond count in `mol.bonds` differs from its "typical" valence
  for common organic elements (H:1, C:4, N:3, O:2, F/Cl/Br/I:1, S:2 or 6,
  P:3 or 5) by more than what the *given* multiplicity's unpaired-electron
  budget can explain. Concretely: sum the "valence deficit" across all atoms
  (typical − actual bonds, floored at 0) and compare against
  `multiplicity - 1`; a mismatch is reported per-atom rather than as one
  aggregate message, so the user knows which atom to fix.
  This table is intentionally small and heuristic (organic elements only);
  an atom outside the table contributes no warning (never a false positive
  for chemistry the heuristic doesn't understand).

### 5.2 `xtb_run.rs` — background process, ported from Tautomer

Structures, ported near-verbatim from `Tautomer/src/qchem_interface/xtbrun.rs`
(§3), adapted to Beavyr's `AsyncComputeTaskPool` idiom (§3) instead of a bare
`Task` field:

```rust
#[derive(Resource, Default)]
pub struct XtbOptimizationTask {
    task: Option<Task<Result<XtbOptimizationOutput, String>>>,
    cancel: Option<Arc<AtomicBool>>,
    child: Option<Arc<Mutex<Option<Child>>>>,
    run_dir: Option<PathBuf>,
}

pub struct XtbOptimizationOutput {
    pub trajectory_text: String,   // contents of xtbopt.log
    pub final_energy_hartree: Option<f64>,
}
```

Functions ported directly (same signatures, same behavior, `PathBuf`/`String`
error style unchanged): `resolve_xtb_executable`, `create_xtb_run_dir` (scratch
dir under Beavyr's own state location — see §5.4), `discard_xtb_run_dir`. The
actual spawn+wait+cancel loop (`run_xtb_optimize_cancellable` in Tautomer) is
ported with one change: it returns the *log* file's contents
(`xtbopt.log`, the trajectory) rather than `xtbopt.xyz` (final-frame-only),
since the trajectory is what Beavyr's player consumes; the final energy is
parsed from that same log's last frame's comment line
(`energy: <value> ...`) for the success message in §4.1, with `xtbopt.xyz`
read as a fallback only if `xtbopt.log` is somehow absent (defensive, since a
successful `--opt` run always writes both, per §3).

A Bevy system (added to the existing `Update` schedule, alongside
`orbitals`'s polling system) polls the task once per frame via
`block_on(future::poll_once(...))`; on `Some(Ok(output))` it parses the
trajectory text with the existing `trajectory::parse_multi_xyz`, populates
`TrajectoryState`, sets `playing = true` (kicking off the existing playback
system unchanged), and writes the success message; on `Some(Err(msg))` it
surfaces `msg` as an error banner and leaves the current molecule/trajectory
untouched.

### 5.3 `config.rs` — Beavyr's first on-disk state

```rust
#[derive(Serialize, Deserialize, Default)]
pub struct XtbConfig {
    pub executable_path: String,   // may be empty: falls back to env/PATH
}

pub fn config_path() -> Option<PathBuf>;   // dirs::config_dir().join("beavyr").join("xtb.ron")
pub fn load_xtb_config() -> XtbConfig;     // missing/corrupt file -> XtbConfig::default(), never errors
pub fn save_xtb_config(cfg: &XtbConfig);   // best-effort; failure logs, never panics/blocks UI
```

Adds the `dirs = "5"` dependency (small, no proc-macros, widely used already
by `rfd` transitively — checked in `Cargo.lock`, so its version footprint is
already present in the dependency graph). `ron` is already a dependency.

### 5.4 Scratch directory

Mirrors Tautomer's `xtb_scratch_dir` reasoning exactly (explicitly NOT
`CARGO_MANIFEST_DIR`, which is a compile-time path): use
`dirs::cache_dir().join("beavyr").join("xtb")` (falls back to a temp dir if
`dirs::cache_dir()` is unavailable), one subdirectory per run. **Confirmed**
by reading Tautomer's own call site (`src/ui.rs:2561`, comment: "The geometry
is in hand, so the scratch files have served their purpose. A failed run
keeps its directory for the logs."): `discard_xtb_run_dir` is called only on
the success path. Beavyr adopts the identical rule — delete the run
directory after a successful parse, keep it on any failure or cancellation so
`xtb.stdout`/`xtb.stderr` remain inspectable.

## 6. UI

New collapsible section in the existing left panel (`src/ui.rs`), placed near
the Molecule Builder section (both operate on "the molecule currently on
screen"):

- Charge: `DragValue<i32>`, default 0.
- Multiplicity: `DragValue<i32>`, range 1.., default 1.
- "Auto-fill from structure" button: sets charge to 0 and multiplicity to 1
  (the defaults) — present mainly so a user who fiddled with the fields has
  an obvious reset, since Beavyr (unlike Tautomer) has no per-atom charge to
  sum.
- xtb executable path: text field, pre-filled per §4.3; changes persist to
  `config.rs` on edit (debounced to on-blur/on-Enter, not every keystroke).
- Validation output: red banner for a hard error (from `ElectronicStateError`)
  with the Optimize button disabled; amber list of `ValenceWarning`s
  (never disabling) shown above the button when present.
- "Optimize" button: disabled while a task is already running or a hard
  error stands; spawns the task.
- "Cancel" button: visible only while running; calls `cancel_now()`.
- Status line: idle / running (with elapsed time, matching Tautomer's pattern
  if it has one — otherwise just "Running…") / the success message from
  §4.1 / the last error message.

## 7. Data flow summary

```
User sets charge/multiplicity (or leaves 0/1)
        │
        ▼
validate_electronic_state(mol, charge, multiplicity)
        │                              │
   Err(hard block)               Ok((uhf, warnings))
        │                              │
  banner, button disabled      warnings shown, button enabled
                                       │
                              user clicks Optimize
                                       │
                    resolve_xtb_executable (config/env/PATH)
                                       │
                     spawn background task (AsyncComputeTaskPool)
                                       │
              xtb <workdir>/input.xyz --opt --chrg <c> [--uhf <u>]
                                       │
                         poll each frame (non-blocking)
                                       │
                    ┌──────────────────┴──────────────────┐
                    ▼                                      ▼
              success: read xtbopt.log              failure: exit code / missing file
                    │                                      │
     parse_multi_xyz -> TrajectoryState              error banner, molecule untouched
     switch to Trajectory tool, play once
```

## 8. Testing plan

- `qchem::valence` — pure unit tests, no process/filesystem:
  - **Checked**: Tautomer's `threed_structure/mod.rs` test module has no
    tests targeting `validate_xtb_electronic_state` specifically (its ~20
    `#[test]` functions cover geometry/topology instead), so there is no
    existing coverage to port — this test set is written fresh, covering:
    parity failure, multiplicity-too-high, zero-electron-after-charge,
    unknown-element rejection, and the exact H2O2 `--uhf 1` scenario from §3
    (18 electrons, multiplicity 2 requested) as the regression case that
    motivated this feature.
  - warning path exercised against `examples/phenyl-OCH3_radical.xyz`
    (already in the repo, confirmed 7 C / 7 H / 1 O — a genuine radical
    example) and a clean closed-shell example (`examples/water.xyz`) to
    confirm no false-positive warning.
- `qchem::xtb_run` trajectory handling — unit test against a **captured real
  fixture**: the `xtbopt.log` generated in §3 while researching this spec is
  saved as `examples/h2o2_xtbopt.log` (or a `tests/fixtures/` location — TBD
  at implementation time) and parsed with the existing
  `trajectory::parse_multi_xyz`, asserting frame count (17) and that the
  final frame's energy matches the value observed in §3
  (-9.054669368748 Eh), rather than a hand-written synthetic stub.
- Process-spawning itself: only exercised as an integration test gated on
  finding a real `xtb` on `PATH` or via `$BEAVYR_XTB_PATH` (skipped
  otherwise, printing why). **Confirmed pattern to match**, from
  `src/orbitals/mod.rs`'s real-file tests: `if !Path::new(EXAMPLE).exists() {
  eprintln!("skipping: {EXAMPLE} not present"); return; }` — the same shape,
  substituting an executable-resolves check for a file-exists check.
- Config round-trip (`config.rs`): save then load in a temp dir, assert
  equality; missing/corrupt file falls back to default without panicking.

## 9. Open questions / risks (flagged, not blocking)

- Fixture location for the captured `xtbopt.log` (§8) — `examples/` next to
  the other sample structures, or a dedicated `tests/fixtures/` — is a minor
  call left to implementation time.
- The "typical valence" table for warnings (§5.1) is a heuristic and will
  need tuning; it is explicitly scoped as best-effort and non-blocking so a
  wrong warning is a false alarm, never a false block.
