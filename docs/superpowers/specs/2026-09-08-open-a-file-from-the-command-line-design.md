# Opening a File From the Command Line — Design

Date: 2026-09-08
Status: implemented 2026-09-08

## Goal

Run `./Beavyr formic_acid.xyz` and have the structure on screen when the window
opens. Run `./Beavyr water.hess` and have the frequency analysis already done,
with the Vibrations window open. Run `./Beavyr orbitals.molden` and have the
surface tool loaded and ready.

Today `main()` ignores its arguments entirely. Everything is loaded through a
file dialog after the window is up, which means a structure a user already has
in hand takes four clicks to look at, and Beavyr cannot be the thing a file
manager or a shell alias opens a file *with*.

## What the user types

```
./Beavyr                          # empty viewport, exactly as now
./Beavyr formic_acid.xyz          # structure loaded
./Beavyr trajectory.xyz           # multi-frame file, trajectory ready to play
./Beavyr water.hess               # frequencies analysed, Vibrations open
./Beavyr mo.molden                # orbitals loaded, Surface open
./Beavyr mol.xyz mo.molden        # several files, each into its own tool
./Beavyr --help                   # what it accepts
```

No flags beyond `--help`. Arguments are paths, in any order.

## Recognising the file

By extension first, because that is what the user thinks in, then by content
where the extension does not settle it. The mapping:

| Extension | Read as | Window opened |
|---|---|---|
| `.xyz` | structure; multi-frame becomes a trajectory | Structure, plus Trajectory if it has more than one frame |
| `.hess` | ORCA Hessian | Vibrations |
| `.molden`, `.molden.input` | orbitals and basis | Surface Tools |
| `.engrad` | ORCA Cartesian gradient | Vibrations |
| `.out`, `.log` | ambiguous — see below | depends |

`.inp` is deliberately **not** accepted. It is an ORCA *input* file, and the
Molden dialog's existing `inp` filter is for the `molden.input` spelling, not
for ORCA inputs. Guessing between the two would be wrong about half the time.

### The ambiguous ones

`.out` and `.log` are written by every program in the field. Beavyr already
sniffs content in one place — `gaussian_log::is_gaussian_log` — and this
extends the same idea. In order:

1. **Gaussian output** (`is_gaussian_log`) → frequencies and modes → Vibrations.
2. **ORCA output containing excited states** (`uvvis::orca::parse_orca_tddft`
   succeeds) → absorption spectrum → UV-Vis.
3. Neither → **report it, do not guess**. The message names the file and says
   what was looked for. An ORCA `Freq` output is the case that lands here: the
   force constants are in the companion `.hess`, not in the `.out`, and telling
   the user that is more useful than silently opening nothing.

A single file is never loaded into two tools. An ORCA output with a `Freq` job
carries no Hessian, so cases 1 and 2 cannot both match the same file.

### `.engrad` on its own

A gradient is only meaningful against a Hessian — it is what makes the
reaction-path projection possible. If an `.engrad` is passed without a `.hess`,
it is read, held, and the user is told it is waiting for a Hessian. It is not
silently dropped, and it does not become the structure on screen.

## Several files at once

Each file is loaded into the tool that owns it, in the order given on the
command line. Three of the five kinds carry a geometry (`.xyz`, `.hess`,
`.molden`), and where more than one does, **the last one on the command line is
the structure that ends up on screen** — the others still load their own data.
This is the only rule that is simple enough to predict without reading
documentation, and the common case is one file, where it does not arise.

## What has to change in the code

Every loader is currently shaped as *pick a path, then do the work*:
`load_hessian_from_dialog`, `load_from_dialog` (UV-Vis), and the Molden load
written inline in the button handler. None of them can be called with a path we
already have.

Each is split in two:

* `load_hessian_from_path(path, ...)` — the work, callable from anywhere.
* `load_hessian_from_dialog(...)` — picks a path, then calls the above.

The dialog versions keep their present behaviour exactly, so nothing the user
does today changes. This is a refactor with no behavioural content, and the
existing tests cover the parsing underneath it.

The Molden load is currently written inline inside the button closure and moves
out to a named function as part of this.

## Where it runs

A Bevy `Startup` system, ordered **before** `center_camera_on_startup` so the
camera frames the molecule that was loaded rather than an empty scene. The
arguments are read once, in that system, via `std::env::args()`. No argument
parsing crate: this is one optional flag and a list of paths, and Beavyr's
dependency list is short on purpose.

`Startup` runs after resources are inserted, so `Molecule`, `UiLayout`,
`XtbFreqPanelState`, `OrbitalState` and `UvVisState` are all available to write
into.

## When it goes wrong

Three failures, all of which must be visible:

* the path does not exist, or cannot be read;
* the extension is not one Beavyr knows;
* the extension is known but the file will not parse.

Each is written to `stderr` — a user who launched from a shell is looking at
one — **and** collected into a new `StartupLoadReport` resource, shown as a
dismissible block at the top of the Structure panel, which is the one window
that is always there. Beavyr still opens; a bad argument is not a reason to
refuse to start, and an empty viewport with no explanation is exactly the
silent failure this project's diagnostics exist to avoid.

`--help` prints the table above and exits before the window opens.

## Testing

The dispatch is a pure function and gets the tests:

```rust
fn classify(path: &Path, head: &str) -> FileKind
```

taking the path and the first few KB of the file, returning the kind. Cases:
each extension maps to its kind; `.out` with a Gaussian banner is Gaussian;
`.out` with an ORCA absorption block is TD-DFT; `.out` with neither is
`Unrecognised`; an unknown extension is `Unrecognised`; case-insensitive
extensions (`.XYZ`, `.HESS`) work, since a file that came off a Windows share
or a cluster often is capitalised.

Fixtures already exist for every one of these under `tests/fixtures/` and
`examples/`, so no new test data is needed.

The loading itself is not unit-tested: it writes into Bevy resources, and the
parsers beneath it are already covered by 628 tests.

## Out of scope

* Watching the file for changes and reloading.
* A file-association / MIME entry so a desktop double-click opens Beavyr. That
  is packaging, and belongs with whatever eventually installs the binary.
* Reading `.gbw` or any other binary ORCA file.
* Any flag that runs a calculation without the window. That is the batch
  question, and the decision on 2026-09-08 was not to build it: batch work
  belongs in Python outside Beavyr.


## What was built, and what changed from this design

Implemented in `src/cli.rs`, wired into `main.rs` as a `Startup` system ordered
before the camera is centred. 26 tests.

One change forced by a real file. The design said the classifier would be given
"the first few KB of the file". That is wrong, and `examples/c164-ts1-2-1001.log`
proves it: the job was restarted with Gaussian's `Link1`, so the file opens with
a continuation header and the word "Gaussian" does not appear until 172 kB in.
An 8 kB peek classified it as an output file with nothing usable in it.

Any fixed peek is a guess about how someone else's program lays out its file.
`.out` and `.log` are now read whole, up to 64 MB. Nothing else is read at all,
because every other extension settles the question by itself. Reading twelve
megabytes once at startup is invisible next to compiling the shaders. There is a
test pinned to that file, and it asserts the banner really is beyond the old
window, so the test cannot quietly stop testing anything.

`FileKind` also gained a variant the design did not have:
`OutputWithNothingUsable`, separate from `Unrecognised`. It is what lets an ORCA
`Freq` output be told "the force constants are in the companion .hess" instead
of "not a file Beavyr reads".
