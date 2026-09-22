# ORCA and PySCF as optimisation and frequency backends

**Date:** 2026-09-22
**Stage:** 1 of 4
**Status:** design, awaiting review

## Why this stage first

Beavyr can currently send a job to xTB, to Behemoth, or to its own DREIDING
force field. The goal across all four stages is to add ORCA, Psi4, PySCF and
Sparrow as well, each wired into every part of Beavyr that calls a
quantum-chemistry program.

This stage does ORCA and PySCF only, and only for geometry optimisation and
frequencies. That pairing is deliberate: ORCA is an ordinary executable, PySCF
is a Python library with no command line at all, so between them they exercise
both ways Beavyr will have to reach a program. Everything built here is the
foundation the remaining three stages sit on. If the shape is wrong, it is
wrong once, on two backends, rather than on six.

## What is in this stage

* ORCA, run as the executable at a path the user browses to.
* PySCF, run by generating a Python script and executing it with `python3`.
* Both offered in the Geometry Optimization panel and the Frequencies panel.
* A per-program statement of what each backend can actually do, so a panel can
  grey out a job rather than start one that is guaranteed to fail.
* A choice, per run, between the program's own optimiser and Beavyr's.

## What is not in this stage

Psi4 and Sparrow (stage 2). Excited states and TD-DFT (stage 3). Conformer
re-ranking and transition-state generation (stage 4). Molpro, which is in the
Smite backend list but was not asked for. The surrogate energy-surface code in
Smite. The xTB force-field retry described in "Known issue left alone" below.

## Background: what already exists and is reused

Two ORCA readers are already in Beavyr, already tested against real files:

* `src/qchem_interfaces/orca_engrad.rs` reads the `.engrad` file: energy,
  Cartesian gradient, geometry, with a tolerance check that the gradient and
  the Hessian belong to the same structure.
* `src/qchem_interfaces/orca_hess.rs` reads the `.hess` file: Hessian,
  geometry, per-atom masses so an isotope survives, ORCA's own frequencies and
  infrared intensities, the multiplicity, the energy and a frequency scale
  factor.

A probe run of ORCA 6.1.1 on water at HF/STO-3G confirmed the `.hess` written
today still carries every section those readers expect: `$hessian`,
`$vibrational_frequencies`, `$normal_modes`, `$atoms`, `$multiplicity`,
`$act_energy`, `$frequency_scale_factor` and `$ir_spectrum`.

So for ORCA this stage writes input and runs the program. It does not write a
new parser for the two files that matter most.

Three further pieces are reused rather than rebuilt:

* `dreiding_run.rs` already shows how to make a backend that does not behave
  like xTB look like one to the panels, by writing its results in the format
  the trajectory player and energy plot already read. Every new backend follows
  that pattern.
* The background-task shape in `xtb_optimize.rs` and `xtb_freq.rs`: a task on
  the compute pool, a child process a cancel flag can kill, a scratch directory
  kept on failure and discarded on success.
* The Eckart projection and eigensolver in `src/normalmode/`. Every backend's
  Hessian goes through it, so mode animation, thermochemistry and infrared
  spectra keep working unchanged.

## Design

### 1. Programs declare what they can do

`QcProgram` gains two variants, `Orca` and `PySCF`, and every variant gains a
statement of its capabilities:

* whether it needs an executable path,
* whether it needs Python,
* whether it can optimise a geometry by itself,
* whether it can produce a Hessian, and whether that Hessian is analytic or
  has to be built from gradients by finite differences.

The panels read this rather than hardcoding which programs appear where. Two
things force it. Sparrow, arriving in stage 2, has gradients and an analytic
Hessian but no geometry optimiser at all, so "can this program optimise?" has
to be a question with an answer rather than an assumption. And PySCF needs no
executable path, which is the same question `runs_in_process` already answers
for DREIDING, generalised.

A test asserts that every program's declared capabilities are consistent: a
program that needs no binary offers no path field, a program that cannot
optimise on its own forces Beavyr's optimiser on.

### 2. How Beavyr reaches Python

Beavyr writes a small Python script into the run directory beside the geometry
and executes it with `python3`, taken from the environment Beavyr itself was
launched in.

**Changed while building, from the JSON originally specified.** The script
writes its results in formats Beavyr already reads, and prints only scalars:
the optimisation path as a multi-frame XYZ with `cycle <n> E = <v> Eh`
comments; the Hessian behind a `$hessian` header, which is Turbomole's layout;
the final geometry as an ordinary XYZ; and the energy on standard output as
`energy <value>`, read by key.

There is no path field for a Python backend and no environment management.
Activating the right environment before launching Beavyr is the user's
responsibility, as agreed.

#### The environment is checked before the run, not after

The trap this design has to avoid is the quiet one: the user activates the
right environment in a terminal, then launches Beavyr from a desktop icon,
which inherits the desktop session's Python instead. The job then fails minutes
later with an import error that looks like a broken installation rather than a
launch-environment mismatch.

So the check happens when the panel is opened, not when the job is started.
`python3 -c "import importlib.util; ..."` using `find_spec` answers "is this
module here?" without importing it, measured at 23 milliseconds against this
machine's PySCF. A real import costs 0.7 seconds, so the cheap form is what
runs. The result is cached for the session and refreshed by a "check again"
button.

Three outcomes, three messages:

* **`python3` is not on the path at all.** The panel says Beavyr runs Python
  backends with the `python3` of the environment it was started in, and that no
  `python3` was found there. The Run button is disabled.
* **`python3` is there, the module is not.** The message names both the module
  and the interpreter that was searched, for example that PySCF is not
  installed in `/home/peter/.venvs/science`, and says to activate the
  environment that has it and restart Beavyr. Naming the interpreter is the
  point: it is what makes a launch-environment mismatch visible instead of
  looking like a broken install. The Run button is disabled.
* **The module is there.** The panel states the module version and the
  environment it was found in, so the right environment is confirmed at a
  glance before a long job is started rather than assumed.

A backend is never silently swapped for another when its environment is
missing. The run is blocked and the reason is named.

One further consequence, stated in the panel rather than discovered: the PySCF
environment and the Psi4 environment on this machine are separate, so one
Beavyr session can reach one or the other, never both.

#### Failures once the script is running

A module that imports but then fails, an SCF that does not converge, a
geometry optimiser that gives up, are reported with the program's own message.
Python's traceback is captured to the run directory, which is kept on failure,
and the last line of the exception is what the panel shows.

Reusing those three formats was chosen over JSON because Beavyr carries no JSON
reader, so JSON would have meant writing one. Reusing what is already parsed
meant writing no new parser at all, and it leaves a run directory a person can
read directly -- which matters, because that directory is deliberately kept
when a run fails. Numbers are written at full double precision.

### 3. Geometry optimisation

Two routes, chosen per run, defaulting to the program's own.

**The program's own optimiser.** ORCA is given `Opt` on its keyword line and
run once. PySCF is driven through its `geometric` solver, confirmed present in
the science environment, with a callback recording each step's geometry and
energy.

**Beavyr's optimiser.** Beavyr's existing redundant-internal minimiser
(`optimizer/minimum_redundant.rs`, falling back to
`optimizer/minimum_cartesian.rs` where internal coordinates cannot be built)
drives the backend one gradient at a time, through the same energy-and-gradient
object the transition-state tool already uses. This is what makes constraints
and relaxed scans work with every program, and from stage 2 it is the only
route Sparrow has.

Whichever route ran, the result is written as a multi-frame XYZ whose comment
lines read `cycle <n> E = <value> Eh`, which is what Beavyr's trajectory
player, energy plot and summary window already read. ORCA's own trajectory
comment reads `Coordinates from ORCA-job <name> E <value>`, confirmed by probe
run, so it is translated rather than parsed in a second format. This keeps
every consumer downstream untouched, exactly as `dreiding_run.rs` does.

### 4. Frequencies

ORCA is given `Freq` and the resulting `.hess` goes through the existing
reader. PySCF's analytic Hessian, confirmed working and symmetric to eleven
decimal places on the probe, is written behind the `$hessian` header.

Both then go through Beavyr's own projection and eigensolver, which is what
gives animated modes, thermochemistry that follows the panel's temperature and
cutoff, and an infrared spectrum. ORCA additionally supplies real infrared
intensities through the existing reader; PySCF's route reports no intensities
in this stage, and the panel shows them as unavailable rather than as zeroes.

Where a geometry optimisation and a frequency job run together, the frequencies
belong to the optimised geometry, and the geometry recorded with the Hessian is
the one used for the projection. The existing tolerance check in the `.engrad`
reader is the model: a Hessian from one structure and a gradient from another
produces plausible-looking nonsense, and is caught rather than analysed.

### 5. Method catalogues

One file per program, as chosen: `orca_method.rs` and `pyscf_method.rs`. Each
is written against that program's own documentation and checked against the
installed copy, the way the Behemoth catalogue was built. Switching program
resets the method to that program's default; the catalogues do not try to
translate between each other.

Both programs take a functional and a basis set, so both reuse the existing
field-visibility and validation machinery rather than inventing a second one.
ORCA additionally takes a free-text keyword line for anything the dropdowns do
not cover, following the `additional` field in the Smite interface.

### 6. Error handling

No silent substitution anywhere. A failed run reports what failed and stops.

Three specific traps, all taken from the Smite ORCA interface and confirmed on
the probe run:

* **ORCA can exit zero and still have failed.** The exit status is not
  sufficient. The output must also be checked for
  `ORCA finished by error termination`, which is what the Smite code does.
* **A restart from previous orbitals can poison a run.** Reading a stale `.gbw`
  fails with `error termination in GUESS` or `No orbitals were found`. The
  Smite code detects exactly this, discards the orbital file and reruns from a
  fresh guess. That is a legitimate retry because it changes only the starting
  guess, not the method, and it is carried over. It is reported when it fires.
* **A missing Python module is not a reason to use another backend.** It is
  reported with the module's name and the interpreter that was searched, before
  the run rather than after it. See "The environment is checked before the run"
  above.

The same rule applies to ORCA's executable, through the path checking Beavyr
already does for xTB and Behemoth: an empty or wrong path disables the Run
button and says so, rather than failing at spawn time. The two cases are the
same question asked of two different kinds of backend, and both are answered
before a job starts.

### 7. Testing

Every parser gets a test against real captured output, following the existing
`tests/fixtures/` convention. The probe runs performed while writing this
design supply them: an ORCA optimisation and frequency job on water at
HF/STO-3G, giving a real `.out`, `.hess`, `.engrad`, `_trj.xyz` and `.xyz`, and
a PySCF optimisation and analytic Hessian on the same molecule.

Command construction is tested without running anything, as
`behemoth.rs` already tests its own: build the command, assert on its
arguments. The generated Python script is tested the same way, by asserting on
the text produced for a given configuration.

## Files

New:

* `src/qchem_interfaces/orca_input.rs` — writes the ORCA input file.
* `src/qchem_interfaces/orca_run.rs` — runs ORCA, reads its results back.
* `src/qchem_interfaces/orca_method.rs` — ORCA's method catalogue.
* `src/qchem_interfaces/pyscf_run.rs` — generates and runs the Python script.
* `src/qchem_interfaces/pyscf_method.rs` — PySCF's method catalogue.
* `src/qchem_interfaces/capabilities.rs` — what each program can do.

Changed:

* `program.rs` — two new variants, capabilities, path and Python rules.
* `xtb_optimize.rs` — optimiser choice, dispatch to the new backends.
* `xtb_freq.rs` — dispatch to the new backends.
* `method_ui.rs` — per-program catalogue selection.
* `mod.rs` — the new modules.

Deliberately untouched: `orca_engrad.rs` and `orca_hess.rs`, which already do
their job.

### The three files already sitting in this folder

`orcarun.rs`, `psi4run.rs` and `pyscfrun.rs` are in
`src/qchem_interfaces/` but are not listed in `mod.rs`, so nothing compiles them
and nothing can call them. They are direct ports of the Smite Python files and
they carry the same Tautomer-era assumptions: a fixed working directory, fixed
file names, and no cancellation.

They are superseded by the new files rather than extended, because the new ones
have to fit Beavyr's background-task and cancellation rules, which these
predate. Keeping them would leave two ORCA drivers in one folder, one of them
dead, with confusingly similar names.

Recommendation: delete them. They are in git history if ever wanted. This is
flagged rather than done, because deleting the user's code is the user's call.

## Known issue left alone

`call_xtb_in` in `xtbrun.rs` silently reruns a failed xTB calculation with the
GFN-FF force field at a raised electronic temperature and loose accuracy, and
returns that energy as though the requested method had produced it. The
conformer re-ranking in `refine.rs` ranks those numbers against genuine GFN2
energies. This is noted here so it is written down, and is deliberately outside
this stage: it is a decision about existing behaviour, not about new backends.

## Open question for review

The frequency panel currently shows infrared intensities whenever the source
supplies them. PySCF can produce them, but it needs a separate dipole-derivative
calculation on top of the Hessian, roughly doubling the cost of a frequency
job. This stage leaves them out and marks them unavailable. Say if you would
rather pay for them.
