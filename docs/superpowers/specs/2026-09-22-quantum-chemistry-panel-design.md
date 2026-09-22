# A general quantum-chemistry panel

**Date:** 2026-09-22
**Status:** built; ORCA wired, the other programs not yet

## What this is

One panel that runs any calculation any installed code can do on the structure
currently on screen. The existing panels stay exactly as they are: the
optimizer, the frequency tool, TS generation and UV-Vis are each shaped around
one job and are quicker for it. This is the general one, for everything those
do not cover and for the combinations they cannot express.

It opens from a new round button in the top-right corner, to the left of the
molecule editor's pencil, and its window opens on the right the way the editor
does.

## The panel, top to bottom

1. **Program.** xTB, Behemoth, DREIDING, ORCA, PySCF. Psi4 and Sparrow join
   when their stage lands.
2. **Charge** and **multiplicity**, validated as the other panels validate them.
3. **Method.** What the chosen program offers at the level above a functional:
   Hartree-Fock, DFT, MP2, RI-MP2, MP2-F12, coupled cluster, or a
   semi-empirical parametrisation. The list is per program.
4. **Functional**, shown only when the method is DFT, grouped by family exactly
   as the optimizer panel now groups them.
5. **Basis set**, shown whenever the method needs one, grouped by family. A
   composite functional hides it, because it brings its own.
6. **Extra keywords**, a free-text line passed through verbatim, for whatever
   the dropdowns do not cover.
7. **Job.** Single-point energy, optimisation, frequencies, optimisation and
   frequencies, transition-state optimisation, intrinsic reaction coordinate,
   or TD-DFT.
   Jobs the chosen program cannot do are shown disabled with the reason, not
   hidden -- a user comparing codes should be able to see that ORCA does an IRC
   and xTB does not.
8. **Excited-state options**, shown only when the job is TD-DFT: how many
   states, and whether to use the Tamm-Dancoff approximation.
9. **Resources**: memory and threads, where the program takes them.
10. **Environment**: for a Python backend, the line already written for the
    other panels, naming the module and the interpreter and disabling the run
    if it cannot import.

## What each program can actually do

Established by running them, not from documentation.

| Job | xTB | Behemoth | DREIDING | ORCA | PySCF |
|---|---|---|---|---|---|
| Single-point energy | yes | yes | yes | yes | yes |
| Optimisation | yes | yes | yes | yes | yes |
| Frequencies | yes | yes | yes | yes | yes |
| Optimisation + frequencies | yes | yes | yes | yes | yes |
| Transition-state optimisation | no | no | no | yes | no |
| Intrinsic reaction coordinate | no | no | no | yes | no |
| TD-DFT | no | sTDA only | no | yes | yes |

ORCA's `OptTS` and `IRC` were both run on water and terminated normally.
`TDDFT` is **not** an ORCA keyword: the calculation is requested by a `%tddft`
block alone, alongside an ordinary DFT keyword line, and `! TDDFT` is an input
error. That was found by running it.

Behemoth's entry says "sTDA only" and the panel will say so too. Its simplified
methods are not TD-DFT and are labelled as what they are.

## Where results go

The panel runs the job and hands the result to whichever existing viewer owns
it, rather than growing its own copies:

* an optimisation, a transition-state optimisation or an IRC becomes a
  trajectory, played by the trajectory tool and plotted by the energy plot,
  through the one `cycle <n> E = <v> Eh` format every backend already writes;
* frequencies go to the Vibrations panel, through Beavyr's own projection and
  eigensolver, so modes animate and thermochemistry follows its controls;
* a TD-DFT run goes to the UV-Vis panel.

That last one needs no new parsing at all. A real ORCA TD-DFT run was checked
against Beavyr's existing ORCA excited-state reader: it writes the
`TD-DFT/TDA EXCITED STATES` heading, `STATE n: E= ... au ... eV ... cm**-1`
records and the `ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS`
table that the reader already keys on.

## Methods per program

Each verified against the installed program.

* **ORCA** — HF, DFT, MP2, RI-MP2, MP2-F12, CCSD(T), DLPNO-CCSD(T),
  CCSD(T)-F12. RI-MP2, the F12 methods and the coupled-cluster methods all need
  an auxiliary basis, and the F12 ones additionally need a CABS set; these are
  supplied automatically, the same way the double hybrids' `/C` set already is.
* **PySCF** — HF, DFT, MP2.
* **xTB** — GFN1, GFN2, GFN-FF.
* **Behemoth** — TASI, GFN1-xTB, HF, DFT, MP2, RI-MP2.
* **DREIDING** — none. It is a force field; the method rows disappear.

## Rules

Nothing is silently substituted. A job a program cannot do is disabled with the
reason stated, never quietly run a different way. A method that needs an
auxiliary basis gets one, and the panel says it did, because that changes the
input file and a person reading the run directory should not have to work out
where the extra keyword came from.

Switching program resets the method, functional and basis to that program's
defaults, for the reason the optimizer panel already does it: the catalogues do
not translate, and ORCA's `M062X` is an input error to Behemoth.

## What is built, and what is not

The panel, every control on it, the capability table and the result routing are
built and tested. ORCA runs through it end to end, for all seven jobs.

**The other four programs are not wired through this panel yet.** They appear in
the selector, their catalogues and job support are correct, and choosing one and
pressing Run reports that it is not wired up and points at the panel that does
do the job. That is deliberate: each backend needs its own input construction
and output reading through this path, and quietly running the wrong thing would
be worse than saying so. They follow in the remaining stages.

## Where a run's results actually go

Each is routed from the job's own `result_targets`, so the panel and the job
description cannot disagree:

* an optimisation, saddle search or IRC loads the trajectory into the player,
  set to play through once, and opens the Trajectory window;
* a frequency job writes its Hessian where the Vibrations panel reads it, and
  opens that window with the modes loaded;
* a TD-DFT run's spectrum goes into the UV-Vis panel and opens it;
* every job, including a single point, opens the run summary.

## The run summary

Shows what the job actually did. A geometry job gets a cycle table with the
energy and the change per cycle in kJ/mol, plus the total change, with the last
cycle's self-consistent field history folded away underneath. A job that never
moved the nuclei gets the field history alone. Both end with the energetics.

The cycle table is read back out of the trajectory rather than parsed from the
program's log a second time, so it cannot disagree with the energy plot.

### One trap worth naming

After a TD-DFT run ORCA's `FINAL SINGLE POINT ENERGY` is **not** the ground
state: it is the self-consistent field energy plus the first root's excitation
energy. Verified on a real run, where the output states `E(SCF) = -76.321074728`,
`DE(CIS) = 0.279708561 (Root 1)` and `E(tot) = -76.041366167`, and calls the
last of those final. The summary shows the ground-state figure as the answer and
says in as many words what ORCA's own final number is, because quoting it would
be wrong by the whole excitation energy and would look entirely plausible.

## Testing

Command and input construction is tested without running anything, as the
existing backends are. The job-support table above is asserted in a test, so a
program gaining or losing a capability has one place to change. The
end-to-end tests that run the real programs are extended to cover a
transition-state optimisation and a TD-DFT run through ORCA.
