# DREIDING Force Field in Beavyr -- Design

Written 2026-09-09. A native DREIDING implementation for instant geometry
cleanup while building, designed so the same engine later drives molecular
dynamics.

Reference: S. L. Mayo, B. D. Olafson, W. A. Goddard III, "DREIDING: A Generic
Force Field for Molecular Simulations", *J. Phys. Chem.* **1990**, *94*,
8897-8909. Local copy:
`~/Dropbox/Research_Leuven/Tautomer/assets/forcefield/papers_to_implement/mayo1990.pdf`

Equation and table numbers below refer to that paper.

## Goal

Three uses, one engine:

1. **Now:** press a button in the Fragment Editor and have a sketched or
   fragment-assembled structure relax to a chemically sensible geometry in a
   fraction of a second, with no external program.
2. **Later:** molecular dynamics on the same energy and forces.
3. **Later:** conformational search, ranking many minimised structures.

Uses 2 and 3 are what shape the architecture, and they pull in the same
direction: the topology is built once and shared, while everything mutable
lives in a per-worker scratch object. Retrofitting that split onto a
cleanup-only implementation would mean rewriting the core, so it is built in
from the start even though only the minimiser ships first.

Conformational search adds one requirement MD does not: it runs **many
independent minimisations from one topology, in parallel**. An engine holding
its neighbour list behind `&mut self` would serialise every worker. That is the
single most expensive thing to retrofit here, and the reason for the
topology/workspace split below.

## Why DREIDING

Beavyr can already drive xTB and Behemoth, but both are external processes with
startup cost measured in seconds. Cleanup while building has to feel
instantaneous, which means in-process. DREIDING is the right fit: it is a
*generic* force field whose parameters are derived from element and
hybridisation rather than tabulated per functional group, so it covers the whole
main group plus common metals from a handful of rules, and it is simple enough
to implement correctly and verify against a published paper.

## Atom typing: vendored, not linked

Typing is the fiddly part of any force field, and getting it wrong is silent --
a mistyped atom gives a plausible-looking but wrong geometry.

We use `dreid-typer` (https://github.com/peter88szabo/dreid-typer, upstream
`caltechmsc/dreid-typer` v0.4.0), MIT licensed,
`Copyright (c) 2025 California Institute of Technology, Materials and Process
Simulation Center`, by Tony Kan and William A. Goddard III -- a co-author of the
DREIDING paper itself. It ships 12,779 lines of tests, including
`tests/cases/dreiding_paper.rs`, which tests against the molecules in the paper.

**It is copied into our source tree, not declared as a dependency.** Beavyr
must build from its own checkout with nothing to fetch or link.

### What it gives us

`assign_topology(&MolecularGraph) -> Result<MolecularTopology, TyperError>`
takes elements plus discrete bond orders and returns:

| field | contents |
|---|---|
| `atoms` | `{id, element, atom_type: String, hybridization}` |
| `bonds` | `{atom_ids, order}` with `Resonant` as a distinct order |
| `angles` | every `end1-center-end2` triple |
| `propers` | every proper dihedral |
| `impropers` | every inversion centre |

So it supplies not only the atom types but the **complete term lists**. That is
the whole topology enumeration -- the part that is tedious to write and easy to
get subtly wrong -- already done and already tested.

### Removing its dependencies

Upstream pulls `serde`, `toml` and `thiserror`. Of Beavyr's tree, only `toml`,
`toml_writer` and `serde_spanned` would be genuinely new (the rest already
arrive via Bevy), but we eliminate all three anyway, because the entanglement is
shallow and confined entirely to *rule loading*:

* `serde` is used **only** to deserialise typing rules from TOML, plus
  `Deserialize` shims for `Element` and `Hybridization` that merely delegate to
  their independent `FromStr`.
* `toml` is used **only** in `parse_rules` and one error variant.
* `thiserror` is six `#[derive(Error)]` uses across two files.

The perception engine -- `aromaticity`, `electrons`, `hybridization`,
`kekulize`, `resonance`, `rings`, roughly 4,600 lines and the valuable part --
uses none of them.

So:

1. The 408-line `resources/default.rules.toml` is translated into
   `typer/rules_data.rs`, a static Rust array of `Rule` literals, replacing the
   `include_str!` + `toml::from_str` path. `assign_topology_with_rules` stays
   public, taking `&[Rule]` built in Rust, so custom rules remain possible
   without a file format.
2. The `Deserialize` impls are dropped; `FromStr` is kept.
3. `thiserror` derives are replaced with hand-written `Display` and
   `std::error::Error` impls.

**Net new dependencies: zero.** Beavyr stays at five.

The TOML translation is the one step where a silent transcription error could
break typing, so it is done by a throwaway generator script (run once, in the
scratchpad, using `toml` there and *not* in Beavyr) that parses the file and
emits the Rust literals. Correctness is then *proved*, not assumed, by porting
the upstream test suite into our tree and running it: if the generated rules
differ from the TOML in any way that matters, `dreiding_paper.rs` fails.

Provenance, the upstream version, and this list of changes are recorded in
`src/forcefield/dreiding/typer/README.md`, with `LICENSE-MIT` alongside it.

## Parameter coverage, and a gap to be explicit about

The typer's default ruleset can emit **56** atom types. The paper parameterises
**39** of them (Table I, plus Ti and Ru from Table VI).

**17 types the typer can emit have no published DREIDING parameters:**

`Ag Au Ba Cd Co Cs Cu Hg K Li Mg Mn Ni Pd Pt Rb Sr Tc`

and two more that are only partly covered:

`S_2` and `S_R` -- van der Waals parameters exist (Table II is indexed by
*element*, so sulfur is covered), but Table I gives a bond radius and bond
angle for `S_3` only. This one bites in practice: **thiophene** types as `S_R`.

**Decision: never invent a parameter.** When a structure contains an
unparameterised type, cleanup refuses and says which atom and which type,
rather than silently substituting a guess and returning a wrong geometry. If
extending coverage is wanted later that is a separate, deliberate change --
adding estimated `S_R`/`S_2` radii would be the obvious first candidate, but it
is a decision to be made in the open, not a default.

## Architecture

The engine is pure numerics with no Bevy or egui types anywhere. That is what
makes it usable from a minimiser, from MD, and from tests.

```
src/forcefield/
  mod.rs              -- public entry points used by the UI
  dreiding/
    mod.rs            -- DreidingSystem: the engine
    params.rs         -- Tables I-V as static data, plus lookup
    terms.rs          -- energy and analytic gradient, one fn per term
    neighbors.rs      -- cell list + Verlet list for the nonbonded loop
    minimize.rs       -- L-BFGS driver (ships now)
    typer/            -- vendored dreid-typer
```

### Topology and workspace

The engine is split in two. This is the central design decision.

**`DreidingTopology`** is immutable and built once: atom types, every term with
its parameters already resolved to numbers, the exclusion list, and masses.
Shared as `Arc<DreidingTopology>` across any number of workers.

```rust
pub struct DreidingTopology {
    natoms: usize,
    masses_amu: Vec<f64>,
    atom_types: Vec<&'static str>,   // reported to the user for diagnosis
    bonds:      Vec<BondTerm>,       // (i, j, k, r_e)
    angles:     Vec<AngleTerm>,      // (i, j, k, c, cos_theta0, linear)
    torsions:   Vec<TorsionTerm>,    // (i, j, k, l, v, n, phi0)
    inversions: Vec<InversionTerm>,  // (i, j, k, l, c, psi0)
    vdw:        Vec<VdwParam>,       // per atom: r0, d0
    hbond:      Vec<HBondTerm>,
    exclusions: Exclusions,          // 1-2 and 1-3 pairs, sorted
    rotatable:  Vec<RotatableBond>,  // for conformational search
}
```

**`Workspace`** is per-worker and mutable: the Verlet neighbour list, the
reference positions it was built at, and the L-BFGS history. Cheap to create,
never shared.

```rust
pub struct Workspace {
    neighbors: NeighborList,
    built_at:  Vec<f64>,
    lbfgs:     LbfgsHistory,
}
```

Chemistry perception and term enumeration happen only in
`DreidingTopology::build`. Nothing in the hot path allocates or re-perceives.

```rust
impl DreidingTopology {
    pub fn build(mol: &Molecule) -> Result<Self, BuildError>;

    /// Writes forces (-dE/dx, kcal/mol/Å) into `forces`; returns the energy.
    /// `&self` -- so N workers share one topology and minimise in parallel.
    pub fn energy_and_forces(
        &self, pos: &[f64], work: &mut Workspace, forces: &mut [f64],
    ) -> f64;

    pub fn energy_breakdown(&self, pos: &[f64], work: &mut Workspace) -> EnergyBreakdown;
    pub fn masses_amu(&self) -> &[f64];
    pub fn atom_types(&self) -> &[&'static str];
}
```

Passing the workspace explicitly is what makes all three drivers fall out of
one engine: the minimiser holds one workspace, an MD run holds one, and a
conformational search holds one per thread over a shared `Arc`.

`EnergyBreakdown` (per-term totals) is reported in the UI after cleanup and is
the main debugging handle -- a wrong term shows up immediately as an absurd
component.

### What conformational search will need, and already has

Recorded here so the first version does not accidentally preclude it:

* **Rotatable bonds and their moving sets.** `DreidingTopology` stores
  `rotatable`, derived from the same bond and ring perception the typer already
  performs -- single, acyclic, both ends non-terminal. Beavyr's
  `src/molecule_builder/rotator.rs` already rotates a substituent about a bond
  and can supply the moving-atom partition.
* **Duplicate rejection.** `src/rmsd/kabsch.rs` already provides aligned RMSD.
* **Ranking and storage.** `src/trajectory.rs` already holds multi-frame
  structures, so a conformer ensemble needs no new container.

So search is a driver over `energy_and_forces` plus three pieces that exist.
Nothing beyond `rotatable` is added to the engine for it.

### Units

kcal/mol and Å, DREIDING's native units, with the MD conversions fixed here so
they are not re-derived later:

* acceleration [Å/fs²] = force [kcal/mol/Å] / mass [amu] × `4.184e-4`
* k_B = `0.0019872041` kcal/mol/K

### Non-periodic

Beavyr targets molecules, not crystals, so there is no periodic boundary
condition and no Ewald sum. The nonbonded loop is nevertheless written so a
minimum-image convention could be inserted in one place if that ever changes.

### No charges in the first version

DREIDING's electrostatics need partial charges, and the paper's own Table V
shows the hydrogen-bond parameters are *conditioned on the charge convention*
(D_hb of 4.0, 7.0 or 9.5 kcal/mol for experimental, Gasteiger, or zero
charges). Beavyr has no charge model. We therefore use the **no-charges**
convention throughout, which the paper supports directly: D_hb = 9.5 kcal/mol
with R_hb = 2.75 Å, reproducing the water dimer at -6.03 kcal/mol and 2.92 Å
against 6.13 and 2.94 experimentally.

This is a stated choice, not an oversight. Adding Gasteiger charges later means
switching D_hb to 7.0 to stay consistent, and that coupling is noted in
`params.rs` next to the constant.

## Energy terms

All six with analytic gradients. Numerical-derivative agreement is a test
requirement, not an afterthought.

**1. Bond stretch** (Eq 4a, 6-9). Harmonic:
`E = ½K(R − R_e)²`, with `R_e = R_I + R_J − 0.01 Å` and `K = 700n` (kcal/mol)/Å²
for bond order `n`. `Resonant` uses `n = 1.5`. The paper's Morse form (Eq 5a) is
deliberately not used: it is `DREIDING/M`, and the paper itself recommends the
harmonic form as the default precisely because a sketched starting geometry may
be far from equilibrium, which is our exact case.

**2. Angle bend** (Eq 10a, 12). Harmonic cosine:
`E = ½C[cosθ − cosθ⁰]²` with `C = K / sin²θ⁰` and `K = 100` (kcal/mol)/rad².
For linear centres (`θ⁰ = 180°`), Eq 10' instead: `E = K[1 + cosθ]`. The paper
prefers the cosine form because the plain harmonic in θ does not give zero slope
as θ approaches 180°.

**3. Torsion** (Eq 13, Table IV, cases (a)-(g)). `E = ½V[1 − cos n(φ − φ⁰)]`,
with `V`, `n`, `φ⁰` selected by the hybridisation of the two central atoms:
sp3-sp3 gives (2.0, 3, 180°); sp2-sp3 (1.0, 6, 0°); sp2=sp2 double (45, 2,
180°); resonant-resonant (25, 2, 180°); sp2-sp2 single (5, 2, 180°); the
exocyclic aromatic-aromatic exception (10, 2, 180°); and zero for sp1,
monovalent, or metal centres. Column-16 atoms with two single bonds (O_3, S_3,
Se3, Te3) take φ⁰ = 90°, which is why the paper discusses them as s²p⁴.

The barrier is divided by the number of torsions sharing the central bond, so
that e.g. ethane's nine H-C-C-H torsions sum to one 2.0 kcal/mol barrier.

**4. Inversion** (Eq 28, 29). `E = ½C[cosψ − cosψ⁰]²`, `K = 40`
(kcal/mol)/rad² for planar centres (X_2, X_R, ψ⁰ = 0°) and for C_31
(ψ⁰ = 54.74°); zero for X_3. Per the paper, all three permutations at a centre
are summed, each weighted 1/3.

**5. van der Waals** (Table II). Lennard-Jones 12-6, geometric-mean
combination: `R⁰_IJ = √(R⁰_I R⁰_J)`, `D⁰_IJ = √(D⁰_I D⁰_J)`. Parameters are
per *element*, not per type. 1-2 and 1-3 pairs excluded; 1-4 kept at full
strength, as DREIDING specifies no 1-4 scaling. The paper's exponential-6 (X6)
alternative is not implemented.

**6. Hydrogen bond** (Eq 25ish, Table V).
`E = D_hb[5(R_hb/R)¹² − 6(R_hb/R)¹⁰] cos⁴θ_DHA`, over `H_HB` donors and
N/O/F acceptors, with `D_hb = 9.5`, `R_hb = 2.75 Å` per the no-charges
convention above. Beavyr already has `src/hbonds.rs` for detection geometry to
draw on.

## Minimiser, and the rest of Behemoth's optimizer

**Not written here.** Behemoth's `src/optimizer/` is imported whole into
`src/optimizer/`: its `Objective` trait, Cartesian minimisation, redundant
internals, constrained optimisation, transition-state search, RDA with the
poor-man's nudged elastic band, IRC following, spin-crossing optimisation,
rigid scans and the genetic conformer search -- about 12,000 lines. A second,
weaker implementation in Beavyr would be duplicated maintenance and diverging
behaviour for no gain.

It arrives BLAS-free. Behemoth wrote the interface to be backend-agnostic --
"any backend that can provide `E(x)` and `∇E(x)` can be plugged in" -- and
nothing in the module touches `ndarray-linalg`. Four redirections cover
everything it reached for outside the optimizer:

| Behemoth | here |
|---|---|
| `numeric::linalg::{Backend, LinAlg}` | `normalmode::linalg_shim`, the pure-Rust Jacobi eigensolver already written for the imported frequency code |
| `numeric::linalg::solve_linear_system_array2` | `optimizer::linalg_solve` -- Behemoth's own Gauss-Jordan, which needed no BLAS |
| `utils::atomic_masses::AtomicMasses` | a forwarder onto Beavyr's mass table |
| `rand` | `optimizer::prng`, a seedable xoshiro256** with the same call surface |

Left behind: `xtb_gfn1` and `tasi_eht` (the quantum-chemistry engine, not
optimiser logic), `scan`'s two `*_method` wrappers that call it directly, and
the conformer search's one-thread Rayon pool, which exists so Rayon inside an
electronic-structure method cannot steal the other workers' threads.

`Objective` works in **bohr and Hartree** where DREIDING works in Å and
kcal/mol, so `dreiding::objective` converts at the boundary and neither side
knows about the other's units.

Cleanup uses looser convergence than a geometry optimisation, and both are
looser than Behemoth's quantum-chemistry defaults: there is no point converging
a generic force field's geometry beyond the accuracy of its own parameters.

## Bridging Beavyr to the typer

The typer wants elements and *discrete* bond orders. Beavyr stores
`bonds: Vec<(usize, usize, f32)>` -- index pairs and a distance -- and estimates
order from geometry in `src/bond_order.rs`, which returns `f64` and already
handles the cases that matter here (the peroxide O-O fix means it correctly
gives 1.0 rather than 1.5).

Mapping: `estimate_bond_order` → `GraphBondOrder` by
1.0 → `Single`, 1.5 → `Aromatic`, 2.0 → `Double`, 3.0 → `Triple`.

Element symbols are `String` in Beavyr and an `Element` enum in the typer, via
its `FromStr`. An unrecognised symbol is a `BuildError`, reported, never
skipped.

Because bond order is perceived from *geometry*, a badly sketched structure can
be mistyped -- a stretched double bond read as single. Cleanup therefore reports
the perceived types alongside the energy breakdown, so a surprising result is
diagnosable rather than mysterious.

## Testing

1. **Ported upstream suite** -- `dreiding_paper.rs`, `amino_acids.rs`,
   `nucleic_acids.rs`, `harness.rs`, `integration_tests.rs`, unchanged in
   substance. This is what preserves the "validated" property through
   vendoring, and what catches any error in the TOML-to-Rust rule translation.
2. **Analytic gradient vs central finite differences**, every term
   independently and all terms together, on molecules chosen to exercise each:
   agreement to 1e-6 relative.
3. **Known geometries** -- after minimisation, benzene stays planar with
   equal C-C bonds; ethane's C-C near 1.53 Å; water's angle near 104.5°;
   ammonia stays pyramidal (this is the same pyramidal-versus-planar concern as
   the recent `FRAG_NH2` fix, now checked by energy rather than by construction).
4. **Refusal** -- a ferrocene-like or thiophene input must fail with a message
   naming the atom and type, not produce a geometry.
5. **Energy conservation** -- once MD exists, a velocity-Verlet run on a small
   molecule conserves total energy over 10 ps. Written now as an ignored test so
   the requirement is recorded in the tree rather than in memory.

## What shipped

* **Fragment Editor "Clean up geometry"** -- the original goal. Runs
  synchronously, since a hand-built molecule relaxes in milliseconds.
* **The optimizer and frequency panels** offer DREIDING as a third program
  alongside xTB and Behemoth. `QcProgram::runs_in_process()` distinguishes a
  built-in backend, which needs no path, no scratch directory and nothing to
  kill on cancel. The optimisation writes its trajectory in Behemoth's
  comment-line format, so the energy plot, trajectory player and summary window
  work unchanged.
* **Frequencies** from central differences of the analytic gradient. **No
  infrared intensities**: an intensity is a dipole derivative, a dipole needs
  charges, and there is no charge model. Stated in the panel rather than
  fabricated.
* **A Conformer Search tab**, driving Behemoth's genetic algorithm over
  DREIDING. It asks only how hard to look -- rotors, population and generations
  all follow from the structure. Results load as a trajectory, save as
  multi-frame XYZ, and carry the citation for the algorithm.

## Still to do

* **Molecular dynamics** -- velocity Verlet, a thermostat, constraints. The
  topology/workspace split and the unit constants are already in place.
* **GFN-FF re-ranking of conformers.** Worth doing, since relative conformer
  energies are where a generic force field is weakest, and cheap -- a few
  hundred xTB calls on the survivors rather than the tens of thousands the
  search itself needs. Blocked on `xtbrun::call_xtb` writing fixed filenames
  into the current working directory, so it cannot be driven in a loop until it
  is scratch-directory aware.
* **Panels for the imported tools** -- transition states, IRC, RDA,
  spin-crossing and scans are all in the tree and unexposed.
* **Widening parameter coverage**, if wanted: estimated `S_R`/`S_2` radii would
  bring thiophene and thiones in, and Ti/Ru could be imported from DREIDING/A.
  Both are deliberate decisions rather than defaults.
