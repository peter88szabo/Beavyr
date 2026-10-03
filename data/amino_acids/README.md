# Embedded amino-acid library

The library contains **62 optimized forms of 22 amino acids**: the 20 standard
residues plus selenocysteine (Sec/U) and pyrrolysine (Pyl/O). Coordinates and
chemical metadata are compiled into
[`amino_acids.rs`](../../src/molecule_builder/amino_acids.rs), selected through
Rust enums. Looking up a template requires no files, parsing, allocation, Python,
or xTB at runtime. The [Peptide Builder](PEPTIDE_BUILDER.md) now uses these
templates to assemble chains and extend compatible terminal groups.

## Included forms

Every side-chain state below has **both** neutral amino/carboxylic-acid termini
and zwitterionic termini (terminal N +1, carboxylate −1).

| Residues | Side-chain states | Templates |
| --- | --- | ---: |
| Gly, Ala, Val, Leu, Ile, Ser, Thr, Met, Pro, Asn, Gln, Phe, Trp, Pyl | Neutral | 28 |
| Asp, Glu, Cys, Tyr, Sec | Neutral; deprotonated (−1) | 20 |
| Lys, Arg | Neutral; protonated (+1) | 8 |
| His | Neutral ND1-H; neutral NE2-H; protonated (+1) | 6 |

`TerminalState::Zwitterionic` describes the **termini**, not the total molecular
charge. A charged side chain still contributes ±1. Proline's secondary amine
is handled explicitly. Hydrogen atoms, formal charges and total charge are
stored for each form. Histidine requires an explicit tautomer selection;
`SideChainState::Neutral` does not silently choose one.

These are complete free amino acids, not already linked peptide residues.
Assembly must remove the appropriate terminal atoms, create peptide bonds,
update hydrogens and charges, and relax the resulting chain. Backbone atom
indices support that process. The library does not yet cover
independently cationic/anionic termini, D amino acids, peptide caps, all possible
protonation microstates, or automatic pH-dependent state selection.

## Geometry preparation and verification

Starting atom names, connectivity, ideal coordinates and stereochemical labels
come from the [wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd).
The 22 exact source records are retained in `ccd/`; their URLs and SHA-256
digests are recorded in `manifest.json`. L stereochemistry is retained, including
the additional stereocentres in threonine, isoleucine and pyrrolysine. Glycine
is achiral. RDKit checks the specified stereocentres before and after optimization.

All 62 geometries were optimized with **xTB 6.7.1 (edcfbbe), GFN2-xTB,
ALPB water**, without geometric constraints. Each calculation uses the form's
explicit charge, a closed-shell singlet, very-tight optimization and a numerical
Hessian:

```sh
xtb input.xyz --gfn 2 --alpb water --uhf 0 --ohess verytight --acc 0.1 --parallel 4 --chrg CHARGE
```

See the xTB documentation for
[optimization](https://xtb-docs.readthedocs.io/en/latest/optimization.html) and
[solvation](https://xtb-docs.readthedocs.io/en/latest/gbsa.html).

An accepted geometry must have:

- Normal xTB termination and converged optimization, with final gradient norm
  ≤ 10⁻⁴ hartree/bohr.
- All 3N−6 vibrational frequencies positive; the six rigid modes are identified
  separately and must each have magnitude below 1 cm⁻¹.
- Finite coordinates, unchanged named atom order, plausible bonded distances,
  no unexpected close contacts, and each hydrogen closest to its assigned
  bonded heavy atom, rejecting proton transfers.
- The expected chemical graph, protonation state and stereochemistry.

The first attempt starts from CCD coordinates. Failed candidates are retried
with deterministic RDKit ETKDGv3 seeds; rejected attempts are recorded. These
are **verified local minima**, not a claim to the global lowest-energy conformer
or the dominant solution species. Energies of different compositions/charges
should not be used directly to rank protonation populations.

## Rust access

```rust
use crate::molecule_builder::amino_acids::{AminoAcid, SideChainState, TerminalState};

let alanine = AminoAcid::Alanine
    .template(TerminalState::Neutral, SideChainState::Neutral)
    .expect("supported alanine form");
let atoms = alanine.atoms; // static slice; coordinates in angstrom
let nitrogen = alanine.backbone.nitrogen;

let histidine = AminoAcid::Histidine
    .template(TerminalState::Zwitterionic, SideChainState::HistidineDelta)
    .expect("supported histidine form");
```

`AminoAcid::ALL`, `from_one_letter`, `from_three_letter` and `templates()` support
residue pickers. Unsupported combinations return `None`. Bonds include
explicit aromatic orders, and atom names allow connector lookup without guessing
indices. Each template carries its charge, formula, multiplicity, stereocentres,
optimization energy, gradient, lowest vibrational frequency and provenance hashes.

## Calculation records and reproduction

- `ccd/`: original source chemical components.
- `optimized/`: actual optimized xTB XYZ files.
- `logs/`: complete accepted xTB logs, compressed with gzip.
- `spectra/`: accepted xTB `vibspectrum` files, including rigid modes.
- `manifest.json`: source graphs, initial and optimized coordinates, all
  frequencies, rejected attempts, method arguments, dependency versions, xTB
  executable hash and artifact hashes.

Generation dependencies are isolated from Beavyr. Using Python 3.12 and the
pinned generation requirements:

```sh
taskset -c 0-3 python3.12 -m venv /tmp/beavyr-amino-venv
taskset -c 0-3 /tmp/beavyr-amino-venv/bin/pip install -r scripts/requirements-amino-acids.txt
taskset -c 0-3 env OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 OMP_STACKSIZE=256M \
  /tmp/beavyr-amino-venv/bin/python scripts/generate_amino_acids.py --xtb /absolute/path/to/xtb
```

Run from the repository root. Use up to four CPU IDs allowed by your machine if
0–3 are unavailable. The generator additionally limits its own affinity to four
cores and runs xTB calculations sequentially. Cached results are tied to the
source model, calculation settings, protocol version and xTB executable.
Use a fresh `--work /tmp/beavyr-amino-fresh` directory to recalculate everything.
Source files are downloaded only when absent. New calculations can differ
slightly across platforms; exact accepted records are retained here.

Audit the published data **without running xTB or accessing the network**:

```sh
taskset -c 0-3 env OPENBLAS_NUM_THREADS=1 \
  /tmp/beavyr-amino-venv/bin/python scripts/generate_amino_acids.py --verify-only
taskset -c 0-3 cargo test -j 4 amino_acids -- --test-threads=1
```

The offline audit checks source graphs, SHA-256 digests, optimization termination,
coordinates, stereochemistry, gradient norms and saved frequency spectra. Rust
tests independently check compositions, state labels, bonds, handedness and
agreement between compiled coordinates and saved xTB XYZ files. File access in
these tests is only for validation; the application uses the compiled arrays.
