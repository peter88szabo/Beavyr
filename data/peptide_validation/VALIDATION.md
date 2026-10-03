# Peptide builder validation

Test sequence: `ACDEFGHIKLMNPQRSTVWY` (20 standard L-amino acids; glycine achiral).
This is an extended starting peptide, not insulin or a predicted protein fold.

- 168 heavy atoms; all 159 hydrogens removed and restored; 327 final atoms.
- Hydrogen restoration leaves every heavy-atom coordinate unchanged and is idempotent.
- Each restored H has one explicit covalent neighbor.
- Template-based hydrogen restoration also passed for all 62 embedded residue forms.
- Cleanup/optimization retain explicit bond orders; optimized bond-length ratios remain
  within the test's 0.65–1.45 range, and retained residue stereocenters do not invert.
- Full regression suite: 1,292 passed, 0 failed, 23 deliberately ignored.
- Explicit 20-residue benchmark: passed, including the continuation and stereochemistry checks.
- CPU affinity: cores 0–3; Cargo jobs: 4; benchmark test threads: 1.

| Stage | Cycles | Initial energy (kcal/mol) | Final energy (kcal/mol) | Largest force component (kcal/mol/Å) | Converged |
|---|---:|---:|---:|---:|---|
| DREIDING cleanup | 400 | 572.091299 | 273.648584 | 0.585966 | No |
| DREIDING optimization | 1000 | 273.643352 | 265.221078 | 0.149836 | No |
| Fresh BFGS continuation | 1 | 265.222934 | 265.222248 | 0.106717 | Yes |

Each stage reports its own initial energy after initializing its force-field workspace;
small differences between adjacent stage endpoints are shown without rounding them away.
The final optimization targets were maximum gradient 1e-4 and RMS gradient 5e-5
Hartree/bohr, with the existing energy and step convergence criteria.

The BFGS matrix update now uses the equivalent quadratic-cost rank-two formula.
An independent dense-formula comparison and secant-condition test pass. DREIDING
optimization now permits one bounded restart (at most two 1000-step passes),
records both passes, and labels output as unconverged if the criteria are unmet.

The optimized coordinates are in [peptide20-dreiding.xyz](peptide20-dreiding.xyz).
DREIDING convergence is convergence within this generic force-field model; it does
not establish a global conformational minimum or a biological fold.

Reproduce the chemistry and optimizer benchmark:

```sh
taskset -c 0-3 cargo test -j4 peptide_twenty_residue_hydrogen_and_dreiding_benchmark -- --ignored --nocapture --test-threads=1
```

The benchmark writes `/tmp/beavyr-peptide-20-dreiding.xyz`. Its continuation has a
4000-step test ceiling to diagnose iteration limits; the tested run required one
step. The application uses a stricter 1000-step continuation ceiling.
