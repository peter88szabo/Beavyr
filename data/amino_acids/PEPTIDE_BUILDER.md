# Peptide Builder

Open **Molecular Editor → Peptide Builder…**. The floating window stays open
independently of the Molecular Editor. It uses the [embedded amino-acid
library](README.md); no external program is needed to build a chain.

## Build a new peptide

1. Choose **New peptide** and click the first amino-acid tile.
2. Choose where the next tile goes: **C terminus** appends a residue;
   **N terminus** prepends it. The displayed sequence always reads **N to C**.
3. Add residues, reorder with the arrows, or remove individual rows. Select a
   row to choose its side-chain protonation state or histidine tautomer.
4. Alternatively expand **Paste sequence** and enter one-letter codes, a single
   FASTA record, or separated three-letter names such as `Ala-Gly-Ser`.
   A bare `ALA` means Ala–Leu–Ala; the parser does not guess the intended format.
5. Choose **Extended**, **α-helix**, or custom φ/ψ starting angles and the free
   terminal state. Preview, drag the preview to rotate it, then **Insert peptide**.

The 22 tiles include Sec/U and Pyl/O. Sequences can contain up to 256 residues.
Changing the recipe or displayed structure invalidates the preview, so an old
preview cannot silently be inserted with different settings.

A new peptide is placed beside any existing structure, without replacing it.
An insertion is one Molecular Editor Undo/Redo step. The newly placed atoms are
automatically assigned to the existing **Newly Placed Fragment** controls;
**Fragment Editor ↗** transfers that target to the unchanged Fragment Editor.

## Extend an existing structure

Choose **Extend selected terminus**, then the end of the **host** to extend:

| Host end | Atom to select in the viewport | Connection | Displayed sequence |
| --- | --- | --- | --- |
| C | Carbon of a free COOH or COO⁻ group | Host C to new sequence N | Host — new sequence |
| N | Free primary/secondary amine nitrogen, with explicit H atoms | New sequence C to host N | New sequence — Host |

Click **Use selected atom** to bind the target, choose residues, preview, and
insert. The window reports how many leaving atoms will be removed. Existing
amides, esters, unsupported atoms and ambiguous attachment sites are rejected.
Changing the host geometry invalidates a bound target. The retained host atoms
stay at exactly their original Cartesian coordinates; freezes and isolation are
remapped to retained atom indices.

This also allows attachment to suitable free amino/carboxyl groups in a loaded
structure that lacks residue annotations. Selecting such a side-chain group
will deliberately build a branch there; the builder does not infer biological
chain membership from an XYZ file.

## What geometry is produced

The builder uses the optimized monomers' explicit atoms and bonds. For every
peptide bond it removes the carboxyl leaving oxygen and its hydrogen, where
present, and removes the required amino hydrogens. Internal peptide nitrogens
are neutral. Proline's internal nitrogen retains no hydrogen.

Backbone construction uses internal coordinates, trans peptide bonds and proper
rigid rotations of each monomer. It preserves side-chain stereochemistry and
ring geometry. Carbonyl oxygens and remaining amide hydrogens are arranged in
the peptide plane. Proline's ring fixes its φ rather than permitting an
independent torsion that would tear the ring apart.

Starting inter-residue geometry uses C–N 1.34 Å, CA–C–N 114°, and C–N–CA 123°
(the proline angle follows its ring). These are approximate construction
parameters, consistent with the
[Biotite peptide assembly example](https://www.biotite-python.org/latest/examples/gallery/structure/protein/peptide_assembly.html).
N–CA and CA–C lengths and the N–CA–C angle come from each embedded monomer.
Presets use φ/ψ = −135°/+135° or −57°/−47°; ω is 180°.

The result is a **starting geometry**, not an optimized peptide or a protein
fold prediction. Close-contact reporting detects severe nonbonded overlaps;
it is not a comprehensive steric-energy assessment. Use the existing Fragment
Editor and geometry optimization tools to adjust and relax the chain.

## Residue editing and persistence

The **Edit residues** tab operates on the displayed, annotated peptide. Choose a
residue from the list or **From selected atom**. Replacement, side-chain fitting,
ACE (N-acetyl) / NME (C-methylamide) caps and disulfide linking first produce a
preview; **Apply preview** commits one Molecular Editor Undo/Redo step.

Replacement preserves existing N, CA, C, O and terminal-group coordinates and
IDs, and fits a template side chain about CA. Proline replacement is refused:
its ring couples the nitrogen and backbone geometry, so change it in the
sequence before building. Cross-linked side chains cannot be replaced or fitted.

**Fit side chains** samples six rotations at each eligible side-chain single
bond in two local passes and selects lower steric-overlap scores. Backbone atoms,
ring bonds and cross-links remain fixed. This is deterministic torsion sampling,
not a Dunbrack probability model, a global packing search or a force-field
minimization. Original xTB ring geometries remain intact. This older local fitting action is
retained alongside the new library-based rotamer controls.

Disulfides require two free cysteine sulfurs already 1.8–2.3 Å apart and acceptable
C–S–S angles. Linking removes thiol H atoms and neutralizes sulfur charges. It
never stretches a distant connection or folds a chain. Automatic cyclic-peptide
closure and folding are not implemented.

Molecules now retain atom names, stable IDs for retained atoms, residue and chain
membership, per-atom formal charges and explicit bond orders. These survive
coordinate transforms, peptide edits, Undo/Redo and `.beavyr` project save/load.
Unknown host charges remain unknown. The builder does not overwrite quantum
chemistry charge/multiplicity settings. Generic edits which invalidate atom
membership discard annotations rather than attach them to the wrong atoms.

**Export PDB…** writes named residues, chains, formal charges and explicit
`CONECT` connectivity. PDB is an interchange export, not a full-fidelity project
format: use `.beavyr` to preserve template states and bond orders. XYZ and Recent
Structures remain coordinate-only formats and do not retain residue annotations.

## Hydrogens and display

**Add missing H** restores annotated peptide hydrogens from the selected residue
forms, including amide N, side-chain N, histidine tautomers and charged variants.
Local heavy-atom frames follow rotated side chains; amide H uses the actual
peptide plane. Existing coordinates are unchanged. For unannotated structures,
the bulk tool supports saturated C, neutral O and neutral amine/amide N, while
ambiguous multiple-bond/protonation cases are skipped.

Annotated peptides do not display hydrogen bonds. A chain of **21 or more amino
acids** defaults to **Backbone Trace** (N–CA–C). Shorter peptides retain the current
representation. These rules use residue annotations; an arbitrary XYZ structure
does not provide reliable protein identities. Other representations remain
selectable by the user.

## Validation and implementation

Assembly lives in `src/molecule_builder/peptide.rs`; the window and detached
preview live in `peptide_builder.rs`. Existing fragment transforms, controls,
selection gizmos and colors are not restyled. Integration adds one editor action,
one launcher button and an independent window draw call.

Tests cover all embedded forms, atom counts after condensation, side-chain
charges, backbone torsions, peptide planes, stereochemistry, proline rings,
extension from either end, unchanged retained host coordinates, stale previews,
Undo/Redo and the existing Fragment Editor's transformations.

```sh
taskset -c 0-3 cargo test -j 4 peptide -- --test-threads=1
```

Additional regressions cover named-H restoration across all 62 template forms,
backbone-preserving residue replacement, caps, display thresholds and project
round trips. Run the explicit 20-residue DREIDING benchmark with:

```sh
taskset -c 0-3 cargo test -j 4 peptide_twenty_residue_hydrogen_and_dreiding_benchmark -- --ignored --nocapture --test-threads=1
```

The benchmark uses `ACDEFGHIKLMNPQRSTVWY`, removes all hydrogens, restores them,
checks the unchanged heavy atoms, then performs DREIDING cleanup and optimization.
It reports convergence and energy, validates retained bond lengths and writes
`/tmp/beavyr-peptide-20-dreiding.xyz`. It is not an insulin fold or a folding test.

See the [recorded 20-residue validation and optimized XYZ](../peptide_validation/VALIDATION.md).


## Open an experimental protein

Use **Structure → Open structure…** for `.pdb`, `.ent`, `.cif` or `.mmcif`, or pass
the filename to Beavyr. Import uses [pdbtbx](https://docs.rs/pdbtbx/0.12.0/pdbtbx/)
for coordinates. The first model is loaded. For alternate locations, one complete
residue alternate is selected by mean occupancy (alphabetical tie-break), sharing
blank-alternate atoms. The selected alternate, source chain, residue number,
insertion code and source residue name survive project saves and Undo/Redo.

Standard residues receive named template connectivity and bond orders. Peptide
bonds require adjacent residue identifiers, the same uninterrupted chain and a
plausible C–N distance. PDB `TER` and missing residue numbers prevent implicit
joining. PDB `SSBOND`, `LINK`, covalent `CONECT` and mmCIF covalent/disulfide
`struct_conn` records preserve explicit cross-links. Symmetry-related links to
unloaded atoms are reported and skipped. No biological assembly or crystal
symmetry expansion is performed.

Unknown ligands and waters remain in the molecule. Their internal bonds may be
distance-inferred; they are not given invented amino-acid identities. Missing
heavy atoms are reported, not reconstructed. Imported histidine/state aliases
are recognized, but coordinates alone do not determine protonation: unlabelled
residue states use neutral defaults, and unspecified formal charges remain
unknown. Review residue states before **Add missing H** and set the calculation's
charge explicitly as needed. This importer is not a pH/protonation predictor.

PDB export retains source residue numbers and insertion codes. Single-character
chain IDs are retained where unique; other chain IDs receive unique PDB letters.
Beavyr projects preserve the original mmCIF identifiers without PDB's limits.

## Backbone conformation

In **Edit residues**, choose a residue, then use **Backbone · φ / ψ / ω**:

- Set φ, ψ and incoming ω, then **Preview backbone**. Undefined terminal angles
  are disabled. The Ramachandran plot shows measured φ/ψ values for the structure;
  clicking sets this residue's proposed φ/ψ. α/β labels are reference positions,
  not residue-specific allowed-region contours or an outlier assessment.
- Choose **α-helix**, **β-strand** or **Type-I turn**, select the last residue with
  **through**, then **Preview range preset**. A type-I turn sets its two central
  residues. Presets change φ/ψ and preserve ω; they do not fold a protein.
- Inspect the detached preview and close contacts, then **Apply preview**.
  Molecular Editor **Undo/Redo** restores the complete edit in one step.

Rotations operate on the covalently connected side of each bond. Bond lengths
and local stereochemistry are preserved, but downstream residues outside a
selected range can move as a rigid continuation. Proline φ, ring closure and
cross-links that prevent an independent rotation are refused. Frozen atoms are
checked both when previewing and applying. The Fragment Editor is unchanged.

## Side-chain rotamers

**Find rotamers** queries the embedded Dunbrack 2010 library through the
[`dunbrack` Rust crate](https://docs.rs/dunbrack/0.1.0/dunbrack/) using the measured
backbone angles. Undefined terminal dimensions are averaged over the periodic
10° grid. The current conformation and up to twelve candidates are shown;
arrows cycle their previews. Each candidate rotates the existing side chain,
preserving its bond lengths, ring geometry and backbone coordinates.

Candidates are ranked by squared nonbonded overlap, then library probability.
The overlap score includes internal side-chain clashes and contacts with the
rest of the structure, excluding bonded and angle-connected pairs. It is a
geometric score, not a force-field energy or a hydrogen-bond score. The shown
prior is a library probability, not a predicted binding probability.

After applying a mutation, **Preview nearby repack** searches neighbors within
2–12 Å of the selected residue. It keeps that residue and all backbones fixed,
accepting only lower-overlap side chains in one sequential pass. It is a local
repair tool, not exhaustive global packing. Unsupported or locked neighbors are
skipped and counted in the preview message.

Gly/Ala have no independent side-chain χ; proline ring fitting and Sec/Pyl library
rotamers are not implemented. Disulfide-linked side chains cannot be rotated
independently. The existing template builder and fitting tools remain available.
All rotamer tables are compiled in; no runtime network or Python is required.

Protein regression fixtures include RCSB 1CRN (PDB and mmCIF) and 1UBQ. Run the
import/conformation regressions with:

```sh
taskset -c 0-3 cargo test -j 4 protein_ -- --test-threads=4
```
