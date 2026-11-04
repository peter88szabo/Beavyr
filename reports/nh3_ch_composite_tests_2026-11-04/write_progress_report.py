#!/usr/bin/env python3
"""Create a dated, explicit snapshot; no quantum calculations are launched."""
import csv
import datetime
import json
import math
from pathlib import Path
import textwrap

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'nh3_ch_com_scan_2026-11-03'
jobs = [j for j in json.loads((OLD/'approved_jobs.json').read_text()) if j['status'] == 'PASS']
groups = {i: {j['job']: j for j in jobs if j['point'] == i} for i in range(8)}
KCAL = 627.5094740631
rows = []
for i, g in groups.items():
    low, high = g['cc_f12_adz'], g['mp2_f12_atz']
    rows.append(dict(R=2+8*i/7, direct=low['energy_hartree'],
        composite=low['energy_hartree']+high['energy_hartree']-low['mp2_f12_hartree'],
        conventional=g.get('cc_atz', {}).get('energy_hartree')))
zero = rows[-1]
for r in rows:
    for name in ('direct', 'composite', 'conventional'):
        r[name+'_relative'] = (r[name]-zero[name])*KCAL if r[name] is not None else None
done = sum(r['conventional'] is not None for r in rows)
f12_time = sum(j['wall_seconds'] for j in jobs if j['job']=='cc_f12_adz')
composite_time = sum(j['wall_seconds'] for j in jobs if j['job']!='cc_atz')
stamp = '2026-11-04 session; progress snapshot'
text = f'''CH + NH3: COMPOSITE ELECTRONIC-STRUCTURE METHODS FOR A CAPTURE PES
Progress report -- {stamp}

PURPOSE AND STATUS
The intended application is a flexible machine-learning potential energy surface
(PES) for quasi-classical trajectory (QCT) capture rates at 10-100 K. The intended
training set contains approximately 8,000-12,000 geometries. The aim is capture
into the entrance complex, rather than the complete product-forming reaction.

The initial rigid scan contains eight COM separations from 2 to 10 Angstrom.
All eight MRCC composite energies are complete. {done}/8 conventional Behemoth
UHF-UCCSD(T)/aug-cc-pVTZ energies are accepted. The 4.285714 Angstrom point
failed CCSD convergence and is deliberately omitted, as accepted by the user.
The 5.428571 Angstrom retry was stopped on the user's instruction when the
scope was narrowed to the five distances below. Neither missing point is being
continued. No unconverged energy is used in these tables.

Three additional composite methods have been authorized at ONLY these distances:
3.142857, 6.571429, 7.714286, 8.857143 and 10.000000 Angstrom.
Their additional components are now being calculated. Results are pending
unless explicitly included in an updated comparison table.

COMMON SETUP
System: NH3 + CH, six atoms, neutral doublet. The user supplied the approximate
NH3 geometry; the CH bond is 1.12 Angstrom. Internal coordinates and mutual
orientation remain fixed along the initial scan. CH points along the scan axis,
with C toward NH3. COM masses: N=14.007, C=12.011, H=1.008 amu.
All accepted initial calculations use UHF references and freeze the C and N 1s
cores. Calculations run sequentially, with eight CPU IDs (0-7) and a 16000 MB
memory allowance. No counterpoise correction or geometry relaxation was applied.
Each curve is independently referenced to its OWN energy at 10 Angstrom.
That finite-distance reference is not an independently evaluated monomer sum.

APPLIED METHODS
Here aDZ = aug-cc-pVDZ and aTZ = aug-cc-pVTZ.

A. Direct explicitly correlated baseline, MRCC:
   E_direct = E_CCSD(F12*)(T)/aDZ
The standard, unscaled (T) result printed by MRCC is selected. The separately
printed (T*) and (T+) alternatives are not substituted. The input keyword
CCSD(T)-F12 must therefore be distinguished from the particular output selected.

B. Original composite, MRCC:
   E_comp = E_CCSD(F12*)(T)/aDZ
          + E_MP2-F12/aTZ - E_MP2-F12/aDZ
The MP2-F12/aDZ term is available within the CCSD(T)-F12/aDZ calculation, so the
composite requires two engine runs per geometry. MP2-F12 total energies include
the CABS singles contribution. This is a finite-basis additive correction, not
a demonstrated CCSD(T) complete-basis-set result.

C. Conventional comparison, Behemoth:
   E_conv = E_UHF-UCCSD(T)/aTZ
The entire accepted conventional curve uses Behemoth, including its 10 Angstrom
zero. Its exact-integral route uses Cholesky factors with threshold 1e-10.
MRCC conventional cross-checks agree with Behemoth to approximately 0.000453
kcal/mol at 10 Angstrom and 0.0000374 kcal/mol at 2 Angstrom.

RESULTS OBTAINED SO FAR
Absolute total energies, Hartree. A dash means no accepted result.
R / A       Direct F12/aDZ       Original composite     Conventional aTZ
'''
for r in rows:
    c=f"{r['conventional']:.12f}" if r['conventional'] is not None else '--'
    text += f"{r['R']:8.6f}    {r['direct']:.12f}    {r['composite']:.12f}    {c}\n"
text += '\nRelative energies, kcal/mol; each method zeroed separately at 10 Angstrom.\nR / A        Direct F12/aDZ      Original composite      Conventional aTZ\n'
for r in rows:
    c=f"{r['conventional_relative']:.7f}" if r['conventional'] is not None else '--'
    text += f"{r['R']:8.6f}      {r['direct_relative']:12.7f}          {r['composite_relative']:12.7f}          {c}\n"
text += f'''
WHAT THESE RESULTS SHOW
At 2 Angstrom the original composite differs from conventional aTZ by 0.96115
kcal/mol, versus 1.12307 kcal/mol for direct F12/aDZ. At 3.142857 Angstrom,
the differences are 0.12652 and 0.02549 kcal/mol, respectively.
At 6.571429, 7.714286 and 8.857143 Angstrom the composite differences are
0.00709, 0.00194 and 0.00116 kcal/mol. Direct F12/aDZ is closer to conventional
aTZ at those three distances too. Thus the additive MP2 basis correction does
not consistently improve agreement with this particular comparison curve.

This is NOT evidence that direct F12/aDZ is necessarily more accurate than the
composite. Conventional aTZ also has basis incompleteness; closer agreement can
reflect cancellation of errors. Tiny differences need tighter numerical checks.
The present orientation is repulsive at all calculated finite separations relative
to 10 Angstrom. It does not validate attractive entrance-channel orientations.

Measured successful-job times (single measurements, seconds):
  Eight direct F12/aDZ jobs: {f12_time:.2f}; mean {f12_time/8:.2f} per geometry.
  Sixteen original-composite component jobs: {composite_time:.2f};
  mean {composite_time/8:.2f} per complete composite geometry.
  Conventional accepted points: listed in approved_jobs.csv.
  New composites report the sum of required component wall times, including
  original timings for reused components; newly incurred time is separate.
  Linear projection for 8,000 composite points: {composite_time/8*8000:.0f} seconds.
  Linear projection for 12,000 composite points: {composite_time/8*12000:.0f} seconds.
These projections exclude retries, force calculations, extra electronic states,
validation, and changes in convergence difficulty over the training domain.

CONVERGENCE AND SPIN DIAGNOSTICS
An initial Behemoth calculation at 2 Angstrom found a different UHF solution,
with <S^2>=1.7015327111. It is excluded. Disabling SCF fallback recovered the
reference matching MRCC, with <S^2>=0.7761074171. For a pure doublet, <S^2>=0.75.
Every accepted Behemoth reference was checked against MRCC's aTZ reference
<S^2>, with a tolerance of 0.001. This is a useful guard, not a complete proof
of electronic-state identity or wavefunction reliability.

At 3.142857 Angstrom, Behemoth's SCF density threshold was changed from 1e-8
to 3e-8 after long stalled attempts; its energy threshold remained 1e-10 Eh.
The long retry changed the SCF energy by less than 0.00002 kcal/mol. The accepted
CCSD calculation still had to satisfy its energy and amplitude thresholds.
CCSD convergence difficulties at 4.285714 and 5.428571 Angstrom are reported
explicitly; stable-looking energies alone are not accepted as convergence.

PROPOSED ADDITIONAL COMPOSITES: FIVE-POINT TEST
1. Conventional small-basis coupled-cluster correction:
   E_1 = E_MP2-F12/aTZ + [E_CCSD(T)/aDZ - E_MP2/aDZ].
This substitutes a conventional small-basis higher-order correlation correction
for the explicitly correlated correction. It may reduce cost, but its timing
and accuracy for this open-shell system must be measured, not assumed.

2. MP2-F12/CBS plus the existing small-basis F12 correction:
   E_2 = E_MP2-F12/CBS
       + [E_CCSD(F12*)(T)/aDZ - E_MP2-F12/aDZ].
Retaining the same small-basis correction isolates the effect of improving the
MP2 basis-limit estimate. The planned MP2-F12 basis pair is cc-pVTZ-F12 and
cc-pVQZ-F12. The HF+CABS contribution is taken from the larger basis; the
correlation part is extrapolated with a documented F12-appropriate exponent.
The correlation exponent is 4.071516; all auxiliary bases and geminal
exponents are recorded with the actual inputs. This is a CBS estimate, not
a guarantee that residual basis error has vanished.
This construction follows the established additive MP2-F12/high-level-correction
strategy for noncovalent benchmarks [1,2]; its suitability for CH+NH3 remains
to be demonstrated.

3. junChS-F12-style MRCC adaptation (explicitly approved by the user):
   E_3 = E_CCSD(F12*)(T)/jun-cc-pVTZ
       + [E_MP2-F12/CBS(jun-T,Q) - E_MP2-F12/jun-cc-pVTZ]
       + [E_MP2-F12(all electrons)/cc-pwCVTZ
          - E_MP2-F12(frozen core)/cc-pwCVTZ].
The published junChS-F12 construction uses the jun-T/Q basis pair and an n^-3
extrapolation, with a core-valence correction at MP2-F12/cc-pwCVTZ [3].
Its published calculations use Molpro F12a/F12b and MP2-F12 approximations.
MRCC instead supplies CCSD(F12*) and its 2B/F+K MP2-F12 approximation.
Therefore this test is an MRCC ADAPTATION, not an exact reproduction of the
published named model. Unlike the other frozen-core candidates, it includes
an explicit core-valence correction. The geometries remain the scan geometries;
the published model's geometry-optimization protocol is not being invoked.

All new components will retain the UHF reference and unscaled (T) convention.
New components use scftol=9 and cctol=8. At the user's explicit request, all
existing CCSD(T)-F12/aVDZ and MP2-F12/aVTZ shared components are reused, retaining
their original thresholds and timings. One 10 Angstrom F12/aVDZ recheck had
already finished before this instruction was applied; it is archived separately
and is not substituted into the selected comparison. Full input files,
raw output, timing, basis provenance and component energies will be retained.

IMPLICATIONS FOR A 10-100 K CAPTURE PES
k_B*T spans approximately 0.0199-0.1987 kcal/mol. A generic 1 kcal/mol accuracy
target is inadequate for the weak long-range region at these temperatures.
A provisional long-range fitting target near 0.001 kcal/mol is reasonable to
test, but it is not a demonstrated reference-method accuracy or a guaranteed
rate error. The rate must be checked against increasingly accurate reference
data and independently fitted potentials.

Recommended PES construction: monomer potentials plus an interaction potential,
with an explicit orientation-dependent long-range electrostatic, induction and
dispersion contribution. Fit the remaining interaction smoothly, using monomer
energies evaluated at corresponding distorted geometries. Extend beyond the
present 10 Angstrom endpoint and test trajectory starting distances, e.g. 20,
30 and 50 Angstrom, until the predicted capture rate is insensitive to them.
Spurious long-range neural-network interactions have been documented [4].

Start with a smaller orientation/distance/vibration sample, then add points from
trajectory exploration and model disagreement. Validate attractive approaches,
repulsive walls, angular dependence, forces, and trajectory energy conservation.
Use suitable rotational populations and vibrational zero-point sampling.
Define an inner capture criterion and test its sensitivity, together with impact
parameter coverage and timestep. Assess CH's electronic-state and spin-orbit
structure before relying on a single-surface description at 10-100 K.
Published CH/CD+NH3 work identifies an attractive, barrierless approach toward
the nitrogen lone pair [5], which the present orientation has not validated.

REFERENCES AND PROVENANCE
[1] Some Observations on Counterpoise Corrections for Explicitly Correlated
    Calculations on Noncovalent Interactions. DOI: 10.1021/ct500513b.
    https://pubs.acs.org/doi/10.1021/ct500513b
[2] S66x8 noncovalent interactions revisited: new benchmark and performance of
    composite localized coupled-cluster methods. DOI: 10.1039/D2CP03938A.
    https://pubs.rsc.org/en/content/articlehtml/2022/cp/d2cp03938a
[3] Lupi, Alessandrini, Puzzarini and Barone, The junChS and junChS-F12 models.
    https://arxiv.org/abs/2110.03075
[4] Neural Network Potential Energy Surface for the low temperature Ring Polymer
    Molecular Dynamics of the H2CO + OH reaction.
    https://arxiv.org/abs/2102.02682
[5] Low-temperature gas-phase formation of methanimine, CD + NH3 study.
    https://www.osti.gov/servlets/purl/2278904

Local evidence:
{OLD}
approved_jobs.json/csv: selected successful jobs and timing.
progress_data.csv: all eight composite energies, available conventional energies.
progress_spin_data.csv: component SCF <S^2> and excess over 0.75.
scan_geometries.xyz: the eight geometries.
runs/: exact input, output, error and time records, including excluded attempts.
{ROOT}
New five-point composite test preparation and this dated report.

No Beavyr or Behemoth application source changes or git commits were made for
these benchmark calculations. This is a progress snapshot, not a final claim
that all proposed methods or the production capture PES have been validated.
'''
name='CH_NH3_composite_methods_progress_2026-11-04'
appendix=ROOT/'final_results_appendix.txt'
if appendix.exists():
    text=text.replace('Progress report --', 'Completed five-point comparison --')
    text=text.replace('Their additional components are now being calculated. Results are pending\nunless explicitly included in an updated comparison table.', 'All three alternatives are complete at the five requested distances. Their\nresults, timings and diagnostics are included in the completed-tests section below.')
    text=text.replace('The planned MP2-F12 basis pair is', 'The MP2-F12 basis pair is')
    text=text.replace('All new components will retain', 'All new components retain')
    text=text.replace('New five-point composite test preparation and this dated report.', 'Completed five-point composite tests and this report.')
    text=text.replace('This is a progress snapshot, not a final claim\nthat all proposed methods or the production capture PES have been validated.', 'The five-point method comparison is complete. It does not establish that the\nproduction capture PES or rate coefficients have been validated.')
    text=text.replace('REFERENCES AND PROVENANCE', appendix.read_text()+'\nREFERENCES AND PROVENANCE')
    name='CH_NH3_composite_methods_results_2026-11-04'

(ROOT/(name+'.txt')).write_text(text)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
wrapped=[]
for line in text.splitlines():
    wrapped.extend(textwrap.wrap(line, width=98, replace_whitespace=False,
        drop_whitespace=False) or [''])
page_lines=65
with PdfPages(ROOT/(name+'.pdf'), metadata={'CreationDate': None, 'Title': 'CH + NH3 composite methods: progress report'}) as pdf:
    total=math.ceil(len(wrapped)/page_lines)
    for page,start in enumerate(range(0,len(wrapped),page_lines),1):
        fig=plt.figure(figsize=(8.27,11.69))
        fig.text(.055,.955,'\n'.join(wrapped[start:start+page_lines]),va='top',
            family='DejaVu Sans Mono',fontsize=8.0,linespacing=1.35)
        fig.text(.5,.025,f'CH + NH3 composite methods | {page}/{total}',ha='center',fontsize=8)
        pdf.savefig(fig);plt.close(fig)
    if appendix.exists():
        for filename in ('composite_comparison.png','spin_contamination.png'):
            picture=ROOT/filename
            if picture.exists():
                fig,ax=plt.subplots(figsize=(11.69,8.27) if filename.startswith('composite') else (8.27,11.69))
                ax.imshow(plt.imread(picture));ax.axis('off');fig.subplots_adjust(left=.02,right=.98,bottom=.02,top=.98)
                pdf.savefig(fig);plt.close(fig)
print(ROOT/(name+'.txt'))
print(ROOT/(name+'.pdf'))
