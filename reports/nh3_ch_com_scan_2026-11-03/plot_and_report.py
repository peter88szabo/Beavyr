#!/usr/bin/env python3
"""Plot the completed scan and retain numerical comparisons in kcal/mol and kJ/mol."""
import json
import math
import csv
import re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
rows = sorted(json.loads((ROOT/'points.json').read_text()), key=lambda r:r['distance_angstrom'])
jobs = json.loads((ROOT/'approved_jobs.json').read_text())
assert len(rows) == 8 and len(jobs) == 24
assert all(r['status'] == 'PASS' for r in jobs)
assert rows[-1]['distance_angstrom'] == 10
for key in ('composite_relative_kcal_mol', 'cc_atz_relative_kcal_mol'):
    assert rows[-1][key] == 0
distances = [r['distance_angstrom'] for r in rows]
composite = [r['composite_relative_kcal_mol'] for r in rows]
conventional = [r['cc_atz_relative_kcal_mol'] for r in rows]
differences = [r['composite_minus_cc_atz_kcal_mol'] for r in rows]
plt.rcParams.update({'font.size':11, 'axes.spines.top':False, 'axes.spines.right':False, 'savefig.facecolor':'white'})
fig = plt.figure(figsize=(11.2,9.1), layout='constrained')
grid = fig.add_gridspec(2,2,height_ratios=(1.15,1))
full = fig.add_subplot(grid[0,:])
tail = fig.add_subplot(grid[1,0])
delta = fig.add_subplot(grid[1,1])
for ax,start in ((full,0),(tail,1)):
    ax.plot(distances[start:],conventional[start:],'o-',color='#1769aa',lw=2,ms=6,label='UHF-UCCSD(T)/aTZ (Behemoth)')
    ax.plot(distances[start:],composite[start:],'s--',color='#d66b16',lw=2,ms=6,mfc='white',label='Composite: CCSD(T)-F12/aDZ + ΔMP2-F12')
    ax.axhline(0,color='#777777',lw=.8)
    ax.set_xlabel('Fragment COM distance (Å)')
    ax.set_ylabel('E(R) − E(10 Å) (kcal/mol)')
    ax.grid(alpha=.18)
full.set_title('Full rigid scan: 8 equally spaced distances',loc='left',fontsize=12)
full.legend(loc='best',frameon=False,fontsize=10)
tail.set_title('Longer-range points',loc='left',fontsize=12)
delta.plot(distances,differences,'D-',color='#673ab7',lw=1.7,ms=5)
delta.axhline(0,color='#777777',lw=.8)
delta.grid(alpha=.18)
delta.set_xlabel('Fragment COM distance (Å)')
delta.set_ylabel('Composite − conventional (kcal/mol)')
delta.set_title('Difference between relative-energy curves',loc='left',fontsize=12)
fig.suptitle('NH₃ + CH: composite versus conventional CCSD(T)/aug-cc-pVTZ\nUHF doublet · frozen core · each curve independently zeroed at 10 Å',fontsize=14)
for extension in ('png','svg','pdf'):
    fig.savefig(ROOT/f'comparison_kcal_mol.{extension}',dpi=190)
plt.close(fig)

max_row = max(rows,key=lambda r:abs(r['composite_minus_cc_atz_kcal_mol']))
rms7 = math.sqrt(sum(d*d for d in differences[:-1])/7)
tailmax = max(abs(d) for d in differences[1:-1])
metadata = json.loads((ROOT/'metadata.json').read_text())
lines = ['NH3 + CH rigid COM scan: composite versus conventional UCCSD(T)/aTZ',
         '===================================================================','',
         'Completed 8 distances: 16 MRCC component jobs and 8 Behemoth conventional jobs, one at a time.',
         'Distance range: 2.0 to 10.0 Angstrom, equally spaced by 8/7 Angstrom.',
         'The 10.0 Angstrom reference was computed first; tables/plots are sorted',
         'by increasing distance. All molecular internal coordinates were fixed.',
         'The original NH3 geometry was supplied by the user. CH bond: 1.12 A.',
         'NH3 COM at origin; CH COM at (0,0,R), C toward NH3 and H away.',
         'Masses for COM: N=14.007, C=12.011, H=1.008 amu.',
         '', 'Electronic structure:',
         'MRCC 26.1.3, executable: '+metadata['executable'],
         'Neutral doublet; UHF reference; frozen C and N 1s core; symmetry off.',
         '8 CPU IDs (affinity 0-7), OMP/MKL/OpenBLAS/Rayon thread variables=8.',
         'Memory allowance: 16000 MB per calculation.',
         'MRCC allows 200 SCF/100 CC iterations. Behemoth input files record the',
         'SCF iteration allowance and mixer for each point. At 3.142857 A only,',
         'the density tolerance was adjusted from 1e-8 to 3e-8 after 400 and',
         '2000-cycle attempts stalled near 2e-8. The energy threshold remains',
         '1e-10 Eh and orbital-gradient threshold 1e-9. The SCF energy drift',
         'during the long failed attempt was below 0.00002 kcal/mol.',
         'Other convergence thresholds retain defaults. SCF fallback is',
         'disabled after the excluded 2 A',
         'attempt. Failed SCF attempts are archived and excluded from energies',
         'and selected-job timings. Fresh directories were used for each retry.',
         'Conventional curve: Behemoth release binary, method = UHF CCSD(T),',
         'spherical aug-cc-pVTZ, chemical frozen core (2 orbitals per spin),',
         'Cholesky threshold 1e-10, CC energy/amplitude tolerances 1e-8/1e-6.',
         'The entire conventional curve uses Behemoth, including its 10 A zero.',
         'The 10 A result agrees with MRCC within 0.000453 kcal/mol.',
         'At 2 A, an initial Behemoth SCF fallback found a different UHF solution',
         'with <S^2>=1.7015327111. That attempt is explicitly excluded. Disabling',
         'fallback preserves the UHF solution corresponding to the MRCC curve.',
         'Each subsequent Behemoth reference <S^2> must match the MRCC aTZ',
         'component within 0.001 before its energy is accepted.',
         'User authorization for 8 cores supersedes the earlier four-core limit.',
         '', 'Definition (aDZ=aug-cc-pVDZ, aTZ=aug-cc-pVTZ):',
         'E_comp(R) = E_CCSD(F12*)(T)/aDZ(R)',
         '          + E_MP2-F12/aTZ(R) - E_MP2-F12/aDZ(R).',
         'E_ref(R)  = E_conventional_CCSD(T)/aTZ(R).',
         'DeltaE_method(R) = E_method(R) - E_method(10 A), separately for each curve.',
         'MP2-F12 total energies include CABS singles; fitting/CABS choices are',
         'the MRCC defaults. MP2-F12/aDZ is taken from the CCSD(T)-F12/aDZ output,',
         'avoiding a duplicate calculation. The standard (T) energy is used;',
         'the separately printed scaled (T*)/(T+) alternatives are not substituted.',
         '', 'No counterpoise correction or geometry relaxation was performed.',
         'The 10 A reference is the requested finite-distance reference, not an',
         'independently computed sum of monomer energies at infinite separation.',
         '', 'Relative energies (kcal/mol)',
         'COM/A       Composite      CCSD(T)/aTZ       Difference']
for r in rows:
    lines.append(f"{r['distance_angstrom']:8.5f}    {r['composite_relative_kcal_mol']:12.7f}    {r['cc_atz_relative_kcal_mol']:12.7f}    {r['composite_minus_cc_atz_kcal_mol']:12.7f}")
lines += ['', 'Comparison statistics:',
          f"Maximum absolute curve difference: {abs(max_row['composite_minus_cc_atz_kcal_mol']):.7f} kcal/mol at {max_row['distance_angstrom']:.7f} A.",
          f'RMS difference over the 7 nonreference points: {rms7:.7f} kcal/mol.',
          f'Maximum absolute difference for the 6 points above 2 A and below 10 A: {tailmax:.7f} kcal/mol.',
          '', 'Timing and memory:',
          f"Total elapsed calculation wall time: {sum(j['wall_seconds'] for j in jobs):.2f} seconds.",
          f"Composite components total: {sum(r['composite_wall_seconds'] for r in rows):.2f} seconds.",
          f"Conventional CCSD(T)/aTZ total: {sum(r['cc_atz_wall_seconds'] for r in rows):.2f} seconds.",
          f"Largest measured process RSS: {max(j['peak_rss_kib'] for j in jobs)/1024**2:.3f} GiB.",
          'Times are GNU time wall times for entire engine invocations. RSS is the',
          'maximum measured process/child resident size, not a sum of simultaneous',
          'process memory. These are single measurements, not repeated timing tests.',
          '', 'Spin diagnostics as printed by the calculation engines:']
for j in jobs:
    lines.append(f"R={j['distance_angstrom']:.7f} A, {j['job']}: {j['spin_diagnostics']}")
lines += ['', 'Interpretation limits:',
          'This compares two approximate finite-basis methods; conventional aTZ',
          'is not a complete-basis reference. Differences do not establish which',
          'curve is more accurate. Eight rigid geometries at one orientation do',
          'not establish transferability to other orientations, reactions, or',
          'optimized stationary points. Lines in the figure connect actual points;',
          'no fitted minimum, barrier, or interpolated calculation is claimed.',
          'MRCC default cctol=6 was retained for F12. Printed decimal places preserve',
          'the numerical output, not a claim of equivalent accuracy. Differences',
          'around 0.001 kcal/mol would require tighter convergence checks.',
          '', 'Artifacts:',
          'scan_data.csv: total energies, components, timings, kcal/mol and kJ/mol.',
          'points.json: the same complete point records.',
          'approved_jobs.csv/approved_jobs.json: selected per-calculation results.',
          'Original jobs.csv/jobs.json: initial MRCC jobs retained as separate evidence.',
          'reference_crosscheck.json: validation against MRCC at 10 A.',
          'spin_data.csv: UHF reference <S^2> and excess above 0.75 for every component.',
          'spin_contamination.png/.svg/.pdf: spin diagnostics along the scan.',
          'comparison_kcal_mol.png/.svg/.pdf: plotted comparison.',
          'scan_geometries.xyz: all 8 geometries, increasing COM distance.',
          'runs/: exact MINP (MRCC) or input.inp (Behemoth), geometry.xyz,',
          'output.out, error.err and time.txt. The Behemoth jobs run input.inp.',
          'metadata.json: geometry source, settings, executable hash, scratch path.',
          'run_approved_scan.py and plot_and_report.py: approved calculation/analysis scripts.',
          'run_scan.py and initial_mrcc_* retain the original MRCC-only attempt.',
          'No application source files, installation settings or git commits changed.']

# Spin contamination is a reference-wavefunction diagnostic, not a composite property.
spin_rows = []
for j in jobs:
    log = (Path(j['evidence'])/'output.out').read_text(errors='replace')
    expression = r'Actual <S\^2>:\s*([-+0-9.Ee]+)' if j['engine'] == 'Behemoth' else r'Expectation value of S\^2:\s*([-+0-9.Ee]+)'
    match = re.search(expression, log)
    assert match, j['evidence']
    s2 = float(match[1])
    spin_rows.append(dict(distance_angstrom=j['distance_angstrom'], job=j['job'], engine=j['engine'], basis=j['basis'], reference='UHF', s2=s2, target_s2=0.75, excess_s2=s2-0.75, source=j['evidence']))
spin_rows.sort(key=lambda r:(r['distance_angstrom'], r['job']))
with (ROOT/'spin_data.csv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=list(spin_rows[0])); writer.writeheader(); writer.writerows(spin_rows)
fig, axes = plt.subplots(2,1,figsize=(10.5,7.5),sharex=True,layout='constrained')
for job,label,color,marker in [('cc_f12_adz','MRCC CCSD(T)-F12/aDZ reference','#d66b16','s'),('mp2_f12_atz','MRCC MP2-F12/aTZ reference','#16867a','^'),('cc_atz','Behemoth UCCSD(T)/aTZ reference','#1769aa','o')]:
    group = [r for r in spin_rows if r['job'] == job]
    for ax,key in zip(axes,('s2','excess_s2')):
        ax.plot([r['distance_angstrom'] for r in group],[r[key] for r in group],marker=marker,color=color,lw=1.8,ms=6,mfc='white',label=label)
axes[0].axhline(.75,color='#555555',ls=':',label='Pure doublet: 0.75')
axes[1].axhline(0,color='#555555',ls=':')
axes[0].set_ylabel('Reference <S^2>')
axes[1].set_ylabel('Spin contamination: <S^2> - 0.75')
axes[1].set_xlabel('Fragment COM distance (Å)')
axes[0].legend(frameon=False,fontsize=9)
for ax in axes:ax.grid(alpha=.18)
fig.suptitle('NH3 + CH: UHF reference spin contamination along the COM scan',fontsize=13)
for extension in ('png','svg','pdf'):fig.savefig(ROOT/f'spin_contamination.{extension}',dpi=190)
plt.close(fig)
lines += ['', 'UHF reference spin contamination (<S^2> minus 0.75):', 'COM/A       Component             <S^2>       Excess']
for r in spin_rows:lines.append(f"{r['distance_angstrom']:8.5f}    {r['job']:<18}  {r['s2']:10.7f}  {r['excess_s2']:10.7f}")
lines += ['The composite is an energy expression and has no single wavefunction;', 'its component reference diagnostics are therefore reported separately.', 'These are SCF reference diagnostics. MRCC projected correlated diagnostics', 'are preserved above but are not mixed into the reference plots.', 'MRCC prints reference <S^2> to four decimal places; CSV preserves that precision.']

lines += ['', 'Numerical controls of the accepted Behemoth jobs:']
for j in sorted((j for j in jobs if j['engine']=='Behemoth'), key=lambda j:j['distance_angstrom']):
    inp=(Path(j['evidence'])/'input.inp').read_text()
    lines.append(f"R={j['distance_angstrom']:.7f} A: " + ' '.join(inp.split()))
lines += ['', 'Excluded attempts (not included in the comparison or selected-job timings):']
for path in sorted(ROOT.glob('*retry*.json')):
    data=json.loads(path.read_text())
    if isinstance(data,dict): data=[data]
    for entry in data:
        lines.append(f"{path.name}: {entry.get('distance_angstrom', '?')} A, {entry.get('job', '?')}, {entry.get('wall_seconds', '?')} seconds, {entry.get('failure','see archived output')}")
lines.append('The excluded high-spin 2 A attempt is documented in excluded_high_spin_attempt.json.')
(ROOT/'report.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[lines.index('Relative energies (kcal/mol)'):lines.index('Timing and memory:')]))
print('Saved report and PNG/SVG/PDF figures to', ROOT)
