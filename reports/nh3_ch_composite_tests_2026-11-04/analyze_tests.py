#!/usr/bin/env python3
"""Audit and plot the completed five-point alternative-composite comparison."""
import csv
import json
import math
from pathlib import Path
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'nh3_ch_com_scan_2026-11-03'
jobs=json.loads((ROOT/'jobs.json').read_text())
rows=sorted(json.loads((ROOT/'points.json').read_text()),key=lambda r:r['distance_angstrom'])
assert len(rows)==5 and len(jobs)==45
assert all(j['status']=='PASS' for j in jobs)
assert sum(not j['reused'] for j in jobs)==35
assert {r['point'] for r in rows}=={1,4,5,6,7}
previous=json.loads((OLD/'approved_jobs.json').read_text())
bh={j['point']:j for j in previous if j['job']=='cc_atz' and j['status']=='PASS'}
KCAL=627.5094740631
METHODS=[('direct_f12_adz','Direct F12/aDZ','#1769aa','o'),
 ('original_composite','Original composite','#d66b16','s'),
 ('conventional_small_basis_correction','Conventional small-basis correction','#16867a','^'),
 ('f12_cbs_composite','MP2-F12/CBS composite','#c43c39','D'),
 ('junchs_f12_mrcc','junChS-style MRCC adaptation','#8550a0','v'),
 ('conventional_atz','Conventional UCCSD(T)/aTZ','#333333','x')]

def parse(log,prefix):
    values=[float(l.split(':',1)[1].split()[0]) for s in log.splitlines() if (l:=' '.join(s.split())).startswith(prefix)]
    return values[-1] if values else None

# Independently check raw results rather than trusting the assembled CSV.
for j in jobs:
    log=(Path(j['evidence'])/'output.out').read_text()
    prefix='Total CCSD(T) energy [au]:' if j['job']=='cc_adz' else ('Total CCSD(F12*)(T) energy [au]:' if j['job'].startswith('cc_') else 'MP2-F12 energy [au]:')
    assert abs(parse(log,prefix)-j['energy_hartree'])<1e-11
    assert 'Normal termination of mrcc.' in log
    actual=float(re.search(r'Expectation value of S\^2:\s*([-+0-9.Ee]+)',log)[1])
    assert actual==j['s2']
    if j['job'].startswith('mp2_f12'):
        assert abs(j['reference_hartree']+j['cabs_hartree']+j['mp2_f12_correlation']-j['energy_hartree'])<5e-11
    if 'jun' in j['job']:
        first=int(re.search(r'Total number of basis functions:\s*(\d+)',log)[1])
        assert first==(134 if j['job']=='cc_f12_juntz' else 262),(j['job'],first)

for r in rows:
    r['conventional_atz']=bh[r['point']]['energy_hartree']
    r['conventional_atz_relative_kcal_mol']=(r['conventional_atz']-bh[7]['energy_hartree'])*KCAL
    r['core_valence_relative_kcal_mol']=(r['core_valence_correction']-rows[-1]['core_valence_correction'])*KCAL
    r['junchs_f12_mrcc_frozen_core_relative_kcal_mol']=r['junchs_f12_mrcc_relative_kcal_mol']-r['core_valence_relative_kcal_mol']
    g={j['job']:j for j in jobs if j['point']==r['point']}
    assert abs(g['mp2_f12_cv_fc']['reference_hartree']-g['mp2_f12_cv_ae']['reference_hartree'])<1e-10
    assert abs(r['junchs_f12_mrcc']-(g['cc_f12_juntz']['energy_hartree']+(64*g['mp2_f12_junqz']['energy_hartree']-27*g['cc_f12_juntz']['mp2_f12_hartree'])/37-g['cc_f12_juntz']['mp2_f12_hartree']+g['mp2_f12_cv_ae']['energy_hartree']-g['mp2_f12_cv_fc']['energy_hartree']))<1e-10
    for key,_,_,_ in METHODS:
        assert abs(r[key+'_relative_kcal_mol']-(r[key]-rows[-1][key])*KCAL)<1e-8 if key!='conventional_atz' else True
    for key,_,_,_ in METHODS[:-1]:
        r[key+'_difference_from_conventional_kcal_mol']=r[key+'_relative_kcal_mol']-r['conventional_atz_relative_kcal_mol']

with (ROOT/'comparison_complete.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
spin=[]
for j in jobs:
    spin.append({k:j[k] for k in ('distance_angstrom','job','basis','core','s2','excess_s2','reused','evidence')})
with (ROOT/'spin_diagnostics.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(spin[0]));w.writeheader();w.writerows(sorted(spin,key=lambda j:(j['distance_angstrom'],j['job'])))

x=[r['distance_angstrom'] for r in rows]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(12,9),layout='constrained')
grid=fig.add_gridspec(2,2)
full=fig.add_subplot(grid[0,:]);tail=fig.add_subplot(grid[1,0]);diff=fig.add_subplot(grid[1,1])
for key,label,color,marker in METHODS:
    y=[r[key+'_relative_kcal_mol'] for r in rows]
    for ax,start in ((full,0),(tail,1)):
        ax.plot(x[start:],y[start:],marker=marker,color=color,lw=1.5,ms=5,label=label)
    if key!='conventional_atz':
        diff.plot(x,[r[key+'_difference_from_conventional_kcal_mol'] for r in rows],marker=marker,color=color,lw=1.5,ms=5)
for ax in (full,tail,diff):
    ax.axhline(0,color='0.6',lw=.6);ax.grid(alpha=.18);ax.set_xlabel('COM distance (Å)')
for ax in (full,tail):ax.set_ylabel('E(R) − E(10 Å), kcal/mol')
diff.set_ylabel('Difference from conventional aTZ, kcal/mol')
full.set_title('Five requested distances',loc='left');tail.set_title('Long-range detail',loc='left');diff.set_title('Comparison, not an absolute error estimate',loc='left')
full.legend(ncol=2,frameon=False,fontsize=9)
fig.suptitle('CH + NH3: alternative composite methods\nUHF doublet · fixed geometry/orientation · independent 10 Å energy zeros',fontsize=14)
for ext in ('png','pdf','svg'):fig.savefig(ROOT/f'composite_comparison.{ext}',dpi=190)
plt.close(fig)

fig,axes=plt.subplots(3,1,figsize=(11,10),sharex=True,layout='constrained')
panels=[('Shared and conventional small-basis references',['cc_adz','cc_f12_adz','mp2_f12_atz']),
 ('MP2-F12/CBS reference components',['mp2_f12_tzf12','mp2_f12_qzf12']),
 ('junChS-style reference components',['cc_f12_juntz','mp2_f12_junqz','mp2_f12_cv_fc','mp2_f12_cv_ae'])]
for ax,(title,keys) in zip(axes,panels):
    for key in keys:
        group=sorted((j for j in jobs if j['job']==key),key=lambda j:j['distance_angstrom'])
        ax.plot([j['distance_angstrom'] for j in group],[j['excess_s2'] for j in group],'o-',ms=4,label=key)
    ax.set_title(title,loc='left');ax.set_ylabel('<S^2> - 0.75');ax.grid(alpha=.18);ax.legend(frameon=False,fontsize=8,ncol=2)
axes[-1].set_xlabel('COM distance (Å)')
fig.suptitle('UHF reference spin contamination\nComposite energies have no single associated wavefunction',fontsize=13)
for ext in ('png','pdf','svg'):fig.savefig(ROOT/f'spin_contamination.{ext}',dpi=180)
plt.close(fig)

lines=['COMPLETED FIVE-POINT ALTERNATIVE-COMPOSITE TESTS','',
 'All five requested geometries have results for all three alternative composites.',
 '35 new selected component calculations; 10 original component results reused.',
 'No additional Behemoth calculation was performed for this five-point test.',
 'The conventional aTZ comparison uses its five existing accepted energies.',
 '', 'Relative energies (kcal/mol), independent 10 Angstrom zero for every column.',
 'R/A       Direct F12     Original     Conv. delta     F12/CBS      jun/MRCC       Conv.aTZ']
for r in rows:
    lines.append(f"{r['distance_angstrom']:7.4f} "+' '.join(f"{r[key+'_relative_kcal_mol']:12.7f}" for key,_,_,_ in METHODS))
lines+=['','Absolute total energies (Hartree):']
for key,label,_,_ in METHODS:
    lines.append(label)
    for r in rows:lines.append(f"  {r['distance_angstrom']:9.6f} A   {r[key]:.12f}")
lines += ['', 'Effect of the jun adaptation core-valence correction on relative energies:', 'R/A          Delta(core-valence) / kcal/mol']
for r in rows:lines.append(f"{r['distance_angstrom']:9.6f}    {r['core_valence_relative_kcal_mol']: .8f}")
lines += ['The jun adaptation includes this correction; the other compared methods',
          'freeze the C and N 1s cores. The CSV also gives the jun frozen-core result.',
          'Agreement between the two basis-limit composites is not independent proof',
          'of accuracy, since both rely on MP2-F12-based additive corrections.']
lines+=['','Differences from conventional aTZ (four nonreference points):']
for key,label,_,_ in METHODS[:-1]:
    d=[r[key+'_difference_from_conventional_kcal_mol'] for r in rows[:-1]]
    lines.append(f'{label}: MAE={sum(abs(v) for v in d)/len(d):.7f}, max={max(abs(v) for v in d):.7f} kcal/mol.')
lines+=['These differences measure agreement with finite-basis conventional aTZ;',
 'they do not establish which model is closest to the exact or CBS interaction.',
 '', 'Complete-method component-time sums (seconds):',
 'R/A       Direct     Original    Conv. delta    F12/CBS    jun/MRCC']
for r in rows:
    lines.append(f"{r['distance_angstrom']:7.4f} "+' '.join(f"{r[k]:11.2f}" for k in ('direct_seconds','original_composite_seconds','conventional_correction_seconds','f12_cbs_seconds','junchs_seconds')))
lines += ['', 'Mean complete-method time over these five geometries (seconds):']
for label,key in [('Direct F12/aDZ','direct_seconds'),('Original composite','original_composite_seconds'),('Conventional small-basis correction','conventional_correction_seconds'),('MP2-F12/CBS composite','f12_cbs_seconds'),('junChS-style MRCC adaptation','junchs_seconds')]:
    lines.append(f'{label}: {sum(r[key] for r in rows)/len(rows):.2f}')
successful_new=sum(j['wall_seconds'] for j in jobs if not j['reused'])
restart_cost=sum(j.get('restart_preparation_seconds',0.0) for j in jobs if not j['reused'])
new=successful_new+restart_cost
lines += ['',f'New component cost including required restart preparation: {new:.2f} seconds.',
 f'Accepted engine invocations: {successful_new:.2f} seconds; restart preparation: {restart_cost:.2f} seconds.',
 f"Peak measured process/child RSS among selected new jobs: {max(j['peak_rss_kib'] for j in jobs if not j['reused'])/1024**2:.3f} GiB.",
 'Complete-method sums include original measured times of reused components.',
 'They estimate the cost of assembling a complete method; they are not the',
 'incremental time newly spent in this batch. Components shared by different',
 'methods must not be counted repeatedly in the batch total.',
 'New and reused components have different convergence thresholds, as documented.',
 'The 7.714286 A QZ-F12 SCF initially failed after 400 iterations. Its restart',
 'used the saved density, scfext=20 and scflshift=off, with unchanged thresholds.',
 'The unsuccessful preparation time is included in that composite timing.',
 'The already-started 10 A shared recheck (17.27 seconds) is excluded from selected',
 'method sums and is preserved in already_started_shared_recheck.json.',
 '', 'UHF reference spin diagnostic ranges across the five distances:']
for key in ('cc_adz','cc_f12_adz','mp2_f12_atz','mp2_f12_tzf12','mp2_f12_qzf12','cc_f12_juntz','mp2_f12_junqz','mp2_f12_cv_fc','mp2_f12_cv_ae'):
    vals=[j['s2'] for j in jobs if j['job']==key]
    lines.append(f'{key}: <S^2>={min(vals):.4f} to {max(vals):.4f}; excess above 0.75={min(vals)-.75:.4f} to {max(vals)-.75:.4f}.')
lines += ['MRCC prints these SCF diagnostics to four decimal places.',
 'The full component data are in spin_diagnostics.csv. These are reference',
 'diagnostics, not correlated CC spin expectation values.',
 '', 'Artifacts:',str(ROOT/'comparison_complete.csv'),str(ROOT/'composite_comparison.png'),
 str(ROOT/'spin_contamination.png'),str(ROOT/'jobs.csv'),str(ROOT/'protocol.json'),
 'runs/: exact inputs, output, timings and custom GENBAS where used.',
 'bases/: copied Psi4 source basis files, converted basis and source hashes.',
 'Raw energies, normal termination, basis counts and component identities were audited.','']
(ROOT/'final_results_appendix.txt').write_text('\n'.join(lines))
print('\n'.join(lines[:12]));print('Audited 45 selected components; saved comparison, spin plots and report appendix.')
