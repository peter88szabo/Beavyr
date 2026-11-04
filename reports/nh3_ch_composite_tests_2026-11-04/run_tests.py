#!/usr/bin/env python3
"""Sequential MRCC five-point composite benchmark, resumable from accepted jobs."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time

ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'nh3_ch_com_scan_2026-11-03'
SCRATCH=ROOT.parents[1]/'target/benchmarks/nh3_ch_composite_tests_2026-11-04'
EXE=Path('/home/peter/Programs/MRCC/binary/dmrcc')
POINTS=[7,1,4,5,6]
ALPHA=4.071516
KCAL=627.5094740631
KJ=2625.4996394799

def cor_aux(z):
    return f'dfbasis_cor=atomtype\nN:aug-cc-pwCV{z}Z-RI\nH:aug-cc-pV{z}Z-RI\nC:aug-cc-pwCV{z}Z-RI\n'

jun_t='dfbasis_scf=aug-cc-pVQZ-RI-JK\n'+cor_aux('Q')+'dfbasis_cab=aug-cc-pVTZ-OPTRI\ngamma=1.0\n'
jun_q='dfbasis_scf=aug-cc-pV5Z-RI-JK\n'+cor_aux('5')+'dfbasis_cab=aug-cc-pVQZ-OPTRI\ngamma=1.0\n'
cv_basis='atomtype\nN:cc-pwCVTZ\nH:cc-pVTZ\nC:cc-pwCVTZ'
cv_aux='dfbasis_scf=aug-cc-pVQZ-RI-JK\n'+cor_aux('Q')+'dfbasis_cab=atomtype\nN:cc-pCVTZ-F12-OPTRI\nH:cc-pVTZ-F12-OPTRI\nC:cc-pCVTZ-F12-OPTRI\ngamma=1.0\n'
# Existing shared components are loaded from the original scan and never rerun.
SETTINGS=[
 ('cc_adz','CCSD(T)','aug-cc-pVDZ','frozen',''),
 ('cc_f12_adz','CCSD(T)-F12','aug-cc-pVDZ','frozen',''),
 ('mp2_f12_atz','MP2-F12','aug-cc-pVTZ','frozen',''),
 ('mp2_f12_tzf12','MP2-F12','cc-pVTZ-F12','frozen',''),
 ('mp2_f12_qzf12','MP2-F12','cc-pVQZ-F12','frozen',''),
 ('cc_f12_juntz','CCSD(T)-F12','jun-cc-pVTZ','frozen',jun_t),
 ('mp2_f12_junqz','MP2-F12','jun-cc-pVQZ','frozen',jun_q),
 ('mp2_f12_cv_fc','MP2-F12',cv_basis,'frozen',cv_aux),
 ('mp2_f12_cv_ae','MP2-F12',cv_basis,'corr',cv_aux),
]
protocol=dict(points=[2+8*i/7 for i in POINTS],cores=8,memory_MB=16000,reference='UHF',
 charge=0,multiplicity=2,scftol=9,cctol=8,scfmaxit=400,ccmaxit=200,
 cbs_basis_pair=['cc-pVTZ-F12','cc-pVQZ-F12'],cbs_correlation_exponent=ALPHA,
 cbs_reference='Larger-basis HF+CABS; extrapolate MP2-F12 correlation only',
 cbs_higher_order_correction_basis='aug-cc-pVDZ',
 reused_components='Original cc_f12_adz and mp2_f12_atz at all five points; original thresholds and timings retained',
 jun_label='MRCC adaptation of junChS-F12, not the published Molpro implementation',
 jun_extrapolation='n^-3 extrapolation of jun-T/Q MP2-F12 totals; subtract jun-T total to form correction',
 jun_core_correction='all-electron minus frozen-core MP2-F12/cc-pwCVTZ (H: cc-pVTZ)',
 f12='MRCC 2B/F+K MP2-F12, CCSD(F12*) with standard unscaled (T)',
 basis_source=str(ROOT/'bases/sources.json'),counterpoise=False,
 executable=str(EXE),sha256=hashlib.sha256(EXE.read_bytes()).hexdigest(),
 sources=['https://pubs.acs.org/doi/10.1021/ct500513b','https://arxiv.org/abs/2110.03075'])
(ROOT/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
SCRATCH.mkdir(parents=True,exist_ok=True)
jobs=json.loads((ROOT/'jobs.json').read_text()) if (ROOT/'jobs.json').exists() else []
child=None

def descendants(pid):
    result=[]; pending=[pid]
    while pending:
        p=pending.pop()
        try: more=list(map(int,Path(f'/proc/{p}/task/{p}/children').read_text().split()))
        except OSError: more=[]
        result.extend(more);pending.extend(more)
    return result

def stop(signum,frame):
    if child and child.poll() is None:
        for pid in reversed(descendants(child.pid)):
            try:os.kill(pid,signal.SIGTERM)
            except ProcessLookupError:pass
        try:os.killpg(child.pid,signal.SIGTERM)
        except ProcessLookupError:pass
        try:child.wait(timeout=5)
        except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
    raise SystemExit(128+signum)
signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)

def dump():
    (ROOT/'jobs.json').write_text(json.dumps(jobs,indent=2)+'\n')
    if jobs:
        with (ROOT/'jobs.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for j in jobs for k in j)));w.writeheader();w.writerows(jobs)

def parse(log,prefix):
    lines=[' '.join(l.split()) for l in log.splitlines()]
    values=[float(l.split(':',1)[1].split()[0].replace('D','E')) for l in lines if l.startswith(prefix)]
    return values[-1] if values else None

def cost(j):
    return j['wall_seconds']+j.get('restart_preparation_seconds',0.0)

def summarize():
    rows=[]
    methods=['direct_f12_adz','original_composite','conventional_small_basis_correction','f12_cbs_composite','junchs_f12_mrcc']
    for i in sorted(POINTS):
        g={j['job']:j for j in jobs if j['point']==i and j['status']=='PASS'}
        low=g.get('cc_f12_adz');high=g.get('mp2_f12_atz');small=g.get('cc_adz');t=g.get('mp2_f12_tzf12');q=g.get('mp2_f12_qzf12')
        jt=g.get('cc_f12_juntz');jq=g.get('mp2_f12_junqz');fc=g.get('mp2_f12_cv_fc');ae=g.get('mp2_f12_cv_ae')
        mp2cbs=(q['reference_hartree']+q['cabs_hartree']+(4**ALPHA*q['mp2_f12_correlation']-3**ALPHA*t['mp2_f12_correlation'])/(4**ALPHA-3**ALPHA)) if t and q else None
        juncbs=(64*jq['energy_hartree']-27*jt['mp2_f12_hartree'])/37 if jt and jq else None
        cv=ae['energy_hartree']-fc['energy_hartree'] if ae and fc else None
        rows.append(dict(point=i,distance_angstrom=2+8*i/7,
            direct_f12_adz=low['energy_hartree'],
            original_composite=low['energy_hartree']+high['energy_hartree']-low['mp2_f12_hartree'],
            conventional_small_basis_correction=high['energy_hartree']+small['energy_hartree']-small['mp2_hartree'] if small else None,
            f12_cbs_composite=mp2cbs+low['energy_hartree']-low['mp2_f12_hartree'] if mp2cbs is not None else None,
            junchs_f12_mrcc=jt['energy_hartree']+juncbs-jt['mp2_f12_hartree']+cv if juncbs is not None and cv is not None else None,
            mp2_f12_cbs=mp2cbs,jun_mp2_cbs=juncbs,core_valence_correction=cv,
            direct_seconds=cost(low),original_composite_seconds=cost(low)+cost(high),
            conventional_correction_seconds=cost(small)+cost(high) if small else None,
            conventional_correction_new_seconds=cost(small) if small else None,
            f12_cbs_seconds=cost(low)+cost(t)+cost(q) if t and q else None,
            f12_cbs_new_seconds=cost(t)+cost(q) if t and q else None,
            junchs_seconds=sum(cost(j) for j in [jt,jq,fc,ae]) if all([jt,jq,fc,ae]) else None,
            all_new_component_seconds=sum(cost(j) for j in g.values() if not j['reused'])))
    far=next(r for r in rows if r['point']==7)
    for row in rows:
        for method in methods:
            delta=row[method]-far[method] if row[method] is not None and far[method] is not None else None
            row[method+'_relative_kcal_mol']=delta*KCAL if delta is not None else None
            row[method+'_relative_kj_mol']=delta*KJ if delta is not None else None
    (ROOT/'points.json').write_text(json.dumps(rows,indent=2)+'\n')
    with (ROOT/'comparison.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return rows

# The old 5.428571 A job must finish before any new engine process is started.
while True:
    previous=json.loads((OLD/'approved_jobs.json').read_text())
    if any(j['point']==3 and j['job']=='cc_atz' for j in previous):break
    print('WAITING for the existing 5.428571 A Behemoth calculation; no concurrent job.',flush=True)
    time.sleep(20)

for i in POINTS:
    geometry=(OLD/'runs'/f'p{i+1:02d}_{2+8*i/7:.6f}A_cc_f12_adz'/'geometry.xyz').read_text()
    for key,calc,basis,core,extra in SETTINGS:
        if any(j['point']==i and j['job']==key for j in jobs):continue
        if key in ('mp2_f12_junqz','mp2_f12_cv_fc','mp2_f12_cv_ae') and any(j['point']==i and j['job']=='cc_f12_juntz' and j['status']!='PASS' for j in jobs):
            print('SKIPPED',i,key,'required CC component failed; no retry',flush=True)
            continue
        if key in ('cc_f12_adz','mp2_f12_atz'):
            raise RuntimeError('A required reused component is missing; do not rerun it')
        if (ROOT/'STOP_AFTER_CURRENT').exists():raise SystemExit(0)
        name=f'p{i+1:02d}_{2+8*i/7:.6f}A_{key}'
        work=SCRATCH/name;evidence=ROOT/'runs'/name
        work.mkdir(exist_ok=False);evidence.mkdir(parents=True,exist_ok=False)
        inp=f'calc={calc}\nbasis={basis}\ncharge=0\nmult=2\nscftype=UHF\nmem=16000MB\ncore={core}\nsymm=off\nunit=angs\nmolden=off\nscftol=9\ncctol=8\nscfmaxit=400\nccmaxit=200\n'+extra+'geom=xyz\n6\n\n'+'\n'.join(geometry.splitlines()[2:])+'\n'
        preparation_seconds=0.0
        if i==5 and key in ('mp2_f12_qzf12','cc_f12_juntz'):
            filename='failed_scf_7p714286_qzf12.json' if key=='mp2_f12_qzf12' else 'failed_scf_7p714286_juntz.json'
            failed=json.loads((ROOT/filename).read_text())[0]
            inp=inp.replace('scfmaxit=400', 'scfmaxit=400\nscfiguess=restart\nscfext=20\nscflshift=off')
            for filename in ('SCFDENSITIES','VARS'):
                shutil.copy2(Path(failed['scratch'])/filename,work/filename)
            preparation_seconds=failed['wall_seconds']
            (evidence/'restart_source.json').write_text(json.dumps(failed,indent=2)+'\n')
        for folder in (work,evidence):
            (folder/'MINP').write_text(inp);(folder/'geometry.xyz').write_text(geometry)
        if 'jun' in key:
            for folder in (work,evidence):shutil.copy2(ROOT/'bases/GENBAS',folder/'GENBAS')
        env=dict(os.environ,OMP_NUM_THREADS='8',MKL_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',RAYON_NUM_THREADS='8')
        env['PATH']=str(EXE.parent)+os.pathsep+env.get('PATH','')
        (ROOT/'active.json').write_text(json.dumps(dict(job=name,point=i,output=str(evidence/'output.out')),indent=2))
        print('START',name,flush=True);start=time.perf_counter()
        with (evidence/'output.out').open('w') as out,(evidence/'error.err').open('w') as err:
            child=subprocess.Popen(['/usr/bin/time','-v','-o',str(evidence/'time.txt'),'/usr/bin/taskset','-c','0-7',str(EXE)],cwd=work,env=env,stdout=out,stderr=err,start_new_session=True)
            while True:
                try:code=child.wait(timeout=20);break
                except subprocess.TimeoutExpired:print('RUNNING',name,f'{time.perf_counter()-start:.1f} seconds',flush=True)
        child=None;log=(evidence/'output.out').read_text(errors='replace')
        energy_prefix='Total CCSD(T) energy [au]:' if key=='cc_adz' else ('Total CCSD(F12*)(T) energy [au]:' if key.startswith('cc_') else 'MP2-F12 energy [au]:')
        energy=parse(log,energy_prefix)
        spin=re.search(r'Expectation value of S\^2:\s*([-+0-9.Ee]+)',log)
        s2=float(spin[1]) if spin else None
        failures=[l.strip() for l in log.splitlines() if any(t in l.lower() for t in ('not converged','did not converge','failed to converge','convergence not achieved','convergence not reached'))]
        ok=code==0 and 'Normal termination of mrcc.' in log and energy is not None and math.isfinite(energy) and not failures and s2 is not None and .74<s2<.9
        mp2=parse(log,'Total MP2 energy [au]:') if key=='cc_adz' else parse(log,'MP2 energy [au]:')
        f12=parse(log,'MP2-F12 energy [au]:')
        if key=='cc_adz' and mp2 is None:ok=False
        if key.startswith('cc_f12') and f12 is None:ok=False
        timing=(evidence/'time.txt').read_text();m=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)',timing);wall=0.
        if m:
            for part in m[1].split(':'):wall=wall*60+float(part)
        else:wall=time.perf_counter()-start
        rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',timing)
        row=dict(point=i,distance_angstrom=2+8*i/7,job=key,calc=calc,basis=basis,core=core,status='PASS' if ok else 'FAILED',energy_hartree=energy,mp2_hartree=mp2,mp2_f12_hartree=f12,mp2_f12_correlation=parse(log,'MP2-F12 correlation energy [au]:'),reference_hartree=parse(log,'Reference energy [au]:'),cabs_hartree=parse(log,'CABS singles correction [au]:'),s2=s2,excess_s2=s2-.75 if s2 is not None else None,wall_seconds=wall,peak_rss_kib=int(rss[1]) if rss else None,exit_code=code,failure='; '.join(failures),evidence=str(evidence),scratch=str(work),reused=False,restart_preparation_seconds=preparation_seconds)
        jobs.append(row);dump();summarize();print('DONE',name,row['status'],f'{wall:.2f} seconds',energy,'<S^2>',s2,flush=True)
        if not ok:
            print('FAILURE', '\n'.join(log.splitlines()[-15:]), (evidence/'error.err').read_text()[-2000:],flush=True)
            continue
    print('POINT_COMPLETE',i,flush=True)
rows=summarize()
(ROOT/'active.json').write_text(json.dumps(dict(status='COMPLETE_WITH_DISCARDED_FAILURES' if any(j['status']!='PASS' for j in jobs) else 'COMPLETE',points=len(rows),jobs=len(jobs)),indent=2))
print('COMPLETE',len(rows),'points',len(jobs),'jobs',flush=True)
