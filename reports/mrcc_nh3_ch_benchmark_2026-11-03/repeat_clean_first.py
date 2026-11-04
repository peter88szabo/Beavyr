#!/usr/bin/env python3
"""Repeat the pre-shutdown first case; retain it as excluded warm-up evidence."""
import csv,json,math,os,re,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
rows=json.loads((ROOT/'results.json').read_text());assert len(rows)==12
original=rows[0];assert original['basis']=='cc-pVTZ' and original['cores']==1
metadata=json.loads((ROOT/'metadata.json').read_text())
work=Path(metadata['scratch'])/'clean_cc-pVTZ_1core';work.mkdir()
evidence=ROOT/'clean_cc-pVTZ_1core';evidence.mkdir()
text=(Path(original['directory'])/'MINP').read_text()
(work/'MINP').write_text(text);(evidence/'MINP').write_text(text)
exe=Path(metadata['executable']);env=dict(os.environ,PATH=str(exe.parent)+os.pathsep+os.environ.get('PATH',''),OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
start=time.perf_counter();started=time.strftime('%Y-%m-%d %H:%M:%S %z');load=os.getloadavg()[0]
print('CLEAN_REPEAT_START',started,flush=True)
with (evidence/'mrcc.out').open('w') as out,(evidence/'mrcc.err').open('w') as err:
 child=subprocess.Popen(['/usr/bin/time','-v','-o',str(evidence/'time.txt'),'/usr/bin/taskset','-c','0',str(exe)],cwd=work,env=env,stdout=out,stderr=err,start_new_session=True)
 def stop(sig,frame):
  os.killpg(child.pid,signal.SIGTERM);child.wait();raise SystemExit(128+sig)
 signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
 while True:
  try:code=child.wait(timeout=30);break
  except subprocess.TimeoutExpired:print('CLEAN_REPEAT_RUNNING',round(time.perf_counter()-start,1),flush=True)
wall=time.perf_counter()-start
log=(evidence/'mrcc.out').read_text(errors='replace');lines=[' '.join(l.split()) for l in log.splitlines()]
energies=[float(l.split(':',1)[1].strip().split()[0].replace('D','E')) for l in lines if l.startswith('Total CCSD(T) energy [au]:')]
failures=[l for l in lines if any(t in l.lower() for t in ['not converged','convergence not achieved','failed to converge','convergence not reached'])]
ok=code==0 and 'Normal termination of mrcc.' in log and energies and math.isfinite(energies[-1]) and not failures
timing=(evidence/'time.txt').read_text()
def val(key):
 m=re.search(r'^\s*'+re.escape(key)+r':\s*([\d.]+)',timing,re.M);return float(m[1]) if m else None
for row in rows:row['included_in_comparison']=True
rows[0]['included_in_comparison']=False;rows[0]['exclusion_reason']='Started while user jobs were still active; replaced by clean repeat.'
row=dict(original,status='PASS' if ok else 'FAILED',wall_seconds=wall,energy_hartree=energies[-1] if ok else None,user_cpu_seconds=val('User time (seconds)'),system_cpu_seconds=val('System time (seconds)'),max_process_rss_kb=val('Maximum resident set size (kbytes)'),started=started,load1_start=load,load1_end=os.getloadavg()[0],exit_code=code,failure='; '.join(failures),directory=str(work),included_in_comparison=True)
row['spin_diagnostics']=' | '.join([l for l in lines if any(s in l.lower() for s in ['s**2','s^2','<s2>'])][-4:])
rows.append(row);(ROOT/'results.json').write_text(json.dumps(rows,indent=2))
keys=list(dict.fromkeys(k for row in rows for k in row))
with (ROOT/'results.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
print('CLEAN_REPEAT_DONE',row['status'],wall,row['energy_hartree'],flush=True)
