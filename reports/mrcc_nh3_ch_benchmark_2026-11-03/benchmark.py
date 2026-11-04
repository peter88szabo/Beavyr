#!/usr/bin/env python3
"""Sequential MRCC NH3+CH timings; generated inputs/output are retained."""
import csv,json,math,os,re,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
SCRATCH=Path('/home/peter/DevCode/Beavyr/target/benchmarks')/('mrcc_nh3_ch_'+time.strftime('%Y%m%d_%H%M%S'))
EXE=Path('/home/peter/Programs/MRCC/binary/dmrcc')
MASS={'H':1.008,'C':12.011,'N':14.007}
NH3=[('N',[-0.00845252,-0.00456017,0.02826060]),('H',[1.00399351,-0.01467638,-0.02173799]),('H',[-0.28913590,0.96859848,-0.01450254]),('H',[-0.34319744,-0.43851852,-0.82454485])]
CH=[('C',[0.,0.,0.]),('H',[0.,0.,1.12])]
def com(fragment):
    total=sum(MASS[a] for a,r in fragment)
    return [sum(MASS[a]*r[k] for a,r in fragment)/total for k in range(3)]
def place(fragment,target):
    center=com(fragment)
    return [(a,[r[k]-center[k]+target[k] for k in range(3)]) for a,r in fragment]
geometry=place(NH3,[0.,0.,0.])+place(CH,[0.,0.,3.])
body=''.join(f'{a} {r[0]:.14f} {r[1]:.14f} {r[2]:.14f}\n' for a,r in geometry)
xyz=f'6\nNH3 + CH; neutral doublet; fragment COM separation 3.0 Angstrom; CH bond 1.12 Angstrom\n{body}'
ROOT.mkdir(parents=True,exist_ok=True);SCRATCH.mkdir(parents=True,exist_ok=False)
(ROOT/'NH3_CH_COM_3A.xyz').write_text(xyz)
separation=math.dist(com(geometry[:4]),com(geometry[4:]))
assert abs(separation-3.)<1e-12
assert sum({'N':7,'C':6,'H':1}[a] for a,r in geometry)==17
metadata={'scratch':str(SCRATCH),'executable':str(EXE),'charge':0,'multiplicity':2,'electrons':17,'COM_distance_angstrom':separation,'masses_amu':MASS,'NH3_source':'user supplied coordinates; translation only','CH_bond_angstrom':1.12,'orientation':'NH3 kept in original orientation; CH parallel to +z, C toward NH3, H away; NH3 COM at origin; CH COM at (0,0,3)','memory_MB':2048,'core':'frozen','repetitions':1,'host':subprocess.check_output(['lscpu'],text=True),'initial_memory':subprocess.check_output(['free','-h'],text=True)}
(ROOT/'metadata.json').write_text(json.dumps(metadata,indent=2))
rows=[];child=None

def stop(signum,frame):
    if child is not None and child.poll() is None:
        os.killpg(child.pid,signal.SIGTERM)
        try:child.wait(timeout=5)
        except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
    raise SystemExit(128+signum)
signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)

def flush():
    (ROOT/'results.json').write_text(json.dumps(rows,indent=2))
    if rows:
        keys=list(dict.fromkeys(k for row in rows for k in row))
        with (ROOT/'results.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)

settings=[('UCCSD(T)','CCSD(T)','cc-pVTZ'),('UCCSD(T)','CCSD(T)','aug-cc-pVTZ'),('UCCSD(T)-F12','CCSD(T)-F12','cc-pVTZ-F12'),('UCCSD(T)-F12','CCSD(T)-F12','aug-cc-pVTZ')]
print('GEOMETRY_READY',str(ROOT/'NH3_CH_COM_3A.xyz'),'COM=',separation,flush=True)
for si,(label,calc,basis) in enumerate(settings):
    for cores in [1,2,4]:
        name=f'{si+1}_{basis}_{cores}cores'
        work=SCRATCH/name;work.mkdir()
        evidence=ROOT/name;evidence.mkdir(exist_ok=False)
        inp=f'calc={calc}\nbasis={basis}\ncharge=0\nmult=2\nscftype=UHF\nmem=2048MB\ncore=frozen\nsymm=off\nunit=angs\nmolden=off\nscfmaxit=200\nccmaxit=100\ngeom=xyz\n6\n\n{body}'
        (work/'MINP').write_text(inp);(evidence/'MINP').write_text(inp)
        env=dict(os.environ,PATH=str(EXE.parent)+os.pathsep+os.environ.get('PATH',''),OMP_NUM_THREADS=str(cores),MKL_NUM_THREADS=str(cores))
        cpu_ids=','.join(str(i) for i in range(cores))
        startload=os.getloadavg();start=time.perf_counter();started=time.strftime('%Y-%m-%d %H:%M:%S %z')
        print('START',name,label,basis,'cores=',cores,'at',started,flush=True)
        timedout=False
        with (evidence/'mrcc.out').open('w') as out,(evidence/'mrcc.err').open('w') as err:
            command=['/usr/bin/time','-v','-o',str(evidence/'time.txt'),'/usr/bin/taskset','-c',cpu_ids,str(EXE)]
            child=subprocess.Popen(command,cwd=work,env=env,stdout=out,stderr=err,start_new_session=True)
            next_report=start+30
            while child.poll() is None:
                time.sleep(.25);now=time.perf_counter()
                if now>=next_report:
                    print('RUNNING',name,'elapsed_s=',round(now-start,1),flush=True);next_report=now+30
                if now-start>3600:
                    timedout=True;os.killpg(child.pid,signal.SIGTERM)
                    try:child.wait(timeout=5)
                    except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
                    break
            code=child.wait();elapsed=time.perf_counter()-start;child=None
        log=(evidence/'mrcc.out').read_text(errors='replace');normalized=[' '.join(line.split()) for line in log.splitlines()]
        prefix='Total CCSD(F12*)(T) energy [au]:' if 'F12' in calc else 'Total CCSD(T) energy [au]:'
        energies=[float(line.split(':',1)[1].strip().split()[0].replace('D','E')) for line in normalized if line.startswith(prefix)]
        failures=[line for line in normalized if any(s in line.lower() for s in ['convergence not achieved','not converged','convergence not reached','failed to converge'])]
        ok=code==0 and 'Normal termination of mrcc.' in log and bool(energies) and math.isfinite(energies[-1]) and not failures and not timedout
        timing=(evidence/'time.txt').read_text(errors='replace')
        def timevalue(key):
            match=re.search(r'^\s*'+re.escape(key)+r':\s*([\d.]+)',timing,re.M)
            return float(match.group(1)) if match else None
        row={'method':label,'basis':basis,'cores':cores,'cpu_ids':cpu_ids,'status':'PASS' if ok else 'FAILED','wall_seconds':elapsed,'energy_hartree':energies[-1] if ok else None,'user_cpu_seconds':timevalue('User time (seconds)'),'system_cpu_seconds':timevalue('System time (seconds)'),'max_process_rss_kb':timevalue('Maximum resident set size (kbytes)'),'started':started,'load1_start':startload[0],'load1_end':os.getloadavg()[0],'exit_code':code,'failure':'; '.join(failures) or ('timeout' if timedout else ''),'directory':str(work)}
        s2=[line for line in normalized if any(s in line.lower() for s in ['s**2','s^2','<s2>'])]
        row['spin_diagnostics']=' | '.join(s2[-4:])
        rows.append(row);flush();print('DONE',name,row['status'],'wall_s=',round(elapsed,3),'energy=',row['energy_hartree'],flush=True)
        if not ok:print('FAILURE_TAIL','\n'.join(log.splitlines()[-20:]),flush=True)
print('ALL_FINISHED',sum(r['status']=='PASS' for r in rows),'/',len(rows),flush=True)
