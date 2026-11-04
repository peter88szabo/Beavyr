#!/usr/bin/env python3
"""Run the user-requested Behemoth ROHF or MRCC F12 comparison, sequentially."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

BASE = Path(__file__).resolve().parent
engine = sys.argv[1]
assert engine in ('behemoth', 'mrcc-f12-dz')
core_counts = [int(value) for value in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1, 2, 4]
assert core_counts and all(value in (1, 2, 4) for value in core_counts)
memory_mb = 16000
geometry = (BASE / 'mrcc_nh3_ch_benchmark_2026-11-03/NH3_CH_COM_3A.xyz').read_text()
if engine == 'behemoth':
    root = BASE / 'behemoth_nh3_ch_benchmark_2026-11-03'
    exe = Path('/home/peter/DevCode/Behemoth/target/release/behemoth')
    method, basis, reference = 'UCCSD(T)', 'aug-cc-pVTZ', 'ROHF'
else:
    root = BASE / 'mrcc_f12_dz_nh3_ch_benchmark_2026-11-03'
    exe = Path('/home/peter/Programs/MRCC/binary/dmrcc')
    method, basis, reference = 'UCCSD(T)-F12', 'aug-cc-pVDZ', 'UHF'
root.mkdir(exist_ok=True)
assert not (root / 'results.json').exists(), 'Do not overwrite existing timings.'
(root / 'NH3_CH_COM_3A.xyz').write_text(geometry)
metadata_path = root / 'metadata.json'
metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
metadata.update(engine=engine, executable=str(exe), executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(), method=method, basis=basis, reference=reference, memory_MB=memory_mb, frozen_core=True, charge=0, multiplicity=2, core_counts=core_counts, geometry='NH3_CH_COM_3A.xyz', timing='GNU time, complete executable invocation, one fresh-directory run per core count', started=time.strftime('%Y-%m-%d %H:%M:%S %z'))
metadata_path.write_text(json.dumps(metadata, indent=2) + '\n')
rows = []
child = None

def stop(signum, frame):
    if child is not None and child.poll() is None:
        os.killpg(child.pid, signal.SIGTERM)
        try:
            child.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()
    raise SystemExit(128 + signum)

signal.signal(signal.SIGINT, stop)
signal.signal(signal.SIGTERM, stop)

for cores in core_counts:
    work = root / f'{cores}cores'
    work.mkdir(exist_ok=True)
    assert not (work / 'output.out').exists()
    (work / 'NH3_CH_COM_3A.xyz').write_text(geometry)
    if engine == 'behemoth':
        assert (work / 'input.inp').exists(), 'Inputs must be prepared before the timing.'
        args = [str(exe), 'input.inp']
    else:
        (work / 'MINP').write_text(f'calc=CCSD(T)-F12\nbasis=aug-cc-pVDZ\ncharge=0\nmult=2\nscftype=UHF\nmem={memory_mb}MB\ncore=frozen\nsymm=off\nunit=angs\nmolden=off\nscfmaxit=200\nccmaxit=100\ngeom=xyz\n6\n\n' + '\n'.join(geometry.splitlines()[2:]) + '\n')
        args = [str(exe)]
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS=str(cores), MKL_NUM_THREADS=str(cores), OPENBLAS_NUM_THREADS=str(cores), RAYON_NUM_THREADS=str(cores), PATH=str(exe.parent) + os.pathsep + env.get('PATH', ''))
    cpu_ids = ','.join(str(i) for i in range(cores))
    started = time.strftime('%Y-%m-%d %H:%M:%S %z')
    print('START', engine, method, basis, reference, cores, 'cores', started, flush=True)
    start = time.perf_counter()
    with (work / 'output.out').open('w') as out, (work / 'error.err').open('w') as err:
        child = subprocess.Popen(['/usr/bin/time', '-v', '-o', str(work / 'time.txt'), '/usr/bin/taskset', '-c', cpu_ids, *args], cwd=work, env=env, stdout=out, stderr=err, start_new_session=True)
        while True:
            try:
                code = child.wait(timeout=30)
                break
            except subprocess.TimeoutExpired:
                print('RUNNING', engine, cores, 'cores', round(time.perf_counter() - start, 1), 'seconds', flush=True)
    child = None
    log = (work / 'output.out').read_text(errors='replace')
    normalized = [' '.join(line.split()) for line in log.splitlines()]
    timing = (work / 'time.txt').read_text()
    match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)', timing)
    wall = 0.0
    assert match, 'GNU time did not record wall time'
    for part in match[1].split(':'):
        wall = wall * 60 + float(part)
    prefix = 'Total UCCSD(T) energy:' if engine == 'behemoth' else 'Total CCSD(F12*)(T) energy [au]:'
    energies = [float(line.split(':', 1)[1].strip().split()[0].replace('D', 'E')) for line in normalized if line.startswith(prefix)]
    reference_ok = 'Reference determinant: ROHF' in normalized if engine == 'behemoth' else 'Normal termination of mrcc.' in log
    failures = [line for line in normalized if any(text in line.lower() for text in ('did not converge', 'not converged', 'convergence not achieved', 'failed to converge', 'convergence not reached'))]
    ok = code == 0 and reference_ok and bool(energies) and math.isfinite(energies[-1]) and not failures
    def metric(key):
        found = re.search(r'^\s*' + re.escape(key) + r':\s*([\d.]+)', timing, re.M)
        return float(found[1]) if found else None
    row = dict(engine=engine, method=method, basis=basis, reference=reference, cores=cores, cpu_ids=cpu_ids, status='PASS' if ok else 'FAILED', wall_seconds=wall, energy_hartree=energies[-1] if energies else None, started=started, exit_code=code, failure='; '.join(failures), user_cpu_seconds=metric('User time (seconds)'), system_cpu_seconds=metric('System time (seconds)'), max_process_rss_kb=metric('Maximum resident set size (kbytes)'), directory=str(work))
    rows.append(row)
    (root / 'results.json').write_text(json.dumps(rows, indent=2) + '\n')
    with (root / 'results.csv').open('w') as file:
        writer = csv.DictWriter(file, fieldnames=list(row))
        writer.writeheader()
        writer.writerows(rows)
    print('DONE', engine, cores, 'cores', row['status'], f'{wall:.2f} seconds', 'energy', row['energy_hartree'], flush=True)
    if not ok:
        print('FAILED_OUTPUT', '\n'.join(log.splitlines()[-20:]), (work / 'error.err').read_text(errors='replace')[-2000:], flush=True)
        break
    if (root / 'STOP_AFTER_CURRENT').exists():
        print('STOP_AFTER_CURRENT requested; ending series.', flush=True)
        break
print('SERIES_FINISHED', engine, len(rows), 'runs', flush=True)
