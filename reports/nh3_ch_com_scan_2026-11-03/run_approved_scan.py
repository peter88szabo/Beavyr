#!/usr/bin/env python3
"""Sequential, resumable MRCC rigid COM scan. No application source changes."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import time

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'mrcc_nh3_ch_benchmark_2026-11-03/NH3_CH_COM_3A.xyz'
EXE = Path('/home/peter/Programs/MRCC/binary/dmrcc')
SCRATCH = ROOT.parents[1] / 'target/benchmarks/nh3_ch_com_scan_2026-11-03'
CORES, MEMORY = 8, 16000
HARTREE_KCAL = 627.5094740631
HARTREE_KJ = 2625.4996394799
MASSES = {'N': 14.007, 'C': 12.011, 'H': 1.008}
DISTANCES = [2.0 + 8.0 * i / 7 for i in range(8)]
SETTINGS = [('cc_f12_adz', 'CCSD(T)-F12', 'aug-cc-pVDZ'),
            ('mp2_f12_atz', 'MP2-F12', 'aug-cc-pVTZ'),
            ('cc_atz', 'CCSD(T)', 'aug-cc-pVTZ')]
atoms = [(p[0], list(map(float, p[1:]))) for line in SOURCE.read_text().splitlines()[2:] if (p := line.split())]
assert len(atoms) == 6
def com(fragment):
    mass = sum(MASSES[a] for a, _ in fragment)
    return [sum(MASSES[a] * xyz[k] for a, xyz in fragment) / mass for k in range(3)]
start_a, start_b = com(atoms[:4]), com(atoms[4:])
assert abs(math.dist(start_a, start_b) - 3.0) < 1e-10
geometries = []
for distance in DISTANCES:
    moved = [(a, [xyz[k] - start_a[k] for k in range(3)]) for a, xyz in atoms[:4]]
    moved += [(a, [xyz[k] - start_b[k] + (distance if k == 2 else 0.0) for k in range(3)]) for a, xyz in atoms[4:]]
    assert abs(math.dist(com(moved[:4]), com(moved[4:])) - distance) < 1e-12
    geometries.append(moved)
def xyz_text(i):
    return f'6\nNH3 + CH; neutral doublet; COM {DISTANCES[i]:.12f} Angstrom; rigid fragments\n' + ''.join(f'{a} {p[0]:.14f} {p[1]:.14f} {p[2]:.14f}\n' for a, p in geometries[i])

ROOT.mkdir(exist_ok=True)
SCRATCH.mkdir(parents=True, exist_ok=True)
(ROOT / 'scan_geometries.xyz').write_text(''.join(xyz_text(i) for i in range(8)))
metadata = dict(conventional_engine='Behemoth UHF UCCSD(T), exact Cholesky threshold 1e-10', behemoth_executable='/home/peter/DevCode/Behemoth/target/release/behemoth', behemoth_sha256=hashlib.sha256(Path('/home/peter/DevCode/Behemoth/target/release/behemoth').read_bytes()).hexdigest(), executable=str(EXE), executable_sha256=hashlib.sha256(EXE.read_bytes()).hexdigest(), source_geometry=str(SOURCE), distances_angstrom=DISTANCES, execution_point_order=[7, 0, 1, 4, 5, 6, 2, 3], cores=CORES, cpu_ids=list(range(CORES)), memory_MB=MEMORY, reference='UHF', charge=0, multiplicity=2, core='frozen', rigid_fragments=True, masses_amu=MASSES, composite='CCSD(F12*)(T)/aDZ + MP2-F12/aTZ - MP2-F12/aDZ', comparison='Conventional CCSD(T)/aTZ', zero='Each curve independently referenced to its own 10.0 Angstrom energy', mp2_adz_source='MP2-F12 total energy, including CABS singles, printed within the CCSD(T)-F12/aDZ job', BSSE='No counterpoise correction; same full-complex prescription throughout', scratch=str(SCRATCH))
(ROOT / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
jobs = json.loads((ROOT / 'approved_jobs.json').read_text())
child = None
current = None

def descendants(pid):
    found = []
    pending = [pid]
    while pending:
        parent = pending.pop()
        try:
            children = Path(f'/proc/{parent}/task/{parent}/children').read_text().split()
        except OSError:
            continue
        found.extend(int(p) for p in children)
        pending.extend(int(p) for p in children)
    return found

def stop(signum, frame):
    if child is not None and child.poll() is None:
        pids = descendants(child.pid)
        for pid in reversed(pids):
            try: os.kill(pid, signal.SIGTERM)
            except ProcessLookupError: pass
        try: os.killpg(child.pid, signal.SIGTERM)
        except ProcessLookupError: pass
        try: child.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()
        for pid in pids:
            try: os.kill(pid, signal.SIGKILL)
            except ProcessLookupError: pass
    print('STOPPED', current, flush=True)
    raise SystemExit(128 + signum)
signal.signal(signal.SIGTERM, stop)
signal.signal(signal.SIGINT, stop)

def flush():
    (ROOT / 'approved_jobs.json').write_text(json.dumps(jobs, indent=2) + '\n')
    if jobs:
        fields = list(dict.fromkeys(key for row in jobs for key in row))
        with (ROOT / 'approved_jobs.csv').open('w') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader(); writer.writerows(jobs)
    points = []
    for i, distance in enumerate(DISTANCES):
        group = {r['job']: r for r in jobs if r['point'] == i and r['status'] == 'PASS'}
        if len(group) != 3: continue
        low, high, conv = group['cc_f12_adz'], group['mp2_f12_atz'], group['cc_atz']
        composite = low['energy_hartree'] + high['energy_hartree'] - low['mp2_f12_hartree']
        points.append(dict(point=i, distance_angstrom=distance, cc_f12_adz_hartree=low['energy_hartree'], mp2_f12_adz_hartree=low['mp2_f12_hartree'], mp2_f12_atz_hartree=high['energy_hartree'], basis_correction_hartree=high['energy_hartree']-low['mp2_f12_hartree'], composite_hartree=composite, cc_atz_hartree=conv['energy_hartree'], composite_wall_seconds=low['wall_seconds']+high['wall_seconds'], cc_atz_wall_seconds=conv['wall_seconds']))
    far = next((r for r in points if r['point'] == 7), None)
    if far:
        for row in points:
            for method in ('composite', 'cc_atz', 'cc_f12_adz'):
                delta = row[method + '_hartree'] - far[method + '_hartree']
                row[method + '_relative_kcal_mol'] = delta * HARTREE_KCAL
                row[method + '_relative_kj_mol'] = delta * HARTREE_KJ
            row['composite_minus_cc_atz_kcal_mol'] = row['composite_relative_kcal_mol'] - row['cc_atz_relative_kcal_mol']
    (ROOT / 'points.json').write_text(json.dumps(points, indent=2) + '\n')
    if points:
        with (ROOT / 'scan_data.csv').open('w') as f:
            writer = csv.DictWriter(f, fieldnames=list(points[0]))
            writer.writeheader(); writer.writerows(points)
    return points

def energy(lines, prefix):
    values = [float(line.split(':', 1)[1].split()[0].replace('D', 'E')) for line in lines if line.startswith(prefix)]
    return values[-1] if values else None

for i in [7, 0, 1, 4, 5, 6, 2, 3]:
    if i == 2:
        print("DEFERRED 4.285714 A Behemoth CCSD: solver convergence unresolved", flush=True)
        continue
    for job, calc, basis in SETTINGS:
        if any(r['point'] == i and r['job'] == job and r['status'] == 'PASS' for r in jobs): continue
        if (ROOT / 'STOP_AFTER_CURRENT').exists():
            print('STOP_AFTER_CURRENT: no next job started.', flush=True)
            raise SystemExit(0)
        current = f'p{i+1:02d}_{DISTANCES[i]:.6f}A_' + ('behemoth_cc_atz' if job == 'cc_atz' else job)
        selected_exe = Path('/home/peter/DevCode/Behemoth/target/release/behemoth') if job == 'cc_atz' else EXE
        work = SCRATCH / current
        evidence = ROOT / 'runs' / current
        work.mkdir(exist_ok=False)
        evidence.mkdir(parents=True, exist_ok=False)
        xyz = xyz_text(i)
        inp = f'calc={calc}\nbasis={basis}\ncharge=0\nmult=2\nscftype=UHF\nmem={MEMORY}MB\ncore=frozen\nsymm=off\nunit=angs\nmolden=off\nscfmaxit=200\nccmaxit=100\ngeom=xyz\n6\n\n' + '\n'.join(xyz.splitlines()[2:]) + '\n'
        (work / 'MINP').write_text(inp)
        (evidence / 'MINP').write_text(inp)
        (evidence / 'geometry.xyz').write_text(xyz)
        (work / 'geometry.xyz').write_text(xyz)
        if job == 'cc_atz':
            bh_input = f'nproc = {CORES}\nmemory = {MEMORY}\nmethod = UHF CCSD(T)\nscf_fallback = false\nbasis = aug-cc-pVTZ\ncharge = 0\nmulti = 2\nxyzfile = "geometry.xyz"\n@scf\n  max_iterations = 2000\n  diis_max = 16\nend\n@ccsd\n  frozen_core = auto\n  max_iterations = 100\nend\n'
            if i == 1:
                bh_input = bh_input.replace('  diis_max = 16', '  diis_max = 16\n  density_threshold = 3e-8').replace('max_iterations = 2000', 'max_iterations = 400')
            if i in (2, 3):
                bh_input = bh_input.replace('  diis_max = 16', '  diis_max = 10')
            bh_input = bh_input.replace('max_iterations = 100', 'max_iterations = 400')
            (work / 'input.inp').write_text(bh_input)
            (evidence / 'input.inp').write_text(bh_input)
        env = dict(os.environ, OMP_NUM_THREADS=str(CORES), MKL_NUM_THREADS=str(CORES), OPENBLAS_NUM_THREADS=str(CORES), RAYON_NUM_THREADS=str(CORES))
        env['PATH'] = str(selected_exe.parent) + os.pathsep + env.get('PATH', '')
        start = time.perf_counter()
        started = time.strftime('%Y-%m-%d %H:%M:%S %z')
        print('START', current, started, flush=True)
        (ROOT / 'active.json').write_text(json.dumps(dict(job=current, point=i, distance_angstrom=DISTANCES[i], output=str(evidence/'output.out'), started=started), indent=2))
        with (evidence / 'output.out').open('w') as out, (evidence / 'error.err').open('w') as err:
            child = subprocess.Popen(['/usr/bin/time', '-v', '-o', str(evidence/'time.txt'), '/usr/bin/taskset', '-c', '0-7', str(selected_exe), *(['input.inp'] if job == 'cc_atz' else [])], cwd=work, env=env, stdout=out, stderr=err, start_new_session=True)
            while True:
                try:
                    code = child.wait(timeout=20)
                    break
                except subprocess.TimeoutExpired:
                    print('RUNNING', current, round(time.perf_counter()-start, 1), 'seconds', flush=True)
        child = None
        log = (evidence / 'output.out').read_text(errors='replace')
        lines = [' '.join(line.split()) for line in log.splitlines()]
        prefix = {'cc_f12_adz': 'Total CCSD(F12*)(T) energy [au]:', 'mp2_f12_atz': 'MP2-F12 energy [au]:', 'cc_atz': 'Total UCCSD(T) energy:'}[job]
        total = energy(lines, prefix)
        mp2 = energy(lines, 'MP2-F12 energy [au]:')
        failures = [line for line in lines if any(s in line.lower() for s in ['not converged', 'did not converge', 'convergence not achieved', 'failed to converge', 'convergence not reached'])]
        ok = code == 0 and (('Reference determinant: UHF' in lines) if job == 'cc_atz' else 'Normal termination of mrcc.' in log) and total is not None and math.isfinite(total) and not failures and (job != 'cc_f12_adz' or mp2 is not None and math.isfinite(mp2))
        if job == 'cc_atz':
            actual = re.search(r'Actual <S\^2>:\s*([-+0-9.Ee]+)', log)
            expected_job = next(r for r in jobs if r['point'] == i and r['job'] == 'mp2_f12_atz' and r['status'] == 'PASS')
            expected_log = (Path(expected_job['evidence'])/'output.out').read_text(errors='replace')
            expected = re.search(r'Expectation value of S\^2:\s*([-+0-9.Ee]+)', expected_log)
            if not actual or not expected or abs(float(actual[1])-float(expected[1])) > 0.001:
                ok = False
                failures.append('UHF reference spin diagnostic differs from MRCC aTZ component; reject mixed electronic states')
        timing = (evidence / 'time.txt').read_text()
        match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)', timing)
        wall = 0.0
        if match:
            for part in match[1].split(':'): wall = wall*60 + float(part)
        else:
            wall = time.perf_counter()-start
        rss = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', timing)
        row = dict(engine='Behemoth' if job == 'cc_atz' else 'MRCC', point=i, distance_angstrom=DISTANCES[i], job=job, calc=calc, basis=basis, status='PASS' if ok else 'FAILED', energy_hartree=total, mp2_f12_hartree=mp2, wall_seconds=wall, peak_rss_kib=int(rss[1]) if rss else None, started=started, exit_code=code, failure='; '.join(failures), spin_diagnostics=' | '.join(line for line in lines if 'S^2' in line or 'S**2' in line), evidence=str(evidence), scratch=str(work))
        jobs.append(row)
        points = flush()
        print('DONE', current, row['status'], f'{wall:.2f} seconds', 'energy', total, flush=True)
        if not ok:
            print('FAILURE_TAIL', '\n'.join(log.splitlines()[-25:]), flush=True)
            raise SystemExit(1)
    points = flush()
    point = next(r for r in points if r['point'] == i)
    print('POINT_COMPLETE', json.dumps(point), flush=True)
(ROOT / 'active.json').write_text(json.dumps({'status': 'COMPLETE' if len(points)==8 else 'INCOMPLETE', 'jobs': len(jobs), 'points': len(points)}))
print('SCAN_COMPLETE' if len(points)==8 else 'SCAN_INCOMPLETE', len(points), 'points', len(jobs), 'jobs', flush=True)
