#!/usr/bin/env python3
"""Summarize the completed, sequential NH3 + CH benchmark from retained logs."""
import csv
import json
import re
from pathlib import Path

root = Path(__file__).resolve().parent
rows = json.loads((root / 'results.json').read_text())
assert len(rows) == 13, 'Finish the original series and clean one-core repeat first.'
for row in rows:
    evidence = root / Path(row['directory']).name
    timing = (evidence / 'time.txt').read_text()
    match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)', timing)
    assert match, evidence
    seconds = 0.0
    for part in match[1].split(':'):
        seconds = seconds * 60 + float(part)
    row.setdefault('launcher_wall_seconds', row['wall_seconds'])
    row['wall_seconds'] = seconds
    row['evidence_directory'] = str(evidence)
    row['timing_source'] = 'GNU /usr/bin/time elapsed wall clock'

(root / 'results.json').write_text(json.dumps(rows, indent=2) + '\n')
keys = list(dict.fromkeys(key for row in rows for key in row))
with (root / 'results.csv').open('w') as file:
    writer = csv.DictWriter(file, fieldnames=keys)
    writer.writeheader()
    writer.writerows(rows)

included = [row for row in rows if row['included_in_comparison']]
assert len(included) == 12
groups = {}
for row in included:
    groups.setdefault((row['method'], row['basis']), {})[row['cores']] = row

lines = [
    'MRCC NH3 + CH single-point timing benchmark',
    '==========================================',
    '',
    'System: neutral NH3 + CH, 17 electrons, doublet (multiplicity 2).',
    'The supplied NH3 geometry was translated without rotation or internal changes.',
    'CH was constructed with a C-H distance of 1.12 Angstrom, along +z, carbon',
    'toward ammonia. Fragment COM separation is 3.00000000000000 Angstrom.',
    'Masses: N 14.007, C 12.011, H 1.008 amu. Geometry: NH3_CH_COM_3A.xyz.',
    'These are fixed-geometry single-point calculations; no geometry optimization.',
    '',
    'Executable: /home/peter/Programs/MRCC/binary/dmrcc, version 26.1.3.',
    'Hardware: AMD Ryzen 9 7900, 12 physical cores, 24 logical CPUs.',
    'Each run used 1, 2 or 4 distinct physical cores: CPU IDs 0, 0-1 or 0-3.',
    'OMP_NUM_THREADS and MKL_NUM_THREADS were both set to the requested core count.',
    'Only one benchmark calculation ran at a time; no concurrent build or test.',
    'Settings: scftype=UHF, charge=0, mult=2, core=frozen, mem=2048MB,',
    'symm=off, unit=angs, molden=off, scfmaxit=200, ccmaxit=100.',
    'Standard MRCC SCF/CC convergence thresholds were retained.',
    '',
    'Timing methodology:',
    'Elapsed seconds are GNU time wall times for the complete dmrcc invocation,',
    'including SCF, integral generation, coupled cluster and triples correction.',
    'One fresh-directory measurement per method/basis/core configuration.',
    'The first cc-pVTZ one-core run began before the user stopped other heavy',
    'jobs. Its timing is excluded and replaced by a clean run after the series.',
    'The user confirmed other heavy jobs stopped before the two-core run began.',
    'Desktop applications remained open; this was not an isolated OS benchmark.',
    'Single measurements do not establish statistical timing uncertainty.',
    '',
    'F12 basis choices:',
    'The unrestricted F12 route uses calc=CCSD(T)-F12 with scftype=UHF.',
    'The nonaugmented F12 test uses cc-pVTZ-F12, the corresponding F12 orbital',
    'basis, rather than conventional cc-pVTZ. The augmented test uses',
    'aug-cc-pVTZ. MRCC documents automatic auxiliary/CABS selection for both',
    'cc-pVXZ-F12 and aug-cc-pVXZ; the installed basis library contains the',
    'matching OPTRI sets. No aug-cc-pVTZ-F12 basis was invented or substituted.',
    'F12 energy reported below is MRCC Total CCSD(F12*)(T) energy [au],',
    'not the separately printed (T*) or (T+) alternatives.',
    '',
    'Elapsed wall times (seconds)',
    'Method          Basis                  1 core      2 cores      4 cores',
]
for (method, basis), group in groups.items():
    assert set(group) == {1, 2, 4}
    values = [f"{group[n]['wall_seconds']:.2f}" if group[n]['status'] == 'PASS' else 'FAILED' for n in (1, 2, 4)]
    lines.append(f'{method:<15} {basis:<19} {values[0]:>11} {values[1]:>12} {values[2]:>12}')
lines += ['', 'Energy consistency and speedup (one-core time / multi-core time):']
for (method, basis), group in groups.items():
    good = [r for r in group.values() if r['status'] == 'PASS']
    lines.append(f'{method} / {basis}: {len(good)}/3 successful runs.')
    for n in (1, 2, 4):
        r = group[n]
        energy = f"{r['energy_hartree']:.12f} Eh" if r['energy_hartree'] is not None else 'unavailable'
        lines.append(f"  {n} cores: {r['status']}, {energy}; wall {r['wall_seconds']:.2f} s; CPU user/system {r['user_cpu_seconds']}/{r['system_cpu_seconds']} s.")
    if len(good) == 3:
        spread = max(r['energy_hartree'] for r in good) - min(r['energy_hartree'] for r in good)
        lines.append(f'  Energy spread: {spread:.3e} Eh.')
        lines.append(f"  Speedup: 2 cores {group[1]['wall_seconds']/group[2]['wall_seconds']:.3f}x; 4 cores {group[1]['wall_seconds']/group[4]['wall_seconds']:.3f}x.")
lines += ['', 'Spin diagnostics (as printed by MRCC):']
for (method, basis), group in groups.items():
    lines.append(f"{method} / {basis}, clean one-core: {group[1]['spin_diagnostics']}")
lines += [
    '',
    'Validation and evidence:',
    'PASS requires exit code zero, normal termination, a finite final energy',
    'at the requested level, and no detected SCF/CC nonconvergence message.',
    'Every evidence subdirectory retains MINP, mrcc.out, mrcc.err and time.txt.',
    'results.csv and results.json contain all 13 runs; filter',
    'included_in_comparison=true to obtain the 12 clean comparison records.',
    'metadata.json records hardware, paths and construction settings.',
    'quiet_machine.json records the user shutdown confirmation.',
    'benchmark.py and repeat_clean_first.py retain the exact launch procedure.',
    'Scratch paths are recorded in results.json. No application source,',
    'installation settings or user environment were changed for this benchmark.',
    '',
    'All times above are seconds. Frozen-core results do not test core correlation.',
]
(root / 'report.txt').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines[lines.index('Elapsed wall times (seconds)'):lines.index('Energy consistency and speedup (one-core time / multi-core time):')]))
print('Report:', root / 'report.txt')
