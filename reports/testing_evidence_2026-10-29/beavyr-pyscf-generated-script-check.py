import json
from pathlib import Path
import subprocess

import numpy as np
from pyscf import dft, gto, lib, scf

lib.num_threads(2)
root = Path('/tmp/beavyr-release-pyscf-generated')
mol = gto.M(atom='O 0 0 0; H 0 0 0.96; H 0.92 0 -0.24', basis='sto-3g', verbose=0)
reference = scf.RHF(mol).run()
assert reference.converged
reference_gradient = reference.nuc_grad_method().kernel()
reference_hessian = reference.Hessian().kernel().transpose(0, 2, 1, 3).reshape(9, 9)
results = []
for name in ['energy', 'gradient', 'hessian', 'dft_gradient', 'optimization_missing_dependency']:
    directory = root / name
    # This process was launched after the user's actual venv activation.
    result = subprocess.run(['python', 'run_pyscf.py'], cwd=directory, capture_output=True, text=True)
    (directory / 'stdout.txt').write_text(result.stdout)
    (directory / 'stderr.txt').write_text(result.stderr)
    if name == 'optimization_missing_dependency':
        assert result.returncode != 0, result.stdout
        assert "No module named 'geometric'" in result.stderr, result.stderr
        assert not (directory / 'final.xyz').exists()
        results.append({'job': name, 'result': 'Expected failure: missing geometric; no final geometry'})
        continue
    assert result.returncode == 0, result.stderr
    energy = float(next(line[7:] for line in result.stdout.splitlines() if line.startswith('energy ')))
    assert np.isfinite(energy)
    if name != 'dft_gradient':
        assert abs(energy - reference.e_tot) < 1e-10, (name, energy, reference.e_tot)
    geometry = (directory / 'final.xyz').read_text().splitlines()
    assert int(geometry[0]) == 3
    xyz = np.array([[float(x) for x in row.split()[1:]] for row in geometry[2:]])
    assert np.max(np.abs(xyz - mol.atom_coords(unit='Angstrom'))) < 1e-8
    if name in ['gradient', 'dft_gradient']:
        rows = [row for row in (directory / 'pyscf.grad').read_text().splitlines()
                if row.strip() and not row.strip().startswith(('$', 'cycle'))]
        assert len(rows) == 6
        gradient = np.array([[float(x) for x in row.split()] for row in rows[3:]])
        expected = reference_gradient
        if name == 'dft_gradient':
            ks = dft.RKS(mol)
            ks.xc = 'b3lyp'
            ks.kernel()
            assert ks.converged
            assert abs(energy - ks.e_tot) < 1e-10
            expected = ks.nuc_grad_method().kernel()
        assert np.max(np.abs(gradient - expected)) < 1e-10
    if name == 'hessian':
        rows = [row for row in (directory / 'hessian.txt').read_text().splitlines()
                if row.strip() and not row.strip().startswith('$')]
        matrix = np.array([[float(x) for x in row.split()] for row in rows])
        assert matrix.shape == (9, 9)
        assert np.max(np.abs(matrix - reference_hessian)) < 1e-10
    results.append({'job': name, 'result': 'PASS', 'energy_hartree': energy})
print(json.dumps({'scope': 'Current production PySCF script generator compiled in a temporary adapter; full app dispatch excluded', 'checks': results}, indent=2))
