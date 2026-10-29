import importlib.util
import json
import sys

import numpy as np
import pyscf
from pyscf import gto, lib, scf

lib.num_threads(2)
mol = gto.M(
    atom="O 0 0 0; H 0 0 0.96; H 0.92 0 -0.24",
    basis="sto-3g", charge=0, spin=0, unit="Angstrom", verbose=0,
)

def solve(molecule):
    mf = scf.RHF(molecule)
    mf.conv_tol = 1e-12
    energy = mf.kernel()
    assert mf.converged and np.isfinite(energy)
    return mf, float(energy)

mf, energy = solve(mol)
gradient = mf.nuc_grad_method().kernel()
hessian = mf.Hessian().kernel().transpose(0, 2, 1, 3).reshape(9, 9)
assert -76 < energy < -74, energy
assert np.isfinite(gradient).all() and np.isfinite(hessian).all()
net_gradient = float(np.max(np.abs(gradient.sum(axis=0))))
asymmetry = float(np.max(np.abs(hessian - hessian.T)))
assert net_gradient < 1e-8, net_gradient
assert asymmetry < 1e-8, asymmetry

coords = mol.atom_coords()
displacement = np.zeros_like(coords)
displacement[2, 0] = 1e-4  # Bohr: the derivatives also use atomic units.
mf_plus, e_plus = solve(mol.copy().set_geom_(coords + displacement, unit="Bohr"))
mf_minus, e_minus = solve(mol.copy().set_geom_(coords - displacement, unit="Bohr"))
fd_gradient = (e_plus - e_minus) / 2e-4
gradient_error = abs(fd_gradient - gradient[2, 0])
fd_column = (mf_plus.nuc_grad_method().kernel() - mf_minus.nuc_grad_method().kernel()).ravel() / 2e-4
hessian_error = float(np.max(np.abs(fd_column - hessian[:, 6])))
assert gradient_error < 1e-6, gradient_error
assert hessian_error < 1e-5, hessian_error

print(json.dumps({
    "scope": "Native PySCF environment and numerical derivatives; not the Beavyr runner",
    "interpreter": sys.executable,
    "pyscf_version": pyscf.__version__,
    "method": "RHF/STO-3G water, charge 0, multiplicity 1",
    "energy_hartree": energy,
    "gradient_components": int(gradient.size),
    "hessian_shape": list(hessian.shape),
    "gradient_finite_difference_error": float(gradient_error),
    "hessian_finite_difference_error": hessian_error,
    "geometric_available": importlib.util.find_spec("geometric") is not None,
    "result": "PASS",
}, indent=2))
