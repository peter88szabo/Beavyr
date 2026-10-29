# Written by Beavyr. Runs one PySCF job and leaves its results beside this file.
import sys
import numpy
from pyscf import gto, lib

lib.num_threads(2)

ATOMS = ["O", "H", "H"]
GEOMETRY = "O 0.0000000000 0.0000000000 0.0000000000\nH 0.0000000000 0.0000000000 0.9600000000\nH 0.9200000000 0.0000000000 -0.2400000000\n"

mol = gto.M(atom=GEOMETRY, basis="sto-3g", charge=0, spin=0,
            unit="Angstrom", verbose=0, output=None)

BOHR = 0.52917721092

def build(molecule):
    """The mean-field object for this molecule, restricted or not by spin."""
    from pyscf import scf
    return scf.UHF(molecule) if molecule.spin else scf.RHF(molecule)


def write_xyz(path, molecule, comment):
    coords = molecule.atom_coords() * BOHR
    with open(path, "w") as handle:
        handle.write("%d\n%s\n" % (len(ATOMS), comment))
        for symbol, xyz in zip(ATOMS, coords):
            handle.write("%-3s %18.10f %18.10f %18.10f\n" % (symbol, xyz[0], xyz[1], xyz[2]))

mf = build(mol)
energy = float(mf.kernel())
if not mf.converged:
    sys.stderr.write("PySCF SCF did not converge.\n")
    sys.exit(2)
if not numpy.isfinite(energy):
    sys.stderr.write("PySCF returned a non-finite energy: %r\n" % (energy,))
    sys.exit(2)
write_xyz("final.xyz", mol, "PySCF geometry")

g = mf.nuc_grad_method().kernel()
coords = mol.atom_coords()
with open("pyscf.grad", "w") as handle:
    handle.write("$grad          cartesian gradients\n")
    handle.write("  cycle =      1    SCF energy =  %20.14f   |dE/dxyz| =  0.000000\n" % energy)
    for symbol, xyz in zip(ATOMS, coords):
        handle.write("%22.14f %22.14f %22.14f      %s\n" % (xyz[0], xyz[1], xyz[2], symbol))
    for row in g:
        handle.write("%22.14e %22.14e %22.14e\n" % (row[0], row[1], row[2]))
    handle.write("$end\n")

print("energy %.12f" % energy)
print("natoms %d" % mol.natm)
