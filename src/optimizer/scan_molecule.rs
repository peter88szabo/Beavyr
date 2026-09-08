//! The slice of Behemoth's `Molecule` that its scan code needs.
//!
//! Behemoth's scan functions take its own `Molecule`, which carries a basis set, a method, an SCF
//! state and much else. They use almost none of it: `natoms`, `xyz`, and the `coord` mirror on
//! each atom. Rather than drag that type across -- it is bound to the quantum-chemistry engine --
//! this provides the same field names over Beavyr's molecule, so the imported scan code compiles
//! unchanged.
//!
//! **Coordinates are in bohr**, as Behemoth keeps them, not Å. [`ScanMolecule::from_molecule`] and
//! [`ScanMolecule::positions_angstrom`] convert at the boundary.
//!
//! One precision note: Beavyr stores positions as `f32`, because that is what the renderer wants,
//! so a geometry entering here carries about seven significant figures however precisely it was
//! written in the file. Everything downstream -- the scan, the optimiser, the force field -- then
//! works in `f64`, so this limits only the starting point, not the convergence. It is well below
//! any chemically meaningful displacement, but it does mean an exact round-trip through this type
//! is not bit-identical.

use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
use crate::molecule::Molecule;

/// One atom, with the `coord` field the scan code assigns through.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanAtom {
    pub symbol: String,
    /// bohr
    pub coord: [f64; 3],
}

/// A geometry the scan code can drive, in bohr.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanMolecule {
    pub natoms: usize,
    /// bohr
    pub xyz: Vec<[f64; 3]>,
    pub atom: Vec<ScanAtom>,
}

impl ScanMolecule {
    /// Converts from Beavyr's molecule, Å to bohr.
    pub fn from_molecule(mol: &Molecule) -> Self {
        let xyz: Vec<[f64; 3]> = mol
            .pos
            .iter()
            .map(|p| {
                [
                    p.x as f64 / BOHR_TO_ANGSTROM,
                    p.y as f64 / BOHR_TO_ANGSTROM,
                    p.z as f64 / BOHR_TO_ANGSTROM,
                ]
            })
            .collect();
        Self {
            natoms: mol.atoms.len(),
            atom: mol
                .atoms
                .iter()
                .zip(&xyz)
                .map(|(symbol, coord)| ScanAtom {
                    symbol: symbol.clone(),
                    coord: *coord,
                })
                .collect(),
            xyz,
        }
    }

    /// Parses an XYZ block, in Å, into bohr. Used by the imported scan tests.
    pub fn from_xyz(xyz: &str) -> Self {
        Self::from_molecule(&Molecule::from_xyz(xyz))
    }

    /// Element symbols, in atom order.
    pub fn symbols(&self) -> Vec<&str> {
        self.atom.iter().map(|a| a.symbol.as_str()).collect()
    }

    /// A flat `[x0, y0, z0, ...]` in Å, which is what the force field and the renderer want.
    pub fn positions_angstrom(&self) -> Vec<f64> {
        self.xyz
            .iter()
            .flat_map(|c| {
                [
                    c[0] * BOHR_TO_ANGSTROM,
                    c[1] * BOHR_TO_ANGSTROM,
                    c[2] * BOHR_TO_ANGSTROM,
                ]
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn converts_angstrom_to_bohr_and_back() {
        let mol = ScanMolecule::from_xyz("2\n\nH 0.000 0.000 0.000\nH 0.740 0.000 0.000\n");
        assert_eq!(mol.natoms, 2);
        // 0.74 Å is about 1.3984 bohr. The tolerance is set by Beavyr storing positions as f32,
        // not by the conversion, which is exact in f64.
        assert!((mol.xyz[1][0] - 0.74 / BOHR_TO_ANGSTROM).abs() < 1.0e-6);

        let back = mol.positions_angstrom();
        assert!((back[3] - 0.74).abs() < 1.0e-6);
        assert_eq!(mol.symbols(), vec!["H", "H"]);
    }

    /// The `coord` mirror the scan code writes through must start consistent with `xyz`.
    #[test]
    fn atom_coordinates_mirror_the_xyz_array() {
        let mol = ScanMolecule::from_xyz("3\n\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n");
        for (atom, xyz) in mol.atom.iter().zip(&mol.xyz) {
            assert_eq!(&atom.coord, xyz);
        }
    }
}
