//! Text exports of the edited structure: XYZ and Z-matrix.

use super::zmat2xyz::xyz_to_zmat;
use crate::molecule::Molecule;

/// A standard XYZ file: atom count, comment line, then one atom per line (Å).
pub fn xyz_text(mol: &Molecule, comment: &str) -> String {
    let mut out = format!("{}\n{}\n", mol.atoms.len(), comment.trim());
    for (symbol, p) in mol.atoms.iter().zip(&mol.pos) {
        out.push_str(&format!(
            "{:<2} {:>14.8} {:>14.8} {:>14.8}\n",
            symbol, p.x as f64, p.y as f64, p.z as f64
        ));
    }
    out
}

/// A Gaussian-style Z-matrix built from the current coordinates: each line is
/// the element, then (reference atom, value) pairs for bond (Å), angle and
/// dihedral (degrees), with 1-based atom numbers.
pub fn zmat_text(mol: &Molecule) -> String {
    let mut out = String::new();
    for row in xyz_to_zmat(&mol.atoms, &mol.pos) {
        out.push_str(&row.symbol);
        if let Some(r) = row.bond_ref {
            out.push_str(&format!("  {r} {:.6}", row.bond_len));
        }
        if let Some(r) = row.angle_ref {
            out.push_str(&format!("  {r} {:.4}", row.angle_deg));
        }
        if let Some(r) = row.dihedral_ref {
            out.push_str(&format!("  {r} {:.4}", row.dihedral_deg));
        }
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule_builder::zmat2xyz::{zmat_to_xyz, ZAtom};

    fn water() -> Molecule {
        Molecule::from_xyz("3\n\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
    }

    #[test]
    fn xyz_has_count_and_comment_lines_and_reads_back() {
        let mol = water();
        let text = xyz_text(&mol, "water");
        let mut lines = text.lines();
        assert_eq!(lines.next(), Some("3"));
        assert_eq!(lines.next(), Some("water"));
        let back = Molecule::from_xyz(&text);
        assert_eq!(back.atoms, mol.atoms);
        for (a, b) in back.pos.iter().zip(&mol.pos) {
            assert!(a.distance(*b) < 1e-6);
        }
    }

    #[test]
    fn zmat_text_rebuilds_the_same_internal_geometry() {
        let mol = water();
        let text = zmat_text(&mol);
        let rows: Vec<ZAtom> = text
            .lines()
            .map(|line| {
                let f: Vec<&str> = line.split_whitespace().collect();
                let r = |i: usize| f.get(i).map(|v| v.parse::<usize>().unwrap());
                let v = |i: usize| f.get(i).map_or(0.0, |v| v.parse::<f64>().unwrap());
                ZAtom {
                    symbol: f[0].to_string(),
                    bond_ref: r(1),
                    bond_len: v(2),
                    angle_ref: r(3),
                    angle_deg: v(4),
                    dihedral_ref: r(5),
                    dihedral_deg: v(6),
                }
            })
            .collect();
        assert_eq!(text.lines().next(), Some("O"));
        let pos = zmat_to_xyz(&rows);
        let d = |p: &[bevy::prelude::Vec3], a: usize, b: usize| p[a].distance(p[b]);
        for (a, b) in [(0, 1), (0, 2), (1, 2)] {
            assert!((d(&pos, a, b) - d(&mol.pos, a, b)).abs() < 1e-4);
        }
    }
}
