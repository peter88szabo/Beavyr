use bevy::prelude::*;

/// Positions are stored in Å (Angstrom). All parsing is Å.
#[derive(Resource, Clone)]
pub struct Molecule {
    pub atoms: Vec<String>,
    pub pos: Vec<Vec3>,                  // positions in Å
    pub bonds: Vec<(usize, usize, f32)>, // (i, j, distance in Å)
    pub atom_entities: Vec<Entity>,
    pub bond_entities: Vec<Entity>,
}

impl Molecule {
    pub fn from_xyz(xyz: &str) -> Self {
        let (_n, atoms, qxyz) = parse_xyz_angstrom(xyz);
        let pos: Vec<Vec3> = qxyz
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();

        let mut mol = Molecule {
            atoms,
            pos,
            bonds: vec![],
            atom_entities: vec![],
            bond_entities: vec![],
        };
        mol.recompute_bonds(2.0);
        mol
    }

    pub fn recompute_bonds(&mut self, thresh_scale: f32) {
        self.bonds.clear();
        let n = self.atoms.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let ri = covalent_radius_angstrom(&self.atoms[i]);
                let rj = covalent_radius_angstrom(&self.atoms[j]);
                let cutoff = (ri + rj) * thresh_scale;
                let d = self.pos[i].distance(self.pos[j]);
                if d > 0.0001 && d <= cutoff {
                    self.bonds.push((i, j, d));
                }
            }
        }
    }
}

// ---- Parsers (Ångström) ----

pub fn parse_xyz_angstrom(xyz: &str) -> (usize, Vec<String>, Vec<Vec<f64>>) {
    let lines: Vec<&str> = xyz.lines().filter(|l| !l.trim().is_empty()).collect();
    let n = lines.len();
    let mut atoms = Vec::new();
    let mut qxyz = Vec::with_capacity(n);

    for line in lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        let atom = parts[0].to_string();
        let x = parts[1].parse::<f64>().unwrap();
        let y = parts[2].parse::<f64>().unwrap();
        let z = parts[3].parse::<f64>().unwrap();
        atoms.push(atom);
        qxyz.push(vec![x, y, z]);
    }

    (n, atoms, qxyz)
}

/// Covalent radii in Å (no Bohr conversion anymore)
pub fn covalent_radius_angstrom(sym: &str) -> f32 {
    match sym {
        "H"  => 0.45,
        "C"  => 0.75,
        "N"  => 0.70,
        "O"  => 0.70,
        "F"  => 0.57,
        "S"  => 1.05,
        "Cl" => 1.02,
        "Br" => 1.20,
        "I"  => 1.39,
        _    => 0.77,
    }
}

