// src/molecule.rs

use bevy::prelude::*;

/// A lightweight snapshot of molecular coordinates (Å).
/// This is intentionally just a clone of `pos` so it’s cheap to create/apply.
pub type MolSnapshot = Vec<Vec3>;

/// Positions are stored in Å (Angstrom). All parsing is Å.
#[derive(Resource, Clone)]
pub struct Molecule {
    pub atoms: Vec<String>,
    pub pos: Vec<Vec3>,                           // positions in Å
    pub bonds: Vec<(usize, usize, f32)>,          // (i, j, distance in Å)
    pub hydrogen_bonds: Vec<(usize, usize, f32)>, // (i, j, distance in Å)
    pub atom_entities: Vec<Entity>,
    pub bond_entities: Vec<Entity>,
}

impl Molecule {
    /// Build from simple XYZ text (Å).
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
            hydrogen_bonds: vec![],
            atom_entities: vec![],
            bond_entities: vec![],
        };
        mol.recompute_bonds(2.0, 3.0);
        mol
    }

    /// Replace atom positions with `new_pos` (Å).
    /// This function **only** updates `pos` — it does not recompute bonds or update transforms.
    /// Callers can decide when to trigger recompute/sync.
    pub fn set_pos(&mut self, new_pos: Vec<Vec3>) {
        self.pos = new_pos;
    }

    /// Create a snapshot of the current coordinates (Å).
    pub fn snapshot(&self) -> MolSnapshot {
        self.pos.clone()
    }

    /// Apply a previously captured coordinate snapshot (Å).
    pub fn apply_snapshot(&mut self, snap: &MolSnapshot) {
        self.pos = snap.clone();
    }

    /// Recompute covalent and hydrogen-bond lists based on current positions and settings.
    pub fn recompute_bonds(&mut self, thresh_scale: f32, hbond_cutoff: f32) {
        // Tunable element sets (extend/trim as you like)
        let is_donor = |s: &str| matches!(s, "O" | "N" | "S" | "F" | "P");
        let is_acceptor = |s: &str| matches!(s, "O" | "N" | "S" | "F" | "Cl" | "Br" | "I" | "P");

        let ha_cut = hbond_cutoff.max(1.2); // H···A cutoff (Å)
        let da_cut = (hbond_cutoff + 1.0).min(4.0); // D···A cutoff (Å), a bit looser

        self.bonds.clear();
        self.hydrogen_bonds.clear();

        let n = self.atoms.len();
        if n == 0 {
            return;
        }

        // 1) Covalent bonds (your rule)
        for i in 0..n {
            for j in (i + 1)..n {
                let ri = covalent_radius_angstrom(&self.atoms[i]);
                let rj = covalent_radius_angstrom(&self.atoms[j]);
                let cutoff = (ri + rj) * thresh_scale;
                let d = self.pos[i].distance(self.pos[j]);
                if d > 0.01 && d <= cutoff {
                    self.bonds.push((i, j, d));
                }
            }
        }

        // Build adjacency from covalent bonds
        let mut bonded_to: Vec<Vec<usize>> = vec![Vec::new(); n];
        for &(i, j, _) in &self.bonds {
            bonded_to[i].push(j);
            bonded_to[j].push(i);
        }

        // 2) H-bonds (distance-only), store (H, A, H···A) and SKIP acceptors that are donor-neighbors
        for h in 0..n {
            if self.atoms[h] != "H" {
                continue;
            }

            // donor D must be a heavy atom covalently bonded to this H
            let mut donor: Option<usize> = None;
            for &nb in &bonded_to[h] {
                if self.atoms[nb] != "H" && is_donor(&self.atoms[nb]) {
                    donor = Some(nb);
                    break;
                }
            }
            let Some(d) = donor else {
                continue;
            };

            for a in 0..n {
                if a == h || a == d {
                    continue;
                }
                if self.atoms[a] == "H" {
                    continue;
                }
                if !is_acceptor(&self.atoms[a]) {
                    continue;
                }

                // exclude nearest-neighbor heavy atoms: A must NOT be covalently bonded to D
                if bonded_to[d].contains(&a) {
                    continue;
                }

                // distance checks
                let ha = self.pos[h].distance(self.pos[a]); // H···A
                if ha > ha_cut {
                    continue;
                }
                let da = self.pos[d].distance(self.pos[a]); // D···A
                if da > da_cut {
                    continue;
                }

                // record H→A (avoid duplicates; keep shortest H···A)
                if let Some(ex) = self
                    .hydrogen_bonds
                    .iter_mut()
                    .find(|(hi, ai, _)| *hi == h && *ai == a)
                {
                    if ha < ex.2 {
                        ex.2 = ha;
                    }
                } else {
                    self.hydrogen_bonds.push((h, a, ha));
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
        "H" => 0.45,
        "C" => 0.75,
        "N" => 0.70,
        "O" => 0.70,
        "F" => 0.57,
        "S" => 1.05,
        "Cl" => 1.02,
        "Br" => 1.20,
        "I" => 1.39,
        _ => 0.77,
    }
}

