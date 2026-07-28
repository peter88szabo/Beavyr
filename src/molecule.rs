// src/molecule.rs

use bevy::prelude::*;
use std::collections::HashMap;

/// A lightweight snapshot of molecular coordinates (Å).
/// This is intentionally just a clone of `pos` so it’s cheap to create/apply.
#[allow(dead_code)]
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
    #[allow(dead_code)]
    pub fn snapshot(&self) -> MolSnapshot {
        self.pos.clone()
    }

    /// Apply a previously captured coordinate snapshot (Å).
    #[allow(dead_code)]
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

        self.recompute_covalent_bonds(thresh_scale);

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

    fn recompute_covalent_bonds(&mut self, thresh_scale: f32) {
        const SPATIAL_GRID_ATOM_THRESHOLD: usize = 256;

        let n = self.atoms.len();
        if n < SPATIAL_GRID_ATOM_THRESHOLD {
            for i in 0..n {
                for j in (i + 1)..n {
                    self.push_bond_if_close(i, j, thresh_scale);
                }
            }
            return;
        }

        let max_radius = self
            .atoms
            .iter()
            .map(|sym| covalent_radius_angstrom(sym))
            .fold(0.0_f32, f32::max);
        let cell_size = (max_radius * 2.0 * thresh_scale).max(0.5);

        let mut cells: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
        for (i, &p) in self.pos.iter().enumerate() {
            cells.entry(grid_cell(p, cell_size)).or_default().push(i);
        }

        for (&cell, indices) in &cells {
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let neighbor_cell = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
                        let Some(other_indices) = cells.get(&neighbor_cell) else {
                            continue;
                        };

                        for &i in indices {
                            for &j in other_indices {
                                if j <= i {
                                    continue;
                                }
                                self.push_bond_if_close(i, j, thresh_scale);
                            }
                        }
                    }
                }
            }
        }
    }

    fn push_bond_if_close(&mut self, i: usize, j: usize, thresh_scale: f32) {
        let ri = covalent_radius_angstrom(&self.atoms[i]);
        let rj = covalent_radius_angstrom(&self.atoms[j]);
        let cutoff = (ri + rj) * thresh_scale;
        let d = self.pos[i].distance(self.pos[j]);
        if d > 0.01 && d <= cutoff {
            self.bonds.push((i, j, d));
        }
    }
}

fn grid_cell(p: Vec3, cell_size: f32) -> (i32, i32, i32) {
    (
        (p.x / cell_size).floor() as i32,
        (p.y / cell_size).floor() as i32,
        (p.z / cell_size).floor() as i32,
    )
}

// ---- Parsers (Ångström) ----

pub fn parse_xyz_angstrom(xyz: &str) -> (usize, Vec<String>, Vec<Vec<f64>>) {
    let mut lines: Vec<&str> = xyz.lines().filter(|l| !l.trim().is_empty()).collect();
    if !lines.is_empty() {
        let first_parts: Vec<&str> = lines[0].split_whitespace().collect();
        if first_parts.len() == 1 && first_parts[0].parse::<usize>().is_ok() {
            // Support standard XYZ header: atom count + comment line.
            lines.remove(0);
            if !lines.is_empty() {
                let second_parts: Vec<&str> = lines[0].split_whitespace().collect();
                if second_parts.len() < 4 {
                    lines.remove(0);
                }
            }
        }
    }
    let n = lines.len();
    let mut atoms = Vec::new();
    let mut qxyz = Vec::with_capacity(n);

    for line in lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let atom = parts[0].to_string();
        let Ok(x) = parts[1].parse::<f64>() else {
            continue;
        };
        let Ok(y) = parts[2].parse::<f64>() else {
            continue;
        };
        let Ok(z) = parts[3].parse::<f64>() else {
            continue;
        };
        atoms.push(atom);
        qxyz.push(vec![x, y, z]);
    }

    (n, atoms, qxyz)
}

/// Parse only the first XYZ frame (Å). Returns (atoms, coords, extra_frames).
pub fn parse_xyz_first_frame_angstrom(xyz: &str) -> (Vec<String>, Vec<Vec<f64>>, bool) {
    let lines: Vec<&str> = xyz.lines().filter(|l| !l.trim().is_empty()).collect();
    if lines.is_empty() {
        return (Vec::new(), Vec::new(), false);
    }

    let first_parts: Vec<&str> = lines[0].split_whitespace().collect();
    let has_count = first_parts.len() == 1 && first_parts[0].parse::<usize>().is_ok();
    if !has_count {
        let (_n, atoms, coords) = parse_xyz_angstrom(xyz);
        return (atoms, coords, false);
    }

    let count = first_parts[0].parse::<usize>().unwrap_or(0);
    let mut idx = 1;
    if idx < lines.len() {
        let second_parts: Vec<&str> = lines[idx].split_whitespace().collect();
        if second_parts.len() < 4 {
            idx += 1;
        }
    }

    let mut atoms = Vec::new();
    let mut qxyz = Vec::new();
    let mut consumed = idx;
    for line in lines.iter().skip(idx) {
        consumed += 1;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let atom = parts[0].to_string();
        let Ok(x) = parts[1].parse::<f64>() else {
            continue;
        };
        let Ok(y) = parts[2].parse::<f64>() else {
            continue;
        };
        let Ok(z) = parts[3].parse::<f64>() else {
            continue;
        };
        atoms.push(atom);
        qxyz.push(vec![x, y, z]);
        if atoms.len() >= count {
            break;
        }
    }

    let extra_frames = lines.iter().skip(consumed).any(|l| !l.trim().is_empty());

    (atoms, qxyz, extra_frames)
}

/// Covalent radii in Å (no Bohr conversion anymore)
pub fn covalent_radius_angstrom(sym: &str) -> f32 {
    match sym {
        "H" => 0.45,
        "He" => 0.28,
        "Li" => 1.28,
        "Be" => 0.96,
        "B" => 0.84,
        "C" => 0.75,
        "N" => 0.70,
        "O" => 0.70,
        "F" => 0.57,
        "Ne" => 0.58,
        "Na" => 1.66,
        "Mg" => 1.41,
        "Al" => 1.21,
        "Si" => 1.11,
        "P" => 1.07,
        "S" => 1.05,
        "Cl" => 1.02,
        "Ar" => 1.06,
        "K" => 2.03,
        "Ca" => 1.76,
        "Sc" => 1.70,
        "Ti" => 1.60,
        "V" => 1.53,
        "Cr" => 1.39,
        "Mn" => 1.39,
        "Fe" => 1.32,
        "Co" => 1.26,
        "Ni" => 1.24,
        "Cu" => 1.32,
        "Zn" => 1.22,
        "Ga" => 1.22,
        "Ge" => 1.20,
        "As" => 1.19,
        "Se" => 1.20,
        "Br" => 1.20,
        "Kr" => 1.16,
        "Rb" => 2.20,
        "Sr" => 1.95,
        "Y" => 1.90,
        "Zr" => 1.75,
        "Nb" => 1.64,
        "Mo" => 1.54,
        "Tc" => 1.47,
        "Ru" => 1.46,
        "Rh" => 1.42,
        "Pd" => 1.39,
        "Ag" => 1.45,
        "Cd" => 1.44,
        "In" => 1.42,
        "Sn" => 1.39,
        "Sb" => 1.39,
        "Te" => 1.38,
        "I" => 1.39,
        "Xe" => 1.40,
        "Cs" => 2.44,
        "Ba" => 2.15,
        "La" => 2.07,
        "Ce" => 2.04,
        "Pr" => 2.03,
        "Nd" => 2.01,
        "Pm" => 1.99,
        "Sm" => 1.98,
        "Eu" => 1.98,
        "Gd" => 1.96,
        "Tb" => 1.94,
        "Dy" => 1.92,
        "Ho" => 1.92,
        "Er" => 1.89,
        "Tm" => 1.90,
        "Yb" => 1.87,
        "Lu" => 1.87,
        "Hf" => 1.75,
        "Ta" => 1.70,
        "W" => 1.62,
        "Re" => 1.51,
        "Os" => 1.44,
        "Ir" => 1.41,
        "Pt" => 1.36,
        "Au" => 1.36,
        "Hg" => 1.32,
        "Tl" => 1.45,
        "Pb" => 1.46,
        "Bi" => 1.48,
        "Po" => 1.40,
        "At" => 1.50,
        "Rn" => 1.50,
        "Fr" => 2.60,
        "Ra" => 2.21,
        "Ac" => 2.15,
        "Th" => 2.06,
        "Pa" => 2.00,
        "U" => 1.96,
        "Np" => 1.90,
        "Pu" => 1.87,
        "Am" => 1.80,
        "Cm" => 1.69,
        "Bk" => 1.68,
        "Cf" => 1.68,
        "Es" => 1.65,
        "Fm" => 1.67,
        "Md" => 1.73,
        "No" => 1.76,
        "Lr" => 1.61,
        _ => 0.77,
    }
}
