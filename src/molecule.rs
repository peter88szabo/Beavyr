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
}

impl Molecule {
    /// No atoms at all -- what the program starts with, so it opens on an
    /// empty viewport rather than a structure nobody asked for.
    pub fn empty() -> Self {
        Molecule {
            atoms: Vec::new(),
            pos: Vec::new(),
            bonds: Vec::new(),
            hydrogen_bonds: Vec::new(),
        }
    }

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

        // 2) H-bonds (distance-only), store (H, A, H···A) and SKIP acceptors that are donor-neighbors.
        //
        // Only acceptor-capable heavy atoms can ever match, so the grid is built over
        // that subset alone.  Querying it per hydrogen turns the old O(n²) sweep into
        // O(n · neighbours), which is what keeps trajectory playback smooth on large
        // systems (the whole bond graph is recomputed per frame when `fixed_bonds`
        // is off).
        let acceptors: Vec<usize> = (0..n)
            .filter(|&a| self.atoms[a] != "H" && is_acceptor(&self.atoms[a]))
            .collect();
        if acceptors.is_empty() {
            return;
        }

        let acceptor_grid = SpatialGrid::build(acceptors.iter().copied(), &self.pos, ha_cut);

        // Reused across hydrogens so the per-atom query allocates nothing.
        let mut candidates: Vec<usize> = Vec::new();

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

            let p_h = self.pos[h];
            let p_d = self.pos[d];

            candidates.clear();
            acceptor_grid.collect_neighbors(p_h, &mut candidates);

            for &a in &candidates {
                if a == h || a == d {
                    continue;
                }

                // exclude nearest-neighbor heavy atoms: A must NOT be covalently bonded to D
                if bonded_to[d].contains(&a) {
                    continue;
                }

                // distance checks
                let ha = p_h.distance(self.pos[a]); // H···A
                if ha > ha_cut {
                    continue;
                }
                let da = p_d.distance(self.pos[a]); // D···A
                if da > da_cut {
                    continue;
                }

                // Each acceptor lives in exactly one cell and the 27 queried cells are
                // distinct, so a given (h, a) pair is produced at most once.
                self.hydrogen_bonds.push((h, a, ha));
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

        let grid = SpatialGrid::build(0..n, &self.pos, cell_size);

        let mut candidates: Vec<usize> = Vec::new();
        for i in 0..n {
            let p = self.pos[i];
            candidates.clear();
            grid.collect_neighbors(p, &mut candidates);
            for &j in &candidates {
                // Visit each unordered pair once.
                if j <= i {
                    continue;
                }
                self.push_bond_if_close(i, j, thresh_scale);
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

/// Uniform spatial hash over a subset of atom positions.
///
/// `cell_size` must be >= the query radius, so that every point within that
/// radius of a query position lands in one of the 27 cells around it.
pub struct SpatialGrid {
    cell_size: f32,
    cells: HashMap<(i32, i32, i32), Vec<usize>>,
}

impl SpatialGrid {
    pub fn build(indices: impl Iterator<Item = usize>, pos: &[Vec3], cell_size: f32) -> Self {
        let cell_size = cell_size.max(0.5);
        let mut cells: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
        for i in indices {
            if let Some(&p) = pos.get(i) {
                cells.entry(grid_cell(p, cell_size)).or_default().push(i);
            }
        }
        Self { cell_size, cells }
    }

    /// Append every indexed point in the 27 cells surrounding `p` into `out`.
    /// `out` is caller-owned so the query allocates nothing in a hot loop.
    pub fn collect_neighbors(&self, p: Vec3, out: &mut Vec<usize>) {
        let (cx, cy, cz) = grid_cell(p, self.cell_size);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(bucket) = self.cells.get(&(cx + dx, cy + dy, cz + dz)) {
                        out.extend_from_slice(bucket);
                    }
                }
            }
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

/// Atomic number, 1-based, for the same element set `covalent_radius_angstrom`
/// covers -- generated from that function's own ordering so the two tables
/// cannot silently drift apart.
const ELEMENT_SYMBOLS: [&str; 103] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"
];

pub fn atomic_number(sym: &str) -> Option<u32> {
    ELEMENT_SYMBOLS
        .iter()
        .position(|&s| s == sym)
        .map(|i| i as u32 + 1)
}

/// The inverse of `atomic_number`: element symbol for a given atomic number
/// (1-based), same table so the two can never drift apart.
pub fn element_symbol(z: u32) -> Option<&'static str> {
    ELEMENT_SYMBOLS.get((z as usize).checked_sub(1)?).copied()
}

/// Standard atomic weight in g/mol (amu), same element set and ordering as
/// `atomic_number` -- generated the same way so the two tables cannot
/// silently drift apart. IUPAC 2021 standard atomic weights; for elements
/// with no stable isotope, the mass of the longest-lived known isotope
/// (the usual convention, e.g. NIST's own tables).
pub fn atomic_mass_amu(sym: &str) -> Option<f64> {
    const MASSES: [f64; 103] = [
        1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.18, 22.99, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.63, 74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 88.906, 91.224, 92.906, 95.95, 97.0, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.6, 126.9, 131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24, 145.0, 150.36, 151.96, 157.25, 158.93, 162.5, 164.93, 167.26, 168.93, 173.05, 174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, 204.38, 207.2, 208.98, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.04, 231.04, 238.03, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0, 258.0, 259.0, 262.0
    ];
    atomic_number(sym).map(|z| MASSES[(z - 1) as usize])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The program opens with no structure at all, so every consumer of the
    /// molecule has to cope with zero atoms from the first frame onward.
    #[test]
    fn an_empty_molecule_has_nothing_in_it() {
        let mol = Molecule::empty();
        assert!(mol.atoms.is_empty());
        assert!(mol.pos.is_empty());
        assert!(mol.bonds.is_empty());
        assert!(mol.hydrogen_bonds.is_empty());
    }

    #[test]
    fn perceiving_bonds_on_an_empty_molecule_is_a_no_op() {
        let mut mol = Molecule::empty();
        mol.recompute_bonds(1.2, 2.5);
        assert!(mol.bonds.is_empty());
    }

    /// Loading empty text is the same as starting empty -- the XYZ box begins
    /// blank, and applying it must not be a special case.
    #[test]
    fn empty_xyz_text_parses_to_an_empty_molecule() {
        let mol = Molecule::from_xyz("");
        assert!(mol.atoms.is_empty());
        assert!(mol.pos.is_empty());
    }

    #[test]
    fn atomic_mass_matches_well_known_values() {
        assert!((atomic_mass_amu("H").unwrap() - 1.008).abs() < 1e-6);
        assert!((atomic_mass_amu("C").unwrap() - 12.011).abs() < 1e-6);
        assert!((atomic_mass_amu("O").unwrap() - 15.999).abs() < 1e-6);
        assert_eq!(atomic_mass_amu("Xx"), None);
    }


    /// Deterministic LCG so the fixtures are reproducible without a rand dep.
    struct Lcg(u64);
    impl Lcg {
        fn next_f32(&mut self) -> f32 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((self.0 >> 33) as f32) / ((1u64 << 31) as f32)
        }
    }

    fn fixture(n: usize, box_size: f32, seed: u64) -> (Vec<String>, Vec<Vec3>) {
        let symbols = ["H", "H", "O", "N", "C", "S", "Cl"];
        let mut rng = Lcg(seed);
        let mut atoms = Vec::with_capacity(n);
        let mut pos = Vec::with_capacity(n);
        for i in 0..n {
            atoms.push(symbols[i % symbols.len()].to_string());
            pos.push(Vec3::new(
                rng.next_f32() * box_size,
                rng.next_f32() * box_size,
                rng.next_f32() * box_size,
            ));
        }
        (atoms, pos)
    }

    /// Brute-force reference for covalent bonds.
    fn covalent_reference(
        atoms: &[String],
        pos: &[Vec3],
        thresh_scale: f32,
    ) -> Vec<(usize, usize)> {
        let n = atoms.len();
        let mut out = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let cutoff = (covalent_radius_angstrom(&atoms[i])
                    + covalent_radius_angstrom(&atoms[j]))
                    * thresh_scale;
                let d = pos[i].distance(pos[j]);
                if d > 0.01 && d <= cutoff {
                    out.push((i, j));
                }
            }
        }
        out.sort_unstable();
        out
    }

    /// Brute-force reference for hydrogen bonds, mirroring the pre-grid logic.
    fn hbond_reference(
        atoms: &[String],
        pos: &[Vec3],
        bonds: &[(usize, usize, f32)],
        hbond_cutoff: f32,
    ) -> Vec<(usize, usize)> {
        let is_donor = |s: &str| matches!(s, "O" | "N" | "S" | "F" | "P");
        let is_acceptor = |s: &str| matches!(s, "O" | "N" | "S" | "F" | "Cl" | "Br" | "I" | "P");
        let ha_cut = hbond_cutoff.max(1.2);
        let da_cut = (hbond_cutoff + 1.0).min(4.0);

        let n = atoms.len();
        let mut bonded_to: Vec<Vec<usize>> = vec![Vec::new(); n];
        for &(i, j, _) in bonds {
            bonded_to[i].push(j);
            bonded_to[j].push(i);
        }

        let mut out = Vec::new();
        for h in 0..n {
            if atoms[h] != "H" {
                continue;
            }
            let Some(&d) = bonded_to[h]
                .iter()
                .find(|&&nb| atoms[nb] != "H" && is_donor(&atoms[nb]))
            else {
                continue;
            };
            for a in 0..n {
                if a == h || a == d || atoms[a] == "H" || !is_acceptor(&atoms[a]) {
                    continue;
                }
                if bonded_to[d].contains(&a) {
                    continue;
                }
                if pos[h].distance(pos[a]) > ha_cut {
                    continue;
                }
                if pos[d].distance(pos[a]) > da_cut {
                    continue;
                }
                out.push((h, a));
            }
        }
        out.sort_unstable();
        out
    }

    fn build(atoms: Vec<String>, pos: Vec<Vec3>) -> Molecule {
        Molecule {
            atoms,
            pos,
            bonds: vec![],
            hydrogen_bonds: vec![],
        }
    }

    #[test]
    fn covalent_grid_matches_brute_force_above_threshold() {
        // 400 atoms forces the spatial-grid path (threshold is 256).
        let (atoms, pos) = fixture(400, 14.0, 0x5EED);
        let thresh = 1.2;
        let expected = covalent_reference(&atoms, &pos, thresh);

        let mut mol = build(atoms, pos);
        mol.recompute_bonds(thresh, 2.5);

        let mut got: Vec<(usize, usize)> = mol.bonds.iter().map(|&(i, j, _)| (i, j)).collect();
        got.sort_unstable();

        assert!(!expected.is_empty(), "fixture produced no bonds to compare");
        assert_eq!(got, expected);
    }

    #[test]
    fn covalent_bonds_are_not_duplicated() {
        let (atoms, pos) = fixture(400, 14.0, 0xC0FFEE);
        let mut mol = build(atoms, pos);
        mol.recompute_bonds(1.2, 2.5);

        let mut pairs: Vec<(usize, usize)> = mol.bonds.iter().map(|&(i, j, _)| (i, j)).collect();
        let total = pairs.len();
        pairs.sort_unstable();
        pairs.dedup();
        assert_eq!(total, pairs.len(), "grid emitted duplicate bonds");
    }

    #[test]
    fn hbond_grid_matches_brute_force() {
        for &(n, box_size, seed) in &[(120usize, 9.0f32, 1u64), (400, 14.0, 2), (700, 17.0, 3)] {
            let (atoms, pos) = fixture(n, box_size, seed);
            let cutoff = 2.5;

            let mut mol = build(atoms.clone(), pos.clone());
            mol.recompute_bonds(1.2, cutoff);

            let expected = hbond_reference(&atoms, &pos, &mol.bonds, cutoff);
            let mut got: Vec<(usize, usize)> =
                mol.hydrogen_bonds.iter().map(|&(h, a, _)| (h, a)).collect();
            got.sort_unstable();

            assert_eq!(got, expected, "mismatch for n={n}");
        }
    }

    #[test]
    fn hbond_grid_finds_a_known_water_dimer_contact() {
        let atoms: Vec<String> = ["O", "H", "H", "O", "H", "H"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let pos = vec![
            // Donor water: one O-H points straight down at the acceptor.
            Vec3::new(0.000, 0.000, 0.000),
            Vec3::new(0.000, -0.960, 0.000), // donor H, 1.94 Å from the acceptor O
            Vec3::new(0.930, 0.240, 0.000),
            // Acceptor water, 2.90 Å below the donor O.
            Vec3::new(0.000, -2.900, 0.000),
            Vec3::new(0.500, -3.400, 0.700),
            Vec3::new(-0.500, -3.400, -0.700),
        ];
        let mut mol = build(atoms, pos);
        mol.recompute_bonds(1.2, 2.5);

        assert!(
            mol.hydrogen_bonds.iter().any(|&(h, a, _)| h == 1 && a == 3),
            "expected an H···O contact into the second water, got {:?}",
            mol.hydrogen_bonds
        );
    }
}
