// src/forcefields/dreiding.rs
#![allow(clippy::too_many_arguments)]

use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;

/// Simple 3D vector with minimal ops (no dependencies).
#[derive(Clone, Copy, Debug, Default)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
impl Vec3 {
    #[inline] pub fn new(x: f64, y: f64, z: f64) -> Self { Self { x, y, z } }
    #[inline] pub fn sub(self, o: Self) -> Self { Self::new(self.x - o.x, self.y - o.y, self.z - o.z) }
    #[inline] pub fn dot(self, o: Self) -> f64 { self.x * o.x + self.y * o.y + self.z * o.z }
    #[inline] pub fn cross(self, o: Self) -> Self {
        Self::new(self.y * o.z - self.z * o.y, self.z * o.x - self.x * o.z, self.x * o.y - self.y * o.x)
    }
    #[inline] pub fn norm(self) -> f64 { self.dot(self).sqrt() }
    #[inline] pub fn normalized(self) -> Self {
        let n = self.norm();
        if n == 0.0 { self } else { Self::new(self.x / n, self.y / n, self.z / n) }
    }
}

/// Topology
#[derive(Clone, Debug)]
pub struct Atom {
    pub typ: String,   // DREIDING type (e.g., "C_R", "N_3", "H")
    pub charge: f64,   // fixed partial charge
    pub pos: Vec3,     // Cartesian coordinates (Å)
}

#[derive(Clone, Copy, Debug)]
pub struct Bond { pub i: usize, pub j: usize }

#[derive(Clone, Copy, Debug)]
pub struct Angle { pub i: usize, pub j: usize, pub k: usize } // j is vertex

#[derive(Clone, Copy, Debug)]
pub struct Dihedral { pub i: usize, pub j: usize, pub k: usize, pub l: usize }

/// Improper: j is the central atom; (i, k, l) define the plane around j.
/// We compute χ as the dihedral(i, j, k, l) and use a harmonic about 0.
#[derive(Clone, Copy, Debug)]
pub struct Improper { pub i: usize, pub j: usize, pub k: usize, pub l: usize }

/// DREIDING parameter records
#[derive(Clone, Copy, Debug)]
pub struct BondParams { pub k_r: f64, pub r0: f64 }              // (kcal/mol)/Å^2, Å
#[derive(Clone, Copy, Debug)]
pub struct AngleParams { pub k_th: f64, pub theta0_deg: f64 }    // (kcal/mol)/rad^2, degrees
#[derive(Clone, Copy, Debug)]
pub struct TorsionParams { pub v: f64, pub n: i32 }              // kcal/mol, multiplicity
#[derive(Clone, Copy, Debug)]
pub struct ImproperParams { pub k_chi: f64 }                     // (kcal/mol)/rad^2
#[derive(Clone, Copy, Debug)]
pub struct VdwParams { pub sigma: f64, pub eps: f64 }            // Å, kcal/mol

/// The parameter database interface you plug your tables into.
pub trait ParameterDB {
    /// Return bond parameters for the (type_i, type_j) pair (order-insensitive).
    fn bond_params(&self, ti: &str, tj: &str, bond: Bond) -> Option<BondParams>;

    /// Return angle params for (type_i, type_j, type_k) with j as vertex.
    fn angle_params(&self, ti: &str, tj: &str, tk: &str, ang: Angle) -> Option<AngleParams>;

    /// Return torsion params for (i, j, k, l) types. DREIDING often uses a single (v, n).
    fn torsion_params(&self, ti: &str, tj: &str, tk: &str, tl: &str, dih: Dihedral) -> Option<TorsionParams>;

    /// Return improper (out-of-plane) params centered at j.
    fn improper_params(&self, ti: &str, tj: &str, tk: &str, tl: &str, imp: Improper) -> Option<ImproperParams>;

    /// Per-atom vdW params (σ, ε). Mixing done here or by caller.
    fn vdw_params(&self, t: &str, idx: usize) -> Option<VdwParams>;
}

/// Settings for exclusions and 1–4 scaling (common in many FFs; keep =1.0 if you don’t want scaling).
#[derive(Clone, Copy, Debug)]
pub struct NonbondedScaling {
    pub scale_14_coul: f64, // e.g., 0.5 in some FFs; DREIDING commonly 1.0
    pub scale_14_vdw:  f64, // e.g., 0.5 in some FFs; DREIDING commonly 1.0
}
impl Default for NonbondedScaling {
    fn default() -> Self { Self { scale_14_coul: 1.0, scale_14_vdw: 1.0 } }
}

/// System container
pub struct DreidingSystem<P: ParameterDB> {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub angles: Vec<Angle>,
    pub dihedrals: Vec<Dihedral>,
    pub impropers: Vec<Improper>,
    pub params: P,
    /// If you already have a precomputed list of 1–4 pairs, provide it here. Otherwise it’s inferred.
    pub pairs_14: Option<HashSet<(usize, usize)>>,
    pub scaling: NonbondedScaling,
}

impl<P: ParameterDB> DreidingSystem<P> {
    pub fn new(atoms: Vec<Atom>, bonds: Vec<Bond>, angles: Vec<Angle>, dihedrals: Vec<Dihedral>, impropers: Vec<Improper>, params: P) -> Self {
        Self {
            atoms,
            bonds,
            angles,
            dihedrals,
            impropers,
            params,
            pairs_14: None,
            scaling: NonbondedScaling::default(),
        }
    }

    /// --- Geometry helpers ---
    #[inline] fn dist(&self, i: usize, j: usize) -> f64 {
        let rij = self.atoms[i].pos.sub(self.atoms[j].pos);
        rij.norm()
    }
    #[inline] fn angle_radians(&self, i: usize, j: usize, k: usize) -> f64 {
        let v1 = self.atoms[i].pos.sub(self.atoms[j].pos);
        let v2 = self.atoms[k].pos.sub(self.atoms[j].pos);
        let c = (v1.dot(v2)) / (v1.norm() * v2.norm());
        c.clamp(-1.0, 1.0).acos()
    }
    #[inline] fn dihedral_radians(&self, i: usize, j: usize, k: usize, l: usize) -> f64 {
        let r1 = self.atoms[j].pos.sub(self.atoms[i].pos);
        let r2 = self.atoms[k].pos.sub(self.atoms[j].pos);
        let r3 = self.atoms[l].pos.sub(self.atoms[k].pos);

        let n1 = r1.cross(r2);
        let n2 = r2.cross(r3);

        let n1u = n1.normalized();
        let n2u = n2.normalized();
        let m1 = n1u.cross(r2.normalized());

        let x = n1u.dot(n2u);
        let y = m1.dot(n2u);
        y.atan2(x) // range (-π, π]
    }

    /// --- Energies ---
    pub fn energy_bond(&self) -> f64 {
        let mut e = 0.0;
        for b in &self.bonds {
            let ti = &self.atoms[b.i].typ;
            let tj = &self.atoms[b.j].typ;
            if let Some(p) = self.params.bond_params(ti, tj, *b) {
                let r = self.dist(b.i, b.j);
                e += 0.5 * p.k_r * (r - p.r0) * (r - p.r0);
            }
        }
        e
    }

    pub fn energy_angle(&self) -> f64 {
        let mut e = 0.0;
        for a in &self.angles {
            let ti = &self.atoms[a.i].typ;
            let tj = &self.atoms[a.j].typ;
            let tk = &self.atoms[a.k].typ;
            if let Some(p) = self.params.angle_params(ti, tj, tk, *a) {
                let th = self.angle_radians(a.i, a.j, a.k);
                let th0 = p.theta0_deg.to_radians();
                e += 0.5 * p.k_th * (th - th0) * (th - th0);
            }
        }
        e
    }

    pub fn energy_torsion(&self) -> f64 {
        let mut e = 0.0;
        for d in &self.dihedrals {
            let ti = &self.atoms[d.i].typ;
            let tj = &self.atoms[d.j].typ;
            let tk = &self.atoms[d.k].typ;
            let tl = &self.atoms[d.l].typ;
            if let Some(p) = self.params.torsion_params(ti, tj, tk, tl, *d) {
                let phi = self.dihedral_radians(d.i, d.j, d.k, d.l);
                e += 0.5 * p.v * (1.0 - (p.n as f64 * phi).cos());
            }
        }
        e
    }

    pub fn energy_improper(&self) -> f64 {
        let mut e = 0.0;
        for q in &self.impropers {
            let ti = &self.atoms[q.i].typ;
            let tj = &self.atoms[q.j].typ;
            let tk = &self.atoms[q.k].typ;
            let tl = &self.atoms[q.l].typ;
            if let Some(p) = self.params.improper_params(ti, tj, tk, tl, *q) {
                let chi = self.dihedral_radians(q.i, q.j, q.k, q.l); // χ ~ 0 for planar
                e += 0.5 * p.k_chi * chi * chi;
            }
        }
        e
    }

    /// Returns (E_LJ, E_Coul) for nonbonded interactions.
    /// - Excludes 1–2 and 1–3.
    /// - Applies optional 1–4 scaling.
    /// - Uses Lorentz–Berthelot mixing (σ_ij = (σi+σj)/2, ε_ij = sqrt(εi εj)).
    pub fn energy_nonbonded(&self, cutoff: Option<f64>) -> (f64, f64) {
        let n = self.atoms.len();

        // Build adjacency for 1–2 and 1–3 exclusions
        let mut bonded: HashMap<usize, HashSet<usize>> = HashMap::new();
        for b in &self.bonds {
            bonded.entry(b.i).or_default().insert(b.j);
            bonded.entry(b.j).or_default().insert(b.i);
        }
        let mut angle13: HashSet<(usize, usize)> = HashSet::new();
        for a in &self.angles {
            let (i, k) = (a.i.min(a.k), a.i.max(a.k));
            angle13.insert((i, k));
        }

        // 1–4 inference if not provided
        let pairs_14 = if let Some(set) = &self.pairs_14 {
            set.clone()
        } else {
            let mut set = HashSet::new();
            for d in &self.dihedrals {
                let (i, l) = (d.i.min(d.l), d.i.max(d.l));
                set.insert((i, l));
            }
            set
        };

        let mut e_lj = 0.0;
        let mut e_c  = 0.0;

        const KE_COULOMB: f64 = 332.06371; // kcal·Å/(mol·e^2)

        for i in 0..n {
            let ti = &self.atoms[i].typ;
            let qi = self.atoms[i].charge;
            let vi = self.params.vdw_params(ti, i).expect("vdW params missing");
            for j in (i + 1)..n {
                // Exclusions
                if bonded.get(&i).map_or(false, |s| s.contains(&j)) { continue; } // 1–2
                if angle13.contains(&(i, j)) { continue; }                         // 1–3

                let r = self.dist(i, j);
                if let Some(rc) = cutoff { if r > rc { continue; } }

                let tj = &self.atoms[j].typ;
                let qj = self.atoms[j].charge;
                let vj = self.params.vdw_params(tj, j).expect("vdW params missing");

                // Lorentz–Berthelot mixing
                let sigma = 0.5 * (vi.sigma + vj.sigma);
                let eps = (vi.eps * vj.eps).sqrt();

                let sr = sigma / r;
                let sr6 = sr.powi(6);
                let sr12 = sr6 * sr6;
                let mut lj = 4.0 * eps * (sr12 - sr6);
                let mut coul = KE_COULOMB * qi * qj / r;

                // 1–4 scaling if applicable
                if pairs_14.contains(&(i, j)) {
                    lj   *= self.scaling.scale_14_vdw;
                    coul *= self.scaling.scale_14_coul;
                }

                e_lj += lj;
                e_c  += coul;
            }
        }
        (e_lj, e_c)
    }

    /// Convenience split + total
    pub fn energies_split(&self, cutoff: Option<f64>) -> Energies {
        let eb = self.energy_bond();
        let ea = self.energy_angle();
        let et = self.energy_torsion();
        let ei = self.energy_improper();
        let (evdw, ecoul) = self.energy_nonbonded(cutoff);
        Energies { bonds: eb, angles: ea, torsions: et, impropers: ei, vdw: evdw, coulomb: ecoul }
    }

    pub fn energy_total(&self, cutoff: Option<f64>) -> f64 {
        let s = self.energies_split(cutoff);
        s.total()
    }
}

/// Struct to return separated components.
#[derive(Clone, Copy, Debug, Default)]
pub struct Energies {
    pub bonds: f64,
    pub angles: f64,
    pub torsions: f64,
    pub impropers: f64,
    pub vdw: f64,
    pub coulomb: f64,
}
impl Energies {
    pub fn total(&self) -> f64 {
        self.bonds + self.angles + self.torsions + self.impropers + self.vdw + self.coulomb
    }
}

/* -------------------------
   Minimal example of a ParameterDB
   Replace with your real tables. This is just a stub that panics if asked.
--------------------------*/

pub struct StubParams;
impl ParameterDB for StubParams {
    fn bond_params(&self, _ti: &str, _tj: &str, _bond: Bond) -> Option<BondParams> { None }
    fn angle_params(&self, _ti: &str, _tj: &str, _tk: &str, _ang: Angle) -> Option<AngleParams> { None }
    fn torsion_params(&self, _ti: &str, _tj: &str, _tk: &str, _tl: &str, _d: Dihedral) -> Option<TorsionParams> { None }
    fn improper_params(&self, _ti: &str, _tj: &str, _tk: &str, _tl: &str, _q: Improper) -> Option<ImproperParams> { None }
    fn vdw_params(&self, _t: &str, _idx: usize) -> Option<VdwParams> {
        // Provide something reasonable if you want to test nonbonded quickly
        // Some defaults for testing (not real DREIDING values):
        Some(VdwParams { sigma: 3.5, eps: 0.1 })
    }
}

/* -------------------------
   How you might use this module:

   let atoms = vec![
       Atom { typ: "C_R".into(), charge: 0.0, pos: Vec3::new(0.0, 0.0, 0.0) },
       Atom { typ: "H".into(),   charge: 0.0, pos: Vec3::new(1.1, 0.0, 0.0) },
       // ...
   ];
   let bonds = vec![ Bond { i: 0, j: 1 } ];
   let angles = vec![];
   let dihedrals = vec![];
   let impropers = vec![];

   let mut sys = DreidingSystem::new(atoms, bonds, angles, dihedrals, impropers, MyParams { ... });
   sys.scaling = NonbondedScaling { scale_14_coul: 1.0, scale_14_vdw: 1.0 };

   let e_components = sys.energies_split(Some(12.0)); // 12 Å cutoff, or None for no cutoff
   println!("Total E = {}", e_components.total());
--------------------------*/

