//! Internal coordinate machinery (bonds/angles/dihedrals).
//!
//! - `StrucAnalyze` (bond/angle/dihedral discovery)
//! - `ComputeInternals!`
//! - `BPGMatrix` (Wilson B matrix and pseudo-inverse metric `invG`)
//! - `DiffInternals` (dihedral wrap handling)
//! - `BestFitofdqtoCart!` (iterative best-fit mapping dq -> dx)

use anyhow::{bail, Result};

use crate::normalmode::linalg_shim::{Backend, LinAlg};
use ndarray::{Array1, Array2};

#[derive(Debug, Clone)]
pub struct InternalCoords {
    pub nat: usize,
    pub bonds: Vec<[usize; 2]>,
    pub angles: Vec<[usize; 3]>,
    pub dihedrals: Vec<[usize; 4]>,
}

impl InternalCoords {
    pub fn nint(&self) -> usize {
        self.bonds.len() + self.angles.len() + self.dihedrals.len()
    }

    pub fn nbonds(&self) -> usize {
        self.bonds.len()
    }

    pub fn nangles(&self) -> usize {
        self.angles.len()
    }

    pub fn ndiheds(&self) -> usize {
        self.dihedrals.len()
    }
}

#[derive(Debug, Clone)]
pub struct ConnectivityModel {
    pub rcov_disp: Vec<f64>, // bohr
    pub kcn: f64,
    pub facmin: f64,
}

#[derive(Debug, Clone)]
pub struct Bpg {
    pub b: Vec<Vec<f64>>,     // (nint x nx)
    pub inv_g: Vec<Vec<f64>>, // (nint x nint)
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm2(a: [f64; 3]) -> f64 {
    dot(a, a)
}

fn norm(a: [f64; 3]) -> f64 {
    norm2(a).sqrt()
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn get_atom(x: &[f64], i: usize) -> [f64; 3] {
    let k = 3 * i;
    [x[k], x[k + 1], x[k + 2]]
}

fn set_atom(x: &mut [f64], i: usize, v: [f64; 3]) {
    let k = 3 * i;
    x[k] = v[0];
    x[k + 1] = v[1];
    x[k + 2] = v[2];
}

fn bond_angle(ba: [f64; 3], bc: [f64; 3]) -> f64 {
    let r1 = norm(ba);
    let r2 = norm(bc);
    if r1 == 0.0 || r2 == 0.0 {
        return 0.0;
    }
    let mut c = dot(ba, bc) / (r1 * r2);
    if c > 1.0 {
        c = 1.0;
    } else if c < -1.0 {
        c = -1.0;
    }
    c.acos()
}

fn dihedral_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let ab = sub(b, a);
    let bc = sub(c, b);
    let cd = sub(d, c);
    let t = cross(ab, bc);
    let u = cross(bc, cd);
    let tu = cross(t, u);
    let rt2 = norm2(t);
    let ru2 = norm2(u);
    let vtvu = (rt2 * ru2).sqrt();
    let rbc = norm(bc);
    if vtvu == 0.0 || rbc == 0.0 {
        return 0.0;
    }
    let cosang = dot(t, u) / vtvu;
    let sinang = dot(bc, tu) / (rbc * vtvu);
    sinang.atan2(cosang)
}

pub fn compute_internals(x: &[f64], ic: &InternalCoords) -> Vec<f64> {
    let nb = ic.nbonds();
    let na = ic.nangles();
    let nd = ic.ndiheds();
    let nint = nb + na + nd;
    let mut q = vec![0.0f64; nint];

    for (i, [a, b]) in ic.bonds.iter().enumerate() {
        let ra = get_atom(x, *a);
        let rb = get_atom(x, *b);
        q[i] = norm(sub(rb, ra));
    }

    for (i, [a, b, c]) in ic.angles.iter().enumerate() {
        let ra = get_atom(x, *a);
        let rb = get_atom(x, *b);
        let rc = get_atom(x, *c);
        let ba = sub(ra, rb);
        let bc = sub(rc, rb);
        q[nb + i] = bond_angle(ba, bc);
    }

    for (i, [a, b, c, d]) in ic.dihedrals.iter().enumerate() {
        let ra = get_atom(x, *a);
        let rb = get_atom(x, *b);
        let rc = get_atom(x, *c);
        let rd = get_atom(x, *d);
        q[nb + na + i] = dihedral_angle(ra, rb, rc, rd);
    }

    q
}

pub fn diff_internals(q2: &[f64], q1: &[f64], ndiheds: usize) -> Vec<f64> {
    let nint = q2.len();
    let mut dq = vec![0.0f64; nint];
    for i in 0..nint {
        dq[i] = q2[i] - q1[i];
    }
    // Wrap dihedrals to (-pi, pi]
    for i in 0..ndiheds {
        let k = nint - ndiheds + i;
        let mut v = q2[k] - q1[k];
        if v > std::f64::consts::PI {
            v -= 2.0 * std::f64::consts::PI;
        } else if v < -std::f64::consts::PI {
            v += 2.0 * std::f64::consts::PI;
        }
        dq[k] = v;
    }
    dq
}

fn dist_idx(nat: usize, i: usize, j: usize) -> usize {
    i * nat + j
}

fn distance_table(x: &[f64], nat: usize) -> Vec<f64> {
    let mut rab = vec![0.0f64; nat * nat];
    for i in 0..nat {
        for j in (i + 1)..nat {
            let ri = get_atom(x, i);
            let rj = get_atom(x, j);
            let r = norm(sub(rj, ri));
            rab[dist_idx(nat, i, j)] = r;
            rab[dist_idx(nat, j, i)] = r;
        }
    }
    rab
}

pub fn analyze_structure(x: &[f64], model: &ConnectivityModel) -> Result<InternalCoords> {
    let nat = model.rcov_disp.len();
    if x.len() != 3 * nat {
        bail!(
            "analyze_structure: expected x length {}, got {}",
            3 * nat,
            x.len()
        );
    }
    let rab = distance_table(x, nat);

    let mut con = vec![0u8; nat * nat];
    let mut bonds: Vec<[usize; 2]> = Vec::new();
    for i in 0..(nat - 1) {
        for j in (i + 1)..nat {
            let rco = model.rcov_disp[i] + model.rcov_disp[j];
            let rij = rab[dist_idx(nat, i, j)];
            if rij == 0.0 {
                continue;
            }
            let fac = 1.0 / (1.0 + (-model.kcn * (rco / rij - 1.0)).exp());
            if fac > model.facmin {
                con[dist_idx(nat, i, j)] = 1;
                con[dist_idx(nat, j, i)] = 1;
                bonds.push([i, j]);
            }
        }
    }

    let mut angles: Vec<[usize; 3]> = Vec::new();
    for i in 0..nat {
        let mut bats: Vec<usize> = Vec::new();
        for j in 0..nat {
            if con[dist_idx(nat, i, j)] != 0 {
                bats.push(j);
            }
        }
        match bats.len() {
            0 | 1 => {}
            _ => {
                // 2..4 neighbors; we generalize to all unique pairs.
                for a in 0..bats.len() {
                    for b in (a + 1)..bats.len() {
                        angles.push([bats[a], i, bats[b]]);
                    }
                }
            }
        }
    }

    let mut dihedrals: Vec<[usize; 4]> = Vec::new();
    for i in 0..(nat - 1) {
        // neighbors of i
        let mut neigh_i: Vec<usize> = Vec::new();
        for k in 0..nat {
            if con[dist_idx(nat, i, k)] != 0 {
                neigh_i.push(k);
            }
        }
        if neigh_i.len() < 2 {
            continue;
        }

        for j in (i + 1)..nat {
            if con[dist_idx(nat, i, j)] == 0 {
                continue;
            }

            let mut neigh_j: Vec<usize> = Vec::new();
            for k in 0..nat {
                if con[dist_idx(nat, j, k)] != 0 {
                    neigh_j.push(k);
                }
            }
            if neigh_j.len() < 2 {
                continue;
            }

            for &ii in &neigh_i {
                if ii == j {
                    continue;
                }
                for &jj in &neigh_j {
                    if jj == i {
                        continue;
                    }
                    dihedrals.push([ii, i, j, jj]);
                }
            }
        }
    }

    Ok(InternalCoords {
        nat,
        bonds,
        angles,
        dihedrals,
    })
}

pub fn internal_coords_from_bonds(nat: usize, bonds_in: &[[usize; 2]]) -> Result<InternalCoords> {
    let mut con = vec![0u8; nat * nat];
    let mut bonds: Vec<[usize; 2]> = Vec::new();
    for &[i, j] in bonds_in {
        if i >= nat || j >= nat {
            bail!("internal_coords_from_bonds: bond index out of range ({i}, {j})");
        }
        if i == j {
            bail!("internal_coords_from_bonds: self-bond at {i}");
        }
        if con[dist_idx(nat, i, j)] == 0 {
            bonds.push([i, j]);
            con[dist_idx(nat, i, j)] = 1;
            con[dist_idx(nat, j, i)] = 1;
        }
    }

    let mut angles: Vec<[usize; 3]> = Vec::new();
    for i in 0..nat {
        let mut bats: Vec<usize> = Vec::new();
        for j in 0..nat {
            if con[dist_idx(nat, i, j)] != 0 {
                bats.push(j);
            }
        }
        if bats.len() >= 2 {
            for a in 0..bats.len() {
                for b in (a + 1)..bats.len() {
                    angles.push([bats[a], i, bats[b]]);
                }
            }
        }
    }

    let mut dihedrals: Vec<[usize; 4]> = Vec::new();
    for i in 0..(nat - 1) {
        let mut neigh_i: Vec<usize> = Vec::new();
        for k in 0..nat {
            if con[dist_idx(nat, i, k)] != 0 {
                neigh_i.push(k);
            }
        }
        if neigh_i.len() < 2 {
            continue;
        }

        for j in (i + 1)..nat {
            if con[dist_idx(nat, i, j)] == 0 {
                continue;
            }

            let mut neigh_j: Vec<usize> = Vec::new();
            for k in 0..nat {
                if con[dist_idx(nat, j, k)] != 0 {
                    neigh_j.push(k);
                }
            }
            if neigh_j.len() < 2 {
                continue;
            }

            for &ii in &neigh_i {
                if ii == j {
                    continue;
                }
                for &jj in &neigh_j {
                    if jj == i {
                        continue;
                    }
                    dihedrals.push([ii, i, j, jj]);
                }
            }
        }
    }

    Ok(InternalCoords {
        nat,
        bonds,
        angles,
        dihedrals,
    })
}

pub fn bpg_matrix(x: &[f64], ic: &InternalCoords) -> Result<Bpg> {
    let nat = ic.nat;
    let nx = 3 * nat;
    if x.len() != nx {
        bail!("bpg_matrix: expected x length {}, got {}", nx, x.len());
    }
    let nb = ic.nbonds();
    let na = ic.nangles();
    let nd = ic.ndiheds();
    let nint = nb + na + nd;

    let mut b = vec![vec![0.0f64; nx]; nint];

    // Bonds
    for (i, [a, c]) in ic.bonds.iter().enumerate() {
        let ra = get_atom(x, *a);
        let rc = get_atom(x, *c);
        let ab = sub(rc, ra);
        let dist = norm(ab);
        if dist == 0.0 {
            continue;
        }
        let g = scale(ab, 1.0 / dist);
        for k in 0..3 {
            b[i][3 * a + k] = -g[k];
            b[i][3 * c + k] = g[k];
        }
    }

    // Angles
    for (ii, [a, b0, c]) in ic.angles.iter().enumerate() {
        let j = nb + ii;
        let ra = get_atom(x, *a);
        let rb = get_atom(x, *b0);
        let rc = get_atom(x, *c);
        let ba = sub(ra, rb);
        let bc = sub(rc, rb);
        let rba2 = norm2(ba);
        let rbc2 = norm2(bc);
        let p = cross(bc, ba);
        let rp = norm(p);
        if rp == 0.0 || rba2 == 0.0 || rbc2 == 0.0 {
            continue;
        }
        let inc_a = scale(cross(ba, p), -1.0 / (rba2 * rp));
        let inc_c = scale(cross(bc, p), 1.0 / (rbc2 * rp));
        for k in 0..3 {
            b[j][3 * a + k] = inc_a[k];
            b[j][3 * c + k] = inc_c[k];
            b[j][3 * b0 + k] = -inc_a[k] - inc_c[k];
        }
    }

    // Dihedrals
    for (ii, [a, b0, c, d]) in ic.dihedrals.iter().enumerate() {
        let j = nb + na + ii;
        let ra = get_atom(x, *a);
        let rb = get_atom(x, *b0);
        let rc = get_atom(x, *c);
        let rd = get_atom(x, *d);

        let ab = sub(rb, ra);
        let ac = sub(rc, ra);
        let bc = sub(rc, rb);
        let bd = sub(rd, rb);
        let cd = sub(rd, rc);

        let t = cross(ab, bc);
        let u = cross(bc, cd);
        let rt2 = norm2(t);
        let ru2 = norm2(u);
        let rbc = norm(bc);
        if rt2 == 0.0 || ru2 == 0.0 || rbc == 0.0 {
            continue;
        }

        let inc_t = scale(cross(t, bc), 1.0 / (rt2 * rbc));
        let inc_u = scale(cross(u, bc), -1.0 / (ru2 * rbc));

        let v1 = cross(inc_t, bc);
        let v2 = cross(ac, inc_t);
        let v3 = cross(inc_u, cd);
        let v4 = cross(inc_t, ab);
        let v5 = cross(bd, inc_u);
        let v6 = cross(inc_u, bc);

        let vb = add(v2, v3);
        let vc = add(v4, v5);
        for k in 0..3 {
            b[j][3 * a + k] = v1[k];
            b[j][3 * b0 + k] = vb[k];
            b[j][3 * c + k] = vc[k];
            b[j][3 * d + k] = v6[k];
        }
    }

    // G = B * B^T
    let mut g = vec![vec![0.0f64; nint]; nint];
    for i in 0..nint {
        for j in i..nint {
            let mut acc = 0.0;
            for k in 0..nx {
                acc += b[i][k] * b[j][k];
            }
            g[i][j] = acc;
            g[j][i] = acc;
        }
    }

    // Pseudoinverse of G via symmetric eigendecomposition.
    let la = LinAlg::new(Backend::Auto)?;
    let mut g_a = Array2::<f64>::zeros((nint, nint));
    for i in 0..nint {
        for j in 0..nint {
            g_a[(i, j)] = g[i][j];
        }
    }
    let (evals, evecs) = la.eigh(&g_a)?;
    let mut inv_g = vec![vec![0.0f64; nint]; nint];
    let cutoff = 1e-14_f64;
    for (k, &lam) in evals.iter().enumerate() {
        if lam <= cutoff {
            continue;
        }
        let inv = 1.0 / lam;
        // Add inv * v_k v_k^T (columns in evecs).
        for i in 0..nint {
            let vi = evecs[(i, k)];
            for j in 0..nint {
                inv_g[i][j] += inv * vi * evecs[(j, k)];
            }
        }
    }

    Ok(Bpg { b, inv_g })
}

pub fn dq_to_dx(bpg: &Bpg, dq: &[f64]) -> Result<Vec<f64>> {
    let nint = bpg.b.len();
    if nint == 0 {
        return Ok(Vec::new());
    }
    if dq.len() != nint {
        bail!(
            "dq_to_dx: dq length mismatch: dq={} nint={}",
            dq.len(),
            nint
        );
    }
    let nx = bpg.b[0].len();

    // tmp = invG * dq
    let mut tmp = vec![0.0f64; nint];
    for i in 0..nint {
        let mut acc = 0.0;
        for j in 0..nint {
            acc += bpg.inv_g[i][j] * dq[j];
        }
        tmp[i] = acc;
    }

    // dx = B^T * tmp
    let mut dx = vec![0.0f64; nx];
    for k in 0..nx {
        let mut acc = 0.0;
        for i in 0..nint {
            acc += bpg.b[i][k] * tmp[i];
        }
        dx[k] = acc;
    }
    Ok(dx)
}

pub fn best_fit_dq_to_cart(
    x: &mut [f64],
    qs: &[f64],
    dq: &[f64],
    ic: &InternalCoords,
    bpg: &Bpg,
    n_iter: usize,
) -> Result<()> {
    let nint = ic.nint();
    if qs.len() != nint || dq.len() != nint {
        bail!(
            "best_fit_dq_to_cart: length mismatch: qs={} dq={} nint={}",
            qs.len(),
            dq.len(),
            nint
        );
    }

    let mut tmpx = x.to_vec();
    let mut tmpqs = vec![0.0f64; nint];

    // Initial dx.
    let dx0 = dq_to_dx(bpg, dq)?;
    for i in 0..tmpx.len() {
        tmpx[i] += dx0[i];
    }

    // Iterative “best fit” correction.
    for _ in 0..n_iter {
        tmpqs.clone_from_slice(&compute_internals(&tmpx, ic));
        let dqq_init = diff_internals(&tmpqs, qs, ic.ndiheds());
        let ddq = diff_internals(dq, &dqq_init, ic.ndiheds());
        let dx = dq_to_dx(bpg, &ddq)?;
        for i in 0..tmpx.len() {
            tmpx[i] += dx[i];
        }
    }

    x.copy_from_slice(&tmpx);
    Ok(())
}

/// Helper to flatten `[f64;3]` coordinates into `x` (bohr) for debugging / interop.
pub fn pack_xyz(coords: &[[f64; 3]]) -> Vec<f64> {
    let mut x = Vec::with_capacity(3 * coords.len());
    for c in coords {
        x.push(c[0]);
        x.push(c[1]);
        x.push(c[2]);
    }
    x
}

/// Helper to unpack `x` (bohr) into `[f64;3]` coordinates.
pub fn unpack_xyz(x: &[f64]) -> Vec<[f64; 3]> {
    let nat = x.len() / 3;
    let mut out = Vec::with_capacity(nat);
    for i in 0..nat {
        out.push(get_atom(x, i));
    }
    out
}
