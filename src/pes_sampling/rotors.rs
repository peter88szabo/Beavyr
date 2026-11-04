//! Pebble-style free internal rotations, after vibrational sampling.
//! Connectivity is determined once from equilibrium geometry, never from a
//! randomly displaced frame. Ring bonds cannot define an independent rotor.
use super::input::{Fragment, BOHR};
use std::collections::{BTreeSet, VecDeque};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FragmentSide {
    Smaller,
    First,
    Second,
    Explicit,
}
impl FragmentSide {
    pub fn label(self) -> &'static str {
        match self {
            Self::Smaller => "Smaller fragment",
            Self::First => "B side",
            Self::Second => "C side",
            Self::Explicit => "Atom list",
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ModeTreatment {
    ZpeOnly,
    Replace,
}
#[derive(Clone, Debug)]
pub struct Rotor {
    pub enabled: bool,
    pub dihedral: bool,
    /// UI atom indices are 1-based; the axis is always B-C (entries 1,2).
    pub atoms: [usize; 4],
    pub side: FragmentSide,
    pub atom_list: String,
    /// 0 means unassigned; other numbers index the positive-mode list, 1-based.
    pub mode: usize,
    pub treatment: ModeTreatment,
    pub bond_hx_angstrom: f64,
    pub bond_xx_angstrom: f64,
}
impl Default for Rotor {
    fn default() -> Self {
        Self {
            enabled: true,
            dihedral: true,
            atoms: [1, 2, 3, 4],
            side: FragmentSide::Smaller,
            atom_list: String::new(),
            mode: 0,
            treatment: ModeTreatment::ZpeOnly,
            bond_hx_angstrom: 1.4,
            bond_xx_angstrom: 2.0,
        }
    }
}
#[derive(Clone, Debug)]
pub struct Prepared {
    pub axis: [usize; 2],
    pub moving: Vec<usize>,
    pub mode: Option<usize>,
    pub treatment: ModeTreatment,
}

fn point(x: &[f64], i: usize) -> [f64; 3] {
    [x[3 * i], x[3 * i + 1], x[3 * i + 2]]
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|k| a[k] - b[k])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|k| a[k] * b[k]).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn component(graph: &[Vec<usize>], start: usize, axis: [usize; 2]) -> BTreeSet<usize> {
    let mut seen = BTreeSet::new();
    let mut queue = VecDeque::from([start]);
    while let Some(i) = queue.pop_front() {
        if !seen.insert(i) {
            continue;
        }
        for &j in &graph[i] {
            if (i == axis[0] && j == axis[1]) || (i == axis[1] && j == axis[0]) {
                continue;
            }
            if !seen.contains(&j) {
                queue.push_back(j);
            }
        }
    }
    seen
}
fn indices(text: &str, n: usize) -> Result<BTreeSet<usize>, String> {
    let number = |s: &str| {
        s.parse::<usize>()
            .ok()
            .filter(|i| *i > 0 && *i <= n)
            .ok_or_else(|| format!("Fragment atom index {s:?} is outside 1..={n}."))
    };
    let mut result = BTreeSet::new();
    for token in text
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|s| !s.is_empty())
    {
        if let Some((a, b)) = token.split_once('-') {
            let (a, b) = (number(a)?, number(b)?);
            if a > b {
                return Err("Fragment atom ranges must run from smaller to larger index.".into());
            }
            result.extend((a..=b).map(|i| i - 1));
        } else {
            result.insert(number(token)? - 1);
        }
    }
    if result.is_empty() {
        return Err("Enter the rotating fragment's atom numbers (for example 4,5,7-10).".into());
    }
    Ok(result)
}

pub fn prepare(f: &Fragment, r: &Rotor) -> Result<Prepared, String> {
    let n = f.atoms.len();
    if n < 4 {
        return Err("An internal rotor needs at least four atoms.".into());
    }
    if [r.bond_hx_angstrom, r.bond_xx_angstrom]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.)
    {
        return Err("Rotor bond-detection distances must be positive and finite.".into());
    }
    let defined = if r.dihedral {
        &r.atoms[..]
    } else {
        &r.atoms[1..3]
    };
    if defined.iter().any(|i| *i == 0 || *i > n)
        || defined.iter().collect::<BTreeSet<_>>().len() != defined.len()
    {
        return Err(format!(
            "Choose distinct rotor atom numbers between 1 and {n}."
        ));
    }
    let axis = [r.atoms[1] - 1, r.atoms[2] - 1];
    let mut graph = vec![vec![]; n];
    for i in 0..n {
        for j in 0..i {
            let ih = f.atoms[i].eq_ignore_ascii_case("H");
            let jh = f.atoms[j].eq_ignore_ascii_case("H");
            if ih && jh {
                continue;
            } // Pebble's fragment connectivity convention.
            let cutoff = if ih || jh {
                r.bond_hx_angstrom
            } else {
                r.bond_xx_angstrom
            } / BOHR;
            let v = sub(point(&f.equilibrium, i), point(&f.equilibrium, j));
            if dot(v, v) <= cutoff * cutoff {
                graph[i].push(j);
                graph[j].push(i);
            }
        }
    }
    if !graph[axis[0]].contains(&axis[1]) {
        return Err(
            "The rotor axis is not a detected bond. Check B/C and the bond-detection distances."
                .into(),
        );
    }
    let vector = sub(
        point(&f.equilibrium, axis[1]),
        point(&f.equilibrium, axis[0]),
    );
    if dot(vector, vector) < 1e-16 {
        return Err("Rotor axis atoms coincide.".into());
    }
    if r.dihedral {
        let a = r.atoms[0] - 1;
        let d = r.atoms[3] - 1;
        if !graph[a].contains(&axis[0]) || !graph[d].contains(&axis[1]) {
            return Err("The dihedral must follow bonded atoms A-B-C-D.".into());
        }
        for (outer, center) in [(a, axis[0]), (d, axis[1])] {
            let c = cross(
                vector,
                sub(point(&f.equilibrium, outer), point(&f.equilibrium, center)),
            );
            if dot(c, c) < 1e-16 {
                return Err(
                    "The dihedral is undefined because three of its atoms are collinear.".into(),
                );
            }
        }
    }
    let first = component(&graph, axis[0], axis);
    if first.contains(&axis[1]) {
        return Err("This bond lies in a ring. Rotating one side would break another bond; choose an independent rotatable bond.".into());
    }
    let second = component(&graph, axis[1], axis);
    let off_axis = |s: &BTreeSet<usize>| {
        s.iter()
            .copied()
            .filter(|i| !axis.contains(i))
            .collect::<BTreeSet<_>>()
    };
    let chosen = match r.side {
        FragmentSide::First => off_axis(&first),
        FragmentSide::Second => off_axis(&second),
        FragmentSide::Smaller => off_axis(if first.len() <= second.len() {
            &first
        } else {
            &second
        }),
        FragmentSide::Explicit => {
            let selected = off_axis(&indices(&r.atom_list, n)?);
            if selected != off_axis(&first) && selected != off_axis(&second) {
                return Err("The atom list must contain a complete fragment on one side of B-C (axis atoms are optional). A partial fragment would break its bonds.".into());
            }
            selected
        }
    };
    if chosen.is_empty()
        || !chosen.iter().any(|i| {
            let c = cross(
                vector,
                sub(point(&f.equilibrium, *i), point(&f.equilibrium, axis[0])),
            );
            dot(c, c) > 1e-16
        })
    {
        return Err("The selected fragment has no atoms away from the rotation axis.".into());
    }
    if r.mode > f.frequencies.len() {
        return Err("The assigned rotor mode does not exist in this Hessian.".into());
    }
    Ok(Prepared {
        axis,
        moving: chosen.into_iter().collect(),
        mode: r.mode.checked_sub(1),
        treatment: r.treatment,
    })
}

pub fn prepare_all(f: &Fragment, rotors: &[Rotor]) -> Result<Vec<Prepared>, String> {
    let mut bonds = BTreeSet::new();
    let mut modes = BTreeSet::new();
    let mut prepared = vec![];
    for (i, r) in rotors.iter().enumerate().filter(|(_, r)| r.enabled) {
        let p = prepare(f, r).map_err(|e| format!("Rotor {}: {e}", i + 1))?;
        let mut axis = p.axis;
        axis.sort();
        if !bonds.insert(axis) {
            return Err("The same bond was defined as a rotor more than once.".into());
        }
        if p.mode.is_some_and(|m| !modes.insert(m)) {
            return Err("A normal mode cannot be assigned to multiple rotors.".into());
        }
        prepared.push(p);
    }
    Ok(prepared)
}

pub fn rotate(x: &mut [f64], rotor: &Prepared, angle: f64) -> Result<(), String> {
    let origin = point(x, rotor.axis[0]);
    let vector = sub(point(x, rotor.axis[1]), origin);
    let norm = dot(vector, vector).sqrt();
    if !norm.is_finite() || norm < 1e-10 {
        return Err("Sampled rotor axis has zero or invalid length. Reduce the vibrational/TS displacement.".into());
    }
    let axis = vector.map(|v| v / norm);
    let (sin, cos) = angle.sin_cos();
    for &i in &rotor.moving {
        let v = sub(point(x, i), origin);
        let c = cross(axis, v);
        let parallel = dot(axis, v);
        for k in 0..3 {
            x[3 * i + k] = origin[k] + v[k] * cos + c[k] * sin + axis[k] * parallel * (1. - cos);
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "rotor_tests.rs"]
mod tests;
