//! Element-preserving atom correspondence and rigid-alignment displacements.
use crate::molecule::Molecule;
use bevy::prelude::*;
use std::{
    collections::BTreeMap,
    hash::{Hash, Hasher},
};
pub fn key(m: &Molecule) -> u64 {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    m.atoms.hash(&mut h);
    for p in &m.pos {
        for v in p.to_array() {
            v.to_bits().hash(&mut h);
        }
    }
    h.finish()
}
pub fn validate(a: &Molecule, b: &Molecule, map: &[usize]) -> Result<(), String> {
    if a.atoms.is_empty() || a.atoms.len() != b.atoms.len() || map.len() != a.atoms.len() {
        return Err("Comparison requires equal, nonzero atom counts.".into());
    }
    let mut seen = vec![false; b.atoms.len()];
    for (i, &j) in map.iter().enumerate() {
        if j >= seen.len() || seen[j] || a.atoms[i] != b.atoms[j] {
            return Err(
                "Mapping must pair each A atom with a unique B atom of the same element.".into(),
            );
        }
        seen[j] = true;
    }
    if a.pos.len() != a.atoms.len()
        || b.pos.len() != b.atoms.len()
        || a.pos.iter().chain(&b.pos).any(|p| !p.is_finite())
    {
        return Err("Invalid comparison coordinates.".into());
    }
    Ok(())
}
pub fn parse(a: &Molecule, b: &Molecule, text: &str) -> Result<Vec<usize>, String> {
    let map = text
        .split(|c: char| c.is_whitespace() || c == ',')
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.parse::<usize>()
                .ok()
                .and_then(|i| i.checked_sub(1))
                .ok_or("Use 1-based B atom indices.".to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    validate(a, b, &map)?;
    Ok(map)
}
pub struct Aligned {
    pub positions: Vec<Vec3>,
    pub distances: Vec<f32>,
    pub rmsd: f64,
}
pub fn aligned(a: &Molecule, b: &Molecule, map: &[usize]) -> Result<Aligned, String> {
    validate(a, b, map)?;
    let ordered: Vec<_> = map.iter().map(|&j| b.pos[j]).collect();
    let t = crate::molecular_alignment::superpose(&ordered, &a.pos)?;
    let positions = t.apply_all(&b.pos);
    let distances: Vec<_> = map
        .iter()
        .enumerate()
        .map(|(i, &j)| positions[j].distance(a.pos[i]))
        .collect();
    if distances.iter().any(|d| !d.is_finite()) {
        return Err("Alignment overflow: check coordinate magnitudes.".into());
    }
    let rmsd = (distances.iter().map(|d| (*d as f64).powi(2)).sum::<f64>()
        / distances.len() as f64)
        .sqrt();
    Ok(Aligned {
        positions,
        distances,
        rmsd,
    })
}
// Hungarian assignment, square finite costs. O(n^3), with a size limit at the caller.
fn assignment(cost: &[Vec<f64>]) -> Vec<usize> {
    let n = cost.len();
    let (mut u, mut v, mut p, mut way) = (
        vec![0.; n + 1],
        vec![0.; n + 1],
        vec![0usize; n + 1],
        vec![0usize; n + 1],
    );
    for i in 1..=n {
        p[0] = i;
        let mut j0 = 0;
        let mut min = vec![f64::INFINITY; n + 1];
        let mut used = vec![false; n + 1];
        loop {
            used[j0] = true;
            let i0 = p[j0];
            let mut delta = f64::INFINITY;
            let mut j1 = 0;
            for j in 1..=n {
                if !used[j] {
                    let c = cost[i0 - 1][j - 1] - u[i0] - v[j];
                    if c < min[j] {
                        min[j] = c;
                        way[j] = j0;
                    }
                    if min[j] < delta {
                        delta = min[j];
                        j1 = j;
                    }
                }
            }
            for j in 0..=n {
                if used[j] {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    min[j] -= delta;
                }
            }
            j0 = j1;
            if p[j0] == 0 {
                break;
            }
        }
        loop {
            let j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
            if j0 == 0 {
                break;
            }
        }
    }
    let mut result = vec![0; n];
    for j in 1..=n {
        result[p[j] - 1] = j - 1;
    }
    result
}
fn groups(m: &Molecule) -> BTreeMap<String, Vec<usize>> {
    let mut g: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, a) in m.atoms.iter().enumerate() {
        g.entry(a.clone()).or_default().push(i);
    }
    g
}
fn descriptors(m: &Molecule, groups: &BTreeMap<String, Vec<usize>>) -> Vec<Vec<f64>> {
    m.pos
        .iter()
        .map(|p| {
            groups
                .values()
                .flat_map(|indices| {
                    let mut d: Vec<_> = indices
                        .iter()
                        .map(|&j| p.as_dvec3().distance(m.pos[j].as_dvec3()))
                        .collect();
                    d.sort_by(f64::total_cmp);
                    d
                })
                .collect()
        })
        .collect()
}
pub fn automatic(a: &Molecule, b: &Molecule) -> Result<Vec<usize>, String> {
    let ga = groups(a);
    let gb = groups(b);
    if a.atoms.is_empty()
        || ga
            .iter()
            .any(|(s, g)| gb.get(s).map(Vec::len) != Some(g.len()))
        || ga.len() != gb.len()
    {
        return Err("The structures must have the same elemental composition.".into());
    }
    if a.atoms.len() > 1000 || ga.values().any(|g| g.len() > 500) {
        return Err("Automatic mapping is limited to 1,000 atoms (500 of one element). Use atom order or enter the correspondence for larger structures.".into());
    }
    if a.pos.len() != a.atoms.len()
        || b.pos.len() != b.atoms.len()
        || a.pos.iter().chain(&b.pos).any(|p| !p.is_finite())
    {
        return Err("Invalid coordinates.".into());
    }
    let da = descriptors(a, &ga);
    let db = descriptors(b, &gb);
    let mut map = vec![0; a.atoms.len()];
    for (s, ai) in &ga {
        let bj = &gb[s];
        let cost: Vec<Vec<_>> = ai
            .iter()
            .map(|&i| {
                bj.iter()
                    .map(|&j| da[i].iter().zip(&db[j]).map(|(x, y)| (x - y).powi(2)).sum())
                    .collect()
            })
            .collect();
        let picked = assignment(&cost);
        for (k, &i) in ai.iter().enumerate() {
            map[i] = bj[picked[k]];
        }
    }
    let mut best = map.clone();
    let mut score = aligned(a, b, &map)?.rmsd;
    for _ in 0..12 {
        let fit = aligned(a, b, &map)?;
        let mut next = map.clone();
        for (s, ai) in &ga {
            let bj = &gb[s];
            let cost: Vec<Vec<_>> = ai
                .iter()
                .map(|&i| {
                    bj.iter()
                        .map(|&j| {
                            fit.positions[j]
                                .as_dvec3()
                                .distance_squared(a.pos[i].as_dvec3())
                        })
                        .collect()
                })
                .collect();
            let picked = assignment(&cost);
            for (k, &i) in ai.iter().enumerate() {
                next[i] = bj[picked[k]];
            }
        }
        if next == map {
            break;
        }
        let value = aligned(a, b, &next)?.rmsd;
        if value < score {
            best = next.clone();
            score = value;
        }
        map = next;
    }
    // Keep known atom order if it gives a better geometric fit.
    let identity: Vec<_> = (0..a.atoms.len()).collect();
    if let Ok(fit) = aligned(a, b, &identity) {
        if fit.rmsd < score {
            best = identity;
        }
    }
    Ok(best)
}
pub fn color(distance: f32, maximum: f32) -> Color {
    let t = (distance / maximum.max(1e-6)).clamp(0., 1.);
    Color::srgb(t, 0.25 + 0.55 * (1. - (2. * t - 1.).abs()), 1. - t)
}
pub fn linked_count(a: usize, b: usize) -> usize {
    match (a > 1, b > 1) {
        (true, true) => a.min(b),
        (true, false) => a,
        (false, true) => b,
        _ => 1,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn global_assignment_is_not_greedy() {
        assert_eq!(assignment(&[vec![1., 2.], vec![1.1, 100.]]), [1, 0]);
    }
    #[test]
    fn mapping_handles_reordering_rotation_and_translation() {
        let a = Molecule::from_xyz("4\nA\nC 0 0 0\nO 1 0 0\nH 0 2 0\nH 0 0 3\n");
        let mut b = a.clone();
        let order = [3, 1, 0, 2];
        b.atoms = order.iter().map(|&i| a.atoms[i].clone()).collect();
        b.pos = order
            .iter()
            .map(|&i| Quat::from_rotation_z(1.2) * a.pos[i] + Vec3::splat(5.))
            .collect();
        let map = automatic(&a, &b).unwrap();
        assert_eq!(map, [2, 1, 3, 0]);
        assert!(aligned(&a, &b, &map).unwrap().rmsd < 1e-5);
        assert!(parse(&a, &b, "3 2 4 4").is_err());
    }
    #[test]
    fn linked_frames_never_run_past_shorter_trajectory() {
        assert_eq!(linked_count(10, 7), 7);
        assert_eq!(linked_count(1, 10), 10);
    }
}
