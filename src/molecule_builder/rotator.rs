// src/molecule_builder/rotator.rs

use std::collections::{HashMap, VecDeque};
use bevy::prelude::*;

/// Which side of the bond should rotate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RotateSide { A, B }

#[inline]
fn dist(a: Vec3, b: Vec3) -> f32 { a.distance(b) }

/// Build adjacency list with simple distance thresholds (in Å).
/// H–H ignored; H–X uses `hx`; X–X uses `xx`.
fn adjacency_list(
    atoms: &[String],
    coords: &[Vec3],
    hx: f32,
    xx: f32,
) -> HashMap<usize, Vec<usize>> {
    let n = atoms.len();
    let mut adj: HashMap<usize, Vec<usize>> = (0..n).map(|i| (i, Vec::new())).collect();

    for i in 0..n {
        for j in (i + 1)..n {
            let hi = atoms[i] == "H";
            let hj = atoms[j] == "H";
            if hi && hj { continue; }
            let thr = if !hi && !hj { xx } else { hx };
            if dist(coords[i], coords[j]) <= thr {
                adj.get_mut(&i).unwrap().push(j);
                adj.get_mut(&j).unwrap().push(i);
            }
        }
    }
    adj
}

/// BFS traversal for one cluster, excluding a specific node.
fn bfs_cluster(
    adj: &HashMap<usize, Vec<usize>>,
    start: usize,
    visited: &mut [bool],
    exclude: usize,
) -> Vec<usize> {
    let mut cluster = Vec::new();
    let mut q = VecDeque::from([start]);
    while let Some(v) = q.pop_front() {
        if visited[v] || v == exclude { continue; }
        visited[v] = true;
        cluster.push(v);
        if let Some(nei) = adj.get(&v) {
            for &u in nei {
                if !visited[u] && u != exclude { q.push_back(u); }
            }
        }
    }
    cluster
}

/// Return the two components you get by cutting bond (a,b).
pub fn split_sides(
    atoms: &[String],
    coords: &[Vec3],
    a: usize,
    b: usize,
    bond_th_hx: f32,
    bond_th_xx: f32,
) -> (Vec<usize>, Vec<usize>) {
    assert!(a < atoms.len() && b < atoms.len());

    let adj = adjacency_list(atoms, coords, bond_th_hx, bond_th_xx);
    let mut visited = vec![false; atoms.len()];
    let side_a = bfs_cluster(&adj, a, &mut visited, b);
    let side_b = bfs_cluster(&adj, b, &mut visited, a);
    (side_a, side_b)
}

/// Rotate the atoms on the chosen side of bond (A,B) by `angle_deg` (degrees).
/// Uses Rodrigues' rotation formula on Vec3.
pub fn rotate_side(
    atoms: &[String],
    coords: &[Vec3],
    a: usize,
    b: usize,
    which: RotateSide,
    angle_deg: f32,
    bond_th_hx: f32,
    bond_th_xx: f32,
) -> Vec<Vec3> {
    assert!(a < atoms.len() && b < atoms.len());

    let (side_a, side_b) = split_sides(atoms, coords, a, b, bond_th_hx, bond_th_xx);
    let (rot_center_idx, rotating_set) = match which {
        RotateSide::A => (a, side_a),
        RotateSide::B => (b, side_b),
    };

    let mut new_coords = coords.to_vec();
    let rot_center = coords[rot_center_idx];
    let axis = (coords[b] - coords[a]).normalize();
    let angle = angle_deg.to_radians();

    // Rodrigues' rotation with Vec3 ops
    for &i in &rotating_set {
        if i == rot_center_idx { continue; }
        let v = coords[i] - rot_center;
        let v_rot = v * angle.cos()
            + axis.cross(v) * angle.sin()
            + axis * (axis.dot(v)) * (1.0 - angle.cos());
        new_coords[i] = rot_center + v_rot;
    }

    new_coords
}

/// Translate the atoms on the chosen side of bond (A,B) along the bond axis by `distance` (Å).
pub fn translate_side(
    atoms: &[String],
    coords: &[Vec3],
    a: usize,
    b: usize,
    which: RotateSide,
    distance: f32,
    bond_th_hx: f32,
    bond_th_xx: f32,
) -> Vec<Vec3> {
    assert!(a < atoms.len() && b < atoms.len());

    let axis_vec = coords[b] - coords[a];
    let axis_len = axis_vec.length();
    if axis_len <= 1e-6 || distance.abs() <= 1e-6 {
        return coords.to_vec();
    }

    let (side_a, side_b) = split_sides(atoms, coords, a, b, bond_th_hx, bond_th_xx);
    let (dir, moving_set) = match which {
        RotateSide::A => (-axis_vec / axis_len, side_a),
        RotateSide::B => (axis_vec / axis_len, side_b),
    };

    let mut new_coords = coords.to_vec();
    let delta = dir * distance;
    for &i in &moving_set {
        new_coords[i] = coords[i] + delta;
    }

    new_coords
}

/// Bend the atoms on the chosen side so the picked atom makes `target_angle_deg`
/// with the bond axis (A→B), using the side's endpoint as the angle vertex.
pub fn bend_side(
    atoms: &[String],
    coords: &[Vec3],
    a: usize,
    b: usize,
    which: RotateSide,
    picked_idx: usize,
    target_angle_deg: f32,
    bond_th_hx: f32,
    bond_th_xx: f32,
) -> Vec<Vec3> {
    assert!(a < atoms.len() && b < atoms.len());
    if picked_idx >= atoms.len() {
        return coords.to_vec();
    }

    let axis_vec = coords[b] - coords[a];
    let axis_len = axis_vec.length();
    if axis_len <= 1e-6 {
        return coords.to_vec();
    }
    let axis = axis_vec / axis_len;

    let (side_a, side_b) = split_sides(atoms, coords, a, b, bond_th_hx, bond_th_xx);
    let (center_idx, rotating_set) = match which {
        RotateSide::A => (a, side_a),
        RotateSide::B => (b, side_b),
    };
    if !rotating_set.contains(&picked_idx) {
        // If the picked atom isn't on the rotating side, don't apply.
        return coords.to_vec();
    }

    let center = coords[center_idx];
    let v = coords[picked_idx] - center;
    let v_len = v.length();
    if v_len <= 1e-6 {
        return coords.to_vec();
    }
    let v_n = v / v_len;

    let mut dot = axis.dot(v_n);
    dot = dot.clamp(-1.0, 1.0);
    let current_angle = dot.acos();
    let target_angle = target_angle_deg.to_radians();
    let delta = target_angle - current_angle;
    if delta.abs() <= 1e-6 {
        return coords.to_vec();
    }

    let mut rot_axis = axis.cross(v_n);
    if rot_axis.length() <= 1e-6 {
        // Pick a stable perpendicular axis when v is (anti)parallel to axis.
        rot_axis = axis.cross(Vec3::X);
        if rot_axis.length() <= 1e-6 {
            rot_axis = axis.cross(Vec3::Y);
        }
        if rot_axis.length() <= 1e-6 {
            return coords.to_vec();
        }
    }
    let rot_axis = rot_axis.normalize();

    let mut new_coords = coords.to_vec();
    for &i in &rotating_set {
        if i == center_idx { continue; }
        let p = coords[i] - center;
        let p_rot = p * delta.cos()
            + rot_axis.cross(p) * delta.sin()
            + rot_axis * (rot_axis.dot(p)) * (1.0 - delta.cos());
        new_coords[i] = center + p_rot;
    }

    new_coords
}
