// src/molecule_builder/rotator.rs

use bevy::prelude::*;
use std::collections::{HashMap, VecDeque};

/// Which side of the bond should rotate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RotateSide {
    A,
    B,
}

/// Return the two disconnected components produced by removing the selected
/// covalent bond. `None` means the endpoints are not a bond or the bond is in
/// a ring, so there is no independently rotatable fragment.
pub fn split_sides_from_bonds(
    atom_count: usize,
    bonds: &[(usize, usize, f32)],
    a: usize,
    b: usize,
) -> Option<(Vec<usize>, Vec<usize>)> {
    if a >= atom_count || b >= atom_count || a == b {
        return None;
    }

    let mut adjacency = vec![Vec::new(); atom_count];
    let mut selected_bond_found = false;
    for &(i, j, _) in bonds {
        if i >= atom_count || j >= atom_count {
            continue;
        }
        if (i == a && j == b) || (i == b && j == a) {
            selected_bond_found = true;
            continue;
        }
        adjacency[i].push(j);
        adjacency[j].push(i);
    }
    if !selected_bond_found {
        return None;
    }

    let component = |start: usize| {
        let mut visited = vec![false; atom_count];
        let mut queue = VecDeque::from([start]);
        let mut result = Vec::new();
        while let Some(index) = queue.pop_front() {
            if visited[index] {
                continue;
            }
            visited[index] = true;
            result.push(index);
            queue.extend(
                adjacency[index]
                    .iter()
                    .copied()
                    .filter(|&next| !visited[next]),
            );
        }
        result
    };

    let side_a = component(a);
    let side_b = component(b);
    if side_a.contains(&b) || side_b.contains(&a) {
        return None;
    }
    Some((side_a, side_b))
}

#[inline]
fn dist(a: Vec3, b: Vec3) -> f32 {
    a.distance(b)
}

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
            if hi && hj {
                continue;
            }
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
        if visited[v] || v == exclude {
            continue;
        }
        visited[v] = true;
        cluster.push(v);
        if let Some(nei) = adj.get(&v) {
            for &u in nei {
                if !visited[u] && u != exclude {
                    q.push_back(u);
                }
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

/// Rotate the atoms on the chosen side of an actual covalent bond.
pub fn rotate_side_with_bonds(
    coords: &[Vec3],
    bonds: &[(usize, usize, f32)],
    a: usize,
    b: usize,
    which: RotateSide,
    angle_deg: f32,
) -> Option<Vec<Vec3>> {
    let axis_vec = *coords.get(b)? - *coords.get(a)?;
    let axis_len = axis_vec.length();
    if axis_len <= 1.0e-6 {
        return None;
    }
    let (side_a, side_b) = split_sides_from_bonds(coords.len(), bonds, a, b)?;
    let (center_index, rotating_set) = match which {
        RotateSide::A => (a, side_a),
        RotateSide::B => (b, side_b),
    };
    let center = coords[center_index];
    let axis = axis_vec / axis_len;
    let angle = angle_deg.to_radians();
    let mut new_coords = coords.to_vec();
    for index in rotating_set {
        if index == center_index {
            continue;
        }
        let offset = coords[index] - center;
        let rotated = offset * angle.cos()
            + axis.cross(offset) * angle.sin()
            + axis * axis.dot(offset) * (1.0 - angle.cos());
        new_coords[index] = center + rotated;
    }
    Some(new_coords)
}

/// Translate exactly the fragment containing the third-picked atom along the
/// selected covalent-bond axis.
pub fn translate_fragment_containing_atom(
    coords: &[Vec3],
    bonds: &[(usize, usize, f32)],
    a: usize,
    b: usize,
    selected_atom: usize,
    distance: f32,
) -> Option<Vec<Vec3>> {
    let axis_vec = *coords.get(b)? - *coords.get(a)?;
    let axis_len = axis_vec.length();
    if axis_len <= 1.0e-6 {
        return None;
    }
    let (side_a, side_b) = split_sides_from_bonds(coords.len(), bonds, a, b)?;
    let (direction, moving_set) = if side_a.contains(&selected_atom) {
        (-axis_vec / axis_len, side_a)
    } else if side_b.contains(&selected_atom) {
        (axis_vec / axis_len, side_b)
    } else {
        return None;
    };

    let mut new_coords = coords.to_vec();
    for index in moving_set {
        new_coords[index] += direction * distance;
    }
    Some(new_coords)
}
/// Bend the chosen covalent-bond side so `picked_idx` reaches the target axis angle.
pub fn bend_side_with_bonds(
    coords: &[Vec3],
    bonds: &[(usize, usize, f32)],
    a: usize,
    b: usize,
    which: RotateSide,
    picked_idx: usize,
    target_angle_deg: f32,
) -> Option<Vec<Vec3>> {
    let axis_vec = *coords.get(b)? - *coords.get(a)?;
    let axis_len = axis_vec.length();
    if axis_len <= 1.0e-6 {
        return None;
    }
    let (side_a, side_b) = split_sides_from_bonds(coords.len(), bonds, a, b)?;
    let (center_index, rotating_set) = match which {
        RotateSide::A => (a, side_a),
        RotateSide::B => (b, side_b),
    };
    if !rotating_set.contains(&picked_idx) {
        return None;
    }

    let center = coords[center_index];
    let axis = axis_vec / axis_len;
    let picked_offset = *coords.get(picked_idx)? - center;
    let picked_len = picked_offset.length();
    if picked_len <= 1.0e-6 {
        return None;
    }
    let picked_direction = picked_offset / picked_len;
    let current_angle = axis.dot(picked_direction).clamp(-1.0, 1.0).acos();
    let delta = target_angle_deg.to_radians() - current_angle;

    let mut rotation_axis = axis.cross(picked_direction);
    if rotation_axis.length() <= 1.0e-6 {
        rotation_axis = axis.cross(Vec3::X);
        if rotation_axis.length() <= 1.0e-6 {
            rotation_axis = axis.cross(Vec3::Y);
        }
        if rotation_axis.length() <= 1.0e-6 {
            return None;
        }
    }
    let rotation_axis = rotation_axis.normalize();
    let mut new_coords = coords.to_vec();
    for index in rotating_set {
        if index == center_index {
            continue;
        }
        let offset = coords[index] - center;
        let rotated = offset * delta.cos()
            + rotation_axis.cross(offset) * delta.sin()
            + rotation_axis * rotation_axis.dot(offset) * (1.0 - delta.cos());
        new_coords[index] = center + rotated;
    }
    Some(new_coords)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cutting_a_chain_bond_returns_the_requested_fragments() {
        let bonds = vec![(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)];
        let (side_a, side_b) = split_sides_from_bonds(4, &bonds, 1, 2).unwrap();
        assert_eq!(side_a, vec![1, 0]);
        assert_eq!(side_b, vec![2, 3]);
    }

    #[test]
    fn a_ring_bond_is_not_presented_as_a_rotatable_fragment() {
        let bonds = vec![(0, 1, 1.0), (1, 2, 1.0), (2, 0, 1.0)];
        assert!(split_sides_from_bonds(3, &bonds, 0, 1).is_none());
    }
    #[test]
    fn translation_moves_the_fragment_containing_the_third_pick() {
        let coords = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(3.0, 0.0, 0.0),
        ];
        let bonds = vec![(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)];

        let moved = translate_fragment_containing_atom(&coords, &bonds, 1, 2, 0, 0.5).unwrap();

        assert_eq!(moved[0], Vec3::new(-0.5, 0.0, 0.0));
        assert_eq!(moved[1], Vec3::new(0.5, 0.0, 0.0));
        assert_eq!(moved[2], coords[2]);
        assert_eq!(moved[3], coords[3]);
    }
}
