//! Orienting, mirroring and flipping the molecule.
//!
//! View operations that change the coordinates rather than the camera, because that is what a
//! figure or an exported file needs: a molecule lying flat in the xy plane stays flat when the
//! file is opened somewhere else, where a camera angle does not travel.
//!
//! # Chirality
//!
//! The distinction runs through this module. [`orient`] and [`flip`] are **proper** rotations and
//! leave a chiral molecule as it was; [`mirror`] is a reflection and gives the enantiomer. A bug
//! that let an improper rotation into [`orient`] would silently swap R for S in a structure on its
//! way to a calculation, so the sign of the rotation is checked rather than assumed, and the tests
//! measure a signed volume rather than trusting it.

use bevy::prelude::Vec3;
use ndarray::Array2;

use crate::molecule::atomic_mass_amu;
use crate::normalmode::linalg_shim::{Backend, LinAlg};

/// Which item of the floating "Orientation" menu is expanded into its choices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MenuItem {
    Orient,
    Mirror,
    Flip,
}

/// A coordinate plane.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Plane {
    #[default]
    Xy,
    Xz,
    Yz,
}

impl Plane {
    pub const ALL: [Plane; 3] = [Plane::Xy, Plane::Xz, Plane::Yz];

    pub fn label(self) -> &'static str {
        match self {
            Self::Xy => "xy",
            Self::Xz => "xz",
            Self::Yz => "yz",
        }
    }

    /// Which axis is perpendicular to the plane: 0 for x, 1 for y, 2 for z.
    fn normal_axis(self) -> usize {
        match self {
            Self::Xy => 2,
            Self::Xz => 1,
            Self::Yz => 0,
        }
    }

    /// The two in-plane axes, longest extent first.
    fn in_plane_axes(self) -> [usize; 2] {
        match self {
            Self::Xy => [0, 1],
            Self::Xz => [0, 2],
            Self::Yz => [1, 2],
        }
    }
}

/// A coordinate axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Axis {
    #[default]
    X,
    Y,
    Z,
}

impl Axis {
    pub const ALL: [Axis; 3] = [Axis::X, Axis::Y, Axis::Z];

    pub fn label(self) -> &'static str {
        match self {
            Self::X => "x",
            Self::Y => "y",
            Self::Z => "z",
        }
    }

    fn index(self) -> usize {
        match self {
            Self::X => 0,
            Self::Y => 1,
            Self::Z => 2,
        }
    }
}

/// The mass-weighted centre of the molecule.
///
/// Mass-weighted rather than geometric, so a heavy atom is not shifted about by a crowd of
/// hydrogens -- and so the axes below are the physical principal axes, the ones a rotational
/// spectrum would refer to.
fn centre_of_mass(atoms: &[String], pos: &[Vec3]) -> Vec3 {
    let mut total = 0.0_f64;
    let mut weighted = [0.0_f64; 3];
    for (symbol, p) in atoms.iter().zip(pos) {
        let mass = atomic_mass_amu(symbol).unwrap_or(12.0);
        total += mass;
        weighted[0] += mass * p.x as f64;
        weighted[1] += mass * p.y as f64;
        weighted[2] += mass * p.z as f64;
    }
    if total <= 0.0 {
        return Vec3::ZERO;
    }
    Vec3::new(
        (weighted[0] / total) as f32,
        (weighted[1] / total) as f32,
        (weighted[2] / total) as f32,
    )
}

/// Principal axes, as columns, ordered by increasing moment of inertia.
///
/// The smallest moment belongs to the axis the molecule is most extended along -- mass sits close
/// to it -- and the largest to the axis it is flattest against. So for a planar molecule the
/// largest-moment axis is the plane normal, which is what makes this the right tool for laying a
/// molecule flat.
fn principal_axes(atoms: &[String], pos: &[Vec3], centre: Vec3) -> Option<[[f64; 3]; 3]> {
    let mut inertia = Array2::<f64>::zeros((3, 3));
    for (symbol, p) in atoms.iter().zip(pos) {
        let mass = atomic_mass_amu(symbol).unwrap_or(12.0);
        let r = [
            (p.x - centre.x) as f64,
            (p.y - centre.y) as f64,
            (p.z - centre.z) as f64,
        ];
        let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        for i in 0..3 {
            for j in 0..3 {
                let kronecker = if i == j { 1.0 } else { 0.0 };
                inertia[(i, j)] += mass * (kronecker * r2 - r[i] * r[j]);
            }
        }
    }

    let solver = LinAlg::new(Backend::Auto).ok()?;
    let (values, vectors) = solver.eigh(&inertia).ok()?;
    // `eigh` returns ascending eigenvalues, so column 0 already has the smallest moment; sorting
    // explicitly costs nothing and does not depend on that promise.
    let mut order: Vec<usize> = (0..3).collect();
    order.sort_by(|a, b| values[*a].partial_cmp(&values[*b]).unwrap());
    Some(std::array::from_fn(|slot| {
        let column = order[slot];
        [
            vectors[(0, column)],
            vectors[(1, column)],
            vectors[(2, column)],
        ]
    }))
}

/// Lays the molecule in `plane`, centred on its centre of mass.
///
/// The axis it is flattest against becomes the plane normal and its longest axis the first
/// in-plane direction, so a planar molecule ends up flat in the plane and a long one runs along
/// it.
///
/// Always a proper rotation. An eigenvector set can come back left-handed, and applying it as-is
/// would mirror the molecule -- turning R into S in a structure on its way to a calculation, with
/// nothing on screen to show it. One axis is negated when that happens.
pub fn orient(atoms: &[String], pos: &[Vec3], plane: Plane) -> Vec<Vec3> {
    let centre = centre_of_mass(atoms, pos);
    let Some(axes) = principal_axes(atoms, pos, centre) else {
        return pos.to_vec();
    };

    // Row `k` of the rotation is the principal axis that should end up along coordinate `k`.
    let [first, second] = plane.in_plane_axes();
    let mut rows = [[0.0_f64; 3]; 3];
    rows[first] = axes[0]; // smallest moment: the molecule's long direction
    rows[second] = axes[1];
    rows[plane.normal_axis()] = axes[2]; // largest moment: the direction it is flattest in

    // A reflection has determinant -1. Negating one row makes it a rotation without changing
    // which axis lies where.
    if determinant(&rows) < 0.0 {
        for value in rows[plane.normal_axis()].iter_mut() {
            *value = -*value;
        }
    }

    pos.iter()
        .map(|p| {
            let r = [
                (p.x - centre.x) as f64,
                (p.y - centre.y) as f64,
                (p.z - centre.z) as f64,
            ];
            Vec3::new(
                dot(rows[0], r) as f32,
                dot(rows[1], r) as f32,
                dot(rows[2], r) as f32,
            )
        })
        .collect()
}

/// Reflects the molecule through `plane`, about its centre of mass.
///
/// This **is** a reflection: a chiral molecule comes back as its enantiomer. That is the point of
/// the operation, and the panel says so.
pub fn mirror(atoms: &[String], pos: &[Vec3], plane: Plane) -> Vec<Vec3> {
    let centre = centre_of_mass(atoms, pos);
    let axis = plane.normal_axis();
    pos.iter()
        .map(|p| {
            let mut r = [p.x - centre.x, p.y - centre.y, p.z - centre.z];
            r[axis] = -r[axis];
            Vec3::new(r[0] + centre.x, r[1] + centre.y, r[2] + centre.z)
        })
        .collect()
}

/// Turns the molecule over: a 180° rotation about `axis`, through its centre of mass.
///
/// A proper rotation, so a chiral molecule is unchanged -- it is the same molecule seen from the
/// other side, which is what distinguishes this from [`mirror`].
pub fn flip(atoms: &[String], pos: &[Vec3], axis: Axis) -> Vec<Vec3> {
    let centre = centre_of_mass(atoms, pos);
    let kept = axis.index();
    pos.iter()
        .map(|p| {
            let mut r = [p.x - centre.x, p.y - centre.y, p.z - centre.z];
            // 180° about an axis keeps that component and negates the other two.
            for (index, value) in r.iter_mut().enumerate() {
                if index != kept {
                    *value = -*value;
                }
            }
            Vec3::new(r[0] + centre.x, r[1] + centre.y, r[2] + centre.z)
        })
        .collect()
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn determinant(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Benzene, flat but lying in an awkward plane so orienting it has work to do.
    fn tilted_benzene() -> (Vec<String>, Vec<Vec3>) {
        let mut atoms = Vec::new();
        let mut pos = Vec::new();
        // A flat hexagon in its own plane, then rotated about x and y.
        let (sx, cx) = (0.6_f32.sin(), 0.6_f32.cos());
        let (sy, cy) = (0.9_f32.sin(), 0.9_f32.cos());
        for k in 0..6 {
            let a = std::f32::consts::TAU * k as f32 / 6.0;
            for (symbol, radius) in [("C", 1.39_f32), ("H", 2.48)] {
                let flat = Vec3::new(radius * a.cos(), radius * a.sin(), 0.0);
                // Rotate about x, then y.
                let after_x = Vec3::new(flat.x, flat.y * cx - flat.z * sx, flat.y * sx + flat.z * cx);
                let after_y = Vec3::new(
                    after_x.x * cy + after_x.z * sy,
                    after_x.y,
                    -after_x.x * sy + after_x.z * cy,
                );
                atoms.push(symbol.to_string());
                pos.push(after_y + Vec3::new(3.0, -2.0, 5.0));
            }
        }
        (atoms, pos)
    }

    /// A chiral centre: four different atoms about a carbon. The signed volume of three bonds
    /// tells R from S, so it is what the chirality tests measure.
    fn chiral_centre() -> (Vec<String>, Vec<Vec3>) {
        let atoms = ["C", "H", "F", "Cl", "Br"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let pos = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.63, 0.63, 0.63),
            Vec3::new(-0.63, -0.63, 0.63),
            Vec3::new(-0.63, 0.63, -0.63),
            Vec3::new(0.63, -0.63, -0.63),
        ];
        (atoms, pos)
    }

    /// The signed volume of the three bonds from atom 0: positive or negative according to
    /// handedness, and it flips under reflection but not under rotation.
    fn handedness(pos: &[Vec3]) -> f32 {
        let u = pos[1] - pos[0];
        let v = pos[2] - pos[0];
        let w = pos[3] - pos[0];
        u.cross(v).dot(w)
    }

    fn spread(pos: &[Vec3], axis: usize) -> f32 {
        let component = |p: &Vec3| match axis {
            0 => p.x,
            1 => p.y,
            _ => p.z,
        };
        let lo = pos.iter().map(component).fold(f32::MAX, f32::min);
        let hi = pos.iter().map(component).fold(f32::MIN, f32::max);
        hi - lo
    }

    fn distances(pos: &[Vec3]) -> Vec<f32> {
        let mut out = Vec::new();
        for i in 0..pos.len() {
            for j in (i + 1)..pos.len() {
                out.push((pos[i] - pos[j]).length());
            }
        }
        out
    }

    #[test]
    fn orienting_a_flat_molecule_lays_it_in_the_plane() {
        let (atoms, pos) = tilted_benzene();
        for plane in Plane::ALL {
            let oriented = orient(&atoms, &pos, plane);
            let normal = spread(&oriented, plane.normal_axis());
            assert!(
                normal < 1.0e-4,
                "{} plane: the molecule still spans {normal} across the normal",
                plane.label()
            );
            // And it really is spread out within the plane, not collapsed.
            for axis in plane.in_plane_axes() {
                assert!(spread(&oriented, axis) > 1.0);
            }
        }
    }

    /// Orienting must not deform the molecule: every internal distance survives.
    #[test]
    fn orienting_is_rigid() {
        let (atoms, pos) = tilted_benzene();
        let before = distances(&pos);
        for plane in Plane::ALL {
            let after = distances(&orient(&atoms, &pos, plane));
            for (a, b) in before.iter().zip(&after) {
                assert!((a - b).abs() < 1.0e-3, "a distance changed from {a} to {b}");
            }
        }
    }

    /// The one that matters most: orienting must not turn R into S. An eigenvector set can come
    /// back left-handed, and using it unchecked would mirror the molecule with nothing on screen
    /// to show it.
    #[test]
    fn orienting_preserves_chirality() {
        let (atoms, pos) = chiral_centre();
        let before = handedness(&pos);
        assert!(before.abs() > 0.1, "the test centre is not chiral");
        for plane in Plane::ALL {
            let after = handedness(&orient(&atoms, &pos, plane));
            assert!(
                before.signum() == after.signum(),
                "{} plane: handedness went from {before} to {after}",
                plane.label()
            );
        }
    }

    /// Mirroring is a reflection, so it must swap handedness. This is the operation's purpose.
    #[test]
    fn mirroring_gives_the_enantiomer() {
        let (atoms, pos) = chiral_centre();
        let before = handedness(&pos);
        for plane in Plane::ALL {
            let after = handedness(&mirror(&atoms, &pos, plane));
            assert!(
                before.signum() != after.signum(),
                "{} plane: handedness stayed {before}",
                plane.label()
            );
        }
    }

    #[test]
    fn mirroring_twice_returns_the_original() {
        let (atoms, pos) = tilted_benzene();
        for plane in Plane::ALL {
            let once = mirror(&atoms, &pos, plane);
            let twice = mirror(&atoms, &once, plane);
            for (a, b) in pos.iter().zip(&twice) {
                assert!((*a - *b).length() < 1.0e-4);
            }
        }
    }

    /// Flipping is a proper rotation, so handedness survives -- that is what distinguishes it from
    /// mirroring, and a chemist choosing between the two is relying on the difference.
    #[test]
    fn flipping_preserves_chirality_and_is_rigid() {
        let (atoms, pos) = chiral_centre();
        let before = handedness(&pos);
        let before_distances = distances(&pos);
        for axis in Axis::ALL {
            let flipped = flip(&atoms, &pos, axis);
            let after = handedness(&flipped);
            assert!(
                before.signum() == after.signum(),
                "{} axis: handedness went from {before} to {after}",
                axis.label()
            );
            for (a, b) in before_distances.iter().zip(&distances(&flipped)) {
                assert!((a - b).abs() < 1.0e-4);
            }
        }
    }

    #[test]
    fn flipping_twice_returns_the_original() {
        let (atoms, pos) = tilted_benzene();
        for axis in Axis::ALL {
            let once = flip(&atoms, &pos, axis);
            let twice = flip(&atoms, &once, axis);
            for (a, b) in pos.iter().zip(&twice) {
                assert!((*a - *b).length() < 1.0e-4);
            }
        }
    }

    /// A flip really does turn the molecule over rather than leaving it alone.
    #[test]
    fn flipping_actually_moves_the_molecule() {
        let (atoms, pos) = tilted_benzene();
        for axis in Axis::ALL {
            let flipped = flip(&atoms, &pos, axis);
            let moved = pos
                .iter()
                .zip(&flipped)
                .map(|(a, b)| (*a - *b).length())
                .fold(0.0_f32, f32::max);
            assert!(moved > 0.5, "{} axis: nothing moved", axis.label());
        }
    }

    /// An empty or single-atom molecule must come back unchanged rather than panicking on a
    /// degenerate inertia tensor.
    #[test]
    fn degenerate_molecules_are_left_alone() {
        assert!(orient(&[], &[], Plane::Xy).is_empty());
        let atoms = vec!["C".to_string()];
        let pos = vec![Vec3::new(1.0, 2.0, 3.0)];
        let oriented = orient(&atoms, &pos, Plane::Xy);
        assert_eq!(oriented.len(), 1);
        // A lone atom lands on the origin, being its own centre of mass.
        assert!(oriented[0].length() < 1.0e-5);
    }

    #[test]
    fn every_plane_and_axis_is_labelled() {
        assert_eq!(Plane::ALL.len(), 3);
        assert_eq!(Axis::ALL.len(), 3);
        for plane in Plane::ALL {
            assert!(!plane.label().is_empty());
            // The normal is not one of the in-plane axes.
            assert!(!plane.in_plane_axes().contains(&plane.normal_axis()));
        }
        for axis in Axis::ALL {
            assert!(!axis.label().is_empty());
        }
    }
}
