//! Ring puckering, so that five- and six-membered rings get their own conformers.
//!
//! A torsion-space search cannot reach them by turning bonds: rotating one ring dihedral prises
//! the ring open rather than giving another shape. The ring's degrees of freedom are instead its
//! **out-of-plane displacements**, and Cremer and Pople showed those decompose into a small,
//! complete set of puckering coordinates (*J. Am. Chem. Soc.* **1975**, *97*, 1354).
//!
//! For a six-membered ring that is a sphere: an amplitude `Q` with polar angles `(θ, φ)`. The
//! poles are the two chairs, the equator holds the boats and twist-boats, and the mid-latitudes
//! the half-chairs and envelopes. For a five-membered ring it is a circle: an amplitude `q` and a
//! single phase `φ` running round the pseudorotation itinerary of envelopes and twists.
//!
//! This is also the language sugar conformations are written in -- ⁴C₁, ¹C₄, ²S₀, E₃ -- so
//! sampling `(θ, φ)` is sampling exactly what a carbohydrate chemist means by a ring conformer.

use std::f64::consts::PI;

/// Out-of-plane amplitude used when generating a six-ring conformer, **in bohr**.
///
/// Cyclohexane's chair sits near 0.63 Å, which is this value. The unit matters: everything in the
/// conformer search carries coordinates in bohr, and a puckering amplitude is a length like any
/// other. Setting it to 0.63 as though it were Å makes every generated ring half as deep as it
/// should be -- which still looks like a pucker, and still relaxes to something, but is not enough
/// to swing a substituent between axial and equatorial, so the two chairs of a substituted ring
/// collapse into one conformer.
///
/// Generated geometries are relaxed immediately afterwards, so this only has to land the ring in
/// the right basin, not at the exact minimum.
pub const SIX_RING_AMPLITUDE: f64 = 1.19;

/// Out-of-plane amplitude for a five-ring, **in bohr**. Cyclopentane's pucker is near 0.42 Å.
pub const FIVE_RING_AMPLITUDE: f64 = 0.79;

/// A ring, as a cycle of atom indices in connected order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    /// Atom indices, each bonded to the next and the last bonded to the first.
    pub atoms: Vec<usize>,
}

impl Ring {
    pub fn size(&self) -> usize {
        self.atoms.len()
    }
}

/// Finds every five- and six-membered ring, each as a cycle in connected order.
///
/// Larger rings are left alone: their conformational space is not a small sphere or circle and
/// wants a different treatment. Smaller ones (three, four) are effectively rigid.
///
/// Fused systems yield each ring separately, which is the useful behaviour -- a decalin or a
/// pyranose fused to something else still puckers ring by ring, and the local optimiser resolves
/// any conflict between two rings sharing a bond.
pub fn find_rings(natoms: usize, bonds: &[[usize; 2]]) -> Vec<Ring> {
    let mut adjacency = vec![Vec::new(); natoms];
    for &[a, b] in bonds {
        if a < natoms && b < natoms {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
    }

    let mut found: Vec<Vec<usize>> = Vec::new();
    let mut path = Vec::new();
    for start in 0..natoms {
        // Each cycle is enumerated once by requiring its lowest-numbered atom to be the start.
        path.clear();
        path.push(start);
        walk(start, start, &adjacency, &mut path, &mut found);
    }

    found
        .into_iter()
        .map(|atoms| Ring { atoms })
        .collect()
}

/// Depth-first walk collecting cycles of length 5 and 6 that return to `start`.
fn walk(
    start: usize,
    current: usize,
    adjacency: &[Vec<usize>],
    path: &mut Vec<usize>,
    found: &mut Vec<Vec<usize>>,
) {
    const MAX_RING: usize = 6;
    for &next in &adjacency[current] {
        if next == start && path.len() >= 5 {
            let mut canonical = canonicalise(path);
            // Both directions round the ring describe the same cycle; keep one.
            if !found.iter().any(|seen| *seen == canonical) {
                let mut reversed = path.clone();
                reversed[1..].reverse();
                reversed = canonicalise(&reversed);
                if !found.iter().any(|seen| *seen == reversed) {
                    found.push(std::mem::take(&mut canonical));
                }
            }
            continue;
        }
        // Only extend through atoms above the start, so each ring is found from its lowest atom.
        if next <= start || path.contains(&next) || path.len() >= MAX_RING {
            continue;
        }
        path.push(next);
        walk(start, next, adjacency, path, found);
        path.pop();
    }
}

/// Rotates a cycle so it begins at its lowest atom, keeping the traversal order.
fn canonicalise(path: &[usize]) -> Vec<usize> {
    let lowest = path
        .iter()
        .enumerate()
        .min_by_key(|(_, atom)| **atom)
        .map(|(index, _)| index)
        .unwrap_or(0);
    let mut out = Vec::with_capacity(path.len());
    out.extend_from_slice(&path[lowest..]);
    out.extend_from_slice(&path[..lowest]);
    out
}

/// A ring conformer, as a point on the Cremer-Pople sphere (six-rings) or circle (five-rings).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pucker {
    /// Total out-of-plane amplitude, Å.
    pub amplitude: f64,
    /// Polar angle in degrees. Meaningless for a five-ring, where the space is a circle.
    pub theta_deg: f64,
    /// Phase angle in degrees, running round the pseudorotation itinerary.
    pub phi_deg: f64,
}

/// The canonical conformers of a ring of this size, as a list to sample from.
///
/// For six-rings this is the classical set: two chairs at the poles, boats and twist-boats
/// alternating round the equator, and half-chairs and envelopes between. For five-rings it is the
/// pseudorotation itinerary, envelopes and twists alternating every 18°.
///
/// A discrete list rather than continuous coordinates because ring conformers *are* discrete
/// basins -- there is no useful geometry halfway between a chair and a boat -- and because it
/// gives the search a gene it can mutate meaningfully.
pub fn canonical_conformers(size: usize) -> Vec<Pucker> {
    match size {
        6 => {
            let mut out = vec![
                // The two chairs. φ is meaningless at a pole.
                Pucker { amplitude: SIX_RING_AMPLITUDE, theta_deg: 0.0, phi_deg: 0.0 },
                Pucker { amplitude: SIX_RING_AMPLITUDE, theta_deg: 180.0, phi_deg: 0.0 },
            ];
            // Equator: boats at φ = 0, 60, ...; twist-boats at 30, 90, ...
            for step in 0..12 {
                out.push(Pucker {
                    amplitude: SIX_RING_AMPLITUDE,
                    theta_deg: 90.0,
                    phi_deg: 30.0 * step as f64,
                });
            }
            // Mid-latitudes: half-chairs and envelopes, on the way between chair and equator.
            for theta_deg in [50.8, 129.2] {
                for step in 0..12 {
                    out.push(Pucker {
                        amplitude: SIX_RING_AMPLITUDE,
                        theta_deg,
                        phi_deg: 30.0 * step as f64,
                    });
                }
            }
            out
        }
        5 => (0..20)
            .map(|step| Pucker {
                amplitude: FIVE_RING_AMPLITUDE,
                // The space is a circle; θ is carried but unused.
                theta_deg: 90.0,
                phi_deg: 18.0 * step as f64,
            })
            .collect(),
        _ => Vec::new(),
    }
}

/// The ring's mean plane: its centre and unit normal, by Cremer and Pople's construction.
fn mean_plane(coords: &[f64], ring: &[usize]) -> ([f64; 3], [f64; 3]) {
    let n = ring.len();
    let mut centre = [0.0; 3];
    for &atom in ring {
        for k in 0..3 {
            centre[k] += coords[3 * atom + k];
        }
    }
    for k in 0..3 {
        centre[k] /= n as f64;
    }

    // Two in-plane vectors built from the ring's own periodicity, so the plane does not depend on
    // which atom happens to be first.
    let mut r1 = [0.0; 3];
    let mut r2 = [0.0; 3];
    for (j, &atom) in ring.iter().enumerate() {
        let angle = 2.0 * PI * j as f64 / n as f64;
        for k in 0..3 {
            let offset = coords[3 * atom + k] - centre[k];
            r1[k] += offset * angle.sin();
            r2[k] += offset * angle.cos();
        }
    }
    let mut normal = [
        r1[1] * r2[2] - r1[2] * r2[1],
        r1[2] * r2[0] - r1[0] * r2[2],
        r1[0] * r2[1] - r1[1] * r2[0],
    ];
    let length = (normal.iter().map(|v| v * v).sum::<f64>()).sqrt();
    if length > 1.0e-12 {
        for v in normal.iter_mut() {
            *v /= length;
        }
    } else {
        // Degenerate: a collapsed ring has no meaningful plane. Any axis will do; the geometry is
        // rejected downstream as unphysical anyway.
        normal = [0.0, 0.0, 1.0];
    }
    (centre, normal)
}

/// Out-of-plane displacement of each ring atom from the mean plane, in ring order.
fn displacements(coords: &[f64], ring: &[usize]) -> Vec<f64> {
    let (centre, normal) = mean_plane(coords, ring);
    ring.iter()
        .map(|&atom| {
            (0..3)
                .map(|k| (coords[3 * atom + k] - centre[k]) * normal[k])
                .sum()
        })
        .collect()
}

/// Measures a ring's Cremer-Pople puckering coordinates from a geometry.
///
/// Returns `None` for a ring size this module does not treat.
pub fn measure(coords: &[f64], ring: &[usize]) -> Option<Pucker> {
    let n = ring.len();
    let z = displacements(coords, ring);
    match n {
        6 => {
            let (mut q2_cos, mut q2_sin, mut q3) = (0.0, 0.0, 0.0);
            for (j, zj) in z.iter().enumerate() {
                let angle = 4.0 * PI * j as f64 / 6.0;
                q2_cos += zj * angle.cos();
                q2_sin -= zj * angle.sin();
                q3 += zj * if j % 2 == 0 { 1.0 } else { -1.0 };
            }
            let scale2 = (2.0 / 6.0_f64).sqrt();
            let scale3 = (1.0 / 6.0_f64).sqrt();
            let (q2_cos, q2_sin, q3) = (scale2 * q2_cos, scale2 * q2_sin, scale3 * q3);
            let q2 = (q2_cos * q2_cos + q2_sin * q2_sin).sqrt();
            Some(Pucker {
                amplitude: (q2 * q2 + q3 * q3).sqrt(),
                theta_deg: q2.atan2(q3).to_degrees(),
                phi_deg: q2_sin.atan2(q2_cos).to_degrees().rem_euclid(360.0),
            })
        }
        5 => {
            let (mut q_cos, mut q_sin) = (0.0, 0.0);
            for (j, zj) in z.iter().enumerate() {
                let angle = 4.0 * PI * j as f64 / 5.0;
                q_cos += zj * angle.cos();
                q_sin -= zj * angle.sin();
            }
            let scale = (2.0 / 5.0_f64).sqrt();
            let (q_cos, q_sin) = (scale * q_cos, scale * q_sin);
            Some(Pucker {
                amplitude: (q_cos * q_cos + q_sin * q_sin).sqrt(),
                theta_deg: 90.0,
                phi_deg: q_sin.atan2(q_cos).to_degrees().rem_euclid(360.0),
            })
        }
        _ => None,
    }
}

/// The out-of-plane displacement each ring atom should have for a given conformer, in ring order.
///
/// This is the inverse of [`measure`]: Cremer and Pople's expansion read forwards.
pub fn target_displacements(size: usize, pucker: &Pucker) -> Option<Vec<f64>> {
    let phi = pucker.phi_deg.to_radians();
    match size {
        6 => {
            let theta = pucker.theta_deg.to_radians();
            let q2 = pucker.amplitude * theta.sin();
            let q3 = pucker.amplitude * theta.cos();
            let scale2 = (2.0 / 6.0_f64).sqrt();
            let scale3 = (1.0 / 6.0_f64).sqrt();
            Some(
                (0..6)
                    .map(|j| {
                        let angle = 4.0 * PI * j as f64 / 6.0;
                        scale2 * q2 * (phi + angle).cos()
                            + scale3 * q3 * if j % 2 == 0 { 1.0 } else { -1.0 }
                    })
                    .collect(),
            )
        }
        5 => {
            let scale = (2.0 / 5.0_f64).sqrt();
            Some(
                (0..5)
                    .map(|j| {
                        let angle = 4.0 * PI * j as f64 / 5.0;
                        scale * pucker.amplitude * (phi + angle).cos()
                    })
                    .collect(),
            )
        }
        _ => None,
    }
}

/// Which canonical conformer a measured pucker is closest to.
///
/// Compared as directions rather than as raw angles: for a six-ring the conformers sit on a
/// sphere, where φ wraps and becomes meaningless near the poles, so the honest distance is the
/// angle between unit vectors. A five-ring's space is a circle, so it is the wrapped difference
/// in φ.
pub fn nearest_conformer(size: usize, measured: &Pucker, conformers: &[Pucker]) -> usize {
    let direction = |p: &Pucker| {
        let (theta, phi) = (p.theta_deg.to_radians(), p.phi_deg.to_radians());
        [theta.sin() * phi.cos(), theta.sin() * phi.sin(), theta.cos()]
    };
    let mut best = 0;
    let mut best_score = f64::NEG_INFINITY;
    for (index, candidate) in conformers.iter().enumerate() {
        let score = if size == 6 {
            let (a, b) = (direction(measured), direction(candidate));
            a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
        } else {
            let difference = (measured.phi_deg - candidate.phi_deg).to_radians();
            difference.cos()
        };
        if score > best_score {
            best_score = score;
            best = index;
        }
    }
    best
}

/// Reshapes a ring in place to the given conformer, carrying its substituents along.
///
/// Each ring atom is swung about the axis through its two ring neighbours -- the classic corner
/// flap -- until it sits at the out-of-plane displacement the conformer calls for, and everything
/// hanging off it turns through the same rotation.
///
/// The rotation matters, and a translation will not do. Inverting a chair has to swap a
/// substituent between axial and equatorial, which is a change of *direction*: shifting the methyl
/// of methylcyclohexane along the ring normal leaves it pointing equatorially, so the geometry
/// relaxes straight back to the conformer it started from and the search finds one shape where
/// there are two. Rotating about the neighbour axis carries the substituent's orientation with the
/// ring, which is what makes axial and equatorial distinct.
///
/// The plane, the current displacements and the targets are all measured once from the incoming
/// geometry, so the atoms move independently rather than each one chasing a plane the previous one
/// has just shifted. That leaves bond lengths and angles slightly strained, which the local
/// optimisation after every generated geometry resolves; the job here is to land in the right
/// basin.
///
/// `branches[j]` holds the atoms to carry with ring atom `ring[j]`, excluding the ring itself.
pub fn apply_pucker(
    coords: &mut [f64],
    ring: &[usize],
    branches: &[Vec<usize>],
    pucker: &Pucker,
) -> bool {
    let n = ring.len();
    let Some(target) = target_displacements(n, pucker) else {
        return false;
    };
    if branches.len() != n {
        return false;
    }
    let (centre, normal) = mean_plane(coords, ring);
    let original = coords.to_vec();

    for j in 0..n {
        let atom = ring[j];
        let previous = ring[(j + n - 1) % n];
        let next = ring[(j + 1) % n];

        // The flap axis runs through the two ring neighbours, so swinging about it moves this
        // atom out of the plane without dragging the rest of the ring.
        let mut axis = [0.0; 3];
        let mut midpoint = [0.0; 3];
        for k in 0..3 {
            axis[k] = original[3 * next + k] - original[3 * previous + k];
            midpoint[k] = 0.5 * (original[3 * next + k] + original[3 * previous + k]);
        }
        let axis_length = (axis.iter().map(|v| v * v).sum::<f64>()).sqrt();
        if axis_length < 1.0e-9 {
            continue;
        }
        for v in axis.iter_mut() {
            *v /= axis_length;
        }

        // Split this atom's offset from the midpoint into the part along the axis, which the
        // rotation leaves alone, and the part it sweeps.
        let offset = [
            original[3 * atom] - midpoint[0],
            original[3 * atom + 1] - midpoint[1],
            original[3 * atom + 2] - midpoint[2],
        ];
        let along = dot(offset, axis);
        let perpendicular = [
            offset[0] - along * axis[0],
            offset[1] - along * axis[1],
            offset[2] - along * axis[2],
        ];
        let swept = cross(axis, perpendicular);

        // Out-of-plane displacement as a function of the rotation angle α:
        //   z(α) = c + a·cos α + b·sin α
        let c = dot(
            [
                midpoint[0] - centre[0],
                midpoint[1] - centre[1],
                midpoint[2] - centre[2],
            ],
            normal,
        ) + along * dot(axis, normal);
        let a = dot(perpendicular, normal);
        let b = dot(swept, normal);
        let radius = (a * a + b * b).sqrt();
        if radius < 1.0e-9 {
            continue;
        }

        // Solve c + R·cos(α − ψ) = target, taking the root that turns the ring least.
        let phase = b.atan2(a);
        let wanted = ((target[j] - c) / radius).clamp(-1.0, 1.0);
        let offset_angle = wanted.acos();
        let candidates = [phase + offset_angle, phase - offset_angle];
        let alpha = candidates
            .into_iter()
            .map(|angle| {
                let wrapped = angle.rem_euclid(2.0 * PI);
                if wrapped > PI {
                    wrapped - 2.0 * PI
                } else {
                    wrapped
                }
            })
            .min_by(|x, y| x.abs().partial_cmp(&y.abs()).unwrap())
            .unwrap_or(0.0);

        rotate_about(coords, &original, atom, midpoint, axis, alpha);
        for &carried in &branches[j] {
            rotate_about(coords, &original, carried, midpoint, axis, alpha);
        }
    }
    true
}

/// Rotates one atom about an axis through `midpoint`, reading its position from `original`.
fn rotate_about(
    coords: &mut [f64],
    original: &[f64],
    atom: usize,
    midpoint: [f64; 3],
    axis: [f64; 3],
    angle: f64,
) {
    let offset = [
        original[3 * atom] - midpoint[0],
        original[3 * atom + 1] - midpoint[1],
        original[3 * atom + 2] - midpoint[2],
    ];
    let along = dot(offset, axis);
    let perpendicular = [
        offset[0] - along * axis[0],
        offset[1] - along * axis[1],
        offset[2] - along * axis[2],
    ];
    let swept = cross(axis, perpendicular);
    let (cos, sin) = (angle.cos(), angle.sin());
    for k in 0..3 {
        coords[3 * atom + k] = midpoint[k]
            + along * axis[k]
            + perpendicular[k] * cos
            + swept[k] * sin;
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// For each ring atom, the atoms that should move with it: everything reachable without passing
/// through another ring atom.
pub fn substituent_branches(
    natoms: usize,
    bonds: &[[usize; 2]],
    ring: &[usize],
) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); natoms];
    for &[a, b] in bonds {
        if a < natoms && b < natoms {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
    }
    let in_ring: Vec<bool> = (0..natoms).map(|a| ring.contains(&a)).collect();

    ring.iter()
        .map(|&anchor| {
            let mut branch = Vec::new();
            let mut seen = in_ring.clone();
            let mut stack = vec![anchor];
            while let Some(current) = stack.pop() {
                for &next in &adjacency[current] {
                    // Never walk into another ring atom: that side belongs to it.
                    if seen[next] {
                        continue;
                    }
                    seen[next] = true;
                    branch.push(next);
                    stack.push(next);
                }
            }
            branch
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A planar hexagon, and the same ring puckered into an ideal chair.
    fn hexagon(pucker: f64) -> (Vec<f64>, Vec<usize>) {
        let mut coords = Vec::new();
        for j in 0..6 {
            let angle = 2.0 * PI * j as f64 / 6.0;
            coords.extend_from_slice(&[
                1.46 * angle.cos(),
                1.46 * angle.sin(),
                if j % 2 == 0 { pucker } else { -pucker },
            ]);
        }
        (coords, vec![0, 1, 2, 3, 4, 5])
    }

    fn distance(coords: &[f64], a: usize, b: usize) -> f64 {
        (0..3)
            .map(|k| (coords[3 * a + k] - coords[3 * b + k]).powi(2))
            .sum::<f64>()
            .sqrt()
    }

    fn ring_bonds(size: usize) -> Vec<[usize; 2]> {
        (0..size).map(|j| [j, (j + 1) % size]).collect()
    }

    #[test]
    fn finds_a_six_ring_and_a_five_ring() {
        let six = find_rings(6, &ring_bonds(6));
        assert_eq!(six.len(), 1);
        assert_eq!(six[0].size(), 6);

        let five = find_rings(5, &ring_bonds(5));
        assert_eq!(five.len(), 1);
        assert_eq!(five[0].size(), 5);
    }

    /// Consecutive atoms in the returned cycle must actually be bonded, or a pucker applied in
    /// that order would be nonsense.
    #[test]
    fn the_returned_cycle_is_in_connected_order() {
        let bonds = ring_bonds(6);
        let rings = find_rings(6, &bonds);
        let atoms = &rings[0].atoms;
        for j in 0..atoms.len() {
            let a = atoms[j];
            let b = atoms[(j + 1) % atoms.len()];
            assert!(
                bonds.iter().any(|&[x, y]| (x == a && y == b) || (x == b && y == a)),
                "{a} and {b} are adjacent in the cycle but not bonded"
            );
        }
    }

    #[test]
    fn a_chain_has_no_rings_and_large_rings_are_ignored() {
        let chain: Vec<[usize; 2]> = (0..5).map(|j| [j, j + 1]).collect();
        assert!(find_rings(6, &chain).is_empty());
        // A seven-ring is outside what this module treats.
        assert!(find_rings(7, &ring_bonds(7)).is_empty());
        // And a four-ring, which is effectively rigid.
        assert!(find_rings(4, &ring_bonds(4)).is_empty());
    }

    /// Fused rings must be reported separately -- decalin puckers ring by ring.
    #[test]
    fn fused_six_rings_are_found_separately() {
        // Two six-rings sharing the 0-1 bond.
        let mut bonds = vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]];
        bonds.extend_from_slice(&[[1, 6], [6, 7], [7, 8], [8, 9], [9, 0]]);
        let rings = find_rings(10, &bonds);
        assert_eq!(rings.len(), 2, "found {rings:?}");
        assert!(rings.iter().all(|r| r.size() == 6));
    }

    #[test]
    fn a_planar_ring_has_no_pucker() {
        let (coords, ring) = hexagon(0.0);
        let measured = measure(&coords, &ring).expect("six-ring");
        assert!(measured.amplitude < 1.0e-9, "amplitude {}", measured.amplitude);
    }

    /// An alternating up-down hexagon is a chair, which sits at a pole of the sphere: θ near 0 or
    /// 180, and all the amplitude in `q3`.
    #[test]
    fn an_alternating_ring_measures_as_a_chair() {
        let (coords, ring) = hexagon(0.25);
        let measured = measure(&coords, &ring).expect("six-ring");
        assert!(measured.amplitude > 0.1, "amplitude {}", measured.amplitude);
        let theta = measured.theta_deg;
        assert!(
            theta < 5.0 || theta > 175.0,
            "a chair must be at a pole, got θ = {theta}"
        );
    }

    /// The expansion and its inverse must agree exactly. Tested by placing the ring atoms at the
    /// displacements [`target_displacements`] asks for and measuring them back, which isolates the
    /// maths from how [`apply_pucker`] chooses to move the atoms.
    #[test]
    fn the_expansion_and_its_inverse_agree() {
        for size in [5usize, 6] {
            for wanted in canonical_conformers(size) {
                let target = target_displacements(size, &wanted).expect("known size");
                let ring: Vec<usize> = (0..size).collect();

                // A flat ring, then each atom lifted to its target displacement **along the
                // plane's own normal**. Cremer and Pople define that normal as R₁ × R₂, which for
                // a ring numbered counter-clockwise in the xy-plane points along −z -- so
                // displacing naively along +z would measure back as the conformer on the opposite
                // side of the sphere. Everything inside this module shares one `mean_plane`, so
                // that convention never leaks out; a test that builds a geometry by hand has to
                // respect it.
                let mut coords = Vec::new();
                for j in 0..size {
                    let angle = 2.0 * PI * j as f64 / size as f64;
                    coords.extend_from_slice(&[1.46 * angle.cos(), 1.46 * angle.sin(), 0.0]);
                }
                let (_, normal) = mean_plane(&coords, &ring);
                for (j, z) in target.iter().enumerate() {
                    for k in 0..3 {
                        coords[3 * j + k] += z * normal[k];
                    }
                }
                let got = measure(&coords, &ring).expect("known size");

                assert!(
                    (got.amplitude - wanted.amplitude).abs() < 1.0e-9,
                    "size {size}: amplitude {} vs {}",
                    got.amplitude,
                    wanted.amplitude
                );
                if size == 6 {
                    assert!(
                        (got.theta_deg - wanted.theta_deg).abs() < 1.0e-6,
                        "size {size}: θ {} vs {}",
                        got.theta_deg,
                        wanted.theta_deg
                    );
                }
                let at_pole = size == 6 && (wanted.theta_deg < 1.0 || wanted.theta_deg > 179.0);
                if !at_pole {
                    let difference = (got.phi_deg - wanted.phi_deg).rem_euclid(360.0);
                    let difference = difference.min(360.0 - difference);
                    assert!(
                        difference < 1.0e-6,
                        "size {size}: φ {} vs {}",
                        got.phi_deg,
                        wanted.phi_deg
                    );
                }
            }
        }
    }

    /// Applying a conformer must land in *its* basin, which is all the search needs -- every
    /// generated geometry is relaxed immediately afterwards.
    ///
    /// It does not land exactly on the requested point, and cannot: each atom swings about its two
    /// ring neighbours as they were before the move, so by the time the last atom has swung the
    /// first one's neighbours have shifted. Requiring exactness here would be requiring the wrong
    /// thing; requiring the right *identity* is the real test.
    #[test]
    fn applying_a_conformer_lands_in_its_basin() {
        for size in [5usize, 6] {
            let conformers = canonical_conformers(size);
            let ring: Vec<usize> = (0..size).collect();
            let branches = vec![Vec::new(); size];

            for (index, wanted) in conformers.iter().enumerate() {
                let mut coords = Vec::new();
                for j in 0..size {
                    let angle = 2.0 * PI * j as f64 / size as f64;
                    coords.extend_from_slice(&[1.46 * angle.cos(), 1.46 * angle.sin(), 0.0]);
                }
                assert!(apply_pucker(&mut coords, &ring, &branches, wanted));
                let got = measure(&coords, &ring).expect("known size");

                // The amplitude comes out a little short of the request, but the ring is properly
                // puckered rather than nearly flat.
                assert!(
                    got.amplitude > 0.7 * wanted.amplitude,
                    "size {size}, conformer {index}: amplitude {:.3} against {:.3} requested",
                    got.amplitude,
                    wanted.amplitude
                );
                // And it is nearest to the conformer that was asked for, not a neighbouring one.
                let nearest = nearest_conformer(size, &got, &conformers);
                assert_eq!(
                    nearest, index,
                    "size {size}: asked for conformer {index} (θ {:.0}°, φ {:.0}°), landed \
                     nearest {nearest} (θ {:.1}°, φ {:.1}°)",
                    wanted.theta_deg, wanted.phi_deg, got.theta_deg, got.phi_deg
                );
            }
        }
    }

    /// The two chairs must be genuinely different geometries, or sampling them is pointless.
    #[test]
    fn the_two_chairs_are_distinct() {
        let conformers = canonical_conformers(6);
        let ring: Vec<usize> = (0..6).collect();
        let branches = vec![Vec::new(); 6];

        let mut shapes = Vec::new();
        for pucker in [conformers[0], conformers[1]] {
            let (mut coords, _) = hexagon(0.0);
            assert!(apply_pucker(&mut coords, &ring, &branches, &pucker));
            shapes.push(displacements(&coords, &ring));
        }
        // One chair is the other inverted, so every displacement flips sign.
        for (a, b) in shapes[0].iter().zip(&shapes[1]) {
            assert!((a + b).abs() < 1.0e-9, "{a} and {b} are not inverses");
            assert!(a.abs() > 0.1, "chair displacement {a} is too small");
        }
    }

    #[test]
    fn there_are_the_expected_numbers_of_canonical_conformers() {
        // Two chairs, twelve on the equator, twenty-four at the mid-latitudes.
        assert_eq!(canonical_conformers(6).len(), 38);
        // The five-ring pseudorotation itinerary, every 18°.
        assert_eq!(canonical_conformers(5).len(), 20);
        assert!(canonical_conformers(7).is_empty());
    }

    /// A substituent must travel with the ring atom it hangs from, rigidly.
    #[test]
    fn substituents_follow_their_ring_atom() {
        let mut bonds = ring_bonds(6);
        bonds.push([0, 6]);
        let ring: Vec<usize> = (0..6).collect();
        let branches = substituent_branches(7, &bonds, &ring);
        assert_eq!(branches[0], vec![6], "atom 6 hangs off ring atom 0");
        assert!(branches[1..].iter().all(|b| b.is_empty()));

        let (mut coords, _) = hexagon(0.0);
        coords.extend_from_slice(&[1.46, 0.0, 1.09]); // the substituent, above atom 0
        let bond_before = distance(&coords, 0, 6);

        let chair = canonical_conformers(6)[0];
        assert!(apply_pucker(&mut coords, &ring, &branches, &chair));

        // The rotation is rigid, so the bond it hangs by is unchanged. Its displacement is *not*
        // the same as the ring atom's -- it sits further from the flap axis and sweeps further --
        // which is exactly what carries axial and equatorial character.
        let bond_after = distance(&coords, 0, 6);
        assert!(
            (bond_after - bond_before).abs() < 1.0e-9,
            "the substituent bond went from {bond_before} to {bond_after}"
        );
    }

    /// The decisive test for substituted rings: the two chairs must hold a substituent at
    /// genuinely different angles to the ring. That is the axial/equatorial distinction -- axial
    /// points along the ring normal, equatorial lies close to the ring plane -- and it is why the
    /// pucker rotates a substituent rather than translating it. A translation leaves the methyl of
    /// methylcyclohexane pointing the same way and the two chairs collapse into one conformer.
    #[test]
    fn the_two_chairs_hold_a_substituent_at_different_angles() {
        let mut bonds = ring_bonds(6);
        bonds.push([0, 6]);
        let ring: Vec<usize> = (0..6).collect();
        let branches = substituent_branches(7, &bonds, &ring);
        let conformers = canonical_conformers(6);

        let mut angles = Vec::new();
        for pucker in [conformers[0], conformers[1]] {
            let (mut coords, _) = hexagon(0.0);
            coords.extend_from_slice(&[2.30, 0.0, 0.55]);
            assert!(apply_pucker(&mut coords, &ring, &branches, &pucker));

            // Angle between the ring-to-substituent bond and the ring normal.
            let (_, normal) = mean_plane(&coords, &ring);
            let bond = [
                coords[3 * 6] - coords[0],
                coords[3 * 6 + 1] - coords[1],
                coords[3 * 6 + 2] - coords[2],
            ];
            let length = (bond.iter().map(|v| v * v).sum::<f64>()).sqrt();
            let cosine = dot(bond, normal) / length;
            angles.push(cosine.clamp(-1.0, 1.0).acos().to_degrees());
        }

        let separation = (angles[0] - angles[1]).abs();
        assert!(
            separation > 20.0,
            "the two chairs hold the substituent at {:.1}° and {:.1}° -- too alike to be axial \
             and equatorial",
            angles[0],
            angles[1]
        );
    }

    /// A branch must stop at the next ring atom, or puckering one ring atom would drag the whole
    /// ring with it.
    #[test]
    fn branches_do_not_leak_through_the_ring() {
        let bonds = ring_bonds(6);
        let ring: Vec<usize> = (0..6).collect();
        let branches = substituent_branches(6, &bonds, &ring);
        assert!(
            branches.iter().all(|b| b.is_empty()),
            "a bare ring has no substituents, got {branches:?}"
        );
    }
}
