//! Reference-frame bookkeeping for attaching a fragment to an existing atom.
//!
//! A fragment is stored as a self-contained Z-matrix, so its first three atoms
//! necessarily have incomplete internal coordinates: atom 1 has no references
//! at all, atom 2 has no angle or dihedral reference, and atom 3 has no
//! dihedral reference.  Splicing that into a host molecule without filling
//! those in leaves them `None`, and `zmat_to_xyz` then falls back to atom 1 of
//! the *whole molecule* -- an arbitrary distant atom.  The fragment's internal
//! geometry gets rebuilt in an unrelated frame, which is what makes an inserted
//! fragment fold back into the structure it was attached to.
//!
//! The fix follows Molden's `AddFrag` (`src/xwin.c`): the missing references
//! are inherited from the attachment atom's own chain, walking back up the host
//! Z-matrix.
//!
//! ```text
//!   Molden, for attachment row `row`:
//!     row == 0 -> angle_ref = 2, dihedral_ref = 3
//!     row == 1 -> angle_ref = 1, dihedral_ref = 3
//!     else     -> angle_ref = iz(row,0)   (the row's own bond reference)
//!                 dihedral_ref = iz(row,1) (the row's own angle reference)
//! ```

use super::zmat2xyz::ZAtom;
use crate::molecule::covalent_radius_angstrom;
use bevy::prelude::*;

/// Angle and dihedral references to use for an atom bonded to `attach`.
///
/// Indices are 0-based; `None` means the host is too small to supply one.
/// Mirrors Molden's rule, including its special cases for the first two rows,
/// which have no chain of their own to inherit.
pub fn host_reference_chain(zmat: &[ZAtom], attach: usize) -> (Option<usize>, Option<usize>) {
    let n = zmat.len();
    let valid = |candidate: usize| (candidate < n && candidate != attach).then_some(candidate);

    match attach {
        0 => (valid(1), valid(2)),
        1 => (valid(0), valid(2)),
        _ => {
            let angle = zmat[attach].bond_ref.and_then(|r| valid(r - 1));
            let dihedral = zmat[attach].angle_ref.and_then(|r| valid(r - 1));
            // Fall back to any earlier atoms rather than leaving a hole; an
            // unfilled reference is what caused the original defect.
            let angle = angle.or_else(|| (0..n).find(|&i| i != attach));
            let dihedral = dihedral.or_else(|| (0..n).find(|&i| Some(i) != angle && i != attach));
            (angle, dihedral)
        }
    }
}

/// Complete the references of a fragment being spliced onto `attach`.
///
/// `atoms` is the fragment's own Z-matrix; on return every reference is a
/// 1-based index into the combined Z-matrix.  `positions` gives each fragment
/// atom's eventual 0-based index in the combined matrix, so a fragment-internal
/// reference can be translated, and `attach` is the 0-based host atom the
/// fragment hangs off.
///
/// `connector_overrides` lets a caller pin the connector's angle and dihedral
/// references, which is what the pick-three-atoms workflow supplies.
pub fn assign_fragment_references(
    atoms: &mut [ZAtom],
    positions: &[usize],
    host: &[ZAtom],
    attach: usize,
    connector_overrides: (Option<usize>, Option<usize>),
) {
    let (chain_angle, chain_dihedral) = host_reference_chain(host, attach);
    let to_one_based = |index: usize| index + 1;

    for (i, atom) in atoms.iter_mut().enumerate() {
        // Translate references that point inside the fragment.
        let translate = |reference: Option<usize>| -> Option<usize> {
            reference
                .and_then(|r| positions.get(r - 1).copied())
                .map(to_one_based)
        };

        match i {
            0 => {
                // The connector hangs off the host atom itself.
                atom.bond_ref = Some(to_one_based(attach));
                atom.angle_ref = connector_overrides.0.or(chain_angle).map(to_one_based);
                atom.dihedral_ref = connector_overrides.1.or(chain_dihedral).map(to_one_based);
            }
            1 => {
                // Bonded to the connector; the missing angle and dihedral
                // references come from the host, per Molden.
                atom.bond_ref = translate(atom.bond_ref).or(Some(to_one_based(positions[0])));
                atom.angle_ref = translate(atom.angle_ref).or(Some(to_one_based(attach)));
                atom.dihedral_ref = translate(atom.dihedral_ref)
                    .or(chain_angle.map(to_one_based))
                    .or(Some(to_one_based(attach)));
            }
            2 => {
                atom.bond_ref = translate(atom.bond_ref).or(Some(to_one_based(positions[0])));
                atom.angle_ref = translate(atom.angle_ref).or(Some(to_one_based(positions[0])));
                atom.dihedral_ref = translate(atom.dihedral_ref).or(Some(to_one_based(attach)));
            }
            _ => {
                // Deeper atoms are fully referenced within the fragment.
                atom.bond_ref = translate(atom.bond_ref).or(Some(to_one_based(positions[0])));
                atom.angle_ref = translate(atom.angle_ref).or(Some(to_one_based(attach)));
                atom.dihedral_ref = translate(atom.dihedral_ref)
                    .or(chain_angle.map(to_one_based))
                    .or(Some(to_one_based(attach)));
            }
        }
    }
}

/// Smallest distance from any atom in `group` to any host atom not in `ignore`.
///
/// Used while searching placements, where the fragment's coordinates exist
/// separately from the host's rather than in one combined array.
pub fn clearance_against(group: &[Vec3], host: &[Vec3], ignore: &[usize]) -> f32 {
    let mut smallest = f32::MAX;
    for &p in group {
        for (index, &h) in host.iter().enumerate() {
            if ignore.contains(&index) {
                continue;
            }
            smallest = smallest.min(p.distance(h));
        }
    }
    if smallest == f32::MAX {
        0.0
    } else {
        smallest
    }
}

/// Host atoms a fragment is legitimately close to: the attachment atom and the
/// two it borrows its reference frame from.
///
/// Scoring a placement without excluding these would penalise a correct
/// geometry simply for being bonded.
pub fn attachment_partners(zmat: &[ZAtom], attach: usize) -> Vec<usize> {
    let mut partners = vec![attach];
    if let Some(atom) = zmat.get(attach) {
        for reference in [atom.bond_ref, atom.angle_ref] {
            if let Some(r) = reference {
                if r >= 1 && r - 1 < zmat.len() {
                    partners.push(r - 1);
                }
            }
        }
    }
    partners.sort_unstable();
    partners.dedup();
    partners
}

/// Smallest distance between any atom of `group` and any atom outside it.
///
/// `ignore` names host atoms that are legitimately close -- the attachment atom
/// and its reference partners -- so a correct placement is not penalised for
/// being bonded.
pub fn min_clearance(coords: &[Vec3], group: &[usize], ignore: &[usize]) -> f32 {
    if coords.is_empty() || group.is_empty() {
        return 0.0;
    }
    let mut in_group = vec![false; coords.len()];
    for &index in group {
        if index < in_group.len() {
            in_group[index] = true;
        }
    }

    let mut smallest = f32::MAX;
    for &index in group {
        if index >= coords.len() {
            continue;
        }
        let position = coords[index];
        for (other, &other_position) in coords.iter().enumerate() {
            if in_group[other] || ignore.contains(&other) {
                continue;
            }
            smallest = smallest.min(position.distance(other_position));
        }
    }
    if smallest == f32::MAX {
        0.0
    } else {
        smallest
    }
}

/// Unit direction of a new bond from `a`, given internal coordinates measured
/// against `b` (angle reference) and `c` (dihedral reference).
///
/// This is the exact inverse of the frame `zmat_to_xyz` builds, so a direction
/// round-trips through [`internal_from_direction`] without drift.  Deriving the
/// pair rather than reusing the fragment's own numbers is the whole point: a
/// fragment's first atoms carry placeholder zeros for coordinates it had no
/// reference to measure, and a 0-degree angle plants an atom straight back
/// along the bond it was attached by.
pub fn direction_from_internal(
    a: Vec3,
    b: Vec3,
    c: Vec3,
    angle_deg: f64,
    dihedral_deg: f64,
) -> Vec3 {
    let e1 = (a - b).normalize_or_zero();
    let e2 = (b - c).normalize_or_zero();
    let mut n = e1.cross(e2);
    n = if n.length() < 1.0e-6 {
        any_perpendicular(e1)
    } else {
        n.normalize()
    };
    let m = n.cross(e1);

    let theta = (angle_deg as f32).to_radians();
    let phi = (dihedral_deg as f32).to_radians();
    (-theta.cos() * e1 + theta.sin() * phi.cos() * m + theta.sin() * phi.sin() * n)
        .normalize_or_zero()
}

/// Angle and dihedral, in degrees, that reproduce `direction` in the frame of
/// `a`, `b`, `c`.  Inverse of [`direction_from_internal`].
pub fn internal_from_direction(a: Vec3, b: Vec3, c: Vec3, direction: Vec3) -> (f64, f64) {
    let e1 = (a - b).normalize_or_zero();
    let e2 = (b - c).normalize_or_zero();
    let mut n = e1.cross(e2);
    n = if n.length() < 1.0e-6 {
        any_perpendicular(e1)
    } else {
        n.normalize()
    };
    let m = n.cross(e1);

    let v = direction.normalize_or_zero();
    let angle = (-v.dot(e1)).clamp(-1.0, 1.0).acos().to_degrees();
    let dihedral = v.dot(n).atan2(v.dot(m)).to_degrees();
    (angle as f64, dihedral.rem_euclid(360.0) as f64)
}

fn any_perpendicular(v: Vec3) -> Vec3 {
    let candidate = if v.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    v.cross(candidate).normalize_or_zero()
}

/// Place a fragment rigidly in world space.
///
/// `local` holds the fragment's own coordinates with the connector at index 0.
/// The connector's free valence -- the direction with no fragment atoms in it --
/// is aimed back at the host, so the fragment body ends up along `direction`.
/// `spin_deg` rotates the fragment about that bond, which is a free internal
/// rotation and therefore the natural parameter to search for clearance.
/// The exact rigid rotation `place_fragment` applies to a fragment's local
/// coordinates, exposed on its own so auxiliary vectors (a plane normal, in
/// particular) can be carried through the identical transform rather than a
/// separately re-derived one -- re-deriving it from a different, fake local
/// shape would generally produce a different "open valence" and therefore a
/// different rotation.
///
/// `open` is the connector's own "missing" direction (see
/// [`fragment_open_valence`]) and must be supplied by the caller, computed
/// from the connector's actual bonded neighbours -- not derived here from
/// every other atom in `local`, which for a multi-atom fragment like
/// `-CH=CH2` would wrongly include atoms that are two bonds away from the
/// connector (the terminal `=CH2`'s own hydrogens) as if they were its direct
/// substituents.
pub fn fragment_rotation(local: &[Vec3], open: Vec3, direction: Vec3, spin_deg: f32) -> Quat {
    if local.is_empty() {
        return Quat::IDENTITY;
    }
    let align = Quat::from_rotation_arc(open, -direction.normalize_or_zero());
    let spin = Quat::from_axis_angle(direction.normalize_or_zero(), spin_deg.to_radians());
    spin * align
}

/// The connector's own "open valence": the direction opposite the sum of its
/// *actual bonded* local neighbours (by covalent-radius distance, since a
/// fragment's local shape carries no Z-matrix bond references of its own).
///
/// Deliberately not "every other atom in the fragment": for a multi-atom
/// fragment such as `-CH=CH2`, the terminal `=CH2` carbon's own two hydrogens
/// are two bonds away from the connector and are not among its substituents,
/// even though they are elsewhere in the same `local` array. Summing them in
/// anyway skews the computed direction by tens of degrees -- exactly the
/// defect that made a vinyl group's own connector land at the wrong angle
/// from its host, before any second fragment was even attached.
/// The direction a connector's remaining bond should point, given the bonds it
/// already has and the angle its element wants between them.
///
/// The naive answer -- negate the sum of the existing bond vectors -- is right
/// for three bonds, and for a planar centre with two, but wrong elsewhere:
///
/// * With **one** bond it gives exactly 180 degrees, which is only correct for
///   a genuine sp centre. On an -OH fragment that put the hydroxyl hydrogen in
///   a straight line with the atom the oxygen bonds to.
/// * With **two** bonds it lies in their plane, which is right for sp2 but
///   flattens a pyramidal centre such as an amine nitrogen.
///
/// `ideal_angle_deg` is the angle wanted between the new bond and each
/// existing one; see `bond_order::ideal_bond_angle_deg`. `max_bond_order` is
/// the highest estimated order among the bonds the connector already has,
/// which is what distinguishes a genuinely sp or sp2 centre -- an alkyne
/// carbon, an aromatic nitrogen -- from a partly-built sp3 one that merely
/// happens to have few bonds so far.
pub fn open_valence_direction(
    neighbor_dirs: &[Vec3],
    ideal_angle_deg: f32,
    max_bond_order: f64,
) -> Option<Vec3> {
    /// Above this, two existing bonds are taken to mark a planar centre, whose
    /// third bond belongs in their plane. Below it the centre is pyramidal.
    const PLANAR_THRESHOLD_DEG: f32 = 115.0;
    /// A bond at least this strong forces the centre flat: a double bond, or
    /// the ~1.5 of an aromatic ring.
    const PLANAR_BOND_ORDER: f64 = 1.25;
    /// A double or triple bond on a centre with one other bond leaves it
    /// linear -- an alkyne or a nitrile has no choice about it.
    const LINEAR_BOND_ORDER: f64 = 2.0;

    let dirs: Vec<Vec3> = neighbor_dirs
        .iter()
        .map(|d| d.normalize_or_zero())
        .filter(|d| d.length() > 0.5)
        .collect();
    let theta = ideal_angle_deg.to_radians();

    match dirs.len() {
        0 => None,
        1 => {
            let n = dirs[0];
            if max_bond_order >= LINEAR_BOND_ORDER {
                // sp: the remaining bond is straight through.
                return Some(-n);
            }
            // Any perpendicular will do: the fragment is rigid-body rotated
            // into place afterwards, so only the angle matters, not which way
            // round the remaining freedom is spent.
            let p = any_perpendicular(n);
            (p.length() > 0.5).then(|| (n * theta.cos() + p * theta.sin()).normalize_or_zero())
        }
        2 => {
            let (a, b) = (dirs[0], dirs[1]);
            let alpha = a.dot(b).clamp(-1.0, 1.0).acos();
            let bisector = (a + b).normalize_or_zero();
            if bisector.length() < 0.5 {
                return None;
            }
            if alpha.to_degrees() >= PLANAR_THRESHOLD_DEG
                || max_bond_order >= PLANAR_BOND_ORDER
            {
                // Planar: the remaining bond opposes the other two, in plane.
                // A narrow angle is not enough to call a centre pyramidal --
                // pyrrole's nitrogen spans only 107 degrees and is flat -- so
                // the bond orders get the final say.
                return Some(-bisector);
            }
            let normal = a.cross(b).normalize_or_zero();
            if normal.length() < 0.5 {
                return None;
            }
            // Out of plane, at `theta` from both. With `s` the bisector and
            // `m` the normal, a direction `-s cos(phi) + m sin(phi)` makes an
            // angle with each existing bond whose cosine is
            // `-cos(phi) cos(alpha/2)`, so solve that for `phi`.
            let cos_phi = (-theta.cos() / (alpha * 0.5).cos()).clamp(-1.0, 1.0);
            let sin_phi = (1.0 - cos_phi * cos_phi).max(0.0).sqrt();
            Some((-bisector * cos_phi + normal * sin_phi).normalize_or_zero())
        }
        _ => {
            // Three or more: their sum already points away from the gap.
            let sum: Vec3 = dirs.iter().copied().sum();
            (sum.length() > 1.0e-4).then(|| -sum.normalize())
        }
    }
}

/// Unit directions making `angle_deg` with `axis`, spread evenly around it.
///
/// A centre with a single existing bond leaves a whole cone of chemically
/// equivalent directions for the next one: the angle is fixed but the rotation
/// about the existing bond is free. Enumerating the cone lets the caller spend
/// that freedom on whatever clears the rest of the molecule best, instead of
/// on an arbitrary perpendicular.
pub fn cone_directions(axis: Vec3, angle_deg: f32, count: usize) -> Vec<Vec3> {
    let axis = axis.normalize_or_zero();
    let perpendicular = any_perpendicular(axis);
    if axis.length() < 0.5 || perpendicular.length() < 0.5 || count == 0 {
        return Vec::new();
    }
    let other = axis.cross(perpendicular).normalize_or_zero();
    let theta = angle_deg.to_radians();
    (0..count)
        .map(|k| {
            let phi = std::f32::consts::TAU * k as f32 / count as f32;
            (axis * theta.cos()
                + (perpendicular * phi.cos() + other * phi.sin()) * theta.sin())
            .normalize_or_zero()
        })
        .collect()
}

pub fn fragment_open_valence(local: &[Vec3], connector_symbol: &str, local_symbols: &[&str]) -> Vec3 {
    if local.len() < 2 {
        return -Vec3::Z;
    }
    let connector = local[0];
    let others: Vec<(Vec3, &str)> = local[1..]
        .iter()
        .copied()
        .zip(local_symbols.iter().copied())
        .collect();
    let bonded = geometric_bonded_neighbors(connector, connector_symbol, &others);
    let dirs: Vec<Vec3> = bonded
        .iter()
        .map(|&p| (p - connector).normalize_or_zero())
        .collect();

    // The strongest bond the connector already has, estimated from its length
    // in the fragment's own geometry. This is what tells an alkyne carbon or
    // an aromatic nitrogen apart from a merely under-coordinated sp3 centre.
    let max_order = bonded
        .iter()
        .filter_map(|&p| {
            let symbol = others
                .iter()
                .find(|(q, _)| q.distance(p) < 1.0e-4)
                .map(|(_, s)| *s)?;
            Some(crate::bond_order::estimate_bond_order(
                connector_symbol,
                symbol,
                connector.distance(p),
            ))
        })
        .fold(0.0f64, f64::max);

    // The angle this element wants once the connector's own bond is made.
    let ideal = crate::bond_order::ideal_bond_angle_deg(connector_symbol, dirs.len() + 1) as f32;
    open_valence_direction(&dirs, ideal, max_order).unwrap_or(-Vec3::Z)
}

pub fn place_fragment(local: &[Vec3], open: Vec3, origin: Vec3, direction: Vec3, spin_deg: f32) -> Vec<Vec3> {
    if local.is_empty() {
        return Vec::new();
    }
    let connector = local[0];
    let rotation = fragment_rotation(local, open, direction, spin_deg);

    local
        .iter()
        .map(|&p| origin + rotation * (p - connector))
        .collect()
}

/// Rewrite a fragment's internal coordinates to match its placed geometry.
///
/// Every value is measured from `coords`, so the Z-matrix reproduces the
/// placement exactly instead of inheriting the placeholders a self-contained
/// fragment necessarily carries.
pub fn derive_internals(coords: &[Vec3], atoms: &mut [ZAtom], positions: &[usize]) {
    // Indexed directly by `positions[k]` rather than zipped against `atoms`'s
    // own iteration order: `atoms` is the *whole* combined Z-matrix here, not
    // a small standalone Vec of just the fragment, so `positions[k]` is
    // almost never `k` itself. Zipping silently pairs the wrong rows together
    // instead of failing loudly -- it updated early host atoms while leaving
    // the actual fragment rows untouched, which looked identical to this
    // function never having run at all.
    for &index in positions {
        let Some(atom) = atoms.get(index) else {
            continue;
        };
        let (Some(r1), Some(r2), Some(r3)) = (atom.bond_ref, atom.angle_ref, atom.dihedral_ref)
        else {
            continue;
        };
        let (Some(&p), Some(&a), Some(&b), Some(&c)) = (
            coords.get(index),
            coords.get(r1 - 1),
            coords.get(r2 - 1),
            coords.get(r3 - 1),
        ) else {
            continue;
        };

        let bond_len = p.distance(a) as f64;
        let (angle, dihedral) = internal_from_direction(a, b, c, p - a);
        let atom = &mut atoms[index];
        atom.bond_len = bond_len;
        atom.angle_deg = angle;
        atom.dihedral_deg = dihedral;
    }
}

/// Exact geometric direction that completes a regular trigonal-planar (two
/// existing neighbours near 120 degrees) or tetrahedral (three existing
/// neighbours near 109.47 degrees) arrangement around `center`.
///
/// Both cases use the same identity: unit vectors spread evenly around a
/// centre sum to zero (three coplanar vectors at 120 degrees; four vectors to
/// a regular tetrahedron's corners), so the one missing direction is exactly
/// the negated sum of the others. This is what keeps a newly attached sp2
/// substituent exactly coplanar with the other two without a separate
/// "detect the plane, snap the dihedral" step: the sum of two coplanar
/// vectors is itself coplanar.
///
/// Returns `None` when the neighbour count or measured angles do not clearly
/// match one of these two patterns. Forcing a rule from an ambiguous count
/// (zero or one neighbours) or an irregular angle would misplace the
/// fragment as confidently as this fixes the regular cases -- callers should
/// fall back to a clearance-driven search in that case.
pub fn ideal_direction_from_neighbors(center: Vec3, neighbors: &[Vec3]) -> Option<Vec3> {
    const SP2_ANGLE_DEG: f32 = 120.0;
    const SP3_ANGLE_DEG: f32 = 109.471;
    const TOLERANCE_DEG: f32 = 15.0;

    let directions: Vec<Vec3> = neighbors
        .iter()
        .map(|&p| (p - center).normalize_or_zero())
        .filter(|d| d.length() > 0.5)
        .collect();

    let angle_between = |a: Vec3, b: Vec3| a.dot(b).clamp(-1.0, 1.0).acos().to_degrees();

    match directions.len() {
        2 => {
            let angle = angle_between(directions[0], directions[1]);
            if (angle - SP2_ANGLE_DEG).abs() > TOLERANCE_DEG {
                return None;
            }
            let sum = directions[0] + directions[1];
            (sum.length() > 1.0e-3).then(|| -sum.normalize())
        }
        3 => {
            for i in 0..3 {
                for j in (i + 1)..3 {
                    if (angle_between(directions[i], directions[j]) - SP3_ANGLE_DEG).abs()
                        > TOLERANCE_DEG
                    {
                        return None;
                    }
                }
            }
            let sum = directions[0] + directions[1] + directions[2];
            (sum.length() > 1.0e-3).then(|| -sum.normalize())
        }
        _ => None,
    }
}

/// The direction completing a tetrahedral centre that already has two bonds.
///
/// Negating the sum of the two existing bond vectors -- the right answer for
/// three bonds, and for a trigonal centre with two -- is wrong here: it lands
/// **in** the plane of the two bonds, at 125.26 degrees from each, not the
/// 109.47 a tetrahedron needs. The two remaining vertices are out of that
/// plane, symmetric about it; this returns one of them.
///
/// Derivation: with `s` the in-plane bisector and `n` the plane normal, the
/// remaining vertices are `-s cos(54.7356) +/- n sin(54.7356)`, whose
/// components are `sqrt(1/3)` and `sqrt(2/3)`.
pub fn tetrahedral_completion(center: Vec3, first: Vec3, second: Vec3) -> Option<Vec3> {
    /// cos(54.7356 deg) = sqrt(1/3): half the tetrahedral angle's complement.
    const IN_PLANE: f32 = 0.577_350_3;
    /// sin(54.7356 deg) = sqrt(2/3).
    const OUT_OF_PLANE: f32 = 0.816_496_6;

    let a = (first - center).normalize_or_zero();
    let b = (second - center).normalize_or_zero();
    if a.length() < 0.5 || b.length() < 0.5 {
        return None;
    }
    let bisector = (a + b).normalize_or_zero();
    let normal = a.cross(b).normalize_or_zero();
    // Collinear bonds leave the normal undefined, so there is no plane to sit
    // out of and no tetrahedron to complete.
    if bisector.length() < 0.5 || normal.length() < 0.5 {
        return None;
    }
    Some((-bisector * IN_PLANE + normal * OUT_OF_PLANE).normalize_or_zero())
}

/// Direction that completes a linear (sp) arrangement, when the single
/// existing neighbour is short enough to be a genuine double or triple bond.
///
/// A lone existing neighbour is geometrically ambiguous on its own: a
/// divalent atom joined only by ordinary single bonds (an ether oxygen, a
/// thiol sulfur, a carbon three-quarters of the way through being built up)
/// is not linear, and this application does not track bond order explicitly.
/// Gating on bond length is what turns "probably sp" into "definitely sp":
/// once a multiple bond is confirmed, a carbon in the middle of it can only
/// ever have two sigma bonds, at 180 degrees, by definition -- so completing
/// it linearly is the only valid geometry, not a guess.
pub fn linear_direction_for_multiple_bond(
    center: Vec3,
    neighbor: Vec3,
    single_bond_length: f32,
) -> Option<Vec3> {
    // Generous enough to catch a real triple bond (ratio ~0.80 for C-C) while
    // safely excluding an ordinary single bond (ratio ~1.0).
    const TRIPLE_BOND_RATIO: f32 = 0.87;

    let direction = (neighbor - center).normalize_or_zero();
    if direction.length() < 0.5 || single_bond_length <= 0.0 {
        return None;
    }
    let observed = center.distance(neighbor);
    (observed / single_bond_length <= TRIPLE_BOND_RATIO).then_some(-direction)
}

/// Distance factor that counts a pair as bonded, matching
/// `Molecule::recompute_bonds`.
const BOND_THRESH_SCALE: f32 = 1.2;

/// Real bonded neighbours of `atom` in the host structure, found by
/// covalent-radius distance -- the same test `Molecule::recompute_bonds` uses
/// to draw the bonds the user actually sees.
///
/// A Z-matrix `bond_ref` chain cannot serve as the bond graph here: for a
/// structure that came in through `xyz_to_zmat`, every reference is assigned
/// from file order (row `i` points at row `i - 1`), so the references
/// routinely name atoms that are not bonded at all. Perceiving the bonds from
/// the coordinates keeps the valence rules correct for loaded structures as
/// well as for ones built here.
///
/// `exclude` is typically the vacancy being filled and the fragment's own
/// (not yet meaningfully placed) atoms.
pub fn host_bonded_neighbors(
    zmat: &[ZAtom],
    coords: &[Vec3],
    atom: usize,
    exclude: &[usize],
) -> Vec<usize> {
    let symbols: Vec<&str> = zmat.iter().map(|a| a.symbol.as_str()).collect();
    bonded_indices(&symbols, coords, atom, exclude)
}

/// The same perception straight from a molecule's own symbols and
/// coordinates, which is all this ever needed: it is the geometry that
/// decides, not the Z-matrix.
pub fn bonded_indices(
    symbols: &[&str],
    coords: &[Vec3],
    atom: usize,
    exclude: &[usize],
) -> Vec<usize> {
    let Some(&center) = coords.get(atom) else {
        return Vec::new();
    };
    let Some(center_symbol) = symbols.get(atom) else {
        return Vec::new();
    };
    let mut found: Vec<usize> = (0..coords.len().min(symbols.len()))
        .filter(|&j| j != atom && !exclude.contains(&j))
        .filter(|&j| {
            let cutoff = (covalent_radius_angstrom(center_symbol)
                + covalent_radius_angstrom(symbols[j]))
                * BOND_THRESH_SCALE;
            let d = center.distance(coords[j]);
            d > 0.01 && d <= cutoff
        })
        .collect();
    // Nearest first, so a caller wanting "the atom this one hangs off" can
    // just take the first.
    found.sort_by(|&a, &b| {
        center
            .distance(coords[a])
            .total_cmp(&center.distance(coords[b]))
    });
    found
}

/// Real bonded neighbours of `center` among `others`, by covalent-radius
/// distance -- the same test `Molecule::recompute_bonds` uses. A fragment's
/// own local Cartesian shape carries no Z-matrix bond references to consult,
/// so its connectivity has to be measured geometrically instead.
pub fn geometric_bonded_neighbors(
    center: Vec3,
    center_symbol: &str,
    others: &[(Vec3, &str)],
) -> Vec<Vec3> {
    others
        .iter()
        .filter_map(|&(p, symbol)| {
            let cutoff =
                (covalent_radius_angstrom(center_symbol) + covalent_radius_angstrom(symbol))
                    * BOND_THRESH_SCALE;
            let d = center.distance(p);
            (d > 0.01 && d <= cutoff).then_some(p)
        })
        .collect()
}


#[cfg(test)]
mod tetrahedral_tests {
    use super::*;

    /// Every direction on the cone makes the requested angle with the axis,
    /// and they are spread around it rather than bunched.
    #[test]
    fn cone_directions_all_make_the_requested_angle() {
        let axis = Vec3::new(0.3, -0.5, 0.81).normalize();
        let dirs = cone_directions(axis, 104.5, 12);
        assert_eq!(dirs.len(), 12);
        for d in &dirs {
            let angle = axis.dot(*d).clamp(-1.0, 1.0).acos().to_degrees();
            assert!((angle - 104.5).abs() < 0.01, "got {angle}");
            assert!((d.length() - 1.0).abs() < 1.0e-5);
        }
        // Spread: the first and the opposite one differ.
        assert!(dirs[0].distance(dirs[6]) > 1.0, "the cone is not degenerate");
    }

    #[test]
    fn a_degenerate_cone_request_yields_nothing() {
        assert!(cone_directions(Vec3::ZERO, 104.5, 8).is_empty());
        assert!(cone_directions(Vec3::X, 104.5, 0).is_empty());
    }

    /// A 180-degree cone collapses onto the reversed axis, which is what an
    /// sp centre wants.
    #[test]
    fn a_straight_cone_is_the_reversed_axis() {
        for d in cone_directions(Vec3::Z, 180.0, 5) {
            assert!((d - -Vec3::Z).length() < 1.0e-5, "{d:?}");
        }
    }

    /// Two bonds at the tetrahedral angle: the completion must sit at 109.47
    /// from both, and out of their plane.
    #[test]
    fn completing_a_tetrahedral_centre_gives_the_tetrahedral_angle() {
        let center = Vec3::ZERO;
        let half = 109.471_f32.to_radians() * 0.5;
        let a = Vec3::new(half.sin(), 0.0, half.cos());
        let b = Vec3::new(-half.sin(), 0.0, half.cos());
        let c = tetrahedral_completion(center, a, b).expect("a completion");

        let angle = |u: Vec3, v: Vec3| u.normalize().dot(v.normalize()).acos().to_degrees();
        assert!((angle(a, c) - 109.471).abs() < 0.01, "{}", angle(a, c));
        assert!((angle(b, c) - 109.471).abs() < 0.01, "{}", angle(b, c));
        assert!(c.y.abs() > 0.5, "it must leave the a-b plane, got y = {}", c.y);
        assert!((c.length() - 1.0).abs() < 1.0e-5);
    }

    /// What the naive negated sum would have given, for contrast: in-plane at
    /// 125.26 degrees. This is the bug the function exists to avoid.
    #[test]
    fn the_negated_sum_would_have_been_planar_and_too_wide() {
        let half = 109.471_f32.to_radians() * 0.5;
        let a = Vec3::new(half.sin(), 0.0, half.cos());
        let b = Vec3::new(-half.sin(), 0.0, half.cos());
        let naive = -(a + b).normalize();
        let angle = a.normalize().dot(naive).acos().to_degrees();
        assert!((angle - 125.26).abs() < 0.1, "{angle}");
        assert!(naive.y.abs() < 1.0e-5, "and it stays in the plane");
    }

    /// Collinear bonds have no plane to sit out of.
    #[test]
    fn collinear_bonds_have_no_tetrahedral_completion() {
        let center = Vec3::ZERO;
        assert!(tetrahedral_completion(center, Vec3::X, -Vec3::X).is_none());
        assert!(tetrahedral_completion(center, Vec3::X, Vec3::X * 2.0).is_none());
    }

    #[test]
    fn a_bond_of_zero_length_is_rejected() {
        assert!(tetrahedral_completion(Vec3::ZERO, Vec3::ZERO, Vec3::X).is_none());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(
        symbol: &str,
        bond: Option<usize>,
        angle: Option<usize>,
        dihedral: Option<usize>,
    ) -> ZAtom {
        ZAtom {
            symbol: symbol.to_string(),
            bond_ref: bond,
            bond_len: 1.5,
            angle_ref: angle,
            angle_deg: 109.471,
            dihedral_ref: dihedral,
            dihedral_deg: 180.0,
        }
    }

    /// A four-atom chain: C1, C2-C1, C3-C2-C1, C4-C3-C2-C1.
    fn chain() -> Vec<ZAtom> {
        vec![
            atom("C", None, None, None),
            atom("C", Some(1), None, None),
            atom("C", Some(2), Some(1), None),
            atom("C", Some(3), Some(2), Some(1)),
        ]
    }

    #[test]
    fn the_chain_is_inherited_from_the_attachment_atom() {
        let host = chain();
        // Attaching at atom 4 (index 3): its own bond_ref is 3 and angle_ref 2,
        // so those become the new atom's angle and dihedral references.
        assert_eq!(host_reference_chain(&host, 3), (Some(2), Some(1)));
        // Attaching at atom 3 (index 2): bond_ref 2, angle_ref 1.
        assert_eq!(host_reference_chain(&host, 2), (Some(1), Some(0)));
    }

    #[test]
    fn the_first_two_rows_use_moldens_special_cases() {
        let host = chain();
        assert_eq!(host_reference_chain(&host, 0), (Some(1), Some(2)));
        assert_eq!(host_reference_chain(&host, 1), (Some(0), Some(2)));
    }

    #[test]
    fn a_reference_never_points_at_the_attachment_atom_itself() {
        let host = chain();
        for attach in 0..host.len() {
            let (angle, dihedral) = host_reference_chain(&host, attach);
            assert_ne!(angle, Some(attach), "angle ref aliased attach {attach}");
            assert_ne!(
                dihedral,
                Some(attach),
                "dihedral ref aliased attach {attach}"
            );
        }
    }

    /// The defect: a spliced fragment must come out with every reference filled.
    /// Leaving one `None` sends `zmat_to_xyz` to atom 1 of the whole molecule.
    #[test]
    fn splicing_leaves_no_unfilled_references() {
        let host = chain();
        // A three-atom fragment as it comes out of `build_fragment_zmat`:
        // incomplete by construction.
        let mut fragment = vec![
            atom("C", None, None, None),
            atom("H", Some(1), None, None),
            atom("H", Some(1), Some(2), None),
        ];
        let positions = vec![4, 5, 6]; // appended after the four host atoms
        assign_fragment_references(&mut fragment, &positions, &host, 3, (None, None));

        for (i, a) in fragment.iter().enumerate() {
            assert!(a.bond_ref.is_some(), "fragment atom {i} has no bond ref");
            assert!(a.angle_ref.is_some(), "fragment atom {i} has no angle ref");
            assert!(
                a.dihedral_ref.is_some(),
                "fragment atom {i} has no dihedral ref"
            );
        }

        // The connector hangs off the attachment atom, and inherits its chain.
        assert_eq!(fragment[0].bond_ref, Some(4));
        assert_eq!(fragment[0].angle_ref, Some(3));
        assert_eq!(fragment[0].dihedral_ref, Some(2));
        // The second fragment atom is bonded to the connector.
        assert_eq!(fragment[1].bond_ref, Some(5));
        assert_eq!(fragment[1].angle_ref, Some(4));
    }

    #[test]
    fn explicit_connector_references_win_over_the_inherited_chain() {
        let host = chain();
        let mut fragment = vec![atom("O", None, None, None)];
        assign_fragment_references(&mut fragment, &[4], &host, 3, (Some(0), Some(1)));
        assert_eq!(fragment[0].angle_ref, Some(1));
        assert_eq!(fragment[0].dihedral_ref, Some(2));
    }

    #[test]
    fn fragment_internal_references_are_translated_not_offset_blindly() {
        let host = chain();
        let mut fragment = vec![
            atom("C", None, None, None),
            atom("C", Some(1), None, None),
            atom("H", Some(2), Some(1), None),
        ];
        // Connector replaces host atom 3; the rest append.
        let positions = vec![2, 4, 5];
        assign_fragment_references(&mut fragment, &positions, &host, 1, (None, None));
        // Fragment atom 3 referenced fragment atoms 2 and 1, which live at
        // combined indices 4 and 2.
        assert_eq!(fragment[2].bond_ref, Some(5));
        assert_eq!(fragment[2].angle_ref, Some(3));
    }

    #[test]
    fn partners_cover_the_attachment_atom_and_its_frame() {
        let host = chain();
        // Atom 4 references atoms 3 and 2, so all three are partners.
        let mut partners = attachment_partners(&host, 3);
        partners.sort_unstable();
        assert_eq!(partners, vec![1, 2, 3]);
        // Atom 1 has no references of its own.
        assert_eq!(attachment_partners(&host, 0), vec![0]);
    }

    #[test]
    fn clearance_ignores_the_attachment_partners() {
        let coords = vec![
            Vec3::ZERO,
            Vec3::new(1.5, 0.0, 0.0),
            Vec3::new(3.0, 0.0, 0.0),
        ];
        // Atom 2 is the group; atom 1 is its bonded partner and excluded.
        let with_ignore = min_clearance(&coords, &[2], &[1]);
        assert!((with_ignore - 3.0).abs() < 1e-6, "got {with_ignore}");
        let without = min_clearance(&coords, &[2], &[]);
        assert!((without - 1.5).abs() < 1e-6, "got {without}");
    }

    /// The conversion must invert exactly, or every derived coordinate is
    /// subtly wrong in a way that looks like a geometry bug.
    #[test]
    fn direction_and_internal_coordinates_round_trip() {
        let a = Vec3::new(0.4, -0.2, 1.1);
        let b = Vec3::new(-1.1, 0.3, 0.2);
        let c = Vec3::new(-0.7, 1.6, -0.9);

        for angle in [30.0, 70.0, 109.471, 150.0] {
            for dihedral in [0.0, 45.0, 120.0, 180.0, 275.0] {
                let direction = direction_from_internal(a, b, c, angle, dihedral);
                let (back_angle, back_dihedral) = internal_from_direction(a, b, c, direction);
                assert!(
                    (back_angle - angle).abs() < 1.0e-3,
                    "angle {angle} came back as {back_angle}"
                );
                let delta = (back_dihedral - dihedral).rem_euclid(360.0);
                let delta = delta.min(360.0 - delta);
                assert!(
                    delta < 1.0e-3,
                    "dihedral {dihedral} came back as {back_dihedral}"
                );
            }
        }
    }

    /// A placed fragment must keep its own shape: rigid motion only.
    #[test]
    fn placement_preserves_internal_distances() {
        let local = vec![
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, 1.07),
            Vec3::new(1.008, 0.0, -0.357),
            Vec3::new(-0.504, -0.874, -0.357),
        ];
        let open = fragment_open_valence(&local, "C", &["H", "H", "H"]);
        let placed = place_fragment(&local, open, Vec3::new(2.0, 1.0, -1.0), Vec3::X, 37.0);

        for i in 0..local.len() {
            for j in (i + 1)..local.len() {
                let before = local[i].distance(local[j]);
                let after = placed[i].distance(placed[j]);
                assert!(
                    (before - after).abs() < 1.0e-4,
                    "distance {i}-{j} changed {before} -> {after}"
                );
            }
        }
    }

    /// The fragment body must end up along the requested direction, with its
    /// open valence facing back at the host rather than into the molecule.
    #[test]
    fn placement_points_the_body_along_the_bond() {
        let local = vec![
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, 1.07),
            Vec3::new(1.008, 0.0, -0.357),
            Vec3::new(-0.504, -0.874, -0.357),
        ];
        let host = Vec3::ZERO;
        let direction = Vec3::new(1.0, 1.0, 0.0).normalize();
        let origin = host + direction * 1.54;
        let open = fragment_open_valence(&local, "C", &["H", "H", "H"]);
        let placed = place_fragment(&local, open, origin, direction, 0.0);

        // Every substituent must lie on the far side of the connector from the
        // host, which is what "not folded back" means quantitatively.
        let to_host = (host - placed[0]).normalize();
        for (i, &p) in placed.iter().enumerate().skip(1) {
            let to_atom = (p - placed[0]).normalize();
            let angle = to_host.dot(to_atom).clamp(-1.0, 1.0).acos().to_degrees();
            assert!(
                angle > 90.0,
                "fragment atom {i} sits {angle:.1} deg from the host bond"
            );
        }
    }

    #[test]
    fn sp2_direction_completes_a_120_degree_arrangement() {
        let center = Vec3::ZERO;
        // Two neighbours at exactly 120 degrees, in the xy-plane.
        let n1 = Vec3::new(1.0, 0.0, 0.0);
        let n2 = Vec3::new((-120.0f32).to_radians().cos(), (-120.0f32).to_radians().sin(), 0.0);
        let direction = ideal_direction_from_neighbors(center, &[n1, n2])
            .expect("two 120-degree neighbours should be recognised as sp2");

        for &n in &[n1, n2] {
            let angle = direction
                .normalize()
                .dot(n.normalize())
                .clamp(-1.0, 1.0)
                .acos()
                .to_degrees();
            assert!((angle - 120.0).abs() < 1e-2, "angle to neighbour was {angle}");
        }
        // Exactly coplanar: the new direction has no out-of-plane (z) component.
        assert!(direction.z.abs() < 1e-5, "direction left the plane: {direction:?}");
    }

    #[test]
    fn sp3_direction_completes_a_tetrahedral_arrangement() {
        // Three of a regular tetrahedron's four vertex directions.
        let a = Vec3::new(1.0, 1.0, 1.0).normalize();
        let b = Vec3::new(1.0, -1.0, -1.0).normalize();
        let c = Vec3::new(-1.0, 1.0, -1.0).normalize();
        let expected_fourth = Vec3::new(-1.0, -1.0, 1.0).normalize();

        let direction = ideal_direction_from_neighbors(Vec3::ZERO, &[a, b, c])
            .expect("three 109.47-degree neighbours should be recognised as sp3");

        assert!(
            direction.dot(expected_fourth) > 0.999,
            "direction {direction:?} does not match the fourth tetrahedral vertex"
        );
    }

    #[test]
    fn irregular_angles_are_not_forced_into_a_rule() {
        // Neighbours at 90 degrees match neither sp2 nor sp3.
        let n1 = Vec3::X;
        let n2 = Vec3::Y;
        assert!(ideal_direction_from_neighbors(Vec3::ZERO, &[n1, n2]).is_none());
    }

    #[test]
    fn one_or_four_neighbours_are_left_to_the_caller() {
        assert!(ideal_direction_from_neighbors(Vec3::ZERO, &[Vec3::X]).is_none());
        assert!(ideal_direction_from_neighbors(Vec3::ZERO, &[]).is_none());
        let tet = [
            Vec3::new(1.0, 1.0, 1.0).normalize(),
            Vec3::new(1.0, -1.0, -1.0).normalize(),
            Vec3::new(-1.0, 1.0, -1.0).normalize(),
            Vec3::new(-1.0, -1.0, 1.0).normalize(),
        ];
        assert!(ideal_direction_from_neighbors(Vec3::ZERO, &tet).is_none());
    }

    #[test]
    fn a_triple_bond_length_neighbour_forces_linearity() {
        // A C#C-length bond (1.19 A) against a C-C single-bond sum (~1.50 A).
        let center = Vec3::ZERO;
        let neighbor = Vec3::new(1.19, 0.0, 0.0);
        let direction = linear_direction_for_multiple_bond(center, neighbor, 1.50)
            .expect("a triple-bond-length neighbour should force linearity");
        assert!(direction.dot(Vec3::NEG_X) > 0.999, "expected the opposite direction");
    }

    #[test]
    fn an_ordinary_single_bond_does_not_force_linearity() {
        let center = Vec3::ZERO;
        let neighbor = Vec3::new(1.50, 0.0, 0.0); // exactly the single-bond length
        assert!(linear_direction_for_multiple_bond(center, neighbor, 1.50).is_none());
    }

    #[test]
    fn a_double_bond_length_neighbour_does_not_force_linearity() {
        // A C=C-length bond (1.34 A): sp2, not sp -- must NOT be treated as linear.
        let center = Vec3::ZERO;
        let neighbor = Vec3::new(1.34, 0.0, 0.0);
        assert!(linear_direction_for_multiple_bond(center, neighbor, 1.50).is_none());
    }

    #[test]
    fn host_bonded_neighbors_perceives_bonds_from_the_coordinates() {
        // A bent three-atom chain whose Z-matrix references deliberately
        // disagree with the real bonding: row 2 chains to row 1 the way
        // `xyz_to_zmat` would, even though it sits bonded to row 0.
        let zmat = vec![
            atom("C", None, None, None),
            atom("C", Some(1), None, None),
            atom("H", Some(2), Some(1), None),
        ];
        let coords = vec![
            Vec3::ZERO,
            Vec3::new(1.335, 0.0, 0.0),
            Vec3::new(-0.544, 0.943, 0.0),
        ];
        let mut neighbors = host_bonded_neighbors(&zmat, &coords, 0, &[]);
        neighbors.sort_unstable();
        assert_eq!(neighbors, vec![1, 2]);
    }

    #[test]
    fn host_bonded_neighbors_excludes_requested_indices() {
        let zmat = vec![
            atom("C", None, None, None),
            atom("C", Some(1), None, None),
            atom("H", Some(2), Some(1), None),
        ];
        let coords = vec![
            Vec3::ZERO,
            Vec3::new(1.335, 0.0, 0.0),
            Vec3::new(-0.544, 0.943, 0.0),
        ];
        assert_eq!(host_bonded_neighbors(&zmat, &coords, 0, &[1]), vec![2]);
    }

    /// A row anchored to the atom that is about to move must be re-pointed at
    /// something stable, and must not be left in a frame whose dihedral is
    /// undefined -- `zmat_to_xyz` resolves those with an arbitrary
    /// perpendicular, so the atom would land somewhere else.
 
    /// When the only atoms available would give a collinear frame, the row is
    /// left exactly as it was: a frame whose dihedral is undefined would put
    /// the atom somewhere else entirely once `zmat_to_xyz` picked its own
    /// arbitrary perpendicular.
 
    /// Rows that never referenced the moving atom are left completely alone.
 
    #[test]
    fn geometric_bonded_neighbors_matches_covalent_radius_distance() {
        let center = Vec3::ZERO;
        let bonded = Vec3::new(1.35, 0.0, 0.0); // within C=C range
        let far = Vec3::new(5.0, 0.0, 0.0);
        let found = geometric_bonded_neighbors(center, "C", &[(bonded, "C"), (far, "C")]);
        assert_eq!(found, vec![bonded]);
    }
}
