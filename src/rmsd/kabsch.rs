//! Root-mean-square deviation and optimal rigid-body superposition.
//!
//! The Kabsch algorithm finds the rotation minimising the RMSD between two
//! equally-sized, index-matched point sets: centre both on their centroids,
//! form the 3x3 cross-covariance `H = P^T Q`, take `H = U S V^T`, and the
//! optimal rotation is `R = V * diag(1, 1, d) * U^T` where
//! `d = sign(det(V U^T))`. That `d` is not cosmetic: without it the
//! decomposition can return the best *reflection* rather than the best
//! rotation, which would superimpose a molecule onto its own mirror image
//! and report a misleadingly small RMSD for two genuinely different
//! enantiomers.
//!
//! Only a 3x3 SVD is ever needed, so `crate::numerics` supplies it directly
//! and no BLAS/LAPACK is involved.

use bevy::prelude::Vec3;

use crate::numerics::mat3::{determinant, mat_vec, matmul, transpose, zeros, Mat3};
use crate::numerics::svd3::svd3;

/// What the panel reports after comparing two structures.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RmsdReport {
    pub atom_count: usize,
    /// RMSD exactly as the two structures sit, with no fitting at all.
    pub raw_rmsd: f64,
    /// RMSD after optimal translation + rotation -- the smallest value any
    /// rigid-body placement can achieve for this atom pairing.
    pub aligned_rmsd: f64,
    /// The largest single-atom deviation after alignment, which localises a
    /// change an averaged number hides.
    pub max_atom_deviation: f64,
}

/// Plain RMSD between two index-matched structures, with no superposition.
pub fn rmsd(a: &[Vec3], b: &[Vec3]) -> Result<f64, String> {
    check_pair(a, b)?;
    let sum: f64 = a
        .iter()
        .zip(b)
        .map(|(p, q)| (*p - *q).length_squared() as f64)
        .sum();
    Ok((sum / a.len() as f64).sqrt())
}

/// The rigid-body transform that best maps `mobile` onto `reference`.
#[derive(Debug, Clone, Copy)]
pub struct Superposition {
    pub rotation: Mat3,
    pub mobile_centroid: Vec3,
    pub reference_centroid: Vec3,
}

impl Superposition {
    /// Applies the transform to a single point: shift to the mobile
    /// centroid, rotate, then shift out to the reference centroid.
    pub fn apply(&self, p: Vec3) -> Vec3 {
        let centred = p - self.mobile_centroid;
        let rotated = mat_vec(
            &self.rotation,
            [centred.x as f64, centred.y as f64, centred.z as f64],
        );
        self.reference_centroid
            + Vec3::new(rotated[0] as f32, rotated[1] as f32, rotated[2] as f32)
    }

    pub fn apply_all(&self, points: &[Vec3]) -> Vec<Vec3> {
        points.iter().map(|&p| self.apply(p)).collect()
    }
}

pub fn centroid(points: &[Vec3]) -> Vec3 {
    if points.is_empty() {
        return Vec3::ZERO;
    }
    points.iter().copied().sum::<Vec3>() / points.len() as f32
}

/// Kabsch superposition of `mobile` onto `reference`.
pub fn superpose(mobile: &[Vec3], reference: &[Vec3]) -> Result<Superposition, String> {
    check_pair(mobile, reference)?;

    let mobile_centroid = centroid(mobile);
    let reference_centroid = centroid(reference);

    // Cross-covariance of the centred coordinates.
    let mut h: Mat3 = zeros();
    for (p, q) in mobile.iter().zip(reference) {
        let pc = *p - mobile_centroid;
        let qc = *q - reference_centroid;
        let pv = [pc.x as f64, pc.y as f64, pc.z as f64];
        let qv = [qc.x as f64, qc.y as f64, qc.z as f64];
        for i in 0..3 {
            for j in 0..3 {
                h[i][j] += pv[i] * qv[j];
            }
        }
    }

    let decomposed = svd3(&h);
    let vt_u = matmul(&decomposed.v, &transpose(&decomposed.u));
    // Guard against the decomposition handing back a reflection instead of
    // a rotation -- see the module note.
    let d = if determinant(&vt_u) < 0.0 { -1.0 } else { 1.0 };
    let mut correction = crate::numerics::mat3::IDENTITY;
    correction[2][2] = d;
    let rotation = matmul(&decomposed.v, &matmul(&correction, &transpose(&decomposed.u)));

    Ok(Superposition {
        rotation,
        mobile_centroid,
        reference_centroid,
    })
}

/// Compares `mobile` against `reference`, reporting both the unfitted and
/// the best-fit RMSD along with the worst single-atom deviation.
pub fn compare(mobile: &[Vec3], reference: &[Vec3]) -> Result<RmsdReport, String> {
    let raw_rmsd = rmsd(mobile, reference)?;
    let superposition = superpose(mobile, reference)?;
    let aligned = superposition.apply_all(mobile);
    let aligned_rmsd = rmsd(&aligned, reference)?;
    let max_atom_deviation = aligned
        .iter()
        .zip(reference)
        .map(|(p, q)| (*p - *q).length() as f64)
        .fold(0.0f64, f64::max);

    Ok(RmsdReport {
        atom_count: mobile.len(),
        raw_rmsd,
        aligned_rmsd,
        max_atom_deviation,
    })
}

fn check_pair(a: &[Vec3], b: &[Vec3]) -> Result<(), String> {
    if a.is_empty() || b.is_empty() {
        return Err("Cannot compare an empty structure.".to_string());
    }
    if a.len() != b.len() {
        return Err(format!(
            "Structures have different atom counts ({} vs {}); they cannot be \
             compared atom by atom.",
            a.len(),
            b.len()
        ));
    }
    Ok(())
}

/// Confirms two structures name the same elements in the same order, which
/// index-matched RMSD silently assumes.
pub fn same_atom_ordering(a: &[String], b: &[String]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x == y)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::FRAC_PI_2;

    fn sample() -> Vec<Vec3> {
        vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.5, 0.0, 0.0),
            Vec3::new(0.0, 1.2, 0.0),
            Vec3::new(0.3, 0.4, 1.1),
        ]
    }

    fn rotate_z(points: &[Vec3], angle: f32) -> Vec<Vec3> {
        let (s, c) = angle.sin_cos();
        points
            .iter()
            .map(|p| Vec3::new(c * p.x - s * p.y, s * p.x + c * p.y, p.z))
            .collect()
    }

    #[test]
    fn rmsd_of_a_structure_with_itself_is_zero() {
        let a = sample();
        assert!(rmsd(&a, &a).unwrap() < 1e-9);
    }

    #[test]
    fn rmsd_matches_a_hand_computed_value() {
        // Every atom displaced by exactly 0.5 A along x -> RMSD = 0.5.
        let a = sample();
        let b: Vec<Vec3> = a.iter().map(|p| *p + Vec3::new(0.5, 0.0, 0.0)).collect();
        assert!((rmsd(&a, &b).unwrap() - 0.5).abs() < 1e-6);
    }

    #[test]
    fn mismatched_atom_counts_are_rejected() {
        let a = sample();
        let b = &a[..3];
        assert!(rmsd(&a, b).is_err());
        assert!(superpose(&a, b).is_err());
    }

    #[test]
    fn empty_structures_are_rejected() {
        assert!(rmsd(&[], &[]).is_err());
    }

    /// The core property: a pure translation + rotation of a structure is
    /// still the *same* structure, so best-fit RMSD must come out at zero
    /// even though the raw RMSD does not.
    #[test]
    fn superposition_removes_a_pure_rotation_and_translation() {
        let reference = sample();
        let moved: Vec<Vec3> = rotate_z(&reference, FRAC_PI_2)
            .iter()
            .map(|p| *p + Vec3::new(3.0, -2.0, 1.0))
            .collect();

        let report = compare(&moved, &reference).unwrap();
        assert!(
            report.raw_rmsd > 1.0,
            "a displaced copy should look far apart before fitting: {}",
            report.raw_rmsd
        );
        assert!(
            report.aligned_rmsd < 1e-5,
            "after fitting it is the same structure: {}",
            report.aligned_rmsd
        );
        assert!(report.max_atom_deviation < 1e-5);
    }

    #[test]
    fn superposition_actually_moves_the_atoms_onto_the_reference() {
        let reference = sample();
        let moved: Vec<Vec3> = rotate_z(&reference, 0.7)
            .iter()
            .map(|p| *p + Vec3::new(-1.0, 4.0, 0.5))
            .collect();
        let fitted = superpose(&moved, &reference).unwrap().apply_all(&moved);
        for (f, r) in fitted.iter().zip(&reference) {
            assert!((*f - *r).length() < 1e-4, "{f:?} vs {r:?}");
        }
    }

    /// A genuine structural difference must survive fitting -- superposition
    /// removes placement, not chemistry.
    #[test]
    fn a_real_deformation_is_not_fitted_away() {
        let reference = sample();
        let mut deformed = reference.clone();
        deformed[3] += Vec3::new(0.0, 0.0, 0.8);
        let report = compare(&deformed, &reference).unwrap();
        assert!(
            report.aligned_rmsd > 0.05,
            "a moved atom must still register: {}",
            report.aligned_rmsd
        );
        assert!(report.max_atom_deviation > report.aligned_rmsd);
    }

    /// Without the determinant correction the fit can collapse a structure
    /// onto its mirror image; a chiral arrangement must not fit onto its
    /// own reflection.
    #[test]
    fn a_mirror_image_is_not_reported_as_identical() {
        let reference = sample();
        let mirrored: Vec<Vec3> = reference.iter().map(|p| Vec3::new(-p.x, p.y, p.z)).collect();
        let report = compare(&mirrored, &reference).unwrap();
        assert!(
            report.aligned_rmsd > 0.1,
            "a reflection is not a rotation: {}",
            report.aligned_rmsd
        );
    }

    #[test]
    fn the_rotation_is_a_proper_rotation() {
        let reference = sample();
        let moved = rotate_z(&reference, 1.1);
        let sup = superpose(&moved, &reference).unwrap();
        // det(R) == +1 for a rotation, -1 would be a reflection.
        assert!((determinant(&sup.rotation) - 1.0).abs() < 1e-6);
        // R^T R == I.
        let should_be_identity = matmul(&transpose(&sup.rotation), &sup.rotation);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((should_be_identity[i][j] - expected).abs() < 1e-6);
            }
        }
    }

    /// A planar molecule makes the cross-covariance rank-deficient, the case
    /// that breaks a naive `A v / s` SVD.
    #[test]
    fn a_planar_structure_still_superposes() {
        let reference = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.4, 0.0, 0.0),
            Vec3::new(0.7, 1.2, 0.0),
        ];
        let moved: Vec<Vec3> = rotate_z(&reference, 0.9)
            .iter()
            .map(|p| *p + Vec3::new(2.0, 1.0, 0.0))
            .collect();
        let report = compare(&moved, &reference).unwrap();
        assert!(report.aligned_rmsd < 1e-4, "{}", report.aligned_rmsd);
    }

    #[test]
    fn same_atom_ordering_detects_a_reordered_structure() {
        let a: Vec<String> = ["C", "H", "H"].iter().map(|s| s.to_string()).collect();
        let b: Vec<String> = ["H", "C", "H"].iter().map(|s| s.to_string()).collect();
        assert!(same_atom_ordering(&a, &a));
        assert!(!same_atom_ordering(&a, &b));
    }
}
