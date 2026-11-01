//! Proper rigid alignment, with explicit rank-zero/rank-one handling for atoms and linear molecules.
use crate::rmsd::kabsch::{centroid, Superposition};
use bevy::prelude::{Quat, Vec3};
fn linear_axis(points: &[Vec3], center: Vec3) -> Option<Vec3> {
    let longest = points
        .iter()
        .map(|p| *p - center)
        .max_by(|a, b| a.length_squared().total_cmp(&b.length_squared()))?;
    let axis = longest.try_normalize()?;
    points
        .iter()
        .all(|p| {
            (*p - center).cross(axis).length_squared() <= 1e-12 * longest.length_squared().max(1.)
        })
        .then_some(axis)
}
fn line_rotation(line: &[Vec3], other: &[Vec3], a: Vec3, b: Vec3, axis: Vec3) -> Quat {
    let target = line
        .iter()
        .zip(other)
        .map(|(p, q)| (*q - b) * (*p - a).dot(axis))
        .sum::<Vec3>();
    target
        .try_normalize()
        .map(|target| Quat::from_rotation_arc(axis, target))
        .unwrap_or(Quat::IDENTITY)
}
pub fn superpose(mobile: &[Vec3], reference: &[Vec3]) -> Result<Superposition, String> {
    if mobile.is_empty()
        || mobile.len() != reference.len()
        || mobile.iter().chain(reference).any(|p| !p.is_finite())
    {
        return Err(
            "Alignment needs nonempty structures with matching atom counts and finite coordinates."
                .into(),
        );
    }
    let a = centroid(mobile);
    let b = centroid(reference);
    let rotation = if mobile.iter().all(|p| p.distance_squared(a) < 1e-20)
        || reference.iter().all(|p| p.distance_squared(b) < 1e-20)
    {
        Quat::IDENTITY
    } else if let Some(axis) = linear_axis(mobile, a) {
        line_rotation(mobile, reference, a, b, axis)
    } else if let Some(axis) = linear_axis(reference, b) {
        line_rotation(reference, mobile, b, a, axis).conjugate()
    } else {
        return crate::rmsd::kabsch::superpose(mobile, reference);
    };
    let x = rotation * Vec3::X;
    let y = rotation * Vec3::Y;
    let z = rotation * Vec3::Z;
    Ok(Superposition {
        rotation: [
            [x.x as f64, y.x as f64, z.x as f64],
            [x.y as f64, y.y as f64, z.y as f64],
            [x.z as f64, y.z as f64, z.z as f64],
        ],
        mobile_centroid: a,
        reference_centroid: b,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rotated_diatomic_and_single_atom_have_zero_aligned_rmsd() {
        for (a, b) in [
            (
                vec![Vec3::ZERO, Vec3::X],
                vec![Vec3::new(2., 2., 0.), Vec3::new(2., 3., 0.)],
            ),
            (vec![Vec3::ONE], vec![Vec3::X]),
        ] {
            let s = superpose(&a, &b).unwrap();
            assert!(crate::rmsd::kabsch::rmsd(&s.apply_all(&a), &b).unwrap() < 1e-6);
        }
    }
    #[test]
    fn linear_alignment_is_a_proper_rotation() {
        let a = [Vec3::ZERO, Vec3::X, Vec3::X * 2.];
        let b = [Vec3::Z * 2., Vec3::Z, Vec3::ZERO];
        let s = superpose(&a, &b).unwrap();
        assert!((crate::numerics::mat3::determinant(&s.rotation) - 1.).abs() < 1e-5);
        assert!(crate::rmsd::kabsch::rmsd(&s.apply_all(&a), &b).unwrap() < 1e-6);
    }
}
