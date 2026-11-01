//! Vibronic analysis against closed-form answers.
//!
//! The quantities here are simple enough to be known exactly for a constructed
//! system, so the tests check numbers rather than plausibility. Two identities
//! carry most of the weight:
//!
//! * the reorganisation energy of a mode is `g_k^2 / (2 w_k^2)`, derived
//!   independently of the `S_k w_k` route the code takes;
//! * the Duschinsky matrix of a state against itself is the identity, and stays
//!   the identity when that state is rigidly rotated and translated -- which is
//!   the only thing that can catch a broken frame alignment.

use ndarray::Array2;

use super::{duschinsky, franck_condon_activity, AMU_TO_ELECTRON_MASS, CM1_TO_AU};

/// Three atoms with deliberately unequal masses, so a mass-weighting error
/// cannot cancel.
const MASSES_AMU: [f64; 3] = [15.9949, 1.00783, 2.01410];

/// A mass-metric-orthonormal mode set, `sum_i m_i L_ik^2 = 1`.
///
/// Built by Gram-Schmidt in the mass metric from a fixed, arbitrary starting
/// set, so it is deterministic and has no accidental symmetry that could hide
/// an index error.
fn orthonormal_modes(count: usize) -> Array2<f64> {
    let dimensions = 3 * MASSES_AMU.len();
    let masses: Vec<f64> = (0..dimensions)
        .map(|i| MASSES_AMU[i / 3] * AMU_TO_ELECTRON_MASS)
        .collect();
    let mut columns: Vec<Vec<f64>> = Vec::new();
    for k in 0..count {
        // An arbitrary but reproducible seed vector.
        let mut vector: Vec<f64> = (0..dimensions)
            .map(|i| ((i * 7 + k * 13) % 11) as f64 - 5.0)
            .collect();
        for existing in &columns {
            let overlap: f64 = (0..dimensions)
                .map(|i| masses[i] * existing[i] * vector[i])
                .sum();
            for i in 0..dimensions {
                vector[i] -= overlap * existing[i];
            }
        }
        let norm: f64 = (0..dimensions)
            .map(|i| masses[i] * vector[i] * vector[i])
            .sum::<f64>()
            .sqrt();
        assert!(norm > 1.0e-8, "the seed vectors are degenerate");
        for value in &mut vector {
            *value /= norm;
        }
        columns.push(vector);
    }
    let mut modes = Array2::<f64>::zeros((dimensions, count));
    for (k, column) in columns.iter().enumerate() {
        for (i, value) in column.iter().enumerate() {
            modes[(i, k)] = *value;
        }
    }
    modes
}

/// Every Franck-Condon quantity, against closed forms.
///
/// The reorganisation energy is the sharpest of these because it is checked
/// against `g_k^2 / (2 w_k^2)`, a route that shares no arithmetic with the
/// `S_k w_k` the code computes: an error in the exponent of the displacement
/// would have to conspire with an error in the Huang-Rhys factor to survive
/// both.
#[test]
fn franck_condon_quantities_match_their_closed_forms() {
    let modes = orthonormal_modes(4);
    let frequencies_cm1 = [1600.0, 1200.0, 800.0, 400.0];
    let gradient: Vec<f64> = (0..9).map(|i| 0.01 * ((i % 5) as f64 - 2.0)).collect();
    let which: Vec<usize> = (0..4).collect();

    let analysis = franck_condon_activity(&modes, &frequencies_cm1, &gradient, &MASSES_AMU, &which)
        .expect("the analysis must run");
    assert_eq!(analysis.modes.len(), 4);

    for mode in &analysis.modes {
        let omega_au = mode.frequency_cm1 * CM1_TO_AU;
        let g = mode.gradient_au;

        // The projection itself.
        let expected_g: f64 = (0..9).map(|i| modes[(i, mode.index)] * gradient[i]).sum();
        assert!((g - expected_g).abs() < 1.0e-14);

        // Displacement: the minimum of `1/2 w^2 Q^2 + g Q` is at `-g/w^2`,
        // and `sqrt(w)` makes it dimensionless.
        let expected_displacement = -g / omega_au.powi(2) * omega_au.sqrt();
        assert!(
            (mode.displacement - expected_displacement).abs()
                < 1.0e-12 * expected_displacement.abs().max(1.0),
            "displacement {} vs {expected_displacement}",
            mode.displacement
        );

        // Huang-Rhys.
        assert!((mode.huang_rhys - 0.5 * mode.displacement.powi(2)).abs() < 1.0e-15);

        // Reorganisation energy, by the independent route.
        let expected_lambda_au = g * g / (2.0 * omega_au * omega_au);
        let reported_lambda_au = mode.reorganisation_cm1 * CM1_TO_AU;
        assert!(
            (reported_lambda_au - expected_lambda_au).abs()
                < 1.0e-10 * expected_lambda_au.abs().max(1.0e-12),
            "mode {}: reorganisation {reported_lambda_au:.6e} vs closed form \
             {expected_lambda_au:.6e}",
            mode.index
        );

        // Resonance-Raman intensity is `w^2 delta^2`, which is also `g^2/w`.
        let expected_intensity = g * g / omega_au;
        assert!(
            (mode.resonance_raman - expected_intensity).abs()
                < 1.0e-10 * expected_intensity.abs().max(1.0e-12)
        );
    }

    // The total is the sum.
    let summed: f64 = analysis.modes.iter().map(|m| m.reorganisation_cm1).sum();
    assert!((analysis.reorganisation_cm1 - summed).abs() < 1.0e-12);
}

/// A complete mode set accounts for the whole gradient; an orthogonal one
/// accounts for none of it.
///
/// `unaccounted_fraction` is the diagnostic that catches a gradient and a
/// Hessian computed at different geometries, so both of its limits are worth
/// pinning. The metric matters: were it computed against the Cartesian norm
/// rather than the mass-weighted one, the first case would report a large
/// unaccounted fraction for a perfectly complete set.
#[test]
fn the_unaccounted_fraction_reaches_both_limits() {
    let dimensions = 9;
    let complete = orthonormal_modes(dimensions);
    let frequencies = vec![1000.0; dimensions];
    let gradient: Vec<f64> = (0..dimensions).map(|i| 0.02 * (i as f64 + 1.0)).collect();
    let all: Vec<usize> = (0..dimensions).collect();

    let full = franck_condon_activity(&complete, &frequencies, &gradient, &MASSES_AMU, &all)
        .expect("the complete analysis must run");
    assert!(
        full.unaccounted_fraction < 1.0e-7,
        "a complete mode set left {:.3e} of the gradient unaccounted for",
        full.unaccounted_fraction
    );

    // A gradient built from one mode, analysed with the others: nothing is
    // recovered, so the whole gradient is unaccounted for.
    let masses_au: Vec<f64> = (0..dimensions)
        .map(|i| MASSES_AMU[i / 3] * AMU_TO_ELECTRON_MASS)
        .collect();
    let along_first: Vec<f64> = (0..dimensions)
        .map(|i| masses_au[i] * complete[(i, 0)])
        .collect();
    let others: Vec<usize> = (1..dimensions).collect();
    let orthogonal =
        franck_condon_activity(&complete, &frequencies, &along_first, &MASSES_AMU, &others)
            .expect("the orthogonal analysis must run");
    assert!(
        orthogonal.unaccounted_fraction > 1.0 - 1.0e-7,
        "a gradient orthogonal to every analysed mode reported {:.6} accounted for",
        1.0 - orthogonal.unaccounted_fraction
    );
    for mode in &orthogonal.modes {
        assert!(
            mode.gradient_au.abs() < 1.0e-12,
            "mode {} picked up {} from an orthogonal gradient",
            mode.index,
            mode.gradient_au
        );
    }
}

/// An imaginary mode is skipped rather than given an invented progression.
#[test]
fn imaginary_modes_are_skipped() {
    let modes = orthonormal_modes(3);
    let frequencies = [-500.0, 1000.0, 1500.0];
    let gradient: Vec<f64> = (0..9).map(|i| 0.01 * (i as f64)).collect();
    let analysis = franck_condon_activity(&modes, &frequencies, &gradient, &MASSES_AMU, &[0, 1, 2])
        .expect("the analysis must run");
    assert_eq!(
        analysis.modes.len(),
        2,
        "the imaginary mode must be dropped"
    );
    assert!(analysis.modes.iter().all(|mode| mode.index != 0));
}

/// A mass list of the wrong length is refused rather than silently misindexed.
///
/// Passing `3N` masses instead of `N` would read the first `N` entries and
/// quietly corrupt the diagnostic, which is the kind of error that survives
/// review.
#[test]
fn a_wrong_length_mass_list_is_refused() {
    let modes = orthonormal_modes(3);
    let frequencies = [1000.0, 1200.0, 1400.0];
    let gradient = vec![0.01; 9];
    let per_coordinate: Vec<f64> = (0..9).map(|i| MASSES_AMU[i / 3]).collect();
    let error =
        franck_condon_activity(&modes, &frequencies, &gradient, &per_coordinate, &[0, 1, 2])
            .expect_err("a per-coordinate mass list must be refused");
    assert!(format!("{error}").contains("one mass per atom"));
}

/// Modes orthogonal to translation and rotation: genuine vibrations.
///
/// The Gram-Schmidt set above is fine wherever only the mode algebra matters,
/// but any test that *displaces the geometry* needs real vibrations. A vector
/// with translational or rotational content moves the molecular frame, the
/// mass-weighted superposition then removes that motion -- correctly -- and
/// the displacement no longer matches what a naive calculation predicts.
fn vibrational_modes() -> Array2<f64> {
    let atoms = MASSES_AMU.len();
    let dimensions = 3 * atoms;
    let masses: Vec<f64> = (0..dimensions)
        .map(|i| MASSES_AMU[i / 3] * AMU_TO_ELECTRON_MASS)
        .collect();
    let coords = geometry();
    let total: f64 = MASSES_AMU.iter().map(|m| m * AMU_TO_ELECTRON_MASS).sum();
    let centre: [f64; 3] = std::array::from_fn(|axis| {
        (0..atoms)
            .map(|atom| masses[3 * atom] * coords[3 * atom + axis])
            .sum::<f64>()
            / total
    });

    let dot =
        |a: &[f64], b: &[f64]| -> f64 { (0..dimensions).map(|i| masses[i] * a[i] * b[i]).sum() };
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let push = |candidate: Vec<f64>, basis: &mut Vec<Vec<f64>>| {
        let mut vector = candidate;
        for existing in basis.iter() {
            let overlap = dot(existing, &vector);
            for i in 0..dimensions {
                vector[i] -= overlap * existing[i];
            }
        }
        let norm = dot(&vector, &vector).sqrt();
        if norm > 1.0e-6 {
            for value in &mut vector {
                *value /= norm;
            }
            basis.push(vector);
            true
        } else {
            false
        }
    };

    // Three translations, then three rotations about the mass centre.
    for axis in 0..3 {
        let translation: Vec<f64> = (0..dimensions)
            .map(|i| if i % 3 == axis { 1.0 } else { 0.0 })
            .collect();
        push(translation, &mut basis);
    }
    for axis in 0..3 {
        let rotation: Vec<f64> = (0..atoms)
            .flat_map(|atom| {
                let r: [f64; 3] = std::array::from_fn(|k| coords[3 * atom + k] - centre[k]);
                let mut unit = [0.0; 3];
                unit[axis] = 1.0;
                // e_axis x r
                [
                    unit[1] * r[2] - unit[2] * r[1],
                    unit[2] * r[0] - unit[0] * r[2],
                    unit[0] * r[1] - unit[1] * r[0],
                ]
                .into_iter()
            })
            .collect();
        push(rotation, &mut basis);
    }
    let rigid = basis.len();

    // Whatever is left is vibrational.
    let mut vibrations: Vec<Vec<f64>> = Vec::new();
    for seed in 0..dimensions {
        let candidate: Vec<f64> = (0..dimensions)
            .map(|i| ((i * 5 + seed * 3) % 7) as f64 - 3.0)
            .collect();
        let mut combined = basis.clone();
        combined.extend(vibrations.iter().cloned());
        let before = combined.len();
        if push(candidate, &mut combined) && combined.len() > before {
            vibrations.push(combined.pop().expect("just pushed"));
        }
    }
    assert_eq!(
        vibrations.len(),
        dimensions - rigid,
        "a non-linear three-atom molecule has three vibrations"
    );

    let mut modes = Array2::<f64>::zeros((dimensions, vibrations.len()));
    for (k, column) in vibrations.iter().enumerate() {
        for (i, value) in column.iter().enumerate() {
            modes[(i, k)] = *value;
        }
    }
    modes
}

/// A geometry, for the Duschinsky tests.
fn geometry() -> Vec<f64> {
    vec![
        0.0, 0.0, 0.2216, //
        0.0, 1.4309, -0.8865, //
        0.0, -1.4309, -0.8865,
    ]
}

/// A state against itself gives the identity and no displacement.
#[test]
fn duschinsky_of_a_state_against_itself_is_the_identity() {
    let modes = orthonormal_modes(3);
    let frequencies = [1600.0, 3700.0, 3800.0];
    let coords = geometry();
    let which: Vec<usize> = (0..3).collect();

    let result = duschinsky(
        &modes,
        &frequencies,
        &coords,
        &modes,
        &frequencies,
        &coords,
        &MASSES_AMU,
        &which,
        &which,
    )
    .expect("the analysis must run");

    for row in 0..3 {
        for column in 0..3 {
            let expected = if row == column { 1.0 } else { 0.0 };
            assert!(
                (result.rotation[(row, column)] - expected).abs() < 1.0e-10,
                "J[{row}][{column}] = {}, expected {expected}",
                result.rotation[(row, column)]
            );
        }
    }
    for (index, value) in result.displacement.iter().enumerate() {
        assert!(
            value.abs() < 1.0e-10,
            "mode {index} reported a displacement of {value} between a state and itself"
        );
    }
    assert!(result.total_huang_rhys() < 1.0e-18);
    assert!((result.mean_participation() - 1.0).abs() < 1.0e-10);
}

/// A rigid rotation and translation of one state leaves the Duschinsky matrix
/// alone.
///
/// This is the test that justifies the frame alignment, and the only one that
/// can fail if it is dropped: without it the rotation between the two frames
/// leaks into `J` and reads as mode mixing that is not there. The rotation
/// below is a deliberately awkward one about all three axes.
#[test]
fn a_rigid_motion_does_not_leak_into_the_duschinsky_matrix() {
    let modes = orthonormal_modes(3);
    let frequencies = [1600.0, 3700.0, 3800.0];
    let coords = geometry();
    let which: Vec<usize> = (0..3).collect();

    // A rotation by 0.7 rad about z then 0.4 rad about x, and a translation.
    let (c1, s1) = (0.7_f64.cos(), 0.7_f64.sin());
    let (c2, s2) = (0.4_f64.cos(), 0.4_f64.sin());
    let rotate = |v: [f64; 3]| -> [f64; 3] {
        let after_z = [c1 * v[0] - s1 * v[1], s1 * v[0] + c1 * v[1], v[2]];
        [
            after_z[0],
            c2 * after_z[1] - s2 * after_z[2],
            s2 * after_z[1] + c2 * after_z[2],
        ]
    };
    let translation = [1.3_f64, -0.8, 2.1];

    let moved_coords: Vec<f64> = (0..3)
        .flat_map(|atom| {
            let point: [f64; 3] = std::array::from_fn(|axis| coords[3 * atom + axis]);
            let turned = rotate(point);
            (0..3).map(move |axis| turned[axis] + translation[axis])
        })
        .collect();
    let mut moved_modes = Array2::<f64>::zeros(modes.dim());
    for column in 0..modes.dim().1 {
        for atom in 0..3 {
            let vector: [f64; 3] = std::array::from_fn(|axis| modes[(3 * atom + axis, column)]);
            let turned = rotate(vector);
            for axis in 0..3 {
                moved_modes[(3 * atom + axis, column)] = turned[axis];
            }
        }
    }

    let result = duschinsky(
        &modes,
        &frequencies,
        &coords,
        &moved_modes,
        &frequencies,
        &moved_coords,
        &MASSES_AMU,
        &which,
        &which,
    )
    .expect("the analysis must run");

    for row in 0..3 {
        for column in 0..3 {
            let expected = if row == column { 1.0 } else { 0.0 };
            assert!(
                (result.rotation[(row, column)] - expected).abs() < 1.0e-8,
                "a rigid motion leaked into J[{row}][{column}] = {}",
                result.rotation[(row, column)]
            );
        }
    }
    for value in &result.displacement {
        assert!(
            value.abs() < 1.0e-8,
            "a rigid motion produced a displacement of {value}"
        );
    }
}

/// The Duschinsky matrix of a complete mode set is orthogonal.
///
/// `J = L_A^T M L_B` with both sets orthonormal in the mass metric, so
/// `J^T J = I` whenever the two sets span the same space. Checking it costs
/// nothing and catches a dropped mass factor, which would leave `J` finite but
/// not orthogonal.
#[test]
fn the_duschinsky_matrix_is_orthogonal_for_complete_mode_sets() {
    let dimensions = 9;
    let set_a = orthonormal_modes(dimensions);
    // A different complete set: reverse the Gram-Schmidt order by seeding
    // differently, so the two span the same space but are not aligned.
    let mut set_b = Array2::<f64>::zeros(set_a.dim());
    for column in 0..dimensions {
        for row in 0..dimensions {
            set_b[(row, column)] = set_a[(row, dimensions - 1 - column)];
        }
    }
    let frequencies = vec![1000.0; dimensions];
    let coords = geometry();
    let which: Vec<usize> = (0..dimensions).collect();

    let result = duschinsky(
        &set_a,
        &frequencies,
        &coords,
        &set_b,
        &frequencies,
        &coords,
        &MASSES_AMU,
        &which,
        &which,
    )
    .expect("the analysis must run");

    for i in 0..dimensions {
        for j in 0..dimensions {
            let product: f64 = (0..dimensions)
                .map(|k| result.rotation[(k, i)] * result.rotation[(k, j)])
                .sum();
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (product - expected).abs() < 1.0e-9,
                "(J^T J)[{i}][{j}] = {product}, expected {expected}"
            );
        }
    }
}

/// A displaced minimum shows up in `K`, with the magnitude the harmonic model
/// predicts.
///
/// Displacing the geometry along one mass-weighted mode by `d` should give
/// `K_k = d sqrt(w_k)` for that mode and nothing for the others. This is the
/// only test of the displacement vector's scale, and `K` is the piece the
/// original implementation did not compute at all.
#[test]
fn the_displacement_vector_reports_a_known_shift() {
    // Genuine vibrations: displacing along a vector with rigid-body content
    // would move the frame, which the superposition then removes.
    let modes = vibrational_modes();
    let frequencies = [1600.0, 3700.0, 3800.0];
    let coords = geometry();
    let which: Vec<usize> = (0..3).collect();

    // Move along mode 1 by a mass-weighted amount `d`, i.e. in Cartesian
    // coordinates by `d * L_i1`.
    let amount = 0.05_f64;
    let displaced: Vec<f64> = (0..9)
        .map(|index| coords[index] + amount * modes[(index, 1)])
        .collect();

    let result = duschinsky(
        &modes,
        &frequencies,
        &coords,
        &modes,
        &frequencies,
        &displaced,
        &MASSES_AMU,
        &which,
        &which,
    )
    .expect("the analysis must run");

    for (index, value) in result.displacement.iter().enumerate() {
        let omega_au = frequencies[index] * CM1_TO_AU;
        let expected = if index == 1 {
            amount * omega_au.sqrt()
        } else {
            0.0
        };
        assert!(
            (value - expected).abs() < 1.0e-6 * expected.abs().max(1.0e-3),
            "K[{index}] = {value}, expected {expected}"
        );
    }
}
