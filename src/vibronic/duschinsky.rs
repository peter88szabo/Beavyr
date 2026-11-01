//! How one electronic state's normal modes are built out of another's.
//!
//! # The two objects
//!
//! A vibronic transition between two states with different geometries *and*
//! different force constants is described by
//!
//! ```text
//! Q_A = J Q_B + d (mass-weighted normal coordinates)
//! K_k = sqrt(omega_Ak) d_k (reported dimensionless displacement)
//! ```
//!
//! `J` is the Duschinsky rotation, saying how much of each mode of one state
//! appears in each mode of the other. `K` is the displacement between the two
//! minima, expressed in normal coordinates. Both are needed for a vibronic
//! band shape; `J` alone answers the narrower and very common question of
//! which mode descends from which.
//!
//! ```text
//! J = L_A^T M L_B
//! K_k = sqrt(omega_Ak) sum_i m_i L^A_ik (x^B_i - x^A_i)
//! ```
//!
//! The mass matrix appears explicitly because the stored columns are
//! Cartesian and normalised against the mass metric, not mass-weighted. `K` is
//! reported in dimensionless normal coordinates, which is what a
//! Franck-Condon progression is written in.
//!
//! # Alignment, and why it is mass-weighted
//!
//! The two geometries are superposed first. A normal mode is defined in the
//! molecule's own frame, and two optimisations have no reason to leave their
//! structures in the same one; without alignment the rotation between frames
//! leaks into `J` and looks like mode mixing.
//!
//! The superposition is **mass-weighted**, which is what the Eckart
//! conditions actually require. An unweighted fit answers a different
//! question -- it treats a hydrogen and a metal alike -- and for a
//! hydride-rich molecule the two rotations differ enough to matter.
//!
//! Only the rotation is applied to a mode. A mode is a displacement, so
//! translating it would be meaningless; the translation belongs to `K`, where
//! it is removed by referring both geometries to their mass centres.

use anyhow::{ensure, Result};
use ndarray::Array2;

use super::AMU_TO_ELECTRON_MASS;

/// One state's modes expressed in another's.
#[derive(Debug, Clone)]
pub struct Duschinsky {
    /// `J[i][j]`: the amount of state-A mode `i` in state-B mode `j`. Rows are
    /// A, columns are B.
    ///
    /// Note this is the transpose of the convention some references use, so
    /// check the index order before comparing against a published matrix.
    pub rotation: Array2<f64>,
    /// `K`, the displacement between the minima along each state-A mode, in
    /// dimensionless normal coordinates.
    pub displacement: Vec<f64>,
    /// State-A mode indices labelling the rows of `rotation` and the entries
    /// of `displacement`.
    pub rows: Vec<usize>,
    /// State-B mode indices labelling the columns.
    pub columns: Vec<usize>,
    pub row_frequencies_cm1: Vec<f64>,
    pub column_frequencies_cm1: Vec<f64>,
}

/// One state-B mode and the state-A mode it most resembles.
#[derive(Debug, Clone, PartialEq)]
pub struct ModeCorrespondence {
    pub column_index: usize,
    pub column_frequency_cm1: f64,
    pub row_index: usize,
    pub row_frequency_cm1: f64,
    /// `J_ij^2` for the best partner, normalised over the block: the fraction
    /// of this mode that is the partner. One is a pure correspondence.
    pub overlap: f64,
    /// Inverse participation ratio over `J^2` down the column: roughly how
    /// many state-A modes this one is built from. One means no mixing.
    pub participation: f64,
    pub shift_cm1: f64,
}

impl Duschinsky {
    /// For each state-B mode, the state-A mode it most resembles.
    pub fn correspondence(&self) -> Vec<ModeCorrespondence> {
        let (rows, columns) = self.rotation.dim();
        let mut out = Vec::with_capacity(columns);
        for column in 0..columns {
            let mut best = 0usize;
            let mut best_square = -1.0_f64;
            let mut sum_square = 0.0_f64;
            let mut sum_fourth = 0.0_f64;
            for row in 0..rows {
                let square = self.rotation[(row, column)] * self.rotation[(row, column)];
                sum_square += square;
                sum_fourth += square * square;
                if square > best_square {
                    best_square = square;
                    best = row;
                }
            }
            // Normalised over the block, so a column spread across a subset of
            // modes still reports a meaningful fraction.
            let overlap = if sum_square > 0.0 {
                best_square / sum_square
            } else {
                0.0
            };
            let participation = if sum_fourth > 0.0 {
                (sum_square * sum_square) / sum_fourth
            } else {
                0.0
            };
            out.push(ModeCorrespondence {
                column_index: self.columns[column],
                column_frequency_cm1: self.column_frequencies_cm1[column],
                row_index: self.rows[best],
                row_frequency_cm1: self.row_frequencies_cm1[best],
                overlap,
                participation,
                shift_cm1: self.column_frequencies_cm1[column] - self.row_frequencies_cm1[best],
            });
        }
        out
    }

    /// Mean participation over all columns: one number for how much the second
    /// state scrambles the first state's modes.
    ///
    /// One is no mixing at all. Larger means the second state's vibrations are
    /// combinations of several of the first's, which is the situation that
    /// opens intramolecular energy-redistribution channels.
    pub fn mean_participation(&self) -> f64 {
        let correspondence = self.correspondence();
        if correspondence.is_empty() {
            return 0.0;
        }
        correspondence
            .iter()
            .map(|mode| mode.participation)
            .sum::<f64>()
            / correspondence.len() as f64
    }

    /// The total Huang-Rhys factor implied by `K` alone, `sum_k K_k^2 / 2`.
    ///
    /// Comparable with the sum from [`super::franck_condon_activity`] when the
    /// two states are related by a small displacement, and a useful check that
    /// the two routes to a Franck-Condon picture agree.
    pub fn total_huang_rhys(&self) -> f64 {
        0.5 * self.displacement.iter().map(|k| k * k).sum::<f64>()
    }
}

/// The Duschinsky rotation and displacement between two mode sets.
///
/// `coords_a` and `coords_b` are flat Cartesian coordinates in bohr, three per
/// atom. `which_a` and `which_b` select the genuine vibrations of each state.
#[allow(clippy::too_many_arguments)]
pub fn duschinsky(
    modes_a: &Array2<f64>,
    frequencies_a_cm1: &[f64],
    coords_a: &[f64],
    modes_b: &Array2<f64>,
    frequencies_b_cm1: &[f64],
    coords_b: &[f64],
    masses_amu: &[f64],
    which_a: &[usize],
    which_b: &[usize],
) -> Result<Duschinsky> {
    let atoms = masses_amu.len();
    let dimensions = 3 * atoms;
    ensure!(
        modes_a.dim().0 == dimensions && modes_b.dim().0 == dimensions,
        "the two mode sets describe different numbers of atoms than the mass list"
    );
    ensure!(
        coords_a.len() == dimensions && coords_b.len() == dimensions,
        "coordinates do not match the mass list"
    );
    ensure!(
        frequencies_a_cm1.len() == modes_a.dim().1 && frequencies_b_cm1.len() == modes_b.dim().1,
        "each mode set needs one frequency per column"
    );

    ensure!(
        !which_a.is_empty() && !which_b.is_empty(),
        "Select at least one mode in each state"
    );
    ensure!(
        which_a.iter().all(|&i| i < modes_a.ncols())
            && which_b.iter().all(|&i| i < modes_b.ncols()),
        "Mode index is out of range"
    );
    ensure!(
        masses_amu.iter().all(|m| m.is_finite() && *m > 0.0),
        "Masses must be finite and positive"
    );
    ensure!(
        coords_a
            .iter()
            .chain(coords_b)
            .chain(modes_a.iter())
            .chain(modes_b.iter())
            .chain(frequencies_a_cm1)
            .chain(frequencies_b_cm1)
            .all(|v| v.is_finite()),
        "Non-finite geometry or normal modes"
    );
    let masses_au: Vec<f64> = masses_amu
        .iter()
        .map(|mass| mass * AMU_TO_ELECTRON_MASS)
        .collect();

    // Mass-weighted superposition of B onto A, and the same rotation applied
    // to B's modes so they stay in step with their geometry.
    let (rotation, shifted_b, centred_a) =
        mass_weighted_superposition(coords_b, coords_a, &masses_au)?;

    let rotate = |vector: [f64; 3]| -> [f64; 3] {
        std::array::from_fn(|row| {
            (0..3)
                .map(|column| rotation[row][column] * vector[column])
                .sum()
        })
    };

    let mut matrix = Array2::<f64>::zeros((which_a.len(), which_b.len()));
    for (row, &a) in which_a.iter().enumerate() {
        for (column, &b) in which_b.iter().enumerate() {
            let mut accumulated = 0.0;
            for atom in 0..atoms {
                let mass = masses_au[atom];
                let mode_a: [f64; 3] = std::array::from_fn(|axis| modes_a[(3 * atom + axis, a)]);
                let mode_b = rotate(std::array::from_fn(|axis| modes_b[(3 * atom + axis, b)]));
                accumulated += mass * (0..3).map(|axis| mode_a[axis] * mode_b[axis]).sum::<f64>();
            }
            matrix[(row, column)] = accumulated;
        }
    }

    // `K_k = sum_i m_i L^A_ik (x^B_i - x^A_i)`, with both geometries referred
    // to their mass centres and B rotated onto A, then scaled to dimensionless
    // normal coordinates by `sqrt(w_k)`.
    let displacement = which_a
        .iter()
        .map(|&a| {
            let mut value = 0.0;
            for atom in 0..atoms {
                let mass = masses_au[atom];
                for axis in 0..3 {
                    let index = 3 * atom + axis;
                    value += mass * modes_a[(index, a)] * (shifted_b[index] - centred_a[index]);
                }
            }
            let frequency_au = frequencies_a_cm1[a] * super::CM1_TO_AU;
            if frequency_au > 0.0 {
                value * frequency_au.sqrt()
            } else {
                0.0
            }
        })
        .collect();

    Ok(Duschinsky {
        rotation: matrix,
        displacement,
        rows: which_a.to_vec(),
        columns: which_b.to_vec(),
        row_frequencies_cm1: which_a.iter().map(|&i| frequencies_a_cm1[i]).collect(),
        column_frequencies_cm1: which_b.iter().map(|&i| frequencies_b_cm1[i]).collect(),
    })
}

/// Mass-weighted Kabsch superposition of `mobile` onto `reference`.
///
/// Returns the rotation, the rotated and centred mobile coordinates, and the
/// centred reference coordinates. Both are referred to their **mass centres**,
/// which is what removes the translation the Eckart conditions also remove.
type Superposition = ([[f64; 3]; 3], Vec<f64>, Vec<f64>);

pub(super) fn mass_weighted_superposition(
    mobile: &[f64],
    reference: &[f64],
    masses_au: &[f64],
) -> Result<Superposition> {
    let atoms = masses_au.len();
    let total_mass: f64 = masses_au.iter().sum();
    ensure!(total_mass > 0.0, "the total mass must be positive");

    let centre = |coords: &[f64]| -> [f64; 3] {
        std::array::from_fn(|axis| {
            (0..atoms)
                .map(|atom| masses_au[atom] * coords[3 * atom + axis])
                .sum::<f64>()
                / total_mass
        })
    };
    let shift = |coords: &[f64], origin: [f64; 3]| -> Vec<f64> {
        (0..3 * atoms)
            .map(|index| coords[index] - origin[index % 3])
            .collect()
    };
    let centred_mobile = shift(mobile, centre(mobile));
    let centred_reference = shift(reference, centre(reference));

    // The mass-weighted correlation matrix, followed by a quaternion fit.
    let mut correlation = [[0.0f64; 3]; 3];
    for atom in 0..atoms {
        let mass = masses_au[atom];
        for row in 0..3 {
            for column in 0..3 {
                correlation[row][column] +=
                    mass * centred_mobile[3 * atom + row] * centred_reference[3 * atom + column];
            }
        }
    }
    let rotation = rotation_from_correlation(&correlation)?;

    let rotated: Vec<f64> = (0..atoms)
        .flat_map(|atom| {
            let point: [f64; 3] = std::array::from_fn(|axis| centred_mobile[3 * atom + axis]);
            let turned: [f64; 3] = std::array::from_fn(|row| {
                (0..3)
                    .map(|column| rotation[row][column] * point[column])
                    .sum()
            });
            turned.into_iter()
        })
        .collect();
    Ok((rotation, rotated, centred_reference))
}

/// The proper rotation that best maps the mobile frame onto the reference,
/// from the correlation matrix, by the quaternion method.
///
/// Builds the symmetric `4x4` key matrix whose largest eigenvector is the
/// optimal quaternion. That route is used rather than an SVD because it
/// cannot produce an improper rotation: a reflection is simply not
/// representable as a quaternion, so the degenerate case that forces an SVD
/// implementation to flip the sign of its smallest singular value does not
/// arise.
fn rotation_from_correlation(correlation: &[[f64; 3]; 3]) -> Result<[[f64; 3]; 3]> {
    let c = correlation;
    let key = [
        [
            c[0][0] + c[1][1] + c[2][2],
            c[1][2] - c[2][1],
            c[2][0] - c[0][2],
            c[0][1] - c[1][0],
        ],
        [
            c[1][2] - c[2][1],
            c[0][0] - c[1][1] - c[2][2],
            c[0][1] + c[1][0],
            c[2][0] + c[0][2],
        ],
        [
            c[2][0] - c[0][2],
            c[0][1] + c[1][0],
            -c[0][0] + c[1][1] - c[2][2],
            c[1][2] + c[2][1],
        ],
        [
            c[0][1] - c[1][0],
            c[2][0] + c[0][2],
            c[1][2] + c[2][1],
            -c[0][0] - c[1][1] + c[2][2],
        ],
    ];

    // Solve all four eigenvectors: a fixed power-iteration seed can be exactly
    // orthogonal to the largest root for a 180-degree rotation.
    use crate::normalmode::linalg_shim::{Backend, LinAlg};
    let matrix = Array2::from_shape_fn((4, 4), |(i, j)| key[i][j]);
    let (values, vectors) = LinAlg::new(Backend::Auto)?.eigh(&matrix)?;
    let largest = (0..4)
        .max_by(|&i, &j| values[i].total_cmp(&values[j]))
        .unwrap();
    let vector: [f64; 4] = std::array::from_fn(|i| vectors[(i, largest)]);

    let [w, x, y, z] = vector;
    Ok([
        [
            w * w + x * x - y * y - z * z,
            2.0 * (x * y - w * z),
            2.0 * (x * z + w * y),
        ],
        [
            2.0 * (x * y + w * z),
            w * w - x * x + y * y - z * z,
            2.0 * (y * z - w * x),
        ],
        [
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            w * w - x * x - y * y + z * z,
        ],
    ])
}
