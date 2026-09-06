//! Minimal fixed-size 3x3 real matrix arithmetic.
//!
//! Small enough to be exact and dependency-free: everything here is written
//! for the 3x3 case specifically, which is all the Kabsch superposition and
//! the inertia-tensor style problems ever need. Nothing in Beavyr links a
//! BLAS, so these are plain loops rather than calls into one.

/// A 3x3 matrix in row-major order: `m[row][col]`.
pub type Mat3 = [[f64; 3]; 3];

pub const IDENTITY: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

pub fn zeros() -> Mat3 {
    [[0.0; 3]; 3]
}

pub fn transpose(a: &Mat3) -> Mat3 {
    let mut out = zeros();
    for i in 0..3 {
        for j in 0..3 {
            out[j][i] = a[i][j];
        }
    }
    out
}

pub fn matmul(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut out = zeros();
    for i in 0..3 {
        for j in 0..3 {
            let mut sum = 0.0;
            for k in 0..3 {
                sum += a[i][k] * b[k][j];
            }
            out[i][j] = sum;
        }
    }
    out
}

/// `a * v`, treating `v` as a column vector.
pub fn mat_vec(a: &Mat3, v: [f64; 3]) -> [f64; 3] {
    let mut out = [0.0; 3];
    for i in 0..3 {
        out[i] = a[i][0] * v[0] + a[i][1] * v[1] + a[i][2] * v[2];
    }
    out
}

/// Closed-form determinant -- exact for 3x3, no pivoting needed.
pub fn determinant(a: &Mat3) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_is_the_multiplicative_identity() {
        let a: Mat3 = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0]];
        assert_eq!(matmul(&a, &IDENTITY), a);
        assert_eq!(matmul(&IDENTITY, &a), a);
    }

    #[test]
    fn transpose_swaps_indices_and_is_its_own_inverse() {
        let a: Mat3 = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let t = transpose(&a);
        assert_eq!(t[0][1], a[1][0]);
        assert_eq!(t[2][0], a[0][2]);
        assert_eq!(transpose(&t), a);
    }

    #[test]
    fn matmul_matches_a_hand_computed_product() {
        let a: Mat3 = [[1.0, 2.0, 0.0], [0.0, 1.0, 3.0], [1.0, 0.0, 1.0]];
        let b: Mat3 = [[2.0, 0.0, 1.0], [1.0, 3.0, 0.0], [0.0, 1.0, 4.0]];
        let c = matmul(&a, &b);
        assert_eq!(c[0], [4.0, 6.0, 1.0]);
        assert_eq!(c[1], [1.0, 6.0, 12.0]);
        assert_eq!(c[2], [2.0, 1.0, 5.0]);
    }

    #[test]
    fn determinant_of_the_identity_is_one_and_of_a_singular_matrix_is_zero() {
        assert!((determinant(&IDENTITY) - 1.0).abs() < 1e-12);
        // Two identical rows -> singular.
        let singular: Mat3 = [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        assert!(determinant(&singular).abs() < 1e-12);
    }

    #[test]
    fn determinant_flips_sign_for_a_reflection() {
        let reflection: Mat3 = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!((determinant(&reflection) + 1.0).abs() < 1e-12);
    }

    #[test]
    fn mat_vec_rotates_a_basis_vector() {
        // 90 degrees about z: x -> y.
        let rot: Mat3 = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let v = mat_vec(&rot, [1.0, 0.0, 0.0]);
        assert!((v[0] - 0.0).abs() < 1e-12);
        assert!((v[1] - 1.0).abs() < 1e-12);
    }
}
