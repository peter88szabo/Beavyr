//! Evaluation of basis functions and molecular orbitals at a point in space.
//!
//! A shell is evaluated once into a scratch buffer holding its Cartesian
//! components in the generator's order, then either transformed to real
//! spherical functions or permuted, in both cases ending in Molden's file
//! order so the coefficients from `[MO]` line up index for index.

use super::basis::{cartesian_component_norm, Shell};
use super::cart2sph::{
    cartesian_powers, molden_cartesian_column, molden_spherical_row, n_cartesian, n_spherical,
    transform, MAX_ANGULAR_MOMENTUM,
};

/// Beyond this the Gaussian is below ~1e-8 of its peak and is skipped.
const EXPONENT_CUTOFF: f64 = 18.4;

/// Radius outside which a shell contributes nothing measurable.
///
/// Driven by the most diffuse primitive, since that is the slowest to decay.
pub fn shell_cutoff_radius(shell: &Shell) -> f64 {
    let smallest = shell
        .exponents
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min);
    if !smallest.is_finite() || smallest <= 0.0 {
        return f64::INFINITY;
    }
    (EXPONENT_CUTOFF / smallest).sqrt()
}

/// Largest number of Cartesian components in a supported shell (g: 15).
pub const MAX_SHELL_COMPONENTS: usize = 15;

/// Everything about a shell's angular part that depends only on `l`.
///
/// Built once per angular momentum and shared.  The ordering tables and the
/// per-component normalisation used to be recomputed for every shell at every
/// grid point, which a single orbital could absorb but a density quadratic
/// form cannot.
struct ShellLayout {
    /// Generator-order Cartesian powers paired with their relative norms.
    powers: Vec<([u8; 3], f64)>,
    /// One sparse row per Molden spherical slot: (Cartesian column, coefficient).
    spherical_rows: Vec<Vec<(usize, f64)>>,
    /// Molden Cartesian slot to generator column.
    cartesian_permutation: Vec<usize>,
}

fn shell_layout(l: usize) -> &'static ShellLayout {
    use std::sync::OnceLock;
    static LAYOUTS: OnceLock<Vec<ShellLayout>> = OnceLock::new();
    &LAYOUTS.get_or_init(|| {
        (0..=MAX_ANGULAR_MOMENTUM)
            .map(|l| {
                let powers = cartesian_powers(l)
                    .iter()
                    .map(|&p| (p, cartesian_component_norm(p)))
                    .collect();

                let n_cart = n_cartesian(l);
                let table = transform(l);
                let spherical_rows = (0..n_spherical(l))
                    .map(|slot| {
                        let row = molden_spherical_row(l, slot);
                        (0..n_cart)
                            .filter_map(|c| {
                                let coefficient = table[row * n_cart + c];
                                (coefficient != 0.0).then_some((c, coefficient))
                            })
                            .collect()
                    })
                    .collect();

                let cartesian_permutation = (0..n_cart)
                    .map(|slot| molden_cartesian_column(l, slot))
                    .collect();

                ShellLayout {
                    powers,
                    spherical_rows,
                    cartesian_permutation,
                }
            })
            .collect()
    })[l]
}

/// Scratch space reused across grid points so evaluation allocates nothing.
#[derive(Default)]
pub struct EvalScratch {
    cartesian: Vec<f64>,
}

/// Values of every basis function in one shell, in Molden file order.
///
/// `d` is the displacement from the shell centre to the sample point, in the
/// same length unit as the exponents (Bohr).  Writes `shell.n_functions()`
/// entries into `out` and returns `false` when the shell is too far away to
/// contribute, leaving `out` untouched.
pub fn shell_values(
    shell: &Shell,
    d: [f64; 3],
    scratch: &mut EvalScratch,
    out: &mut [f64],
) -> bool {
    let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];

    let mut radial = 0.0;
    for (&a, &c) in shell.exponents.iter().zip(&shell.coefficients) {
        let arg = a * r2;
        if arg < EXPONENT_CUTOFF {
            radial += c * (-arg).exp();
        }
    }
    if radial == 0.0 {
        return false;
    }

    let l = shell.l as usize;
    let layout = shell_layout(l);

    scratch.cartesian.clear();
    scratch
        .cartesian
        .extend(layout.powers.iter().map(|&(p, norm)| {
            let monomial = d[0].powi(p[0] as i32) * d[1].powi(p[1] as i32) * d[2].powi(p[2] as i32);
            radial * norm * monomial
        }));

    if shell.spherical {
        for (slot, row) in layout.spherical_rows.iter().enumerate() {
            out[slot] = row
                .iter()
                .map(|&(c, coefficient)| coefficient * scratch.cartesian[c])
                .sum();
        }
    } else {
        for (slot, &column) in layout.cartesian_permutation.iter().enumerate() {
            out[slot] = scratch.cartesian[column];
        }
    }
    true
}

/// Add one shell's contribution to a contracted value.
///
/// `weights` are the coefficients multiplying this shell's basis functions,
/// which for a molecular orbital is the matching slice of its coefficient
/// vector.
pub fn accumulate_shell(
    shell: &Shell,
    d: [f64; 3],
    weights: &[f64],
    scratch: &mut EvalScratch,
    out: &mut f64,
) {
    let mut values = [0.0f64; MAX_SHELL_COMPONENTS];
    let n = shell.n_functions();
    if !shell_values(shell, d, scratch, &mut values[..n]) {
        return;
    }
    for (&weight, &value) in weights.iter().zip(&values[..n]) {
        *out += weight * value;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn s_shell(alpha: f64) -> Shell {
        let mut s = Shell {
            center: 0,
            l: 0,
            spherical: true,
            exponents: vec![alpha],
            raw_coefficients: vec![1.0],
            coefficients: vec![],
        };
        s.normalize();
        s
    }

    /// A single normalised s primitive is (2a/pi)^(3/4) exp(-a r^2).
    #[test]
    fn s_function_matches_the_analytic_gaussian() {
        let alpha = 0.8;
        let shell = s_shell(alpha);
        let mut scratch = EvalScratch::default();

        for r in [0.0, 0.3, 1.0, 2.5] {
            let mut value = 0.0;
            accumulate_shell(&shell, [r, 0.0, 0.0], &[1.0], &mut scratch, &mut value);
            let expected = (2.0 * alpha / PI).powf(0.75) * (-alpha * r * r).exp();
            assert!(
                (value - expected).abs() < 1e-12,
                "r={r}: {value} != {expected}"
            );
        }
    }

    /// An s function integrates to 1 over all space: check on a fine grid.
    #[test]
    fn s_function_is_normalised_on_a_numerical_grid() {
        let alpha = 1.1;
        let shell = s_shell(alpha);
        let mut scratch = EvalScratch::default();

        let n = 90;
        let extent = 4.0;
        let h = 2.0 * extent / n as f64;
        let mut integral = 0.0;
        for ix in 0..n {
            let x = -extent + (ix as f64 + 0.5) * h;
            for iy in 0..n {
                let y = -extent + (iy as f64 + 0.5) * h;
                for iz in 0..n {
                    let z = -extent + (iz as f64 + 0.5) * h;
                    let mut v = 0.0;
                    accumulate_shell(&shell, [x, y, z], &[1.0], &mut scratch, &mut v);
                    integral += v * v;
                }
            }
        }
        integral *= h * h * h;
        assert!(
            (integral - 1.0).abs() < 1e-3,
            "integral of psi^2 was {integral}, expected 1"
        );
    }

    /// Every component of a d shell must be individually normalised, which is
    /// the convention the spherical transform assumes on its input.
    #[test]
    fn every_cartesian_d_component_is_normalised() {
        let alpha = 0.9;
        let mut shell = Shell {
            center: 0,
            l: 2,
            spherical: false,
            exponents: vec![alpha],
            raw_coefficients: vec![1.0],
            coefficients: vec![],
        };
        shell.normalize();
        let mut scratch = EvalScratch::default();

        let n = 70;
        let extent = 4.0;
        let h = 2.0 * extent / n as f64;

        for slot in 0..6 {
            let mut weights = vec![0.0; 6];
            weights[slot] = 1.0;
            let mut integral = 0.0;
            for ix in 0..n {
                let x = -extent + (ix as f64 + 0.5) * h;
                for iy in 0..n {
                    let y = -extent + (iy as f64 + 0.5) * h;
                    for iz in 0..n {
                        let z = -extent + (iz as f64 + 0.5) * h;
                        let mut v = 0.0;
                        accumulate_shell(&shell, [x, y, z], &weights, &mut scratch, &mut v);
                        integral += v * v;
                    }
                }
            }
            integral *= h * h * h;
            assert!(
                (integral - 1.0).abs() < 5e-3,
                "cartesian d component {slot} had norm {integral}"
            );
        }
    }

    /// The five pure d functions must also be individually normalised, and
    /// mutually orthogonal.  This is the check that catches a wrong Cartesian
    /// normalisation feeding the Schlegel transform.
    #[test]
    fn pure_d_functions_are_orthonormal() {
        let alpha = 0.9;
        let mut shell = Shell {
            center: 0,
            l: 2,
            spherical: true,
            exponents: vec![alpha],
            raw_coefficients: vec![1.0],
            coefficients: vec![],
        };
        shell.normalize();
        let mut scratch = EvalScratch::default();

        let n = 70;
        let extent = 4.0;
        let h = 2.0 * extent / n as f64;

        let mut overlap = [[0.0f64; 5]; 5];
        for ix in 0..n {
            let x = -extent + (ix as f64 + 0.5) * h;
            for iy in 0..n {
                let y = -extent + (iy as f64 + 0.5) * h;
                for iz in 0..n {
                    let z = -extent + (iz as f64 + 0.5) * h;
                    let mut values = [0.0f64; 5];
                    for (slot, value) in values.iter_mut().enumerate() {
                        let mut weights = vec![0.0; 5];
                        weights[slot] = 1.0;
                        accumulate_shell(&shell, [x, y, z], &weights, &mut scratch, value);
                    }
                    for a in 0..5 {
                        for b in 0..5 {
                            overlap[a][b] += values[a] * values[b];
                        }
                    }
                }
            }
        }

        for a in 0..5 {
            for b in 0..5 {
                let s = overlap[a][b] * h * h * h;
                let expected = if a == b { 1.0 } else { 0.0 };
                assert!(
                    (s - expected).abs() < 5e-3,
                    "pure d overlap[{a}][{b}] = {s}, expected {expected}"
                );
            }
        }
    }
}
