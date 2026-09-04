//! Cartesian-to-real-spherical transformation for Gaussian shells.
//!
//! The coefficient generator is adapted from Behemoth's
//! `src/integrals/cartesian_to_spherical.rs`, which implements equation 15 of
//! Schlegel and Frisch, *Int. J. Quantum Chem.* **54**, 83-87 (1995) with the
//! Condon-Shortley phase.  Only the coefficient table is needed here; the
//! integral-block machinery around it in Behemoth is not.
//!
//! Two ordering conventions meet in this file and must not be confused:
//!
//! | | Cartesian order | Spherical order |
//! |---|---|---|
//! | Schlegel generator (below) | `lx` descending, then `ly` descending | `m = -l ..= +l` |
//! | Molden file format | `xx, yy, zz, xy, xz, yz` (d) | `0, +1, -1, +2, -2` |
//!
//! [`molden_spherical_row`] and [`molden_cartesian_column`] are the only places
//! that know about Molden's ordering, so a convention bug has exactly one home.

/// Molden's `[9G]` is the highest angular momentum the format defines.
pub const MAX_ANGULAR_MOMENTUM: usize = 4;

/// Number of Cartesian components in a shell.
#[inline]
pub const fn n_cartesian(l: usize) -> usize {
    (l + 1) * (l + 2) / 2
}

/// Number of real spherical components in a shell.
#[inline]
pub const fn n_spherical(l: usize) -> usize {
    2 * l + 1
}

/// Cartesian exponents `(lx, ly, lz)` in the generator's own order.
///
/// `lx` descends outermost, then `ly`.  For d this yields
/// `xx, xy, xz, yy, yz, zz`.
pub fn cartesian_powers(l: usize) -> &'static [[u8; 3]] {
    use std::sync::OnceLock;
    static TABLES: OnceLock<Vec<Vec<[u8; 3]>>> = OnceLock::new();
    &TABLES.get_or_init(|| {
        (0..=MAX_ANGULAR_MOMENTUM)
            .map(|l| {
                let mut out = Vec::with_capacity(n_cartesian(l));
                for lx in (0..=l).rev() {
                    for ly in (0..=l - lx).rev() {
                        out.push([lx as u8, ly as u8, (l - lx - ly) as u8]);
                    }
                }
                out
            })
            .collect()
    })[l]
}

/// Position of spherical component `m` within a Molden `[5D]/[7F]/[9G]` block.
///
/// Molden lists `m = 0`, then `+1, -1, +2, -2, ...`, so slot `p` holds
/// `m = 0` for `p == 0`, `m = +(p+1)/2` for odd `p`, and `m = -(p/2)` for even
/// `p > 0`.  Returns the row index into the generator's `m = -l ..= +l` table.
pub fn molden_spherical_row(l: usize, slot: usize) -> usize {
    debug_assert!(slot < n_spherical(l));
    let m: isize = if slot == 0 {
        0
    } else if slot % 2 == 1 {
        ((slot + 1) / 2) as isize
    } else {
        -((slot / 2) as isize)
    };
    (m + l as isize) as usize
}

/// Column index in the generator's Cartesian order for Molden slot `slot`.
///
/// Molden's Cartesian orders are fixed tables rather than a loop nest, so they
/// are written out literally and mapped back onto `cartesian_powers`.
pub fn molden_cartesian_column(l: usize, slot: usize) -> usize {
    use std::sync::OnceLock;
    static TABLES: OnceLock<Vec<Vec<usize>>> = OnceLock::new();
    TABLES.get_or_init(|| {
        (0..=MAX_ANGULAR_MOMENTUM)
            .map(|l| {
                molden_cartesian_powers(l)
                    .iter()
                    .map(|target| {
                        cartesian_powers(l)
                            .iter()
                            .position(|p| p == target)
                            .expect("every Molden component exists in generator order")
                    })
                    .collect()
            })
            .collect()
    })[l][slot]
}

/// Cartesian exponents in Molden's documented order.
pub fn molden_cartesian_powers(l: usize) -> &'static [[u8; 3]] {
    use std::sync::OnceLock;
    static TABLES: OnceLock<Vec<Vec<[u8; 3]>>> = OnceLock::new();
    &TABLES.get_or_init(|| {
        (0..=MAX_ANGULAR_MOMENTUM)
            .map(molden_cartesian_powers_uncached)
            .collect()
    })[l]
}

fn molden_cartesian_powers_uncached(l: usize) -> Vec<[u8; 3]> {
    // Molden writes Cartesian shells in these fixed sequences.
    const D: [[u8; 3]; 6] = [
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 2],
        [1, 1, 0],
        [1, 0, 1],
        [0, 1, 1],
    ];
    const F: [[u8; 3]; 10] = [
        [3, 0, 0],
        [0, 3, 0],
        [0, 0, 3],
        [1, 2, 0],
        [2, 1, 0],
        [2, 0, 1],
        [1, 0, 2],
        [0, 1, 2],
        [0, 2, 1],
        [1, 1, 1],
    ];
    const G: [[u8; 3]; 15] = [
        [4, 0, 0],
        [0, 4, 0],
        [0, 0, 4],
        [3, 1, 0],
        [3, 0, 1],
        [1, 3, 0],
        [0, 3, 1],
        [1, 0, 3],
        [0, 1, 3],
        [2, 2, 0],
        [2, 0, 2],
        [0, 2, 2],
        [2, 1, 1],
        [1, 2, 1],
        [1, 1, 2],
    ];
    match l {
        0 => vec![[0, 0, 0]],
        1 => vec![[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        2 => D.to_vec(),
        3 => F.to_vec(),
        4 => G.to_vec(),
        _ => panic!("Molden format defines Cartesian shells only through g"),
    }
}

/// Row-major table with spherical rows (`m = -l ..= +l`) and Cartesian columns
/// (generator order).  Cached per angular momentum.
pub fn transform(l: usize) -> &'static [f64] {
    use std::sync::OnceLock;
    static TABLES: OnceLock<Vec<Vec<f64>>> = OnceLock::new();
    &TABLES.get_or_init(|| (0..=MAX_ANGULAR_MOMENTUM).map(build_transform).collect())[l]
}

fn build_transform(l: usize) -> Vec<f64> {
    let n_cart = n_cartesian(l);
    let n_sph = n_spherical(l);
    let mut matrix = vec![0.0; n_sph * n_cart];

    // s and p are already pure; the generator's order is the conventional one.
    if l <= 1 {
        for c in 0..n_cart {
            matrix[c * n_cart + c] = 1.0;
        }
        return matrix;
    }

    let powers = cartesian_powers(l);
    for sph in 0..n_sph {
        let m = sph as isize - l as isize;
        for (cart, p) in powers.iter().enumerate() {
            let c = normalized_real_coefficient(p[0] as usize, p[1] as usize, p[2] as usize, m);
            matrix[sph * n_cart + cart] = if c.abs() < 1.0e-15 { 0.0 } else { c };
        }
    }
    matrix
}

#[derive(Clone, Copy, Debug, Default)]
struct Complex {
    real: f64,
    imaginary: f64,
}

impl Complex {
    fn scale(self, f: f64) -> Self {
        Self {
            real: f * self.real,
            imaginary: f * self.imaginary,
        }
    }
    fn add(self, rhs: Self) -> Self {
        Self {
            real: self.real + rhs.real,
            imaginary: self.imaginary + rhs.imaginary,
        }
    }
    fn sub(self, rhs: Self) -> Self {
        Self {
            real: self.real - rhs.real,
            imaginary: self.imaginary - rhs.imaginary,
        }
    }
}

/// Coefficient of one normalized Cartesian Gaussian in a normalized complex
/// pure spherical Gaussian: Schlegel and Frisch equation 15.
fn normalized_complex_coefficient(lx: usize, ly: usize, lz: usize, m: isize) -> Complex {
    let l = lx + ly + lz;
    let abs_m = m.unsigned_abs();
    if (lx + ly + abs_m) % 2 != 0 {
        return Complex::default();
    }

    let j = ((lx + ly) as isize - abs_m as isize) / 2;

    let numerator = factorial(2 * lx)
        * factorial(2 * ly)
        * factorial(2 * lz)
        * factorial(l)
        * factorial(l - abs_m);
    let denominator =
        factorial(2 * l) * factorial(l + abs_m) * factorial(lx) * factorial(ly) * factorial(lz);
    let mut leading = (numerator / denominator).sqrt() * 0.5_f64.powi(l as i32) / factorial(l);

    let max_i = (l - abs_m) / 2;
    let mut angular_sum = 0.0;
    for i in 0..=max_i {
        if (i as isize) < j {
            continue;
        }
        let term =
            binomial(l as isize, i as isize) * binomial(i as isize, j) * factorial(2 * l - 2 * i)
                / factorial(l - abs_m - 2 * i);
        if i % 2 == 0 {
            angular_sum += term;
        } else {
            angular_sum -= term;
        }
    }
    leading *= angular_sum;

    let mut azimuthal = Complex::default();
    for k in 0..=l {
        let cartesian_power = lx as isize - 2 * k as isize;
        let term = binomial(j, k as isize) * binomial(abs_m as isize, cartesian_power);
        if term == 0.0 {
            continue;
        }
        match (abs_m as isize - lx as isize + 2 * k as isize).rem_euclid(4) {
            0 => azimuthal.real += term,
            1 => azimuthal.imaginary += term,
            2 => azimuthal.real -= term,
            3 => azimuthal.imaginary -= term,
            _ => unreachable!(),
        }
    }

    let mut result = azimuthal.scale(leading);
    if m >= 0 {
        if m % 2 != 0 {
            result = result.scale(-1.0);
        }
    } else {
        result.imaginary = -result.imaginary;
    }
    result
}

/// Real spherical coefficient, rows ordered `m = -l ..= +l`.
fn normalized_real_coefficient(lx: usize, ly: usize, lz: usize, m: isize) -> f64 {
    if m == 0 {
        return normalized_complex_coefficient(lx, ly, lz, 0).real;
    }

    let abs_m = m.unsigned_abs() as isize;
    let negative = normalized_complex_coefficient(lx, ly, lz, -abs_m);
    let positive = normalized_complex_coefficient(lx, ly, lz, abs_m);
    let condon = if abs_m % 2 == 0 { 1.0 } else { -1.0 };
    let inv_sqrt2 = std::f64::consts::FRAC_1_SQRT_2;

    if m > 0 {
        negative.add(positive.scale(condon)).scale(inv_sqrt2).real
    } else {
        -negative
            .sub(positive.scale(condon))
            .scale(inv_sqrt2)
            .imaginary
    }
}

fn factorial(n: usize) -> f64 {
    (1..=n).map(|k| k as f64).product()
}

fn binomial(n: isize, k: isize) -> f64 {
    if n < 0 || k < 0 || k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut acc = 1.0;
    for i in 0..k {
        acc = acc * (n - i) as f64 / (i + 1) as f64;
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn component_counts_match_the_format() {
        assert_eq!((n_cartesian(2), n_spherical(2)), (6, 5)); // 6D / 5D
        assert_eq!((n_cartesian(3), n_spherical(3)), (10, 7)); // 10F / 7F
        assert_eq!((n_cartesian(4), n_spherical(4)), (15, 9)); // 15G / 9G
    }

    #[test]
    fn molden_spherical_slots_follow_zero_plus_minus() {
        // [5D] is documented as D 0, D+1, D-1, D+2, D-2.
        let l = 2;
        let rows: Vec<isize> = (0..5)
            .map(|slot| molden_spherical_row(l, slot) as isize - l as isize)
            .collect();
        assert_eq!(rows, vec![0, 1, -1, 2, -2]);

        let l = 3;
        let rows: Vec<isize> = (0..7)
            .map(|slot| molden_spherical_row(l, slot) as isize - l as isize)
            .collect();
        assert_eq!(rows, vec![0, 1, -1, 2, -2, 3, -3]);
    }

    #[test]
    fn molden_cartesian_slots_map_onto_generator_columns() {
        // Molden d order xx, yy, zz, xy, xz, yz against generator order
        // xx, xy, xz, yy, yz, zz.
        let cols: Vec<usize> = (0..6).map(|s| molden_cartesian_column(2, s)).collect();
        assert_eq!(cols, vec![0, 3, 5, 1, 2, 4]);
    }

    #[test]
    fn every_molden_slot_maps_to_a_distinct_column() {
        for l in 0..=MAX_ANGULAR_MOMENTUM {
            let mut seen: Vec<usize> = (0..n_cartesian(l))
                .map(|s| molden_cartesian_column(l, s))
                .collect();
            let total = seen.len();
            seen.sort_unstable();
            seen.dedup();
            assert_eq!(seen.len(), total, "duplicate column mapping for l={l}");
        }
    }

    /// Overlap of two normalised Cartesian Gaussians on the same centre with
    /// the same exponent.  Cartesian components are *not* mutually orthogonal,
    /// so this metric is needed to state the transform's correctness.
    fn cartesian_overlap(p: [u8; 3], q: [u8; 3], a: f64) -> f64 {
        // int x^n exp(-2a x^2) dx = (n-1)!! / (4a)^(n/2) * sqrt(pi / 2a), n even.
        fn axis(n: usize, a: f64) -> f64 {
            if n % 2 == 1 {
                return 0.0;
            }
            double_factorial_odd(n / 2) / (4.0 * a).powf(n as f64 / 2.0)
                * (std::f64::consts::PI / (2.0 * a)).sqrt()
        }
        fn double_factorial_odd(n: usize) -> f64 {
            let mut acc = 1.0;
            let mut k = 2 * n as isize - 1;
            while k > 1 {
                acc *= k as f64;
                k -= 2;
            }
            acc
        }
        fn norm(p: [u8; 3], a: f64) -> f64 {
            let l = (p[0] + p[1] + p[2]) as usize;
            let d = double_factorial_odd(p[0] as usize)
                * double_factorial_odd(p[1] as usize)
                * double_factorial_odd(p[2] as usize);
            (2.0 * a / std::f64::consts::PI).powf(0.75) * (4.0 * a).powf(l as f64 / 2.0) / d.sqrt()
        }

        norm(p, a)
            * norm(q, a)
            * axis((p[0] + q[0]) as usize, a)
            * axis((p[1] + q[1]) as usize, a)
            * axis((p[2] + q[2]) as usize, a)
    }

    /// The transform maps normalised Cartesian functions onto orthonormal pure
    /// functions, so `T S Tᵀ = I` where `S` is the Cartesian overlap metric.
    ///
    /// Plain `T Tᵀ = I` does *not* hold: Cartesian Gaussians overlap each other
    /// (`<xx|yy>` is nonzero), which is exactly what `S` accounts for.
    #[test]
    fn transform_is_orthonormal_under_the_cartesian_metric() {
        let a = 1.0;
        for l in 2..=MAX_ANGULAR_MOMENTUM {
            let t = transform(l);
            let nc = n_cartesian(l);
            let ns = n_spherical(l);
            let powers = cartesian_powers(l);

            let s: Vec<f64> = (0..nc)
                .flat_map(|i| {
                    (0..nc)
                        .map(|j| cartesian_overlap(powers[i], powers[j], a))
                        .collect::<Vec<_>>()
                })
                .collect();

            for x in 0..ns {
                for y in 0..ns {
                    let mut acc = 0.0;
                    for i in 0..nc {
                        for j in 0..nc {
                            acc += t[x * nc + i] * s[i * nc + j] * t[y * nc + j];
                        }
                    }
                    let expected = if x == y { 1.0 } else { 0.0 };
                    assert!(
                        (acc - expected).abs() < 1.0e-10,
                        "l={l} pure {x},{y}: {acc} != {expected}"
                    );
                }
            }
        }
    }

    #[test]
    fn s_and_p_shells_are_the_identity() {
        assert_eq!(transform(0), &[1.0]);
        assert_eq!(transform(1), &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
    }
}
