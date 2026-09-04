//! Gaussian basis model and normalisation.
//!
//! A Molden `[GTO]` shell is a contraction of primitive Gaussians sharing a
//! centre and an angular momentum.  The contraction coefficients as written in
//! the file are *not* directly usable: they assume the primitives carry their
//! own normalisation, and most producers additionally expect the contracted
//! radial function to be renormalised to unit norm.
//!
//! Convention used here, matching Molden's `pnorm` treatment in `gaussian.f`:
//! every Cartesian component is individually normalised.  For a d shell that
//! means `xx` carries `1/sqrt(3)` relative to the shell factor while `xy`
//! carries `1`, i.e. `xy` is larger by `sqrt(3)` -- precisely Molden's
//! `pnorm(8..10) = sqrt(3)`.  This is also the convention the Schlegel
//! transform in [`super::cart2sph`] expects on its input.

use super::cart2sph::{n_cartesian, n_spherical};

/// One contracted Gaussian shell.
#[derive(Debug, Clone)]
pub struct Shell {
    /// Index of the atom this shell sits on.
    pub center: usize,
    /// Angular momentum: 0 = s, 1 = p, 2 = d, 3 = f, 4 = g.
    pub l: u8,
    /// True when the file declared this angular momentum pure (`[5D]` etc).
    pub spherical: bool,
    pub exponents: Vec<f64>,
    /// Contraction coefficients exactly as read from the file.
    pub raw_coefficients: Vec<f64>,
    /// Coefficients with primitive normalisation folded in and the contraction
    /// renormalised.  Filled by [`Shell::normalize`].
    pub coefficients: Vec<f64>,
}

impl Shell {
    /// Number of basis functions this shell contributes.
    pub fn n_functions(&self) -> usize {
        if self.spherical {
            n_spherical(self.l as usize)
        } else {
            n_cartesian(self.l as usize)
        }
    }

    /// Fold primitive normalisation into the coefficients and renormalise the
    /// contracted radial function to unit norm.
    ///
    /// The primitive factor is the one shared by the whole shell,
    /// `(2a/pi)^(3/4) (4a)^(l/2) / sqrt((2l-1)!!)`; the per-component
    /// difference between `xx` and `xy` is applied later by
    /// [`cartesian_component_norm`], so that a single radial evaluation serves
    /// every component of the shell.
    pub fn normalize(&mut self) {
        let l = self.l as usize;
        let df = double_factorial_odd(l);

        // Self-overlap of the contraction. S_ij is the overlap of *normalised*
        // primitives, so the sum uses the raw coefficients; folding the
        // primitive factor in first would count it twice.
        //   S_ij = (2 sqrt(ai aj) / (ai + aj))^(l + 3/2)
        let mut norm_squared = 0.0;
        for (i, &ai) in self.exponents.iter().enumerate() {
            for (j, &aj) in self.exponents.iter().enumerate() {
                let s = (2.0 * (ai * aj).sqrt() / (ai + aj)).powf(l as f64 + 1.5);
                norm_squared += self.raw_coefficients[i] * self.raw_coefficients[j] * s;
            }
        }
        let contraction_scale = if norm_squared > 0.0 {
            1.0 / norm_squared.sqrt()
        } else {
            1.0
        };

        self.coefficients = self
            .exponents
            .iter()
            .zip(&self.raw_coefficients)
            .map(|(&a, &c)| {
                let n = (2.0 * a / std::f64::consts::PI).powf(0.75)
                    * (4.0 * a).powf(l as f64 / 2.0)
                    / df.sqrt();
                c * n * contraction_scale
            })
            .collect();
    }
}

/// Relative normalisation of one Cartesian component within a shell.
///
/// This is Molden's `pnorm`: `sqrt((2l-1)!! / ((2lx-1)!! (2ly-1)!! (2lz-1)!!))`.
/// It is `1` for `x^l` and grows for mixed components -- `sqrt(3)` for `xy`,
/// `sqrt(15)` for `xyz` at f, and so on.
pub fn cartesian_component_norm(powers: [u8; 3]) -> f64 {
    let l = (powers[0] + powers[1] + powers[2]) as usize;
    let denominator = double_factorial_odd(powers[0] as usize)
        * double_factorial_odd(powers[1] as usize)
        * double_factorial_odd(powers[2] as usize);
    (double_factorial_odd(l) / denominator).sqrt()
}

/// `(2n-1)!!`, with `(-1)!! = 1`.
pub fn double_factorial_odd(n: usize) -> f64 {
    let mut acc = 1.0;
    let mut k = 2 * n as isize - 1;
    while k > 1 {
        acc *= k as f64;
        k -= 2;
    }
    acc
}

/// Basis-function layout: which shell each AO index belongs to.
#[derive(Debug, Clone, Default)]
pub struct BasisLayout {
    /// First basis-function index of each shell.
    pub shell_offsets: Vec<usize>,
    pub n_basis: usize,
}

impl BasisLayout {
    pub fn build(shells: &[Shell]) -> Self {
        let mut shell_offsets = Vec::with_capacity(shells.len());
        let mut n_basis = 0;
        for s in shells {
            shell_offsets.push(n_basis);
            n_basis += s.n_functions();
        }
        Self {
            shell_offsets,
            n_basis,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn double_factorial_matches_known_values() {
        assert_eq!(double_factorial_odd(0), 1.0); // (-1)!!
        assert_eq!(double_factorial_odd(1), 1.0); // 1!!
        assert_eq!(double_factorial_odd(2), 3.0); // 3!!
        assert_eq!(double_factorial_odd(3), 15.0); // 5!!
        assert_eq!(double_factorial_odd(4), 105.0); // 7!!
    }

    #[test]
    fn component_norms_reproduce_moldens_pnorm() {
        // Molden pnorm: 1 for x^l, sqrt(3) for xy/xz/yz at d.
        assert!((cartesian_component_norm([2, 0, 0]) - 1.0).abs() < 1e-12);
        assert!((cartesian_component_norm([1, 1, 0]) - 3.0_f64.sqrt()).abs() < 1e-12);
        // f: xyz carries sqrt(15).
        assert!((cartesian_component_norm([1, 1, 1]) - 15.0_f64.sqrt()).abs() < 1e-12);
        // f: x^2 y carries sqrt(5).
        assert!((cartesian_component_norm([2, 1, 0]) - 5.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn a_single_primitive_s_shell_normalises_to_unit_norm() {
        let mut shell = Shell {
            center: 0,
            l: 0,
            spherical: true,
            exponents: vec![1.3],
            raw_coefficients: vec![1.0],
            coefficients: vec![],
        };
        shell.normalize();
        // For one primitive the contraction renormalisation is a no-op, so the
        // coefficient is exactly the primitive normalisation constant.
        let expected = (2.0 * 1.3 / std::f64::consts::PI).powf(0.75);
        assert!((shell.coefficients[0] - expected).abs() < 1e-12);
    }

    #[test]
    fn a_contracted_shell_is_renormalised_to_unit_norm() {
        let mut shell = Shell {
            center: 0,
            l: 1,
            spherical: true,
            exponents: vec![3.0, 0.7, 0.2],
            raw_coefficients: vec![0.15, 0.55, 0.45],
            coefficients: vec![],
        };
        shell.normalize();

        // Rebuild the norm from the finished coefficients by dividing the
        // primitive factor back out, which is what the evaluator effectively
        // does when it sums c_i exp(-a_i r^2).
        let l = 1.0_f64;
        let df = double_factorial_odd(1);
        let unit: Vec<f64> = shell
            .exponents
            .iter()
            .zip(&shell.coefficients)
            .map(|(&a, &c)| {
                let n =
                    (2.0 * a / std::f64::consts::PI).powf(0.75) * (4.0 * a).powf(0.5) / df.sqrt();
                c / n
            })
            .collect();
        let mut norm2 = 0.0;
        for (i, &ai) in shell.exponents.iter().enumerate() {
            for (j, &aj) in shell.exponents.iter().enumerate() {
                let s = (2.0 * (ai * aj).sqrt() / (ai + aj)).powf(l + 1.5);
                norm2 += unit[i] * unit[j] * s;
            }
        }
        assert!((norm2 - 1.0).abs() < 1e-12, "contraction norm was {norm2}");
    }

    #[test]
    fn function_counts_respect_the_spherical_flag() {
        let mk = |l, spherical| Shell {
            center: 0,
            l,
            spherical,
            exponents: vec![1.0],
            raw_coefficients: vec![1.0],
            coefficients: vec![],
        };
        assert_eq!(mk(2, true).n_functions(), 5); // [5D]
        assert_eq!(mk(2, false).n_functions(), 6); // 6D
        assert_eq!(mk(3, true).n_functions(), 7); // [7F]
        assert_eq!(mk(4, true).n_functions(), 9); // [9G]
    }

    #[test]
    fn layout_offsets_accumulate() {
        let mk = |l, spherical| Shell {
            center: 0,
            l,
            spherical,
            exponents: vec![1.0],
            raw_coefficients: vec![1.0],
            coefficients: vec![],
        };
        let shells = vec![mk(0, true), mk(1, true), mk(2, true)];
        let layout = BasisLayout::build(&shells);
        assert_eq!(layout.shell_offsets, vec![0, 1, 4]);
        assert_eq!(layout.n_basis, 9);
    }
}
