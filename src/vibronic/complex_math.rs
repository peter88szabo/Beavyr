//! Small complex arithmetic and pivoted dense solves for harmonic Gaussian integrals.
//! No BLAS or additional dependency is required.
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct C {
    pub re: f64,
    pub im: f64,
}
impl C {
    pub const ZERO: Self = Self { re: 0., im: 0. };
    pub const ONE: Self = Self { re: 1., im: 0. };
    pub const I: Self = Self { re: 0., im: 1. };
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    pub fn abs2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }
    pub fn exp(self) -> Self {
        let r = self.re.exp();
        Self::new(r * self.im.cos(), r * self.im.sin())
    }
    pub fn ln(self) -> Self {
        Self::new(self.abs().ln(), self.im.atan2(self.re))
    }
    pub fn finite(self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
}
impl From<f64> for C {
    fn from(v: f64) -> Self {
        Self::new(v, 0.)
    }
}
impl Add for C {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self::new(self.re + b.re, self.im + b.im)
    }
}
impl Sub for C {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self::new(self.re - b.re, self.im - b.im)
    }
}
impl Neg for C {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}
impl Mul for C {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self::new(
            self.re * b.re - self.im * b.im,
            self.re * b.im + self.im * b.re,
        )
    }
}
impl Div for C {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        let a = b.abs2();
        Self::new(
            (self.re * b.re + self.im * b.im) / a,
            (self.im * b.re - self.re * b.im) / a,
        )
    }
}
impl Mul<f64> for C {
    type Output = Self;
    fn mul(self, b: f64) -> Self {
        Self::new(self.re * b, self.im * b)
    }
}
impl Div<f64> for C {
    type Output = Self;
    fn div(self, b: f64) -> Self {
        self * (1. / b)
    }
}
impl AddAssign for C {
    fn add_assign(&mut self, b: Self) {
        *self = *self + b;
    }
}
impl SubAssign for C {
    fn sub_assign(&mut self, b: Self) {
        *self = *self - b;
    }
}
impl MulAssign for C {
    fn mul_assign(&mut self, b: Self) {
        *self = *self * b;
    }
}
/// An LU decomposition with row pivoting. Complex symmetric matrices are not Hermitian.
pub struct Lu {
    a: Vec<C>,
    pivots: Vec<usize>,
    n: usize,
    pub logdet: C,
}
impl Lu {
    pub fn new(mut a: Vec<C>, n: usize) -> Result<Self, String> {
        if n == 0 || a.len() != n * n || a.iter().any(|z| !z.finite()) {
            return Err("Invalid Gaussian matrix".into());
        }
        let mut pivots = vec![];
        let mut logdet = C::ZERO;
        for k in 0..n {
            let p = (k..n)
                .max_by(|&i, &j| a[i * n + k].abs2().total_cmp(&a[j * n + k].abs2()))
                .unwrap();
            if a[p * n + k].abs() < 1e-280 {
                return Err("Singular harmonic Gaussian integral".into());
            }
            pivots.push(p);
            if p != k {
                for j in 0..n {
                    a.swap(k * n + j, p * n + j);
                }
                logdet += C::new(0., std::f64::consts::PI);
            }
            let diag = a[k * n + k];
            logdet += diag.ln();
            for i in k + 1..n {
                a[i * n + k] = a[i * n + k] / diag;
                let l = a[i * n + k];
                for j in k + 1..n {
                    let v = l * a[k * n + j];
                    a[i * n + j] -= v;
                }
            }
        }
        Ok(Self {
            a,
            pivots,
            n,
            logdet,
        })
    }
    pub fn solve(&self, b: &[C]) -> Vec<C> {
        let n = self.n;
        let mut x = b.to_vec();
        for k in 0..n {
            x.swap(k, self.pivots[k]);
        }
        for i in 0..n {
            for j in 0..i {
                let v = self.a[i * n + j] * x[j];
                x[i] -= v;
            }
        }
        for i in (0..n).rev() {
            for j in i + 1..n {
                let v = self.a[i * n + j] * x[j];
                x[i] -= v;
            }
            x[i] = x[i] / self.a[i * n + i];
        }
        x
    }
    pub fn inverse(&self) -> Vec<C> {
        let n = self.n;
        let mut out = vec![C::ZERO; n * n];
        for j in 0..n {
            let mut b = vec![C::ZERO; n];
            b[j] = C::ONE;
            let x = self.solve(&b);
            for i in 0..n {
                out[i * n + j] = x[i];
            }
        }
        out
    }
}
/// Inverse DFT, with 1/N normalization. Input length must be a power of two.
pub fn inverse_fft(values: &mut [C]) {
    let n = values.len();
    assert!(n.is_power_of_two());
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            values.swap(i, j);
        }
    }
    let mut width = 2;
    while width <= n {
        let step = C::new(0., std::f64::consts::TAU / width as f64).exp();
        for start in (0..n).step_by(width) {
            let mut w = C::ONE;
            for k in 0..width / 2 {
                let a = values[start + k];
                let b = w * values[start + k + width / 2];
                values[start + k] = a + b;
                values[start + k + width / 2] = a - b;
                w *= step;
            }
        }
        width *= 2;
    }
    for v in values {
        *v = *v / n as f64;
    }
}
