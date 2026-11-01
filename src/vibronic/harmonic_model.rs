//! Two harmonic electronic surfaces in a common mass-weighted coordinate frame.
use super::{
    complex_math::{Lu, C},
    duschinsky,
    input::{self, HessianSource, StateModes},
    CM1_TO_AU,
};
use ndarray::Array2;
#[derive(Clone, Debug)]
pub struct HarmonicModel {
    pub initial: Vec<f64>,
    pub final_: Vec<f64>,
    /// Q_initial = J Q_final + displacement; all Q are mass-weighted atomic units.
    pub j: Array2<f64>,
    pub displacement: Vec<f64>,
    /// Electronic minimum-to-minimum energy gap, Hartree (final minus initial).
    pub electronic_gap: f64,
    pub initial_labels: Vec<usize>,
    pub final_labels: Vec<usize>,
}
impl HarmonicModel {
    pub fn n(&self) -> usize {
        self.initial.len()
    }
    pub fn zero_zero(&self) -> f64 {
        self.electronic_gap + self.zpe_difference()
    }
    pub fn zpe_difference(&self) -> f64 {
        0.5 * (self.final_.iter().sum::<f64>() - self.initial.iter().sum::<f64>())
    }
    pub fn validate(&self) -> Result<(), String> {
        let n = self.n();
        if n == 0
            || self.final_.len() != n
            || self.j.dim() != (n, n)
            || self.displacement.len() != n
            || self
                .initial
                .iter()
                .chain(&self.final_)
                .any(|v| !v.is_finite() || *v <= 0.)
            || self
                .j
                .iter()
                .chain(&self.displacement)
                .any(|v| !v.is_finite())
            || !self.electronic_gap.is_finite()
        {
            return Err("Invalid two-state harmonic model".into());
        }
        Ok(())
    }
    pub fn inverse_rotation(&self) -> Result<(Vec<f64>, f64), String> {
        self.validate()?;
        let lu = Lu::new(self.j.iter().map(|v| C::from(*v)).collect(), self.n())?;
        let inverse: Vec<_> = lu.inverse().iter().map(|v| v.re).collect();
        if inverse.iter().any(|v| !v.is_finite() || v.abs() > 1e6) {
            return Err(
                "The mode mapping is nearly singular; inspect the two geometries and mode sets"
                    .into(),
            );
        }
        Ok((inverse, -lu.logdet.re))
    }
    pub fn reversed(&self) -> Result<Self, String> {
        let (inv, _) = self.inverse_rotation()?;
        let n = self.n();
        Ok(Self {
            initial: self.final_.clone(),
            final_: self.initial.clone(),
            j: Array2::from_shape_vec((n, n), inv.clone()).unwrap(),
            displacement: (0..n)
                .map(|i| {
                    -(0..n)
                        .map(|k| inv[i * n + k] * self.displacement[k])
                        .sum::<f64>()
                })
                .collect(),
            electronic_gap: -self.electronic_gap,
            initial_labels: self.final_labels.clone(),
            final_labels: self.initial_labels.clone(),
        })
    }
}
pub fn from_states(a: &StateModes, b: &StateModes, gap: f64) -> Result<HarmonicModel, String> {
    input::same_atoms(&a.atoms, &b.atoms)?;
    if a.imaginary > 0 || b.imaginary > 0 {
        return Err("A two-minimum vibronic spectrum requires stable minima: one Hessian has imaginary modes. Optimise that state before proceeding.".into());
    }
    if a.selected.len() != b.selected.len() {
        return Err("Both states must have the same number of physical vibrational modes".into());
    }
    if a.masses
        .iter()
        .zip(&b.masses)
        .any(|(x, y)| (x - y).abs() > 1e-5 * x.abs().max(1.))
    {
        return Err("The two states must use identical isotope masses".into());
    }
    let d = duschinsky(
        &a.modes,
        &a.frequencies,
        &a.coords,
        &b.modes,
        &b.frequencies,
        &b.coords,
        &a.masses,
        &a.selected,
        &b.selected,
    )
    .map_err(|e| e.to_string())?;
    let initial: Vec<_> = a
        .selected
        .iter()
        .map(|&i| a.frequencies[i] * CM1_TO_AU)
        .collect();
    let model = HarmonicModel {
        displacement: d
            .displacement
            .iter()
            .zip(&initial)
            .map(|(k, w)| k / w.sqrt())
            .collect(),
        initial,
        final_: b
            .selected
            .iter()
            .map(|&i| b.frequencies[i] * CM1_TO_AU)
            .collect(),
        j: d.rotation,
        electronic_gap: gap,
        initial_labels: a.selected.clone(),
        final_labels: b.selected.clone(),
    };
    model.validate()?;
    Ok(model)
}
pub fn load(
    a: &HessianSource,
    b: &HessianSource,
    gap: f64,
) -> Result<(HarmonicModel, StateModes, StateModes), String> {
    let a = input::load(a, 0.)?;
    let b = input::load(b, 0.)?;
    let model = from_states(&a, &b, gap)?;
    Ok((model, a, b))
}
