//! Bounded configurational sampling of the single unstable saddle direction.
//! An imaginary mode is not a thermal oscillator and contributes no ZPE here.
use super::input::{Fragment, BOHR};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    Both,
    Positive,
    Negative,
}
impl Side {
    pub fn label(self) -> &'static str {
        match self {
            Self::Both => "Both sides",
            Self::Positive => "Positive only",
            Self::Negative => "Negative only",
        }
    }
}
#[derive(Clone, Debug)]
pub struct Settings {
    pub enabled: bool,
    pub max_displacement_angstrom: f64,
    pub side: Side,
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            enabled: true,
            max_displacement_angstrom: 0.2,
            side: Side::Both,
        }
    }
}
pub fn validate(f: &Fragment, s: &Settings) -> Result<(), String> {
    if s.enabled && f.imaginary_cm1.len() == 1 {
        if !s.max_displacement_angstrom.is_finite() || s.max_displacement_angstrom <= 0. {
            return Err("TS-mode maximum displacement must be positive and finite (Å).".into());
        }
        if f.imaginary_modes.dim() != (f.equilibrium.len(), 1) || peak(f) <= 0. {
            return Err("The projected TS eigenvector is missing or invalid.".into());
        }
    }
    Ok(())
}
fn peak(f: &Fragment) -> f64 {
    (0..f.atoms.len())
        .map(|a| {
            (0..3)
                .map(|k| f.imaginary_modes[(3 * a + k, 0)].powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .fold(0., f64::max)
}
pub fn apply(x: &mut [f64], f: &Fragment, s: &Settings, uniform: f64) -> Result<(), String> {
    if !s.enabled || f.imaginary_cm1.len() != 1 {
        return Ok(());
    }
    validate(f, s)?;
    let signed = match s.side {
        Side::Both => 2. * uniform - 1.,
        Side::Positive => uniform,
        Side::Negative => -uniform,
    };
    let factor = signed * s.max_displacement_angstrom / BOHR / peak(f);
    for (i, value) in x.iter_mut().enumerate() {
        *value += factor * f.imaginary_modes[(i, 0)];
    }
    Ok(())
}
