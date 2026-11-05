//! Chemist-facing summaries: the per-mode overview shown before a run and
//! histograms of geometric quantities and energies after it.
use super::{
    input::{Fragment, BOHR, CM_TO_AU},
    output::SamplingOutput,
    rotors::ModeTreatment,
    sampling::KB,
    state::*,
};

/// One row of the mode table.
#[derive(Clone, Debug)]
pub struct ModeRow {
    pub index: usize,
    pub cm1: f64,
    /// Mode energy in hartree: the largest allowed energy, or the mean for
    /// Smite thermal sampling (which has no hard cutoff).
    pub energy: f64,
    /// False when `energy` is a mean rather than a maximum.
    pub is_maximum: bool,
    /// Largest atom displacement at the turning point, Å.
    pub max_displacement: f64,
    pub note: Option<&'static str>,
}

/// Deterministic energy of one mode for the overview (hartree).
pub fn mode_energy(e: &Excitation, w: f64) -> (f64, bool) {
    match e.population {
        Population::Zpe => (0.5 * w, true),
        Population::Quantum => (w * (e.quantum as f64 + 0.5), true),
        Population::Energy => (e.energy_kcal / 627.51, true),
        Population::Thermal => match e.temperature_rule {
            TemperatureRule::BoltzmannTail => (
                w * (super::sampling::tail_quantum(w, e.temperature, e.retained_population) + 0.5),
                true,
            ),
            TemperatureRule::ZpePlusKt => (0.5 * w + KB * e.temperature, true),
            TemperatureRule::SmiteThermal => {
                let x = w / (KB * e.temperature);
                let occupation = if e.temperature > 0. && x < 700. {
                    1. / x.exp_m1()
                } else {
                    0.
                };
                (w * (occupation + 0.5), false)
            }
        },
    }
}

/// Largest atom displacement (Å) of mode `j` per unit normal coordinate.
fn atom_scale(f: &Fragment, j: usize) -> f64 {
    (0..f.atoms.len())
        .map(|i| {
            (0..3)
                .map(|k| f.modes[(3 * i + k, j)].powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .fold(0., f64::max)
        * BOHR
}

/// Per-mode energies and turning-point displacements for a partner.
pub fn mode_overview(f: &Fragment, p: &Partner) -> Vec<ModeRow> {
    f.frequencies
        .iter()
        .enumerate()
        .map(|(j, &w)| {
            let e = if p.use_overrides {
                p.overrides.get(j).unwrap_or(&p.excitation)
            } else {
                &p.excitation
            };
            let rotor = p
                .rotors
                .iter()
                .find(|r| r.enabled && r.mode == j + 1)
                .map(|r| r.treatment);
            let ((energy, is_maximum), note) = match rotor {
                Some(ModeTreatment::ZpeOnly) => ((0.5 * w, true), Some("rotor: ZPE only")),
                Some(ModeTreatment::Replace) => ((0., true), Some("replaced by rotor")),
                None => (mode_energy(e, w), None),
            };
            ModeRow {
                index: j,
                cm1: w / CM_TO_AU,
                energy,
                is_maximum,
                max_displacement: (2. * energy).sqrt() / w * atom_scale(f, j),
                note,
            }
        })
        .collect()
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Quantity {
    #[default]
    Distance,
    Angle,
    Energy,
}
impl Quantity {
    pub fn label(self) -> &'static str {
        match self {
            Self::Distance => "Distance",
            Self::Angle => "Angle",
            Self::Energy => "Harmonic energy",
        }
    }
    pub fn atoms(self) -> usize {
        match self {
            Self::Distance => 2,
            Self::Angle => 3,
            Self::Energy => 0,
        }
    }
}

/// Histogram choice in the panel; atom numbers are 1-based.
#[derive(Clone, Debug)]
pub struct HistogramView {
    pub quantity: Quantity,
    pub atoms: [usize; 3],
    pub bins: usize,
}
impl Default for HistogramView {
    fn default() -> Self {
        Self {
            quantity: Quantity::Distance,
            atoms: [1, 2, 3],
            bins: 30,
        }
    }
}

fn point(x: &[f64], i: usize) -> [f64; 3] {
    [x[3 * i], x[3 * i + 1], x[3 * i + 2]]
}
fn delta(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Values of the chosen quantity for every sample, with its unit label.
pub fn values(
    out: &SamplingOutput,
    view: &HistogramView,
    unit: EnergyUnit,
) -> Result<(Vec<f64>, String), String> {
    let n = view.quantity.atoms();
    let idx: Vec<usize> = view.atoms[..n].iter().map(|a| a.wrapping_sub(1)).collect();
    if idx.iter().any(|&i| i >= out.atoms.len()) {
        return Err(format!("Atom numbers must lie between 1 and {}.", out.atoms.len()));
    }
    if (1..n).any(|k| idx[..k].contains(&idx[k])) {
        return Err("Choose different atoms.".into());
    }
    let samples = out.samples.iter();
    Ok(match view.quantity {
        Quantity::Distance => (
            samples
                .map(|s| norm(delta(point(&s.coordinates, idx[0]), point(&s.coordinates, idx[1]))))
                .collect(),
            "Å".into(),
        ),
        Quantity::Angle => (
            samples
                .map(|s| {
                    let b = point(&s.coordinates, idx[1]);
                    let u = delta(point(&s.coordinates, idx[0]), b);
                    let v = delta(point(&s.coordinates, idx[2]), b);
                    let c = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (norm(u) * norm(v));
                    c.clamp(-1., 1.).acos().to_degrees()
                })
                .collect(),
            "°".into(),
        ),
        Quantity::Energy => (
            samples.map(|s| unit.from_hartree(s.energy)).collect(),
            unit.label().into(),
        ),
    })
}

/// Equal-width bins between the smallest and largest value.
pub fn histogram(values: &[f64], bins: usize) -> (f64, f64, Vec<usize>) {
    let bins = bins.max(1);
    let lo = values.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut counts = vec![0; bins];
    if values.is_empty() {
        return (0., 0., counts);
    }
    let width = (hi - lo).max(1e-12);
    for v in values {
        let k = (((v - lo) / width) * bins as f64) as usize;
        counts[k.min(bins - 1)] += 1;
    }
    (lo, hi, counts)
}
