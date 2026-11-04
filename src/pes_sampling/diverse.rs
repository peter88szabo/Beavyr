//! Farthest-point selection of a well-spread subset of sampled geometries.
//!
//! Descriptor: all inverse interatomic distances 1/r_ij (Å⁻¹). It is
//! invariant to translation and rotation and emphasises short contacts, where
//! the PES changes fastest. It is not invariant to permuting identical atoms.
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

#[derive(Clone, Debug)]
pub struct Diversity {
    pub enabled: bool,
    /// Candidates generated per kept sample (≥ 1).
    pub oversampling: f64,
}

impl Default for Diversity {
    fn default() -> Self {
        Self {
            enabled: true,
            oversampling: 3.,
        }
    }
}

/// Upper bound on per-run selection work (kept × candidates × descriptor length).
pub const MAX_WORK: f64 = 2e10;
/// Upper bound on stored descriptor values (f32) for one block.
pub const MAX_DESCRIPTOR_VALUES: usize = 50_000_000;

impl Diversity {
    pub fn validate(&self) -> Result<(), String> {
        if self.enabled && (!self.oversampling.is_finite() || !(1.0..=100.).contains(&self.oversampling)) {
            return Err("Diverse-subset oversampling must lie between 1 and 100.".into());
        }
        Ok(())
    }
    /// Candidates to generate for `kept` final samples.
    pub fn candidates(&self, kept: usize) -> usize {
        if self.enabled {
            ((kept as f64 * self.oversampling).ceil() as usize).max(kept)
        } else {
            kept
        }
    }
}

pub fn descriptor_len(atoms: usize) -> usize {
    atoms * atoms.saturating_sub(1) / 2
}

pub fn descriptor(x: &[f64]) -> Vec<f32> {
    let n = x.len() / 3;
    let mut d = Vec::with_capacity(descriptor_len(n));
    for i in 0..n {
        for j in 0..i {
            let r2: f64 = (0..3).map(|k| (x[3 * i + k] - x[3 * j + k]).powi(2)).sum();
            d.push((1. / r2.sqrt().max(1e-6)) as f32);
        }
    }
    d
}

/// Indices of `kept` candidates in selection order: start with the candidate
/// closest to the candidate mean, then repeatedly take the candidate farthest
/// (Euclidean descriptor distance) from everything already chosen.
pub fn farthest_point(
    candidates: &[Vec<f64>],
    kept: usize,
    cancel: &AtomicBool,
    progress: &AtomicUsize,
) -> Result<Vec<usize>, String> {
    let n = candidates.len();
    if kept >= n {
        return Ok((0..n).collect());
    }
    let descriptors: Vec<Vec<f32>> = candidates.iter().map(|x| descriptor(x)).collect();
    let len = descriptors.first().map_or(0, Vec::len);
    if len == 0 {
        return Ok((0..kept).collect());
    }
    let mut mean = vec![0f64; len];
    for d in &descriptors {
        for (m, v) in mean.iter_mut().zip(d) {
            *m += *v as f64 / n as f64;
        }
    }
    let distance = |a: &[f32], b: &[f32]| -> f64 {
        a.iter().zip(b).map(|(x, y)| ((x - y) as f64).powi(2)).sum()
    };
    let first = (0..n)
        .min_by(|&a, &b| {
            let da: f64 = descriptors[a].iter().zip(&mean).map(|(x, m)| (*x as f64 - m).powi(2)).sum();
            let db: f64 = descriptors[b].iter().zip(&mean).map(|(x, m)| (*x as f64 - m).powi(2)).sum();
            da.total_cmp(&db)
        })
        .unwrap();
    let mut nearest = vec![f64::INFINITY; n];
    let mut chosen = vec![false; n];
    let mut order = Vec::with_capacity(kept);
    let mut current = first;
    loop {
        chosen[current] = true;
        order.push(current);
        progress.fetch_add(1, Ordering::Relaxed);
        if order.len() == kept {
            return Ok(order);
        }
        if cancel.load(Ordering::Relaxed) {
            return Err("Cancelled. Previous results were retained.".into());
        }
        let mut best = (f64::NEG_INFINITY, usize::MAX);
        for i in 0..n {
            if chosen[i] {
                continue;
            }
            let d = distance(&descriptors[i], &descriptors[current]);
            if d < nearest[i] {
                nearest[i] = d;
            }
            if nearest[i] > best.0 {
                best = (nearest[i], i);
            }
        }
        current = best.1;
    }
}
