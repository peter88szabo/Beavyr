//! Behemoth's `utils::atomic_masses::AtomicMasses` surface, backed by Beavyr's own mass table.
//!
//! The imported optimizer code asks for masses through this type. Rather than carry a second
//! table -- which would be one more place for the two projects to disagree about an isotope --
//! this forwards to `molecule::atomic_mass_amu`.

use crate::molecule::atomic_mass_amu;

/// Standard atomic masses, in amu.
#[derive(Debug, Default, Clone)]
pub struct AtomicMasses;

impl AtomicMasses {
    pub fn new() -> Self {
        Self
    }

    /// Mass of `element` in amu.
    ///
    /// Falls back to carbon's mass for an unrecognised symbol, matching Behemoth's behaviour of
    /// returning a finite number rather than panicking deep inside an optimisation.
    pub fn get_mass(&self, element: &str) -> f64 {
        atomic_mass_amu(element).unwrap_or(12.0)
    }

    pub fn get_mass_vector(&self, elements: Vec<&str>) -> Vec<f64> {
        elements.iter().map(|e| self.get_mass(e)).collect()
    }
}
