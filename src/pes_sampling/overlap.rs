//! Pebble's atom-pair overlap cutoffs, in Angstrom rather than bohr.
#[derive(Clone)]
pub struct OverlapFilter {
    pub enabled: bool,
    pub hydrogen_hydrogen: f64,
    pub hydrogen_other: f64,
    pub other_other: f64,
}

impl Default for OverlapFilter {
    fn default() -> Self {
        Self {
            enabled: true,
            // Pebble/src/atomic_overlap.py: H_H, H_X and X_X thresholds.
            hydrogen_hydrogen: 0.4,
            hydrogen_other: 0.5,
            other_other: 0.8,
        }
    }
}

impl OverlapFilter {
    pub fn validate(&self) -> Result<(), String> {
        if [
            self.hydrogen_hydrogen,
            self.hydrogen_other,
            self.other_other,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.)
        {
            return Err("Overlap cutoffs must be finite and nonnegative (Å).".into());
        }
        Ok(())
    }

    /// Check all pairs, including pairs spanning two collision partners.
    pub fn accepts(&self, atoms: &[String], coordinates: &[f64]) -> bool {
        if !self.enabled {
            return true;
        }
        for i in 0..atoms.len() {
            for j in 0..i {
                let cutoff = match (atoms[i] == "H", atoms[j] == "H") {
                    (true, true) => self.hydrogen_hydrogen,
                    (true, false) | (false, true) => self.hydrogen_other,
                    (false, false) => self.other_other,
                };
                let distance_squared: f64 = (0..3)
                    .map(|k| (coordinates[3 * i + k] - coordinates[3 * j + k]).powi(2))
                    .sum();
                // Pebble rejects distances strictly below the threshold.
                if distance_squared < cutoff * cutoff {
                    return false;
                }
            }
        }
        true
    }
}
