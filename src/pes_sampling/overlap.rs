//! Pebble's atom-pair overlap cutoffs, in Angstrom rather than bohr, plus an
//! element-aware covalent-radius clash rule.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct OverlapFilter {
    pub enabled: bool,
    pub hydrogen_hydrogen: f64,
    pub hydrogen_other: f64,
    pub other_other: f64,
    /// Also reject any pair closer than `covalent_factor` × (rA + rB), with
    /// covalent radii; element-aware, so heavy atoms are covered too.
    pub covalent: bool,
    pub covalent_factor: f64,
}

impl Default for OverlapFilter {
    fn default() -> Self {
        Self {
            enabled: true,
            // Pebble/src/atomic_overlap.py: H_H, H_X and X_X thresholds.
            hydrogen_hydrogen: 0.4,
            hydrogen_other: 0.5,
            other_other: 0.8,
            covalent: true,
            covalent_factor: 0.5,
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
            || !self.covalent_factor.is_finite()
            || !(0.0..=1.5).contains(&self.covalent_factor)
        {
            return Err("Overlap cutoffs must be finite and nonnegative (Å); the covalent-radius factor must lie between 0 and 1.5.".into());
        }
        Ok(())
    }

    /// Check all pairs, including pairs spanning two collision partners.
    pub fn accepts(&self, atoms: &[String], coordinates: &[f64]) -> bool {
        if !self.enabled {
            return true;
        }
        let radii: Vec<f64> = if self.covalent {
            atoms
                .iter()
                .map(|a| crate::molecule::covalent_radius_angstrom(a) as f64)
                .collect()
        } else {
            vec![]
        };
        for i in 0..atoms.len() {
            for j in 0..i {
                let mut cutoff = match (atoms[i] == "H", atoms[j] == "H") {
                    (true, true) => self.hydrogen_hydrogen,
                    (true, false) | (false, true) => self.hydrogen_other,
                    (false, false) => self.other_other,
                };
                if self.covalent {
                    cutoff = cutoff.max(self.covalent_factor * (radii[i] + radii[j]));
                }
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
