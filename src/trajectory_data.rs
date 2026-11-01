//! Optional frame provenance and energies; absent energies remain absent.
use crate::trajectory::TrajectoryFrame;
use serde::{Deserialize, Serialize};
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Metadata {
    pub comment: String,
    pub energy_hartree: Option<f64>,
    #[serde(default)]
    pub energy_is_relative: bool,
}
impl Metadata {
    pub fn comment(comment: String) -> Self {
        Self {
            energy_hartree: parse_energy(&comment),
            energy_is_relative: comment.to_ascii_lowercase().contains("relative_energy"),
            comment,
        }
    }
    pub fn energy(energy: f64, relative: bool, comment: String) -> Self {
        Self {
            comment,
            energy_hartree: energy.is_finite().then_some(energy),
            energy_is_relative: relative,
        }
    }
}
pub fn parse_energy(comment: &str) -> Option<f64> {
    let lower = comment.to_ascii_lowercase();
    let normalized = lower.replace(['=', ':', ','], " ");
    let words: Vec<_> = normalized.split_whitespace().collect();
    let at = words.iter().position(|w| {
        matches!(
            *w,
            "e" | "energy" | "relative_energy" | "relative_energy_kcal"
        )
    })?;
    let value = words.get(at + 1)?.replace('d', "e").parse::<f64>().ok()?;
    if !value.is_finite() {
        return None;
    }
    let unit = words.get(at + 2).copied();
    let factor = if words[at] == "relative_energy_kcal" {
        1. / 627.5094740631
    } else {
        match unit {
            Some("ev") => 1. / 27.211386245988,
            Some("kj/mol") => 1. / 2625.499639479163,
            Some("kcal/mol") => 1. / 627.5094740631,
            Some("eh" | "hartree" | "hartrees" | "au" | "a.u.") => 1.,
            // xTB's documented trajectory comments use energy: <Hartree> gnorm: …
            _ if lower.trim_start().starts_with("energy:") => 1.,
            _ => return None,
        }
    };
    Some(value * factor)
}
pub fn fingerprint(frames: &[TrajectoryFrame]) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::collections::hash_map::DefaultHasher::new();
    frames.len().hash(&mut h);
    for f in frames {
        f.atoms.hash(&mut h);
        for p in &f.pos {
            for c in p.to_array() {
                c.to_bits().hash(&mut h);
            }
        }
        f.metadata.comment.hash(&mut h);
        f.metadata.energy_hartree.map(f64::to_bits).hash(&mut h);
        f.metadata.energy_is_relative.hash(&mut h);
    }
    h.finish()
}
pub fn xyz(frames: &[TrajectoryFrame], indices: impl IntoIterator<Item = usize>) -> String {
    use std::fmt::Write;
    let mut out = String::new();
    for i in indices {
        let f = &frames[i];
        let mut comment = f.metadata.comment.replace(['\n', '\r'], " ");
        if comment.is_empty() {
            comment = format!("Frame {}", i + 1);
            if let Some(e) = f.metadata.energy_hartree {
                comment.push_str(&format!(
                    " {}= {e:.14} Eh",
                    if f.metadata.energy_is_relative {
                        "Relative_energy"
                    } else {
                        "E"
                    }
                ));
            }
        }
        writeln!(out, "{}\n{comment}", f.atoms.len()).unwrap();
        for (a, p) in f.atoms.iter().zip(&f.pos) {
            writeln!(out, "{a:<3} {:16.8} {:16.8} {:16.8}", p.x, p.y, p.z).unwrap();
        }
    }
    out
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn energy_units_are_explicit_and_frame_numbers_are_not_energies() {
        assert_eq!(parse_energy("cycle 7 E = -5.25 Eh"), Some(-5.25));
        assert_eq!(parse_energy("energy: -3.4 gnorm: 0.01"), Some(-3.4));
        assert!((parse_energy("Energy = 27.211386245988 eV").unwrap() - 1.).abs() < 1e-12);
        for c in [
            "frame 42",
            "PES_sample=8 seed=12",
            "E=NaN Eh",
            "energy=9 unknown",
            "100",
        ] {
            assert_eq!(parse_energy(c), None);
        }
    }
    #[test]
    fn mixed_energy_availability_and_comments_survive_xyz_and_old_projects() {
        let text = "1\ncycle 1 E = -2 Eh\nH 0 0 0\n1\nPES_sample=2 seed=1\nH 1 0 0\n";
        let f = crate::trajectory::parse_multi_xyz(text).unwrap();
        assert_eq!(f[0].metadata.energy_hartree, Some(-2.));
        assert!(f[1].metadata.energy_hartree.is_none());
        let copy = xyz(&f, 0..2);
        assert!(copy.contains("PES_sample=2 seed=1"));
        let old: TrajectoryFrame = ron::from_str("(atoms:[\"H\"],pos:[(0.0,0.0,0.0)])").unwrap();
        assert!(old.metadata.comment.is_empty());
    }
}
