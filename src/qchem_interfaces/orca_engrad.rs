//! Reading ORCA's `.engrad` file: the Cartesian gradient at one geometry.
//!
//! The gradient is not in the `.hess`. That file has twelve sections and none
//! of them is a gradient -- ORCA's only `$gradients` keyword belongs to its MD
//! trajectory format, not to the Hessian file. So a reaction-path projection,
//! which needs the gradient as well as the Hessian, needs a second file.
//!
//! `.engrad` is what ORCA writes for an `EnGrad` run, and its layout is fixed:
//! `#`-delimited comment banners, then the atom count, the energy, `3N`
//! gradient components in Eh/bohr, and finally the geometry as atomic number
//! and position in bohr.
//!
//! ```text
//! #
//! # Number of atoms
//! #
//!  3
//! #
//! # The current total energy in Eh
//! #
//!     -76.372268232289
//! #
//! # The current gradient in Eh/bohr
//! #
//!       -0.000053014439
//!        ...
//! #
//! # The atomic numbers and current coordinates in Bohr
//! #
//!    8    -0.0488064    0.0000000   -0.0344820
//! ```
//!
//! The banners are read rather than the blank lines counted, because the
//! comment text is what ORCA's own reader keys on and it is more stable than
//! a line offset.

use crate::molecule::element_symbol;

/// How far apart two geometries may be and still be taken as the same one.
///
/// The `.engrad` prints coordinates to seven decimals where the `.hess` prints
/// twelve, so an honest pair differs by rounding alone -- 6.7e-8 bohr on the
/// 43-atom example. Anything above this is a different structure, which is the
/// mistake worth catching: a gradient from one geometry with a Hessian from
/// another produces frequencies that look plausible and are meaningless.
pub const GEOMETRY_TOLERANCE_BOHR: f64 = 1.0e-4;

#[derive(Debug, Clone, PartialEq)]
pub struct OrcaGradient {
    pub energy_hartree: f64,
    /// `3N` components in Eh/bohr, x, y, z per atom.
    pub gradient_bohr: Vec<f64>,
    /// Element symbols, converted from the atomic numbers the file gives.
    pub atoms: Vec<String>,
    /// `3N` coordinates in bohr, for checking this really is the geometry the
    /// Hessian was computed at.
    pub coords_bohr: Vec<f64>,
}

impl OrcaGradient {
    pub fn natoms(&self) -> usize {
        self.atoms.len()
    }

    /// The gradient's Euclidean norm, in Eh/bohr.
    ///
    /// What decides whether a reaction-path projection means anything: at a
    /// stationary point this is essentially zero and there is no path
    /// direction to remove.
    pub fn norm(&self) -> f64 {
        self.gradient_bohr.iter().map(|g| g * g).sum::<f64>().sqrt()
    }

    /// Checks this gradient belongs to the geometry the Hessian was computed
    /// at, by element and by position.
    pub fn check_matches(&self, atoms: &[String], coords_bohr: &[f64]) -> Result<(), String> {
        if self.atoms.len() != atoms.len() {
            return Err(format!(
                "The gradient file has {} atoms but the Hessian has {}.",
                self.atoms.len(),
                atoms.len()
            ));
        }
        for (i, (mine, theirs)) in self.atoms.iter().zip(atoms).enumerate() {
            if !mine.eq_ignore_ascii_case(theirs) {
                return Err(format!(
                    "Atom {} is {theirs} in the Hessian but {mine} in the gradient file.",
                    i + 1
                ));
            }
        }
        if self.coords_bohr.len() != coords_bohr.len() {
            return Err("The two files disagree about the number of coordinates.".to_string());
        }
        let mut worst = 0.0_f64;
        let mut worst_atom = 0usize;
        for i in 0..self.atoms.len() {
            let d2: f64 = (0..3)
                .map(|k| {
                    let d = self.coords_bohr[3 * i + k] - coords_bohr[3 * i + k];
                    d * d
                })
                .sum();
            let d = d2.sqrt();
            if d > worst {
                worst = d;
                worst_atom = i;
            }
        }
        if worst > GEOMETRY_TOLERANCE_BOHR {
            return Err(format!(
                "The gradient was computed at a different geometry: atom {} ({}) differs by \
                 {worst:.3e} bohr. A gradient and a Hessian from different structures give \
                 meaningless frequencies, so this pair is refused.",
                worst_atom + 1,
                self.atoms[worst_atom]
            ));
        }
        Ok(())
    }
}

/// The line index just after the banner whose text contains `needle`.
///
/// ORCA writes each field under a three-line `#` banner, so the value starts
/// on the line after the banner's closing `#`.
fn value_start(lines: &[&str], needle: &str) -> Option<usize> {
    let banner = lines
        .iter()
        .position(|l| l.trim_start().starts_with('#') && l.to_lowercase().contains(needle))?;
    // Skip the banner line itself and the closing `#` beneath it.
    let mut i = banner + 1;
    while i < lines.len() && lines[i].trim_start().starts_with('#') {
        i += 1;
    }
    Some(i)
}

/// Reads the numbers following a banner, stopping at the next banner.
fn numbers_after<'a>(lines: &[&'a str], start: usize) -> Vec<&'a str> {
    let mut out = Vec::new();
    for line in &lines[start..] {
        let trimmed = line.trim();
        if trimmed.starts_with('#') {
            break;
        }
        if trimmed.is_empty() {
            continue;
        }
        out.push(*line);
    }
    out
}

pub fn parse_orca_engrad(text: &str) -> Result<OrcaGradient, String> {
    let lines: Vec<&str> = text.lines().collect();

    let natoms_at = value_start(&lines, "number of atoms")
        .ok_or_else(|| "Not an ORCA .engrad file: no atom-count section.".to_string())?;
    let natoms: usize = numbers_after(&lines, natoms_at)
        .first()
        .and_then(|l| l.trim().split_whitespace().next())
        .and_then(|t| t.parse().ok())
        .ok_or_else(|| "Could not read the atom count.".to_string())?;
    if natoms == 0 {
        return Err("The gradient file reports zero atoms.".to_string());
    }

    let energy_at = value_start(&lines, "total energy")
        .ok_or_else(|| "No total-energy section in the .engrad file.".to_string())?;
    let energy_hartree: f64 = numbers_after(&lines, energy_at)
        .first()
        .and_then(|l| l.trim().split_whitespace().next())
        .and_then(|t| t.parse().ok())
        .ok_or_else(|| "Could not read the total energy.".to_string())?;

    let grad_at = value_start(&lines, "current gradient")
        .ok_or_else(|| "No gradient section in the .engrad file.".to_string())?;
    let mut gradient_bohr = Vec::with_capacity(3 * natoms);
    for line in numbers_after(&lines, grad_at) {
        for token in line.split_whitespace() {
            gradient_bohr.push(
                token
                    .parse::<f64>()
                    .map_err(|_| format!("Could not parse the gradient value {token:?}."))?,
            );
        }
    }
    if gradient_bohr.len() != 3 * natoms {
        return Err(format!(
            "The gradient has {} values; expected {} for {natoms} atoms.",
            gradient_bohr.len(),
            3 * natoms
        ));
    }

    // The geometry block is optional in principle; without it the pair cannot
    // be checked, which is worth refusing rather than trusting.
    let geom_at = value_start(&lines, "atomic numbers").ok_or_else(|| {
        "The .engrad file carries no geometry, so it cannot be checked against the Hessian."
            .to_string()
    })?;
    let mut atoms = Vec::with_capacity(natoms);
    let mut coords_bohr = Vec::with_capacity(3 * natoms);
    for line in numbers_after(&lines, geom_at).into_iter().take(natoms) {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 4 {
            return Err(format!("Malformed geometry line in the .engrad file: {line:?}"));
        }
        let z: u32 = f[0]
            .parse()
            .map_err(|_| format!("Could not parse the atomic number {:?}.", f[0]))?;
        atoms.push(
            element_symbol(z)
                .ok_or_else(|| format!("Unknown atomic number {z} in the gradient file."))?
                .to_string(),
        );
        for k in 1..4 {
            coords_bohr.push(
                f[k].parse::<f64>()
                    .map_err(|_| format!("Could not parse the coordinate {:?}.", f[k]))?,
            );
        }
    }
    if atoms.len() != natoms {
        return Err(format!(
            "The geometry lists {} atoms; expected {natoms}.",
            atoms.len()
        ));
    }

    Ok(OrcaGradient {
        energy_hartree,
        gradient_bohr,
        atoms,
        coords_bohr,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn example(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join(name)
    }

    fn real() -> OrcaGradient {
        let text =
            std::fs::read_to_string(example("Exc_root_OnlyFreq_from_optimized.engrad")).unwrap();
        parse_orca_engrad(&text).expect("the example .engrad should parse")
    }

    /// The real 43-atom pair: an excited-state frequency job run at the
    /// ground-state geometry, so the gradient is genuinely non-zero and the
    /// reaction-path projection means something.
    #[test]
    fn the_example_gradient_is_read_whole() {
        let g = real();
        assert_eq!(g.natoms(), 43);
        assert_eq!(g.gradient_bohr.len(), 129);
        assert_eq!(g.coords_bohr.len(), 129);
        assert!((g.energy_hartree - -2299.235455890749).abs() < 1e-9);
        assert!((g.gradient_bohr[0] - -0.002827772183).abs() < 1e-12);
        assert!((g.gradient_bohr[1] - 0.001385281473).abs() < 1e-12);
        assert_eq!(g.atoms[0], "Fe", "atomic number 26 is iron");
        assert_eq!(g.atoms[1], "O");
    }

    /// Well above the floor below which a reaction path has no direction.
    #[test]
    fn the_example_is_not_a_stationary_point() {
        let norm = real().norm();
        assert!((norm - 5.6828e-3).abs() < 1e-6, "norm was {norm:e}");
        assert!(norm > 1.0e-4, "otherwise there is no path to project out");
    }

    /// The pair has to describe one geometry, and this one does -- to within
    /// the .engrad's seven printed decimals.
    #[test]
    fn the_example_pair_agrees_on_the_geometry() {
        let hess_text =
            std::fs::read_to_string(example("Exc_root_OnlyFreq_from_optimized.hess")).unwrap();
        let hess = super::super::orca_hess::parse_orca_hess(&hess_text).unwrap();
        real()
            .check_matches(&hess.atoms, &hess.coords_bohr)
            .expect("the example pair is one geometry");
    }

    /// The mistake worth catching: a gradient from a different structure gives
    /// frequencies that look plausible and mean nothing.
    #[test]
    fn a_gradient_from_another_geometry_is_refused() {
        let g = real();
        let mut moved = g.coords_bohr.clone();
        moved[0] += 0.05;
        let err = g.check_matches(&g.atoms, &moved).expect_err("different geometry");
        assert!(err.contains("different geometry"), "{err}");
        assert!(err.contains("atom 1"), "the message names the atom: {err}");
    }

    #[test]
    fn a_mismatched_atom_count_is_refused() {
        let g = real();
        let err = g
            .check_matches(&g.atoms[..5].to_vec(), &g.coords_bohr[..15])
            .expect_err("wrong atom count");
        assert!(err.contains("43") && err.contains('5'), "{err}");
    }

    /// Same count, different elements -- a gradient for the wrong molecule.
    #[test]
    fn a_different_element_is_refused() {
        let g = real();
        let mut atoms = g.atoms.clone();
        atoms[3] = "C".to_string();
        let err = g.check_matches(&atoms, &g.coords_bohr).expect_err("wrong element");
        assert!(err.contains("Atom 4"), "{err}");
    }

    /// Rounding between the two files' precisions must not trip the check.
    #[test]
    fn seven_decimal_rounding_still_counts_as_the_same_geometry() {
        let g = real();
        let nudged: Vec<f64> = g.coords_bohr.iter().map(|c| c + 5.0e-8).collect();
        assert!(g.check_matches(&g.atoms, &nudged).is_ok());
    }

    #[test]
    fn a_file_that_is_not_an_engrad_says_so() {
        let err = parse_orca_engrad("just some text\n").expect_err("not an engrad");
        assert!(err.contains("Not an ORCA .engrad"), "{err}");
    }

    /// A truncated gradient is an error, not a silently short vector.
    #[test]
    fn a_truncated_gradient_is_refused() {
        let text = "#\n# Number of atoms\n#\n 3\n#\n# The current total energy in Eh\n#\n -76.0\n\
                    #\n# The current gradient in Eh/bohr\n#\n 0.1\n 0.2\n\
                    #\n# The atomic numbers and current coordinates in Bohr\n#\n\
                    8 0.0 0.0 0.0\n1 0.0 0.0 1.8\n1 1.8 0.0 0.0\n";
        let err = parse_orca_engrad(text).expect_err("short gradient");
        assert!(err.contains("expected 9"), "{err}");
    }

    /// Without a geometry the pair cannot be checked, which is worth refusing
    /// rather than trusting.
    #[test]
    fn an_engrad_with_no_geometry_is_refused() {
        let text = "#\n# Number of atoms\n#\n 1\n#\n# The current total energy in Eh\n#\n -1.0\n\
                    #\n# The current gradient in Eh/bohr\n#\n 0.1\n 0.2\n 0.3\n";
        let err = parse_orca_engrad(text).expect_err("no geometry");
        assert!(err.contains("cannot be checked"), "{err}");
    }
}
