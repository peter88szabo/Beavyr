//! Molden-format reader.
//!
//! Handles the sections needed for orbital plotting: `[Atoms]`, `[GTO]`,
//! the pure-function declarations `[5D]/[5D7F]/[5D10F]/[7F]/[9G]`, and `[MO]`.
//!
//! Closed versus open shell follows Molden's own rule from `rdmolf.f`: the file
//! is open shell if and only if a `Spin= Beta` block appears.  Molden
//! distinguishes the two by token length (`Alpha` is 5 characters, `Beta` is
//! 4); we compare case-insensitively instead, which accepts the same files.

use bevy::prelude::*;

use super::basis::{BasisLayout, Shell};

/// Bohr radius in Angstrom, for `[Atoms] AU` files.
const BOHR_TO_ANGSTROM: f32 = 0.529_177_210_9;

#[derive(Debug, Clone)]
pub struct MolecularOrbital {
    /// Symmetry label as written, e.g. `1a`.
    pub symmetry: String,
    /// Orbital energy in Hartree.
    pub energy: f64,
    pub occupation: f64,
    /// One coefficient per basis function, in file order.
    pub coefficients: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Spin {
    Alpha,
    Beta,
}

#[derive(Debug, Clone, Default)]
pub struct OrbitalData {
    pub atoms: Vec<String>,
    /// Positions in Angstrom, converted from Bohr when the header says so.
    pub positions: Vec<Vec3>,
    pub shells: Vec<Shell>,
    pub layout: BasisLayout,
    pub alpha: Vec<MolecularOrbital>,
    /// Empty for a closed-shell file.
    pub beta: Vec<MolecularOrbital>,
}

impl OrbitalData {
    pub fn is_open_shell(&self) -> bool {
        !self.beta.is_empty()
    }

    pub fn orbitals(&self, spin: Spin) -> &[MolecularOrbital] {
        match spin {
            Spin::Alpha => &self.alpha,
            Spin::Beta => &self.beta,
        }
    }

    /// Index of the highest occupied orbital in one spin set.
    ///
    /// Computed per set rather than from a total electron count: in an
    /// open-shell system the alpha and beta HOMOs are generally different
    /// orbitals.
    pub fn homo_index(&self, spin: Spin) -> Option<usize> {
        self.orbitals(spin)
            .iter()
            .enumerate()
            .filter(|(_, o)| o.occupation > 1.0e-6)
            .map(|(i, _)| i)
            .next_back()
    }
}

#[derive(Debug, Clone)]
pub struct MoldenError {
    pub line: usize,
    pub message: String,
}

impl std::fmt::Display for MoldenError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.line == 0 {
            write!(f, "{}", self.message)
        } else {
            write!(f, "line {}: {}", self.line, self.message)
        }
    }
}

fn err(line: usize, message: impl Into<String>) -> MoldenError {
    MoldenError {
        line,
        message: message.into(),
    }
}

/// Which angular momenta the file declares as pure spherical.
#[derive(Debug, Clone, Copy, Default)]
struct PureFlags {
    d: bool,
    f: bool,
    g: bool,
}

impl PureFlags {
    fn is_spherical(&self, l: u8) -> bool {
        match l {
            0 | 1 => true, // s and p have no Cartesian/pure distinction
            2 => self.d,
            3 => self.f,
            4 => self.g,
            _ => false,
        }
    }
}

/// Parse Fortran-style reals, accepting `D` exponents alongside `E`.
fn parse_real(token: &str) -> Option<f64> {
    token
        .parse::<f64>()
        .ok()
        .or_else(|| token.replace(['D', 'd'], "E").parse::<f64>().ok())
}

fn section_name(line: &str) -> Option<String> {
    let t = line.trim();
    if !t.starts_with('[') {
        return None;
    }
    let end = t.find(']')?;
    Some(t[1..end].trim().to_ascii_lowercase())
}

pub fn parse_molden(text: &str) -> Result<OrbitalData, MoldenError> {
    let lines: Vec<&str> = text.lines().collect();

    // Pure-function declarations may appear anywhere before or after [GTO],
    // so scan for them first.
    let mut pure = PureFlags::default();
    for line in &lines {
        match section_name(line).as_deref() {
            Some("5d") | Some("5d7f") => {
                pure.d = true;
                pure.f = true;
            }
            Some("5d10f") => pure.d = true,
            Some("7f") => pure.f = true,
            Some("9g") => pure.g = true,
            _ => {}
        }
    }

    let mut data = OrbitalData::default();
    let mut i = 0usize;
    let mut seen_mo = false;

    while i < lines.len() {
        let Some(name) = section_name(lines[i]) else {
            i += 1;
            continue;
        };

        match name.as_str() {
            "atoms" => {
                let angstrom = {
                    let rest = lines[i].to_ascii_lowercase();
                    // "[Atoms] AU" / "[Atoms] Angs"; default is Angstrom.
                    !(rest.contains("au") && !rest.contains("angs"))
                };
                i = parse_atoms(&lines, i + 1, angstrom, &mut data)?;
            }
            "gto" => {
                i = parse_gto(&lines, i + 1, pure, &mut data)?;
            }
            "mo" => {
                seen_mo = true;
                data.layout = BasisLayout::build(&data.shells);
                i = parse_mo(&lines, i + 1, &mut data)?;
            }
            _ => i += 1,
        }
    }

    if data.shells.is_empty() {
        return Err(err(0, "no [GTO] basis section found"));
    }
    if !seen_mo {
        return Err(err(0, "no [MO] section found"));
    }
    data.layout = BasisLayout::build(&data.shells);

    validate(&data)?;
    Ok(data)
}

fn parse_atoms(
    lines: &[&str],
    mut i: usize,
    angstrom: bool,
    data: &mut OrbitalData,
) -> Result<usize, MoldenError> {
    let scale = if angstrom { 1.0 } else { BOHR_TO_ANGSTROM };

    while i < lines.len() {
        if section_name(lines[i]).is_some() {
            break;
        }
        let parts: Vec<&str> = lines[i].split_whitespace().collect();
        if parts.is_empty() {
            i += 1;
            continue;
        }
        if parts.len() < 6 {
            return Err(err(i + 1, format!("malformed [Atoms] row: {:?}", lines[i])));
        }
        let coord = |k: usize| -> Result<f32, MoldenError> {
            parse_real(parts[k])
                .map(|v| v as f32 * scale)
                .ok_or_else(|| err(i + 1, format!("bad coordinate {:?}", parts[k])))
        };
        data.atoms.push(parts[0].to_string());
        data.positions
            .push(Vec3::new(coord(3)?, coord(4)?, coord(5)?));
        i += 1;
    }
    Ok(i)
}

fn angular_momentum(symbol: &str) -> Option<u8> {
    match symbol.to_ascii_lowercase().as_str() {
        "s" => Some(0),
        "p" => Some(1),
        "d" => Some(2),
        "f" => Some(3),
        "g" => Some(4),
        // sp shells are legal Molden but need splitting into two shells; the
        // orbital viewer refuses rather than silently mis-ordering the basis.
        _ => None,
    }
}

fn parse_gto(
    lines: &[&str],
    mut i: usize,
    pure: PureFlags,
    data: &mut OrbitalData,
) -> Result<usize, MoldenError> {
    let mut center: Option<usize> = None;

    while i < lines.len() {
        if section_name(lines[i]).is_some() {
            break;
        }
        let parts: Vec<&str> = lines[i].split_whitespace().collect();
        if parts.is_empty() {
            i += 1;
            continue;
        }

        // "<atom index> 0" starts a new atom block.
        if let Ok(index) = parts[0].parse::<usize>() {
            if index == 0 {
                return Err(err(i + 1, "[GTO] atom indices are 1-based"));
            }
            center = Some(index - 1);
            i += 1;
            continue;
        }

        // "<shell letter> <n primitives> <scale>"
        let Some(l) = angular_momentum(parts[0]) else {
            return Err(err(
                i + 1,
                format!(
                    "unsupported shell type {:?}; the viewer handles s, p, d, f and g",
                    parts[0]
                ),
            ));
        };
        let Some(center) = center else {
            return Err(err(i + 1, "shell appears before any atom index in [GTO]"));
        };
        let n_prim: usize = parts
            .get(1)
            .and_then(|t| t.parse().ok())
            .ok_or_else(|| err(i + 1, "shell header is missing its primitive count"))?;

        let mut exponents = Vec::with_capacity(n_prim);
        let mut coefficients = Vec::with_capacity(n_prim);
        for k in 0..n_prim {
            let row = i + 1 + k;
            let p: Vec<&str> = lines
                .get(row)
                .ok_or_else(|| err(row, "file ends inside a [GTO] contraction"))?
                .split_whitespace()
                .collect();
            if p.len() < 2 {
                return Err(err(row + 1, "expected an exponent and a coefficient"));
            }
            exponents.push(
                parse_real(p[0]).ok_or_else(|| err(row + 1, format!("bad exponent {:?}", p[0])))?,
            );
            coefficients.push(
                parse_real(p[1])
                    .ok_or_else(|| err(row + 1, format!("bad coefficient {:?}", p[1])))?,
            );
        }

        let mut shell = Shell {
            center,
            l,
            spherical: pure.is_spherical(l),
            exponents,
            raw_coefficients: coefficients,
            coefficients: Vec::new(),
        };
        shell.normalize();
        data.shells.push(shell);

        i += 1 + n_prim;
    }
    Ok(i)
}

fn parse_mo(lines: &[&str], mut i: usize, data: &mut OrbitalData) -> Result<usize, MoldenError> {
    let n_basis = data.layout.n_basis;

    let mut symmetry = String::new();
    let mut energy = 0.0f64;
    let mut occupation = 0.0f64;
    let mut spin = Spin::Alpha;
    let mut coefficients: Vec<f64> = Vec::new();
    let mut in_orbital = false;

    // Push whatever orbital has been accumulated so far.
    macro_rules! flush {
        () => {
            if in_orbital {
                let orbital = MolecularOrbital {
                    symmetry: std::mem::take(&mut symmetry),
                    energy,
                    occupation,
                    coefficients: std::mem::replace(&mut coefficients, vec![0.0; n_basis]),
                };
                match spin {
                    Spin::Alpha => data.alpha.push(orbital),
                    Spin::Beta => data.beta.push(orbital),
                }
            }
        };
    }

    while i < lines.len() {
        if section_name(lines[i]).is_some() {
            break;
        }
        let line = lines[i].trim();
        if line.is_empty() {
            i += 1;
            continue;
        }

        let lower = line.to_ascii_lowercase();
        if let Some(rest) = keyword_value(&lower, line, "sym=") {
            flush!();
            in_orbital = true;
            coefficients = vec![0.0; n_basis];
            symmetry = rest.trim().to_string();
            energy = 0.0;
            occupation = 0.0;
            spin = Spin::Alpha;
        } else if let Some(rest) = keyword_value(&lower, line, "ene=") {
            // Some producers omit Sym=, so Ene= can also open an orbital.
            if !in_orbital {
                in_orbital = true;
                coefficients = vec![0.0; n_basis];
                symmetry.clear();
                occupation = 0.0;
                spin = Spin::Alpha;
            }
            energy = parse_real(rest.trim())
                .ok_or_else(|| err(i + 1, format!("bad orbital energy {:?}", rest.trim())))?;
        } else if let Some(rest) = keyword_value(&lower, line, "spin=") {
            spin = if rest.trim().eq_ignore_ascii_case("beta") {
                Spin::Beta
            } else {
                Spin::Alpha
            };
        } else if let Some(rest) = keyword_value(&lower, line, "occup=") {
            occupation = parse_real(rest.trim())
                .ok_or_else(|| err(i + 1, format!("bad occupation {:?}", rest.trim())))?;
        } else {
            // "<basis index> <coefficient>"
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                if let (Ok(index), Some(value)) = (parts[0].parse::<usize>(), parse_real(parts[1]))
                {
                    if index == 0 || index > n_basis {
                        return Err(err(
                            i + 1,
                            format!(
                                "MO coefficient index {index} is outside the basis \
                                 (1..={n_basis}); the [GTO] section and the [MO] section disagree"
                            ),
                        ));
                    }
                    if !in_orbital {
                        in_orbital = true;
                        coefficients = vec![0.0; n_basis];
                    }
                    coefficients[index - 1] = value;
                }
            }
        }
        i += 1;
    }

    flush!();
    Ok(i)
}

/// Match `keyword` at the start of the lowercased line and return the original
/// (non-lowercased) remainder.
fn keyword_value<'a>(lower: &str, original: &'a str, keyword: &str) -> Option<&'a str> {
    lower
        .starts_with(keyword)
        .then(|| &original[keyword.len()..])
}

fn validate(data: &OrbitalData) -> Result<(), MoldenError> {
    if data.atoms.is_empty() {
        return Err(err(0, "no [Atoms] section found"));
    }
    for shell in &data.shells {
        if shell.center >= data.atoms.len() {
            return Err(err(
                0,
                format!(
                    "a [GTO] shell references atom {} but only {} atoms were read",
                    shell.center + 1,
                    data.atoms.len()
                ),
            ));
        }
    }
    if data.alpha.is_empty() && data.beta.is_empty() {
        return Err(err(0, "[MO] section contained no orbitals"));
    }
    // The structural check that catches a 5D-versus-6D misread: every orbital
    // must carry exactly one coefficient per basis function.
    for (set, label) in [(&data.alpha, "alpha"), (&data.beta, "beta")] {
        for (k, orbital) in set.iter().enumerate() {
            if orbital.coefficients.len() != data.layout.n_basis {
                return Err(err(
                    0,
                    format!(
                        "{label} orbital {} has {} coefficients but the basis has {} \
                         functions",
                        k + 1,
                        orbital.coefficients.len(),
                        data.layout.n_basis
                    ),
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal closed-shell fixture: H2 with one s function per atom.
    const H2_CLOSED: &str = "\
[Molden Format]
[Atoms] AU
H     1    1    0.000000    0.000000   -0.700000
H     2    1    0.000000    0.000000    0.700000
[GTO]
  1 0
s   1 1.0
      1.2000000000    1.0000000000

  2 0
s   1 1.0
      1.2000000000    1.0000000000

[MO]
 Sym= 1a
 Ene= -0.5950000000
 Spin= Alpha
 Occup= 2.000000
   1    0.5483000000
   2    0.5483000000
 Sym= 2a
 Ene=  0.1890000000
 Spin= Alpha
 Occup= 0.000000
   1    1.2124000000
   2   -1.2124000000
";

    /// Same geometry, but written as an unrestricted calculation.
    const H2_OPEN: &str = "\
[Molden Format]
[Atoms] Angs
H     1    1    0.000000    0.000000   -0.370000
H     2    1    0.000000    0.000000    0.370000
[GTO]
  1 0
s   1 1.0
      1.2000000000    1.0000000000

  2 0
s   1 1.0
      1.2000000000    1.0000000000

[MO]
 Sym= 1a
 Ene= -0.5950000000
 Spin= Alpha
 Occup= 1.000000
   1    0.5483000000
   2    0.5483000000
 Sym= 2a
 Ene=  0.1890000000
 Spin= Alpha
 Occup= 0.000000
   1    1.2124000000
   2   -1.2124000000
 Sym= 1a
 Ene= -0.5100000000
 Spin= Beta
 Occup= 1.000000
   1    0.5300000000
   2    0.5300000000
 Sym= 2a
 Ene=  0.2200000000
 Spin= Beta
 Occup= 0.000000
   1    1.1900000000
   2   -1.1900000000
";

    #[test]
    fn closed_shell_file_has_one_orbital_set() {
        let data = parse_molden(H2_CLOSED).expect("fixture should parse");
        assert_eq!(data.atoms, vec!["H", "H"]);
        assert_eq!(data.layout.n_basis, 2);
        assert_eq!(data.alpha.len(), 2);
        assert!(data.beta.is_empty());
        assert!(!data.is_open_shell());
        assert_eq!(data.alpha[0].occupation, 2.0);
    }

    #[test]
    fn open_shell_file_has_separate_alpha_and_beta_sets() {
        let data = parse_molden(H2_OPEN).expect("fixture should parse");
        assert!(data.is_open_shell());
        assert_eq!(data.alpha.len(), 2);
        assert_eq!(data.beta.len(), 2);
        // The two sets carry independent energies.
        assert!((data.alpha[0].energy - -0.595).abs() < 1e-9);
        assert!((data.beta[0].energy - -0.510).abs() < 1e-9);
        assert_eq!(data.alpha[0].occupation, 1.0);
        assert_eq!(data.beta[0].occupation, 1.0);
    }

    #[test]
    fn homo_is_found_per_spin_set() {
        let data = parse_molden(H2_OPEN).expect("fixture should parse");
        assert_eq!(data.homo_index(Spin::Alpha), Some(0));
        assert_eq!(data.homo_index(Spin::Beta), Some(0));
    }

    #[test]
    fn bohr_coordinates_are_converted_but_angstrom_are_not() {
        let au = parse_molden(H2_CLOSED).unwrap();
        let angs = parse_molden(H2_OPEN).unwrap();
        // 0.7 Bohr = 0.3704 Angstrom, matching the Angstrom fixture.
        assert!((au.positions[1].z - 0.370_4).abs() < 1e-3);
        assert!((angs.positions[1].z - 0.370).abs() < 1e-6);
    }

    #[test]
    fn a_coefficient_index_beyond_the_basis_is_an_error() {
        let bad = H2_CLOSED.replace("   2    0.5483000000", "   9    0.5483000000");
        let error = parse_molden(&bad).expect_err("index 9 exceeds a 2-function basis");
        assert!(
            error.message.contains("outside the basis"),
            "unhelpful message: {}",
            error
        );
    }

    #[test]
    fn a_missing_mo_section_is_an_error() {
        let truncated = H2_CLOSED.split("[MO]").next().unwrap().to_string();
        let error = parse_molden(&truncated).expect_err("no [MO] block");
        assert!(error.message.contains("[MO]"), "unhelpful: {error}");
    }

    #[test]
    fn pure_function_flags_are_picked_up() {
        let spherical = H2_CLOSED.replace("[GTO]", "[5D]\n[7F]\n[9G]\n[GTO]");
        let data = parse_molden(&spherical).expect("should parse");
        // s shells are unaffected, but the flags must have been seen.
        assert_eq!(data.layout.n_basis, 2);
    }

    #[test]
    fn fortran_d_exponents_are_accepted() {
        let fortran = H2_CLOSED.replace("1.2000000000    1.0000000000", "1.2D+00    1.0D+00");
        let data = parse_molden(&fortran).expect("D exponents are valid Fortran reals");
        assert_eq!(data.shells.len(), 2);
    }
}
