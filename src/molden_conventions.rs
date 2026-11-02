//! Reconcile primitive-normalized exports with the viewer's normalized-primitive convention.
//! ORCA also reverses the real f(|m|=3), g(|m|=3,4) phases.
//! References: Molden renorm.f / rdmolf.f; IOData's documented ORCA conventions.
//! Anonymous files are changed only when contraction norms identify this convention
//! and the resulting molecular orbitals pass an independent overlap-norm check.
use super::OrbitalData;
use crate::orbitals::{
    basis::{double_factorial_odd, Shell},
    overlap::Overlap,
};
fn primitive(exponent: f64, l: u8) -> f64 {
    (2. * exponent / std::f64::consts::PI).powf(0.75) * (4. * exponent).powf(l as f64 / 2.)
        / double_factorial_odd(l as usize).sqrt()
}
fn included_norm(shell: &Shell) -> f64 {
    let mut value = 0.;
    for (i, &a) in shell.exponents.iter().enumerate() {
        for (j, &b) in shell.exponents.iter().enumerate() {
            value += shell.raw_coefficients[i] / primitive(a, shell.l) * shell.raw_coefficients[j]
                / primitive(b, shell.l)
                * (2. * (a * b).sqrt() / (a + b)).powf(shell.l as f64 + 1.5);
        }
    }
    value
}
fn convert(data: &mut OrbitalData, orca: bool) {
    for shell in &mut data.shells {
        for (c, &a) in shell.raw_coefficients.iter_mut().zip(&shell.exponents) {
            *c /= primitive(a, shell.l);
        }
        shell.normalize();
    }
    if orca {
        for (i, shell) in data.shells.iter().enumerate() {
            if shell.spherical && shell.l >= 3 {
                let offset = data.layout.shell_offsets[i];
                for slot in 5..shell.n_functions() {
                    for o in data.alpha.iter_mut().chain(&mut data.beta) {
                        o.coefficients[offset + slot] *= -1.;
                    }
                }
            }
        }
    }
}
fn norm_error(data: &OrbitalData) -> Result<f64, String> {
    let overlap = Overlap::build(data)?;
    let n = overlap.n;
    let mut error = 0.0_f64;
    // Both spin sets and occupied/virtual orbitals: a few localized MOs alone
    // cannot distinguish all basis conventions.
    for o in data.alpha.iter().chain(&data.beta) {
        let mut norm = 0.;
        for i in 0..n {
            norm += o.coefficients[i]
                * overlap.matrix[i * n..(i + 1) * n]
                    .iter()
                    .zip(&o.coefficients)
                    .map(|(s, c)| s * c)
                    .sum::<f64>();
        }
        if !norm.is_finite() {
            return Err("Nonfinite Molden orbital norm".into());
        }
        error = error.max((norm - 1.).abs());
    }
    Ok(error)
}
pub fn apply(data: &mut OrbitalData, text: &str) -> Result<(), String> {
    let lower = text.to_ascii_lowercase();
    let orca = lower.contains("molden file created by orca_2mkl");
    let xtb = lower.contains("beavyr source: xtb");
    let sp: Vec<_> = data.shells.iter().filter(|s| s.l <= 1).collect();
    let includes = !sp.is_empty() && sp.iter().all(|s| (included_norm(s) - 1.).abs() < 2e-4);
    if xtb || (orca && includes) {
        convert(data, orca);
        return Ok(());
    }
    // Unlabelled xTB exports have no producer name. A contracted basis whose
    // printed (unnormalized) primitives are already normalized is a candidate,
    // but the MO norms must also prove that the conversion is appropriate.
    if !includes || !sp.iter().any(|s| s.exponents.len() > 1) {
        return Ok(());
    }
    let mut candidate = data.clone();
    convert(&mut candidate, false);
    if norm_error(&candidate)? < 2e-4 && norm_error(data)? > 0.01 {
        *data = candidate;
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn primitive_embedded_and_standard_files_give_same_orbitals() {
        let raw = [0.4, 0.6];
        let exp = [0.7, 0.2];
        let mut shell = Shell {
            center: 0,
            l: 0,
            spherical: true,
            exponents: exp.to_vec(),
            raw_coefficients: raw.to_vec(),
            coefficients: vec![],
        };
        shell.normalize();
        let make = |title: &str, c: &[f64]| {
            format!("[Molden Format]\n[Title]\n{title}\n[Atoms] Angs\nH 1 1 0 0 0\n[GTO]\n1 0\ns 2 1.0\n0.7 {:.15}\n0.2 {:.15}\n\n[MO]\nSym= 1\nEne= -0.5\nSpin= Alpha\nOccup= 1.0\n1 1.0\n",c[0],c[1])
        };
        let standard = super::super::parse_molden(&make("Standard", &raw)).unwrap();
        // A one-AO example cannot distinguish conventions by MO norms: both
        // alternatives normalize that single AO. Use the explicit producer
        // marker here; the multi-AO native fixture below tests anonymous files.
        let embedded = super::super::parse_molden(&make(
            "Beavyr source: xTB (primitive-normalized coefficients)",
            &shell.coefficients,
        ))
        .unwrap();
        for (a, b) in standard.shells[0]
            .coefficients
            .iter()
            .zip(&embedded.shells[0].coefficients)
        {
            assert!((a - b).abs() < 1e-10);
        }
        let ambiguous = super::super::parse_molden(&make("", &shell.coefficients)).unwrap();
        assert_eq!(ambiguous.shells[0].raw_coefficients.len(), 2);
        for (actual, printed) in ambiguous.shells[0]
            .raw_coefficients
            .iter()
            .zip(&shell.coefficients)
        {
            assert!((actual - printed).abs() < 1e-14);
        }
    }
}

#[cfg(test)]
mod native_exports {
    use super::*;
    #[test]
    fn actual_orca_and_anonymous_xtb_exports_have_normalized_orbitals() {
        for (text, count) in [
            (
                include_str!("../tests/fixtures/orca_water_embedded_primitive_norm.molden"),
                10.,
            ),
            (
                include_str!("../tests/fixtures/xtb_water_embedded_primitive_norm.molden"),
                8.,
            ),
        ] {
            let data = super::super::parse_molden(text).unwrap();
            assert!(norm_error(&data).unwrap() < 1e-5);
            let p = Overlap::build(&data).unwrap().population(&data).unwrap();
            assert!((p.iter().sum::<f64>() - count).abs() < 1e-5);
        }
    }
}
