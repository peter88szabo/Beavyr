//! Read final, method-specific MRCC results; never substitute SCF for a failed CC job.
use super::*;
use crate::uvvis::types::{ExcitedState, TddftResult, EV_NM, HARTREE_TO_EV};
fn number(s: &str) -> Option<f64> {
    s.replace(['D', 'd'], "E")
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
}
fn value(line: &str) -> Option<f64> {
    line.split_once(':')?.1.split_whitespace().find_map(number)
}
pub fn completed(text: &str) -> Result<(), String> {
    if !text.contains("Normal termination of mrcc.") {
        return Err(
            "MRCC did not terminate normally. See the retained output for the failed calculation."
                .into(),
        );
    }
    if [
        "not converged",
        "convergence not achieved",
        "convergence not reached",
        "failed to converge",
    ]
    .iter()
    .any(|s| text.to_ascii_lowercase().contains(s))
    {
        return Err("MRCC reports an unconverged calculation.".into());
    }
    Ok(())
}
pub fn parse_energy(text: &str, m: QcMethod) -> Result<f64, String> {
    completed(text)?;
    let prefixes: &[&str] = match m {
        QcMethod::Hf => &["***FINAL HARTREE-FOCK ENERGY:"],
        QcMethod::Dft => &["***FINAL KOHN-SHAM ENERGY:"],
        QcMethod::Mp2 => &["Total MP2 energy [au]:"],
        QcMethod::RiMp2 => &["DF-MP2 energy [au]:"],
        QcMethod::Ccsd => &["Total CCSD energy [au]:"],
        QcMethod::CcsdT | QcMethod::UccsdT => &["Total CCSD(T) energy [au]:"],
        QcMethod::Ccsdt => &["Total CCSDT energy [au]:"],
        QcMethod::CcsdtQ => &["Total CCSDT(Q) energy [au]:"],
        QcMethod::Ccsdtq => &["Total CCSDTQ energy [au]:"],
        QcMethod::Casscf => &["***FINAL MCSCF ENERGY:"],
        QcMethod::Nevpt2 => &["Total SC-NEVPT2 energy [au]:"],
        QcMethod::Drpa => &["DF-dRPA energy [au]:", "Total dRPA energy [au]:"],
        QcMethod::Adc2 | QcMethod::ScsAdc2 | QcMethod::SosAdc2 => &["Total MP2 energy:"],
        QcMethod::Mp2F12 => &[
            "MP2-F12 energy [au]:",
            "DF-MP2-F12 energy [au]:",
            "Total MP2-F12 energy [au]:",
        ],
        QcMethod::CcsdF12 => &[
            "CCSD-F12 energy [au]:",
            "Total CCSD-F12 energy [au]:",
            "Total CCSD(F12*) energy [au]:",
        ],
        QcMethod::CcsdTF12 => &[
            "CCSD(T)-F12 energy [au]:",
            "Total CCSD(T)-F12 energy [au]:",
            "Total CCSD(F12*)(T) energy [au]:",
        ],
        _ => return Err("Unsupported MRCC energy method.".into()),
    };
    text.lines().rev().find_map(|l| prefixes.iter().any(|p|l.split_whitespace().collect::<Vec<_>>().join(" ").starts_with(p)).then(||value(l)).flatten())
        .ok_or_else(||format!("MRCC did not report a finite final {} energy. No lower-level energy was substituted.",m.label()))
}
pub fn parse_gradient(text: &str, n: usize) -> Result<Vec<f64>, String> {
    completed(text)?;
    let tail = text
        .rsplit_once("Cartesian gradient in original orientation [au]:")
        .ok_or("MRCC returned no analytical gradient in the input orientation.")?
        .1;
    let mut g = Vec::with_capacity(3 * n);
    for (i, line) in tail
        .lines()
        .filter(|l| !l.trim().is_empty())
        .take(n)
        .enumerate()
    {
        let c: Vec<_> = line.split_whitespace().collect();
        if c.len() != 5 || c[0].parse::<usize>().ok() != Some(i + 1) {
            return Err("Incomplete or unordered MRCC Cartesian gradient.".into());
        }
        for s in &c[2..] {
            g.push(number(s).ok_or("Non-finite MRCC gradient component")?)
        }
    }
    if g.len() != 3 * n {
        return Err("MRCC Cartesian gradient has the wrong atom count.".into());
    }
    Ok(g)
}
pub fn parse_spectrum(text: &str, source: &Path) -> Result<TddftResult, String> {
    completed(text)?;
    let lines: Vec<_> = text.lines().collect();
    // The last final energy table replaces preliminary CIS/TDHF roots in ADC output.
    let begin = lines
        .iter()
        .rposition(|l| l.contains("Final results of") && l.contains("excited states"))
        .ok_or("MRCC returned no final excitation-energy table.")?;
    let is_adc = lines[begin].contains("ADC(2)");
    let mut states = Vec::new();
    for l in lines.iter().skip(begin + 1) {
        let c: Vec<_> = l.split_whitespace().collect();
        if c.first().and_then(|s| s.parse::<usize>().ok()).is_none() {
            if !states.is_empty() {
                break;
            } else {
                continue;
            }
        }
        let root = c[0].parse::<usize>().unwrap();
        let offset = if is_adc { 1 } else { 2 };
        if c.len() != offset + 4 || root != states.len() + 1 {
            return Err("Malformed MRCC excitation-energy table.".into());
        }
        let au = number(c[offset])
            .filter(|v| *v > 0.)
            .ok_or("Invalid MRCC excitation energy")?;
        let ev = au * HARTREE_TO_EV;
        states.push(ExcitedState {
            root,
            energy_au: au,
            energy_ev: ev,
            energy_cm1: ev * 8065.54393734921,
            wavelength_nm: EV_NM / ev,
            s2: None,
            multiplicity: Some(1),
            oscillator_strength: None,
            d2_au2: None,
            transition_dipole_au: None,
            excitations: vec![],
        });
    }
    if states.is_empty() {
        return Err("MRCC excitation table is empty.".into());
    }
    if let Some(begin) = lines
        .iter()
        .enumerate()
        .skip(begin + 1)
        .filter(|(_, l)| l.contains("dipole moment in length gauge [au]"))
        .map(|(index, _)| index)
        .last()
    {
        let mut count = 0;
        for l in lines.iter().skip(begin + 1) {
            let c: Vec<_> = l.split_whitespace().collect();
            let Some(root) = c.first().and_then(|s| s.parse::<usize>().ok()) else {
                if count > 0 {
                    break;
                } else {
                    continue;
                }
            };
            if c.len() < 6 || root == 0 || root > states.len() {
                return Err("Malformed MRCC oscillator-strength table.".into());
            }
            // ADC adds an irrep column before x,y,z; take the last five columns.
            let k = c.len() - 5;
            let d = [
                number(c[k]).ok_or("Invalid MRCC dipole")?,
                number(c[k + 1]).ok_or("Invalid MRCC dipole")?,
                number(c[k + 2]).ok_or("Invalid MRCC dipole")?,
            ];
            let f = number(c[k + 4])
                .filter(|v| *v >= 0.)
                .ok_or("Invalid MRCC oscillator strength")?;
            let s = &mut states[root - 1];
            s.oscillator_strength = Some(f);
            s.transition_dipole_au = Some(d);
            s.d2_au2 = number(c[k + 3]);
            count += 1;
        }
        if count != states.len() {
            return Err("Incomplete MRCC oscillator-strength table.".into());
        }
    }
    Ok(TddftResult {
        source: source.to_owned(),
        program: "MRCC".into(),
        method: lines[begin]
            .trim()
            .trim_end_matches(':')
            .replace("Final results of the ", "")
            .replace(" calculations for the excited states", ""),
        unrestricted: false,
        ground_state_s2: None,
        states,
        notes: vec![],
    })
}
