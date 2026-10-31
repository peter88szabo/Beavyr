//! Read both Behemoth's compact semiempirical and atom-labelled AO gradient tables.

pub fn parse_gradient(stdout: &str, natoms: usize) -> Result<(f64, Vec<f64>), String> {
    if stdout.contains("SCF did not converge") {
        return Err("Behemoth SCF did not converge.".into());
    }
    let lines: Vec<_> = stdout.lines().collect();
    let heading = lines
        .iter()
        .rposition(|line| line.to_ascii_lowercase().contains("nuclear gradient"))
        .ok_or("Behemoth printed no nuclear gradient")?;
    let mut energy = None;
    let mut gradient = Vec::with_capacity(3 * natoms);
    for line in &lines[heading + 1..] {
        let trimmed = line.trim();
        if let Some(value) = trimmed.strip_prefix("Energy:") {
            energy = first_number(value);
            continue;
        }
        let fields: Vec<_> = trimmed.split_whitespace().collect();
        if matches!(fields.len(), 4 | 5) {
            if let Ok(index) = fields[0].parse::<usize>() {
                let start = fields.len() - 3;
                if let Ok(values) = fields[start..]
                    .iter()
                    .map(|s| s.parse::<f64>())
                    .collect::<Result<Vec<_>, _>>()
                {
                    if index != gradient.len() / 3 + 1 {
                        return Err("Behemoth gradient atom order is invalid.".into());
                    }
                    gradient.extend(values);
                    if gradient.len() == 3 * natoms {
                        break;
                    }
                    continue;
                }
            }
        }
        if !gradient.is_empty() {
            break;
        }
    }
    if energy.is_none() {
        let before = &lines[..heading];
        let title = lines[heading].to_ascii_lowercase();
        // Correlation summaries precede a generic "Nuclear gradient" heading.
        // Never use the reference SCF energy in place of a missing correlated total.
        let mp2_section = before
            .iter()
            .rev()
            .find(|l| l.contains("analytical gradient"));
        let label = if title == "nuclear gradient"
            && mp2_section.is_some_and(|l| l.contains("RI-MP2"))
        {
            "RHF plus RI-MP2 energy"
        } else if title == "nuclear gradient" && mp2_section.is_some_and(|l| l.contains("MP2")) {
            "MP2 total energy"
        } else if title.contains("rks") || title.contains("uks") || title.contains("corrected") {
            "Total energy (all terms):"
        } else if title.contains("rhf") || title.contains("uhf") || title.contains("rohf") {
            "SCF total energy:"
        } else {
            ""
        };
        if !label.is_empty() {
            energy = before.iter().rev().find_map(|line| {
                line.trim()
                    .strip_prefix(label)
                    .and_then(|s| first_number(s.trim_start_matches(':')))
            });
        }
    }
    let energy = energy.ok_or("Behemoth printed no gradient energy")?;
    if gradient.len() != 3 * natoms || natoms == 0 {
        return Err(format!(
            "Behemoth printed {} gradient components for {natoms} atoms",
            gradient.len()
        ));
    }
    if !energy.is_finite() || gradient.iter().any(|g| !g.is_finite()) {
        return Err("Behemoth returned non-finite energy or gradient.".into());
    }
    Ok((energy, gradient))
}
fn first_number(text: &str) -> Option<f64> {
    text.split_whitespace().next()?.parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn scan_behemoth_reads_ao_and_correlated_totals_without_reference_fallback() {
        for (heading, summary) in [
            (
                "Analytical RHF nuclear gradient",
                "SCF total energy: -75.1 Hartree",
            ),
            (
                "Analytical RKS (PBE) nuclear gradient",
                "Total energy (all terms): -75.1 Hartree",
            ),
            (
                "Nuclear gradient",
                "Closed-shell all-electron MP2 analytical gradient\nMP2 total energy -75.1 Hartree",
            ),
            (
                "Nuclear gradient",
                "Frozen-core RI-MP2 analytical gradient\nRHF plus RI-MP2 energy -75.1 Hartree",
            ),
        ] {
            let log=format!("SCF total energy: -74.0 Hartree\n{summary}\n{heading}\nUnits: Hartree/bohr\nAtom Element X Y Z\n1 O 0.1 0.2 0.3\n2 H -0.1 -0.2 -0.3\n");
            assert_eq!(
                parse_gradient(&log, 2).unwrap(),
                (-75.1, vec![0.1, 0.2, 0.3, -0.1, -0.2, -0.3])
            );
            assert!(parse_gradient(&log.replace("-75.1", "NaN"), 2).is_err());
            if heading == "Nuclear gradient" {
                assert!(parse_gradient(
                    &log.lines()
                        .filter(|l| !l.contains("-75.1"))
                        .collect::<Vec<_>>()
                        .join("\n"),
                    2
                )
                .is_err());
            }
        }
    }
}
