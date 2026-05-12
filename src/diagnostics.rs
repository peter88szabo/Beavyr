use bevy::prelude::*;

use crate::molecule::{covalent_radius_angstrom, Molecule};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiagnosticSeverity {
    Warning,
    Info,
}

#[derive(Clone, Debug)]
pub struct MoleculeDiagnostic {
    pub severity: DiagnosticSeverity,
    pub atom: Option<usize>,
    pub message: String,
}

#[derive(Clone, Debug, Default)]
pub struct DiagnosticsReport {
    pub items: Vec<MoleculeDiagnostic>,
}

impl DiagnosticsReport {
    pub fn warning_count(&self) -> usize {
        self.items
            .iter()
            .filter(|d| d.severity == DiagnosticSeverity::Warning)
            .count()
    }

    pub fn info_count(&self) -> usize {
        self.items
            .iter()
            .filter(|d| d.severity == DiagnosticSeverity::Info)
            .count()
    }
}

pub fn analyze_molecule(mol: &Molecule) -> DiagnosticsReport {
    let n = mol.atoms.len().min(mol.pos.len());
    let mut report = DiagnosticsReport::default();
    if n == 0 {
        return report;
    }

    let mut adjacency = vec![Vec::new(); n];
    for &(i, j, _) in &mol.bonds {
        if i < n && j < n {
            adjacency[i].push(j);
            adjacency[j].push(i);
        }
    }

    for i in 0..n {
        let sym = mol.atoms[i].as_str();
        let degree = adjacency[i].len();
        let h_neighbors = adjacency[i]
            .iter()
            .filter(|&&j| mol.atoms.get(j).map(|s| s.as_str()) == Some("H"))
            .count();

        if degree == 0 {
            report.items.push(MoleculeDiagnostic {
                severity: DiagnosticSeverity::Warning,
                atom: Some(i),
                message: format!("{}{} is isolated", sym, i + 1),
            });
        }

        if let Some(max_valence) = max_common_valence(sym) {
            if degree > max_valence {
                report.items.push(MoleculeDiagnostic {
                    severity: DiagnosticSeverity::Warning,
                    atom: Some(i),
                    message: format!(
                        "{}{} has valence {} > {}",
                        sym,
                        i + 1,
                        degree,
                        max_valence
                    ),
                });
            }
        }

        if let Some(neutral_valence) = neutral_single_bond_valence(sym) {
            if degree > neutral_valence {
                let hint = match sym {
                    "N" if degree == 4 => "possible formal positive nitrogen",
                    "O" if degree >= 3 => "possible oxonium or over-coordinated oxygen",
                    "S" if degree > 2 => "hypervalent sulfur candidate",
                    "P" if degree > 3 => "hypervalent phosphorus candidate",
                    _ => "formal-charge or bond-order check recommended",
                };
                report.items.push(MoleculeDiagnostic {
                    severity: DiagnosticSeverity::Info,
                    atom: Some(i),
                    message: format!("{}{}: {}", sym, i + 1, hint),
                });
            }
        }

        if let Some(target) = typical_saturated_valence(sym) {
            if should_flag_missing_hydrogen(sym, degree, h_neighbors, target) {
                report.items.push(MoleculeDiagnostic {
                    severity: DiagnosticSeverity::Info,
                    atom: Some(i),
                    message: format!(
                        "{}{} may be missing H or an explicit multiple bond (valence {}/{})",
                        sym,
                        i + 1,
                        degree,
                        target
                    ),
                });
            }
        }
    }

    add_close_contact_diagnostics(mol, n, &adjacency, &mut report);

    report
}

fn max_common_valence(sym: &str) -> Option<usize> {
    match sym {
        "H" => Some(1),
        "C" => Some(4),
        "N" => Some(4),
        "O" => Some(2),
        "F" | "Cl" | "Br" | "I" => Some(1),
        _ => None,
    }
}

fn neutral_single_bond_valence(sym: &str) -> Option<usize> {
    match sym {
        "C" => Some(4),
        "N" => Some(3),
        "O" => Some(2),
        "F" | "Cl" | "Br" | "I" => Some(1),
        "P" => Some(3),
        "S" => Some(2),
        _ => None,
    }
}

fn typical_saturated_valence(sym: &str) -> Option<usize> {
    match sym {
        "C" => Some(4),
        "N" => Some(3),
        "O" => Some(2),
        "S" => Some(2),
        "P" => Some(3),
        _ => None,
    }
}

fn should_flag_missing_hydrogen(sym: &str, degree: usize, h_neighbors: usize, target: usize) -> bool {
    if degree == 0 || degree >= target {
        return false;
    }

    match sym {
        "C" => degree <= 3,
        "N" => degree <= 2,
        "O" => degree == 1 && h_neighbors == 0,
        "S" => degree == 1 && h_neighbors == 0,
        "P" => degree <= 2,
        _ => false,
    }
}

fn add_close_contact_diagnostics(
    mol: &Molecule,
    n: usize,
    adjacency: &[Vec<usize>],
    report: &mut DiagnosticsReport,
) {
    const MAX_CLOSE_CONTACTS: usize = 20;
    let mut contacts = Vec::new();

    for i in 0..n {
        for j in (i + 1)..n {
            if adjacency[i].contains(&j) {
                continue;
            }
            let d = mol.pos[i].distance(mol.pos[j]);
            let limit = close_contact_limit(&mol.atoms[i], &mol.atoms[j]);
            if d < limit {
                contacts.push((i, j, d, limit));
            }
        }
    }

    contacts.sort_by(|a, b| {
        let severity_a = a.2 / a.3;
        let severity_b = b.2 / b.3;
        severity_a
            .partial_cmp(&severity_b)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    for (i, j, d, limit) in contacts.into_iter().take(MAX_CLOSE_CONTACTS) {
        report.items.push(MoleculeDiagnostic {
            severity: DiagnosticSeverity::Warning,
            atom: Some(i),
            message: format!(
                "Close contact: {}{}-{}{} = {:.2} Å (limit {:.2} Å)",
                mol.atoms[i],
                i + 1,
                mol.atoms[j],
                j + 1,
                d,
                limit
            ),
        });
    }
}

fn close_contact_limit(a: &str, b: &str) -> f32 {
    let ri = vdw_radius_angstrom(a).unwrap_or_else(|| covalent_radius_angstrom(a) + 0.8);
    let rj = vdw_radius_angstrom(b).unwrap_or_else(|| covalent_radius_angstrom(b) + 0.8);
    (ri + rj) * 0.55
}

fn vdw_radius_angstrom(sym: &str) -> Option<f32> {
    match sym {
        "H" => Some(1.20),
        "C" => Some(1.70),
        "N" => Some(1.55),
        "O" => Some(1.52),
        "F" => Some(1.47),
        "P" => Some(1.80),
        "S" => Some(1.80),
        "Cl" => Some(1.75),
        "Br" => Some(1.85),
        "I" => Some(1.98),
        _ => None,
    }
}
