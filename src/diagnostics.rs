use bevy::prelude::*;
use std::collections::HashSet;

use crate::molecule::{covalent_radius_angstrom, Molecule, SpatialGrid};

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

/// Cached diagnostics, refreshed only after the molecule changes.
#[derive(Resource, Default)]
pub struct DiagnosticsCache {
    pub report: DiagnosticsReport,
    /// Set by the UI while the Diagnostics panel is expanded.  The analysis is
    /// far too expensive to run per frame during trajectory playback for a
    /// panel nobody is looking at.
    pub wanted: bool,
}

/// Keep the expensive diagnostics pass out of the egui frame loop.
///
/// Runs only while the Diagnostics panel is open, and then only when the
/// molecule actually changed (plus once on the frame the panel is opened, so
/// it never shows a stale or empty report).
pub fn refresh_diagnostics_cache(
    mut changes: MessageReader<crate::events::MoleculeChanged>,
    mol: Res<Molecule>,
    mut cache: ResMut<DiagnosticsCache>,
    mut was_wanted: Local<bool>,
) {
    let molecule_changed = changes.read().next().is_some();
    let just_opened = cache.wanted && !*was_wanted;
    *was_wanted = cache.wanted;

    if !cache.wanted {
        return;
    }
    if molecule_changed || just_opened {
        cache.report = analyze_molecule(&mol);
    }
}

pub fn analyze_molecule(mol: &Molecule) -> DiagnosticsReport {
    let n = mol.atoms.len().min(mol.pos.len());
    let mut report = DiagnosticsReport::default();
    if n == 0 {
        return report;
    }

    let mut adjacency = vec![Vec::new(); n];
    // Membership set for the close-contact pass, which would otherwise do a
    // linear scan of `adjacency[i]` inside an already quadratic loop.
    let mut bonded_pairs: HashSet<(usize, usize)> = HashSet::with_capacity(mol.bonds.len() * 2);
    for &(i, j, _) in &mol.bonds {
        if i < n && j < n {
            adjacency[i].push(j);
            adjacency[j].push(i);
            bonded_pairs.insert((i.min(j), i.max(j)));
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
                    message: format!("{}{} has valence {} > {}", sym, i + 1, degree, max_valence),
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

    add_close_contact_diagnostics(mol, n, &bonded_pairs, &mut report);

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

fn should_flag_missing_hydrogen(
    sym: &str,
    degree: usize,
    h_neighbors: usize,
    target: usize,
) -> bool {
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
    bonded_pairs: &HashSet<(usize, usize)>,
    report: &mut DiagnosticsReport,
) {
    const MAX_CLOSE_CONTACTS: usize = 20;
    let mut contacts = Vec::new();

    // A contact can only be flagged when d < (r_i + r_j) * 0.55, so no pair further
    // apart than `max_limit` matters.  Sizing the grid to that bound turns the
    // all-pairs sweep into a neighbour query.
    let max_effective_radius = mol.atoms[..n]
        .iter()
        .map(|sym| effective_contact_radius(sym))
        .fold(0.0_f32, f32::max);
    let max_limit = (max_effective_radius * 2.0 * 0.55).max(0.5);

    let grid = SpatialGrid::build(0..n, &mol.pos, max_limit);
    let mut candidates: Vec<usize> = Vec::new();

    for i in 0..n {
        candidates.clear();
        grid.collect_neighbors(mol.pos[i], &mut candidates);
        for &j in &candidates {
            // Each unordered pair once.
            if j <= i {
                continue;
            }
            if bonded_pairs.contains(&(i, j)) {
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

fn effective_contact_radius(sym: &str) -> f32 {
    vdw_radius_angstrom(sym).unwrap_or_else(|| covalent_radius_angstrom(sym) + 0.8)
}

fn close_contact_limit(a: &str, b: &str) -> f32 {
    (effective_contact_radius(a) + effective_contact_radius(b)) * 0.55
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

