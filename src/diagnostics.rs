use bevy::prelude::*;
use std::collections::HashSet;

use crate::bond_order::{classify_valence, total_bond_orders, ValenceVerdict};

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

    // Bond *orders*, not neighbour counts. A neighbour count calls every sp2
    // carbon under-coordinated (three neighbours, but one of them a double
    // bond) and every transferring hydrogen over-coordinated (two neighbours,
    // but each only half a bond). Summing estimated orders gets both right,
    // which is what makes a missing hydrogen worth reporting at all.
    let totals = total_bond_orders(mol);

    for i in 0..n {
        let sym = mol.atoms[i].as_str();
        let degree = adjacency[i].len();
        let total = totals.get(i).copied().unwrap_or(0.0);

        if degree == 0 {
            report.items.push(MoleculeDiagnostic {
                severity: DiagnosticSeverity::Warning,
                atom: Some(i),
                message: format!("{}{} is isolated", sym, i + 1),
            });
            continue;
        }

        let (verdict, reference) = classify_valence(sym, total);
        let label = format!("{}{}", sym, i + 1);
        match verdict {
            ValenceVerdict::Satisfied => {}
            ValenceVerdict::Short => {
                let target = reference.unwrap_or(total);
                let missing = target - total;
                // A whole bond short of the usual valence is the case worth
                // naming: on a closed-shell structure that is a hydrogen
                // nobody added. Geometry cannot tell that from a radical, so
                // both readings are offered rather than one asserted.
                let hydrogens = (missing + 0.25).floor().max(1.0) as usize;
                report.items.push(MoleculeDiagnostic {
                    severity: DiagnosticSeverity::Info,
                    atom: Some(i),
                    message: format!(
                        "{label}: bond order {:.1} of {:.0} \u{2014} {} H may be missing, \
                         or this is a radical",
                        total, target, hydrogens
                    ),
                });
            }
            ValenceVerdict::Excess => {
                let target = reference.unwrap_or(total);
                report.items.push(MoleculeDiagnostic {
                    severity: DiagnosticSeverity::Warning,
                    atom: Some(i),
                    message: format!(
                        "{label}: bond order {:.1} exceeds {:.0}, the highest common \
                         valence for {sym}",
                        total, target
                    ),
                });
            }
            ValenceVerdict::Unusual => {
                report.items.push(MoleculeDiagnostic {
                    severity: DiagnosticSeverity::Info,
                    atom: Some(i),
                    message: format!(
                        "{label}: bond order {:.1} sits between {sym}'s common valences",
                        total
                    ),
                });
            }
        }
    }

    add_close_contact_diagnostics(mol, n, &bonded_pairs, &mut report);

    report
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


#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::Molecule;
    use bevy::prelude::Vec3;

    fn from_xyz_text(text: &str) -> Molecule {
        let mut mol = Molecule::from_xyz(text);
        mol.recompute_bonds(1.2, 2.5);
        mol
    }

    fn from_example(name: &str) -> Molecule {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join(name);
        from_xyz_text(&std::fs::read_to_string(&path).unwrap())
    }

    fn messages(mol: &Molecule) -> Vec<String> {
        analyze_molecule(mol)
            .items
            .iter()
            .map(|d| d.message.clone())
            .collect()
    }

    /// The transition state from examples/: its transferring hydrogen sits
    /// between donor and acceptor and used to be reported as divalent.
    fn transition_state() -> Molecule {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("c164-ts1-2-1001.log");
        let text = std::fs::read_to_string(&path).unwrap();
        let g = crate::qchem_interfaces::gaussian_log::parse_gaussian_log(&text).unwrap();
        let mut mol = Molecule::empty();
        mol.atoms = g.atoms.clone();
        mol.pos = g
            .coords_angstrom
            .chunks(3)
            .map(|c| Vec3::new(c[0] as f32, c[1] as f32, c[2] as f32))
            .collect();
        mol.recompute_bonds(1.2, 2.5);
        mol
    }

    /// The whole point: a transferring hydrogen has two half bonds, which is
    /// one bond, not two.
    #[test]
    fn a_transferring_hydrogen_is_not_reported_as_divalent() {
        let found = messages(&transition_state());
        assert!(
            !found.iter().any(|m| m.starts_with("H18")),
            "the transferring H is still flagged: {found:?}"
        );
    }

    /// An sp2 carbon has three neighbours but a full valence of four. A real
    /// aromatic molecule must produce nothing at all.
    #[test]
    fn an_aromatic_molecule_is_clean() {
        let found = messages(&from_example("phenyl-OCH3.xyz"));
        assert!(found.is_empty(), "a plain aromatic ether should be clean: {found:?}");
    }

    /// The feature the panel is for: a genuinely absent hydrogen.
    #[test]
    fn a_missing_hydrogen_is_reported_with_how_many() {
        // Methane with one H removed.
        let ch3 = from_xyz_text(
            "C  0.000  0.000  0.000\n\
             H  0.629  0.629  0.629\n\
             H -0.629 -0.629  0.629\n\
             H -0.629  0.629 -0.629\n",
        );
        let found = messages(&ch3);
        assert_eq!(found.len(), 1, "{found:?}");
        assert!(found[0].starts_with("C1:"), "{found:?}");
        assert!(found[0].contains("1 H may be missing"), "{found:?}");
        assert!(found[0].contains("radical"), "both readings offered: {found:?}");
    }

    #[test]
    fn two_missing_hydrogens_are_counted_as_two() {
        let ch2 = from_xyz_text(
            "C  0.000  0.000  0.000\n\
             H  0.629  0.629  0.629\n\
             H -0.629 -0.629  0.629\n",
        );
        let found = messages(&ch2);
        assert_eq!(found.len(), 1, "{found:?}");
        assert!(found[0].contains("2 H may be missing"), "{found:?}");
    }

    /// Complete small molecules must be silent.
    #[test]
    fn complete_molecules_produce_nothing() {
        let water = from_xyz_text(
            "O  0.000  0.000  0.119\nH  0.000  0.763 -0.477\nH  0.000 -0.763 -0.477\n",
        );
        assert!(messages(&water).is_empty(), "{:?}", messages(&water));

        // Ethylene: each carbon has two H and one C=C, summing to four.
        let ethylene = from_xyz_text(
            "C  0.000  0.000  0.667\n\
             C  0.000  0.000 -0.667\n\
             H  0.000  0.923  1.238\n\
             H  0.000 -0.923  1.238\n\
             H  0.000  0.923 -1.238\n\
             H  0.000 -0.923 -1.238\n",
        );
        assert!(messages(&ethylene).is_empty(), "{:?}", messages(&ethylene));
    }

    /// An isolated atom is still worth saying, but only once and not as a
    /// valence complaint.
    #[test]
    fn a_lone_atom_is_reported_as_isolated_only() {
        let lone = from_xyz_text("C  0.000  0.000  0.000\nO 10.0 10.0 10.0\n");
        let found = messages(&lone);
        assert_eq!(found.len(), 2, "{found:?}");
        assert!(found.iter().all(|m| m.contains("isolated")), "{found:?}");
    }

    #[test]
    fn an_empty_molecule_produces_nothing() {
        assert!(messages(&Molecule::empty()).is_empty());
    }
}
