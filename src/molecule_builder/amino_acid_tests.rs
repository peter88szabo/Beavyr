use super::*;
use std::collections::HashSet;

fn neighbours(t: &AminoAcidTemplate, i: usize) -> Vec<usize> {
    t.bonds
        .iter()
        .filter_map(|b| {
            if b.atoms[0] == i {
                Some(b.atoms[1])
            } else if b.atoms[1] == i {
                Some(b.atoms[0])
            } else {
                None
            }
        })
        .collect()
}

fn hydrogen_count(t: &AminoAcidTemplate, name: &str) -> usize {
    neighbours(t, t.atom_index(name).unwrap())
        .iter()
        .filter(|&&i| t.atoms[i].element == "H")
        .count()
}

fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[test]
fn all_residues_and_supported_states_are_unambiguous() {
    assert_eq!(AminoAcid::ALL.len(), 22);
    assert_eq!(TEMPLATES.len(), 62);
    let mut ids = HashSet::new();
    let mut letters = HashSet::new();
    for residue in AminoAcid::ALL {
        assert!(letters.insert(residue.one_letter()));
        assert_eq!(
            AminoAcid::from_one_letter(residue.one_letter().to_ascii_lowercase()),
            Some(residue)
        );
        assert_eq!(
            AminoAcid::from_three_letter(residue.three_letter()),
            Some(residue)
        );
        let expected = match residue {
            AminoAcid::Histidine => 6,
            AminoAcid::AsparticAcid
            | AminoAcid::GlutamicAcid
            | AminoAcid::Cysteine
            | AminoAcid::Tyrosine
            | AminoAcid::Selenocysteine
            | AminoAcid::Arginine
            | AminoAcid::Lysine => 4,
            _ => 2,
        };
        assert_eq!(residue.templates().len(), expected);
        let mut states = HashSet::new();
        for t in residue.templates() {
            assert_eq!(t.amino_acid, residue);
            assert!(ids.insert(t.id), "duplicate {}", t.id);
            assert!(states.insert((t.termini, t.side_chain)));
            assert!(std::ptr::eq(
                residue.template(t.termini, t.side_chain).unwrap(),
                t
            ));
        }
    }
    assert!(AminoAcid::from_one_letter('X').is_none());
    assert!(AminoAcid::from_three_letter("unknown").is_none());
    assert!(AminoAcid::Alanine
        .template(TerminalState::Neutral, SideChainState::Protonated)
        .is_none());
    assert!(AminoAcid::Histidine
        .template(TerminalState::Neutral, SideChainState::Neutral)
        .is_none());
}

#[test]
fn elemental_compositions_and_electron_counts_match_the_named_forms() {
    // Independent neutral free-monomer formulas: C,H,N,O,S,Se.
    let formulas = [
        [2, 5, 1, 2, 0, 0],
        [3, 7, 1, 2, 0, 0],
        [5, 11, 1, 2, 0, 0],
        [6, 13, 1, 2, 0, 0],
        [6, 13, 1, 2, 0, 0],
        [3, 7, 1, 3, 0, 0],
        [4, 9, 1, 3, 0, 0],
        [3, 7, 1, 2, 1, 0],
        [5, 11, 1, 2, 1, 0],
        [5, 9, 1, 2, 0, 0],
        [4, 7, 1, 4, 0, 0],
        [5, 9, 1, 4, 0, 0],
        [4, 8, 2, 3, 0, 0],
        [5, 10, 2, 3, 0, 0],
        [6, 9, 3, 2, 0, 0],
        [6, 14, 2, 2, 0, 0],
        [6, 14, 4, 2, 0, 0],
        [9, 11, 1, 2, 0, 0],
        [9, 11, 1, 3, 0, 0],
        [11, 12, 2, 2, 0, 0],
        [3, 7, 1, 2, 0, 1],
        [12, 21, 3, 3, 0, 0],
    ];
    for (residue, base) in AminoAcid::ALL.into_iter().zip(formulas) {
        for t in residue.templates() {
            let charge = match t.side_chain {
                SideChainState::Protonated => 1,
                SideChainState::Deprotonated => -1,
                _ => 0,
            };
            let mut expected = base;
            expected[1] += charge;
            let mut actual = [0; 6];
            for atom in t.atoms {
                let i = ["C", "H", "N", "O", "S", "Se"]
                    .iter()
                    .position(|&s| s == atom.element)
                    .unwrap();
                actual[i] += 1;
            }
            assert_eq!(actual, expected, "{}", t.id);
            assert_eq!(t.charge as i32, charge, "{}", t.id);
            assert_eq!(
                t.atoms.iter().map(|a| a.formal_charge as i32).sum::<i32>(),
                charge
            );
            let electrons = actual
                .iter()
                .zip([6, 1, 7, 8, 16, 34])
                .map(|(n, z)| n * z)
                .sum::<i32>()
                - charge;
            assert_eq!(electrons % 2, 0);
            assert_eq!(t.multiplicity, 1);
        }
    }
}

#[test]
fn termini_and_side_chain_hydrogens_match_the_explicit_protonation_labels() {
    for t in &TEMPLATES {
        let backbone = t.backbone;
        for (index, name) in [
            (backbone.nitrogen, "N"),
            (backbone.alpha_carbon, "CA"),
            (backbone.carbonyl_carbon, "C"),
            (backbone.carbonyl_oxygen, "O"),
            (backbone.terminal_oxygen, "OXT"),
        ] {
            assert_eq!(t.atoms[index].name, name);
        }
        let zwitterionic = t.termini == TerminalState::Zwitterionic;
        assert_eq!(
            t.atoms[backbone.nitrogen].formal_charge,
            i8::from(zwitterionic)
        );
        assert_eq!(
            t.atoms[backbone.terminal_oxygen].formal_charge,
            -i8::from(zwitterionic)
        );
        assert_eq!(
            hydrogen_count(t, "N"),
            if t.amino_acid == AminoAcid::Proline {
                1
            } else {
                2
            } + usize::from(zwitterionic)
        );
        assert_eq!(hydrogen_count(t, "OXT"), usize::from(!zwitterionic));
        match t.amino_acid {
            AminoAcid::Histidine => {
                let expected = match t.side_chain {
                    SideChainState::HistidineDelta => (1, 0),
                    SideChainState::HistidineEpsilon => (0, 1),
                    SideChainState::Protonated => (1, 1),
                    _ => panic!("ambiguous histidine state"),
                };
                assert_eq!(
                    (hydrogen_count(t, "ND1"), hydrogen_count(t, "NE2")),
                    expected
                );
            }
            AminoAcid::Lysine => assert_eq!(
                hydrogen_count(t, "NZ"),
                2 + usize::from(t.side_chain == SideChainState::Protonated)
            ),
            AminoAcid::Arginine => assert_eq!(
                hydrogen_count(t, "NH2"),
                1 + usize::from(t.side_chain == SideChainState::Protonated)
            ),
            _ => {}
        }
        let acidic = match t.amino_acid {
            AminoAcid::AsparticAcid => Some("OD2"),
            AminoAcid::GlutamicAcid => Some("OE2"),
            AminoAcid::Cysteine => Some("SG"),
            AminoAcid::Tyrosine => Some("OH"),
            AminoAcid::Selenocysteine => Some("SE"),
            _ => None,
        };
        if let Some(name) = acidic {
            assert_eq!(
                hydrogen_count(t, name),
                usize::from(t.side_chain == SideChainState::Neutral)
            );
        }
    }
}

#[test]
fn all_graphs_are_connected_with_unique_names_bonds_and_attached_hydrogens() {
    for t in &TEMPLATES {
        let names: HashSet<_> = t.atoms.iter().map(|a| a.name).collect();
        assert_eq!(names.len(), t.atoms.len(), "{}", t.id);
        let mut bonds = HashSet::new();
        for b in t.bonds {
            let [a, b] = b.atoms;
            assert!(a != b && a < t.atoms.len() && b < t.atoms.len());
            assert!(bonds.insert((a.min(b), a.max(b))));
            let d = difference(t.atoms[a].position_angstrom, t.atoms[b].position_angstrom);
            let distance = d.iter().map(|v| v * v).sum::<f64>().sqrt();
            assert!(
                (0.65..2.4).contains(&distance),
                "{} bond {a}-{b}: {distance}",
                t.id
            );
        }
        let mut seen = HashSet::from([0]);
        let mut stack = vec![0];
        while let Some(i) = stack.pop() {
            for j in neighbours(t, i) {
                if seen.insert(j) {
                    stack.push(j);
                }
            }
        }
        assert_eq!(seen.len(), t.atoms.len());
        for (i, a) in t.atoms.iter().enumerate() {
            assert!(a.position_angstrom.iter().all(|v| v.is_finite()));
            if a.element == "H" {
                let adjacent = neighbours(t, i);
                assert_eq!(adjacent.len(), 1);
                assert_ne!(t.atoms[adjacent[0]].element, "H");
            }
        }
    }
}

#[test]
fn optimized_stereocenters_keep_the_curated_handedness_including_threonine_and_isoleucine() {
    for t in &TEMPLATES {
        let count = match t.amino_acid {
            AminoAcid::Glycine => 0,
            AminoAcid::Isoleucine | AminoAcid::Threonine => 2,
            AminoAcid::Pyrrolysine => 3,
            _ => 1,
        };
        assert_eq!(t.stereocenters.len(), count);
        for center in t.stereocenters {
            let p = center.neighbours.map(|i| t.atoms[i].position_angstrom);
            let a = difference(p[0], p[3]);
            let b = difference(p[1], p[3]);
            let c = difference(p[2], p[3]);
            let volume = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]);
            assert!(volume.abs() > 0.1, "{} flat stereocenter", t.id);
            assert_eq!(
                volume.signum() as i8,
                center.reference_handedness,
                "{}",
                t.id
            );
            let expected = match (t.amino_acid, t.atoms[center.atom].name) {
                (AminoAcid::Cysteine | AminoAcid::Selenocysteine, "CA") => Configuration::R,
                (AminoAcid::Threonine, "CB") => Configuration::R,
                (AminoAcid::Pyrrolysine, "CA2" | "CG2") => Configuration::R,
                _ => Configuration::S,
            };
            assert_eq!(center.configuration, expected);
        }
    }
}

#[test]
fn compiled_coordinates_and_energies_match_the_actual_xtb_output_files() {
    let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("data/amino_acids/optimized");
    for t in &TEMPLATES {
        // Test-only file access. Production templates are entirely static data.
        let text = std::fs::read_to_string(base.join(format!("{}.xyz", t.id))).unwrap();
        let mut lines = text.lines();
        assert_eq!(
            lines.next().unwrap().parse::<usize>().unwrap(),
            t.atoms.len()
        );
        let header: Vec<_> = lines.next().unwrap().split_whitespace().collect();
        assert!((header[1].parse::<f64>().unwrap() - t.optimization.energy_hartree).abs() < 1e-11);
        let coords: Vec<_> = lines.collect();
        assert_eq!(coords.len(), t.atoms.len());
        for (line, atom) in coords.into_iter().zip(t.atoms) {
            let fields: Vec<_> = line.split_whitespace().collect();
            assert_eq!(fields[0], atom.element);
            for axis in 0..3 {
                assert!(
                    (fields[axis + 1].parse::<f64>().unwrap() - atom.position_angstrom[axis]).abs()
                        < 1e-11
                );
            }
        }
        let o = t.optimization;
        assert!(o.energy_hartree.is_finite());
        assert!((0.0..=1e-4).contains(&o.gradient_norm_hartree_per_bohr));
        assert!(o.lowest_vibrational_wavenumber_cm1 > 0.0);
        assert_eq!(o.vibrational_mode_count, 3 * t.atoms.len() - 6);
        for hash in [o.ccd_sha256, o.xtb_log_sha256] {
            assert_eq!(hash.len(), 64);
            assert!(hash.bytes().all(|b| b.is_ascii_hexdigit()));
        }
    }
}
