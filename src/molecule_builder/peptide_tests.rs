use super::*;
use std::collections::BTreeMap;

fn recipe(sequence: &str) -> Recipe {
    Recipe {
        residues: parse_sequence(sequence).unwrap(),
        ..Default::default()
    }
}

fn counts(mol: &Molecule) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for a in &mol.atoms {
        *counts.entry(a.clone()).or_default() += 1;
    }
    counts
}

fn close(a: f32, b: f32) {
    assert!((a - b).abs() < 0.002, "{a} != {b}");
}
fn torsion_close(a: f32, b: f32) {
    close((a - b + 180.0).rem_euclid(360.0) - 180.0, 0.0);
}

#[test]
fn sequence_parsing_is_explicit_and_rejects_invalid_or_excessive_input() {
    assert_eq!(recipe("ala").sequence(), "ALA");
    assert_eq!(recipe("Ala-Gly-Sec-Pyl").sequence(), "AGUO");
    assert_eq!(recipe(">peptide\nAG\nSK").sequence(), "AGSK");
    assert!(parse_sequence("ABZ").is_err());
    assert!(parse_sequence("Ala-Missing").is_err());
    assert!(parse_sequence(">empty").is_err());
    assert!(parse_sequence(&"A".repeat(MAX_RESIDUES + 1)).is_err());
}

#[test]
fn condensation_removes_water_and_preserves_side_chain_charge_in_all_residues() {
    for acid in AminoAcid::ALL {
        for template in acid.templates() {
            let mut r = recipe("GGG");
            r.termini = template.termini;
            r.residues[1] = Residue {
                amino: acid,
                side_chain: template.side_chain,
            };
            let built = build(&r, None).unwrap();
            let mut expected = BTreeMap::<String, usize>::new();
            for residue in &r.residues {
                let t = residue
                    .amino
                    .template(r.termini, residue.side_chain)
                    .unwrap();
                for a in t.atoms {
                    *expected.entry(a.element.into()).or_default() += 1;
                }
            }
            *expected.get_mut("H").unwrap() -= 4;
            *expected.get_mut("O").unwrap() -= 2;
            assert_eq!(
                counts(&built.molecule),
                expected,
                "{} {:?}",
                acid.name(),
                template.side_chain
            );
            assert_eq!(built.charge, Some(template.charge as i32));
            assert!(built.molecule.pos.iter().all(|p| p.is_finite()));
            let adj = adjacency(&built.molecule);
            for (i, a) in built.molecule.atoms.iter().enumerate() {
                if a == "H" {
                    assert_eq!(adj[i].len(), 1);
                }
            }
        }
    }
}

#[test]
fn backbone_angles_trans_peptide_planes_and_proline_ring_are_preserved() {
    for shape in [Shape::Extended, Shape::Alpha, Shape::Custom] {
        let mut r = recipe("AGPITG");
        r.shape = shape;
        r.phi = -80.0;
        r.psi = 100.0;
        let parts = assemble(&r).unwrap();
        for pair in parts.windows(2) {
            let (a, b) = (&pair[0], &pair[1]);
            let (an, aca, ac, ao) = (a.pos[a.n()], a.pos[a.ca()], a.pos[a.c()], a.pos[a.o()]);
            let (bn, bca, bc) = (b.pos[b.n()], b.pos[b.ca()], b.pos[b.c()]);
            close(ac.distance(bn), PEPTIDE_BOND);
            torsion_close(torsion(an, aca, ac, bn), r.angles().1);
            torsion_close(torsion(aca, ac, bn, bca), 180.0);
            if b.template.amino_acid != AminoAcid::Proline {
                torsion_close(torsion(ac, bn, bca, bc), r.angles().0);
                let h = b
                    .neighbors(b.n())
                    .into_iter()
                    .find(|&i| b.template.atoms[i].element == "H" && !b.removed.contains(&i))
                    .unwrap();
                close(
                    (b.pos[h] - bn).dot((ac - bn).cross(bca - bn).normalize()),
                    0.0,
                );
            } else {
                let cd = b.pos[b.template.atom_index("CD").unwrap()];
                close((cd - bn).dot((ac - bn).cross(bca - bn).normalize()), 0.0);
            }
            close((ao - ac).dot((aca - ac).cross(bn - ac).normalize()), 0.0);
        }
    }
}

#[test]
fn every_source_stereocenter_and_proline_ring_survives_assembly() {
    let r = Recipe {
        residues: AminoAcid::ALL.into_iter().map(Residue::new).collect(),
        ..Default::default()
    };
    for part in assemble(&r).unwrap() {
        for center in part.template.stereocenters {
            let p = center.neighbours.map(|i| part.pos[i]);
            let handedness = (p[0] - p[3]).dot((p[1] - p[3]).cross(p[2] - p[3])).signum() as i8;
            assert_eq!(
                handedness, center.reference_handedness,
                "{}",
                part.template.id
            );
        }
        if part.template.amino_acid == AminoAcid::Proline {
            for names in [
                ("N", "CA"),
                ("CA", "CB"),
                ("CB", "CG"),
                ("CG", "CD"),
                ("CD", "N"),
            ] {
                let a = part.template.atom_index(names.0).unwrap();
                let b = part.template.atom_index(names.1).unwrap();
                let source_a =
                    Vec3::from_array(part.template.atoms[a].position_angstrom.map(|x| x as f32));
                let source_b =
                    Vec3::from_array(part.template.atoms[b].position_angstrom.map(|x| x as f32));
                close(
                    part.pos[a].distance(part.pos[b]),
                    source_a.distance(source_b),
                );
            }
        }
    }
}

#[test]
fn attachment_at_either_end_removes_leaving_atoms_and_keeps_host_coordinates_exact() {
    for code in ["A", "P", "G"] {
        for termini in [TerminalState::Neutral, TerminalState::Zwitterionic] {
            let mut r = recipe(code);
            r.termini = termini;
            let mut host = build(&r, None).unwrap().molecule;
            host.recompute_bonds(1.2, 3.0);
            let template = r.residues[0]
                .amino
                .template(termini, SideChainState::Neutral)
                .unwrap();
            for end in [End::C, End::N] {
                let index = if end == End::C {
                    template.backbone.carbonyl_carbon
                } else {
                    template.backbone.nitrogen
                };
                let anchor = attachment(&host, index, end).unwrap();
                for sequence in ["SG", "PG", "GP"] {
                    let added = build(&recipe(sequence), Some((&host, &anchor))).unwrap();
                    for (i, mapped) in added.host_mapping.iter().enumerate() {
                        if let Some(j) = mapped {
                            assert_eq!(host.pos[i], added.molecule.pos[*j]);
                            assert_eq!(host.atoms[i], added.molecule.atoms[*j]);
                        } else {
                            assert!(anchor.leaving.contains(&i));
                        }
                    }
                    let [a, b, _] = added.picks;
                    close(
                        added.molecule.pos[a].distance(added.molecule.pos[b]),
                        PEPTIDE_BOND,
                    );
                    assert!(!added.moving.contains(&a));
                    assert!(added.moving.contains(&b));
                    assert_eq!(added.charge, Some(if termini == TerminalState::Neutral { 0 } else if end == End::C { 1 } else { -1 }));
                    let expected = build(&recipe(&format!("{code}{sequence}")), None).unwrap();
                    // Free termini in mixed neutral/zwitterionic extension can
                    // change hydrogen count, but never the heavy composition.
                    for element in ["C", "N", "O", "S"] {
                        assert_eq!(
                            counts(&added.molecule).get(element),
                            counts(&expected.molecule).get(element)
                        );
                    }
                    // An existing peptide nitrogen cannot be acylated again.
                    let linked_n = if end == End::C { b } else { a };
                    assert!(attachment(&added.molecule, linked_n, End::N).is_err());
                }
            }
        }
    }
}

#[test]
fn invalid_anchors_fail_and_new_peptide_never_replaces_host() {
    let host = build(&recipe("AS"), None).unwrap().molecule;
    assert!(attachment(&host, host.atoms.len(), End::C).is_err());
    assert!(attachment(&host, 0, End::C).is_err());
    let combined = beside(build(&recipe("GK"), None).unwrap(), &host);
    assert_eq!(&combined.molecule.atoms[..host.atoms.len()], &host.atoms);
    assert_eq!(&combined.molecule.pos[..host.pos.len()], &host.pos);
    assert_eq!(combined.moving[0], host.atoms.len());
    assert!(combined
        .host_mapping
        .iter()
        .enumerate()
        .all(|(i, j)| *j == Some(i)));
}
