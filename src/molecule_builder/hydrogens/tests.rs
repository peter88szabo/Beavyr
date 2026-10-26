use super::*;

fn molecule(atoms: &[&str], pos: &[Vec3]) -> Molecule {
    Molecule {
        atoms: atoms.iter().map(|s| s.to_string()).collect(),
        pos: pos.to_vec(),
        bonds: vec![],
        hydrogen_bonds: vec![],
        topology: None,
    }
}

fn angle(mol: &Molecule, center: usize, a: usize, b: usize) -> f32 {
    (mol.pos[a] - mol.pos[center])
        .normalize()
        .dot((mol.pos[b] - mol.pos[center]).normalize())
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees()
}

#[test]
fn bare_carbon_becomes_tetrahedral_methane_and_repeated_click_is_a_noop() {
    let origin = Vec3::new(3.0, -2.0, 1.0);
    let mut mol = molecule(&["C"], &[origin]);
    assert_eq!(complete(&mut mol).unwrap().added, 4);
    assert_eq!(mol.pos[0], origin);
    for a in 1..5 {
        assert!((mol.pos[a].distance(origin) - single_bond_length("C", "H")).abs() < 1.0e-5);
        for b in a + 1..5 {
            assert!((angle(&mol, 0, a, b) - 109.471).abs() < 0.02);
        }
    }
    let positions = mol.pos.clone();
    assert_eq!(complete(&mut mol).unwrap().added, 0);
    assert_eq!(positions, mol.pos);
}

#[test]
fn isolated_oxygen_and_existing_oh_form_bent_water() {
    for mut mol in [
        molecule(&["O"], &[Vec3::ZERO]),
        molecule(&["O", "H"], &[Vec3::ZERO, Vec3::X * 0.97]),
    ] {
        let before = mol.pos.clone();
        assert_eq!(complete(&mut mol).unwrap().added, 3 - before.len());
        assert_eq!(mol.atoms.len(), 3);
        assert_eq!(&mol.pos[..before.len()], &before);
        assert!((angle(&mol, 0, 1, 2) - 104.5).abs() < 0.02);
        for h in 1..3 {
            assert!((mol.pos[h].length() - 0.97).abs() < 1.0e-5);
        }
    }
}

#[test]
fn methanol_skeleton_fills_both_carbon_and_oxygen_without_moving_them() {
    let mut mol = molecule(&["C", "O"], &[Vec3::ZERO, Vec3::X * 1.43]);
    let before = mol.pos.clone();
    assert_eq!(complete(&mut mol).unwrap().added, 4);
    assert_eq!(&mol.pos[..2], &before);
    mol.recompute_bonds(1.2, 2.5);
    assert_eq!(mol.bonds.len(), 5);
    let oh = mol
        .bonds
        .iter()
        .find_map(|&(a, b, _)| (a == 1 && b > 1).then_some(b))
        .unwrap();
    assert!((angle(&mol, 1, 0, oh) - 104.5).abs() < 0.02);
    for h in 2..5 {
        assert!((angle(&mol, 0, 1, h) - 109.471).abs() < 0.02);
    }
    assert_eq!(complete(&mut mol).unwrap().added, 0);
}

#[test]
fn ethane_has_no_extra_cross_bonds_even_with_an_inflated_display_threshold() {
    let mut mol = molecule(&["C", "C"], &[Vec3::ZERO, Vec3::X * 1.54]);
    mol.recompute_bonds(3.0, 4.0);
    assert_eq!(complete(&mut mol).unwrap().added, 6);
    mol.recompute_bonds(1.2, 2.5);
    assert_eq!(mol.bonds.len(), 7);
    for i in 2..8 {
        assert_eq!(
            mol.bonds
                .iter()
                .filter(|&&(a, b, _)| a == i || b == i)
                .count(),
            1
        );
    }
    assert_eq!(complete(&mut mol).unwrap().added, 0);
}

#[test]
fn existing_hydrogens_and_non_target_elements_are_preserved() {
    let mut mol = molecule(
        &["Ne", "C", "H"],
        &[Vec3::new(20.0, 0.0, 0.0), Vec3::ZERO, Vec3::X * 1.07],
    );
    let before = mol.clone();
    assert_eq!(complete(&mut mol).unwrap().added, 3);
    assert_eq!(&mol.atoms[..3], &before.atoms);
    assert_eq!(&mol.pos[..3], &before.pos);
    assert!(mol.pos[3..].iter().all(|p| p.distance(mol.pos[0]) > 10.0));
}

#[test]
fn unsaturated_centres_get_only_their_open_valence() {
    // C=O becomes formaldehyde and C=C ethylene: the double bond counts
    // toward the valence, so each carbon gets two H and the oxygen none.
    for (symbols, distance, added) in [(["C", "O"], 1.22, 2), (["C", "C"], 1.34, 4)] {
        let mut mol = molecule(&symbols, &[Vec3::ZERO, Vec3::X * distance]);
        let result = complete(&mut mol).unwrap();
        assert_eq!(result.added, added, "{symbols:?} at {distance}");
        assert_eq!(result.skipped, 0);
    }
}

#[test]
fn ambiguous_lone_aromatic_and_partial_bonds_are_not_saturated() {
    for (symbols, distance) in [(["C", "C"], 1.40), (["C", "C"], 1.70)] {
        let mut mol = molecule(&symbols, &[Vec3::ZERO, Vec3::X * distance]);
        let result = complete(&mut mol).unwrap();
        assert_eq!(result.added, 0, "{symbols:?} at {distance}");
        assert_eq!(result.skipped, 2);
        assert_eq!(mol.atoms.len(), 2);
    }
}

#[test]
fn divalent_ether_oxygen_does_not_get_an_extra_hydrogen() {
    let theta = 110_f32.to_radians();
    let mut mol = molecule(
        &["O", "C", "C"],
        &[
            Vec3::ZERO,
            Vec3::X * 1.43,
            Vec3::new(theta.cos(), theta.sin(), 0.0) * 1.43,
        ],
    );
    assert_eq!(complete(&mut mol).unwrap().added, 6);
    mol.recompute_bonds(1.2, 2.5);
    assert_eq!(
        mol.bonds
            .iter()
            .filter(|&&(a, b, _)| a == 0 || b == 0)
            .count(),
        2
    );
}

#[test]
fn malformed_coordinates_and_overlapping_atoms_are_not_modified() {
    for mut mol in [
        molecule(&["C"], &[]),
        molecule(&["O"], &[Vec3::splat(f32::INFINITY)]),
    ] {
        let before = mol.clone();
        assert!(complete(&mut mol).is_err());
        assert_eq!(mol.atoms, before.atoms);
        assert_eq!(mol.pos, before.pos);
    }
    let mut mol = molecule(&["C", "C"], &[Vec3::ZERO, Vec3::ZERO]);
    assert_eq!(complete(&mut mol).unwrap().added, 0);
    assert_eq!(mol.atoms.len(), 2);
}

#[test]
fn smart_hydrogen_refuses_a_full_valence() {
    let mut mol = molecule(&["C"], &[Vec3::ZERO]);
    complete(&mut mol).unwrap();
    let before = mol.pos.clone();
    assert!(append_smart_hydrogen(&mut mol, 0)
        .unwrap_err()
        .contains("full valence"));
    assert_eq!(mol.pos, before);
}

#[test]
fn bulk_and_repeated_smart_additions_produce_identical_coordinates() {
    let initial = molecule(
        &["C", "O"],
        &[Vec3::new(-4.0, 2.0, 0.5), Vec3::new(-2.57, 2.0, 0.5)],
    );
    let mut bulk = initial.clone();
    let mut manual = initial;
    complete(&mut bulk).unwrap();
    for host in [0, 0, 0, 1] {
        append_smart_hydrogen(&mut manual, host).unwrap();
    }
    assert_eq!(manual.atoms, bulk.atoms);
    assert_eq!(manual.pos, bulk.pos);
}

#[test]
fn undo_restores_exact_original_and_survives_a_repeated_noop_click() {
    let mut mol = molecule(&["C", "O"], &[Vec3::ZERO, Vec3::X * 1.43]);
    let before = mol.clone();
    let mut state = HydrogenState::default();
    let mut settings = MolSettings::default();
    assert!(state.add(&mut mol, &mut settings));
    assert!(settings.geometry_dirty && settings.bond_topology_dirty);
    assert!(!state.add(&mut mol, &mut settings));
    assert!(state.undo(&mut mol, &mut settings));
    assert_eq!(mol.atoms, before.atoms);
    assert_eq!(mol.pos, before.pos);
    assert!(!state.undo(&mut mol, &mut settings));
}

#[test]
fn undo_cannot_overwrite_later_edits_or_another_loaded_structure() {
    for replace in [false, true] {
        let mut mol = molecule(&["C"], &[Vec3::ZERO]);
        let mut state = HydrogenState::default();
        let mut settings = MolSettings::default();
        state.add(&mut mol, &mut settings);
        if replace {
            mol = molecule(&["O"], &[Vec3::ZERO]);
        } else {
            mol.pos[1].x += 0.05;
        }
        let expected = mol.clone();
        assert!(!state.undo(&mut mol, &mut settings));
        assert_eq!(mol.atoms, expected.atoms);
        assert_eq!(mol.pos, expected.pos);
    }
}

#[test]
fn new_hydrogen_editor_rows_reference_their_hosts_and_preserve_geometry() {
    for (mut mol, initial) in [
        (molecule(&["C"], &[Vec3::ZERO]), 1),
        (molecule(&["C", "O"], &[Vec3::ZERO, Vec3::X * 1.43]), 2),
    ] {
        complete(&mut mol).unwrap();
        let zmat = builder_zmat(&mol, initial);
        for row in &zmat[initial..] {
            let host = row.bond_ref.unwrap() - 1;
            assert!(host < initial);
            assert!(
                (row.bond_len - single_bond_length(&mol.atoms[host], "H") as f64).abs() < 1.0e-5
            );
        }
        let rebuilt = zmat2xyz::zmat_to_xyz(&zmat);
        for a in 0..mol.atoms.len() {
            for b in a + 1..mol.atoms.len() {
                assert!(
                    (rebuilt[a].distance(rebuilt[b]) - mol.pos[a].distance(mol.pos[b])).abs()
                        < 1.0e-4
                );
            }
        }
    }
}

#[test]
fn nitrogen_completion_adds_ammonia_and_planar_amide_hydrogens() {
    let mut ammonia=molecule(&["N"],&[Vec3::ZERO]);assert_eq!(complete(&mut ammonia).unwrap().added,3);assert_eq!(complete(&mut ammonia).unwrap().added,0);
    // Formamide heavy atoms. The short amide C-N must not be treated as a double bond.
    let mut amide=molecule(&["C","O","N"],&[Vec3::ZERO,Vec3::X*1.23,Vec3::new(-0.67,1.16,0.0)]);
    complete(&mut amide).unwrap();amide.recompute_bonds(1.2,2.5);
    let hs:Vec<_>=amide.bonds.iter().filter_map(|&(a,b,_)| if a==2 && amide.atoms[b]=="H" {Some(b)} else if b==2 && amide.atoms[a]=="H" {Some(a)} else {None}).collect();
    assert_eq!(hs.len(),2);for h in hs {assert!(amide.pos[h].z.abs()<1e-4,"amide H must be planar");}
    assert_eq!(complete(&mut amide).unwrap().added,0);
}

#[test]
fn open_sites_on_aromatic_alkene_and_carbonyl_carbons_get_one_h() {
    // Benzene with one H removed (C6H5).
    let r = 1.39_f32;
    let rh = r + 1.09;
    let mut atoms = vec![];
    let mut pos = vec![];
    for k in 0..6 {
        let a = (k as f32) * std::f32::consts::PI / 3.0;
        atoms.push("C");
        pos.push(Vec3::new(r * a.cos(), r * a.sin(), 0.0));
    }
    for k in 1..6 {
        let a = (k as f32) * std::f32::consts::PI / 3.0;
        atoms.push("H");
        pos.push(Vec3::new(rh * a.cos(), rh * a.sin(), 0.0));
    }
    let mut phenyl = molecule(&atoms, &pos);
    // Vinyl, H2C=CH.
    let mut vinyl = molecule(
        &["C", "C", "H", "H", "H"],
        &[
            Vec3::ZERO,
            Vec3::new(1.33, 0.0, 0.0),
            Vec3::new(-0.55, 0.94, 0.0),
            Vec3::new(-0.55, -0.94, 0.0),
            Vec3::new(1.88, 0.94, 0.0),
        ],
    );
    // Formyl, H-C=O.
    let mut formyl = molecule(
        &["C", "O", "H"],
        &[Vec3::ZERO, Vec3::new(1.21, 0.0, 0.0), Vec3::new(-0.55, 0.94, 0.0)],
    );
    for (name, mol, host) in [
        ("phenyl", &mut phenyl, 0),
        ("vinyl", &mut vinyl, 1),
        ("formyl", &mut formyl, 0),
    ] {
        let before = mol.pos.clone();
        let result = complete(mol).unwrap();
        assert_eq!((result.added, result.skipped), (1, 0), "{name}");
        assert_eq!(mol.pos[..before.len()], before[..], "{name}: existing atoms moved");
        let h = *mol.pos.last().unwrap();
        let length = h.distance(mol.pos[host]);
        assert!((length - single_bond_length("C", "H")).abs() < 1.0e-4, "{name}");
        // The new H lies in the plane of the sp2 centre.
        assert!(h.z.abs() < 1.0e-3, "{name}: H out of plane");
        assert_eq!(complete(mol).unwrap().added, 0, "{name}: second click adds nothing");
    }
}

#[test]
fn conjugated_ring_carbon_with_double_bond_lengths_on_both_sides_gets_its_h() {
    // N-methyl-2-nitropyrrole as drawn in the editor, missing the H on C6.
    // C6 has C=N (1.30 Å) and C=C (1.34 Å) lengths at once; summing those
    // orders reads as a full valence, which is impossible for a ring carbon.
    let rows = [
        ("C", [0.0, 0.0, 0.0]),
        ("H", [1.08899999, 0.0, 0.0]),
        ("H", [-0.36300099, 1.02671909, 0.0]),
        ("H", [-0.36300302, -0.51335996, -0.88916808]),
        ("N", [-0.48999652, -0.69296622, 1.20025051]),
        ("C", [0.27723429, -1.40935135, 1.96885228]),
        ("C", [-1.73024261, -0.64553106, 1.59030211]),
        ("C", [-1.83181119, -1.37666774, 2.70266390]),
        ("C", [-0.48010388, -1.89174235, 2.95868254]),
        ("H", [-0.16203259, -2.53181744, 3.78162527]),
        ("N", [-2.82193136, 0.08140799, 0.92646819]),
        ("O", [-3.86514425, 0.19068474, 1.52968085]),
        ("O", [-2.60509348, 0.52240282, -0.17924687]),
        ("H", [-2.71235228, -1.55082357, 3.28508735]),
    ];
    let atoms: Vec<&str> = rows.iter().map(|r| r.0).collect();
    let pos: Vec<Vec3> = rows.iter().map(|r| Vec3::from_array(r.1)).collect();
    let mut mol = molecule(&atoms, &pos);
    let result = complete(&mut mol).unwrap();
    assert_eq!(result.added, 1, "{result:?}");
    let h = mol.pos[14];
    assert!((h.distance(mol.pos[5]) - single_bond_length("C", "H")).abs() < 1.0e-4);
    // The new H sits in the ring plane, pointing away from the ring.
    let ring_centre = [4, 5, 6, 7, 8].iter().map(|&i| mol.pos[i]).sum::<Vec3>() / 5.0;
    assert!((h - mol.pos[5]).dot(mol.pos[5] - ring_centre) > 0.0);
    // The single-atom tool agrees with the bulk one.
    let mut smart = molecule(&atoms, &pos);
    append_smart_hydrogen(&mut smart, 5).unwrap();
    assert_eq!(smart.pos[14], h);
    assert_eq!(complete(&mut mol).unwrap().added, 0);
}
