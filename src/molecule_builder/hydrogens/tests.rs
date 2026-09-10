use super::*;

fn molecule(atoms: &[&str], pos: &[Vec3]) -> Molecule {
    Molecule {
        atoms: atoms.iter().map(|s| s.to_string()).collect(),
        pos: pos.to_vec(),
        bonds: vec![],
        hydrogen_bonds: vec![],
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
        &["N", "C", "H"],
        &[Vec3::new(20.0, 0.0, 0.0), Vec3::ZERO, Vec3::X * 1.07],
    );
    let before = mol.clone();
    assert_eq!(complete(&mut mol).unwrap().added, 3);
    assert_eq!(&mol.atoms[..3], &before.atoms);
    assert_eq!(&mol.pos[..3], &before.pos);
    assert!(mol.pos[3..].iter().all(|p| p.distance(mol.pos[0]) > 10.0));
}

#[test]
fn carbonyl_alkene_aromatic_and_partial_bonds_are_not_saturated() {
    for (symbols, distance) in [
        (["C", "O"], 1.22),
        (["C", "C"], 1.34),
        (["C", "C"], 1.40),
        (["C", "C"], 1.70),
    ] {
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
