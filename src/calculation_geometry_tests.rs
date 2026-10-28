use super::*;
use crate::conformer::{constraints::Constraint, ConformerSettings, SearchEngine, Thoroughness};
use crate::molecule_builder::{
    builder_ui::{EditorRotateState, ZMatrixBuilderState},
    history::BuilderStep,
    topology::Topology,
    zmat2xyz::xyz_to_zmat,
};
use crate::optimizer::constrained::ConstraintCoordinate;
use bevy::prelude::Vec3;
use std::sync::atomic::AtomicBool;

fn water() -> Molecule {
    Molecule::from_xyz("5\n\nX 0 0 0\nO 0 0 0\nH 1.1 0 0\nXx 3 0 0\nH -0.2 1.05 0\n")
}

#[test]
fn dummy_removal_preserves_real_coordinates_and_explicit_chemistry() {
    let mut mol = water();
    let mut topology = Topology::from_molecule(&mol);
    topology.bonds = vec![(1, 2, 1), (1, 4, 1), (0, 1, 1)];
    for atom in &mut topology.atoms {
        atom.charge = Some(0);
    }
    mol.topology = Some(topology);
    let (real, mapping) = prepare(&mol).unwrap();
    assert_eq!(mapping, [None, Some(0), Some(1), None, Some(2)]);
    assert_eq!(real.atoms, ["O", "H", "H"]);
    assert_eq!(real.pos, [mol.pos[1], mol.pos[2], mol.pos[4]]);
    let topology = real.topology.unwrap();
    assert!(topology.valid_for(&real.atoms));
    assert_eq!(topology.bonds, [(0, 1, 1), (0, 2, 1)]);
    assert_eq!(
        topology.atoms.iter().map(|a| a.id).collect::<Vec<_>>(),
        [2, 3, 5]
    );
    assert_eq!(topology.charge(), Some(0));
}

#[test]
fn removal_updates_display_indices_and_has_an_undo_snapshot() {
    let ctx = egui::Context::default();
    let mut mol = water();
    let original = mol.clone();
    let mut builder = ZMatrixBuilderState::default();
    builder.zmat = xyz_to_zmat(&mol.atoms, &mol.pos);
    builder.selection.sync(&mol);
    builder.selection.frozen.insert(4);
    builder.selection.hidden.insert(2);
    let mut editor = EditorRotateState::default();
    let mut settings = crate::settings::MolSettings::default();
    let mut traj = crate::trajectory::TrajectoryState::default();
    let mut history = crate::structure_history::StructureHistory::detached();
    let mut conformer = crate::conformer::ConformerRun::default();
    conformer.settings.constraints = vec![Constraint {
        coordinate: ConstraintCoordinate::Bond([1, 4]),
        value: 2.0,
    }];
    conformer.trail.active = vec![crate::conformer::afir::ActiveBond { a: 1, b: 2 }];
    prepare_display(&ctx, &mut mol).unwrap();
    assert_eq!(mol.atoms, ["O", "H", "H"]);
    assert!(finish_display_removal(
        &ctx,
        &mol,
        &mut builder,
        &mut editor,
        &mut settings,
        &mut traj,
        &mut history,
        &mut conformer
    ));
    assert_eq!(builder.selection.frozen, [2].into());
    assert_eq!(builder.selection.hidden, [1].into());
    assert_eq!(
        conformer.settings.constraints[0].coordinate,
        ConstraintCoordinate::Bond([0, 2])
    );
    assert_eq!(
        conformer.trail.active[0],
        crate::conformer::afir::ActiveBond { a: 0, b: 1 }
    );
    assert!(builder.skip_next_sync);
    let before = builder
        .history
        .undo(BuilderStep {
            name: String::new(),
            atoms: mol.atoms.clone(),
            pos: mol.pos.clone(),
            topology: mol.topology.clone(),
            zmat: builder.zmat.clone(),
        })
        .unwrap();
    assert_eq!(before.atoms, original.atoms);
    assert_eq!(before.pos, original.pos);
    assert!(!finish_display_removal(
        &ctx,
        &mol,
        &mut builder,
        &mut editor,
        &mut settings,
        &mut traj,
        &mut history,
        &mut conformer
    ));
}

#[test]
fn malformed_and_all_dummy_inputs_are_rejected_before_editing() {
    let ctx = egui::Context::default();
    let mut mol = Molecule::from_xyz("1\n\nX 0 0 0\n");
    assert!(prepare_display(&ctx, &mut mol)
        .unwrap_err()
        .contains("no real atoms"));
    assert_eq!(mol.atoms, ["X"]);
    mol = water();
    mol.pos.pop();
    assert!(prepare(&mol).is_err());
    assert!(real_coordinates(&["C".into()], &[0.0, 0.0]).is_err());
    let (atoms, positions) = real_coordinates(
        &["X".into(), "C".into()],
        &[0.0, 0.0, 0.0, 1.234567890123, 2.0, 3.0],
    )
    .unwrap();
    assert_eq!(atoms, ["C"]);
    assert_eq!(positions, [1.234567890123, 2.0, 3.0]);
}

#[test]
fn ferrocene_cleanup_and_search_explain_the_same_dreiding_limitation() {
    let mol = Molecule::from_xyz(include_str!("../tests/fixtures/ferrocene_builder.xyz"));
    let cleanup = crate::molecule_builder::cleanup::relax(
        &mol,
        &crate::molecule_builder::selection::Selection::default(),
        crate::molecule_builder::cleanup::Scope::Whole,
        &AtomicBool::new(false),
    )
    .err()
    .unwrap()
    .to_string();
    let search = crate::conformer::search(&mol, ConformerSettings::default())
        .unwrap_err()
        .to_string();
    for text in [&cleanup, &search] {
        for needed in [
            "DREIDING cannot handle",
            "Quantum Chemistry",
            "xTB",
            "geometry optimization",
            "total charge",
            "spin multiplicity",
            "valence electrons not defined for element Fe",
        ] {
            assert!(text.contains(needed), "Missing {needed}: {text}");
        }
        assert!(!text.contains("not a recognised element"));
    }
    assert_eq!(cleanup, search);
}

#[test]
fn unsupported_other_elements_get_the_same_guidance() {
    let mol = Molecule::from_xyz("2\n\nLi 0 0 0\nF 1.564 0 0\n");
    let text = crate::conformer::search(&mol, ConformerSettings::default())
        .unwrap_err()
        .to_string();
    assert!(
        text.contains("xTB") && text.contains("spin multiplicity"),
        "{text}"
    );
}

#[test]
fn organic_conformer_search_works_with_interleaved_dummies() {
    let mut mol = Molecule::from_xyz("14\n\nC -1.92 0.26 0\nC -0.64 -0.57 0\nC 0.64 0.26 0\nC 1.92 -0.57 0\nH -2.82 -0.36 0\nH -1.96 0.9 0.89\nH -1.96 0.9 -0.89\nH -0.62 -1.22 0.89\nH -0.62 -1.22 -0.89\nH 0.62 0.91 0.89\nH 0.62 0.91 -0.89\nH 2.82 0.05 0\nH 1.96 -1.21 0.89\nH 1.96 -1.21 -0.89\n");
    mol.atoms.insert(0, "X".into());
    mol.pos.insert(0, Vec3::ZERO);
    mol.atoms.insert(4, "Xx".into());
    mol.pos.insert(4, Vec3::ZERO);
    mol.topology = None;
    mol.recompute_bonds(1.2, 2.5);
    let outcome = crate::conformer::search(
        &mol,
        ConformerSettings {
            thoroughness: Thoroughness::Quick,
            reproducible: true,
            ncore: 1,
            engine: SearchEngine::Dreiding,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(outcome.atoms.len(), 14);
    assert!(!outcome.atoms.iter().any(|a| crate::molecule::is_dummy(a)));
    assert!(!outcome.conformers.is_empty());
    assert!(outcome
        .conformers
        .iter()
        .all(|c| c.positions_angstrom.len() == 42));
}

#[test]
fn dummy_constraint_is_rejected_instead_of_becoming_a_different_atom() {
    let mol = water();
    let error = crate::conformer::search(
        &mol,
        ConformerSettings {
            constraints: vec![Constraint {
                coordinate: ConstraintCoordinate::Bond([0, 2]),
                value: 2.0,
            }],
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(matches!(
        error,
        crate::conformer::ConformerError::Constraints(_)
    ));
    assert!(error.to_string().contains("dummy atom"));
}
