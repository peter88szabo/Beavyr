use super::*;
use crate::trajectory::TrajectoryFrame;

struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        Self(
            crate::qchem_interfaces::xtb_optimize::create_xtb_run_dir(
                &std::env::temp_dir().join("beavyr_history_tests"),
            )
            .unwrap(),
        )
    }
    fn history(&self) -> StructureHistory {
        StructureHistory::open(Some(self.0.clone()))
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn molecule(distance: f32) -> Molecule {
    Molecule {
        atoms: vec!["C".into(), "H".into()],
        pos: vec![Vec3::ZERO, Vec3::new(distance, 1.2345678, -0.000000123)],
        bonds: vec![],
        hydrogen_bonds: vec![],
    }
}
fn remember(history: &mut StructureHistory, distance: f32, name: &str) {
    history.remember(Snapshot::capture(&molecule(distance)).unwrap(), name);
    assert!(!history.is_error, "{:?}", history.message);
}
#[test]
fn coordinates_and_names_survive_a_new_session() {
    let scratch = Scratch::new();
    let mut first = scratch.history();
    remember(&mut first, 1.1234567, "Drawn α structure");
    drop(first);
    let next_day = scratch.history();
    assert_eq!(next_day.entries.len(), 1);
    assert_eq!(next_day.entries[0].name, "Drawn α structure");
    assert_eq!(
        next_day.entries[0].snapshot,
        Snapshot::capture(&molecule(1.1234567)).unwrap()
    );
}
#[test]
fn separate_windows_merge_saved_structures() {
    let scratch = Scratch::new();
    let mut first = scratch.history();
    let mut second = scratch.history();
    remember(&mut first, 1.0, "window A");
    remember(&mut second, 2.0, "window B");
    first.refresh();
    assert_eq!(first.entries.len(), 2);
    assert_eq!(first.entries[0].name, "window B");
    assert_eq!(first.entries[1].name, "window A");
}
#[test]
fn history_keeps_only_twenty_distinct_structures_and_promotes_reloads() {
    let scratch = Scratch::new();
    let mut history = scratch.history();
    for i in 0..25 {
        remember(&mut history, i as f32, &format!("structure {i}"));
    }
    assert_eq!(history.entries.len(), LIMIT);
    assert_eq!(history.entries[0].name, "structure 24");
    assert_eq!(history.entries[19].name, "structure 5");
    remember(&mut history, 10.0, "revisited");
    assert_eq!(history.entries.len(), LIMIT);
    assert_eq!(history.entries[0].name, "revisited");
    assert_eq!(read_entries(&scratch.0).unwrap().0.len(), LIMIT);
    assert_eq!(
        history
            .entries
            .iter()
            .filter(|e| e.snapshot == Snapshot::capture(&molecule(10.0)).unwrap())
            .count(),
        1
    );
}
#[test]
fn replacing_or_clearing_preserves_the_latest_drawn_geometry() {
    let scratch = Scratch::new();
    let mut history = scratch.history();
    let traj = TrajectoryState::default();
    let mut current = molecule(1.0);
    // Arbitrary editing creates no entries until a boundary or explicit save.
    current.pos[1].x = 1.7;
    assert!(history.entries.is_empty());
    let previous = history.before(&current, &traj);
    current = molecule(2.0);
    history.replaced(previous, &current, &traj, "loaded.xyz");
    assert_eq!(history.entries.len(), 2);
    assert_eq!(history.entries[1].snapshot.positions[1].x, 1.7);
    current.pos[1].x = 2.8;
    let previous = history.before(&current, &traj);
    current = Molecule::empty();
    history.replaced(previous, &current, &traj, "Drawn structure");
    assert_eq!(history.entries.len(), 3);
    assert_eq!(history.entries[0].snapshot.positions[1].x, 2.8);
    assert!(history.entries.iter().all(|e| !e.snapshot.atoms.is_empty()));
}
#[test]
fn trajectory_loading_and_playback_do_not_fill_the_history() {
    let scratch = Scratch::new();
    let mut history = scratch.history();
    let mut traj = TrajectoryState::default();
    let mut current = molecule(1.0);
    let previous = history.before(&current, &traj);
    traj.frames = (0..100)
        .map(|i| TrajectoryFrame {
            atoms: current.atoms.clone(),
            pos: molecule(i as f32).pos,
        })
        .collect();
    current = molecule(0.0);
    history.replaced(previous, &current, &traj, "trajectory.xyz");
    assert_eq!(history.entries.len(), 1); // outgoing standalone structure
    for i in 0..100 {
        current.pos.clone_from(&traj.frames[i].pos);
        assert!(history.before(&current, &traj).is_none());
    }
    let previous = history.before(&current, &traj);
    traj.clear_for_structure();
    current = molecule(200.0);
    history.replaced(previous, &current, &traj, "next.xyz");
    assert_eq!(history.entries.len(), 2);
}
#[test]
fn manual_snapshot_is_independent_and_restore_stops_old_playback() {
    let scratch = Scratch::new();
    let mut history = scratch.history();
    let mut mol = molecule(1.0);
    history.save_name = "My edited structure".into();
    history.save_current(&mol);
    let saved = history.entries[0].clone();
    mol.pos[1].x = 4.0;
    assert_eq!(saved.snapshot.positions[1].x, 1.0);
    let mut traj = TrajectoryState {
        playing: true,
        overlay_enabled: true,
        fixed_bonds: true,
        ..Default::default()
    };
    history.restore(saved, &mut mol, &mut traj);
    assert_eq!(mol.pos[1].x, 1.0);
    assert!(!traj.playing);
    assert!(!traj.overlay_enabled);
    assert!(!traj.fixed_bonds);
    assert!(traj.frames.is_empty());
}
#[test]
fn partial_or_corrupt_files_do_not_hide_valid_history() {
    let scratch = Scratch::new();
    let mut history = scratch.history();
    remember(&mut history, 1.0, "good");
    fs::write(scratch.0.join("structure-incomplete.tmp"), "2\npartial").unwrap();
    fs::write(scratch.0.join("structure-bad.xyz"), "not an XYZ").unwrap();
    let reloaded = scratch.history();
    assert_eq!(reloaded.entries.len(), 1);
    assert!(reloaded.message.unwrap().contains("unreadable"));
    assert!(scratch.0.join("structure-bad.xyz").exists());
}
#[test]
fn empty_or_nonfinite_structures_are_not_saved() {
    let scratch = Scratch::new();
    let mut history = scratch.history();
    history.save_current(&Molecule::empty());
    assert!(history.is_error);
    let mut invalid = molecule(1.0);
    invalid.pos[0].x = f32::NAN;
    history.save_current(&invalid);
    assert!(history.is_error);
    assert!(history.entries.is_empty());
}

#[test]
fn simultaneous_windows_preserve_every_save_below_the_limit() {
    let scratch = Scratch::new();
    let mut workers = Vec::new();
    for window in 0..4 {
        let directory = scratch.0.clone();
        workers.push(std::thread::spawn(move || {
            let mut history = StructureHistory::open(Some(directory));
            for i in 0..4 {
                remember(
                    &mut history,
                    (window * 10 + i) as f32,
                    &format!("window {window} entry {i}"),
                );
            }
        }));
    }
    for worker in workers {
        worker.join().unwrap();
    }
    let history = scratch.history();
    assert_eq!(history.entries.len(), 16);
    assert_eq!(read_entries(&scratch.0).unwrap().0.len(), 16);
    for window in 0..4 {
        assert!(history
            .entries
            .iter()
            .any(|entry| entry.name == format!("window {window} entry 3")));
    }
}

#[test]
fn standalone_xyz_loading_replaces_old_playback_in_the_real_loader() {
    use crate::events::MoleculeChanged;
    use crate::settings::MolSettings;
    use bevy::ecs::system::SystemState;
    let scratch = Scratch::new();
    let path = scratch.0.join("single.xyz");
    fs::write(&path, "1\ncarbon\nC 3 2 1\n").unwrap();
    let mut world = World::new();
    world.init_resource::<Messages<MoleculeChanged>>();
    let mut system: SystemState<MessageWriter<MoleculeChanged>> = SystemState::new(&mut world);
    let mut writer = system.get_mut(&mut world).unwrap();
    let old = molecule(1.0);
    let mut traj = TrajectoryState {
        frames: vec![
            TrajectoryFrame {
                atoms: old.atoms.clone(),
                pos: old.pos.clone()
            };
            2
        ],
        playing: true,
        overlay_enabled: true,
        ..Default::default()
    };
    let mut mol = old;
    let mut buffer = crate::ui::XyzBuffer::default();
    crate::ui::load_xyz_from_path(
        &path,
        &mut buffer,
        &mut traj,
        &mut mol,
        &mut MolSettings::default(),
        &mut crate::camera::OrbitCamera::default(),
        &mut writer,
    )
    .unwrap();
    assert_eq!(mol.atoms, ["C"]);
    assert!(!traj.playing);
    assert!(!traj.overlay_enabled);
    assert!(traj.frames.is_empty());
}
