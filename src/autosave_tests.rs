use super::*;
use bevy::ecs::system::RunSystemOnce;
struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "beavyr-recovery-test-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn snapshot() -> Snapshot {
    let mut world = super::super::tests::world();
    world.insert_resource(crate::molecule::Molecule::from_xyz(
        "2\nH2\nH 0 0 0\nH 0 0 0.74\n",
    ));
    world
        .resource_mut::<crate::measurements::Measurements>()
        .new_pair(0, 1, 0.74);
    Snapshot::capture(&world, Workspace::capture(&egui::Context::default())).unwrap()
}
#[test]
fn recovery_retains_three_generations_and_excludes_live_sessions() {
    let root = Directory::new();
    let Completed::Initialized(session, _, _) = initialize(&root.0).unwrap() else {
        panic!()
    };
    for generation in 1..=5 {
        save(&session.directory, generation, snapshot()).unwrap();
    }
    assert!(discover(&root.0).unwrap().is_empty());
    let directory = session.directory.clone();
    drop(session);
    let entries = discover(&root.0).unwrap();
    assert_eq!(entries.len(), 3);
    assert!(!directory.join("recovery-00000002.beavyr").exists());
    for entry in entries {
        let s: Snapshot = ron::from_str(&fs::read_to_string(entry.path).unwrap()).unwrap();
        s.validate().unwrap();
    }
}
#[test]
fn corrupted_recovery_never_replaces_editor_and_valid_restore_is_undoable() {
    let root = Directory::new();
    let path = root.0.join("restore.beavyr");
    let mut world = super::super::tests::world();
    world.insert_resource(crate::molecule::Molecule::from_xyz(
        "1\noriginal\nO 1 2 3\n",
    ));
    fs::write(&path, "invalid").unwrap();
    assert!(recover(
        &mut world,
        &path,
        Workspace::capture(&egui::Context::default())
    )
    .is_err());
    assert_eq!(world.resource::<crate::molecule::Molecule>().atoms, ["O"]);
    fs::write(&path, ron::ser::to_string(&snapshot()).unwrap()).unwrap();
    recover(
        &mut world,
        &path,
        Workspace::capture(&egui::Context::default()),
    )
    .unwrap();
    world
        .run_system_once(super::super::restore_measurements)
        .unwrap();
    assert_eq!(
        world
            .resource::<crate::measurements::Measurements>()
            .pairs
            .len(),
        1
    );
    assert_eq!(
        world.resource::<crate::molecule::Molecule>().atoms,
        ["H", "H"]
    );
    assert!(world.resource::<ProjectState>().path.is_none());
    assert_eq!(
        world
            .resource::<ProjectState>()
            .previous
            .as_ref()
            .unwrap()
            .atoms,
        ["O"]
    );
    assert!(!world
        .resource::<crate::qchem_panel::run::QcRunTask>()
        .is_running());
}
#[test]
fn autosave_does_not_touch_manual_projects() {
    let root = Directory::new();
    let manual = root.0.join("manual.beavyr");
    fs::write(&manual, b"keep").unwrap();
    let Completed::Initialized(session, _, _) = initialize(&root.0).unwrap() else {
        panic!()
    };
    save(&session.directory, 1, snapshot()).unwrap();
    assert_eq!(fs::read(manual).unwrap(), b"keep");
}

#[test]
fn recovery_preserves_completed_calculation_and_atomic_analysis_without_starting_jobs() {
    let mut world = super::super::tests::world();
    let text = include_str!("../tests/fixtures/behemoth_ch2o_canonicalMO.molden");
    let data = crate::orbitals::molden::parse_molden(text).unwrap();
    let frame = crate::trajectory::TrajectoryFrame {
        atoms: data.atoms.clone(),
        pos: data.positions.clone(),
        metadata: Default::default(),
    };
    let polar_values = crate::polarizability::data::compute(&data).unwrap();
    let expected = polar_values.alpha.clone();
    let mut polar = crate::polarizability::State::default();
    polar.results = Some(std::sync::Arc::new(crate::polarizability::Results {
        frames: vec![frame.clone()],
        values: vec![polar_values],
        method: "test".into(),
    }));
    world.insert_resource(polar);
    let f = crate::reactivity::data::frontier(
        frame,
        crate::reactivity::data::Wave::Embedded(text.into()),
        crate::orbitals::molden::Spin::Alpha,
    )
    .unwrap();
    let mut react = crate::reactivity::State::default();
    react.results = Some(std::sync::Arc::new(crate::reactivity::data::Results {
        version: 1,
        method: "test".into(),
        frames: vec![f],
    }));
    world.insert_resource(react);
    world
        .resource_mut::<crate::qchem_panel::run::QcRunTask>()
        .last = Some(crate::qchem_panel::run::QcRunOutput {
        scan: None,
        job: crate::qchem_interfaces::job::JobType::SinglePoint,
        program: crate::qchem_interfaces::program::QcProgram::PySCF,
        log: "completed output".into(),
        trajectory_text: None,
        hessian_path: None,
        molden_path: None,
        saveable: vec![],
        kept_dir: None,
        energetics: Default::default(),
        scf: vec![],
        converged: Some(true),
    });
    let snapshot =
        Snapshot::capture(&world, Workspace::capture(&egui::Context::default())).unwrap();
    let serialized = ron::ser::to_string(&snapshot).unwrap();
    let saved: Snapshot = ron::from_str(&serialized).unwrap();
    let orbitals = saved.validate().unwrap();
    saved.apply(&mut world, orbitals);
    assert_eq!(
        world
            .resource::<crate::polarizability::State>()
            .results
            .as_ref()
            .unwrap()
            .values[0]
            .alpha,
        expected
    );
    assert_eq!(
        world
            .resource::<crate::reactivity::State>()
            .results
            .as_ref()
            .unwrap()
            .frames
            .len(),
        1
    );
    assert!(!world.resource::<crate::reactivity::State>().is_running());
    let qc = world.resource::<crate::qchem_panel::run::QcRunTask>();
    assert!(!qc.is_running());
    assert_eq!(qc.last.as_ref().unwrap().log, "completed output");
}
