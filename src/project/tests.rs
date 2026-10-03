use super::*;

fn world() -> World {
    let mut world = World::new();
    world.insert_resource(Molecule::empty());
    world.init_resource::<XyzBuffer>();
    world.init_resource::<MolSettings>();
    world.init_resource::<OrbitCamera>();
    world.init_resource::<TrajectoryState>();
    world.init_resource::<QcPanelState>();
    world.init_resource::<crate::qchem_interfaces::xtb_optimize::XtbPanelState>();
    world.init_resource::<XtbFrequencyTask>();
    world.init_resource::<XtbFreqPanelState>();
    world.init_resource::<UvVisState>();
    world.init_resource::<OrbitalState>();
    world.init_resource::<UiLayout>();
    world.init_resource::<ProjectState>();
    world.init_resource::<crate::qchem_panel::run::QcRunTask>();
    world.init_resource::<crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask>();
    world.init_resource::<crate::uvvis::run::SpectrumTask>();
    world.init_resource::<crate::conformer::ConformerRun>();
    world.init_resource::<crate::ts_generation::TsGeneration>();
    world.init_resource::<crate::measurements::Measurements>();
    world.init_resource::<Messages<crate::events::MoleculeChanged>>();
    world
}

#[test]
fn portable_project_restores_embedded_data_settings_and_camera() {
    let mut world = world();
    let text = include_str!("../../tests/fixtures/behemoth_ch2o_canonicalMO.molden");
    let data = crate::orbitals::molden::parse_molden(text).unwrap();
    world.insert_resource(Molecule {
        atoms: data.atoms.clone(),
        pos: data.positions.clone(),
        bonds: vec![],
        hydrogen_bonds: vec![],
    });
    let mut orbitals = OrbitalState::default();
    orbitals.adopt(
        data,
        Some("canonicalMO.molden".into()),
        None,
        Some(text.into()),
    );
    orbitals.selected = Some(3);
    orbitals.isovalue = 0.11;
    world
        .resource_mut::<UvVisState>()
        .orbital_link
        .bind(&orbitals);
    world.insert_resource(orbitals);
    world.resource_mut::<UvVisState>().width_nm = 7.0;
    world.resource_mut::<QcPanelState>().advanced_text = "manual input stays editable".into();
    world.resource_mut::<QcPanelState>().config.nproc = 32;
    world.resource_mut::<OrbitCamera>().target = Vec3::new(1.0, 2.0, 3.0);
    let workspace = Workspace::capture(&egui::Context::default());
    let saved = Snapshot::capture(&world, workspace).unwrap();
    let text = ron::ser::to_string(&saved).unwrap();
    let loaded: Snapshot = ron::from_str(&text).unwrap();
    let orbitals = loaded.validate().unwrap();
    loaded.apply(&mut world, orbitals);
    assert_eq!(world.resource::<OrbitalState>().selected, Some(3));
    assert_eq!(world.resource::<OrbitalState>().isovalue, 0.11);
    assert!(world
        .resource::<UvVisState>()
        .orbital_link
        .matches(world.resource::<OrbitalState>()));
    assert_eq!(world.resource::<UvVisState>().width_nm, 7.0);
    assert_eq!(world.resource::<QcPanelState>().config.nproc, 4);
    assert_eq!(
        world.resource::<QcPanelState>().advanced_text,
        "manual input stays editable"
    );
    assert_eq!(
        world.resource::<OrbitCamera>().target,
        Vec3::new(1.0, 2.0, 3.0)
    );
}

#[test]
fn frequency_data_survives_without_source_files() {
    let mut world = world();
    let mut panel = XtbFreqPanelState::default();
    let mut task = XtbFrequencyTask::default();
    let geometry = include_str!("../../tests/fixtures/h2o2_hessian_geometry.xyz");
    let (_, atoms, coords) = crate::molecule::parse_xyz_angstrom(geometry);
    let coords: Vec<f64> = coords.into_iter().flatten().collect();
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/h2o2_hessian");
    assert!(
        crate::qchem_interfaces::xtb_freq::load_bare_hessian_from_path(
            &path, &atoms, &coords, &mut panel, &mut task
        )
    );
    let frequencies = task.result.as_ref().unwrap().frequencies_cm1.clone();
    let archive = task.archive();
    let text = ron::ser::to_string(&archive).unwrap();
    let archive: FrequencyArchive = ron::from_str(&text).unwrap();
    archive.validate().unwrap();
    world
        .resource_mut::<XtbFrequencyTask>()
        .restore_archive(archive);
    assert_eq!(
        world
            .resource::<XtbFrequencyTask>()
            .result
            .as_ref()
            .unwrap()
            .frequencies_cm1,
        frequencies
    );
}

#[test]
fn corrupt_or_future_projects_are_rejected_before_changes() {
    let world = world();
    let mut snapshot =
        Snapshot::capture(&world, Workspace::capture(&egui::Context::default())).unwrap();
    snapshot.version += 1;
    assert!(snapshot.validate().is_err());
    snapshot.version = VERSION;
    snapshot.atoms.push("C".into());
    assert!(snapshot.validate().is_err());
    assert!(ron::from_str::<Snapshot>("broken project").is_err());
    assert!(world.resource::<Molecule>().atoms.is_empty());
}

#[test]
fn failed_atomic_save_preserves_existing_destination() {
    let path = std::env::temp_dir().join(format!("beavyr-project-test-{}", std::process::id()));
    fs::create_dir_all(&path).unwrap();
    let file = path.join("project.beavyr");
    atomic_write(&file, b"first").unwrap();
    atomic_write(&file, b"second").unwrap();
    assert_eq!(fs::read(&file).unwrap(), b"second");
    assert!(atomic_write(&path, b"cannot replace a directory").is_err());
    assert_eq!(fs::read(&file).unwrap(), b"second");
    fs::remove_dir_all(path).unwrap();
}

#[test]
fn opening_invalid_data_keeps_the_current_session_and_valid_open_is_reversible() {
    let mut world = world();
    world.insert_resource(Molecule::from_xyz("1\noriginal\nC 1 2 3\n"));
    let ctx = egui::Context::default();
    let path =
        std::env::temp_dir().join(format!("beavyr-open-project-{}.beavyr", std::process::id()));
    fs::write(&path, "not a project").unwrap();
    world.resource_mut::<ProjectState>().request =
        Some(Request::Open(path.clone(), Workspace::capture(&ctx)));
    process(&mut world);
    assert!(world.resource::<ProjectState>().error);
    assert_eq!(world.resource::<Molecule>().atoms, ["C"]);
    let mut replacement = Snapshot::capture(&world, Workspace::capture(&ctx)).unwrap();
    replacement.atoms = vec!["O".into()];
    replacement.positions = vec![Vec3::ZERO];
    fs::write(&path, ron::ser::to_string(&replacement).unwrap()).unwrap();
    world.resource_mut::<ProjectState>().request =
        Some(Request::Open(path.clone(), Workspace::capture(&ctx)));
    process(&mut world);
    assert!(
        !world.resource::<ProjectState>().error,
        "{:?}",
        world.resource::<ProjectState>().message
    );
    assert_eq!(world.resource::<Molecule>().atoms, ["O"]);
    world.resource_mut::<ProjectState>().request = Some(Request::Restore(Workspace::capture(&ctx)));
    process(&mut world);
    assert_eq!(world.resource::<Molecule>().atoms, ["C"]);
    assert_eq!(world.resource::<Molecule>().pos, [Vec3::new(1.0, 2.0, 3.0)]);
    fs::remove_file(path).unwrap();
}

#[test]
fn workspace_round_trip_preserves_window_position_and_docking() {
    let ctx = egui::Context::default();
    let mut window_id = egui::Id::NULL;
    for _ in 0..2 {
        let mut output = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1000.0, 800.0),
                )),
                ..Default::default()
            },
            |ui| {
                let response = egui::Window::new("Project layout test")
                    .default_pos([120.0, 180.0])
                    .show(ui.ctx(), |ui| {
                        ui.label("Content");
                    })
                    .unwrap();
                window_id = response.response.layer_id.id;
            },
        );
        output.textures_delta.clear();
    }
    crate::molecule_builder::builder_ui::restore_editor_docking(&ctx, true);
    let position = ctx.memory(|m| m.area_rect(window_id));
    assert!(position.is_some());
    let workspace: Workspace =
        ron::from_str(&ron::ser::to_string(&Workspace::capture(&ctx)).unwrap()).unwrap();
    let restored = egui::Context::default();
    workspace.apply(&restored);
    // egui recomputes area extents on its next layout pass; position is persistent.
    assert_eq!(
        restored.memory(|m| m.area_rect(window_id)).map(|r| r.min),
        position.map(|r| r.min)
    );
    assert!(crate::molecule_builder::builder_ui::editor_docked(
        &restored
    ));
    let saved = ron::ser::to_string(&Workspace::capture(&ctx)).unwrap();
    let mut workspace: Option<Workspace> = Some(ron::from_str(&saved).unwrap());
    let mut output = restored.run_ui(egui::RawInput::default(), |ui| {
        if let Some(workspace) = workspace.take() {
            workspace.apply(ui.ctx());
        }
        let mut text = String::new();
        ui.text_edit_singleline(&mut text).request_focus();
    });
    output.textures_delta.clear();
}
