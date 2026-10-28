use super::*;

#[test]
fn unsupported_ferrocene_jobs_queue_the_quantum_chemistry_popup() {
    let pool = AsyncComputeTaskPool::get_or_init(bevy::tasks::TaskPool::default);
    let mol = crate::molecule::Molecule::from_xyz(include_str!("../../tests/fixtures/ferrocene_builder.xyz"));
    for job in [JobType::SinglePoint, JobType::Gradient, JobType::Optimize,
        JobType::Hessian, JobType::Frequencies, JobType::OptimizeThenFrequencies] {
        let atoms = mol.atoms.clone();
        let positions: Vec<_> = mol.pos.iter().flat_map(|p| [p.x as f64,p.y as f64,p.z as f64]).collect();
        let dir = std::env::temp_dir().join(format!("beavyr-dreiding-popup-{}-{}",std::process::id(),job.label()));
        let workdir = dir.clone();
        let worker = pool.spawn(async move {
            run_cancellable(QcProgram::Dreiding,None,job,Path::new(""),&workdir,&atoms,&positions,
                0,1,&MethodConfig::default(),ExcitedStateOptions::default(),false,None,
                &AtomicBool::new(false),&Mutex::new(None))
        });
        let start = Instant::now();
        while !worker.is_finished() && start.elapsed() < Duration::from_secs(10) {
            std::thread::sleep(Duration::from_millis(1));
        }
        assert!(worker.is_finished());
        let task = QcRunTask { task:Some(worker), run_dir:Some(dir.clone()), ..Default::default() };
        let mut app = App::new();
        app.insert_resource(task).insert_resource(mol.clone())
            .init_resource::<crate::qchem_panel::QcPanelState>()
            .init_resource::<crate::ui_layout::UiLayout>()
            .init_resource::<crate::trajectory::TrajectoryState>()
            .init_resource::<crate::settings::MolSettings>()
            .init_resource::<crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask>()
            .init_resource::<crate::qchem_interfaces::xtb_optimize::XtbPanelState>()
            .init_resource::<crate::qchem_interfaces::xtb_freq::XtbFreqPanelState>()
            .init_resource::<crate::qchem_interfaces::xtb_freq::XtbFrequencyTask>()
            .init_resource::<crate::uvvis::UvVisState>()
            .init_resource::<crate::orbitals::OrbitalState>()
            .add_message::<crate::events::MoleculeChanged>()
            .add_systems(Update,poll_qc_run);
        app.update();
        let task = app.world().resource::<QcRunTask>();
        let warning = task.dreiding_warning.as_deref().expect("unsupported DREIDING must queue a popup");
        assert!(warning.contains("xTB") && warning.contains("spin multiplicity"),"{warning}");
        assert!(warning.contains("Fe"),"{warning}");
        assert!(!app.world().resource::<crate::qchem_panel::QcPanelState>().summary_open);
        assert!(task.last.is_none());
        let _ = fs::remove_dir_all(dir);
    }
}

#[test]
fn calculator_energy_and_frequencies_exclude_dummy_atoms_consistently() {
    let path = std::env::temp_dir().join(format!("beavyr-dummy-calculator-{}", std::process::id()));
    let atoms = ["X".into(), "O".into(), "H".into(), "Xx".into(), "H".into()];
    let positions = [0.0,0.0,0.0, 0.0,0.0,0.0, 1.1,0.0,0.0, 3.0,0.0,0.0, -0.2,1.05,0.0];
    for job in [JobType::SinglePoint, JobType::Frequencies] {
        let output = run_cancellable(QcProgram::Dreiding, None, job, Path::new(""),
            &path, &atoms, &positions, 0, 1, &MethodConfig::default(),
            ExcitedStateOptions::default(), false, None, &AtomicBool::new(false), &Mutex::new(None)).unwrap();
        if job == JobType::SinglePoint {
            assert!(output.energetics.final_energy_hartree.unwrap().is_finite());
        } else {
            let text = fs::read_to_string(output.hessian_path.unwrap()).unwrap();
            assert_eq!(text.lines().filter(|line| !line.starts_with('$')).count(), 9);
        }
    }
    let _ = fs::remove_dir_all(path);
}

#[test]
fn dreiding_dispatch_forwards_the_cancel_request() {
    let path = std::env::temp_dir().join(format!(
        "beavyr-dreiding-dispatch-cancel-{}",
        std::process::id()
    ));
    let atoms = ["O".into(), "H".into(), "H".into()];
    let positions = [0.0, 0.0, 0.0, 1.3, 0.0, 0.0, 0.0, 1.3, 0.0];
    let result = run_cancellable(
        QcProgram::Dreiding,
        None,
        JobType::OptimizeThenFrequencies,
        Path::new(""),
        &path,
        &atoms,
        &positions,
        0,
        1,
        &MethodConfig::default(),
        ExcitedStateOptions::default(),
        false,
        None,
        &AtomicBool::new(true),
        &Mutex::new(None),
    );
    assert_eq!(result.err().as_deref(), Some("Cancelled."));
    let _ = fs::remove_dir_all(path);
}

#[test]
fn cancelling_a_completed_dreiding_worker_prevents_result_application() {
    use crate::molecule::Molecule;
    use crate::trajectory::TrajectoryState;
    let pool = AsyncComputeTaskPool::get_or_init(bevy::tasks::TaskPool::default);
    let output = QcRunOutput {
        job: JobType::Optimize,
        program: QcProgram::Dreiding,
        log: "completed before Cancel".into(),
        trajectory_text: Some("3\ncycle 1 E = 0.01 Eh\nO 0 0 0\nH 1 0 0\nH 0 1 0\n".into()),
        hessian_path: None,
        molden_path: None,
        saveable: vec![],
        kept_dir: None,
        energetics: Energetics {
            final_energy_hartree: Some(0.01),
            ..Default::default()
        },
        scf: vec![],
        converged: Some(true),
    };
    let worker = pool.spawn(async move { Ok(output) });
    let start = Instant::now();
    while !worker.is_finished() && start.elapsed() < Duration::from_secs(10) {
        std::thread::sleep(Duration::from_millis(1));
    }
    assert!(worker.is_finished());

    let mut task = QcRunTask {
        task: Some(worker),
        cancel: Some(Arc::new(AtomicBool::new(false))),
        ..Default::default()
    };
    task.cancel();
    let original = Molecule::from_xyz("1\ncurrent structure\nC 8 9 10\n");
    let mut app = App::new();
    app.insert_resource(task)
        .insert_resource(original.clone())
        .init_resource::<crate::qchem_panel::QcPanelState>()
        .init_resource::<crate::ui_layout::UiLayout>()
        .init_resource::<TrajectoryState>()
        .init_resource::<crate::settings::MolSettings>()
        .init_resource::<crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask>()
        .init_resource::<crate::qchem_interfaces::xtb_optimize::XtbPanelState>()
        .init_resource::<crate::qchem_interfaces::xtb_freq::XtbFreqPanelState>()
        .init_resource::<crate::qchem_interfaces::xtb_freq::XtbFrequencyTask>()
        .init_resource::<crate::uvvis::UvVisState>()
        .init_resource::<crate::orbitals::OrbitalState>()
        .add_message::<crate::events::MoleculeChanged>()
        .add_systems(Update, poll_qc_run);
    app.update();
    let world = app.world();
    assert_eq!(world.resource::<Molecule>().atoms, original.atoms);
    assert_eq!(world.resource::<Molecule>().pos, original.pos);
    assert!(world.resource::<TrajectoryState>().frames.is_empty());
    assert!(
        !world
            .resource::<crate::qchem_panel::QcPanelState>()
            .summary_open
    );
    let task = world.resource::<QcRunTask>();
    assert!(!task.is_running());
    assert!(task.last.is_none());
    assert_eq!(task.last_message.as_deref(), Some("Cancelled."));
}
