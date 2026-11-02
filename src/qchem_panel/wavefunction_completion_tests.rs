use super::*;
use crate::{
    molecule::Molecule,
    orbitals::{molden::parse_molden, OrbitalState},
    trajectory::TrajectoryState,
};
const MOLDEN: &str = include_str!("../../tests/fixtures/tracking_water_a.molden");
struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "beavyr-qc-wave-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&p).unwrap();
        Self(p)
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn geometry() -> Molecule {
    let d = parse_molden(MOLDEN).unwrap();
    Molecule {
        atoms: d.atoms,
        pos: d.positions,
        ..Molecule::empty()
    }
}
fn output(dir: &Path, job: JobType) -> QcRunOutput {
    let path = dir.join("wave.molden");
    fs::write(&path, MOLDEN).unwrap();
    QcRunOutput {
        scan: None,
        job,
        program: QcProgram::PySCF,
        log: "finished".into(),
        trajectory_text: None,
        hessian_path: None,
        molden_path: Some(path),
        saveable: vec![],
        kept_dir: Some(dir.to_owned()),
        energetics: Default::default(),
        scf: vec![],
        converged: Some(true),
    }
}
fn complete(output: QcRunOutput, mol: Molecule, requested: bool, cancelled: bool) -> App {
    let pool = AsyncComputeTaskPool::get_or_init(bevy::tasks::TaskPool::default);
    let worker = pool.spawn(async move { Ok(output) });
    let start = Instant::now();
    while !worker.is_finished() && start.elapsed() < Duration::from_secs(5) {
        std::thread::sleep(Duration::from_millis(1));
    }
    assert!(worker.is_finished());
    let task = QcRunTask {
        task: Some(worker),
        wavefunction_requested: requested,
        cancel: Some(Arc::new(AtomicBool::new(cancelled))),
        ..Default::default()
    };
    let mut app = App::new();
    app.insert_resource(task)
        .insert_resource(mol)
        .init_resource::<crate::qchem_panel::QcPanelState>()
        .init_resource::<crate::ui_layout::UiLayout>()
        .init_resource::<TrajectoryState>()
        .init_resource::<crate::settings::MolSettings>()
        .init_resource::<crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask>()
        .init_resource::<crate::qchem_interfaces::xtb_optimize::XtbPanelState>()
        .init_resource::<crate::qchem_interfaces::xtb_freq::XtbFreqPanelState>()
        .init_resource::<crate::qchem_interfaces::xtb_freq::XtbFrequencyTask>()
        .init_resource::<crate::uvvis::UvVisState>()
        .init_resource::<OrbitalState>()
        .init_resource::<crate::trajectory_wavefunctions::State>()
        .add_message::<crate::events::MoleculeChanged>()
        .add_systems(Update, poll_qc_run);
    app.world_mut()
        .resource_mut::<crate::trajectory_wavefunctions::State>()
        .follow = true;
    app.world_mut().resource_mut::<OrbitalState>().show_surface = false;
    app.update();
    app
}
#[test]
fn single_point_opens_surface_and_replaces_old_trajectory_follow() {
    let dir = Directory::new();
    let app = complete(
        output(&dir.0, JobType::SinglePoint),
        geometry(),
        true,
        false,
    );
    let w = app.world();
    let orbital = w.resource::<OrbitalState>();
    assert!(orbital.data.is_some() && orbital.show_surface && orbital.selected.is_some());
    assert!(
        w.resource::<crate::ui_layout::UiLayout>().open[crate::ui_layout::Tab::Surface.index()]
    );
    assert!(
        !w.resource::<crate::trajectory_wavefunctions::State>()
            .follow
    );
    assert!(!w.resource::<QcRunTask>().last_is_error);
}
#[test]
fn optimization_wavefunction_keeps_the_final_matching_geometry() {
    for job in [JobType::Optimize, JobType::OptimizeThenFrequencies] {
        let dir = Directory::new();
        let final_mol = geometry();
        let mut initial = final_mol.clone();
        initial.pos[1].x += 0.2;
        let frame = |m: &Molecule| crate::trajectory::TrajectoryFrame {
            atoms: m.atoms.clone(),
            pos: m.pos.clone(),
            metadata: Default::default(),
        };
        let frames = vec![frame(&initial), frame(&final_mol)];
        let mut out = output(&dir.0, job);
        out.trajectory_text = Some(crate::trajectory_data::xyz(&frames, [0, 1]));
        let app = complete(out, initial, true, false);
        let w = app.world();
        assert_eq!(w.resource::<TrajectoryState>().current_frame, 1);
        assert!(!w.resource::<TrajectoryState>().playing);
        assert!(wavefunction_output::matches(
            w.resource::<OrbitalState>().data.as_ref().unwrap(),
            &w.resource::<Molecule>().atoms,
            &w.resource::<Molecule>().pos
        ));
    }
}
#[test]
fn missing_or_invalid_wavefunctions_are_reported_and_changed_structure_is_preserved() {
    let dir = Directory::new();
    let mut out = output(&dir.0, JobType::SinglePoint);
    out.molden_path = None;
    let app = complete(out, geometry(), true, false);
    assert!(app.world().resource::<QcRunTask>().last_is_error);
    let out = output(&dir.0, JobType::SinglePoint);
    fs::write(out.molden_path.as_ref().unwrap(), "broken").unwrap();
    let app = complete(out, geometry(), true, false);
    assert!(app.world().resource::<QcRunTask>().last_is_error);
    let mut edited = geometry();
    edited.pos[0].x += 1.;
    let app = complete(
        output(&dir.0, JobType::SinglePoint),
        edited.clone(),
        true,
        false,
    );
    assert_eq!(app.world().resource::<Molecule>().pos, edited.pos);
    assert!(app.world().resource::<OrbitalState>().data.is_none());
}
#[test]
fn cancelled_completion_does_not_load_wavefunction() {
    let dir = Directory::new();
    let app = complete(output(&dir.0, JobType::SinglePoint), geometry(), true, true);
    assert!(app.world().resource::<OrbitalState>().data.is_none());
    assert!(app.world().resource::<QcRunTask>().last.is_none());
}
#[test]
fn xtb_and_behemoth_request_their_native_molden_flag() {
    for (p, flag) in [
        (QcProgram::Xtb, "--molden"),
        (QcProgram::Behemoth, "--wf2molden"),
    ] {
        let mut c = std::process::Command::new("unused");
        wavefunction_output::request(&mut c, p);
        assert_eq!(
            c.get_args().collect::<Vec<_>>(),
            [std::ffi::OsStr::new(flag)]
        );
    }
}
#[test]
#[ignore = "requires BEAVYR_TEST_XTB pointing to an installed executable"]
fn native_xtb_wavefunctions_load_for_energy_and_optimization() {
    let binary = PathBuf::from(std::env::var("BEAVYR_TEST_XTB").unwrap());
    for job in [JobType::SinglePoint, JobType::Optimize] {
        let dir = Directory::new();
        let m = geometry();
        let coords = m
            .pos
            .iter()
            .flat_map(|p| p.to_array().map(f64::from))
            .collect::<Vec<_>>();
        let mut config = MethodConfig::default();
        config.nproc = 1;
        let out = run_cancellable(
            QcProgram::Xtb,
            None,
            job,
            &binary,
            &dir.0,
            &m.atoms,
            &coords,
            0,
            1,
            &config,
            Default::default(),
            true,
            None,
            &AtomicBool::new(false),
            &Mutex::new(None),
        )
        .unwrap();
        assert!(out.molden_path.as_ref().unwrap().is_file());
        let app = complete(out, m, true, false);
        assert!(
            !app.world().resource::<QcRunTask>().last_is_error,
            "{:?}",
            app.world().resource::<QcRunTask>().last_message
        );
        assert!(app.world().resource::<OrbitalState>().data.is_some());
        assert!(!app.world().resource::<TrajectoryState>().playing);
    }
}
