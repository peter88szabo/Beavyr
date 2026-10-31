use super::*;
#[test]
fn geometry_search_defaults_match_behemoth_cli() {
    use crate::optimizer::transition_state::*;
    let o = search::Options::default();
    assert_eq!(
        o.ts.coordinate_mode,
        TransitionStateCoordinateMode::Internal
    );
    assert_eq!(
        o.ts.ts_initial_hessian,
        TransitionStateInitialHessian::Exact
    );
    assert_eq!(
        o.ts.mode_tracking_coordinates,
        TransitionStateModeTracking::Internal
    );
    assert_eq!(o.ts.hessian_recalc_interval, 0);
    assert_eq!(o.ts.max_cycles, 100);
    assert!(o.ts.project_eckart && o.ts.final_hessian);
    assert_eq!(o.minimum.coordinates, Coordinates::Primitive);
    assert_eq!(o.minimum.max_cycles, 300);
}
fn scratch(label: &str) -> PathBuf {
    let p = std::env::temp_dir().join(format!(
        "beavyr-geometry-search-{label}-{}",
        std::process::id()
    ));
    fs::create_dir_all(&p).unwrap();
    p
}
fn model(n: usize) -> ConnectivityModel {
    ConnectivityModel {
        rcov_disp: vec![1.; n],
        kcn: 7.5,
        facmin: 0.6,
    }
}
struct Quadratic(f64);
impl Objective for Quadratic {
    fn energy(&mut self, x: &[f64]) -> anyhow::Result<f64> {
        Ok(0.5 * (self.0 * x[0] * x[0] + x[1..].iter().map(|v| v * v).sum::<f64>()))
    }
    fn gradient(&mut self, x: &[f64], g: &mut [f64]) -> anyhow::Result<()> {
        g.copy_from_slice(x);
        g[0] *= self.0;
        Ok(())
    }
}
#[test]
fn geometry_search_saddle_verification_and_step_limit() {
    use crate::optimizer::transition_state::*;
    let mut settings = search::Options::default();
    settings.ts.coordinate_mode = TransitionStateCoordinateMode::Cartesian;
    settings.ts.project_eckart = false;
    let dir = scratch("analytic-ts");
    let atoms = vec!["He".into()];
    let r = compute(
        &mut Quadratic(-1.),
        vec![0.2, 0.1, -0.1],
        model(1),
        &atoms,
        &settings,
        JobType::TransitionState,
        &dir,
    )
    .unwrap();
    assert!(r.converged, "{}", r.message);
    assert!(r.message.contains("Negative modes: 1"));
    let o = output(QcProgram::PySCF, JobType::TransitionState, &dir, &atoms, r).unwrap();
    assert!(o
        .saveable
        .iter()
        .any(|(n, p)| *n == "TS Hessian" && p.exists()));
    let r = compute(
        &mut Quadratic(1.),
        vec![0., 0., 0.],
        model(1),
        &atoms,
        &settings,
        JobType::TransitionState,
        &dir,
    )
    .unwrap();
    assert!(!r.converged);
    settings.ts.max_cycles = 1;
    let r = compute(
        &mut Quadratic(-1.),
        vec![0.5, 0.2, 0.1],
        model(1),
        &atoms,
        &settings,
        JobType::TransitionState,
        &dir,
    )
    .unwrap();
    assert!(!r.converged);
}
#[test]
fn geometry_search_cartesian_minimum_algorithms() {
    for algorithm in [
        MinimumAlgorithm::Bfgs,
        MinimumAlgorithm::SteepestDescent,
        MinimumAlgorithm::Bb1,
        MinimumAlgorithm::Bb2,
    ] {
        let options = search::Minimum {
            coordinates: Coordinates::Cartesian,
            algorithm,
            ..Default::default()
        };
        let r =
            minimum::compute(&mut Quadratic(1.), vec![0.2, 0.1, -0.1], model(1), &options).unwrap();
        assert!(r.converged, "{algorithm:?}: {}", r.message);
        assert!(r.energy.abs() < 1e-6);
    }
}
#[test]
fn geometry_search_dummy_references_remap_and_reject_dummy_selection() {
    let atoms = vec!["X".into(), "H".into(), "H".into()];
    let mut o = search::Options::default();
    o.reaction = search::Reaction::Bond;
    o.reaction_atoms = "2,3".into();
    o.extra_bonds = "2,3".into();
    let mapped = search::remap_dummies(&o, &atoms).unwrap();
    assert_eq!(mapped.reaction_atoms, "1,2");
    assert_eq!(mapped.extra_bonds, "1,2");
    o.reaction_atoms = "1,2".into();
    assert!(search::remap_dummies(&o, &atoms).is_err());
}
#[test]
fn geometry_search_native_behemoth_command_and_convergence() {
    let mut config = MethodConfig::default();
    config.search.ts.max_cycles = 7;
    config.search.ts.hessian_recalc_interval = 3;
    let c = native::command(
        Path::new("/user chosen/behemoth"),
        Path::new("/tmp"),
        0,
        1,
        &config,
    )
    .unwrap();
    let args: Vec<_> = c
        .get_args()
        .map(|a| a.to_string_lossy().into_owned())
        .collect();
    assert!(args.windows(2).any(|a| a == ["--ts-max-cycles", "7"]));
    assert!(args
        .windows(2)
        .any(|a| a == ["--ts-hessian-recalc-interval", "3"]));
    let log="Transition-state result\n converged : true\n cycles : 8\n E(final) : -5.0 Eh\n negative_modes : 1\n";
    let xyz = "1\nfinal\nHe 0 0 0\n";
    let history = "1\nstep= 1 Ene = -4.0 GradMax = 0.1\nX 0.1 0 0\n";
    let atoms = vec!["He".into()];
    assert!(
        native::read_result(log, xyz, history, &atoms)
            .unwrap()
            .converged
    );
    assert!(
        !native::read_result(
            &log.replace("negative_modes : 1", "negative_modes : 0"),
            xyz,
            history,
            &atoms
        )
        .unwrap()
        .converged
    );
    assert!(
        !native::read_result(&log.replace("negative_modes : 1", ""), xyz, history, &atoms)
            .unwrap()
            .converged
    );
    assert!(native::read_result(&log.replace("-5.0", "NaN"), xyz, history, &atoms).is_err());
}
#[test]
fn geometry_search_orca_hessian_interval_and_cancellation() {
    let mut config = MethodConfig::default();
    config.search.orca_recalc_hessian = 9;
    let input = |c: &MethodConfig| {
        orca_run::panel_input_text(
            JobType::TransitionState,
            QcMethod::Hf,
            0,
            1,
            c,
            Default::default(),
            "",
            "H 0 0 0\nH 0 0 0.74\n",
        )
    };
    assert!(input(&config).contains("Recalc_Hess 9"));
    config.search.orca_recalc_hessian = 0;
    assert!(!input(&config).contains("Recalc_Hess"));
    assert!(run(
        QcProgram::Xtb,
        Some(QcMethod::Gfn2),
        JobType::TransitionState,
        Path::new(""),
        Path::new("/tmp"),
        &[],
        &[],
        0,
        1,
        &config,
        &AtomicBool::new(true),
        &Mutex::new(None)
    )
    .err()
    .unwrap()
    .contains("Cancelled"));
}
#[test]
#[ignore = "requires installed backends and isolated XDG_CONFIG_HOME"]
fn geometry_search_installed_minima() {
    assert!(Path::new(&std::env::var("XDG_CONFIG_HOME").unwrap()).starts_with(std::env::temp_dir()));
    let atoms = vec!["O".into(), "H".into(), "H".into()];
    let xyz = [0., 0., 0., 0.8, 0., 0.6, -0.8, 0., 0.6];
    for (program, env, method) in [
        (QcProgram::PySCF, "BEAVYR_TEST_PYSCF", QcMethod::Hf),
        (QcProgram::Orca, "BEAVYR_TEST_ORCA", QcMethod::Hf),
        (QcProgram::Xtb, "BEAVYR_TEST_XTB", QcMethod::Gfn2),
    ] {
        let binary = std::env::var(env).unwrap();
        if program == QcProgram::PySCF {
            crate::qchem_interfaces::python_env::launch::save(
                "pyscf",
                crate::qchem_interfaces::python_env::launch::Settings {
                    interpreter: binary.clone(),
                    ..Default::default()
                },
            )
            .unwrap();
        }
        for coords in [
            Coordinates::Primitive,
            Coordinates::Pulay,
            Coordinates::Cartesian,
        ] {
            let mut config = MethodConfig::default();
            config.basis = "sto-3g".into();
            config.nproc = 1;
            config.search.driver = search::Driver::Beavyr;
            config.search.minimum.coordinates = coords;
            let dir = scratch(&format!("{program:?}-{coords:?}"));
            let out = run(
                program,
                Some(method),
                if program == QcProgram::PySCF && coords == Coordinates::Primitive {
                    JobType::OptimizeThenFrequencies
                } else {
                    JobType::Optimize
                },
                Path::new(&binary),
                &dir,
                &atoms,
                &xyz,
                0,
                1,
                &config,
                &AtomicBool::new(false),
                &Mutex::new(None),
            )
            .unwrap();
            assert_eq!(out.converged, Some(true), "{}: {}", dir.display(), out.log);
            assert!(out.trajectory_text.is_some());
            println!("{}: {}", dir.display(), out.log);
            if program == QcProgram::PySCF && coords == Coordinates::Primitive {
                assert!(out.hessian_path.as_ref().is_some_and(|p| p.exists()));
                let out = write_orbitals(
                    out,
                    Some(method),
                    Path::new(&binary),
                    &dir,
                    0,
                    1,
                    &config,
                    &AtomicBool::new(false),
                    &Mutex::new(None),
                )
                .unwrap();
                assert!(out.molden_path.as_ref().is_some_and(|p| p.exists()));
            }
        }
    }
}
#[test]
#[ignore = "requires installed backends and isolated XDG_CONFIG_HOME"]
fn geometry_search_installed_ammonia_ts() {
    assert!(Path::new(&std::env::var("XDG_CONFIG_HOME").unwrap()).starts_with(std::env::temp_dir()));
    let atoms = vec!["N".into(), "H".into(), "H".into(), "H".into()];
    let xyz = [
        0., 0., 0., 1., 0., 0., -0.5, 0.8660254, 0., -0.5, -0.8660254, 0.,
    ];
    for (program, env, method) in [
        (QcProgram::PySCF, "BEAVYR_TEST_PYSCF", QcMethod::Hf),
        (QcProgram::Xtb, "BEAVYR_TEST_XTB", QcMethod::Gfn1),
        (QcProgram::Xtb, "BEAVYR_TEST_XTB", QcMethod::Gfn2),
        (QcProgram::Behemoth, "BEAVYR_TEST_BEHEMOTH", QcMethod::Gfn1),
    ] {
        let binary = std::env::var(env).unwrap();
        if program == QcProgram::PySCF {
            crate::qchem_interfaces::python_env::launch::save(
                "pyscf",
                crate::qchem_interfaces::python_env::launch::Settings {
                    interpreter: binary.clone(),
                    ..Default::default()
                },
            )
            .unwrap();
        }
        let mut config = MethodConfig::default();
        config.basis = "sto-3g".into();
        config.nproc = 1;
        // Cartesian coordinates retain the inversion coordinate at the planar guess.
        config.search.ts.coordinate_mode =
            crate::optimizer::transition_state::TransitionStateCoordinateMode::Cartesian;
        config.search.ts.max_cycles = 30;
        let dir = scratch(&format!("{program:?}-{method:?}-ammonia-ts"));
        let out = run(
            program,
            Some(method),
            JobType::TransitionState,
            Path::new(&binary),
            &dir,
            &atoms,
            &xyz,
            0,
            1,
            &config,
            &AtomicBool::new(false),
            &Mutex::new(None),
        )
        .unwrap();
        if program == QcProgram::Xtb && method == QcMethod::Gfn2 {
            // The installed GFN2 build gives zero negative projected modes
            // here; a stationary structure alone must not pass TS validation.
            assert_eq!(out.converged, Some(false));
            assert!(out.log.contains("Negative modes: 0"));
        } else {
            assert_eq!(out.converged, Some(true), "{}: {}", dir.display(), out.log);
        }
        println!("{}: {}", dir.display(), out.log);
    }
}
struct BondMaximum;
impl Objective for BondMaximum {
    fn energy(&mut self, x: &[f64]) -> anyhow::Result<f64> {
        let r = (0..3)
            .map(|i| (x[i + 3] - x[i]).powi(2))
            .sum::<f64>()
            .sqrt();
        Ok(-0.5 * (r - 1.5).powi(2))
    }
    fn gradient(&mut self, x: &[f64], g: &mut [f64]) -> anyhow::Result<()> {
        let r = (0..3)
            .map(|i| (x[i + 3] - x[i]).powi(2))
            .sum::<f64>()
            .sqrt();
        for i in 0..3 {
            g[i] = (r - 1.5) * (x[i + 3] - x[i]) / r;
            g[i + 3] = -g[i];
        }
        Ok(())
    }
}
#[test]
fn geometry_search_ts_all_coordinate_systems() {
    use crate::optimizer::transition_state::*;
    for coords in [
        TransitionStateCoordinateMode::Cartesian,
        TransitionStateCoordinateMode::Internal,
        TransitionStateCoordinateMode::PulayRedundant,
    ] {
        let mut settings = search::Options::default();
        settings.ts.coordinate_mode = coords;
        settings.reaction = search::Reaction::Bond;
        settings.reaction_atoms = "1,2".into();
        let dir = scratch(&format!("bond-ts-{coords:?}"));
        let r = compute(
            &mut BondMaximum,
            vec![0., 0., 0., 1.6, 0., 0.],
            model(2),
            &["H".into(), "H".into()],
            &settings,
            JobType::TransitionState,
            &dir,
        )
        .unwrap();
        assert!(r.converged, "{coords:?}: {}", r.message);
        assert!(r.energy.abs() < 1e-7);
    }
}
#[test]
fn geometry_search_ui_exposes_native_and_local_controls() {
    use bevy_egui::egui;
    let ctx = egui::Context::default();
    let task = QcRunTask::default();
    for (program, job, driver) in [
        (
            QcProgram::Orca,
            JobType::TransitionState,
            search::Driver::Native,
        ),
        (
            QcProgram::Orca,
            JobType::TransitionState,
            search::Driver::Beavyr,
        ),
        (QcProgram::PySCF, JobType::Optimize, search::Driver::Native),
        (QcProgram::PySCF, JobType::Optimize, search::Driver::Beavyr),
        (
            QcProgram::Behemoth,
            JobType::TransitionState,
            search::Driver::Native,
        ),
    ] {
        let mut panel = crate::qchem_panel::QcPanelState::default();
        panel.program = program;
        panel.job = job;
        panel.config.search.driver = driver;
        let mut labels = vec![];
        for _ in 0..2 {
            let mut frame = ctx.run_ui(Default::default(), |ui| {
                crate::qchem_panel::geometry_search_controls::settings(ui, &mut panel, &task)
            });
            frame.textures_delta.clear();
            labels = frame
                .shapes
                .into_iter()
                .filter_map(|s| {
                    if let egui::Shape::Text(t) = s.shape {
                        Some(t.galley.text().to_owned())
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
        }
        if job == JobType::TransitionState {
            assert!(labels
                .iter()
                .any(|s| s == "Recalculate Hessian every N steps"));
        }
        if program == QcProgram::PySCF {
            assert!(labels.iter().any(|s| s == "PySCF / geomeTRIC"));
            assert!(labels.iter().any(|s| s == "Beavyr"));
        }
        if program == QcProgram::Behemoth || driver == search::Driver::Beavyr {
            assert!(labels.iter().any(|s| s == "Primitive internal"));
        }
    }
}
