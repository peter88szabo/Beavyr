use super::*;
use crate::qchem_interfaces::job::scan::{Coordinate, FixedCoordinate, Options};

fn scratch(label: &str) -> PathBuf {
    let id = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("beavyr-scan-{label}-{id}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}
fn water() -> (Vec<String>, Vec<f64>) {
    (
        vec!["O".into(), "H".into(), "H".into()],
        vec![0., 0., 0., 0.96, 0., 0., -0.24, 0.93, 0.],
    )
}
fn report(atoms: Vec<String>, options: Options) -> ScanReport {
    ScanReport {
        program: "test".into(),
        method: "test".into(),
        driver: "Beavyr".into(),
        options,
        atoms,
        points: vec![],
        complete: false,
        message: String::new(),
    }
}
#[test]
fn scan_validation_remaps_dummy_atoms_and_rejects_invalid_constraints() {
    let atoms = vec!["X".into(), "O".into(), "H".into(), "X".into(), "H".into()];
    let mut o = Options {
        atoms: "2 3".into(),
        constraints: vec![FixedCoordinate {
            coordinate: Coordinate::Angle,
            atoms: "3 2 5".into(),
            value: 104.5,
        }],
        ..Default::default()
    };
    let clean = o.remap_dummies(&atoms).unwrap();
    assert_eq!(clean.atoms, "1 2");
    assert_eq!(clean.constraints[0].atoms, "2 1 3");
    o.atoms = "1 3".into();
    assert!(o.prepare(&atoms).unwrap_err().contains("dummy"));
    o.atoms = "2 3".into();
    o.constraints[0] = FixedCoordinate {
        coordinate: Coordinate::Bond,
        atoms: "3 2".into(),
        value: 1.,
    };
    assert!(o
        .prepare(&atoms)
        .unwrap_err()
        .contains("both scanned and fixed"));
    o.constraints.clear();
    o.points = 1;
    assert!(o.prepare(&atoms).is_err());
    o.points = 3;
    o.end = f64::NAN;
    assert!(o.prepare(&atoms).is_err());
}
#[test]
fn scan_ranges_include_both_ends_and_old_settings_load() {
    let o = Options {
        start: 2.,
        end: 1.,
        points: 3,
        ..Default::default()
    };
    assert_eq!(o.values(), vec![2., 1.5, 1.]);
    let mut config = MethodConfig::default();
    config.scan = o.clone();
    let saved = ron::to_string(&config).unwrap();
    assert_eq!(ron::from_str::<MethodConfig>(&saved).unwrap(), config);
    let field = format!("scan:{},", ron::to_string(&o).unwrap());
    assert!(saved.contains(&field));
    let old = saved.replace(&field, "");
    assert_eq!(
        ron::from_str::<MethodConfig>(&old).unwrap().scan,
        Options::default()
    );
}
#[test]
fn scan_rigid_bonds_angles_and_wrapped_dihedrals_reach_target() {
    let (atoms, x) = water();
    for (c, selection, target) in [
        (Coordinate::Bond, "1 2", 1.3),
        (Coordinate::Angle, "2 1 3", 120.),
    ] {
        let o = Options {
            mode: Mode::Rigid,
            coordinate: c,
            atoms: selection.into(),
            move_fragment: false,
            ..Default::default()
        };
        let p = o.prepare(&atoms).unwrap();
        let moved = geometry::rigid(&atoms, &x, &o, &p, target).unwrap();
        assert!(
            geometry::residual(c, geometry::measure(c, &p.atoms, &moved).unwrap(), target).abs()
                < 1e-6
        );
        for i in 0..atoms.len() {
            if i != *p.atoms.last().unwrap() {
                assert_eq!(moved[3 * i..3 * i + 3], x[3 * i..3 * i + 3]);
            }
        }
    }
    let atoms = vec!["C".into(); 4];
    let x = vec![0., 1., 0., 0., 0., 0., 1.5, 0., 0., 1.5, 1., 1.];
    let o = Options {
        coordinate: Coordinate::Dihedral,
        atoms: "1 2 3 4".into(),
        start: -170.,
        end: 170.,
        move_fragment: false,
        ..Default::default()
    };
    let p = o.prepare(&atoms).unwrap();
    for target in [-180., -170., -45., 0., 45., 170., 180., 270.] {
        let moved = geometry::rigid(&atoms, &x, &o, &p, target).unwrap();
        assert!(
            geometry::residual(
                o.coordinate,
                geometry::measure(o.coordinate, &p.atoms, &moved).unwrap(),
                target
            )
            .abs()
                < 1e-6
        );
        let ic = geometry::internals(&atoms, &moved, &o, &p).unwrap();
        let bohr: Vec<_> = moved.iter().map(|x| x / BOHR_TO_ANGSTROM).collect();
        let q = crate::optimizer::internal_coords::compute_internals(&bohr, &ic);
        let index = ic
            .dihedrals
            .iter()
            .position(|a| *a == [0, 1, 2, 3] || *a == [3, 2, 1, 0])
            .unwrap();
        assert!(
            geometry::residual(
                o.coordinate,
                q[ic.bonds.len() + ic.angles.len() + index].to_degrees(),
                target
            )
            .abs()
                < 1e-6
        );
    }
}
#[test]
fn scan_fragment_moves_attached_atoms_and_rejects_rings() {
    let atoms = vec!["C".into(), "C".into(), "H".into()];
    let x = vec![0., 0., 0., 1.5, 0., 0., 2.5, 0., 0.];
    let o = Options {
        mode: Mode::Rigid,
        ..Default::default()
    };
    let p = o.prepare(&atoms).unwrap();
    let moved = geometry::rigid(&atoms, &x, &o, &p, 2.).unwrap();
    assert!((moved[6] - moved[3] - 1.).abs() < 1e-8);
    assert!((moved[6] - 3.).abs() < 1e-8);
    let atoms = vec!["C".into(); 3];
    let x = vec![0., 0., 0., 1.5, 0., 0., 0.75, 1.3, 0.];
    assert!(geometry::rigid(&atoms, &x, &o, &p, 2.)
        .unwrap_err()
        .contains("ring"));
}
struct PairPotential;
impl Objective for PairPotential {
    fn energy(&mut self, x: &[f64]) -> anyhow::Result<f64> {
        let mut g = vec![0.; x.len()];
        self.energy_gradient(x, &mut g)
    }
    fn gradient(&mut self, x: &[f64], g: &mut [f64]) -> anyhow::Result<()> {
        self.energy_gradient(x, g)?;
        Ok(())
    }
    fn energy_gradient(&mut self, x: &[f64], g: &mut [f64]) -> anyhow::Result<f64> {
        g.fill(0.);
        let mut e = 0.;
        for (a, b, r0) in [(0, 1, 1.8), (0, 2, 1.8), (1, 2, 2.9)] {
            let d: Vec<_> = (0..3).map(|k| x[3 * b + k] - x[3 * a + k]).collect();
            let r = d.iter().map(|v| v * v).sum::<f64>().sqrt();
            e += 0.5 * (r - r0).powi(2);
            for k in 0..3 {
                let v = (r - r0) * d[k] / r;
                g[3 * b + k] += v;
                g[3 * a + k] -= v;
            }
        }
        Ok(e)
    }
}
#[test]
fn scan_relaxed_preserves_fixed_coordinates_and_saves_matching_energies() {
    let (atoms, x) = water();
    let o = Options {
        start: 0.95,
        end: 1.05,
        points: 3,
        constraints: vec![FixedCoordinate {
            coordinate: Coordinate::Angle,
            atoms: "2 1 3".into(),
            value: 104.5,
        }],
        ..Default::default()
    };
    for mode in [
        scan::RelaxationCoordinates::Primitive,
        scan::RelaxationCoordinates::Pulay,
    ] {
        let mut o = o.clone();
        o.relaxation_coordinates = mode;
        let mut r = report(atoms.clone(), o);
        let dir = scratch("constraints");
        compute(
            &mut PairPotential,
            &x,
            &mut r,
            &dir,
            &AtomicBool::new(false),
        )
        .unwrap();
        assert_eq!(r.points.len(), 3);
        for p in &r.points {
            assert_eq!(p.converged, Some(true), "{p:?}");
            assert!((p.actual - p.target).abs() < 1e-3);
            assert!(
                (geometry::measure(Coordinate::Angle, &[1, 0, 2], &p.coordinates).unwrap() - 104.5)
                    .abs()
                    < 1e-3
            );
            let x: Vec<_> = p.coordinates.iter().map(|v| v / BOHR_TO_ANGSTROM).collect();
            assert!((PairPotential.energy(&x).unwrap() - p.energy).abs() < 1e-10);
        }
        let restored: ScanReport =
            ron::from_str(&fs::read_to_string(dir.join("scan_results.ron")).unwrap()).unwrap();
        assert_eq!(restored.points.len(), 3);
        assert_eq!(
            crate::trajectory::parse_multi_xyz(&r.xyz()).unwrap().len(),
            3
        );
        assert_eq!(r.csv().lines().count(), 4);
    }
}
#[test]
fn scan_nonconvergence_is_reported_and_stop_option_is_honoured() {
    let (atoms, x) = water();
    for stop in [true, false] {
        let o = Options {
            start: 1.2,
            end: 1.3,
            points: 3,
            max_cycles: 1,
            stop_on_failure: stop,
            ..Default::default()
        };
        let mut r = report(atoms.clone(), o);
        compute(
            &mut PairPotential,
            &x,
            &mut r,
            &scratch("nonconvergence"),
            &AtomicBool::new(false),
        )
        .unwrap();
        assert_eq!(r.points.len(), if stop { 1 } else { 3 });
        assert!(r.points.iter().any(|p| p.converged == Some(false)));
    }
}
#[test]
fn scan_cancellation_keeps_completed_points_and_does_not_add_cancelled_point() {
    struct Cancel<'a> {
        flag: &'a AtomicBool,
        calls: usize,
    }
    impl Objective for Cancel<'_> {
        fn energy(&mut self, _: &[f64]) -> anyhow::Result<f64> {
            self.calls += 1;
            if self.calls == 2 {
                self.flag.store(true, Ordering::SeqCst);
            }
            Ok(-1.)
        }
        fn gradient(&mut self, _: &[f64], _: &mut [f64]) -> anyhow::Result<()> {
            unreachable!()
        }
    }
    let flag = AtomicBool::new(false);
    let mut objective = Cancel {
        flag: &flag,
        calls: 0,
    };
    let (atoms, x) = water();
    let mut r = report(
        atoms,
        Options {
            mode: Mode::Rigid,
            points: 3,
            start: 0.9,
            end: 1.1,
            ..Default::default()
        },
    );
    let dir = scratch("cancel");
    assert_eq!(
        compute(&mut objective, &x, &mut r, &dir, &flag).unwrap_err(),
        "Cancelled."
    );
    assert_eq!(r.points.len(), 1);
    assert_eq!(
        crate::trajectory::parse_multi_xyz(&fs::read_to_string(dir.join("scan.xyz")).unwrap())
            .unwrap()
            .len(),
        1
    );
}
#[test]
fn scan_limits_correlated_relaxation_and_spectrum_only_backends() {
    for p in [QcProgram::Orca, QcProgram::PySCF, QcProgram::Behemoth] {
        assert!(scan::method_problem(p, Some(QcMethod::Casscf), Mode::Rigid).is_none());
        assert!(scan::method_problem(p, Some(QcMethod::Casscf), Mode::Relaxed).is_some());
    }
    assert!(scan::method_problem(QcProgram::StdaXtb, None, Mode::Rigid).is_some());
}
#[test]
fn scan_pyscf_mp2_uses_correlated_energy_and_gradient() {
    let (atoms, x) = water();
    let text = pyscf_script(
        JobType::Gradient,
        QcMethod::Mp2,
        0,
        1,
        &MethodConfig::default(),
        &atoms,
        &x,
    );
    assert!(text.contains("correlation = mp.MP2(mf)"));
    assert!(text.contains("g = correlation.nuc_grad_method().kernel()"));
    assert!(!text.contains("mf = dft.RKS"));
}

#[test]
#[ignore = "requires BEAVYR_TEST_ORCA and BEAVYR_TEST_BEHEMOTH executables"]
fn scan_installed_orca_and_behemoth() {
    let (atoms, x) = water();
    for (program, env, method, drivers) in [
        (
            QcProgram::Orca,
            "BEAVYR_TEST_ORCA",
            QcMethod::Hf,
            vec![Driver::Orca, Driver::Beavyr],
        ),
        (
            QcProgram::Behemoth,
            "BEAVYR_TEST_BEHEMOTH",
            QcMethod::Hf,
            vec![Driver::Beavyr],
        ),
        (
            QcProgram::Xtb,
            "BEAVYR_TEST_XTB",
            QcMethod::Gfn2,
            vec![Driver::Beavyr],
        ),
    ] {
        let path = std::env::var(env).unwrap();
        for driver in drivers {
            for mode in [Mode::Rigid, Mode::Relaxed] {
                let dir = scratch(&format!("native-{program:?}-{driver:?}-{mode:?}"));
                let mut config = MethodConfig::default();
                config.basis = "sto-3g".into();
                config.nproc = 1;
                config.scan = Options {
                    driver,
                    mode,
                    start: 0.95,
                    end: 1.0,
                    points: 3,
                    ..Default::default()
                };
                if mode == Mode::Relaxed {
                    config.scan.constraints.push(FixedCoordinate {
                        coordinate: Coordinate::Angle,
                        atoms: "2 1 3".into(),
                        value: 104.5,
                    });
                }
                let output = run(
                    program,
                    Some(method),
                    Path::new(&path),
                    &dir,
                    &atoms,
                    &x,
                    0,
                    1,
                    &config,
                    &AtomicBool::new(false),
                    &Mutex::new(None),
                )
                .unwrap();
                let result = output.scan.unwrap();
                assert!(result.complete, "{}: {}", dir.display(), result.message);
                assert_eq!(result.points.len(), 3);
                for p in &result.points {
                    assert!((p.target - p.actual).abs() < 1e-3);
                    assert!(p.energy.is_finite());
                }
                println!(
                    "{program:?} {driver:?} {mode:?}: {:?}; {}",
                    result.points.iter().map(|p| p.energy).collect::<Vec<_>>(),
                    dir.display()
                );
            }
        }
    }
}
#[test]
#[ignore = "requires isolated XDG_CONFIG_HOME and BEAVYR_TEST_PYSCF / BEAVYR_TEST_PSI4 interpreters"]
fn scan_installed_python_backends() {
    // Never change the user's saved environment from this opt-in integration test.
    let config_root =
        std::env::var("XDG_CONFIG_HOME").expect("Use an isolated test config directory");
    assert!(Path::new(&config_root).starts_with(std::env::temp_dir()));
    let (atoms, x) = water();
    for (program, module, env) in [
        (QcProgram::PySCF, "pyscf", "BEAVYR_TEST_PYSCF"),
        (QcProgram::Psi4Py, "psi4", "BEAVYR_TEST_PSI4"),
    ] {
        crate::qchem_interfaces::python_env::launch::save(
            module,
            crate::qchem_interfaces::python_env::launch::Settings {
                interpreter: std::env::var(env).unwrap(),
                ..Default::default()
            },
        )
        .unwrap();
        for mode in [Mode::Rigid, Mode::Relaxed] {
            let mut config = MethodConfig::default();
            config.basis = "sto-3g".into();
            config.nproc = 1;
            config.scan = Options {
                mode,
                start: 0.95,
                end: 1.0,
                points: 3,
                ..Default::default()
            };
            let dir = scratch(&format!("native-{program:?}-{mode:?}"));
            let output = run(
                program,
                Some(QcMethod::Hf),
                Path::new(""),
                &dir,
                &atoms,
                &x,
                0,
                1,
                &config,
                &AtomicBool::new(false),
                &Mutex::new(None),
            )
            .unwrap();
            let r = output.scan.unwrap();
            assert!(r.complete, "{}: {}", dir.display(), r.message);
            assert_eq!(r.points.len(), 3);
            println!(
                "{program:?} {mode:?}: {:?}; {}",
                r.points.iter().map(|p| p.energy).collect::<Vec<_>>(),
                dir.display()
            );
        }
    }
}

#[cfg(unix)]
#[test]
fn scan_native_uses_exact_executable_and_never_reports_bad_exit_as_success() {
    use std::os::unix::fs::PermissionsExt;
    let dir = scratch("selected executable");
    let binary = dir.join("chosen orca");
    fs::write(&binary,"#!/bin/sh\n[ \"$#\" -eq 1 ] || exit 21\n[ \"$1\" = orca_scan.inp ] || exit 22\nprintf '1 -1.1\\n2 -1.2\\n' > orca_scan.xyzact.dat\nprintf 'FINAL SINGLE POINT ENERGY -1.1\\nFINAL SINGLE POINT ENERGY -1.2\\nORCA TERMINATED NORMALLY\\n'\nexit 9\n").unwrap();
    fs::set_permissions(&binary, fs::Permissions::from_mode(0o700)).unwrap();
    let (atoms, x) = water();
    let mut config = MethodConfig::default();
    config.scan = Options {
        driver: Driver::Orca,
        mode: Mode::Rigid,
        points: 2,
        start: 0.9,
        end: 1.,
        ..Default::default()
    };
    let output = run(
        QcProgram::Orca,
        Some(QcMethod::Hf),
        &binary,
        &dir,
        &atoms,
        &x,
        0,
        1,
        &config,
        &AtomicBool::new(false),
        &Mutex::new(None),
    )
    .unwrap();
    let result = output.scan.unwrap();
    assert!(!result.complete);
    assert_eq!(result.points.len(), 2);
    assert!(result.message.contains("exited with"));
    assert!(!fs::read_to_string(dir.join("orca_scan.inp"))
        .unwrap()
        .contains("mpirun"));
}
#[test]
fn scan_dreiding_removes_dummies_and_handles_both_modes() {
    let (mut atoms, mut x) = water();
    atoms.insert(0, "X".into());
    x.splice(0..0, [9., 9., 9.]);
    for mode in [Mode::Rigid, Mode::Relaxed] {
        let mut config = MethodConfig::default();
        config.scan = Options {
            mode,
            atoms: "2 3".into(),
            points: 3,
            start: 0.95,
            end: 1.,
            ..Default::default()
        };
        let output = run(
            QcProgram::Dreiding,
            None,
            Path::new(""),
            &scratch("dreiding"),
            &atoms,
            &x,
            0,
            1,
            &config,
            &AtomicBool::new(false),
            &Mutex::new(None),
        )
        .unwrap();
        let r = output.scan.unwrap();
        assert!(r.complete, "{}", r.message);
        assert_eq!(r.atoms, vec!["O", "H", "H"]);
        assert_eq!(r.options.atoms, "1 2");
        assert_eq!(r.points.len(), 3);
    }
}

#[test]
#[ignore = "requires BEAVYR_TEST_ORCA executable"]
fn scan_installed_orca_angles_and_dihedrals() {
    let path = std::env::var("BEAVYR_TEST_ORCA").unwrap();
    let (water_atoms, water_x) = water();
    for (coordinate, selected, atoms, x, start, end) in [
        (Coordinate::Angle, "2 1 3", water_atoms, water_x, 100., 110.),
        (
            Coordinate::Dihedral,
            "1 2 3 4",
            vec!["H".into(), "O".into(), "O".into(), "H".into()],
            vec![-0.3, 0.9, 0., 0., 0., 0., 1.45, 0., 0., 1.75, 0.2, 0.9],
            90.,
            110.,
        ),
    ] {
        for driver in [Driver::Orca, Driver::Beavyr] {
            let mut config = MethodConfig::default();
            config.nproc = 1;
            config.basis = "sto-3g".into();
            config.scan = Options {
                driver,
                coordinate,
                atoms: selected.into(),
                start,
                end,
                points: 3,
                ..Default::default()
            };
            let dir = scratch(&format!("angular-{driver:?}-{coordinate:?}"));
            let r = run(
                QcProgram::Orca,
                Some(QcMethod::Hf),
                Path::new(&path),
                &dir,
                &atoms,
                &x,
                0,
                1,
                &config,
                &AtomicBool::new(false),
                &Mutex::new(None),
            )
            .unwrap()
            .scan
            .unwrap();
            assert!(r.complete, "{}: {}", dir.display(), r.message);
            for p in &r.points {
                assert!(geometry::residual(coordinate, p.actual, p.target).abs() < 1e-3);
            }
            println!(
                "ORCA {driver:?} {coordinate:?}: {:?}; {}",
                r.points.iter().map(|p| p.energy).collect::<Vec<_>>(),
                dir.display()
            );
        }
    }
}
#[test]
#[ignore = "requires isolated XDG_CONFIG_HOME and BEAVYR_TEST_BEHEMOTH / BEAVYR_TEST_PYSCF"]
fn scan_installed_mp2_energy_and_gradients() {
    let root = std::env::var("XDG_CONFIG_HOME").unwrap();
    assert!(Path::new(&root).starts_with(std::env::temp_dir()));
    crate::qchem_interfaces::python_env::launch::save(
        "pyscf",
        crate::qchem_interfaces::python_env::launch::Settings {
            interpreter: std::env::var("BEAVYR_TEST_PYSCF").unwrap(),
            ..Default::default()
        },
    )
    .unwrap();
    let (atoms, x) = water();
    for (program, method) in [
        (QcProgram::Behemoth, QcMethod::Mp2),
        (QcProgram::Behemoth, QcMethod::RiMp2),
        (QcProgram::PySCF, QcMethod::Mp2),
    ] {
        let binary = if program == QcProgram::Behemoth {
            std::env::var("BEAVYR_TEST_BEHEMOTH").unwrap()
        } else {
            String::new()
        };
        let mut config = MethodConfig::default();
        config.nproc = 1;
        config.basis = "def2-svp".into();
        config.scan = Options {
            start: 0.96,
            end: 0.98,
            points: 2,
            ..Default::default()
        };
        let dir = scratch(&format!("mp2-{program:?}-{method:?}"));
        let r = run(
            program,
            Some(method),
            Path::new(&binary),
            &dir,
            &atoms,
            &x,
            0,
            1,
            &config,
            &AtomicBool::new(false),
            &Mutex::new(None),
        )
        .unwrap()
        .scan
        .unwrap();
        assert!(r.complete, "{}: {}", dir.display(), r.message);
        assert!(r.points.iter().all(|p| p.energy < -75.9));
        println!(
            "{program:?} {method:?}: {:?}; {}",
            r.points.iter().map(|p| p.energy).collect::<Vec<_>>(),
            dir.display()
        );
    }
}

#[test]
#[ignore = "requires BEAVYR_TEST_BEHEMOTH executable"]
fn scan_installed_behemoth_builtin() {
    let binary = std::env::var("BEAVYR_TEST_BEHEMOTH").unwrap();
    let (atoms, x) = water();
    for (coordinate, selected, start, end) in [
        (Coordinate::Bond, "1 2", 0.95, 1.0),
        (Coordinate::Angle, "2 1 3", 100., 110.),
    ] {
        for mode in [Mode::Rigid, Mode::Relaxed] {
            let mut config = MethodConfig::default();
            config.nproc = 1;
            config.basis = "sto-3g".into();
            config.scan = Options {
                driver: Driver::Behemoth,
                coordinate,
                atoms: selected.into(),
                start,
                end,
                points: 3,
                mode,
                move_fragment: false,
                ..Default::default()
            };
            if mode == Mode::Relaxed {
                config
                    .scan
                    .constraints
                    .push(if coordinate == Coordinate::Bond {
                        FixedCoordinate {
                            coordinate: Coordinate::Angle,
                            atoms: "2 1 3".into(),
                            value: 104.5,
                        }
                    } else {
                        FixedCoordinate {
                            coordinate: Coordinate::Bond,
                            atoms: "1 2".into(),
                            value: 0.96,
                        }
                    });
            }
            let dir = scratch(&format!("builtin-behemoth-{coordinate:?}-{mode:?}"));
            let r = run(
                QcProgram::Behemoth,
                Some(QcMethod::Hf),
                Path::new(&binary),
                &dir,
                &atoms,
                &x,
                0,
                1,
                &config,
                &AtomicBool::new(false),
                &Mutex::new(None),
            )
            .unwrap()
            .scan
            .unwrap();
            assert!(r.complete, "{}: {}", dir.display(), r.message);
            assert_eq!(r.points.len(), 3);
            println!(
                "Behemoth built-in {coordinate:?} {mode:?}: {:?}; {}",
                r.points.iter().map(|p| p.energy).collect::<Vec<_>>(),
                dir.display()
            );
        }
    }
}

#[test]
#[ignore = "requires rebuilt BEAVYR_TEST_BEHEMOTH with extended scan CLI"]
fn scan_installed_behemoth_extended_controls() {
    let binary = std::env::var("BEAVYR_TEST_BEHEMOTH").unwrap();
    let (atoms, x) = water();
    for stop in [true, false] {
        let mut config = MethodConfig::default();
        config.nproc = 1;
        config.basis = "sto-3g".into();
        config.scan = Options {
            driver: Driver::Behemoth,
            start: 1.2,
            end: 1.3,
            points: 3,
            max_cycles: 1,
            stop_on_failure: stop,
            ..Default::default()
        };
        let dir = scratch("native-behemoth-step-limit");
        let r = run(
            QcProgram::Behemoth,
            Some(QcMethod::Hf),
            Path::new(&binary),
            &dir,
            &atoms,
            &x,
            0,
            1,
            &config,
            &AtomicBool::new(false),
            &Mutex::new(None),
        )
        .unwrap()
        .scan
        .unwrap();
        assert!(!r.complete);
        assert_eq!(r.points.len(), if stop { 1 } else { 3 });
        assert!(r.points.iter().any(|p| p.converged == Some(false)));
        println!(
            "Behemoth stop={stop}, {} saved points; {}",
            r.points.len(),
            dir.display()
        );
    }
    let mut config = MethodConfig::default();
    config.nproc = 1;
    config.basis = "sto-3g".into();
    config.scan = Options {
        driver: Driver::Behemoth,
        start: 0.7,
        end: 0.8,
        points: 2,
        ..Default::default()
    };
    let dir = scratch("native-behemoth-diatomic");
    let r = run(
        QcProgram::Behemoth,
        Some(QcMethod::Hf),
        Path::new(&binary),
        &dir,
        &["H".into(), "H".into()],
        &[0., 0., 0., 0.74, 0., 0.],
        0,
        1,
        &config,
        &AtomicBool::new(false),
        &Mutex::new(None),
    )
    .unwrap()
    .scan
    .unwrap();
    assert!(r.complete, "{}", r.message);
    assert!(r.points.iter().all(|p| p.steps == 0));
    config.scan = Options {
        driver: Driver::Behemoth,
        mode: Mode::Rigid,
        coordinate: Coordinate::Dihedral,
        atoms: "1 2 3 4".into(),
        start: 90.,
        end: 110.,
        points: 3,
        move_fragment: true,
        ..Default::default()
    };
    let dir = scratch("native-behemoth-rigid-dihedral");
    let r = run(
        QcProgram::Behemoth,
        Some(QcMethod::Hf),
        Path::new(&binary),
        &dir,
        &["H".into(), "O".into(), "O".into(), "H".into()],
        &[-0.3, 0.9, 0., 0., 0., 0., 1.45, 0., 0., 1.75, 0.2, 0.9],
        0,
        1,
        &config,
        &AtomicBool::new(false),
        &Mutex::new(None),
    )
    .unwrap()
    .scan
    .unwrap();
    assert!(r.complete, "{}", r.message);
    assert_eq!(r.points.len(), 3);
    println!("Behemoth rigid fragment dihedral passed; {}", dir.display());
}

#[test]
fn scan_bond_motion_preserves_fixed_side_or_midpoint() {
    use crate::qchem_interfaces::job::scan::BondMotion;
    let (atoms, x) = water();
    for motion in [BondMotion::Second, BondMotion::First, BondMotion::Symmetric] {
        for fragment in [false, true] {
            let o = Options {
                mode: Mode::Rigid,
                move_fragment: fragment,
                bond_motion: motion,
                ..Default::default()
            };
            let moved = geometry::rigid(&atoms, &x, &o, &o.prepare(&atoms).unwrap(), 1.16).unwrap();
            let expected = match motion {
                BondMotion::Second => (0., 1.16),
                BondMotion::First => (-0.2, 0.96),
                BondMotion::Symmetric => (-0.1, 1.06),
            };
            assert!((moved[0] - expected.0).abs() < 1e-8);
            assert!((moved[3] - expected.1).abs() < 1e-8);
            let expected_third = if fragment { x[6] + expected.0 } else { x[6] };
            assert!((moved[6] - expected_third).abs() < 1e-8);
        }
    }
}
