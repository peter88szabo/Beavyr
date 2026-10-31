use super::*;
use crate::optimizer::traits::Objective;
use crate::qchem_interfaces::job::mecp::{self, Coordinates};

struct Surface(f64);
impl Objective for Surface {
    fn energy(&mut self, x: &[f64]) -> anyhow::Result<f64> {
        Ok(0.5 * (x[0] - self.0).powi(2) + 0.5 * x[1..].iter().map(|v| v * v).sum::<f64>())
    }
    fn gradient(&mut self, x: &[f64], g: &mut [f64]) -> anyhow::Result<()> {
        g.copy_from_slice(x);
        g[0] -= self.0;
        Ok(())
    }
}
#[test]
fn mecp_two_surfaces_converge_to_crossing_and_export_both_energies() {
    let result = optimize_spin_crossing(
        vec![0.05, 0.02, 0.],
        None,
        &mut Surface(1.),
        &mut Surface(-1.),
        SpinCrossingOptions {
            coordinate_mode: crate::optimizer::spin_crossing::SpinCrossingCoordinateMode::Cartesian,
            trajectory_file: None,
            verbosity: 0,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(result.converged);
    assert!((result.energy_a - result.energy_b).abs() < 5e-5);
    let (xyz, csv) = exports(&result, &["He".into()]);
    assert!((result.energy_a - 0.5).abs() < 1e-4);
    assert!(xyz.contains("Ea=") && xyz.contains("Eb="));
    assert!(csv.contains("state_a_energy_Eh,state_b_energy_Eh,gap_Eh"));
    assert!(crate::trajectory::parse_multi_xyz(&xyz).unwrap().len() >= 2);
}
#[test]
fn mecp_step_limit_is_not_convergence() {
    let result = optimize_spin_crossing(
        vec![0.3, 0.1, 0.],
        None,
        &mut Surface(1.),
        &mut Surface(-1.),
        SpinCrossingOptions {
            coordinate_mode: crate::optimizer::spin_crossing::SpinCrossingCoordinateMode::Cartesian,
            max_cycles: 1,
            trajectory_file: None,
            verbosity: 0,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(!result.converged);
}
#[test]
fn mecp_validates_both_multiplicities_and_excludes_spin_independent_models() {
    let mol = crate::molecule::Molecule::from_xyz("2\nH2\nH 0 0 0\nH 0.74 0 0\n");
    let mut options = mecp::Options::default();
    assert!(mecp::validate(QcProgram::Orca, Some(QcMethod::Hf), &options, &mol, 0, 1).is_ok());
    options.multiplicity_b = 2;
    assert!(
        mecp::validate(QcProgram::Orca, Some(QcMethod::Hf), &options, &mol, 0, 1)
            .unwrap_err()
            .contains("State B")
    );
    options.multiplicity_b = 1;
    assert!(mecp::validate(QcProgram::Orca, Some(QcMethod::Hf), &options, &mol, 0, 1).is_err());
    options.multiplicity_b = 3;
    assert!(mecp::validate(QcProgram::Dreiding, None, &options, &mol, 0, 1).is_err());
    assert!(mecp::validate(
        QcProgram::PySCF,
        Some(QcMethod::CcsdT),
        &options,
        &mol,
        0,
        1
    )
    .is_err());
    options.coordinates = Coordinates::Cartesian;
    let mut config = MethodConfig::default();
    config.mecp = options;
    assert_eq!(
        ron::from_str::<MethodConfig>(&ron::to_string(&config).unwrap()).unwrap(),
        config
    );
}
#[test]
fn mecp_cancelled_run_does_not_launch_backend() {
    assert_eq!(
        run(
            QcProgram::Orca,
            Some(QcMethod::Hf),
            Path::new("does-not-exist"),
            Path::new("/tmp/unused-mecp"),
            &["H".into(), "H".into()],
            &[0., 0., 0., 0.74, 0., 0.],
            0,
            1,
            &MethodConfig::default(),
            &AtomicBool::new(true),
            &Mutex::new(None)
        )
        .err()
        .unwrap(),
        "Cancelled."
    );
}

#[test]
#[ignore = "Requires installed backend paths in BEAVYR_TEST_ORCA, BEAVYR_TEST_BEHEMOTH and Python environments"]
fn mecp_installed_backends_evaluate_distinct_spin_surfaces() {
    assert!(
        std::env::var("XDG_CONFIG_HOME")
            .unwrap_or_default()
            .starts_with("/tmp/"),
        "Use an isolated test config"
    );
    let atoms = vec!["C".into(), "H".into(), "H".into()];
    let x = vec![0., 0., 0., 1.05, 0., 0., -0.3, 1.0, 0.];
    for (program, variable) in [
        (QcProgram::Orca, "BEAVYR_TEST_ORCA"),
        (QcProgram::Behemoth, "BEAVYR_TEST_BEHEMOTH"),
        (QcProgram::PySCF, "BEAVYR_TEST_PYSCF"),
        (QcProgram::Psi4Py, "BEAVYR_TEST_PSI4"),
    ] {
        let binary = PathBuf::from(std::env::var(variable).unwrap());
        if let Some(module) = program.python_module() {
            crate::qchem_interfaces::python_env::launch::save(
                module,
                crate::qchem_interfaces::python_env::launch::Settings {
                    interpreter: binary.to_string_lossy().into_owned(),
                    ..Default::default()
                },
            )
            .unwrap();
            crate::qchem_interfaces::python_env::forget();
        }
        let mut config = MethodConfig::default();
        config.nproc = 1;
        config.basis = "sto-3g".into();
        config.mecp.max_cycles = 1;
        config.mecp.coordinates = Coordinates::Cartesian;
        let dir =
            std::env::temp_dir().join(format!("beavyr-mecp-{:?}-{}", program, std::process::id()));
        let result = run(
            program,
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
        assert_eq!(result.converged, Some(false));
        let csv = fs::read_to_string(dir.join("mecp.csv")).unwrap();
        let row: Vec<_> = csv.lines().nth(1).unwrap().split(',').collect();
        let gap = row[3].parse::<f64>().unwrap();
        assert!(gap.abs() > 1e-4, "spin surfaces must differ: {csv}");
        assert!(dir.join("state_a").is_dir() && dir.join("state_b").is_dir());
        println!(
            "{program:?} MECP distinct surfaces: gap={gap}, {}",
            dir.display()
        );
    }
}

#[test]
#[ignore = "Requires BEAVYR_TEST_ORCA executable"]
fn mecp_installed_orca_native_convergence_and_step_limit() {
    let binary = PathBuf::from(std::env::var("BEAVYR_TEST_ORCA").unwrap());
    let atoms = vec!["C".into(), "H".into(), "H".into()];
    let x = vec![0., 0., 0., 1.05, 0., 0., -0.3, 1., 0.];
    for (limit, expected) in [(1, false), (50, true)] {
        let mut config = MethodConfig::default();
        config.nproc = 1;
        config.basis = "sto-3g".into();
        config.mecp.driver = mecp::Driver::Orca;
        config.mecp.max_cycles = limit;
        let dir = std::env::temp_dir().join(format!(
            "beavyr-mecp-native-orca-{limit}-{}",
            std::process::id()
        ));
        let result = run(
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
        assert_eq!(result.converged, Some(expected), "{}", result.log);
        assert!(dir.join("mecp.csv").is_file());
        println!("ORCA native MECP limit={limit}: {}", result.log);
    }
}

#[test]
#[ignore = "Requires BEAVYR_TEST_ORCA executable"]
fn mecp_installed_beavyr_search_reaches_crossing() {
    let binary = PathBuf::from(std::env::var("BEAVYR_TEST_ORCA").unwrap());
    let mut config = MethodConfig::default();
    config.nproc = 1;
    config.basis = "sto-3g".into();
    config.mecp.max_cycles = 100;
    let dir = std::env::temp_dir().join(format!("beavyr-mecp-full-search-{}", std::process::id()));
    let result = run(
        QcProgram::Orca,
        Some(QcMethod::Hf),
        &binary,
        &dir,
        &["C".into(), "H".into(), "H".into()],
        &[0., 0., 0., 1.05, 0., 0., -0.3, 1., 0.],
        0,
        1,
        &config,
        &AtomicBool::new(false),
        &Mutex::new(None),
    )
    .unwrap();
    println!("Beavyr MECP with ORCA: {}\n{}", result.log, dir.display());
    assert_eq!(result.converged, Some(true), "{}", result.log);
}

#[test]
#[ignore = "Requires BEAVYR_TEST_ORCA and original SpinCross CH3O+ fixture"]
fn mecp_original_fortran_ch3o_reference_case() {
    let path =
        Path::new("/home/peter/Dropbox/Research_Leuven/SpinCross_code/Test_CH3Op/Initial_Geom.xyz");
    let text = fs::read_to_string(path).unwrap();
    let (_, atoms, coords) = crate::molecule::parse_xyz_angstrom(&text);
    let x: Vec<_> = coords.into_iter().flatten().collect();
    let mut config = MethodConfig::default();
    config.nproc = 1;
    config.memory_mb = 1000;
    config.functional = "B3LYP".into();
    config.basis = "6-31G".into();
    config.mecp.max_cycles = 40;
    let dir = std::env::temp_dir().join(format!("beavyr-mecp-fortran-ch3o-{}", std::process::id()));
    let binary = PathBuf::from(std::env::var("BEAVYR_TEST_ORCA").unwrap());
    let result = run(
        QcProgram::Orca,
        Some(QcMethod::Dft),
        &binary,
        &dir,
        &atoms,
        &x,
        1,
        1,
        &config,
        &AtomicBool::new(false),
        &Mutex::new(None),
    )
    .unwrap();
    println!(
        "Original SpinCross CH3O+ case: {}\n{}",
        result.log,
        dir.display()
    );
    assert_eq!(result.converged, Some(true), "{}", result.log);
    let csv = fs::read_to_string(dir.join("mecp.csv")).unwrap();
    let row: Vec<_> = csv.lines().last().unwrap().split(',').collect();
    assert!((row[1].parse::<f64>().unwrap() + 114.530635).abs() < 5e-5);
    assert!((row[2].parse::<f64>().unwrap() + 114.530632).abs() < 5e-5);
}
