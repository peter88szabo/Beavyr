use super::*;
use std::os::unix::fs::PermissionsExt;

#[test]
#[ignore = "requires the configured Python environment with PySCF and geometric"]
fn real_pyscf_preserves_actual_optimization_convergence() {
    use crate::qchem_interfaces::pyscf_run;

    let config = MethodConfig {
        functional: "hf".into(),
        basis: "sto-3g".into(),
        nproc: 2,
        ..Default::default()
    };
    let atoms = vec!["O".into(), "H".into(), "H".into()];
    let positions = [0.0, 0.0, 0.0, 0.0, 0.0, 0.96, 0.92, 0.0, -0.24];
    let script = pyscf_run::script_text(
        JobType::Optimize,
        0,
        1,
        &config,
        ExcitedStateOptions::default(),
        false,
        &atoms,
        &positions,
    );
    let directory = std::env::temp_dir().join(format!(
        "beavyr-real-pyscf-convergence-{}",
        std::process::id()
    ));
    for expected in [false, true] {
        let workdir = directory.join(if expected {
            "converged"
        } else {
            "iteration-limited"
        });
        fs::create_dir_all(&workdir).unwrap();
        let input = if expected {
            script.clone()
        } else {
            let input = script.replace(
                "kernel(mf, callback=record)",
                "kernel(mf, callback=record, maxsteps=1)",
            );
            assert_ne!(input, script, "the test must limit the real optimizer");
            input
        };
        let output = run_pyscf(
            JobType::Optimize,
            &workdir,
            &atoms,
            &positions,
            0,
            1,
            &config,
            ExcitedStateOptions::default(),
            false,
            Some(&input),
            &AtomicBool::new(false),
            &Mutex::new(None),
        )
        .unwrap_or_else(|error| panic!("{error}; inspect {}", workdir.display()));
        assert_eq!(output.converged, Some(expected));
        assert!(output.energetics.final_energy_hartree.unwrap().is_finite());
        assert!(!output.trajectory_text.unwrap().is_empty());
        let final_geometry = fs::read_to_string(pyscf_run::job_files(&workdir).geometry).unwrap();
        assert!(final_geometry.contains(if expected {
            "optimised geometry"
        } else {
            "unconverged geometry"
        }));
        println!(
            "Real PySCF geometry convergence: {expected}; energy {:.12} Hartree",
            output.energetics.final_energy_hartree.unwrap()
        );
    }
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn pyscf_optimization_preserves_false_true_and_missing_geometry_status() {
    const CHILD: &str = "BEAVYR_PYSCF_COMPLETION_TEST";
    let directory =
        std::env::temp_dir().join(format!("beavyr-pyscf-completion-{}", std::process::id()));
    if std::env::var_os(CHILD).is_none() {
        fs::create_dir_all(&directory).unwrap();
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args(["--exact", "qchem_panel::run::pyscf_completion_tests::pyscf_optimization_preserves_false_true_and_missing_geometry_status", "--nocapture", "--test-threads=1"])
            .env(CHILD, "1")
            .env("XDG_CONFIG_HOME", directory.join("config"))
            .output().unwrap();
        fs::remove_dir_all(&directory).unwrap();
        assert!(
            output.status.success(),
            "{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        return;
    }

    fs::create_dir_all(&directory).unwrap();
    let binary = directory.join("selected python");
    let config = crate::qchem_interfaces::config::config_path_for(QcProgram::Xtb)
        .unwrap()
        .with_file_name("python_pyscf.ron");
    fs::create_dir_all(config.parent().unwrap()).unwrap();
    fs::write(
        config,
        format!(
            "(interpreter: {:?}, activation: \"\", shell: \"/bin/bash\")",
            binary.to_str().unwrap()
        ),
    )
    .unwrap();
    let run = |job, workdir: &Path| {
        fs::create_dir_all(workdir).unwrap();
        run_pyscf(
            job,
            workdir,
            &["H".into(), "H".into()],
            &[0.0, 0.0, 0.0, 0.0, 0.0, 0.74],
            0,
            1,
            &MethodConfig::default(),
            ExcitedStateOptions::default(),
            false,
            None,
            &AtomicBool::new(false),
            &Mutex::new(None),
        )
    };
    for (marker, expected) in [
        ("0", Some(false)),
        ("1", Some(true)),
        ("", None),
        ("invalid", None),
    ] {
        let status = if marker.is_empty() {
            String::new()
        } else {
            format!("printf 'BEAVYR_GEOMETRY_CONVERGED {marker}\\n'\n")
        };
        fs::write(&binary, format!("#!/bin/sh\nprintf '2\\ncycle 1 E = -1.1 Eh\\nH 0 0 0\\nH 0 0 0.74\\n' > trajectory.xyz\n{status}printf 'SCF converged\\nenergy -1.1\\n'\n")).unwrap();
        fs::set_permissions(&binary, fs::Permissions::from_mode(0o700)).unwrap();
        let workdir = directory.join(format!("run-{marker}"));
        for job in [JobType::Optimize, JobType::OptimizeThenFrequencies] {
            let output = run(job, &workdir);
            if let Some(converged) = expected {
                let output = output.unwrap_or_else(|e| panic!("{e}"));
                assert_eq!(output.converged, Some(converged));
                assert_eq!(output.energetics.final_energy_hartree, Some(-1.1));
                assert!(output.trajectory_text.is_some());
            } else {
                assert!(output
                    .err()
                    .unwrap()
                    .contains("without a valid geometry-optimization convergence status"));
            }
        }
        let output = run(JobType::SinglePoint, &workdir).unwrap();
        assert_eq!(
            output.converged, None,
            "SCF-only jobs have no geometry convergence status"
        );
    }
    fs::remove_dir_all(directory).unwrap();
}
