use super::*;
use std::os::unix::fs::PermissionsExt;

fn run_input(binary: &Path, directory: &Path, input: &str) -> Result<QcRunOutput, String> {
    run_cancellable(
        QcProgram::Orca,
        Some(QcMethod::Hf),
        JobType::SinglePoint,
        binary,
        directory,
        &["H".into(), "H".into()],
        &[0.0, 0.0, 0.0, 0.0, 0.0, 0.74],
        0,
        1,
        &MethodConfig {
            nproc: 4,
            ..Default::default()
        },
        ExcitedStateOptions::default(),
        false,
        Some(input),
        &AtomicBool::new(false),
        &Mutex::new(None),
    )
}

#[test]
fn configured_orca_runs_without_a_path_mpi_precheck() {
    const CHILD: &str = "BEAVYR_ORCA_DIRECT_LAUNCH_TEST";
    if std::env::var_os(CHILD).is_none() {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "qchem_panel::run::orca_launch_tests::configured_orca_runs_without_a_path_mpi_precheck",
                "--nocapture",
                "--test-threads=1",
            ])
            .env(CHILD, "1")
            .env("PATH", "/no-mpi-or-other-executables")
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        return;
    }

    let directory =
        std::env::temp_dir().join(format!("beavyr configured orca {}", std::process::id()));
    fs::create_dir_all(&directory).unwrap();
    let binary = directory.join("selected executable");
    fs::write(
        &binary,
        "#!/bin/sh\n[ \"$#\" -eq 1 ] || exit 21\n[ \"$1\" = orca_job.inp ] || exit 22\n/bin/cp \"$1\" received.inp || exit 23\nprintf 'FINAL SINGLE POINT ENERGY -1.1\\nORCA TERMINATED NORMALLY\\n'\n",
    )
    .unwrap();
    fs::set_permissions(&binary, fs::Permissions::from_mode(0o700)).unwrap();
    let input = "! RHF def2-SVP\n%pal nprocs 4 end\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n";
    let output = run_input(&binary, &directory, input)
        .unwrap_or_else(|error| panic!("The configured executable must be invoked: {error}"));
    assert_eq!(
        fs::read_to_string(directory.join("received.inp")).unwrap(),
        input
    );
    assert_eq!(output.energetics.final_energy_hartree, Some(-1.1));
    fs::remove_dir_all(directory).unwrap();
}

#[test]
#[ignore = "Runs the supplied ORCA executable on BEAVYR_ORCA_INPUT"]
fn real_orca_runs_the_supplied_input_unchanged() {
    let binary = std::env::var_os("BEAVYR_ORCA").expect("Set BEAVYR_ORCA to the executable");
    let input_path = std::env::var_os("BEAVYR_ORCA_INPUT").expect("Set BEAVYR_ORCA_INPUT");
    let input = fs::read_to_string(input_path).unwrap();
    let directory =
        std::env::temp_dir().join(format!("beavyr-orca-supplied-input-{}", std::process::id()));
    let output = run_input(Path::new(&binary), &directory, &input)
        .unwrap_or_else(|error| panic!("{error}; inspect {}", directory.display()));
    assert_eq!(
        fs::read_to_string(orca_run::job_files(&directory).input).unwrap(),
        input
    );
    assert!(output.log.contains("4 parallel MPI-processes"));
    assert!(output.energetics.final_energy_hartree.unwrap().is_finite());
    println!(
        "Unchanged four-core ORCA input completed through Beavyr: {}",
        directory.display()
    );
}
