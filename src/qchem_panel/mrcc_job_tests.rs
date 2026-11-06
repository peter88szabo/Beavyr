use super::*;
#[test]
#[ignore = "Native MRCC: set BEAVYR_MRCC_EXECUTABLE; runs on at most four assigned cores"]
fn mrcc_native_energy_spectrum_gradient_and_optimization() {
    let binary = PathBuf::from(
        std::env::var_os("BEAVYR_MRCC_EXECUTABLE").expect("Set BEAVYR_MRCC_EXECUTABLE"),
    );
    let root = std::env::temp_dir().join(format!("beavyr-mrcc-native-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    println!("MRCC native evidence: {}", root.display());
    let atoms: Vec<String> = ["O", "H", "H"].iter().map(|s| s.to_string()).collect();
    let xyz = [0., 0., 0., 0.758, 0., 0.59, -0.758, 0., 0.59];
    let mut config = MethodConfig {
        basis: "STO-3G".into(),
        functional: "PBE0".into(),
        nproc: 2,
        memory_mb: 768,
        ..Default::default()
    };
    config.active_space.electrons = 4;
    config.active_space.orbitals = 3;
    let exc = ExcitedStateOptions {
        states: 2,
        ..Default::default()
    };
    let cancel = AtomicBool::new(false);
    let child = Mutex::new(None);
    for (index, &method) in mrcc::METHODS.iter().enumerate() {
        config.basis = if method.is_f12() {
            "cc-pVDZ-F12"
        } else if matches!(
            method,
            QcMethod::RiMp2
                | QcMethod::Drpa
                | QcMethod::Adc2
                | QcMethod::ScsAdc2
                | QcMethod::SosAdc2
        ) {
            "cc-pVDZ"
        } else {
            "STO-3G"
        }
        .into();
        let job = if mrcc::adc(method) {
            JobType::TdDft
        } else {
            JobType::SinglePoint
        };
        let result = run_cancellable(
            QcProgram::Mrcc,
            Some(method),
            job,
            &binary,
            &root.join(format!("method_{index}")),
            &atoms,
            &xyz,
            0,
            1,
            &config,
            exc,
            true,
            None,
            &cancel,
            &child,
        )
        .unwrap_or_else(|e| panic!("{}: {e}", method.label()));
        assert!(result.energetics.final_energy_hartree.unwrap().is_finite());
        assert!(result.molden_path.as_ref().unwrap().is_file());
        println!(
            "PASS {}: {:?}",
            method.label(),
            result.energetics.final_energy_hartree
        );
    }
    config.basis = "STO-3G".into();
    for tda in [true, false] {
        let result = run_cancellable(
            QcProgram::Mrcc,
            Some(QcMethod::Dft),
            JobType::TdDft,
            &binary,
            &root.join(format!("tddft_{tda}")),
            &atoms,
            &xyz,
            0,
            1,
            &config,
            ExcitedStateOptions { tda, ..exc },
            false,
            None,
            &cancel,
            &child,
        )
        .unwrap();
        assert_eq!(
            mrcc::parse_spectrum(&result.log, Path::new("spectrum"))
                .unwrap()
                .states
                .len(),
            2
        );
        println!("PASS TDDFT tda={tda}");
    }
    for method in [
        QcMethod::Hf,
        QcMethod::Dft,
        QcMethod::RiMp2,
        QcMethod::Ccsd,
        QcMethod::Ccsdt,
        QcMethod::Ccsdtq,
    ] {
        config.basis = if method == QcMethod::RiMp2 {
            "cc-pVDZ"
        } else {
            "STO-3G"
        }
        .into();
        let dir = root.join(format!("gradient_{method:?}"));
        let result = run_cancellable(
            QcProgram::Mrcc,
            Some(method),
            JobType::Gradient,
            &binary,
            &dir,
            &atoms,
            &xyz,
            0,
            1,
            &config,
            exc,
            false,
            None,
            &cancel,
            &child,
        )
        .unwrap_or_else(|e| panic!("{method:?}: {e}"));
        let (energy, grad) =
            crate::qchem_interfaces::gradient::read(QcProgram::Mrcc, &dir, 3, &result.log).unwrap();
        assert_eq!(grad.len(), 9);
        assert_eq!(energy, result.energetics.final_energy_hartree.unwrap());
        println!("PASS analytical gradient {}", method.label());
    }
    config.basis = "STO-3G".into();
    let radical = vec!["O".into(), "H".into()];
    let result = run_cancellable(
        QcProgram::Mrcc,
        Some(QcMethod::CcsdT),
        JobType::SinglePoint,
        &binary,
        &root.join("open_shell_uccsd_t"),
        &radical,
        &[0., 0., 0., 0., 0., 0.97],
        0,
        2,
        &config,
        exc,
        false,
        None,
        &cancel,
        &child,
    )
    .unwrap();
    println!(
        "PASS open-shell UCCSD(T) {:?}",
        result.energetics.final_energy_hartree
    );
    // Check a derivative of an arbitrarily oriented molecule against energy differences.
    let rotated = [0.1, 0.2, 0.3, 0.49, 0.77, 0.85, -0.32, 0.88, 0.81];
    let g = run_cancellable(
        QcProgram::Mrcc,
        Some(QcMethod::Hf),
        JobType::Gradient,
        &binary,
        &root.join("fd_gradient"),
        &atoms,
        &rotated,
        0,
        1,
        &config,
        exc,
        false,
        None,
        &cancel,
        &child,
    )
    .unwrap();
    let gradient = mrcc::parse_gradient(&g.log, 3).unwrap();
    let h = 1e-4;
    let mut energies = vec![];
    for sign in [-1., 1.] {
        let mut displaced = rotated;
        displaced[5] += sign * h * 0.529177210903;
        let e = run_cancellable(
            QcProgram::Mrcc,
            Some(QcMethod::Hf),
            JobType::SinglePoint,
            &binary,
            &root.join(format!("fd_{sign}")),
            &atoms,
            &displaced,
            0,
            1,
            &config,
            exc,
            false,
            None,
            &cancel,
            &child,
        )
        .unwrap();
        energies.push(e.energetics.final_energy_hartree.unwrap());
    }
    let numerical = (energies[1] - energies[0]) / (2. * h);
    assert!(
        (numerical - gradient[5]).abs() < 2e-5,
        "numerical {numerical}, analytic {}",
        gradient[5]
    );
    println!(
        "PASS original-frame analytical derivative: analytic {}, finite difference {numerical}",
        gradient[5]
    );
    config.basis = "cc-pCVDZ".into();
    config.correlated_electrons = Some(10);
    let result = run_cancellable(
        QcProgram::Mrcc,
        Some(QcMethod::CcsdT),
        JobType::SinglePoint,
        &binary,
        &root.join("core_ccsd_t"),
        &atoms,
        &xyz,
        0,
        1,
        &config,
        exc,
        false,
        None,
        &cancel,
        &child,
    )
    .unwrap();
    assert!(result.log.lines().any(
        |l| l.contains("Number of occupied/virtual orbitals:") && l.split_once(':').is_some_and(|(_,counts)|counts.split_whitespace().next()==Some("5"))
    ));
    println!(
        "PASS core correlation {:?}",
        result.energetics.final_energy_hartree
    );
    config.basis = "STO-3G".into();
    config.correlated_electrons = None;
    config.search.minimum.max_cycles = 40;
    let atoms = vec!["H".into(), "H".into()];
    let xyz = [0., 0., 0., 0.3, 0.4, 0.8];
    let result = run_cancellable(
        QcProgram::Mrcc,
        Some(QcMethod::Hf),
        JobType::Optimize,
        &binary,
        &root.join("optimization"),
        &atoms,
        &xyz,
        0,
        1,
        &config,
        exc,
        false,
        None,
        &cancel,
        &child,
    )
    .unwrap();
    assert_eq!(result.converged, Some(true));
    assert!(result.trajectory_text.is_some());
    println!(
        "PASS Beavyr/MRCC optimization {:?}",
        result.energetics.final_energy_hartree
    );
}

#[cfg(unix)]
#[test]
fn mrcc_cancellation_kills_companion_processes() {
    let dir = std::env::temp_dir().join(format!(
        "beavyr-mrcc-cancel-{}-{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    fs::create_dir_all(&dir).unwrap();
    let cancel = Arc::new(AtomicBool::new(false));
    let flag = cancel.clone();
    let thread = std::thread::spawn(move || {
        std::thread::sleep(Duration::from_millis(100));
        flag.store(true, Ordering::SeqCst);
    });
    let mut cmd = std::process::Command::new("/bin/sh");
    cmd.args([
        "-c",
        "(sleep 0.4; printf survived > companion-survived) & wait",
    ]);
    let result = process::run(cmd, &dir, "cancel-test", &cancel);
    thread.join().unwrap();
    assert!(result.unwrap_err().contains("Cancelled"));
    std::thread::sleep(Duration::from_millis(500));
    assert!(!dir.join("companion-survived").exists());
    fs::remove_dir_all(dir).unwrap();
}
