use super::*;

const TABLE: &str = "excitation energies, transition moments and TDA amplitudes\nstate    eV      nm       fL         Rv\n    1    4.000   310.0    0.0000    0.0\n    2    5.000   248.0    0.2500    0.0\nsum rule f_L 0.25\n";

#[test]
fn grimme_reader_uses_the_absorption_columns_and_keeps_dark_states() {
    let parsed = parse_stda(TABLE, Path::new("spectrum.out")).unwrap();
    assert_eq!(parsed.program, "sTDA (Grimme)");
    assert_eq!(parsed.states.len(), 2);
    assert_eq!(parsed.states[0].oscillator_strength, Some(0.0));
    assert_eq!(parsed.states[1].oscillator_strength, Some(0.25));
    assert_eq!(parsed.states[1].energy_ev, 5.0);
    assert!((parsed.states[1].energy_cm1 - 40327.71968674605).abs() < 1e-6);
    assert!(parse_stda("sTDA calculation failed", Path::new("error.out")).is_err());
    assert!(parse_stda(&TABLE.replace("0.2500", "NaN"), Path::new("bad.out")).is_err());
}

#[test]
fn legacy_excited_state_settings_load_with_spectrum_defaults() {
    let options: ExcitedStateOptions = ron::from_str("(states:12,tda:false)").unwrap();
    assert_eq!(options.states, 12);
    assert!(!options.tda);
    assert_eq!(options.energy_window_ev, 10);
    assert_eq!(options.spectrum_method, ExcitedStateMethod::StdaXtbGfn1);
}

#[cfg(unix)]
mod process {
    use super::*;
    use std::os::unix::fs::PermissionsExt;
    struct Scratch(PathBuf);
    impl Scratch {
        fn new() -> Self {
            Self(
                create_xtb_run_dir(&std::env::temp_dir().join("beavyr-spectrum-adapter-tests"))
                    .unwrap(),
            )
        }
        fn executable(&self, name: &str, script: &str) -> PathBuf {
            let path = self.0.join(name);
            fs::write(&path, format!("#!/bin/sh\nset -eu\n{script}\n")).unwrap();
            fs::set_permissions(&path, fs::Permissions::from_mode(0o700)).unwrap();
            path
        }
        fn run(
            &self,
            binary: &Path,
            generator: &Path,
            program: QcProgram,
            options: ExcitedStateOptions,
            cancel: &AtomicBool,
            child: &Mutex<Option<Child>>,
        ) -> Result<QcRunOutput, String> {
            let atoms = ["C", "O", "H", "H"].map(str::to_owned);
            run_with_generator(
                program,
                JobType::TdDft,
                binary,
                &self.0,
                &atoms,
                &[0., 0., 0., 1.2, 0., 0., -0.6, 0.9, 0., -0.6, -0.9, 0.],
                -2,
                1,
                &MethodConfig::default(),
                options,
                cancel,
                child,
                generator.to_str().unwrap(),
            )
        }
    }
    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn standalone_runs_user_paths_in_sequence_with_charge_and_no_geometry_flags() {
        let scratch = Scratch::new();
        let generator = scratch.executable(
            "generator with spaces",
            r#"
[ "$#" = 1 ]
[ "$1" = input.xyz ]
[ "$(cat .CHRG)" = -2 ]
[ "$(cat .UHF)" = 0 ]
[ "$OMP_NUM_THREADS" = 1 ]
[ -s input.xyz ]
printf 'wavefunction' > wfn.xtb
"#,
        );
        fs::write(scratch.0.join("expected-spectrum"), TABLE).unwrap();
        let binary = scratch.executable(
            "stda with spaces",
            r#"
[ -s wfn.xtb ]
[ "$#" = 3 ]
[ "$1" = -xtb ]
[ "$2" = -e ]
[ "$3" = 12 ]
cat expected-spectrum
printf 'spectrum' > tda.dat
"#,
        );
        let output = scratch
            .run(
                &binary,
                &generator,
                QcProgram::StdaXtb,
                ExcitedStateOptions {
                    energy_window_ev: 12,
                    ..Default::default()
                },
                &AtomicBool::new(false),
                &Mutex::new(None),
            )
            .unwrap();
        assert!(output.trajectory_text.is_none());
        assert!(output.hessian_path.is_none());
        assert!(output.converged.is_none());
        assert_eq!(
            parse_stda(&output.log, Path::new("run"))
                .unwrap()
                .states
                .len(),
            2
        );
    }

    #[test]
    fn failed_generator_never_starts_stda_and_retains_diagnostics() {
        let scratch = Scratch::new();
        let generator = scratch.executable("bad-generator", "echo missing-parameters >&2\nexit 7");
        let binary = scratch.executable("stda", "touch should-not-run");
        let error = scratch
            .run(
                &binary,
                &generator,
                QcProgram::StdaXtb,
                Default::default(),
                &AtomicBool::new(false),
                &Mutex::new(None),
            )
            .err()
            .unwrap();
        assert!(error.contains("missing-parameters"), "{error}");
        assert!(!scratch.0.join("should-not-run").exists());
        assert!(scratch.0.join("xtb4stda.stderr").exists());
    }

    #[test]
    fn cancellation_stops_the_active_generator_and_never_starts_stda() {
        let scratch = Scratch::new();
        let generator = scratch.executable("slow-generator", "exec sleep 30");
        let binary = scratch.executable("stda", "touch should-not-run");
        let cancel = AtomicBool::new(false);
        let child = Mutex::new(None);
        std::thread::scope(|scope| {
            let result = scope.spawn(|| {
                scratch.run(
                    &binary,
                    &generator,
                    QcProgram::StdaXtb,
                    Default::default(),
                    &cancel,
                    &child,
                )
            });
            let deadline = Instant::now() + Duration::from_secs(5);
            while child.lock().unwrap().is_none() && Instant::now() < deadline {
                std::thread::sleep(Duration::from_millis(10));
            }
            cancel.store(true, Ordering::SeqCst);
            assert_eq!(result.join().unwrap().err().as_deref(), Some("Cancelled."));
        });
        assert!(child.lock().unwrap().is_none());
        assert!(!scratch.0.join("should-not-run").exists());
    }

    #[test]
    fn behemoth_spectrum_routes_use_the_requested_reference_and_response() {
        for (method, tda, reference, response, fixture) in [
            (
                ExcitedStateMethod::StdaXtbGfn1,
                true,
                "xtb",
                "--stda",
                include_str!("../../tests/fixtures/behemoth_ch2o_stda.log"),
            ),
            (
                ExcitedStateMethod::StdaTasi,
                true,
                "tasi",
                "--stda",
                include_str!("../../tests/fixtures/behemoth_ch2o_stda_tasi.log"),
            ),
            (
                ExcitedStateMethod::StdaDft,
                true,
                "rks",
                "--stda",
                include_str!("../../tests/fixtures/behemoth_ch2o_stda.log"),
            ),
            (
                ExcitedStateMethod::StdaDft,
                false,
                "rks",
                "--stddft",
                include_str!("../../tests/fixtures/behemoth_ch2o_stddft.log"),
            ),
        ] {
            let scratch = Scratch::new();
            fs::write(scratch.0.join("fixture.out"), fixture).unwrap();
            let binary =
                scratch.executable("behemoth", "printf '%s\\n' \"$@\" > args\ncat fixture.out");
            let output = scratch
                .run(
                    &binary,
                    Path::new("unused"),
                    QcProgram::Behemoth,
                    ExcitedStateOptions {
                        spectrum_method: method,
                        tda,
                        ..Default::default()
                    },
                    &AtomicBool::new(false),
                    &Mutex::new(None),
                )
                .unwrap();
            let args = fs::read_to_string(scratch.0.join("args")).unwrap();
            assert!(args.contains(&format!("--method\n{reference}\n")), "{args}");
            assert!(args.lines().any(|arg| arg == response), "{args}");
            assert!(!args
                .lines()
                .any(|arg| ["--opt", "--hessian"].contains(&arg)));
            if method == ExcitedStateMethod::StdaXtbGfn1 {
                assert!(args.contains("--stda-orbitals\ngfn1"));
            }
            assert!(crate::uvvis::detect_and_parse(&output.log, Path::new("run")).is_ok());
        }
    }
}

#[test]
fn behemoth_tasi_ignores_dummies_but_rejects_unparameterized_elements() {
    let options = ExcitedStateOptions {
        spectrum_method: ExcitedStateMethod::StdaTasi,
        ..Default::default()
    };
    assert!(validate(
        QcProgram::Behemoth,
        options,
        &MethodConfig::default(),
        1,
        &["C".into(), "X".into()]
    )
    .is_ok());
    assert!(validate(
        QcProgram::Behemoth,
        options,
        &MethodConfig::default(),
        1,
        &["Fe".into()]
    )
    .unwrap_err()
    .contains("Fe"));
}

#[test]
#[ignore = "requires an installed Behemoth; set BEAVYR_TEST_BEHEMOTH"]
fn real_behemoth_spectrum_through_general_calculator() {
    let binary =
        PathBuf::from(std::env::var("BEAVYR_TEST_BEHEMOTH").expect("set BEAVYR_TEST_BEHEMOTH"));
    let atoms = ["C", "O", "H", "H"].map(str::to_owned);
    for (method, tda) in [
        (ExcitedStateMethod::StdaXtbGfn1, true),
        (ExcitedStateMethod::StdaTasi, true),
        (ExcitedStateMethod::StdaDft, true),
        (ExcitedStateMethod::StdaDft, false),
    ] {
        let dir =
            create_xtb_run_dir(&std::env::temp_dir().join("beavyr-real-spectrum-tests")).unwrap();
        let result = super::super::run_cancellable(
            QcProgram::Behemoth,
            None,
            JobType::TdDft,
            &binary,
            &dir,
            &atoms,
            &[0., 0., 0., 1.2, 0., 0., -0.6, 0.9, 0., -0.6, -0.9, 0.],
            0,
            1,
            &MethodConfig {
                nproc: 2,
                ..Default::default()
            },
            ExcitedStateOptions {
                spectrum_method: method,
                tda,
                energy_window_ev: 15,
                states: 5,
                ..Default::default()
            },
            false,
            None,
            &AtomicBool::new(false),
            &Mutex::new(None),
        );
        match result {
            Ok(output) => {
                let spectrum =
                    crate::uvvis::detect_and_parse(&output.log, &dir.join("behemoth.stdout"))
                        .unwrap();
                assert!(!spectrum.states.is_empty());
                eprintln!("{method:?}, tda={tda}: {} roots", spectrum.states.len());
            }
            Err(error)
                if method == ExcitedStateMethod::StdaXtbGfn1 && error.contains("unstable") =>
            {
                eprintln!(
                    "{method:?}: engine rejected unstable reference, propagated correctly: {error}"
                );
            }
            Err(error) => panic!("{method:?}, tda={tda}: {error}. Files: {}", dir.display()),
        }
        fs::remove_dir_all(dir).unwrap();
    }
}
