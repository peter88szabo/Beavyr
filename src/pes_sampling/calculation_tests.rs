use super::*;

fn output(program: QcProgram, job: JobType, path: std::path::PathBuf) -> QcRunOutput {
    QcRunOutput {
        program,
        job,
        hessian_path: None,
        saveable: vec![("Hessian ($hessian)", path)],
        scan: None,
        log: String::new(),
        trajectory_text: None,
        molden_path: None,
        kept_dir: None,
        energetics: Default::default(),
        scf: vec![],
        converged: None,
    }
}

#[test]
fn pes_calculation_handoff_uses_each_backends_hessian_geometry() {
    let dir = std::env::temp_dir().join(format!("beavyr-pes-handoff-{}", std::process::id()));
    fs::create_dir_all(dir.join("hess")).unwrap();
    // Stretch along y. Using the unrelated x-oriented trajectory below would
    // remove this stretch as a rotation, so it catches mismatched geometries.
    let mut h = ndarray::Array2::<f64>::zeros((6, 6));
    h[(1, 1)] = 0.2;
    h[(4, 4)] = 0.2;
    h[(1, 4)] = -0.2;
    h[(4, 1)] = -0.2;
    let mut text = String::from("$hessian\n");
    for row in h.rows() {
        for x in row {
            text.push_str(&format!(" {x}"));
        }
        text.push('\n');
    }
    text.push_str("$end\n");
    let path = dir.join("hessian.txt");
    fs::write(&path, text).unwrap();
    let xyz = "2\nHessian geometry\nH 10 10 0\nH 10 10.74 0\n";
    for name in [
        "hess/input.xyz",
        "hessian_geometry.xyz",
        "final.xyz",
        "geometry.xyz",
    ] {
        fs::write(dir.join(name), xyz).unwrap();
    }
    fs::write(
        dir.join("hess/g98.out"),
        "Standard orientation:\n1 1 0 10 10 0\n2 1 0 10 10.74 0\n---\n",
    )
    .unwrap();
    for program in [
        QcProgram::Xtb,
        QcProgram::Behemoth,
        QcProgram::Dreiding,
        QcProgram::PySCF,
        QcProgram::Psi4,
        QcProgram::Psi4Py,
        QcProgram::SparrowPy,
    ] {
        for job in [
            JobType::Frequencies,
            JobType::Hessian,
            JobType::OptimizeThenFrequencies,
        ] {
            let mut out = output(program, job, path.clone());
            out.trajectory_text = Some("2\nunrelated input frame\nH 0 0 0\nH 0.74 0 0\n".into());
            if job != JobType::Hessian {
                out.hessian_path = Some(path.clone());
            }
            let f = load_calculation(&out).unwrap();
            assert_eq!(f.frequencies.len(), 1, "{program:?}, {job:?}");
            assert!((f.equilibrium[4] - f.equilibrium[1] - 0.74 / input::BOHR).abs() < 1e-12);
            let mut state = PesSampling::default();
            state.settings.workflow = Workflow::Bimolecular;
            state.settings.partners[0].rigid = true;
            state.settings.partners[0].use_overrides = true;
            adopt(&mut state, f);
            assert_eq!(state.settings.workflow, Workflow::Unimolecular);
            let p = &state.settings.partners[0];
            assert_eq!(p.source, Source::Hessian);
            assert!(!p.rigid && !p.use_overrides);
            assert_eq!(p.overrides.len(), 1);
            assert!(input::prepare(p).is_ok());
        }
    }
    fs::remove_file(dir.join("hess/g98.out")).unwrap();
    assert!(
        load_calculation(&output(QcProgram::Xtb, JobType::Frequencies, path))
            .unwrap_err()
            .contains("g98.out")
    );
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn pes_calculation_handoff_loads_self_contained_orca_hessian() {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/TS_Gamma-3-6_ZZ-S_11_Compound_1.hess");
    let f = load_calculation(&output(
        QcProgram::Orca,
        JobType::OptimizeThenFrequencies,
        path,
    ))
    .unwrap();
    assert_eq!(f.atoms.len(), 19);
    assert_eq!(f.frequencies.len(), 50);
    assert_eq!(f.excluded.len(), 1);
}

#[test]
fn pes_calculation_button_is_only_available_for_hessian_jobs() {
    for job in JobType::ALL {
        assert_eq!(
            supports_sampling(job),
            matches!(
                job,
                JobType::Frequencies | JobType::Hessian | JobType::OptimizeThenFrequencies
            )
        );
    }
    let mut out = output(QcProgram::Orca, JobType::Hessian, "missing.hess".into());
    assert!(load_calculation(&out).unwrap_err().contains("Cannot read"));
    out.saveable.clear();
    assert!(load_calculation(&out)
        .unwrap_err()
        .contains("did not produce a Hessian"));
}
