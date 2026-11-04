use super::*;
fn water() -> Vec<String> {
    vec!["O".into(), "H".into(), "H".into()]
}
fn config() -> MethodConfig {
    MethodConfig {
        basis: "cc-pVDZ".into(),
        functional: "PBE".into(),
        ..Default::default()
    }
}
fn fixture(n: &str) -> String {
    std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/mrcc")
            .join(format!("{n}.out")),
    )
    .unwrap()
}
#[test]
fn mrcc_final_energies_are_method_specific() {
    for (name, m, expected) in [
        ("hf", QcMethod::Hf, -1.1108503890136419),
        ("mp2", QcMethod::Mp2, -74.999194270367),
        ("rimp2_grad", QcMethod::RiMp2, -76.228465000238),
        ("ccsd", QcMethod::Ccsd, -75.013238313373),
        ("ccsd_t", QcMethod::CcsdT, -75.013305880922),
        ("ccsdt", QcMethod::Ccsdt, -75.013331793225),
        ("ccsdt_q", QcMethod::CcsdtQ, -75.013348110741),
        ("ccsdtq", QcMethod::Ccsdtq, -75.013355981003),
        ("casscf", QcMethod::Casscf, -74.975482245087),
        ("nevpt2", QcMethod::Nevpt2, -75.003881253632),
        ("drpa", QcMethod::Drpa, -76.327589402624),
        ("mp2_f12", QcMethod::Mp2F12, -76.357502719615),
        ("ccsd_f12", QcMethod::CcsdF12, -76.351989318982),
        ("ccsd_t_f12", QcMethod::CcsdTF12, -76.358776277629),
    ] {
        let e = parse_energy(&fixture(name), m).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!((e - expected).abs() < 1e-10, "{name}: {e}");
    }
    // Higher-order output ends with an MP2 summary: it must never replace CCSDTQ.
    assert!(parse_energy(&fixture("mp2"), QcMethod::CcsdT).is_err());
    assert!(parse_energy(
        &fixture("ccsd").replace("Normal termination of mrcc.", "aborted"),
        QcMethod::Ccsd
    )
    .is_err());
    assert!(parse_energy(
        "Total CCSD energy [au]: NaN\nNormal termination of mrcc.",
        QcMethod::Ccsd
    )
    .is_err());
}
#[test]
fn mrcc_adc_uses_final_correlated_roots_and_length_gauge_strengths() {
    let r = parse_spectrum(&fixture("adc2"), Path::new("adc.out")).unwrap();
    assert_eq!(r.states.len(), 2);
    assert!((r.states[0].energy_au - 0.29533911).abs() < 1e-10);
    assert_eq!(r.states[0].oscillator_strength, Some(0.0269575));
    assert_eq!(r.states[0].transition_dipole_au, Some([-0.3700195, 0., 0.]));
    assert_eq!(r.states[1].oscillator_strength, Some(0.));
    for name in ["scs_adc2", "sos_adc2", "td_dft"] {
        assert_eq!(
            parse_spectrum(&fixture(name), Path::new(name))
                .unwrap()
                .states
                .len(),
            2
        );
    }
    assert!(parse_spectrum(
        &fixture("adc2").replace("Normal termination of mrcc.", ""),
        Path::new("x")
    )
    .is_err());
}
#[test]
fn mrcc_gradient_requires_original_frame_and_complete_finite_rows() {
    let text = fixture("hf");
    let g = parse_gradient(&text, 2).unwrap();
    assert_eq!(g.len(), 6);
    assert!((g[2] + 0.0737829661).abs() < 1e-12);
    assert!(parse_gradient(&text, 3).is_err());
    assert!(parse_gradient(&text.replace("in original orientation ", ""), 2).is_err());
    assert!(parse_gradient(&text.replace("-0.0737829661", "NaN"), 2).is_err());
}
#[test]
fn mrcc_capabilities_restrict_jobs_and_references() {
    for m in METHODS {
        assert_eq!(
            job_problem(Some(*m), JobType::Optimize).is_none(),
            analytical_gradient(*m)
        );
        assert!(job_problem(Some(*m), JobType::OptimizeThenFrequencies).is_some());
    }
    let c = config();
    assert!(validate(JobType::Gradient, QcMethod::RiMp2, 0, 1, &c, &water()).is_ok());
    assert!(validate(JobType::Gradient, QcMethod::RiMp2, 0, 3, &c, &water()).is_err());
    assert!(validate(JobType::TdDft, QcMethod::Adc2, 0, 3, &c, &water()).is_err());
    assert!(validate(JobType::SinglePoint, QcMethod::CcsdT, 0, 2, &c, &water()).is_err());
}
#[test]
fn mrcc_core_correlation_cannot_use_a_valence_basis() {
    let mut c = config();
    c.mrcc.correlate_core = true;
    assert!(
        validate(JobType::SinglePoint, QcMethod::CcsdT, 0, 1, &c, &water())
            .unwrap_err()
            .contains("core-correlation basis")
    );
    for basis in ["cc-pCVDZ", "cc-pwCVTZ", "aug-cc-pCVQZ"] {
        c.basis = basis.into();
        assert!(validate(JobType::SinglePoint, QcMethod::CcsdT, 0, 1, &c, &water()).is_ok());
    }
    for basis in ["def2-TZVPP", "cc-pVDZ", "cc-pCVTZ-bad", "fake-cc-pCVTZ"] {
        assert!(!core_basis(basis));
    }
}
#[test]
fn mrcc_inputs_define_charge_spin_active_space_and_requested_roots() {
    let mut c = config();
    let atoms = water();
    let xyz = [0., 0., 0., 0.75, 0., 0.59, -0.75, 0., 0.59];
    c.active_space.electrons = 4;
    c.active_space.orbitals = 3;
    let inp = input(
        JobType::SinglePoint,
        QcMethod::Nevpt2,
        0,
        1,
        &c,
        Default::default(),
        false,
        &atoms,
        &xyz,
    )
    .unwrap();
    for k in [
        "calc=MCSCF",
        "docc=3",
        "mact=3",
        "nevpt=on",
        "symm=off",
        "unit=angs",
        "geom=xyz\n3\n\n",
    ] {
        assert!(inp.contains(k), "{k}")
    }
    let mut exc = ExcitedStateOptions::default();
    exc.states = 2;
    let inp = input(
        JobType::TdDft,
        QcMethod::Dft,
        0,
        1,
        &c,
        exc,
        false,
        &atoms,
        &xyz,
    )
    .unwrap();
    assert!(inp.contains("calc=TDA"));
    assert!(inp.contains("nstate=3"));
    let inp = input(
        JobType::SinglePoint,
        QcMethod::UccsdT,
        0,
        1,
        &c,
        exc,
        false,
        &atoms,
        &xyz,
    )
    .unwrap();
    assert!(inp.contains("calc=CCSD(T)"));
    assert!(inp.contains("scftype=UHF"));
    let inp = input(
        JobType::Gradient,
        QcMethod::Ccsd,
        0,
        1,
        &c,
        exc,
        false,
        &atoms,
        &xyz,
    )
    .unwrap();
    assert!(inp.contains("dens=2"));
    assert!(inp.contains("ccprog=mrcc"));
    assert!(validate_custom(&inp, &inp).is_ok());
    assert!(validate_custom(&inp, &format!("scfmaxit=500\n{inp}")).is_ok());
    assert!(validate_custom(&inp, &inp.replace("core=frozen", "core=corr")).is_err());
    assert!(validate_custom(&inp, &format!("geoopt=on\n{inp}")).is_err());
}
#[test]
fn mrcc_unconverged_energy_is_rejected_even_with_normal_termination() {
    let good = fixture("ccsd");
    let bad = format!("Convergence not achieved in 2 iterations!\n{good}");
    assert!(parse_energy(&bad, QcMethod::Ccsd).is_err());
}
#[test]
fn mrcc_dummy_atoms_are_omitted_before_input_and_validation() {
    let atoms = vec!["X".into(), "H".into(), "H".into()];
    let text = input(
        JobType::SinglePoint,
        QcMethod::Hf,
        0,
        1,
        &config(),
        Default::default(),
        false,
        &atoms,
        &[0., 0., 0., 0., 0., 0., 0., 0., 0.8],
    )
    .unwrap();
    assert!(text.contains("geom=xyz\n2\n\nH "));
    assert!(!text.lines().any(|l| l.starts_with("X ")));
}
#[test]
fn mrcc_launch_uses_the_selected_executable_without_shell_arguments() {
    let dir = std::env::temp_dir();
    let selected = std::env::current_exe().unwrap();
    let cmd = command(&selected, &dir, 2).unwrap();
    assert_eq!(cmd.get_program(), std::fs::canonicalize(selected).unwrap());
    assert_eq!(cmd.get_args().count(), 0);
    assert_eq!(cmd.get_current_dir(), Some(dir.as_path()));
    assert!(cmd
        .get_envs()
        .any(|(key, val)| key == "OMP_NUM_THREADS" && val == Some(std::ffi::OsStr::new("2"))));
}

#[test]
fn mrcc_missing_correlated_strengths_do_not_reuse_cis_strengths() {
    let text = fixture("adc2");
    let start = text
        .rfind("Excited               dipole moment in length gauge")
        .unwrap();
    let without = format!("{}\nNormal termination of mrcc.\n", &text[..start]);
    let result = parse_spectrum(&without, Path::new("no-final-strengths")).unwrap();
    assert!(result
        .states
        .iter()
        .all(|s| s.oscillator_strength.is_none()));
}
#[test]
fn mrcc_core_basis_uses_correct_hydrogen_member() {
    let mut c = config();
    c.basis = "cc-pCVDZ".into();
    c.mrcc.correlate_core = true;
    let text = input(
        JobType::SinglePoint,
        QcMethod::CcsdT,
        0,
        1,
        &c,
        Default::default(),
        false,
        &water(),
        &[0., 0., 0., 0.75, 0., 0.59, -0.75, 0., 0.59],
    )
    .unwrap();
    assert!(text.contains("basis=atomtype\nO:cc-pCVDZ\nH:cc-pvdz\n"));
    assert!(text.contains("core=corr"));
    assert!(validate_custom(&text, &text).is_ok());
    assert!(validate_custom(&text, &text.replace("O:cc-pCVDZ", "O:cc-pVDZ")).is_err());
}
