use super::*;
use crate::qchem_interfaces::{job, psi4_run, pyscf_run};

fn water() -> (Vec<String>, Vec<f64>, MethodConfig) {
    let mut config = MethodConfig::default();
    config.basis = "sto-3g".into();
    config.nproc = 2;
    config.memory_mb = 1000;
    config.active_space = correlated::ActiveSpace {
        electrons: 4,
        orbitals: 3,
    };
    (
        vec!["O".into(), "H".into(), "H".into()],
        vec![0., 0., 0., 0.758602, 0., 0.504284, -0.758602, 0., 0.504284],
        config,
    )
}

#[test]
fn correlated_choices_are_available_and_new_jobs_are_energy_only() {
    for program in [QcProgram::PySCF, QcProgram::Orca, QcProgram::Behemoth] {
        for method in [QcMethod::CcsdT, QcMethod::Casscf, QcMethod::Nevpt2] {
            assert!(job::methods_for(program).contains(&method));
            assert!(
                correlated::unsupported_job(program, Some(method), JobType::SinglePoint).is_none()
            );
            if program != QcProgram::Orca || method != QcMethod::CcsdT {
                for job in [
                    JobType::Gradient,
                    JobType::Optimize,
                    JobType::Frequencies,
                    JobType::TdDft,
                ] {
                    assert!(correlated::unsupported_job(program, Some(method), job).is_some());
                }
            }
        }
    }
    for program in [QcProgram::Psi4, QcProgram::Psi4Py] {
        assert!(job::methods_for(program).contains(&QcMethod::CcsdT));
        assert!(!job::methods_for(program).contains(&QcMethod::Casscf));
    }
}

#[test]
fn correlated_settings_preserve_loading_of_older_projects() {
    let config = MethodConfig::default();
    let saved = ron::to_string(&config).unwrap();
    let old = saved.replace("active_space:(electrons:0,orbitals:0),", "");
    assert_ne!(old, saved);
    assert_eq!(ron::from_str::<MethodConfig>(&old).unwrap(), config);
    let (_, _, with_space) = water();
    assert_eq!(ron::from_str::<MethodConfig>(&ron::to_string(&with_space).unwrap()).unwrap(), with_space);
}

#[test]
fn correlated_active_space_rejects_bad_electron_and_spin_counts() {
    let (mut atoms, _, config) = water();
    atoms.push("X".into());
    let validate = |active: &correlated::ActiveSpace, mult| {
        correlated::validate_active_space(QcMethod::Casscf, active, &atoms, 0, mult)
    };
    assert!(validate(&config.active_space, 1).is_ok());
    assert!(validate(&Default::default(), 1).is_err());
    assert!(validate(&config.active_space, 2).is_err());
    assert!(validate(
        &correlated::ActiveSpace {
            electrons: 12,
            orbitals: 7
        },
        1
    )
    .is_err());
    assert!(validate(
        &correlated::ActiveSpace {
            electrons: 4,
            orbitals: 1
        },
        1
    )
    .is_err());
    assert!(validate(
        &correlated::ActiveSpace {
            electrons: 3,
            orbitals: 3
        },
        2
    )
    .is_err());
}

#[test]
fn correlated_energy_parser_never_substitutes_a_reference_energy() {
    for method in [QcMethod::CcsdT, QcMethod::Casscf, QcMethod::Nevpt2] {
        assert!(parse_energy(method, "Final HF energy: -75.0\n").is_err());
    }
    assert_eq!(
        parse_energy(QcMethod::CcsdT, "Total CCSD(T) energy: -75.1 Eh\n").unwrap(),
        -75.1
    );
    assert_eq!(
        parse_energy(QcMethod::CcsdT, "Total UCCSD(T) energy: -75.2 Eh\n").unwrap(),
        -75.2
    );
    assert_eq!(
        parse_energy(QcMethod::Casscf, "CASSCF summary\n total energy -75.3 Eh\n").unwrap(),
        -75.3
    );
    assert!(parse_energy(
        QcMethod::Casscf,
        "CASSCF summary\n total energy -75.3 Eh\n NOT converged"
    )
    .is_err());
    assert_eq!(parse_energy(QcMethod::Nevpt2,
        "CASSCF summary\n total energy -75.3 Eh\nSC-NEVPT2  root 0\n E(2) -0.1\n E(reference) -75.3\n E(total) -75.4\n").unwrap(), -75.4);
    assert!(parse_energy(QcMethod::Nevpt2, "SC-NEVPT2 root 0\n E(total) NaN\n").is_err());
}

#[test]
fn correlated_inputs_use_selected_method_and_space() {
    let (atoms, coords, config) = water();
    for method in [QcMethod::CcsdT, QcMethod::Casscf, QcMethod::Nevpt2] {
        let script = pyscf_run::panel_script_text(
            JobType::SinglePoint,
            method,
            0,
            1,
            &config,
            Default::default(),
            true,
            &atoms,
            &coords,
        );
        assert_eq!(
            script.matches("print(\"energy %.12f\" % energy)").count(),
            1
        );
        assert!(!script.contains("dft.RKS"));
        if method == QcMethod::CcsdT {
            assert!(script.contains("correlation.ccsd_t()"));
        } else {
            assert!(script.contains("active_electrons = 4"));
            assert!(script.contains("active_orbitals = 3"));
            assert!(script.contains("if not cas.converged:"));
        }
        if method == QcMethod::Nevpt2 {
            assert!(script.contains("mrpt.NEVPT(cas).kernel()"));
        }
        let text = input_text(method, 0, 1, &config, false);
        assert!(text.contains(if method == QcMethod::CcsdT {
            "HF CCSD(T)"
        } else {
            "HF CASSCF[4e,3o]"
        }));
    }
    let orca = orca_run::panel_input_text(
        JobType::SinglePoint,
        QcMethod::Nevpt2,
        0,
        1,
        &config,
        Default::default(),
        "",
        &orca_run::xyz_body(&atoms, &coords),
    );
    assert!(orca.contains("nel 4\n  norb 3\n  mult 1\n  nroots 1\n  PTMethod SC_NEVPT2"));
    let psi4 = psi4_run::script_text(
        JobType::SinglePoint,
        QcMethod::CcsdT,
        0,
        1,
        &config,
        Default::default(),
        false,
        &atoms,
        &coords,
    );
    assert!(psi4.to_ascii_lowercase().contains("ccsd(t)"));
}

fn scratch(label: &str) -> PathBuf {
    let id = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!("beavyr-correlated-{label}-{id}"));
    fs::create_dir_all(&path).unwrap();
    path
}

fn execute(binary: &str, input: &str, directory: &Path) -> String {
    let output = std::process::Command::new(binary)
        .arg(input)
        .current_dir(directory)
        .env("OMP_NUM_THREADS", "2")
        .env("OPENBLAS_NUM_THREADS", "2")
        .env("MKL_NUM_THREADS", "2")
        .output()
        .unwrap();
    fs::write(directory.join("stdout.log"), &output.stdout).unwrap();
    fs::write(directory.join("stderr.log"), &output.stderr).unwrap();
    assert!(
        output.status.success(),
        "Calculation failed: {}\n{}\n{}",
        directory.display(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

#[test]
#[ignore = "requires BEAVYR_TEST_PYSCF interpreter with PySCF"]
fn correlated_native_pyscf_energies() {
    let binary = std::env::var("BEAVYR_TEST_PYSCF").unwrap();
    let (atoms, coords, config) = water();
    for method in [QcMethod::CcsdT, QcMethod::Casscf, QcMethod::Nevpt2] {
        let dir = scratch("pyscf");
        fs::write(
            dir.join("job.py"),
            pyscf_run::panel_script_text(
                JobType::SinglePoint,
                method,
                0,
                1,
                &config,
                Default::default(),
                true,
                &atoms,
                &coords,
            ),
        )
        .unwrap();
        let log = execute(&binary, "job.py", &dir);
        let energy = pyscf_run::parse_energy(&log).unwrap();
        assert!((-76.0..-74.0).contains(&energy), "{energy}");
        assert!(dir.join("orbitals.molden.input").is_file());
        println!(
            "PySCF {}: {energy:.12} Eh; {}",
            method.label(),
            dir.display()
        );
    }
}

#[test]
#[ignore = "requires BEAVYR_TEST_BEHEMOTH executable"]
fn correlated_native_behemoth_energies() {
    let binary = std::env::var("BEAVYR_TEST_BEHEMOTH").unwrap();
    let (atoms, coords, config) = water();
    for method in [QcMethod::CcsdT, QcMethod::Casscf, QcMethod::Nevpt2] {
        let dir = scratch("behemoth");
        let result = run_behemoth(
            method,
            Path::new(&binary),
            &dir,
            &atoms,
            &coords,
            0,
            1,
            &config,
            true,
            &AtomicBool::new(false),
            &Mutex::new(None),
        )
        .unwrap();
        let energy = result.energetics.final_energy_hartree.unwrap();
        assert!((-76.0..-74.0).contains(&energy), "{energy}");
        assert!(
            result.molden_path.is_some(),
            "No Molden result: {}",
            dir.display()
        );
        println!(
            "Behemoth {}: {energy:.12} Eh; {}",
            method.label(),
            dir.display()
        );
    }
}

#[test]
#[ignore = "requires BEAVYR_TEST_ORCA executable"]
fn correlated_native_orca_energies() {
    let binary = std::env::var("BEAVYR_TEST_ORCA").unwrap();
    let (atoms, coords, mut config) = water();
    config.nproc = 1;
    config.basis = "def2-SVP".into();
    for method in [QcMethod::CcsdT, QcMethod::Casscf, QcMethod::Nevpt2] {
        let dir = scratch("orca");
        fs::write(
            dir.join("job.inp"),
            orca_run::panel_input_text(
                JobType::SinglePoint,
                method,
                0,
                1,
                &config,
                Default::default(),
                "",
                &orca_run::xyz_body(&atoms, &coords),
            ),
        )
        .unwrap();
        let log = execute(&binary, "job.inp", &dir);
        assert!(log.contains("ORCA TERMINATED NORMALLY"));
        let energy = orca_run::parse_final_energy(&log).unwrap();
        assert!((-77.0..-74.0).contains(&energy), "{energy}");
        println!(
            "ORCA {}: {energy:.12} Eh; {}",
            method.label(),
            dir.display()
        );
    }
}

#[test]
fn ccsd_t_reference_is_chosen_per_spin_state_in_every_program() {
    use correlated::Reference::*;
    // Defaults: RHF closed shell, UHF open shell (ORCA's own default).
    assert_eq!(correlated::reference(None, 1), Rhf);
    assert_eq!(correlated::reference(None, 2), Uhf);
    assert_eq!(correlated::reference(Some(Rohf), 3), Rohf);
    assert_eq!(correlated::reference(Some(Uhf), 1), Uhf);
    // A choice that no longer fits the spin state falls back to the default.
    assert_eq!(correlated::reference(Some(Rohf), 1), Rhf);
    assert_eq!(correlated::reference(Some(Rhf), 2), Uhf);
    assert_eq!(correlated::Reference::options(1), [Rhf, Uhf]);
    assert_eq!(correlated::Reference::options(2), [Uhf, Rohf]);
    for program in [QcProgram::Orca, QcProgram::Mrcc, QcProgram::Behemoth] {
        assert!(correlated::reference_choice(program, Some(QcMethod::CcsdT)));
        assert!(!correlated::reference_choice(program, Some(QcMethod::Mp2)));
    }
    assert!(!job::methods_for(QcProgram::Mrcc).contains(&QcMethod::UccsdT));

    let radical: Vec<String> = vec!["O".into(), "H".into()];
    let xyz = [0., 0., 0., 0., 0., 0.97];
    let mut config = MethodConfig { basis: "cc-pVDZ".into(), ..Default::default() };
    for (choice, orca, mrcc, behemoth) in [
        (None, "! UHF CCSD(T)", "scftype=UHF", "method = UHF CCSD(T)"),
        (Some(Rohf), "! ROHF CCSD(T)", "scftype=ROHF", "method = ROHF CCSD(T)"),
    ] {
        config.cc_reference = choice;
        let o = orca_run::panel_input_text(JobType::SinglePoint, QcMethod::CcsdT, 0, 2, &config,
            Default::default(), "", &orca_run::xyz_body(&radical, &xyz));
        assert!(o.contains(orca), "{o}");
        let m = crate::qchem_interfaces::mrcc::input(JobType::SinglePoint, QcMethod::CcsdT, 0, 2,
            &config, Default::default(), false, &radical, &xyz).unwrap();
        assert!(m.contains(mrcc), "{m}");
        assert!(input_text(QcMethod::CcsdT, 0, 2, &config, false).contains(behemoth));
    }
    // Closed-shell UHF on request.
    config.cc_reference = Some(Uhf);
    assert!(input_text(QcMethod::CcsdT, 0, 1, &config, false).contains("method = UHF CCSD(T)"));
}

#[test]
fn correlated_electrons_set_the_frozen_core_in_orca_and_mrcc() {
    let (atoms, coords, mut config) = water();
    // Default: the programs' own frozen core, nothing written.
    assert_eq!(correlated::frozen_electrons(None, &atoms, 0, 1), Ok(None));
    let body = orca_run::xyz_body(&atoms, &coords);
    let orca = |c: &MethodConfig| orca_run::panel_input_text(JobType::SinglePoint, QcMethod::CcsdT, 0, 1, c, Default::default(), "", &body);
    assert!(!orca(&config).contains("FrozenCore"));
    config.basis = "cc-pCVDZ".into();
    let mrcc = |c: &MethodConfig| crate::qchem_interfaces::mrcc::input(JobType::SinglePoint, QcMethod::CcsdT, 0, 1, c, Default::default(), false, &atoms, &coords).unwrap();
    assert!(mrcc(&config).contains("core=frozen\n"));
    // All 10 electrons: no frozen core.
    config.correlated_electrons = Some(10);
    assert!(orca(&config).contains("%method\n  FrozenCore FC_NONE\nend\n"));
    assert!(mrcc(&config).contains("core=0\n"));
    // 8 correlated: 2 electrons (one spatial orbital) frozen.
    config.correlated_electrons = Some(8);
    assert!(orca(&config).contains("%method\n  FrozenCore -2\nend\n"));
    assert!(mrcc(&config).contains("core=1\n"));
    // Impossible counts are refused.
    assert!(correlated::frozen_electrons(Some(9), &atoms, 0, 1).is_err());
    assert!(correlated::frozen_electrons(Some(12), &atoms, 0, 1).is_err());
    assert!(correlated::frozen_electrons(Some(0), &atoms, 0, 1).is_err());
    assert_eq!(correlated::frozen_electrons(Some(9), &atoms, 1, 2), Ok(Some(0)));
    // Only correlated methods of ORCA and MRCC offer it.
    assert!(correlated::correlated_electrons_choice(QcProgram::Orca, Some(QcMethod::DlpnoCcsdT)));
    assert!(correlated::correlated_electrons_choice(QcProgram::Mrcc, Some(QcMethod::Ccsdt)));
    assert!(!correlated::correlated_electrons_choice(QcProgram::Orca, Some(QcMethod::Dft)));
    assert!(!correlated::correlated_electrons_choice(QcProgram::Mrcc, Some(QcMethod::Hf)));
    assert!(!correlated::correlated_electrons_choice(QcProgram::PySCF, Some(QcMethod::CcsdT)));
    // A stale count is ignored for a method that has no frozen core.
    let hf = crate::qchem_interfaces::mrcc::input(JobType::SinglePoint, QcMethod::Hf, 0, 1, &MethodConfig { correlated_electrons: Some(9), ..config.clone() }, Default::default(), false, &atoms, &coords).unwrap();
    assert!(hf.contains("core=frozen\n"));
}

#[test]
fn mrcc_offers_the_manuals_functionals() {
    use crate::qchem_interfaces::mrcc::FUNCTIONALS;
    let names: Vec<&str> = FUNCTIONALS.iter().map(|(n, _)| *n).collect();
    for f in ["B3LYP", "wB97M-V", "SCAN0", "DSDPBEP86", "RS-PBE-P86", "B3LYP-D3", "dRPA75"] {
        assert!(names.contains(&f), "{f}");
    }
    let mut unique = names.clone();
    unique.sort();
    unique.dedup();
    assert_eq!(unique.len(), names.len());
    assert!(!names.contains(&"M06"), "M06 is not in the MRCC manual");
}
