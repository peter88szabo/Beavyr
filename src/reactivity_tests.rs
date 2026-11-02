use super::data::*;
fn h2(occupations: [f64; 2]) -> String {
    let s = (-0.5_f64).exp();
    let a = 1. / (2. * (1. + s)).sqrt();
    let b = 1. / (2. * (1. - s)).sqrt();
    format!("[Molden Format]\n[Atoms] AU\nH 1 1 0 0 0\nH 2 1 1 0 0\n[GTO]\n1 0\ns 1 1.0\n1.0 1.0\n\n2 0\ns 1 1.0\n1.0 1.0\n\n[MO]\nSym= 1\nEne= -0.5\nSpin= Alpha\nOccup= {}\n1 {a:.15}\n2 {a:.15}\nSym= 2\nEne= 0.2\nSpin= Alpha\nOccup= {}\n1 {b:.15}\n2 {:.15}\n",occupations[0],occupations[1],-b)
}
fn geometry() -> crate::trajectory::TrajectoryFrame {
    let d = crate::orbitals::molden::parse_molden(&h2([2., 0.])).unwrap();
    crate::trajectory::TrajectoryFrame {
        atoms: d.atoms,
        pos: d.positions,
        metadata: Default::default(),
    }
}
#[test]
fn frozen_orbital_finite_difference_matches_frontier_and_normalizes() {
    let finite = finite(
        geometry(),
        Sources::Finite {
            reference: Wave::Embedded(h2([2., 0.])),
            added: Wave::Embedded(h2([2., 1.])),
            removed: Wave::Embedded(h2([1., 0.])),
        },
        [-1., -1.1, -0.5],
    )
    .unwrap();
    let frontier = frontier(
        geometry(),
        Wave::Embedded(h2([2., 0.])),
        crate::orbitals::molden::Spin::Alpha,
    )
    .unwrap();
    for values in [&finite.plus, &finite.minus, &frontier.plus, &frontier.minus] {
        assert!((values.iter().sum::<f64>() - 1.).abs() < 1e-7);
        assert!((values[0] - 0.5).abs() < 1e-7);
    }
    let a = field(&finite, Descriptor::Dual, 0.4).unwrap();
    let b = field(&frontier, Descriptor::Dual, 0.4).unwrap();
    assert!(a.max_abs() > 0.01);
    assert!(a
        .values
        .iter()
        .zip(&b.values)
        .all(|(a, b)| (a - b).abs() < 1e-6));
    assert!((finite.global.eta - 0.2 * HARTREE_EV).abs() < 1e-12);
    assert!(Descriptor::Dual.value(&finite, &[0, 1]).unwrap().abs() < 1e-7);
}
#[test]
fn finite_difference_rejects_wrong_electron_counts() {
    assert!(finite(
        geometry(),
        Sources::Finite {
            reference: Wave::Embedded(h2([2., 0.])),
            added: Wave::Embedded(h2([2., 0.])),
            removed: Wave::Embedded(h2([1., 0.]))
        },
        [-1., -1.1, -0.5]
    )
    .err()
    .unwrap()
    .contains("occupations"));
}
#[test]
fn frontier_rejects_fractional_occupations() {
    assert!(frontier(
        geometry(),
        Wave::Embedded(h2([1.5, 0.5])),
        crate::orbitals::molden::Spin::Alpha
    )
    .is_err());
}
#[test]
fn saved_results_roundtrip_keeps_all_descriptors() {
    let frame = frontier(
        geometry(),
        Wave::Embedded(h2([2., 0.])),
        crate::orbitals::molden::Spin::Alpha,
    )
    .unwrap();
    let r = Results {
        version: 1,
        method: "test".into(),
        frames: vec![frame],
    };
    let text = ron::ser::to_string(&r).unwrap();
    let read: Results = ron::from_str(&text).unwrap();
    read.validate().unwrap();
    assert!((read.frames[0].global.eta - 0.35 * HARTREE_EV).abs() < 1e-12);
}

#[test]
fn real_behemoth_molden_frontiers_have_unit_mulliken_sums() {
    let text = include_str!("../tests/fixtures/behemoth_ch2o_canonicalMO.molden");
    let d = crate::orbitals::molden::parse_molden(text).unwrap();
    let geometry = crate::trajectory::TrajectoryFrame {
        atoms: d.atoms,
        pos: d.positions,
        metadata: Default::default(),
    };
    let f = frontier(
        geometry,
        Wave::Embedded(text.into()),
        crate::orbitals::molden::Spin::Alpha,
    )
    .unwrap();
    assert!((f.plus.iter().sum::<f64>() - 1.).abs() < 1e-5);
    assert!((f.minus.iter().sum::<f64>() - 1.).abs() < 1e-5);
}
#[test]
#[ignore = "requires BEAVYR_TEST_XTB; performs three real one-core xTB calculations"]
fn real_xtb_finite_difference() {
    let binary = std::path::PathBuf::from(std::env::var("BEAVYR_TEST_XTB").unwrap());
    let mol =
        crate::molecule::Molecule::from_xyz("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n");
    let frame = crate::trajectory::TrajectoryFrame {
        atoms: mol.atoms,
        pos: mol.pos,
        metadata: Default::default(),
    };
    let root = std::env::temp_dir().join(format!("beavyr-real-fukui-{}", std::process::id()));
    std::fs::create_dir(&root).unwrap();
    let result = super::job::finite(
        crate::trajectory_wavefunctions::job::Engine::Xtb {
            binary,
            gfn: 2,
            charge: 0,
            multiplicity: 1,
            cores: 1,
        },
        vec![(0, frame)],
        0,
        [1, 2, 2],
        &root,
        Default::default(),
        std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0)),
    )
    .unwrap();
    assert_eq!(result.frames.len(), 1);
    let f = &result.frames[0];
    assert!((f.plus.iter().sum::<f64>() - 1.).abs() < 1e-4);
    assert!((f.minus.iter().sum::<f64>() - 1.).abs() < 1e-4);
}

#[test]
fn behemoth_molden_matches_independent_pyscf_populations_and_spatial_values() {
    // Reference: PySCF 2.14.0 molden.load + int1e_ovlp / GTOval_sph, 2026-11-02.
    let text = include_str!("../tests/fixtures/behemoth_ch2o_canonicalMO.molden");
    let d = crate::orbitals::molden::parse_molden(text).unwrap();
    let s = overlap::Overlap::build(&d).unwrap();
    let expected = [
        [
            0.0484408127600372,
            0.68790711694004,
            0.1318260351499348,
            0.13182603514993393,
        ],
        [
            0.6363815924474174,
            0.3493685741692037,
            0.007124916691545829,
            0.007124916691545884,
        ],
    ];
    for (k, i) in [7, 8].into_iter().enumerate() {
        let pop = s.orbital(&d.alpha[i], 4).unwrap();
        for (a, b) in pop.iter().zip(expected[k]) {
            assert!((a - b).abs() < 1e-6, "{a} != {b}");
        }
    }
    let points = [[0.3, 0.2, 0.1], [1., 0.5, -1.2], [-0.6, 0.3, 2.]];
    let expected = [
        [0.04041971044449583, 0.0032449044203226303],
        [-0.04127140649784849, 0.21830566683660865],
        [0.1027031310171779, 0.17178732642093958],
    ];
    for (p, values) in points.iter().zip(expected) {
        for (k, i) in [7, 8].into_iter().enumerate() {
            let mut value = 0.;
            for (s, shell) in d.shells.iter().enumerate() {
                let c = d.positions[shell.center].as_dvec3() / crate::orbitals::grid::BOHR;
                let delta = [p[0] - c.x, p[1] - c.y, p[2] - c.z];
                let start = d.layout.shell_offsets[s];
                crate::orbitals::evaluate::accumulate_shell(
                    shell,
                    delta,
                    &d.alpha[i].coefficients[start..start + shell.n_functions()],
                    &mut Default::default(),
                    &mut value,
                );
            }
            assert!((value - values[k]).abs() < 1e-6, "{value} != {}", values[k]);
        }
    }
}

#[test]
#[ignore = "requires BEAVYR_TEST_ORCA; performs three real one-core HF calculations"]
fn real_orca_finite_difference() {
    let binary = std::path::PathBuf::from(std::env::var("BEAVYR_TEST_ORCA").unwrap());
    let mol =
        crate::molecule::Molecule::from_xyz("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n");
    let frame = crate::trajectory::TrajectoryFrame {
        atoms: mol.atoms,
        pos: mol.pos,
        metadata: Default::default(),
    };
    let root = std::env::temp_dir().join(format!("beavyr-real-orca-fukui-{}", std::process::id()));
    std::fs::create_dir(&root).unwrap();
    let mut config = crate::qchem_interfaces::method::MethodConfig::default();
    config.basis = "def2-SVP".into();
    config.nproc = 1;
    let result = super::job::finite(
        crate::trajectory_wavefunctions::job::Engine::QuantumChemistry {
            binary,
            program: crate::qchem_interfaces::program::QcProgram::Orca,
            method: Some(crate::qchem_interfaces::job::QcMethod::Hf),
            config,
            charge: 0,
            multiplicity: 1,
        },
        vec![(0, frame)],
        0,
        [1, 2, 2],
        &root,
        Default::default(),
        std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0)),
    )
    .unwrap();
    assert_eq!(result.frames.len(), 1);
    let f = &result.frames[0];
    assert!((f.plus.iter().sum::<f64>() - 1.).abs() < 1e-4);
    assert!((f.minus.iter().sum::<f64>() - 1.).abs() < 1e-4);
    assert!(f.global.ionization > 0.);
    let density = field(f, Descriptor::Dual, 0.3).unwrap();
    assert!(density.max_abs() > 0.0001);
}
