use super::*;

fn hydrogen() -> HessianData {
    let mut h = Array2::zeros((6, 6));
    h[(0, 0)] = 0.3;
    h[(3, 3)] = 0.3;
    h[(0, 3)] = -0.3;
    h[(3, 0)] = -0.3 + 2e-6;
    from_geometry(
        vec!["H".into(), "H".into()],
        vec![0., 0., 0., 1.4, 0., 0.],
        h,
    )
    .unwrap()
}

#[test]
fn hessian_standard_preserves_matrix_geometry_isotopes_and_metadata() {
    let mut data = hydrogen();
    data.masses_amu[1] = 2.01410177812;
    data.frequency_scale_factor = Some(0.98);
    data.ir_intensities = vec![0., 0., 0., 0., 0., 17.3];
    let text = standard_text(&data).unwrap();
    let loaded = parse(&text, None).unwrap();
    assert_eq!(data.hessian, loaded.hessian);
    assert_eq!(data.coords_bohr, loaded.coords_bohr);
    assert_eq!(data.masses_amu, loaded.masses_amu);
    assert_eq!(data.ir_intensities, loaded.ir_intensities);
    assert_eq!(data.frequency_scale_factor, loaded.frequency_scale_factor);
    // Smite/numpy.loadtxt ignores comment metadata and receives exactly 3N rows.
    let numeric: Vec<_> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(numeric.len(), 6);
    assert!(numeric.iter().all(|l| l.split_whitespace().count() == 6));
}

#[test]
fn hessian_internal_pyscf_smite_and_legacy_formats_share_one_matrix() {
    let data = hydrogen();
    let xyz = (
        data.atoms.clone(),
        data.coords_bohr.iter().map(|v| v * BOHR).collect(),
    );
    let numeric = standard_text(&data)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect::<Vec<_>>()
        .join("\n");
    for text in [
        numeric.clone(),
        format!("$hessian\n{}\n$end\n", numeric.replace('e', "D")),
        format!("6\n{}", numeric.replace('\n', " ")),
    ] {
        assert_eq!(parse(&text, Some(&xyz)).unwrap().hessian, data.hessian);
        assert!(parse(&text, None).unwrap_err().contains("no geometry"));
    }
}

#[test]
fn hessian_small_asymmetry_is_averaged_but_corruption_is_rejected() {
    let mut h = hydrogen().hessian;
    h[(0, 3)] = -0.300002;
    h[(3, 0)] = -0.3;
    let s = symmetric(&h).unwrap();
    assert_eq!(s, s.t());
    assert!((s[(0, 3)] + 0.300001).abs() < 1e-15);
    h[(3, 0)] = 0.3;
    assert!(symmetric(&h).unwrap_err().contains("large asymmetry"));
    h[(0, 0)] = f64::NAN;
    assert!(symmetric(&h).is_err());
    assert!(symmetric(&Array2::zeros((5, 5))).is_err());
}

#[test]
fn hessian_real_orca_and_behemoth_output_export_and_reload() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let orca = load(
        &root.join("examples/TS_Gamma-3-6_ZZ-S_11_Compound_1.hess"),
        None,
    )
    .unwrap();
    let copy = parse(&standard_text(&orca).unwrap(), None).unwrap();
    assert_eq!(copy.hessian, orca.hessian);
    assert_eq!(copy.ir_intensities, orca.ir_intensities);
    let text =
        std::fs::read_to_string(root.join("tests/fixtures/behemoth_water_hessian.log")).unwrap();
    let xyz = (
        vec!["O".into(), "H".into(), "H".into()],
        vec![0., 0., 0., 0., 0.76, 0.59, 0., -0.76, 0.59],
    );
    let behemoth = parse(&text, Some(&xyz)).unwrap();
    assert_eq!(behemoth.hessian, behemoth.hessian.t());
    assert_eq!(
        parse(&standard_text(&behemoth).unwrap(), None)
            .unwrap()
            .hessian,
        behemoth.hessian
    );
}

#[test]
fn hessian_real_pyscf_pbe_water_exports_as_smite_matrix() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let xyz = geometry(&std::fs::read_to_string(root.join("pyscf_pbe_water_hessian.xyz")).unwrap())
        .unwrap();
    let data = load(&root.join("pyscf_pbe_water_hessian.txt"), Some(&xyz)).unwrap();
    assert_eq!(data.hessian.dim(), (9, 9));
    assert_eq!(data.hessian, data.hessian.t());
    let text = standard_text(&data).unwrap();
    let imported = parse(&text, None).unwrap();
    assert_eq!(imported.hessian, data.hessian);
    assert_eq!(imported.coords_bohr, data.coords_bohr);
    // Optional evidence for an independent numpy/Smite load, not a runtime setting.
    if let Some(path) = std::env::var_os("BEAVYR_TEST_HESSIAN_EXPORT") {
        std::fs::write(path, text).unwrap();
    }
}
