use super::*;
use crate::qchem_interfaces::job::JobType;

#[test]
fn hessian_export_preserves_native_and_uses_calculated_pyscf_geometry() {
    let dir = std::env::temp_dir().join(format!(
        "beavyr-pyscf-hessian-export-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let native = dir.join("hessian.txt");
    let mut h = ndarray::Array2::<f64>::eye(6);
    h[(0, 3)] = -0.2;
    h[(3, 0)] = -0.200002;
    let text = format!(
        "$hessian\n{}\n$end\n",
        h.rows()
            .into_iter()
            .map(|row| row
                .iter()
                .map(|v| format!("{v:.14e}"))
                .collect::<Vec<_>>()
                .join(" "))
            .collect::<Vec<_>>()
            .join("\n")
    );
    std::fs::write(&native, &text).unwrap();
    std::fs::write(
        dir.join("final.xyz"),
        "2\ncalculated geometry\nH 4 0 0\nH 4.74 0 0\n",
    )
    .unwrap();
    let mut output = QcRunOutput {
        program: QcProgram::PySCF,
        job: JobType::OptimizeThenFrequencies,
        hessian_path: Some(native.clone()),
        saveable: vec![("Hessian ($hessian)", native.clone())],
        scan: None,
        log: String::new(),
        trajectory_text: Some("2\nwrong frame\nH 0 0 0\nH 0 1 0\n".into()),
        molden_path: None,
        kept_dir: Some(dir.clone()),
        energetics: Default::default(),
        scf: vec![],
        converged: None,
    };
    prepare(&mut output).unwrap();
    assert_eq!(std::fs::read_to_string(&native).unwrap(), text);
    let first_path = output.hessian_path.clone();
    prepare(&mut output).unwrap();
    assert_eq!(
        output.hessian_path, first_path,
        "preparation must be idempotent"
    );
    let exported = output.hessian_path.as_ref().unwrap();
    assert_eq!(&output.saveable[0].1, exported);
    let loaded = hessian_io::load(exported, None).unwrap();
    assert!((loaded.coords_bohr[0] - 4. / hessian_io::BOHR).abs() < 1e-12);
    assert!((loaded.hessian[(0, 3)] + 0.200001).abs() < 1e-12);
    assert_eq!(loaded.hessian, loaded.hessian.t());
    std::fs::remove_dir_all(dir).unwrap();
}
