use super::*;
use crate::optimizer::traits::Objective;

fn ferrocene() -> Molecule {
    let mol = Molecule::from_xyz(include_str!("../../tests/fixtures/ferrocene_builder.xyz"));
    crate::calculation_geometry::prepare(&mol).unwrap().0
}

fn settings(binary: PathBuf, engine: SearchEngine) -> ConformerSettings {
    ConformerSettings {
        engine,
        xtb_binary: Some(binary),
        charge: 0,
        multiplicity: 1,
        thoroughness: Thoroughness::Quick,
        ncore: 1,
        reproducible: true,
        ..Default::default()
    }
}

#[cfg(unix)]
struct FakeXtb {
    root: PathBuf,
    binary: PathBuf,
}

#[cfg(unix)]
impl FakeXtb {
    fn new() -> Self {
        use std::os::unix::fs::PermissionsExt;
        let root = std::env::temp_dir().join(format!(
            "beavyr-conformer-engine-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&root).unwrap();
        let binary = root.join("xtb");
        let script = r##"#!/bin/sh
printf '%s\n' "$*" >> 'CALLS_FILE'
read -r n < "$1"
{
    printf '$grad\ncycle 1\n'
    tail -n +3 "$1"
    i=0
    while [ "$i" -lt "$n" ]; do
        printf '0.0 0.0 0.0\n'
        i=$((i+1))
    done
    printf '$end\n'
} > gradient
printf 'TOTAL ENERGY -1.0\n'
"##
        .replace("CALLS_FILE", root.join("calls").to_str().unwrap());
        std::fs::write(&binary, script).unwrap();
        std::fs::set_permissions(&binary, std::fs::Permissions::from_mode(0o700)).unwrap();
        Self { root, binary }
    }
    fn calls(&self) -> String {
        std::fs::read_to_string(self.root.join("calls")).unwrap()
    }
}

#[cfg(unix)]
impl Drop for FakeXtb {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.root);
    }
}

#[cfg(unix)]
#[test]
fn ferrocene_search_uses_xtb_without_dreiding_typing() {
    let mol = ferrocene();
    assert!(DreidingTopology::build(&mol).is_err());
    let fake = FakeXtb::new();
    let outcome = search(&mol, settings(fake.binary.clone(), SearchEngine::Gfn2)).unwrap();
    assert_eq!(outcome.atoms, mol.atoms);
    assert!(!outcome.conformers.is_empty());
    assert!(outcome.radical_warning.is_none());
    for call in fake.calls().lines() {
        assert!(call.contains("--gfn 2"), "{call}");
        assert!(call.contains("-c 0") && !call.contains("--uhf"), "{call}");
    }
}

#[cfg(unix)]
#[test]
fn single_conformer_fallback_keeps_each_selected_xtb_method() {
    let mol = ferrocene();
    let coords: Vec<_> = mol
        .pos
        .iter()
        .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
        .map(|x| x / BOHR_TO_ANGSTROM)
        .collect();
    for (engine, args) in [
        (SearchEngine::GfnFf, "--gfnff"),
        (SearchEngine::Gfn1, "--gfn 1"),
        (SearchEngine::Gfn2, "--gfn 2"),
    ] {
        let fake = FakeXtb::new();
        let outcome = single_external_conformer(
            &settings(fake.binary.clone(), engine),
            &coords,
            &mol.atoms,
            Duration::ZERO,
            1,
            1,
            0,
        )
        .unwrap();
        assert_eq!(outcome.conformers.len(), 1);
        assert!(fake.calls().lines().all(|call| call.contains(args)));
    }
}

#[test]
#[ignore = "Requires a real xTB executable supplied in BEAVYR_TEST_XTB"]
fn real_xtb_evaluates_the_ferrocene_derivative() {
    let binary = std::env::var_os("BEAVYR_TEST_XTB")
        .map(PathBuf::from)
        .expect("Set BEAVYR_TEST_XTB to the xTB executable");
    let mol = ferrocene();
    let coords: Vec<_> = mol
        .pos
        .iter()
        .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
        .map(|x| x / BOHR_TO_ANGSTROM)
        .collect();
    let workdir =
        std::env::temp_dir().join(format!("beavyr-real-xtb-ferrocene-{}", std::process::id()));
    let selected = settings(binary, SearchEngine::Gfn2);
    let SearchObjective::Xtb(mut objective) =
        objective(&selected, None, &mol.atoms, workdir.clone()).unwrap()
    else {
        panic!("The selected engine must be xTB");
    };
    let mut gradient = vec![0.0; coords.len()];
    let energy = objective.energy_gradient(&coords, &mut gradient).unwrap();
    assert!(energy.is_finite() && gradient.iter().all(|g| g.is_finite()));
    assert!(
        !objective.fell_back,
        "GFN2 should run without a force-field fallback"
    );
    println!(
        "GFN2-xTB ferrocene derivative: {energy:.12} Hartree; {} gradient components",
        gradient.len()
    );
    let _ = std::fs::remove_dir_all(workdir);
}
