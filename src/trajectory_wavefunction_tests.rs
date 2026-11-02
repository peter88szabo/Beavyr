use super::*;
struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "beavyr-wave-test-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn frame() -> TrajectoryFrame {
    TrajectoryFrame {
        atoms: vec!["H".into()],
        pos: vec![Vec3::ZERO],
        metadata: Default::default(),
    }
}
pub const HYDROGEN:&str="[Molden Format]\n[Atoms] Angs\nH 1 1 0 0 0\n[GTO]\n1 0\ns 1 1.0\n1.0 1.0\n\n[MO]\nSym= 1\nEne= -0.5\nSpin= Alpha\nOccup= 1.0\n1 1.0\n";
#[test]
fn wavefunction_input_removes_dummies_and_checks_ranges() {
    let mut f = frame();
    f.atoms.push("X".into());
    f.pos.push(Vec3::ONE);
    let selected = prepared(&[f], 1, 1, 1).unwrap();
    assert_eq!(selected[0].1.atoms, ["H"]);
    assert!(prepared(&[], 1, 1, 1).is_err());
    assert!(prepared(&[frame()], 1, 1, 0).is_err());
}
#[test]
fn wavefunction_manifest_roundtrip_rejects_wrong_geometry_and_parent_paths() {
    let root = Directory::new();
    fs::write(root.0.join("h.molden"), HYDROGEN).unwrap();
    let mut manifest = Manifest {
        version: 1,
        method: "test".into(),
        charge: 0,
        multiplicity: 2,
        frames: vec![Frame {
            source_index: 0,
            geometry: frame(),
            molden: "h.molden".into(),
            spin_resolved: false,
            checksum: None,
        }],
        requested_frames: 1,
        complete: true,
    };
    let path = root.0.join(MANIFEST);
    fs::write(&path, ron::ser::to_string(&manifest).unwrap()).unwrap();
    let series = read(&path).unwrap();
    let (d, _) = load(&series, 0).unwrap();
    assert_eq!(d.atoms, ["H"]);
    let mut moved = frame();
    moved.pos[0].x = 1.;
    assert!(check_geometry(&d, &moved).is_err());
    manifest.frames[0].molden = "../outside.molden".into();
    fs::write(&path, ron::ser::to_string(&manifest).unwrap()).unwrap();
    assert!(read(&path).is_err());
}
#[cfg(unix)]
#[test]
fn wavefunction_cancel_kills_process_and_keeps_no_success_manifest() {
    use std::os::unix::fs::PermissionsExt;
    let root = Directory::new();
    let binary = root.0.join("slow-xtb");
    fs::write(&binary, "#!/bin/sh\nexec sleep 30\n").unwrap();
    fs::set_permissions(&binary, fs::Permissions::from_mode(0o700)).unwrap();
    let run = root.0.join("run");
    fs::create_dir(&run).unwrap();
    let control = Control::default();
    let c = control.clone();
    let r = run.clone();
    let task = std::thread::spawn(move || {
        calculate(
            Engine::Xtb {
                binary,
                gfn: 2,
                charge: 0,
                multiplicity: 2,
                cores: 1,
            },
            vec![(0, frame())],
            r,
            c,
        )
    });
    let deadline = std::time::Instant::now() + Duration::from_secs(5);
    while control.child.lock().unwrap().is_none() {
        assert!(std::time::Instant::now() < deadline);
        std::thread::sleep(Duration::from_millis(5));
    }
    control.cancel();
    let error = task.join().unwrap().err().unwrap();
    assert!(error.contains("cancelled"));
    assert!(!run.join(MANIFEST).exists());
    assert!(control.child.lock().unwrap().is_none());
}
#[test]
#[ignore = "requires BEAVYR_TEST_XTB pointing to an installed executable"]
fn wavefunction_real_gfn2_trajectory() {
    let binary = PathBuf::from(std::env::var("BEAVYR_TEST_XTB").expect("set BEAVYR_TEST_XTB"));
    let root = Directory::new();
    let mol =
        crate::molecule::Molecule::from_xyz("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n");
    let f = TrajectoryFrame {
        atoms: mol.atoms,
        pos: mol.pos,
        metadata: Default::default(),
    };
    let mut second = f.clone();
    second.pos[1].x += 0.02;
    let series = calculate(
        Engine::Xtb {
            binary,
            gfn: 2,
            charge: 0,
            multiplicity: 1,
            cores: 1,
        },
        vec![(0, f), (1, second)],
        root.0.clone(),
        Control::default(),
    )
    .unwrap();
    assert!(series.manifest.complete);
    let(d,_)=load(&series,0).unwrap();
    let spin=crate::orbitals::molden::Spin::Alpha;
    let track=crate::orbital_tracking::calculate(&series,0,d.homo_index(spin).unwrap(),spin,0.5,0.1,&Control::default()).unwrap();
    assert!(track.choices.iter().all(Option::is_some));assert!(track.issues.is_empty());
    let control=Control::default();let resumed=resume::run(root.0.clone(),control.clone()).unwrap();
    assert_eq!(control.reused.load(Ordering::Relaxed),2);assert!(resumed.manifest.complete);
    for i in 0..2 {
        let (d, _) = load(&series, i).unwrap();
        let v = crate::polarizability::data::compute(&d).unwrap();
        assert_eq!(v.alpha.len(), 3);
        assert!(v.alpha.iter().all(|a| a.is_finite() && *a > 0.));
        assert!(
            series.manifest.frames[i]
                .geometry
                .metadata
                .energy_hartree
                .unwrap()
                < -1.
        );
    }
}

#[test]
fn resume_verifies_saved_files_and_rejects_changed_inputs() {
    let root=Directory::new();
    let engine=Engine::Xtb{binary:"unused".into(),gfn:2,charge:0,multiplicity:2,cores:1};
    let selected=vec![(3,frame()),(5,frame())];
    resume::save(&root.0,&engine,&selected).unwrap();
    let plan=resume::read_plan(&root.0).unwrap();
    let mut done=frame();done.metadata.energy_hartree=Some(-0.5);
    fs::write(root.0.join("h.molden"),HYDROGEN).unwrap();
    let mut m=Manifest{version:1,method:engine.label(),charge:0,multiplicity:2,frames:vec![Frame{source_index:3,geometry:done,molden:"h.molden".into(),spin_resolved:false,checksum:Some(resume::checksum(HYDROGEN.as_bytes()))}],requested_frames:2,complete:false};
    let save=|m:&Manifest|fs::write(root.0.join(MANIFEST),ron::ser::to_string(m).unwrap()).unwrap();save(&m);
    assert_eq!(resume::verified(&root.0,&plan,&Control::default()).unwrap().len(),1);
    fs::write(root.0.join("h.molden"),format!("{HYDROGEN}\n")).unwrap();
    assert!(resume::verified(&root.0,&plan,&Control::default()).err().unwrap().contains("modified"));
    fs::write(root.0.join("h.molden"),HYDROGEN).unwrap();m.frames[0].source_index=4;save(&m);
    assert!(resume::verified(&root.0,&plan,&Control::default()).is_err());
    assert!(resume::save(&root.0,&engine,&selected).is_err());
    let plan_text=fs::read_to_string(root.0.join(resume::PLAN)).unwrap();fs::write(root.0.join(resume::PLAN),format!("{plan_text} ")).unwrap();assert!(resume::read_plan(&root.0).is_err());
    let lock=resume::lock(&root.0).unwrap();assert!(resume::lock(&root.0).is_err());drop(lock);assert!(resume::lock(&root.0).is_ok());
}

#[cfg(unix)]
#[test]
fn resume_retries_failed_frame_without_recomputing_or_overwriting_completed_frames(){
    use std::os::unix::fs::PermissionsExt;
    let root=Directory::new();let binary=root.0.join("mock-xtb");
    let script=format!("#!/bin/sh\ncase \"$PWD\" in *frame-000002) exit 1;; esac\nprintf 'TOTAL ENERGY -0.5 Eh\\n'\ncat > molden.input <<'MOLDEN'\n{HYDROGEN}\nMOLDEN\n");
    fs::write(&binary,script).unwrap();fs::set_permissions(&binary,fs::Permissions::from_mode(0o700)).unwrap();
    let engine=Engine::Xtb{binary,gfn:2,charge:0,multiplicity:2,cores:1};
    assert!(calculate(engine,vec![(0,frame()),(1,frame())],root.0.clone(),Control::default()).is_err());
    let before=fs::read(root.0.join("frame-000001/wavefunction.molden")).unwrap();
    let control=Control::default();let series=resume::run(root.0.clone(),control.clone()).unwrap();
    assert_eq!(control.reused.load(Ordering::Relaxed),1);assert_eq!(series.manifest.frames.len(),2);assert!(series.manifest.complete);
    assert_eq!(before,fs::read(root.0.join("frame-000001/wavefunction.molden")).unwrap());
    assert!(root.0.join("frame-000002").is_dir());assert!(root.0.join("frame-000002-retry-1/wavefunction.molden").is_file());
    assert!(!root.0.join("frame-000001-retry-1").exists());
}
