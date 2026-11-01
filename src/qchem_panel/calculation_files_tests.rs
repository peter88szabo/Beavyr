//! File lifecycle tests: incremental output, complete exports, and real panel submissions.
use super::*;
use crate::qchem_panel::{calculation_export as export, calculation_record as record, live_output};
use std::io::Write;
use std::sync::atomic::AtomicU64;

struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        Self(
            create_xtb_run_dir(&std::env::temp_dir().join("beavyr-calculation-files-tests"))
                .unwrap(),
        )
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn output() -> QcRunOutput {
    QcRunOutput {
        scan: None,
        job: JobType::Optimize,
        program: QcProgram::Dreiding,
        log: "Energy = 0.125 Hartree\n".into(),
        trajectory_text: Some("1\nfirst\nC 0 0 0\n1\nlast\nC 1 2 3\n".into()),
        hessian_path: None,
        molden_path: None,
        saveable: vec![],
        kept_dir: None,
        energetics: Energetics {
            final_energy_hartree: Some(0.125),
            ..Default::default()
        },
        scf: vec![],
        converged: Some(false),
    }
}

fn begin(root: &Path) {
    record::begin(
        root,
        &record::Request {
            version: env!("CARGO_PKG_VERSION"),
            program: QcProgram::Dreiding,
            method: None,
            job: JobType::Optimize,
            executable: Path::new(""),
            atoms: &["C".into(), "X".into()],
            coordinates_angstrom: &[0., 0., 0., 9., 9., 9.],
            charge: 0,
            multiplicity: 1,
            config: &MethodConfig::default(),
            excited: Default::default(),
            write_orbitals: false,
            custom_input: Some("exact input\n  keep whitespace\n"),
        },
    )
    .unwrap();
}

#[test]
fn live_output_reads_appended_nested_output_and_handles_truncation() {
    let dir = Directory::new();
    let nested = dir.0.join("evaluation_000001/energy");
    fs::create_dir_all(&nested).unwrap();
    let file = nested.join("program.stdout");
    fs::write(&file, "SCF iteration 1\n").unwrap();
    let first = live_output::read(&dir.0, None).unwrap();
    assert_eq!(first.text, "SCF iteration 1\n");
    let mut writer = fs::OpenOptions::new().append(true).open(&file).unwrap();
    writeln!(writer, "SCF iteration 2").unwrap();
    assert!(live_output::read(&dir.0, first.selected.as_deref())
        .unwrap()
        .text
        .contains("iteration 2"));
    fs::write(&file, "restarted\n").unwrap();
    assert_eq!(
        live_output::read(&dir.0, first.selected.as_deref())
            .unwrap()
            .text,
        "restarted\n"
    );
    fs::write(nested.join("program.stderr"), "error details\n").unwrap();
    assert_eq!(
        live_output::read(
            &dir.0,
            Some(Path::new("evaluation_000001/energy/program.stderr"))
        )
        .unwrap()
        .text,
        "error details\n"
    );
}

#[test]
fn live_output_is_bounded_utf8_safe_and_export_keeps_every_byte() {
    let source = Directory::new();
    let dest = Directory::new();
    let contents = format!("{}\nFinal energy: -10.5\n", "λ".repeat(400_000));
    fs::write(source.0.join("orca_job.out"), contents.as_bytes()).unwrap();
    let snapshot = live_output::read(&source.0, None).unwrap();
    assert!(snapshot.truncated);
    assert!(snapshot.text.len() <= 512 * 1024);
    assert!(!snapshot.text.contains('\u{fffd}'));
    assert!(snapshot.text.ends_with("Final energy: -10.5\n"));
    let saved = export::save(&source.0, &dest.0, &AtomicU64::new(0)).unwrap();
    assert_eq!(
        fs::read(saved.directory.join("orca_job.out")).unwrap(),
        contents.as_bytes()
    );
}

#[test]
fn live_output_rejects_paths_outside_the_calculation() {
    let dir = Directory::new();
    assert!(live_output::read(&dir.0, Some(Path::new("../outside.out"))).is_err());
    assert!(live_output::read(&dir.0, Some(Path::new("/tmp/outside.out"))).is_err());
    #[cfg(unix)]
    {
        let outside = Directory::new();
        fs::write(outside.0.join("outside.out"), "outside").unwrap();
        std::os::unix::fs::symlink(outside.0.join("outside.out"), dir.0.join("link.out")).unwrap();
        assert!(live_output::read(&dir.0, None).unwrap().files.is_empty());
        assert!(live_output::read(&dir.0, Some(Path::new("link.out"))).is_err());
    }
}

#[test]
fn complete_export_preserves_binary_files_nested_files_empty_directories_and_old_exports() {
    let source = Directory::new();
    let dest = Directory::new();
    fs::create_dir_all(source.0.join("opt/empty")).unwrap();
    for (name, content) in [
        ("job.inp", &b"! RHF\n"[..]),
        ("opt/orbitals.gbw", &[0, 255, 1, 128][..]),
        ("opt/empty.out", &b""[..]),
        ("beavyr_export_inventory_0.ron", &b"native file"[..]),
    ] {
        fs::write(source.0.join(name), content).unwrap();
    }
    let bytes = AtomicU64::new(0);
    let first = export::save(&source.0, &dest.0, &bytes).unwrap();
    let second = export::save(&source.0, &dest.0, &AtomicU64::new(0)).unwrap();
    assert_ne!(first.directory, second.directory);
    assert_eq!(first.files, 4);
    assert_eq!(first.bytes, bytes.load(Ordering::Relaxed));
    for target in [&first.directory, &second.directory] {
        assert_eq!(
            fs::read(target.join("opt/orbitals.gbw")).unwrap(),
            vec![0, 255, 1, 128]
        );
        assert!(target.join("opt/empty").is_dir());
        assert_eq!(
            fs::read_to_string(target.join("beavyr_export_inventory_0.ron")).unwrap(),
            "native file"
        );
        let inventory: Vec<(PathBuf, u64)> = ron::from_str(
            &fs::read_to_string(target.join("beavyr_export_inventory_1.ron")).unwrap(),
        )
        .unwrap();
        assert_eq!(inventory.len(), first.files);
    }
    assert_eq!(
        fs::read(source.0.join("opt/orbitals.gbw")).unwrap(),
        vec![0, 255, 1, 128]
    );
    assert!(export::save(&source.0, &source.0, &AtomicU64::new(0))
        .unwrap_err()
        .contains("outside"));
    assert!(
        export::save(&source.0, &source.0.join("opt"), &AtomicU64::new(0))
            .unwrap_err()
            .contains("outside")
    );
}

#[cfg(unix)]
#[test]
fn complete_export_marks_incomplete_copies_without_following_links() {
    let source = Directory::new();
    let dest = Directory::new();
    fs::write(source.0.join("original"), "unchanged").unwrap();
    std::os::unix::fs::symlink("original", source.0.join("z_link")).unwrap();
    let error = export::save(&source.0, &dest.0, &AtomicU64::new(0)).unwrap_err();
    assert!(error.contains("symbolic link") && error.contains("incomplete copy"));
    let target = fs::read_dir(&dest.0)
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .path();
    assert!(target.join("INCOMPLETE_EXPORT_0.txt").exists());
    assert_eq!(
        fs::read_to_string(source.0.join("original")).unwrap(),
        "unchanged"
    );
}

#[test]
fn calculation_record_preserves_settings_geometry_trajectory_and_nonconvergence() {
    let dir = Directory::new();
    begin(&dir.0);
    let record = dir.0.join(record::RECORD);
    assert!(fs::read_to_string(record.join("initial.xyz"))
        .unwrap()
        .starts_with("1\n"));
    assert_eq!(
        fs::read_to_string(record.join("custom_input.txt")).unwrap(),
        "exact input\n  keep whitespace\n"
    );
    assert!(fs::read_to_string(record.join("settings.ron"))
        .unwrap()
        .contains("multiplicity: 1"));
    let mut result = Ok(output());
    record::finish(&dir.0, &mut result, Duration::from_secs(2)).unwrap();
    assert_eq!(result.unwrap().kept_dir.as_ref(), Some(&dir.0));
    assert!(fs::read_to_string(record.join("status.txt"))
        .unwrap()
        .starts_with("Finished without convergence"));
    assert_eq!(
        fs::read_to_string(record.join("final.xyz")).unwrap().trim(),
        "1\nlast\nC 1 2 3"
    );
    assert!(fs::read_to_string(record.join("trajectory.xyz"))
        .unwrap()
        .contains("first"));
}

#[test]
fn failed_cancelled_and_incomplete_scan_records_keep_their_distinct_statuses() {
    for (error, status) in [
        ("backend failed\nSCF details", "Failed"),
        ("Cancelled.", "Cancelled"),
    ] {
        let dir = Directory::new();
        begin(&dir.0);
        record::finish(&dir.0, &mut Err(error.into()), Duration::ZERO).unwrap();
        let record = dir.0.join(record::RECORD);
        assert!(fs::read_to_string(record.join("status.txt"))
            .unwrap()
            .starts_with(status));
        assert_eq!(fs::read_to_string(record.join("error.txt")).unwrap(), error);
        assert!(record.join("initial.xyz").exists());
    }
    let dir = Directory::new();
    begin(&dir.0);
    let mut out = output();
    out.scan = Some(Arc::new(crate::qchem_panel::scan_results::ScanReport {
        program: "DREIDING".into(),
        method: "DREIDING".into(),
        driver: "Beavyr".into(),
        options: Default::default(),
        atoms: vec!["C".into()],
        points: vec![],
        complete: false,
        message: "Stopped after a failed point".into(),
    }));
    record::finish(&dir.0, &mut Ok(out), Duration::ZERO).unwrap();
    assert!(
        fs::read_to_string(dir.0.join(record::RECORD).join("status.txt"))
            .unwrap()
            .starts_with("Incomplete scan")
    );
    record::mark_cancelled(&dir.0);
    assert!(
        fs::read_to_string(dir.0.join(record::RECORD).join("status.txt"))
            .unwrap()
            .starts_with("Cancelled")
    );
}

#[test]
fn complete_export_retains_hessian_and_orbitals_from_outside_the_run() {
    let dir = Directory::new();
    let external = Directory::new();
    let destination = Directory::new();
    begin(&dir.0);
    let mut out = output();
    let hessian = dir.0.join("hessian.txt");
    let molden = external.0.join("orbitals.molden");
    fs::write(&hessian, "$hessian\n1\n$end\n").unwrap();
    fs::write(&molden, "[Molden Format]\n").unwrap();
    out.hessian_path = Some(hessian);
    out.molden_path = Some(molden);
    record::finish(&dir.0, &mut Ok(out), Duration::ZERO).unwrap();
    let saved = export::save(&dir.0, &destination.0, &AtomicU64::new(0)).unwrap();
    assert_eq!(
        fs::read_to_string(
            saved
                .directory
                .join("beavyr_record/artifact_001_orbitals.molden")
        )
        .unwrap(),
        "[Molden Format]\n"
    );
    assert!(saved.directory.join("hessian.txt").exists());
}

#[test]
fn panel_submissions_retain_energy_only_and_failed_runs_after_polling() {
    AsyncComputeTaskPool::get_or_init(bevy::tasks::TaskPool::default);
    for fails in [false, true] {
        let molecule = crate::molecule::Molecule::from_xyz(if fails {
            include_str!("../../tests/fixtures/ferrocene_builder.xyz")
        } else {
            "1\ncarbon\nC 0 0 0\n"
        });
        let mut task = QcRunTask::default();
        task.start(
            QcProgram::Dreiding,
            None,
            JobType::SinglePoint,
            PathBuf::new(),
            molecule.atoms,
            molecule.pos.iter().flat_map(|p| [p.x as f64, p.y as f64, p.z as f64]).collect(),
            0,
            1,
            MethodConfig::default(),
            Default::default(),
            false,
            None,
        );
        let dir = Directory(task.live_output.root.clone().unwrap());
        let mut app = App::new();
        app.insert_resource(task)
            .insert_resource(crate::molecule::Molecule::from_xyz("1\ncurrent\nC 0 0 0\n"))
            .init_resource::<crate::qchem_panel::QcPanelState>()
            .init_resource::<crate::ui_layout::UiLayout>()
            .init_resource::<crate::trajectory::TrajectoryState>()
            .init_resource::<crate::settings::MolSettings>()
            .init_resource::<crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask>()
            .init_resource::<crate::qchem_interfaces::xtb_optimize::XtbPanelState>()
            .init_resource::<crate::qchem_interfaces::xtb_freq::XtbFreqPanelState>()
            .init_resource::<crate::qchem_interfaces::xtb_freq::XtbFrequencyTask>()
            .init_resource::<crate::uvvis::UvVisState>()
            .init_resource::<crate::orbitals::OrbitalState>()
            .add_message::<crate::events::MoleculeChanged>()
            .add_systems(Update, poll_qc_run);
        let start = Instant::now();
        loop {
            app.update();
            if !app.world().resource::<QcRunTask>().is_running() {
                break;
            }
            assert!(start.elapsed() < Duration::from_secs(10));
            std::thread::sleep(Duration::from_millis(5));
        }
        let task = app.world().resource::<QcRunTask>();
        assert_eq!(task.last_is_error, fails, "{:?}", task.last_message);
        assert_eq!(task.live_output.root.as_ref(), Some(&dir.0));
        assert!(dir.0.join(record::RECORD).join("settings.ron").is_file());
        let status = fs::read_to_string(dir.0.join(record::RECORD).join("status.txt")).unwrap();
        assert!(
            status.starts_with(if fails { "Failed" } else { "Completed" }),
            "{status}"
        );
        if !fails {
            assert_eq!(task.last.as_ref().unwrap().kept_dir.as_ref(), Some(&dir.0));
        }
        let dest = Directory::new();
        assert!(
            export::save(&dir.0, &dest.0, &AtomicU64::new(0))
                .unwrap()
                .files
                >= 5
        );
    }
}
