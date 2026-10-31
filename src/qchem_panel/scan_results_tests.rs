use super::*;
use crate::qchem_interfaces::{job::JobType, program::QcProgram};
use crate::qchem_panel::run::QcRunOutput;
use crate::{molecule::Molecule, trajectory::TrajectoryState};
pub(super) fn output() -> QcRunOutput {
    QcRunOutput {
        job: JobType::Scan,
        program: QcProgram::Dreiding,
        log: String::new(),
        trajectory_text: None,
        hessian_path: None,
        molden_path: None,
        saveable: vec![],
        kept_dir: None,
        energetics: Default::default(),
        scf: vec![],
        converged: None,
        scan: Some(Arc::new(ScanReport {
            program: "DREIDING".into(),
            method: "DREIDING".into(),
            driver: "Beavyr".into(),
            options: Default::default(),
            atoms: vec!["H".into(), "H".into()],
            points: (0..3)
                .map(|i| ScanPoint {
                    target: 1. + i as f64 * 0.1,
                    actual: 1. + i as f64 * 0.1,
                    energy: -1. + i as f64 * 0.01,
                    converged: None,
                    steps: 0,
                    coordinates: vec![0., 0., 0., 1. + i as f64 * 0.1, 0., 0.],
                })
                .collect(),
            complete: true,
            message: "3 scan points completed.".into(),
        })),
    }
}
#[test]
fn scan_results_open_once_per_result_and_survive_other_jobs() {
    let mut task = QcRunTask::default();
    let mut s = WindowState::default();
    s.sync(&task);
    assert!(!s.open);
    task.last = Some(output());
    s.sync(&task);
    assert!(s.open);
    s.open = false;
    s.sync(&task);
    assert!(!s.open);
    task.last = None;
    s.sync(&task);
    assert!(s.report.is_some());
    assert!(!s.open);
    task.last = Some(output());
    s.sync(&task);
    assert!(s.open);
}
struct Fixture {
    ctx: egui::Context,
    task: QcRunTask,
    mol: Molecule,
    traj: TrajectoryState,
    actions: Actions,
    copied: Option<String>,
}
impl Fixture {
    fn frame(&mut self, events: Vec<egui::Event>) -> egui::FullOutput {
        let mut output = self.ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1200., 900.),
                )),
                events,
                ..Default::default()
            },
            |ui| {
                let (_, a) = window(ui.ctx(), &self.task, &mut self.mol, &mut self.traj);
                self.actions.loaded |= a.loaded;
                self.actions.trajectory |= a.trajectory;
                egui::Window::new("Setup")
                    .default_pos(egui::Pos2::ZERO)
                    .show(ui.ctx(), |ui| reopen_button(ui, &self.task));
            },
        );
        for command in &output.platform_output.commands {
            if let egui::OutputCommand::CopyText(text) = command {
                self.copied = Some(text.clone());
            }
        }
        output.textures_delta.clear();
        output
    }
    fn click(&mut self, label: &str) {
        self.frame(vec![]);
        let out = self.frame(vec![]);
        let pos = out
            .shapes
            .iter()
            .find_map(|s| match &s.shape {
                egui::Shape::Text(t) if t.galley.text() == label => {
                    Some(t.pos + t.galley.size() / 2.)
                }
                _ => None,
            })
            .unwrap_or_else(|| panic!("Missing {label}"));
        for pressed in [true, false] {
            self.frame(vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ]);
        }
    }
}
#[test]
fn scan_results_window_loads_selection_and_trajectory_and_reopens() {
    let mut task = QcRunTask::default();
    task.last = Some(output());
    let mut f = Fixture {
        ctx: Default::default(),
        task,
        mol: Molecule::from_xyz("1\noriginal\nHe 0 0 0\n"),
        traj: Default::default(),
        actions: Default::default(),
        copied: None,
    };
    f.frame(vec![]);
    assert!(is_open(&f.ctx, &f.task));
    assert!(!f.actions.loaded);
    assert!(state(&f.ctx).plot_open, "the plot must open automatically");
    assert!(matches!(state(&f.ctx).view.unit, EnergyUnit::Kcal));
    let mut closed_plot = state(&f.ctx);
    closed_plot.plot_open = false;
    store(&f.ctx, closed_plot);
    f.click("Show Plot");
    assert!(state(&f.ctx).plot_open);
    f.click("Reverse order");
    assert!(state(&f.ctx).view.reversed);
    let displayed = f.frame(vec![]);
    let row_y = |label: &str| {
        displayed
            .shapes
            .iter()
            .find_map(|s| match &s.shape {
                egui::Shape::Text(t) if t.galley.text() == label => Some(t.pos.y),
                _ => None,
            })
            .unwrap_or_else(|| panic!("Missing numerical table value {label}"))
    };
    assert!(row_y("1.20000") < row_y("1.10000"));
    assert!(row_y("1.10000") < row_y("1.00000"));
    f.click("kJ/mol");
    f.click("Last point");
    f.click("Copy data");
    let text = f
        .copied
        .as_ref()
        .expect("Copy data must reach clipboard output");
    assert!(text.contains("relative_energy_kJ_mol"));
    let rows: Vec<_> = text.lines().skip(1).collect();
    assert!(rows[0].starts_with("3,"));
    assert!(rows[2].starts_with("1,"));
    assert_eq!(rows[2].split(',').nth(4), Some("0.00000000"));
    let s = state(&f.ctx);
    let xyz = s.report.as_ref().unwrap().xyz_ordered(s.view);
    let first = xyz.find("cycle 3 ").unwrap();
    let last = xyz.find("cycle 1 ").unwrap();
    assert!(
        first < last,
        "saved XYZ must follow the displayed table order"
    );

    f.click("Show selected");
    assert!(f.actions.loaded);
    assert_eq!(f.mol.atoms, vec!["H", "H"]);
    f.click("Load as Trajectory");
    assert!(f.actions.trajectory);
    assert_eq!(f.traj.frames.len(), 3);
    assert!((f.traj.frames[0].pos[1].x - 1.2).abs() < 1e-6);
    assert!((f.traj.frames[2].pos[1].x - 1.0).abs() < 1e-6);
    assert!(!f.traj.playing);
    assert!(matches!(f.traj.mode, crate::trajectory::PlaybackMode::Once));
    let mut s = state(&f.ctx);
    s.open = false;
    store(&f.ctx, s);
    f.frame(vec![]);
    assert!(!is_open(&f.ctx, &f.task));
    f.click("Open scan results");
    assert!(is_open(&f.ctx, &f.task));
}
