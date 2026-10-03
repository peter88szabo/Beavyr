use super::*;

struct Fixture {
    ctx: egui::Context,
    panel: QcPanelState,
    task: QcRunTask,
    opt: XtbPanelState,
    mol: Molecule,
}

impl Fixture {
    fn new() -> Self {
        let ctx = egui::Context::default();
        ctx.style_mut_of(ctx.theme(), |style| {
            style.text_styles = [
                (egui::TextStyle::Heading, egui::FontId::proportional(22.0)),
                (egui::TextStyle::Body, egui::FontId::proportional(16.0)),
                (egui::TextStyle::Monospace, egui::FontId::monospace(15.0)),
                (egui::TextStyle::Button, egui::FontId::proportional(16.0)),
                (egui::TextStyle::Small, egui::FontId::proportional(13.0)),
            ]
            .into();
        });
        store(
            &ctx,
            Presentation {
                modern: true,
                ..Default::default()
            },
        );
        Self {
            ctx,
            panel: QcPanelState {
                open: true,
                ..Default::default()
            },
            task: Default::default(),
            opt: XtbPanelState {
                paths: Vec::new(),
                ..Default::default()
            },
            mol: Molecule::from_xyz("1\ncarbon\nC 0 0 0\n"),
        }
    }

    fn frame(&mut self, events: Vec<egui::Event>, height: f32) -> egui::FullOutput {
        let mut output = self.ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1280.0, height),
                )),
                events,
                ..Default::default()
            },
            |ui| {
                let ctx = ui.ctx().clone();
                if self.panel.open {
                    qc_panel(
                        &ctx,
                        ui,
                        &mut self.panel,
                        &mut self.task,
                        &mut self.opt,
                        &self.mol,
                    );
                }
                input_editor_window(&ctx, &mut self.panel, &mut self.task, &self.opt, &self.mol);
            },
        );
        output.textures_delta.clear();
        output
    }

    fn click(&mut self, label: &str) {
        self.frame(vec![], 720.0);
        let output = self.frame(vec![], 720.0);
        let position = label_rect(&output, label)
            .unwrap_or_else(|| {
                let state = presentation(&self.ctx);
                let labels: Vec<_> = output
                    .shapes
                    .iter()
                    .filter_map(|shape| match &shape.shape {
                        egui::Shape::Text(text) => Some(text.galley.text()),
                        _ => None,
                    })
                    .collect();
                panic!(
                    "Missing control {label:?}; search={:?}; labels={labels:?}",
                    state.functional_search
                )
            })
            .center();
        for pressed in [true, false] {
            self.frame(
                vec![
                    egui::Event::PointerMoved(position),
                    egui::Event::PointerButton {
                        pos: position,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::NONE,
                    },
                ],
                720.0,
            );
        }
    }

    fn search(&mut self, query: &str) {
        self.click("Search by name or keyword");
        self.frame(vec![egui::Event::Text(query.to_owned())], 720.0);
        self.frame(vec![], 720.0);
    }
}

fn label_rect(output: &egui::FullOutput, label: &str) -> Option<egui::Rect> {
    fn find(shape: &egui::Shape, label: &str) -> Option<egui::Rect> {
        match shape {
            egui::Shape::Text(text) if text.galley.text() == label => {
                Some(egui::Rect::from_min_size(text.pos, text.galley.size()))
            }
            egui::Shape::Vec(shapes) => shapes.iter().find_map(|shape| find(shape, label)),
            _ => None,
        }
    }
    output
        .shapes
        .iter()
        .find_map(|shape| find(&shape.shape, label))
}

#[test]
fn prominent_code_chooser_uses_existing_program_defaults_without_changing_geometry() {
    let mut fixture = Fixture::new();
    fixture.panel.job = JobType::TdDft;
    let geometry = fixture.mol.pos.clone();
    let optimizer_program = fixture.opt.program;
    fixture.click("ORCA");
    fixture.click("DREIDING");
    assert_eq!(fixture.panel.program, QcProgram::Dreiding);
    assert_eq!(fixture.panel.method, None);
    assert_eq!(fixture.panel.job, JobType::Optimize);
    assert_eq!(fixture.opt.program, optimizer_program);
    assert_eq!(fixture.mol.pos, geometry);
    assert!(!fixture.task.is_running());
}

#[test]
fn switching_classic_and_docking_preserves_the_shared_setup_and_last_run() {
    let mut fixture = Fixture::new();
    fixture.panel.config.basis = "def2-TZVP".into();
    fixture.panel.config.nproc = 3;
    fixture.panel.job = JobType::Gradient;
    fixture.panel.advanced = true;
    fixture.panel.advanced_text = "! HF def2-TZVP\n* xyz 0 1\nC 0 0 0\n*\n".into();
    fixture.task.last_message = Some("Previous calculation finished".into());
    let input = fixture.panel.advanced_text.clone();
    assert!(!is_docked(&fixture.ctx));
    fixture.click("Dock");
    assert!(is_docked(&fixture.ctx));
    fixture.click("Classic");
    assert!(!is_modern(&fixture.ctx));
    fixture.frame(vec![], 720.0);
    let classic = fixture.frame(vec![], 720.0);
    assert!(
        label_rect(&classic, "Program").is_some(),
        "Classic window open={}; painted labels={:?}",
        fixture.panel.open,
        classic
            .shapes
            .iter()
            .filter_map(|shape| match &shape.shape {
                egui::Shape::Text(text) => Some(text.galley.text()),
                _ => None,
            })
            .collect::<Vec<_>>()
    );
    fixture.click("Modern");
    assert!(is_modern(&fixture.ctx));
    assert!(is_docked(&fixture.ctx));
    fixture.click("Float");
    assert!(!is_docked(&fixture.ctx));
    assert_eq!(fixture.panel.config.basis, "def2-TZVP");
    assert_eq!(fixture.panel.config.nproc, 3);
    assert_eq!(fixture.panel.job, JobType::Gradient);
    assert_eq!(fixture.panel.advanced_text, input);
    assert!(fixture.panel.advanced);
    assert_eq!(
        fixture.task.last_message.as_deref(),
        Some("Previous calculation finished")
    );
}

#[test]
fn unsupported_jobs_remain_disabled_in_the_modern_picker() {
    let mut fixture = Fixture::new();
    fixture.panel.set_program(QcProgram::Dreiding);
    fixture.click("More…");
    fixture.click("TD-DFT");
    assert_eq!(fixture.panel.job, JobType::Optimize);
    assert!(!fixture.task.is_running());
}

#[test]
fn searching_for_a_composite_functional_uses_the_catalogue_keyword_and_hides_basis() {
    let mut fixture = Fixture::new();
    fixture.click("HF");
    fixture.click("DFT");
    fixture.click("B3LYP");
    fixture.search("r2scan-3c");
    fixture.click("r2SCAN-3c");
    assert_eq!(fixture.panel.config.functional, "r2SCAN-3c");
    assert!(!fixture.panel.needs_basis());
    let text = fixture
        .panel
        .generated_input(&fixture.mol.atoms, &[0.0, 0.0, 0.0]);
    assert!(text.contains("r2SCAN-3c"));
    assert!(!text.contains("def2-SVP"));
}

#[test]
fn code_chooser_and_action_stay_visible_with_long_input_at_laptop_height() {
    for docked in [false, true] {
        let mut fixture = Fixture::new();
        fixture.panel.advanced = true;
        fixture.panel.advanced_text = "! HF def2-SVP\n".repeat(40);
        let mut state = presentation(&fixture.ctx);
        state.docked = docked;
        state.configure = true;
        store(&fixture.ctx, state);
        for height in [720.0, 600.0] {
            fixture.frame(vec![], height);
            let output = fixture.frame(vec![], height);
            for label in [
                "Quantum chemistry code",
                "ORCA",
                "Configure…",
                "CPU cores",
                "Memory (MB)",
                "Start optimization",
            ] {
                let rect = label_rect(&output, label)
                    .unwrap_or_else(|| panic!("Missing {label} at {height}"));
                assert!(
                    rect.top() >= 0.0 && rect.bottom() <= height,
                    "{label} is offscreen: {rect:?}"
                );
                assert!(
                    rect.left() >= 0.0 && rect.right() <= 1280.0,
                    "{label} overflows: {rect:?}"
                );
            }
        }
    }
}

#[test]
fn invalid_electronic_state_disables_launch_before_any_calculation_starts() {
    let mut fixture = Fixture::new();
    fixture.panel.set_program(QcProgram::Dreiding);
    fixture.panel.multiplicity = 2; // Carbon has an even electron count.
    let code = code_status(
        &mut presentation(&fixture.ctx),
        &fixture.panel,
        &fixture.opt,
    );
    assert!(blocking_reason(&fixture.panel, &fixture.mol, &code).is_some());
    fixture.click("Start optimization");
    assert!(!fixture.task.is_running());
    assert!(
        !fixture.panel.show_warnings,
        "the disabled action must not enter the launch path"
    );
}

#[test]
fn unavailable_code_is_blocked_and_custom_input_is_ignored_for_non_file_backends() {
    let mut fixture = Fixture::new();
    let code = code_status(
        &mut presentation(&fixture.ctx),
        &fixture.panel,
        &fixture.opt,
    );
    assert!(!code.ready);
    assert!(blocking_reason(&fixture.panel, &fixture.mol, &code)
        .unwrap()
        .contains("executable path"));
    fixture.panel.advanced = true;
    fixture.panel.advanced_text.clear();
    fixture.panel.set_program(QcProgram::Dreiding);
    let code = code_status(
        &mut presentation(&fixture.ctx),
        &fixture.panel,
        &fixture.opt,
    );
    assert!(code.ready);
    assert!(blocking_reason(&fixture.panel, &fixture.mol, &code).is_none());
}

#[test]
fn resources_remain_within_the_four_core_limit_in_both_layouts() {
    let mut fixture = Fixture::new();
    fixture.panel.config.nproc = 64;
    fixture.frame(vec![], 720.0);
    assert_eq!(fixture.panel.config.nproc, 4);
    let mut state = presentation(&fixture.ctx);
    state.modern = false;
    store(&fixture.ctx, state);
    fixture.panel.config.nproc = 64;
    fixture.frame(vec![], 720.0);
    assert_eq!(fixture.panel.config.nproc, 4);
}

#[test]
fn floating_input_editor_keeps_the_draft_after_close_and_syncs_current_settings() {
    let mut fixture = Fixture::new();
    fixture.click("Edit input…");
    assert!(fixture.panel.advanced);
    fixture.frame(vec![egui::Event::Text("# manual change\n".into())], 720.0);
    assert!(fixture.panel.advanced_text.contains("# manual change"));
    let draft = fixture.panel.advanced_text.clone();
    fixture.click("Close editor");
    fixture.click("Edit input…");
    assert_eq!(fixture.panel.advanced_text, draft);
    fixture.panel.config.basis = "def2-TZVP".into();
    fixture.panel.config.nproc = 3;
    fixture.mol.pos[0] = Vec3::new(1.0, 2.0, 3.0);
    fixture.click("Sync from settings");
    assert_eq!(
        fixture.panel.advanced_text,
        fixture
            .panel
            .generated_input(&fixture.mol.atoms, &[1.0, 2.0, 3.0])
    );
    assert!(!fixture.panel.advanced_text.contains("# manual change"));
    assert!(fixture.panel.advanced);
}

#[test]
fn input_editor_works_with_setup_closed_and_does_not_launch_an_unavailable_code() {
    let mut fixture = Fixture::new();
    fixture.click("Edit input…");
    fixture.panel.open = false;
    fixture.frame(vec![], 600.0);
    let output = fixture.frame(vec![], 600.0);
    for label in [
        "Save input…",
        "Sync from settings",
        "Run input",
        "Close editor",
    ] {
        let rect = label_rect(&output, label).unwrap_or_else(|| panic!("Missing {label}"));
        assert!(rect.top() >= 0.0 && rect.bottom() <= 600.0);
        assert!(rect.left() >= 0.0 && rect.right() <= 1280.0);
    }
    fixture.click("Run input");
    assert!(!fixture.task.is_running());
    assert!(!fixture.panel.show_warnings);
    fixture.click("Close editor");
    assert!(
        fixture.panel.advanced,
        "closing keeps the manual input active"
    );
}

#[test]
fn changing_code_requires_explicit_input_sync_without_losing_the_old_draft() {
    let mut fixture = Fixture::new();
    fixture.click("Edit input…");
    let draft = fixture.panel.advanced_text.clone();
    fixture.panel.set_program(QcProgram::Psi4);
    assert!(input_editor::program_mismatch(&fixture.ctx, &fixture.panel).is_some());
    assert_eq!(fixture.panel.advanced_text, draft);
    fixture.click("Sync from settings");
    assert!(input_editor::program_mismatch(&fixture.ctx, &fixture.panel).is_none());
    assert_eq!(
        fixture.panel.advanced_text,
        fixture
            .panel
            .generated_input(&fixture.mol.atoms, &[0.0, 0.0, 0.0])
    );
}
