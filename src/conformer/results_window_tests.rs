use super::*;

struct Fixture {
    ctx: egui::Context,
    run: ConformerRun,
    mol: Molecule,
    traj: TrajectoryState,
    settings: MolSettings,
    open_trajectory: bool,
}

impl Fixture {
    fn new() -> Self {
        Self {
            ctx: egui::Context::default(),
            run: ConformerRun {
                outcome: Some(super::super::tests::outcome()),
                started: Some(Instant::now()),
                ..Default::default()
            },
            mol: Molecule::empty(),
            traj: TrajectoryState::default(),
            settings: MolSettings::default(),
            open_trajectory: false,
        }
    }

    fn frame(
        &mut self,
        events: Vec<egui::Event>,
        reopen: bool,
    ) -> (egui::FullOutput, Option<egui::Rect>, bool) {
        let mut result = (None, ResultActions::default());
        let mut output = self.ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1200.0, 900.0),
                )),
                events,
                ..Default::default()
            },
            |ui| {
                if reopen {
                    reopen_button(ui, &self.run);
                }
                result = window(
                    ui.ctx(),
                    &mut self.run,
                    &mut self.mol,
                    &mut self.traj,
                    &self.settings,
                );
            },
        );
        output.textures_delta.clear();
        self.open_trajectory = result.1.open_trajectory;
        (output, result.0, result.1.loaded)
    }

    fn close(&mut self) {
        let mut value = state(&self.ctx);
        value.open = false;
        store(&self.ctx, value);
    }

    fn click(&mut self, position: egui::Pos2, reopen: bool) -> bool {
        let mut loaded = false;
        for pressed in [true, false] {
            loaded |= self
                .frame(
                    vec![
                        egui::Event::PointerMoved(position),
                        egui::Event::PointerButton {
                            pos: position,
                            button: egui::PointerButton::Primary,
                            pressed,
                            modifiers: egui::Modifiers::NONE,
                        },
                    ],
                    reopen,
                )
                .2;
        }
        loaded
    }
}

fn text_position(output: &egui::FullOutput, label: &str) -> egui::Pos2 {
    output
        .shapes
        .iter()
        .find_map(|shape| match &shape.shape {
            egui::Shape::Text(text) if text.galley.text() == label => {
                Some(text.pos + text.galley.size() / 2.0)
            }
            _ => None,
        })
        .unwrap_or_else(|| panic!("Missing {label}"))
}

#[test]
fn completion_opens_results_without_drawing_the_setup_panel() {
    let mut fixture = Fixture::new();
    fixture.frame(vec![], false);
    let (output, rect, _) = fixture.frame(vec![], false);
    assert!(rect.is_some());
    text_position(&output, "Conformer Search Results");
    text_position(&output, "Population");
}

#[test]
fn dismissed_results_stay_closed_and_the_reopen_button_works() {
    let mut fixture = Fixture::new();
    fixture.frame(vec![], false);
    fixture.close();
    assert!(fixture.frame(vec![], false).1.is_none());
    let (output, rect, _) = fixture.frame(vec![], true);
    assert!(rect.is_none());
    let position = text_position(&output, "Open conformer results");
    fixture.click(position, true);
    assert!(fixture.frame(vec![], false).1.is_some());
}

#[test]
fn a_new_finished_search_opens_again_but_refinement_keeps_a_closed_window_closed() {
    let mut fixture = Fixture::new();
    fixture.frame(vec![], false);
    fixture.close();
    fixture.run.outcome.as_mut().unwrap().refined_with_xtb = true;
    assert!(fixture.frame(vec![], false).1.is_none());
    fixture.run.started = fixture
        .run
        .started
        .map(|start| start + std::time::Duration::from_secs(1));
    assert!(fixture.frame(vec![], false).1.is_some());
}

#[test]
fn viewing_a_result_restores_its_atoms_and_reports_a_geometry_change() {
    let mut fixture = Fixture::new();
    fixture.frame(vec![], false);
    let (output, _, _) = fixture.frame(vec![], false);
    let position = text_position(&output, "View");
    assert!(fixture.click(position, false));
    assert_eq!(fixture.mol.atoms, ["C", "H"]);
    assert_eq!(fixture.run.previewing, Some(0));
    assert!(!fixture.open_trajectory);
    assert!((fixture.mol.pos[1].x - 1.1).abs() < 1.0e-5);
}

#[test]
fn loading_results_as_a_trajectory_works_from_the_separate_window() {
    let mut fixture = Fixture::new();
    fixture.frame(vec![], false);
    let (output, _, _) = fixture.frame(vec![], false);
    let position = text_position(&output, "Load as trajectory");
    assert!(fixture.click(position, false));
    assert_eq!(fixture.traj.frames.len(), 2);
    assert!(fixture.open_trajectory);
    assert_eq!(fixture.mol.atoms, ["C", "H"]);
}
