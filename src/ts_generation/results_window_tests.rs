use super::*;
use crate::optimizer::rda::RdaGuessResult;
use crate::ts_generation::GenerationOutput;

fn output(algorithm: Algorithm) -> GenerationOutput {
    let atoms = vec!["H".to_owned(), "H".to_owned()];
    GenerationOutput {
        algorithm,
        method: "test energy model".into(),
        guess: RdaGuessResult {
            atoms: atoms.clone(),
            q: vec![0., 0., 0., 3., 0., 0.],
            reaction_direction: None,
            points: vec![],
            direction: String::new(),
            bracketed: true,
            message: "Test TS guess".into(),
            rda_trajectory_file: None,
            rda_chain_file: None,
            poormans_neb_file: None,
            poormans_neb_profile_file: None,
        },
        energy_hartree: -1.0,
        guess_xyz: "2\nguess\nH 0 0 0\nH 1.5 0 0\n".into(),
        path_xyz: String::new(),
        search_xyz: String::new(),
        profile_dat: None,
        frames: [1.0, 2.0]
            .into_iter()
            .map(|x| TrajectoryFrame {
                atoms: atoms.clone(),
                pos: vec![Vec3::ZERO, Vec3::X * x],
            })
            .collect(),
    }
}

#[test]
fn both_completed_methods_open_once_and_a_new_run_reopens() {
    for algorithm in [Algorithm::Rda, Algorithm::PoorMansNeb] {
        let mut generation = TsGeneration::default();
        let mut window = WindowState::default();
        window.sync(&generation);
        assert!(!window.open);
        generation.output = Some(output(algorithm));
        window.sync(&generation);
        assert!(window.open);
        window.open = false;
        window.sync(&generation);
        assert!(!window.open, "closing must persist across frames");
        // A fast new run may finish between UI frames, without showing an empty result.
        generation.control = Arc::new(RunControl::default());
        generation.output = Some(output(algorithm));
        window.sync(&generation);
        assert!(window.open);
        generation.output = None;
        generation.is_error = true;
        window.sync(&generation);
        assert!(
            !window.open,
            "a failed or cancelled run has no result window"
        );
    }
}

struct Fixture {
    ctx: egui::Context,
    generation: TsGeneration,
    mol: Molecule,
    traj: TrajectoryState,
    displayed: bool,
}
impl Fixture {
    fn new(algorithm: Algorithm) -> Self {
        let mut generation = TsGeneration::default();
        generation.output = Some(output(algorithm));
        Self {
            ctx: Default::default(),
            generation,
            mol: Molecule::from_xyz("1\noriginal\nHe 0 0 0\n"),
            traj: Default::default(),
            displayed: false,
        }
    }
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
                // The setup panel is deliberately absent.
                let (_, displayed) = window(
                    ui.ctx(),
                    &mut self.generation,
                    &mut self.mol,
                    &mut self.traj,
                );
                self.displayed |= displayed;
                egui::Window::new("Reopen test")
                    .default_pos(egui::pos2(0., 0.))
                    .show(ui.ctx(), |ui| {
                        reopen_button(ui, &self.generation);
                    });
            },
        );
        output.textures_delta.clear();
        output
    }
    fn click(&mut self, label: &str) {
        self.frame(vec![]);
        let output = self.frame(vec![]);
        let position = output
            .shapes
            .iter()
            .find_map(|shape| match &shape.shape {
                egui::Shape::Text(text) if text.galley.text() == label => {
                    Some(text.pos + text.galley.size() / 2.0)
                }
                _ => None,
            })
            .unwrap_or_else(|| panic!("Missing button {label}"));
        for pressed in [true, false] {
            self.frame(vec![
                egui::Event::PointerMoved(position),
                egui::Event::PointerButton {
                    pos: position,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ]);
        }
    }
}

#[test]
fn independent_window_keeps_guess_and_trajectory_actions_and_can_be_reopened() {
    for algorithm in [Algorithm::Rda, Algorithm::PoorMansNeb] {
        let mut fixture = Fixture::new(algorithm);
        fixture.frame(vec![]);
        assert!(is_open(&fixture.ctx, &fixture.generation));
        assert!(
            !fixture.displayed,
            "opening results must not change the molecule"
        );
        fixture.click("Show TS guess");
        assert!(fixture.displayed);
        assert_eq!(fixture.mol.atoms, ["H", "H"]);
        assert!((fixture.mol.pos[1].x - 3.0 * BOHR_TO_ANGSTROM as f32).abs() < 1e-6);
        assert_eq!(fixture.traj.frames.len(), 1);
        fixture.displayed = false;
        fixture.click("Show path in Trajectory");
        assert!(fixture.displayed);
        assert_eq!(fixture.traj.frames.len(), 2);
        assert!(!fixture.traj.playing);
        assert!(matches!(
            fixture.traj.mode,
            crate::trajectory::PlaybackMode::Once
        ));
        let mut closed = state(&fixture.ctx);
        closed.open = false;
        store(&fixture.ctx, closed);
        fixture.frame(vec![]);
        assert!(!is_open(&fixture.ctx, &fixture.generation));
        fixture.click("Open TS results");
        assert!(is_open(&fixture.ctx, &fixture.generation));
    }
}
