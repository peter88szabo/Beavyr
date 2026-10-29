use super::*;
use bevy::ecs::system::SystemState;

struct SurfaceWindow {
    ctx: egui::Context,
    state: OrbitalState,
    mol: Molecule,
    settings: MolSettings,
    world: World,
    messages: SystemState<MessageWriter<'static, MoleculeChanged>>,
    open: bool,
}

impl SurfaceWindow {
    fn new() -> Self {
        let mut world = World::new();
        world.init_resource::<Messages<MoleculeChanged>>();
        let messages = SystemState::new(&mut world);
        let mut state = OrbitalState::default();
        let text = include_str!("../../tests/fixtures/behemoth_ch2o_canonicalMO.molden");
        state.adopt(parse_molden(text).unwrap(), None, None, Some(text.into()));
        Self {
            ctx: egui::Context::default(),
            state,
            mol: Molecule::empty(),
            settings: MolSettings::default(),
            world,
            messages,
            open: true,
        }
    }

    fn frame(&mut self, events: Vec<egui::Event>) -> egui::Rect {
        let mut rects = Vec::new();
        let mut writer = self.messages.get_mut(&mut self.world).unwrap();
        let mut output = self.ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1400.0, 1200.0),
                )),
                events,
                ..Default::default()
            },
            |ui| {
                crate::ui_layout::editor_section(
                    ui,
                    true,
                    &mut self.open,
                    &mut rects,
                    crate::ui_layout::Tab::Surface.default_size(),
                    "Surface Tools",
                    |ui| {
                        orbital_panel(
                            ui,
                            true,
                            &mut self.state,
                            &mut self.mol,
                            &mut self.settings,
                            &mut writer,
                        );
                    },
                );
            },
        );
        output.textures_delta.clear();
        rects[0]
    }

    fn drag(&mut self, start: egui::Pos2, end: egui::Pos2) -> egui::Rect {
        self.frame(vec![egui::Event::PointerMoved(start)]);
        self.frame(vec![egui::Event::PointerButton {
            pos: start,
            button: egui::PointerButton::Primary,
            pressed: true,
            modifiers: egui::Modifiers::NONE,
        }]);
        self.frame(vec![egui::Event::PointerMoved(end)]);
        self.frame(vec![egui::Event::PointerButton {
            pos: end,
            button: egui::PointerButton::Primary,
            pressed: false,
            modifiers: egui::Modifiers::NONE,
        }]);
        self.frame(vec![]);
        self.frame(vec![])
    }
}

#[test]
fn surface_window_keeps_a_user_selected_height_after_vertical_resizing() {
    let mut window = SurfaceWindow::new();
    window.frame(vec![]);
    window.frame(vec![]);
    let before = window.frame(vec![]);
    let corner = before.right_bottom() - egui::vec2(1.0, 1.0);
    let taller = window.drag(corner, corner + egui::vec2(0.0, 140.0));
    assert!(
        taller.height() > before.height() + 100.0,
        "vertical resize did not persist: {before:?} -> {taller:?}"
    );
    assert!((taller.width() - before.width()).abs() < 3.0);
    let corner = taller.right_bottom() - egui::vec2(1.0, 1.0);
    let shorter = window.drag(corner, corner - egui::vec2(0.0, 80.0));
    assert!(
        shorter.height() < taller.height() - 50.0,
        "window cannot shrink vertically: {taller:?} -> {shorter:?}"
    );
    assert!((shorter.width() - taller.width()).abs() < 3.0);
}
