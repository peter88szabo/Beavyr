use super::*;

fn frame(ctx: &egui::Context, state: &mut State, events: Vec<egui::Event>) -> egui::FullOutput {
    let mut output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1280., 800.),
            )),
            events,
            ..Default::default()
        },
        |ui| {
            window(ui.ctx(), state, false);
        },
    );
    // Headless frames do not upload textures to a GPU renderer.
    output.textures_delta.clear();
    output
}

fn label(output: &egui::FullOutput, text: &str) -> egui::Rect {
    fn find(shape: &egui::Shape, text: &str) -> Option<egui::Rect> {
        match shape {
            egui::Shape::Text(s) if s.galley.text() == text => {
                Some(egui::Rect::from_min_size(s.pos, s.galley.size()))
            }
            egui::Shape::Vec(shapes) => shapes.iter().find_map(|s| find(s, text)),
            _ => None,
        }
    }
    output
        .shapes
        .iter()
        .find_map(|s| find(&s.shape, text))
        .unwrap_or_else(|| panic!("Missing {text}"))
}

#[test]
fn output_window_opens_independently_and_search_results_can_be_copied() {
    let ctx = egui::Context::default();
    let mut state = State {
        open: true,
        search: "energy".into(),
        matching_only: true,
        snapshot: Snapshot {
            text: "iteration 1\nEnergy = -1.0\nconverged\n".into(),
            ..Default::default()
        },
        ..Default::default()
    };
    frame(&ctx, &mut state, vec![]);
    let output = frame(&ctx, &mut state, vec![]);
    label(&output, "Live output");
    label(&output, "1 matching lines in the displayed output");
    label(&output, "Energy = -1.0");
    let position = label(&output, "Copy displayed output").center();
    let mut copied = None;
    for pressed in [true, false] {
        let output = frame(
            &ctx,
            &mut state,
            vec![
                egui::Event::PointerMoved(position),
                egui::Event::PointerButton {
                    pos: position,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ],
        );
        for command in output.platform_output.commands {
            if let egui::OutputCommand::CopyText(text) = command {
                copied = Some(text);
            }
        }
    }
    assert_eq!(copied.as_deref(), Some("Energy = -1.0"));
    state.open = false;
    let mut output = ctx.run_ui(Default::default(), |ui| {
        assert!(window(ui.ctx(), &mut state, false).is_none())
    });
    output.textures_delta.clear();
}
