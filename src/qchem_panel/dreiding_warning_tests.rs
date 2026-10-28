use super::*;

fn frame(
    ctx: &egui::Context,
    task: &mut QcRunTask,
    events: Vec<egui::Event>,
    size: egui::Vec2,
) -> (egui::FullOutput, Option<egui::Rect>) {
    let mut rect = None;
    let mut output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, size)),
            events,
            ..Default::default()
        },
        |ui| {
            rect = dreiding_warning_window(ui.ctx(), task);
        },
    );
    output.textures_delta.clear();
    (output, rect)
}

#[test]
fn unsupported_system_popup_is_centered_and_acknowledged() {
    for size in [egui::vec2(1280.0, 720.0), egui::vec2(800.0, 600.0)] {
        let ctx = egui::Context::default();
        let mut task = QcRunTask::default();
        let cause = crate::forcefield::dreiding::BuildError::Typing(
            "valence electrons not defined for element Fe".into(),
        );
        task.dreiding_warning = Some(crate::forcefield::dreiding::diagnostics::message(&cause));
        frame(&ctx, &mut task, vec![], size);
        let (output, rect) = frame(&ctx, &mut task, vec![], size);
        let rect = rect.expect("the pending failure must open a popup");
        assert!(
            (rect.center() - egui::pos2(size.x / 2.0, size.y / 2.0)).length() < 2.0,
            "{rect:?}"
        );
        let texts: Vec<_> = output
            .shapes
            .iter()
            .filter_map(|s| match &s.shape {
                egui::Shape::Text(t) => Some(t.galley.text()),
                _ => None,
            })
            .collect();
        for needed in [
            "DREIDING unavailable",
            "Quantum Chemistry",
            "xTB",
            "total charge",
            "spin multiplicity",
        ] {
            assert!(
                texts.iter().any(|t| t.contains(needed)),
                "Missing {needed}: {texts:?}"
            );
        }
        let ok = output
            .shapes
            .iter()
            .find_map(|s| match &s.shape {
                egui::Shape::Text(t) if t.galley.text() == "OK" => {
                    Some(t.pos + t.galley.size() / 2.0)
                }
                _ => None,
            })
            .unwrap();
        for pressed in [true, false] {
            frame(
                &ctx,
                &mut task,
                vec![
                    egui::Event::PointerMoved(ok),
                    egui::Event::PointerButton {
                        pos: ok,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::NONE,
                    },
                ],
                size,
            );
        }
        assert!(task.dreiding_warning.is_none());
        assert!(frame(&ctx, &mut task, vec![], size).1.is_none());
    }
}

#[test]
fn ordinary_errors_and_success_do_not_open_a_dreiding_popup() {
    let ctx = egui::Context::default();
    let mut task = QcRunTask::default();
    task.last_is_error = true;
    task.last_message = Some("Cancelled.".into());
    assert!(frame(&ctx, &mut task, vec![], egui::vec2(800.0, 600.0))
        .1
        .is_none());
}
