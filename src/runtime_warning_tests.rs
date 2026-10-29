use super::*;

fn frame(ctx: &egui::Context, events: Vec<egui::Event>) -> (egui::FullOutput, Option<egui::Rect>) {
    let mut rect = None;
    let mut output = ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1000.0, 700.0),
            )),
            events,
            ..Default::default()
        },
        |ui| {
            rect = window(ui.ctx());
        },
    );
    output.textures_delta.clear();
    (output, rect)
}

#[test]
fn missing_python_environment_opens_a_centered_actionable_popup_once() {
    let ctx = egui::Context::default();
    let environment = crate::qchem_interfaces::python_env::PythonEnv::ModuleMissing {
        module: "pyscf",
        executable: "/usr/bin/python3".into(),
        prefix: "/usr".into(),
    };
    crate::qchem_interfaces::python_env::notify_missing(&ctx, &environment);
    crate::qchem_interfaces::python_env::notify_missing(&ctx, &environment);
    frame(&ctx, vec![]);
    let (output, rect) = frame(&ctx, vec![]);
    assert!((rect.unwrap().center() - egui::pos2(500.0, 350.0)).length() < 2.0);
    let texts: Vec<_> = output
        .shapes
        .iter()
        .filter_map(|s| match &s.shape {
            egui::Shape::Text(t) => Some(t.galley.text()),
            _ => None,
        })
        .collect();
    for expected in [
        "Python environment setup required",
        "Conda",
        "Environment activation",
        "/usr/bin/python3",
        "pyscf",
    ] {
        assert!(
            texts.iter().any(|t| t.contains(expected)),
            "Missing {expected}: {texts:?}"
        );
    }
    let ok = output
        .shapes
        .iter()
        .find_map(|s| match &s.shape {
            egui::Shape::Text(t) if t.galley.text() == "OK" => Some(t.pos + t.galley.size() / 2.0),
            _ => None,
        })
        .unwrap();
    for pressed in [true, false] {
        frame(
            &ctx,
            vec![
                egui::Event::PointerMoved(ok),
                egui::Event::PointerButton {
                    pos: ok,
                    button: egui::PointerButton::Primary,
                    pressed,
                    modifiers: egui::Modifiers::NONE,
                },
            ],
        );
    }
    crate::qchem_interfaces::python_env::notify_missing(&ctx, &environment);
    assert!(
        frame(&ctx, vec![]).1.is_none(),
        "Acknowledged warnings must not reopen each frame"
    );
}

#[test]
fn ready_python_environments_do_not_open_a_warning() {
    let ctx = egui::Context::default();
    let environment = crate::qchem_interfaces::python_env::PythonEnv::Ready {
        module: "pyscf",
        version: "test".into(),
        executable: "/env/bin/python3".into(),
        prefix: "/env".into(),
    };
    crate::qchem_interfaces::python_env::notify_missing(&ctx, &environment);
    assert!(frame(&ctx, vec![]).1.is_none());
}

#[test]
fn no_python_interpreter_message_explains_saved_environment_settings() {
    let environment = crate::qchem_interfaces::python_env::PythonEnv::NoInterpreter {
        detail: "not found".into(),
    };
    let text = environment.detail();
    assert!(
        text.contains("Conda") && text.contains("micromamba") && text.contains("Save & Check"),
        "{text}"
    );
}

#[test]
fn copying_an_error_copies_all_diagnostics_and_keeps_the_popup_open() {
    let ctx = egui::Context::default();
    let message = "PySCF optimization failed.\n\nTraceback:\nModuleNotFoundError: No module named 'geometric'\nOutput kept in /tmp/example.";
    observe(&ctx, "calculator", "Calculation failed", Some(message));
    frame(&ctx, vec![]);
    let (output, _) = frame(&ctx, vec![]);
    let position = output
        .shapes
        .iter()
        .find_map(|shape| match &shape.shape {
            egui::Shape::Text(text) if text.galley.text() == "Copy error" => {
                Some(text.pos + text.galley.size() / 2.0)
            }
            _ => None,
        })
        .expect("The error window needs a visible Copy error button");
    let mut copied = false;
    for pressed in [true, false] {
        let (output, rect) = frame(
            &ctx,
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
        copied |= output.platform_output.commands.iter().any(|command| {
            matches!(command,
            egui::OutputCommand::CopyText(text) if text == message)
        });
        assert!(rect.is_some());
    }
    assert!(copied, "The full error text must be sent to the clipboard");
}

#[test]
fn stored_errors_do_not_reopen_every_frame_and_can_report_a_later_failure() {
    let ctx = egui::Context::default();
    observe(&ctx, "calculator", "Calculation failed", Some("Failed"));
    // Acknowledge the queue directly; GUI dismissal is covered above.
    ctx.data_mut(|data| {
        data.get_temp_mut_or_default::<Notices>(egui::Id::new("runtime_setup_notices"))
            .0
            .clear()
    });
    observe(&ctx, "calculator", "Calculation failed", Some("Failed"));
    assert!(frame(&ctx, vec![]).1.is_none());
    observe(&ctx, "calculator", "Calculation failed", None);
    observe(&ctx, "calculator", "Calculation failed", Some("Failed"));
    assert!(frame(&ctx, vec![]).1.is_some());
}

#[test]
fn successful_and_cancelled_calculations_do_not_raise_an_error_popup() {
    let ctx = egui::Context::default();
    observe(&ctx, "calculator", "Calculation failed", None);
    observe(&ctx, "optimizer", "Optimization failed", Some("Cancelled."));
    assert!(frame(&ctx, vec![]).1.is_none());
}
