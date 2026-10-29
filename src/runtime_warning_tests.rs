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
