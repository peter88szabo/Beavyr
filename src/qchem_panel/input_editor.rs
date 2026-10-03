//! A floating editor for the exact input used by the shared calculation runner.

use crate::ui_style::ResponseExt;

use super::*;

const STATE_ID: &str = "qc_input_editor_state_v1";

#[derive(Clone, Default)]
struct EditorState {
    open: bool,
    focus: bool,
    program: Option<QcProgram>,
    file_name: String,
    save_status: Option<(String, bool)>,
}

fn state(ctx: &egui::Context) -> EditorState {
    ctx.data_mut(|data| {
        data.get_temp_mut_or_default::<EditorState>(egui::Id::new(STATE_ID))
            .clone()
    })
}

fn store_state(ctx: &egui::Context, state: EditorState) {
    ctx.data_mut(|data| data.insert_temp(egui::Id::new(STATE_ID), state));
}

fn sync(state: &mut EditorState, panel: &mut QcPanelState, mol: &Molecule) {
    let positions: Vec<_> = mol
        .pos
        .iter()
        .flat_map(|position| [position.x as f64, position.y as f64, position.z as f64])
        .collect();
    panel.advanced_text = panel.generated_input(&mol.atoms, &positions);
    panel.advanced = true;
    state.program = Some(panel.program);
    state.file_name = panel.input_file_name().to_owned();
    state.save_status = None;
}

pub(super) fn entry(ui: &mut egui::Ui, panel: &mut QcPanelState, mol: &Molecule, running: bool) {
    if !panel.has_input_file() {
        return;
    }
    let mut state = state(ui.ctx());
    ui.horizontal_wrapped(|ui| {
        if ui
            .add_enabled(!running, egui::Button::new("Edit input…"))
            .on_hover_text("Open the input in a large, resizable editor")
            .clicked_once()
        {
            if panel.advanced_text.trim().is_empty() {
                sync(&mut state, panel, mol);
            } else if state.program.is_none() {
                state.program = Some(panel.program);
                state.file_name = panel.input_file_name().to_owned();
            }
            panel.advanced = true;
            state.open = true;
            state.focus = true;
        }
        if panel.advanced
            && ui
                .add_enabled(!running, egui::Button::new("Use settings"))
                .on_hover_text("Run input generated from the controls; keep the manual draft")
                .clicked_once()
        {
            panel.advanced = false;
            state.open = false;
        }
    });
    store_state(ui.ctx(), state);
}

pub(super) fn program_mismatch(ctx: &egui::Context, panel: &QcPanelState) -> Option<String> {
    if !panel.advanced || !panel.has_input_file() {
        return None;
    }
    let original = state(ctx).program?;
    (original != panel.program).then(|| {
        format!(
            "Input is for {}. Sync from settings to use {}.",
            original.label(),
            panel.program.label()
        )
    })
}

fn save(state: &mut EditorState, text: &str) {
    if let Some(target) = rfd::FileDialog::new()
        .set_file_name(&state.file_name)
        .save_file()
    {
        state.save_status = Some(match std::fs::write(&target, text) {
            Ok(()) => (format!("Saved {}", target.display()), false),
            Err(error) => (format!("Could not save input: {error}"), true),
        });
    }
}

pub(crate) fn show(
    ctx: &egui::Context,
    panel: &mut QcPanelState,
    task: &mut QcRunTask,
    opt_panel: &XtbPanelState,
    mol: &Molecule,
) -> Option<egui::Rect> {
    let mut state = state(ctx);
    if !state.open {
        return None;
    }
    let mut open = true;
    let mut close = false;
    let viewport = ctx.viewport_rect();
    let title = format!(
        "{} Input Editor",
        state.program.unwrap_or(panel.program).label()
    );
    let result = egui::Window::new(title)
        .id(egui::Id::new("qc_input_editor_window_v1"))
        .open(&mut open)
        .resizable(true)
        .default_size([880.0, (viewport.height() - 96.0).max(360.0)])
        .min_size([600.0, 360.0])
        .default_pos([
            (viewport.center().x - 440.0).max(viewport.left() + 16.0),
            viewport.top() + 32.0,
        ])
        .show(ctx, |ui| {
            style(ui);
            let running = task.is_running();
            let mut presentation = presentation(ctx);
            let code = code_status(&mut presentation, panel, opt_panel);
            let reason = if !panel.has_input_file() {
                Some("Choose a code with an editable input file to run this input".to_owned())
            } else {
                program_mismatch(ctx, panel).or_else(|| blocking_reason(panel, mol, &code))
            };
            ui.horizontal_wrapped(|ui| {
                if ui.button("Save input…").clicked_once() {
                    save(&mut state, &panel.advanced_text);
                }
                if ui
                    .add_enabled(
                        !running && panel.has_input_file(),
                        egui::Button::new("Sync from settings"),
                    )
                    .on_hover_text("Replace the draft with input from the current controls and molecule")
                    .clicked_once()
                {
                    sync(&mut state, panel, mol);
                }
                if running {
                    ui.spinner();
                    if ui.button("Cancel calculation").clicked_once() {
                        task.cancel();
                    }
                } else {
                    let run = ui.add_enabled(
                        reason.is_none(),
                        egui::Button::new(
                            egui::RichText::new("Run input").color(egui::Color32::WHITE),
                        )
                        .fill(egui::Color32::from_rgb(38, 96, 159)),
                    );
                    if let Some(reason) = &reason {
                        run.clone().on_disabled_hover_text(reason);
                    }
                    if run.clicked_once() {
                        start_run(panel, task, opt_panel, mol);
                    }
                }
                close = ui.button("Close editor").clicked_once();
            });
            ui.add(egui::Label::new(
                egui::RichText::new("This text is used exactly as written. Sync replaces it; closing keeps your draft.")
                    .small().color(WARNING),
            ).truncate());
            if let Some((message, error)) = &state.save_status {
                ui.add(egui::Label::new(egui::RichText::new(message).small()
                    .color(if *error { ERROR } else { ACCENT })).truncate())
                    .on_hover_text(message);
            } else if let Some(reason) = &reason {
                ui.add(egui::Label::new(egui::RichText::new(reason).small().color(ERROR)).truncate())
                    .on_hover_text(reason);
            }
            ui.separator();
            let width = ui.available_width();
            let height = ui.available_height().max(180.0);
            let rows = (height / ui.text_style_height(&egui::TextStyle::Monospace))
                .floor().max(8.0) as usize;
            egui::ScrollArea::both()
                .id_salt("qc_input_editor_scroll_v1")
                .auto_shrink([false, false])
                .max_height(height)
                .show(ui, |ui| {
                    let text = ui.add_enabled(
                        !running,
                        egui::TextEdit::multiline(&mut panel.advanced_text)
                            .id(egui::Id::new("qc_input_editor_text_v1"))
                            .code_editor()
                            .desired_width(width)
                            .desired_rows(rows),
                    );
                    if state.focus {
                        text.request_focus();
                        state.focus = false;
                    }
                    if text.changed() {
                        state.save_status = None;
                    }
                });
        });
    state.open = open && !close;
    store_state(ctx, state);
    result.map(|window| window.response.rect)
}
