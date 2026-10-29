//! Actionable setup notices shared by calculation panels.

use crate::ui_style::ResponseExt;
use bevy_egui::egui;

#[derive(Clone, Default)]
struct Notices(Vec<(String, String)>);

/// Observe stored errors even while their setup panel is closed. Acknowledged
/// errors stay closed until the source clears or reports a different failure.
pub(crate) fn observe(ctx: &egui::Context, source: &str, title: &str, message: Option<&str>) {
    let message = message.filter(|m| !m.trim().is_empty() && !m.starts_with("Cancelled."));
    let id = egui::Id::new(("runtime_error_source", source));
    let changed = ctx.data_mut(|data| {
        let old = data.get_temp::<String>(id);
        match message {
            Some(message) => {
                data.insert_temp(id, message.to_owned());
                old.as_deref() != Some(message)
            }
            None => {
                data.remove::<String>(id);
                false
            }
        }
    });
    if changed {
        push(ctx, title, message.unwrap());
    }
}

pub(crate) fn push(ctx: &egui::Context, title: &str, message: &str) {
    ctx.data_mut(|data| {
        let pending = &mut data
            .get_temp_mut_or_default::<Notices>(egui::Id::new("runtime_setup_notices"))
            .0;
        let notice = (title.to_owned(), message.to_owned());
        if !pending.contains(&notice) {
            pending.push(notice);
        }
    });
}

/// A missing Python environment should be announced once, rather than every frame.
pub(crate) fn once(ctx: &egui::Context, title: &str, message: &str) {
    let first = ctx.data_mut(|data| {
        data.get_temp_mut_or_default::<std::collections::HashSet<String>>(egui::Id::new(
            "runtime_setup_notices_seen",
        ))
        .insert(format!("{title}\n{message}"))
    });
    if first {
        push(ctx, title, message);
    }
}

pub(crate) fn window(ctx: &egui::Context) -> Option<egui::Rect> {
    let (title, message) = ctx
        .data_mut(|data| data.get_temp::<Notices>(egui::Id::new("runtime_setup_notices")))?
        .0
        .first()?
        .clone();
    let response = egui::Modal::new(egui::Id::new("runtime_setup_modal"))
        .show(ctx, |ui| message_body(ui, &title, &message));
    if response.inner || response.should_close() {
        ctx.data_mut(|data| {
            let pending = &mut data
                .get_temp_mut_or_default::<Notices>(egui::Id::new("runtime_setup_notices"))
                .0;
            if !pending.is_empty() {
                pending.remove(0);
            }
        });
    }
    Some(response.response.rect)
}

/// Shared by setup errors and the DREIDING acknowledgement popup.
pub(crate) fn message_body(ui: &mut egui::Ui, title: &str, message: &str) -> bool {
    ui.set_width(540.0_f32.min((ui.ctx().content_rect().width() - 48.0).max(160.0)));
    ui.heading(title);
    ui.add_space(8.0);
    egui::ScrollArea::vertical()
        .max_height((ui.ctx().content_rect().height() - 200.0).max(100.0))
        .show(ui, |ui| {
            ui.add(egui::Label::new(message).selectable(true).wrap());
        });
    ui.add_space(12.0);
    ui.horizontal(|ui| {
        if ui.button("Copy error").clicked_once() {
            ui.ctx().copy_text(message.to_owned());
        }
        ui.button("OK").clicked_once()
    })
    .inner
}

#[cfg(test)]
#[path = "runtime_warning_tests.rs"]
mod tests;
