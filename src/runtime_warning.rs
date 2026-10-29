//! Actionable setup notices shared by calculation panels.

use bevy_egui::egui;

#[derive(Clone, Default)]
struct Notices(Vec<(String, String)>);

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
    let response = egui::Modal::new(egui::Id::new("runtime_setup_modal")).show(ctx, |ui| {
        ui.set_width(480.0_f32.min((ctx.content_rect().width() - 48.0).max(160.0)));
        ui.heading(title);
        ui.add_space(8.0);
        egui::ScrollArea::vertical()
            .max_height((ctx.content_rect().height() - 160.0).max(100.0))
            .show(ui, |ui| {
                ui.label(message);
            });
        ui.add_space(12.0);
        ui.button("OK").clicked()
    });
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

#[cfg(test)]
#[path = "runtime_warning_tests.rs"]
mod tests;
