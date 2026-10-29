//! Compact controls for trajectory playback.

use super::{PlaybackMode, TrajectoryState};
use crate::ui_style::ResponseExt;
use bevy_egui::egui;

pub(crate) fn play_button(ui: &mut egui::Ui, direction: i32) -> egui::Response {
    let label = if direction < 0 {
        "Play backward"
    } else {
        "Play forward"
    };
    let button = if direction > 0 {
        crate::ui_style::primary("")
    } else {
        egui::Button::new("")
    };
    let response = ui
        .add(button.min_size(egui::vec2(36.0, 30.0)))
        .on_hover_text(label);
    response.widget_info(|| {
        egui::WidgetInfo::labeled(egui::WidgetType::Button, ui.is_enabled(), label)
    });
    if ui.is_rect_visible(response.rect) {
        let center = response.rect.center();
        let sign = if direction < 0 { -1.0 } else { 1.0 };
        let mut points = vec![
            center + egui::vec2(-sign * 5.0, -7.0),
            center + egui::vec2(sign * 6.0, 0.0),
            center + egui::vec2(-sign * 5.0, 7.0),
        ];
        if direction < 0 {
            points.reverse();
        }
        let color = if direction > 0 {
            egui::Color32::WHITE
        } else {
            ui.style().interact(&response).fg_stroke.color
        };
        ui.painter().add(egui::Shape::convex_polygon(
            points,
            color,
            egui::Stroke::NONE,
        ));
    }
    response
}

pub(crate) fn mode_row(ui: &mut egui::Ui, traj: &mut TrajectoryState) {
    ui.scope(|ui| {
        ui.spacing_mut().button_padding.x = 6.0;
        ui.spacing_mut().item_spacing.x = 6.0;
        ui.horizontal(|ui| {
            if ui
                .add(egui::Button::new("Fixed bonds").selected(traj.fixed_bonds))
                .on_hover_text("Keep the current bond topology during playback.")
                .clicked_once()
            {
                traj.fixed_bonds = !traj.fixed_bonds;
                traj.last_applied = None;
                if traj.fixed_bonds {
                    traj.overlay_dirty = true;
                }
            }
            for (label, mode) in [
                ("Loop", PlaybackMode::Loop),
                ("Once", PlaybackMode::Once),
                ("PingPong", PlaybackMode::PingPong),
            ] {
                if ui
                    .add(egui::Button::new(label).selected(traj.mode == mode))
                    .clicked_once()
                {
                    traj.mode = mode;
                }
            }
        });
    });
}

pub(crate) fn stop_when_panel_closes(traj: &mut TrajectoryState, was_open: bool, is_open: bool) {
    if was_open && !is_open {
        traj.playing = false;
        traj.accum = 0.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn closing_the_panel_stops_at_the_current_frame_and_reopening_does_not_resume() {
        let mut traj = TrajectoryState {
            playing: true,
            current_frame: 4,
            accum: 0.04,
            ..Default::default()
        };
        stop_when_panel_closes(&mut traj, true, true);
        assert!(traj.playing);
        stop_when_panel_closes(&mut traj, true, false);
        assert!(!traj.playing);
        assert_eq!(traj.current_frame, 4);
        assert_eq!(traj.accum, 0.0);
        stop_when_panel_closes(&mut traj, false, true);
        assert!(!traj.playing);
    }

    #[test]
    fn a_closed_trajectory_panel_does_not_stop_a_frequency_animation() {
        let mut traj = TrajectoryState {
            playing: true,
            mode: PlaybackMode::Loop,
            ..Default::default()
        };
        stop_when_panel_closes(&mut traj, false, false);
        assert!(traj.playing);
    }
}
