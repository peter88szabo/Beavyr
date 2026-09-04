use bevy::input::mouse::MouseMotion;
use bevy::prelude::*;
use bevy_egui::input::EguiWantsInput;

use crate::events::{AtomPicked, ToolKind, ViewportClicked};
use crate::measurements::Measurements;
use crate::molecule::Molecule;
use crate::molecule_builder::builder_ui::EditorRotateState;

pub mod screen;

pub struct PickerPlugin;

impl Plugin for PickerPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<ClickTracker>()
            .add_systems(Update, pick_and_emit_atom);
    }
}

/// How far the pointer may travel between press and release and still count as
/// a click rather than a camera orbit. Left-drag orbits the view, so without
/// this every rotation would land on whatever atom the drag started over -- or,
/// starting on empty space, would clear the selection.
const CLICK_TRAVEL_PX: f32 = 5.0;

/// Press-to-release state for one potential click.
///
/// The hit test runs at press time, because the camera grabs and hides the
/// cursor for the duration of a left drag, leaving no trustworthy cursor
/// position at release. The *decision* still waits for the release, so a press
/// that turns into an orbit reports nothing at all.
#[derive(Resource, Default)]
pub struct ClickTracker {
    /// What was under the cursor when the button went down. The outer option
    /// is "is a press in flight"; the inner one is the atom hit, where `None`
    /// means background -- itself a meaningful click, since it deselects.
    candidate: Option<Option<usize>>,
    travel: f32,
}

impl ClickTracker {
    fn press(&mut self, hit: Option<usize>) {
        self.candidate = Some(hit);
        self.travel = 0.0;
    }

    fn moved(&mut self, distance: f32) {
        self.travel += distance;
    }

    /// The click to report, or `None` if no press was in flight or the pointer
    /// travelled far enough to count as a drag.
    fn release(&mut self) -> Option<Option<usize>> {
        let candidate = self.candidate.take();
        let travelled = std::mem::take(&mut self.travel);
        candidate.filter(|_| travelled <= CLICK_TRAVEL_PX)
    }

    fn cancel(&mut self) {
        self.candidate = None;
        self.travel = 0.0;
    }
}

/// Centralized screen-space picking:
/// - Ignores clicks over egui panels
/// - Converts window cursor to PHYSICAL px
/// - Emits `ViewportClicked` for every deliberate click, atom or background
/// - Emits `AtomPicked` on top of that for each armed tool (Measurements /
///   builder rotate tool), which need several ordered picks
fn pick_and_emit_atom(
    egui_wants_input: Res<EguiWantsInput>,
    buttons: Res<ButtonInput<MouseButton>>,
    windows: Query<&Window>,
    q_cam: Query<(&Camera, &GlobalTransform), (With<Camera3d>, With<crate::scene::MainCamera>)>,
    mol: Option<Res<Molecule>>,
    measurements: Option<Res<Measurements>>,
    builder: Option<Res<EditorRotateState>>,
    panel_regions: Res<crate::ui::UiPanelRegions>,
    mut mouse_motion: MessageReader<MouseMotion>,
    mut tracker: ResMut<ClickTracker>,
    mut ev_pick: MessageWriter<AtomPicked>,
    mut ev_click: MessageWriter<ViewportClicked>,
) {
    // Drain the motion messages every frame, whether or not they are needed,
    // so a stale backlog can never be charged against a later press.
    let travel: f32 = mouse_motion.read().map(|m| m.delta.length()).sum();

    let Some(mol) = mol else {
        tracker.cancel();
        return;
    };
    if mol.pos.is_empty() {
        tracker.cancel();
        return;
    }

    if buttons.just_pressed(MouseButton::Left) {
        tracker.cancel();
        if let Some(hit) = press_target(
            &egui_wants_input,
            &windows,
            &q_cam,
            &mol,
            &panel_regions,
        ) {
            tracker.press(hit);
        }
    } else if buttons.pressed(MouseButton::Left) {
        // Not on the press frame itself: that frame's motion is the movement
        // that brought the pointer to the atom, not a drag away from it.
        tracker.moved(travel);
    }

    if !buttons.just_released(MouseButton::Left) {
        return;
    }
    let Some(hit) = tracker.release() else {
        return;
    };

    ev_click.write(ViewportClicked { hit });

    // Tools that need several atoms in a particular order stay explicitly
    // armed; a background click carries no atom for them.
    let Some(hit_idx) = hit else {
        return;
    };
    let want_measure = measurements
        .is_some_and(|m| m.is_active || m.angle_active || m.dihedral_active);
    let want_builder = builder.is_some_and(|b| b.active);
    if want_measure {
        ev_pick.write(AtomPicked {
            index: hit_idx,
            tool: ToolKind::Measurements,
        });
    }
    if want_builder {
        ev_pick.write(AtomPicked {
            index: hit_idx,
            tool: ToolKind::Builder,
        });
    }
}

/// The atom under the cursor at press time, or `Some(None)` for a press on
/// empty background. `None` means the press was not aimed at the 3D view at
/// all (it landed on a panel, or there is nothing to pick against) and so
/// should not start a click.
fn press_target(
    egui_wants_input: &EguiWantsInput,
    windows: &Query<&Window>,
    q_cam: &Query<(&Camera, &GlobalTransform), (With<Camera3d>, With<crate::scene::MainCamera>)>,
    mol: &Molecule,
    panel_regions: &crate::ui::UiPanelRegions,
) -> Option<Option<usize>> {
    // Skip while egui has captured the pointer or a popup is open.
    if egui_wants_input.is_using_pointer() || egui_wants_input.is_popup_open() {
        return None;
    }

    let window = windows.single().ok()?;

    // Same blind spot as the orbit camera: egui does not know the hand-built
    // control panels occupy the screen, and `is_using_pointer` is still false
    // on the frame a button is first pressed.  Without this, clicking a panel
    // control also picks whatever atom happens to sit behind it.
    let mouse_logical = window.cursor_position()?;
    if panel_regions.covers(mouse_logical) {
        return None;
    }
    let (cam, cam_xf) = q_cam.single().ok()?;

    // Convert logical -> physical px for world_to_viewport
    let mouse_px = mouse_logical * window.scale_factor() as f32;

    Some(
        screen::find_nearest_atom_screen_space(cam, cam_xf, &mol.pos, mouse_px, 18.0)
            .map(|(hit_idx, _)| hit_idx),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_stationary_press_and_release_is_a_click() {
        let mut tracker = ClickTracker::default();
        tracker.press(Some(3));
        tracker.moved(1.0);
        tracker.moved(1.5);
        assert_eq!(tracker.release(), Some(Some(3)));
    }

    #[test]
    fn a_press_that_drags_away_reports_nothing() {
        let mut tracker = ClickTracker::default();
        tracker.press(Some(3));
        tracker.moved(CLICK_TRAVEL_PX + 0.1);
        assert_eq!(
            tracker.release(),
            None,
            "a drag orbits the camera and must leave the selection alone"
        );
    }

    #[test]
    fn small_wobbles_accumulate_into_a_drag() {
        let mut tracker = ClickTracker::default();
        tracker.press(None);
        for _ in 0..10 {
            tracker.moved(1.0);
        }
        assert_eq!(tracker.release(), None);
    }

    #[test]
    fn a_background_click_is_reported_as_a_click_with_no_atom() {
        let mut tracker = ClickTracker::default();
        tracker.press(None);
        assert_eq!(
            tracker.release(),
            Some(None),
            "clicking empty space is what clears the selection"
        );
    }

    #[test]
    fn releasing_without_a_press_reports_nothing() {
        let mut tracker = ClickTracker::default();
        assert_eq!(tracker.release(), None);
    }

    #[test]
    fn a_cancelled_press_does_not_survive_into_the_next_release() {
        let mut tracker = ClickTracker::default();
        tracker.press(Some(1));
        tracker.cancel();
        assert_eq!(tracker.release(), None);
    }

    #[test]
    fn travel_does_not_carry_over_from_an_earlier_drag() {
        let mut tracker = ClickTracker::default();
        tracker.press(Some(0));
        tracker.moved(50.0);
        assert_eq!(tracker.release(), None);

        tracker.press(Some(2));
        assert_eq!(
            tracker.release(),
            Some(Some(2)),
            "the previous drag's travel must not disqualify this click"
        );
    }
}
