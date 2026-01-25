use bevy::prelude::*;
use bevy_egui::EguiContexts;

use crate::events::{AtomPicked, ToolKind};
use crate::molecule::Molecule;
use crate::measurements::Measurements;
use crate::molecule_builder::builder_ui::{EditorRotateState, ZMatrixBuilderState};

pub mod screen;

pub struct PickerPlugin;

impl Plugin for PickerPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Update, pick_and_emit_atom);
    }
}

/// Centralized screen-space picking:
/// - Ignores clicks over egui panels
/// - Converts window cursor to PHYSICAL px
/// - Emits AtomPicked for each active tool (Measurements / Builder)
fn pick_and_emit_atom(
    mut contexts: EguiContexts,
    buttons: Res<ButtonInput<MouseButton>>,
    windows: Query<&Window>,
    q_cam: Query<(&Camera, &GlobalTransform), (With<Camera3d>, With<crate::scene::MainCamera>)>,
    mol: Option<Res<Molecule>>,
    measurements: Option<Res<Measurements>>,
    builder: Option<Res<EditorRotateState>>,
    zmat_builder: Option<Res<ZMatrixBuilderState>>,
    mut ev_pick: EventWriter<AtomPicked>,
) {
    let Some(mol) = mol else { return; };
    if mol.pos.is_empty() { return; }

    // Which tools want picking right now?
    let mut want_measure = false;
    let mut want_builder = false;
    if let Some(m) = &measurements {
        want_measure = m.is_active || m.angle_active || m.dihedral_active;
    }
    if let Some(b) = &builder {
        want_builder = b.active;
    }
    if let Some(z) = &zmat_builder {
        want_builder = want_builder || z.pick_active || z.frag_pick_active;
    }
    if !(want_measure || want_builder) { return; }

    // Must be a fresh LMB click
    if !buttons.just_pressed(MouseButton::Left) { return; }

    // Skip if pointer is over any egui area (panels/overlays)
    let Ok(ctx) = contexts.ctx_mut() else { return; };
    if ctx.is_pointer_over_area() { return; }

    let Ok(window) = windows.single() else { return; };
    let Ok((cam, cam_xf)) = q_cam.single() else { return; };
    let Some(mouse_logical) = window.cursor_position() else { return; };

    // Convert logical -> physical px for world_to_viewport
    let mouse_px = mouse_logical * window.scale_factor() as f32;

    if let Some((hit_idx, _)) = screen::find_nearest_atom_screen_space(
        cam, cam_xf, &mol.pos, mouse_px, 18.0,
    ) {
        if want_measure {
            ev_pick.write(AtomPicked { index: hit_idx, tool: ToolKind::Measurements });
        }
        if want_builder {
            ev_pick.write(AtomPicked { index: hit_idx, tool: ToolKind::Builder });
        }
    }
}
