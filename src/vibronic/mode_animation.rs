//! Explicit mode animation requests, consumed by the viewer after drawing UI.
use super::input::StateModes;
use crate::{
    molecule::Molecule,
    trajectory::{TrajectoryFrame, TrajectoryState},
};
use bevy::prelude::Vec3;
use bevy_egui::egui;
#[derive(Default)]
pub struct Animation {
    pub excited: bool,
    pub mode: usize,
    pub amplitude: f32,
    pub pending: Option<(TrajectoryFrame, Vec<TrajectoryFrame>)>,
    pub stop: bool,
    owned: Option<(TrajectoryFrame, TrajectoryFrame)>,
}
fn same(a: &TrajectoryFrame, b: &TrajectoryFrame) -> bool {
    a.atoms == b.atoms && a.pos == b.pos
}
impl Animation {
    /// True while this animation owns the displayed trajectory.
    pub fn is_running(&self) -> bool {
        self.owned.is_some()
    }
    pub fn needs_apply(&self, open: bool) -> bool {
        self.pending.is_some() || self.stop || (!open && self.owned.is_some())
    }
    pub fn apply(&mut self, open: bool, mol: &mut Molecule, traj: &mut TrajectoryState) -> bool {
        // A different trajectory/structure revokes ownership; closing this window
        // must never replace a subsequently loaded user's geometry.
        if let Some((_, first)) = &self.owned {
            if !traj.frames.first().is_some_and(|f| same(f, first)) {
                self.owned = None;
            }
        }
        if let Some((equilibrium, frames)) = self.pending.take() {
            let first = frames[0].clone();
            traj.clear_for_structure();
            traj.frames = frames;
            traj.playing = true;
            traj.mode = crate::trajectory::PlaybackMode::Loop;
            traj.fps = 30.;
            mol.atoms = first.atoms.clone();
            mol.pos = first.pos.clone();
            mol.topology = None;
            self.owned = Some((equilibrium, first));
            self.stop = false;
            return true;
        }
        if self.stop || !open {
            self.stop = false;
            if let Some((equilibrium, _)) = self.owned.take() {
                crate::qchem_interfaces::xtb_freq::stop_mode_animation(
                    traj,
                    &equilibrium.atoms,
                    &equilibrium.pos,
                );
                mol.atoms = equilibrium.atoms;
                mol.pos = equilibrium.pos;
                mol.topology = None;
                return true;
            }
        }
        false
    }
}
pub fn controls(ui: &mut egui::Ui, a: &mut Animation, ground: &StateModes, excited: &StateModes) {
    use crate::ui_style::ResponseExt;
    ui.horizontal_wrapped(|ui| {
        ui.label("Animate mode");
        let ground_clicked = ui
            .selectable_value(&mut a.excited, false, "Ground")
            .clicked();
        let excited_clicked = ui
            .selectable_value(&mut a.excited, true, "Excited")
            .clicked();
        if ground_clicked || excited_clicked {
            a.mode = 0;
        }
        let state = if a.excited { excited } else { ground };
        a.mode = a.mode.min(state.selected.len().saturating_sub(1));
        egui::ComboBox::from_id_salt("vibronic_animate_mode")
            .selected_text(format!(
                "{} · {:.1} cm⁻¹",
                state.selected[a.mode] + 1,
                state.frequencies[state.selected[a.mode]]
            ))
            .show_ui(ui, |ui| {
                for (k, &i) in state.selected.iter().enumerate() {
                    ui.selectable_value(
                        &mut a.mode,
                        k,
                        format!("{} · {:.1} cm⁻¹", i + 1, state.frequencies[i]),
                    );
                }
            });
        if a.amplitude <= 0. {
            a.amplitude = 0.3;
        }
        ui.add(
            egui::DragValue::new(&mut a.amplitude)
                .range(0.01..=2.)
                .speed(0.01)
                .suffix(" Å"),
        );
        if ui.button("Animate").clicked_once() {
            let pos = state
                .coords
                .chunks_exact(3)
                .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32) * 0.5291772)
                .collect::<Vec<_>>();
            let frames = crate::qchem_interfaces::xtb_freq::animate_mode(
                &state.atoms,
                &pos,
                &state.modes,
                state.selected[a.mode],
                a.amplitude,
                60,
            );
            a.pending = Some((
                TrajectoryFrame { metadata: Default::default(),
                    atoms: state.atoms.clone(),
                    pos,
                },
                frames,
            ));
        }
        if a.owned.is_some() && ui.button("Stop").clicked_once() {
            a.stop = true;
        }
    });
}
