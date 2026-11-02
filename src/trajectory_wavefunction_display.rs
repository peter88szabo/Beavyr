//! Adopt only the current frame's wavefunction; keep at most one pending disk read.
use super::*;
use crate::trajectory_wavefunctions::{job, State};
#[derive(Resource, Default)]
struct Preview {
    task: Option<Task<Result<(OrbitalData, String), String>>>,
    requested: Option<(PathBuf, usize)>,
    adopted: Option<(PathBuf, usize)>,
    choice: Option<(Quantity, Spin, Option<usize>, f32)>,
    tracking_key: (bool,usize),
    tracked_selection: Option<(Spin,usize)>,
}
pub struct Plugin;
impl bevy::prelude::Plugin for Plugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<Preview>()
            .add_systems(
                PostUpdate,
                update
                    .after(bevy_egui::EguiPostUpdateSet::EndPass)
                    .before(clear_orbitals_on_structure_change),
            )
            .add_systems(PostUpdate, ready.after(rebuild_orbital_surface));
    }
}
fn update(
    mut preview: ResMut<Preview>,
    mut state: ResMut<State>,
    mut orbital: ResMut<OrbitalState>,
    mut traj: ResMut<crate::trajectory::TrajectoryState>,
    mol: Res<crate::molecule::Molecule>,
) {
    let Some(series) = state
        .series
        .clone()
        .filter(|s| state.follow && s.matches(&traj))
    else {
        preview.task = None;
        preview.requested = None;
        preview.adopted = None;
        state.ready = None;
        return;
    };
    let key = (series.root.clone(), traj.current_frame);
    if state.tracking.enabled && preview.adopted.as_ref()==Some(&key) {
        if let Some((spin,index))=preview.tracked_selection {
            if orbital.quantity!=Quantity::Orbital || orbital.spin!=spin || orbital.selected!=Some(index) {
                state.tracking.enabled=false;state.tracking.message="Orbital selection changed. Track this selection to choose a new anchor.".into();
            }
        }
    }
    let tracking_key=(state.tracking.enabled,state.tracking.track.as_ref().map_or(0,|t|std::sync::Arc::as_ptr(t)as usize));
    if preview.tracking_key!=tracking_key {preview.tracking_key=tracking_key;preview.adopted=None;preview.requested=None;preview.task=None;preview.tracked_selection=None;}

    if orbital.data.is_some() {
        preview.choice = Some((
            orbital.quantity,
            orbital.spin,
            orbital.selected,
            orbital.isovalue,
        ));
    }
    if preview.adopted.as_ref() == Some(&key) {
        if orbital.data.is_none() {
            state.follow = false;
            state.ready = None;
        }
        return;
    }
    state.ready = None;
    if preview.requested.as_ref() != Some(&key) {
        let index = traj.current_frame;
        preview.task =
            Some(AsyncComputeTaskPool::get().spawn(async move { job::load(&series, index) }));
        preview.requested = Some(key.clone());
    }
    let Some(task) = preview.task.as_mut() else {
        return;
    };
    if let Some(result) = block_on(future::poll_once(task)) {
        preview.task = None;
        match result {
            Ok((mut data, text)) if describes(&data, &mol) => {
                let choice = preview.choice;
                if choice.is_some_and(|(q, _, _, _)| q == Quantity::SpinDensity)
                    && !data.is_open_shell()
                {
                    state.error = Some("This frame's Molden file has no separate alpha/beta orbitals. Choose total density or an orbital, or compute spin-resolved wavefunctions with another backend.".into());
                    state.follow = false;
                    return;
                }
                let tracked=if state.tracking.enabled {
                    let track=state.tracking.track.as_ref().unwrap();
                    match crate::orbital_tracking::apply(&mut data,track,key.1) {
                        Ok(c)=>Some((track.spin,c.index)),
                        Err(e)=>{traj.playing=false;state.follow=false;state.error=Some(e);return;}
                    }
                }else{None};
                orbital.adopt(
                    data,
                    Some(format!("Frame {}", key.1 + 1)),
                    Some("Trajectory wavefunction".into()),
                    Some(text),
                );
                if let Some((quantity, spin, selected, isovalue)) = choice {
                    orbital.quantity = quantity;
                    orbital.spin = spin;
                    orbital.isovalue = isovalue;
                    orbital.selected = selected
                        .filter(|i| *i < orbital.data.as_ref().unwrap().orbitals(spin).len());
                }
                if let Some((spin,index))=tracked {orbital.quantity=Quantity::Orbital;orbital.spin=spin;orbital.selected=Some(index);}
                preview.tracked_selection=tracked;
                if orbital.quantity == Quantity::Orbital && orbital.selected.is_none() {
                    state.follow = false;
                    state.error=Some("The selected orbital index is unavailable in this frame. Choose an available orbital before resuming.".into());
                }
                preview.adopted = Some(key);
            }
            Ok(_) => {
                state.follow = false;
                state.error =
                    Some("The displayed geometry changed while loading its wavefunction.".into());
            }
            Err(e) => {
                state.follow = false;
                state.error = Some(e);
            }
        }
    }
}
fn ready(
    mut state: ResMut<State>,
    orbital: Res<OrbitalState>,
    traj: Res<crate::trajectory::TrajectoryState>,
    mol: Res<crate::molecule::Molecule>,
) {
    if state.follow
        && state.series.as_ref().is_some_and(|s| s.matches(&traj))
        && orbital.data.as_ref().is_some_and(|d| describes(d, &mol))
        && !orbital.is_busy()
        && !orbital.surface_dirty
        && (!orbital.show_surface || orbital.field.is_some())
    {
        state.ready = Some(traj.current_frame);
    } else {
        state.ready = None;
    }
}
