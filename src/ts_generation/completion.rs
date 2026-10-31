//! Receive completed TS searches and retain their results.
use super::*;
use bevy::tasks::{block_on, futures_lite::future};

pub fn poll_generation(mut state: ResMut<TsGeneration>) {
    let Some(task) = state.task.as_mut() else {
        return;
    };
    let Some(result) = block_on(future::poll_once(task)) else {
        return;
    };
    state.task = None;
    state.duration = state.started.take().map(|t| t.elapsed());
    match result {
        Ok(output) => {
            state.message = Some(format!("{} TS guess generated.", output.algorithm.label()));
            // The energy profile opens by itself, as an optimisation's
            // energy plot does: it is the first thing worth looking at, and
            // it should not wait behind a button. Only when there is one --
            // not every method writes a profile.
            state.energy_plot_open = output.profile_dat.is_some();
            state.output = Some(output);
            state.is_error = false;
        }
        Err(error) => {
            state.message = Some(error);
            state.is_error = true;
        }
    }
}
