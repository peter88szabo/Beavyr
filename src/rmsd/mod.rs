//! Structure comparison: RMSD and optimal (Kabsch) superposition between
//! trajectory frames or loaded structures.
//!
//! This file deliberately contains declarations and Bevy wiring only --
//! every algorithm lives in `kabsch`, and the panel lives in `ui`.

pub mod kabsch;
pub mod ui;

use bevy::prelude::*;

/// Which frame everything is compared against, plus the last computed
/// numbers for the panel to show.
#[derive(Resource)]
pub struct RmsdState {
    /// `None` means "the first frame", the default reference. A frame the
    /// user explicitly pins with the panel's own button is stored here so
    /// scrubbing elsewhere does not silently move the reference.
    pub reference_frame: Option<usize>,
    /// Reference geometry captured at the moment it was pinned, so it stays
    /// valid even if the trajectory frames are later realigned in place.
    pub reference_atoms: Vec<String>,
    pub reference_pos: Vec<Vec3>,
    pub last_report: Option<kabsch::RmsdReport>,
    pub last_error: Option<String>,
    pub last_message: Option<String>,
}

impl Default for RmsdState {
    fn default() -> Self {
        Self {
            reference_frame: None,
            reference_atoms: Vec::new(),
            reference_pos: Vec::new(),
            last_report: None,
            last_error: None,
            last_message: None,
        }
    }
}

impl RmsdState {
    /// The frame index currently serving as the reference: whatever the user
    /// pinned, or frame 0 by default.
    pub fn effective_reference(&self) -> usize {
        self.reference_frame.unwrap_or(0)
    }
}
