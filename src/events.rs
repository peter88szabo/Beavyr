use bevy::prelude::*;

/// Which tool requested the pick.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ToolKind {
    Measurements,
    Builder,
}

/// Emitted when the shared picker reports a picked atom for a specific tool.
#[derive(Event, Debug, Clone, Copy)]
pub struct AtomPicked {
    pub index: usize,
    pub tool: ToolKind,
}

/// Reasons why the molecule changed (useful for UI/state integration).
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MoleculeChangeReason {
    /// Parsed/loaded text XYZ; recenter camera?
    ParseXyz { recenter: bool },
    /// Programmatic position change via `Molecule::set_pos`.
    SetPos,
    /// Fragment rotation applied in builder.
    BuilderRotate,
    /// Undo/Redo action (if you add stacks later).
    Undo,
    Redo,
}

/// Emitted whenever `Molecule` data changes and views should refresh.
#[derive(Event, Debug, Clone, Copy)]
pub struct MoleculeChanged {
    pub reason: MoleculeChangeReason,
}

impl MoleculeChanged {
    pub fn parse_xyz(recenter: bool) -> Self {
        Self { reason: MoleculeChangeReason::ParseXyz { recenter } }
    }
    pub fn set_pos() -> Self {
        Self { reason: MoleculeChangeReason::SetPos }
    }
    #[allow(dead_code)]
    pub fn builder_rotate() -> Self {
        Self { reason: MoleculeChangeReason::BuilderRotate }
    }
}
