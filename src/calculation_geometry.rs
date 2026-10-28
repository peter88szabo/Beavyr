//! Remove drawing placeholders before a calculation and retain an undo snapshot.

use crate::molecule::{is_dummy, Molecule};
use bevy_egui::egui;

const REMOVAL_ID: &str = "calculation_dummy_removal";

#[cfg(test)]
#[path = "calculation_geometry_tests.rs"]
mod tests;

/// Keep backend atom counts and coordinate arrays in agreement for every job type.
pub(crate) fn real_coordinates(
    atoms: &[String],
    positions: &[f64],
) -> Result<(Vec<String>, Vec<f64>), String> {
    if positions.len() != atoms.len() * 3 || positions.iter().any(|p| !p.is_finite()) {
        return Err("The geometry must have three finite coordinates for each atom.".into());
    }
    let mut real_atoms = Vec::new();
    let mut real_positions = Vec::new();
    for (symbol, xyz) in atoms.iter().zip(positions.chunks_exact(3)) {
        if !is_dummy(symbol) {
            real_atoms.push(symbol.clone());
            real_positions.extend_from_slice(xyz);
        }
    }
    if real_atoms.is_empty() {
        return Err(
            "There are no real atoms to calculate. Dummy atoms are drawing placeholders.".into(),
        );
    }
    Ok((real_atoms, real_positions))
}

#[derive(Clone)]
struct Removal {
    before: Molecule,
    mapping: Vec<Option<usize>>,
}

/// Prepare a calculation copy while retaining explicit bonds, charges and residues.
pub(crate) fn prepare(mol: &Molecule) -> Result<(Molecule, Vec<Option<usize>>), String> {
    if mol.atoms.len() != mol.pos.len() || mol.pos.iter().any(|p| !p.is_finite()) {
        return Err("The geometry must have one finite position for each atom.".into());
    }
    let removed = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, symbol)| is_dummy(symbol))
        .map(|(i, _)| i)
        .collect();
    let mut real = mol.clone();
    if real
        .topology
        .as_ref()
        .is_some_and(|t| !t.valid_for(&real.atoms))
    {
        real.topology = None;
    }
    let mapping = crate::molecule_builder::peptide_edit::remove_atoms(&mut real, &removed);
    if real.atoms.is_empty() {
        return Err(
            "There are no real atoms to calculate. Dummy atoms are drawing placeholders.".into(),
        );
    }
    if mol.atoms.len() != real.atoms.len() {
        real.recompute_bonds(1.2, 2.5);
    }
    Ok((real, mapping))
}

/// Called only by a calculation's Run button, before its input is captured.
/// The main UI consumes the snapshot in the same frame to update history and picks.
pub(crate) fn prepare_display(ctx: &egui::Context, mol: &mut Molecule) -> Result<(), String> {
    let (real, mapping) = prepare(mol)?;
    if real.atoms.len() != mol.atoms.len() {
        let removal = Removal {
            before: mol.clone(),
            mapping,
        };
        ctx.data_mut(|data| data.insert_temp(egui::Id::new(REMOVAL_ID), removal));
        *mol = real;
    }
    Ok(())
}

/// Finish the UI bookkeeping once all panels have released their molecule borrows.
#[allow(clippy::too_many_arguments)]
pub(crate) fn finish_display_removal(
    ctx: &egui::Context,
    mol: &Molecule,
    builder: &mut crate::molecule_builder::builder_ui::ZMatrixBuilderState,
    editor: &mut crate::molecule_builder::builder_ui::EditorRotateState,
    settings: &mut crate::settings::MolSettings,
    traj: &mut crate::trajectory::TrajectoryState,
    history: &mut crate::structure_history::StructureHistory,
    conformer: &mut crate::conformer::ConformerRun,
) -> bool {
    use crate::molecule_builder::{builder_ui::set_selected_atom, zmat2xyz::xyz_to_zmat};
    let Some(removal) = ctx.data_mut(|data| {
        let id = egui::Id::new(REMOVAL_ID);
        let removal = data.get_temp::<Removal>(id);
        data.remove::<Removal>(id);
        removal
    }) else {
        return false;
    };
    let before = &removal.before;
    let previous = history.before(before, traj);
    builder.history.record(
        "removed dummy atoms before calculation",
        &before.atoms,
        &before.pos,
        &xyz_to_zmat(&before.atoms, &before.pos),
    );
    builder.history.set_last_topology(before.topology.clone());
    builder.skip_next_sync = true;
    let remap = |indices: &std::collections::BTreeSet<usize>| {
        indices
            .iter()
            .filter_map(|&i| removal.mapping.get(i).copied().flatten())
            .collect()
    };
    let frozen = remap(&builder.selection.frozen);
    let hidden = remap(&builder.selection.hidden);
    builder.selection = Default::default();
    builder.selection.sync(mol);
    builder.selection.frozen = frozen;
    builder.selection.hidden = hidden;
    builder.zmat = xyz_to_zmat(&mol.atoms, &mol.pos);
    builder.edit_refresh = true;
    set_selected_atom(builder, &mol.atoms, None);
    builder.add_atom_active = false;
    builder.add_atom_auto = false;
    builder.add_atom_picks.clear();
    builder.add_atom_status = None;
    builder.edit_preview.clear();
    builder.last_click = None;
    builder.last_frag_snapshot = None;
    builder.frag_undo_visible = false;
    builder.last_remove_snapshot = None;
    builder.redo_remove_visible = false;
    builder.hydrogen_state = Default::default();
    builder.placed_editor = Default::default();
    let window_open = editor.window_open;
    *editor = Default::default();
    editor.window_open = window_open;
    crate::conformer::preparation::remap_display_settings(conformer, &removal.mapping);
    traj.clear_for_structure();
    history.replaced(previous, mol, traj, "Dummy atoms removed");
    settings.geometry_dirty = true;
    settings.bond_topology_dirty = true;
    true
}
