// src/sync.rs
use bevy::prelude::*;

use crate::molecule::Molecule;
use crate::scene::AtomMarker;
use crate::events::{MoleculeChanged, MoleculeChangeReason};

/// System: sync `Molecule.pos` to atom entity `Transform`s whenever a
/// `MoleculeChanged` event fires.
pub fn sync_atoms_after_change(
    mut evr_changed: EventReader<MoleculeChanged>,
    mol: Option<Res<Molecule>>,
    mut q_atoms: Query<&mut Transform, With<AtomMarker>>,
) {
    let Some(mol) = mol else { return; };
    if evr_changed.is_empty() {
        return;
    }

    // Process all events; we only need the last one for actual sync.
    let last_event = evr_changed.read().last().cloned();
    if last_event.is_none() {
        return;
    }

    let positions = &mol.pos;
    for (i, &p) in positions.iter().enumerate() {
        if let Some(ent) = mol.atom_entities.get(i) {
            if let Ok(mut tf) = q_atoms.get_mut(*ent) {
                tf.translation = p;
            }
        }
    }
}

/// OPTIONAL strict feature:
/// After sync, check that every atom entity matches Molecule.pos.
/// If not, emit `AtomTransformsOutOfSync`.
#[cfg(feature = "strict_sync")]
pub fn validate_atom_sync(
    mol: Option<Res<Molecule>>,
    q_atoms: Query<&Transform, With<AtomMarker>>,
    mut evw: EventWriter<crate::events::AtomTransformsOutOfSync>,
) {
    let Some(mol) = mol else { return; };

    for (i, &p) in mol.pos.iter().enumerate() {
        if let Some(ent) = mol.atom_entities.get(i) {
            if let Ok(tf) = q_atoms.get(ent) {
                if (tf.translation - p).length() > 1e-5 {
                    evw.send(crate::events::AtomTransformsOutOfSync);
                    return;
                }
            }
        }
    }
}

