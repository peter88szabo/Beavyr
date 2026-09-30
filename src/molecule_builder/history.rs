//! One undo history for the whole builder.
//!
//! The builder used to have six Undo buttons -- one beside removing an atom,
//! two beside adding a fragment, and three in the fragment editor. Each
//! reversed only its own last action, so pressing one did not reverse what had
//! actually been done last. That is what made undo feel unreliable: it was
//! reliable, but about the wrong thing.
//!
//! This is a plain snapshot stack. Every action records the structure as it
//! stood *before* it ran, together with a name a chemist can read, and Undo
//! walks back through them in the order they happened. Snapshots rather than
//! inverse operations, because an inverse has to be written and kept correct
//! for every action separately, which is the same mistake in a new shape.

use bevy::prelude::Vec3;

use super::zmat2xyz::ZAtom;

/// How many steps are kept.
///
/// Far more than a building session produces, and small enough that the memory
/// is irrelevant: a step is three short vectors, so even a thousand-atom
/// structure costs tens of kilobytes a step. Past this the oldest go.
pub const MAX_STEPS: usize = 200;

/// The structure as it stood at one moment, and what was about to be done.
#[derive(Clone, Debug)]
pub struct BuilderStep {
    /// What the action was called, in the words the button used.
    pub name: String,
    pub atoms: Vec<String>,
    pub pos: Vec<Vec3>,
    pub zmat: Vec<ZAtom>,
}

/// The undo and redo stacks.
#[derive(Clone, Default, Debug)]
pub struct BuilderHistory {
    /// States from before each action, oldest first.
    past: Vec<BuilderStep>,
    /// States undone out of `past`, newest last.
    future: Vec<BuilderStep>,
}

impl BuilderHistory {
    /// Record the structure as it stands, before an action changes it.
    ///
    /// Recording a new action discards the redo stack: once the history has
    /// branched, the old forward path is not somewhere the user can get back
    /// to, and offering Redo for it would be a lie.
    pub fn record(&mut self, name: impl Into<String>, atoms: &[String], pos: &[Vec3], zmat: &[ZAtom]) {
        self.future.clear();
        self.past.push(BuilderStep {
            name: name.into(),
            atoms: atoms.to_vec(),
            pos: pos.to_vec(),
            zmat: zmat.to_vec(),
        });
        if self.past.len() > MAX_STEPS {
            self.past.remove(0);
        }
    }

    /// What Undo will take back, named, so it can be read before it is pressed.
    pub fn next_undo(&self) -> Option<&str> {
        self.past.last().map(|step| step.name.as_str())
    }

    /// What Redo will put back.
    pub fn next_redo(&self) -> Option<&str> {
        self.future.last().map(|step| step.name.as_str())
    }

    /// Step back one action. `current` is the structure as it stands now, kept
    /// so Redo can return to it.
    pub fn undo(&mut self, current: BuilderStep) -> Option<BuilderStep> {
        let step = self.past.pop()?;
        // The name travels with the action, not with the state: what Redo
        // will put back is the action just undone.
        let name = step.name.clone();
        self.future.push(BuilderStep { name, ..current });
        Some(step)
    }

    /// Step forward again.
    pub fn redo(&mut self, current: BuilderStep) -> Option<BuilderStep> {
        let step = self.future.pop()?;
        let name = step.name.clone();
        self.past.push(BuilderStep { name, ..current });
        Some(step)
    }

    /// Forget everything.
    ///
    /// Loading a different structure clears the history, because undoing
    /// across a structure change would restore atoms that no longer belong to
    /// anything.
    pub fn clear(&mut self) {
        self.past.clear();
        self.future.clear();
    }

    /// How many steps can be undone. For tests and for the panel's count.
    pub fn depth(&self) -> usize {
        self.past.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(symbol: &str) -> ZAtom {
        ZAtom {
            symbol: symbol.to_string(),
            bond_ref: None,
            bond_len: 0.0,
            angle_ref: None,
            angle_deg: 0.0,
            dihedral_ref: None,
            dihedral_deg: 0.0,
        }
    }

    fn state(n: usize) -> BuilderStep {
        BuilderStep {
            name: format!("step {n}"),
            atoms: vec!["C".to_string(); n],
            pos: vec![Vec3::splat(n as f32); n],
            zmat: vec![atom("C"); n],
        }
    }

    /// The property that makes a history worth trusting: undo then redo puts
    /// the structure back exactly, whatever the action was.
    #[test]
    fn undo_then_redo_returns_the_structure_unchanged() {
        let mut history = BuilderHistory::default();
        let before = state(2);
        history.record("added CH3", &before.atoms, &before.pos, &before.zmat);
        let after = state(5);

        let undone = history.undo(after.clone()).expect("something to undo");
        assert_eq!(undone.atoms, before.atoms);
        assert_eq!(undone.pos, before.pos);

        let redone = history.redo(undone).expect("something to redo");
        assert_eq!(redone.atoms, after.atoms);
        assert_eq!(redone.pos, after.pos);
        assert_eq!(redone.zmat.len(), after.zmat.len());
    }

    /// Undo must reverse what was last done, not what some particular button
    /// last did. Three actions of different kinds, undone in reverse order.
    #[test]
    fn undo_walks_back_in_the_order_things_happened() {
        let mut history = BuilderHistory::default();
        let s1 = state(1);
        let s2 = state(2);
        let s3 = state(3);
        history.record("added CH3", &s1.atoms, &s1.pos, &s1.zmat);
        history.record("deleted atom 4", &s2.atoms, &s2.pos, &s2.zmat);
        history.record("rotated about 3-7", &s3.atoms, &s3.pos, &s3.zmat);

        assert_eq!(history.next_undo(), Some("rotated about 3-7"));
        let back = history.undo(state(9)).unwrap();
        assert_eq!(back.atoms.len(), 3);
        assert_eq!(history.next_undo(), Some("deleted atom 4"));
        let back = history.undo(back).unwrap();
        assert_eq!(back.atoms.len(), 2);
        assert_eq!(history.next_undo(), Some("added CH3"));
        let back = history.undo(back).unwrap();
        assert_eq!(back.atoms.len(), 1);
        assert!(history.next_undo().is_none());
        assert!(history.undo(back).is_none());
    }

    /// Recording after undoing discards the redo stack. The forward path is
    /// somewhere the user can no longer reach, and offering Redo would lie.
    #[test]
    fn a_new_action_discards_the_redo_stack() {
        let mut history = BuilderHistory::default();
        let s1 = state(1);
        history.record("added CH3", &s1.atoms, &s1.pos, &s1.zmat);
        let back = history.undo(state(4)).unwrap();
        assert!(history.next_redo().is_some());

        history.record("added OH", &back.atoms, &back.pos, &back.zmat);
        assert!(history.next_redo().is_none(), "the old forward path is gone");
        assert_eq!(history.next_redo(), None);
    }

    /// The stack is bounded. A long session must not grow without limit, and
    /// the oldest steps are the ones nobody returns to.
    #[test]
    fn the_history_is_bounded_and_drops_the_oldest() {
        let mut history = BuilderHistory::default();
        for i in 0..(MAX_STEPS + 50) {
            let s = state(1);
            history.record(format!("step {i}"), &s.atoms, &s.pos, &s.zmat);
        }
        assert_eq!(history.depth(), MAX_STEPS);
        // The newest survived; the oldest did not.
        assert_eq!(history.next_undo(), Some(format!("step {}", MAX_STEPS + 49)).as_deref());
    }

    /// Loading a different structure clears it: undoing across that would
    /// restore atoms that no longer belong to anything.
    #[test]
    fn clearing_forgets_both_directions() {
        let mut history = BuilderHistory::default();
        let s = state(2);
        history.record("added CH3", &s.atoms, &s.pos, &s.zmat);
        history.undo(state(3));
        history.clear();
        assert!(history.next_undo().is_none() && history.next_redo().is_none());
        assert_eq!(history.depth(), 0);
    }
}
