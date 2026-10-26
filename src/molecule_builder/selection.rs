//! Multiple selection is independent of the Fragment Editor's ordered picks.
use std::collections::BTreeSet;

use bevy::prelude::*;

use super::builder_ui::ZMatrixBuilderState;
use crate::{molecule::Molecule, settings::MolSettings};

#[derive(Clone, Default)]
pub struct Selection {
    pub atoms: BTreeSet<usize>,
    pub frozen: BTreeSet<usize>,
    pub hidden: BTreeSet<usize>,
    symbols: Vec<String>,
}

impl Selection {
    pub fn sync(&mut self, mol: &Molecule) -> bool {
        if self.symbols == mol.atoms {
            return false;
        }
        if mol.atoms.starts_with(&self.symbols) {
            self.symbols = mol.atoms.clone();
            return false;
        }
        let hidden = !self.hidden.is_empty();
        *self = Self {
            symbols: mol.atoms.clone(),
            ..Default::default()
        };
        hidden
    }

    pub fn click(&mut self, hit: Option<usize>, toggle: bool) {
        if !toggle {
            self.atoms.clear();
        }
        if let Some(i) = hit {
            if !self.hidden.contains(&i) && !self.atoms.insert(i) && toggle {
                self.atoms.remove(&i);
            }
        }
    }

    pub fn active_atom(&self, preferred: Option<usize>, previous: Option<usize>) -> Option<usize> {
        preferred
            .filter(|i| self.atoms.contains(i))
            .or_else(|| previous.filter(|i| self.atoms.contains(i)))
            .or_else(|| self.atoms.first().copied())
    }
}

pub fn sync(
    mut builder: ResMut<ZMatrixBuilderState>,
    mol: Res<Molecule>,
    mut settings: ResMut<MolSettings>,
) {
    if builder.selection.sync(&mol) {
        settings.geometry_dirty = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn molecule() -> Molecule {
        let mut mol = Molecule::from_xyz("5\n\nC 0 0 0\nC 1.5 0 0\nH 2.5 0 0\nO 8 0 0\nH 9 0 0\n");
        mol.bonds = vec![(0, 1, 1.5), (1, 2, 1.0), (3, 4, 1.0)];
        mol
    }
    #[test]
    fn toggle_and_background_clicks() {
        let mut s = Selection::default();
        s.click(Some(0), false);
        s.click(Some(1), true);
        s.click(Some(0), true);
        assert_eq!(s.atoms, [1].into());
        s.click(Some(2), true);
        s.click(None, true);
        assert_eq!(s.atoms.len(), 2);
        s.click(None, false);
        assert!(s.atoms.is_empty());
    }
    #[test]
    fn new_atoms_clear_stale_selection_freezing_and_visibility() {
        let mut s = Selection::default();
        let mut m = molecule();
        s.sync(&m);
        s.atoms.insert(2);
        s.frozen.insert(1);
        s.hidden.insert(0);
        m.pos[0].x += 0.1;
        assert!(!s.sync(&m));
        assert_eq!(s.frozen.len(), 1);
        m.atoms.remove(0);
        assert!(s.sync(&m));
        assert!(s.atoms.is_empty() && s.frozen.is_empty() && s.hidden.is_empty());
    }

    #[test]
    fn appending_a_fragment_preserves_existing_frozen_atoms() {
        let mut s = Selection::default();
        let mut m = molecule();
        s.sync(&m);
        s.frozen = [0, 1].into();
        m.atoms.push("H".into());
        m.pos.push(Vec3::ZERO);
        s.sync(&m);
        assert_eq!(s.frozen, [0, 1].into());
    }

    #[test]
    fn toggling_off_an_atom_also_removes_its_active_highlight() {
        let mut s = Selection::default();
        s.click(Some(0), false);
        s.click(Some(1), true);
        assert_eq!(s.active_atom(Some(1), Some(0)), Some(1));
        s.click(Some(1), true);
        assert_eq!(s.active_atom(Some(1), Some(1)), Some(0));
        s.click(Some(0), true);
        assert_eq!(s.active_atom(Some(0), Some(0)), None);
    }
}
