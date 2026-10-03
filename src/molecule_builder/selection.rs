//! Multiple selection is independent of the Fragment Editor's ordered picks.
use std::collections::BTreeSet;

use bevy::prelude::*;
use bevy_egui::egui;

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

    pub fn expand(&mut self, mol: &Molecule, connected: bool) {
        let mut adjacency = vec![Vec::new(); mol.atoms.len()];
        for &(a, b, _) in &mol.bonds {
            if a < adjacency.len() && b < adjacency.len() {
                adjacency[a].push(b);
                adjacency[b].push(a);
            }
        }
        let mut queue: Vec<_> = self.atoms.iter().copied().collect();
        let initial = queue.len();
        let mut cursor = 0;
        while cursor < queue.len() && (connected || cursor < initial) {
            let i = queue[cursor];
            cursor += 1;
            if let Some(neighbours) = adjacency.get(i) {
                for &j in neighbours {
                    if !self.hidden.contains(&j) && self.atoms.insert(j) {
                        queue.push(j);
                    }
                }
            }
        }
    }

    pub fn invert(&mut self, count: usize) {
        self.atoms = (0..count)
            .filter(|i| !self.hidden.contains(i) && !self.atoms.contains(i))
            .collect();
    }

    /// Compact contextual tools. Nothing here changes coordinates or atom indices.
    pub fn controls(
        &mut self,
        ui: &mut egui::Ui,
        mol: &Molecule,
        settings: &mut MolSettings,
        placed: Option<&[usize]>,
    ) {
        let visible: BTreeSet<_> = (0..mol.atoms.len())
            .filter(|i| {
                !self.hidden.contains(i) && crate::scene::element_visible(&mol.atoms[*i], settings)
            })
            .collect();
        self.atoms.retain(|i| visible.contains(i));
        ui.horizontal_wrapped(|ui| {
            ui.menu_button(format!("Select · {}", self.atoms.len()), |ui| {
                ui.weak("Click replaces · Shift-click toggles");
                if let Some(placed) = placed {
                    if ui.button("Last placed fragment").clicked() {
                        self.atoms = placed
                            .iter()
                            .copied()
                            .filter(|i| *i < mol.atoms.len() && !self.hidden.contains(i))
                            .collect();
                        ui.close();
                    }
                }
                ui.menu_button("By element", |ui| {
                    for symbol in mol.atoms.iter().collect::<BTreeSet<_>>() {
                        if ui.button(symbol).clicked() {
                            self.atoms = mol
                                .atoms
                                .iter()
                                .enumerate()
                                .filter(|(i, s)| *s == symbol && !self.hidden.contains(i))
                                .map(|(i, _)| i)
                                .collect();
                            ui.close();
                        }
                    }
                });
                ui.add_enabled_ui(!self.atoms.is_empty(), |ui| {
                    if ui.button("Connected fragment").clicked() {
                        self.expand(mol, true);
                        ui.close();
                    }
                    if ui.button("Bonded neighbours").clicked() {
                        self.expand(mol, false);
                        ui.close();
                    }
                });
                if ui.button("All visible").clicked() {
                    self.atoms = (0..mol.atoms.len())
                        .filter(|i| !self.hidden.contains(i))
                        .collect();
                    ui.close();
                }
                if ui.button("Invert").clicked() {
                    self.invert(mol.atoms.len());
                    ui.close();
                }
                if ui.button("Clear selection").clicked() {
                    self.atoms.clear();
                    ui.close();
                }
            })
            .response
            .on_hover_text(
                "Select atoms for cleanup or freezing; Fragment Editor picks remain independent.",
            );
            self.atoms.retain(|i| visible.contains(i));
            if !self.atoms.is_empty() {
                if ui
                    .small_button("Freeze")
                    .on_hover_text("Keep selected atoms fixed during quick cleanup")
                    .clicked()
                {
                    self.frozen.extend(&self.atoms);
                }
                if ui
                    .small_button("Isolate")
                    .on_hover_text(
                        "Hide unselected atoms; calculations still include the entire structure",
                    )
                    .clicked()
                {
                    self.hidden = (0..mol.atoms.len())
                        .filter(|i| !self.atoms.contains(i))
                        .collect();
                    settings.geometry_dirty = true;
                }
            }
            if !self.hidden.is_empty() && ui.small_button("Show all").clicked() {
                self.hidden.clear();
                settings.geometry_dirty = true;
            }
            if !self.frozen.is_empty() {
                ui.menu_button(format!("Frozen · {}", self.frozen.len()), |ui| {
                    if ui.button("Unfreeze selected").clicked() {
                        self.frozen.retain(|i| !self.atoms.contains(i));
                        ui.close();
                    }
                    if ui.button("Unfreeze all").clicked() {
                        self.frozen.clear();
                        ui.close();
                    }
                });
            }
        });
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
    fn neighbours_and_connected_components_are_distinct() {
        let mut s = Selection::default();
        s.click(Some(0), false);
        s.expand(&molecule(), false);
        assert_eq!(s.atoms, [0, 1].into());
        s.expand(&molecule(), true);
        assert_eq!(s.atoms, [0, 1, 2].into());
    }
    #[test]
    fn toggle_invert_and_background_respect_visibility() {
        let mut s = Selection::default();
        s.click(Some(0), false);
        s.click(Some(1), true);
        s.click(Some(0), true);
        assert_eq!(s.atoms, [1].into());
        s.hidden.insert(3);
        s.invert(4);
        assert_eq!(s.atoms, [0, 2].into());
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
