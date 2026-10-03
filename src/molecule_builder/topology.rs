//! Explicit chemistry for structures assembled by the peptide builder.
//! Indices follow the molecule; IDs survive residue edits and atom compaction.
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AtomInfo {
    pub id: u64,
    pub residue: Option<u64>,
    pub name: String,
    pub charge: Option<i8>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ResidueInfo {
    pub id: u64,
    pub chain: u64,
    pub template: String,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Topology {
    pub elements: Vec<String>,
    pub atoms: Vec<AtomInfo>,
    pub residues: Vec<ResidueInfo>,
    /// Bond order: 1/2/3, or 4 for aromatic; 0 for an untyped host bond.
    pub bonds: Vec<(usize, usize, u8)>,
}

impl Topology {
    pub fn valid_for(&self, elements: &[String]) -> bool {
        use std::collections::HashSet;
        let ids: HashSet<_> = self.atoms.iter().map(|a| a.id).collect();
        let residues: HashSet<_> = self.residues.iter().map(|r| r.id).collect();
        self.elements == elements
            && self.atoms.len() == elements.len()
            && ids.len() == self.atoms.len()
            && residues.len() == self.residues.len()
            && self
                .atoms
                .iter()
                .all(|a| a.residue.is_none_or(|r| residues.contains(&r)))
            && self
                .bonds
                .iter()
                .all(|&(a, b, o)| a < elements.len() && b < elements.len() && a != b && o <= 4)
    }

    pub fn peptide_residue_count(&self) -> usize {
        self.residues
            .iter()
            .filter(|r| super::peptide_edit::template(&r.template).is_some())
            .count()
    }
    pub fn longest_chain(&self) -> usize {
        let mut counts = std::collections::HashMap::new();
        for r in &self.residues {
            if super::peptide_edit::template(&r.template).is_some() {
                *counts.entry(r.chain).or_insert(0) += 1;
            }
        }
        counts.into_values().max().unwrap_or(0)
    }
    pub fn backbone_atom(&self, i: usize) -> bool {
        self.atoms
            .get(i)
            .is_some_and(|a| a.residue.is_none() || matches!(a.name.as_str(), "N" | "CA" | "C"))
    }
    pub fn display_defaults(&self, settings: &mut crate::settings::MolSettings) {
        if self.peptide_residue_count() > 0 {
            settings.show_hbonds = false;
        }
        if self.longest_chain() > 20 {
            settings.representation = crate::settings::RepresentationMode::BackboneTrace;
        }
    }

    pub fn charge(&self) -> Option<i32> {
        self.atoms.iter().map(|a| a.charge.map(i32::from)).sum()
    }

    pub fn from_molecule(mol: &crate::molecule::Molecule) -> Self {
        mol.topology
            .as_ref()
            .filter(|t| t.valid_for(&mol.atoms))
            .cloned()
            .unwrap_or_else(|| Self {
                elements: mol.atoms.clone(),
                atoms: (0..mol.atoms.len())
                    .map(|i| AtomInfo {
                        id: i as u64 + 1,
                        residue: None,
                        name: String::new(),
                        charge: None,
                    })
                    .collect(),
                residues: vec![],
                bonds: mol.bonds.iter().map(|&(a, b, _)| (a, b, 0)).collect(),
            })
    }

    pub fn retain(&mut self, mapping: &[Option<usize>]) {
        self.elements = self
            .elements
            .iter()
            .zip(mapping)
            .filter(|(_, m)| m.is_some())
            .map(|(a, _)| a.clone())
            .collect();
        self.atoms = self
            .atoms
            .iter()
            .zip(mapping)
            .filter(|(_, m)| m.is_some())
            .map(|(a, _)| a.clone())
            .collect();
        self.bonds = self
            .bonds
            .iter()
            .filter_map(|&(a, b, o)| Some((mapping[a]?, mapping[b]?, o)))
            .collect();
        self.residues
            .retain(|r| self.atoms.iter().any(|a| a.residue == Some(r.id)));
    }
}
