//! Hydrogen completion shared by the one-atom smart tool and the bulk C/O tool.
//! Work in the displayed Cartesian coordinates: existing atoms never move.

use bevy::prelude::Vec3;
use bevy_egui::egui;

use super::attach;
use super::zmat2xyz::{self, ZAtom};
use crate::bond_order::{
    common_valences, estimate_bond_order, ideal_bond_angle_deg, single_bond_length,
};
use crate::molecule::{covalent_radius_angstrom, Molecule};
use crate::settings::MolSettings;

/// Place one H using perceived connectivity, bond order, VSEPR and the space
/// around the host. Both manual smart addition and bulk completion call this.
pub(super) fn append_smart_hydrogen(mol: &mut Molecule, host: usize) -> Result<(), String> {
    if mol.atoms.len() != mol.pos.len() || mol.pos.iter().any(|p| !p.is_finite()) {
        return Err("The structure contains missing or non-finite coordinates.".into());
    }
    let Some(symbol) = mol.atoms.get(host) else {
        return Err("That atom no longer exists.".into());
    };
    let symbols: Vec<&str> = mol.atoms.iter().map(String::as_str).collect();
    let neighbors = attach::bonded_indices(&symbols, &mol.pos, host, &[]);
    let center = mol.pos[host];
    let orders: Vec<f64> = neighbors
        .iter()
        .map(|&i| estimate_bond_order(symbol, symbols[i], center.distance(mol.pos[i])))
        .collect();
    let total: f64 = orders.iter().sum();
    if let Some(valences) = common_valences(symbol) {
        if !valences.iter().any(|&v| total + 1.0 <= v + 1.0e-6) {
            return Err(format!(
                "{symbol}{} already has a full valence ({total:.1}).",
                host + 1
            ));
        }
    }
    if mol
        .pos
        .iter()
        .enumerate()
        .any(|(i, p)| i != host && center.distance(*p) < 0.01)
    {
        return Err(format!(
            "{symbol}{} overlaps another atom; separate them first.",
            host + 1
        ));
    }

    let dirs: Vec<Vec3> = neighbors
        .iter()
        .map(|&i| (mol.pos[i] - center).normalize())
        .collect();
    let max_order = orders.iter().copied().fold(1.0_f64, f64::max);
    let angle = if max_order >= 2.5 {
        180.0
    } else if max_order >= 1.25 {
        120.0
    } else {
        ideal_bond_angle_deg(symbol, dirs.len() + 1) as f32
    };
    let candidates = match dirs.as_slice() {
        [] => vec![Vec3::X, Vec3::Y, Vec3::Z, -Vec3::X, -Vec3::Y, -Vec3::Z],
        [axis] => attach::cone_directions(*axis, angle, 24),
        [a, b] => {
            let mut candidates = Vec::new();
            if let Some(d) = attach::open_valence_direction(&dirs, angle, max_order) {
                candidates.push(d);
                let normal = a.cross(*b).normalize_or_zero();
                // The two tetrahedral vacancies have opposite handedness.
                candidates.push((d - 2.0 * d.dot(normal) * normal).normalize_or_zero());
            } else {
                // Opposite existing bonds: choose in the plane perpendicular
                // to them, without collapsing onto either existing atom.
                candidates.extend(attach::cone_directions(*a, 90.0, 24));
            }
            candidates
        }
        _ => {
            if let Some(d) = attach::open_valence_direction(&dirs, angle, max_order) {
                vec![d]
            } else {
                let normal = dirs[0].cross(dirs[1]).normalize_or_zero();
                vec![normal, -normal]
            }
        }
    };
    let length = single_bond_length(symbol, "H");
    // Never place H where it would be perceived as bonded to a second atom.
    // Among equivalent rotations choose the largest minimum clearance.
    let best = candidates
        .into_iter()
        .filter_map(|d| {
            if !d.is_finite() || d.length_squared() < 0.5 {
                return None;
            }
            let p = center + length * d.normalize();
            let clearance = mol
                .pos
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != host)
                .map(|(i, q)| {
                    p.distance(*q)
                        / (1.2
                            * (covalent_radius_angstrom("H")
                                + covalent_radius_angstrom(symbols[i])))
                })
                .fold(f32::INFINITY, f32::min);
            (clearance > 1.0).then_some((p, clearance))
        })
        .max_by(|a, b| a.1.total_cmp(&b.1));
    let Some((position, _)) = best else {
        return Err(format!(
            "No clear H position at {symbol}{}; adjust or clean up the geometry.",
            host + 1
        ));
    };
    mol.atoms.push("H".into());
    mol.pos.push(position);
    Ok(())
}

/// Synchronize the editor after adding H, retaining the actual host as each
/// new row's bond reference rather than the preceding hydrogen in file order.
pub(super) fn builder_zmat(mol: &Molecule, first_added: usize) -> Vec<ZAtom> {
    let mut zmat = zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos);
    let symbols: Vec<&str> = mol.atoms.iter().map(String::as_str).collect();
    for index in first_added..mol.atoms.len() {
        if symbols[index] != "H" {
            continue;
        }
        let neighbors = attach::bonded_indices(&symbols, &mol.pos, index, &[]);
        let [host] = neighbors.as_slice() else {
            continue;
        };
        if *host >= index {
            continue;
        }
        let (angle_ref, dihedral_ref) = attach::host_reference_chain(&zmat[..index], *host);
        let direction = mol.pos[index] - mol.pos[*host];
        let (angle_deg, dihedral_deg) = match (angle_ref, dihedral_ref) {
            (Some(a), Some(d)) => {
                attach::internal_from_direction(mol.pos[*host], mol.pos[a], mol.pos[d], direction)
            }
            (Some(a), None) => (
                direction
                    .normalize()
                    .dot((mol.pos[a] - mol.pos[*host]).normalize())
                    .clamp(-1.0, 1.0)
                    .acos()
                    .to_degrees() as f64,
                0.0,
            ),
            _ => (0.0, 0.0),
        };
        zmat[index] = ZAtom {
            symbol: "H".into(),
            bond_ref: Some(host + 1),
            bond_len: direction.length() as f64,
            angle_ref: angle_ref.map(|i| i + 1),
            angle_deg,
            dihedral_ref: dihedral_ref.map(|i| i + 1),
            dihedral_deg,
        };
    }
    zmat
}

#[derive(Default, Debug)]
struct Completion {
    added: usize,
    skipped: usize,
}

fn complete(mol: &mut Molecule) -> Result<Completion, String> {
    if mol.atoms.len() != mol.pos.len() || mol.pos.iter().any(|p| !p.is_finite()) {
        return Err("The structure contains missing or non-finite coordinates.".into());
    }
    let original_count = mol.atoms.len();
    let mut result = Completion::default();
    for host in 0..original_count {
        let target = match mol.atoms[host].as_str() {
            "C" => 4,
            "O" => 2,
            _ => continue,
        };
        let symbols: Vec<&str> = mol.atoms.iter().map(String::as_str).collect();
        let neighbors = attach::bonded_indices(&symbols, &mol.pos, host, &[]);
        if neighbors.len() >= target {
            continue;
        }
        // This action assumes saturated C and neutral divalent O. A shorter
        // or stretched bond is evidence against that assumption, not another
        // missing H. In particular leave carbonyl and aromatic centres alone.
        if neighbors.iter().any(|&i| {
            estimate_bond_order(
                symbols[host],
                symbols[i],
                mol.pos[host].distance(mol.pos[i]),
            ) != 1.0
        }) {
            result.skipped += 1;
            continue;
        }
        for _ in neighbors.len()..target {
            if append_smart_hydrogen(mol, host).is_err() {
                result.skipped += 1;
                break;
            }
            result.added += 1;
        }
    }
    Ok(result)
}

#[derive(Clone)]
struct Undo {
    before: Molecule,
    after: Molecule,
}

#[derive(Clone, Default)]
pub struct HydrogenState {
    report: Option<Result<String, String>>,
    undo: Option<Undo>,
}

impl HydrogenState {
    pub(super) fn invalidate_if_changed(&mut self, mol: &Molecule) {
        if self
            .undo
            .as_ref()
            .is_some_and(|u| u.after.atoms != mol.atoms || u.after.pos != mol.pos)
        {
            self.undo = None;
            self.report = None;
        }
    }

    fn add(&mut self, mol: &mut Molecule, settings: &mut MolSettings) -> bool {
        self.invalidate_if_changed(mol);
        let before = mol.clone();
        match complete(mol) {
            Ok(result) => {
                self.report = Some(Ok(format!(
                    "Added {} H atom(s).{}",
                    result.added,
                    if result.skipped == 0 {
                        String::new()
                    } else {
                        format!(" Skipped {} C/O centre(s) with multiple/partial bonds or unsuitable geometry.", result.skipped)
                    }
                )));
                if result.added == 0 {
                    return false;
                }
                refresh(mol, settings);
                self.undo = Some(Undo {
                    before,
                    after: mol.clone(),
                });
                true
            }
            Err(message) => {
                self.report = Some(Err(message));
                false
            }
        }
    }

    fn undo(&mut self, mol: &mut Molecule, settings: &mut MolSettings) -> bool {
        self.invalidate_if_changed(mol);
        let Some(undo) = self.undo.take() else {
            return false;
        };
        *mol = undo.before;
        refresh(mol, settings);
        self.report = Some(Ok("Hydrogen addition undone.".into()));
        true
    }
}

fn refresh(mol: &mut Molecule, settings: &mut MolSettings) {
    mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
    settings.geometry_dirty = true;
    settings.bond_topology_dirty = true;
}

/// Returns true for a successful edit or undo, so the caller can synchronize
/// the builder and detach animation of the previous atom list.
pub(super) fn controls(
    ui: &mut egui::Ui,
    state: &mut HydrogenState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
) -> bool {
    state.invalidate_if_changed(mol);
    let mut changed = false;
    ui.horizontal(|ui| {
        if ui.add_enabled(!mol.atoms.is_empty(), egui::Button::new("Add missing H"))
            .on_hover_text("Uses the same chemically aware placement as Add Atom (auto): fills saturated sp3 carbon to four bonds and oxygen to two. Existing atoms stay fixed. Detected multiple/partial bonds are skipped.")
            .clicked() {
            changed = state.add(mol, settings);
        }
        if ui.add_enabled(state.undo.is_some(), egui::Button::new("Undo H addition")).clicked() {
            changed |= state.undo(mol, settings);
        }
    });
    ui.weak("Assumes sp3 C and neutral divalent O; other elements are unchanged.");
    match &state.report {
        Some(Ok(message)) => {
            ui.label(egui::RichText::new(message).small());
        }
        Some(Err(message)) => {
            ui.colored_label(egui::Color32::from_rgb(220, 120, 90), message);
        }
        None => {}
    }
    changed
}

#[cfg(test)]
mod tests;
