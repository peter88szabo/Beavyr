//! Backbone torsions are graph rotations: no stretching, ring tearing or hidden relaxation.
use super::{
    peptide::{self, Build},
    peptide_edit::{atom, finish, metadata},
};
use crate::molecule::Molecule;
use bevy::prelude::*;
use std::collections::HashSet;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Phi,
    Psi,
    Omega,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Preset {
    Helix,
    Strand,
    TurnI,
}

pub fn quartet(mol: &Molecule, id: u64, kind: Kind) -> Result<[usize; 4], String> {
    let t = metadata(mol)?;
    let n = atom(mol, id, "N")?;
    let ca = atom(mol, id, "CA")?;
    let c = atom(mol, id, "C")?;
    let adj = peptide::adjacency(mol);
    let neighbor = |center: usize, name: &str| -> Result<usize, String> {
        let ns: Vec<_> = adj[center]
            .iter()
            .copied()
            .filter(|&i| t.atoms[i].residue != Some(id) && t.atoms[i].name == name)
            .collect();
        if ns.len() == 1 {
            Ok(ns[0])
        } else {
            Err("This torsion is undefined at a terminus, chain break or ambiguous branch.".into())
        }
    };
    let ids = match kind {
        Kind::Phi => [neighbor(n, "C")?, n, ca, c],
        Kind::Psi => [n, ca, c, neighbor(c, "N")?],
        Kind::Omega => {
            let pc = neighbor(n, "C")?;
            let prev = t.atoms[pc].residue.ok_or("Previous residue is unknown.")?;
            [atom(mol, prev, "CA")?, pc, n, ca]
        }
    };
    if !ids.windows(2).all(|p| adj[p[0]].contains(&p[1])) {
        return Err("Backbone connectivity is incomplete.".into());
    }
    Ok(ids)
}
pub fn measured(mol: &Molecule, id: u64, kind: Kind) -> Option<f32> {
    let [a, b, c, d] = quartet(mol, id, kind).ok()?;
    if (mol.pos[a] - mol.pos[b])
        .cross(mol.pos[c] - mol.pos[b])
        .length_squared()
        < 1e-8
        || (mol.pos[d] - mol.pos[c])
            .cross(mol.pos[b] - mol.pos[c])
            .length_squared()
            < 1e-8
    {
        return None;
    }
    let v = peptide::torsion(mol.pos[a], mol.pos[b], mol.pos[c], mol.pos[d]);
    v.is_finite().then_some(v)
}
pub fn wrap(x: f32) -> f32 {
    (x + 180.0).rem_euclid(360.0) - 180.0
}

pub fn moving_side(mol: &Molecule, b: usize, c: usize) -> Result<Vec<usize>, String> {
    let adj = peptide::adjacency(mol);
    if !adj.get(b).is_some_and(|ns| ns.contains(&c)) {
        return Err("The torsion axis is not a covalent bond.".into());
    }
    let mut side = HashSet::from([c]);
    let mut stack = vec![c];
    while let Some(i) = stack.pop() {
        for &j in &adj[i] {
            if (i == b && j == c) || (i == c && j == b) {
                continue;
            }
            if side.insert(j) {
                stack.push(j);
            }
        }
    }
    if side.contains(&b) {
        return Err("A ring or cross-link blocks this rotation. Proline φ and cyclic backbones require a closure solver.".into());
    }
    let mut side: Vec<_> = side.into_iter().collect();
    side.sort_unstable();
    Ok(side)
}
pub fn set_torsion(mol: &mut Molecule, ids: [usize; 4], target: f32) -> Result<(), String> {
    if !target.is_finite() {
        return Err("Torsion angles must be finite.".into());
    }
    let [a, b, c, d] = ids;
    if ids.iter().any(|&i| i >= mol.pos.len()) {
        return Err("Torsion atom no longer exists.".into());
    }
    let current = peptide::torsion(mol.pos[a], mol.pos[b], mol.pos[c], mol.pos[d]);
    if !current.is_finite()
        || (mol.pos[a] - mol.pos[b])
            .cross(mol.pos[c] - mol.pos[b])
            .length_squared()
            < 1e-8
        || (mol.pos[d] - mol.pos[c])
            .cross(mol.pos[b] - mol.pos[c])
            .length_squared()
            < 1e-8
    {
        return Err("Torsion is undefined for collinear atoms.".into());
    }
    let delta = wrap(target - current);
    if delta.abs() < 1e-4 {
        return Ok(());
    }
    let side = moving_side(mol, b, c)?;
    if !side.contains(&d) || side.contains(&a) {
        return Err("The torsion is not an independent four-atom path.".into());
    }
    let origin = mol.pos[b];
    let axis = (mol.pos[c] - origin).normalize();
    let rotation = Quat::from_axis_angle(axis, delta.to_radians());
    for i in side {
        if i != c {
            mol.pos[i] = origin + rotation * (mol.pos[i] - origin);
        }
    }
    let after = peptide::torsion(mol.pos[a], mol.pos[b], mol.pos[c], mol.pos[d]);
    if wrap(after - target).abs() > 0.05 {
        return Err("Could not reach the requested torsion.".into());
    }
    Ok(())
}
pub fn build_edit(
    original: &Molecule,
    mut result: Molecule,
    id: u64,
    message: String,
) -> Result<Build, String> {
    let picks = [
        atom(&result, id, "N")?,
        atom(&result, id, "CA")?,
        atom(&result, id, "C")?,
    ];
    result.recompute_bonds(1.2, 2.5);
    let moving = (0..original.atoms.len())
        .filter(|&i| result.pos[i] != original.pos[i])
        .collect();
    Ok(finish(
        result,
        (0..original.atoms.len()).map(Some).collect(),
        moving,
        picks,
        message,
    ))
}
pub fn edit(
    mol: &Molecule,
    id: u64,
    phi: Option<f32>,
    psi: Option<f32>,
    omega: Option<f32>,
) -> Result<Build, String> {
    let mut result = mol.clone();
    for (k, value) in [(Kind::Omega, omega), (Kind::Phi, phi), (Kind::Psi, psi)] {
        if let Some(v) = value {
            let ids = quartet(&result, id, k)?;
            set_torsion(&mut result, ids, v)?;
        }
    }
    build_edit(mol, result, id, "Backbone torsion edit".into())
}
pub fn preset(mol: &Molecule, first: u64, last: u64, preset: Preset) -> Result<Build, String> {
    let t = metadata(mol)?;
    let a = t
        .residues
        .iter()
        .position(|r| r.id == first)
        .ok_or("Select the first residue.")?;
    let b = t
        .residues
        .iter()
        .position(|r| r.id == last)
        .ok_or("Select the last residue.")?;
    if a > b {
        return Err("Choose the range in N-to-C order.".into());
    }
    let range = &t.residues[a..=b];
    if range.iter().any(|r| r.chain != range[0].chain) {
        return Err("A backbone preset cannot cross chains.".into());
    }
    if preset == Preset::TurnI && range.len() != 2 {
        return Err("A type-I turn preset acts on its two central residues; select exactly two adjacent residues.".into());
    }
    let adj = peptide::adjacency(mol);
    for r in range {
        let n = atom(mol, r.id, "N")?;
        let ca = atom(mol, r.id, "CA")?;
        let c = atom(mol, r.id, "C")?;
        if !adj[n].contains(&ca) || !adj[ca].contains(&c) {
            return Err("Backbone connectivity is incomplete.".into());
        }
    }
    for pair in range.windows(2) {
        let c = atom(mol, pair[0].id, "C")?;
        let n = atom(mol, pair[1].id, "N")?;
        if !mol
            .bonds
            .iter()
            .any(|&(x, y, _)| (x == c && y == n) || (y == c && x == n))
        {
            return Err("The range crosses a chain break.".into());
        }
    }
    let mut result = mol.clone();
    for (i, r) in range.iter().enumerate() {
        let (phi, psi) = match preset {
            Preset::Helix => (-57.0, -47.0),
            Preset::Strand => (-135.0, 135.0),
            Preset::TurnI => {
                if i == 0 {
                    (-60.0, -30.0)
                } else {
                    (-90.0, 0.0)
                }
            }
        };
        for (k, v) in [(Kind::Phi, phi), (Kind::Psi, psi)] {
            let center = atom(&result, r.id, if k == Kind::Phi { "N" } else { "C" })?;
            let other = if k == Kind::Phi { "C" } else { "N" };
            if adj[center]
                .iter()
                .any(|&j| t.atoms[j].residue != Some(r.id) && t.atoms[j].name == other)
            {
                let ids = quartet(&result, r.id, k)?;
                set_torsion(&mut result, ids, v)?;
            }
        }
    }
    build_edit(mol, result, first, format!("{preset:?} backbone preset"))
}
