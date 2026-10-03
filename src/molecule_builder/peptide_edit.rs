//! Residue edits retain the existing backbone and explicit chemical topology.
use super::{
    amino_acids::{AminoAcid, AminoAcidTemplate, BondOrder, TEMPLATES},
    peptide::{self, Build, Residue},
    topology::{AtomInfo, Topology},
};
use crate::molecule::Molecule;
use bevy::prelude::*;
use std::collections::HashSet;

pub fn template(id: &str) -> Option<&'static AminoAcidTemplate> {
    TEMPLATES.iter().find(|t| t.id == id)
}
pub fn metadata(mol: &Molecule) -> Result<&Topology, String> {
    mol.topology.as_ref().filter(|t| t.valid_for(&mol.atoms)).ok_or_else(|| "This structure has no residue annotations. Build a peptide or open its Beavyr project first.".into())
}
pub fn atom(mol: &Molecule, residue: u64, name: &str) -> Result<usize, String> {
    metadata(mol)?
        .atoms
        .iter()
        .position(|a| a.residue == Some(residue) && a.name == name)
        .ok_or_else(|| format!("Residue {residue} has no {name} atom."))
}
pub fn remove_atoms(mol: &mut Molecule, removed: &HashSet<usize>) -> Vec<Option<usize>> {
    let mut next = 0;
    let mapping: Vec<_> = (0..mol.atoms.len())
        .map(|i| {
            if removed.contains(&i) {
                None
            } else {
                let n = next;
                next += 1;
                Some(n)
            }
        })
        .collect();
    if let Some(t) = mol.topology.as_mut() {
        t.retain(&mapping);
    }
    mol.atoms = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(i, _)| !removed.contains(i))
        .map(|(_, a)| a.clone())
        .collect();
    mol.pos = mol
        .pos
        .iter()
        .enumerate()
        .filter(|(i, _)| !removed.contains(i))
        .map(|(_, p)| *p)
        .collect();
    mol.bonds = mol
        .bonds
        .iter()
        .filter_map(|&(a, b, d)| Some((mapping[a]?, mapping[b]?, d)))
        .collect();
    mol.hydrogen_bonds.clear();
    mapping
}
fn finish(
    mut mol: Molecule,
    mapping: Vec<Option<usize>>,
    moving: Vec<usize>,
    picks: [usize; 3],
    message: String,
) -> Build {
    mol.recompute_bonds(1.2, 2.5);
    Build {
        charge: mol.topology.as_ref().and_then(Topology::charge),
        clashes: peptide::clashes(&mol, 0),
        removed_host_atoms: mapping.iter().filter(|i| i.is_none()).count(),
        sequence: message,
        host_mapping: mapping,
        molecule: mol,
        moving,
        picks,
    }
}
fn order(o: BondOrder) -> u8 {
    match o {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
    }
}

pub fn replace(mol: &Molecule, id: u64, residue: &Residue, pack: bool) -> Result<Build, String> {
    let t = metadata(mol)?;
    let info = t
        .residues
        .iter()
        .find(|r| r.id == id)
        .ok_or("Select a residue.")?;
    let old = template(&info.template).ok_or("Caps cannot be replaced with amino acids.")?;
    if old.amino_acid == AminoAcid::Proline || residue.amino == AminoAcid::Proline {
        return Err("Proline replacement requires rebuilding the ring and backbone. Change it in the sequence before building.".into());
    }
    let new = residue
        .amino
        .template(old.termini, residue.side_chain)
        .ok_or("Unsupported residue state.")?;
    let n = atom(mol, id, "N")?;
    let ca = atom(mol, id, "CA")?;
    let c = atom(mol, id, "C")?;
    let adj = peptide::adjacency(mol);
    // Preserve the peptide plane, terminal groups and their H. Replace CA-H too (Gly has two).
    let keep: HashSet<_> = t
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| {
            a.residue == Some(id) && matches!(a.name.as_str(), "N" | "CA" | "C" | "O" | "OXT")
        })
        .map(|(i, _)| i)
        .collect();
    let removed: HashSet<_> = t
        .atoms
        .iter()
        .enumerate()
        .filter(|(i, a)| {
            a.residue == Some(id)
                && !keep.contains(i)
                && !(mol.atoms[*i] == "H" && adj[*i].iter().any(|j| keep.contains(j) && *j != ca))
        })
        .map(|(i, _)| i)
        .collect();
    if t.bonds.iter().any(|&(a, b, _)| {
        (removed.contains(&a) && t.atoms[b].residue != Some(id))
            || (removed.contains(&b) && t.atoms[a].residue != Some(id))
    }) {
        return Err("Remove this residue's side-chain cross-link before replacing it.".into());
    }
    let mut positions: Vec<_> = new
        .atoms
        .iter()
        .map(|a| Vec3::from_array(a.position_angstrom.map(|x| x as f32)))
        .collect();
    let from = [
        positions[new.backbone.alpha_carbon],
        positions[new.backbone.nitrogen],
        positions[new.backbone.carbonyl_carbon],
    ];
    peptide::transform(&mut positions, from, [mol.pos[ca], mol.pos[n], mol.pos[c]]);
    let mut result = mol.clone();
    let mapping = remove_atoms(&mut result, &removed);
    let topology = result.topology.as_mut().unwrap();
    let mut next = t.atoms.iter().map(|a| a.id).max().unwrap_or(0) + 1;
    let mut local = vec![None; new.atoms.len()];
    let mut moving = vec![];
    for (i, a) in new.atoms.iter().enumerate() {
        if let Some(index) = topology
            .atoms
            .iter()
            .position(|x| x.residue == Some(id) && x.name == a.name)
        {
            local[i] = Some(index);
            continue;
        }
        let hydrogen_on_backbone = a.element == "H"
            && new.bonds.iter().any(|b| {
                b.atoms.contains(&i)
                    && b.atoms
                        .iter()
                        .any(|&j| j != i && matches!(new.atoms[j].name, "N" | "OXT"))
            });
        if matches!(a.name, "N" | "CA" | "C" | "O" | "OXT") || hydrogen_on_backbone {
            continue;
        }
        let index = result.atoms.len();
        local[i] = Some(index);
        moving.push(index);
        result.atoms.push(a.element.into());
        result.pos.push(positions[i]);
        topology.atoms.push(AtomInfo {
            id: next,
            residue: Some(id),
            name: a.name.into(),
            charge: Some(a.formal_charge),
        });
        next += 1;
    }
    for b in new.bonds {
        if let (Some(a), Some(bi)) = (local[b.atoms[0]], local[b.atoms[1]]) {
            if moving.contains(&a) || moving.contains(&bi) {
                topology.bonds.push((a, bi, order(b.order)));
            }
        }
    }
    topology
        .residues
        .iter_mut()
        .find(|r| r.id == id)
        .unwrap()
        .template = new.id.into();
    topology.elements = result.atoms.clone();
    result.recompute_bonds(1.2, 2.5);
    if pack {
        pack_sidechains(&mut result, Some(id))?;
    }
    let picks = [
        mapping[n].unwrap(),
        mapping[ca].unwrap(),
        mapping[c].unwrap(),
    ];
    Ok(finish(
        result,
        mapping,
        moving,
        picks,
        format!("replace residue {id} with {}", residue.amino.three_letter()),
    ))
}

/// Deterministic torsion sampling. Ring bonds and bonds crossing residue boundaries never rotate.
/// This is local steric fitting, not a statistical rotamer library or an energy minimizer.
pub fn pack_sidechains(mol: &mut Molecule, only: Option<u64>) -> Result<usize, String> {
    let t = metadata(mol)?.clone();
    let adj = peptide::adjacency(mol);
    let mut rotations = vec![];
    for &(a, b, o) in &t.bonds {
        if o != 1 || mol.atoms[a] == "H" || mol.atoms[b] == "H" {
            continue;
        }
        let Some(r) = t.atoms[a].residue else {
            continue;
        };
        if t.atoms[b].residue != Some(r) || only.is_some_and(|id| id != r) {
            continue;
        }
        if t.bonds.iter().any(|&(i, j, _)| {
            [(i, j), (j, i)].iter().any(|&(x, y)| {
                t.atoms[x].residue == Some(r)
                    && t.atoms[y].residue != Some(r)
                    && !matches!(t.atoms[x].name.as_str(), "N" | "C")
            })
        }) {
            continue;
        }
        if [a, b]
            .iter()
            .any(|&i| matches!(t.atoms[i].name.as_str(), "N" | "C" | "O" | "OXT"))
        {
            continue;
        }
        let mut side = HashSet::from([b]);
        let mut stack = vec![b];
        while let Some(i) = stack.pop() {
            for &j in &adj[i] {
                if (i == b && j == a) || (i == a && j == b) {
                    continue;
                }
                if side.insert(j) {
                    stack.push(j);
                }
            }
        }
        if side.contains(&a) {
            continue;
        } // ring
        let side: Vec<_> = if side
            .iter()
            .any(|&i| t.atoms[i].residue != Some(r) || t.atoms[i].name == "CA")
        {
            (0..mol.atoms.len())
                .filter(|i| !side.contains(i) && t.atoms[*i].residue == Some(r))
                .collect()
        } else {
            side.into_iter().collect()
        };
        if side
            .iter()
            .any(|&i| matches!(t.atoms[i].name.as_str(), "N" | "CA" | "C" | "O" | "OXT"))
        {
            continue;
        }
        if !side
            .iter()
            .any(|&i| i != a && i != b && mol.atoms[i] != "H")
        {
            continue;
        }
        rotations.push((a, b, side));
    }
    let mut changed = 0;
    for _ in 0..2 {
        for (a, b, side) in &rotations {
            let center = mol.pos[*a];
            let axis = (mol.pos[*b] - center).normalize_or_zero();
            if axis.length_squared() < 0.5 {
                continue;
            }
            let original: Vec<_> = side.iter().map(|&i| mol.pos[i]).collect();
            let mut best = original.clone();
            let mut score = steric_score(mol, side, &adj);
            let initial = score;
            for step in 1..6 {
                let q = Quat::from_axis_angle(axis, step as f32 * std::f32::consts::PI / 3.0);
                for (&i, &p) in side.iter().zip(&original) {
                    mol.pos[i] = center + q * (p - center);
                }
                let candidate = steric_score(mol, side, &adj);
                if candidate + 1e-5 < score {
                    score = candidate;
                    best = side.iter().map(|&i| mol.pos[i]).collect();
                }
            }
            for (&i, &p) in side.iter().zip(&best) {
                mol.pos[i] = p;
            }
            if score + 1e-5 < initial {
                changed += 1;
            }
        }
    }
    mol.recompute_bonds(1.2, 2.5);
    Ok(changed)
}
fn steric_score(mol: &Molecule, moving: &[usize], adj: &[Vec<usize>]) -> f32 {
    let mut score = 0.0;
    for &i in moving {
        for j in 0..mol.atoms.len() {
            if moving.contains(&j)
                || adj[i].contains(&j)
                || adj[i].iter().any(|&k| adj[k].contains(&j))
            {
                continue;
            }
            let radius = |s: &str| match s {
                "H" => 1.1,
                "N" => 1.55,
                "O" => 1.52,
                "S" | "Se" => 1.8,
                _ => 1.7,
            };
            let overlap = (0.8 * (radius(&mol.atoms[i]) + radius(&mol.atoms[j]))
                - mol.pos[i].distance(mol.pos[j]))
            .max(0.0);
            score += overlap * overlap;
        }
    }
    score
}
pub fn pack(mol: &Molecule, id: Option<u64>) -> Result<Build, String> {
    let mut result = mol.clone();
    let count = pack_sidechains(&mut result, id)?;
    let t = metadata(&result)?;
    let r = id
        .or_else(|| t.residues.first().map(|r| r.id))
        .ok_or("No residues.")?;
    let picks = [
        atom(&result, r, "N")?,
        atom(&result, r, "CA")?,
        atom(&result, r, "C")?,
    ];
    let moving = (0..mol.atoms.len())
        .filter(|&i| mol.pos[i] != result.pos[i])
        .collect();
    Ok(finish(
        result,
        (0..mol.atoms.len()).map(Some).collect(),
        moving,
        picks,
        format!("fit side chains ({count} torsion changes)"),
    ))
}

/// Restore named template hydrogens using a local heavy-atom frame. Peptide amides
/// use the actual adjacent carbonyl, not the free monomer nitrogen geometry.
pub fn restore_hydrogens(mol: &mut Molecule) -> Result<usize, String> {
    let initial = mol.atoms.len();
    let residues = metadata(mol)?.residues.clone();
    for r in residues {
        let Some(t) = template(&r.template) else {
            continue;
        };
        let positions: Vec<_> = t
            .atoms
            .iter()
            .map(|a| Vec3::from_array(a.position_angstrom.map(|x| x as f32)))
            .collect();
        for (hi, h) in t.atoms.iter().enumerate().filter(|(_, a)| a.element == "H") {
            if atom(mol, r.id, h.name).is_ok() {
                continue;
            }
            let Some(bond) = t.bonds.iter().find(|b| b.atoms.contains(&hi)) else {
                continue;
            };
            let ti = if bond.atoms[0] == hi {
                bond.atoms[1]
            } else {
                bond.atoms[0]
            };
            let Ok(host) = atom(mol, r.id, t.atoms[ti].name) else {
                continue;
            };
            let adj = peptide::adjacency(mol);
            let external = adj[host].iter().copied().find(|&i| {
                mol.topology.as_ref().unwrap().atoms[i].residue != Some(r.id) && mol.atoms[i] != "H"
            });
            // Condensed N has one H (zero for Pro); cross-linked sulfur has none.
            if external.is_some() && t.atoms[ti].name == "SG" {
                continue;
            }
            let amide = t.atoms[ti].name == "N" && external.is_some();
            if amide
                && (t.amino_acid == AminoAcid::Proline
                    || adj[host].iter().any(|&i| mol.atoms[i] == "H"))
            {
                continue;
            }
            let expected = t
                .bonds
                .iter()
                .filter(|b| {
                    b.atoms.contains(&ti) && b.atoms.iter().any(|&i| t.atoms[i].element == "H")
                })
                .count();
            if !amide && adj[host].iter().filter(|&&i| mol.atoms[i] == "H").count() >= expected {
                continue;
            }
            let p = if amide {
                let ca = atom(mol, r.id, "CA")?;
                mol.pos[host]
                    - 1.01
                        * ((mol.pos[external.unwrap()] - mol.pos[host]).normalize()
                            + (mol.pos[ca] - mol.pos[host]).normalize())
                        .normalize()
            } else {
                let mut anchors: Vec<_> = t
                    .atoms
                    .iter()
                    .enumerate()
                    .filter(|(i, a)| *i != ti && a.element != "H")
                    .filter_map(|(i, a)| atom(mol, r.id, a.name).ok().map(|j| (i, j)))
                    .collect();
                anchors.sort_by(|a, b| {
                    positions[a.0]
                        .distance_squared(positions[ti])
                        .total_cmp(&positions[b.0].distance_squared(positions[ti]))
                });
                let Some(&(ai, aj)) = anchors.first() else {
                    continue;
                };
                let Some(&(bi, bj)) = anchors.iter().skip(1).find(|(bi, bj)| {
                    (positions[ai] - positions[ti])
                        .cross(positions[*bi] - positions[ti])
                        .length_squared()
                        > 0.01
                        && (mol.pos[aj] - mol.pos[host])
                            .cross(mol.pos[*bj] - mol.pos[host])
                            .length_squared()
                            > 0.01
                }) else {
                    continue;
                };
                let mut point = [positions[hi]];
                peptide::transform(
                    &mut point,
                    [positions[ti], positions[ai], positions[bi]],
                    [mol.pos[host], mol.pos[aj], mol.pos[bj]],
                );
                point[0]
            };
            if !p.is_finite() {
                return Err("Cannot restore H at a degenerate residue geometry.".into());
            }
            let index = mol.atoms.len();
            let topo = mol.topology.as_mut().unwrap();
            let id = topo.atoms.iter().map(|a| a.id).max().unwrap_or(0) + 1;
            topo.atoms.push(AtomInfo {
                id,
                residue: Some(r.id),
                name: h.name.into(),
                charge: Some(0),
            });
            topo.elements.push("H".into());
            topo.bonds.push((host, index, 1));
            mol.atoms.push("H".into());
            mol.pos.push(p);
            mol.bonds.push((host, index, p.distance(mol.pos[host])));
        }
    }
    Ok(mol.atoms.len() - initial)
}

pub fn pdb(mol: &Molecule) -> Result<String, String> {
    use std::fmt::Write;
    let t = metadata(mol)?;
    if mol.atoms.len() > 99999
        || mol
            .pos
            .iter()
            .any(|p| !p.is_finite() || p.min_element() < -999.999 || p.max_element() > 9999.999)
    {
        return Err("Structure exceeds PDB coordinate or atom limits.".into());
    }
    let mut chains: Vec<_> = t.residues.iter().map(|r| r.chain).collect();
    chains.sort_unstable();
    chains.dedup();
    if chains.len() > 62 || t.residues.len() > 9999 {
        return Err("Structure exceeds PDB chain or residue limits.".into());
    }
    let alphabet = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    let mut text = String::from("REMARK   Generated by Beavyr; explicit residue connectivity\n");
    for (i, a) in t.atoms.iter().enumerate() {
        let info = a
            .residue
            .and_then(|id| t.residues.iter().enumerate().find(|(_, r)| r.id == id));
        let (res, name, chain) = if let Some((ri, r)) = info {
            (
                ri + 1,
                template(&r.template)
                    .map(|x| x.amino_acid.three_letter())
                    .unwrap_or(&r.template),
                alphabet[chains.iter().position(|c| *c == r.chain).unwrap()] as char,
            )
        } else {
            (1, "UNK", ' ')
        };
        let aname = if a.name.is_empty() {
            mol.atoms[i].as_str()
        } else {
            a.name.as_str()
        };
        let record = if info.is_some_and(|(_, r)| template(&r.template).is_some()) {
            "ATOM  "
        } else {
            "HETATM"
        };
        let aligned = if aname.len() < 4 && mol.atoms[i].len() == 1 {
            format!(" {:<3}", aname)
        } else {
            format!("{:<4}", aname)
        };
        let charge = match a.charge {
            Some(q) if q != 0 => format!("{}{}", q.unsigned_abs(), if q > 0 { '+' } else { '-' }),
            _ => "  ".into(),
        };
        let p = mol.pos[i];
        writeln!(text,"{record}{:>5} {:4} {:>3} {}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}{:>2}",i+1,aligned,name,chain,res,p.x,p.y,p.z,1.0,0.0,mol.atoms[i],charge).unwrap();
    }
    for &(a, b, _) in &t.bonds {
        writeln!(text, "CONECT{:>5}{:>5}", a + 1, b + 1).unwrap();
    }
    text.push_str("END\n");
    Ok(text)
}

/// Cap a free terminus without moving any retained host atom.
pub fn cap(mol: &Molecule, id: u64, end: peptide::End) -> Result<Build, String> {
    let index = atom(mol, id, if end == peptide::End::N { "N" } else { "C" })?;
    let anchor = peptide::attachment(mol, index, end)?;
    let source = metadata(mol)?;
    let chain = source
        .residues
        .iter()
        .find(|r| r.id == id)
        .ok_or("Missing residue")?
        .chain;
    let mut result = mol.clone();
    let removed = anchor.leaving.iter().copied().collect();
    let mapping = remove_atoms(&mut result, &removed);
    let host = mapping[index].unwrap();
    let neighbor = mapping[anchor.neighbor].unwrap();
    let p = result.pos[host];
    let q = result.pos[neighbor];
    let connector = p + anchor.direction * 1.34;
    let mut coords = vec![connector];
    let mut symbols = vec![];
    let mut names = vec![];
    let mut edges = vec![];
    let name = if end == peptide::End::N {
        let methyl = peptide::place(q, p, connector, 1.51, 114.0, 180.0);
        let oxygen = connector
            - 1.23 * ((p - connector).normalize() + (methyl - connector).normalize()).normalize();
        coords.extend([oxygen, methyl]);
        symbols.extend(["C", "O", "C"]);
        names.extend(["C", "O", "CH3"]);
        edges.extend([(0, 1, 2), (0, 2, 1)]);
        "ACE"
    } else {
        let methyl = peptide::place(q, p, connector, 1.46, 123.0, 180.0);
        let hydrogen = connector
            - 1.01 * ((p - connector).normalize() + (methyl - connector).normalize()).normalize();
        coords.extend([hydrogen, methyl]);
        symbols.extend(["N", "H", "C"]);
        names.extend(["N", "H", "CH3"]);
        edges.extend([(0, 1, 1), (0, 2, 1)]);
        "NME"
    };
    for (k, torsion) in [60.0, 180.0, 300.0].into_iter().enumerate() {
        coords.push(peptide::place(
            p, connector, coords[2], 1.09, 109.471, torsion,
        ));
        symbols.push("H");
        names.push(["HH31", "HH32", "HH33"][k]);
        edges.push((2, 3 + k, 1));
    }
    let offset = result.atoms.len();
    let t = result.topology.as_mut().unwrap();
    t.atoms[host].charge = Some(0);
    let rid = t.residues.iter().map(|r| r.id).max().unwrap_or(0) + 1;
    let next = source.atoms.iter().map(|a| a.id).max().unwrap_or(0) + 1;
    let at =
        t.residues.iter().position(|r| r.id == id).unwrap() + usize::from(end == peptide::End::C);
    t.residues.insert(
        at,
        super::topology::ResidueInfo {
            id: rid,
            chain,
            template: name.into(),
        },
    );
    for i in 0..coords.len() {
        result.atoms.push(symbols[i].into());
        result.pos.push(coords[i]);
        t.atoms.push(AtomInfo {
            id: next + i as u64,
            residue: Some(rid),
            name: names[i].into(),
            charge: Some(0),
        });
    }
    t.elements = result.atoms.clone();
    t.bonds.push((host, offset, 1));
    t.bonds.extend(
        edges
            .into_iter()
            .map(|(a, b, o)| (a + offset, b + offset, o)),
    );
    let moving = (offset..result.atoms.len()).collect();
    Ok(finish(
        result,
        mapping,
        moving,
        [host, offset, offset + 2],
        format!("add {name} cap"),
    ))
}

/// Connect only an already positioned pair. Do not silently stretch a remote bond.
pub fn disulfide(mol: &Molecule, first: u64, second: u64) -> Result<Build, String> {
    if first == second {
        return Err("Choose two different cysteine residues.".into());
    }
    let a = atom(mol, first, "SG")?;
    let b = atom(mol, second, "SG")?;
    let adj = peptide::adjacency(mol);
    if !(1.8..=2.3).contains(&mol.pos[a].distance(mol.pos[b])) {
        return Err(format!("Sulfur atoms are {:.2} Å apart. Position them at 1.8–2.3 Å before linking; this action does not fold the chain.",mol.pos[a].distance(mol.pos[b])));
    }
    for (i, j) in [(a, b), (b, a)] {
        let heavy: Vec<_> = adj[i]
            .iter()
            .copied()
            .filter(|&k| mol.atoms[k] != "H")
            .collect();
        if heavy.len() != 1 || mol.atoms[heavy[0]] != "C" {
            return Err("Choose two free cysteine thiols/thiolates.".into());
        }
        let angle = (mol.pos[heavy[0]] - mol.pos[i])
            .normalize()
            .dot((mol.pos[j] - mol.pos[i]).normalize())
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees();
        if !(80.0..=125.0).contains(&angle) {
            return Err(
                "The C–S–S angle is unsuitable. Adjust the side chains before linking.".into(),
            );
        }
    }
    let removed = adj[a]
        .iter()
        .chain(&adj[b])
        .copied()
        .filter(|&i| mol.atoms[i] == "H")
        .collect();
    let mut result = mol.clone();
    let mapping = remove_atoms(&mut result, &removed);
    let a = mapping[a].unwrap();
    let b = mapping[b].unwrap();
    let t = result.topology.as_mut().unwrap();
    t.atoms[a].charge = Some(0);
    t.atoms[b].charge = Some(0);
    t.bonds.push((a, b, 1));
    let ca = atom(&result, first, "CB")?;
    Ok(finish(
        result,
        mapping,
        vec![],
        [ca, a, b],
        "form disulfide".into(),
    ))
}

#[cfg(test)]
mod tests {
    use super::super::peptide::{parse_sequence, Recipe};
    use super::*;
    fn chain(sequence: &str) -> Molecule {
        peptide::build(
            &Recipe {
                residues: parse_sequence(sequence).unwrap(),
                ..Default::default()
            },
            None,
        )
        .unwrap()
        .molecule
    }
    fn strip_h(mol: &mut Molecule) {
        let removed = mol
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.as_str() == "H")
            .map(|(i, _)| i)
            .collect();
        remove_atoms(mol, &removed);
    }
    #[test]
    fn peptide_explicit_bonds_survive_recompute_and_charge_is_persistent() {
        let mut mol = chain("AKD");
        let before = mol.topology.clone();
        let bonds = mol.bonds.len();
        mol.recompute_bonds(3.0, 3.0);
        assert_eq!(mol.bonds.len(), bonds);
        assert_eq!(mol.topology, before);
        let text = ron::to_string(mol.topology.as_ref().unwrap()).unwrap();
        let t: Topology = ron::from_str(&text).unwrap();
        assert!(t.valid_for(&mol.atoms));
        assert_eq!(t.charge(), Some(0));
    }
    #[test]
    fn peptide_replacement_retains_backbone_ids_and_coordinates() {
        let mol = chain("AGK");
        let replaced = replace(&mol, 2, &Residue::new(AminoAcid::Tryptophan), true).unwrap();
        for name in ["N", "CA", "C", "O"] {
            let a = atom(&mol, 2, name).unwrap();
            let b = atom(&replaced.molecule, 2, name).unwrap();
            assert_eq!(mol.pos[a], replaced.molecule.pos[b]);
            assert_eq!(
                mol.topology.as_ref().unwrap().atoms[a].id,
                replaced.molecule.topology.as_ref().unwrap().atoms[b].id
            );
        }
        assert!(replaced
            .molecule
            .topology
            .as_ref()
            .unwrap()
            .valid_for(&replaced.molecule.atoms));
        for (i, m) in replaced.host_mapping.iter().enumerate() {
            if let Some(j) = m {
                if mol.topology.as_ref().unwrap().atoms[i].residue != Some(2) {
                    assert_eq!(mol.pos[i], replaced.molecule.pos[*j]);
                }
            }
        }
        assert!(replace(&mol, 2, &Residue::new(AminoAcid::Proline), true).is_err());
    }
    #[test]
    fn peptide_hydrogen_restoration_covers_every_template_and_preserves_heavy_atoms() {
        for t in &TEMPLATES {
            let recipe = Recipe {
                residues: vec![
                    Residue::new(AminoAcid::Glycine),
                    Residue {
                        amino: t.amino_acid,
                        side_chain: t.side_chain,
                    },
                    Residue::new(AminoAcid::Glycine),
                ],
                termini: t.termini,
                ..Default::default()
            };
            let mut mol = peptide::build(&recipe, None).unwrap().molecule;
            let count = mol.atoms.len();
            let q = mol.topology.as_ref().unwrap().charge();
            pack_sidechains(&mut mol, None).unwrap();
            strip_h(&mut mol);
            let before = mol.pos.clone();
            super::super::hydrogens::complete(&mut mol).unwrap();
            assert_eq!(mol.atoms.len(), count, "{}", t.id);
            assert_eq!(&mol.pos[..before.len()], &before);
            assert_eq!(mol.topology.as_ref().unwrap().charge(), q);
            assert_eq!(
                super::super::hydrogens::complete(&mut mol).unwrap().added,
                0
            );
            for &(a, b, _) in &mol.bonds {
                if mol.atoms[a] == "H" || mol.atoms[b] == "H" {
                    let d = mol.pos[a].distance(mol.pos[b]);
                    assert!((0.8..1.65).contains(&d), "{} H bond {d}", t.id);
                }
            }
        }
    }
    #[test]
    fn peptide_caps_have_expected_composition_and_no_dangling_hydrogens() {
        for termini in [
            super::super::amino_acids::TerminalState::Neutral,
            super::super::amino_acids::TerminalState::Zwitterionic,
        ] {
            let mol = peptide::build(
                &Recipe {
                    residues: parse_sequence("AG").unwrap(),
                    termini,
                    ..Default::default()
                },
                None,
            )
            .unwrap()
            .molecule;
            let capped = cap(
                &cap(&mol, 1, peptide::End::N).unwrap().molecule,
                2,
                peptide::End::C,
            )
            .unwrap()
            .molecule;
            assert_eq!(capped.topology.as_ref().unwrap().charge(), Some(0));
            let adj = peptide::adjacency(&capped);
            for (i, a) in capped.atoms.iter().enumerate() {
                if a == "H" {
                    assert_eq!(adj[i].len(), 1);
                }
            }
            assert!(cap(&capped, 1, peptide::End::N).is_err());
            assert!(cap(&capped, 2, peptide::End::C).is_err());
            let pdb = pdb(&capped).unwrap();
            assert!(pdb.contains("ACE"));
            assert!(pdb.contains("NME"));
            assert_eq!(
                pdb.lines()
                    .filter(|l| l.starts_with("HETATM") || l.starts_with("ATOM  "))
                    .count(),
                capped.atoms.len()
            );
            for l in pdb
                .lines()
                .filter(|l| l.starts_with("HETATM") || l.starts_with("ATOM  "))
            {
                assert_eq!(l.len(), 80);
            }
        }
    }
    #[test]
    fn peptide_display_switches_only_above_twenty_residues() {
        for n in [1, 20, 21] {
            let mol = chain(&"A".repeat(n));
            let t = mol.topology.as_ref().unwrap();
            let mut s = crate::settings::MolSettings::default();
            t.display_defaults(&mut s);
            assert!(!s.show_hbonds);
            assert_eq!(
                s.representation == crate::settings::RepresentationMode::BackboneTrace,
                n > 20
            );
            assert!(!t.backbone_atom(atom(&mol, 1, "CB").unwrap()));
            assert!(t.backbone_atom(atom(&mol, 1, "CA").unwrap()));
        }
    }
    #[test]
    fn peptide_disulfide_rejects_remote_pairs_without_mutation() {
        let mol = chain("CAAAAC");
        let before = mol.pos.clone();
        assert!(disulfide(&mol, 1, 6).is_err());
        assert_eq!(mol.pos, before);
    }
    #[test]
    fn peptide_positioned_disulfide_and_cap_hydrogens_are_idempotent() {
        let mut mol = peptide::beside(
            peptide::build(
                &Recipe {
                    residues: parse_sequence("C").unwrap(),
                    ..Default::default()
                },
                None,
            )
            .unwrap(),
            &chain("C"),
        )
        .molecule;
        for (id, origin) in [(1, Vec3::ZERO), (2, Vec3::X * 2.05)] {
            let from = [
                mol.pos[atom(&mol, id, "SG").unwrap()],
                mol.pos[atom(&mol, id, "CB").unwrap()],
                mol.pos[atom(&mol, id, "CA").unwrap()],
            ];
            let indices: Vec<_> = mol
                .topology
                .as_ref()
                .unwrap()
                .atoms
                .iter()
                .enumerate()
                .filter(|(_, a)| a.residue == Some(id))
                .map(|(i, _)| i)
                .collect();
            let mut pos: Vec<_> = indices.iter().map(|&i| mol.pos[i]).collect();
            peptide::transform(
                &mut pos,
                from,
                [origin, origin + Vec3::Y, origin + Vec3::Y + Vec3::X],
            );
            for (&i, &p) in indices.iter().zip(&pos) {
                mol.pos[i] = p;
            }
        }
        let linked = disulfide(&mol, 1, 2).unwrap();
        assert_eq!(linked.removed_host_atoms, 2);
        assert_eq!(linked.charge, Some(0));
        let mut result = linked.molecule;
        assert_eq!(restore_hydrogens(&mut result).unwrap(), 0);
        assert!(disulfide(&result, 1, 2).is_err());
        let mut capped = cap(
            &cap(&chain("AG"), 1, peptide::End::N).unwrap().molecule,
            2,
            peptide::End::C,
        )
        .unwrap()
        .molecule;
        let count = capped.atoms.len();
        strip_h(&mut capped);
        let result = super::super::hydrogens::complete(&mut capped).unwrap();
        assert_eq!(result.skipped, 0);
        assert_eq!(capped.atoms.len(), count);
        assert_eq!(
            super::super::hydrogens::complete(&mut capped)
                .unwrap()
                .added,
            0
        );
    }

    #[test]
    #[ignore = "20-residue DREIDING benchmark; run explicitly with --ignored --nocapture"]
    fn peptide_twenty_residue_hydrogen_and_dreiding_benchmark() {
        use crate::forcefield::dreiding::{
            objective::{cleanup_options, relax_angstrom},
            DreidingTopology,
        };
        let mut mol = chain("ACDEFGHIKLMNPQRSTVWY");
        assert_eq!(mol.topology.as_ref().unwrap().peptide_residue_count(), 20);
        pack_sidechains(&mut mol, None).unwrap();
        let total = mol.atoms.len();
        strip_h(&mut mol);
        let heavy = mol.atoms.len();
        let fixed = mol.pos.clone();
        let h = super::super::hydrogens::complete(&mut mol).unwrap();
        assert_eq!(mol.atoms.len(), total);
        assert_eq!(h.added, total - heavy);
        assert_eq!(h.skipped, 0);
        assert_eq!(&mol.pos[..heavy], &fixed);
        mol.recompute_bonds(1.2, 2.5);
        let adj = peptide::adjacency(&mol);
        for (i, a) in mol.atoms.iter().enumerate() {
            if a == "H" {
                assert_eq!(adj[i].len(), 1);
            }
        }
        let topology = DreidingTopology::build(&mol).expect("all standard amino acids must type");
        let start: Vec<_> = mol
            .pos
            .iter()
            .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
            .collect();
        let (clean, report) = relax_angstrom(&topology, &start, cleanup_options()).unwrap();
        println!("PEPTIDE CLEANUP: {heavy} heavy atoms, {} H restored, {total} total; cycles {}; E {:.6} -> {:.6} kcal/mol; force {:.6}; converged {}",h.added,report.cycles,report.initial_energy,report.final_energy,report.max_force,report.converged);
        assert!(report.final_energy.is_finite() && report.final_energy < report.initial_energy);
        assert!(clean.iter().all(|x| x.is_finite()));
        let mut options = cleanup_options();
        options.max_cycles = 1000;
        options.max_gradient = 1e-4;
        options.rms_gradient = 5e-5;
        let (optimized, report) = relax_angstrom(&topology, &clean, options).unwrap();
        println!(
            "PEPTIDE OPTIMIZATION: cycles {}; E {:.6} -> {:.6} kcal/mol; force {:.6}; converged {}",
            report.cycles,
            report.initial_energy,
            report.final_energy,
            report.max_force,
            report.converged
        );
        let (optimized, report) = if report.converged {
            (optimized, report)
        } else {
            let mut continuation = cleanup_options();
            continuation.max_cycles = 4000;
            continuation.max_gradient = 1e-4;
            continuation.rms_gradient = 5e-5;
            let result = relax_angstrom(&topology, &optimized, continuation).unwrap();
            println!("PEPTIDE CONTINUATION: cycles {}; E {:.6} -> {:.6} kcal/mol; force {:.6}; converged {}",result.1.cycles,result.1.initial_energy,result.1.final_energy,result.1.max_force,result.1.converged);
            result
        };
        assert!(
            report.final_energy.is_finite() && report.final_energy <= report.initial_energy + 1e-6
        );
        for &(a, b, _) in &mol.bonds {
            let p = |i| {
                Vec3::new(
                    optimized[3 * i] as f32,
                    optimized[3 * i + 1] as f32,
                    optimized[3 * i + 2] as f32,
                )
            };
            let ratio = p(a).distance(p(b)) / mol.pos[a].distance(mol.pos[b]);
            assert!(
                (0.65..1.45).contains(&ratio),
                "bond {a}-{b} changed excessively: {ratio}"
            );
        }
        for residue in &mol.topology.as_ref().unwrap().residues {
            let Some(t) = template(&residue.template) else {
                continue;
            };
            for stereo in t.stereocenters {
                let ids: Option<Vec<_>> = stereo
                    .neighbours
                    .iter()
                    .map(|&i| atom(&mol, residue.id, t.atoms[i].name).ok())
                    .collect();
                if let Some(ids) = ids {
                    let handedness =
                        |p: [Vec3; 4]| (p[0] - p[3]).dot((p[1] - p[3]).cross(p[2] - p[3]));
                    let before = handedness(std::array::from_fn(|i| mol.pos[ids[i]]));
                    let after = handedness(std::array::from_fn(|i| {
                        Vec3::new(
                            optimized[3 * ids[i]] as f32,
                            optimized[3 * ids[i] + 1] as f32,
                            optimized[3 * ids[i] + 2] as f32,
                        )
                    }));
                    assert!(
                        before * after > 0.0,
                        "Residue {} changed stereochemistry",
                        residue.id
                    );
                }
            }
        }
        let mut xyz = format!(
            "{total}\nACDEFGHIKLMNPQRSTVWY; DREIDING converged={} E={} kcal/mol\n",
            report.converged, report.final_energy
        );
        for (a, p) in mol.atoms.iter().zip(optimized.chunks_exact(3)) {
            xyz.push_str(&format!("{a} {:.8} {:.8} {:.8}\n", p[0], p[1], p[2]));
        }
        std::fs::write("/tmp/beavyr-peptide-20-dreiding.xyz", xyz).unwrap();
    }
}
