//! Dunbrack candidates fitted by rigid bond rotations to the existing side chain.
use super::{
    amino_acids::AminoAcid,
    peptide::{self, Build},
    peptide_conformation::{self as conf, Kind},
    peptide_edit::{self, atom, metadata},
};
use crate::molecule::Molecule;
use bevy::prelude::*;
use dunbrack::Residue as _;
use std::collections::{BTreeMap, HashSet};

#[derive(Clone)]
pub struct Candidate {
    pub positions: Vec<Vec3>,
    pub probability: Option<f32>,
    pub chi: Vec<f32>,
    pub score: f32,
    pub current: bool,
}
#[derive(Clone)]
pub struct Choices {
    pub candidates: Vec<Candidate>,
    pub marginal: bool,
    pub residue: u64,
    source: Vec<Vec3>,
    topology: Option<super::topology::Topology>,
}
fn definitions(a: AminoAcid) -> Vec<[&'static str; 4]> {
    use AminoAcid::*;
    let path: &[&str] = match a {
        Valine => &["N", "CA", "CB", "CG1"],
        Isoleucine => &["N", "CA", "CB", "CG1", "CD1"],
        Leucine => &["N", "CA", "CB", "CG", "CD1"],
        Serine => &["N", "CA", "CB", "OG"],
        Threonine => &["N", "CA", "CB", "OG1"],
        Cysteine => &["N", "CA", "CB", "SG"],
        AsparticAcid | Asparagine => &["N", "CA", "CB", "CG", "OD1"],
        GlutamicAcid | Glutamine => &["N", "CA", "CB", "CG", "CD", "OE1"],
        Histidine => &["N", "CA", "CB", "CG", "ND1"],
        Phenylalanine | Tyrosine | Tryptophan => &["N", "CA", "CB", "CG", "CD1"],
        Methionine => &["N", "CA", "CB", "CG", "SD", "CE"],
        Lysine => &["N", "CA", "CB", "CG", "CD", "CE", "NZ"],
        Arginine => &["N", "CA", "CB", "CG", "CD", "NE", "CZ"],
        _ => &[],
    };
    path.windows(4).map(|p| [p[0], p[1], p[2], p[3]]).collect()
}
fn raw(a: AminoAcid, phi: f32, psi: f32) -> Vec<(Vec<u8>, f32, Vec<f32>)> {
    macro_rules! query {
        ($t:ident) => {
            dunbrack::$t::rotamers(phi, psi)
                .map(|r| (r.r.to_vec(), r.prob, r.chi_mean.to_vec()))
                .collect()
        };
    }
    use AminoAcid::*;
    match a {
        Valine => query!(Val),
        Isoleucine => query!(Ile),
        Leucine => query!(Leu),
        Serine => query!(Ser),
        Threonine => query!(Thr),
        Cysteine => query!(Cyh),
        AsparticAcid => query!(Asp),
        Asparagine => query!(Asn),
        GlutamicAcid => query!(Glu),
        Glutamine => query!(Gln),
        Histidine => query!(His),
        Phenylalanine => query!(Phe),
        Tyrosine => query!(Tyr),
        Tryptophan => query!(Trp),
        Methionine => query!(Met),
        Lysine => query!(Lys),
        Arginine => query!(Arg),
        _ => vec![],
    }
}
fn library(a: AminoAcid, phi: Option<f32>, psi: Option<f32>) -> Vec<(f32, Vec<f32>)> {
    // At termini average the unknown backbone dimension over the periodic grid,
    // rather than silently pretending that every terminus is alpha-helical.
    let grid = || {
        (0..36)
            .map(|i| -180.0 + i as f32 * 10.0)
            .collect::<Vec<_>>()
    };
    let phis = phi.map(|v| vec![v]).unwrap_or_else(grid);
    let psis = psi.map(|v| vec![v]).unwrap_or_else(grid);
    let n = (phis.len() * psis.len()) as f32;
    let mut sums: BTreeMap<Vec<u8>, (f32, Vec<(f32, f32)>)> = BTreeMap::new();
    for &p in &phis {
        for &q in &psis {
            for (key, prob, chi) in raw(a, p, q) {
                let entry = sums
                    .entry(key)
                    .or_insert_with(|| (0.0, vec![(0.0, 0.0); chi.len()]));
                entry.0 += prob / n;
                for (i, v) in chi.iter().enumerate() {
                    let (s, c) = v.to_radians().sin_cos();
                    entry.1[i].0 += prob * s;
                    entry.1[i].1 += prob * c;
                }
            }
        }
    }
    sums.into_values()
        .map(|(p, v)| {
            (
                p,
                v.into_iter()
                    .map(|(s, c)| s.atan2(c).to_degrees())
                    .collect(),
            )
        })
        .collect()
}
// Count each affected nonbonded pair once, including clashes between two
// moving side-chain atoms after several chi rotations.
fn overlap_score(mol: &Molecule, moving: &[usize], adj: &[Vec<usize>]) -> f32 {
    let moving: HashSet<_> = moving.iter().copied().collect();
    let radius = |s: &str| match s {
        "H" => 1.1,
        "N" => 1.55,
        "O" => 1.52,
        "S" | "Se" => 1.8,
        _ => 1.7,
    };
    let mut score = 0.0;
    let mut ordered: Vec<_> = moving.iter().copied().collect();
    ordered.sort_unstable();
    for i in ordered {
        for j in 0..mol.atoms.len() {
            if i == j
                || (moving.contains(&j) && j < i)
                || adj[i].contains(&j)
                || adj[i].iter().any(|&k| adj[k].contains(&j))
            {
                continue;
            }
            let overlap = (0.8 * (radius(&mol.atoms[i]) + radius(&mol.atoms[j]))
                - mol.pos[i].distance(mol.pos[j]))
            .max(0.0);
            score += overlap * overlap;
        }
    }
    score
}
pub fn choices(mol: &Molecule, id: u64) -> Result<Choices, String> {
    let t = metadata(mol)?;
    let r = t
        .residues
        .iter()
        .find(|r| r.id == id)
        .ok_or("Select a residue.")?;
    let a = peptide_edit::template(&r.template)
        .ok_or("This residue has no rotamer library.")?
        .amino_acid;
    let defs = definitions(a);
    if defs.is_empty() {
        return Err("No independent library rotamers for this residue. Gly/Ala have no χ; proline requires ring fitting; Sec/Pyl are not covered by this library.".into());
    }
    let mut ids = vec![];
    let mut moving = HashSet::new();
    for names in defs {
        let mut q = [0; 4];
        for (i, name) in names.iter().enumerate() {
            q[i] = atom(mol, id, name)?;
        }
        let side = conf::moving_side(mol, q[1], q[2])?;
        if side.iter().any(|&i| t.atoms[i].residue != Some(id)) {
            return Err("A side-chain cross-link blocks independent rotamers.".into());
        }
        moving.extend(side);
        ids.push(q);
    }
    let moving: Vec<_> = moving.into_iter().collect();
    let adj = peptide::adjacency(mol);
    let current = Candidate {
        positions: mol.pos.clone(),
        probability: None,
        chi: ids
            .iter()
            .map(|&[a, b, c, d]| peptide::torsion(mol.pos[a], mol.pos[b], mol.pos[c], mol.pos[d]))
            .collect(),
        score: overlap_score(mol, &moving, &adj),
        current: true,
    };
    let phi = conf::measured(mol, id, Kind::Phi);
    let psi = conf::measured(mol, id, Kind::Psi);
    let mut candidates = vec![];
    for (probability, chi) in library(a, phi, psi) {
        if !probability.is_finite() || probability < 1e-5 {
            continue;
        }
        let mut candidate = mol.clone();
        for (&q, &v) in ids.iter().zip(&chi) {
            conf::set_torsion(&mut candidate, q, v)?;
        }
        let score = overlap_score(&candidate, &moving, &adj);
        candidates.push(Candidate {
            positions: candidate.pos,
            probability: Some(probability),
            chi,
            score,
            current: false,
        });
        candidates.sort_by(|a, b| {
            a.score.total_cmp(&b.score).then_with(|| {
                b.probability
                    .unwrap_or(0.0)
                    .total_cmp(&a.probability.unwrap_or(0.0))
            })
        });
        candidates.truncate(12);
    }
    candidates.insert(0, current);
    Ok(Choices {
        candidates,
        marginal: phi.is_none() || psi.is_none(),
        residue: id,
        source: mol.pos.clone(),
        topology: mol.topology.clone(),
    })
}
pub fn preview(mol: &Molecule, choices: &Choices, index: usize) -> Result<Build, String> {
    let c = choices.candidates.get(index).ok_or("Choose a rotamer.")?;
    if choices.source != mol.pos
        || choices.topology != mol.topology
        || c.positions.len() != mol.pos.len()
    {
        return Err("Structure changed; generate rotamers again.".into());
    }
    let mut result = mol.clone();
    result.pos = c.positions.clone();
    conf::build_edit(
        mol,
        result,
        choices.residue,
        if c.current {
            "Current side chain".into()
        } else {
            format!("Rotamer · local overlap {:.3}", c.score)
        },
    )
}
pub fn repack_neighbors(mol: &Molecule, id: u64, radius: f32) -> Result<Build, String> {
    if !radius.is_finite() || !(2.0..=12.0).contains(&radius) {
        return Err("Choose a neighborhood radius from 2 to 12 Å.".into());
    }
    let t = metadata(mol)?;
    let target: Vec<_> = t
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.residue == Some(id))
        .map(|(i, _)| mol.pos[i])
        .collect();
    if target.is_empty() {
        return Err("Select a residue.".into());
    }
    let neighbors: Vec<_> = t
        .residues
        .iter()
        .filter(|r| {
            r.id != id
                && t.atoms.iter().enumerate().any(|(i, a)| {
                    a.residue == Some(r.id)
                        && mol.atoms[i] != "H"
                        && target
                            .iter()
                            .any(|p| p.distance_squared(mol.pos[i]) <= radius * radius)
                })
        })
        .map(|r| r.id)
        .collect();
    let mut result = mol.clone();
    let mut changed = 0;
    let mut skipped = 0;
    for r in neighbors {
        match choices(&result, r) {
            Ok(options) => {
                let best = options
                    .candidates
                    .iter()
                    .skip(1)
                    .find(|c| c.score + 1e-5 < options.candidates[0].score);
                if let Some(best) = best {
                    result.pos = best.positions.clone();
                    changed += 1;
                }
            }
            Err(_) => skipped += 1,
        }
    }
    conf::build_edit(
        mol,
        result,
        id,
        format!("Repacked {changed} neighboring side chains; {skipped} unavailable/locked"),
    )
}
