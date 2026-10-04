//! Peptide assembly from the embedded, stereochemically checked monomers.
//! Produces starting geometries, not folded or reoptimized proteins.
use bevy::prelude::*;
use std::collections::HashSet;

use super::amino_acids::{AminoAcid, AminoAcidTemplate, SideChainState, TerminalState};
use super::topology::{AtomInfo, ResidueInfo, Topology};
use crate::molecule::{covalent_radius_angstrom, Molecule};

pub const MAX_RESIDUES: usize = 256;
const PEPTIDE_BOND: f32 = 1.34;
const CA_C_N: f32 = 114.0;
const C_N_CA: f32 = 123.0;

#[derive(Clone, Debug, PartialEq)]
pub struct Residue {
    pub amino: AminoAcid,
    pub side_chain: SideChainState,
}

impl Residue {
    pub fn new(amino: AminoAcid) -> Self {
        Self {
            amino,
            side_chain: if amino == AminoAcid::Histidine {
                SideChainState::HistidineEpsilon
            } else {
                SideChainState::Neutral
            },
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum End {
    C,
    N,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Shape {
    Extended,
    Alpha,
    Custom,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Recipe {
    /// Always in chemical N-to-C order, irrespective of the tile insertion end.
    pub residues: Vec<Residue>,
    pub termini: TerminalState,
    pub shape: Shape,
    pub phi: f32,
    pub psi: f32,
}

impl Default for Recipe {
    fn default() -> Self {
        Self {
            residues: Vec::new(),
            termini: TerminalState::Neutral,
            shape: Shape::Extended,
            phi: -135.0,
            psi: 135.0,
        }
    }
}

impl Recipe {
    pub fn angles(&self) -> (f32, f32) {
        match self.shape {
            Shape::Extended => (-135.0, 135.0),
            Shape::Alpha => (-57.0, -47.0),
            Shape::Custom => (self.phi, self.psi),
        }
    }
    pub fn sequence(&self) -> String {
        self.residues.iter().map(|r| r.amino.one_letter()).collect()
    }
}

pub fn parse_sequence(text: &str) -> Result<Vec<Residue>, String> {
    let trimmed = text.trim();
    if trimmed.is_empty() {
        return Err("Enter a one-letter sequence or separated three-letter names.".into());
    }
    if trimmed
        .lines()
        .filter(|line| line.trim_start().starts_with('>'))
        .count()
        > 1
    {
        return Err("Paste one FASTA sequence at a time.".into());
    }
    // Three-letter names are an explicit separated format; a bare ALA is the
    // one-letter sequence Ala-Leu-Ala, never guessed to mean something else.
    let three_letter = trimmed.contains(|c: char| c == '-' || c == ',' || c.is_whitespace())
        && trimmed
            .split(|c: char| c == '-' || c == ',' || c.is_whitespace())
            .filter(|s| !s.is_empty())
            .all(|s| AminoAcid::from_three_letter(s).is_some());
    let amino = if three_letter {
        trimmed
            .split(|c: char| c == '-' || c == ',' || c.is_whitespace())
            .filter(|s| !s.is_empty())
            .map(|s| AminoAcid::from_three_letter(s).ok_or_else(|| format!("Unknown residue: {s}")))
            .collect::<Result<Vec<_>, _>>()?
    } else {
        let mut result = Vec::new();
        for line in trimmed
            .lines()
            .filter(|line| !line.trim_start().starts_with('>'))
        {
            for c in line.chars().filter(|c| !c.is_whitespace()) {
                result.push(AminoAcid::from_one_letter(c)
                    .ok_or_else(|| format!("Unknown residue '{c}'. Use one-letter codes, including U and O, or names such as Ala-Gly-Ser."))?);
            }
        }
        result
    };
    if amino.is_empty() || amino.len() > MAX_RESIDUES {
        return Err(format!("Choose 1–{MAX_RESIDUES} residues."));
    }
    Ok(amino.into_iter().map(Residue::new).collect())
}

#[derive(Clone)]
struct Part {
    template: &'static AminoAcidTemplate,
    pos: Vec<Vec3>,
    removed: HashSet<usize>,
    charges: Vec<i8>,
}

impl Part {
    fn new(residue: &Residue, termini: TerminalState) -> Result<Self, String> {
        let template = residue
            .amino
            .template(termini, residue.side_chain)
            .ok_or("This residue does not have the requested protonation state.")?;
        Ok(Self {
            template,
            pos: template
                .atoms
                .iter()
                .map(|a| Vec3::from_array(a.position_angstrom.map(|x| x as f32)))
                .collect(),
            removed: HashSet::new(),
            charges: template.atoms.iter().map(|a| a.formal_charge).collect(),
        })
    }
    fn n(&self) -> usize {
        self.template.backbone.nitrogen
    }
    fn ca(&self) -> usize {
        self.template.backbone.alpha_carbon
    }
    fn c(&self) -> usize {
        self.template.backbone.carbonyl_carbon
    }
    fn o(&self) -> usize {
        self.template.backbone.carbonyl_oxygen
    }
    fn neighbors(&self, index: usize) -> Vec<usize> {
        self.template
            .bonds
            .iter()
            .filter_map(|b| {
                if b.atoms[0] == index {
                    Some(b.atoms[1])
                } else if b.atoms[1] == index {
                    Some(b.atoms[0])
                } else {
                    None
                }
            })
            .collect()
    }
    fn acid_leaving_group(&mut self) {
        let oxt = self.template.backbone.terminal_oxygen;
        self.removed.insert(oxt);
        for i in self.neighbors(oxt) {
            if self.template.atoms[i].element == "H" {
                self.removed.insert(i);
            }
        }
    }
    fn amide(&mut self, previous_c: Vec3) {
        let n = self.n();
        let keep = usize::from(self.template.amino_acid != AminoAcid::Proline);
        let hydrogens: Vec<_> = self
            .neighbors(n)
            .into_iter()
            .filter(|&i| self.template.atoms[i].element == "H")
            .collect();
        for &i in hydrogens.iter().skip(keep) {
            self.removed.insert(i);
        }
        self.charges[n] = 0;
        if keep == 1 {
            let direction = -((previous_c - self.pos[n]).normalize()
                + (self.pos[self.ca()] - self.pos[n]).normalize())
            .normalize();
            self.pos[hydrogens[0]] = self.pos[n] + 1.01 * direction;
        }
    }
    fn carbonyl(&mut self, next_n: Vec3) {
        let c = self.pos[self.c()];
        let direction =
            -((self.pos[self.ca()] - c).normalize() + (next_n - c).normalize()).normalize();
        let o = self.o();
        self.pos[o] = c + 1.23 * direction;
    }
    fn incoming_geometry(&self, requested_phi: f32) -> (f32, f32) {
        if self.template.amino_acid != AminoAcid::Proline {
            return (requested_phi, C_N_CA);
        }
        // Proline's N–CA bond is in a ring: never twist it independently.
        let n = self.pos[self.n()];
        let ca = self.pos[self.ca()];
        let cd = self.pos[self.template.atom_index("CD").unwrap()];
        let previous_c =
            n - ((ca - n).normalize() + (cd - n).normalize()).normalize() * PEPTIDE_BOND;
        (
            torsion(previous_c, n, ca, self.pos[self.c()]),
            angle(previous_c, n, ca),
        )
    }
}

fn angle(a: Vec3, b: Vec3, c: Vec3) -> f32 {
    (a - b)
        .normalize()
        .dot((c - b).normalize())
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees()
}

pub(super) fn torsion(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> f32 {
    let axis = (c - b).normalize();
    let v = a - b - (a - b).dot(axis) * axis;
    let w = d - c - (d - c).dot(axis) * axis;
    axis.cross(v).dot(w).atan2(v.dot(w)).to_degrees()
}

/// Internal-coordinate placement: angle B-C-D and torsion A-B-C-D.
pub(super) fn place(a: Vec3, b: Vec3, c: Vec3, length: f32, theta: f32, dihedral: f32) -> Vec3 {
    let axis = (c - b).normalize();
    let v = (a - b - (a - b).dot(axis) * axis).normalize();
    let w = axis.cross(v);
    let (st, ct) = theta.to_radians().sin_cos();
    let (sd, cd) = dihedral.to_radians().sin_cos();
    c + length * (-ct * axis + st * (cd * v + sd * w))
}

fn frame(p: [Vec3; 3]) -> Mat3 {
    let x = (p[1] - p[0]).normalize();
    let z = x.cross(p[2] - p[0]).normalize();
    Mat3::from_cols(x, z.cross(x), z)
}

pub(super) fn transform(positions: &mut [Vec3], from: [Vec3; 3], to: [Vec3; 3]) {
    let rotation = frame(to) * frame(from).transpose();
    for p in positions {
        *p = to[0] + rotation * (*p - from[0]);
    }
}

fn assemble(recipe: &Recipe) -> Result<Vec<Part>, String> {
    if recipe.residues.is_empty() || recipe.residues.len() > MAX_RESIDUES {
        return Err(format!("Choose 1–{MAX_RESIDUES} residues."));
    }
    let (phi, psi) = recipe.angles();
    if !phi.is_finite() || !psi.is_finite() {
        return Err("Backbone angles must be finite.".into());
    }
    let mut parts = Vec::<Part>::new();
    for residue in &recipe.residues {
        let mut part = Part::new(residue, recipe.termini)?;
        if let Some(previous) = parts.last_mut() {
            let prev_n = previous.pos[previous.n()];
            let prev_ca = previous.pos[previous.ca()];
            let prev_c = previous.pos[previous.c()];
            let source = [part.pos[part.n()], part.pos[part.ca()], part.pos[part.c()]];
            let (effective_phi, n_angle) = part.incoming_geometry(phi);
            let n = place(prev_n, prev_ca, prev_c, PEPTIDE_BOND, CA_C_N, psi);
            let ca = place(
                prev_ca,
                prev_c,
                n,
                source[0].distance(source[1]),
                n_angle,
                180.0,
            );
            let c = place(
                prev_c,
                n,
                ca,
                source[1].distance(source[2]),
                angle(source[0], source[1], source[2]),
                effective_phi,
            );
            transform(&mut part.pos, source, [n, ca, c]);
            previous.acid_leaving_group();
            previous.carbonyl(n);
            part.amide(prev_c);
        }
        parts.push(part);
    }
    Ok(parts)
}

#[derive(Clone, Debug)]
pub struct Attachment {
    pub end: End,
    pub atom: usize,
    pub neighbor: usize,
    pub leaving: Vec<usize>,
    // Direction of the new bond at the host atom, determined before deletion.
    pub direction: Vec3,
}

pub(super) fn adjacency(mol: &Molecule) -> Vec<Vec<usize>> {
    let mut result = vec![Vec::new(); mol.atoms.len()];
    for &(a, b, _) in &mol.bonds {
        if a < result.len() && b < result.len() {
            result[a].push(b);
            result[b].push(a);
        }
    }
    result
}

/// Only free carboxyls and saturated primary/secondary amines can condense.
/// Arbitrary selected atoms and already linked peptide nitrogens are rejected.
pub fn attachment(mol: &Molecule, atom: usize, end: End) -> Result<Attachment, String> {
    if mol.atoms.len() != mol.pos.len() || atom >= mol.atoms.len() {
        return Err("Select a terminal atom in the structure.".into());
    }
    let adj = adjacency(mol);
    let heavy: Vec<_> = adj[atom]
        .iter()
        .copied()
        .filter(|&i| mol.atoms[i] != "H")
        .collect();
    let hydrogens: Vec<_> = adj[atom]
        .iter()
        .copied()
        .filter(|&i| mol.atoms[i] == "H")
        .collect();
    match end {
        End::C => {
            if mol.atoms[atom] != "C" || heavy.len() != 3 || !hydrogens.is_empty() {
                return Err(
                    "For C-terminal growth, select the carbon of a free COOH or COO⁻ group.".into(),
                );
            }
            let oxygens: Vec<_> = heavy
                .iter()
                .copied()
                .filter(|&i| mol.atoms[i] == "O")
                .collect();
            let carbon = heavy.iter().copied().find(|&i| mol.atoms[i] == "C");
            if oxygens.len() != 2 || carbon.is_none() {
                return Err("The selected carbon is not a free carboxyl carbon.".into());
            }
            if oxygens
                .iter()
                .any(|&o| adj[o].iter().any(|&j| j != atom && mol.atoms[j] != "H"))
            {
                return Err(
                    "This carboxyl group is esterified; choose a free COOH or COO⁻ terminus."
                        .into(),
                );
            }
            let protonated: Vec<_> = oxygens
                .iter()
                .copied()
                .filter(|&o| adj[o].iter().any(|&i| mol.atoms[i] == "H"))
                .collect();
            if protonated.len() > 1 {
                return Err("Ambiguous carboxyl group.".into());
            }
            let named_oxt = mol
                .topology
                .as_ref()
                .filter(|t| t.valid_for(&mol.atoms))
                .and_then(|t| oxygens.iter().copied().find(|&i| t.atoms[i].name == "OXT"));
            let leaving_o = named_oxt
                .or_else(|| protonated.first().copied())
                .unwrap_or_else(|| {
                    *oxygens
                        .iter()
                        .max_by(|&&a, &&b| {
                            mol.pos[a]
                                .distance(mol.pos[atom])
                                .total_cmp(&mol.pos[b].distance(mol.pos[atom]))
                        })
                        .unwrap()
                });
            let retained_o = *oxygens.iter().find(|&&i| i != leaving_o).unwrap();
            let neighbor = carbon.unwrap();
            let direction = -((mol.pos[neighbor] - mol.pos[atom]).normalize()
                + (mol.pos[retained_o] - mol.pos[atom]).normalize())
            .normalize();
            let mut leaving = vec![leaving_o];
            leaving.extend(
                adj[leaving_o]
                    .iter()
                    .copied()
                    .filter(|&i| mol.atoms[i] == "H"),
            );
            Ok(Attachment {
                end,
                atom,
                neighbor,
                leaving,
                direction,
            })
        }
        End::N => {
            if mol.atoms[atom] != "N"
                || !(1..=2).contains(&heavy.len())
                || heavy.iter().any(|&i| mol.atoms[i] != "C")
                || hydrogens.len() < 3 - heavy.len()
                || hydrogens.len() > 4 - heavy.len()
            {
                return Err("For N-terminal growth, select a free primary or secondary amine N with explicit hydrogens.".into());
            }
            if heavy.iter().any(|&c| {
                mol.pos[c].distance(mol.pos[atom]) < 1.36
                    || adj[c]
                        .iter()
                        .any(|&o| mol.atoms[o] == "O" && mol.pos[o].distance(mol.pos[c]) < 1.32)
            }) {
                return Err("This nitrogen is already acylated or unsaturated; choose a free amino terminus.".into());
            }
            // Prefer an alpha carbon in cyclic amines, without requiring residue metadata.
            let neighbor = heavy
                .iter()
                .copied()
                .find(|&c| {
                    adj[c].iter().any(|&j| {
                        mol.atoms[j] == "C" && adj[j].iter().any(|&o| mol.atoms[o] == "O")
                    })
                })
                .unwrap_or(heavy[0]);
            let keep_h = 2 - heavy.len();
            let other = if keep_h == 1 {
                hydrogens[0]
            } else {
                *heavy.iter().find(|&&i| i != neighbor).unwrap()
            };
            let direction = -((mol.pos[neighbor] - mol.pos[atom]).normalize()
                + (mol.pos[other] - mol.pos[atom]).normalize())
            .normalize();
            Ok(Attachment {
                end,
                atom,
                neighbor,
                leaving: hydrogens.into_iter().skip(keep_h).collect(),
                direction,
            })
        }
    }
}

#[derive(Clone)]
pub struct Build {
    pub molecule: Molecule,
    pub moving: Vec<usize>,
    pub picks: [usize; 3],
    pub charge: Option<i32>,
    pub clashes: usize,
    pub removed_host_atoms: usize,
    pub sequence: String,
    /// Original host index -> retained index. Used to preserve freezes/isolation.
    pub host_mapping: Vec<Option<usize>>,
}

pub fn build(recipe: &Recipe, host: Option<(&Molecule, &Attachment)>) -> Result<Build, String> {
    let mut parts = assemble(recipe)?;
    if let Some((mol, anchor)) = host {
        // Revalidate before indexing: callers can retain a selection across edits.
        let fresh = attachment(mol, anchor.atom, anchor.end)?;
        if fresh.leaving != anchor.leaving || fresh.neighbor != anchor.neighbor {
            return Err("The attachment site changed. Pick it again.".into());
        }
        let host_atom = mol.pos[fresh.atom];
        let host_neighbor = mol.pos[fresh.neighbor];
        let connector = host_atom + fresh.direction * PEPTIDE_BOND;
        let (from, to) = match anchor.end {
            End::C => {
                let first = &parts[0];
                let n = first.pos[first.n()];
                let ca = first.pos[first.ca()];
                let c = first.pos[first.c()];
                let (phi, theta) = first.incoming_geometry(recipe.angles().0);
                let virtual_c = place(c, ca, n, PEPTIDE_BOND, theta, phi);
                let new_ca = place(
                    host_neighbor,
                    host_atom,
                    connector,
                    n.distance(ca),
                    theta,
                    180.0,
                );
                ([n, virtual_c, ca], [connector, host_atom, new_ca])
            }
            End::N => {
                let last = parts.last().unwrap();
                let c = last.pos[last.c()];
                let ca = last.pos[last.ca()];
                let o = last.pos[last.o()];
                let virtual_n =
                    c - ((ca - c).normalize() + (o - c).normalize()).normalize() * PEPTIDE_BOND;
                let theta = angle(ca, c, virtual_n);
                let new_ca = place(
                    host_neighbor,
                    host_atom,
                    connector,
                    c.distance(ca),
                    theta,
                    180.0,
                );
                ([c, virtual_n, ca], [connector, host_atom, new_ca])
            }
        };
        for part in &mut parts {
            transform(&mut part.pos, from, to);
        }
        match anchor.end {
            End::C => parts[0].amide(host_atom),
            End::N => {
                let last = parts.last_mut().unwrap();
                last.acid_leaving_group();
                last.carbonyl(host_atom);
            }
        }
    }
    let mut molecule = Molecule::empty();
    let mut host_mapping = Vec::new();
    let mut edges = Vec::new();
    if let Some((mol, anchor)) = host {
        host_mapping.resize(mol.atoms.len(), None);
        for (i, entry) in host_mapping.iter_mut().enumerate() {
            if !anchor.leaving.contains(&i) {
                *entry = Some(molecule.atoms.len());
                molecule.atoms.push(mol.atoms[i].clone());
                molecule.pos.push(mol.pos[i]);
            }
        }
        for &(a, b, _) in &mol.bonds {
            if let (Some(a), Some(b)) = (host_mapping[a], host_mapping[b]) {
                edges.push((a, b));
            }
        }
    }
    let mut topology = if let Some((mol, anchor)) = host {
        let mut topology = Topology::from_molecule(mol);
        topology.retain(&host_mapping);
        if anchor.end == End::N {
            topology.atoms[host_mapping[anchor.atom].unwrap()].charge = Some(0);
        }
        topology
    } else {
        Topology::default()
    };
    let mut next_atom = topology.atoms.iter().map(|a| a.id).max().unwrap_or(0) + 1;
    let mut next_residue = topology.residues.iter().map(|r| r.id).max().unwrap_or(0) + 1;
    let first_residue_id = next_residue;
    let chain = host
        .and_then(|(_, anchor)| host_mapping[anchor.atom])
        .and_then(|i| topology.atoms[i].residue)
        .and_then(|id| {
            topology
                .residues
                .iter()
                .find(|r| r.id == id)
                .map(|r| r.chain)
        })
        .unwrap_or_else(|| topology.residues.iter().map(|r| r.chain).max().unwrap_or(0) + 1);
    let start = molecule.atoms.len();
    let mut mappings = Vec::new();

    for part in &parts {
        let residue_id = next_residue;
        next_residue += 1;
        topology.residues.push(ResidueInfo {
            origin: None,
            id: residue_id,
            chain,
            template: part.template.id.into(),
        });
        let mut mapping = vec![None; part.pos.len()];
        for (i, entry) in mapping.iter_mut().enumerate() {
            if !part.removed.contains(&i) {
                *entry = Some(molecule.atoms.len());
                molecule.atoms.push(part.template.atoms[i].element.into());
                molecule.pos.push(part.pos[i]);
                topology.atoms.push(AtomInfo {
                    id: next_atom,
                    residue: Some(residue_id),
                    name: part.template.atoms[i].name.into(),
                    charge: Some(part.charges[i]),
                });
                next_atom += 1;
            }
        }
        for bond in part.template.bonds {
            if let (Some(a), Some(b)) = (mapping[bond.atoms[0]], mapping[bond.atoms[1]]) {
                edges.push((a, b));
                let order = match bond.order {
                    super::amino_acids::BondOrder::Single => 1,
                    super::amino_acids::BondOrder::Double => 2,
                    super::amino_acids::BondOrder::Triple => 3,
                    super::amino_acids::BondOrder::Aromatic => 4,
                };
                topology.bonds.push((a, b, order));
            }
        }
        mappings.push(mapping);
    }
    for i in 1..parts.len() {
        edges.push((
            mappings[i - 1][parts[i - 1].c()].unwrap(),
            mappings[i][parts[i].n()].unwrap(),
        ));
    }
    let first_n = mappings[0][parts[0].n()].unwrap();
    let first_ca = mappings[0][parts[0].ca()].unwrap();
    let last = parts.len() - 1;
    let picks = if let Some((_, anchor)) = host {
        let host_index = host_mapping[anchor.atom].unwrap();
        let (connector, side) = match anchor.end {
            End::C => (first_n, first_ca),
            End::N => (
                mappings[last][parts[last].c()].unwrap(),
                mappings[last][parts[last].ca()].unwrap(),
            ),
        };
        edges.push((host_index, connector));
        [host_index, connector, side]
    } else {
        [first_n, first_ca, mappings[0][parts[0].c()].unwrap()]
    };
    if molecule.pos.iter().any(|p| !p.is_finite()) {
        return Err("Degenerate attachment geometry; select a different site.".into());
    }
    molecule.bonds = edges
        .into_iter()
        .map(|(a, b)| (a, b, molecule.pos[a].distance(molecule.pos[b])))
        .collect();
    for &(a, b, _) in &molecule.bonds {
        if !topology
            .bonds
            .iter()
            .any(|&(x, y, _)| (x == a && y == b) || (x == b && y == a))
        {
            topology.bonds.push((a, b, 1));
        }
    }
    if let Some((_, anchor)) = host {
        if anchor.end == End::N {
            let host_residue = topology.atoms[host_mapping[anchor.atom].unwrap()].residue;
            if let Some(at) = topology
                .residues
                .iter()
                .position(|r| Some(r.id) == host_residue)
            {
                let added: Vec<_> = topology
                    .residues
                    .iter()
                    .filter(|r| r.id >= first_residue_id)
                    .cloned()
                    .collect();
                topology.residues.retain(|r| r.id < first_residue_id);
                topology.residues.splice(at..at, added);
            }
        }
    }
    topology.elements = molecule.atoms.clone();
    let charge = topology.charge();
    molecule.topology = Some(topology);
    let clashes = clashes(&molecule, start);
    Ok(Build {
        moving: (start..molecule.atoms.len()).collect(),
        molecule,
        picks,
        charge,
        clashes,
        removed_host_atoms: host.map_or(0, |(_, a)| a.leaving.len()),
        sequence: recipe.sequence(),
        host_mapping,
    })
}

pub(super) fn clashes(mol: &Molecule, first_new: usize) -> usize {
    let adj = adjacency(mol);
    let radii: Vec<_> = mol
        .atoms
        .iter()
        .map(|a| covalent_radius_angstrom(a))
        .collect();
    let mut count = 0;
    for i in first_new..mol.atoms.len() {
        for j in 0..i {
            if adj[i].contains(&j) || adj[i].iter().any(|&k| adj[k].contains(&j)) {
                continue;
            }
            if mol.pos[i].distance_squared(mol.pos[j]) < (0.85 * (radii[i] + radii[j])).powi(2) {
                count += 1;
            }
        }
    }
    count
}

/// A new peptide can coexist with the displayed structure without replacing it.
pub fn beside(mut build: Build, host: &Molecule) -> Build {
    if host.atoms.is_empty() {
        return build;
    }
    let offset = host.atoms.len();
    let host_max = host
        .pos
        .iter()
        .map(|p| p.x)
        .fold(f32::NEG_INFINITY, f32::max);
    let new_min = build
        .molecule
        .pos
        .iter()
        .map(|p| p.x)
        .fold(f32::INFINITY, f32::min);
    let shift = Vec3::X * (host_max - new_min + 5.0);
    let mut topology = Topology::from_molecule(host);
    let new = build.molecule.topology.as_ref().unwrap();
    let atom_offset = topology.atoms.iter().map(|a| a.id).max().unwrap_or(0);
    let residue_offset = topology.residues.iter().map(|r| r.id).max().unwrap_or(0);
    let chain_offset = topology.residues.iter().map(|r| r.chain).max().unwrap_or(0);
    topology
        .atoms
        .extend(new.atoms.iter().cloned().map(|mut a| {
            a.id += atom_offset;
            a.residue = a.residue.map(|r| r + residue_offset);
            a
        }));
    topology
        .residues
        .extend(new.residues.iter().cloned().map(|mut r| {
            r.id += residue_offset;
            r.chain += chain_offset;
            r
        }));
    topology.bonds.extend(
        new.bonds
            .iter()
            .map(|&(a, b, o)| (a + offset, b + offset, o)),
    );
    let mut molecule = host.clone();
    molecule.atoms.extend(build.molecule.atoms);
    molecule
        .pos
        .extend(build.molecule.pos.into_iter().map(|p| p + shift));
    molecule.bonds.extend(
        build
            .molecule
            .bonds
            .into_iter()
            .map(|(a, b, d)| (a + offset, b + offset, d)),
    );
    topology.elements = molecule.atoms.clone();
    build.charge = topology.charge();
    molecule.topology = Some(topology);
    build.molecule = molecule;
    build.moving.iter_mut().for_each(|i| *i += offset);
    build.picks.iter_mut().for_each(|i| *i += offset);
    build.host_mapping = (0..offset).map(Some).collect();

    build
}

#[cfg(test)]
#[path = "peptide_tests.rs"]
mod tests;
