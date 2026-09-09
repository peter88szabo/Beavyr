use std::collections::{HashSet, VecDeque};

use anyhow::{bail, ensure, Result};

use super::rings::{self, Pucker};
use crate::optimizer::internal_coords::{analyze_structure, ConnectivityModel};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TorsionKind {
    Rotatable,
    CisTrans,
    /// A five- or six-membered ring's pucker.
    ///
    /// Not a torsion at all, but it is a degree of freedom the search samples the same way -- one
    /// gene, mutated and crossed over like any other -- so it rides in the same list. Turning a
    /// single ring bond would prise the ring open; a ring changes shape by moving its atoms out of
    /// their mean plane, which is what [`super::rings`] does.
    RingPucker,
}

#[derive(Debug, Clone)]
pub struct PreparedTorsion {
    pub atoms: [usize; 4],
    pub central_bond: [usize; 2],
    pub moving_atoms: Vec<usize>,
    pub kind: TorsionKind,
    /// Set only for [`TorsionKind::RingPucker`]: the ring this gene reshapes.
    pub ring: Option<PreparedRing>,
    pub user_selected: bool,
}

/// A ring the search can reshape, with everything needed to do it.
#[derive(Debug, Clone)]
pub struct PreparedRing {
    /// Ring atoms in connected order.
    pub atoms: Vec<usize>,
    /// For each ring atom, the substituent atoms that move with it.
    pub branches: Vec<Vec<usize>>,
    /// The canonical conformers this ring can take, which the gene indexes.
    pub conformers: Vec<Pucker>,
}

fn adjacency(natoms: usize, bonds: &[[usize; 2]]) -> Vec<Vec<usize>> {
    let mut graph = vec![Vec::new(); natoms];
    for &[a, b] in bonds {
        graph[a].push(b);
        graph[b].push(a);
    }
    for neighbours in &mut graph {
        neighbours.sort_unstable();
    }
    graph
}

fn has_edge(graph: &[Vec<usize>], a: usize, b: usize) -> bool {
    graph[a].binary_search(&b).is_ok()
}

fn component_without_edge(graph: &[Vec<usize>], start: usize, blocked: [usize; 2]) -> Vec<usize> {
    let mut seen = vec![false; graph.len()];
    let mut queue = VecDeque::from([start]);
    seen[start] = true;
    while let Some(atom) = queue.pop_front() {
        for &next in &graph[atom] {
            if (atom == blocked[0] && next == blocked[1])
                || (atom == blocked[1] && next == blocked[0])
            {
                continue;
            }
            if !seen[next] {
                seen[next] = true;
                queue.push_back(next);
            }
        }
    }
    seen.into_iter()
        .enumerate()
        .filter_map(|(atom, present)| present.then_some(atom))
        .collect()
}

fn edge_is_in_ring(graph: &[Vec<usize>], a: usize, b: usize) -> bool {
    component_without_edge(graph, a, [a, b]).contains(&b)
}

fn methyl_endpoint(
    elements: &[String],
    graph: &[Vec<usize>],
    endpoint: usize,
    other: usize,
) -> bool {
    elements[endpoint].eq_ignore_ascii_case("C")
        && graph[endpoint].len() == 4
        && graph[endpoint]
            .iter()
            .filter(|&&atom| atom != other)
            .all(|&atom| elements[atom].eq_ignore_ascii_case("H"))
}

fn representative_neighbour(
    graph: &[Vec<usize>],
    elements: &[String],
    center: usize,
    excluded: usize,
) -> Option<usize> {
    graph[center]
        .iter()
        .copied()
        .filter(|&atom| atom != excluded)
        .min_by_key(|&atom| (elements[atom].eq_ignore_ascii_case("H"), atom))
}

fn prepare_explicit(
    atoms: [usize; 4],
    kind: TorsionKind,
    natoms: usize,
    graph: &[Vec<usize>],
) -> Result<PreparedTorsion> {
    ensure!(
        atoms.iter().all(|&atom| atom < natoms),
        "conformer-search torsion {:?} contains an atom index outside 0..{}",
        atoms,
        natoms.saturating_sub(1)
    );
    ensure!(
        atoms.iter().copied().collect::<HashSet<_>>().len() == 4,
        "conformer-search torsion {:?} must contain four distinct atoms",
        atoms
    );
    for pair in atoms.windows(2) {
        ensure!(
            has_edge(graph, pair[0], pair[1]),
            "conformer-search torsion {:?} is not a bonded atom sequence",
            atoms
        );
    }
    let central = [atoms[1], atoms[2]];
    ensure!(
        !edge_is_in_ring(graph, central[0], central[1]),
        "conformer-search torsion {:?} rotates a ring bond; ring-bond torsions are not supported",
        atoms
    );
    let moving_atoms = component_without_edge(graph, central[1], central);
    ensure!(
        !moving_atoms.contains(&central[0]),
        "conformer-search torsion {:?} does not divide the molecular graph",
        atoms
    );
    Ok(PreparedTorsion {
        atoms,
        central_bond: central,
        moving_atoms,
        ring: None,
        kind,
        user_selected: true,
    })
}

/// Discover acyclic, non-terminal torsions or validate an explicit user list.
///
/// XYZ input does not carry bond orders. Consequently automatic discovery is
/// deliberately geometric and should be overridden for amides, conjugated
/// bonds, or any chemically curated search space.
pub fn discover_torsions(
    coordinates_bohr: &[f64],
    elements: &[String],
    connectivity: &ConnectivityModel,
    explicit_rotatable: &[[usize; 4]],
    explicit_cis_trans: &[[usize; 4]],
    exclude_methyl_rotors: bool,
) -> Result<(Vec<[usize; 2]>, Vec<PreparedTorsion>)> {
    let natoms = elements.len();
    ensure!(
        coordinates_bohr.len() == 3 * natoms,
        "conformer-search topology expected {} Cartesian coordinates, got {}",
        3 * natoms,
        coordinates_bohr.len()
    );
    ensure!(
        connectivity.rcov_disp.len() == natoms,
        "conformer-search connectivity has {} radii for {natoms} atoms",
        connectivity.rcov_disp.len()
    );
    if natoms < 4 {
        bail!("conformer search requires at least four atoms and one torsional degree of freedom");
    }

    let internals = analyze_structure(coordinates_bohr, connectivity)?;
    let graph = adjacency(natoms, &internals.bonds);
    let mut torsions = Vec::new();
    let mut central_bonds = HashSet::new();

    for &atoms in explicit_cis_trans {
        let prepared = prepare_explicit(atoms, TorsionKind::CisTrans, natoms, &graph)?;
        let key = sorted_bond(prepared.central_bond);
        ensure!(
            central_bonds.insert(key),
            "central bond {}-{} is listed more than once in conformer-search torsions",
            key[0] + 1,
            key[1] + 1
        );
        torsions.push(prepared);
    }

    if !explicit_rotatable.is_empty() {
        for &atoms in explicit_rotatable {
            let prepared = prepare_explicit(atoms, TorsionKind::Rotatable, natoms, &graph)?;
            let key = sorted_bond(prepared.central_bond);
            ensure!(
                central_bonds.insert(key),
                "central bond {}-{} is listed more than once in conformer-search torsions",
                key[0] + 1,
                key[1] + 1
            );
            torsions.push(prepared);
        }
    } else {
        for &[left, right] in &internals.bonds {
            let key = sorted_bond([left, right]);
            if central_bonds.contains(&key)
                || graph[left].len() < 2
                || graph[right].len() < 2
                || edge_is_in_ring(&graph, left, right)
                || (exclude_methyl_rotors
                    && (methyl_endpoint(elements, &graph, left, right)
                        || methyl_endpoint(elements, &graph, right, left)))
            {
                continue;
            }
            let Some(a) = representative_neighbour(&graph, elements, left, right) else {
                continue;
            };
            let Some(d) = representative_neighbour(&graph, elements, right, left) else {
                continue;
            };
            torsions.push(PreparedTorsion {
                atoms: [a, left, right, d],
                central_bond: [left, right],
                moving_atoms: component_without_edge(&graph, right, [left, right]),
                kind: TorsionKind::Rotatable,
                ring: None,
                user_selected: false,
            });
            central_bonds.insert(key);
        }
    }

    // Rings are degrees of freedom too, and ones no rotatable bond can reach. Appended after the
    // torsions so a genome reads torsions-then-rings.
    for ring in rings::find_rings(natoms, &internals.bonds) {
        let conformers = rings::canonical_conformers(ring.size());
        if conformers.is_empty() {
            continue;
        }
        let branches = rings::substituent_branches(natoms, &internals.bonds, &ring.atoms);
        torsions.push(PreparedTorsion {
            // A ring has no single dihedral or central bond; these fields go unused for this kind
            // and are filled from the ring so they stay in range.
            atoms: [ring.atoms[0], ring.atoms[1], ring.atoms[2], ring.atoms[3]],
            central_bond: [ring.atoms[0], ring.atoms[1]],
            moving_atoms: Vec::new(),
            kind: TorsionKind::RingPucker,
            ring: Some(PreparedRing {
                atoms: ring.atoms,
                branches,
                conformers,
            }),
            user_selected: false,
        });
    }

    if torsions.is_empty() {
        bail!("no rotatable torsion was found; specify `torsions` explicitly if the automatic geometric selection omitted the intended bond");
    }
    Ok((internals.bonds, torsions))
}

fn sorted_bond(mut bond: [usize; 2]) -> [usize; 2] {
    if bond[1] < bond[0] {
        bond.swap(0, 1);
    }
    bond
}
