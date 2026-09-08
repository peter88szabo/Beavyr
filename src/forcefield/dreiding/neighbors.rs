//! Verlet neighbour list for the nonbonded terms, built from a cell list.
//!
//! The bonded terms are a fixed, short list, but the nonbonded ones are O(N²) if computed
//! naively. For a one-off cleanup that hardly matters; for molecular dynamics, where the same
//! pairs are needed every step for millions of steps, it is the whole cost. So the list is built
//! once and reused until an atom has moved far enough that a pair could have entered the cutoff
//! unnoticed.
//!
//! The list lives in the per-worker [`super::Workspace`], not in the shared topology, which is
//! what lets a conformational search minimise many structures in parallel against one topology.

/// Nonbonded cutoff (Å). DREIDING does not prescribe one; this is the usual range for a 12-6
/// potential, by which point the well is ~1e-4 of its depth.
pub const CUTOFF: f64 = 11.0;

/// Extra margin (Å) included when building, so the list survives some motion before a rebuild.
pub const SKIN: f64 = 2.0;

const CUTOFF_WITH_SKIN: f64 = CUTOFF + SKIN;

/// Pairs within the cutoff plus skin, excluding bonded exclusions.
#[derive(Debug, Default, Clone)]
pub struct NeighborList {
    pairs: Vec<(u32, u32)>,
    /// Positions the list was built at, to decide when it is stale.
    reference: Vec<f64>,
}

impl NeighborList {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn pairs(&self) -> &[(u32, u32)] {
        &self.pairs
    }

    /// Rebuilds if any atom has moved more than half the skin since the last build.
    ///
    /// Half, not the whole skin, because *two* atoms can move toward each other -- each by up to
    /// half -- and close the margin between them.
    pub fn refresh<F>(&mut self, pos: &[f64], is_excluded: F)
    where
        F: Fn(usize, usize) -> bool,
    {
        if self.is_stale(pos) {
            self.rebuild(pos, is_excluded);
        }
    }

    fn is_stale(&self, pos: &[f64]) -> bool {
        if self.reference.len() != pos.len() {
            return true;
        }
        let limit = (SKIN * 0.5) * (SKIN * 0.5);
        pos.chunks_exact(3)
            .zip(self.reference.chunks_exact(3))
            .any(|(now, then)| {
                let (dx, dy, dz) = (now[0] - then[0], now[1] - then[1], now[2] - then[2]);
                dx * dx + dy * dy + dz * dz > limit
            })
    }

    fn rebuild<F>(&mut self, pos: &[f64], is_excluded: F)
    where
        F: Fn(usize, usize) -> bool,
    {
        let natoms = pos.len() / 3;
        self.pairs.clear();
        self.reference.clear();
        self.reference.extend_from_slice(pos);

        // Below this count the cell-list bookkeeping costs more than the pairs it saves.
        const DIRECT_LOOP_LIMIT: usize = 64;
        let cut2 = CUTOFF_WITH_SKIN * CUTOFF_WITH_SKIN;

        if natoms <= DIRECT_LOOP_LIMIT {
            for i in 0..natoms {
                for j in (i + 1)..natoms {
                    if !is_excluded(i, j) && distance_squared(pos, i, j) <= cut2 {
                        self.pairs.push((i as u32, j as u32));
                    }
                }
            }
            return;
        }

        let cells = CellList::build(pos, CUTOFF_WITH_SKIN);
        let mut candidates = Vec::new();
        for i in 0..natoms {
            candidates.clear();
            cells.neighbors_of(pos, i, &mut candidates);
            for &j in &candidates {
                // Each pair is visited from both atoms; keep it once.
                if j > i && !is_excluded(i, j) && distance_squared(pos, i, j) <= cut2 {
                    self.pairs.push((i as u32, j as u32));
                }
            }
        }
    }
}

#[inline]
fn distance_squared(pos: &[f64], i: usize, j: usize) -> f64 {
    let dx = pos[3 * i] - pos[3 * j];
    let dy = pos[3 * i + 1] - pos[3 * j + 1];
    let dz = pos[3 * i + 2] - pos[3 * j + 2];
    dx * dx + dy * dy + dz * dz
}

/// A uniform spatial grid, so each atom only tests atoms in the 27 cells around it.
///
/// Beavyr targets molecules rather than crystals, so there is no periodic wrapping here. A
/// minimum-image convention would go in [`Self::neighbors_of`] if that ever changed.
struct CellList {
    /// Atom indices, grouped by cell.
    contents: Vec<Vec<usize>>,
    dims: [usize; 3],
    origin: [f64; 3],
    cell_size: f64,
}

impl CellList {
    fn build(pos: &[f64], cell_size: f64) -> Self {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for atom in pos.chunks_exact(3) {
            for axis in 0..3 {
                lo[axis] = lo[axis].min(atom[axis]);
                hi[axis] = hi[axis].max(atom[axis]);
            }
        }
        let dims = std::array::from_fn(|axis| {
            (((hi[axis] - lo[axis]) / cell_size).floor() as usize + 1).max(1)
        });
        let mut list = Self {
            contents: vec![Vec::new(); dims[0] * dims[1] * dims[2]],
            dims,
            origin: lo,
            cell_size,
        };
        for i in 0..pos.len() / 3 {
            let idx = list.flat_index(list.cell_of(pos, i));
            list.contents[idx].push(i);
        }
        list
    }

    fn cell_of(&self, pos: &[f64], i: usize) -> [usize; 3] {
        std::array::from_fn(|axis| {
            let raw = ((pos[3 * i + axis] - self.origin[axis]) / self.cell_size).floor();
            (raw.max(0.0) as usize).min(self.dims[axis] - 1)
        })
    }

    fn flat_index(&self, c: [usize; 3]) -> usize {
        (c[2] * self.dims[1] + c[1]) * self.dims[0] + c[0]
    }

    fn neighbors_of(&self, pos: &[f64], i: usize, out: &mut Vec<usize>) {
        let c = self.cell_of(pos, i);
        for dz in -1i64..=1 {
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let n: Option<[usize; 3]> = (|| {
                        let mut n = [0usize; 3];
                        for (axis, d) in [dx, dy, dz].into_iter().enumerate() {
                            let v = c[axis] as i64 + d;
                            if v < 0 || v >= self.dims[axis] as i64 {
                                return None;
                            }
                            n[axis] = v as usize;
                        }
                        Some(n)
                    })();
                    if let Some(n) = n {
                        out.extend_from_slice(&self.contents[self.flat_index(n)]);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two atoms in a line, spaced so only some pairs fall inside the cutoff.
    fn chain(natoms: usize, spacing: f64) -> Vec<f64> {
        (0..natoms)
            .flat_map(|i| [i as f64 * spacing, 0.0, 0.0])
            .collect()
    }

    #[test]
    fn finds_every_pair_within_the_cutoff() {
        let pos = chain(10, 2.0);
        let mut list = NeighborList::new();
        list.refresh(&pos, |_, _| false);
        // Spacing 2 Å over 10 atoms spans 18 Å, so distant pairs must be dropped.
        let cut2 = (CUTOFF + SKIN) * (CUTOFF + SKIN);
        let expected = (0..10)
            .flat_map(|i| ((i + 1)..10).map(move |j| (i, j)))
            .filter(|&(i, j)| distance_squared(&pos, i, j) <= cut2)
            .count();
        assert_eq!(list.pairs().len(), expected);
        assert!(list.pairs().len() < 45, "nothing was actually excluded");
    }

    /// The cell-list path and the direct path must agree, or behaviour would change silently at
    /// whatever size the switch happens.
    #[test]
    fn cell_list_and_direct_loop_agree() {
        // Two blobs, one either side of the size threshold, with the same geometry.
        let natoms = 200;
        let pos: Vec<f64> = (0..natoms)
            .flat_map(|i| {
                let f = i as f64;
                [f % 7.0 * 1.7, (f / 7.0) % 7.0 * 1.9, (f / 49.0) * 2.1]
            })
            .collect();

        let mut via_cells = NeighborList::new();
        via_cells.refresh(&pos, |_, _| false);

        let cut2 = (CUTOFF + SKIN) * (CUTOFF + SKIN);
        let mut direct: Vec<(u32, u32)> = (0..natoms)
            .flat_map(|i| ((i + 1)..natoms).map(move |j| (i, j)))
            .filter(|&(i, j)| distance_squared(&pos, i, j) <= cut2)
            .map(|(i, j)| (i as u32, j as u32))
            .collect();

        let mut from_cells = via_cells.pairs().to_vec();
        from_cells.sort_unstable();
        direct.sort_unstable();
        assert_eq!(from_cells, direct);
    }

    #[test]
    fn honours_exclusions() {
        let pos = chain(4, 1.5);
        let mut list = NeighborList::new();
        // Exclude everything adjacent, as 1-2 exclusion would.
        list.refresh(&pos, |i, j| j == i + 1);
        assert!(list.pairs().iter().all(|&(i, j)| j != i + 1));
        assert_eq!(list.pairs().len(), 3); // (0,2) (0,3) (1,3)
    }

    #[test]
    fn does_not_rebuild_for_small_motion_but_does_for_large() {
        let pos = chain(8, 2.0);
        let mut list = NeighborList::new();
        list.refresh(&pos, |_, _| false);

        // A nudge well under half the skin leaves the reference untouched.
        let mut nudged = pos.clone();
        nudged[0] += SKIN * 0.1;
        list.refresh(&nudged, |_, _| false);
        assert_eq!(list.reference[0], pos[0], "rebuilt when it did not need to");

        // A move past half the skin must rebuild, and the rebuilt list must be right for the new
        // geometry -- shoving an end atom 2 Å along this chain really does pull a distant pair
        // inside the cutoff, so the count is expected to change.
        let mut shoved = pos.clone();
        shoved[0] += SKIN;
        list.refresh(&shoved, |_, _| false);
        assert_eq!(list.reference[0], shoved[0], "failed to rebuild");

        let cut2 = (CUTOFF + SKIN) * (CUTOFF + SKIN);
        let natoms = shoved.len() / 3;
        let expected = (0..natoms)
            .flat_map(|i| ((i + 1)..natoms).map(move |j| (i, j)))
            .filter(|&(i, j)| distance_squared(&shoved, i, j) <= cut2)
            .count();
        assert_eq!(list.pairs().len(), expected);
    }
}
