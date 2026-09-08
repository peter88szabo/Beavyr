//! The DREIDING generic force field.
//!
//! S. L. Mayo, B. D. Olafson and W. A. Goddard III, "DREIDING: A Generic Force Field for
//! Molecular Simulations", *J. Phys. Chem.* **1990**, *94*, 8897-8909. Equation and table
//! numbers in this module refer to that paper.
//!
//! DREIDING derives its parameters from element and hybridization rather than tabulating them per
//! functional group, so a small number of rules covers the whole main group and several metals.
//! That generality is what makes it the right fit for cleaning up structures a user has just
//! sketched, where the chemistry is not known in advance.
//!
//! # Layout
//!
//! * [`typer`] -- atom typing and topology perception (vendored; see its `README.md`)
//! * [`params`] -- Tables I-V, transcribed
//! * [`terms`] -- the six energy terms and their analytic gradients
//! * [`neighbors`] -- Verlet neighbour list for the nonbonded terms
//!
//! # Design
//!
//! See `docs/superpowers/specs/2026-09-09-dreiding-force-field-design.md`. In outline: an
//! immutable `DreidingTopology`, built once and shareable across threads, plus a per-worker
//! mutable `Workspace` holding the neighbour list and minimiser history. Geometry cleanup,
//! molecular dynamics and conformational search are then all drivers over one
//! `energy_and_forces`, rather than three engines.

pub mod typer;

/// DREIDING parameters, transcribed from the paper.
pub mod params;

/// Energy terms and their analytic gradients.
pub mod terms;

/// Verlet neighbour list for the nonbonded terms.
pub mod neighbors;

use std::collections::HashSet;
use std::fmt;

use crate::bond_order::estimate_bond_order;
use crate::molecule::{atomic_mass_amu, Molecule};

use self::typer::{
    assign_topology, Element, GraphBondOrder, MolecularGraph, MolecularTopology,
    TopologyBondOrder,
};

/// Harmonic bond stretch, Eq 4a.
#[derive(Debug, Clone, Copy)]
struct BondTerm {
    i: usize,
    j: usize,
    /// Force constant, kcal/mol/Å².
    k: f64,
    /// Equilibrium length, Å.
    r_e: f64,
}

/// Harmonic-cosine angle bend, Eq 10a, or Eq 10' when `linear`.
#[derive(Debug, Clone, Copy)]
struct AngleTerm {
    i: usize,
    j: usize,
    k: usize,
    /// `K/sin²θ⁰` for Eq 10a, or plain `K` for Eq 10'.
    c: f64,
    cos_theta0: f64,
    linear: bool,
}

/// Proper torsion, Eq 13. `v` already carries the division by the number of (I,L) pairs.
#[derive(Debug, Clone, Copy)]
struct TorsionTerm {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    v: f64,
    n: f64,
    phi0: f64,
}

/// Inversion, Eq 28. One per permutation, so each carries its 1/3 weight in `c`.
#[derive(Debug, Clone, Copy)]
struct InversionTerm {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    c: f64,
    cos_psi0: f64,
    planar: bool,
}

/// Per-atom Lennard-Jones parameters, Table II.
#[derive(Debug, Clone, Copy)]
struct VdwParam {
    r0: f64,
    d0: f64,
}

/// Explicit hydrogen bond, Table V.
#[derive(Debug, Clone, Copy)]
struct HBondTerm {
    donor: usize,
    hydrogen: usize,
    acceptor: usize,
}

/// A rotatable single bond and the atoms that move when it turns.
///
/// Not used by the minimiser. It is here because conformational search needs it and it falls out
/// of perception the topology already does; deriving it later would mean re-perceiving rings.
#[derive(Debug, Clone)]
pub struct RotatableBond {
    pub j: usize,
    pub k: usize,
    /// Atoms on `k`'s side of the bond, which a rotation about `j-k` moves.
    pub moving: Vec<usize>,
}

/// 1-2 and 1-3 pairs, which DREIDING excludes from the nonbonded terms.
///
/// 1-4 pairs are deliberately *not* here: DREIDING specifies no 1-4 scaling, so those interact at
/// full strength.
#[derive(Debug, Default, Clone)]
struct Exclusions {
    pairs: HashSet<(u32, u32)>,
}

impl Exclusions {
    fn insert(&mut self, i: usize, j: usize) {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.pairs.insert((a as u32, b as u32));
    }

    fn contains(&self, i: usize, j: usize) -> bool {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.pairs.contains(&(a as u32, b as u32))
    }
}

/// Why a molecule could not be given a DREIDING energy.
#[derive(Debug)]
pub enum BuildError {
    /// An element symbol Beavyr holds as text is not one the typer recognises.
    UnknownElement { index: usize, symbol: String },
    /// Chemical perception or rule-based typing failed.
    Typing(String),
    /// The typer assigned a type the paper does not parameterise.
    ///
    /// Reported rather than approximated: an invented bond radius yields a confident, wrong
    /// geometry, which is worse than no geometry. See [`params::is_parameterized`].
    Unparameterized {
        index: usize,
        symbol: String,
        atom_type: String,
    },
    /// The structure has no atoms, or no bonds to work with.
    Empty,
    /// The bond list is not chemically possible, so typing it would be meaningless.
    ///
    /// This is nearly always the display bond threshold being too generous rather than a bad
    /// structure: at a wide enough threshold the two hydrogens of water are "bonded" to each
    /// other, which makes a three-membered ring that perception quite reasonably reads as
    /// aromatic, and every atom is then mistyped. Caught here because the alternative is a
    /// confident, wrong geometry with nothing to indicate why.
    ImplausibleConnectivity { index: usize, symbol: String, degree: usize },
}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownElement { index, symbol } => write!(
                f,
                "atom {} is '{symbol}', which is not a recognised element",
                index + 1
            ),
            Self::Typing(message) => write!(f, "could not perceive the chemistry: {message}"),
            Self::Unparameterized {
                index,
                symbol,
                atom_type,
            } => write!(
                f,
                "atom {} ({symbol}) types as {atom_type}, which DREIDING does not parameterise. \
                 The published parameters cover 34 of the 56 types the typer can assign; the gaps \
                 are the metals outside Na, Ca, Fe and Zn, plus S_2 and S_R (so thiophene, \
                 thiones and sulfoxides are not covered).",
                index + 1
            ),
            Self::Empty => write!(f, "there is no structure to work on"),
            Self::ImplausibleConnectivity {
                index,
                symbol,
                degree,
            } => write!(
                f,
                "atom {} ({symbol}) is shown with {degree} bonds, which is not chemically \
                 possible. The bond threshold in Appearance is probably too high -- lower it \
                 until the drawn bonds are right, then try again.",
                index + 1
            ),
        }
    }
}

impl std::error::Error for BuildError {}

/// Energy split by term, in kcal/mol.
///
/// Reported after a cleanup, and the main handle for diagnosing a surprising result: a term that
/// has gone wrong shows up here as an absurd component long before the geometry explains itself.
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct EnergyBreakdown {
    pub bond: f64,
    pub angle: f64,
    pub torsion: f64,
    pub inversion: f64,
    pub vdw: f64,
    pub hbond: f64,
}

impl EnergyBreakdown {
    pub fn total(&self) -> f64 {
        self.bond + self.angle + self.torsion + self.inversion + self.vdw + self.hbond
    }
}

/// Per-worker mutable scratch: the neighbour list and the positions it was built at.
///
/// Separate from [`DreidingTopology`] so that many workers can share one topology and minimise
/// different structures in parallel -- which is what conformational search needs, and what an
/// engine holding its neighbour list behind `&mut self` would prevent.
#[derive(Debug, Default, Clone)]
pub struct Workspace {
    neighbors: neighbors::NeighborList,
}

impl Workspace {
    pub fn new() -> Self {
        Self::default()
    }
}

/// An immutable DREIDING topology: atom types and every term with its parameters resolved.
///
/// Built once by [`Self::build`]; positions then vary freely. Nothing in the hot path allocates or
/// re-perceives chemistry, so a million dynamics steps re-run only [`Self::energy_and_forces`].
#[derive(Debug, Clone)]
pub struct DreidingTopology {
    natoms: usize,
    masses_amu: Vec<f64>,
    atom_types: Vec<String>,
    bonds: Vec<BondTerm>,
    angles: Vec<AngleTerm>,
    torsions: Vec<TorsionTerm>,
    inversions: Vec<InversionTerm>,
    vdw: Vec<VdwParam>,
    hbonds: Vec<HBondTerm>,
    exclusions: Exclusions,
    rotatable: Vec<RotatableBond>,
}

impl DreidingTopology {
    pub fn natoms(&self) -> usize {
        self.natoms
    }

    pub fn masses_amu(&self) -> &[f64] {
        &self.masses_amu
    }

    /// The perceived DREIDING type of each atom.
    ///
    /// Worth showing the user alongside an energy: bond orders are perceived from *geometry*, so a
    /// badly sketched structure can be mistyped -- a stretched double bond read as single -- and
    /// seeing the types makes a surprising result diagnosable instead of mysterious.
    pub fn atom_types(&self) -> &[String] {
        &self.atom_types
    }

    /// Rotatable bonds and the atoms each one moves. For conformational search.
    pub fn rotatable_bonds(&self) -> &[RotatableBond] {
        &self.rotatable
    }

    /// Writes `-dE/dx` (kcal/mol/Å) into `forces` and returns the energy (kcal/mol).
    ///
    /// Takes `&self`, so N workers can share one topology.
    pub fn energy_and_forces(
        &self,
        pos: &[f64],
        work: &mut Workspace,
        forces: &mut [f64],
    ) -> f64 {
        let breakdown = self.accumulate(pos, work, forces);
        for f in forces.iter_mut() {
            *f = -*f;
        }
        breakdown.total()
    }

    /// Energy split by term, without keeping the gradient.
    pub fn energy_breakdown(&self, pos: &[f64], work: &mut Workspace) -> EnergyBreakdown {
        let mut scratch = vec![0.0; self.natoms * 3];
        self.accumulate(pos, work, &mut scratch)
    }

    /// Sums every term, leaving `dE/dx` (not forces) in `grad`.
    fn accumulate(
        &self,
        pos: &[f64],
        work: &mut Workspace,
        grad: &mut [f64],
    ) -> EnergyBreakdown {
        grad.iter_mut().for_each(|g| *g = 0.0);
        let mut e = EnergyBreakdown::default();

        for b in &self.bonds {
            e.bond += terms::bond(pos, grad, b.i, b.j, b.k, b.r_e);
        }
        for a in &self.angles {
            e.angle += terms::angle(pos, grad, a.i, a.j, a.k, a.c, a.cos_theta0, a.linear);
        }
        for t in &self.torsions {
            e.torsion += terms::torsion(pos, grad, t.i, t.j, t.k, t.l, t.v, t.n, t.phi0);
        }
        for v in &self.inversions {
            e.inversion +=
                terms::inversion(pos, grad, v.i, v.j, v.k, v.l, v.c, v.cos_psi0, v.planar);
        }

        work.neighbors
            .refresh(pos, |i, j| self.exclusions.contains(i, j));
        for &(i, j) in work.neighbors.pairs() {
            let (i, j) = (i as usize, j as usize);
            // Geometric-mean combination, per the paper.
            let r0 = (self.vdw[i].r0 * self.vdw[j].r0).sqrt();
            let d0 = (self.vdw[i].d0 * self.vdw[j].d0).sqrt();
            e.vdw += terms::vdw_pair(pos, grad, i, j, r0, d0);
        }

        for h in &self.hbonds {
            e.hbond += terms::hbond(
                pos,
                grad,
                h.donor,
                h.hydrogen,
                h.acceptor,
                params::D_HB,
                params::R_HB,
            );
        }
        e
    }
}

/// Where an atom type sits for the purposes of choosing a torsion barrier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TorsionClass {
    Sp3,
    Sp2,
    Resonant,
    Sp1,
    /// Monovalent atoms and metals. The paper gives these no torsional barrier at all (Eq 20).
    Terminal,
}

/// Reads the hybridisation off the type-name suffix, which is how DREIDING encodes it.
///
/// The main-group types outside the first row use a bare digit (`Si3`, `Te3`) where the first row
/// uses an underscore (`C_3`, `O_3`), so the sp3 test is on the trailing `3` either way.
fn torsion_class(atom_type: &str) -> TorsionClass {
    if atom_type.ends_with("_1") {
        TorsionClass::Sp1
    } else if atom_type.ends_with("_2") {
        TorsionClass::Sp2
    } else if atom_type.ends_with("_R") {
        TorsionClass::Resonant
    } else if atom_type.ends_with('3') {
        TorsionClass::Sp3
    } else {
        TorsionClass::Terminal
    }
}

/// Elements of the oxygen column, which the paper treats separately for torsions.
///
/// Its reasoning: these are written `X_3` when they carry two single bonds, but their barriers are
/// better understood from the s²p⁴ configuration, with a lone pair perpendicular to the two bonds.
/// That is why HSSH prefers a ~90° torsion, and crystalline HOOH likewise.
fn is_column_16(element: Element) -> bool {
    matches!(element, Element::O | Element::S | Element::Se | Element::Te)
}

/// The DREIDING torsion barrier for one dihedral, before division by the (I,L) count.
///
/// Returns `None` where the paper prescribes no barrier. Implements Eqs 14-23, which are keyed on
/// the two central atoms and -- for case (j) -- on the terminal atom too.
fn torsion_parameters(
    class_i: TorsionClass,
    class_j: TorsionClass,
    class_k: TorsionClass,
    col16_j: bool,
    col16_k: bool,
    order: TopologyBondOrder,
) -> Option<(f64, f64, f64)> {
    use TorsionClass::{Resonant, Sp1, Sp2, Sp3, Terminal};

    // (g), Eq 20: sp1 centres, monovalent atoms and metals get no barrier.
    if matches!(class_j, Sp1 | Terminal) || matches!(class_k, Sp1 | Terminal) {
        return None;
    }

    let sp2ish = |c: TorsionClass| matches!(c, Sp2 | Resonant);

    match (class_j, class_k) {
        // Both sp3.
        (Sp3, Sp3) => {
            if col16_j && col16_k {
                // (h), Eq 21: two oxygen-column sp3 atoms, e.g. HOOH and HSSH.
                Some((2.0, 2.0, 90.0))
            } else {
                // (a), Eq 14. The paper notes explicitly that an oxygen-column sp3 atom bonded to
                // an sp3 atom of another column falls back to this.
                Some((2.0, 3.0, 180.0))
            }
        }

        // One sp3, one sp2 or resonant.
        (Sp3, k) if sp2ish(k) => single_bond_sp3_sp2(class_i, col16_j),
        (j, Sp3) if sp2ish(j) => single_bond_sp3_sp2(class_i, col16_k),

        // Both sp2 or resonant: the bond order decides.
        (j, k) if sp2ish(j) && sp2ish(k) => match order {
            // (c), Eq 16.
            TopologyBondOrder::Double => Some((45.0, 2.0, 180.0)),
            // (d), Eq 17: a resonance bond of order 1.5.
            TopologyBondOrder::Resonant => Some((25.0, 2.0, 180.0)),
            // (e), Eq 18, e.g. the middle bond of butadiene. Case (f) -- the exocyclic
            // aromatic-aromatic exception of Eq 19 -- is applied by the caller, which is where
            // ring membership is known.
            TopologyBondOrder::Single => Some((5.0, 2.0, 180.0)),
            TopologyBondOrder::Triple => None,
        },

        _ => None,
    }
}

/// Cases (b)/(i)/(j) -- a single bond joining an sp3 centre to an sp2 or resonant one.
///
/// `col16_sp3` says whether the sp3 end is an oxygen-column atom.
fn single_bond_sp3_sp2(class_i: TorsionClass, col16_sp3: bool) -> Option<(f64, f64, f64)> {
    use TorsionClass::{Resonant, Sp2};
    if col16_sp3 {
        // (i), Eq 22: the oxygen-like lone pair prefers to overlap the sp2 orbitals, giving a
        // planar preference rather than the 6-fold barrier of Eq 15.
        return Some((2.0, 2.0, 180.0));
    }
    if matches!(class_i, Sp2 | Resonant) {
        // (b), Eq 15: 6-fold, as for the C-C bond of acetate.
        Some((1.0, 6.0, 0.0))
    } else {
        // (j), Eq 23: propene's 3-fold barrier, with the sp3 centre eclipsing the double bond.
        // The distinction from (b) is the terminal atom I, not the central pair.
        Some((2.0, 3.0, 180.0))
    }
}

impl DreidingTopology {
    /// Types a molecule and resolves every DREIDING term for it.
    ///
    /// # Errors
    ///
    /// [`BuildError::Unparameterized`] if any atom types to something the paper does not cover.
    /// The alternative -- guessing a bond radius -- would return a confident wrong geometry, so it
    /// is reported instead and the caller refuses.
    pub fn build(mol: &Molecule) -> Result<Self, BuildError> {
        let natoms = mol.atoms.len();
        if natoms == 0 {
            return Err(BuildError::Empty);
        }

        // --- Beavyr's molecule to the typer's graph -------------------------------------------
        //
        // Beavyr stores bonds as index pairs plus a distance and perceives order from geometry, so
        // the discrete order the typer wants comes from `bond_order::estimate_bond_order`.
        let mut graph = MolecularGraph::new();
        for (index, symbol) in mol.atoms.iter().enumerate() {
            let element = symbol.parse::<Element>().map_err(|_| {
                BuildError::UnknownElement {
                    index,
                    symbol: symbol.clone(),
                }
            })?;
            graph.add_atom(element);
        }
        for &(i, j, distance) in &mol.bonds {
            let order = estimate_bond_order(&mol.atoms[i], &mol.atoms[j], distance);
            graph
                .add_bond(i, j, discrete_order(order))
                .map_err(|e| BuildError::Typing(e.to_string()))?;
        }

        check_connectivity_is_plausible(mol)?;

        let topology =
            assign_topology(&graph).map_err(|e| BuildError::Typing(describe_typer_error(&e)))?;

        // --- Refuse before doing any work if coverage is missing ------------------------------
        for atom in &topology.atoms {
            if !params::is_parameterized(&atom.atom_type) {
                return Err(BuildError::Unparameterized {
                    index: atom.id,
                    symbol: mol.atoms[atom.id].clone(),
                    atom_type: atom.atom_type.clone(),
                });
            }
        }

        let atom_types: Vec<String> = topology
            .atoms
            .iter()
            .map(|a| a.atom_type.clone())
            .collect();
        // An unknown mass would only arise for an element the typer already rejected; 12.0 keeps
        // dynamics finite rather than propagating a NaN if that ever changes.
        let masses_amu = mol
            .atoms
            .iter()
            .map(|s| atomic_mass_amu(s).unwrap_or(12.0))
            .collect();

        let adjacency = build_adjacency(natoms, &topology);
        let mut me = Self {
            natoms,
            masses_amu,
            atom_types,
            bonds: Vec::new(),
            angles: Vec::new(),
            torsions: Vec::new(),
            inversions: Vec::new(),
            vdw: Vec::new(),
            hbonds: Vec::new(),
            exclusions: Exclusions::default(),
            rotatable: Vec::new(),
        };
        me.resolve_bonds(&topology);
        me.resolve_angles(&topology);
        me.resolve_torsions(&topology, &adjacency);
        me.resolve_inversions(&topology);
        me.resolve_vdw(&topology);
        me.resolve_hbonds(&topology, &adjacency);
        me.rotatable = rotatable_bonds(&topology, &adjacency);
        Ok(me)
    }

    /// Eqs 6-9: `Rₑ = R⁰_I + R⁰_J − 0.01 Å` and `K(n) = n × 700`.
    fn resolve_bonds(&mut self, topology: &MolecularTopology) {
        for bond in &topology.bonds {
            let (i, j) = bond.atom_ids;
            let (r_i, _) = params::geometry(&self.atom_types[i]).expect("coverage checked");
            let (r_j, _) = params::geometry(&self.atom_types[j]).expect("coverage checked");
            let n = order_multiplicity(bond.order);
            self.bonds.push(BondTerm {
                i,
                j,
                k: params::BOND_K_PER_ORDER * n,
                r_e: r_i + r_j - params::BOND_LENGTH_CORRECTION,
            });
            self.exclusions.insert(i, j);
        }
    }

    /// Eq 10a, or Eq 10' at a linear centre. `θ⁰` comes from the *central* atom's type.
    fn resolve_angles(&mut self, topology: &MolecularTopology) {
        for angle in &topology.angles {
            let (i, j, k) = angle.atom_ids;
            let (_, theta0_deg) = params::geometry(&self.atom_types[j]).expect("coverage checked");
            let theta0 = theta0_deg.to_radians();
            // A linear centre makes K/sin²θ⁰ diverge, so Eq 10' is used there instead.
            let linear = theta0_deg >= 179.999;
            let c = if linear {
                params::ANGLE_K
            } else {
                params::ANGLE_K / (theta0.sin() * theta0.sin())
            };
            self.angles.push(AngleTerm {
                i,
                j,
                k,
                c,
                cos_theta0: theta0.cos(),
                linear,
            });
            self.exclusions.insert(i, k);
        }
    }

    /// Eq 13 with the barriers of Eqs 14-23, each divided by the number of (I,L) pairs sharing the
    /// central bond -- the paper's "2/9 for each of the nine possibilities of I and L".
    fn resolve_torsions(&mut self, topology: &MolecularTopology, adjacency: &[Vec<usize>]) {
        use std::collections::HashMap;

        // How many dihedrals share each central bond.
        let mut share_count: HashMap<(usize, usize), usize> = HashMap::new();
        for proper in &topology.propers {
            let (_, j, k, _) = proper.atom_ids;
            *share_count.entry(ordered(j, k)).or_insert(0) += 1;
        }

        let order_of: HashMap<(usize, usize), TopologyBondOrder> = topology
            .bonds
            .iter()
            .map(|b| (ordered(b.atom_ids.0, b.atom_ids.1), b.order))
            .collect();

        for proper in &topology.propers {
            let (i, j, k, l) = proper.atom_ids;
            let Some(&order) = order_of.get(&ordered(j, k)) else {
                continue;
            };

            let (class_i, class_j, class_k) = (
                torsion_class(&self.atom_types[i]),
                torsion_class(&self.atom_types[j]),
                torsion_class(&self.atom_types[k]),
            );
            let col16_j = is_column_16(topology.atoms[j].element);
            let col16_k = is_column_16(topology.atoms[k].element);

            let Some((mut v, n, mut phi0_deg)) =
                torsion_parameters(class_i, class_j, class_k, col16_j, col16_k, order)
            else {
                continue;
            };

            // (f), Eq 19: an *exocyclic* single bond between two aromatic atoms -- biphenyl's
            // central bond, or a phenyl ester -- doubles the barrier of Eq 18. Ring membership is
            // only known here, which is why this correction is not in `torsion_parameters`.
            if class_j == TorsionClass::Resonant
                && class_k == TorsionClass::Resonant
                && order == TopologyBondOrder::Single
                && !bond_is_in_ring(adjacency, j, k)
            {
                v = 10.0;
                phi0_deg = 180.0;
            }

            let shared = share_count.get(&ordered(j, k)).copied().unwrap_or(1).max(1);
            self.torsions.push(TorsionTerm {
                i,
                j,
                k,
                l,
                v: v / shared as f64,
                n,
                phi0: phi0_deg.to_radians(),
            });
        }
    }

    /// Eqs 28-29, expanded to the three permutations the paper asks for, each weighted 1/3.
    ///
    /// The typer emits impropers only at three-coordinate sp2 and resonant centres, which is
    /// exactly where DREIDING wants them: `X_3` has `K = 0`, and the paper notes that the
    /// inversion barriers of NH3 and PH3 come out well from the angle terms alone. The
    /// `C_31` case (ψ⁰ = 54.74°) belongs to the united-atom types, which never arise here because
    /// Beavyr always has explicit hydrogens.
    fn resolve_inversions(&mut self, topology: &MolecularTopology) {
        for improper in &topology.impropers {
            let (p1, p2, centre, p3) = improper.atom_ids;
            for (j, k, l) in [(p1, p2, p3), (p2, p3, p1), (p3, p1, p2)] {
                self.inversions.push(InversionTerm {
                    i: centre,
                    j,
                    k,
                    l,
                    c: params::INVERSION_K / 3.0,
                    cos_psi0: 1.0,
                    planar: true,
                });
            }
        }
    }

    fn resolve_vdw(&mut self, topology: &MolecularTopology) {
        self.vdw = topology
            .atoms
            .iter()
            .map(|atom| {
                let (r0, d0) = params::vdw(&atom.atom_type).expect("coverage checked");
                VdwParam { r0, d0 }
            })
            .collect();
    }

    /// Hydrogen bonds: an `H_HB` hydrogen, the heavy atom it is bonded to, and every N/O/F that
    /// could accept from it.
    ///
    /// Every candidate triple is enumerated once at build time rather than searched each step. The
    /// 12-10 term decays fast and the angular factor cuts off past 90°, so distant or badly
    /// oriented candidates contribute exactly zero.
    fn resolve_hbonds(&mut self, topology: &MolecularTopology, adjacency: &[Vec<usize>]) {
        let acceptors: Vec<usize> = topology
            .atoms
            .iter()
            .filter(|a| matches!(a.element, Element::N | Element::O | Element::F))
            .map(|a| a.id)
            .collect();

        for atom in &topology.atoms {
            if atom.atom_type != "H_HB" {
                continue;
            }
            let Some(&donor) = adjacency[atom.id].first() else {
                continue;
            };
            for &acceptor in &acceptors {
                // Its own donor, and anything already bonded to the hydrogen, are not acceptors.
                if acceptor == donor || adjacency[atom.id].contains(&acceptor) {
                    continue;
                }
                self.hbonds.push(HBondTerm {
                    donor,
                    hydrogen: atom.id,
                    acceptor,
                });
            }
        }
    }
}

/// Rejects a bond list that cannot be right, before it can quietly mistype every atom.
///
/// Only hydrogen is checked, because hydrogen is where a slightly generous distance threshold
/// does damage first and most visibly: two hydrogens 1.5 Å apart are a normal H-H contact in
/// water or methane, not a bond. The one legitimate two-coordinate hydrogen is DREIDING's `H_b`,
/// the bridging hydrogen of diborane, so a hydrogen bonded to two borons is allowed through.
fn check_connectivity_is_plausible(mol: &Molecule) -> Result<(), BuildError> {
    let mut degree = vec![0usize; mol.atoms.len()];
    let mut boron_neighbours = vec![0usize; mol.atoms.len()];
    for &(i, j, _) in &mol.bonds {
        degree[i] += 1;
        degree[j] += 1;
        if mol.atoms[j] == "B" {
            boron_neighbours[i] += 1;
        }
        if mol.atoms[i] == "B" {
            boron_neighbours[j] += 1;
        }
    }
    for (index, symbol) in mol.atoms.iter().enumerate() {
        if symbol != "H" || degree[index] <= 1 {
            continue;
        }
        // Diborane's bridging hydrogen genuinely has two bonds, both to boron.
        if degree[index] == 2 && boron_neighbours[index] == 2 {
            continue;
        }
        return Err(BuildError::ImplausibleConnectivity {
            index,
            symbol: symbol.clone(),
            degree: degree[index],
        });
    }
    Ok(())
}

/// Maps Beavyr's continuous perceived bond order onto the typer's discrete one.
///
/// A partial order (0.5, from `bond_order::PARTIAL_BOND_ORDER`) has no discrete equivalent and is
/// treated as single, which is the conservative reading.
fn discrete_order(order: f64) -> GraphBondOrder {
    if order >= 2.5 {
        GraphBondOrder::Triple
    } else if order >= 1.75 {
        GraphBondOrder::Double
    } else if order >= 1.25 {
        GraphBondOrder::Aromatic
    } else {
        GraphBondOrder::Single
    }
}

/// Eq 9: the multiplier on the single-bond force constant and dissociation energy.
fn order_multiplicity(order: TopologyBondOrder) -> f64 {
    match order {
        TopologyBondOrder::Single => 1.0,
        TopologyBondOrder::Resonant => 1.5,
        TopologyBondOrder::Double => 2.0,
        TopologyBondOrder::Triple => 3.0,
    }
}

fn ordered(i: usize, j: usize) -> (usize, usize) {
    if i < j {
        (i, j)
    } else {
        (j, i)
    }
}

fn build_adjacency(natoms: usize, topology: &MolecularTopology) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); natoms];
    for bond in &topology.bonds {
        let (i, j) = bond.atom_ids;
        adjacency[i].push(j);
        adjacency[j].push(i);
    }
    adjacency
}

/// Whether `j` and `k` are still connected without using the `j-k` bond -- that is, whether the
/// bond lies in a ring.
fn bond_is_in_ring(adjacency: &[Vec<usize>], j: usize, k: usize) -> bool {
    let mut seen = vec![false; adjacency.len()];
    let mut stack = vec![j];
    seen[j] = true;
    while let Some(current) = stack.pop() {
        for &next in &adjacency[current] {
            // Walking j->k directly would prove nothing; that is the bond in question.
            if (current == j && next == k) || (current == k && next == j) {
                continue;
            }
            if next == k {
                return true;
            }
            if !seen[next] {
                seen[next] = true;
                stack.push(next);
            }
        }
    }
    false
}

/// Single, acyclic bonds with a non-terminal atom at each end, plus the atoms each rotation moves.
///
/// For conformational search. Bonds of order 1.5 are excluded along with double and triple bonds:
/// a resonant bond is not free to rotate.
fn rotatable_bonds(
    topology: &MolecularTopology,
    adjacency: &[Vec<usize>],
) -> Vec<RotatableBond> {
    let mut out = Vec::new();
    for bond in &topology.bonds {
        if bond.order != TopologyBondOrder::Single {
            continue;
        }
        let (j, k) = bond.atom_ids;
        if adjacency[j].len() < 2 || adjacency[k].len() < 2 {
            continue;
        }
        if bond_is_in_ring(adjacency, j, k) {
            continue;
        }
        // Everything reachable from k without crossing back over the j-k bond.
        let mut seen = vec![false; adjacency.len()];
        seen[j] = true;
        seen[k] = true;
        let mut stack = vec![k];
        let mut moving = vec![k];
        while let Some(current) = stack.pop() {
            for &next in &adjacency[current] {
                if !seen[next] {
                    seen[next] = true;
                    moving.push(next);
                    stack.push(next);
                }
            }
        }
        out.push(RotatableBond { j, k, moving });
    }
    out
}

/// Flattens a typer error into one line, including its source, which carries the detail.
fn describe_typer_error(error: &typer::TyperError) -> String {
    use std::error::Error;
    match error.source() {
        Some(source) => format!("{error}: {source}"),
        None => error.to_string(),
    }
}


/// Bridge to the imported optimizer, including unit conversion.
pub mod objective;

#[cfg(test)]
mod tests {
    use super::objective::{cleanup_options, relax_angstrom};
    use super::*;

    /// Standard geometries, so bond perception sees sensible distances. Perception works from
    /// geometry, so a nonsense input would be mistyped before the force field ever ran.
    const WATER: &str = "3\n\nO 0.000 0.000 0.000\nH 0.960 0.000 0.000\nH -0.240 0.929 0.000\n";

    const ETHANE: &str = "8\n\n\
        C  0.000  0.000  0.000\n\
        C  0.000  0.000  1.530\n\
        H  1.028  0.000 -0.363\n\
        H -0.514  0.890 -0.363\n\
        H -0.514 -0.890 -0.363\n\
        H -1.028  0.000  1.893\n\
        H  0.514 -0.890  1.893\n\
        H  0.514  0.890  1.893\n";

    const AMMONIA: &str = "4\n\n\
        N  0.000  0.000  0.117\n\
        H  0.000  0.938 -0.274\n\
        H  0.812 -0.469 -0.274\n\
        H -0.812 -0.469 -0.274\n";

    /// Experimental thiophene: S-C 1.714, C2-C3 1.370, C3-C4 1.423, C-S-C 92.2 deg.
    const THIOPHENE: &str = "9\n\n\
        S  0.0000  0.0000  0.0000\n\
        C  1.2350  1.1885  0.0000\n\
        C  0.7130  2.4551  0.0000\n\
        C -0.7130  2.4551  0.0000\n\
        C -1.2350  1.1885  0.0000\n\
        H  2.3140  1.1405  0.0000\n\
        H  1.3403  3.3342  0.0000\n\
        H -1.3403  3.3342  0.0000\n\
        H -2.3140  1.1405  0.0000\n";

    /// A planar hexagon of carbons with radial hydrogens.
    fn benzene(pucker: f64) -> String {
        let mut xyz = String::from("12\n\n");
        for k in 0..6 {
            let a = (k as f64) * std::f64::consts::FRAC_PI_3;
            // Alternate the pucker so the ring starts non-planar.
            let z = if k % 2 == 0 { pucker } else { -pucker };
            xyz += &format!("C {:.4} {:.4} {:.4}\n", 1.39 * a.cos(), 1.39 * a.sin(), z);
        }
        for k in 0..6 {
            let a = (k as f64) * std::f64::consts::FRAC_PI_3;
            let z = if k % 2 == 0 { pucker } else { -pucker };
            xyz += &format!("H {:.4} {:.4} {:.4}\n", 2.48 * a.cos(), 2.48 * a.sin(), z);
        }
        xyz
    }

    /// `Molecule::from_xyz` hardcodes a very generous bond threshold; the running application
    /// uses `settings.bond_thresh_scale`, which defaults to 1.2. Tests use the app's value, since
    /// that is the connectivity a cleanup will actually be handed.
    fn molecule(xyz: &str) -> Molecule {
        let mut mol = Molecule::from_xyz(xyz);
        mol.recompute_bonds(1.2, 2.5);
        mol
    }

    fn flat(mol: &Molecule) -> Vec<f64> {
        mol.pos
            .iter()
            .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
            .collect()
    }

    fn distance(pos: &[f64], i: usize, j: usize) -> f64 {
        let d: Vec<f64> = (0..3).map(|c| pos[3 * i + c] - pos[3 * j + c]).collect();
        d.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    fn angle_deg(pos: &[f64], i: usize, j: usize, k: usize) -> f64 {
        let u: Vec<f64> = (0..3).map(|c| pos[3 * i + c] - pos[3 * j + c]).collect();
        let v: Vec<f64> = (0..3).map(|c| pos[3 * k + c] - pos[3 * j + c]).collect();
        let dot: f64 = u.iter().zip(&v).map(|(a, b)| a * b).sum();
        let nu = u.iter().map(|a| a * a).sum::<f64>().sqrt();
        let nv = v.iter().map(|a| a * a).sum::<f64>().sqrt();
        (dot / (nu * nv)).clamp(-1.0, 1.0).acos().to_degrees()
    }

    #[test]
    fn types_water_and_ethane_as_the_paper_does() {
        let water = DreidingTopology::build(&molecule(WATER)).unwrap();
        assert_eq!(water.atom_types()[0], "O_3");
        // Hydrogen on an atom with lone pairs is a hydrogen-bond donor.
        assert_eq!(water.atom_types()[1], "H_HB");

        let ethane = DreidingTopology::build(&molecule(ETHANE)).unwrap();
        assert_eq!(ethane.atom_types()[0], "C_3");
        assert_eq!(ethane.atom_types()[2], "H_");
    }

    #[test]
    fn enumerates_the_expected_terms_for_ethane() {
        let ethane = DreidingTopology::build(&molecule(ETHANE)).unwrap();
        assert_eq!(ethane.bonds.len(), 7); // 1 C-C + 6 C-H
        assert_eq!(ethane.angles.len(), 12); // 6 H-C-H + 6 H-C-C
        assert_eq!(ethane.torsions.len(), 9); // 9 H-C-C-H
        assert!(ethane.inversions.is_empty(), "sp3 centres take no inversion");
    }

    /// The paper divides the barrier among the dihedrals sharing a central bond: "2/9 for each of
    /// the nine possibilities of I and L". Without it ethane's barrier is nine times too large.
    #[test]
    fn ethane_torsion_barrier_is_divided_among_its_nine_dihedrals() {
        let ethane = DreidingTopology::build(&molecule(ETHANE)).unwrap();
        for t in &ethane.torsions {
            assert!(
                (t.v - 2.0 / 9.0).abs() < 1.0e-12,
                "barrier is {}, want 2/9",
                t.v
            );
            assert_eq!(t.n, 3.0);
        }
        // Summed over the nine, the barrier is the paper's 2 kcal/mol, not 18.
        let total: f64 = ethane.torsions.iter().map(|t| t.v).sum();
        assert!((total - 2.0).abs() < 1.0e-12);
    }

    /// The whole assembled gradient, not just the individual terms, against finite differences.
    /// This is what catches a mistake in how terms are combined rather than in a term itself.
    #[test]
    fn assembled_gradient_matches_finite_differences() {
        for (name, xyz) in [("water", WATER.to_string()), ("ethane", ETHANE.to_string()), ("benzene", benzene(0.05))] {
            let mol = molecule(&xyz);
            let topology = DreidingTopology::build(&mol).unwrap();
            let pos = flat(&mol);
            let mut work = Workspace::new();

            let mut forces = vec![0.0; pos.len()];
            topology.energy_and_forces(&pos, &mut work, &mut forces);

            const H: f64 = 1.0e-6;
            for c in 0..pos.len() {
                let mut plus = pos.clone();
                let mut minus = pos.clone();
                plus[c] += H;
                minus[c] -= H;
                let mut scratch = vec![0.0; pos.len()];
                let e_plus = topology.energy_and_forces(&plus, &mut work, &mut scratch);
                let e_minus = topology.energy_and_forces(&minus, &mut work, &mut scratch);
                // forces = -dE/dx
                let numeric = -(e_plus - e_minus) / (2.0 * H);
                let scale = forces[c].abs().max(numeric.abs()).max(1.0);
                assert!(
                    (forces[c] - numeric).abs() / scale < 1.0e-5,
                    "{name}: coordinate {c} analytic {} vs numeric {numeric}",
                    forces[c]
                );
            }
        }
    }

    /// Table I predicts these directly: Rₑ = R⁰_I + R⁰_J − 0.01, and θ⁰ from the central atom.
    /// Water has one angle and no unexcluded nonbonded pair, so it must relax to exactly those.
    #[test]
    fn water_relaxes_to_its_table_i_geometry() {
        let mol = molecule(WATER);
        let topology = DreidingTopology::build(&mol).unwrap();
        let (relaxed, report) =
            relax_angstrom(&topology, &flat(&mol), cleanup_options()).unwrap();

        // O_3 0.660 + H_HB 0.330 − 0.01
        let expected_oh = 0.660 + 0.330 - 0.01;
        assert!(
            (distance(&relaxed, 0, 1) - expected_oh).abs() < 2.0e-3,
            "O-H is {}, want {expected_oh}",
            distance(&relaxed, 0, 1)
        );
        assert!(
            (angle_deg(&relaxed, 1, 0, 2) - 104.51).abs() < 0.5,
            "H-O-H is {}, want 104.51",
            angle_deg(&relaxed, 1, 0, 2)
        );
        assert!(report.converged, "{}", report.message);
        assert!(report.final_energy <= report.initial_energy);
    }

    #[test]
    fn ethane_relaxes_to_the_expected_bond_lengths() {
        let mol = molecule(ETHANE);
        let topology = DreidingTopology::build(&mol).unwrap();
        let (relaxed, report) =
            relax_angstrom(&topology, &flat(&mol), cleanup_options()).unwrap();

        // C_3 0.770 twice, less 0.01: the paper's own worked value of 1.53 Å.
        assert!(
            (distance(&relaxed, 0, 1) - 1.53).abs() < 0.02,
            "C-C is {}, want 1.53",
            distance(&relaxed, 0, 1)
        );
        // C_3 0.770 + H_ 0.330 − 0.01
        assert!(
            (distance(&relaxed, 0, 2) - 1.09).abs() < 0.02,
            "C-H is {}, want 1.09",
            distance(&relaxed, 0, 2)
        );
        assert!(report.converged, "{}", report.message);
    }

    /// Ammonia must stay pyramidal. DREIDING gives `X_3` centres no inversion term at all, and the
    /// paper claims the angle terms alone describe the NH3 inversion barrier well -- so this tests
    /// that claim holds in the implementation.
    #[test]
    fn ammonia_stays_pyramidal() {
        let mol = molecule(AMMONIA);
        let topology = DreidingTopology::build(&mol).unwrap();
        assert_eq!(topology.atom_types()[0], "N_3");
        let (relaxed, _) = relax_angstrom(&topology, &flat(&mol), cleanup_options()).unwrap();

        // N_3's θ⁰ is 106.7°; a planar centre would be 120°.
        for (i, j) in [(1, 2), (2, 3), (1, 3)] {
            let a = angle_deg(&relaxed, i, 0, j);
            assert!(
                (a - 106.7).abs() < 1.5,
                "H-N-H is {a}, want 106.7 (flattening would give 120)"
            );
        }
    }

    /// Benzene must flatten and equalise. This is the inversion term's job: the paper points out
    /// that angle terms alone "will not in general lead to the proper restoring force toward the
    /// planar configuration", which is why an explicit four-body term exists.
    #[test]
    fn puckered_benzene_flattens_with_equal_bonds() {
        let mol = molecule(&benzene(0.10));
        let topology = DreidingTopology::build(&mol).unwrap();
        assert_eq!(topology.atom_types()[0], "C_R");
        assert_eq!(topology.inversions.len(), 18, "6 centres x 3 permutations");

        let start = flat(&mol);
        let start_spread = (0..6)
            .map(|k| start[3 * k + 2].abs())
            .fold(0.0_f64, f64::max);
        assert!(start_spread > 0.05, "the test input was not actually puckered");

        let (relaxed, report) = relax_angstrom(&topology, &start, cleanup_options()).unwrap();

        let end_spread = (0..6)
            .map(|k| relaxed[3 * k + 2].abs())
            .fold(0.0_f64, f64::max);
        assert!(
            end_spread < start_spread * 0.5,
            "ring did not flatten: out-of-plane spread went from {start_spread} to {end_spread}"
        );

        // C_R 0.700 twice, less 0.01.
        for k in 0..6 {
            let d = distance(&relaxed, k, (k + 1) % 6);
            assert!((d - 1.39).abs() < 0.03, "C-C {k} is {d}, want 1.39");
        }
        assert!(report.final_energy < report.initial_energy);
    }

    /// Aromatic sulfur is the coverage gap that bites in practice. It must refuse by name rather
    /// than invent a radius and hand back a confident, wrong geometry.
    #[test]
    fn thiophene_is_refused_by_name() {
        let mol = molecule(THIOPHENE);
        let error = DreidingTopology::build(&mol).expect_err("S_R has no published parameters");
        match &error {
            BuildError::Unparameterized {
                symbol, atom_type, ..
            } => {
                assert_eq!(symbol, "S");
                assert_eq!(atom_type, "S_R");
            }
            other => panic!("expected an unparameterised-type refusal, got {other:?}"),
        }
        // The message must name the atom and say what is missing, not just fail.
        let text = error.to_string();
        assert!(text.contains("S_R"), "{text}");
        assert!(text.contains("thiophene"), "{text}");
    }

    #[test]
    fn an_unparameterised_metal_is_refused() {
        // Lithium fluoride: the typer types Li, the paper does not parameterise it.
        let mol = molecule("2\n\nLi 0.000 0.000 0.000\nF 0.000 0.000 1.564\n");
        let error = DreidingTopology::build(&mol).expect_err("Li has no published parameters");
        assert!(matches!(error, BuildError::Unparameterized { .. }), "{error:?}");
    }

    #[test]
    fn an_unknown_element_is_reported_not_skipped() {
        let mol = molecule("2\n\nXx 0.000 0.000 0.000\nH 0.000 0.000 1.000\n");
        match DreidingTopology::build(&mol) {
            Err(BuildError::UnknownElement { symbol, .. }) => assert_eq!(symbol, "Xx"),
            other => panic!("expected an unknown-element error, got {other:?}"),
        }
    }

    #[test]
    fn an_empty_structure_is_refused() {
        assert!(matches!(
            DreidingTopology::build(&Molecule::empty()),
            Err(BuildError::Empty)
        ));
    }

    /// Rotatable bonds drive conformational search. Ethane's C-C is the only one; its C-H bonds
    /// end in a terminal atom and a ring bond would be excluded.
    #[test]
    fn finds_ethanes_single_rotatable_bond() {
        let ethane = DreidingTopology::build(&molecule(ETHANE)).unwrap();
        assert_eq!(ethane.rotatable_bonds().len(), 1);
        let bond = &ethane.rotatable_bonds()[0];
        assert_eq!(bond.moving.len(), 4, "one carbon and its three hydrogens");

        // Benzene's ring bonds are all cyclic, so none rotate.
        let benzene = DreidingTopology::build(&molecule(&benzene(0.0))).unwrap();
        assert!(benzene.rotatable_bonds().is_empty());
    }
}
