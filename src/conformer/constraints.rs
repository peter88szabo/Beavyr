//! Geometric constraints, for searching the conformers of a transition state.
//!
//! A transition state is a saddle point, so it cannot simply be relaxed: let the optimiser go and
//! it slides downhill to reactant or product. The usual way to explore its conformers is to hold
//! the reacting part fixed and search everything else -- fix the forming and breaking bond
//! lengths, and let the rest of the molecule find its shapes.
//!
//! # What to fix
//!
//! The motif is `A-X-B`: an atom `X` between the partner it is leaving and the one it is joining.
//! A hydrogen in flight, the carbon of an SN2 displacement, a metal centre exchanging a ligand.
//! Three levels of restriction are offered, because how much needs holding depends on how far the
//! rest of the molecule can wander before the reacting part follows it:
//!
//! 1. [`ActiveSiteLevel::Bonds`] -- the two lengths `A-X` and `X-B`. The least that works.
//! 2. [`ActiveSiteLevel::BondsAndAngle`] -- and the `A-X-B` angle. Needed once the geometry can
//!    bend at `X`, which it will for anything but a near-linear transfer.
//! 3. [`ActiveSiteLevel::BondsAngleAndDihedrals`] -- and the dihedrals that orient the `A-X-B`
//!    unit against the rest of the molecule, so the reacting group cannot rotate away from the
//!    face it is attacking.
//!
//! Several sites can be constrained at once -- a concerted reaction has more than one -- and each
//! is added on its own, so the user builds up exactly the set they want.
//!
//! # The value is read from the geometry
//!
//! A constraint fixes a coordinate at whatever the current structure has, rather than at a number
//! typed in. That is what makes it usable: the structure came from a transition-state search, so
//! its active bond lengths are the ones worth holding, and asking a chemist to retype them to four
//! decimal places would only introduce error.

use std::f64::consts::PI;

use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
use crate::optimizer::constrained::{ConstraintCoordinate, ConstraintTarget};

/// How much of an `A-X-B` active site to hold fixed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ActiveSiteLevel {
    /// The two bond lengths, `A-X` and `X-B`.
    #[default]
    Bonds,
    /// The two bond lengths and the `A-X-B` angle.
    BondsAndAngle,
    /// The bond lengths, the angle, and the dihedrals orienting the site.
    BondsAngleAndDihedrals,
}

impl ActiveSiteLevel {
    pub const ALL: [ActiveSiteLevel; 3] = [
        ActiveSiteLevel::Bonds,
        ActiveSiteLevel::BondsAndAngle,
        ActiveSiteLevel::BondsAngleAndDihedrals,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Self::Bonds => "Bond lengths",
            Self::BondsAndAngle => "Lengths + angle",
            Self::BondsAngleAndDihedrals => "Lengths + angle + dihedrals",
        }
    }

    pub fn description(self) -> &'static str {
        match self {
            Self::Bonds => "Holds A-X and X-B. The least that keeps a saddle point from sliding.",
            Self::BondsAndAngle => {
                "Also holds the A-X-B angle, which anything but a near-linear transfer needs."
            }
            Self::BondsAngleAndDihedrals => {
                "Also holds the site's orientation, so the reacting group cannot rotate off the \
                 face it is attacking."
            }
        }
    }
}

/// One fixed coordinate, with the value taken from the geometry and a description for the user.
#[derive(Debug, Clone, PartialEq)]
pub struct Constraint {
    pub coordinate: ConstraintCoordinate,
    /// bohr for a length, radians for an angle or dihedral -- the units the optimiser works in.
    pub value: f64,
}

impl Constraint {
    /// The atoms this holds, in the order the coordinate names them.
    pub fn atoms(&self) -> Vec<usize> {
        match self.coordinate {
            ConstraintCoordinate::Bond(a) => a.to_vec(),
            ConstraintCoordinate::Angle(a) => a.to_vec(),
            ConstraintCoordinate::Dihedral(a) => a.to_vec(),
        }
    }

    /// A line for the user: what is held, between which atoms, at what value.
    ///
    /// Lengths in Å and angles in degrees, whatever the optimiser uses internally -- nobody reads
    /// a bond length in bohr.
    pub fn describe(&self, symbols: &[String]) -> String {
        let name = |index: usize| match symbols.get(index) {
            Some(symbol) => format!("{symbol}{}", index + 1),
            None => format!("atom {}", index + 1),
        };
        match self.coordinate {
            ConstraintCoordinate::Bond([a, b]) => format!(
                "{} - {} at {:.3} Å",
                name(a),
                name(b),
                self.value * BOHR_TO_ANGSTROM
            ),
            ConstraintCoordinate::Angle([a, b, c]) => format!(
                "{} - {} - {} at {:.1}°",
                name(a),
                name(b),
                name(c),
                self.value.to_degrees()
            ),
            ConstraintCoordinate::Dihedral([a, b, c, d]) => format!(
                "{} - {} - {} - {} at {:.1}°",
                name(a),
                name(b),
                name(c),
                name(d),
                self.value.to_degrees()
            ),
        }
    }

    /// The form the constrained optimiser wants.
    pub fn target(&self) -> ConstraintTarget {
        ConstraintTarget {
            coordinate: self.coordinate,
            target: self.value,
        }
    }
}

/// Reads a bond length from a geometry, in bohr.
pub fn measure_bond(coords_bohr: &[f64], a: usize, b: usize) -> f64 {
    (0..3)
        .map(|k| (coords_bohr[3 * a + k] - coords_bohr[3 * b + k]).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Reads an angle from a geometry, in radians.
pub fn measure_angle(coords_bohr: &[f64], a: usize, b: usize, c: usize) -> f64 {
    let u: Vec<f64> = (0..3)
        .map(|k| coords_bohr[3 * a + k] - coords_bohr[3 * b + k])
        .collect();
    let v: Vec<f64> = (0..3)
        .map(|k| coords_bohr[3 * c + k] - coords_bohr[3 * b + k])
        .collect();
    let dot: f64 = u.iter().zip(&v).map(|(x, y)| x * y).sum();
    let nu = u.iter().map(|x| x * x).sum::<f64>().sqrt();
    let nv = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if nu < 1.0e-12 || nv < 1.0e-12 {
        return 0.0;
    }
    (dot / (nu * nv)).clamp(-1.0, 1.0).acos()
}

/// Reads a dihedral from a geometry, in radians, signed.
pub fn measure_dihedral(coords_bohr: &[f64], a: usize, b: usize, c: usize, d: usize) -> f64 {
    let at = |i: usize| {
        [
            coords_bohr[3 * i],
            coords_bohr[3 * i + 1],
            coords_bohr[3 * i + 2],
        ]
    };
    let sub = |p: [f64; 3], q: [f64; 3]| [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
    let cross = |p: [f64; 3], q: [f64; 3]| {
        [
            p[1] * q[2] - p[2] * q[1],
            p[2] * q[0] - p[0] * q[2],
            p[0] * q[1] - p[1] * q[0],
        ]
    };
    let dot = |p: [f64; 3], q: [f64; 3]| p[0] * q[0] + p[1] * q[1] + p[2] * q[2];

    let b1 = sub(at(b), at(a));
    let b2 = sub(at(c), at(b));
    let b3 = sub(at(d), at(c));
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let lb2 = dot(b2, b2).sqrt();
    if lb2 < 1.0e-12 {
        return 0.0;
    }
    let m = dot(cross(n1, n2), [b2[0] / lb2, b2[1] / lb2, b2[2] / lb2]);
    m.atan2(dot(n1, n2))
}

/// Builds the constraints for one `A-X-B` active site at the requested level.
///
/// `adjacency` is needed only for [`ActiveSiteLevel::BondsAngleAndDihedrals`], to find a
/// substituent on each of `A` and `B` to hang the orienting dihedrals from. A site whose ends
/// carry nothing -- a bare atom or ion -- simply gets no dihedral, rather than an error: the
/// orientation of a lone atom is not a thing to fix.
pub fn active_site(
    coords_bohr: &[f64],
    a: usize,
    x: usize,
    b: usize,
    level: ActiveSiteLevel,
    adjacency: &[Vec<usize>],
) -> Vec<Constraint> {
    let mut out = vec![
        Constraint {
            coordinate: ConstraintCoordinate::Bond([a, x]),
            value: measure_bond(coords_bohr, a, x),
        },
        Constraint {
            coordinate: ConstraintCoordinate::Bond([x, b]),
            value: measure_bond(coords_bohr, x, b),
        },
    ];
    if matches!(
        level,
        ActiveSiteLevel::BondsAndAngle | ActiveSiteLevel::BondsAngleAndDihedrals
    ) {
        out.push(Constraint {
            coordinate: ConstraintCoordinate::Angle([a, x, b]),
            value: measure_angle(coords_bohr, a, x, b),
        });
    }
    if level == ActiveSiteLevel::BondsAngleAndDihedrals {
        // One dihedral per end, each hanging off a substituent of that end. Together they pin the
        // site's orientation without over-determining it.
        if let Some(w) = substituent(adjacency, a, &[x, b]) {
            out.push(Constraint {
                coordinate: ConstraintCoordinate::Dihedral([w, a, x, b]),
                value: measure_dihedral(coords_bohr, w, a, x, b),
            });
        }
        if let Some(z) = substituent(adjacency, b, &[x, a]) {
            out.push(Constraint {
                coordinate: ConstraintCoordinate::Dihedral([a, x, b, z]),
                value: measure_dihedral(coords_bohr, a, x, b, z),
            });
        }
    }
    out
}

/// A neighbour of `atom` that is not one of `avoid`, preferring a heavy atom.
///
/// Heavy first because a dihedral to a hydrogen is a weaker handle on orientation -- a methyl's
/// hydrogens spin freely, so pinning one of them pins the methyl's rotation as much as the site's
/// orientation.
fn substituent(adjacency: &[Vec<usize>], atom: usize, avoid: &[usize]) -> Option<usize> {
    let candidates = adjacency.get(atom)?;
    candidates
        .iter()
        .copied()
        .find(|n| !avoid.contains(n) && adjacency.get(*n).is_some_and(|a| a.len() > 1))
        .or_else(|| candidates.iter().copied().find(|n| !avoid.contains(n)))
}

/// Whether a constraint and a search degree of freedom would fight over the same coordinate.
///
/// A fixed dihedral must not also be sampled, and a bond held at a length must not have the search
/// twisting the torsion about it: the constrained optimiser would pull the geometry back every
/// time, wasting the whole generation. The clash is reported so the caller can drop the degree of
/// freedom instead.
pub fn clashes_with_central_bond(constraint: &Constraint, j: usize, k: usize) -> bool {
    let same_pair = |p: usize, q: usize| (p == j && q == k) || (p == k && q == j);
    match constraint.coordinate {
        ConstraintCoordinate::Bond([a, b]) => same_pair(a, b),
        // An angle pins the bond pair at its centre, so a torsion about either is constrained too.
        ConstraintCoordinate::Angle([a, b, c]) => same_pair(a, b) || same_pair(b, c),
        // A dihedral's own central bond is exactly the torsion it fixes.
        ConstraintCoordinate::Dihedral([_, b, c, _]) => same_pair(b, c),
    }
}

/// Whether a set of constraints is self-consistent enough to optimise against.
///
/// The constrained optimiser needs its coordinates to be independent; the same coordinate listed
/// twice makes the system singular. Duplicates are the easy mistake -- adding the same active site
/// twice, or a bond that a previous site already held -- so they are caught here with a message
/// rather than in a linear solve.
pub fn check(constraints: &[Constraint], symbols: &[String]) -> Result<(), String> {
    for (index, constraint) in constraints.iter().enumerate() {
        for other in &constraints[index + 1..] {
            if same_coordinate(constraint, other) {
                return Err(format!(
                    "{} is fixed twice. Remove one of them.",
                    constraint.describe(symbols)
                ));
            }
        }
        let atoms = constraint.atoms();
        let mut sorted = atoms.clone();
        sorted.sort_unstable();
        sorted.dedup();
        if sorted.len() != atoms.len() {
            return Err(format!(
                "{} names the same atom twice.",
                constraint.describe(symbols)
            ));
        }
    }
    Ok(())
}

/// Whether two constraints hold the same coordinate, in either direction.
fn same_coordinate(left: &Constraint, right: &Constraint) -> bool {
    use ConstraintCoordinate::{Angle, Bond, Dihedral};
    match (left.coordinate, right.coordinate) {
        (Bond(a), Bond(b)) => a == b || a == [b[1], b[0]],
        (Angle(a), Angle(b)) => a == b || a == [b[2], b[1], b[0]],
        (Dihedral(a), Dihedral(b)) => a == b || a == [b[3], b[2], b[1], b[0]],
        _ => false,
    }
}

/// Radians, wrapped into (-π, π].
pub fn wrap(angle: f64) -> f64 {
    let wrapped = angle.rem_euclid(2.0 * PI);
    if wrapped > PI {
        wrapped - 2.0 * PI
    } else {
        wrapped
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A bent A-X-B site: H in flight between two oxygens, with a substituent on each end.
    /// Coordinates in bohr.
    fn transfer_site() -> (Vec<f64>, Vec<Vec<usize>>, Vec<String>) {
        // 0: C (substituent on A)   1: O (A)   2: H (X)   3: O (B)   4: C (substituent on B)
        let angstrom: Vec<f64> = vec![
            -1.40, 0.60, 0.00, // C
            -0.60, 0.00, 0.00, // O
            0.00, 0.55, 0.00, // H, off the O-O axis so the angle is not 180
            0.75, 0.00, 0.00, // O
            1.60, 0.70, 0.30, // C
        ];
        let bohr = angstrom.iter().map(|a| a / BOHR_TO_ANGSTROM).collect();
        let adjacency = vec![vec![1], vec![0, 2], vec![1, 3], vec![2, 4], vec![3]];
        let symbols = ["C", "O", "H", "O", "C"].iter().map(|s| s.to_string()).collect();
        (bohr, adjacency, symbols)
    }

    #[test]
    fn the_first_level_fixes_the_two_bond_lengths() {
        let (coords, adjacency, _) = transfer_site();
        let set = active_site(&coords, 1, 2, 3, ActiveSiteLevel::Bonds, &adjacency);
        assert_eq!(set.len(), 2);
        assert!(matches!(set[0].coordinate, ConstraintCoordinate::Bond([1, 2])));
        assert!(matches!(set[1].coordinate, ConstraintCoordinate::Bond([2, 3])));

        // And at the lengths the structure actually has.
        let expected = measure_bond(&coords, 1, 2);
        assert!((set[0].value - expected).abs() < 1.0e-12);
    }

    #[test]
    fn the_second_level_adds_the_angle() {
        let (coords, adjacency, _) = transfer_site();
        let set = active_site(&coords, 1, 2, 3, ActiveSiteLevel::BondsAndAngle, &adjacency);
        assert_eq!(set.len(), 3);
        assert!(matches!(
            set[2].coordinate,
            ConstraintCoordinate::Angle([1, 2, 3])
        ));
        // The site is bent, so the angle must not come back as 180 degrees.
        let degrees = set[2].value.to_degrees();
        assert!(
            degrees > 60.0 && degrees < 175.0,
            "angle came out {degrees}, which does not match a bent site"
        );
    }

    #[test]
    fn the_third_level_adds_one_dihedral_at_each_end() {
        let (coords, adjacency, _) = transfer_site();
        let set = active_site(
            &coords,
            1,
            2,
            3,
            ActiveSiteLevel::BondsAngleAndDihedrals,
            &adjacency,
        );
        assert_eq!(set.len(), 5, "two bonds, one angle, two dihedrals");
        let dihedrals: Vec<_> = set
            .iter()
            .filter(|c| matches!(c.coordinate, ConstraintCoordinate::Dihedral(_)))
            .collect();
        assert_eq!(dihedrals.len(), 2);
        // Each hangs off a substituent, so it reaches outside the A-X-B unit.
        assert!(matches!(
            dihedrals[0].coordinate,
            ConstraintCoordinate::Dihedral([0, 1, 2, 3])
        ));
        assert!(matches!(
            dihedrals[1].coordinate,
            ConstraintCoordinate::Dihedral([1, 2, 3, 4])
        ));
    }

    /// A site whose ends carry nothing gets no dihedral rather than an error: the orientation of a
    /// lone atom is not a thing to fix.
    #[test]
    fn a_bare_site_gets_no_dihedrals() {
        // Three atoms, nothing hanging off either end.
        let coords: Vec<f64> = vec![-1.4, 0.0, 0.0, 0.0, 0.3, 0.0, 1.4, 0.0, 0.0];
        let adjacency = vec![vec![1], vec![0, 2], vec![1]];
        let set = active_site(
            &coords,
            0,
            1,
            2,
            ActiveSiteLevel::BondsAngleAndDihedrals,
            &adjacency,
        );
        assert_eq!(set.len(), 3, "two bonds and the angle, no dihedral: {set:?}");
    }

    /// Measured values must agree with the geometry, since everything downstream holds them.
    #[test]
    fn the_measurements_are_right() {
        // A right angle at the origin, in bohr.
        let coords: Vec<f64> = vec![2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 3.0, 4.0];
        assert!((measure_bond(&coords, 0, 1) - 2.0).abs() < 1.0e-12);
        assert!((measure_angle(&coords, 0, 1, 2).to_degrees() - 90.0).abs() < 1.0e-9);
        // 0-1-2-3 with atom 3 lifted straight out of the 0-1-2 plane: a 90 degree dihedral.
        assert!((measure_dihedral(&coords, 0, 1, 2, 3).to_degrees().abs() - 90.0).abs() < 1.0e-9);
    }

    #[test]
    fn a_duplicated_constraint_is_reported() {
        let (coords, adjacency, symbols) = transfer_site();
        let mut set = active_site(&coords, 1, 2, 3, ActiveSiteLevel::Bonds, &adjacency);
        assert!(check(&set, &symbols).is_ok());

        // Adding the same site again duplicates every coordinate.
        set.extend(active_site(&coords, 1, 2, 3, ActiveSiteLevel::Bonds, &adjacency));
        let problem = check(&set, &symbols).expect_err("a duplicate must be caught");
        assert!(problem.contains("fixed twice"), "{problem}");

        // The reversed bond is the same coordinate, not a new one.
        let reversed = vec![
            Constraint {
                coordinate: ConstraintCoordinate::Bond([1, 2]),
                value: 1.0,
            },
            Constraint {
                coordinate: ConstraintCoordinate::Bond([2, 1]),
                value: 1.0,
            },
        ];
        assert!(check(&reversed, &symbols).is_err());
    }

    /// A constrained coordinate must remove the matching search degree of freedom, or the search
    /// twists a torsion the optimiser then pulls straight back.
    #[test]
    fn a_constraint_claims_the_torsion_it_fixes() {
        let bond = Constraint {
            coordinate: ConstraintCoordinate::Bond([1, 2]),
            value: 1.0,
        };
        assert!(clashes_with_central_bond(&bond, 1, 2));
        assert!(clashes_with_central_bond(&bond, 2, 1), "either direction");
        assert!(!clashes_with_central_bond(&bond, 3, 4));

        let angle = Constraint {
            coordinate: ConstraintCoordinate::Angle([1, 2, 3]),
            value: 1.0,
        };
        // Both bonds of the angle are claimed.
        assert!(clashes_with_central_bond(&angle, 1, 2));
        assert!(clashes_with_central_bond(&angle, 2, 3));
        assert!(!clashes_with_central_bond(&angle, 1, 3));

        let dihedral = Constraint {
            coordinate: ConstraintCoordinate::Dihedral([0, 1, 2, 3]),
            value: 1.0,
        };
        // Its central bond is the torsion it fixes.
        assert!(clashes_with_central_bond(&dihedral, 1, 2));
        assert!(!clashes_with_central_bond(&dihedral, 0, 1));
    }

    #[test]
    fn every_level_is_labelled_and_described() {
        assert_eq!(ActiveSiteLevel::ALL.len(), 3);
        for level in ActiveSiteLevel::ALL {
            assert!(!level.label().is_empty());
            assert!(!level.description().is_empty());
        }
        assert_eq!(ActiveSiteLevel::default(), ActiveSiteLevel::Bonds);
    }

    #[test]
    fn a_constraint_describes_itself_in_chemists_units() {
        let (coords, adjacency, symbols) = transfer_site();
        let set = active_site(&coords, 1, 2, 3, ActiveSiteLevel::BondsAndAngle, &adjacency);
        let bond = set[0].describe(&symbols);
        // Å not bohr, and one-based atom numbers.
        assert!(bond.contains("O2"), "{bond}");
        assert!(bond.contains("H3"), "{bond}");
        assert!(bond.contains('Å'), "{bond}");
        let angle = set[2].describe(&symbols);
        assert!(angle.contains('°'), "{angle}");
    }
}
