//! Scan coordinates, fragment displacement, and internal-coordinate preparation.
use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
use crate::optimizer::constrained::ConstraintCoordinate;
use crate::optimizer::internal_coords::{
    analyze_structure, internal_coords_from_bonds, ConnectivityModel, InternalCoords,
};
use crate::qchem_interfaces::job::scan::{BondMotion, Coordinate, Options, Prepared};
use bevy::math::{DQuat, DVec3};

fn v(x: &[f64], i: usize) -> DVec3 {
    DVec3::new(x[3 * i], x[3 * i + 1], x[3 * i + 2])
}
fn put(x: &mut [f64], i: usize, p: DVec3) {
    x[3 * i..3 * i + 3].copy_from_slice(&p.to_array());
}
pub(crate) fn measure(c: Coordinate, a: &[usize], x: &[f64]) -> Result<f64, String> {
    if a.len() != c.count() || a.iter().any(|&i| i * 3 + 2 >= x.len()) {
        return Err("Invalid scan atom indices.".into());
    }
    let value = match c {
        Coordinate::Bond => v(x, a[0]).distance(v(x, a[1])),
        Coordinate::Angle => {
            let u = v(x, a[0]) - v(x, a[1]);
            let w = v(x, a[2]) - v(x, a[1]);
            if u.length() < 1e-12 || w.length() < 1e-12 {
                return Err("Angle contains coincident atoms.".into());
            }
            (u.dot(w) / (u.length() * w.length()))
                .clamp(-1., 1.)
                .acos()
                .to_degrees()
        }
        Coordinate::Dihedral => {
            let u = v(x, a[1]) - v(x, a[0]);
            let w = v(x, a[2]) - v(x, a[1]);
            let z = v(x, a[3]) - v(x, a[2]);
            let n = u.cross(w);
            let m = w.cross(z);
            if n.length() < 1e-12 || m.length() < 1e-12 || w.length() < 1e-12 {
                return Err(
                    "A dihedral cannot be defined by collinear or coincident atoms.".into(),
                );
            }
            (w.normalize().dot(n.cross(m))).atan2(n.dot(m)).to_degrees()
        }
    };
    if value.is_finite() {
        Ok(value)
    } else {
        Err("Non-finite scan coordinate.".into())
    }
}
pub(super) fn residual(c: Coordinate, actual: f64, target: f64) -> f64 {
    if c == Coordinate::Dihedral {
        (actual - target + 180.).rem_euclid(360.) - 180.
    } else {
        actual - target
    }
}
fn model(atoms: &[String]) -> ConnectivityModel {
    ConnectivityModel {
        rcov_disp: atoms
            .iter()
            .map(|a| {
                crate::molecule::covalent_radius_angstrom(a) as f64 * (4. / 3.) / BOHR_TO_ANGSTROM
            })
            .collect(),
        kcn: 7.5,
        facmin: 0.6,
    }
}
pub(super) fn internals(
    atoms: &[String],
    x_ang: &[f64],
    options: &Options,
    prepared: &Prepared,
) -> Result<InternalCoords, String> {
    let x: Vec<f64> = x_ang.iter().map(|v| v / BOHR_TO_ANGSTROM).collect();
    let original = analyze_structure(&x, &model(atoms)).map_err(|e| e.to_string())?;
    let mut bonds = original.bonds;
    let targets = std::iter::once(options.coordinate.target(&prepared.atoms, options.start))
        .chain(prepared.fixed.iter().cloned());
    for target in targets {
        let a = match target.coordinate {
            ConstraintCoordinate::Bond(a) => a.to_vec(),
            ConstraintCoordinate::Angle(a) => a.to_vec(),
            ConstraintCoordinate::Dihedral(a) => a.to_vec(),
        };
        for pair in a.windows(2) {
            if !bonds
                .iter()
                .any(|b| *b == [pair[0], pair[1]] || *b == [pair[1], pair[0]])
            {
                bonds.push([pair[0], pair[1]]);
            }
        }
    }
    internal_coords_from_bonds(atoms.len(), &bonds).map_err(|e| e.to_string())
}
fn moving_atoms(
    atoms: &[String],
    x: &[f64],
    c: Coordinate,
    selected: &[usize],
    fragment: bool,
) -> Result<Vec<usize>, String> {
    if !fragment {
        return Ok(vec![*selected.last().unwrap()]);
    }
    let (fixed, moving) = match c {
        Coordinate::Bond => (selected[0], selected[1]),
        Coordinate::Angle => (selected[1], selected[2]),
        Coordinate::Dihedral => (selected[1], selected[2]),
    };
    let bohr: Vec<_> = x.iter().map(|v| v / BOHR_TO_ANGSTROM).collect();
    let ic = analyze_structure(&bohr, &model(atoms)).map_err(|e| e.to_string())?;
    let mut found = vec![moving];
    let mut cursor = 0;
    while cursor < found.len() {
        let atom = found[cursor];
        cursor += 1;
        for &[a, b] in &ic.bonds {
            if (a == fixed && b == moving) || (b == fixed && a == moving) {
                continue;
            }
            let next = if a == atom {
                b
            } else if b == atom {
                a
            } else {
                continue;
            };
            if !found.contains(&next) {
                found.push(next);
            }
        }
    }
    if found.contains(&fixed)
        || selected[..selected.len() - 1]
            .iter()
            .any(|a| *a != moving && found.contains(a))
    {
        return Err("The selected bond belongs to a ring or reconnecting fragment. Disable ‘Move connected fragment’, or use a relaxed scan.".into());
    }
    if !found.contains(selected.last().unwrap()) {
        return Err("The final selected atom is not in the moving fragment.".into());
    }
    Ok(found)
}
pub(super) fn rigid(
    atoms: &[String],
    original: &[f64],
    options: &Options,
    prepared: &Prepared,
    target: f64,
) -> Result<Vec<f64>, String> {
    let a = &prepared.atoms;
    let current = measure(options.coordinate, a, original)?;
    let moving = moving_atoms(
        atoms,
        original,
        options.coordinate,
        a,
        options.move_fragment,
    )?;
    let mut x = original.to_vec();
    if options.coordinate == Coordinate::Bond {
        let axis = v(original, a[1]) - v(original, a[0]);
        if axis.length() < 1e-12 {
            return Err("Cannot scan a bond between coincident atoms.".into());
        }
        let shift = axis.normalize() * (target - current);
        let fraction = if options.bond_motion == BondMotion::Symmetric {
            0.5
        } else {
            1.0
        };
        if options.bond_motion != BondMotion::First {
            for &i in &moving {
                put(&mut x, i, v(original, i) + shift * fraction);
            }
        }
        if options.bond_motion != BondMotion::Second {
            let first = moving_atoms(
                atoms,
                original,
                Coordinate::Bond,
                &[a[1], a[0]],
                options.move_fragment,
            )?;
            for i in first {
                put(&mut x, i, v(original, i) - shift * fraction);
            }
        }
    } else {
        let origin = v(original, a[1]);
        let axis = if options.coordinate == Coordinate::Angle {
            (v(original, a[0]) - origin).cross(v(original, a[2]) - origin)
        } else {
            v(original, a[2]) - origin
        };
        if axis.length() < 1e-12 {
            return Err("Cannot rotate about an undefined scan axis.".into());
        }
        let delta = residual(options.coordinate, target, current).to_radians();
        for sign in [1., -1.] {
            let rotation = DQuat::from_axis_angle(axis.normalize(), sign * delta);
            for &i in &moving {
                put(&mut x, i, origin + rotation * (v(original, i) - origin));
            }
            if residual(
                options.coordinate,
                measure(options.coordinate, a, &x)?,
                target,
            )
            .abs()
                < 1e-6
            {
                break;
            }
        }
    }
    if residual(
        options.coordinate,
        measure(options.coordinate, a, &x)?,
        target,
    )
    .abs()
        > 1e-5
    {
        return Err("Could not construct the requested scan geometry.".into());
    }
    Ok(x)
}
