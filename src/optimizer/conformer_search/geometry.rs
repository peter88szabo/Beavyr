use anyhow::{bail, ensure, Result};

use super::topology::{PreparedTorsion, TorsionKind};

const ANGSTROM_TO_BOHR: f64 = 1.0 / 0.529_177_210_903;
const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

fn atom(x: &[f64], index: usize) -> [f64; 3] {
    let offset = 3 * index;
    [x[offset], x[offset + 1], x[offset + 2]]
}

fn set_atom(x: &mut [f64], index: usize, value: [f64; 3]) {
    let offset = 3 * index;
    x[offset..offset + 3].copy_from_slice(&value);
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}

pub(crate) fn wrap_angle(mut angle: f64) -> f64 {
    while angle <= -std::f64::consts::PI {
        angle += TWO_PI;
    }
    while angle > std::f64::consts::PI {
        angle -= TWO_PI;
    }
    angle
}

pub(crate) fn dihedral_angle(x: &[f64], indices: [usize; 4]) -> f64 {
    let [a, b, c, d] = indices.map(|index| atom(x, index));
    let ab = sub(b, a);
    let bc = sub(c, b);
    let cd = sub(d, c);
    let left = cross(ab, bc);
    let right = cross(bc, cd);
    let denominator = norm(left) * norm(right);
    if denominator <= 1.0e-15 || norm(bc) <= 1.0e-15 {
        return 0.0;
    }
    let cosine = (dot(left, right) / denominator).clamp(-1.0, 1.0);
    let sine = dot(bc, cross(left, right)) / (norm(bc) * denominator);
    sine.atan2(cosine)
}

fn rotate_atoms_about_bond(
    x: &mut [f64],
    bond: [usize; 2],
    moving_atoms: &[usize],
    angle: f64,
) -> Result<()> {
    let origin = atom(x, bond[0]);
    let axis_vector = sub(atom(x, bond[1]), origin);
    let axis_norm = norm(axis_vector);
    ensure!(
        axis_norm > 1.0e-12,
        "cannot rotate around a zero-length bond"
    );
    let axis = axis_vector.map(|value| value / axis_norm);
    let cosine = angle.cos();
    let sine = angle.sin();
    for &index in moving_atoms {
        let relative = sub(atom(x, index), origin);
        let axial = dot(axis, relative);
        let perpendicular = cross(axis, relative);
        let rotated = [
            relative[0] * cosine + perpendicular[0] * sine + axis[0] * axial * (1.0 - cosine),
            relative[1] * cosine + perpendicular[1] * sine + axis[1] * axial * (1.0 - cosine),
            relative[2] * cosine + perpendicular[2] * sine + axis[2] * axial * (1.0 - cosine),
        ];
        set_atom(
            x,
            index,
            [
                origin[0] + rotated[0],
                origin[1] + rotated[1],
                origin[2] + rotated[2],
            ],
        );
    }
    Ok(())
}

fn set_dihedral(x: &mut [f64], torsion: &PreparedTorsion, target: f64) -> Result<()> {
    let current = dihedral_angle(x, torsion.atoms);
    let desired_change = wrap_angle(target - current);
    if desired_change.abs() < 1.0e-12 {
        return Ok(());
    }

    // Determine the sign implied by the atom ordering instead of assuming one
    // convention for every torsion representation.
    let probe = 1.0e-5;
    let mut tested = x.to_vec();
    rotate_atoms_about_bond(
        &mut tested,
        torsion.central_bond,
        &torsion.moving_atoms,
        probe,
    )?;
    let response = wrap_angle(dihedral_angle(&tested, torsion.atoms) - current);
    if response.abs() < 1.0e-10 {
        bail!("torsion {:?} is geometrically singular", torsion.atoms);
    }
    let signed_rotation = desired_change * probe / response;
    rotate_atoms_about_bond(
        x,
        torsion.central_bond,
        &torsion.moving_atoms,
        signed_rotation,
    )
}

pub(crate) fn coordinates_from_genome(
    template: &[f64],
    torsions: &[PreparedTorsion],
    genome: &[f64],
) -> Result<Vec<f64>> {
    ensure!(
        torsions.len() == genome.len(),
        "conformer genome has {} values for {} torsions",
        genome.len(),
        torsions.len()
    );
    let mut coordinates = template.to_vec();
    for (torsion, &target) in torsions.iter().zip(genome) {
        set_dihedral(&mut coordinates, torsion, target)?;
    }
    Ok(coordinates)
}

pub(crate) fn genome_from_coordinates(
    coordinates: &[f64],
    torsions: &[PreparedTorsion],
) -> Vec<f64> {
    torsions
        .iter()
        .map(|torsion| match torsion.kind {
            TorsionKind::Rotatable => dihedral_angle(coordinates, torsion.atoms),
            TorsionKind::CisTrans => {
                let angle = dihedral_angle(coordinates, torsion.atoms);
                if angle.abs() <= std::f64::consts::FRAC_PI_2 {
                    0.0
                } else {
                    std::f64::consts::PI
                }
            }
        })
        .collect()
}

fn distance(x: &[f64], a: usize, b: usize) -> f64 {
    norm(sub(atom(x, a), atom(x, b)))
}

pub(crate) fn sensible_geometry(
    coordinates: &[f64],
    natoms: usize,
    bonds: &[[usize; 2]],
    minimum_nonbonded_angstrom: f64,
    maximum_bonded_angstrom: f64,
) -> bool {
    if coordinates.len() != 3 * natoms || coordinates.iter().any(|value| !value.is_finite()) {
        return false;
    }
    let minimum_nonbonded = minimum_nonbonded_angstrom * ANGSTROM_TO_BOHR;
    let maximum_bonded = maximum_bonded_angstrom * ANGSTROM_TO_BOHR;
    let mut bonded = vec![false; natoms * natoms];
    for &[a, b] in bonds {
        if distance(coordinates, a, b) > maximum_bonded {
            return false;
        }
        bonded[a * natoms + b] = true;
        bonded[b * natoms + a] = true;
    }
    for a in 0..natoms {
        for b in (a + 1)..natoms {
            if !bonded[a * natoms + b] && distance(coordinates, a, b) < minimum_nonbonded {
                return false;
            }
        }
    }
    true
}
