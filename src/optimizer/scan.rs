//! Rigid bond and angle scans, imported from Behemoth.
//!
//! The scan drivers are generic over an energy callback -- `F: FnMut(&Molecule) -> Result<f64>` --
//! so they carry across unchanged and work with any energy source, the DREIDING force field
//! included. Behemoth also ships `rigid_bond_scan_method` and `rigid_angle_scan_method`, thin
//! wrappers that call its own quantum chemistry through `Molecule::method(...)`; those are the
//! only quantum-chemistry-bound part and are not imported. Pass a callback instead.
//!
//! `Molecule` here is [`crate::optimizer::scan_molecule::ScanMolecule`], which presents the few
//! fields the scan code touches. Coordinates are in **bohr**.

//! Rigid coordinate scans.
//!
//! This module intentionally implements only the simplest scan needed for
//! diagnostics: set one bond distance or one bond angle to a sequence of values
//! and evaluate an energy callback at each generated geometry.  These are not
//! relaxed scans; no geometry optimization or constraint projection is performed.

use anyhow::{bail, Result};

use crate::optimizer::scan_molecule::ScanMolecule as Molecule;

// Matches `xyzparser::parse_xyz`'s ANGSTROM_TO_BOHR deliberately (see the
// note there) rather than the more precise CODATA 2018 value, so a scan step
// given in Angstrom lands on the same bohr geometry the rest of the pipeline
// would parse it to.
const ANGSTROM_TO_BOHR: f64 = 1.8897259886;
const BOHR_TO_ANGSTROM: f64 = 1.0 / ANGSTROM_TO_BOHR;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondScanMoveMode {
    /// Keep the first atom fixed and move only the second atom along the
    /// original bond axis.  This is the default rigid scan mode.
    MoveSecondAtom,
    /// Keep the second atom fixed and move only the first atom along the
    /// original bond axis.
    MoveFirstAtom,
    /// Move both atoms symmetrically around their original midpoint.
    SymmetricPair,
}

#[derive(Debug, Clone)]
pub struct BondScanOptions {
    /// Human-facing atom numbers, one-based, e.g. `[1, 2]`.
    pub atom_pair: [usize; 2],
    /// First scan distance in Angstrom.
    pub start_angstrom: f64,
    /// Last scan distance in Angstrom.
    pub end_angstrom: f64,
    /// Number of scan points including start and end.
    pub n_points: usize,
    /// Rigid coordinate update mode.
    pub move_mode: BondScanMoveMode,
}

impl BondScanOptions {
    pub fn new(
        atom_pair: [usize; 2],
        start_angstrom: f64,
        end_angstrom: f64,
        n_points: usize,
    ) -> Self {
        Self {
            atom_pair,
            start_angstrom,
            end_angstrom,
            n_points,
            move_mode: BondScanMoveMode::MoveSecondAtom,
        }
    }
}

#[derive(Debug, Clone)]
pub struct BondScanPoint {
    pub index: usize,
    pub distance_angstrom: f64,
    pub energy_hartree: f64,
    pub coords_bohr: Vec<[f64; 3]>,
}

#[derive(Debug, Clone)]
pub struct BondScanResult {
    pub atom_pair: [usize; 2],
    pub move_mode: BondScanMoveMode,
    pub points: Vec<BondScanPoint>,
}

impl BondScanResult {
    pub fn print_table(&self) {
        println!("{}", self.format_table());
    }

    pub fn format_table(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!(
            "Rigid bond scan for atoms {}-{} (no relaxation)\n",
            self.atom_pair[0], self.atom_pair[1]
        ));
        out.push_str(&format!("{:<6} {:>14} {:>18}\n", "point", "R/A", "E/Ha"));
        for p in &self.points {
            out.push_str(&format!(
                "{:<6} {:>14.6} {:>18.10}\n",
                p.index + 1,
                p.distance_angstrom,
                p.energy_hartree
            ));
        }
        out
    }
}

/// Run a rigid bond scan with a caller-supplied energy method.
///
/// Atom indices in `options.atom_pair` are one-based.  The callback receives a
/// complete `Molecule` with coordinates already updated to the current scan
/// point, and must return the energy in Hartree.
pub fn rigid_bond_scan<F>(
    molecule: &Molecule,
    options: BondScanOptions,
    mut energy_method: F,
) -> Result<BondScanResult>
where
    F: FnMut(&Molecule) -> Result<f64>,
{
    validate_options(molecule, &options)?;
    let atom_a = options.atom_pair[0] - 1;
    let atom_b = options.atom_pair[1] - 1;
    let axis = initial_bond_axis(molecule, atom_a, atom_b)?;
    let distances = scan_distances_angstrom(
        options.start_angstrom,
        options.end_angstrom,
        options.n_points,
    )?;

    let mut points = Vec::with_capacity(distances.len());
    for (index, distance_angstrom) in distances.into_iter().enumerate() {
        let mut scan_mol = molecule.clone();
        set_bond_distance(
            &mut scan_mol,
            atom_a,
            atom_b,
            &axis,
            distance_angstrom * ANGSTROM_TO_BOHR,
            options.move_mode,
        )?;
        let energy_hartree = energy_method(&scan_mol)?;
        points.push(BondScanPoint {
            index,
            distance_angstrom,
            energy_hartree,
            coords_bohr: scan_mol.xyz.clone(),
        });
    }

    Ok(BondScanResult {
        atom_pair: options.atom_pair,
        move_mode: options.move_mode,
        points,
    })
}

/// Run a rigid bond scan using the crate's normal `Molecule::method(...)`
/// dispatch.
///
/// This is the convenience entry point for ordinary method scans.  For
/// post-processing methods such as xTB+Frozen-APSG, use `rigid_bond_scan` and

pub fn bond_distance_angstrom(molecule: &Molecule, atom_pair: [usize; 2]) -> Result<f64> {
    if atom_pair[0] == 0 || atom_pair[1] == 0 {
        bail!("atom_pair uses one-based atom numbers; zero is invalid");
    }
    let a = atom_pair[0] - 1;
    let b = atom_pair[1] - 1;
    if a >= molecule.natoms || b >= molecule.natoms {
        bail!(
            "atom_pair {:?} is outside molecule with {} atoms",
            atom_pair,
            molecule.natoms
        );
    }
    Ok(distance_bohr(molecule.xyz[a], molecule.xyz[b]) * BOHR_TO_ANGSTROM)
}

#[derive(Debug, Clone)]
pub struct AngleScanOptions {
    /// Human-facing atom numbers, one-based.  `[A, B, C]` means angle A-B-C,
    /// with atom B as the vertex.
    pub atom_triple: [usize; 3],
    /// First scan angle in degrees.
    pub start_degrees: f64,
    /// Last scan angle in degrees.
    pub end_degrees: f64,
    /// Number of scan points including start and end.
    pub n_points: usize,
}

impl AngleScanOptions {
    pub fn new(
        atom_triple: [usize; 3],
        start_degrees: f64,
        end_degrees: f64,
        n_points: usize,
    ) -> Self {
        Self {
            atom_triple,
            start_degrees,
            end_degrees,
            n_points,
        }
    }
}

#[derive(Debug, Clone)]
pub struct AngleScanPoint {
    pub index: usize,
    pub angle_degrees: f64,
    pub energy_hartree: f64,
    pub coords_bohr: Vec<[f64; 3]>,
}

#[derive(Debug, Clone)]
pub struct AngleScanResult {
    pub atom_triple: [usize; 3],
    pub points: Vec<AngleScanPoint>,
}

impl AngleScanResult {
    pub fn print_table(&self) {
        println!("{}", self.format_table());
    }

    pub fn format_table(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!(
            "Rigid angle scan for atoms {}-{}-{} (no relaxation)\n",
            self.atom_triple[0], self.atom_triple[1], self.atom_triple[2]
        ));
        out.push_str(&format!(
            "{:<6} {:>14} {:>18}\n",
            "point", "angle/deg", "E/Ha"
        ));
        for p in &self.points {
            out.push_str(&format!(
                "{:<6} {:>14.6} {:>18.10}\n",
                p.index + 1,
                p.angle_degrees,
                p.energy_hartree
            ));
        }
        out
    }
}

/// Run a rigid angle scan with a caller-supplied energy method.
///
/// Atom indices in `options.atom_triple` are one-based.  `[A, B, C]` defines
/// angle A-B-C, with B as the vertex.  The scan keeps A and B fixed, preserves
/// the B-C distance from the input geometry, and rotates atom C in the initial
/// A-B-C plane.  No relaxation is performed.
pub fn rigid_angle_scan<F>(
    molecule: &Molecule,
    options: AngleScanOptions,
    mut energy_method: F,
) -> Result<AngleScanResult>
where
    F: FnMut(&Molecule) -> Result<f64>,
{
    validate_angle_options(molecule, &options)?;
    let atom_a = options.atom_triple[0] - 1;
    let atom_b = options.atom_triple[1] - 1;
    let atom_c = options.atom_triple[2] - 1;
    let frame = initial_angle_frame(molecule, atom_a, atom_b, atom_c)?;
    let angles = scan_values(options.start_degrees, options.end_degrees, options.n_points)?;

    let mut points = Vec::with_capacity(angles.len());
    for (index, angle_degrees) in angles.into_iter().enumerate() {
        let mut scan_mol = molecule.clone();
        set_bond_angle(
            &mut scan_mol,
            atom_b,
            atom_c,
            &frame,
            angle_degrees.to_radians(),
        )?;
        let energy_hartree = energy_method(&scan_mol)?;
        points.push(AngleScanPoint {
            index,
            angle_degrees,
            energy_hartree,
            coords_bohr: scan_mol.xyz.clone(),
        });
    }

    Ok(AngleScanResult {
        atom_triple: options.atom_triple,
        points,
    })
}

/// Run a rigid angle scan using the crate's normal `Molecule::method(...)`

pub fn bond_angle_degrees(molecule: &Molecule, atom_triple: [usize; 3]) -> Result<f64> {
    validate_atom_triple(molecule, atom_triple)?;
    let a = atom_triple[0] - 1;
    let b = atom_triple[1] - 1;
    let c = atom_triple[2] - 1;
    let ba = sub(molecule.xyz[a], molecule.xyz[b]);
    let bc = sub(molecule.xyz[c], molecule.xyz[b]);
    let nba = norm(ba);
    let nbc = norm(bc);
    if nba <= 1.0e-14 || nbc <= 1.0e-14 {
        bail!("cannot compute angle with coincident vertex/outer atom");
    }
    let cosang = (dot(ba, bc) / (nba * nbc)).clamp(-1.0, 1.0);
    Ok(cosang.acos().to_degrees())
}

#[derive(Debug, Clone, Copy)]
struct AngleFrame {
    vertex_to_a: [f64; 3],
    in_plane_perp: [f64; 3],
    vertex_to_c_distance: f64,
}

fn validate_angle_options(molecule: &Molecule, options: &AngleScanOptions) -> Result<()> {
    validate_atom_triple(molecule, options.atom_triple)?;
    if options.n_points == 0 {
        bail!("angle scan requires at least one point");
    }
    validate_angle_degrees(options.start_degrees)?;
    validate_angle_degrees(options.end_degrees)?;
    Ok(())
}

fn validate_atom_triple(molecule: &Molecule, atom_triple: [usize; 3]) -> Result<()> {
    if molecule.natoms == 0 {
        bail!("cannot scan an empty molecule");
    }
    if atom_triple.iter().any(|&i| i == 0) {
        bail!("atom_triple uses one-based atom numbers; zero is invalid");
    }
    if atom_triple[0] == atom_triple[1]
        || atom_triple[1] == atom_triple[2]
        || atom_triple[0] == atom_triple[2]
    {
        bail!("angle scan atom triple must contain three different atoms");
    }
    if atom_triple.iter().any(|&i| i > molecule.natoms) {
        bail!(
            "atom_triple {:?} is outside molecule with {} atoms",
            atom_triple,
            molecule.natoms
        );
    }
    Ok(())
}

fn validate_angle_degrees(angle: f64) -> Result<()> {
    if !(0.0..180.0).contains(&angle) {
        bail!("angle scan values must be greater than 0 and less than 180 degrees");
    }
    Ok(())
}

fn initial_angle_frame(
    molecule: &Molecule,
    atom_a: usize,
    atom_b: usize,
    atom_c: usize,
) -> Result<AngleFrame> {
    let vertex = molecule.xyz[atom_b];
    let ba = sub(molecule.xyz[atom_a], vertex);
    let bc = sub(molecule.xyz[atom_c], vertex);
    let ba_len = norm(ba);
    let bc_len = norm(bc);
    if ba_len <= 1.0e-14 || bc_len <= 1.0e-14 {
        bail!("cannot define angle scan frame with coincident vertex/outer atom");
    }
    let e1 = scale(ba, 1.0 / ba_len);
    let bc_unit = scale(bc, 1.0 / bc_len);
    let projected = sub(bc_unit, scale(e1, dot(bc_unit, e1)));
    let perp_len = norm(projected);
    let in_plane_perp = if perp_len > 1.0e-12 {
        scale(projected, 1.0 / perp_len)
    } else {
        arbitrary_perpendicular(e1)
    };
    Ok(AngleFrame {
        vertex_to_a: e1,
        in_plane_perp,
        vertex_to_c_distance: bc_len,
    })
}

fn set_bond_angle(
    molecule: &mut Molecule,
    atom_b: usize,
    atom_c: usize,
    frame: &AngleFrame,
    angle_radians: f64,
) -> Result<()> {
    let vertex = molecule.xyz[atom_b];
    let direction = add(
        scale(frame.vertex_to_a, angle_radians.cos()),
        scale(frame.in_plane_perp, angle_radians.sin()),
    );
    let rc = add(vertex, scale(direction, frame.vertex_to_c_distance));
    set_atom_coord(molecule, atom_c, rc)
}

fn scan_values(start: f64, end: f64, n_points: usize) -> Result<Vec<f64>> {
    if n_points == 0 {
        bail!("scan requires at least one point");
    }
    if n_points == 1 {
        return Ok(vec![start]);
    }
    let step = (end - start) / (n_points as f64 - 1.0);
    Ok((0..n_points).map(|i| start + step * i as f64).collect())
}

fn validate_options(molecule: &Molecule, options: &BondScanOptions) -> Result<()> {
    if molecule.natoms == 0 {
        bail!("cannot scan an empty molecule");
    }
    if options.atom_pair[0] == 0 || options.atom_pair[1] == 0 {
        bail!("atom_pair uses one-based atom numbers; zero is invalid");
    }
    if options.atom_pair[0] == options.atom_pair[1] {
        bail!("bond scan atom pair must contain two different atoms");
    }
    if options.atom_pair[0] > molecule.natoms || options.atom_pair[1] > molecule.natoms {
        bail!(
            "atom_pair {:?} is outside molecule with {} atoms",
            options.atom_pair,
            molecule.natoms
        );
    }
    if options.n_points == 0 {
        bail!("bond scan requires at least one point");
    }
    if options.start_angstrom <= 0.0 || options.end_angstrom <= 0.0 {
        bail!("bond scan distances must be positive Angstrom values");
    }
    Ok(())
}

fn scan_distances_angstrom(start: f64, end: f64, n_points: usize) -> Result<Vec<f64>> {
    scan_values(start, end, n_points)
}

fn initial_bond_axis(molecule: &Molecule, atom_a: usize, atom_b: usize) -> Result<[f64; 3]> {
    let ra = molecule.xyz[atom_a];
    let rb = molecule.xyz[atom_b];
    let dx = [rb[0] - ra[0], rb[1] - ra[1], rb[2] - ra[2]];
    let norm = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
    if norm <= 1.0e-14 {
        bail!(
            "cannot define bond scan axis: atoms {} and {} are coincident",
            atom_a + 1,
            atom_b + 1
        );
    }
    Ok([dx[0] / norm, dx[1] / norm, dx[2] / norm])
}

fn set_bond_distance(
    molecule: &mut Molecule,
    atom_a: usize,
    atom_b: usize,
    axis_ab: &[f64; 3],
    distance_bohr: f64,
    move_mode: BondScanMoveMode,
) -> Result<()> {
    match move_mode {
        BondScanMoveMode::MoveSecondAtom => {
            let ra = molecule.xyz[atom_a];
            let rb = [
                ra[0] + distance_bohr * axis_ab[0],
                ra[1] + distance_bohr * axis_ab[1],
                ra[2] + distance_bohr * axis_ab[2],
            ];
            set_atom_coord(molecule, atom_b, rb)?;
        }
        BondScanMoveMode::MoveFirstAtom => {
            let rb = molecule.xyz[atom_b];
            let ra = [
                rb[0] - distance_bohr * axis_ab[0],
                rb[1] - distance_bohr * axis_ab[1],
                rb[2] - distance_bohr * axis_ab[2],
            ];
            set_atom_coord(molecule, atom_a, ra)?;
        }
        BondScanMoveMode::SymmetricPair => {
            let ra0 = molecule.xyz[atom_a];
            let rb0 = molecule.xyz[atom_b];
            let mid = [
                0.5 * (ra0[0] + rb0[0]),
                0.5 * (ra0[1] + rb0[1]),
                0.5 * (ra0[2] + rb0[2]),
            ];
            let half = 0.5 * distance_bohr;
            let ra = [
                mid[0] - half * axis_ab[0],
                mid[1] - half * axis_ab[1],
                mid[2] - half * axis_ab[2],
            ];
            let rb = [
                mid[0] + half * axis_ab[0],
                mid[1] + half * axis_ab[1],
                mid[2] + half * axis_ab[2],
            ];
            set_atom_coord(molecule, atom_a, ra)?;
            set_atom_coord(molecule, atom_b, rb)?;
        }
    }
    Ok(())
}

fn set_atom_coord(molecule: &mut Molecule, atom: usize, coord: [f64; 3]) -> Result<()> {
    if atom >= molecule.natoms {
        bail!(
            "atom index {} outside molecule with {} atoms",
            atom + 1,
            molecule.natoms
        );
    }
    molecule.xyz[atom] = coord;
    molecule.atom[atom].coord = coord;
    Ok(())
}

fn distance_bohr(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn arbitrary_perpendicular(a: [f64; 3]) -> [f64; 3] {
    let trial = if a[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let projected = sub(trial, scale(a, dot(trial, a)));
    let n = norm(projected);
    if n <= 1.0e-14 {
        [0.0, 0.0, 1.0]
    } else {
        scale(projected, 1.0 / n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rigid_bond_scan_sets_requested_distances() -> Result<()> {
        let mol = Molecule::from_xyz("2\n\nH 0.0 0.0 0.0\nH 0.7 0.0 0.0\n");
        let options = BondScanOptions::new([1, 2], 0.8, 1.2, 3);
        let result = rigid_bond_scan(&mol, options, |m| bond_distance_angstrom(m, [1, 2]))?;

        assert_eq!(result.points.len(), 3);
        assert!((result.points[0].distance_angstrom - 0.8).abs() < 1.0e-12);
        assert!((result.points[1].distance_angstrom - 1.0).abs() < 1.0e-12);
        assert!((result.points[2].distance_angstrom - 1.2).abs() < 1.0e-12);
        assert!((result.points[0].energy_hartree - 0.8).abs() < 1.0e-10);
        assert!((result.points[1].energy_hartree - 1.0).abs() < 1.0e-10);
        assert!((result.points[2].energy_hartree - 1.2).abs() < 1.0e-10);
        Ok(())
    }

    #[test]
    fn rigid_bond_scan_rejects_zero_based_atoms() {
        let mol = Molecule::from_xyz("2\n\nH 0.0 0.0 0.0\nH 0.7 0.0 0.0\n");
        let options = BondScanOptions::new([0, 1], 0.8, 1.2, 3);
        assert!(rigid_bond_scan(&mol, options, |_| Ok(0.0)).is_err());
    }

    #[test]
    fn rigid_angle_scan_sets_requested_angles() -> Result<()> {
        let mol = Molecule::from_xyz("3\n\nH 1.0 0.0 0.0\nO 0.0 0.0 0.0\nH 0.0 1.0 0.0\n");
        let options = AngleScanOptions::new([1, 2, 3], 80.0, 120.0, 3);
        let result = rigid_angle_scan(&mol, options, |m| bond_angle_degrees(m, [1, 2, 3]))?;

        assert_eq!(result.points.len(), 3);
        assert!((result.points[0].energy_hartree - 80.0).abs() < 1.0e-10);
        assert!((result.points[1].energy_hartree - 100.0).abs() < 1.0e-10);
        assert!((result.points[2].energy_hartree - 120.0).abs() < 1.0e-10);
        Ok(())
    }

    #[test]
    fn rigid_angle_scan_rejects_duplicate_atoms() {
        let mol = Molecule::from_xyz("3\n\nH 1.0 0.0 0.0\nO 0.0 0.0 0.0\nH 0.0 1.0 0.0\n");
        let options = AngleScanOptions::new([1, 2, 1], 80.0, 120.0, 3);
        assert!(rigid_angle_scan(&mol, options, |_| Ok(0.0)).is_err());
    }
}
