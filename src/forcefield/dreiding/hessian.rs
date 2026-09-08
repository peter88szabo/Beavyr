//! Second derivatives of the DREIDING energy, for vibrational analysis.
//!
//! The Hessian is built by central differences of the **analytic** gradient, which is a different
//! proposition from differencing energies: 6N gradient evaluations rather than O(N²) energy ones,
//! each exact to machine precision rather than to the convergence of an SCF. For a force field
//! answering in microseconds, a full Hessian is effectively free -- a 60-atom molecule needs 360
//! gradient calls and finishes in well under a second.
//!
//! # What this cannot give you
//!
//! **Infrared intensities.** An intensity is the square of the dipole derivative along a normal
//! mode, and a dipole needs partial charges. DREIDING as implemented here has no charge model --
//! see `params::D_HB`, where the paper's own hydrogen-bond parameter is conditioned on that
//! choice -- so there is nothing to differentiate. Frequencies and normal modes are available;
//! the intensity column is not, and is left empty rather than filled with a fabricated number.
//!
//! Frequencies from a generic force field are also worth reading with care. DREIDING's own paper
//! reports it against equilibrium *structures*, not vibrational spectra, and its single bond force
//! constant of 700 kcal/mol/Å² is deliberately the same for every element pair. Treat these as a
//! quick check that a geometry really is a minimum, and as a starting point for mode assignment --
//! not as a substitute for a quantum-chemical Hessian.

use anyhow::{anyhow, Result};
use ndarray::Array2;

use crate::normalmode::fd_hessian::numerical_hessian_from_gradients;
use crate::optimizer::traits::Objective;

use super::objective::DreidingObjective;
use super::DreidingTopology;

/// Displacement for the central difference, in bohr.
///
/// Larger than a quantum-chemistry code would use, and deliberately so. There the step fights
/// SCF noise from below and anharmonicity from above; here the gradient is analytic and noiseless,
/// so the only error term is the O(h²) truncation, and this size keeps that far below the
/// accuracy the force field itself claims.
pub const STEP_BOHR: f64 = 5.0e-3;

/// The mass-independent Cartesian Hessian in Hartree/bohr², plus the gradient at the input
/// geometry in Hartree/bohr.
///
/// Coordinates in and out are bohr, matching [`crate::optimizer::traits::Objective`] and the
/// frequency code, not Beavyr's usual Å.
pub fn hessian_and_gradient(
    topology: &DreidingTopology,
    coords_bohr: &[f64],
) -> Result<(Array2<f64>, Vec<f64>)> {
    if coords_bohr.len() != topology.natoms() * 3 {
        return Err(anyhow!(
            "DREIDING Hessian expected {} coordinates for {} atoms, got {}",
            topology.natoms() * 3,
            topology.natoms(),
            coords_bohr.len()
        ));
    }

    let mut objective = DreidingObjective::new(topology);
    let mut gradient = vec![0.0; coords_bohr.len()];
    objective
        .gradient(coords_bohr, &mut gradient)
        .map_err(|e| anyhow!("DREIDING gradient failed: {e}"))?;

    let hessian = numerical_hessian_from_gradients(coords_bohr, STEP_BOHR, 0, |displaced| {
        let mut g = vec![0.0; displaced.len()];
        objective.gradient(displaced, &mut g)?;
        Ok(g)
    })?;

    Ok((hessian, gradient))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
    use crate::molecule::Molecule;

    fn molecule(xyz: &str) -> Molecule {
        let mut mol = Molecule::from_xyz(xyz);
        mol.recompute_bonds(1.2, 2.5);
        mol
    }

    fn coords_bohr(mol: &Molecule) -> Vec<f64> {
        mol.pos
            .iter()
            .flat_map(|p| {
                [
                    p.x as f64 / BOHR_TO_ANGSTROM,
                    p.y as f64 / BOHR_TO_ANGSTROM,
                    p.z as f64 / BOHR_TO_ANGSTROM,
                ]
            })
            .collect()
    }

    const WATER: &str = "3\n\nO 0.000 0.000 0.000\nH 0.980 0.000 0.000\nH -0.245 0.949 0.000\n";

    #[test]
    fn the_hessian_is_symmetric_and_the_right_shape() {
        let mol = molecule(WATER);
        let topology = DreidingTopology::build(&mol).unwrap();
        let (hessian, gradient) = hessian_and_gradient(&topology, &coords_bohr(&mol)).unwrap();

        assert_eq!(hessian.dim(), (9, 9));
        assert_eq!(gradient.len(), 9);
        for i in 0..9 {
            for j in 0..9 {
                assert!(
                    (hessian[(i, j)] - hessian[(j, i)]).abs() < 1.0e-10,
                    "asymmetric at ({i}, {j})"
                );
            }
        }
    }

    /// A Hessian must have three zero eigenvalues for translation and, for a non-linear molecule,
    /// three more for rotation. Their absence is the classic sign of a gradient or Hessian bug,
    /// so this is the strongest cheap check available.
    ///
    /// It has to be done **at a relaxed geometry**. Translational invariance holds everywhere, but
    /// the rotational eigenvalues only vanish at a stationary point: away from one they pick up a
    /// contribution proportional to the residual gradient, since rotating a molecule that is
    /// still being pulled on does change its energy.
    #[test]
    fn six_eigenvalues_vanish_for_translation_and_rotation() {
        use crate::forcefield::dreiding::objective::{cleanup_options, relax_angstrom};
        use crate::normalmode::linalg_shim::{Backend, LinAlg};

        let mol = molecule(WATER);
        let topology = DreidingTopology::build(&mol).unwrap();
        let start: Vec<f64> = mol
            .pos
            .iter()
            .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
            .collect();
        let (relaxed, report) = relax_angstrom(&topology, &start, cleanup_options()).unwrap();
        assert!(report.converged, "{}", report.message);
        let relaxed_bohr: Vec<f64> = relaxed.iter().map(|a| a / BOHR_TO_ANGSTROM).collect();

        let (hessian, _) = hessian_and_gradient(&topology, &relaxed_bohr).unwrap();
        let solver = LinAlg::new(Backend::Auto).unwrap();
        let (values, _) = solver.eigh(&hessian).unwrap();
        let mut magnitudes: Vec<f64> = values.iter().map(|v| v.abs()).collect();
        magnitudes.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Water: 9 degrees of freedom, 3 vibrations, so 6 must vanish. The threshold allows for
        // the central difference's O(h²) truncation and the residual force the minimiser stops
        // at; the gap to the first real vibration is four orders of magnitude, so there is no
        // ambiguity about which is which.
        for (index, value) in magnitudes.iter().take(6).enumerate() {
            assert!(
                *value < 1.0e-4,
                "eigenvalue {index} is {value}, expected ~0 for a rigid-body mode: {magnitudes:?}"
            );
        }
        // And the three genuine vibrations must be clearly non-zero.
        assert!(
            magnitudes[6] > 1.0e-2,
            "no vibrational eigenvalue survived: {magnitudes:?}"
        );
    }

    /// At a relaxed geometry the gradient must vanish, which is what makes the Hessian's
    /// eigenvalues meaningful as force constants.
    #[test]
    fn the_gradient_vanishes_at_a_relaxed_geometry() {
        use crate::forcefield::dreiding::objective::{cleanup_options, relax_angstrom};

        let mol = molecule(WATER);
        let topology = DreidingTopology::build(&mol).unwrap();
        let start: Vec<f64> = mol
            .pos
            .iter()
            .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
            .collect();
        let (relaxed, _) = relax_angstrom(&topology, &start, cleanup_options()).unwrap();

        let relaxed_bohr: Vec<f64> = relaxed.iter().map(|a| a / BOHR_TO_ANGSTROM).collect();
        let (_, gradient) = hessian_and_gradient(&topology, &relaxed_bohr).unwrap();
        let largest = gradient.iter().fold(0.0_f64, |acc, g| acc.max(g.abs()));
        assert!(largest < 1.0e-4, "largest gradient component is {largest}");
    }

    /// A relaxed minimum must have no imaginary frequency: every vibrational eigenvalue positive.
    #[test]
    fn a_relaxed_minimum_has_no_negative_curvature() {
        use crate::forcefield::dreiding::objective::{cleanup_options, relax_angstrom};
        use crate::normalmode::linalg_shim::{Backend, LinAlg};

        let mol = molecule(WATER);
        let topology = DreidingTopology::build(&mol).unwrap();
        let start: Vec<f64> = mol
            .pos
            .iter()
            .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
            .collect();
        let (relaxed, _) = relax_angstrom(&topology, &start, cleanup_options()).unwrap();
        let relaxed_bohr: Vec<f64> = relaxed.iter().map(|a| a / BOHR_TO_ANGSTROM).collect();

        let (hessian, _) = hessian_and_gradient(&topology, &relaxed_bohr).unwrap();
        let solver = LinAlg::new(Backend::Auto).unwrap();
        let (values, _) = solver.eigh(&hessian).unwrap();
        let mut sorted: Vec<f64> = values.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // The three smallest may sit slightly either side of zero; the vibrations must not.
        for value in sorted.iter().skip(6) {
            assert!(*value > 0.0, "negative curvature at a minimum: {sorted:?}");
        }
    }

    #[test]
    fn a_coordinate_count_mismatch_is_reported() {
        let mol = molecule(WATER);
        let topology = DreidingTopology::build(&mol).unwrap();
        assert!(hessian_and_gradient(&topology, &[0.0; 6]).is_err());
    }
}
