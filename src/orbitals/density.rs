//! Total electron density and spin density.
//!
//! Both are quadratic forms in the atomic-orbital basis:
//!
//! ```text
//!   rho(r) = sum_{mu,nu} P_{mu,nu} phi_mu(r) phi_nu(r)
//! ```
//!
//! with `P` the density matrix.  Going through `P` rather than summing
//! `n_i |psi_i|^2` over occupied orbitals matters at scale: the orbital sum
//! costs one dot product per occupied orbital at every grid point (185 of them
//! for a medium open-shell system), whereas the quadratic form is a single pass
//! over the basis functions that are actually significant there.  It also
//! handles fractional occupations without a special case.

use bevy::prelude::*;
use bevy::tasks::{ComputeTaskPool, TaskPool};

use super::evaluate::{shell_cutoff_radius, shell_values, EvalScratch, MAX_SHELL_COMPONENTS};
use super::grid::{GridSpec, ScalarField, BOHR};
use super::molden::OrbitalData;

/// Occupations below this are treated as empty.
const OCCUPATION_FLOOR: f64 = 1.0e-8;
/// Basis functions smaller than this are dropped from the quadratic form.
const VALUE_FLOOR: f64 = 1.0e-8;

/// Which scalar field to build.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum Quantity {
    /// A single molecular orbital, psi_i.
    #[default]
    Orbital,
    /// Total electron density, rho_alpha + rho_beta.
    Density,
    /// Spin density, rho_alpha - rho_beta.
    SpinDensity,
}

impl Quantity {
    /// Whether the field takes both signs, and so needs two lobes drawn.
    ///
    /// Total density is non-negative everywhere; an orbital and a spin density
    /// are not.
    pub fn is_signed(self) -> bool {
        !matches!(self, Quantity::Density)
    }

    /// Sensible isovalue slider bounds.
    ///
    /// The ranges differ by orders of magnitude, so a single range would be
    /// useless in at least one mode.
    pub fn isovalue_range(self) -> (f32, f32) {
        match self {
            Quantity::Orbital => (0.001, 0.5),
            Quantity::Density => (0.001, 1.0),
            Quantity::SpinDensity => (0.0005, 0.5),
        }
    }

    /// Isovalue to adopt when switching into this mode.
    pub fn default_isovalue(self) -> f32 {
        match self {
            Quantity::Orbital => 0.10,
            Quantity::Density => 0.05,
            Quantity::SpinDensity => 0.01,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Quantity::Orbital => "Orbitals",
            Quantity::Density => "Electron density",
            Quantity::SpinDensity => "Spin density",
        }
    }

    /// Spin density is only meaningful when the two spin sets differ.
    pub fn is_available(self, data: &OrbitalData) -> bool {
        !matches!(self, Quantity::SpinDensity) || data.is_open_shell()
    }
}

/// Density matrices in the AO basis, row-major `n_basis x n_basis`.
#[derive(Debug)]
pub struct DensityMatrices {
    pub n_basis: usize,
    /// `P_alpha + P_beta`.
    pub total: Vec<f64>,
    /// `P_alpha - P_beta`; `None` for a closed-shell file.
    pub spin: Option<Vec<f64>>,
}

impl DensityMatrices {
    /// Accumulate `sum_i occ_i c_i c_i^T` over one spin set.
    fn accumulate(matrix: &mut [f64], n: usize, orbitals: &[super::molden::MolecularOrbital]) {
        for orbital in orbitals {
            if orbital.occupation <= OCCUPATION_FLOOR {
                continue;
            }
            let c = &orbital.coefficients;
            // Only the basis functions this orbital actually uses; a localised
            // orbital leaves most of a large basis untouched.
            let used: Vec<usize> = (0..n).filter(|&i| c[i].abs() > VALUE_FLOOR).collect();
            for &i in &used {
                let weighted = orbital.occupation * c[i];
                let row = i * n;
                for &j in &used {
                    matrix[row + j] += weighted * c[j];
                }
            }
        }
    }

    pub fn build(data: &OrbitalData) -> Self {
        let n = data.layout.n_basis;

        if data.is_open_shell() {
            let mut alpha = vec![0.0; n * n];
            let mut beta = vec![0.0; n * n];
            Self::accumulate(&mut alpha, n, &data.alpha);
            Self::accumulate(&mut beta, n, &data.beta);

            let total = alpha
                .iter()
                .zip(&beta)
                .map(|(a, b)| a + b)
                .collect::<Vec<f64>>();
            let spin = alpha
                .iter()
                .zip(&beta)
                .map(|(a, b)| a - b)
                .collect::<Vec<f64>>();

            Self {
                n_basis: n,
                total,
                spin: Some(spin),
            }
        } else {
            // A closed-shell file carries occupation 2.0 on its single set, so
            // accumulating it directly already gives the total density.
            let mut total = vec![0.0; n * n];
            Self::accumulate(&mut total, n, &data.alpha);
            Self {
                n_basis: n,
                total,
                spin: None,
            }
        }
    }

    fn matrix_for(&self, quantity: Quantity) -> Option<&[f64]> {
        match quantity {
            Quantity::Density => Some(&self.total),
            Quantity::SpinDensity => self.spin.as_deref(),
            Quantity::Orbital => None,
        }
    }
}

/// Sample a density-like quantity onto `spec`.
///
/// Returns an all-zero field if the requested quantity is unavailable, which
/// is only reachable for a spin density on a closed-shell file.
pub fn sample_density(
    data: &OrbitalData,
    matrices: &DensityMatrices,
    quantity: Quantity,
    spec: &GridSpec,
) -> ScalarField {
    let mut values = vec![0.0f32; spec.point_count()];
    let Some(matrix) = matrices.matrix_for(quantity) else {
        return ScalarField {
            spec: spec.clone(),
            values,
        };
    };

    let n = matrices.n_basis;
    let centers: Vec<[f64; 3]> = data
        .positions
        .iter()
        .map(|p| [p.x as f64 / BOHR, p.y as f64 / BOHR, p.z as f64 / BOHR])
        .collect();

    let [nx, ny, _nz] = spec.dims;
    let step = spec.step as f64 / BOHR;
    let origin = [
        spec.origin.x as f64 / BOHR,
        spec.origin.y as f64 / BOHR,
        spec.origin.z as f64 / BOHR,
    ];
    let slab = nx * ny;

    // `get_or_init` so this is callable from a plain unit test.
    ComputeTaskPool::get_or_init(TaskPool::default).scope(|scope| {
        for (iz, plane) in values.chunks_mut(slab).enumerate() {
            let centers = &centers;
            scope.spawn(async move {
                let mut scratch = EvalScratch::default();
                let mut components = [0.0f64; MAX_SHELL_COMPONENTS];
                // Significant basis functions at the current point.
                let mut indices: Vec<usize> = Vec::with_capacity(256);
                let mut amplitudes: Vec<f64> = Vec::with_capacity(256);

                let z = origin[2] + iz as f64 * step;
                for iy in 0..ny {
                    let y = origin[1] + iy as f64 * step;
                    for ix in 0..nx {
                        let x = origin[0] + ix as f64 * step;

                        indices.clear();
                        amplitudes.clear();
                        for (index, shell) in data.shells.iter().enumerate() {
                            let c = centers[shell.center];
                            let d = [x - c[0], y - c[1], z - c[2]];
                            let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                            let cutoff = shell_cutoff_radius(shell);
                            if r2 > cutoff * cutoff {
                                continue;
                            }
                            let count = shell.n_functions();
                            if !shell_values(shell, d, &mut scratch, &mut components[..count]) {
                                continue;
                            }
                            let start = data.layout.shell_offsets[index];
                            for (slot, &value) in components[..count].iter().enumerate() {
                                if value.abs() > VALUE_FLOOR {
                                    indices.push(start + slot);
                                    amplitudes.push(value);
                                }
                            }
                        }

                        // Symmetric quadratic form: diagonal once, off-diagonal
                        // twice, so only the upper triangle is visited.
                        let mut total = 0.0;
                        for a in 0..indices.len() {
                            let row = indices[a] * n;
                            let va = amplitudes[a];
                            total += matrix[row + indices[a]] * va * va;
                            let mut cross = 0.0;
                            for b in (a + 1)..indices.len() {
                                cross += matrix[row + indices[b]] * amplitudes[b];
                            }
                            total += 2.0 * va * cross;
                        }

                        plane[iy * nx + ix] = total as f32;
                    }
                }
            });
        }
    });

    ScalarField {
        spec: spec.clone(),
        values,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::orbitals::evaluate::{accumulate_shell, EvalScratch};
    use crate::orbitals::grid::plan_grid;
    use crate::orbitals::molden::parse_molden;

    /// H2 with one s function per atom, closed shell, two electrons.
    ///
    /// The MO coefficient is the analytic normalisation for this geometry and
    /// exponent, not an invented number:
    ///
    /// ```text
    ///   S12 = exp(-0.6 R^2) = 0.309341616  at R = 0.74 A = 1.398577 a0
    ///   c   = 1 / sqrt(2 (1 + S12))        = 0.617957370
    /// ```
    ///
    /// so the orbital is exactly normalised and the density integrates to 2.
    /// Getting this wrong is easy and looks like a bug in the density code: a
    /// coefficient of 0.5483 gives an orbital of norm 0.787 and an integral of
    /// 1.5745, which is correct arithmetic on an unnormalised wavefunction.
    const H2: &str = "\
[Molden Format]
[Atoms] Angs
H     1    1    0.000000    0.000000   -0.370000
H     2    1    0.000000    0.000000    0.370000
[GTO]
  1 0
s   1 1.0
      1.2000000000    1.0000000000

  2 0
s   1 1.0
      1.2000000000    1.0000000000

[MO]
 Sym= 1a
 Ene= -0.5950000000
 Spin= Alpha
 Occup= 2.000000
   1    0.6179573700
   2    0.6179573700
";

    /// Evaluate one molecular orbital directly, for cross-checking the
    /// density-matrix contraction against the already-validated orbital path.
    fn orbital_value(data: &OrbitalData, orbital_index: usize, spin: bool, p: [f64; 3]) -> f64 {
        let orbitals = if spin { &data.beta } else { &data.alpha };
        let coefficients = &orbitals[orbital_index].coefficients;
        let mut scratch = EvalScratch::default();
        let mut total = 0.0;
        for (index, shell) in data.shells.iter().enumerate() {
            let c = data.positions[shell.center];
            let center = [c.x as f64 / BOHR, c.y as f64 / BOHR, c.z as f64 / BOHR];
            let d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
            let start = data.layout.shell_offsets[index];
            let weights = &coefficients[start..start + shell.n_functions()];
            accumulate_shell(shell, d, weights, &mut scratch, &mut total);
        }
        total
    }

    /// Density from the matrix, at one point, without a grid.
    fn density_at(
        data: &OrbitalData,
        matrices: &DensityMatrices,
        quantity: Quantity,
        p: [f64; 3],
    ) -> f64 {
        let n = matrices.n_basis;
        let matrix = matrices
            .matrix_for(quantity)
            .expect("quantity has a matrix");
        let mut scratch = EvalScratch::default();
        let mut phi = vec![0.0f64; n];
        for (index, shell) in data.shells.iter().enumerate() {
            let c = data.positions[shell.center];
            let center = [c.x as f64 / BOHR, c.y as f64 / BOHR, c.z as f64 / BOHR];
            let d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
            let start = data.layout.shell_offsets[index];
            let count = shell.n_functions();
            for slot in 0..count {
                let mut weights = vec![0.0; count];
                weights[slot] = 1.0;
                let mut value = 0.0;
                accumulate_shell(shell, d, &weights, &mut scratch, &mut value);
                phi[start + slot] = value;
            }
        }
        let mut total = 0.0;
        for i in 0..n {
            for j in 0..n {
                total += matrix[i * n + j] * phi[i] * phi[j];
            }
        }
        total
    }

    /// The absolute check: integrating the density must recover the electron
    /// count.  A two-electron system with a single smooth s primitive per atom
    /// has no nuclear cusp sharp enough to defeat a uniform grid, so this can
    /// be done honestly.
    #[test]
    fn integrating_the_density_recovers_the_electron_count() {
        let data = parse_molden(H2).expect("fixture parses");
        let matrices = DensityMatrices::build(&data);
        assert!(matrices.spin.is_none(), "closed shell has no spin density");

        let spec = plan_grid(&data.positions, 5.0, 0.08, 20_000_000);
        let field = sample_density(&data, &matrices, Quantity::Density, &spec);

        let cell = (spec.step as f64).powi(3) / BOHR.powi(3);
        let electrons: f64 = field.values.iter().map(|&v| v as f64).sum::<f64>() * cell;

        assert!(
            (electrons - 2.0).abs() < 0.01,
            "integrated density was {electrons} electrons, expected 2"
        );
    }

    /// The density matrix must reproduce the occupation-weighted sum of squared
    /// orbitals pointwise.  This checks the matrix construction and the
    /// quadratic form against the orbital path that is already validated, with
    /// no grid-integration error in the way.
    #[test]
    fn the_density_matrix_matches_the_orbital_sum() {
        let data = parse_molden(H2).expect("fixture parses");
        let matrices = DensityMatrices::build(&data);

        for point in [
            [0.0, 0.0, 0.0],
            [0.3, -0.2, 0.7],
            [1.4, 0.9, -1.1],
            [-2.0, 0.4, 0.3],
        ] {
            let from_matrix = density_at(&data, &matrices, Quantity::Density, point);

            let mut from_orbitals = 0.0;
            for (index, orbital) in data.alpha.iter().enumerate() {
                if orbital.occupation <= OCCUPATION_FLOOR {
                    continue;
                }
                let psi = orbital_value(&data, index, false, point);
                from_orbitals += orbital.occupation * psi * psi;
            }

            assert!(
                (from_matrix - from_orbitals).abs() < 1.0e-10,
                "at {point:?}: matrix gave {from_matrix}, orbital sum gave {from_orbitals}"
            );
        }
    }

    #[test]
    fn total_density_is_never_negative() {
        let data = parse_molden(H2).expect("fixture parses");
        let matrices = DensityMatrices::build(&data);
        let spec = plan_grid(&data.positions, 3.0, 0.15, 5_000_000);
        let field = sample_density(&data, &matrices, Quantity::Density, &spec);

        let minimum = field.values.iter().cloned().fold(f32::MAX, f32::min);
        assert!(minimum >= -1.0e-6, "density went negative: {minimum}");
    }

    #[test]
    fn spin_density_is_unavailable_for_a_closed_shell_file() {
        let data = parse_molden(H2).expect("fixture parses");
        assert!(!Quantity::SpinDensity.is_available(&data));
        assert!(Quantity::Density.is_available(&data));
        assert!(Quantity::Orbital.is_available(&data));
    }

    #[test]
    fn only_the_total_density_is_unsigned() {
        assert!(!Quantity::Density.is_signed());
        assert!(Quantity::SpinDensity.is_signed());
        assert!(Quantity::Orbital.is_signed());
    }
}
