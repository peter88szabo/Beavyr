//! Sampling a molecular orbital onto a uniform scalar field.
//!
//! The loop is over shells rather than grid points.  Each shell touches only
//! the points inside its own cutoff sphere, so cost scales with the basis's
//! spatial extent rather than with `n_shells * n_points`.  Combined with
//! skipping shells whose MO coefficients are all negligible -- most of them,
//! for a localised orbital -- this is what makes a 546-function system usable.

use bevy::prelude::*;
use bevy::tasks::{ComputeTaskPool, TaskPool};

use super::evaluate::{accumulate_shell, shell_cutoff_radius, EvalScratch};
use super::molden::{OrbitalData, Spin};

/// Angstrom per Bohr.
pub const BOHR: f64 = 0.529_177_210_9;

/// Coefficients below this contribute nothing visible and are skipped.
const COEFFICIENT_FLOOR: f64 = 1.0e-7;

#[derive(Debug, Clone)]
pub struct GridSpec {
    /// Corner of the first sample, in Angstrom.
    pub origin: Vec3,
    /// Uniform spacing in Angstrom.
    pub step: f32,
    pub dims: [usize; 3],
}

impl GridSpec {
    pub fn point_count(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    pub fn position(&self, ix: usize, iy: usize, iz: usize) -> Vec3 {
        self.origin + Vec3::new(ix as f32, iy as f32, iz as f32) * self.step
    }
}

/// Sampled orbital values on a uniform grid.
#[derive(Debug, Clone)]
pub struct ScalarField {
    pub spec: GridSpec,
    pub values: Vec<f32>,
}

impl ScalarField {
    #[inline]
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        (iz * self.spec.dims[1] + iy) * self.spec.dims[0] + ix
    }

    #[inline]
    pub fn at(&self, ix: usize, iy: usize, iz: usize) -> f32 {
        self.values[self.index(ix, iy, iz)]
    }

    /// Largest magnitude in the field, used to offer a sensible isovalue.
    pub fn max_abs(&self) -> f32 {
        self.values.iter().fold(0.0f32, |m, v| m.max(v.abs()))
    }
}

/// Choose a grid that covers the molecule with padding, capped in total size.
///
/// The step is relaxed rather than the box truncated: clipping the box would
/// slice through orbital lobes, whereas a coarser step only softens them.
pub fn plan_grid(
    positions: &[Vec3],
    padding: f32,
    preferred_step: f32,
    max_points: usize,
) -> GridSpec {
    let mut min = Vec3::splat(f32::MAX);
    let mut max = Vec3::splat(f32::MIN);
    for &p in positions {
        min = min.min(p);
        max = max.max(p);
    }
    if positions.is_empty() {
        min = Vec3::ZERO;
        max = Vec3::ZERO;
    }
    min -= Vec3::splat(padding);
    max += Vec3::splat(padding);

    let extent = max - min;
    let mut step = preferred_step.max(0.01);

    // Grow the step until the point budget is met.
    loop {
        let dims = [
            ((extent.x / step).ceil() as usize + 1).max(2),
            ((extent.y / step).ceil() as usize + 1).max(2),
            ((extent.z / step).ceil() as usize + 1).max(2),
        ];
        if dims[0] * dims[1] * dims[2] <= max_points {
            return GridSpec {
                origin: min,
                step,
                dims,
            };
        }
        step *= 1.15;
    }
}

/// Sample one molecular orbital onto `spec`.
///
/// Runs on the compute task pool, split into z-slabs.
pub fn sample_orbital(
    data: &OrbitalData,
    spin: Spin,
    orbital_index: usize,
    spec: &GridSpec,
) -> ScalarField {
    let orbitals = data.orbitals(spin);
    let coefficients = &orbitals[orbital_index].coefficients;

    // Everything below works in Bohr, matching the Gaussian exponents.
    let centers: Vec<[f64; 3]> = data
        .positions
        .iter()
        .map(|p| [p.x as f64 / BOHR, p.y as f64 / BOHR, p.z as f64 / BOHR])
        .collect();

    // Shells that actually matter for this orbital.
    let active: Vec<usize> = data
        .shells
        .iter()
        .enumerate()
        .filter(|(index, shell)| {
            let start = data.layout.shell_offsets[*index];
            coefficients[start..start + shell.n_functions()]
                .iter()
                .any(|c| c.abs() > COEFFICIENT_FLOOR)
        })
        .map(|(index, _)| index)
        .collect();

    // `nz` is implied by the chunking below, which walks one z-slab per task.
    let [nx, ny, _nz] = spec.dims;
    let step = spec.step as f64 / BOHR;
    let origin = [
        spec.origin.x as f64 / BOHR,
        spec.origin.y as f64 / BOHR,
        spec.origin.z as f64 / BOHR,
    ];

    let mut values = vec![0.0f32; spec.point_count()];
    let slab = nx * ny;

    ComputeTaskPool::get_or_init(TaskPool::default).scope(|scope| {
        for (iz, plane) in values.chunks_mut(slab).enumerate() {
            let centers = &centers;
            let active = &active;
            scope.spawn(async move {
                let mut scratch = EvalScratch::default();
                let z = origin[2] + iz as f64 * step;
                for iy in 0..ny {
                    let y = origin[1] + iy as f64 * step;
                    for ix in 0..nx {
                        let x = origin[0] + ix as f64 * step;
                        let mut value = 0.0;
                        for &index in active {
                            let shell = &data.shells[index];
                            let c = centers[shell.center];
                            let d = [x - c[0], y - c[1], z - c[2]];
                            let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                            let cut = shell_cutoff_radius(shell);
                            if r2 > cut * cut {
                                continue;
                            }
                            let start = data.layout.shell_offsets[index];
                            let weights = &coefficients[start..start + shell.n_functions()];
                            accumulate_shell(shell, d, weights, &mut scratch, &mut value);
                        }
                        plane[iy * nx + ix] = value as f32;
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

    #[test]
    fn grid_covers_the_molecule_with_padding() {
        let positions = vec![Vec3::new(-1.0, 0.0, 0.0), Vec3::new(1.0, 2.0, 0.0)];
        let spec = plan_grid(&positions, 3.0, 0.5, 10_000_000);
        assert!(spec.origin.x <= -4.0 + 1e-6);
        assert!(spec.origin.y <= -3.0 + 1e-6);
        let far = spec.position(spec.dims[0] - 1, spec.dims[1] - 1, spec.dims[2] - 1);
        assert!(far.x >= 4.0 - 1e-6, "far corner {far:?}");
        assert!(far.y >= 5.0 - 1e-6, "far corner {far:?}");
    }

    #[test]
    fn the_point_budget_is_respected_by_relaxing_the_step() {
        let positions = vec![Vec3::splat(-10.0), Vec3::splat(10.0)];
        let fine = plan_grid(&positions, 3.0, 0.05, 200_000);
        assert!(fine.point_count() <= 200_000, "{}", fine.point_count());
        assert!(fine.step > 0.05, "step should have been relaxed");
    }

    #[test]
    fn an_empty_molecule_still_produces_a_usable_grid() {
        let spec = plan_grid(&[], 2.0, 0.5, 1_000_000);
        assert!(spec.dims.iter().all(|&d| d >= 2));
        assert!(spec.point_count() > 0);
    }
}
