//! Isosurface extraction from a sampled scalar field.
//!
//! Marching *tetrahedra* rather than marching cubes.  Each cell is split into
//! six tetrahedra along a body diagonal, and a tetrahedron has only three
//! distinct crossing cases (one corner inside, two inside, three inside),
//! which can be enumerated by inspection instead of transcribed from a
//! 256-entry table.  It costs roughly twice the triangles of marching cubes
//! and, unlike marching cubes, has no ambiguous-face cases to resolve.
//!
//! Vertex normals come from the field gradient rather than from triangle
//! faces, which also hides the tetrahedral tessellation pattern.

use bevy::prelude::*;
use bevy::tasks::{ComputeTaskPool, TaskPool};

use super::grid::ScalarField;

/// Triangle soup ready to become a Bevy mesh.
#[derive(Debug, Default, Clone)]
pub struct SurfaceMesh {
    pub positions: Vec<[f32; 3]>,
    pub normals: Vec<[f32; 3]>,
}

impl SurfaceMesh {
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    pub fn triangle_count(&self) -> usize {
        self.positions.len() / 3
    }
}

/// Corner offsets, indexed so that bit 0 is x, bit 1 is y, bit 2 is z.
const CORNERS: [[usize; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
    [0, 1, 1],
];

/// Six tetrahedra sharing the 0-6 body diagonal.  This decomposition tiles the
/// cube exactly and matches across neighbouring cells, so no cracks appear.
const TETRAHEDRA: [[usize; 4]; 6] = [
    [0, 1, 2, 6],
    [0, 2, 3, 6],
    [0, 3, 7, 6],
    [0, 7, 4, 6],
    [0, 4, 5, 6],
    [0, 5, 1, 6],
];

/// Extract the surface where the field equals `isovalue`.
///
/// A positive `isovalue` yields the lobe where psi is positive; pass a negative
/// value for the other lobe.  Normals always point away from the enclosed
/// volume, so both lobes are lit consistently.
pub fn extract(field: &ScalarField, isovalue: f32) -> SurfaceMesh {
    let [nx, ny, nz] = field.spec.dims;
    let mut mesh = SurfaceMesh::default();
    if nx < 2 || ny < 2 || nz < 2 || isovalue == 0.0 {
        return mesh;
    }

    // For a negative isovalue the inside/outside sense flips, so the gradient
    // has to be negated to keep normals pointing outward.
    let normal_sign = if isovalue > 0.0 { -1.0 } else { 1.0 };

    let mut values = [0.0f32; 8];
    let mut points = [Vec3::ZERO; 8];
    let mut gradients = [Vec3::ZERO; 8];

    for iz in 0..nz - 1 {
        for iy in 0..ny - 1 {
            for ix in 0..nx - 1 {
                for (c, offset) in CORNERS.iter().enumerate() {
                    let (x, y, z) = (ix + offset[0], iy + offset[1], iz + offset[2]);
                    values[c] = field.at(x, y, z);
                    points[c] = field.spec.position(x, y, z);
                    gradients[c] = gradient(field, x, y, z);
                }

                // Cheap rejection: no crossing if every corner is on one side.
                let mut any_inside = false;
                let mut any_outside = false;
                for &v in &values {
                    if inside(v, isovalue) {
                        any_inside = true;
                    } else {
                        any_outside = true;
                    }
                }
                if !(any_inside && any_outside) {
                    continue;
                }

                for tet in &TETRAHEDRA {
                    emit_tetrahedron(
                        tet,
                        &values,
                        &points,
                        &gradients,
                        isovalue,
                        normal_sign,
                        &mut mesh,
                    );
                }
            }
        }
    }

    mesh
}

/// A point is "inside" the surface when it is further from zero than the
/// isovalue, on the isovalue's own side.
#[inline]
fn inside(value: f32, isovalue: f32) -> bool {
    if isovalue > 0.0 {
        value >= isovalue
    } else {
        value <= isovalue
    }
}

/// Central-difference gradient, falling back to one-sided at the boundary.
fn gradient(field: &ScalarField, ix: usize, iy: usize, iz: usize) -> Vec3 {
    let [nx, ny, nz] = field.spec.dims;
    let h = field.spec.step;

    let axis = |lo: usize, hi: usize, n: usize, sample: &dyn Fn(usize) -> f32| -> f32 {
        let _ = n;
        (sample(hi) - sample(lo)) / ((hi - lo) as f32 * h)
    };

    let gx = {
        let lo = ix.saturating_sub(1);
        let hi = (ix + 1).min(nx - 1);
        if hi == lo {
            0.0
        } else {
            axis(lo, hi, nx, &|i| field.at(i, iy, iz))
        }
    };
    let gy = {
        let lo = iy.saturating_sub(1);
        let hi = (iy + 1).min(ny - 1);
        if hi == lo {
            0.0
        } else {
            axis(lo, hi, ny, &|i| field.at(ix, i, iz))
        }
    };
    let gz = {
        let lo = iz.saturating_sub(1);
        let hi = (iz + 1).min(nz - 1);
        if hi == lo {
            0.0
        } else {
            axis(lo, hi, nz, &|i| field.at(ix, iy, i))
        }
    };

    Vec3::new(gx, gy, gz)
}

#[allow(clippy::too_many_arguments)]
fn emit_tetrahedron(
    tet: &[usize; 4],
    values: &[f32; 8],
    points: &[Vec3; 8],
    gradients: &[Vec3; 8],
    isovalue: f32,
    normal_sign: f32,
    mesh: &mut SurfaceMesh,
) {
    // Partition the four corners.
    let mut inside_corners = [0usize; 4];
    let mut outside_corners = [0usize; 4];
    let (mut n_in, mut n_out) = (0usize, 0usize);
    for &c in tet {
        if inside(values[c], isovalue) {
            inside_corners[n_in] = c;
            n_in += 1;
        } else {
            outside_corners[n_out] = c;
            n_out += 1;
        }
    }

    let interp = |a: usize, b: usize| -> ([f32; 3], [f32; 3]) {
        let va = values[a];
        let vb = values[b];
        let denominator = vb - va;
        let t = if denominator.abs() < f32::EPSILON {
            0.5
        } else {
            ((isovalue - va) / denominator).clamp(0.0, 1.0)
        };
        let position = points[a].lerp(points[b], t);
        let gradient = gradients[a].lerp(gradients[b], t) * normal_sign;
        let normal = gradient.normalize_or_zero();
        (position.into(), normal.into())
    };

    let mut push = |a: ([f32; 3], [f32; 3]), b: ([f32; 3], [f32; 3]), c: ([f32; 3], [f32; 3])| {
        mesh.positions.extend_from_slice(&[a.0, b.0, c.0]);
        mesh.normals.extend_from_slice(&[a.1, b.1, c.1]);
    };

    match (n_in, n_out) {
        // Entirely on one side: nothing to draw.
        (0, _) | (_, 0) => {}

        // One corner inside: a single triangle cutting off that corner.
        (1, 3) => {
            let a = inside_corners[0];
            push(
                interp(a, outside_corners[0]),
                interp(a, outside_corners[1]),
                interp(a, outside_corners[2]),
            );
        }

        // Three inside: the mirror case, one triangle around the lone outsider.
        (3, 1) => {
            let b = outside_corners[0];
            push(
                interp(b, inside_corners[0]),
                interp(b, inside_corners[1]),
                interp(b, inside_corners[2]),
            );
        }

        // Two inside, two outside: the crossing is a quad, split into two
        // triangles sharing the diagonal.
        (2, 2) => {
            let (a, b) = (inside_corners[0], inside_corners[1]);
            let (c, d) = (outside_corners[0], outside_corners[1]);
            let ac = interp(a, c);
            let ad = interp(a, d);
            let bc = interp(b, c);
            let bd = interp(b, d);
            push(ac, ad, bd);
            push(ac, bd, bc);
        }

        _ => unreachable!("a tetrahedron has exactly four corners"),
    }
}

/// One line segment of the wireframe net.
pub type Segment = [Vec3; 2];

/// Never place planes closer than this fraction of a grid cell.  Below it the
/// net gains no detail the sampling can support, and the cost keeps climbing.
const MIN_PLANE_STRIDE: f32 = 0.25;
/// Bound on planes per axis, so a very fine setting cannot stall a frame.
const MAX_PLANES_PER_AXIS: usize = 240;

/// Iso-contours drawn on evenly spaced planes: a chicken-wire net.
///
/// The triangle mesh cannot serve as an adjustable wireframe -- its edge count
/// is fixed by the sampling resolution, and thinning it out leaves a
/// random-looking scatter.  Contouring planes at a chosen world-space
/// `plane_spacing` instead gives a regular net whose density is a free
/// parameter, independent of how finely the orbital was sampled.
///
/// `plane_spacing` may be *smaller* than the grid step: planes are positioned
/// at fractional grid coordinates and the field is interpolated along the
/// plane normal, so the net can be finer than the sampling lattice.
///
/// Each plane is contoured with marching squares, which needs no case table:
/// a cell has a crossing on an edge whenever its two ends straddle the
/// isovalue, and those crossings are joined in pairs.
pub fn contour_net(field: &ScalarField, isovalue: f32, plane_spacing: f32) -> Vec<Segment> {
    let dims = field.spec.dims;
    let mut out = Vec::new();
    if dims.iter().any(|&d| d < 2) || isovalue == 0.0 || plane_spacing <= 0.0 {
        return out;
    }

    // Convert the world-space spacing into grid cells.
    let stride = (plane_spacing / field.spec.step).max(MIN_PLANE_STRIDE);

    // Enumerate the planes first, then contour them in parallel: each plane is
    // independent, and at the finer settings there are several hundred of them.
    // Extraction runs synchronously inside the surface rebuild, so it has to
    // fit a frame rather than merely be fast on average.
    let mut planes: Vec<(usize, f32)> = Vec::new();
    for axis in 0..3 {
        let span = (dims[axis] - 1) as f32;
        let mut placed = 0usize;
        let mut offset = stride * 0.5; // keep planes off the box faces
        while offset < span && placed < MAX_PLANES_PER_AXIS {
            planes.push((axis, offset));
            offset += stride;
            placed += 1;
        }
    }
    if planes.is_empty() {
        return out;
    }

    // `get_or_init` rather than `get` so unit tests need no Bevy app.
    let pool = ComputeTaskPool::get_or_init(TaskPool::default);
    let chunk = planes.len().div_ceil(pool.thread_num().max(1) * 4).max(1);

    let batches = pool.scope(|scope| {
        for group in planes.chunks(chunk) {
            scope.spawn(async move {
                let mut scratch = Vec::new();
                let mut local = Vec::new();
                for &(axis, offset) in group {
                    contour_plane(field, isovalue, axis, offset, &mut scratch, &mut local);
                }
                local
            });
        }
    });

    for batch in batches {
        out.extend(batch);
    }
    out
}

/// Grid index triple for `axis` held at `along`, with the two in-plane
/// coordinates `u` and `v`.
#[inline]
fn plane_index(axis: usize, along: usize, u: usize, v: usize) -> (usize, usize, usize) {
    match axis {
        0 => (along, u, v),
        1 => (u, along, v),
        _ => (u, v, along),
    }
}

/// Marching squares across one plane of constant `axis` coordinate.
///
/// `offset` is a fractional grid coordinate along `axis`, so the plane need not
/// coincide with a sampling plane; values are interpolated between the two that
/// bracket it.
///
/// The plane's values are gathered once into `scratch` and cells then read from
/// it, rather than each cell re-interpolating its four corners: neighbouring
/// cells share corners, so the naive form does roughly four times the
/// interpolation work.  Corner *positions* are built only inside the crossing
/// branch, because the overwhelming majority of cells lie entirely inside or
/// outside the surface and need no geometry at all.
fn contour_plane(
    field: &ScalarField,
    isovalue: f32,
    axis: usize,
    offset: f32,
    scratch: &mut Vec<f32>,
    out: &mut Vec<Segment>,
) {
    let dims = field.spec.dims;
    let n_along = dims[axis];
    let (n_u, n_v) = match axis {
        0 => (dims[1], dims[2]),
        1 => (dims[0], dims[2]),
        _ => (dims[0], dims[1]),
    };
    if n_u < 2 || n_v < 2 {
        return;
    }

    let lower = offset.floor().max(0.0) as usize;
    let upper = (lower + 1).min(n_along - 1);
    let blend = offset - lower as f32;

    // Gather the plane once.
    scratch.clear();
    scratch.reserve(n_u * n_v);
    for v in 0..n_v {
        for u in 0..n_u {
            let (ax, ay, az) = plane_index(axis, lower, u, v);
            let (bx, by, bz) = plane_index(axis, upper, u, v);
            let a = field.at(ax, ay, az);
            let b = field.at(bx, by, bz);
            scratch.push(a + (b - a) * blend);
        }
    }

    // In-plane basis, so corner positions are two multiply-adds.
    let step = field.spec.step;
    let plane_coordinate = field.spec.origin[axis] + offset * step;
    let (u_axis, v_axis) = match axis {
        0 => (1usize, 2usize),
        1 => (0, 2),
        _ => (0, 1),
    };
    let corner_position = |u: usize, v: usize| -> Vec3 {
        let mut p = Vec3::ZERO;
        p[axis] = plane_coordinate;
        p[u_axis] = field.spec.origin[u_axis] + u as f32 * step;
        p[v_axis] = field.spec.origin[v_axis] + v as f32 * step;
        p
    };

    // Corner order around the cell, so consecutive entries share an edge.
    const RING: [(usize, usize); 4] = [(0, 0), (1, 0), (1, 1), (0, 1)];

    for v in 0..n_v - 1 {
        let row = v * n_u;
        let row_above = row + n_u;
        for u in 0..n_u - 1 {
            let values = [
                scratch[row + u],
                scratch[row + u + 1],
                scratch[row_above + u + 1],
                scratch[row_above + u],
            ];

            // Reject the common case before touching any geometry.
            let first = inside(values[0], isovalue);
            if values[1..]
                .iter()
                .all(|&value| inside(value, isovalue) == first)
            {
                continue;
            }

            let mut crossings = [Vec3::ZERO; 4];
            let mut found = 0usize;
            for corner in 0..4 {
                let next = (corner + 1) % 4;
                let (a, b) = (values[corner], values[next]);
                if inside(a, isovalue) == inside(b, isovalue) {
                    continue;
                }
                let denominator = b - a;
                let t = if denominator.abs() < f32::EPSILON {
                    0.5
                } else {
                    ((isovalue - a) / denominator).clamp(0.0, 1.0)
                };
                let (du0, dv0) = RING[corner];
                let (du1, dv1) = RING[next];
                let p0 = corner_position(u + du0, v + dv0);
                let p1 = corner_position(u + du1, v + dv1);
                crossings[found] = p0.lerp(p1, t);
                found += 1;
            }

            match found {
                2 => out.push([crossings[0], crossings[1]]),
                // Saddle: two opposite corners inside.  Either pairing is a
                // valid contour; picking one consistently avoids flicker.
                4 => {
                    out.push([crossings[0], crossings[1]]);
                    out.push([crossings[2], crossings[3]]);
                }
                _ => {}
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::orbitals::grid::GridSpec;

    /// A cone-shaped lobe peaking at the origin: `f(p) = peak - |p|`.
    ///
    /// Modelled on a real orbital lobe, which is largest at its centre and
    /// falls off outward, so the isovalue semantics match production use:
    /// the surface at isovalue `v` is a sphere of radius `peak - v`, and
    /// raising the isovalue shrinks it.
    ///
    /// `sign` flips the whole field to produce a negative lobe.
    fn lobe_field(n: usize, extent: f32, peak: f32, sign: f32) -> ScalarField {
        let step = 2.0 * extent / (n - 1) as f32;
        let spec = GridSpec {
            origin: Vec3::splat(-extent),
            step,
            dims: [n, n, n],
        };
        let mut values = Vec::with_capacity(spec.point_count());
        for iz in 0..n {
            for iy in 0..n {
                for ix in 0..n {
                    let p = spec.position(ix, iy, iz);
                    values.push(sign * (peak - p.length()));
                }
            }
        }
        ScalarField { spec, values }
    }

    fn positive_lobe(n: usize) -> ScalarField {
        lobe_field(n, 4.0, 4.0, 1.0)
    }

    fn negative_lobe(n: usize) -> ScalarField {
        lobe_field(n, 4.0, 4.0, -1.0)
    }

    fn surface_area(mesh: &SurfaceMesh) -> f32 {
        mesh.positions
            .chunks_exact(3)
            .map(|t| {
                let a = Vec3::from(t[0]);
                let b = Vec3::from(t[1]);
                let c = Vec3::from(t[2]);
                (b - a).cross(c - a).length() * 0.5
            })
            .sum()
    }

    /// The self-validating test: a sphere of radius R must come out with area
    /// 4 pi R^2, which no bookkeeping error in the case table survives.
    /// The self-validating test: a sphere of radius R must come out with area
    /// 4 pi R^2, which no bookkeeping error in the case analysis survives.
    #[test]
    fn a_sphere_has_the_analytic_surface_area() {
        let field = positive_lobe(65);
        let radius = 2.0f32;
        let mesh = extract(&field, 4.0 - radius);

        assert!(!mesh.is_empty(), "sphere produced no triangles");
        let area = surface_area(&mesh);
        let expected = 4.0 * std::f32::consts::PI * radius * radius;
        let error = (area - expected).abs() / expected;
        assert!(
            error < 0.02,
            "area {area} vs analytic {expected} ({:.1}% off)",
            error * 100.0
        );
    }

    /// The negative lobe must mesh identically to the positive one.  Both are
    /// drawn for every orbital, so a sign error here would show as one lobe
    /// missing or inside out.
    #[test]
    fn the_negative_lobe_matches_the_positive_one() {
        let radius = 2.0f32;
        let positive = surface_area(&extract(&positive_lobe(65), 4.0 - radius));
        let negative = surface_area(&extract(&negative_lobe(65), -(4.0 - radius)));
        assert!(
            (positive - negative).abs() / positive < 1e-4,
            "positive lobe area {positive} != negative lobe area {negative}"
        );
    }

    /// Normals must point away from the enclosed volume, for both signs.
    #[test]
    fn normals_point_outward_for_both_lobes() {
        for (field, isovalue, label) in [
            (positive_lobe(49), 2.0f32, "positive"),
            (negative_lobe(49), -2.0f32, "negative"),
        ] {
            let mesh = extract(&field, isovalue);
            let mut checked = 0;
            for (position, normal) in mesh.positions.iter().zip(&mesh.normals) {
                let p = Vec3::from(*position);
                let n = Vec3::from(*normal);
                if n.length() < 0.5 || p.length() < 1e-3 {
                    continue;
                }
                assert!(
                    p.normalize().dot(n) > 0.8,
                    "{label} lobe: normal {n:?} is not outward at {p:?}"
                );
                checked += 1;
            }
            assert!(
                checked > 100,
                "{label} lobe: only checked {checked} vertices"
            );
        }
    }

    #[test]
    fn an_isovalue_beyond_the_field_range_yields_nothing() {
        let field = positive_lobe(33);
        assert!(extract(&field, 99.0).is_empty());
        assert!(extract(&field, -99.0).is_empty());
    }

    #[test]
    fn a_larger_isovalue_encloses_a_smaller_surface() {
        let field = positive_lobe(65);
        // isovalue 1.0 -> radius 3;  isovalue 2.5 -> radius 1.5
        let big = surface_area(&extract(&field, 1.0));
        let small = surface_area(&extract(&field, 2.5));
        assert!(
            big > small * 2.0,
            "radius 3 area {big} should dwarf radius 1.5 area {small}"
        );
    }

    /// The net must lie on the surface: every endpoint should sit at the
    /// isovalue, i.e. at the sphere's radius.
    #[test]
    fn contour_net_endpoints_lie_on_the_surface() {
        let field = positive_lobe(65);
        let radius = 2.0f32;
        let net = contour_net(&field, 4.0 - radius, 0.5);

        assert!(!net.is_empty(), "no contour segments produced");
        for segment in &net {
            for point in segment {
                let r = point.length();
                assert!(
                    (r - radius).abs() < 0.15,
                    "segment endpoint at radius {r}, expected {radius}"
                );
            }
        }
    }

    /// Density is the whole point of the feature: closer planes must draw
    /// more lines.
    #[test]
    fn a_finer_spacing_gives_a_denser_net() {
        let field = positive_lobe(65);
        let dense = contour_net(&field, 2.0, 0.2).len();
        let sparse = contour_net(&field, 2.0, 1.0).len();
        assert!(
            dense > sparse * 2,
            "0.2 spacing gave {dense} segments, 1.0 gave {sparse}"
        );
    }

    /// The net must be able to go finer than the sampling grid, which is what
    /// an integer plane stride could not express.
    #[test]
    fn spacing_below_the_grid_step_still_refines_the_net() {
        let field = positive_lobe(65);
        let step = field.spec.step;
        let at_grid = contour_net(&field, 2.0, step).len();
        let sub_grid = contour_net(&field, 2.0, step * 0.4).len();
        assert!(
            sub_grid > at_grid,
            "sub-grid spacing gave {sub_grid} segments, grid spacing gave {at_grid}"
        );
    }

    #[test]
    fn the_net_works_for_the_negative_lobe_too() {
        let positive = contour_net(&positive_lobe(49), 2.0, 0.5).len();
        let negative = contour_net(&negative_lobe(49), -2.0, 0.5).len();
        assert_eq!(positive, negative, "lobes should mirror exactly");
    }

    #[test]
    fn an_isovalue_beyond_the_range_produces_no_net() {
        assert!(contour_net(&positive_lobe(33), 99.0, 0.5).is_empty());
    }
}
