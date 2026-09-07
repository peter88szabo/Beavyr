//! Molecular orbital viewer: Molden import, basis evaluation, isosurfaces.
//!
//! Flow: a `.molden` file is parsed into [`molden::OrbitalData`]; selecting an
//! orbital spawns a background sampling job producing a [`grid::ScalarField`];
//! the field is meshed at the current isovalue into two lobe entities.
//!
//! The field is retained for the currently displayed orbital, so moving the
//! isovalue slider only re-runs the mesher.

pub mod basis;
pub mod cart2sph;
pub mod density;
pub mod evaluate;
pub mod grid;
pub mod isosurface;
pub mod molden;
pub mod ui;

use std::path::PathBuf;
use std::sync::Arc;

use bevy::asset::RenderAssetUsages;
use bevy::camera::visibility::RenderLayers;
use bevy::mesh::PrimitiveTopology;
use bevy::prelude::*;
use bevy::tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task};

use crate::scene::LAYER_MAIN;
use density::{sample_density, DensityMatrices, Quantity};
use grid::{plan_grid, sample_orbital, ScalarField};
use molden::{OrbitalData, Spin};

/// Padding around the molecule's bounding box, in Angstrom.
const GRID_PADDING: f32 = 3.0;
/// Upper bound on grid points, to keep the worst case bounded.
const MAX_GRID_POINTS: usize = 2_000_000;

/// Grid spacing in Angstrom for each surface-resolution setting.
///
/// Higher setting means finer sampling and a smoother surface, at a cost that
/// grows with the cube of the refinement.  Index 0 is unused so the slider can
/// read 1..=6 like the atom-resolution controls.
const RESOLUTION_SPACING: [f32; 8] = [0.0, 0.40, 0.30, 0.25, 0.20, 0.15, 0.12, 0.09];

/// Lowest and highest selectable surface resolution.
pub const MIN_RESOLUTION: u32 = 1;
pub const MAX_RESOLUTION: u32 = 7;
/// Matches the 0.20 A spacing the viewer used before the control existed.
const DEFAULT_RESOLUTION: u32 = 4;

/// Plane spacing for the wireframe net is `MESH_SPACING_BASE / density`, so a
/// higher slider value means a denser net -- the direction people expect.
/// The top of the range lands well below a typical grid step, which an integer
/// plane stride could not reach.
const MESH_SPACING_BASE: f32 = 2.0;
pub const MIN_MESH_DENSITY: u32 = 10;
pub const MAX_MESH_DENSITY: u32 = 30;
const DEFAULT_MESH_DENSITY: u32 = 15;

/// Identifies a sampled field.  The resolution is part of the identity, so
/// changing it invalidates the cache and triggers a resample rather than a
/// bare remesh.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct FieldKey {
    quantity: Quantity,
    spin: Spin,
    orbital: usize,
    resolution: u32,
}

/// What a sampling job produces.
///
/// A density job also returns the matrices it had to build, so the next one
/// reuses them instead of paying for construction again.
struct SampleOutput {
    field: ScalarField,
    density: Option<Arc<DensityMatrices>>,
}

/// How the isosurface is drawn.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum SurfaceStyle {
    /// Shaded triangles.
    #[default]
    Solid,
    /// Chicken-wire net of iso-contours on evenly spaced grid planes.
    Mesh,
}

/// Gizmo group for the wireframe net.
///
/// Gizmos rather than a line mesh, so the net gets a real width in screen
/// pixels; a `LineList` mesh is stuck at one pixel.  The segments are cached
/// and only recomputed when the field, isovalue or spacing change, so redrawing
/// them each frame costs no more than the bond lines already do.
#[derive(Default, Reflect, GizmoConfigGroup)]
pub struct OrbitalMeshGizmos;

/// Which lobe an entity draws.
#[derive(Component, Clone, Copy, PartialEq, Eq, Debug)]
pub enum OrbitalLobe {
    Positive,
    Negative,
}

#[derive(Resource)]
pub struct OrbitalState {
    pub data: Option<Arc<OrbitalData>>,
    pub file_name: Option<String>,
    pub last_dir: Option<PathBuf>,
    pub warning: Option<String>,
    /// Where this orbital set came from, when it was not a file the user
    /// opened: the method line of the excited-state run that produced it.
    /// `None` for a loaded file, whose name is shown instead.
    pub provenance: Option<String>,
    /// The Molden text behind a computed set, kept so it can be saved. A run's
    /// scratch directory is discarded, so this is the only copy.
    pub molden_text: Option<String>,

    /// Which scalar field the panel is showing.
    pub quantity: Quantity,
    pub spin: Spin,
    /// Orbital chosen in the UI, as an index into the current spin set.
    pub selected: Option<usize>,
    /// AO-basis density matrices, built on first use and then reused.
    density: Option<Arc<DensityMatrices>>,

    pub isovalue: f32,
    /// Sampling density of the scalar field; see [`RESOLUTION_SPACING`].
    pub resolution: u32,
    pub style: SurfaceStyle,
    /// Wireframe density; higher draws more contour planes.
    pub mesh_density: u32,
    pub mesh_width: f32,
    /// Cached wireframe segments for each lobe.
    mesh_positive: Vec<isosurface::Segment>,
    mesh_negative: Vec<isosurface::Segment>,
    pub show_surface: bool,
    pub opacity: f32,
    pub positive_color: Color,
    pub negative_color: Color,

    /// Sampled field for `field_key`, kept so the isovalue slider only remeshes.
    field: Option<ScalarField>,
    field_key: Option<FieldKey>,
    job: Option<Task<SampleOutput>>,
    job_key: Option<FieldKey>,
    surface_dirty: bool,
    /// Triangles (solid) or segments (mesh) last extracted, shown in the panel.
    pub triangle_count: usize,
}

impl Default for OrbitalState {
    fn default() -> Self {
        Self {
            data: None,
            file_name: None,
            last_dir: None,
            warning: None,
            provenance: None,
            molden_text: None,
            quantity: Quantity::Orbital,
            spin: Spin::Alpha,
            selected: None,
            density: None,
            isovalue: 0.10,
            resolution: DEFAULT_RESOLUTION,
            style: SurfaceStyle::Solid,
            mesh_density: DEFAULT_MESH_DENSITY,
            mesh_width: 1.5,
            mesh_positive: Vec::new(),
            mesh_negative: Vec::new(),
            show_surface: true,
            opacity: 1.0,
            positive_color: Color::srgb(0.20, 0.45, 0.95),
            negative_color: Color::srgb(0.90, 0.25, 0.25),
            field: None,
            field_key: None,
            job: None,
            job_key: None,
            surface_dirty: false,
            triangle_count: 0,
        }
    }
}

impl OrbitalState {
    /// True while a sampling job is in flight.
    pub fn is_busy(&self) -> bool {
        self.job.is_some()
    }

    /// Peak magnitude of the current field, for isovalue guidance.
    pub fn field_peak(&self) -> Option<f32> {
        self.field.as_ref().map(|f| f.max_abs())
    }

    /// Request a remesh without resampling.
    pub fn mark_surface_dirty(&mut self) {
        self.surface_dirty = true;
    }

    /// World-space distance between wireframe contour planes.
    ///
    /// Inverted from the slider so that a higher setting is a denser net.
    pub fn mesh_plane_spacing(&self) -> f32 {
        let density = self.mesh_density.clamp(MIN_MESH_DENSITY, MAX_MESH_DENSITY);
        MESH_SPACING_BASE / density as f32
    }

    /// Grid spacing in Angstrom for the current resolution setting.
    pub fn grid_spacing(&self) -> f32 {
        let index = self.resolution.clamp(MIN_RESOLUTION, MAX_RESOLUTION) as usize;
        RESOLUTION_SPACING[index]
    }

    /// The grid that would be used for `data` at the current resolution.
    ///
    /// Exposed so the panel can show the real cost, including any relaxation
    /// the point budget forces.
    pub fn planned_grid(&self, data: &OrbitalData) -> grid::GridSpec {
        plan_grid(
            &data.positions,
            GRID_PADDING,
            self.grid_spacing(),
            MAX_GRID_POINTS,
        )
    }

    /// Drop the loaded file and everything derived from it.
    /// Takes on a parsed orbital set, however it arrived.
    ///
    /// One place for the several fields that have to move together -- the
    /// data, what to show first, and the isovalue that suits it -- so a set
    /// loaded from a file and one computed by the spectrum tool cannot end up
    /// in different states.
    ///
    /// Opens on the HOMO, which is what one usually wants to look at first and
    /// is also the orbital most excitations start from.
    pub fn adopt(
        &mut self,
        data: OrbitalData,
        file_name: Option<String>,
        provenance: Option<String>,
        molden_text: Option<String>,
    ) {
        let spin = Spin::Alpha;
        self.selected = data.homo_index(spin);
        self.spin = spin;
        self.quantity = Quantity::Orbital;
        self.isovalue = Quantity::Orbital.default_isovalue();
        self.file_name = file_name;
        self.provenance = provenance;
        self.molden_text = molden_text;
        self.warning = None;
        // A previous set's cached field and density describe different
        // orbitals entirely.
        self.density = None;
        self.field = None;
        self.field_key = None;
        self.job = None;
        self.job_key = None;
        self.data = Some(Arc::new(data));
        self.mark_surface_dirty();
    }

    pub fn clear(&mut self) {
        self.data = None;
        self.file_name = None;
        self.provenance = None;
        self.molden_text = None;
        self.selected = None;
        self.density = None;
        self.field = None;
        self.field_key = None;
        self.job = None;
        self.job_key = None;
        self.surface_dirty = true;
    }
}

pub struct OrbitalPlugin;

impl Plugin for OrbitalPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<OrbitalState>()
            .init_gizmo_group::<OrbitalMeshGizmos>()
            .add_systems(
                Update,
                (start_orbital_job, poll_orbital_job, rebuild_orbital_surface).chain(),
            )
            .add_systems(Update, (configure_orbital_mesh_gizmos, draw_orbital_mesh));
    }
}

/// Spawn a sampling job when the selection differs from what is cached or
/// already running.
fn start_orbital_job(mut state: ResMut<OrbitalState>) {
    let Some(data) = state.data.clone() else {
        return;
    };
    let quantity = state.quantity;
    if !quantity.is_available(&data) {
        return;
    }

    // Only an orbital needs a selection; the density fields are properties of
    // the whole wavefunction.
    let selected = match quantity {
        Quantity::Orbital => match state.selected {
            Some(index) if index < data.orbitals(state.spin).len() => index,
            _ => return,
        },
        _ => 0,
    };

    let key = FieldKey {
        quantity,
        spin: state.spin,
        orbital: selected,
        resolution: state.resolution,
    };
    if state.field_key == Some(key) || state.job_key == Some(key) {
        return;
    }

    let spec = state.planned_grid(&data);
    let spin = state.spin;
    let existing = state.density.clone();

    let task = AsyncComputeTaskPool::get().spawn(async move {
        match quantity {
            Quantity::Orbital => SampleOutput {
                field: sample_orbital(&data, spin, selected, &spec),
                density: existing,
            },
            _ => {
                // Building the matrices is the expensive part of a first
                // density request, so it happens here rather than on the main
                // thread at load time.
                let matrices = match existing {
                    Some(matrices) => matrices,
                    None => Arc::new(DensityMatrices::build(&data)),
                };
                let field = sample_density(&data, &matrices, quantity, &spec);
                SampleOutput {
                    field,
                    density: Some(matrices),
                }
            }
        }
    });

    state.job = Some(task);
    state.job_key = Some(key);
}

/// Adopt a finished field.  Non-blocking: the viewport keeps rendering while
/// the job runs.
fn poll_orbital_job(mut state: ResMut<OrbitalState>) {
    let Some(task) = state.job.as_mut() else {
        return;
    };
    let Some(output) = block_on(future::poll_once(task)) else {
        return;
    };
    state.job = None;
    state.field_key = state.job_key.take();
    state.field = Some(output.field);
    if output.density.is_some() {
        state.density = output.density;
    }
    state.surface_dirty = true;
}

/// Re-extract the two lobes.  Cheap enough to run synchronously, which is what
/// makes the isovalue slider live.
fn rebuild_orbital_surface(
    mut commands: Commands,
    mut state: ResMut<OrbitalState>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    existing: Query<Entity, With<OrbitalLobe>>,
) {
    if !state.surface_dirty {
        return;
    }
    state.surface_dirty = false;

    for entity in existing.iter() {
        commands.entity(entity).despawn();
    }

    // Read the settings out before borrowing the field, so the accumulator
    // below does not fight the immutable borrow of `state.field`.
    let isovalue = state.isovalue;
    let opacity = state.opacity;
    let colors = [state.positive_color, state.negative_color];
    let signed = state.quantity.is_signed();
    let draw =
        state.show_surface && (state.quantity != Quantity::Orbital || state.selected.is_some());

    let style = state.style;
    let mesh_spacing = state.mesh_plane_spacing();

    // Wireframe segments are rebuilt here too, so both styles share one
    // invalidation path.
    let mut net_positive = Vec::new();
    let mut net_negative = Vec::new();
    let mut triangles = 0usize;

    if let (Some(field), true, SurfaceStyle::Mesh) = (state.field.as_ref(), draw, style) {
        net_positive = isosurface::contour_net(field, isovalue, mesh_spacing);
        if signed {
            net_negative = isosurface::contour_net(field, -isovalue, mesh_spacing);
        }
        triangles = net_positive.len() + net_negative.len();
    }

    if let (Some(field), true, SurfaceStyle::Solid) = (state.field.as_ref(), draw, style) {
        let alpha_mode = if opacity >= 0.999 {
            AlphaMode::Opaque
        } else {
            AlphaMode::Blend
        };

        for (lobe, iso, color) in [
            (OrbitalLobe::Positive, isovalue, colors[0]),
            (OrbitalLobe::Negative, -isovalue, colors[1]),
        ] {
            // Total density is non-negative, so it has no second lobe.
            if !signed && lobe == OrbitalLobe::Negative {
                continue;
            }
            let surface = isosurface::extract(field, iso);
            if surface.is_empty() {
                continue;
            }
            triangles += surface.triangle_count();

            let mut mesh = Mesh::new(
                PrimitiveTopology::TriangleList,
                RenderAssetUsages::default(),
            );
            mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, surface.positions);
            mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, surface.normals);

            let base = color.to_srgba();
            let material = materials.add(StandardMaterial {
                base_color: Color::srgba(base.red, base.green, base.blue, opacity),
                perceptual_roughness: 0.45,
                metallic: 0.0,
                alpha_mode,
                // A lobe clipped by the grid box is an open surface, so drawing
                // both sides avoids holes when looking into the cut.
                cull_mode: None,
                ..default()
            });

            commands.spawn((
                Mesh3d(meshes.add(mesh)),
                MeshMaterial3d(material),
                Transform::default(),
                lobe,
                RenderLayers::layer(LAYER_MAIN),
            ));
        }
    }

    state.mesh_positive = net_positive;
    state.mesh_negative = net_negative;
    state.triangle_count = triangles;
}

fn configure_orbital_mesh_gizmos(
    mut config_store: ResMut<GizmoConfigStore>,
    state: Res<OrbitalState>,
) {
    let (config, _) = config_store.config_mut::<OrbitalMeshGizmos>();
    config.enabled = matches!(state.style, SurfaceStyle::Mesh) && state.show_surface;
    config.line.width = state.mesh_width.clamp(0.5, 20.0);
    // Constant screen-space width: a wireframe should stay legible at the back
    // of the molecule, not taper away with distance.
    config.line.perspective = false;
}

fn draw_orbital_mesh(mut gizmos: Gizmos<OrbitalMeshGizmos>, state: Res<OrbitalState>) {
    if !matches!(state.style, SurfaceStyle::Mesh) || !state.show_surface {
        return;
    }
    for segment in &state.mesh_positive {
        gizmos.line(segment[0], segment[1], state.positive_color);
    }
    for segment in &state.mesh_negative {
        gizmos.line(segment[0], segment[1], state.negative_color);
    }
}

#[cfg(test)]
mod adopt_tests {
    use super::*;

    fn behemoth_orbitals() -> OrbitalData {
        let text = include_str!("../../tests/fixtures/behemoth_ch2o_canonicalMO.molden");
        crate::orbitals::molden::parse_molden(text).expect("Behemoth's Molden parses")
    }

    /// A computed set arrives with a provenance line and the text to save,
    /// where a loaded file has a name and neither -- it is already on disk.
    #[test]
    fn a_computed_set_carries_its_provenance_and_its_text() {
        let mut state = OrbitalState::default();
        let text = "molden text".to_string();
        state.adopt(
            behemoth_orbitals(),
            None,
            Some("From this sTDA / Behemoth RKS PBE0 run".to_string()),
            Some(text.clone()),
        );
        assert!(state.data.is_some());
        assert!(state.file_name.is_none(), "there is no file behind it yet");
        assert_eq!(state.molden_text.as_deref(), Some(text.as_str()));
        assert!(state.provenance.as_deref().unwrap().contains("PBE0"));
    }

    #[test]
    fn a_loaded_file_has_a_name_and_no_provenance() {
        let mut state = OrbitalState::default();
        state.adopt(behemoth_orbitals(), Some("mine.molden".to_string()), None, None);
        assert_eq!(state.file_name.as_deref(), Some("mine.molden"));
        assert!(state.provenance.is_none(), "its name is its provenance");
        assert!(state.molden_text.is_none(), "already on disk, nothing to save");
    }

    /// Adopting opens on the HOMO, which is where most excitations start.
    #[test]
    fn adopting_opens_on_the_homo_as_an_orbital() {
        let data = behemoth_orbitals();
        let homo = data.homo_index(Spin::Alpha);
        let mut state = OrbitalState::default();
        state.adopt(data, None, None, None);
        assert_eq!(state.selected, homo);
        assert_eq!(state.quantity, Quantity::Orbital);
        assert_eq!(state.spin, Spin::Alpha);
    }

    /// A second set must not be drawn with the first one's cached field, which
    /// describes different orbitals entirely.
    #[test]
    fn adopting_discards_the_previous_sets_caches() {
        let mut state = OrbitalState::default();
        state.adopt(behemoth_orbitals(), None, None, None);
        state.warning = Some("stale".to_string());
        state.adopt(behemoth_orbitals(), None, Some("run 2".to_string()), None);
        assert!(state.warning.is_none(), "a fresh set clears the old warning");
        assert!(state.surface_dirty, "the surface has to be rebuilt");
    }

    /// Closing a computed set must forget the text and the provenance too, or
    /// the panel would offer to save orbitals that are no longer shown.
    #[test]
    fn clearing_forgets_a_computed_sets_provenance_and_text() {
        let mut state = OrbitalState::default();
        state.adopt(
            behemoth_orbitals(),
            None,
            Some("run".to_string()),
            Some("text".to_string()),
        );
        state.clear();
        assert!(state.data.is_none());
        assert!(state.provenance.is_none());
        assert!(state.molden_text.is_none());
    }
}

#[cfg(test)]
mod real_file_tests {
    use super::basis::Shell;
    use super::evaluate::{accumulate_shell, shell_cutoff_radius, EvalScratch};
    use super::molden::{parse_molden, Spin};
    use std::path::Path;

    const EXAMPLE: &str = "examples/Exc_TPSSh-D4-def2-TZVP_SCMwater_16.molden";

    fn load() -> Option<super::molden::OrbitalData> {
        if !Path::new(EXAMPLE).exists() {
            return None;
        }
        Some(parse_molden(&std::fs::read_to_string(EXAMPLE).unwrap()).expect("example must parse"))
    }

    #[test]
    fn the_example_file_parses_as_open_shell_with_matching_sets() {
        let Some(data) = load() else {
            eprintln!("skipping: {EXAMPLE} not present");
            return;
        };
        assert_eq!(data.atoms.len(), 43);
        assert!(data.is_open_shell(), "file declares Spin= Beta blocks");
        assert_eq!(data.alpha.len(), 546);
        assert_eq!(data.beta.len(), 546);
        // The structural check: n_basis from [GTO] must equal the MO length.
        assert_eq!(data.layout.n_basis, 546);
        assert_eq!(data.alpha[0].coefficients.len(), 546);

        // Alpha and beta HOMOs are computed independently.
        assert!(data.homo_index(Spin::Alpha).is_some());
        assert!(data.homo_index(Spin::Beta).is_some());
    }

    #[test]
    fn the_example_uses_spherical_d_f_and_g_shells() {
        let Some(data) = load() else { return };
        let mut seen = [false; 5];
        for s in &data.shells {
            seen[s.l as usize] = true;
            assert!(
                s.spherical,
                "file declares [5D][7F][9G]; l={} should be pure",
                s.l
            );
        }
        assert!(seen[2] && seen[3], "def2-TZVP on Fe has d and f shells");
    }

    /// Integrate each basis function of a shell over a box scaled to that
    /// shell's own extent.  Every function must have unit norm.
    ///
    /// Only meaningful for shells whose primitives span a narrow exponent
    /// range: a uniform grid sized to the most diffuse primitive cannot
    /// resolve the cusp of a tight one.  Contracted shells are covered by
    /// [`radial_norm`] instead.
    fn shell_function_norms(shell: &Shell) -> Vec<f64> {
        let n = 48usize;
        let extent = shell_cutoff_radius(shell).min(12.0);
        let h = 2.0 * extent / n as f64;
        let n_func = shell.n_functions();
        let mut scratch = EvalScratch::default();
        let mut norms = vec![0.0; n_func];

        for slot in 0..n_func {
            let mut weights = vec![0.0; n_func];
            weights[slot] = 1.0;
            let mut acc = 0.0;
            for ix in 0..n {
                let x = -extent + (ix as f64 + 0.5) * h;
                for iy in 0..n {
                    let y = -extent + (iy as f64 + 0.5) * h;
                    for iz in 0..n {
                        let z = -extent + (iz as f64 + 0.5) * h;
                        let mut v = 0.0;
                        accumulate_shell(shell, [x, y, z], &weights, &mut scratch, &mut v);
                        acc += v * v;
                    }
                }
            }
            norms[slot] = acc * h * h * h;
        }
        norms
    }

    /// Norm of a shell's radial function by 1D quadrature.
    ///
    /// Separating the angular part makes this exact for arbitrarily tight
    /// contractions, which a uniform 3D grid cannot manage:
    ///
    ///   |phi|^2 = [4 pi (2l-1)!! / (2l+1)!!] * integral R(r)^2 r^(2l+2) dr
    ///
    /// for the `x^l` component, whose Cartesian factor is 1.
    fn radial_norm(shell: &Shell) -> f64 {
        let l = shell.l as usize;
        let cutoff = shell_cutoff_radius(shell).min(60.0);
        let steps = 400_000usize;
        let h = cutoff / steps as f64;

        let mut integral = 0.0;
        for k in 0..steps {
            let r = (k as f64 + 0.5) * h;
            let r2 = r * r;
            let mut radial = 0.0;
            for (&a, &c) in shell.exponents.iter().zip(&shell.coefficients) {
                radial += c * (-a * r2).exp();
            }
            integral += radial * radial * r.powi(2 * l as i32 + 2);
        }
        integral *= h;

        fn dfo(n: usize) -> f64 {
            let mut acc = 1.0;
            let mut k = 2 * n as isize - 1;
            while k > 1 {
                acc *= k as f64;
                k -= 2;
            }
            acc
        }
        // integral over the sphere of (x/r)^(2l) = 4 pi (2l-1)!! / (2l+1)!!
        let angular = 4.0 * std::f64::consts::PI * dfo(l) / dfo(l + 1);
        integral * angular
    }

    /// Every contracted shell in the real basis must have a normalised radial
    /// function.  This exercises the `[GTO]` parse and the contraction
    /// renormalisation across all 43 atoms, including the tight Fe core.
    #[test]
    fn every_shell_in_the_real_file_has_a_normalised_radial_function() {
        let Some(data) = load() else {
            eprintln!("skipping: {EXAMPLE} not present");
            return;
        };

        let mut worst = (0.0f64, 0usize, 0u8);
        for (index, shell) in data.shells.iter().enumerate() {
            let norm = radial_norm(shell);
            let error = (norm - 1.0).abs();
            if error > worst.0 {
                worst = (error, index, shell.l);
            }
        }
        assert!(
            worst.0 < 0.01,
            "shell {} (l={}) is off by {:.4} from unit norm",
            worst.1,
            worst.2,
            worst.0
        );
    }

    /// Single-primitive shells can be integrated on a uniform 3D grid, which
    /// also exercises the angular part: the pure d, f and g transforms.
    #[test]
    fn single_primitive_functions_from_the_real_file_are_normalised() {
        let Some(data) = load() else { return };

        let mut checked = [false; 5];
        for shell in &data.shells {
            let l = shell.l as usize;
            if checked[l] || shell.exponents.len() != 1 {
                continue;
            }
            checked[l] = true;

            for (slot, norm) in shell_function_norms(shell).iter().enumerate() {
                assert!(
                    (norm - 1.0).abs() < 0.02,
                    "l={l} component {slot} on atom {} has norm {norm}, expected 1",
                    shell.center
                );
            }
        }
        assert!(
            checked[0] && checked[1] && checked[2],
            "expected single-primitive s, p and d shells in def2-TZVP; got {checked:?}"
        );
    }

    /// End-to-end: parse the real file, sample its alpha HOMO onto a grid, and
    /// mesh both lobes.  This is the only test that exercises the whole chain
    /// on production data, including pure d/f/g shells and 546 coefficients.
    #[test]
    fn the_real_homo_samples_into_a_two_lobed_surface() {
        use super::grid::{plan_grid, sample_orbital};
        use super::isosurface;
        use bevy::tasks::{ComputeTaskPool, TaskPool};

        let Some(data) = load() else {
            eprintln!("skipping: {EXAMPLE} not present");
            return;
        };
        ComputeTaskPool::get_or_init(TaskPool::default);

        let spin = Spin::Alpha;
        let homo = data.homo_index(spin).expect("an occupied orbital exists");

        let started = std::time::Instant::now();
        let spec = plan_grid(&data.positions, 3.0, 0.25, 2_000_000);
        let field = sample_orbital(&data, spin, homo, &spec);
        let elapsed = started.elapsed();

        eprintln!(
            "sampled HOMO ({homo}) on {}x{}x{} = {} points in {:.2?}",
            spec.dims[0],
            spec.dims[1],
            spec.dims[2],
            spec.point_count(),
            elapsed
        );

        let peak = field.max_abs();
        assert!(
            (0.01..5.0).contains(&peak),
            "peak |psi| of {peak} is not physically plausible for a valence orbital"
        );

        // A molecular orbital has both signs; finding only one would mean the
        // spherical transform or a coefficient sign had been lost.
        let most_positive = field.values.iter().cloned().fold(f32::MIN, f32::max);
        let most_negative = field.values.iter().cloned().fold(f32::MAX, f32::min);
        assert!(
            most_positive > 0.0 && most_negative < 0.0,
            "orbital is single-signed: max {most_positive}, min {most_negative}"
        );

        let isovalue = peak * 0.25;
        let positive = isosurface::extract(&field, isovalue);
        let negative = isosurface::extract(&field, -isovalue);
        assert!(
            !positive.is_empty() && !negative.is_empty(),
            "expected both lobes at isovalue {isovalue}: {} and {} triangles",
            positive.triangle_count(),
            negative.triangle_count()
        );
        eprintln!(
            "lobes: {} + {} triangles at isovalue {isovalue:.4}",
            positive.triangle_count(),
            negative.triangle_count()
        );
    }

    /// The wireframe net must work on a real orbital, not just analytic
    /// spheres, and the spacing control must actually thin it out.
    #[test]
    fn the_real_homo_produces_a_wireframe_net() {
        use super::grid::sample_orbital;
        use super::isosurface;
        use bevy::tasks::{ComputeTaskPool, TaskPool};

        let Some(data) = load() else { return };
        ComputeTaskPool::get_or_init(TaskPool::default);

        let spin = Spin::Alpha;
        let homo = data.homo_index(spin).unwrap();
        let spec = super::plan_grid(&data.positions, 3.0, 0.25, 2_000_000);
        let field = sample_orbital(&data, spin, homo, &spec);
        let isovalue = field.max_abs() * 0.25;

        let dense = isosurface::contour_net(&field, isovalue, 0.10).len();
        let coarse = isosurface::contour_net(&field, isovalue, 0.80).len();
        eprintln!("net segments: 0.10 A spacing -> {dense}, 0.80 A -> {coarse}");

        assert!(dense > 0, "no wireframe produced for a real orbital");
        assert!(coarse < dense, "spacing did not thin the net");
    }

    /// Timing probe for the wireframe: the density slider rebuilds the net on
    /// every change, so extraction has to stay inside a frame budget.
    #[test]
    fn wireframe_extraction_stays_within_a_frame_budget() {
        use super::grid::sample_orbital;
        use super::isosurface;
        use bevy::tasks::{ComputeTaskPool, TaskPool};

        let Some(data) = load() else { return };
        ComputeTaskPool::get_or_init(TaskPool::default);

        let spin = Spin::Alpha;
        let homo = data.homo_index(spin).unwrap();
        // The grid the app actually uses at the default resolution.
        let spec = super::plan_grid(&data.positions, 3.0, 0.20, 2_000_000);
        let field = sample_orbital(&data, spin, homo, &spec);
        let isovalue = 0.10;

        for density in [
            super::MIN_MESH_DENSITY,
            super::DEFAULT_MESH_DENSITY,
            super::MAX_MESH_DENSITY,
        ] {
            let spacing = super::MESH_SPACING_BASE / density as f32;
            let started = std::time::Instant::now();
            let positive = isosurface::contour_net(&field, isovalue, spacing);
            let negative = isosurface::contour_net(&field, -isovalue, spacing);
            let elapsed = started.elapsed();
            let segments = positive.len() + negative.len();
            eprintln!("density {density} ({spacing:.3} A): {segments} segments in {elapsed:.1?}");

            assert!(segments > 0, "density {density} produced no net");
            // Generous enough not to flake on a slow machine, tight enough to
            // catch the 290 ms behaviour this path had before it was gathered
            // per plane and parallelised.
            assert!(
                elapsed.as_millis() < 100,
                "density {density} took {elapsed:.1?}; the wireframe rebuild runs \
                 synchronously in the surface rebuild and must fit a frame"
            );
        }
    }

    /// Total density and spin density on the real open-shell Fe complex.
    ///
    /// Checks the two properties that any correct result must have -- a
    /// non-negative total density, and a spin density that takes both signs on
    /// a high-spin system -- and reports the cost, since the estimate that
    /// justified building this at all was that the density-matrix route lands
    /// within a few seconds rather than the ~40 s an orbital-by-orbital sum
    /// would need.
    #[test]
    fn the_real_file_yields_a_density_and_a_spin_density() {
        use super::density::{sample_density, DensityMatrices, Quantity};
        use bevy::tasks::{ComputeTaskPool, TaskPool};

        let Some(data) = load() else {
            eprintln!("skipping: {EXAMPLE} not present");
            return;
        };
        ComputeTaskPool::get_or_init(TaskPool::default);

        let started = std::time::Instant::now();
        let matrices = DensityMatrices::build(&data);
        eprintln!("density matrices built in {:.2?}", started.elapsed());
        assert!(
            matrices.spin.is_some(),
            "open-shell file must give a spin density"
        );

        let spec = super::plan_grid(&data.positions, 3.0, 0.20, 2_000_000);
        eprintln!("grid: {} points", spec.point_count());

        let started = std::time::Instant::now();
        let total = sample_density(&data, &matrices, Quantity::Density, &spec);
        eprintln!("total density sampled in {:.2?}", started.elapsed());

        let started = std::time::Instant::now();
        let spin = sample_density(&data, &matrices, Quantity::SpinDensity, &spec);
        eprintln!("spin density sampled in {:.2?}", started.elapsed());

        // A total density is non-negative by construction.
        let minimum = total.values.iter().cloned().fold(f32::MAX, f32::min);
        assert!(
            minimum > -1.0e-4,
            "total density went negative ({minimum}), so the matrix or the \
             quadratic form is wrong"
        );
        assert!(total.max_abs() > 0.1, "total density is implausibly flat");

        // Alpha and beta differ on a high-spin metal centre, so the spin
        // density must be signed rather than uniformly positive.
        let most_positive = spin.values.iter().cloned().fold(f32::MIN, f32::max);
        let most_negative = spin.values.iter().cloned().fold(f32::MAX, f32::min);
        eprintln!("spin density range: {most_negative:.4} .. {most_positive:.4}");
        assert!(
            most_positive > 0.0 && most_negative < 0.0,
            "spin density is single-signed"
        );
    }
}
