//! Immutable export settings and per-frame geometry, independent of the editor.
use crate::{
    orbitals::{
        density::{DensityMatrices, Quantity},
        molden::{OrbitalData, Spin},
        OrbitalState, SurfaceStyle,
    },
    trajectory::TrajectoryFrame,
};
use bevy::prelude::*;
use std::sync::Arc;
#[derive(Clone, Default)]
pub struct Options {
    pub surface: bool,
    pub measurements: bool,
    pub labels: bool,
    pub caption: bool,
    pub halos: bool,
    pub reactivity: bool,
}
#[derive(Clone)]
pub struct Surface {
    quantity: Quantity,
    spin: Spin,
    index: usize,
    spacing: f32,
    isovalue: f32,
    style: SurfaceStyle,
    mesh_spacing: f32,
    opacity: f32,
    colors: [Color; 2],
}
#[derive(Clone)]
pub struct Reactivity {
    pub results: Arc<crate::reactivity::data::Results>,
    pub descriptor: crate::reactivity::data::Descriptor,
    pub color_scale: crate::property_legend::Scale,
    pub surface: bool,
    pub halos: bool,
    pub labels: bool,
    pub scale: f32,
    pub opacity: f32,
    pub iso: f32,
    pub spacing: f32,
    pub colors: [Color; 2],
}
#[derive(Clone, Default)]
pub struct Extras {
    pub reactivity: Option<Reactivity>,
    pub options: Options,
    pub measurements: Option<crate::measurements::Measurements>,
    pub labels: Option<crate::saved_views::Library>,
    surface: Option<Surface>,
    wavefunctions: Option<Arc<crate::trajectory_wavefunctions::job::Series>>,
    tracking: Option<Arc<crate::orbital_tracking::Track>>,
    static_data: Option<Arc<OrbitalData>>,
    pub polar: Option<Arc<crate::polarizability::Results>>,
    pub polar_scale: f32,
    pub polar_colors: crate::property_legend::Scale,
    pub polar_opacity: f32,
    pub polar_labels: bool,
    pub polar_angstrom: bool,
}
pub struct Geometry {
    pub meshes: Vec<(Mesh, StandardMaterial, Transform)>,
}
impl Extras {
    pub fn snapshot(
        options: Options,
        orbital: &OrbitalState,
        wave: &crate::trajectory_wavefunctions::State,
        polar: &crate::polarizability::State,
        reactivity: &crate::reactivity::State,
        measurements: &crate::measurements::Measurements,
        labels: &crate::saved_views::Library,
        frames: &[TrajectoryFrame],
        rotation: bool,
    ) -> Result<Self, String> {
        let surface = options.surface.then(|| Surface {
            quantity: orbital.quantity,
            spin: orbital.spin,
            index: orbital.selected.unwrap_or(0),
            spacing: orbital.grid_spacing(),
            isovalue: orbital.isovalue,
            style: orbital.style,
            mesh_spacing: orbital.mesh_plane_spacing(),
            opacity: orbital.opacity,
            colors: [orbital.positive_color, orbital.negative_color],
        });
        let react = if options.reactivity {
            let results = reactivity
                .results
                .clone()
                .ok_or("Compute reactivity before exporting it.")?;
            let mut settings = OrbitalState::default();
            settings.resolution = reactivity.resolution;
            Some(Reactivity {
                results,
                descriptor: reactivity.descriptor,
                surface: reactivity.surface,
                halos: reactivity.halos,
                labels: reactivity.labels,
                scale: reactivity.scale,
                opacity: reactivity.opacity,
                iso: reactivity.isovalue,
                spacing: settings.grid_spacing(),
                colors: [reactivity.positive, reactivity.negative],
                color_scale: reactivity.colors.clone(),
            })
        } else {
            None
        };
        let result = Self {
            reactivity: react,
            measurements: options.measurements.then(|| measurements.clone()),
            labels: options.labels.then(|| labels.clone()),
            surface,
            wavefunctions: wave.series.clone(),
            static_data: if rotation { orbital.data.clone() } else { None },
            polar: if options.halos {
                polar.results.clone()
            } else {
                None
            },
            tracking: if wave.tracking.enabled {wave.tracking.track.clone()} else {None},
            polar_scale: polar.scale,
            polar_colors: polar.colors.clone(),
            polar_opacity: polar.opacity,
            polar_labels: polar.labels,
            polar_angstrom: polar.angstrom,
            options,
        };
        for frame in frames {
            if let Some(r) = &result.reactivity {
                if !r.results.frames.iter().any(|f| {
                    crate::trajectory_wavefunctions::job::geometry_matches(&f.geometry, frame)
                }) {
                    return Err("Reactivity data are missing for an exported geometry.".into());
                }
            }
            if result.surface.is_some() {
                if rotation {
                    crate::trajectory_wavefunctions::job::check_geometry(
                        result.static_data.as_deref().ok_or(
                            "Load or compute a wavefunction before exporting its surface.",
                        )?,
                        frame,
                    )?;
                } else if !result.wavefunctions.as_ref().is_some_and(|s| {
                    s.manifest.frames.iter().any(|f| {
                        crate::trajectory_wavefunctions::job::geometry_matches(&f.geometry, frame)
                    })
                }) {
                    return Err("A matching wavefunction is needed for every exported geometry. Compute or load a wavefunction trajectory in Surface Tools first.".into());
                }
            }
            if result.options.halos
                && !result.polar.as_ref().is_some_and(|p| {
                    p.frames
                        .iter()
                        .any(|f| crate::trajectory_wavefunctions::job::geometry_matches(f, frame))
                })
            {
                return Err("Compute atomic polarizabilities for every exported geometry in Surface Tools first.".into());
            }
        }
        Ok(result)
    }
    pub fn values(&self, frame: &TrajectoryFrame) -> Option<&crate::polarizability::data::Values> {
        let r = self.polar.as_ref()?;
        let i = r
            .frames
            .iter()
            .position(|f| crate::trajectory_wavefunctions::job::geometry_matches(f, frame))?;
        r.values.get(i)
    }
    pub fn geometry(
        &self,
        frame: &TrajectoryFrame,
        hidden: &std::collections::BTreeSet<usize>,
    ) -> Result<Geometry, String> {
        use crate::orbitals::{density, grid, isosurface};
        use bevy::{asset::RenderAssetUsages, mesh::PrimitiveTopology};
        let mut meshes = vec![];
        if let Some(s) = &self.surface {
            let mut index=s.index;
            let data = if let Some(d) = &self.static_data {
                d.clone()
            } else {
                let series = self
                    .wavefunctions
                    .as_ref()
                    .ok_or("No trajectory wavefunctions")?;
                let i = series
                    .manifest
                    .frames
                    .iter()
                    .position(|f| {
                        crate::trajectory_wavefunctions::job::geometry_matches(&f.geometry, frame)
                    })
                    .ok_or("No matching wavefunction")?;
                let mut d=crate::trajectory_wavefunctions::job::load(series,i)?.0;
                if let Some(t)=&self.tracking {if s.quantity==Quantity::Orbital {
                    if t.root!=series.root || t.spin!=s.spin {return Err("Orbital tracking does not match the export spin or wavefunction series.".into());}
                    index=crate::orbital_tracking::apply(&mut d,t,i)?.index;
                }}
                Arc::new(d)
            };
            crate::trajectory_wavefunctions::job::check_geometry(&data, frame)?;
            if !s.quantity.is_available(&data) {
                return Err(
                    "Spin density requires separate alpha/beta orbitals in every exported frame."
                        .into(),
                );
            }
            let grid = grid::plan_grid(&data.positions, 3., s.spacing, 2_000_000);
            let field = if s.quantity == Quantity::Orbital {
                if index >= data.orbitals(s.spin).len() {
                    return Err("The selected orbital is missing from an exported frame.".into());
                }
                grid::sample_orbital(&data, s.spin, index, &grid)
            } else {
                density::sample_density(&data, &DensityMatrices::build(&data), s.quantity, &grid)
            };
            for (i, color) in s.colors.iter().enumerate() {
                if i == 1 && !s.quantity.is_signed() {
                    continue;
                }
                let iso = if i == 0 { s.isovalue } else { -s.isovalue };
                let mut mesh;
                if s.style == SurfaceStyle::Solid {
                    let surface = isosurface::extract(&field, iso);
                    if surface.is_empty() {
                        continue;
                    }
                    mesh = Mesh::new(
                        PrimitiveTopology::TriangleList,
                        RenderAssetUsages::default(),
                    );
                    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, surface.positions);
                    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, surface.normals);
                } else {
                    let lines = isosurface::contour_net(&field, iso, s.mesh_spacing);
                    if lines.is_empty() {
                        continue;
                    }
                    let pos: Vec<[f32; 3]> = lines
                        .into_iter()
                        .flat_map(|l| l.map(|p| p.to_array()))
                        .collect();
                    let normal = vec![[0., 1., 0.]; pos.len()];
                    mesh = Mesh::new(PrimitiveTopology::LineList, RenderAssetUsages::default());
                    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, pos);
                    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normal);
                }
                let c = color.to_srgba();
                meshes.push((
                    mesh,
                    StandardMaterial {
                        base_color: Color::srgba(c.red, c.green, c.blue, s.opacity),
                        alpha_mode: if s.opacity >= 0.999 {
                            AlphaMode::Opaque
                        } else {
                            AlphaMode::Blend
                        },
                        perceptual_roughness: 0.45,
                        cull_mode: None,
                        unlit: s.style == SurfaceStyle::Mesh,
                        ..default()
                    },
                    Transform::default(),
                ));
            }
        }
        if let Some(v) = self.values(frame) {
            for (i, &a) in v.alpha.iter().enumerate() {
                if hidden.contains(&i) {
                    continue;
                }
                let radius = crate::polarizability::data::radius(a, self.polar_scale);
                if radius <= 0. {
                    continue;
                }
                meshes.push((
                    Sphere::new(radius).mesh().uv(32, 20),
                    StandardMaterial {
                        base_color: self.polar_colors.color(a,Color::srgb(0.25,0.75,1.),Color::srgb(0.25,0.75,1.),self.polar_opacity),
                        alpha_mode: AlphaMode::Blend,
                        perceptual_roughness: 0.65,
                        ..default()
                    },
                    Transform::from_translation(frame.pos[i]),
                ));
            }
        }
        if let Some(r) = &self.reactivity {
            let f = r
                .results
                .frames
                .iter()
                .find(|f| {
                    crate::trajectory_wavefunctions::job::geometry_matches(&f.geometry, frame)
                })
                .ok_or("Missing reactivity frame")?;
            if r.surface && r.descriptor.local() {
                let field = crate::reactivity::data::field(f, r.descriptor, r.spacing)?;
                for (i, iso) in [r.iso, -r.iso].into_iter().enumerate() {
                    let surface = isosurface::extract(&field, iso);
                    if surface.is_empty() {
                        continue;
                    }
                    let mut mesh = Mesh::new(
                        PrimitiveTopology::TriangleList,
                        RenderAssetUsages::default(),
                    );
                    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, surface.positions);
                    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, surface.normals);
                    let c = r.colors[i].to_srgba();
                    meshes.push((
                        mesh,
                        StandardMaterial {
                            base_color: Color::srgba(c.red, c.green, c.blue, r.opacity),
                            alpha_mode: AlphaMode::Blend,
                            cull_mode: None,
                            ..default()
                        },
                        Transform::default(),
                    ));
                }
            }
            if r.halos {
                if let Some(values) = r.descriptor.atoms(f) {
                    for (i, value) in values.into_iter().enumerate() {
                        if hidden.contains(&i) || value.abs() < 1e-12 {
                            continue;
                        }
                        let c = r.color_scale.color(value,r.colors[0],r.colors[1],r.opacity).to_srgba();
                        meshes.push((
                            Sphere::new(r.scale * value.abs().cbrt() as f32)
                                .mesh()
                                .uv(32, 20),
                            StandardMaterial {
                                base_color: Color::srgba(c.red, c.green, c.blue, r.opacity),
                                alpha_mode: AlphaMode::Blend,
                                ..default()
                            },
                            Transform::from_translation(frame.pos[i]),
                        ));
                    }
                }
            }
        }
        Ok(Geometry { meshes })
    }
    pub fn needs_geometry(&self) -> bool {
        self.surface.is_some() || self.polar.is_some() || self.reactivity.is_some()
    }
}
