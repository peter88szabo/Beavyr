//! Off-screen molecular scenes with independent cameras and entity ownership.
//! They never replace the editor's molecule or use its atom/bond marker components.
use crate::{
    camera::OrbitCamera,
    molecule::Molecule,
    settings::{BondColorMode, MolSettings, RepresentationMode as R},
};
use bevy::render::render_resource::{TextureFormat, TextureUsages};
use bevy::{
    camera::{visibility::RenderLayers, RenderTarget},
    prelude::*,
};
use std::collections::{BTreeSet, HashMap};

pub struct Scene {
    pub image: Handle<Image>,
    pub camera: Entity,
    entities: Vec<Entity>,
    geometry: Vec<Entity>,
    sphere: Handle<Mesh>,
    cylinder: Handle<Mesh>,
    materials: HashMap<String, Handle<StandardMaterial>>,
    layer: usize,
}
impl Scene {
    pub fn new(
        commands: &mut Commands,
        images: &mut Assets<Image>,
        meshes: &mut Assets<Mesh>,
        settings: &MolSettings,
        orbit: &OrbitCamera,
        projection: &Projection,
        size: UVec2,
        layer: usize,
    ) -> Self {
        let mut image =
            Image::new_target_texture(size.x, size.y, TextureFormat::Bgra8UnormSrgb, None);
        image.texture_descriptor.usage |= TextureUsages::COPY_SRC;
        let image = images.add(image);
        let mut projection = projection.clone();
        match &mut projection {
            Projection::Perspective(p) => p.aspect_ratio = size.x as f32 / size.y as f32,
            Projection::Orthographic(p) => {
                let height = p.area.height();
                p.area = Rect::from_center_size(
                    Vec2::ZERO,
                    Vec2::new(height * size.x as f32 / size.y as f32, height),
                );
            }
            _ => {}
        }
        let camera = commands
            .spawn((
                Camera3d::default(),
                Camera {
                    clear_color: ClearColorConfig::Custom(settings.bg_color),
                    ..default()
                },
                projection,
                RenderTarget::Image(image.clone().into()),
                camera_transform(orbit),
                RenderLayers::layer(layer),
            ))
            .id();
        let mut entities = vec![camera];
        let lighting = [
            (
                Vec3::splat(settings.light_distance),
                settings.light_intensity,
            ),
            (
                Vec3::new(
                    -settings.fill_distance,
                    settings.fill_distance * 0.5,
                    settings.fill_distance,
                ),
                settings.fill_intensity,
            ),
            (
                Vec3::new(
                    settings.rim_distance,
                    settings.rim_distance * 0.5,
                    -settings.rim_distance,
                ),
                settings.rim_intensity,
            ),
        ];
        for (index, (offset, intensity)) in lighting.into_iter().enumerate() {
            if index > 0 && settings.lighting_mode == crate::settings::LightingMode::SinglePoint {
                continue;
            }
            entities.push(
                commands
                    .spawn((
                        PointLight {
                            intensity,
                            range: 200.,
                            shadow_maps_enabled: false,
                            ..default()
                        },
                        Transform::from_translation(orbit.target + orbit.orientation * offset),
                        RenderLayers::layer(layer),
                    ))
                    .id(),
            );
        }
        Self {
            image,
            camera,
            entities,
            geometry: vec![],
            sphere: meshes.add(Sphere::new(1.).mesh().uv(32, 20)),
            cylinder: meshes.add(Cylinder::new(1., 1.)),
            materials: HashMap::new(),
            layer,
        }
    }
    fn material(
        &mut self,
        symbol: &str,
        settings: &MolSettings,
        materials: &mut Assets<StandardMaterial>,
    ) -> Handle<StandardMaterial> {
        self.materials
            .entry(symbol.to_string())
            .or_insert_with(|| {
                let color = if symbol == "__bond" {
                    settings.uniform_bond_color
                } else {
                    settings
                        .element_colors
                        .get(symbol)
                        .copied()
                        .unwrap_or_else(|| crate::color_schemes::color_for(symbol, settings.scheme))
                };
                materials.add(StandardMaterial {
                    base_color: color,
                    metallic: settings.metallic,
                    perceptual_roughness: settings.roughness,
                    reflectance: settings.reflectance,
                    emissive: if settings.use_emissive {
                        settings.emissive_color.to_linear() * settings.emissive_strength
                    } else {
                        LinearRgba::BLACK
                    },
                    ..default()
                })
            })
            .clone()
    }
    fn piece(
        &mut self,
        commands: &mut Commands,
        mesh: Handle<Mesh>,
        material: Handle<StandardMaterial>,
        transform: Transform,
    ) {
        self.geometry.push(
            commands
                .spawn((
                    Mesh3d(mesh),
                    MeshMaterial3d(material),
                    transform,
                    RenderLayers::layer(self.layer),
                ))
                .id(),
        );
    }
    pub fn molecule(
        &mut self,
        commands: &mut Commands,
        materials: &mut Assets<StandardMaterial>,
        molecule: &Molecule,
        settings: &MolSettings,
        hidden: &BTreeSet<usize>,
    ) {
        self.molecule_with_colors(commands, materials, molecule, settings, hidden, None);
    }
    fn atom_material(&mut self, index: usize, symbol: &str, settings: &MolSettings, materials: &mut Assets<StandardMaterial>, colors: Option<&[Color]>) -> Handle<StandardMaterial> {
        let base = self.material(symbol, settings, materials);
        let Some(color) = colors.and_then(|c| c.get(index)).copied() else { return base; };
        let key = format!("displacement-{index}");
        if let Some(handle) = self.materials.get(&key) {
            if let Some(mut material) = materials.get_mut(handle) { material.base_color = color; }
            return handle.clone();
        }
        let mut material = materials.get(&base).cloned().unwrap_or_default();
        material.base_color = color;
        let handle = materials.add(material);
        self.materials.insert(key, handle.clone());
        handle
    }
    pub fn molecule_with_colors(
        &mut self, commands: &mut Commands, materials: &mut Assets<StandardMaterial>,
        molecule: &Molecule, settings: &MolSettings, hidden: &BTreeSet<usize>, colors: Option<&[Color]>,
    ) {
        for e in self.geometry.drain(..) {
            commands.entity(e).despawn();
        }
        let visible = |i: usize| {
            !hidden.contains(&i)
                && crate::scene::element_visible(&molecule.atoms[i], settings)
                && (settings.representation != R::BackboneTrace
                    || (molecule.atoms[i] != "H"
                        && molecule
                            .topology
                            .as_ref()
                            .is_none_or(|t| t.backbone_atom(i))))
        };
        let atoms = matches!(
            settings.representation,
            R::BallAndStick | R::SpaceFilling | R::LowResBallsAndLines
        ) || (settings.representation == R::BackboneTrace && settings.trace_show_atoms);
        if atoms {
            for (i, position) in molecule.pos.iter().enumerate() {
                if visible(i) {
                    let material = self.atom_material(i, &molecule.atoms[i], settings, materials, colors);
                    self.piece(
                        commands,
                        self.sphere.clone(),
                        material,
                        Transform::from_translation(*position).with_scale(Vec3::splat(
                            crate::scene::atom_display_radius(&molecule.atoms[i], settings),
                        )),
                    );
                }
            }
        }
        if settings.representation == R::SpaceFilling {
            return;
        }
        for &(i, j, _) in &molecule.bonds {
            if !visible(i)
                || !visible(j)
                || crate::molecule::is_dummy(&molecule.atoms[i])
                || crate::molecule::is_dummy(&molecule.atoms[j])
            {
                continue;
            }
            let a = molecule.pos[i];
            let b = molecule.pos[j];
            let length = a.distance(b);
            if length < 1e-5 {
                continue;
            }
            let lines = matches!(
                settings.representation,
                R::LinesOnly | R::LowResBallsAndLines | R::BackboneTrace
            );
            let radius = if lines {
                0.008 * settings.line_bond_thickness
            } else if settings.representation == R::SticksRounded {
                settings.stick_radius
            } else {
                (crate::molecule::covalent_radius_angstrom(&molecule.atoms[i]).min(
                    crate::molecule::covalent_radius_angstrom(&molecule.atoms[j]),
                ) * settings.bond_radius_pct)
                    .max(0.01)
            };
            let split = lines
                || settings.bond_color_mode == BondColorMode::AtomSplit
                || settings.representation == R::SticksRounded;
            for half in 0..if split { 2 } else { 1 } {
                let (start, end) = if !split {
                    (a, b)
                } else if half == 0 {
                    (a, (a + b) * 0.5)
                } else {
                    ((a + b) * 0.5, b)
                };
                let mat = self.material(
                    if !split {
                        "__bond"
                    } else if half == 0 {
                        &molecule.atoms[i]
                    } else {
                        &molecule.atoms[j]
                    },
                    settings,
                    materials,
                );
                self.piece(
                    commands,
                    self.cylinder.clone(),
                    mat.clone(),
                    Transform::from_translation((start + end) * 0.5)
                        .with_rotation(Quat::from_rotation_arc(Vec3::Y, (end - start).normalize()))
                        .with_scale(Vec3::new(radius, start.distance(end), radius)),
                );
                if settings.representation == R::SticksRounded {
                    self.piece(
                        commands,
                        self.sphere.clone(),
                        mat,
                        Transform::from_translation(if half == 0 { a } else { b })
                            .with_scale(Vec3::splat(radius)),
                    );
                }
            }
        }
    }
    /// Additional geometry owned by this private scene and cleared with each frame.
    pub fn overlay_mesh(&mut self, commands: &mut Commands, meshes: &mut Assets<Mesh>, materials: &mut Assets<StandardMaterial>, mesh: Mesh, material: StandardMaterial, transform: Transform) {
        self.geometry.push(commands.spawn((Mesh3d(meshes.add(mesh)), MeshMaterial3d(materials.add(material)), transform, RenderLayers::layer(self.layer))).id());
    }
    pub fn dispose(self, commands: &mut Commands, images: &mut Assets<Image>) {
        for entity in self.entities.into_iter().chain(self.geometry) {
            commands.entity(entity).despawn();
        }
        images.remove(self.image.id());
    }
}
pub fn camera_transform(orbit: &OrbitCamera) -> Transform {
    Transform::from_translation(orbit.target + orbit.orientation * Vec3::Z * orbit.radius)
        .with_rotation(orbit.orientation)
}
