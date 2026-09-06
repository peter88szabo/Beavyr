// src/scene.rs

use bevy::camera::visibility::RenderLayers;
use bevy::camera::Viewport;
use bevy::math::primitives::{Cone, Cylinder, Sphere};
use bevy::prelude::*;
use bevy_egui::PrimaryEguiContext;
use std::collections::HashMap;

use crate::camera::OrbitCamera;
use crate::color_schemes::color_for;
use crate::molecule::{covalent_radius_angstrom, Molecule};
use crate::settings::{BondColorMode, LightingMode, MolSettings, RepresentationMode};

// NEW: react to coordinate-change events (only to mark dirty when topology might have changed)
use crate::events::MoleculeChangeReason;
use crate::events::MoleculeChanged;

// default molecule (Å)
//pub const DEFAULT_WATER: &str = r#"
//O 0.00000    0.00000    0.12008
//H 0.00000    0.71604   -0.48061
//H 0.00000   -0.71604   -0.48061
//"#;

pub const DEFAULT_WATER: &str = r#"
 C  3.40  -0.28 -0.90
 C  2.32  -0.59  0.09
 C  0.97  -0.78 -0.39
 O  0.04  -1.30  0.19
 C  2.64  -0.82  1.40
 O  1.65  -0.78  2.33
 O  2.11  -1.34  3.55
 C  2.99  -3.12  1.30
 O  2.92  -3.49  2.59
 O  4.01  -2.93  3.31
 H  3.64  -0.67  1.79
 H  0.81  -0.43 -1.43
 H  3.17   0.62 -1.45
 H  4.36  -0.14 -0.41
 H  3.51  -1.08 -1.63
 H  3.98  -3.11  0.85
 H  2.15  -3.50  0.74
 H  3.52  -2.59  4.10
 H  0.37  -1.95  0.79
"#;

#[derive(Component)]
pub struct AtomMarker {
    pub index: usize,
}

#[derive(Clone, Copy)]
pub enum BondVisual {
    Uniform {
        base_len: f32,
    },
    Split {
        start_frac: f32,
        end_frac: f32,
        base_len: f32,
        color_atom: usize,
    },
    Cap {
        atom: usize,
    },
}

#[derive(Component)]
pub struct BondMarker {
    pub i: usize,
    pub j: usize,
    pub visual: BondVisual,
}
#[derive(Component)]
pub struct MainCamera;
#[derive(Component)]
pub struct AxisCamera;

#[derive(Component)]
pub struct KeyLight;
#[derive(Component)]
pub struct FillLight;
#[derive(Component)]
pub struct RimLight;

#[derive(Default, Reflect, GizmoConfigGroup)]
pub struct BondLineGizmos;

#[derive(Default, Resource)]
pub struct GeometryCache {
    sphere_meshes: HashMap<SphereMeshKey, Handle<Mesh>>,
    cylinder_meshes: HashMap<CylinderMeshKey, Handle<Mesh>>,
    materials: HashMap<MaterialKey, Handle<StandardMaterial>>,
}

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
struct SphereMeshKey {
    radius_bits: u32,
    resolution: u32,
}

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
struct CylinderMeshKey {
    radius_bits: u32,
    height_bits: u32,
}

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
struct MaterialKey {
    base_color: u32,
    metallic_milli: u16,
    roughness_milli: u16,
    reflectance_milli: u16,
    emissive_milli: [u16; 4],
}

pub const LAYER_MAIN: usize = 0;
pub const LAYER_AXES: usize = 1;

const AXIS_VP_SIZE: u32 = 140;
const AXIS_VP_MARGIN: u32 = 12;
/// Distance from the gizmo's origin to an arrow tip, in its own world units
/// (`axis_len + head_len` as spawned in `setup`).
const AXIS_TIP_LEN: f32 = 1.4 + 0.28;
/// The axis camera's fixed orbit radius (its spawn position's length);
/// `sync_axis_camera_to_main` rotates it but preserves this distance.
const AXIS_CAM_RADIUS: f32 = 4.330_127;
/// Bevy's default perspective field of view.
const AXIS_CAM_FOV: f32 = std::f32::consts::FRAC_PI_4;
/// How far past the arrow tip each letter sits, as a multiple of the tip's
/// own projected distance. Greater than one, so the letter clears the arrow
/// head instead of overlapping it.
const AXIS_LABEL_OVERSHOOT: f32 = 1.22;
/// Keep labels this many physical pixels clear of the window edge.
const AXIS_LABEL_PADDING: f32 = 12.0;

/// Screen positions, in logical egui points, where the orientation gizmo's
/// X/Y/Z arrow tips appear for the given camera orientation.
///
/// The gizmo camera mirrors the main camera's rotation, so an axis direction
/// projects onto the little viewport simply by its components along the
/// camera's right and up vectors -- no separate render-world query needed.
pub fn axis_gizmo_label_positions(
    window_physical: UVec2,
    scale_factor: f32,
    cam_xform: &GlobalTransform,
) -> [(&'static str, Vec2, Color); 3] {
    let half = AXIS_VP_SIZE as f32 * 0.5;
    let centre_physical = Vec2::new(
        window_physical.x as f32 - AXIS_VP_MARGIN as f32 - half,
        window_physical.y as f32 - AXIS_VP_MARGIN as f32 - half,
    );
    let centre = centre_physical / scale_factor.max(0.001);

    // How far a fully sideways axis reaches, as a fraction of the half
    // viewport, under this camera's perspective.
    let reach = (AXIS_TIP_LEN / AXIS_CAM_RADIUS) / (AXIS_CAM_FOV * 0.5).tan();
    // Placed *past* the arrow head rather than on it, so the letter reads
    // clearly instead of sitting on top of the arrow it names.
    let radius = reach * half * AXIS_LABEL_OVERSHOOT / scale_factor.max(0.001);

    let right = cam_xform.right().as_vec3();
    let up = cam_xform.up().as_vec3();
    let project = |dir: Vec3| -> Vec2 {
        // Screen y grows downward, hence the negated up component.
        let raw = centre + Vec2::new(dir.dot(right), -dir.dot(up)) * radius;
        // Overshooting can push a label off the window edge, since the gizmo
        // sits in the corner -- keep them just inside it. Clamped against the
        // physical window, like every other dimension of the gizmo, so the
        // letters do not drift when the display scale changes.
        let max = (window_physical.as_vec2() - Vec2::splat(AXIS_LABEL_PADDING))
            / scale_factor.max(0.001);
        let min = Vec2::splat(AXIS_LABEL_PADDING / scale_factor.max(0.001));
        raw.clamp(min, max)
    };

    [
        ("X", project(Vec3::X), Color::srgb(1.0, 0.1, 0.1)),
        ("Y", project(Vec3::Y), Color::srgb(0.1, 1.0, 0.1)),
        ("Z", project(Vec3::Z), Color::srgb(0.1, 0.4, 1.0)),
    ]
}
const LARGE_MOLECULE_ATOM_THRESHOLD: usize = 1000;

fn element_visible(sym: &str, settings: &MolSettings) -> bool {
    settings
        .element_visibility
        .get(sym)
        .copied()
        .unwrap_or(true)
}

/// Sphere tessellation level for the current representation.
fn atom_display_resolution(settings: &MolSettings) -> u32 {
    match settings.representation {
        RepresentationMode::BallAndStick => settings.atom_resolution,
        RepresentationMode::LowResBallsAndLines => settings.low_res_atom_resolution,
        RepresentationMode::BackboneTrace => settings.low_res_atom_resolution,
        RepresentationMode::SpaceFilling => settings.cpk_atom_resolution,
        _ => settings.atom_resolution,
    }
}

/// Extra per-representation multiplier applied on top of `atom_scale`.
fn representation_atom_scale(settings: &MolSettings) -> f32 {
    match settings.representation {
        RepresentationMode::LowResBallsAndLines => settings.low_res_atom_scale,
        RepresentationMode::BackboneTrace => settings.trace_atom_scale,
        RepresentationMode::SpaceFilling => settings.cpk_atom_scale,
        _ => 1.0,
    }
}

/// Radius of the drawn sphere for an atom.
fn atom_display_radius(sym: &str, settings: &MolSettings) -> f32 {
    covalent_radius_angstrom(sym) * settings.atom_scale * representation_atom_scale(settings)
}

/// Radius the bond geometry treats each atom's sphere as occupying.  Bonds only
/// exist in representations where `representation_atom_scale` is 1.0, so this
/// deliberately omits that factor.
fn bond_endpoint_sphere_radius(sym: &str, settings: &MolSettings) -> f32 {
    covalent_radius_angstrom(sym) * settings.atom_scale
}

/// Cylinder radius for the bond between atoms `i` and `j`.
fn bond_display_radius(mol: &Molecule, settings: &MolSettings, i: usize, j: usize) -> f32 {
    let ri = bond_endpoint_sphere_radius(&mol.atoms[i], settings);
    let rj = bond_endpoint_sphere_radius(&mol.atoms[j], settings);
    let base = ri.min(rj);
    if matches!(settings.representation, RepresentationMode::SticksRounded) {
        settings.stick_radius.clamp(0.01, 1.0)
    } else {
        (base * settings.bond_radius_pct).clamp(0.01, base)
    }
}

/// Fractions along the i->j axis covered by the visible bond cylinder.
///
/// In ball-and-stick the cylinder starts inside atom i's sphere and stops inside
/// atom j's, overlapping slightly so no seam shows.  Rounded sticks span the
/// full axis and rely on cap spheres instead.
fn bond_extent_fracs(
    mol: &Molecule,
    settings: &MolSettings,
    i: usize,
    j: usize,
    len: f32,
) -> (f32, f32) {
    if matches!(settings.representation, RepresentationMode::SticksRounded) || len <= 1.0e-5 {
        return (0.0, 1.0);
    }
    let radius = bond_display_radius(mol, settings, i, j);
    let overlap = radius * 0.40; // 40% of bond radius
    let ri = bond_endpoint_sphere_radius(&mol.atoms[i], settings);
    let rj = bond_endpoint_sphere_radius(&mol.atoms[j], settings);
    ((ri - overlap) / len, (len - (rj - overlap)) / len)
}

fn apply_large_molecule_representation(mol: &Molecule, settings: &mut MolSettings) {
    if mol.atoms.len() <= LARGE_MOLECULE_ATOM_THRESHOLD {
        return;
    }

    if matches!(
        settings.representation,
        RepresentationMode::BallAndStick
            | RepresentationMode::SticksRounded
            | RepresentationMode::SpaceFilling
    ) {
        settings.representation = RepresentationMode::LowResBallsAndLines;
    }
}

fn color_key(color: Color) -> u32 {
    let s = color.to_srgba();
    let r = (s.red.clamp(0.0, 1.0) * 255.0).round() as u32;
    let g = (s.green.clamp(0.0, 1.0) * 255.0).round() as u32;
    let b = (s.blue.clamp(0.0, 1.0) * 255.0).round() as u32;
    let a = (s.alpha.clamp(0.0, 1.0) * 255.0).round() as u32;
    r | (g << 8) | (b << 16) | (a << 24)
}

fn quantize_unit_milli(v: f32) -> u16 {
    (v.clamp(0.0, 1.0) * 1000.0).round() as u16
}

fn quantize_nonnegative_milli(v: f32) -> u16 {
    (v.max(0.0).min(65.535) * 1000.0).round() as u16
}

fn emissive_key_milli(emissive: LinearRgba) -> [u16; 4] {
    [
        quantize_nonnegative_milli(emissive.red),
        quantize_nonnegative_milli(emissive.green),
        quantize_nonnegative_milli(emissive.blue),
        quantize_unit_milli(emissive.alpha),
    ]
}

fn cached_sphere_mesh(
    cache: &mut GeometryCache,
    meshes: &mut Assets<Mesh>,
    radius: f32,
    resolution: u32,
) -> Handle<Mesh> {
    let key = SphereMeshKey {
        radius_bits: radius.to_bits(),
        resolution,
    };
    cache
        .sphere_meshes
        .entry(key)
        .or_insert_with(|| meshes.add(Sphere::new(radius).mesh().ico(resolution).unwrap()))
        .clone()
}

fn cached_cylinder_mesh(
    cache: &mut GeometryCache,
    meshes: &mut Assets<Mesh>,
    radius: f32,
    height: f32,
) -> Handle<Mesh> {
    let key = CylinderMeshKey {
        radius_bits: radius.to_bits(),
        height_bits: height.to_bits(),
    };
    cache
        .cylinder_meshes
        .entry(key)
        .or_insert_with(|| meshes.add(Mesh::from(Cylinder::new(radius, height))))
        .clone()
}

fn cached_material(
    cache: &mut GeometryCache,
    materials: &mut Assets<StandardMaterial>,
    base_color: Color,
    metallic: f32,
    roughness: f32,
    reflectance: f32,
    emissive: LinearRgba,
) -> Handle<StandardMaterial> {
    let key = MaterialKey {
        base_color: color_key(base_color),
        metallic_milli: quantize_unit_milli(metallic),
        roughness_milli: quantize_unit_milli(roughness),
        reflectance_milli: quantize_unit_milli(reflectance),
        emissive_milli: emissive_key_milli(emissive),
    };
    cache
        .materials
        .entry(key)
        .or_insert_with(|| {
            materials.add(StandardMaterial {
                base_color,
                metallic,
                perceptual_roughness: roughness,
                reflectance,
                emissive,
                ..default()
            })
        })
        .clone()
}

pub fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    settings: Res<MolSettings>,
) {
    commands.insert_resource(ClearColor(settings.bg_color));

    // Ambient (resource)
    commands.insert_resource(GlobalAmbientLight {
        brightness: settings.ambient_brightness,
        color: settings.ambient_color,
        affects_lightmapped_meshes: true,
    });

    // Main camera
    commands.spawn((
        Camera3d::default(),
        PrimaryEguiContext,
        Transform::from_xyz(0.0, 7.0, 14.0).looking_at(Vec3::new(0.0, 1.0, 0.0), Vec3::Y),
        MainCamera,
        RenderLayers::layer(LAYER_MAIN),
    ));

    // Key light (front-right-up)
    commands.spawn((
        PointLight {
            shadow_maps_enabled: false,
            intensity: settings.light_intensity,
            range: 200.0,
            ..default()
        },
        Transform::from_xyz(
            settings.light_distance,
            settings.light_distance,
            settings.light_distance,
        ),
        KeyLight,
        RenderLayers::layer(LAYER_MAIN),
    ));

    // Fill (front-left)
    commands.spawn((
        PointLight {
            shadow_maps_enabled: false,
            intensity: settings.fill_intensity,
            range: 200.0,
            ..default()
        },
        Transform::from_xyz(
            -settings.fill_distance,
            settings.fill_distance * 0.5,
            settings.fill_distance,
        ),
        FillLight,
        RenderLayers::layer(LAYER_MAIN),
    ));

    // Rim (back-right)
    commands.spawn((
        PointLight {
            shadow_maps_enabled: false,
            intensity: settings.rim_intensity,
            range: 200.0,
            ..default()
        },
        Transform::from_xyz(
            settings.rim_distance,
            settings.rim_distance * 0.5,
            -settings.rim_distance,
        ),
        RimLight,
        RenderLayers::layer(LAYER_MAIN),
    ));

    // ---- Mini XYZ axes (separate layer + camera) ----
    let axis_radius = 0.06;
    let axis_len = 1.4;
    let head_len = 0.28;

    let shaft = meshes.add(Mesh::from(Cylinder::new(axis_radius, axis_len)));
    let head = meshes.add(Mesh::from(Cone {
        radius: axis_radius * 2.3,
        height: head_len,
    }));

    let red = materials.add(StandardMaterial {
        base_color: Color::srgb(1.0, 0.1, 0.1),
        emissive: Color::srgb(1.0, 0.1, 0.1).into(),
        ..default()
    });
    let green = materials.add(StandardMaterial {
        base_color: Color::srgb(0.1, 1.0, 0.1),
        emissive: Color::srgb(0.1, 1.0, 0.1).into(),
        ..default()
    });
    let blue = materials.add(StandardMaterial {
        base_color: Color::srgb(0.1, 0.4, 1.0),
        emissive: Color::srgb(0.1, 0.4, 1.0).into(),
        ..default()
    });

    let axes_root = commands
        .spawn((
            Transform::default(),
            Visibility::Visible,
            RenderLayers::layer(LAYER_AXES),
        ))
        .id();

    let to_x = Quat::from_rotation_arc(Vec3::Y, Vec3::X);
    let to_z = Quat::from_rotation_arc(Vec3::Y, Vec3::Z);

    let x_shaft = commands
        .spawn((
            Mesh3d(shaft.clone()),
            MeshMaterial3d(red.clone()),
            Transform {
                rotation: to_x,
                translation: Vec3::X * (axis_len * 0.5),
                ..default()
            },
            RenderLayers::layer(LAYER_AXES),
        ))
        .id();
    commands.entity(x_shaft).insert(ChildOf(axes_root));
    let x_head = commands
        .spawn((
            Mesh3d(head.clone()),
            MeshMaterial3d(red),
            Transform {
                rotation: to_x,
                translation: Vec3::X * (axis_len + head_len * 0.5),
                ..default()
            },
            RenderLayers::layer(LAYER_AXES),
        ))
        .id();
    commands.entity(x_head).insert(ChildOf(axes_root));

    let y_shaft = commands
        .spawn((
            Mesh3d(shaft.clone()),
            MeshMaterial3d(green.clone()),
            Transform {
                translation: Vec3::Y * (axis_len * 0.5),
                ..default()
            },
            RenderLayers::layer(LAYER_AXES),
        ))
        .id();
    commands.entity(y_shaft).insert(ChildOf(axes_root));
    let y_head = commands
        .spawn((
            Mesh3d(head.clone()),
            MeshMaterial3d(green),
            Transform {
                translation: Vec3::Y * (axis_len + head_len * 0.5),
                ..default()
            },
            RenderLayers::layer(LAYER_AXES),
        ))
        .id();
    commands.entity(y_head).insert(ChildOf(axes_root));

    let z_shaft = commands
        .spawn((
            Mesh3d(shaft),
            MeshMaterial3d(blue.clone()),
            Transform {
                rotation: to_z,
                translation: Vec3::Z * (axis_len * 0.5),
                ..default()
            },
            RenderLayers::layer(LAYER_AXES),
        ))
        .id();
    commands.entity(z_shaft).insert(ChildOf(axes_root));
    let z_head = commands
        .spawn((
            Mesh3d(head),
            MeshMaterial3d(blue),
            Transform {
                rotation: to_z,
                translation: Vec3::Z * (axis_len + head_len * 0.5),
                ..default()
            },
            RenderLayers::layer(LAYER_AXES),
        ))
        .id();
    commands.entity(z_head).insert(ChildOf(axes_root));

    // Axis camera (overlay)
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(2.5, 2.5, 2.5).looking_at(Vec3::ZERO, Vec3::Y),
        AxisCamera,
        RenderLayers::layer(LAYER_AXES),
        Camera {
            order: 1,
            clear_color: ClearColorConfig::None,
            viewport: None,
            ..default()
        },
    ));
}

/// Center the orbit camera on the initial molecule when the app starts.
pub fn center_camera_on_startup(mol: Res<Molecule>, mut cam: ResMut<OrbitCamera>) {
    auto_fit_camera(&mol.pos, &mut cam);
}

/// React to molecule coordinate changes (event).  
/// If the change may affect bond topology/meshes, mark the scene dirty so
/// `rebuild_if_dirty` will rebuild geometry next frame.
pub fn react_to_molecule_changed_mark_dirty(
    mut evr: MessageReader<MoleculeChanged>,
    mut settings: ResMut<MolSettings>,
    mol: Res<Molecule>,
    mut cam: ResMut<OrbitCamera>,
) {
    for ev in evr.read() {
        match ev.reason {
            MoleculeChangeReason::ParseXyz { recenter } => {
                if recenter {
                    auto_fit_camera(&mol.pos, &mut cam);
                }
                // Suggest a cheap representation once, when the molecule is loaded.
                // Doing this inside `rebuild_if_dirty` used to silently undo the
                // user's choice on every subsequent rebuild.
                apply_large_molecule_representation(&mol, &mut settings);
                settings.geometry_dirty = true;
                settings.bond_topology_dirty = true;
            }
            MoleculeChangeReason::SetPos
            | MoleculeChangeReason::BuilderRotate
            | MoleculeChangeReason::Undo
            | MoleculeChangeReason::Redo => {
                settings.coords_dirty = true;
            }
        }
    }
}

pub fn sync_coordinates_if_dirty(
    mut settings: ResMut<MolSettings>,
    mol: Res<Molecule>,
    mut q_atoms: Query<(&AtomMarker, &mut Transform)>,
    mut q_bonds: Query<(&BondMarker, &mut Transform), Without<AtomMarker>>,
) {
    if !settings.coords_dirty || settings.geometry_dirty || settings.bond_topology_dirty {
        return;
    }

    for (marker, mut tf) in q_atoms.iter_mut() {
        let Some(&pos) = mol.pos.get(marker.index) else {
            settings.geometry_dirty = true;
            settings.coords_dirty = false;
            return;
        };
        tf.translation = pos;
    }

    for (marker, mut tf) in q_bonds.iter_mut() {
        let (Some(&p0), Some(&p1)) = (mol.pos.get(marker.i), mol.pos.get(marker.j)) else {
            settings.geometry_dirty = true;
            settings.coords_dirty = false;
            return;
        };
        if !update_bond_transform(&mut tf, p0, p1, marker.visual) {
            settings.geometry_dirty = true;
            settings.coords_dirty = false;
            return;
        }
    }

    settings.coords_dirty = false;
}

fn update_bond_transform(tf: &mut Transform, p0: Vec3, p1: Vec3, visual: BondVisual) -> bool {
    let dir = p1 - p0;
    let len = dir.length();
    if len <= 1.0e-5 {
        return false;
    }

    let dir_n = dir / len;
    let rot = Quat::from_rotation_arc(Vec3::Y, dir_n);

    match visual {
        BondVisual::Uniform { base_len } => {
            if base_len <= 1.0e-5 {
                return false;
            }
            tf.translation = (p0 + p1) * 0.5;
            tf.rotation = rot;
            tf.scale = Vec3::new(1.0, len / base_len, 1.0);
        }
        BondVisual::Split {
            start_frac,
            end_frac,
            base_len,
            ..
        } => {
            if base_len <= 1.0e-5 {
                return false;
            }
            let segment_len = len * (end_frac - start_frac).max(0.0);
            if segment_len <= 1.0e-5 {
                return false;
            }
            tf.translation = p0 + dir_n * (len * (start_frac + end_frac) * 0.5);
            tf.rotation = rot;
            tf.scale = Vec3::new(1.0, segment_len / base_len, 1.0);
        }
        BondVisual::Cap { atom } => {
            tf.translation = if atom == 0 { p0 } else { p1 };
            tf.rotation = Quat::IDENTITY;
            tf.scale = Vec3::ONE;
        }
    }

    true
}

fn material_params(settings: &MolSettings) -> (f32, f32, f32, LinearRgba) {
    let metallic = settings.metallic.clamp(0.0, 1.0);
    let rough = settings.roughness.clamp(0.02, 1.0);
    let refl = settings.reflectance.clamp(0.0, 1.0);
    let emissive = if settings.use_emissive {
        let lr = settings.emissive_color.to_linear();
        LinearRgba::new(
            lr.red * settings.emissive_strength.max(0.0),
            lr.green * settings.emissive_strength.max(0.0),
            lr.blue * settings.emissive_strength.max(0.0),
            1.0,
        )
    } else {
        LinearRgba::BLACK
    };
    (metallic, rough, refl, emissive)
}

pub fn update_lighting_if_dirty(
    mut settings: ResMut<MolSettings>,
    mut clear: ResMut<ClearColor>,
    mut q_key: Query<&mut PointLight, (With<KeyLight>, Without<FillLight>, Without<RimLight>)>,
    mut q_fill: Query<&mut PointLight, (With<FillLight>, Without<KeyLight>, Without<RimLight>)>,
    mut q_rim: Query<&mut PointLight, (With<RimLight>, Without<KeyLight>, Without<FillLight>)>,
    mut amb: ResMut<GlobalAmbientLight>,
) {
    if !settings.lighting_dirty {
        return;
    }
    settings.lighting_dirty = false;

    *clear = ClearColor(settings.bg_color);
    amb.color = settings.ambient_color;
    amb.brightness = if settings.use_ambient {
        settings.ambient_brightness
    } else {
        0.0
    };

    if let Ok(mut l) = q_key.single_mut() {
        l.intensity = settings.light_intensity;
    }
    if let Ok(mut l) = q_fill.single_mut() {
        l.intensity = if matches!(settings.lighting_mode, LightingMode::ThreePoint) {
            settings.fill_intensity
        } else {
            0.0
        };
    }
    if let Ok(mut l) = q_rim.single_mut() {
        l.intensity = if matches!(settings.lighting_mode, LightingMode::ThreePoint) {
            settings.rim_intensity
        } else {
            0.0
        };
    }
}

pub fn update_materials_if_dirty(
    mut settings: ResMut<MolSettings>,
    mut cache: ResMut<GeometryCache>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mol: Res<Molecule>,
    mut q_atoms: Query<(&AtomMarker, &mut MeshMaterial3d<StandardMaterial>)>,
    mut q_bonds: Query<(&BondMarker, &mut MeshMaterial3d<StandardMaterial>), Without<AtomMarker>>,
) {
    if !settings.materials_dirty || settings.geometry_dirty || settings.bond_topology_dirty {
        return;
    }
    settings.materials_dirty = false;

    let (metallic, rough, refl, emissive) = material_params(&settings);

    for (marker, mut mat) in q_atoms.iter_mut() {
        let Some(sym) = mol.atoms.get(marker.index) else {
            settings.geometry_dirty = true;
            return;
        };
        let color = settings
            .element_colors
            .get(sym)
            .copied()
            .unwrap_or_else(|| color_for(sym, settings.scheme));
        *mat = MeshMaterial3d(cached_material(
            &mut cache,
            &mut materials,
            color,
            metallic,
            rough,
            refl,
            emissive,
        ));
    }

    for (marker, mut mat) in q_bonds.iter_mut() {
        let color = match marker.visual {
            BondVisual::Uniform { .. } => settings.uniform_bond_color,
            BondVisual::Split { color_atom, .. } => {
                let Some(sym) = mol.atoms.get(color_atom) else {
                    settings.geometry_dirty = true;
                    return;
                };
                settings
                    .element_colors
                    .get(sym)
                    .copied()
                    .unwrap_or_else(|| color_for(sym, settings.scheme))
            }
            BondVisual::Cap { atom } => {
                let idx = if atom == 0 { marker.i } else { marker.j };
                let Some(sym) = mol.atoms.get(idx) else {
                    settings.geometry_dirty = true;
                    return;
                };
                settings
                    .element_colors
                    .get(sym)
                    .copied()
                    .unwrap_or_else(|| color_for(sym, settings.scheme))
            }
        };
        *mat = MeshMaterial3d(cached_material(
            &mut cache,
            &mut materials,
            color,
            metallic,
            rough,
            refl,
            emissive,
        ));
    }
}

/// Re-point `Mesh3d` handles when only sizes/tessellation changed.
///
/// Atom scale, atom resolution, bond radius and stick thickness all leave the
/// entity set untouched — every atom and bond still exists, it just wants a
/// different cached mesh.  Handling that here instead of via `geometry_dirty`
/// is what keeps those sliders from despawning and respawning the whole scene
/// on every frame of a drag.
pub fn update_meshes_if_dirty(
    mut settings: ResMut<MolSettings>,
    mut cache: ResMut<GeometryCache>,
    mut meshes: ResMut<Assets<Mesh>>,
    mol: Res<Molecule>,
    mut q_atoms: Query<(&AtomMarker, &mut Mesh3d)>,
    mut q_bonds: Query<(&mut BondMarker, &mut Mesh3d, &mut Transform), Without<AtomMarker>>,
) {
    if !settings.meshes_dirty || settings.geometry_dirty || settings.bond_topology_dirty {
        return;
    }
    settings.meshes_dirty = false;

    let atom_resolution = atom_display_resolution(&settings);

    for (marker, mut mesh) in q_atoms.iter_mut() {
        let Some(sym) = mol.atoms.get(marker.index) else {
            settings.geometry_dirty = true;
            return;
        };
        let radius = atom_display_radius(sym, &settings);
        *mesh = Mesh3d(cached_sphere_mesh(
            &mut cache,
            &mut meshes,
            radius,
            atom_resolution,
        ));
    }

    let n_atoms = mol.atoms.len();
    for (mut marker, mut mesh, mut tf) in q_bonds.iter_mut() {
        let (Some(&p0), Some(&p1)) = (mol.pos.get(marker.i), mol.pos.get(marker.j)) else {
            settings.geometry_dirty = true;
            return;
        };
        if marker.i >= n_atoms || marker.j >= n_atoms {
            settings.geometry_dirty = true;
            return;
        }

        let radius = bond_display_radius(&mol, &settings, marker.i, marker.j);
        let len = (p1 - p0).length();

        match marker.visual {
            BondVisual::Uniform { .. } => {
                *mesh = Mesh3d(cached_cylinder_mesh(&mut cache, &mut meshes, radius, 1.0));
            }
            BondVisual::Split { color_atom, .. } => {
                let (start, end) = bond_extent_fracs(&mol, &settings, marker.i, marker.j, len);
                let mid = (start + end) * 0.5;
                // The entity for atom i owns the first half of the span, j the second.
                let (start_frac, end_frac) = if color_atom == marker.i {
                    (start, mid)
                } else {
                    (mid, end)
                };
                marker.visual = BondVisual::Split {
                    start_frac,
                    end_frac,
                    base_len: 1.0,
                    color_atom,
                };
                *mesh = Mesh3d(cached_cylinder_mesh(&mut cache, &mut meshes, radius, 1.0));
            }
            BondVisual::Cap { .. } => {
                *mesh = Mesh3d(cached_sphere_mesh(&mut cache, &mut meshes, radius, 2));
            }
        }

        if !update_bond_transform(&mut tf, p0, p1, marker.visual) {
            // The spheres grew until the visible cylinder vanished.  Shrink it away
            // rather than forcing a full respawn on every frame of the drag.
            tf.translation = (p0 + p1) * 0.5;
            tf.scale = Vec3::splat(1.0e-4);
        }
    }
}

pub fn rebuild_if_dirty(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut cache: ResMut<GeometryCache>,
    mut mol: ResMut<Molecule>,
    mut settings: ResMut<MolSettings>,
    q_atoms: Query<Entity, With<AtomMarker>>,
    q_bonds: Query<Entity, With<BondMarker>>,
) {
    if !settings.geometry_dirty && !settings.bond_topology_dirty {
        return;
    }
    let topology_dirty = settings.bond_topology_dirty;
    settings.geometry_dirty = false;
    settings.bond_topology_dirty = false;
    settings.materials_dirty = false;
    settings.meshes_dirty = false;
    settings.coords_dirty = false;

    // Despawn old geometry
    for e in q_atoms.iter() {
        commands.entity(e).despawn();
    }
    for e in q_bonds.iter() {
        commands.entity(e).despawn();
    }

    // Recomputing connectivity is the expensive half of a rebuild, so only do it
    // when something that actually affects topology changed.  A representation
    // switch or an element-visibility toggle respawns entities over the bond
    // list we already have.
    if topology_dirty {
        mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
    }

    let (metallic, rough, refl, emissive) = material_params(&settings);

    // Atoms
    let draw_atoms = matches!(
        settings.representation,
        RepresentationMode::BallAndStick
            | RepresentationMode::LowResBallsAndLines
            | RepresentationMode::SpaceFilling
    ) || (matches!(settings.representation, RepresentationMode::BackboneTrace)
        && settings.trace_show_atoms);
    if draw_atoms {
        let atom_resolution = atom_display_resolution(&settings);

        for i in 0..mol.atoms.len() {
            let sym = &mol.atoms[i];
            if matches!(settings.representation, RepresentationMode::BackboneTrace) && sym == "H" {
                continue;
            }
            if !element_visible(sym, &settings) {
                continue;
            }
            let r_cov = atom_display_radius(sym, &settings);

            let color = settings
                .element_colors
                .get(sym)
                .copied()
                .unwrap_or_else(|| color_for(sym, settings.scheme));

            let sphere_mesh = cached_sphere_mesh(&mut cache, &mut meshes, r_cov, atom_resolution);
            let mat = cached_material(
                &mut cache,
                &mut materials,
                color,
                metallic,
                rough,
                refl,
                emissive,
            );

            let pos = mol.pos[i];
            commands.spawn((
                Mesh3d(sphere_mesh),
                MeshMaterial3d(mat),
                Transform::from_translation(pos),
                AtomMarker { index: i },
                RenderLayers::layer(LAYER_MAIN),
            ));
        }
    }

    if matches!(
        settings.representation,
        RepresentationMode::LowResBallsAndLines
            | RepresentationMode::LinesOnly
            | RepresentationMode::BackboneTrace
            | RepresentationMode::SpaceFilling
    ) {
        return;
    }

    // Bonds
    let bonds = mol.bonds.clone();
    let mut degree: Vec<usize> = Vec::new();
    if matches!(settings.representation, RepresentationMode::SticksRounded) {
        degree = vec![0; mol.atoms.len()];
        for &(i, j, _d) in &bonds {
            if i < degree.len() {
                degree[i] += 1;
            }
            if j < degree.len() {
                degree[j] += 1;
            }
        }
    }
    for (i, j, len_cc) in bonds {
        if i >= mol.atoms.len() || j >= mol.atoms.len() {
            continue;
        }
        if !element_visible(&mol.atoms[i], &settings) || !element_visible(&mol.atoms[j], &settings)
        {
            continue;
        }

        let p0 = mol.pos[i];
        let p1 = mol.pos[j];
        if len_cc <= 0.0001 {
            continue;
        }

        let dir = p1 - p0;
        let dir_n = dir.normalize();

        // Bond cylinder/capsule radius based on smaller atom radius
        let radius = bond_display_radius(&mol, &settings, i, j);

        // Common rotation (align local +Y with bond direction)
        let rot = Quat::from_rotation_arc(Vec3::Y, dir_n);

        let bond_color_mode =
            if matches!(settings.representation, RepresentationMode::SticksRounded) {
                BondColorMode::AtomSplit
            } else {
                settings.bond_color_mode
            };

        match bond_color_mode {
            BondColorMode::Uniform => {
                // Single cached unit cylinder scaled to the current bond length.
                let center = (p0 + p1) * 0.5;
                let bond_mesh = cached_cylinder_mesh(&mut cache, &mut meshes, radius, 1.0);

                let mat = cached_material(
                    &mut cache,
                    &mut materials,
                    settings.uniform_bond_color,
                    metallic,
                    rough,
                    refl,
                    emissive,
                );

                commands.spawn((
                    Mesh3d(bond_mesh),
                    MeshMaterial3d(mat),
                    Transform {
                        translation: center,
                        rotation: rot,
                        scale: Vec3::new(1.0, len_cc, 1.0),
                        ..default()
                    },
                    BondMarker {
                        i,
                        j,
                        visual: BondVisual::Uniform { base_len: 1.0 },
                    },
                    RenderLayers::layer(LAYER_MAIN),
                ));
            }

            // AtomSplit: two segments with a clean joint at the midpoint,
            // optionally extended to overlap for rounded-stick mode.
            BondColorMode::AtomSplit => {
                let (start_frac, end_frac) = bond_extent_fracs(&mol, &settings, i, j, len_cc);
                let start = p0 + dir_n * (start_frac * len_cc);
                let end = p0 + dir_n * (end_frac * len_cc);

                let visible_len = end_frac - start_frac;
                if visible_len <= 0.0 {
                    continue;
                }
                let visible_len = visible_len * len_cc;

                let half_len = visible_len * 0.5;

                let mid = (start + end) * 0.5;
                let c0 = mid - dir_n * (visible_len * 0.25);
                let c1 = mid + dir_n * (visible_len * 0.25);
                let cyl_half = cached_cylinder_mesh(&mut cache, &mut meshes, radius, 1.0);
                let (mesh_i, mesh_j) = (cyl_half.clone(), cyl_half);

                let col_i = settings
                    .element_colors
                    .get(&mol.atoms[i])
                    .copied()
                    .unwrap_or_else(|| color_for(&mol.atoms[i], settings.scheme));
                let col_j = settings
                    .element_colors
                    .get(&mol.atoms[j])
                    .copied()
                    .unwrap_or_else(|| color_for(&mol.atoms[j], settings.scheme));

                let mat_i = cached_material(
                    &mut cache,
                    &mut materials,
                    col_i,
                    metallic,
                    rough,
                    refl,
                    emissive,
                );
                let mat_j = cached_material(
                    &mut cache,
                    &mut materials,
                    col_j,
                    metallic,
                    rough,
                    refl,
                    emissive,
                );

                commands.spawn((
                    Mesh3d(mesh_i),
                    MeshMaterial3d(mat_i.clone()),
                    Transform {
                        translation: c0,
                        rotation: rot,
                        scale: Vec3::new(1.0, half_len, 1.0),
                        ..default()
                    },
                    BondMarker {
                        i,
                        j,
                        visual: BondVisual::Split {
                            start_frac,
                            end_frac: (start_frac + end_frac) * 0.5,
                            base_len: 1.0,
                            color_atom: i,
                        },
                    },
                    RenderLayers::layer(LAYER_MAIN),
                ));

                commands.spawn((
                    Mesh3d(mesh_j),
                    MeshMaterial3d(mat_j.clone()),
                    Transform {
                        translation: c1,
                        rotation: rot,
                        scale: Vec3::new(1.0, half_len, 1.0),
                        ..default()
                    },
                    BondMarker {
                        i,
                        j,
                        visual: BondVisual::Split {
                            start_frac: (start_frac + end_frac) * 0.5,
                            end_frac,
                            base_len: 1.0,
                            color_atom: j,
                        },
                    },
                    RenderLayers::layer(LAYER_MAIN),
                ));

                if matches!(settings.representation, RepresentationMode::SticksRounded) {
                    let cap_mesh = cached_sphere_mesh(&mut cache, &mut meshes, radius, 2);
                    if i < degree.len() {
                        commands.spawn((
                            Mesh3d(cap_mesh.clone()),
                            MeshMaterial3d(mat_i.clone()),
                            Transform::from_translation(p0),
                            BondMarker {
                                i,
                                j,
                                visual: BondVisual::Cap { atom: 0 },
                            },
                            RenderLayers::layer(LAYER_MAIN),
                        ));
                    }
                    if j < degree.len() {
                        commands.spawn((
                            Mesh3d(cap_mesh.clone()),
                            MeshMaterial3d(mat_j.clone()),
                            Transform::from_translation(p1),
                            BondMarker {
                                i,
                                j,
                                visual: BondVisual::Cap { atom: 1 },
                            },
                            RenderLayers::layer(LAYER_MAIN),
                        ));
                    }
                }
            }
        }
    } // end bonds
}

/// Draw bonds as simple line segments (used in line-based representations).
pub fn draw_bond_lines(
    mut gizmos: Gizmos<BondLineGizmos>,
    mol: Res<Molecule>,
    settings: Res<MolSettings>,
) {
    if !matches!(
        settings.representation,
        RepresentationMode::LowResBallsAndLines
            | RepresentationMode::LinesOnly
            | RepresentationMode::BackboneTrace
    ) {
        return;
    }

    let n_pos = mol.pos.len();
    let n_atoms = mol.atoms.len();
    if n_pos == 0 || n_atoms == 0 {
        return;
    }

    for &(i, j, _len_cc) in &mol.bonds {
        if i >= n_pos || j >= n_pos || i >= n_atoms || j >= n_atoms {
            continue;
        }
        if !element_visible(&mol.atoms[i], &settings) || !element_visible(&mol.atoms[j], &settings)
        {
            continue;
        }
        if matches!(settings.representation, RepresentationMode::BackboneTrace) {
            if mol.atoms[i] == "H" || mol.atoms[j] == "H" {
                continue;
            }
        }

        let p0 = mol.pos[i];
        let p1 = mol.pos[j];
        let d = p1 - p0;
        if d.length() <= 1e-4 {
            continue;
        }

        let a = p0;
        let b = p1;
        let mid = (a + b) * 0.5;
        let col_i = settings
            .element_colors
            .get(&mol.atoms[i])
            .copied()
            .unwrap_or_else(|| color_for(&mol.atoms[i], settings.scheme));
        let col_j = settings
            .element_colors
            .get(&mol.atoms[j])
            .copied()
            .unwrap_or_else(|| color_for(&mol.atoms[j], settings.scheme));
        gizmos.line(a, mid, col_i);
        gizmos.line(mid, b, col_j);
    }
}

pub fn configure_bond_line_gizmos(
    mut cfg_store: ResMut<GizmoConfigStore>,
    settings: Res<MolSettings>,
) {
    let (cfg, _) = cfg_store.config_mut::<BondLineGizmos>();
    cfg.enabled = matches!(
        settings.representation,
        RepresentationMode::LowResBallsAndLines
            | RepresentationMode::LinesOnly
            | RepresentationMode::BackboneTrace
    );
    cfg.line.width = settings.line_bond_thickness.max(1.0).min(24.0);
    cfg.line.perspective = true;
}

/// Keep the light rig centered on the current camera target (molecule center).
pub fn update_light_positions(
    settings: Res<MolSettings>,
    cam: Res<crate::camera::OrbitCamera>,
    q_cam: Query<&GlobalTransform, (With<MainCamera>, Without<AxisCamera>)>,
    mut q_key: Query<&mut Transform, (With<KeyLight>, Without<FillLight>, Without<RimLight>)>,
    mut q_fill: Query<&mut Transform, (With<FillLight>, Without<KeyLight>, Without<RimLight>)>,
    mut q_rim: Query<&mut Transform, (With<RimLight>, Without<KeyLight>, Without<FillLight>)>,
) {
    let Ok(cam_tf) = q_cam.single() else {
        return;
    };
    let rot = cam_tf.rotation();
    let t = cam.target;

    if let Ok(mut tf) = q_key.single_mut() {
        tf.translation = t + rot
            * Vec3::new(
                settings.light_distance,
                settings.light_distance,
                settings.light_distance,
            );
    }
    if let Ok(mut tf) = q_fill.single_mut() {
        tf.translation = t + rot
            * Vec3::new(
                -settings.fill_distance,
                settings.fill_distance * 0.5,
                settings.fill_distance,
            );
    }
    if let Ok(mut tf) = q_rim.single_mut() {
        tf.translation = t + rot
            * Vec3::new(
                settings.rim_distance,
                settings.rim_distance * 0.5,
                -settings.rim_distance,
            );
    }
}

// ---------------------- Axis viewport & sync ----------------------

pub fn update_axis_viewport_on_resize(
    windows: Query<&Window>,
    mut q_axis_cam: Query<&mut Camera, With<AxisCamera>>,
) {
    let Ok(window) = windows.single() else {
        return;
    };
    let Ok(mut cam) = q_axis_cam.single_mut() else {
        return;
    };

    let w = window.physical_width();
    let h = window.physical_height();

    if w == 0 || h == 0 {
        cam.viewport = None;
        return;
    }

    let vp_w = AXIS_VP_SIZE;
    let vp_h = AXIS_VP_SIZE;
    // Bottom-right: the left edge belongs to the tool rail in the
    // floating-window layout, so the orientation gizmo would sit under it.
    let x = w.saturating_sub(vp_w + AXIS_VP_MARGIN);
    let y = h.saturating_sub(vp_h + AXIS_VP_MARGIN);

    let wanted = Viewport {
        physical_position: UVec2::new(x, y),
        physical_size: UVec2::new(vp_w, vp_h),
        ..default()
    };

    // This system runs every frame.  Assigning unconditionally would mark the
    // axis Camera as changed each frame and force it back through render-world
    // extraction for nothing.
    let unchanged = cam.viewport.as_ref().is_some_and(|current| {
        current.physical_position == wanted.physical_position
            && current.physical_size == wanted.physical_size
            && current.depth == wanted.depth
    });
    if unchanged {
        return;
    }

    cam.viewport = Some(wanted);
}

pub fn sync_axis_camera_to_main(
    q_main_cam: Query<&Transform, (With<Camera3d>, With<MainCamera>, Without<AxisCamera>)>,
    mut q_axis_cam: Query<&mut Transform, (With<Camera3d>, With<AxisCamera>, Without<MainCamera>)>,
) {
    let Ok(main_tf) = q_main_cam.single() else {
        return;
    };
    let Ok(mut axis_tf) = q_axis_cam.single_mut() else {
        return;
    };

    let radius = axis_tf.translation.length();
    let forward: Vec3 = main_tf.forward().as_vec3();
    let up: Vec3 = main_tf.up().as_vec3();
    let right: Vec3 = main_tf.right().as_vec3();

    axis_tf.translation = -forward.normalize() * radius;
    *axis_tf = Transform {
        translation: axis_tf.translation,
        rotation: Quat::from_mat3(&Mat3::from_cols(right, up, -forward)),
        ..*axis_tf
    };
}

fn compute_centroid(points: &[Vec3]) -> Option<Vec3> {
    if points.is_empty() {
        return None;
    }
    let mut acc = Vec3::ZERO;
    for &p in points {
        acc += p;
    }
    Some(acc / (points.len() as f32))
}

fn auto_fit_camera(points: &[Vec3], cam: &mut OrbitCamera) {
    let Some(c) = compute_centroid(points) else {
        return;
    };
    cam.target = c;

    let mut max_r = 0.0_f32;
    for &p in points {
        let r = p.distance(c);
        if r > max_r {
            max_r = r;
        }
    }

    // Keep a comfortable framing margin and clamp to usable bounds.
    let desired = (max_r * 3.0).max(2.0);
    cam.radius = desired.clamp(2.0, 200.0);
}

#[cfg(test)]
mod axis_gizmo_tests {
    use super::*;

    /// A camera looking down -Z with +Y up: the classic front view, where X
    /// runs right, Y runs up the screen, and Z points at the viewer.
    fn front_view() -> GlobalTransform {
        GlobalTransform::from(Transform::from_xyz(0.0, 0.0, 10.0).looking_at(Vec3::ZERO, Vec3::Y))
    }

    fn labels() -> [(&'static str, Vec2, Color); 3] {
        axis_gizmo_label_positions(UVec2::new(1000, 800), 1.0, &front_view())
    }

    /// The letters must clear the arrow heads (that was the complaint: they
    /// overlapped the arrows) while still staying on screen.
    #[test]
    fn labels_clear_the_arrow_tips_without_leaving_the_window() {
        let half = AXIS_VP_SIZE as f32 * 0.5;
        let centre = Vec2::new(
            1000.0 - AXIS_VP_MARGIN as f32 - half,
            800.0 - AXIS_VP_MARGIN as f32 - half,
        );
        let tip = (AXIS_TIP_LEN / AXIS_CAM_RADIUS) / (AXIS_CAM_FOV * 0.5).tan() * half;
        for (name, pos, _) in labels() {
            // Axes pointing at the viewer project to the centre and have no
            // tip on screen to clear.
            let offset = (pos - centre).length();
            if offset > 1.0 {
                assert!(
                    offset > tip,
                    "{name} label sits at {offset:.1} pt, inside the {tip:.1} pt arrow tip -- it \
                     would overlap the arrow"
                );
            }
            assert!(pos.x >= 0.0 && pos.x <= 1000.0, "{name} off screen: {pos:?}");
            assert!(pos.y >= 0.0 && pos.y <= 800.0, "{name} off screen: {pos:?}");
        }
    }

    /// A label that would land past the window edge is pulled back inside
    /// rather than clipped away.
    #[test]
    fn overshooting_labels_are_clamped_into_the_window() {
        // A tiny window makes the corner gizmo's overshoot run off the edge.
        let labels = axis_gizmo_label_positions(UVec2::new(200, 160), 1.0, &front_view());
        for (name, pos, _) in labels {
            assert!(
                pos.x >= AXIS_LABEL_PADDING - 0.01 && pos.x <= 200.0 - AXIS_LABEL_PADDING + 0.01,
                "{name} x out of bounds: {pos:?}"
            );
            assert!(
                pos.y >= AXIS_LABEL_PADDING - 0.01 && pos.y <= 160.0 - AXIS_LABEL_PADDING + 0.01,
                "{name} y out of bounds: {pos:?}"
            );
        }
    }

    #[test]
    fn x_is_right_and_y_is_up_in_a_front_view() {
        let half = AXIS_VP_SIZE as f32 * 0.5;
        let centre = Vec2::new(
            1000.0 - AXIS_VP_MARGIN as f32 - half,
            800.0 - AXIS_VP_MARGIN as f32 - half,
        );
        let found = labels();
        let x = found[0].1 - centre;
        let y = found[1].1 - centre;
        assert!(x.x > 1.0, "X should project to the right: {x:?}");
        assert!(x.y.abs() < 1.0, "X should not drift vertically: {x:?}");
        // Screen y grows downward, so "up" is negative.
        assert!(y.y < -1.0, "Y should project upward: {y:?}");
        assert!(y.x.abs() < 1.0, "Y should not drift horizontally: {y:?}");
    }

    #[test]
    fn an_axis_pointing_at_the_viewer_collapses_to_the_centre() {
        let half = AXIS_VP_SIZE as f32 * 0.5;
        let centre = Vec2::new(
            1000.0 - AXIS_VP_MARGIN as f32 - half,
            800.0 - AXIS_VP_MARGIN as f32 - half,
        );
        // Z points straight at the camera in a front view, so it has no
        // sideways component to project.
        let z = labels()[2].1 - centre;
        assert!(z.length() < 1.0, "Z should sit on the origin: {z:?}");
    }

    #[test]
    fn rotating_the_camera_moves_the_labels() {
        let side = GlobalTransform::from(
            Transform::from_xyz(10.0, 0.0, 0.0).looking_at(Vec3::ZERO, Vec3::Y),
        );
        let front_x = labels()[0].1;
        let side_x = axis_gizmo_label_positions(UVec2::new(1000, 800), 1.0, &side)[0].1;
        assert!(
            (front_x - side_x).length() > 5.0,
            "the labels must follow the camera, not stay put"
        );
    }

    /// The gizmo's viewport is a fixed size in *physical* pixels, so for one
    /// and the same physical window the labels must land on the same physical
    /// pixels whatever the scale factor -- which means their logical
    /// coordinates scale, rather than staying put.
    #[test]
    fn the_scale_factor_only_changes_logical_coordinates_not_physical_ones() {
        let physical = UVec2::new(2000, 1600);
        let at_1x = axis_gizmo_label_positions(physical, 1.0, &front_view())[0].1;
        let at_2x = axis_gizmo_label_positions(physical, 2.0, &front_view())[0].1;
        assert!(
            (at_1x - at_2x * 2.0).length() < 0.01,
            "same physical window should give the same physical position: \
             {at_1x:?} vs {:?}",
            at_2x * 2.0
        );
        assert!(
            (at_1x - at_2x).length() > 1.0,
            "the logical coordinates themselves must differ between scale factors"
        );
    }
}
