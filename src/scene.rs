use bevy::math::primitives::{Capsule3d, Cone, Cylinder, Sphere};
use bevy::prelude::*;
use bevy::render::camera::Viewport;
use bevy::render::view::RenderLayers;

use crate::color_schemes::color_for;
use crate::molecule::{covalent_radius_angstrom, Molecule};
use crate::settings::{LightingMode, MolSettings};

// default molecule (Å)
pub const DEFAULT_WATER: &str = r#"
O 0.00000   0.00000   0.22700
H 0.00000   1.35300  -0.90800
H 0.00000  -1.35300  -0.90800
"#;

#[derive(Component)] pub struct AtomMarker;
#[derive(Component)] pub struct BondMarker;
#[derive(Component)] pub struct MainCamera;
#[derive(Component)] pub struct AxisCamera;

#[derive(Component)] pub struct KeyLight;
#[derive(Component)] pub struct FillLight;
#[derive(Component)] pub struct RimLight;

pub const LAYER_MAIN: usize = 0;
pub const LAYER_AXES: usize = 1;

const AXIS_VP_SIZE: u32 = 140;
const AXIS_VP_MARGIN: u32 = 12;

pub fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    settings: Res<MolSettings>,
) {
    commands.insert_resource(ClearColor(settings.bg_color));

    // Ambient (resource)
    commands.insert_resource(AmbientLight {
        brightness: settings.ambient_brightness,
        color: settings.ambient_color,
        affects_lightmapped_meshes: true,
    });

    // Main camera
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(0.0, 7.0, 14.0).looking_at(Vec3::new(0.0, 1.0, 0.0), Vec3::Y),
        MainCamera,
        RenderLayers::layer(LAYER_MAIN),
    ));

    // Key light (front-right-up)
    commands.spawn((
        PointLight {
            shadows_enabled: false,
            intensity: settings.light_intensity,
            range: 200.0,
            ..default()
        },
        Transform::from_xyz(settings.light_distance, settings.light_distance, settings.light_distance),
        KeyLight,
        RenderLayers::layer(LAYER_MAIN),
    ));

    // Fill (front-left)
    commands.spawn((
        PointLight {
            shadows_enabled: false,
            intensity: settings.fill_intensity,
            range: 200.0,
            ..default()
        },
        Transform::from_xyz(-settings.fill_distance, settings.fill_distance * 0.5, settings.fill_distance),
        FillLight,
        RenderLayers::layer(LAYER_MAIN),
    ));

    // Rim (back-right)
    commands.spawn((
        PointLight {
            shadows_enabled: false,
            intensity: settings.rim_intensity,
            range: 200.0,
            ..default()
        },
        Transform::from_xyz(settings.rim_distance, settings.rim_distance * 0.5, -settings.rim_distance),
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
        .spawn((Transform::default(), Visibility::Visible, RenderLayers::layer(LAYER_AXES)))
        .id();

    let to_x = Quat::from_rotation_arc(Vec3::Y, Vec3::X);
    let to_z = Quat::from_rotation_arc(Vec3::Y, Vec3::Z);

    let x_shaft = commands
        .spawn((
            Mesh3d(shaft.clone()),
            MeshMaterial3d(red.clone()),
            Transform { rotation: to_x, translation: Vec3::X * (axis_len * 0.5), ..default() },
            RenderLayers::layer(LAYER_AXES),
        ))
        .id();
    commands.entity(x_shaft).insert(ChildOf(axes_root));
    let x_head = commands
        .spawn((
            Mesh3d(head.clone()),
            MeshMaterial3d(red),
            Transform { rotation: to_x, translation: Vec3::X * (axis_len + head_len * 0.5), ..default() },
            RenderLayers::layer(LAYER_AXES),
        ))
        .id();
    commands.entity(x_head).insert(ChildOf(axes_root));

    let y_shaft = commands
        .spawn((
            Mesh3d(shaft.clone()),
            MeshMaterial3d(green.clone()),
            Transform { translation: Vec3::Y * (axis_len * 0.5), ..default() },
            RenderLayers::layer(LAYER_AXES),
        ))
        .id();
    commands.entity(y_shaft).insert(ChildOf(axes_root));
    let y_head = commands
        .spawn((
            Mesh3d(head.clone()),
            MeshMaterial3d(green),
            Transform { translation: Vec3::Y * (axis_len + head_len * 0.5), ..default() },
            RenderLayers::layer(LAYER_AXES),
        ))
        .id();
    commands.entity(y_head).insert(ChildOf(axes_root));

    let z_shaft = commands
        .spawn((
            Mesh3d(shaft),
            MeshMaterial3d(blue.clone()),
            Transform { rotation: to_z, translation: Vec3::Z * (axis_len * 0.5), ..default() },
            RenderLayers::layer(LAYER_AXES),
        ))
        .id();
    commands.entity(z_shaft).insert(ChildOf(axes_root));
    let z_head = commands
        .spawn((
            Mesh3d(head),
            MeshMaterial3d(blue),
            Transform { rotation: to_z, translation: Vec3::Z * (axis_len + head_len * 0.5), ..default() },
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
        Camera { order: 1, clear_color: ClearColorConfig::None, viewport: None, ..default() },
    ));
}

pub fn rebuild_if_dirty(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut clear: ResMut<ClearColor>,
    mut mol: ResMut<Molecule>,
    mut settings: ResMut<MolSettings>,
    q_atoms: Query<Entity, With<AtomMarker>>,
    q_bonds: Query<Entity, With<BondMarker>>,
    mut q_key: Query<(&mut PointLight, &mut Transform), (With<KeyLight>, Without<FillLight>, Without<RimLight>)>,
    mut q_fill: Query<(&mut PointLight, &mut Transform), (With<FillLight>, Without<KeyLight>, Without<RimLight>)>,
    mut q_rim:  Query<(&mut PointLight, &mut Transform), (With<RimLight>,  Without<KeyLight>, Without<FillLight>)>,
    mut amb: ResMut<AmbientLight>,
) {
    if !settings.dirty { return; }
    settings.dirty = false;

    // Background
    *clear = ClearColor(settings.bg_color);

    // Ambient
    amb.color = settings.ambient_color;
    amb.brightness = if settings.use_ambient { settings.ambient_brightness } else { 0.0 };

    // Update lights to current settings
    if let Ok((mut l, mut tf)) = q_key.get_single_mut() {
        l.intensity = settings.light_intensity;
        tf.translation = Vec3::new(settings.light_distance, settings.light_distance, settings.light_distance);
    }
    if let Ok((mut l, mut tf)) = q_fill.get_single_mut() {
        l.intensity = if matches!(settings.lighting_mode, LightingMode::ThreePoint) { settings.fill_intensity } else { 0.0 };
        tf.translation = Vec3::new(-settings.fill_distance, settings.fill_distance * 0.5, settings.fill_distance);
    }
    if let Ok((mut l, mut tf)) = q_rim.get_single_mut() {
        l.intensity = if matches!(settings.lighting_mode, LightingMode::ThreePoint) { settings.rim_intensity } else { 0.0 };
        tf.translation = Vec3::new(settings.rim_distance, settings.rim_distance * 0.5, -settings.rim_distance);
    }

    // Despawn old geometry
    for e in q_atoms.iter() { commands.entity(e).despawn(); }
    for e in q_bonds.iter() { commands.entity(e).despawn(); }
    mol.atom_entities.clear();
    mol.bond_entities.clear();

    mol.recompute_bonds(settings.bond_thresh_scale);

    // material parameters from settings
    let metallic = settings.metallic.clamp(0.0, 1.0);
    let rough    = settings.roughness.clamp(0.02, 1.0);
    let refl     = settings.reflectance.clamp(0.0, 1.0);

    let emissive = if settings.use_emissive {
        let lr = settings.emissive_color.to_linear();
        LinearRgba::new(
            lr.red   * settings.emissive_strength.max(0.0),
            lr.green * settings.emissive_strength.max(0.0),
            lr.blue  * settings.emissive_strength.max(0.0),
            1.0,
        )
    } else {
        LinearRgba::BLACK
    };

    // Atoms
    for i in 0..mol.atoms.len() {
        let sym = &mol.atoms[i];
        let r_cov = covalent_radius_angstrom(sym) * settings.atom_scale;

        // If you prefer the palette function, use it; otherwise the overrides map:
        let color = settings
            .element_colors
            .get(sym)
            .copied()
            .unwrap_or_else(|| color_for(sym, settings.scheme));

        let sphere_mesh = Sphere::new(r_cov).mesh().ico(5).unwrap();

        let mat = materials.add(StandardMaterial {
            base_color: color,
            metallic,
            perceptual_roughness: rough,
            reflectance: refl,
            emissive,
            ..default()
        });

        let pos = mol.pos[i];
        let ent = commands
            .spawn((
                Mesh3d(meshes.add(sphere_mesh)),
                MeshMaterial3d(mat),
                Transform::from_translation(pos),
                AtomMarker,
                RenderLayers::layer(LAYER_MAIN),
            ))
            .id();
        mol.atom_entities.push(ent);
    }

    // Bonds
    let bonds = mol.bonds.clone();
    for (i, j, d) in bonds {
        let p0 = mol.pos[i];
        let p1 = mol.pos[j];
        let dir = p1 - p0;
        let len = d;
        if len <= 0.0001 { continue; }

        let center = (p0 + p1) * 0.5;
        let rot = Quat::from_rotation_arc(Vec3::Y, dir.normalize());

        // radius based on smaller atom radius
        let ri = covalent_radius_angstrom(&mol.atoms[i]) * settings.atom_scale;
        let rj = covalent_radius_angstrom(&mol.atoms[j]) * settings.atom_scale;
        let base = ri.min(rj);
        let radius = (base * settings.bond_radius_pct).clamp(0.01, base);

        let half_length = (len * 0.5 - radius).max(0.0);
        let cap = Mesh::from(Capsule3d { radius, half_length });

        let mat = materials.add(StandardMaterial {
            base_color: Color::srgb(0.85, 0.85, 0.88),
            metallic,
            perceptual_roughness: rough,
            reflectance: refl,
            emissive,
            ..default()
        });

        let ent = commands
            .spawn((
                Mesh3d(meshes.add(cap)),
                MeshMaterial3d(mat),
                Transform { translation: center, rotation: rot, ..default() },
                BondMarker,
                RenderLayers::layer(LAYER_MAIN),
            ))
            .id();
        mol.bond_entities.push(ent);
    }
}

// ---------------------- Axis viewport & sync ----------------------

pub fn update_axis_viewport_on_resize(
    windows: Query<&Window>,
    mut q_axis_cam: Query<&mut Camera, With<AxisCamera>>,
) {
    let Ok(window) = windows.single() else { return; };
    let Ok(mut cam) = q_axis_cam.single_mut() else { return; };

    let w = window.physical_width();
    let h = window.physical_height();

    if w == 0 || h == 0 {
        cam.viewport = None;
        return;
    }

    let vp_w = AXIS_VP_SIZE;
    let vp_h = AXIS_VP_SIZE;
    let x = AXIS_VP_MARGIN;
    let y = h.saturating_sub(vp_h + AXIS_VP_MARGIN);

    cam.viewport = Some(Viewport {
        physical_position: UVec2::new(x, y),
        physical_size: UVec2::new(vp_w, vp_h),
        ..default()
    });
}

pub fn sync_axis_camera_to_main(
    q_main_cam: Query<&Transform, (With<Camera3d>, With<MainCamera>, Without<AxisCamera>)>,
    mut q_axis_cam: Query<&mut Transform, (With<Camera3d>, With<AxisCamera>, Without<MainCamera>)>,
) {
    let Ok(main_tf) = q_main_cam.single() else { return; };
    let Ok(mut axis_tf) = q_axis_cam.single_mut() else { return; };

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

