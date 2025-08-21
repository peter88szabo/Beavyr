use bevy::math::primitives::{Capsule3d, Cone, Cylinder, Sphere};
use bevy::prelude::*;
use bevy::render::camera::Viewport;
use bevy::render::view::RenderLayers;

use crate::color_schemes::color_for;
use crate::molecule::{covalent_radius_angstrom, Molecule};
use crate::settings::{BondColorMode, LightingMode, MolSettings};

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

    mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);

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

        let color = settings
            .element_colors
            .get(sym)
            .copied()
            .unwrap_or_else(|| color_for(sym, settings.scheme));

        let sphere_mesh = Sphere::new(r_cov).mesh().ico(settings.atom_resolution).unwrap();
        //let sphere_mesh = Sphere::new(r_cov).mesh().ico(6).unwrap();

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
    for (i, j, len_cc) in bonds {
        let p0 = mol.pos[i];
        let p1 = mol.pos[j];
        if len_cc <= 0.0001 { continue; }

        let dir = p1 - p0;
        let dir_n = dir.normalize();

        // Visual sphere radii (atom spheres)
        let ri_sphere = covalent_radius_angstrom(&mol.atoms[i]) * settings.atom_scale;
        let rj_sphere = covalent_radius_angstrom(&mol.atoms[j]) * settings.atom_scale;

        // Bond cylinder/capsule radius based on smaller atom radius
        let base = ri_sphere.min(rj_sphere);
        let radius = (base * settings.bond_radius_pct).clamp(0.01, base);

        // Common rotation (align local +Y with bond direction)
        let rot = Quat::from_rotation_arc(Vec3::Y, dir_n);

        match settings.bond_color_mode {
            BondColorMode::Uniform => {
                // Single capsule centered between atom centers, like the original
                let center = (p0 + p1) * 0.5;
                // This is the original half-length formula
                // half_length = (center-to-center distance)/2 - radius
                let half_length = (len_cc * 0.5 - radius).max(0.0);
                let cap = Mesh::from(Capsule3d { radius, half_length });

                let mat = materials.add(StandardMaterial {
                    base_color: settings.uniform_bond_color,
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

            // flat cylinders, slight overlap into spheres) ===
            BondColorMode::AtomSplit => {
                // Push the visible segment slightly inside each atom sphere to guarantee no seam.
                // Using a generous overlap fraction looks best on bright backgrounds.
                let overlap = radius * 0.25; // 25% of bond radius
                let start = p0 + dir_n * (ri_sphere - overlap);
                let end   = p1 - dir_n * (rj_sphere - overlap);

                let visible_len = (end - start).length();
                if visible_len <= 0.0 { continue; }

                let half_len = visible_len * 0.5;

                // Centers for each half (flat ends meet at the exact midpoint)
                let mid = (start + end) * 0.5;
                let c0 = mid - dir_n * (half_len * 0.5);
                let c1 = mid + dir_n * (half_len * 0.5);

                let cyl_half = Mesh::from(Cylinder::new(radius, half_len));

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

                let mat_i = materials.add(StandardMaterial {
                    base_color: col_i,
                    metallic,
                    perceptual_roughness: rough,
                    reflectance: refl,
                    emissive,
                    ..default()
                });
                let mat_j = materials.add(StandardMaterial {
                    base_color: col_j,
                    metallic,
                    perceptual_roughness: rough,
                    reflectance: refl,
                    emissive,
                    ..default()
                });

                let e0 = commands
                    .spawn((
                        Mesh3d(meshes.add(cyl_half.clone())),
                        MeshMaterial3d(mat_i),
                        Transform { translation: c0, rotation: rot, ..default() },
                        BondMarker,
                        RenderLayers::layer(LAYER_MAIN),
                    ))
                    .id();
                mol.bond_entities.push(e0);

                let e1 = commands
                    .spawn((
                        Mesh3d(meshes.add(cyl_half)),
                        MeshMaterial3d(mat_j),
                        Transform { translation: c1, rotation: rot, ..default() },
                        BondMarker,
                        RenderLayers::layer(LAYER_MAIN),
                    ))
                    .id();
                mol.bond_entities.push(e1);
            }
        }
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

