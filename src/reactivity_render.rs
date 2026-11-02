//! Independent signed reactivity surfaces and atomic halos. Stale jobs cannot replace a new frame.
use super::*;
use bevy::{asset::RenderAssetUsages, camera::visibility::RenderLayers, mesh::PrimitiveTopology};
#[derive(Component)]
struct Surface;
#[derive(Component)]
struct Halo;
#[derive(Resource, Default)]
struct Cache {
    key: Option<(usize, usize, u32, u32, data::Descriptor)>,
    task: Option<Task<Result<[crate::orbitals::isosurface::SurfaceMesh; 2], String>>>,
    halo_key: Vec<u64>,
    sphere: Option<Handle<Mesh>>,
    materials: Vec<Handle<StandardMaterial>>,
    appearance: Vec<u32>,
}
pub struct Plugin;
impl bevy::prelude::Plugin for Plugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<Cache>()
            .add_systems(
                PostUpdate,
                update
                    .after(bevy_egui::EguiPostUpdateSet::EndPass)
                    .before(bevy::transform::TransformSystems::Propagate),
            )
            .add_systems(EguiPrimaryContextPass, labels.after(super::ui));
    }
}
#[allow(clippy::too_many_arguments)]
fn update(
    mut commands: Commands,
    mut state: ResMut<State>,
    traj: Res<TrajectoryState>,
    mol: Res<crate::molecule::Molecule>,
    mut cache: ResMut<Cache>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    surfaces: Query<Entity, With<Surface>>,
    halos: Query<Entity, With<Halo>>,
    builder: Res<crate::molecule_builder::builder_ui::ZMatrixBuilderState>,
) {
    let current = state.current_index(&mol, &traj);
    let result = state.results.clone();
    if cache.materials.is_empty() {
        for _ in 0..2 {
            cache.materials.push(materials.add(StandardMaterial {
                alpha_mode: AlphaMode::Blend,
                cull_mode: None,
                perceptual_roughness: 0.45,
                ..default()
            }));
        }
    }
    let appearance = vec![
        state.opacity.to_bits(),
        state.positive.to_srgba().red.to_bits(),
        state.positive.to_srgba().green.to_bits(),
        state.positive.to_srgba().blue.to_bits(),
        state.negative.to_srgba().red.to_bits(),
        state.negative.to_srgba().green.to_bits(),
        state.negative.to_srgba().blue.to_bits(),
    ];
    if cache.appearance != appearance {
        cache.appearance = appearance;
        for (handle, color) in cache.materials.iter().zip([state.positive, state.negative]) {
            if let Some(mut m) = materials.get_mut(handle) {
                let c = color.to_srgba();
                m.base_color = Color::srgba(c.red, c.green, c.blue, state.opacity);
                m.alpha_mode = if state.opacity >= 0.999 {
                    AlphaMode::Opaque
                } else {
                    AlphaMode::Blend
                };
            }
        }
    }
    let key = current.and_then(|i| {
        let r = result.as_ref()?;
        if !state.surface || !state.descriptor.local() {
            return None;
        }
        Some((
            Arc::as_ptr(r) as usize,
            i,
            state.isovalue.to_bits(),
            state.resolution,
            state.descriptor,
        ))
    });
    if cache.key != key {
        cache.key = key;
        cache.task = None;
        state.ready = None;
        for e in &surfaces {
            commands.entity(e).despawn();
        }
        if let Some((_, i, _, _, descriptor)) = key {
            let r = result.clone().unwrap();
            let iso = state.isovalue;
            let mut settings = crate::orbitals::OrbitalState::default();
            settings.resolution = state.resolution;
            let spacing = settings.grid_spacing();
            cache.task = Some(AsyncComputeTaskPool::get().spawn(async move {
                let field = data::field(&r.frames[i], descriptor, spacing)?;
                Ok([
                    crate::orbitals::isosurface::extract(&field, iso),
                    crate::orbitals::isosurface::extract(&field, -iso),
                ])
            }));
        }
    }
    if let Some(task) = &mut cache.task {
        if let Some(result) = block_on(future::poll_once(task)) {
            cache.task = None;
            match result {
                Ok(lobes) => {
                    for (i, lobe) in lobes.into_iter().enumerate() {
                        if lobe.is_empty() {
                            continue;
                        }
                        let mut mesh = Mesh::new(
                            PrimitiveTopology::TriangleList,
                            RenderAssetUsages::default(),
                        );
                        mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, lobe.positions);
                        mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, lobe.normals);
                        commands.spawn((
                            Surface,
                            Mesh3d(meshes.add(mesh)),
                            MeshMaterial3d(cache.materials[i].clone()),
                            Transform::default(),
                            RenderLayers::layer(crate::scene::LAYER_MAIN),
                        ));
                    }
                    state.ready = current;
                }
                Err(e) => {
                    state.error = Some(e);
                    state.surface = false;
                    state.ready = None;
                }
            }
        }
    }
    let values = current
        .and_then(|i| state.descriptor.atoms(&result.as_ref()?.frames[i]))
        .filter(|_| state.halos);
    let mut halo_key = vec![state.scale.to_bits() as u64];
    if let Some(v) = &values {
        for (i, (&x, p)) in v.iter().zip(&mol.pos).enumerate() {
            halo_key.extend([
                x.to_bits(),
                p.x.to_bits() as u64,
                p.y.to_bits() as u64,
                p.z.to_bits() as u64,
                builder.selection.hidden.contains(&i) as u64,
            ]);
        }
    }
    if halo_key != cache.halo_key {
        cache.halo_key = halo_key;
        for e in &halos {
            commands.entity(e).despawn();
        }
        if let Some(values) = values {
            let sphere = cache
                .sphere
                .get_or_insert_with(|| meshes.add(Sphere::new(1.).mesh().uv(32, 20)))
                .clone();
            for (i, (&value, &position)) in values.iter().zip(&mol.pos).enumerate() {
                if builder.selection.hidden.contains(&i) || value.abs() < 1e-12 {
                    continue;
                }
                let radius = state.scale * value.abs().cbrt() as f32;
                commands.spawn((
                    Halo,
                    Mesh3d(sphere.clone()),
                    MeshMaterial3d(cache.materials[usize::from(value < 0.)].clone()),
                    Transform::from_translation(position).with_scale(Vec3::splat(radius)),
                    RenderLayers::layer(crate::scene::LAYER_MAIN),
                ));
            }
        }
    }
}
fn labels(
    mut contexts: EguiContexts,
    state: Res<State>,
    mol: Res<crate::molecule::Molecule>,
    traj: Res<TrajectoryState>,
    camera: Query<(&Camera, &GlobalTransform), With<crate::scene::MainCamera>>,
    windows: Query<&Window, With<bevy::window::PrimaryWindow>>,
    regions: Res<crate::ui::UiPanelRegions>,
    builder: Res<crate::molecule_builder::builder_ui::ZMatrixBuilderState>,
) {
    if !state.labels {
        return;
    }
    let Some(i) = state.current_index(&mol, &traj) else {
        return;
    };
    let Some(values) = state
        .descriptor
        .atoms(&state.results.as_ref().unwrap().frames[i])
    else {
        return;
    };
    let (Ok(ctx), Ok((camera, transform)), Ok(window)) =
        (contexts.ctx_mut(), camera.single(), windows.single())
    else {
        return;
    };
    let painter = ctx.layer_painter(egui::LayerId::new(
        egui::Order::Background,
        egui::Id::new("reactivity_labels"),
    ));
    for (i, (&value, &position)) in values.iter().zip(&mol.pos).enumerate() {
        if builder.selection.hidden.contains(&i) {
            continue;
        }
        if let Ok(p) = camera.world_to_viewport(transform, position) {
            let p = p / window.scale_factor() as f32;
            let p = egui::pos2(p.x + 8., p.y + 8.);
            if regions.covers(Vec2::new(p.x, p.y)) {
                continue;
            }
            let text = format!("{value:+.4}");
            for (offset, color) in [
                (egui::vec2(1., 1.), egui::Color32::BLACK),
                (egui::Vec2::ZERO, egui::Color32::WHITE),
            ] {
                painter.text(
                    p + offset,
                    egui::Align2::LEFT_TOP,
                    &text,
                    egui::FontId::proportional(14.),
                    color,
                );
            }
        }
    }
}
