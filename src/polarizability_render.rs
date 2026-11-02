//! Atom-centred translucent CPA spheres and optional numeric overlays.
use super::*;
use bevy::camera::visibility::RenderLayers;
#[derive(Component)]
struct Halo;
#[derive(Resource, Default)]
struct Cache {
    key: Vec<u64>,
    sphere: Option<Handle<Mesh>>,
    materials: Vec<Handle<StandardMaterial>>,
}
pub struct Plugin;
impl bevy::prelude::Plugin for Plugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<Cache>()
            .add_systems(
                PostUpdate,
                halos
                    .after(bevy_egui::EguiPostUpdateSet::EndPass)
                    .before(bevy::transform::TransformSystems::Propagate),
            )
            .add_systems(EguiPrimaryContextPass, labels.after(super::ui));
    }
}
fn halos(
    mut commands: Commands,
    state: Res<State>,
    mol: Res<crate::molecule::Molecule>,
    traj: Res<TrajectoryState>,
    mut cache: ResMut<Cache>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    existing: Query<Entity, With<Halo>>,
    builder: Res<crate::molecule_builder::builder_ui::ZMatrixBuilderState>,
) {
    let values = state.current(&mol, &traj).filter(|_| state.show);
    let mut key = vec![state.scale.to_bits() as u64, state.opacity.to_bits() as u64];
    key.extend(state.colors.bounds.map(f64::to_bits));
    if let Some(v) = values {
        for (i, (&a, p)) in v.alpha.iter().zip(&mol.pos).enumerate() {
            key.extend([
                a.to_bits(),
                p.x.to_bits() as u64,
                p.y.to_bits() as u64,
                p.z.to_bits() as u64,
                builder.selection.hidden.contains(&i) as u64,
            ]);
        }
    }
    if cache.key == key {
        return;
    }
    cache.key = key;
    for e in &existing {
        commands.entity(e).despawn();
    }
    let Some(values) = values else {
        return;
    };
    let sphere = cache
        .sphere
        .get_or_insert_with(|| meshes.add(Sphere::new(1.).mesh().uv(32, 20)))
        .clone();
    cache.materials.clear();
    for (i, (&alpha, &pos)) in values.alpha.iter().zip(&mol.pos).enumerate() {
        if builder.selection.hidden.contains(&i) {
            continue;
        }
        let radius = data::radius(alpha, state.scale);
        if radius > 0. && radius.is_finite() {
            let material=materials.add(StandardMaterial{
                base_color:state.colors.color(alpha,Color::srgb(0.25,0.75,1.),Color::srgb(0.25,0.75,1.),state.opacity),
                alpha_mode:AlphaMode::Blend,perceptual_roughness:0.65,..default()
            });
            cache.materials.push(material.clone());
            commands.spawn((
                Halo,
                Mesh3d(sphere.clone()),
                MeshMaterial3d(material.clone()),
                Transform::from_translation(pos).with_scale(Vec3::splat(radius)),
                RenderLayers::layer(crate::scene::LAYER_MAIN),
            ));
        }
    }
}
fn labels(
    mut contexts: EguiContexts,
    state: Res<State>,
    mol: Res<crate::molecule::Molecule>,
    traj: Res<TrajectoryState>,
    cameras: Query<(&Camera, &GlobalTransform), With<crate::scene::MainCamera>>,
    windows: Query<&Window, With<bevy::window::PrimaryWindow>>,
    regions: Res<crate::ui::UiPanelRegions>,
    builder: Res<crate::molecule_builder::builder_ui::ZMatrixBuilderState>,
) {
    if !state.labels {
        return;
    }
    let Some(values) = state.current(&mol, &traj) else {
        return;
    };
    let (Ok(ctx), Ok((camera, transform)), Ok(window)) =
        (contexts.ctx_mut(), cameras.single(), windows.single())
    else {
        return;
    };
    let painter = ctx.layer_painter(egui::LayerId::new(
        egui::Order::Background,
        egui::Id::new("atomic_polarizability_labels"),
    ));
    for (i, (&alpha, &position)) in values.alpha.iter().zip(&mol.pos).enumerate() {
        if builder.selection.hidden.contains(&i) {
            continue;
        }
        if let Ok(p) = camera.world_to_viewport(transform, position) {
            let p = p / window.scale_factor() as f32;
            let p = egui::pos2(p.x + 8., p.y - 8.);
            if regions.covers(Vec2::new(p.x, p.y)) {
                continue;
            }
            let text = format!("{:.3} {}", alpha * state.factor(), state.unit());
            for (offset, color) in [
                (egui::vec2(1., 1.), egui::Color32::BLACK),
                (egui::Vec2::ZERO, egui::Color32::WHITE),
            ] {
                painter.text(
                    p + offset,
                    egui::Align2::LEFT_BOTTOM,
                    &text,
                    egui::FontId::proportional(14.),
                    color,
                );
            }
        }
    }
}
