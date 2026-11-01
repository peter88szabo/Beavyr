//! Two isolated molecular views. Comparison transforms never touch the main molecule.
use crate::{
    camera::OrbitCamera,
    isolated_scene::Scene,
    molecule::Molecule,
    settings::MolSettings,
    trajectory::{TrajectoryFrame, TrajectoryState},
    ui_style::ResponseExt,
};
use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts, EguiPrimaryContextPass, EguiTextureHandle};
const OPEN: &str = "open_side_by_side";
pub struct Plugin;
impl bevy::prelude::Plugin for Plugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<State>()
            .add_systems(EguiPrimaryContextPass, ui.after(crate::ui::ui_panel));
    }
}
struct Side {
    molecule: Option<Molecule>,
    frames: Vec<TrajectoryFrame>,
    frame: usize,
    title: String,
    orbit: OrbitCamera,
    scene: Option<Scene>,
    texture: Option<egui::TextureId>,
    dirty: bool,
    size: UVec2,
}
impl Default for Side {
    fn default() -> Self {
        Self {
            molecule: None,
            frames: vec![],
            frame: 0,
            title: String::new(),
            orbit: Default::default(),
            scene: None,
            texture: None,
            dirty: true,
            size: UVec2::new(600, 600),
        }
    }
}
#[derive(Resource)]
struct State {
    open: bool,
    sync: bool,
    sides: [Side; 2],
    message: String,
    alignment_undo: Option<(Vec<Vec3>, OrbitCamera)>,
}
impl Default for State {
    fn default() -> Self {
        Self {
            open: false,
            sync: true,
            sides: Default::default(),
            message: String::new(),
            alignment_undo: None,
        }
    }
}
pub fn button(ui: &mut egui::Ui) {
    if ui.button("Side by side…").clicked_once() {
        ui.ctx()
            .data_mut(|d| d.insert_temp(egui::Id::new(OPEN), true));
    }
}
fn fit(orbit: &mut OrbitCamera, mol: &Molecule) {
    if mol.pos.is_empty() {
        return;
    }
    orbit.target = crate::rmsd::kabsch::centroid(&mol.pos);
    orbit.radius = (mol
        .pos
        .iter()
        .map(|p| p.distance(orbit.target))
        .fold(0., f32::max)
        * 2.8)
        .max(2.);
    orbit.rotating = false;
    orbit.panning = false;
}
fn motion(orbit: &mut OrbitCamera, rotation: Vec2, pan: Vec2, zoom: f32) {
    let yaw = Quat::from_rotation_y(-rotation.x * 0.006);
    let pitch = Quat::from_axis_angle(orbit.orientation * Vec3::X, -rotation.y * 0.006);
    orbit.orientation = (yaw * pitch * orbit.orientation).normalize();
    orbit.radius = (orbit.radius * zoom).clamp(0.1, 100000.);
    orbit.target += orbit.orientation * Vec3::new(-pan.x, pan.y, 0.) * orbit.radius * 0.002;
}
fn load(
    path: &std::path::Path,
    settings: &MolSettings,
) -> Result<(Molecule, Vec<TrajectoryFrame>), String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    if crate::molecule_builder::protein_import::is_protein_path(path) {
        let cif = path
            .extension()
            .and_then(|s| s.to_str())
            .is_some_and(|s| matches!(s.to_ascii_lowercase().as_str(), "cif" | "mmcif"));
        let mut mol = crate::molecule_builder::protein_import::parse(&text, cif)?.molecule;
        mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
        return Ok((mol, vec![]));
    }
    let frames = crate::trajectory::parse_multi_xyz(&text)?;
    let first = frames.first().ok_or("No frames found.")?;
    let mut mol = Molecule::empty();
    mol.atoms = first.atoms.clone();
    mol.pos = first.pos.clone();
    mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
    Ok((mol, frames))
}
#[allow(clippy::too_many_arguments)]
fn ui(
    mut commands: Commands,
    mut contexts: EguiContexts,
    mut state: ResMut<State>,
    mol: Res<Molecule>,
    traj: Res<TrajectoryState>,
    settings: Res<MolSettings>,
    orbit: Res<OrbitCamera>,
    camera: Query<&Projection, With<crate::scene::MainCamera>>,
    mut images: ResMut<Assets<Image>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut regions: ResMut<crate::ui::UiPanelRegions>,
) {
    let Ok(ctx) = contexts.ctx_mut().map(|c| c.clone()) else {
        return;
    };
    if ctx
        .data_mut(|d| d.remove_temp::<bool>(egui::Id::new(OPEN)))
        .unwrap_or(false)
    {
        state.open = true;
        for side in &mut state.sides {
            if side.molecule.is_none() && !mol.atoms.is_empty() {
                side.molecule = Some(mol.clone());
                side.orbit = orbit.clone();
                side.title = "Current structure".into();
                side.dirty = true;
            }
        }
    }
    if !state.open {
        for side in &mut state.sides {
            if let Some(scene) = side.scene.take() {
                contexts.remove_image(scene.image.id());
                scene.dispose(&mut commands, &mut images);
                side.texture = None;
                side.dirty = true;
            }
        }
        return;
    }
    let Ok(projection) = camera.single() else {
        return;
    };
    for (index, side) in state.sides.iter_mut().enumerate() {
        if side.scene.is_none() && side.molecule.is_some() {
            let scene = Scene::new(
                &mut commands,
                &mut images,
                &mut meshes,
                &settings,
                &side.orbit,
                projection,
                side.size,
                21 + index,
            );
            side.texture = Some(contexts.add_image(EguiTextureHandle::Strong(scene.image.clone())));
            side.scene = Some(scene);
            side.dirty = true;
        }
        if side.dirty {
            if let (Some(scene), Some(molecule)) = (&mut side.scene, &side.molecule) {
                scene.molecule(
                    &mut commands,
                    &mut materials,
                    molecule,
                    &settings,
                    &Default::default(),
                );
                side.dirty = false;
            }
        }
    }
    let mut open = true;
    let mut actions = vec![];
    let mut movements = vec![];
    let mut sync_changed = false;
    let mut align = false;
    let mut restore = false;
    let window = egui::Window::new("Side-by-side molecule views")
        .open(&mut open)
        .resizable(true)
        .default_size([1050., 620.])
        .min_size([480., 300.])
        .show(&ctx, |ui| {
            crate::ui_style::modern(ui);
            ui.horizontal_wrapped(|ui| {
                sync_changed = ui
                    .checkbox(&mut state.sync, "Synchronize rotation and zoom")
                    .changed();
                if ui.button("Align B to A").clicked_once() {
                    align = true;
                }
                if ui
                    .add_enabled(
                        state.alignment_undo.is_some(),
                        egui::Button::new("Undo alignment"),
                    )
                    .clicked_once()
                {
                    restore = true;
                }
            });
            if !state.message.is_empty() {
                ui.label(&state.message);
            }
            let height = (ui.available_height() - 85.).max(160.);
            ui.columns(2, |columns| {
                for (i, ui) in columns.iter_mut().enumerate() {
                    ui.push_id(i, |ui| {
                        let side = &mut state.sides[i];
                        ui.horizontal_wrapped(|ui| {
                            ui.strong(if i == 0 { "A" } else { "B" });
                            if ui.button("Capture current").clicked_once() {
                                actions.push((i, 0));
                            }
                            if ui.button("Open…").clicked_once() {
                                actions.push((i, 1));
                            }
                            if ui
                                .add_enabled(
                                    !traj.frames.is_empty(),
                                    egui::Button::new("Use trajectory"),
                                )
                                .clicked_once()
                            {
                                actions.push((i, 2));
                            }
                            if ui.button("Fit").clicked_once() {
                                actions.push((i, 3));
                            }
                        });
                        ui.add(egui::Label::new(&side.title).truncate());
                        if side.frames.len() > 1 {
                            let mut frame = side.frame + 1;
                            if ui
                                .add(
                                    egui::Slider::new(&mut frame, 1..=side.frames.len())
                                        .text("Frame"),
                                )
                                .changed()
                            {
                                side.frame = frame - 1;
                                actions.push((i, 4));
                            }
                        }
                        let size = egui::vec2(ui.available_width(), height);
                        let physical = UVec2::new(
                            (size.x * ctx.pixels_per_point()).round().clamp(64., 2048.) as u32,
                            (size.y * ctx.pixels_per_point()).round().clamp(64., 2048.) as u32,
                        );
                        if physical != side.size {
                            side.size = physical;
                            if let Some(scene) = &side.scene {
                                if let Some(mut image) = images.get_mut(scene.image.id()) {
                                    image.resize(bevy::render::render_resource::Extent3d {
                                        width: physical.x,
                                        height: physical.y,
                                        depth_or_array_layers: 1,
                                    });
                                }
                            }
                        }
                        if let Some(texture) = side.texture {
                            let response = ui.add(
                                egui::Image::new((texture, size))
                                    .fit_to_exact_size(size)
                                    .sense(egui::Sense::drag()),
                            );
                            let delta = ui.input(|input| input.pointer.delta());
                            let rotation = if response.dragged_by(egui::PointerButton::Primary) {
                                Vec2::new(delta.x, delta.y)
                            } else {
                                Vec2::ZERO
                            };
                            let pan = if response.dragged_by(egui::PointerButton::Secondary) {
                                Vec2::new(delta.x, delta.y)
                            } else {
                                Vec2::ZERO
                            };
                            let zoom = if response.hovered() {
                                (-ui.input(|input| input.smooth_scroll_delta.y) * 0.002).exp()
                            } else {
                                1.
                            };
                            if rotation != Vec2::ZERO || pan != Vec2::ZERO || zoom != 1. {
                                movements.push((i, rotation, pan, zoom));
                            }
                        } else {
                            ui.allocate_space(size);
                        }
                    });
                }
            });
        });
    state.open = open;
    if let Some(w) = window {
        let r = w.response.rect;
        regions
            .windows
            .push(Rect::new(r.min.x, r.min.y, r.max.x, r.max.y));
    }
    if sync_changed && state.sync {
        state.sides[1].orbit.orientation = state.sides[0].orbit.orientation;
    }
    for (i, rotation, pan, zoom) in movements {
        let sync = state.sync;
        for (j, side) in state.sides.iter_mut().enumerate() {
            if j == i || sync {
                motion(
                    &mut side.orbit,
                    rotation,
                    if j == i { pan } else { Vec2::ZERO },
                    zoom,
                );
            }
        }
    }
    for (i, action) in actions {
        state.message.clear();
        state.alignment_undo = None;
        let side = &mut state.sides[i];
        let result = (|| -> Result<(), String> {
            match action {
                0 => {
                    side.molecule = Some(mol.clone());
                    side.frames.clear();
                    side.title = "Current structure".into();
                    side.orbit = orbit.clone();
                }
                1 => {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("Structures", &["xyz", "pdb", "ent", "cif", "mmcif"])
                        .pick_file()
                    {
                        let (molecule, frames) = load(&path, &settings)?;
                        fit(&mut side.orbit, &molecule);
                        side.molecule = Some(molecule);
                        side.frames = frames;
                        side.frame = 0;
                        side.title = path
                            .file_name()
                            .unwrap_or_default()
                            .to_string_lossy()
                            .into();
                    }
                }
                2 => {
                    side.frames = traj.frames.clone();
                    side.frame = traj.current_frame.min(side.frames.len().saturating_sub(1));
                    side.title = "Viewer trajectory".into();
                    let f = &side.frames[side.frame];
                    let mut m = mol.clone();
                    if m.atoms != f.atoms {
                        m.topology = None;
                    }
                    m.atoms = f.atoms.clone();
                    m.pos = f.pos.clone();
                    m.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                    side.molecule = Some(m);
                }
                3 => {
                    if let Some(m) = &side.molecule {
                        fit(&mut side.orbit, m);
                    }
                }
                4 => {
                    let f = &side.frames[side.frame];
                    if let Some(m) = &mut side.molecule {
                        if m.atoms != f.atoms {
                            m.topology = None;
                        }
                        m.atoms = f.atoms.clone();
                        m.pos = f.pos.clone();
                        m.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                    }
                }
                _ => {}
            }
            side.dirty = true;
            Ok(())
        })();
        if let Err(e) = result {
            crate::runtime_warning::push(&ctx, "Comparison view", &e);
        }
    }
    if align {
        let result = (|| -> Result<(), String> {
            let a = state.sides[0]
                .molecule
                .as_ref()
                .ok_or("Load molecule A first.")?
                .clone();
            let b = state.sides[1]
                .molecule
                .as_ref()
                .ok_or("Load molecule B first.")?;
            if a.atoms != b.atoms {
                return Err("Alignment requires the same elements in the same atom order.".into());
            }
            let transform = crate::molecular_alignment::superpose(&b.pos, &a.pos)?;
            let positions = transform.apply_all(&b.pos);
            let rmsd = crate::rmsd::kabsch::rmsd(&positions, &a.pos)?;
            state.alignment_undo = Some((b.pos.clone(), state.sides[1].orbit.clone()));
            state.sides[1].molecule.as_mut().unwrap().pos = positions;
            state.sides[1].dirty = true;
            state.sides[1].orbit.target = state.sides[0].orbit.target;
            state.message = format!("Aligned RMSD: {rmsd:.6} Å (matching atom order)");
            Ok(())
        })();
        if let Err(e) = result {
            crate::runtime_warning::push(&ctx, "Comparison alignment", &e);
        }
    }
    if restore {
        if let Some((positions, camera)) = state.alignment_undo.take() {
            if let Some(m) = &mut state.sides[1].molecule {
                m.pos = positions;
                state.sides[1].orbit = camera;
                state.sides[1].dirty = true;
                state.message.clear();
            }
        }
    }
    for side in &state.sides {
        if let Some(scene) = &side.scene {
            commands
                .entity(scene.camera)
                .insert(crate::isolated_scene::camera_transform(&side.orbit));
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn comparison_motion_operates_on_an_independent_camera() {
        let main = OrbitCamera::default();
        let mut side = main.clone();
        motion(&mut side, Vec2::new(12., 5.), Vec2::ZERO, 0.8);
        assert_eq!(main.radius, 12.);
        assert!((side.radius - 9.6).abs() < 1e-5);
        assert_ne!(side.orientation, main.orientation);
        assert!((side.orientation.length() - 1.).abs() < 1e-6);
    }
    #[test]
    fn fitting_centers_without_modifying_geometry() {
        let molecule = Molecule::from_xyz("2\ntest\nH 10 0 0\nH 12 0 0\n");
        let old = molecule.pos.clone();
        let mut camera = OrbitCamera::default();
        fit(&mut camera, &molecule);
        assert_eq!(camera.target, Vec3::new(11., 0., 0.));
        assert_eq!(molecule.pos, old);
    }
}
