//! Independent off-screen animation capture, one encoded frame at a time.
use crate::{
    camera::OrbitCamera,
    isolated_scene::Scene,
    molecule::Molecule,
    settings::MolSettings,
    trajectory::{TrajectoryFrame, TrajectoryState},
    ui_style::ResponseExt,
};
use bevy::{
    prelude::*,
    render::view::screenshot::{Screenshot, ScreenshotCaptured},
    tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task},
};
use bevy_egui::{egui, EguiContexts, EguiPrimaryContextPass};
use std::{
    collections::BTreeSet,
    path::PathBuf,
    process::Child,
    sync::{Arc, Mutex},
};
#[path = "exporting/animation_encoder.rs"]
mod encoder;
use encoder::{Encoder, Format};
#[cfg(test)]
#[path = "exporting/animation_render_tests.rs"]
mod render_tests;
const OPEN: &str = "open_animation_export";
pub struct Plugin;
impl bevy::prelude::Plugin for Plugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<State>()
            .add_systems(EguiPrimaryContextPass, ui.after(crate::ui::ui_panel))
            .add_systems(Update, tick.after(crate::UpdateSet::Extras));
    }
}
#[derive(Clone, Copy, PartialEq, Eq)]
enum Source {
    Trajectory,
    Rotation,
}
#[derive(Resource)]
struct State {
    open: bool,
    format: Format,
    source: Source,
    width: u16,
    height: u16,
    fps: u16,
    first: usize,
    last: usize,
    stride: usize,
    rotation_frames: usize,
    repeat: bool,
    ffmpeg: String,
    request: Option<Request>,
    job: Option<Job>,
    message: Option<String>,
    error: Option<String>,
}
impl Default for State {
    fn default() -> Self {
        Self {
            open: false,
            format: Format::Gif,
            source: Source::Trajectory,
            width: 800,
            height: 600,
            fps: 20,
            first: 1,
            last: 1,
            stride: 1,
            rotation_frames: 120,
            repeat: true,
            ffmpeg: "ffmpeg".into(),
            request: None,
            job: None,
            message: None,
            error: None,
        }
    }
}
struct Request {
    frames: Vec<TrajectoryFrame>,
    count: usize,
    rotation: bool,
    template: Molecule,
    settings: MolSettings,
    hidden: BTreeSet<usize>,
    orbit: OrbitCamera,
    projection: Projection,
    path: PathBuf,
    fixed_bonds: bool,
}
struct Job {
    request: Request,
    scene: Scene,
    encoder: Arc<Mutex<Encoder>>,
    child: Option<Arc<Mutex<Child>>>,
    slot: Arc<Mutex<Option<Image>>>,
    capture: Option<Entity>,
    capture_started: Option<std::time::Instant>,
    cancel_flag: Arc<std::sync::atomic::AtomicBool>,
    worker: Option<Task<Result<(), String>>>,
    index: usize,
    settle: u8,
    finalizing: bool,
    cancel: bool,
}
pub fn button(ui: &mut egui::Ui) {
    if ui.button("Export animation…").clicked_once() {
        ui.ctx()
            .data_mut(|d| d.insert_temp(egui::Id::new(OPEN), true));
    }
}

fn ui(
    mut contexts: EguiContexts,
    mut state: ResMut<State>,
    traj: Res<TrajectoryState>,
    mol: Res<Molecule>,
    settings: Res<MolSettings>,
    orbit: Res<OrbitCamera>,
    builder: Res<crate::molecule_builder::builder_ui::ZMatrixBuilderState>,
    camera: Query<&Projection, With<crate::scene::MainCamera>>,
    mut regions: ResMut<crate::ui::UiPanelRegions>,
) {
    let Ok(ctx) = contexts.ctx_mut() else {
        return;
    };
    if ctx
        .data_mut(|d| d.remove_temp::<bool>(egui::Id::new(OPEN)))
        .unwrap_or(false)
    {
        state.open = true;
        state.last = traj.frames.len().max(1);
        state.repeat = matches!(traj.mode, crate::trajectory::PlaybackMode::Loop);
        if traj.frames.is_empty() {
            state.source = Source::Rotation;
        }
    }
    if let Some(error) = state.error.take() {
        crate::runtime_warning::push(ctx, "Animation export", &error);
    }
    if !state.open {
        return;
    }
    let mut open = true;
    let window = egui::Window::new("Export animation")
        .open(&mut open)
        .resizable(true)
        .default_size([450., 380.])
        .show(ctx, |ui| {
            crate::ui_style::modern(ui);
            let busy = state.job.is_some() || state.request.is_some();
            ui.add_enabled_ui(!busy, |ui| {
                ui.horizontal(|ui| {
                    ui.selectable_value(&mut state.format, Format::Gif, "GIF");
                    ui.selectable_value(&mut state.format, Format::Mp4, "MP4");
                });
                ui.horizontal_wrapped(|ui| {
                    ui.selectable_value(
                        &mut state.source,
                        Source::Trajectory,
                        "Trajectory / vibration",
                    );
                    ui.selectable_value(
                        &mut state.source,
                        Source::Rotation,
                        "Rotate current structure",
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Size (px)");
                    ui.add(egui::DragValue::new(&mut state.width).range(2..=4096));
                    ui.label("×");
                    ui.add(egui::DragValue::new(&mut state.height).range(2..=4096));
                });
                ui.horizontal(|ui| {
                    ui.label("Frames/s");
                    ui.add(egui::DragValue::new(&mut state.fps).range(1..=100));
                    if state.format == Format::Gif {
                        ui.checkbox(&mut state.repeat, "Loop GIF");
                    }
                });
                if state.source == Source::Trajectory {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Frames");
                        ui.add(
                            egui::DragValue::new(&mut state.first)
                                .range(1..=traj.frames.len().max(1)),
                        );
                        ui.label("to");
                        ui.add(
                            egui::DragValue::new(&mut state.last)
                                .range(1..=traj.frames.len().max(1)),
                        );
                        ui.label("Stride");
                        ui.add(egui::DragValue::new(&mut state.stride).range(1..=100000));
                    });
                } else {
                    ui.horizontal(|ui| {
                        ui.label("Frames per turn");
                        ui.add(egui::DragValue::new(&mut state.rotation_frames).range(2..=10000));
                    });
                }
                if state.format == Format::Mp4 {
                    ui.label("FFmpeg executable");
                    ui.horizontal(|ui| {
                        ui.text_edit_singleline(&mut state.ffmpeg);
                        if ui.button("Browse…").clicked_once() {
                            if let Some(p) = rfd::FileDialog::new().pick_file() {
                                state.ffmpeg = p.display().to_string();
                            }
                        }
                    });
                }
                if ui
                    .add_enabled(!mol.atoms.is_empty(), egui::Button::new("Export…"))
                    .clicked_once()
                {
                    let result = (|| -> Result<(), String> {
                        let rotation = state.source == Source::Rotation;
                        let frames = if rotation {
                            vec![TrajectoryFrame {
                                metadata: Default::default(),
                                atoms: mol.atoms.clone(),
                                pos: mol.pos.clone(),
                            }]
                        } else {
                            selected_frames(&traj.frames, state.first, state.last, state.stride)?
                        };
                        if frames.iter().any(|f| {
                            f.atoms.len() != f.pos.len() || f.pos.iter().any(|p| !p.is_finite())
                        }) {
                            return Err(
                                "Animation frames must contain finite coordinates for every atom."
                                    .into(),
                            );
                        }
                        let projection = camera.single().map_err(|e| e.to_string())?.clone();
                        let Some(mut path) = rfd::FileDialog::new()
                            .add_filter("Animation", &[state.format.ext()])
                            .set_file_name(format!("beavyr-animation.{}", state.format.ext()))
                            .save_file()
                        else {
                            return Ok(());
                        };
                        path.set_extension(state.format.ext());
                        state.request = Some(Request {
                            count: if rotation {
                                state.rotation_frames
                            } else {
                                frames.len()
                            },
                            frames,
                            rotation,
                            template: mol.clone(),
                            settings: settings.clone(),
                            hidden: builder.selection.hidden.clone(),
                            orbit: orbit.clone(),
                            projection,
                            path,
                            fixed_bonds: traj.fixed_bonds,
                        });
                        state.message = None;
                        Ok(())
                    })();
                    if let Err(e) = result {
                        state.error = Some(e);
                    }
                }
            });
            if let Some(job) = &mut state.job {
                ui.add(
                    egui::ProgressBar::new(job.index as f32 / job.request.count as f32)
                        .text(format!("{} / {} frames", job.index, job.request.count)),
                );
                if ui.button("Cancel export").clicked_once() {
                    job.cancel = true;
                }
                ctx.request_repaint();
            }
            if let Some(message) = &state.message {
                ui.label(message);
            }
        });
    state.open = open;
    if let Some(window) = window {
        let r = window.response.rect;
        regions
            .windows
            .push(Rect::new(r.min.x, r.min.y, r.max.x, r.max.y));
    }
}
fn selected_frames(
    frames: &[TrajectoryFrame],
    first: usize,
    last: usize,
    stride: usize,
) -> Result<Vec<TrajectoryFrame>, String> {
    if first == 0 || last < first || last > frames.len() || stride == 0 {
        return Err("Choose a valid frame range from the loaded trajectory.".into());
    }
    Ok((first - 1..last)
        .step_by(stride)
        .map(|i| frames[i].clone())
        .collect())
}
fn draw(job: &mut Job, commands: &mut Commands, materials: &mut Assets<StandardMaterial>) {
    let frame = &job.request.frames[if job.request.rotation { 0 } else { job.index }];
    let mut mol = job.request.template.clone();
    if mol.atoms != frame.atoms {
        mol.topology = None;
    }
    mol.atoms = frame.atoms.clone();
    mol.pos = frame.pos.clone();
    if !job.request.fixed_bonds || mol.atoms != job.request.template.atoms {
        mol.recompute_bonds(
            job.request.settings.bond_thresh_scale,
            job.request.settings.hbond_cutoff,
        );
    }
    job.scene.molecule(
        commands,
        materials,
        &mol,
        &job.request.settings,
        &job.request.hidden,
    );
    if job.request.rotation {
        let mut orbit = job.request.orbit.clone();
        orbit.orientation = Quat::from_rotation_y(
            std::f32::consts::TAU * job.index as f32 / job.request.count as f32,
        ) * orbit.orientation;
        commands
            .entity(job.scene.camera)
            .insert(crate::isolated_scene::camera_transform(&orbit));
    }
    job.settle = 3;
}
fn tick(
    mut commands: Commands,
    mut state: ResMut<State>,
    mut images: ResMut<Assets<Image>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    if let Some(request) = state.request.take() {
        match Encoder::new(
            &request.path,
            state.format,
            state.width,
            state.height,
            state.fps,
            state.repeat,
            &state.ffmpeg,
        ) {
            Err(e) => state.error = Some(e),
            Ok(encoder) => {
                let scene = Scene::new(
                    &mut commands,
                    &mut images,
                    &mut meshes,
                    &request.settings,
                    &request.orbit,
                    &request.projection,
                    UVec2::new(state.width as u32, state.height as u32),
                    20,
                );
                let child = encoder.child.clone();
                let cancel_flag = encoder.cancel.clone();
                let mut job = Job {
                    request,
                    scene,
                    encoder: Arc::new(Mutex::new(encoder)),
                    child,
                    slot: Default::default(),
                    capture: None,
                    capture_started: None,
                    cancel_flag,
                    worker: None,
                    index: 0,
                    settle: 0,
                    finalizing: false,
                    cancel: false,
                };
                draw(&mut job, &mut commands, &mut materials);
                state.job = Some(job);
            }
        }
    }
    let Some(mut job) = state.job.take() else {
        return;
    };
    let mut error = None;
    let mut done = false;
    if job.cancel {
        job.cancel_flag
            .store(true, std::sync::atomic::Ordering::SeqCst);
        if let Some(child) = &job.child {
            if let Ok(mut child) = child.lock() {
                let _ = child.kill();
            }
        }
    }
    if let Some(worker) = &mut job.worker {
        if let Some(result) = block_on(future::poll_once(worker)) {
            job.worker = None;
            match result {
                Err(e) => error = Some(e),
                Ok(()) => {
                    if job.finalizing {
                        done = true;
                    } else {
                        job.index += 1;
                        if job.index < job.request.count && !job.cancel {
                            draw(&mut job, &mut commands, &mut materials);
                        }
                    }
                }
            }
        }
    }
    if job
        .capture_started
        .is_some_and(|t| t.elapsed() > std::time::Duration::from_secs(60))
        && job.capture.is_some()
    {
        error = Some("Timed out waiting for the graphics renderer to capture a frame.".into());
    }
    if job.worker.is_none() && !job.cancel && !done && error.is_none() {
        if job.index == job.request.count {
            let encoder = job.encoder.clone();
            job.finalizing = true;
            job.worker = Some(
                AsyncComputeTaskPool::get()
                    .spawn(async move { encoder.lock().map_err(|e| e.to_string())?.finish() }),
            );
        } else {
            let image = job.slot.lock().ok().and_then(|mut slot| slot.take());
            if let Some(image) = image {
                if let Some(entity) = job.capture.take() {
                    commands.entity(entity).try_despawn();
                }
                let encoder = job.encoder.clone();
                job.worker = Some(AsyncComputeTaskPool::get().spawn(async move {
                    let rgba = image
                        .try_into_dynamic()
                        .map_err(|e| e.to_string())?
                        .to_rgba8();
                    encoder
                        .lock()
                        .map_err(|e| e.to_string())?
                        .frame(rgba.as_raw())
                }));
            } else if job.capture.is_none() {
                if job.settle > 0 {
                    job.settle -= 1;
                } else {
                    job.capture_started = Some(std::time::Instant::now());
                    let slot = job.slot.clone();
                    job.capture = Some(
                        commands
                            .spawn(Screenshot::image(job.scene.image.clone()))
                            .observe(move |event: On<ScreenshotCaptured>| {
                                if let Ok(mut slot) = slot.lock() {
                                    *slot = Some(event.image.clone());
                                }
                            })
                            .id(),
                    );
                }
            }
        }
    }
    if (job.cancel && job.worker.is_none()) || done || error.is_some() {
        if let Some(entity) = job.capture.take() {
            commands.entity(entity).try_despawn();
        }
        state.message = Some(if done {
            format!("Saved {}", job.request.path.display())
        } else if job.cancel {
            "Export cancelled.".into()
        } else {
            "Export failed.".into()
        });
        if !job.cancel {
            state.error = error;
        }
        job.scene.dispose(&mut commands, &mut images);
    } else {
        state.job = Some(job);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn frame_ranges_are_inclusive_and_stride_keeps_order() {
        let f: Vec<_> = (0..5)
            .map(|i| TrajectoryFrame {
                metadata: Default::default(),
                atoms: vec!["C".into()],
                pos: vec![Vec3::new(i as f32, 0., 0.)],
            })
            .collect();
        let out = selected_frames(&f, 2, 5, 2).unwrap();
        assert_eq!(
            out.iter().map(|f| f.pos[0].x).collect::<Vec<_>>(),
            vec![1., 3.]
        );
        assert!(selected_frames(&f, 0, 5, 1).is_err());
        assert!(selected_frames(&f, 1, 5, 0).is_err());
    }
}
