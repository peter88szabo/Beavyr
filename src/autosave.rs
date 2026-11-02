//! Per-process recovery generations. Manual project paths are never autosave targets.
use super::{ProjectState, Snapshot, Workspace};
use crate::ui_style::ResponseExt;
use bevy::{
    prelude::*,
    tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task},
};
use bevy_egui::{egui, EguiContexts, EguiPrimaryContextPass};
use serde::{Deserialize, Serialize};
use std::{
    fs::{self, File},
    path::{Path, PathBuf},
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};

pub struct Plugin;
impl bevy::prelude::Plugin for Plugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<State>()
            .add_systems(
                Update,
                process.after(super::process).before(crate::UpdateSet::Core),
            )
            .add_systems(Update, capture.after(crate::UpdateSet::Extras))
            .add_systems(EguiPrimaryContextPass, ui.after(crate::ui::ui_panel));
    }
}
#[derive(Clone, Serialize, Deserialize)]
struct Settings {
    enabled: bool,
    seconds: u64,
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            enabled: true,
            seconds: 60,
        }
    }
}
struct Session {
    directory: PathBuf,
    _lock: File,
}
#[derive(Clone)]
struct Entry {
    path: PathBuf,
    modified: SystemTime,
    bytes: u64,
}
enum Completed {
    Initialized(Session, Vec<Entry>, Settings),
    Saved(PathBuf),
    Refreshed(Vec<Entry>),
}
#[derive(Resource)]
struct State {
    initialized: bool,
    root: Option<PathBuf>,
    session: Option<Session>,
    task: Option<Task<Result<Completed, String>>>,
    settings: Settings,
    entries: Vec<Entry>,
    open: bool,
    deadline: Instant,
    capture: bool,
    workspace: Option<Workspace>,
    recover: Option<(PathBuf, Workspace)>,
    message: String,
    error: Option<String>,
    generation: u64,
    last_saved: Option<PathBuf>,
}
impl Default for State {
    fn default() -> Self {
        Self {
            initialized: false,
            root: None,
            session: None,
            task: None,
            settings: Default::default(),
            entries: vec![],
            open: false,
            deadline: Instant::now() + Duration::from_secs(60),
            capture: false,
            workspace: None,
            recover: None,
            message: String::new(),
            error: None,
            generation: 0,
            last_saved: None,
        }
    }
}
fn root() -> Result<PathBuf, String> {
    #[cfg(target_os = "windows")]
    let base = std::env::var_os("LOCALAPPDATA")
        .map(PathBuf::from)
        .map(|p| p.join("Beavyr"));
    #[cfg(target_os = "macos")]
    let base = std::env::var_os("HOME")
        .map(|p| PathBuf::from(p).join("Library/Application Support/Beavyr"));
    #[cfg(not(any(target_os = "windows", target_os = "macos")))]
    let base = std::env::var_os("XDG_STATE_HOME")
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|p| PathBuf::from(p).join(".local/state")))
        .map(|p| p.join("beavyr"));
    base.map(|p| p.join("recovery"))
        .ok_or("Cannot determine the recovery directory.".into())
}
fn settings(root: &Path) -> Settings {
    let mut s = fs::read_to_string(root.join("settings.ron"))
        .ok()
        .and_then(|s| ron::from_str::<Settings>(&s).ok())
        .unwrap_or_default();
    s.seconds = s.seconds.clamp(10, 3600);
    s
}
fn initialize(root: &Path) -> Result<Completed, String> {
    fs::create_dir_all(root).map_err(|e| e.to_string())?;
    let entries = discover(root)?;
    let id = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    let directory = root.join(format!("session-{id}-{}", std::process::id()));
    fs::create_dir(&directory).map_err(|e| e.to_string())?;
    let lock = File::options()
        .create_new(true)
        .read(true)
        .write(true)
        .open(directory.join("active.lock"))
        .map_err(|e| e.to_string())?;
    lock.try_lock().map_err(|e| e.to_string())?;
    Ok(Completed::Initialized(
        Session {
            directory,
            _lock: lock,
        },
        entries,
        settings(root),
    ))
}
fn discover(root: &Path) -> Result<Vec<Entry>, String> {
    let mut entries = vec![];
    for item in fs::read_dir(root).map_err(|e| e.to_string())? {
        let item = item.map_err(|e| e.to_string())?;
        if !item.file_type().map_err(|e| e.to_string())?.is_dir()
            || !item.file_name().to_string_lossy().starts_with("session-")
        {
            continue;
        }
        let Ok(lock) = File::options()
            .read(true)
            .write(true)
            .open(item.path().join("active.lock"))
        else {
            continue;
        };
        if lock.try_lock().is_err() {
            continue;
        } // A different Beavyr still owns this session.
        for file in fs::read_dir(item.path()).map_err(|e| e.to_string())? {
            let file = file.map_err(|e| e.to_string())?;
            if !file.file_type().map_err(|e| e.to_string())?.is_file()
                || file.path().extension().and_then(|e| e.to_str()) != Some("beavyr")
            {
                continue;
            }
            let m = file.metadata().map_err(|e| e.to_string())?;
            entries.push(Entry {
                path: file.path(),
                modified: m.modified().unwrap_or(UNIX_EPOCH),
                bytes: m.len(),
            });
        }
    }
    entries.sort_by(|a, b| b.modified.cmp(&a.modified));
    Ok(entries)
}
fn save(directory: &Path, generation: u64, snapshot: Snapshot) -> Result<Completed, String> {
    snapshot.validate().map_err(|e| e.to_string())?;
    let text = ron::ser::to_string(&snapshot).map_err(|e| e.to_string())?;
    let path = directory.join(format!("recovery-{generation:08}.beavyr"));
    super::atomic_write(&path, text.as_bytes()).map_err(|e| e.to_string())?;
    // Rotate only this process's recovery generations, after the new file is durable.
    if generation > 3 {
        let _ = fs::remove_file(directory.join(format!("recovery-{:08}.beavyr", generation - 3)));
    }
    Ok(Completed::Saved(path))
}
fn recover(world: &mut World, path: &Path, workspace: Workspace) -> Result<(), String> {
    if super::busy(world) {
        return Err("Stop running calculations before restoring a session.".into());
    }
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let snapshot: Snapshot =
        ron::from_str(&text).map_err(|e| format!("Invalid recovery snapshot: {e}"))?;
    let orbitals = snapshot.validate().map_err(|e| e.to_string())?;
    let previous = Snapshot::capture(world, workspace).map_err(|e| e.to_string())?;
    snapshot.apply(world, orbitals);
    let mut project = world.resource_mut::<ProjectState>();
    project.previous = Some(previous);
    project.path = None; // Save project must never silently overwrite a recovery generation.
    project.message =
        Some("Session recovered. Interrupted calculations were not restarted.".into());
    project.error = false;
    Ok(())
}
fn process(world: &mut World) {
    world.resource_scope(|world, mut state: Mut<State>| {
        if !state.initialized {
            state.initialized = true;
            match root() {
                Ok(root) => {
                    state.root = Some(root.clone());
                    state.task =
                        Some(AsyncComputeTaskPool::get().spawn(async move { initialize(&root) }));
                }
                Err(e) => state.error = Some(e),
            }
        }
        if let Some(task) = &mut state.task {
            if let Some(result) = block_on(future::poll_once(task)) {
                state.task = None;
                match result {
                    Ok(Completed::Initialized(session, entries, settings)) => {
                        state.session = Some(session);
                        state.settings = settings;
                        state.open = !entries.is_empty();
                        state.entries = entries;
                        state.deadline =
                            Instant::now() + Duration::from_secs(state.settings.seconds);
                    }
                    Ok(Completed::Saved(path)) => {
                        state.message = "Recovery snapshot saved.".into();
                        state.last_saved = Some(path);
                    }
                    Ok(Completed::Refreshed(entries)) => state.entries = entries,
                    Err(e) => state.error = Some(e),
                }
            }
        }
        if let Some((path, workspace)) = state.recover.take() {
            match recover(world, &path, workspace) {
                Ok(()) => {
                    state.message =
                        "Session recovered. Use Project → Restore previous session to undo.".into();
                    state.open = false;
                }
                Err(e) => state.error = Some(e),
            }
        }
        if state.settings.enabled
            && Instant::now() >= state.deadline
            && state.task.is_none()
            && state.session.is_some()
        {
            state.capture = true;
        }
    });
}
fn capture(world: &mut World) {
    world.resource_scope(|world, mut state: Mut<State>| {
        if let Some(workspace) = state.workspace.take() {
            state.capture = false;
            state.deadline = Instant::now() + Duration::from_secs(state.settings.seconds);
            if state.task.is_some() {
                return;
            }
            let Some(directory) = state.session.as_ref().map(|s| s.directory.clone()) else {
                return;
            };
            if world
                .resource::<crate::molecule::Molecule>()
                .atoms
                .is_empty()
                && world
                    .resource::<crate::ui::XyzBuffer>()
                    .text
                    .trim()
                    .is_empty()
            {
                return;
            }
            match Snapshot::capture(world, workspace) {
                Ok(snapshot) => {
                    state.generation += 1;
                    let generation = state.generation;
                    state.task = Some(
                        AsyncComputeTaskPool::get()
                            .spawn(async move { save(&directory, generation, snapshot) }),
                    );
                }
                Err(e) => state.error = Some(e.to_string()),
            }
        }
    });
}
fn ui(
    mut contexts: EguiContexts,
    mut state: ResMut<State>,
    mut regions: ResMut<crate::ui::UiPanelRegions>,
) {
    let Ok(ctx) = contexts.ctx_mut() else {
        return;
    };
    if state.capture {
        state.workspace = Some(Workspace::capture(ctx));
    }
    if ctx
        .data_mut(|d| d.remove_temp::<bool>(egui::Id::new("open_autosave")))
        .unwrap_or(false)
    {
        state.open = true;
    }
    if let Some(error) = state.error.take() {
        state.message = format!("Autosave/recovery failed: {error}");
        crate::runtime_warning::push(ctx, "Autosave / recovery", &state.message);
    }
    if !state.open {
        return;
    }
    let mut open = true;
    let window=egui::Window::new("Autosave / recovery").open(&mut open).resizable(true).default_size([600.,400.]).show(ctx,|ui| {
        crate::ui_style::modern(ui);
        let mut changed=ui.checkbox(&mut state.settings.enabled,"Autosave").changed();
        ui.horizontal(|ui|{ui.label("Interval (seconds)");changed|=ui.add(egui::DragValue::new(&mut state.settings.seconds).range(10..=3600)).changed();});
        if changed {
            state.deadline=Instant::now()+Duration::from_secs(state.settings.seconds);
            if let Some(root)=&state.root {
                let result=ron::ser::to_string(&state.settings).map_err(|e|e.to_string()).and_then(|s|super::atomic_write(&root.join("settings.ron"),s.as_bytes()).map_err(|e|e.to_string()));
                if let Err(e)=result {state.error=Some(e);}
            }
        }
        ui.horizontal(|ui|{
            if ui.add_enabled(state.task.is_none()&&state.session.is_some(),egui::Button::new("Save recovery now")).clicked_once(){state.capture=true;}
            if ui.add_enabled(state.task.is_none(),egui::Button::new("Refresh")).clicked_once(){if let Some(root)=state.root.clone(){state.task=Some(AsyncComputeTaskPool::get().spawn(async move{discover(&root).map(Completed::Refreshed)}));}}
        });
        ui.label("Recovery copies are separate from saved projects. Running calculations are not restarted.");
        if !state.message.is_empty(){ui.label(&state.message);}
        if state.task.is_some(){ui.spinner();ctx.request_repaint();}
        if let Some(path)=&state.last_saved {ui.label(format!("Latest: {}",path.display()));}
        ui.separator();
        let mut restore=None;
        egui::ScrollArea::vertical().show(ui,|ui|{
            if state.entries.is_empty(){ui.label("No previous sessions available.");}
            for e in &state.entries {
                ui.push_id(&e.path,|ui|{
                    ui.horizontal(|ui|{
                        let age=SystemTime::now().duration_since(e.modified).unwrap_or_default().as_secs()/60;
                        ui.label(format!("{age} min ago · {:.1} MB",e.bytes as f64/1_000_000.));
                        if ui.button("Restore").clicked_once(){restore=Some(e.path.clone());}
                        if ui.button("Copy path").clicked_once(){ctx.copy_text(e.path.display().to_string());}
                    });
                    ui.small(e.path.parent().and_then(|p|p.file_name()).unwrap_or_default().to_string_lossy());
                });
            }
        });
        if let Some(path)=restore{state.recover=Some((path,Workspace::capture(ctx)));}
    });
    state.open = open;
    if let Some(w) = window {
        let r = w.response.rect;
        regions
            .windows
            .push(Rect::new(r.min.x, r.min.y, r.max.x, r.max.y));
    }
}

#[cfg(test)]
#[path = "autosave_tests.rs"]
mod tests;
