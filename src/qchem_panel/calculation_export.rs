//! Copies a completed run into a new directory; sources and previous exports are untouched.
use crate::ui_style::ResponseExt;
use bevy::tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task};
use bevy_egui::egui;
use std::{
    fs,
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicU64, Ordering},
        Arc,
    },
    time::{SystemTime, UNIX_EPOCH},
};
const REQUEST: &str = "qc_complete_export_request";
const BUSY: &str = "qc_complete_export_busy";
#[derive(Default)]
pub struct State {
    task: Option<Task<Result<Saved, String>>>,
    pub open: bool,
    bytes: Arc<AtomicU64>,
    result: Option<Result<Saved, String>>,
}
#[derive(Clone, Debug)]
pub struct Saved {
    pub directory: PathBuf,
    pub files: usize,
    pub bytes: u64,
}
#[derive(Clone, Default)]
struct Request {
    source: PathBuf,
    parent: PathBuf,
}
pub fn button(ui: &mut egui::Ui, root: Option<&Path>, enabled: bool) {
    let busy = ui
        .ctx()
        .data(|d| d.get_temp::<bool>(egui::Id::new(BUSY)))
        .unwrap_or(false);
    if ui
        .add_enabled(
            enabled && root.is_some() && !busy,
            egui::Button::new("Save complete calculation…"),
        )
        .clicked_once()
    {
        if let (Some(source), Some(parent)) = (
            root,
            rfd::FileDialog::new()
                .set_title("Choose a folder for the complete calculation")
                .pick_folder(),
        ) {
            ui.ctx().data_mut(|d| {
                d.insert_temp(
                    egui::Id::new(REQUEST),
                    Request {
                        source: source.to_owned(),
                        parent,
                    },
                )
            });
        }
    }
}
fn copy_tree(
    source: &Path,
    target: &Path,
    relative: &Path,
    inventory: &mut Vec<(PathBuf, u64)>,
    bytes: &AtomicU64,
    depth: usize,
) -> Result<(), String> {
    if depth > 64 {
        return Err("Calculation directory nesting exceeds 64 levels".into());
    }
    let mut entries = fs::read_dir(source)
        .map_err(|e| format!("Cannot read {}: {e}", source.display()))?
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| e.to_string())?;
    entries.sort_by_key(|e| e.file_name());
    for entry in entries {
        let path = entry.path();
        let dest = target.join(entry.file_name());
        let relative = relative.join(entry.file_name());
        let kind = entry.file_type().map_err(|e| e.to_string())?;
        if kind.is_symlink() {
            return Err(format!("Cannot create a self-contained export from symbolic link {}. Original calculation files are unchanged.",path.display()));
        }
        if kind.is_dir() {
            fs::create_dir(&dest).map_err(|e| e.to_string())?;
            copy_tree(&path, &dest, &relative, inventory, bytes, depth + 1)?;
        } else if kind.is_file() {
            let before = entry.metadata().map_err(|e| e.to_string())?;
            let copied = fs::copy(&path, &dest)
                .map_err(|e| format!("Cannot copy {}: {e}", path.display()))?;
            let after = fs::metadata(&path).map_err(|e| e.to_string())?;
            if copied != before.len()
                || after.len() != before.len()
                || after.modified().ok() != before.modified().ok()
            {
                return Err(format!("{} changed while saving. Wait until all calculation processes finish and retry.",path.display()));
            }
            bytes.fetch_add(copied, Ordering::Relaxed);
            inventory.push((relative, copied));
        } else {
            return Err(format!("Unsupported calculation file: {}", path.display()));
        }
    }
    Ok(())
}
pub fn save(source: &Path, parent: &Path, bytes: &AtomicU64) -> Result<Saved, String> {
    let source = source
        .canonicalize()
        .map_err(|e| format!("Cannot find calculation files: {e}"))?;
    let parent = parent
        .canonicalize()
        .map_err(|e| format!("Cannot open destination: {e}"))?;
    if !source.is_dir() || !parent.is_dir() {
        return Err("Calculation and destination must be directories".into());
    }
    if parent.starts_with(&source) {
        return Err("Choose a destination outside the calculation directory.".into());
    }
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    let mut attempt = 0;
    let target = loop {
        let name = format!(
            "Beavyr-calculation-{stamp}-{}-{attempt}",
            std::process::id()
        );
        let path = parent.join(name);
        match fs::create_dir(&path) {
            Ok(()) => break path,
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => attempt += 1,
            Err(e) => return Err(format!("Cannot create destination: {e}")),
        }
    };
    let mut inventory = vec![];
    let copied = (|| {
        copy_tree(&source, &target, Path::new(""), &mut inventory, bytes, 0)?;
        let total = inventory.iter().map(|(_, size)| size).sum();
        let count = inventory.len();
        let manifest = ron::ser::to_string_pretty(&inventory, ron::ser::PrettyConfig::default())
            .map_err(|e| e.to_string())?;
        write_unique(
            &target,
            "beavyr_export_inventory",
            "ron",
            manifest.as_bytes(),
        )
        .map_err(|e| format!("Cannot write export inventory: {e}"))?;
        Ok(Saved {
            directory: target.clone(),
            files: count,
            bytes: total,
        })
    })();
    copied.map_err(|error: String| {
        // Never overwrite even an unusually named native artifact in a partial copy.
        let _ = write_unique(&target, "INCOMPLETE_EXPORT", "txt", error.as_bytes());
        format!(
            "{error}\nThe incomplete copy is in {}. No source files were changed.",
            target.display()
        )
    })
}

fn write_unique(dir: &Path, stem: &str, extension: &str, content: &[u8]) -> std::io::Result<()> {
    use std::io::Write;
    for number in 0u64.. {
        match fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(dir.join(format!("{stem}_{number}.{extension}")))
        {
            Ok(mut file) => return file.write_all(content),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e),
        }
    }
    unreachable!()
}
pub fn window(ctx: &egui::Context, state: &mut State) -> Option<egui::Rect> {
    if let Some(request) = ctx.data_mut(|d| d.remove_temp::<Request>(egui::Id::new(REQUEST))) {
        if state.task.is_none() {
            state.bytes = Default::default();
            let bytes = state.bytes.clone();
            state.result = None;
            state.open = true;
            state.task = Some(
                AsyncComputeTaskPool::get()
                    .spawn(async move { save(&request.source, &request.parent, &bytes) }),
            );
        }
    }
    if let Some(task) = &mut state.task {
        if let Some(result) = block_on(future::poll_once(task)) {
            if let Err(e) = &result {
                crate::runtime_warning::push(ctx, "Save complete calculation", e);
            }
            state.result = Some(result);
            state.task = None;
        } else {
            ctx.request_repaint_after(std::time::Duration::from_millis(100));
        }
    }
    ctx.data_mut(|d| d.insert_temp(egui::Id::new(BUSY), state.task.is_some()));
    if !state.open {
        return None;
    }
    let mut open = true;
    let response = egui::Window::new("Save complete calculation")
        .id(egui::Id::new("qc_complete_export_window"))
        .open(&mut open)
        .resizable(true)
        .default_size([480., 180.])
        .show(ctx, |ui| {
            if state.task.is_some() {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(format!(
                        "Saving… {:.1} MiB copied",
                        state.bytes.load(Ordering::Relaxed) as f64 / 1048576.
                    ));
                });
            }
            if let Some(result) = &state.result {
                match result {
                    Ok(saved) => {
                        ui.label(format!(
                            "Saved {} files ({:.1} MiB)",
                            saved.files,
                            saved.bytes as f64 / 1048576.
                        ));
                        ui.add(
                            egui::Label::new(saved.directory.display().to_string())
                                .selectable(true),
                        );
                        if ui.button("Copy folder path").clicked_once() {
                            ctx.copy_text(saved.directory.display().to_string());
                        }
                    }
                    Err(e) => {
                        ui.label("The calculation could not be completely saved.");
                        if ui.button("Show error").clicked_once() {
                            crate::runtime_warning::push(ctx, "Save complete calculation", e);
                        }
                    }
                }
            }
        });
    state.open = open;
    response.map(|r| r.response.rect)
}
