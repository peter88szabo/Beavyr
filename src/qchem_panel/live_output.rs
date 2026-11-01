//! Background, bounded reads of calculation output. No process command is modified.
use crate::ui_style::ResponseExt;
use bevy::tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task};
use bevy_egui::egui;
use std::{
    fs,
    io::{Read, Seek, SeekFrom},
    path::{Path, PathBuf},
    time::{Duration, Instant, SystemTime},
};
const TAIL_BYTES: u64 = 512 * 1024;
#[derive(Clone, Debug)]
pub struct LogFile {
    pub relative: PathBuf,
    pub bytes: u64,
    pub modified: Option<SystemTime>,
}
#[derive(Default)]
pub struct Snapshot {
    pub files: Vec<LogFile>,
    pub selected: Option<PathBuf>,
    pub text: String,
    pub truncated: bool,
    pub omitted_files: usize,
}
pub struct State {
    pub open: bool,
    pub root: Option<PathBuf>,
    choice: Option<PathBuf>,
    snapshot: Snapshot,
    task: Option<Task<Result<Snapshot, String>>>,
    last_read: Option<Instant>,
    error: Option<String>,
    follow: bool,
    search: String,
    matching_only: bool,
}
impl Default for State {
    fn default() -> Self {
        Self {
            open: false,
            root: None,
            choice: None,
            snapshot: Default::default(),
            task: None,
            last_read: None,
            error: None,
            follow: true,
            search: String::new(),
            matching_only: false,
        }
    }
}
impl State {
    pub fn start(&mut self, root: PathBuf) {
        self.root = Some(root);
        self.choice = None;
        self.snapshot = Default::default();
        self.task = None;
        self.last_read = None;
        self.error = None;
        self.follow = true;
    }
}
fn is_log(path: &Path) -> bool {
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    matches!(
        path.extension().and_then(|s| s.to_str()),
        Some("out" | "err" | "log" | "stdout" | "stderr")
    ) || matches!(
        name.as_ref(),
        "status.txt"
            | "error.txt"
            | "summary.txt"
            | "scan_progress.txt"
            | "search_summary.txt"
            | "mecp_summary.txt"
    )
}
fn discover(
    root: &Path,
    relative: &Path,
    files: &mut Vec<LogFile>,
    depth: usize,
) -> Result<(), String> {
    if depth > 32 {
        return Ok(());
    }
    let here = root.join(relative);
    for entry in fs::read_dir(&here)
        .map_err(|e| format!("Cannot read output directory {}: {e}", here.display()))?
    {
        let entry = entry.map_err(|e| e.to_string())?;
        let kind = entry.file_type().map_err(|e| e.to_string())?;
        let path = relative.join(entry.file_name());
        if kind.is_dir() {
            discover(root, &path, files, depth + 1)?;
        } else if kind.is_file() && is_log(&path) {
            let metadata = entry.metadata().map_err(|e| e.to_string())?;
            files.push(LogFile {
                relative: path,
                bytes: metadata.len(),
                modified: metadata.modified().ok(),
            });
        }
    }
    Ok(())
}
pub fn read(root: &Path, choice: Option<&Path>) -> Result<Snapshot, String> {
    let mut files = vec![];
    discover(root, Path::new(""), &mut files, 0)?;
    files.sort_by(|a, b| {
        b.modified
            .cmp(&a.modified)
            .then_with(|| a.relative.cmp(&b.relative))
    });
    let selected = choice.map(Path::to_owned).or_else(|| {
        files
            .iter()
            .find(|f| f.bytes > 0 && !f.relative.starts_with(super::calculation_record::RECORD))
            .or_else(|| files.iter().find(|f| f.bytes > 0))
            .map(|f| f.relative.clone())
    });
    let (mut text, mut truncated) = (String::new(), false);
    if let Some(relative) = &selected {
        if relative.is_absolute()
            || relative
                .components()
                .any(|c| !matches!(c, std::path::Component::Normal(_)))
        {
            return Err("Output path must stay within the calculation directory".into());
        }
        let path = root.join(relative);
        let canonical = path.canonicalize().map_err(|e| e.to_string())?;
        if !canonical.starts_with(root.canonicalize().map_err(|e| e.to_string())?) {
            return Err("Output file is outside the calculation directory".into());
        }
        let mut file =
            fs::File::open(&path).map_err(|e| format!("Cannot read {}: {e}", path.display()))?;
        let len = file.metadata().map_err(|e| e.to_string())?.len();
        let offset = len.saturating_sub(TAIL_BYTES);
        truncated = offset > 0;
        file.seek(SeekFrom::Start(offset))
            .map_err(|e| e.to_string())?;
        let mut bytes = vec![];
        file.take(TAIL_BYTES)
            .read_to_end(&mut bytes)
            .map_err(|e| e.to_string())?;
        let start = if truncated {
            bytes
                .iter()
                .take(3)
                .take_while(|b| (**b & 0xc0) == 0x80)
                .count()
        } else {
            0
        };
        text = String::from_utf8_lossy(&bytes[start..]).into_owned();
    }
    let omitted_files = files.len().saturating_sub(200);
    if files.len() > 200 {
        let chosen = selected
            .as_ref()
            .and_then(|p| files.iter().find(|f| &f.relative == p))
            .cloned();
        files.truncate(200);
        if let Some(f) = chosen {
            if !files.iter().any(|p| p.relative == f.relative) {
                files.push(f);
            }
        }
    }
    Ok(Snapshot {
        files,
        selected,
        text,
        truncated,
        omitted_files,
    })
}
pub fn window(ctx: &egui::Context, state: &mut State, running: bool) -> Option<egui::Rect> {
    if !state.open {
        return None;
    }
    if let Some(task) = &mut state.task {
        if let Some(result) = block_on(future::poll_once(task)) {
            match result {
                Ok(snapshot) => {
                    state.snapshot = snapshot;
                    state.error = None;
                }
                Err(e) => state.error = Some(e),
            }
            state.task = None;
            state.last_read = Some(Instant::now());
        }
    }
    let refresh = state
        .last_read
        .is_none_or(|t| t.elapsed() >= Duration::from_millis(if running { 350 } else { 1500 }));
    if state.task.is_none() && refresh {
        if let Some(root) = state.root.clone() {
            let choice = state.choice.clone();
            state.task = Some(
                AsyncComputeTaskPool::get().spawn(async move { read(&root, choice.as_deref()) }),
            );
        }
    }
    ctx.request_repaint_after(Duration::from_millis(150));
    let mut open = true;
    let response = egui::Window::new("Live output")
        .id(egui::Id::new("qc_live_output_window"))
        .frame(crate::ui_style::window_frame(ctx))
        .open(&mut open)
        .resizable(true)
        .default_size([850., 550.])
        .min_size([360., 220.])
        .show(ctx, |ui| {
            crate::ui_style::modern(ui);
            ui.horizontal_wrapped(|ui| {
                if running {
                    ui.spinner();
                    ui.label("Running");
                } else {
                    ui.label("Finished");
                }
                ui.checkbox(&mut state.follow, "Follow output");
                if ui.button("Copy displayed output").clicked_once() {
                    let text = if state.matching_only && !state.search.is_empty() {
                        let query = state.search.to_lowercase();
                        state.snapshot.text.lines()
                            .filter(|line| line.to_lowercase().contains(&query))
                            .collect::<Vec<_>>().join("\n")
                    } else {
                        state.snapshot.text.clone()
                    };
                    ctx.copy_text(text);
                }
            });
            ui.horizontal_wrapped(|ui| {
                let old = state.choice.clone();
                egui::ComboBox::from_id_salt("qc_live_file")
                    .width(ui.available_width().min(400.))
                    .selected_text(state.choice.as_ref().map(|p| p.display().to_string())
                        .unwrap_or_else(|| "Latest output".into()))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut state.choice, None, "Latest output");
                        for file in &state.snapshot.files {
                            ui.selectable_value(&mut state.choice,
                                Some(file.relative.clone()), file.relative.display().to_string());
                        }
                    });
                if old != state.choice {
                    state.task = None;
                    state.last_read = None;
                    state.follow = true;
                }
                super::calculation_export::button(ui, state.root.as_deref(), !running);
            });
            ui.horizontal_wrapped(|ui| {
                ui.label("Find");
                ui.add(egui::TextEdit::singleline(&mut state.search).desired_width(220.));
                ui.checkbox(&mut state.matching_only, "Matching lines only");
            });
            if let Some(error) = &state.error {
                ui.label(error);
                if ui.button("Copy error").clicked_once() {
                    ctx.copy_text(error.clone());
                }
            }
            if state.snapshot.truncated {
                ui.small("Showing the last 512 KiB. Save complete calculation includes the full output file.");
            }
            if state.snapshot.omitted_files > 0 {
                ui.small(format!("Showing the latest 200 output files; {} older files remain available in the complete calculation export.", state.snapshot.omitted_files));
            }
            if let Some(selected) = &state.snapshot.selected {
                ui.small(selected.display().to_string());
            }
            if state.snapshot.text.is_empty() {
                ui.label(if running { "Waiting for the program to write output…" }
                         else { "No output has been written." });
            } else {
                let query = state.search.to_lowercase();
                let lines: Vec<_> = state.snapshot.text.lines()
                    .filter(|line| !state.matching_only || query.is_empty()
                        || line.to_lowercase().contains(&query)).collect();
                if !query.is_empty() {
                    let matches = state.snapshot.text.lines()
                        .filter(|line| line.to_lowercase().contains(&query)).count();
                    ui.small(format!("{matches} matching lines in the displayed output"));
                }
                let row_height = ui.text_style_height(&egui::TextStyle::Monospace);
                egui::ScrollArea::both()
                    .id_salt(("qc_live_text", &state.snapshot.selected))
                    .auto_shrink([false, false])
                    .stick_to_bottom(state.follow)
                    .show_rows(ui, row_height, lines.len(), |ui, range| {
                        for i in range {
                            let line = lines[i];
                            let mut text = egui::RichText::new(line).monospace();
                            if !query.is_empty() && line.to_lowercase().contains(&query) {
                                text = text.background_color(ui.visuals().selection.bg_fill);
                            }
                            ui.add(egui::Label::new(text)
                                .wrap_mode(egui::TextWrapMode::Extend).selectable(true));
                        }
                    });
            }
        });
    state.open = open;
    response.map(|r| r.response.rect)
}

#[cfg(test)]
#[path = "live_output_tests.rs"]
mod tests;
