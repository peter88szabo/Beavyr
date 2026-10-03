//! Persistent history of explicit structure loads and user snapshots.
//! Each entry is an atomically written XYZ file, shared between viewer windows.
#[cfg(test)]
mod tests;

use crate::ui_style::ResponseExt;

use crate::molecule::{composition, Molecule};
use crate::trajectory::{parse_multi_xyz, TrajectoryState};
use anyhow::{bail, Context, Result};
use bevy::prelude::*;
use bevy_egui::egui;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};
const LIMIT: usize = 20;
const HEADER: &str = "Beavyr history v1";
static SERIAL: AtomicU64 = AtomicU64::new(0);

#[derive(Clone, Debug, PartialEq)]
pub struct Snapshot {
    pub atoms: Vec<String>,
    pub positions: Vec<Vec3>,
}
impl Snapshot {
    pub fn capture(mol: &Molecule) -> Option<Self> {
        if mol.atoms.is_empty()
            || mol.atoms.len() != mol.pos.len()
            || mol.pos.iter().any(|p| !p.is_finite())
            || mol
                .atoms
                .iter()
                .any(|a| a.is_empty() || a.chars().any(char::is_whitespace))
        {
            return None;
        }
        Some(Self {
            atoms: mol.atoms.clone(),
            positions: mol.pos.clone(),
        })
    }
    fn xyz(&self, name: &str, saved: u128) -> String {
        let name: String = name
            .chars()
            .map(|c| if c.is_control() { ' ' } else { c })
            .collect();
        let mut text = format!("{}\n{HEADER}\t{saved}\t{name}\n", self.atoms.len());
        for (atom, p) in self.atoms.iter().zip(&self.positions) {
            text.push_str(&format!("{atom} {} {} {}\n", p.x, p.y, p.z));
        }
        text
    }
}
#[derive(Clone)]
pub struct Entry {
    path: PathBuf,
    saved: u128,
    pub name: String,
    pub snapshot: Snapshot,
}
fn history_directory() -> Option<PathBuf> {
    #[cfg(target_os = "windows")]
    if let Some(base) = std::env::var_os("APPDATA") {
        return Some(PathBuf::from(base).join("Beavyr/structure_history"));
    }
    #[cfg(target_os = "macos")]
    if let Some(home) = std::env::var_os("HOME") {
        return Some(
            PathBuf::from(home).join("Library/Application Support/Beavyr/structure_history"),
        );
    }
    std::env::var_os("XDG_CONFIG_HOME")
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|home| PathBuf::from(home).join(".config")))
        .map(|base| base.join("beavyr/structure_history"))
}
fn read_entry(path: PathBuf) -> Result<Entry> {
    let text = fs::read_to_string(&path)?;
    let mut comment = text
        .lines()
        .nth(1)
        .context("Missing history metadata")?
        .splitn(3, '\t');
    if comment.next() != Some(HEADER) {
        bail!("Unrecognized history entry");
    }
    let saved = comment
        .next()
        .context("Missing save time")?
        .parse::<u128>()?;
    let name = comment
        .next()
        .context("Missing structure name")?
        .to_string();
    let mut frames = parse_multi_xyz(&text).map_err(anyhow::Error::msg)?;
    if frames.len() != 1 {
        bail!("A history entry must contain one structure");
    }
    let frame = frames.remove(0);
    let snapshot = Snapshot::capture(&Molecule {
        atoms: frame.atoms,
        pos: frame.pos,
        bonds: vec![],
        hydrogen_bonds: vec![],
    })
    .context("Invalid structure coordinates")?;
    Ok(Entry {
        path,
        saved,
        name,
        snapshot,
    })
}
fn read_entries(directory: &Path) -> Result<(Vec<Entry>, usize)> {
    let files = match fs::read_dir(directory) {
        Ok(files) => files,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok((Vec::new(), 0)),
        Err(error) => return Err(error.into()),
    };
    let mut entries = Vec::new();
    let mut unreadable = 0;
    for file in files {
        let path = file?.path();
        if path.extension().is_none_or(|ext| ext != "xyz")
            || !path
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .starts_with("structure-")
        {
            continue;
        }
        match read_entry(path.clone()) {
            Ok(entry) => entries.push(entry),
            Err(_) if !path.exists() => {} // concurrently pruned by another window
            Err(_) => unreadable += 1,
        }
    }
    entries.sort_by(|a, b| b.saved.cmp(&a.saved).then_with(|| b.path.cmp(&a.path)));
    Ok((entries, unreadable))
}
fn unique_entries(entries: Vec<Entry>) -> (Vec<Entry>, Vec<PathBuf>) {
    let mut keep: Vec<Entry> = Vec::new();
    let mut discard = Vec::new();
    for entry in entries {
        if keep.len() >= LIMIT || keep.iter().any(|kept| kept.snapshot == entry.snapshot) {
            discard.push(entry.path);
        } else {
            keep.push(entry);
        }
    }
    (keep, discard)
}
#[derive(Resource)]
pub struct StructureHistory {
    directory: Option<PathBuf>,
    pub entries: Vec<Entry>,
    pub viewer_name: String,
    pub save_name: String,
    pub message: Option<String>,
    pub is_error: bool,
    refreshed: Instant,
}
impl Default for StructureHistory {
    fn default() -> Self {
        Self::open(history_directory())
    }
}
impl StructureHistory {
    /// A history with nowhere to save, for tests: nothing they do may land in
    /// the user's real Recent Structures.
    #[cfg(test)]
    pub(crate) fn detached() -> Self {
        Self::open(None)
    }

    fn open(directory: Option<PathBuf>) -> Self {
        let mut history = Self {
            directory,
            entries: Vec::new(),
            viewer_name: "Drawn structure".into(),
            save_name: String::new(),
            message: None,
            is_error: false,
            refreshed: Instant::now(),
        };
        history.refresh();
        history
    }
    fn refresh(&mut self) {
        self.refreshed = Instant::now();
        let Some(directory) = &self.directory else {
            self.message = Some(
                "No configuration directory is available for persistent structure history.".into(),
            );
            self.is_error = true;
            return;
        };
        match read_entries(directory) {
            Ok((entries, unreadable)) => {
                self.entries = unique_entries(entries).0;
                if unreadable > 0 {
                    self.message = Some(format!("Skipped {unreadable} unreadable history file(s); the files were kept on disk."));
                    self.is_error = true;
                }
            }
            Err(error) => {
                self.message = Some(format!("Could not read structure history: {error}"));
                self.is_error = true;
            }
        }
    }
    fn write_snapshot(&mut self, snapshot: Snapshot, name: &str) -> Result<()> {
        let directory = self
            .directory
            .as_ref()
            .context("No persistent history directory available")?;
        fs::create_dir_all(directory)?;
        let saved = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
        let serial = SERIAL.fetch_add(1, Ordering::Relaxed);
        let stem = format!("structure-{saved}-{}-{serial}", std::process::id());
        let pending = directory.join(format!("{stem}.tmp"));
        let finished = directory.join(format!("{stem}.xyz"));
        let result = (|| -> Result<()> {
            let mut file = fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&pending)?;
            file.write_all(snapshot.xyz(name, saved).as_bytes())?;
            file.sync_all()?;
            drop(file);
            fs::rename(&pending, &finished)?;
            Ok(())
        })();
        if result.is_err() {
            let _ = fs::remove_file(&pending);
        }
        result?;
        let (entries, _) = read_entries(directory)?;
        let (keep, discard) = unique_entries(entries);
        for path in discard {
            match fs::remove_file(path) {
                Ok(()) => {}
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
                Err(error) => return Err(error.into()),
            }
        }
        self.entries = keep;
        self.refreshed = Instant::now();
        Ok(())
    }
    fn remember(&mut self, snapshot: Snapshot, name: &str) {
        match self.write_snapshot(snapshot, name) {
            Ok(()) => {
                self.message = Some(format!("Saved {name} to history."));
                self.is_error = false;
            }
            Err(error) => {
                self.message = Some(format!("Could not save structure history: {error:#}"));
                self.is_error = true;
            }
        }
    }
    /// Snapshot the outgoing structure only at an explicit replacement boundary.
    pub fn before(&self, mol: &Molecule, traj: &TrajectoryState) -> Option<(Snapshot, String)> {
        if traj.frames.len() > 1 {
            return None;
        }
        Snapshot::capture(mol).map(|snapshot| (snapshot, self.viewer_name.clone()))
    }
    pub fn replaced(
        &mut self,
        previous: Option<(Snapshot, String)>,
        mol: &Molecule,
        traj: &TrajectoryState,
        name: &str,
    ) {
        if let Some((snapshot, label)) = previous {
            self.remember(snapshot, &label);
        }
        self.viewer_name = name.into();
        if traj.frames.len() <= 1 {
            if let Some(snapshot) = Snapshot::capture(mol) {
                self.remember(snapshot, name);
            }
        }
    }
    pub fn save_current(&mut self, mol: &Molecule) {
        let Some(snapshot) = Snapshot::capture(mol) else {
            self.message = Some("There is no valid structure to save.".into());
            self.is_error = true;
            return;
        };
        let name = if self.save_name.trim().is_empty() {
            self.viewer_name.clone()
        } else {
            self.save_name.trim().to_string()
        };
        self.remember(snapshot, &name);
    }
    pub fn restore(&mut self, entry: Entry, mol: &mut Molecule, traj: &mut TrajectoryState) {
        let previous = self.before(mol, traj);
        mol.atoms = entry.snapshot.atoms;
        mol.pos = entry.snapshot.positions;
        mol.bonds.clear();
        mol.hydrogen_bonds.clear();
        traj.clear_for_structure();
        self.replaced(previous, mol, traj, &entry.name);
    }
}
fn status(ui: &mut egui::Ui, history: &StructureHistory) {
    if let Some(message) = &history.message {
        if history.is_error {
            ui.colored_label(egui::Color32::from_rgb(225, 110, 90), message);
        } else {
            ui.weak(message);
        }
    }
}
pub fn panel(
    ui: &mut egui::Ui,
    history: &mut StructureHistory,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
) -> bool {
    if history.refreshed.elapsed() >= Duration::from_secs(2) {
        history.refresh();
    }
    let mut selected = None;
    ui.scope(|ui| {
        crate::ui_style::modern(ui);
        crate::ui_style::heading(ui, "Your saved structures", "Newest first · shared across windows and sessions");
        ui.horizontal(|ui| {
            ui.small(format!("{} of 20 saved", history.entries.len()));
            ui.menu_button("How history works", |ui| {
                ui.set_max_width(320.0);
                ui.label("Loading another structure or clearing the display preserves the previous one.");
                ui.label("Use Save to history in Molecule Editor to keep an edited version. Edits and trajectory playback do not add entries automatically.");
            });
        });
        let filter_id = egui::Id::new("recent_structures_search_v1");
        let mut filter = ui.ctx().data_mut(|data| data.get_temp::<String>(filter_id).unwrap_or_default());
        ui.add(egui::TextEdit::singleline(&mut filter)
            .desired_width(f32::INFINITY).hint_text("Search by name or formula"));
        ui.ctx().data_mut(|data| data.insert_temp(filter_id, filter.clone()));
        status(ui, history);
        let query = filter.trim().to_lowercase();
        let mut matches = 0;
        for (index, entry) in history.entries.iter().enumerate() {
            let comp = composition(&entry.snapshot.atoms);
            if !entry.name.to_lowercase().contains(&query) && !comp.formula.to_lowercase().contains(&query) {
                continue;
            }
            matches += 1;
            ui.push_id(&entry.path, |ui| {
                crate::ui_style::group(ui, "", |ui| {
                    ui.horizontal(|ui| {
                        ui.label(egui::RichText::new(format!("{:02}", index + 1)).small().color(crate::ui_style::ACCENT));
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            if ui.add(crate::ui_style::primary("Load")).clicked_once() {
                                selected = Some(entry.clone());
                            }
                            ui.add(egui::Label::new(egui::RichText::new(&entry.name).strong()).truncate())
                                .on_hover_text(&entry.name);
                        });
                    });
                    ui.horizontal_wrapped(|ui| {
                        ui.label(egui::RichText::new(&comp.formula).color(crate::ui_style::ACCENT));
                        ui.small(format!("{} atoms", comp.natoms));
                        let age = SystemTime::now()
                            .duration_since(UNIX_EPOCH).unwrap_or_default().as_nanos()
                            .saturating_sub(entry.saved) / 1_000_000_000;
                        ui.weak(if age < 60 { "Just now".into() }
                            else if age < 3600 { format!("{} min ago", age / 60) }
                            else if age < 86400 { format!("{} h ago", age / 3600) }
                            else { format!("{} days ago", age / 86400) });
                    });
                });
            });
        }
        if history.entries.is_empty() {
            crate::ui_style::group(ui, "No saved structures yet", |ui| {
                ui.small("Save a version from Molecule Editor to keep it here.");
            });
        } else if matches == 0 {
            ui.weak("No structures match your search.");
        }
    });
    ui.ctx().request_repaint_after(Duration::from_secs(2));
    if let Some(entry) = selected {
        history.restore(entry, mol, traj);
        return true;
    }
    false
}
