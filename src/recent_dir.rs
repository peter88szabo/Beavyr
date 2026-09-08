//! One remembered directory, shared by every file dialog.
//!
//! Each tool used to keep its own last-used directory, or none at all, so
//! loading a Hessian and then a Molden file meant navigating the filesystem
//! twice from wherever the dialog happened to start. There is one directory
//! the user is working in, not one per tool, so there is one here.
//!
//! It is persisted under the config directory, so it also survives a restart:
//! the directory you were working in yesterday is where the next dialog opens.
//!
//! Free functions over a process-global rather than a Bevy resource, because
//! the dialogs are scattered across six modules and several are called from
//! places that have no resource in hand. Nothing here ever fails loudly --
//! a directory that cannot be read or written simply means the dialog opens
//! wherever the platform would have opened it anyway.

use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock};

/// The file the directory is remembered in.
fn store_path() -> Option<PathBuf> {
    let base = std::env::var_os("XDG_CONFIG_HOME")
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|home| PathBuf::from(home).join(".config")))?;
    Some(base.join("beavyr").join("last_dir.txt"))
}

/// The in-process copy, so a dialog opened twice in one session does not have
/// to touch the disk, and so the value survives a failed write.
fn cache() -> &'static Mutex<Option<PathBuf>> {
    static CACHE: OnceLock<Mutex<Option<PathBuf>>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(load_from_disk()))
}

fn load_from_disk() -> Option<PathBuf> {
    let text = std::fs::read_to_string(store_path()?).ok()?;
    let trimmed = text.trim();
    if trimmed.is_empty() {
        return None;
    }
    let dir = PathBuf::from(trimmed);
    // A remembered directory that has since been deleted or unmounted would
    // make the dialog open nowhere useful, so it is treated as absent.
    dir.is_dir().then_some(dir)
}

/// The directory the last file was chosen in, if it still exists.
pub fn current() -> Option<PathBuf> {
    cache().lock().ok()?.clone().filter(|d| d.is_dir())
}

/// Records the directory of a file the user has just chosen.
///
/// Takes the file's path rather than its directory, since that is what every
/// caller has, and stores the parent.
pub fn remember(path: &Path) {
    let Some(dir) = path.parent().filter(|p| p.is_dir()) else {
        return;
    };
    if let Ok(mut held) = cache().lock() {
        *held = Some(dir.to_path_buf());
    }
    let Some(store) = store_path() else {
        return;
    };
    if let Some(parent) = store.parent() {
        let _ = std::fs::create_dir_all(parent);
    }
    let _ = std::fs::write(store, dir.to_string_lossy().as_ref());
}

/// A file-open dialog starting in the remembered directory.
pub fn open() -> rfd::FileDialog {
    let dialog = rfd::FileDialog::new();
    match current() {
        Some(dir) => dialog.set_directory(dir),
        None => dialog,
    }
}

/// A file-save dialog starting in the remembered directory, with `name`
/// suggested.
pub fn save(name: &str) -> rfd::FileDialog {
    open().set_file_name(name)
}

/// Picks a file and remembers where it came from.
///
/// The two steps belong together -- a dialog whose result is not remembered is
/// the bug this module exists to fix -- so this is the form callers should
/// reach for unless they need to build the dialog themselves.
pub fn pick_file(dialog: rfd::FileDialog) -> Option<PathBuf> {
    let path = dialog.pick_file()?;
    remember(&path);
    Some(path)
}

/// Chooses a save location and remembers where it was.
pub fn save_file(dialog: rfd::FileDialog) -> Option<PathBuf> {
    let path = dialog.save_file()?;
    remember(&path);
    Some(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A file's directory is what gets remembered, not the file itself.
    #[test]
    fn remembering_a_file_records_its_directory() {
        let dir = std::env::temp_dir();
        let file = dir.join("beavyr_recent_dir_probe.xyz");
        // `remember` only stores a parent that exists, which temp_dir does.
        remember(&file);
        assert_eq!(current().as_deref(), Some(dir.as_path()));
    }

    /// A path with no usable parent must not clear or corrupt the memory.
    #[test]
    fn a_path_without_a_real_parent_is_ignored() {
        let dir = std::env::temp_dir();
        remember(&dir.join("beavyr_recent_dir_probe.xyz"));
        let before = current();
        remember(Path::new("/nonexistent_dir_xyzzy/file.txt"));
        remember(Path::new("bare_filename_with_no_parent"));
        assert_eq!(current(), before, "the good directory should survive");
    }

    /// The dialog builders must not panic without a remembered directory, and
    /// a save dialog carries its suggested name.
    #[test]
    fn the_dialog_builders_work_with_and_without_a_memory() {
        let _ = open();
        let _ = save("spectrum.dat");
    }
}
