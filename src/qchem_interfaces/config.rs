//! Beavyr's first piece of on-disk state: just the xTB executable path.
//!
//! The design considered a RON config file (matching `Tautomer`, which
//! already depends on `serde`/`ron`), but Beavyr's own `Cargo.toml` carries
//! neither dependency, and one string does not justify adding them. A plain
//! text file -- the path, and nothing else -- is the whole format.

use std::path::PathBuf;

/// Where the xTB path is remembered, or `None` if there is no sensible
/// location (`$HOME` unset) -- callers treat that exactly like an empty file.
pub fn config_path() -> Option<PathBuf> {
    let base = std::env::var_os("XDG_CONFIG_HOME")
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|home| PathBuf::from(home).join(".config")))?;
    Some(base.join("beavyr").join("xtb_path.txt"))
}

/// The remembered xTB path, or an empty string if nothing is saved yet, the
/// file is missing, or it cannot be read. Never errors: a missing config is
/// simply "nothing configured", not a failure.
pub fn load_xtb_path() -> String {
    let Some(path) = config_path() else {
        return String::new();
    };
    std::fs::read_to_string(path)
        .unwrap_or_default()
        .trim()
        .to_string()
}

/// Best-effort save. A failure (unwritable config directory, for instance)
/// is not surfaced to the caller: forgetting the path across restarts is a
/// minor inconvenience, never a reason to interrupt the user's session.
pub fn save_xtb_path(path: &str) {
    let Some(config_file) = config_path() else {
        return;
    };
    if let Some(parent) = config_file.parent() {
        let _ = std::fs::create_dir_all(parent);
    }
    let _ = std::fs::write(config_file, path.trim());
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Exercises the real save/load round trip against a temp `$HOME`,
    /// rather than mocking the filesystem -- this is a thin enough wrapper
    /// that the real calls are the whole thing worth testing.
    #[test]
    fn save_then_load_round_trips() {
        let dir = std::env::temp_dir().join(format!(
            "beavyr_config_test_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        // SAFETY: this test does not run concurrently with anything else
        // that reads HOME/XDG_CONFIG_HOME in this process (the whole test
        // binary runs single-threaded here via --test-threads=1, and no
        // other test in this crate touches either variable).
        unsafe {
            std::env::remove_var("XDG_CONFIG_HOME");
            std::env::set_var("HOME", &dir);
        }

        assert_eq!(load_xtb_path(), "", "nothing saved yet");
        save_xtb_path("/usr/local/bin/xtb");
        assert_eq!(load_xtb_path(), "/usr/local/bin/xtb");

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_missing_file_loads_as_empty_rather_than_erroring() {
        let dir = std::env::temp_dir().join(format!(
            "beavyr_config_test_missing_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        // Deliberately do not create `dir`.
        unsafe {
            std::env::remove_var("XDG_CONFIG_HOME");
            std::env::set_var("HOME", &dir);
        }
        assert_eq!(load_xtb_path(), "");
    }
}
