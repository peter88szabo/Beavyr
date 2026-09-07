//! Saving and loading the user's own element colour schemes.
//!
//! One file per scheme under `~/.config/beavyr/color_schemes/`, plus a
//! one-line file naming the scheme to start with. The format is a plain text
//! table -- `SYMBOL RRGGBB`, one element per line -- following the same
//! reasoning as the xTB path file: Beavyr depends on neither `serde` nor
//! `ron`, and a colour table does not justify adding them.
//!
//! Nothing here ever fails loudly. A scheme that cannot be written is worth
//! telling the user about, so `save` reports it, but a missing or unreadable
//! scheme is simply "not saved yet".

use std::collections::HashMap;
use std::path::PathBuf;

use bevy::prelude::Color;

/// Where saved schemes live, or `None` when there is no config directory to
/// put them in (`$HOME` and `$XDG_CONFIG_HOME` both unset).
pub fn schemes_dir() -> Option<PathBuf> {
    let base = std::env::var_os("XDG_CONFIG_HOME")
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|home| PathBuf::from(home).join(".config")))?;
    Some(base.join("beavyr").join("color_schemes"))
}

/// The file naming the scheme to load at startup.
fn default_marker_path() -> Option<PathBuf> {
    Some(schemes_dir()?.join("default.txt"))
}

/// A scheme name reduced to something safe to use as a file name.
///
/// Anything that is not alphanumeric, a dash, a dot or a space becomes an
/// underscore, and leading dots go: a name is user text, and it must not be
/// able to escape the schemes directory or hide the file it writes.
pub fn sanitize_name(name: &str) -> String {
    let cleaned: String = name
        .trim()
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' || c == '.' || c == ' ' {
                c
            } else {
                '_'
            }
        })
        .collect();
    let cleaned = cleaned.trim_matches(|c: char| c == '.' || c == ' ').to_string();
    if cleaned.is_empty() {
        "unnamed".to_string()
    } else {
        cleaned
    }
}

/// The scheme as its file contents: one `SYMBOL RRGGBB` line per element,
/// sorted so that saving the same scheme twice produces the same file.
pub fn format_scheme(colors: &HashMap<String, Color>) -> String {
    let mut symbols: Vec<&String> = colors.keys().collect();
    symbols.sort();
    let mut out = String::new();
    for symbol in symbols {
        let srgba = colors[symbol].to_srgba();
        let byte = |v: f32| (v.clamp(0.0, 1.0) * 255.0).round() as u8;
        out.push_str(&format!(
            "{symbol} {:02X}{:02X}{:02X}\n",
            byte(srgba.red),
            byte(srgba.green),
            byte(srgba.blue)
        ));
    }
    out
}

/// Reads a scheme file. Unparsable lines are skipped rather than failing the
/// whole scheme: one bad line should not cost the user every other colour.
pub fn parse_scheme(text: &str) -> HashMap<String, Color> {
    let mut colors = HashMap::new();
    for line in text.lines() {
        let mut fields = line.split_whitespace();
        let (Some(symbol), Some(hex)) = (fields.next(), fields.next()) else {
            continue;
        };
        if hex.len() != 6 {
            continue;
        }
        let channel = |at: usize| u8::from_str_radix(&hex[at..at + 2], 16).ok();
        let (Some(r), Some(g), Some(b)) = (channel(0), channel(2), channel(4)) else {
            continue;
        };
        colors.insert(
            symbol.to_string(),
            Color::srgb(r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0),
        );
    }
    colors
}

/// Every saved scheme's name, alphabetically.
pub fn list_saved() -> Vec<String> {
    let Some(dir) = schemes_dir() else {
        return Vec::new();
    };
    let Ok(entries) = std::fs::read_dir(dir) else {
        return Vec::new();
    };
    let mut names: Vec<String> = entries
        .filter_map(|entry| {
            let path = entry.ok()?.path();
            if path.extension()? != "scheme" {
                return None;
            }
            Some(path.file_stem()?.to_string_lossy().to_string())
        })
        .collect();
    names.sort();
    names
}

/// Writes a scheme, returning why it could not be written.
///
/// Unlike the xTB path, this one reports failure: the user pressed a button
/// labelled "Save", so silently not saving would be a lie.
pub fn save(name: &str, colors: &HashMap<String, Color>) -> Result<String, String> {
    let safe = sanitize_name(name);
    let Some(dir) = schemes_dir() else {
        return Err("No config directory to save into ($HOME is unset).".to_string());
    };
    std::fs::create_dir_all(&dir)
        .map_err(|e| format!("Cannot create {}: {e}", dir.display()))?;
    let path = dir.join(format!("{safe}.scheme"));
    std::fs::write(&path, format_scheme(colors))
        .map_err(|e| format!("Cannot write {}: {e}", path.display()))?;
    Ok(safe)
}

/// Loads a saved scheme, or `None` if there is no such scheme.
pub fn load(name: &str) -> Option<HashMap<String, Color>> {
    let path = schemes_dir()?.join(format!("{}.scheme", sanitize_name(name)));
    let text = std::fs::read_to_string(path).ok()?;
    let colors = parse_scheme(&text);
    (!colors.is_empty()).then_some(colors)
}

/// Records which scheme to load at startup.
pub fn set_default(name: &str) -> Result<(), String> {
    let Some(path) = default_marker_path() else {
        return Err("No config directory to save into ($HOME is unset).".to_string());
    };
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|e| format!("Cannot create {}: {e}", parent.display()))?;
    }
    std::fs::write(&path, sanitize_name(name))
        .map_err(|e| format!("Cannot write {}: {e}", path.display()))
}

/// The scheme to load at startup, if one was chosen and still exists.
pub fn default_name() -> Option<String> {
    let name = std::fs::read_to_string(default_marker_path()?)
        .ok()?
        .trim()
        .to_string();
    (!name.is_empty()).then_some(name)
}

/// The default scheme's colours, if a default is set and readable. What
/// `MolSettings::default` consults before falling back to a built-in.
pub fn load_default() -> Option<(String, HashMap<String, Color>)> {
    let name = default_name()?;
    let colors = load(&name)?;
    Some((name, colors))
}

/// Forgets a saved scheme, and the default marker if it named this one.
pub fn delete(name: &str) -> Result<(), String> {
    let safe = sanitize_name(name);
    let Some(dir) = schemes_dir() else {
        return Err("No config directory.".to_string());
    };
    let path = dir.join(format!("{safe}.scheme"));
    std::fs::remove_file(&path).map_err(|e| format!("Cannot remove {}: {e}", path.display()))?;
    if default_name().as_deref() == Some(safe.as_str()) {
        if let Some(marker) = default_marker_path() {
            let _ = std::fs::remove_file(marker);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> HashMap<String, Color> {
        let mut colors = HashMap::new();
        colors.insert("H".to_string(), Color::srgb(1.0, 1.0, 1.0));
        colors.insert("C".to_string(), Color::srgb(0.0, 0.0, 0.0));
        colors.insert("O".to_string(), Color::srgb(1.0, 0.0, 0.0));
        colors
    }

    #[test]
    fn a_scheme_round_trips_through_its_file_format() {
        let original = sample();
        let recovered = parse_scheme(&format_scheme(&original));
        assert_eq!(recovered.len(), original.len());
        for (symbol, color) in &original {
            let got = recovered[symbol].to_srgba();
            let want = color.to_srgba();
            // Eight bits per channel, so a round trip is exact to 1/255.
            assert!((got.red - want.red).abs() < 0.005, "{symbol} red");
            assert!((got.green - want.green).abs() < 0.005, "{symbol} green");
            assert!((got.blue - want.blue).abs() < 0.005, "{symbol} blue");
        }
    }

    #[test]
    fn the_file_format_is_symbol_and_hex_one_element_per_line() {
        let text = format_scheme(&sample());
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 3);
        // Sorted, so the output is stable across saves.
        assert_eq!(lines[0], "C 000000");
        assert_eq!(lines[1], "H FFFFFF");
        assert_eq!(lines[2], "O FF0000");
    }

    /// One malformed line must not cost the user the rest of the scheme.
    #[test]
    fn unparsable_lines_are_skipped_not_fatal() {
        let colors = parse_scheme(
            "H FFFFFF\n\
             this is not a colour\n\
             C 00\n\
             O ZZZZZZ\n\
             N 3050F8\n\
             \n",
        );
        assert_eq!(colors.len(), 2, "H and N survive: {colors:?}");
        assert!(colors.contains_key("H"));
        assert!(colors.contains_key("N"));
    }

    #[test]
    fn an_empty_file_yields_no_colours() {
        assert!(parse_scheme("").is_empty());
    }

    /// A scheme name is user text, so it must not be able to escape the
    /// schemes directory or write a hidden file.
    #[test]
    fn names_cannot_escape_the_schemes_directory() {
        assert_eq!(sanitize_name("../../etc/passwd"), "_.._etc_passwd");
        assert_eq!(sanitize_name("/absolute"), "_absolute");
        assert!(!sanitize_name("...hidden").starts_with('.'));
        assert!(!sanitize_name("a/b").contains('/'));
        assert!(!sanitize_name("a\\b").contains('\\'));
    }

    #[test]
    fn ordinary_names_survive_sanitising() {
        assert_eq!(sanitize_name("my scheme"), "my scheme");
        assert_eq!(sanitize_name("Jmol-tweaked"), "Jmol-tweaked");
        assert_eq!(sanitize_name("  padded  "), "padded");
        assert_eq!(sanitize_name("v1.2"), "v1.2");
    }

    #[test]
    fn an_empty_name_becomes_something_usable() {
        assert_eq!(sanitize_name(""), "unnamed");
        assert_eq!(sanitize_name("   "), "unnamed");
        assert_eq!(sanitize_name("..."), "unnamed");
    }

    /// Saving, listing, loading, defaulting and deleting, against a real
    /// directory rather than a mock.
    #[test]
    fn a_saved_scheme_can_be_listed_loaded_defaulted_and_deleted() {
        let temp = std::env::temp_dir().join(format!(
            "beavyr_scheme_test_{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&temp);
        std::fs::create_dir_all(&temp).unwrap();
        // SAFETY: single-threaded test run (`--test-threads=1`), and the
        // variable is restored below.
        let previous = std::env::var_os("XDG_CONFIG_HOME");
        unsafe { std::env::set_var("XDG_CONFIG_HOME", &temp) };

        let saved = save("My Scheme", &sample()).expect("save");
        assert_eq!(saved, "My Scheme");
        assert_eq!(list_saved(), vec!["My Scheme".to_string()]);

        let loaded = load("My Scheme").expect("load");
        assert_eq!(loaded.len(), 3);

        assert!(default_name().is_none(), "no default until one is set");
        set_default("My Scheme").expect("set default");
        assert_eq!(default_name().as_deref(), Some("My Scheme"));
        let (name, colors) = load_default().expect("load default");
        assert_eq!(name, "My Scheme");
        assert_eq!(colors.len(), 3);

        delete("My Scheme").expect("delete");
        assert!(list_saved().is_empty());
        assert!(
            default_name().is_none(),
            "deleting the default scheme clears the marker"
        );

        match previous {
            Some(value) => unsafe { std::env::set_var("XDG_CONFIG_HOME", value) },
            None => unsafe { std::env::remove_var("XDG_CONFIG_HOME") },
        }
        let _ = std::fs::remove_dir_all(&temp);
    }

    #[test]
    fn loading_a_scheme_that_was_never_saved_is_none() {
        assert!(load("definitely-not-saved-9f3a").is_none());
    }
}
