//! Saved Python launch settings. Activation affects only the calculation's child process.

use std::ffi::{OsStr, OsString};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Mutex;

#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub(crate) struct Settings {
    pub interpreter: String,
    pub activation: String,
    pub shell: String,
}

impl Default for Settings {
    fn default() -> Self {
        let fallback = if cfg!(target_os = "macos") {
            "/bin/zsh"
        } else {
            "/bin/bash"
        };
        let shell = if cfg!(windows) {
            "powershell.exe".to_owned()
        } else {
            std::env::var("SHELL")
                .ok()
                .filter(|s| {
                    matches!(
                        Path::new(s).file_name().and_then(OsStr::to_str),
                        Some("bash" | "zsh")
                    )
                })
                .unwrap_or_else(|| fallback.to_owned())
        };
        Self {
            interpreter: String::new(),
            activation: String::new(),
            shell,
        }
    }
}

fn expanded_path(value: &str) -> OsString {
    let value = value.trim();
    if let Some(relative) = value.strip_prefix("~/") {
        if let Some(home) = std::env::var_os("HOME") {
            return PathBuf::from(home).join(relative).into_os_string();
        }
    }
    value.into()
}

fn ps_literal(value: &OsStr) -> String {
    format!("'{}'", value.to_string_lossy().replace('\'', "''"))
}

impl Settings {
    pub fn interpreter(&self) -> OsString {
        if self.interpreter.trim().is_empty() {
            super::PYTHON.into()
        } else {
            expanded_path(&self.interpreter)
        }
    }

    /// Arguments are separate from the activation text, so paths and scripts remain literal.
    pub fn command(&self, arguments: impl IntoIterator<Item = impl AsRef<OsStr>>) -> Command {
        let arguments: Vec<OsString> = arguments
            .into_iter()
            .map(|a| a.as_ref().to_owned())
            .collect();
        let interpreter = self.interpreter();
        if self.activation.trim().is_empty() {
            let mut command = Command::new(interpreter);
            command.args(arguments);
            return command;
        }
        let shell = expanded_path(&self.shell);
        let shell_name = self
            .shell
            .rsplit(['/', '\\'])
            .next()
            .unwrap_or("")
            .to_ascii_lowercase();
        let mut command = Command::new(shell);
        if matches!(
            shell_name.as_str(),
            "powershell" | "powershell.exe" | "pwsh" | "pwsh.exe"
        ) {
            let call = std::iter::once(interpreter.as_os_str())
                .chain(arguments.iter().map(OsString::as_os_str))
                .map(ps_literal)
                .collect::<Vec<_>>()
                .join(" ");
            let script = format!(
                "$ErrorActionPreference = 'Stop'\ntry {{\n{}\nif (-not $?) {{ throw 'Python environment activation failed.' }}\n$global:LASTEXITCODE = 0\n& {call}\nif (-not $?) {{ if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}; exit 1 }}\nexit $LASTEXITCODE\n}} catch {{ [Console]::Error.WriteLine($_); exit 1 }}",
                self.activation
            );
            command.args(["-NoLogo", "-NonInteractive", "-Command", &script]);
        } else {
            // Interactive mode loads the user's shell hooks, including conda/micromamba.
            // There is no terminal UI, and the calculation replaces this shell with exec.
            let script = format!(
                "set -e\ntrap 'echo \"Python environment activation failed.\" >&2' ERR\n{}\nexec \"$@\"",
                self.activation
            );
            command
                .args(["-i", "-c", &script, "beavyr-python"])
                .arg(interpreter)
                .args(arguments);
        }
        command
    }
}

static SETTINGS: Mutex<Vec<(&'static str, Settings)>> = Mutex::new(Vec::new());

fn config_path(module: &str) -> Option<PathBuf> {
    if !matches!(module, "pyscf" | "psi4" | "scine_sparrow") {
        return None;
    }
    super::super::config::config_path_for(super::super::program::QcProgram::Xtb)
        .map(|p| p.with_file_name(format!("python_{module}.ron")))
}

pub(crate) fn load_from(path: &Path) -> Result<Settings, String> {
    match std::fs::read_to_string(path) {
        Ok(text) => ron::from_str(&text)
            .map_err(|e| format!("Cannot read Python settings {}: {e}", path.display())),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(Settings::default()),
        Err(e) => Err(format!(
            "Cannot read Python settings {}: {e}",
            path.display()
        )),
    }
}

pub(crate) fn save_to(path: &Path, settings: &Settings) -> Result<(), String> {
    let text = ron::ser::to_string_pretty(settings, ron::ser::PrettyConfig::default())
        .map_err(|e| e.to_string())?;
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    std::fs::write(path, text)
        .map_err(|e| format!("Cannot save Python settings {}: {e}", path.display()))
}

pub fn settings(module: &'static str) -> Result<Settings, String> {
    let mut saved = SETTINGS
        .lock()
        .map_err(|_| "Python settings are unavailable".to_owned())?;
    if let Some((_, value)) = saved.iter().find(|(name, _)| *name == module) {
        return Ok(value.clone());
    }
    let value = match config_path(module) {
        Some(path) => load_from(&path)?,
        None => Settings::default(),
    };
    saved.push((module, value.clone()));
    Ok(value)
}

pub fn save(module: &'static str, value: Settings) -> Result<(), String> {
    let path = config_path(module)
        .ok_or_else(|| "Cannot locate the Python settings directory".to_owned())?;
    save_to(&path, &value)?;
    let mut saved = SETTINGS
        .lock()
        .map_err(|_| "Python settings are unavailable".to_owned())?;
    saved.retain(|(name, _)| *name != module);
    saved.push((module, value));
    super::forget();
    Ok(())
}

pub fn command(
    module: &'static str,
    arguments: impl IntoIterator<Item = impl AsRef<OsStr>>,
) -> Command {
    match settings(module) {
        Ok(value) => value.command(arguments),
        // An unreadable configuration must never silently select another environment.
        Err(_) => Command::new(""),
    }
}

#[cfg(test)]
#[path = "python_launch_tests.rs"]
mod tests;
