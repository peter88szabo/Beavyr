//! Is there a Python here that can import what this backend needs?
//!
//! Asked when a panel is opened, not when a job is started. The failure this
//! exists to prevent is the quiet one: the environment is activated in a
//! terminal, Beavyr is launched from a desktop icon, and Beavyr inherits the
//! desktop session's Python instead. The job then dies minutes later with an
//! import error that reads like a broken installation rather than the wrong
//! environment.
//!
//! So every message names the interpreter it actually used. That is the whole
//! point: "PySCF is not installed" is ambiguous, "PySCF is not installed in
//! /home/peter/.venvs/science" is not.
//!
//! The check is cheap because it does not import anything. `find_spec` answers
//! "is this module here?" by looking, measured at 48 ms against this machine's
//! PySCF where a real import costs 0.7 s. Cheap enough to run on opening a
//! panel; the answer is cached until the user saves new environment settings
//! or requests another check.
//!
//! Beavyr never falls back to another backend when an environment is missing.
//! The run is blocked and the reason is named.

use std::sync::Mutex;

#[path = "python_launch.rs"]
pub(crate) mod launch;
#[path = "python_env_ui.rs"]
mod environment_ui;
pub(crate) use environment_ui::environment_row;

/// Default Python name, resolved after any configured activation command.
pub const PYTHON: &str = if cfg!(windows) { "python" } else { "python3" };

/// What was found when we went looking.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PythonEnv {
    /// The interpreter, shell, or environment activation could not be run.
    NoInterpreter {
        /// Why it could not be started, as the operating system put it.
        detail: String,
    },
    /// There is a Python, but it cannot import the module.
    ModuleMissing {
        module: &'static str,
        /// The interpreter that was searched. Naming it is what makes a
        /// launch-environment mismatch visible.
        executable: String,
        prefix: String,
    },
    /// Everything is in place.
    Ready {
        module: &'static str,
        /// Empty when the module is importable but declares no version.
        version: String,
        executable: String,
        prefix: String,
    },
}

impl PythonEnv {
    /// Whether a run may start.
    pub fn is_ready(&self) -> bool {
        matches!(self, PythonEnv::Ready { .. })
    }

    /// The short line the panel shows inline.
    ///
    /// Short on purpose. This sits in the middle of a form, between the charge
    /// and the level of theory, and a paragraph there pushes everything else
    /// off screen and reads as an error even when it is only a note. What it
    /// has to carry is the module and where it looked; why, and what to do, is
    /// in [`PythonEnv::detail`] on hover.
    pub fn short(&self) -> String {
        match self {
            PythonEnv::NoInterpreter { .. } => "Python environment unavailable".into(),
            PythonEnv::ModuleMissing { module, prefix, .. } => {
                format!("{module} not found in {}", shorten(prefix))
            }
            PythonEnv::Ready { module, version, prefix, .. } => {
                if version.is_empty() {
                    format!("{module} in {}", shorten(prefix))
                } else {
                    format!("{module} {version} in {}", shorten(prefix))
                }
            }
        }
    }

    /// The full explanation, for a tooltip.
    ///
    /// Everything the short line leaves out: which interpreter was searched,
    /// how to change the interpreter or activation command in Configure.
    pub fn detail(&self) -> String {
        match self {
            PythonEnv::NoInterpreter { detail } => format!(
                "The selected Python environment could not be started. Open Configure and set \
                 the Python interpreter or an Environment activation command for your venv, \
                 Conda or micromamba environment. Save & Check to try again; no restart is needed.\n\n{detail}"
            ),
            PythonEnv::ModuleMissing { module, executable, prefix } => format!(
                "{module} was not found in {prefix}.\n\nBeavyr searched {executable}. Open Configure \
                 and select the Python interpreter or Environment activation command for the venv, \
                 Conda or micromamba environment containing {module}. Save & Check applies the \
                 setting immediately; no restart is needed."
            ),
            PythonEnv::Ready { module, version, executable, prefix } => {
                let version = if version.is_empty() {
                    String::new()
                } else {
                    format!(" {version}")
                };
                format!("{module}{version}, found in {prefix} via {executable}.")
            }
        }
    }
}

pub(crate) fn notify_missing(ctx: &bevy_egui::egui::Context, environment: &PythonEnv) {
    if !environment.is_ready() {
        crate::runtime_warning::once(ctx, "Python environment setup required", &environment.detail());
    }
}

/// A home-relative path where that is shorter, so a long prefix does not push
/// the rest of the row off the panel.
fn shorten(path: &str) -> String {
    match std::env::var("HOME") {
        Ok(home) if !home.is_empty() && path.starts_with(&home) => {
            format!("~{}", &path[home.len()..])
        }
        _ => path.to_string(),
    }
}

/// The script that does the looking.
///
/// Plain keyed lines rather than JSON: Beavyr carries no JSON reader, the
/// producer and the consumer are both ours, and a run directory that is kept
/// on failure is one a person reads. Every value is on one line after its key.
fn probe_source(module: &str) -> String {
    format!(
        r#"import sys
found = False
version = ""
try:
    import importlib.util
    found = importlib.util.find_spec({module:?}) is not None
except Exception:
    found = False
if found:
    try:
        import importlib.metadata
        version = importlib.metadata.version({module:?})
    except Exception:
        version = ""
print("executable", sys.executable)
print("prefix", sys.prefix)
print("found", "yes" if found else "no")
print("version", version)
"#
    )
}

/// The value after `key` on its own line, or an empty string.
fn field<'a>(text: &'a str, key: &str) -> &'a str {
    text.lines()
        .find_map(|line| line.strip_prefix(key))
        .map(|rest| rest.trim())
        .unwrap_or("")
}

/// Runs a fresh probe using the backend's saved settings.
#[cfg(test)]
pub fn probe(module: &'static str) -> PythonEnv {
    match launch::settings(module) {
        Ok(settings) => probe_with_settings(module, &settings),
        Err(detail) => PythonEnv::NoInterpreter { detail },
    }
}

pub(crate) fn probe_with_settings(module: &'static str, settings: &launch::Settings) -> PythonEnv {
    let output = settings.command(["-c", &probe_source(module)]).output();

    let output = match output {
        Ok(output) => output,
        Err(err) => {
            return PythonEnv::NoInterpreter {
                detail: format!("{}: {err}", settings.interpreter().to_string_lossy()),
            }
        }
    };

    let stdout = String::from_utf8_lossy(&output.stdout);
    let executable = field(&stdout, "executable").to_string();
    // No interpreter line means the interpreter ran but produced nothing we
    // recognise, which is not a working Python for our purposes either.
    if !output.status.success() || executable.is_empty() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        let detail = stderr
            .lines()
            .rev()
            .take(8)
            .collect::<Vec<_>>().into_iter().rev().collect::<Vec<_>>().join("\n");
        let detail = if detail.trim().is_empty() { "Python or environment activation failed without diagnostic output".into() } else { detail };
        return PythonEnv::NoInterpreter { detail };
    }

    let prefix = field(&stdout, "prefix").to_string();
    if field(&stdout, "found") == "yes" {
        PythonEnv::Ready {
            module,
            version: field(&stdout, "version").to_string(),
            executable,
            prefix,
        }
    } else {
        PythonEnv::ModuleMissing {
            module,
            executable,
            prefix,
        }
    }
}

/// Keep a result associated with the exact configuration that was checked.
static CACHE: Mutex<Vec<(&'static str, launch::Settings, PythonEnv)>> = Mutex::new(Vec::new());

/// The cached answer for `module`, probing the first time it is asked for.
pub fn cached(module: &'static str) -> PythonEnv {
    let settings = match launch::settings(module) {
        Ok(settings) => settings,
        Err(detail) => return PythonEnv::NoInterpreter { detail },
    };
    if let Ok(cache) = CACHE.lock() {
        if let Some((_, _, env)) = cache.iter().find(|(name, checked, _)| *name == module && *checked == settings) {
            return env.clone();
        }
    }
    let env = probe_with_settings(module, &settings);
    if let Ok(mut cache) = CACHE.lock() {
        cache.retain(|(name, _, _)| *name != module);
        cache.push((module, settings, env.clone()));
    }
    env
}

/// Saving settings also clears probes, including when packages were installed
/// into the same environment since its previous check.
pub fn forget() {
    if let Ok(mut cache) = CACHE.lock() {
        cache.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The probe reports a module that is genuinely importable here, together
    /// with the interpreter it used. Run against the real environment rather
    /// than a mock: what is being tested is whether we can actually tell.
    #[test]
    fn a_module_that_is_present_is_found() {
        // `json` is in every Python standard library, so this passes wherever
        // there is a Python at all, and is skipped where there is not.
        match probe("json") {
            PythonEnv::Ready { module, executable, prefix, .. } => {
                assert_eq!(module, "json");
                assert!(!executable.is_empty(), "the interpreter must be named");
                assert!(!prefix.is_empty());
            }
            PythonEnv::NoInterpreter { .. } => {
                // No Python on this machine: nothing to assert.
            }
            other => panic!("json should be importable: {other:?}"),
        }
    }

    /// A module that does not exist is reported as missing, and the message
    /// names both it and the interpreter that was searched -- without which a
    /// wrong-environment mistake is indistinguishable from a broken install.
    #[test]
    fn a_missing_module_names_itself_and_the_interpreter() {
        match probe("beavyr_module_that_does_not_exist") {
            PythonEnv::ModuleMissing { module, executable, prefix } => {
                assert_eq!(module, "beavyr_module_that_does_not_exist");
                assert!(!executable.is_empty());
                let env = PythonEnv::ModuleMissing {
                    module,
                    executable: executable.clone(),
                    prefix: prefix.clone(),
                };
                // The inline line is short: the module and where we looked,
                // and nothing else. It sits in the middle of a form.
                let short = env.short();
                assert!(short.contains(module), "{short}");
                assert!(
                    short.len() < 90,
                    "the inline line must stay short, got {} chars: {short}",
                    short.len()
                );
                assert!(
                    !short.contains("restart"),
                    "advice belongs in the tooltip, not the form: {short}"
                );

                // The tooltip carries everything the short line leaves out.
                let detail = env.detail();
                assert!(detail.contains(module), "{detail}");
                assert!(detail.contains(&prefix), "{detail}");
                assert!(detail.contains(&executable), "{detail}");
                assert!(detail.contains("Save & Check"), "{detail}");
            }
            PythonEnv::NoInterpreter { .. } => {}
            other => panic!("that module cannot exist: {other:?}"),
        }
    }

    /// Only a ready environment lets a run start.
    #[test]
    fn only_a_ready_environment_is_ready() {
        assert!(!PythonEnv::NoInterpreter { detail: "x".into() }.is_ready());
        assert!(!PythonEnv::ModuleMissing {
            module: "pyscf",
            executable: "/usr/bin/python3".into(),
            prefix: "/usr".into(),
        }
        .is_ready());
        assert!(PythonEnv::Ready {
            module: "pyscf",
            version: "2.12.1".into(),
            executable: "/usr/bin/python3".into(),
            prefix: "/usr".into(),
        }
        .is_ready());
    }

    /// A ready line states the version and where it came from, so the right
    /// environment is confirmed at a glance before a long job is started.
    #[test]
    fn a_ready_message_states_the_version_and_the_environment() {
        let env = PythonEnv::Ready {
            module: "pyscf",
            version: "2.12.1".into(),
            executable: "/home/p/.venvs/science/bin/python3".into(),
            prefix: "/home/p/.venvs/science".into(),
        };
        let short = env.short();
        assert!(short.contains("pyscf 2.12.1"), "{short}");
        assert!(short.contains(".venvs/science"), "{short}");
        assert!(short.len() < 90, "{} chars: {short}", short.len());
        // The interpreter path belongs in the tooltip, not the form.
        assert!(!short.contains("bin/python3"), "{short}");
        assert!(env.detail().contains("bin/python3"), "{}", env.detail());
    }

    /// A path under the home directory shortens, because the full one is long
    /// enough to crowd the row it shares.
    #[test]
    fn a_home_path_is_shortened() {
        let Ok(home) = std::env::var("HOME") else {
            return;
        };
        if home.is_empty() {
            return;
        }
        let shortened = shorten(&format!("{home}/miniforge3"));
        assert_eq!(shortened, "~/miniforge3");
        // Anything outside home is left alone.
        assert_eq!(shorten("/usr/local"), "/usr/local");
    }

    /// Keyed lines are read by key, not by position, so an extra line printed
    /// by a future probe cannot shift the values.
    #[test]
    fn fields_are_read_by_key() {
        let text = "executable /a/b/python3\nprefix /a/b\nfound yes\nversion 1.2.3\n";
        assert_eq!(field(text, "executable"), "/a/b/python3");
        assert_eq!(field(text, "prefix"), "/a/b");
        assert_eq!(field(text, "found"), "yes");
        assert_eq!(field(text, "version"), "1.2.3");
        assert_eq!(field(text, "absent"), "");
    }

    /// The probe script names the module it was asked about, quoted, so a name
    /// with awkward characters cannot break out of the string.
    #[test]
    fn the_probe_script_quotes_the_module_name() {
        let source = probe_source("pyscf");
        assert!(source.contains("find_spec(\"pyscf\")"), "{source}");
        assert!(source.contains("version(\"pyscf\")"), "{source}");
    }
}
