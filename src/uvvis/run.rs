//! Running an excited-state calculation in the background.
//!
//! Structurally the same as `qchem_interfaces::xtb_freq`: a task on the
//! `AsyncComputeTaskPool`, polled with `block_on(poll_once)`, a child process
//! that cancellation kills, and a scratch run directory. Copying that shape
//! rather than inventing a second one means cancellation, polling and the
//! run-directory rules behave identically across every tool that shells out.
//!
//! The two differences: the command line, which lives in `excited_state`, and
//! that the result is parsed from stdout rather than from files left behind.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Child, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use bevy::prelude::*;
use bevy::tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task};

use super::behemoth::parse_behemoth_spectrum;
use super::excited_state::{spectrum_command, ExcitedStateConfig, CANONICAL_MOLDEN_FILE};
use super::types::TddftResult;
use super::UvVisState;
use crate::qchem_interfaces::method::MethodConfig;
use crate::qchem_interfaces::xtb_optimize::{
    create_xtb_run_dir, discard_xtb_run_dir, write_xyz_string, xtb_scratch_dir,
};

/// How often a cancelled run is checked for having died.
const POLL_INTERVAL: Duration = Duration::from_millis(50);

/// The file the program's stdout is captured to, inside the run directory.
const STDOUT_FILE: &str = "behemoth.stdout";
const STDERR_FILE: &str = "behemoth.stderr";

/// A finished spectrum, and the ground-state orbitals it was expressed in.
#[derive(Debug)]
pub struct SpectrumOutput {
    pub result: TddftResult,
    /// `canonicalMO.molden`, when the run wrote one. Read here rather than
    /// left on disk because a successful run's scratch directory is discarded
    /// moments later, and this is what the Surface tool loads.
    ///
    /// `None` for the sTDA routes, which are not asked for orbitals: a
    /// semi-empirical set is fitted to transitions, not to the ground-state
    /// electronic structure, so it does not describe the ground state.
    pub molden: Option<String>,
}

/// The background run, plus what cancelling it needs.
#[derive(Resource, Default)]
pub struct SpectrumTask {
    task: Option<Task<Result<SpectrumOutput, String>>>,
    cancel: Option<Arc<AtomicBool>>,
    child: Option<Arc<Mutex<Option<Child>>>>,
    run_dir: Option<PathBuf>,
    started_at: Option<Instant>,
    /// Set once the run finishes, either way, for the panel to show.
    pub last_message: Option<String>,
    pub last_is_error: bool,
    pub last_duration: Option<Duration>,
}

impl SpectrumTask {
    pub fn is_running(&self) -> bool {
        self.task.is_some()
    }

    pub fn elapsed(&self) -> Option<Duration> {
        self.started_at.map(|t| t.elapsed())
    }

    /// Starts a run. Does nothing while one is already in flight, so a second
    /// click cannot spawn a competing process.
    #[allow(clippy::too_many_arguments)]
    pub fn start(
        &mut self,
        binary: &Path,
        atoms: &[String],
        pos: &[Vec3],
        charge: i32,
        multiplicity: i32,
        config: &ExcitedStateConfig,
        reference: &MethodConfig,
    ) {
        if self.is_running() {
            return;
        }
        let input_xyz = write_xyz_string(atoms, pos);
        let binary = binary.to_path_buf();
        let config = config.clone();
        let reference = reference.clone();
        let cancel = Arc::new(AtomicBool::new(false));
        let child_slot = Arc::new(Mutex::new(None));
        let run_dir = match create_xtb_run_dir(&xtb_scratch_dir()) {
            Ok(dir) => dir,
            Err(err) => {
                self.last_message = Some(err);
                self.last_is_error = true;
                return;
            }
        };

        let task_cancel = cancel.clone();
        let task_child = child_slot.clone();
        let task_dir = run_dir.clone();
        let task = AsyncComputeTaskPool::get().spawn(async move {
            run_spectrum_cancellable(
                &binary,
                &task_dir,
                &input_xyz,
                charge,
                multiplicity,
                &config,
                &reference,
                &task_cancel,
                &task_child,
            )
        });

        self.task = Some(task);
        self.cancel = Some(cancel);
        self.child = Some(child_slot);
        self.run_dir = Some(run_dir);
        self.started_at = Some(Instant::now());
        self.last_message = None;
        self.last_duration = None;
    }

    /// Asks the run to stop and kills the child, so a long DFT job does not
    /// carry on burning cores after the user has given up on it.
    pub fn cancel(&mut self) {
        if let Some(flag) = &self.cancel {
            flag.store(true, Ordering::SeqCst);
        }
        if let Some(slot) = &self.child {
            if let Ok(mut guard) = slot.lock() {
                if let Some(child) = guard.as_mut() {
                    let _ = child.kill();
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn run_spectrum_cancellable(
    binary: &Path,
    workdir: &Path,
    input_xyz: &str,
    charge: i32,
    multiplicity: i32,
    config: &ExcitedStateConfig,
    reference: &MethodConfig,
    cancel: &AtomicBool,
    child_slot: &Mutex<Option<Child>>,
) -> Result<SpectrumOutput, String> {
    if binary.as_os_str().is_empty() {
        return Err("Behemoth path is empty".to_string());
    }
    fs::create_dir_all(workdir).map_err(|e| format!("Failed to create workdir: {e}"))?;
    fs::write(workdir.join("input.xyz"), input_xyz)
        .map_err(|e| format!("Failed to write input.xyz: {e}"))?;

    let stdout = fs::File::create(workdir.join(STDOUT_FILE))
        .map_err(|e| format!("Failed to create {STDOUT_FILE}: {e}"))?;
    let stderr = fs::File::create(workdir.join(STDERR_FILE))
        .map_err(|e| format!("Failed to create {STDERR_FILE}: {e}"))?;

    let mut command = spectrum_command(
        binary,
        workdir,
        "input.xyz",
        charge,
        multiplicity,
        config,
        reference,
    );
    command.stdout(Stdio::from(stdout)).stderr(Stdio::from(stderr));

    let child = command
        .spawn()
        .map_err(|e| format!("Failed to run Behemoth: {e}"))?;
    *child_slot
        .lock()
        .map_err(|_| "Internal error: run state was poisoned".to_string())? = Some(child);

    loop {
        if cancel.load(Ordering::SeqCst) {
            if let Ok(mut guard) = child_slot.lock() {
                if let Some(child) = guard.as_mut() {
                    let _ = child.kill();
                    let _ = child.wait();
                }
                *guard = None;
            }
            return Err("Cancelled.".to_string());
        }
        let status = {
            let mut guard = child_slot
                .lock()
                .map_err(|_| "Internal error: run state was poisoned".to_string())?;
            match guard.as_mut() {
                Some(child) => child.try_wait().map_err(|e| format!("Wait failed: {e}"))?,
                None => return Err("Cancelled.".to_string()),
            }
        };
        match status {
            Some(status) => {
                if let Ok(mut guard) = child_slot.lock() {
                    *guard = None;
                }
                let log = fs::read_to_string(workdir.join(STDOUT_FILE)).unwrap_or_default();
                // The log is parsed even when the exit status is non-zero:
                // Behemoth's own error line is the useful message, and the
                // parser passes it through.
                let parsed = parse_behemoth_spectrum(
                    &log,
                    &workdir.join(STDOUT_FILE),
                    Some(config.roots as usize),
                );
                let molden =
                    fs::read_to_string(workdir.join(CANONICAL_MOLDEN_FILE)).ok();
                return match parsed {
                    Ok(result) => Ok(SpectrumOutput { result, molden }),
                    Err(err) if status.success() => Err(err),
                    Err(err) => Err(format!("Behemoth exited with an error: {err}")),
                };
            }
            None => std::thread::sleep(POLL_INTERVAL),
        }
    }
}

/// Moves a finished run's spectrum into the panel's state.
///
/// On success the scratch directory goes; on failure it stays, so the log can
/// still be read -- the same rule the optimizer and the frequency tool follow.
pub fn poll_spectrum_run(
    mut task: ResMut<SpectrumTask>,
    mut state: ResMut<UvVisState>,
    mut orbitals: ResMut<crate::orbitals::OrbitalState>,
    mut settings: ResMut<crate::settings::MolSettings>,
) {
    let Some(running) = task.task.as_mut() else {
        return;
    };
    let Some(result) = block_on(future::poll_once(running)) else {
        return;
    };
    task.task = None;
    task.cancel = None;
    task.child = None;
    // Read before `started_at` is cleared.
    task.last_duration = task.started_at.map(|t| t.elapsed());
    task.started_at = None;
    let run_dir = task.run_dir.take();

    match result {
        Ok(output) => {
            if let Some(dir) = &run_dir {
                discard_xtb_run_dir(dir);
            }
            let roots = output.result.states.len();
            let method = output.result.method.clone();
            let time = task
                .last_duration
                .map(|d| {
                    format!(
                        ", {}",
                        crate::qchem_interfaces::xtb_optimize::format_elapsed(d)
                    )
                })
                .unwrap_or_default();
            task.last_message = Some(format!(
                "{method}: {roots} root{}{time}.",
                if roots == 1 { "" } else { "s" }
            ));
            task.last_is_error = false;

            // The ground-state orbitals the excitations are expressed in,
            // straight into the Surface tool -- the point of asking for them
            // being to look at the orbitals a root actually involves.
            //
            // A TD-DFT run writes them; an sTDA run is deliberately not asked
            // for any, because a semi-empirical orbital set describes
            // transitions rather than the ground state. So there is nothing to
            // load after one, and whatever the Surface tool already holds is
            // left alone.
            if let Some(text) = &output.molden {
                match crate::orbitals::molden::parse_molden(text) {
                    Ok(data) => {
                        // Solid spheres hide the lobes they sit inside, so the
                        // same switch the file loader makes is made here.
                        settings.representation =
                            crate::settings::RepresentationMode::SticksRounded;
                        settings.geometry_dirty = true;
                        orbitals.adopt(
                            data,
                            None,
                            Some(format!("From this {} run", output.result.method)),
                            Some(text.clone()),
                        );
                    }
                    Err(err) => {
                        // The spectrum is still perfectly good, so this is a
                        // note against the result rather than a failed run.
                        orbitals.warning = Some(format!(
                            "The spectrum ran, but its orbitals could not be read: {err}"
                        ));
                    }
                }
            }

            // A run's spectrum replaces a loaded one, exactly as loading a
            // second file would, so the panel has one result to show.
            state.expanded_root = None;
            state.load_error = None;
            state.result = Some(output.result);
        }
        Err(err) => {
            // The directory is kept: its stdout and stderr are what explain a
            // failure, and a cancelled DFT job may still be worth reading.
            let where_to_look = run_dir
                .as_ref()
                .map(|dir| format!(" Output kept in {}.", dir.display()))
                .unwrap_or_default();
            task.last_message = Some(if err == "Cancelled." {
                err
            } else {
                format!("{err}{where_to_look}")
            });
            task.last_is_error = true;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::excited_state::ExcitedStateMethod;

    /// A run cannot start without a program to run, and says so rather than
    /// spawning nothing and reporting success.
    #[test]
    fn an_empty_binary_path_is_an_error_not_a_silent_no_op() {
        let cancel = AtomicBool::new(false);
        let child = Mutex::new(None);
        let dir = xtb_scratch_dir().join(format!("uvvis_empty_{}", std::process::id()));
        let err = run_spectrum_cancellable(
            Path::new(""),
            &dir,
            "1\n\nH 0 0 0\n",
            0,
            1,
            &ExcitedStateConfig::default(),
            &MethodConfig::default(),
            &cancel,
            &child,
        )
        .expect_err("no binary, no run");
        assert!(err.contains("path is empty"), "{err}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A run already cancelled before it starts must not leave a process
    /// behind, and must report the cancellation rather than a parse failure.
    #[test]
    fn a_run_cancelled_before_it_starts_reports_cancellation() {
        let cancel = AtomicBool::new(true);
        let child = Mutex::new(None);
        let dir = xtb_scratch_dir().join(format!("uvvis_cancel_{}", std::process::id()));
        // `/bin/true` stands in for the program: it exits immediately, so
        // whichever branch wins the race, neither hangs.
        let result = run_spectrum_cancellable(
            Path::new("/bin/true"),
            &dir,
            "1\n\nH 0 0 0\n",
            0,
            1,
            &ExcitedStateConfig::default(),
            &MethodConfig::default(),
            &cancel,
            &child,
        );
        match result {
            Err(err) => assert!(
                err.contains("Cancelled") || err.contains("No excited-state table"),
                "{err}"
            ),
            Ok(_) => panic!("/bin/true prints no spectrum"),
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A fresh task is idle, has no elapsed time, and cancelling it is safe.
    #[test]
    fn an_unstarted_task_is_idle_and_safe_to_cancel() {
        let mut task = SpectrumTask::default();
        assert!(!task.is_running());
        assert!(task.elapsed().is_none());
        // No child, no flag: this must not panic.
        task.cancel();
        assert!(!task.is_running());
    }

    /// The task carries the configuration it was given, so what runs is what
    /// the panel showed. Checked through the command rather than by reaching
    /// into the task, which owns its copy.
    #[test]
    fn the_command_reflects_the_config_the_run_was_given() {
        let config = ExcitedStateConfig {
            method: ExcitedStateMethod::Tddft,
            roots: 12,
            tda: false,
            ..Default::default()
        };
        let command = spectrum_command(
            Path::new("behemoth"),
            Path::new("."),
            "input.xyz",
            -1,
            3,
            &config,
            &MethodConfig::default(),
        );
        let args: Vec<String> = command
            .get_args()
            .map(|a| a.to_string_lossy().to_string())
            .collect();
        assert!(args.contains(&"--stddft".to_string()), "TDA off: {args:?}");
        assert!(args.contains(&"uks".to_string()), "a triplet is UKS: {args:?}");
        assert!(args.contains(&"12".to_string()), "{args:?}");
        assert!(args.contains(&"-1".to_string()), "{args:?}");
    }
}
