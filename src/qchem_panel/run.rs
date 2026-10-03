//! Running a job from the general panel, and handing its results on.
//!
//! The process handling is the shape every other backend in Beavyr uses: a task
//! on the compute pool, a child process a cancel flag can kill, a scratch
//! directory kept on failure and discarded on success. What differs is that one
//! task here covers every program and every job, so the branching is in what
//! gets written before the run and what gets read after it.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Child, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use bevy::prelude::*;
use bevy::tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task};

use crate::qchem_interfaces::job::{ExcitedStateOptions, JobType, QcMethod};
use crate::qchem_interfaces::method::MethodConfig;
use crate::qchem_interfaces::orca_run::{self, Energetics, ScfIteration};
use crate::qchem_interfaces::program::QcProgram;
use crate::qchem_interfaces::xtb_optimize::{
    create_xtb_run_dir, discard_xtb_run_dir, xtb_scratch_dir,
};

/// How often a cancelled run is checked for having died.
const POLL_INTERVAL: Duration = Duration::from_millis(50);

/// Everything a finished run produced.
///
/// Each field is `None` where the job did not produce that thing, rather than
/// empty: "no trajectory" and "a trajectory with no frames" are different, and
/// only the second is worth complaining about.
pub struct QcRunOutput {
    pub job: JobType,
    pub program: QcProgram,
    /// The program's own output, kept because a successful run's scratch
    /// directory is discarded moments later and the summary reads from this.
    pub log: String,
    /// A multi-frame XYZ in Beavyr's own `cycle <n> E = <v> Eh` format.
    pub trajectory_text: Option<String>,
    /// Where the Hessian file was left, for the Vibrations panel to load.
    pub hessian_path: Option<PathBuf>,
    /// The Molden file, when one was asked for. Loaded into the Surface tool
    /// and offered to save.
    pub molden_path: Option<PathBuf>,
    /// Files this run produced that are worth keeping, as (what it is, where
    /// it is). The summary offers each one a Save button.
    ///
    /// The scratch directory is discarded when a run succeeds, so anything
    /// listed here keeps it alive until the panel is done with it.
    pub saveable: Vec<(&'static str, PathBuf)>,
    /// Whether the run directory survived, which it does when a Hessian still
    /// has to be read out of it.
    pub kept_dir: Option<PathBuf>,
    pub energetics: Energetics,
    pub scf: Vec<ScfIteration>,
    /// `None` when the job never optimised anything.
    pub converged: Option<bool>,
}

/// A run in flight, and what cancelling it needs.
#[derive(Resource, Default)]
pub struct QcRunTask {
    task: Option<Task<Result<QcRunOutput, String>>>,
    cancel: Option<Arc<AtomicBool>>,
    child: Option<Arc<Mutex<Option<Child>>>>,
    run_dir: Option<PathBuf>,
    started_at: Option<Instant>,
    /// The last finished run, for the summary window.
    pub last: Option<QcRunOutput>,
    pub last_message: Option<String>,
    pub last_is_error: bool,
    pub last_duration: Option<Duration>,
}

impl QcRunTask {
    pub fn is_running(&self) -> bool {
        self.task.is_some()
    }

    pub fn elapsed(&self) -> Option<Duration> {
        self.started_at.map(|t| t.elapsed())
    }

    /// Asks the run to stop and kills the child, so a long job does not carry
    /// on burning cores after the user has given up on it.
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

    /// Starts a run. Does nothing while one is in flight, so a second click
    /// cannot spawn a competing process.
    #[allow(clippy::too_many_arguments)]
    pub fn start(
        &mut self,
        program: QcProgram,
        method: Option<QcMethod>,
        job: JobType,
        binary: PathBuf,
        atoms: Vec<String>,
        positions_angstrom: Vec<f64>,
        charge: i32,
        multiplicity: i32,
        config: MethodConfig,
        excited: ExcitedStateOptions,
        wavefunction: bool,
        // The input as the user wrote it, when they chose to write it. Sent
        // exactly as it stands, in place of anything generated.
        hand_written: Option<String>,
    ) {
        if self.is_running() {
            return;
        }
        let run_dir = match create_xtb_run_dir(&xtb_scratch_dir()) {
            Ok(dir) => dir,
            Err(err) => {
                self.last_message = Some(err);
                self.last_is_error = true;
                return;
            }
        };

        let cancel = Arc::new(AtomicBool::new(false));
        let child_slot = Arc::new(Mutex::new(None));
        let task_cancel = cancel.clone();
        let task_child = child_slot.clone();
        let task_dir = run_dir.clone();

        let task = AsyncComputeTaskPool::get().spawn(async move {
            run_cancellable(
                program,
                method,
                job,
                &binary,
                &task_dir,
                &atoms,
                &positions_angstrom,
                charge,
                multiplicity,
                &config,
                excited,
                wavefunction,
                hand_written.as_deref(),
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
}

/// Writes the input, runs the program, reads the results.
#[allow(clippy::too_many_arguments)]
fn run_cancellable(
    program: QcProgram,
    method: Option<QcMethod>,
    job: JobType,
    binary: &Path,
    workdir: &Path,
    atoms: &[String],
    positions_angstrom: &[f64],
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    excited: ExcitedStateOptions,
    wavefunction: bool,
    hand_written: Option<&str>,
    cancel: &AtomicBool,
    child_slot: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    if program.needs_binary_path() && binary.as_os_str().is_empty() {
        return Err(format!("Set the {} executable path first.", program.label()));
    }
    fs::create_dir_all(workdir).map_err(|e| format!("Failed to create the run directory: {e}"))?;

    // DREIDING first, and before the method check: a force field has no level of
    // theory to choose, so `method` is rightly `None` for it and demanding one
    // would refuse the only backend that needs nothing chosen. It also runs in
    // this process, so none of the machinery below applies.
    if program == QcProgram::Dreiding {
        return run_dreiding(job, workdir, atoms, positions_angstrom);
    }
    // xTB and Behemoth take their level of theory from the settings, not from
    // a method choice, so they are dispatched before the method is demanded.
    if matches!(program, QcProgram::Xtb | QcProgram::Behemoth) {
        return run_xtb_family(
            program, job, binary, workdir, atoms, positions_angstrom, charge, multiplicity,
            config, cancel, child_slot,
        );
    }

    let method = method.ok_or_else(|| "Choose a method first.".to_string())?;

    match program {
        QcProgram::Orca => {}
        QcProgram::Psi4 | QcProgram::Psi4Py => {
            return run_psi4(
                program, method, job, binary, workdir, atoms, positions_angstrom, charge,
                multiplicity, config, excited, wavefunction, hand_written, cancel, child_slot,
            );
        }
        QcProgram::SparrowPy => {
            return run_sparrow(
                job, workdir, atoms, positions_angstrom, charge, multiplicity, config,
                wavefunction, hand_written, cancel, child_slot,
            );
        }
        QcProgram::PySCF => {
            return run_pyscf(
                job, workdir, atoms, positions_angstrom, charge, multiplicity, config, excited,
                wavefunction, hand_written, cancel, child_slot,
            );
        }
        // Every program is handled above; this is only reached by a program
        // added to the list without a runner, which should say so plainly.
        _ => {
            return Err(format!(
                "{} has no runner in this panel yet.",
                program.label()
            ))
        }
    }

    let files = orca_run::job_files(workdir);
    // A hand-written input is sent exactly as it stands. Nothing is appended
    // and nothing is checked: the point of the editor is that it reaches
    // things no set of controls does, so second-guessing it would defeat it.
    let input = match hand_written {
        Some(text) => text.to_string(),
        None => orca_run::panel_input_text(
            job,
            method,
            charge,
            multiplicity,
            config,
            excited,
            &config.orca_extra,
            &orca_run::xyz_body(atoms, positions_angstrom),
        ),
    };
    fs::write(&files.input, &input).map_err(|e| format!("Failed to write the ORCA input: {e}"))?;

    let stdout = fs::File::create(&files.output)
        .map_err(|e| format!("Failed to create the output file: {e}"))?;
    let stderr = fs::File::create(workdir.join("orca_job.err"))
        .map_err(|e| format!("Failed to create the error file: {e}"))?;

    let mut command = orca_run::run_command(binary, workdir);
    command.stdout(Stdio::from(stdout)).stderr(Stdio::from(stderr));

    let child = command
        .spawn()
        .map_err(|e| format!("Could not start ORCA: {e}"))?;
    *child_slot
        .lock()
        .map_err(|_| "Internal error: the run state was poisoned".to_string())? = Some(child);

    let status = loop {
        if cancel.load(Ordering::SeqCst) {
            if let Ok(mut guard) = child_slot.lock() {
                if let Some(mut child) = guard.take() {
                    let _ = child.kill();
                    let _ = child.wait();
                }
            }
            return Err("Cancelled.".to_string());
        }
        let finished = {
            let mut guard = child_slot
                .lock()
                .map_err(|_| "Internal error: the run state was poisoned".to_string())?;
            match guard.as_mut() {
                Some(child) => child
                    .try_wait()
                    .map_err(|e| format!("Failed while waiting for ORCA: {e}"))?,
                None => return Err("Cancelled.".to_string()),
            }
        };
        match finished {
            Some(status) => {
                if let Ok(mut guard) = child_slot.lock() {
                    *guard = None;
                }
                break status;
            }
            None => std::thread::sleep(POLL_INTERVAL),
        }
    };

    let log = fs::read_to_string(&files.output).unwrap_or_default();
    // ORCA can exit zero having failed, so its own output decides.
    orca_run::check_output(&log, status.success())?;

    // The optimisation path, translated into the one trajectory format every
    // consumer in Beavyr already reads.
    let trajectory_text = job.produces_trajectory().then(|| {
        fs::read_to_string(&files.trajectory)
            .or_else(|_| fs::read_to_string(&files.final_geometry))
            .map(|text| orca_run::trajectory_to_beavyr(&text))
            .unwrap_or_default()
    });

    let hessian_path = job
        .produces_frequencies()
        .then(|| files.hessian.clone())
        .filter(|path| path.is_file());

    // Anything worth keeping. A frequency job's matrix is offered too: it was
    // computed either way, and wanting to save it is not a different job.
    let mut saveable: Vec<(&'static str, PathBuf)> = Vec::new();
    if matches!(job, JobType::Gradient) && files.engrad.is_file() {
        saveable.push(("Gradient (.engrad)", files.engrad.clone()));
    }
    if matches!(
        job,
        JobType::Hessian | JobType::Frequencies | JobType::OptimizeThenFrequencies
    ) && files.hessian.is_file()
    {
        saveable.push(("Hessian (.hess)", files.hessian.clone()));
    }
    if job.produces_trajectory() && files.final_geometry.is_file() {
        saveable.push(("Final geometry (.xyz)", files.final_geometry.clone()));
    }
    saveable.push(("ORCA output (.out)", files.output.clone()));
    // The binary orbital file, which is what ORCA restarts from and what its
    // own tools read. Worth keeping even when no Molden was asked for.
    let gbw = workdir.join(format!("{}.gbw", crate::qchem_interfaces::orca_run::JOB_BASENAME));
    if gbw.is_file() {
        saveable.push(("ORCA orbitals (.gbw)", gbw));
    }

    // Molden, when asked for. ORCA writes none itself: its converter has to be
    // run afterwards on the binary orbital file, which is why this is a second
    // process rather than a flag on the first.
    let molden_path = if wavefunction {
        crate::qchem_interfaces::orca_run::convert_to_molden(binary, workdir)
            .inspect_err(|_| ())
            .ok()
    } else {
        None
    };
    if let Some(path) = &molden_path {
        saveable.push(("Orbitals (Molden)", path.clone()));
    }

    Ok(QcRunOutput {
        job,
        program,
        energetics: orca_run::parse_energetics(&log),
        scf: orca_run::parse_scf_iterations(&log),
        converged: job
            .produces_trajectory()
            .then(|| orca_run::converged(&log)),
        molden_path,
        trajectory_text: trajectory_text.filter(|text| !text.trim().is_empty()),
        // The directory has to survive while anything still has to be read or
        // copied out of it.
        kept_dir: (hessian_path.is_some() || !saveable.is_empty())
            .then(|| workdir.to_path_buf()),
        hessian_path,
        saveable,
        log,
    })
}

/// Runs a Sparrow job.
///
/// Sparrow computes an energy, a gradient and a Hessian, and nothing else. It
/// has no optimiser and no frequency analysis, so an optimisation here is
/// Beavyr's own optimiser calling the script once per step, and a frequency job
/// is Beavyr's own projection and eigensolver on the matrix Sparrow hands over.
#[allow(clippy::too_many_arguments)]
fn run_sparrow(
    job: JobType,
    workdir: &Path,
    atoms: &[String],
    positions_angstrom: &[f64],
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    wavefunction: bool,
    hand_written: Option<&str>,
    cancel: &AtomicBool,
    child_slot: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    use crate::qchem_interfaces::sparrow_run;

    if matches!(job, JobType::TransitionState | JobType::Irc | JobType::TdDft) {
        return Err(format!(
            "Sparrow gives an energy, a gradient and a Hessian, and nothing else. It cannot do \
             {}. Use ORCA or Psi4 for that.",
            job.label()
        ));
    }

    let files = sparrow_run::job_files(workdir);

    // An optimisation is Beavyr's own, driven one gradient at a time, because
    // that is the only route Sparrow has.
    let mut trajectory_text = None;
    let mut geometry = positions_angstrom.to_vec();
    if job.produces_trajectory() {
        let (path, relaxed) = optimise_with_beavyr(
            &files, workdir, atoms, &geometry, charge, multiplicity, config, cancel, child_slot,
        )?;
        fs::write(workdir.join(crate::qchem_interfaces::behemoth::TRAJECTORY_FILE), &path)
            .map_err(|e| format!("Failed to write the trajectory: {e}"))?;
        trajectory_text = Some(path);
        geometry = relaxed;
    }

    // Then whatever property the job asked for, at the geometry we ended on.
    let property_job = if job == JobType::Optimize {
        JobType::SinglePoint
    } else if job == JobType::OptimizeThenFrequencies {
        JobType::Hessian
    } else {
        job
    };
    let log = run_sparrow_once(
        property_job, &files, workdir, atoms, &geometry, charge, multiplicity, config,
        wavefunction, hand_written, cancel, child_slot,
    )?;

    let hessian_path = job
        .produces_frequencies()
        .then(|| files.hessian.clone())
        .filter(|path| path.is_file());

    let mut saveable: Vec<(&'static str, PathBuf)> = Vec::new();
    if files.gradient.is_file() {
        saveable.push(("Gradient (Turbomole $grad)", files.gradient.clone()));
    }
    if files.hessian.is_file() {
        saveable.push(("Hessian ($hessian)", files.hessian.clone()));
    }
    if workdir.join("sparrow.out").is_file() {
        saveable.push(("Sparrow output", workdir.join("sparrow.out")));
    }
    let molden_path = files.molden.is_file().then(|| files.molden.clone());
    if let Some(path) = &molden_path {
        saveable.push(("Orbitals (Molden)", path.clone()));
    }

    Ok(QcRunOutput {
        molden_path,
        job,
        program: QcProgram::SparrowPy,
        energetics: Energetics {
            final_energy_hartree: sparrow_run::parse_energy(&log),
            ..Default::default()
        },
        // A semi-empirical calculator here reports no iteration history.
        scf: Vec::new(),
        converged: job.produces_trajectory().then_some(true),
        trajectory_text,
        kept_dir: (hessian_path.is_some() || !saveable.is_empty())
            .then(|| workdir.to_path_buf()),
        hessian_path,
        saveable,
        log,
    })
}

/// One Sparrow calculation at one geometry. Returns what it printed.
#[allow(clippy::too_many_arguments)]
fn run_sparrow_once(
    job: JobType,
    files: &crate::qchem_interfaces::sparrow_run::JobFiles,
    workdir: &Path,
    atoms: &[String],
    positions_angstrom: &[f64],
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    wavefunction: bool,
    hand_written: Option<&str>,
    cancel: &AtomicBool,
    child_slot: &Mutex<Option<Child>>,
) -> Result<String, String> {
    use crate::qchem_interfaces::sparrow_run;

    fs::write(
        &files.geometry,
        sparrow_run::geometry_text(atoms, positions_angstrom),
    )
    .map_err(|e| format!("Failed to write the Sparrow geometry: {e}"))?;
    fs::write(
        &files.script,
        hand_written.map(str::to_string).unwrap_or_else(|| {
            sparrow_run::script_text(job, charge, multiplicity, config, wavefunction)
        }),
    )
    .map_err(|e| format!("Failed to write the Sparrow script: {e}"))?;

    let stdout = fs::File::create(workdir.join("sparrow.out"))
        .map_err(|e| format!("Failed to create the output file: {e}"))?;
    let stderr = fs::File::create(workdir.join("sparrow.err"))
        .map_err(|e| format!("Failed to create the error file: {e}"))?;
    let mut command = sparrow_run::run_command(workdir);
    command.stdout(Stdio::from(stdout)).stderr(Stdio::from(stderr));

    let status = spawn_and_wait(command, "Sparrow", cancel, child_slot)?;
    let log = fs::read_to_string(workdir.join("sparrow.out")).unwrap_or_default();
    if !status.success() {
        let stderr_text = fs::read_to_string(workdir.join("sparrow.err")).unwrap_or_default();
        return Err(format!(
            "Sparrow stopped: {}",
            sparrow_run::failure_detail(&stderr_text)
        ));
    }
    Ok(log)
}

/// Beavyr's own optimiser, driving Sparrow one gradient at a time.
///
/// Returns the trajectory in the shared format and the geometry it ended on.
#[allow(clippy::too_many_arguments)]
fn optimise_with_beavyr(
    files: &crate::qchem_interfaces::sparrow_run::JobFiles,
    workdir: &Path,
    atoms: &[String],
    start_angstrom: &[f64],
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    cancel: &AtomicBool,
    child_slot: &Mutex<Option<Child>>,
) -> Result<(String, Vec<f64>), String> {
    use crate::optimizer::minimum_cartesian::{
        minimize_cartesian, CartesianMinimumMethod, CartesianMinimumOptions,
    };
    use crate::optimizer::traits::Objective;
    use crate::qchem_interfaces::hessian_file::parse_turbomole_gradient;

    const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;

    /// One Sparrow gradient call, presented as something the optimiser can use.
    struct SparrowObjective<'a> {
        files: &'a crate::qchem_interfaces::sparrow_run::JobFiles,
        workdir: &'a Path,
        atoms: &'a [String],
        charge: i32,
        multiplicity: i32,
        config: &'a MethodConfig,
        cancel: &'a AtomicBool,
        child_slot: &'a Mutex<Option<Child>>,
        /// Every geometry visited, with its energy, for the trajectory.
        frames: Vec<(Vec<f64>, f64)>,
    }

    impl SparrowObjective<'_> {
        fn evaluate(&mut self, x_bohr: &[f64]) -> anyhow::Result<(f64, Vec<f64>)> {
            let angstrom: Vec<f64> = x_bohr.iter().map(|c| c * BOHR_TO_ANGSTROM).collect();
            let log = run_sparrow_once(
                JobType::Gradient,
                self.files,
                self.workdir,
                self.atoms,
                &angstrom,
                self.charge,
                self.multiplicity,
                self.config,
                // Orbitals from an intermediate optimiser step are not wanted;
                // only the final geometry's, written by the call after the
                // optimisation finishes.
                false,
                // Each step is generated: a hand-written script would compute
                // the same geometry every time and the optimiser would never
                // move.
                None,
                self.cancel,
                self.child_slot,
            )
            .map_err(anyhow::Error::msg)?;
            let energy = crate::qchem_interfaces::sparrow_run::parse_energy(&log)
                .ok_or_else(|| anyhow::anyhow!("Sparrow printed no energy"))?;
            let text = fs::read_to_string(&self.files.gradient)
                .map_err(|e| anyhow::anyhow!("Sparrow wrote no gradient: {e}"))?;
            let gradient = parse_turbomole_gradient(&text, self.atoms.len())
                .map_err(anyhow::Error::msg)?;
            self.frames.push((x_bohr.to_vec(), energy));
            Ok((energy, gradient))
        }
    }

    impl Objective for SparrowObjective<'_> {
        fn energy(&mut self, x: &[f64]) -> anyhow::Result<f64> {
            Ok(self.evaluate(x)?.0)
        }

        fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> anyhow::Result<()> {
            let (_, g) = self.evaluate(x)?;
            grad.copy_from_slice(&g);
            Ok(())
        }

        fn energy_gradient(&mut self, x: &[f64], grad: &mut [f64]) -> anyhow::Result<f64> {
            let (energy, g) = self.evaluate(x)?;
            grad.copy_from_slice(&g);
            Ok(energy)
        }
    }

    let start_bohr: Vec<f64> = start_angstrom
        .iter()
        .map(|c| c / BOHR_TO_ANGSTROM)
        .collect();
    let mut objective = SparrowObjective {
        files,
        workdir,
        atoms,
        charge,
        multiplicity,
        config,
        cancel,
        child_slot,
        frames: Vec::new(),
    };
    let options = CartesianMinimumOptions {
        method: CartesianMinimumMethod::Bfgs,
        max_cycles: 300,
        // Quiet: this runs behind a panel, not in a terminal.
        verbosity: 0,
        ..CartesianMinimumOptions::default()
    };
    let result = minimize_cartesian(start_bohr, &mut objective, options)
        .map_err(|e| format!("Beavyr's optimizer stopped: {e}"))?;

    // The shared trajectory format, so the player and the energy plot need no
    // Sparrow branch.
    let mut text = String::new();
    for (index, (coords, energy)) in objective.frames.iter().enumerate() {
        text.push_str(&format!("{}\ncycle {} E = {energy:.10} Eh\n", atoms.len(), index + 1));
        for (atom, xyz) in atoms.iter().zip(coords.chunks_exact(3)) {
            text.push_str(&format!(
                "{atom:<3} {:18.10} {:18.10} {:18.10}\n",
                xyz[0] * BOHR_TO_ANGSTROM,
                xyz[1] * BOHR_TO_ANGSTROM,
                xyz[2] * BOHR_TO_ANGSTROM
            ));
        }
    }
    let relaxed: Vec<f64> = result.x.iter().map(|c| c * BOHR_TO_ANGSTROM).collect();
    Ok((text, relaxed))
}

/// Runs a PySCF job, from a generated script.
///
/// Shaped like the Psi4 route: one script written, one process watched, the
/// results read from files in formats Beavyr already takes.
#[allow(clippy::too_many_arguments)]
fn run_pyscf(
    job: JobType,
    workdir: &Path,
    atoms: &[String],
    positions_angstrom: &[f64],
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    excited: ExcitedStateOptions,
    wavefunction: bool,
    hand_written: Option<&str>,
    cancel: &AtomicBool,
    child_slot: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    use crate::qchem_interfaces::pyscf_run;

    // PySCF has no saddle-point search and no reaction-path follower. The
    // panel disables both, so reaching here would be a bug; it is said plainly
    // either way rather than quietly running something else.
    if matches!(job, JobType::TransitionState | JobType::Irc) {
        return Err(format!(
            "PySCF has no {}. Use ORCA or Psi4 for that.",
            job.label()
        ));
    }

    let files = pyscf_run::job_files(workdir);
    fs::write(
        &files.script,
        hand_written.map(str::to_string).unwrap_or_else(|| {
            pyscf_run::script_text(
                job, charge, multiplicity, config, excited, wavefunction, atoms,
                positions_angstrom,
            )
        }),
    )
    .map_err(|e| format!("Failed to write the PySCF script: {e}"))?;

    let stdout = fs::File::create(workdir.join("pyscf.out"))
        .map_err(|e| format!("Failed to create the output file: {e}"))?;
    let stderr = fs::File::create(workdir.join("pyscf.err"))
        .map_err(|e| format!("Failed to create the error file: {e}"))?;
    let mut command = pyscf_run::run_command(workdir);
    command.stdout(Stdio::from(stdout)).stderr(Stdio::from(stderr));

    let status = spawn_and_wait(command, "PySCF", cancel, child_slot)?;
    let log = fs::read_to_string(workdir.join("pyscf.out")).unwrap_or_default();
    if !status.success() {
        let stderr_text = fs::read_to_string(workdir.join("pyscf.err")).unwrap_or_default();
        return Err(format!(
            "PySCF stopped: {}",
            pyscf_run::failure_detail(&stderr_text)
        ));
    }

    let trajectory_text = job.produces_trajectory().then(|| {
        fs::read_to_string(&files.trajectory)
            .or_else(|_| fs::read_to_string(&files.geometry))
            .unwrap_or_default()
    });
    let hessian_path = job
        .produces_frequencies()
        .then(|| files.hessian.clone())
        .filter(|path| path.is_file());

    let mut saveable: Vec<(&'static str, PathBuf)> = Vec::new();
    if files.gradient.is_file() {
        saveable.push(("Gradient (Turbomole $grad)", files.gradient.clone()));
    }
    if files.hessian.is_file() {
        saveable.push(("Hessian ($hessian)", files.hessian.clone()));
    }
    if files.geometry.is_file() {
        saveable.push(("Final geometry (.xyz)", files.geometry.clone()));
    }
    if files.spectrum.is_file() {
        saveable.push(("Excited states", files.spectrum.clone()));
    }
    if workdir.join("pyscf.out").is_file() {
        saveable.push(("PySCF output", workdir.join("pyscf.out")));
    }
    let molden_path = files.molden.is_file().then(|| files.molden.clone());
    if let Some(path) = &molden_path {
        saveable.push(("Orbitals (Molden)", path.clone()));
    }

    Ok(QcRunOutput {
        molden_path,
        job,
        program: QcProgram::PySCF,
        energetics: Energetics {
            final_energy_hartree: pyscf_run::parse_energy(&log),
            ..Default::default()
        },
        scf: Vec::new(),
        converged: job.produces_trajectory().then_some(true),
        trajectory_text: trajectory_text.filter(|t| !t.trim().is_empty()),
        kept_dir: (hessian_path.is_some() || !saveable.is_empty())
            .then(|| workdir.to_path_buf()),
        hessian_path,
        saveable,
        log,
    })
}

/// Spawns a prepared command and waits for it, killing it if cancelled.
///
/// One copy for every backend that shells out from this panel: the polling,
/// the cancel check and the kill are identical, and three of them would drift.
fn spawn_and_wait(
    mut command: std::process::Command,
    label: &str,
    cancel: &AtomicBool,
    child_slot: &Mutex<Option<Child>>,
) -> Result<std::process::ExitStatus, String> {
    let child = command
        .spawn()
        .map_err(|e| format!("Could not start {label}: {e}"))?;
    *child_slot
        .lock()
        .map_err(|_| "Internal error: the run state was poisoned".to_string())? = Some(child);

    loop {
        if cancel.load(Ordering::SeqCst) {
            if let Ok(mut guard) = child_slot.lock() {
                if let Some(mut child) = guard.take() {
                    let _ = child.kill();
                    let _ = child.wait();
                }
            }
            return Err("Cancelled.".to_string());
        }
        let finished = {
            let mut guard = child_slot
                .lock()
                .map_err(|_| "Internal error: the run state was poisoned".to_string())?;
            match guard.as_mut() {
                Some(child) => child
                    .try_wait()
                    .map_err(|e| format!("Failed while waiting for {label}: {e}"))?,
                None => return Err("Cancelled.".to_string()),
            }
        };
        match finished {
            Some(status) => {
                if let Ok(mut guard) = child_slot.lock() {
                    *guard = None;
                }
                return Ok(status);
            }
            None => std::thread::sleep(POLL_INTERVAL),
        }
    }
}

/// Runs a DREIDING job, here in this process.
///
/// No executable, no scratch process, nothing to cancel -- it answers in
/// milliseconds. Results are still written into the run directory in the same
/// formats every other backend uses, so the trajectory player, the Vibrations
/// panel and the Save buttons need no force-field branch.
fn run_dreiding(
    job: JobType,
    workdir: &Path,
    atoms: &[String],
    positions_angstrom: &[f64],
) -> Result<QcRunOutput, String> {
    use crate::forcefield::dreiding::hessian::hessian_and_gradient;
    use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
    use crate::forcefield::dreiding::DreidingTopology;

    // A saddle search, a reaction path and excited states are all beyond a
    // generic force field. The panel disables them, so reaching here would be
    // a bug rather than a user's mistake -- but it is said plainly either way.
    if matches!(
        job,
        JobType::TransitionState | JobType::Irc | JobType::TdDft
    ) {
        return Err(format!(
            "DREIDING is a force field and cannot do {}. Choose a quantum-chemistry program \
             for that.",
            job.label()
        ));
    }

    let input_xyz = crate::qchem_interfaces::xtb_optimize::write_xyz_string(
        atoms,
        &positions_angstrom
            .chunks_exact(3)
            .map(|p| bevy::prelude::Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32))
            .collect::<Vec<_>>(),
    );

    let mut trajectory_text = None;
    let mut final_energy = None;
    let mut saveable: Vec<(&'static str, PathBuf)> = Vec::new();

    if job.produces_trajectory() {
        let output = crate::qchem_interfaces::dreiding_run::optimize(workdir, &input_xyz)?;
        final_energy = output.final_energy_hartree;
        trajectory_text = Some(output.trajectory_text);
    }

    // The Hessian, at whichever geometry we ended up with.
    let hessian_path = if job.produces_frequencies() || job == JobType::Hessian {
        let geometry = trajectory_text
            .as_deref()
            .and_then(last_frame_xyz)
            .unwrap_or_else(|| input_xyz.clone());
        let molecule = crate::qchem_interfaces::dreiding_run::molecule_for_force_field(&geometry);
        let topology = DreidingTopology::build(&molecule).map_err(|e| e.to_string())?;
        let coords_bohr: Vec<f64> = molecule
            .pos
            .iter()
            .flat_map(|p| {
                [
                    p.x as f64 / BOHR_TO_ANGSTROM,
                    p.y as f64 / BOHR_TO_ANGSTROM,
                    p.z as f64 / BOHR_TO_ANGSTROM,
                ]
            })
            .collect();
        let (hessian, _gradient) =
            hessian_and_gradient(&topology, &coords_bohr).map_err(|e| e.to_string())?;

        // Behind a `$hessian` header, the format every other backend's matrix
        // is written in and the one Beavyr's reader takes.
        let mut text = String::from("$hessian\n");
        for row in hessian.rows() {
            for value in row {
                text.push_str(&format!(" {value:22.14e}"));
            }
            text.push('\n');
        }
        text.push_str("$end\n");
        let path = workdir.join("hessian.txt");
        fs::write(&path, text).map_err(|e| format!("Failed to write the Hessian: {e}"))?;
        saveable.push(("Hessian ($hessian)", path.clone()));
        Some(path)
    } else {
        None
    };

    Ok(QcRunOutput {
        // A force field has no orbitals.
        molden_path: None,
        job,
        program: QcProgram::Dreiding,
        log: String::new(),
        energetics: Energetics {
            final_energy_hartree: final_energy,
            ..Default::default()
        },
        // A force field has no self-consistent field to iterate.
        scf: Vec::new(),
        converged: job.produces_trajectory().then_some(true),
        trajectory_text,
        kept_dir: (hessian_path.is_some() || !saveable.is_empty())
            .then(|| workdir.to_path_buf()),
        hessian_path: job.produces_frequencies().then(|| hessian_path.clone()).flatten(),
        saveable,
    })
}

/// The last frame of a multi-frame XYZ, as its own single-frame XYZ.
/// Runs an xTB or Behemoth job.
///
/// Both have long had their own optimiser, Hessian and gradient code in
/// Beavyr, written for the panels that used to run them. This drives that same
/// code, so a job here is the same calculation it always was: the optimiser's
/// run, the Hessian runner's matrix, and the gradient module's energy.
#[allow(clippy::too_many_arguments)]
fn run_xtb_family(
    program: QcProgram,
    job: JobType,
    binary: &Path,
    workdir: &Path,
    atoms: &[String],
    positions_angstrom: &[f64],
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    cancel: &AtomicBool,
    child_slot: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
    use crate::qchem_interfaces::{gradient, xtb_freq, xtb_optimize};

    if matches!(job, JobType::TransitionState | JobType::Irc | JobType::TdDft) {
        return Err(format!("{} cannot do {} from this panel.", program.label(), job.label()));
    }

    let positions: Vec<bevy::prelude::Vec3> = positions_angstrom
        .chunks_exact(3)
        .map(|p| bevy::prelude::Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32))
        .collect();
    // Dummies are placeholders and never reach the program.
    let (kept_atoms, kept_positions) = crate::molecule::without_dummies(atoms, &positions);
    let input_xyz = xtb_optimize::write_xyz_string(&kept_atoms, &kept_positions);
    let uhf = (multiplicity - 1).max(0);

    let mut log = String::new();
    let mut final_energy = None;
    let mut trajectory_text = None;
    let mut saveable: Vec<(&'static str, PathBuf)> = Vec::new();

    // One energy, with its gradient: the gradient module's own run.
    if matches!(job, JobType::SinglePoint | JobType::Gradient) {
        let dir = workdir.join("energy");
        let coords_bohr: Vec<f64> = kept_positions
            .iter()
            .flat_map(|p| {
                [
                    p.x as f64 / BOHR_TO_ANGSTROM,
                    p.y as f64 / BOHR_TO_ANGSTROM,
                    p.z as f64 / BOHR_TO_ANGSTROM,
                ]
            })
            .collect();
        gradient::prepare(program, &dir, &kept_atoms, &coords_bohr, charge, multiplicity, config)?;
        let mut command = gradient::command(program, binary, &dir, charge, multiplicity, config)?;
        let stdout_path = dir.join("program.stdout");
        let stdout = fs::File::create(&stdout_path)
            .map_err(|e| format!("Failed to create the output file: {e}"))?;
        let stderr = fs::File::create(dir.join("program.stderr"))
            .map_err(|e| format!("Failed to create the error file: {e}"))?;
        command.stdout(Stdio::from(stdout)).stderr(Stdio::from(stderr));
        let status = spawn_and_wait(command, program.label(), cancel, child_slot)?;
        log = fs::read_to_string(&stdout_path).unwrap_or_default();
        if !status.success() {
            return Err(format!(
                "{} stopped with an error ({status}). The run directory {} has its output.",
                program.label(),
                dir.display()
            ));
        }
        let (energy, _gradient) = gradient::read(program, &dir, kept_atoms.len(), &log)?;
        final_energy = Some(energy);
        if job == JobType::Gradient {
            let path = gradient::gradient_path(program, &dir);
            if path.is_file() {
                saveable.push(("Gradient ($grad)", path));
            }
        }
    }

    // An optimisation: the optimiser's own run, its path translated into
    // Beavyr's trajectory format so the player and the energy plot read it.
    if job.produces_trajectory() {
        let output = xtb_optimize::run_optimize_cancellable(
            program,
            binary,
            &workdir.join("opt"),
            &input_xyz,
            charge,
            uhf,
            multiplicity,
            config,
            cancel,
            child_slot,
        )?;
        final_energy = output.final_energy_hartree;
        log = output.log;
        trajectory_text = Some(beavyr_trajectory(
            &output.trajectory_text,
            &output
                .energy_history
                .iter()
                .map(|point| point.energy_hartree)
                .collect::<Vec<_>>(),
        ));
    }

    // The Hessian, at whichever geometry we ended up with, written as a bare
    // `$hessian` matrix like every other backend's.
    let hessian_path = if job.produces_frequencies() || job == JobType::Hessian {
        let geometry = trajectory_text
            .as_deref()
            .and_then(last_frame_xyz)
            .unwrap_or_else(|| input_xyz.clone());
        let (_n, symbols, _coords) = crate::molecule::parse_xyz_angstrom(&geometry);
        let raw = xtb_freq::run_hess_cancellable(
            program,
            binary,
            &workdir.join("hess"),
            &geometry,
            &symbols,
            charge,
            uhf,
            multiplicity,
            false,
            config,
            cancel,
            child_slot,
        )?;
        if log.is_empty() {
            log = fs::read_to_string(workdir.join("hess").join("xtb.stdout")).unwrap_or_default();
        }
        let mut text = String::from("$hessian\n");
        for row in raw.hessian.rows() {
            for value in row {
                text.push_str(&format!(" {value:22.14e}"));
            }
            text.push('\n');
        }
        text.push_str("$end\n");
        let path = workdir.join("hessian.txt");
        fs::write(&path, text).map_err(|e| format!("Failed to write the Hessian: {e}"))?;
        saveable.push(("Hessian ($hessian)", path.clone()));
        Some(path)
    } else {
        None
    };

    Ok(QcRunOutput {
        molden_path: None,
        job,
        program,
        log,
        energetics: Energetics {
            final_energy_hartree: final_energy,
            ..Default::default()
        },
        scf: Vec::new(),
        converged: job.produces_trajectory().then_some(true),
        trajectory_text,
        kept_dir: (hessian_path.is_some() || !saveable.is_empty())
            .then(|| workdir.to_path_buf()),
        hessian_path: job.produces_frequencies().then(|| hessian_path.clone()).flatten(),
        saveable,
    })
}

/// A multi-frame XYZ with each frame's comment rewritten as
/// `cycle <n> E = <value> Eh`, the one format the trajectory player and the
/// energy plot both read.
///
/// xTB writes its own comments (` energy: -5.07 gnorm: ...`), which nothing in
/// the panel's reader understands, so an optimisation's energy plot came out
/// empty. The energies are taken from the optimiser's own history where it
/// has one for the frame, and from the comment's `energy:` otherwise. A frame
/// already in Beavyr's format is left alone.
fn beavyr_trajectory(text: &str, energies: &[f64]) -> String {
    let lines: Vec<&str> = text.lines().collect();
    let mut out = String::new();
    let mut index = 0;
    let mut frame = 0;
    while index < lines.len() {
        let Ok(natoms) = lines[index].trim().parse::<usize>() else {
            index += 1;
            continue;
        };
        // A frame is the count line, a comment, then one line per atom; a
        // truncated last frame is dropped rather than written short.
        if index + natoms + 2 > lines.len() {
            break;
        }
        let comment = lines.get(index + 1).copied().unwrap_or("");
        let energy = energies.get(frame).copied().or_else(|| {
            let after = comment.split("energy:").nth(1)?;
            after.split_whitespace().next()?.parse::<f64>().ok()
        });
        out.push_str(&format!("{natoms}\n"));
        match energy {
            _ if comment.trim_start().starts_with("cycle") => {
                out.push_str(comment);
            }
            Some(value) => out.push_str(&format!("cycle {} E = {value:.12} Eh", frame + 1)),
            None => out.push_str(comment),
        }
        out.push('\n');
        for line in lines.iter().skip(index + 2).take(natoms) {
            out.push_str(line);
            out.push('\n');
        }
        index += natoms + 2;
        frame += 1;
    }
    out
}

fn last_frame_xyz(trajectory: &str) -> Option<String> {
    let lines: Vec<&str> = trajectory.lines().collect();
    let natoms: usize = lines.first()?.trim().parse().ok()?;
    // Frames are the count line, a comment, then one line per atom.
    let per_frame = natoms + 2;
    let frames = lines.len() / per_frame;
    if frames == 0 {
        return None;
    }
    let start = (frames - 1) * per_frame;
    Some(format!("{}\n", lines[start..start + per_frame].join("\n")))
}

/// Runs a Psi4 job, by whichever of its two routes was chosen.
///
/// The two differ only in what file is written and what is launched; everything
/// read back afterwards is the same, because the generated body is the same.
#[allow(clippy::too_many_arguments)]
fn run_psi4(
    program: QcProgram,
    method: QcMethod,
    job: JobType,
    binary: &Path,
    workdir: &Path,
    atoms: &[String],
    positions_angstrom: &[f64],
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    excited: ExcitedStateOptions,
    wavefunction: bool,
    hand_written: Option<&str>,
    cancel: &AtomicBool,
    child_slot: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    use crate::qchem_interfaces::psi4_run::{self, Route};

    let files = psi4_run::job_files(workdir);
    let route = crate::qchem_interfaces::xtb_optimize::psi4_route(program);
    // Sent exactly as written when the user wrote it.
    let (path, text) = match route {
        Route::Executable => (
            files.input.clone(),
            hand_written.map(str::to_string).unwrap_or_else(|| {
                psi4_run::input_text(
                    job, method, charge, multiplicity, config, excited, wavefunction, atoms,
                    positions_angstrom,
                )
            }),
        ),
        Route::Python => (
            files.script.clone(),
            hand_written.map(str::to_string).unwrap_or_else(|| {
                psi4_run::script_text(
                    job, method, charge, multiplicity, config, excited, wavefunction, atoms,
                    positions_angstrom,
                )
            }),
        ),
    };
    fs::write(&path, text).map_err(|e| format!("Failed to write the Psi4 job: {e}"))?;

    let stderr = fs::File::create(workdir.join("psi4.err"))
        .map_err(|e| format!("Failed to create the error file: {e}"))?;
    let mut command = psi4_run::run_command(route, binary, workdir, config.nproc);
    command.stderr(Stdio::from(stderr));
    // The executable writes its own output file; the Python route prints, so
    // that is captured into the same name.
    if route == Route::Python {
        let stdout = fs::File::create(workdir.join("psi4.stdout"))
            .map_err(|e| format!("Failed to create the output file: {e}"))?;
        command.stdout(Stdio::from(stdout));
    }

    let child = command
        .spawn()
        .map_err(|e| format!("Could not start Psi4: {e}"))?;
    *child_slot
        .lock()
        .map_err(|_| "Internal error: the run state was poisoned".to_string())? = Some(child);

    let status = loop {
        if cancel.load(Ordering::SeqCst) {
            if let Ok(mut guard) = child_slot.lock() {
                if let Some(mut child) = guard.take() {
                    let _ = child.kill();
                    let _ = child.wait();
                }
            }
            return Err("Cancelled.".to_string());
        }
        let finished = {
            let mut guard = child_slot
                .lock()
                .map_err(|_| "Internal error: the run state was poisoned".to_string())?;
            match guard.as_mut() {
                Some(child) => child
                    .try_wait()
                    .map_err(|e| format!("Failed while waiting for Psi4: {e}"))?,
                None => return Err("Cancelled.".to_string()),
            }
        };
        match finished {
            Some(status) => {
                if let Ok(mut guard) = child_slot.lock() {
                    *guard = None;
                }
                break status;
            }
            None => std::thread::sleep(POLL_INTERVAL),
        }
    };

    let log = fs::read_to_string(&files.output)
        .or_else(|_| fs::read_to_string(workdir.join("psi4.stdout")))
        .unwrap_or_default();
    if !status.success() {
        let stderr_text = fs::read_to_string(workdir.join("psi4.err")).unwrap_or_default();
        return Err(format!(
            "Psi4 stopped: {}",
            psi4_run::failure_detail(&stderr_text, &log)
        ));
    }

    let trajectory_text = job.produces_trajectory().then(|| {
        fs::read_to_string(&files.trajectory)
            .or_else(|_| fs::read_to_string(&files.geometry))
            .unwrap_or_default()
    });
    let hessian_path = job
        .produces_frequencies()
        .then(|| files.hessian.clone())
        .filter(|path| path.is_file());

    let mut saveable: Vec<(&'static str, PathBuf)> = Vec::new();
    if files.gradient.is_file() {
        saveable.push(("Gradient (Turbomole $grad)", files.gradient.clone()));
    }
    if files.hessian.is_file() {
        saveable.push(("Hessian ($hessian)", files.hessian.clone()));
    }
    if files.geometry.is_file() {
        saveable.push(("Final geometry (.xyz)", files.geometry.clone()));
    }
    if files.spectrum.is_file() {
        saveable.push(("Excited states", files.spectrum.clone()));
    }
    if files.output.is_file() {
        saveable.push(("Psi4 output", files.output.clone()));
    }
    let molden_path = files.molden.is_file().then(|| files.molden.clone());
    if let Some(path) = &molden_path {
        saveable.push(("Orbitals (Molden)", path.clone()));
    }

    Ok(QcRunOutput {
        molden_path,
        job,
        program,
        energetics: Energetics {
            final_energy_hartree: psi4_run::parse_energy(&log),
            ..Default::default()
        },
        // Psi4's iteration table has a layout of its own; the energy and the
        // files are what this panel shows for it rather than a half-read table.
        scf: Vec::new(),
        converged: job.produces_trajectory().then_some(status.success()),
        trajectory_text: trajectory_text.filter(|t| !t.trim().is_empty()),
        kept_dir: (hessian_path.is_some() || !saveable.is_empty())
            .then(|| workdir.to_path_buf()),
        hessian_path,
        saveable,
        log,
    })
}

#[cfg(test)]
mod tests {

    /// Reported: xTB could not be run from this panel at all. It was sent to
    /// "not wired into this panel yet -- use the Geometry Optimization panel",
    /// and that panel no longer exists. A run now reaches xTB itself: with an
    /// executable that is not there, the failure is starting xTB, not being
    /// turned away.
    #[test]
    fn xtb_is_dispatched_rather_than_turned_away() {
        let dir = std::env::temp_dir().join("beavyr_xtb_dispatch_test");
        let _ = fs::remove_dir_all(&dir);
        let cancel = AtomicBool::new(false);
        let child = Mutex::new(None);
        for job in [JobType::SinglePoint, JobType::Optimize, JobType::Frequencies] {
            let result = run_cancellable(
                QcProgram::Xtb,
                None,
                job,
                Path::new("/nonexistent/beavyr-test/xtb"),
                &dir,
                &["O".to_string(), "H".to_string(), "H".to_string()],
                &[0.0, 0.0, 0.0, 0.96, 0.0, 0.0, -0.24, 0.93, 0.0],
                0,
                1,
                &MethodConfig::default(),
                ExcitedStateOptions::default(),
                false,
                None,
                &cancel,
                &child,
            );
            let error = result.err().unwrap_or_default();
            assert!(
                !error.contains("not wired") && !error.contains("no runner"),
                "{job:?}: {error}"
            );
            assert!(!error.is_empty(), "{job:?}: a missing executable must be an error");
        }
        let _ = fs::remove_dir_all(&dir);
    }

    /// xTB's own trajectory comments (` energy: ... gnorm: ...`) mean nothing
    /// to the energy plot, which reads `cycle <n> E = <v> Eh`. Rewritten, the
    /// energies come through -- from the optimiser's history where it has
    /// one, from the comment where it does not.
    #[test]
    fn xtb_trajectory_comments_are_rewritten_for_the_energy_plot() {
        let xtbopt = "\
3
 energy: -5.070000000000 gnorm: 0.020000000000 xtb: 6.5.1
O 0.0 0.0 0.0
H 0.96 0.0 0.0
H -0.24 0.93 0.0
3
 energy: -5.080000000000 gnorm: 0.001000000000 xtb: 6.5.1
O 0.0 0.0 0.0
H 0.95 0.0 0.0
H -0.24 0.92 0.0
";
        let from_history = beavyr_trajectory(xtbopt, &[-5.07, -5.08]);
        let energies = crate::qchem_interfaces::behemoth::parse_trajectory_energies(&from_history);
        assert_eq!(energies, vec![-5.07, -5.08]);

        let from_comments = beavyr_trajectory(xtbopt, &[]);
        let energies = crate::qchem_interfaces::behemoth::parse_trajectory_energies(&from_comments);
        assert_eq!(energies, vec![-5.07, -5.08]);

        // Geometry untouched, frame count unchanged.
        assert_eq!(from_history.lines().count(), xtbopt.lines().count());
        assert!(from_history.contains("H 0.95 0.0 0.0"));
    }

    use super::*;

    /// A force field has no level of theory to choose, so its method is rightly
    /// `None`. Demanding one refused the single backend that needs nothing
    /// chosen, with a message about xTB that had nothing to do with it.
    #[test]
    fn dreiding_runs_without_a_method() {
        let dir = std::env::temp_dir().join(format!("beavyr_dreiding_panel_{}", std::process::id()));
        let _ = fs::create_dir_all(&dir);

        let atoms: Vec<String> = ["O", "H", "H"].iter().map(|s| s.to_string()).collect();
        let water = [0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692];

        let result = run_cancellable(
            QcProgram::Dreiding,
            // The panel passes `None` here, because `methods_for(Dreiding)` is
            // empty. That must not be an error.
            None,
            JobType::Optimize,
            Path::new(""),
            &dir,
            &atoms,
            &water,
            0,
            1,
            &MethodConfig::default(),
            ExcitedStateOptions::default(),
            // A force field has no orbitals to write.
            false,
            // Built from the controls, not written by hand.
            None,
            &AtomicBool::new(false),
            &Mutex::new(None),
        );

        match result {
            Ok(output) => {
                assert_eq!(output.program, QcProgram::Dreiding);
                assert!(
                    output.trajectory_text.is_some(),
                    "an optimisation should produce a path"
                );
            }
            Err(err) => {
                assert!(
                    !err.contains("Choose a method"),
                    "a force field needs no method: {err}"
                );
                // Perceiving water's atom types can legitimately fail on a
                // structure the force field does not understand; what must not
                // happen is being refused for want of a method.
            }
        }

        let _ = fs::remove_dir_all(&dir);
    }

    /// The jobs a force field genuinely cannot do are refused by name, not by
    /// being quietly run as something else.
    #[test]
    fn dreiding_refuses_what_a_force_field_cannot_do() {
        let dir = std::env::temp_dir().join(format!("beavyr_dreiding_no_{}", std::process::id()));
        let atoms: Vec<String> = vec!["O".to_string()];
        for job in [JobType::TransitionState, JobType::Irc, JobType::TdDft] {
            let err = match run_dreiding(job, &dir, &atoms, &[0.0, 0.0, 0.0]) {
                Err(err) => err,
                Ok(_) => panic!("a force field cannot do {}", job.label()),
            };
            assert!(err.contains("force field"), "{err}");
            assert!(err.contains(job.label()), "{err}");
        }
    }

    /// The last frame of a trajectory is what a follow-on Hessian is taken at.
    #[test]
    fn the_last_frame_is_taken_from_a_trajectory() {
        let trajectory = concat!(
            "3\ncycle 1 E = -1.0 Eh\nO 0 0 0\nH 0 0 1\nH 0 1 0\n",
            "3\ncycle 2 E = -2.0 Eh\nO 0 0 9\nH 0 0 8\nH 0 8 0\n"
        );
        let last = last_frame_xyz(trajectory).expect("two frames are here");
        assert!(last.contains("cycle 2"), "{last}");
        assert!(last.contains("O 0 0 9"), "{last}");
        assert!(!last.contains("cycle 1"), "{last}");
        assert_eq!(last.lines().count(), 5, "count, comment, three atoms");
        assert!(last_frame_xyz("").is_none());
    }
}

/// Moves a finished run's results into the panels that show them.
///
/// The routing is the job's own, from `JobType::result_targets`, so the panel
/// and the job description cannot disagree about where a result belongs.
#[allow(clippy::too_many_arguments)]
pub fn poll_qc_run(
    mut task: ResMut<QcRunTask>,
    mut panel: ResMut<super::QcPanelState>,
    mut layout: ResMut<crate::ui_layout::UiLayout>,
    mut traj: ResMut<crate::trajectory::TrajectoryState>,
    mut mol: ResMut<crate::molecule::Molecule>,
    mut settings: ResMut<crate::settings::MolSettings>,
    mut cam: Query<&mut crate::camera::OrbitCamera>,
    mut ev_changed: MessageWriter<crate::events::MoleculeChanged>,
    mut opt_task: ResMut<crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask>,
    mut opt_panel: ResMut<crate::qchem_interfaces::xtb_optimize::XtbPanelState>,
    mut freq_panel: ResMut<crate::qchem_interfaces::xtb_freq::XtbFreqPanelState>,
    mut freq_task: ResMut<crate::qchem_interfaces::xtb_freq::XtbFrequencyTask>,
    mut uvvis: ResMut<crate::uvvis::UvVisState>,
    mut orbitals: ResMut<crate::orbitals::OrbitalState>,
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
    task.last_duration = task.started_at.map(|t| t.elapsed());
    task.started_at = None;
    let run_dir = task.run_dir.take();

    match result {
        Ok(output) => {
            let job = output.job;
            let elapsed = task
                .last_duration
                .map(|d| {
                    format!(
                        " in {}",
                        crate::qchem_interfaces::xtb_optimize::format_elapsed(d)
                    )
                })
                .unwrap_or_default();
            task.last_message = Some(format!("{}{elapsed}.", job.label()));
            task.last_is_error = false;

            // Everything the optimizer panel's summary and energy plot read,
            // so a run started here gets the same two windows as one started
            // there. They are the same calculation; a second pair of windows
            // built to look almost the same would drift.
            opt_task.last_log = Some(output.log.clone());
            opt_task.last_duration = task.last_duration;
            opt_task.last_program = Some(output.program);
            opt_task.last_is_error = false;
            opt_task.last_energy_history = output
                .trajectory_text
                .as_deref()
                .map(|text| {
                    crate::qchem_interfaces::behemoth::parse_trajectory_energies(text)
                        .into_iter()
                        .enumerate()
                        .map(|(index, energy)| {
                            crate::qchem_interfaces::xtb_optimize::EnergyHistoryPoint {
                                iteration: index as u32 + 1,
                                energy_hartree: energy,
                            }
                        })
                        .collect()
                })
                .unwrap_or_default();
            // An optimisation that has just finished opens its energy plot:
            // the curve is the first thing worth looking at, and it should
            // not wait behind a button. Closing it is the user's call; it can
            // be reopened from the summary.
            if opt_task.last_energy_history.len() > 1 {
                opt_panel.energy_plot_open = true;
            }

            // The trajectory, into the player, which the energy plot reads too.
            // Forced to play once rather than loop: watching an optimisation
            // ping-pong back and forth makes no physical sense, which is the
            // same choice the optimizer panel already makes.
            if let Some(text) = &output.trajectory_text {
                if let Ok(frames) = crate::trajectory::parse_multi_xyz(text) {
                    if !frames.is_empty() {
                        traj.frames = frames;
                        traj.current_frame = 0;
                        traj.accum = 0.0;
                        traj.last_applied = None;
                        traj.overlay_dirty = true;
                        let mut camera = cam.iter_mut().next();
                        crate::trajectory::apply_current_frame(
                            &mut traj,
                            &mut mol,
                            &mut settings,
                            &mut ev_changed,
                            true,
                            camera.as_deref_mut(),
                        );
                        ev_changed.write(crate::events::MoleculeChanged::parse_xyz(false));
                        traj.mode = crate::trajectory::PlaybackMode::Once;
                        traj.direction = 1;
                        traj.playing = true;
                        layout.open[crate::ui_layout::Tab::Trajectory.index()] = true;
                    }
                }
            }

            // The Hessian, into the Vibrations panel, which runs it through
            // Beavyr's own projection and eigensolver -- so the modes animate
            // and the thermochemistry follows that panel's own controls.
            if let Some(path) = &output.hessian_path {
                // ORCA's `.hess` carries its own geometry. Everything else is a
                // bare `$hessian` matrix, loaded against the structure it was
                // computed at: the last frame of an optimisation, or the
                // structure that was sent.
                let bare = fs::read_to_string(path)
                    .map(|text| !text.contains("$atoms"))
                    .unwrap_or(false);
                let loaded = if bare {
                    let (atoms, coords) = match output
                        .trajectory_text
                        .as_deref()
                        .and_then(last_frame_xyz)
                    {
                        Some(frame) => {
                            let (_n, symbols, coords) =
                                crate::molecule::parse_xyz_angstrom(&frame);
                            (symbols, coords.into_iter().flatten().collect::<Vec<f64>>())
                        }
                        None => {
                            let (symbols, positions) =
                                crate::molecule::without_dummies(&mol.atoms, &mol.pos);
                            let coords = positions
                                .iter()
                                .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
                                .collect::<Vec<f64>>();
                            (symbols, coords)
                        }
                    };
                    crate::qchem_interfaces::xtb_freq::load_bare_hessian_from_path(
                        path,
                        &atoms,
                        &coords,
                        &mut freq_panel,
                        &mut freq_task,
                    )
                } else {
                    crate::qchem_interfaces::xtb_freq::load_hessian_from_path(
                        path,
                        &mut freq_panel,
                        &mut freq_task,
                    )
                };
                if loaded {
                    layout.open[crate::ui_layout::Tab::Vibrations.index()] = true;
                }
            }

            // The orbitals, into the Surface tool. Asking for them is almost
            // always so they can be looked at, so they are loaded rather than
            // merely offered to save.
            let mut loaded_matching_orbitals = false;
            if let Some(path) = &output.molden_path {
                match fs::read_to_string(path) {
                    Ok(text) => match crate::orbitals::molden::parse_molden(&text) {
                        Ok(data) => {
                            // Solid spheres hide the lobes they sit inside, so
                            // the same switch the file loader makes is made here.
                            settings.representation =
                                crate::settings::RepresentationMode::SticksRounded;
                            settings.geometry_dirty = true;
                            orbitals.adopt(
                                data,
                                None,
                                Some(format!("From this {} run", output.program.label())),
                                Some(text),
                            );
                            loaded_matching_orbitals = true;
                            layout.open[crate::ui_layout::Tab::Surface.index()] = true;
                        }
                        Err(err) => {
                            // The calculation is still perfectly good, so this
                            // is a note against it rather than a failed run.
                            orbitals.warning = Some(format!(
                                "The calculation ran, but its orbitals could not be read: {err}"
                            ));
                        }
                    },
                    Err(err) => {
                        orbitals.warning =
                            Some(format!("The orbital file could not be opened: {err}"));
                    }
                }
            }

            // The excited states, into the UV-Vis panel. ORCA's TD-DFT output
            // is exactly what that reader already takes, so nothing is parsed
            // twice.
            if job.is_excited_state() {
                // Whichever program ran: the dispatcher identifies it from the
                // output itself, so this one branch serves all of them.
                let label = std::path::PathBuf::from(format!("{} run", output.program.label()));
                match crate::uvvis::detect_and_parse(&output.log, &label) {
                    Ok(result) => {
                        uvvis.expanded_root = None;
                        uvvis.load_error = None;
                        uvvis.orbital_link = Default::default();
                        if loaded_matching_orbitals { uvvis.orbital_link.bind(&orbitals); }
                        uvvis.result = Some(result);
                        layout.open[crate::ui_layout::Tab::UvVis.index()] = true;
                    }
                    Err(err) => {
                        task.last_message =
                            Some(format!("{} ran, but its spectrum could not be read: {err}", job.label()));
                        task.last_is_error = true;
                    }
                }
            }

            // Discard the scratch directory unless something still has to be
            // read out of it.
            if output.kept_dir.is_none() {
                if let Some(dir) = &run_dir {
                    discard_xtb_run_dir(dir);
                }
            }

            panel.summary_open = true;
            task.last = Some(output);
        }
        Err(err) => {
            // Kept: its output is what explains a failure, and a cancelled job
            // may still be worth reading.
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
