//! RDA and Poor Man's NEB: endpoint capture, background generation, and results.
//!
//! The algorithms remain in optimizer::rda, imported from Behemoth. This module
//! connects them to the viewer and retains their XYZ outputs for explicit saving.

mod backend;
pub mod help;
pub mod parameters;
pub mod plot;
pub mod ui;

#[cfg(test)]
mod tests;

use std::path::Path;
use std::sync::Arc;
use std::time::{Duration, Instant};

use anyhow::{bail, Context, Result};
use bevy::prelude::*;
use bevy::tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task};

use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
use crate::molecule::Molecule;
use crate::optimizer::rda::{generate_rda_ts_guess, RdaGuessResult, RdaOptions};
use crate::optimizer::traits::Objective;
use crate::qchem_interfaces::xtb_optimize::{create_xtb_run_dir, resolve_program_executable};
use crate::qchem_interfaces::{config, method, method_ui, program::QcProgram, valence};
use crate::trajectory::{parse_multi_xyz, PlaybackMode, TrajectoryFrame, TrajectoryState};

use backend::{BackendConfig, EnergyGradient, RunControl};
use parameters::{Algorithm, Parameters};

#[derive(Clone)]
pub struct Endpoint {
    pub name: String,
    pub atoms: Vec<String>,
    pub positions: Vec<Vec3>,
}

impl Endpoint {
    pub fn capture(mol: &Molecule, name: &str) -> Result<Self> {
        let endpoint = Self {
            name: name.into(),
            atoms: mol.atoms.clone(),
            positions: mol.pos.clone(),
        };
        endpoint.validate()?;
        Ok(endpoint)
    }

    pub fn from_path(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("Could not read {}", path.display()))?;
        let mut frames = parse_multi_xyz(&text).map_err(anyhow::Error::msg)?;
        if frames.len() != 1 {
            bail!("Choose a single-structure XYZ, or display the desired trajectory frame and use 'Use current as A/B'.");
        }
        let frame = frames.remove(0);
        let endpoint = Self {
            name: path
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .into_owned(),
            atoms: frame.atoms,
            positions: frame.pos,
        };
        endpoint.validate()?;
        Ok(endpoint)
    }

    fn validate(&self) -> Result<()> {
        if self.atoms.len() < 2 || self.positions.len() != self.atoms.len() {
            bail!("An endpoint must contain at least two atoms with matching coordinates.");
        }
        if self.positions.iter().any(|p| !p.is_finite()) {
            bail!("Endpoint coordinates must be finite.");
        }
        Ok(())
    }

    pub fn coordinates_bohr(&self) -> Vec<f64> {
        self.positions
            .iter()
            .flat_map(|p| p.to_array().map(|v| v as f64 / BOHR_TO_ANGSTROM))
            .collect()
    }
}

fn validate_endpoints(reactant: &Endpoint, product: &Endpoint) -> Result<()> {
    reactant.validate()?;
    product.validate()?;
    if reactant.atoms != product.atoms {
        bail!("Reactant and product must have the same atoms in the same order. Atom numbers refer to that shared order.");
    }
    if reactant
        .positions
        .iter()
        .zip(&product.positions)
        .all(|(a, b)| a.distance(*b) < 1.0e-6)
    {
        bail!("Reactant and product are identical. Define the changed product geometry first.");
    }
    Ok(())
}

pub struct GenerationOutput {
    pub algorithm: Algorithm,
    pub method: String,
    pub guess: RdaGuessResult,
    pub energy_hartree: f64,
    pub guess_xyz: String,
    pub path_xyz: String,
    pub search_xyz: String,
    pub profile_dat: Option<String>,
    pub frames: Vec<TrajectoryFrame>,
}

#[derive(Resource)]
pub struct TsGeneration {
    pub reactant: Option<Endpoint>,
    pub product: Option<Endpoint>,
    pub parameters: Parameters,
    pub program: QcProgram,
    pub method: method::MethodConfig,
    pub charge: i32,
    pub multiplicity: i32,
    pub xtb_path: String,
    pub behemoth_path: String,
    pub output: Option<GenerationOutput>,
    pub energy_plot_open: bool,
    pub help_open: bool,
    pub message: Option<String>,
    pub is_error: bool,
    task: Option<Task<Result<GenerationOutput, String>>>,
    pub control: Arc<RunControl>,
    started: Option<Instant>,
    pub duration: Option<Duration>,
}

impl Default for TsGeneration {
    fn default() -> Self {
        let path = |program, fallback: &str| {
            let saved = config::load_path(program);
            if saved.is_empty() {
                fallback.into()
            } else {
                saved
            }
        };
        Self {
            reactant: None,
            product: None,
            parameters: Parameters::default(),
            program: QcProgram::Xtb,
            method: method::MethodConfig::default(),
            charge: 0,
            multiplicity: 1,
            xtb_path: path(QcProgram::Xtb, "xtb"),
            behemoth_path: path(QcProgram::Behemoth, "behemoth"),
            output: None,
            energy_plot_open: false,
            help_open: false,
            message: None,
            is_error: false,
            task: None,
            control: Arc::new(RunControl::default()),
            started: None,
            duration: None,
        }
    }
}

impl Drop for TsGeneration {
    fn drop(&mut self) {
        self.control.cancel();
    }
}

impl TsGeneration {
    pub fn is_running(&self) -> bool {
        self.task.is_some()
    }

    pub fn elapsed(&self) -> Option<Duration> {
        self.started.map(|t| t.elapsed()).or(self.duration)
    }

    pub fn input_options(&self) -> Result<RdaOptions> {
        let reactant = self
            .reactant
            .as_ref()
            .context("Capture or load the reactant structure.")?;
        let product = self
            .product
            .as_ref()
            .context("Capture or load the product structure.")?;
        validate_endpoints(reactant, product)?;
        if self.program == QcProgram::Dreiding {
            bail!("Choose xTB or Behemoth for TS generation.");
        }
        let mol = Molecule {
            atoms: reactant.atoms.clone(),
            pos: reactant.positions.clone(),
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        valence::validate_electronic_state(&mol, self.charge, self.multiplicity)
            .map_err(|error| anyhow::anyhow!("{error}"))?;
        let issues = method::validate(
            self.program,
            &self.method,
            self.multiplicity,
            &reactant.atoms,
        );
        if let Some(issue) = issues
            .iter()
            .find(|i| i.severity == method::Severity::Block)
        {
            bail!("{}", issue.message);
        }
        if self.method.nproc == 0 {
            bail!("Use at least one calculation thread.");
        }
        self.parameters.rda_options(reactant.atoms.len())
    }

    pub fn start(&mut self) -> Result<()> {
        if self.is_running() {
            bail!("A TS generation is already running.");
        }
        let options = self.input_options()?;
        let configured = match self.program {
            QcProgram::Xtb => &self.xtb_path,
            _ => &self.behemoth_path,
        };
        let binary =
            resolve_program_executable(self.program, configured).map_err(anyhow::Error::msg)?;
        // Subprocesses change directory; resolve relative executable paths first.
        let binary =
            std::fs::canonicalize(binary).context("Could not resolve the executable path")?;
        let backend = BackendConfig {
            program: self.program,
            binary,
            method: self.method.clone(),
            charge: self.charge,
            multiplicity: self.multiplicity,
        };
        config::save_path(self.program, configured);
        let reactant = self.reactant.clone().unwrap();
        let product = self.product.clone().unwrap();
        let algorithm = self.parameters.algorithm;
        let method = format!(
            "{} / {}",
            self.program.label(),
            method_ui::method_summary(self.program, &self.method)
        );
        let run_dir = create_xtb_run_dir(&std::env::temp_dir().join("beavyr_ts"))
            .map_err(anyhow::Error::msg)?;
        self.control = Arc::new(RunControl::default());
        let control = self.control.clone();
        self.task = Some(AsyncComputeTaskPool::get().spawn(async move {
            let mut objective = EnergyGradient::new(
                backend,
                reactant.atoms.clone(),
                run_dir.join("energy"),
                control.clone(),
            );
            let result = generate(
                &mut objective,
                &reactant,
                &product,
                options,
                algorithm,
                method,
                &run_dir,
            );
            if control.cancelled.load(std::sync::atomic::Ordering::Relaxed) {
                return Err(format!(
                    "TS generation cancelled. Scratch files: {}",
                    run_dir.display()
                ));
            }
            match result {
                Ok(output) => {
                    let _ = std::fs::remove_dir_all(&run_dir);
                    Ok(output)
                }
                Err(error) => Err(format!("{error:#}\nScratch files: {}", run_dir.display())),
            }
        }));
        self.started = Some(Instant::now());
        self.duration = None;
        self.output = None;
        self.energy_plot_open = false;
        self.message = None;
        self.is_error = false;
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn generate<O: Objective>(
    objective: &mut O,
    reactant: &Endpoint,
    product: &Endpoint,
    mut options: RdaOptions,
    algorithm: Algorithm,
    method: String,
    run_dir: &Path,
) -> Result<GenerationOutput> {
    validate_endpoints(reactant, product)?;
    std::fs::create_dir_all(run_dir)?;
    let path = |name: &str| run_dir.join(name).to_string_lossy().into_owned();
    options.rda_trajectory_file = Some(path("rda_search.xyz"));
    options.rda_chain_file = Some(path("rda_chain.xyz"));
    options.poormans_neb = algorithm == Algorithm::PoorMansNeb;
    options.poormans_neb_file = options.poormans_neb.then(|| path("poormans_neb.xyz"));
    options.poormans_neb_profile_file = options.poormans_neb.then(|| path("poormans_neb.dat"));
    let guess = generate_rda_ts_guess(
        objective,
        reactant.atoms.clone(),
        reactant.coordinates_bohr(),
        product.coordinates_bohr(),
        options,
    )?;
    let energy_hartree = objective.energy(&guess.q)?;
    if guess.q.iter().any(|v| !v.is_finite()) || !energy_hartree.is_finite() {
        bail!("TS generation returned invalid coordinates or energy.");
    }
    let guess_xyz = backend::xyz(
        &guess.atoms,
        &guess.q,
        &format!(
            "{} TS guess; {method}; E={energy_hartree:.12} Eh",
            algorithm.label()
        ),
    );
    let path_xyz = std::fs::read_to_string(run_dir.join(match algorithm {
        Algorithm::Rda => "rda_chain.xyz",
        Algorithm::PoorMansNeb => "poormans_neb.xyz",
    }))?;
    let frames = parse_multi_xyz(&path_xyz).map_err(anyhow::Error::msg)?;
    let search_xyz = std::fs::read_to_string(run_dir.join("rda_search.xyz"))?;
    let profile_dat = if algorithm == Algorithm::PoorMansNeb {
        Some(std::fs::read_to_string(run_dir.join("poormans_neb.dat"))?)
    } else {
        None
    };
    Ok(GenerationOutput {
        algorithm,
        method,
        guess,
        energy_hartree,
        guess_xyz,
        path_xyz,
        search_xyz,
        profile_dat,
        frames,
    })
}

pub fn poll_generation(mut state: ResMut<TsGeneration>) {
    let Some(task) = state.task.as_mut() else {
        return;
    };
    let Some(result) = block_on(future::poll_once(task)) else {
        return;
    };
    state.task = None;
    state.duration = state.started.take().map(|t| t.elapsed());
    match result {
        Ok(output) => {
            state.message = Some(format!(
                "{} TS guess generated. Use the buttons below to view or save it.",
                output.algorithm.label()
            ));
            state.output = Some(output);
            state.is_error = false;
        }
        Err(error) => {
            state.message = Some(error);
            state.is_error = true;
        }
    }
}

pub fn show_frames(
    frames: Vec<TrajectoryFrame>,
    mol: &mut Molecule,
    traj: &mut TrajectoryState,
) -> bool {
    let Some(first) = frames.first() else {
        return false;
    };
    mol.atoms.clone_from(&first.atoms);
    mol.pos.clone_from(&first.pos);
    traj.frames = frames;
    traj.current_frame = 0;
    traj.playing = false;
    traj.mode = PlaybackMode::Once;
    traj.direction = 1;
    traj.accum = 0.0;
    traj.last_applied = None;
    traj.overlay_enabled = false;
    traj.overlay_dirty = true;
    traj.fixed_bonds = false;
    traj.current_file = None;
    true
}
