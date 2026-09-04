//! Runs `xtb --hess` (optionally `--grad` too, for the reaction-path
//! projection) in the background, then runs Beavyr's imported
//! (from Behemoth) normal-mode analysis on the result: frequencies, mode
//! shapes, thermochemistry, and an IR spectrum built from xTB's own
//! `vibspectrum` intensities. Structurally mirrors `xtb_optimize.rs`
//! (background task via `AsyncComputeTaskPool`, same cancellation and
//! run-directory rules) rather than introducing a second pattern.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use bevy::prelude::*;
use bevy::tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task};
use bevy_egui::egui;
use ndarray::Array2;

use crate::molecule::{atomic_mass_amu, Molecule};
#[cfg(test)]
use crate::molecule::parse_xyz_angstrom;
use crate::normalmode::thermofuncs::{self, ThermoResults};
use crate::normalmode::{normal_modes_with_projection, EckartMode};
use crate::trajectory::{PlaybackMode, TrajectoryFrame, TrajectoryState};

use super::hessian_file::{parse_turbomole_gradient, parse_turbomole_hessian, parse_vibspectrum};
use super::xtb_optimize::{create_xtb_run_dir, discard_xtb_run_dir, xtb_scratch_dir};

const ANGSTROM_TO_BOHR: f64 = 1.0 / 0.529_177_210_903;
const BOHR_TO_ANGSTROM: f32 = 0.529_177_210_903;
/// Below this gradient norm (Eh/bohr), a reaction-path projection is treated
/// as meaningless (the usual case right after an optimization) and the run
/// falls back to the ordinary translation/rotation-projected Hessian
/// instead of failing outright -- Behemoth's own
/// `eckart_reactionpath_transform` would otherwise `bail!` on a zero
/// gradient, per the design spec.
const REACTION_PATH_GRADIENT_FLOOR: f64 = 1.0e-4;

/// The finished analysis: everything the panel needs to list modes, animate
/// one, and show thermochemistry/spectrum windows.
pub struct FrequencyResult {
    pub atoms: Vec<String>,
    /// Equilibrium geometry in the frame the Hessian was actually computed
    /// in (xTB's own `xtbhess.xyz`, not necessarily the on-screen frame --
    /// see the module-level note on why that distinction matters).
    pub coords_angstrom: Vec<Vec3>,
    /// Columns are Cartesian displacement eigenvectors, ordered ascending
    /// by eigenvalue -- same order as `frequencies_cm1` and `ir_intensities`.
    pub modes: Array2<f64>,
    pub frequencies_cm1: Vec<f64>,
    pub ir_intensities_km_mol: Vec<f64>,
    pub negative_indices: Vec<usize>,
    pub zero_indices: Vec<usize>,
    pub positive_indices: Vec<usize>,
    pub thermo: ThermoResults,
    pub thermo_temp_k: f64,
    pub thermo_freq_cutoff_cm1: f64,
    pub eckart_used: EckartMode,
    /// Set when a reaction-path request was silently downgraded to the
    /// ordinary translation/rotation projection
    /// because the gradient was essentially zero.
    pub reaction_path_fallback_note: Option<String>,
}

struct RawFrequencyOutput {
    atoms: Vec<String>,
    coords_bohr: Vec<f64>,
    hessian: Array2<f64>,
    vib_lines: Vec<super::hessian_file::VibSpectrumLine>,
    gradient_bohr: Option<Vec<f64>>,
}

#[derive(Resource, Default)]
pub struct XtbFrequencyTask {
    task: Option<Task<Result<RawFrequencyOutput, String>>>,
    cancel: Option<Arc<AtomicBool>>,
    child: Option<Arc<Mutex<Option<Child>>>>,
    run_dir: Option<PathBuf>,
    started_at: Option<Instant>,
    /// The multiplicity this run was started with, captured at click time so
    /// `poll_xtb_frequencies` does not need to reach back into the
    /// optimizer's panel state (which the user could have since changed).
    multiplicity_used: i32,
    pub last_message: Option<String>,
    pub last_is_error: bool,
    pub result: Option<FrequencyResult>,
}

impl XtbFrequencyTask {
    pub fn is_running(&self) -> bool {
        self.task.is_some()
    }

    pub fn elapsed(&self) -> Option<Duration> {
        self.started_at.map(|t| t.elapsed())
    }

    pub fn cancel_now(&self) {
        if let Some(cancel) = &self.cancel {
            cancel.store(true, Ordering::Relaxed);
        }
        if let Some(child) = &self.child {
            if let Ok(mut child) = child.lock() {
                if let Some(mut process) = child.take() {
                    let _ = process.kill();
                    let _ = process.wait();
                }
            }
        }
    }

    pub fn start(
        &mut self,
        xtb_path: &Path,
        atoms: &[String],
        pos: &[Vec3],
        charge: i32,
        uhf: i32,
        multiplicity: i32,
        need_gradient: bool,
    ) {
        if self.is_running() {
            return;
        }
        let input_xyz = super::xtb_optimize::write_xyz_string(atoms, pos);
        let xtb_path = xtb_path.to_path_buf();
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
            run_xtb_hess_cancellable(
                &xtb_path,
                &task_dir,
                &input_xyz,
                charge,
                uhf,
                need_gradient,
                &task_cancel,
                &task_child,
            )
        });

        self.task = Some(task);
        self.cancel = Some(cancel);
        self.child = Some(child_slot);
        self.run_dir = Some(run_dir);
        self.started_at = Some(Instant::now());
        self.multiplicity_used = multiplicity;
        self.last_message = None;
    }
}

fn run_xtb_hess_cancellable(
    xtb_path: &Path,
    workdir: &Path,
    input_xyz: &str,
    charge: i32,
    uhf: i32,
    need_gradient: bool,
    cancel: &AtomicBool,
    child_slot: &Mutex<Option<Child>>,
) -> Result<RawFrequencyOutput, String> {
    if xtb_path.as_os_str().is_empty() {
        return Err("XTB path is empty".to_string());
    }
    fs::create_dir_all(workdir).map_err(|e| format!("Failed to create workdir: {e}"))?;
    fs::write(workdir.join("input.xyz"), input_xyz)
        .map_err(|e| format!("Failed to write input.xyz: {e}"))?;

    let stdout = fs::File::create(workdir.join("xtb.stdout"))
        .map_err(|e| format!("Failed to create xtb.stdout: {e}"))?;
    let stderr = fs::File::create(workdir.join("xtb.stderr"))
        .map_err(|e| format!("Failed to create xtb.stderr: {e}"))?;
    let mut cmd = Command::new(xtb_path);
    cmd.current_dir(workdir)
        .arg("input.xyz")
        .arg("--hess")
        .arg("--chrg")
        .arg(charge.to_string())
        .stdout(Stdio::from(stdout))
        .stderr(Stdio::from(stderr));
    if uhf > 0 {
        cmd.arg("--uhf").arg(uhf.to_string());
    }
    if need_gradient {
        cmd.arg("--grad");
    }

    let child = cmd.spawn().map_err(|e| format!("Failed to run xTB: {e}"))?;
    *child_slot
        .lock()
        .map_err(|_| "xTB child-process lock was poisoned".to_string())? = Some(child);
    let status = loop {
        let mut child = child_slot
            .lock()
            .map_err(|_| "xTB child-process lock was poisoned".to_string())?;
        if cancel.load(Ordering::Relaxed) {
            if let Some(mut process) = child.take() {
                let _ = process.kill();
                let _ = process.wait();
            }
            return Err("xTB frequency calculation cancelled".to_string());
        }
        let Some(process) = child.as_mut() else {
            return Err("xTB frequency calculation cancelled".to_string());
        };
        if let Some(status) = process
            .try_wait()
            .map_err(|e| format!("Failed while waiting for xTB: {e}"))?
        {
            child.take();
            break status;
        }
        drop(child);
        std::thread::sleep(Duration::from_millis(50));
    };
    if !status.success() {
        return Err(
            "xTB frequency calculation failed (see xtb.stderr in the run directory)".to_string(),
        );
    }

    // The geometry xTB actually built the Hessian in: `--hess` does not move
    // the structure chemically, but xTB reprints it re-oriented to its own
    // centre-of-mass/principal-axis frame. Using the on-screen geometry
    // instead of this one silently corrupts the Eckart projection --
    // measured directly while building this feature: several frequencies
    // came out up to 12% wrong from exactly this mismatch.
    //
    // `xtbhess.xyz` was tried first for this and turned out to be written
    // only for some point groups (confirmed: a plain `--hess` run on a
    // larger, Cs-symmetry molecule produces no `xtbhess.xyz` at all, while
    // the identical command on a small, symmetric one does) -- reported
    // directly by the user as "xTB did not write xtbhess.xyz" on their own
    // real structure. `g98.out`'s "Standard orientation" block is written
    // in every case tested and is the geometry used instead.
    let g98_text = fs::read_to_string(workdir.join("g98.out"))
        .map_err(|_| "xTB did not write g98.out".to_string())?;
    let (atoms, coords_angstrom_flat) =
        super::hessian_file::parse_g98_geometry(&g98_text)?;
    let coords_bohr: Vec<f64> = coords_angstrom_flat
        .into_iter()
        .map(|c| c * ANGSTROM_TO_BOHR)
        .collect();

    let hessian_text = fs::read_to_string(workdir.join("hessian"))
        .map_err(|_| "xTB did not write a hessian file".to_string())?;
    let hessian = parse_turbomole_hessian(&hessian_text, atoms.len())?;

    let vibspectrum_text = fs::read_to_string(workdir.join("vibspectrum"))
        .map_err(|_| "xTB did not write vibspectrum".to_string())?;
    let vib_lines = parse_vibspectrum(&vibspectrum_text)?;

    let gradient_bohr = if need_gradient {
        let text = fs::read_to_string(workdir.join("gradient"))
            .map_err(|_| "xTB did not write a gradient file".to_string())?;
        Some(parse_turbomole_gradient(&text, atoms.len())?)
    } else {
        None
    };

    Ok(RawFrequencyOutput {
        atoms,
        coords_bohr,
        hessian,
        vib_lines,
        gradient_bohr,
    })
}

/// True when every atom lies on one line through the first two atoms, to
/// within a small tolerance -- linear molecules have 5 low modes instead of
/// the usual 6 (no rotation about the molecular axis to project out).
fn is_linear(coords_bohr: &[f64], natoms: usize) -> bool {
    if natoms < 3 {
        return true;
    }
    let at = |i: usize| Vec3::new(
        coords_bohr[3 * i] as f32,
        coords_bohr[3 * i + 1] as f32,
        coords_bohr[3 * i + 2] as f32,
    );
    let origin = at(0);
    let axis = (at(1) - origin).normalize_or_zero();
    if axis.length() < 1.0e-6 {
        return true;
    }
    (2..natoms).all(|i| {
        let v = at(i) - origin;
        let v_len = v.length();
        v_len < 1.0e-6 || axis.cross(v).length() / v_len < 1.0e-3
    })
}

/// Turns a raw xTB run into the full analysis: normal modes, thermochemistry,
/// and IR intensities matched from `vibspectrum`. Kept separate from the
/// background task so it can be unit tested directly against fixtures.
fn analyze(
    raw: RawFrequencyOutput,
    multiplicity: i32,
    eckart_requested: EckartMode,
    thermo_temp_k: f64,
    thermo_freq_cutoff_cm1: f64,
) -> Result<FrequencyResult, String> {
    let natoms = raw.atoms.len();
    let mass: Vec<f64> = raw
        .atoms
        .iter()
        .map(|s| atomic_mass_amu(s).ok_or_else(|| format!("Unknown element {s:?}")))
        .collect::<Result<_, _>>()?;
    let linear = is_linear(&raw.coords_bohr, natoms);

    let mut eckart_used = eckart_requested;
    let mut reaction_path_fallback_note = None;
    let grad_cart = if eckart_requested == EckartMode::ReactionPath {
        match &raw.gradient_bohr {
            Some(g) => {
                let norm: f64 = g.iter().map(|v| v * v).sum::<f64>().sqrt();
                if norm < REACTION_PATH_GRADIENT_FLOOR {
                    eckart_used = EckartMode::VibRot;
                    reaction_path_fallback_note = Some(format!(
                        "The gradient norm ({norm:.2e} Eh/bohr) is essentially zero -- this \\
                         geometry is already a stationary point, so there is no reaction path \\
                         to project out. Falling back to the ordinary Eckart projection \
                         (translations and rotations removed, vibrations only)."
                    ));
                    None
                } else {
                    Some(g.as_slice())
                }
            }
            None => {
                eckart_used = EckartMode::VibRot;
                reaction_path_fallback_note = Some(
                    "No gradient was computed for this run; falling back to the ordinary \
                     Eckart projection (translations and rotations removed, vibrations only)."
                        .to_string(),
                );
                None
            }
        }
    } else {
        None
    };

    let result = normal_modes_with_projection(
        &mass,
        &raw.hessian,
        linear,
        Some(&raw.coords_bohr),
        eckart_used,
        grad_cart,
    )
    .map_err(|e| e.to_string())?;

    let frequencies_cm1 = thermofuncs::freqs_au_to_cm1(&result.frequencies_au);

    // xTB's own `vibspectrum` lists exactly the same 3N modes in the same
    // ascending order our own eigensolver produces (both are "lowest to
    // highest"), so intensities line up by position.
    let ir_intensities_km_mol: Vec<f64> = if raw.vib_lines.len() == frequencies_cm1.len() {
        raw.vib_lines.iter().map(|l| l.ir_intensity_km_mol).collect()
    } else {
        vec![0.0; frequencies_cm1.len()]
    };

    let freqs_for_thermo: Vec<f64> = result
        .positive_indices
        .iter()
        .map(|&i| frequencies_cm1[i])
        .collect();
    let coords_for_brot: Vec<[f64; 3]> = (0..natoms)
        .map(|i| {
            [
                raw.coords_bohr[3 * i],
                raw.coords_bohr[3 * i + 1],
                raw.coords_bohr[3 * i + 2],
            ]
        })
        .collect();
    let brot_cm1 = thermofuncs::brot_from_coords(&coords_for_brot, &mass);
    let mass_total: f64 = mass.iter().sum();
    const ONE_ATM_PASCAL: f64 = 101_325.0;
    let thermo = thermofuncs::eval_thermo(
        &freqs_for_thermo,
        &brot_cm1,
        mass_total,
        multiplicity as f64,
        thermo_temp_k,
        ONE_ATM_PASCAL,
        thermo_freq_cutoff_cm1,
    );
    let coords_angstrom: Vec<Vec3> = (0..natoms)
        .map(|i| {
            Vec3::new(
                raw.coords_bohr[3 * i] as f32,
                raw.coords_bohr[3 * i + 1] as f32,
                raw.coords_bohr[3 * i + 2] as f32,
            ) * BOHR_TO_ANGSTROM
        })
        .collect();

    Ok(FrequencyResult {
        atoms: raw.atoms,
        coords_angstrom,
        modes: result.modes,
        frequencies_cm1,
        ir_intensities_km_mol,
        negative_indices: result.negative_indices,
        zero_indices: result.zero_indices,
        positive_indices: result.positive_indices,
        thermo,
        thermo_temp_k,
        thermo_freq_cutoff_cm1,
        eckart_used,
        reaction_path_fallback_note,
    })
}

/// Builds one period of animation frames for mode `mode_index`, following
/// the same `q_eq + L * A * cos(theta)` construction `Smite`'s
/// `normalmode/nmodeprint.py::print_single_mode` uses. `amplitude_angstrom`
/// scales the mode's own eigenvector after normalizing it so its largest
/// single-atom Cartesian displacement equals 1 -- a predictable, visible
/// amplitude regardless of the eigenvector's own (physically meaningful but
/// visually arbitrary) normalization.
pub fn animate_mode(
    atoms: &[String],
    coords_angstrom: &[Vec3],
    modes: &Array2<f64>,
    mode_index: usize,
    amplitude_angstrom: f32,
    frame_count: usize,
) -> Vec<TrajectoryFrame> {
    let natoms = atoms.len();
    let mut max_disp = 0.0f32;
    let mut displacement = vec![Vec3::ZERO; natoms];
    for i in 0..natoms {
        let d = Vec3::new(
            modes[(3 * i, mode_index)] as f32,
            modes[(3 * i + 1, mode_index)] as f32,
            modes[(3 * i + 2, mode_index)] as f32,
        );
        displacement[i] = d;
        max_disp = max_disp.max(d.length());
    }
    if max_disp < 1.0e-12 {
        max_disp = 1.0;
    }

    let frame_count = frame_count.max(2);
    (0..frame_count)
        .map(|k| {
            let theta = 2.0 * std::f32::consts::PI * (k as f32) / (frame_count as f32);
            let scale = amplitude_angstrom * theta.cos() / max_disp;
            let pos = (0..natoms)
                .map(|i| coords_angstrom[i] + displacement[i] * scale)
                .collect();
            TrajectoryFrame {
                atoms: atoms.to_vec(),
                pos,
            }
        })
        .collect()
}

/// Loads a mode's animation into the trajectory player, looping (a
/// vibration genuinely repeats, unlike an optimization path) and
/// deliberately overwriting whatever trajectory was already loaded -- the
/// caller is responsible for having told the user this will happen.
pub fn load_mode_animation(
    traj: &mut TrajectoryState,
    atoms: &[String],
    coords_angstrom: &[Vec3],
    modes: &Array2<f64>,
    mode_index: usize,
    amplitude_angstrom: f32,
    speed_fps: f32,
) {
    traj.frames = animate_mode(atoms, coords_angstrom, modes, mode_index, amplitude_angstrom, 60);
    traj.current_frame = 0;
    traj.mode = PlaybackMode::Loop;
    traj.direction = 1;
    traj.playing = true;
    traj.fps = speed_fps;
    traj.accum = 0.0;
    traj.last_applied = None;
    traj.overlay_dirty = true;
}

pub fn poll_xtb_frequencies(
    mut freq_task: ResMut<XtbFrequencyTask>,
    panel_state: Res<XtbFreqPanelState>,
) {
    let Some(task) = freq_task.task.as_mut() else {
        return;
    };
    let Some(result) = block_on(future::poll_once(task)) else {
        return;
    };
    freq_task.task = None;
    freq_task.cancel = None;
    freq_task.child = None;
    freq_task.started_at = None;
    let run_dir = freq_task.run_dir.take();
    let multiplicity = freq_task.multiplicity_used;

    match result {
        Ok(raw) => {
            if let Some(dir) = &run_dir {
                discard_xtb_run_dir(dir);
            }
            match analyze(
                raw,
                multiplicity,
                panel_state.eckart_mode,
                panel_state.thermo_temp_k,
                panel_state.thermo_freq_cutoff_cm1,
            ) {
                Ok(analyzed) => {
                    let note = analyzed
                        .reaction_path_fallback_note
                        .clone()
                        .map(|n| format!(" {n}"))
                        .unwrap_or_default();
                    freq_task.last_message = Some(format!(
                        "xTB frequency calculation succeeded ({} modes).{note}",
                        analyzed.frequencies_cm1.len()
                    ));
                    freq_task.last_is_error = false;
                    freq_task.result = Some(analyzed);
                }
                Err(err) => {
                    freq_task.last_message = Some(format!("Frequency analysis failed: {err}"));
                    freq_task.last_is_error = true;
                }
            }
        }
        Err(err) => {
            freq_task.last_message = Some(match &run_dir {
                Some(dir) => format!("{err} (see {})", dir.display()),
                None => err,
            });
            freq_task.last_is_error = true;
        }
    }
}

#[derive(Resource)]
pub struct XtbFreqPanelState {
    pub eckart_mode: EckartMode,
    pub thermo_temp_k: f64,
    pub thermo_freq_cutoff_cm1: f64,
    /// The mode currently loaded into the trajectory player, whether it is
    /// actively playing or paused -- this is what a row's Animate/Stop
    /// button and the amplitude/speed controls act on.
    pub selected_mode: Option<usize>,
    pub mode_amplitude_angstrom: f32,
    pub mode_speed_fps: f32,
    pub thermo_window_open: bool,
    pub spectrum_window_open: bool,
    pub spectrum_broadening: BroadeningKind,
    pub spectrum_width_cm1: f64,
}

impl Default for XtbFreqPanelState {
    fn default() -> Self {
        Self {
            eckart_mode: EckartMode::VibRot,
            thermo_temp_k: 298.15,
            thermo_freq_cutoff_cm1: 100.0,
            selected_mode: None,
            mode_amplitude_angstrom: 0.12,
            mode_speed_fps: 60.0,
            thermo_window_open: false,
            spectrum_window_open: false,
            spectrum_broadening: BroadeningKind::Gaussian,
            spectrum_width_cm1: 20.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BroadeningKind {
    Gaussian,
    Lorentzian,
    Voigt,
    /// No broadening at all: a vertical line at each predicted frequency,
    /// height equal to its raw IR intensity -- the theoretical stick
    /// spectrum, as opposed to a simulated experimental one.
    Sticks,
}

impl BroadeningKind {
    fn is_sticks(self) -> bool {
        matches!(self, BroadeningKind::Sticks)
    }
}

/// A single line-shape value at `x`, centred on `x0`, with intensity
/// `height` and characteristic `width` (interpreted per shape: Gaussian
/// uses it as sigma, Lorentzian as the half-width-at-half-maximum gamma,
/// Voigt as both via the pseudo-Voigt mixing below).
fn lineshape(kind: BroadeningKind, x: f64, x0: f64, height: f64, width: f64) -> f64 {
    let width = width.max(1.0e-6);
    match kind {
        BroadeningKind::Gaussian => {
            let z = (x - x0) / width;
            height * (-0.5 * z * z).exp()
        }
        BroadeningKind::Lorentzian => {
            let z = (x - x0) / width;
            height / (1.0 + z * z)
        }
        BroadeningKind::Voigt => {
            // Pseudo-Voigt: a linear mix of Gaussian and Lorentzian with the
            // same width, which is what most spectroscopy software calls a
            // practical "Voigt" profile without the cost of the true
            // Voigt convolution integral.
            0.5 * lineshape(BroadeningKind::Gaussian, x, x0, height, width)
                + 0.5 * lineshape(BroadeningKind::Lorentzian, x, x0, height, width)
        }
        // Sticks have no continuous curve -- callers branch on
        // `BroadeningKind::is_sticks` before ever reaching this function.
        BroadeningKind::Sticks => 0.0,
    }
}

/// Samples the broadened IR spectrum at `n_points` evenly spaced points
/// across `[x_min, x_max]`, summing every mode's line shape.
/// Formats a broadened spectrum as two space-separated columns (frequency,
/// intensity), one point per line, sampled at a fixed resolution rather than
/// the display's own point count -- so the exported `.dat` does not depend
/// on how large the plot window happened to be when exported.
pub fn format_spectrum_dat(
    frequencies_cm1: &[f64],
    intensities_km_mol: &[f64],
    kind: BroadeningKind,
    width_cm1: f64,
    x_min: f64,
    x_max: f64,
    resolution_cm1: f64,
) -> String {
    let resolution_cm1 = resolution_cm1.max(1.0e-6);
    let n_points = (((x_max - x_min) / resolution_cm1).round() as usize).max(1) + 1;
    let curve = sample_spectrum(
        frequencies_cm1,
        intensities_km_mol,
        kind,
        width_cm1,
        x_min,
        x_max,
        n_points,
    );
    let mut out = String::new();
    for (x, y) in curve {
        out.push_str(&format!("{x:.4} {y:.8}\n"));
    }
    out
}

/// Exports the theoretical stick spectrum as-is: one line per predicted
/// mode (frequency, IR intensity), the same pairs the mode list and the
/// stick plot itself show -- no resampling grid, since there is no
/// continuous curve to resample.
pub fn format_spectrum_sticks(frequencies_cm1: &[f64], intensities_km_mol: &[f64]) -> String {
    let mut out = String::new();
    for (&f, &i) in frequencies_cm1.iter().zip(intensities_km_mol) {
        out.push_str(&format!("{f:.4} {i:.8}\n"));
    }
    out
}

pub fn sample_spectrum(
    frequencies_cm1: &[f64],
    intensities_km_mol: &[f64],
    kind: BroadeningKind,
    width_cm1: f64,
    x_min: f64,
    x_max: f64,
    n_points: usize,
) -> Vec<(f64, f64)> {
    let n_points = n_points.max(2);
    (0..n_points)
        .map(|i| {
            let x = x_min + (x_max - x_min) * i as f64 / (n_points - 1) as f64;
            let y = frequencies_cm1
                .iter()
                .zip(intensities_km_mol)
                .map(|(&f, &inten)| lineshape(kind, x, f, inten, width_cm1))
                .sum();
            (x, y)
        })
        .collect()
}

/// Draws the "Frequency Analysis (xTB Hessian)" section: Eckart-mode choice,
/// temperature/cutoff for thermochemistry, the Run/Cancel controls, the mode
/// list with an animate button per row, and the Show Thermochemistry / Show
/// IR Spectrum buttons. Called from `ui.rs`, mirroring
/// `xtb_optimize::xtb_optimization_panel`.
pub fn xtb_frequency_panel(
    ui: &mut egui::Ui,
    freq_panel: &mut XtbFreqPanelState,
    freq_task: &mut XtbFrequencyTask,
    opt_panel: &super::xtb_optimize::XtbPanelState,
    mol: &Molecule,
    traj: &mut TrajectoryState,
) {
    let running = freq_task.is_running();

    // The Eckart projection -- which removes translations and rotations
    // from the Hessian, leaving vibrations only -- is always applied; that
    // is not an optional choice for a normal molecule, so there is nothing
    // to toggle for it. (Behemoth's own enum for this is still named
    // `EckartMode::VibRot`, imported as-is.) The one real decision
    // is whether to *additionally* project out the reaction-coordinate
    // direction too, which only makes sense at a transition state or other
    // non-minimum extremum, since it needs an actual (non-zero) gradient.
    ui.horizontal(|ui| {
        let mut reaction_path = freq_panel.eckart_mode == EckartMode::ReactionPath;
        if ui
            .add_enabled(
                !running,
                egui::Checkbox::new(
                    &mut reaction_path,
                    "Reaction-path projection (transition states / extrema only)",
                ),
            )
            .changed()
        {
            freq_panel.eckart_mode = if reaction_path {
                EckartMode::ReactionPath
            } else {
                EckartMode::VibRot
            };
        }
    });

    ui.horizontal(|ui| {
        ui.label("Temperature (K)");
        ui.add_enabled(
            !running,
            egui::DragValue::new(&mut freq_panel.thermo_temp_k).range(1.0..=2000.0),
        );
        ui.label("Cutoff (cm⁻¹)");
        ui.add_enabled(
            !running,
            egui::DragValue::new(&mut freq_panel.thermo_freq_cutoff_cm1).range(0.0..=500.0),
        );
    });

    ui.horizontal(|ui| {
        let can_run = !running && !mol.atoms.is_empty();
        let mut button = ui.add_enabled(can_run, egui::Button::new("Run Frequencies"));
        if mol.atoms.is_empty() {
            button = button.on_disabled_hover_text("There is no structure on screen to analyze.");
        }
        if button.clicked() {
            match super::xtb_optimize::resolve_xtb_executable(&opt_panel.xtb_path) {
                Ok(xtb_path) => {
                    let need_gradient = freq_panel.eckart_mode == EckartMode::ReactionPath;
                    freq_task.start(
                        &xtb_path,
                        &mol.atoms,
                        &mol.pos,
                        opt_panel.charge,
                        (opt_panel.multiplicity - 1).max(0),
                        opt_panel.multiplicity,
                        need_gradient,
                    );
                }
                Err(err) => {
                    freq_task.last_message = Some(err);
                    freq_task.last_is_error = true;
                }
            }
        }
        if running && ui.button("Cancel").clicked() {
            freq_task.cancel_now();
        }
    });

    if running {
        let elapsed = freq_task.elapsed().unwrap_or_default();
        ui.weak(format!(
            "Running… {} (a Hessian calculation has no intermediate progress \
             to report; this is one single-point-like step, not a series of \
             iterations)",
            super::xtb_optimize::format_elapsed(elapsed)
        ));
    }

    if let Some(message) = &freq_task.last_message {
        let color = if freq_task.last_is_error {
            egui::Color32::from_rgb(220, 80, 80)
        } else {
            egui::Color32::from_rgb(80, 180, 100)
        };
        ui.colored_label(color, message);
    }

    let Some(result) = &freq_task.result else {
        return;
    };

    ui.separator();
    let eckart_label = match result.eckart_used {
        EckartMode::Off => "no projection",
        EckartMode::VibRot => "translations/rotations projected out",
        EckartMode::ReactionPath => "reaction-path projection",
    };
    ui.label(format!(
        "{} modes ({} low, {} imaginary) — {eckart_label}",
        result.frequencies_cm1.len(),
        result.zero_indices.len(),
        result.negative_indices.len()
    ));

    egui::ScrollArea::vertical()
        .max_height(180.0)
        .show(ui, |ui| {
            egui::Grid::new("xtb_freq_mode_list")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("Mode");
                    ui.strong("Frequency (cm⁻¹)");
                    ui.strong("IR (km/mol)");
                    ui.strong("");
                    ui.end_row();

                    for i in 0..result.frequencies_cm1.len() {
                        let freq = result.frequencies_cm1[i];
                        let is_imaginary = result.negative_indices.contains(&i);
                        let is_low = result.zero_indices.contains(&i);
                        let label = if is_imaginary {
                            format!("{:.1} i", freq.abs())
                        } else {
                            format!("{freq:.1}")
                        };
                        let color = if is_imaginary {
                            egui::Color32::from_rgb(220, 80, 80)
                        } else if is_low {
                            ui.visuals().weak_text_color()
                        } else {
                            ui.visuals().text_color()
                        };
                        ui.label(format!("{}", i + 1));
                        ui.colored_label(color, label);
                        ui.label(format!("{:.2}", result.ir_intensities_km_mol[i]));
                        // Animate/Stop is a single toggle, not two separate
                        // buttons: at most one mode animates at a time (they
                        // all share the one trajectory player), so the
                        // button's own label is exactly "is this the row
                        // that's currently playing?"
                        let is_animating_this_row =
                            freq_panel.selected_mode == Some(i) && traj.playing;
                        let toggle_label = if is_animating_this_row { "Stop" } else { "Animate" };
                        if ui.button(toggle_label).clicked() {
                            if is_animating_this_row {
                                traj.playing = false;
                            } else {
                                freq_panel.selected_mode = Some(i);
                                load_mode_animation(
                                    traj,
                                    &result.atoms,
                                    &result.coords_angstrom,
                                    &result.modes,
                                    i,
                                    freq_panel.mode_amplitude_angstrom,
                                    freq_panel.mode_speed_fps,
                                );
                            }
                        }
                        ui.end_row();
                    }
                });
        });

    if freq_panel.selected_mode.is_some() {
        // These apply to whichever mode is currently loaded, whether it is
        // playing or paused (Stop leaves the frames in place), so both
        // controls stay live: nudging amplitude or speed while stopped sets
        // up the next Animate press; nudging while playing updates it
        // immediately.
        ui.horizontal(|ui| {
            ui.label("Amplitude (Å)");
            let amplitude_changed = ui
                .add(egui::DragValue::new(&mut freq_panel.mode_amplitude_angstrom)
                    .range(0.01..=50.0)
                    .speed(0.01))
                .changed();
            ui.label("Speed (fps)");
            let speed_changed = ui
                .add(egui::DragValue::new(&mut freq_panel.mode_speed_fps)
                    .range(1.0..=200.0)
                    .speed(0.5))
                .changed();
            if (amplitude_changed || speed_changed) && traj.playing {
                if let Some(i) = freq_panel.selected_mode {
                    load_mode_animation(
                        traj,
                        &result.atoms,
                        &result.coords_angstrom,
                        &result.modes,
                        i,
                        freq_panel.mode_amplitude_angstrom,
                        freq_panel.mode_speed_fps,
                    );
                }
            } else if speed_changed {
                // Paused: just remember the new speed for the next Animate
                // press rather than restarting playback on its own.
                traj.fps = freq_panel.mode_speed_fps;
            }
        });
        ui.weak(
            "Animating a mode replaces whatever trajectory was loaded (e.g. an \
             optimization path) and loops, since a vibration genuinely repeats.",
        );
    }

    ui.horizontal(|ui| {
        let label = if freq_panel.thermo_window_open {
            "Hide Thermochemistry"
        } else {
            "Show Thermochemistry"
        };
        if ui.button(label).clicked() {
            freq_panel.thermo_window_open = !freq_panel.thermo_window_open;
        }
        let label = if freq_panel.spectrum_window_open {
            "Hide IR Spectrum"
        } else {
            "Show IR Spectrum"
        };
        if ui.button(label).clicked() {
            freq_panel.spectrum_window_open = !freq_panel.spectrum_window_open;
        }
    });

    let ctx = ui.ctx().clone();
    thermochemistry_window(&ctx, &mut freq_panel.thermo_window_open, result);
    ir_spectrum_window(&ctx, freq_panel, result);
}

fn thermochemistry_window(ctx: &egui::Context, open: &mut bool, result: &FrequencyResult) {
    if !*open {
        return;
    }
    egui::Window::new("Thermochemistry")
        .open(open)
        .resizable(true)
        .default_size([460.0, 420.0])
        .show(ctx, |ui| {
            egui::ScrollArea::both().show(ui, |ui| {
                let text = thermofuncs::format_thermo(
                    &result.thermo,
                    result.thermo_temp_k,
                    result.thermo_freq_cutoff_cm1,
                    None,
                );
                ui.add(egui::Label::new(
                    egui::RichText::new(text).font(egui::FontId::monospace(13.0)),
                ));
            });
        });
}

fn ir_spectrum_window(
    ctx: &egui::Context,
    freq_panel: &mut XtbFreqPanelState,
    result: &FrequencyResult,
) {
    if !freq_panel.spectrum_window_open {
        return;
    }
    egui::Window::new("IR Spectrum")
        .open(&mut freq_panel.spectrum_window_open)
        .resizable(true)
        .default_size([500.0, 340.0])
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Line shape");
                ui.selectable_value(
                    &mut freq_panel.spectrum_broadening,
                    BroadeningKind::Gaussian,
                    "Gaussian",
                );
                ui.selectable_value(
                    &mut freq_panel.spectrum_broadening,
                    BroadeningKind::Lorentzian,
                    "Lorentzian",
                );
                ui.selectable_value(
                    &mut freq_panel.spectrum_broadening,
                    BroadeningKind::Voigt,
                    "Voigt",
                );
                ui.selectable_value(
                    &mut freq_panel.spectrum_broadening,
                    BroadeningKind::Sticks,
                    "Sticks",
                );
                // A stick spectrum has no width to speak of -- it is drawn
                // (and exported) at the exact predicted frequencies.
                ui.add_enabled_ui(!freq_panel.spectrum_broadening.is_sticks(), |ui| {
                    ui.label("Width (cm⁻¹)");
                    ui.add(
                        egui::DragValue::new(&mut freq_panel.spectrum_width_cm1)
                            .range(0.5..=200.0),
                    );
                });
            });

            let vib_freqs: Vec<f64> = result
                .positive_indices
                .iter()
                .map(|&i| result.frequencies_cm1[i])
                .collect();
            let vib_intens: Vec<f64> = result
                .positive_indices
                .iter()
                .map(|&i| result.ir_intensities_km_mol[i])
                .collect();
            if vib_freqs.is_empty() {
                ui.weak("No vibrational modes to plot.");
                return;
            }
            // Fixed axis convention for an IR spectrum, per the user's own
            // choice: 0 up to the highest predicted frequency plus a fixed
            // 100 cm^-1 margin, rather than a window that tracks the current
            // broadening width and would shift as it is adjusted.
            let x_min = 0.0f64;
            let x_max = vib_freqs.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + 100.0;

            if ui.button("Export .dat…").clicked() {
                // Sticks export the raw predicted (frequency, intensity)
                // pairs as-is -- there is no continuous curve to resample,
                // so the fixed 0.5 cm^-1 grid used for the broadened shapes
                // does not apply here.
                let text = if freq_panel.spectrum_broadening.is_sticks() {
                    format_spectrum_sticks(&vib_freqs, &vib_intens)
                } else {
                    const EXPORT_RESOLUTION_CM1: f64 = 0.5;
                    format_spectrum_dat(
                        &vib_freqs,
                        &vib_intens,
                        freq_panel.spectrum_broadening,
                        freq_panel.spectrum_width_cm1,
                        x_min,
                        x_max,
                        EXPORT_RESOLUTION_CM1,
                    )
                };
                if let Some(path) = rfd::FileDialog::new()
                    .set_file_name("ir_spectrum.dat")
                    .add_filter("Data file", &["dat", "txt"])
                    .save_file()
                {
                    if let Err(err) = std::fs::write(&path, text) {
                        eprintln!("Failed to write {}: {err}", path.display());
                    }
                }
            }

            let is_sticks = freq_panel.spectrum_broadening.is_sticks();
            // Sticks have no continuous curve to sample; the plot's y-scale
            // then comes straight from the raw intensities instead of a
            // broadened curve's peak height.
            let curve = if is_sticks {
                Vec::new()
            } else {
                sample_spectrum(
                    &vib_freqs,
                    &vib_intens,
                    freq_panel.spectrum_broadening,
                    freq_panel.spectrum_width_cm1,
                    x_min,
                    x_max,
                    400,
                )
            };

            const AXIS_LABEL_SIZE: f32 = 13.0;
            const AXIS_TITLE_SIZE: f32 = 16.0;
            let max_y = if is_sticks {
                vib_intens.iter().cloned().fold(0.0f64, f64::max).max(1.0e-9)
            } else {
                curve
                    .iter()
                    .map(|&(_, y)| y)
                    .fold(0.0f64, f64::max)
                    .max(1.0e-9)
            };

            let desired_size = egui::vec2(
                ui.available_width().max(200.0),
                ui.available_height().max(120.0),
            );
            let (rect, response) = ui.allocate_exact_size(desired_size, egui::Sense::hover());
            let painter = ui.painter_at(rect);
            // Fixed white background regardless of the app's own theme, with
            // dark text/gridlines to stay legible against it -- a
            // consistent look independent of light/dark mode.
            painter.rect_filled(rect, 0.0, egui::Color32::WHITE);
            let plot_color = egui::Color32::from_rgb(30, 30, 30);
            let grid_color = plot_color.gamma_multiply(0.35);
            let line_color = egui::Color32::from_rgb(40, 100, 170);

            let plot_rect = egui::Rect::from_min_max(
                rect.min + egui::vec2(50.0, 10.0),
                rect.max - egui::vec2(10.0, 46.0),
            );

            let to_screen = |x: f64, y: f64| -> egui::Pos2 {
                let sx = plot_rect.left()
                    + ((x - x_min) / (x_max - x_min)) as f32 * plot_rect.width();
                let sy = plot_rect.bottom() - (y / max_y) as f32 * plot_rect.height();
                egui::pos2(sx, sy)
            };

            // Ticks at fixed 500 cm^-1 intervals (500, 1000, 1500, ...) up
            // to the top of the axis, rather than a fixed count of evenly
            // spaced but arbitrarily-valued ticks.
            const TICK_STEP_CM1: f64 = 500.0;
            let mut tick = TICK_STEP_CM1;
            while tick <= x_max {
                let sx = plot_rect.left()
                    + ((tick - x_min) / (x_max - x_min)) as f32 * plot_rect.width();
                painter.line_segment(
                    [egui::pos2(sx, plot_rect.top()), egui::pos2(sx, plot_rect.bottom())],
                    egui::Stroke::new(1.0, grid_color),
                );
                painter.text(
                    egui::pos2(sx, plot_rect.bottom() + 6.0),
                    egui::Align2::CENTER_TOP,
                    format!("{tick:.0}"),
                    egui::FontId::monospace(AXIS_LABEL_SIZE),
                    plot_color,
                );
                tick += TICK_STEP_CM1;
            }
            painter.rect_stroke(
                plot_rect,
                0.0,
                egui::Stroke::new(1.5, plot_color.gamma_multiply(0.6)),
                egui::StrokeKind::Outside,
            );

            let peak_points: Vec<egui::Pos2>;
            if is_sticks {
                // One vertical bar per predicted mode, from the baseline up
                // to its own raw intensity -- the theoretical stick
                // spectrum, not a simulated continuous one.
                peak_points = vib_freqs
                    .iter()
                    .zip(&vib_intens)
                    .map(|(&f, &i)| to_screen(f, i))
                    .collect();
                let baseline_y = plot_rect.bottom();
                for &p in &peak_points {
                    painter.line_segment(
                        [egui::pos2(p.x, baseline_y), p],
                        egui::Stroke::new(2.0, line_color),
                    );
                }
            } else {
                let points: Vec<egui::Pos2> =
                    curve.iter().map(|&(x, y)| to_screen(x, y)).collect();
                if points.len() >= 2 {
                    painter.add(egui::Shape::line(points, egui::Stroke::new(2.0, line_color)));
                }

                // Peak markers, one per vibrational mode, drawn at the
                // curve's own height at that frequency (not the mode's raw
                // intensity), so a hover near a peak lands on what is
                // actually visible.
                peak_points = vib_freqs
                    .iter()
                    .map(|&f| {
                        let y = vib_freqs
                            .iter()
                            .zip(&vib_intens)
                            .map(|(&f2, &i2)| {
                                lineshape(
                                    freq_panel.spectrum_broadening,
                                    f,
                                    f2,
                                    i2,
                                    freq_panel.spectrum_width_cm1,
                                )
                            })
                            .sum();
                        to_screen(f, y)
                    })
                    .collect();
                for &p in &peak_points {
                    painter.circle_filled(p, 2.5, line_color);
                }
            }

            if let Some(hover_pos) = response.hover_pos() {
                if let Some((idx, dist)) = peak_points
                    .iter()
                    .enumerate()
                    .map(|(i, p)| (i, p.distance(hover_pos)))
                    .min_by(|a, b| a.1.total_cmp(&b.1))
                {
                    const HOVER_RADIUS: f32 = 14.0;
                    if dist <= HOVER_RADIUS {
                        // Fixed dark colour, matching the plot's own fixed
                        // white background rather than the app's theme --
                        // the theme's own "strong text" is often light,
                        // which would vanish against white in dark mode.
                        let highlight_color = plot_color;
                        painter.circle_stroke(
                            peak_points[idx],
                            5.0,
                            egui::Stroke::new(2.0, highlight_color),
                        );
                        let label = format!(
                            "{:.1} cm⁻¹\n{:.2} km/mol",
                            vib_freqs[idx], vib_intens[idx]
                        );
                        let galley = painter.layout_no_wrap(
                            label,
                            egui::FontId::proportional(AXIS_LABEL_SIZE),
                            highlight_color,
                        );
                        let padding = egui::vec2(6.0, 4.0);
                        let box_size = galley.size() + padding * 2.0;
                        let mut box_pos = peak_points[idx] + egui::vec2(12.0, -box_size.y - 12.0);
                        box_pos.x = box_pos.x.clamp(rect.left(), rect.right() - box_size.x);
                        box_pos.y = box_pos.y.clamp(rect.top(), rect.bottom() - box_size.y);
                        let box_rect = egui::Rect::from_min_size(box_pos, box_size);
                        painter.rect_filled(box_rect, 4.0, egui::Color32::from_rgb(240, 240, 240));
                        painter.rect_stroke(
                            box_rect,
                            4.0,
                            egui::Stroke::new(1.0, plot_color.gamma_multiply(0.5)),
                            egui::StrokeKind::Outside,
                        );
                        painter.galley(box_rect.min + padding, galley, highlight_color);
                    }
                }
            }

            painter.text(
                egui::pos2(plot_rect.center().x, rect.bottom() - 4.0),
                egui::Align2::CENTER_BOTTOM,
                "Wavenumber (cm⁻¹)",
                egui::FontId::proportional(AXIS_TITLE_SIZE),
                plot_color,
            );
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join(name);
        std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("cannot read fixture {}: {e}", path.display()))
    }

    fn fixture_example(name: &str) -> String {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join(name);
        std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("cannot read example {}: {e}", path.display()))
    }

    fn h2o2_raw() -> RawFrequencyOutput {
        let (_n, atoms, coords_angstrom) =
            parse_xyz_angstrom(&fixture("h2o2_hessian_geometry.xyz"));
        let coords_bohr: Vec<f64> = coords_angstrom
            .into_iter()
            .flatten()
            .map(|c| c * ANGSTROM_TO_BOHR)
            .collect();
        let hessian =
            parse_turbomole_hessian(&fixture("h2o2_hessian"), atoms.len()).unwrap();
        let vib_lines = parse_vibspectrum(&fixture("h2o2_vibspectrum")).unwrap();
        RawFrequencyOutput {
            atoms,
            coords_bohr,
            hessian,
            vib_lines,
            gradient_bohr: None,
        }
    }

    /// The full `analyze` pipeline -- normal modes, IR-intensity matching,
    /// and thermochemistry all at once -- against the real captured
    /// fixtures, cross-checked against xTB's own independently-computed
    /// frequencies and zero-point energy.
    #[test]
    fn analyze_reproduces_xtbs_frequencies_and_zpe() {
        let result = analyze(h2o2_raw(), 1, EckartMode::VibRot, 298.15, 100.0).unwrap();

        assert_eq!(result.zero_indices.len(), 6);
        assert_eq!(result.positive_indices.len(), 6);
        assert!(result.negative_indices.is_empty());
        assert_eq!(result.eckart_used, EckartMode::VibRot);
        assert!(result.reaction_path_fallback_note.is_none());

        let mut computed: Vec<f64> = result
            .positive_indices
            .iter()
            .map(|&i| result.frequencies_cm1[i])
            .collect();
        computed.sort_by(|a, b| a.total_cmp(b));
        let expected = [259.56, 1108.86, 1158.96, 1355.45, 3531.45, 3534.68];
        for (c, e) in computed.iter().zip(&expected) {
            assert!((c - e).abs() / e < 0.02, "{c} vs {e}");
        }

        // IR intensities travel with their own frequency, matched by index
        // position -- the loudest computed mode should be the loudest real one.
        let loudest_idx = result
            .positive_indices
            .iter()
            .copied()
            .max_by(|&a, &b| {
                result.ir_intensities_km_mol[a].total_cmp(&result.ir_intensities_km_mol[b])
            })
            .unwrap();
        assert!(
            (result.frequencies_cm1[loudest_idx] - 259.56).abs() < 5.0,
            "expected the 259.56 cm^-1 mode (413 km/mol) to be the loudest, got {} cm^-1",
            result.frequencies_cm1[loudest_idx]
        );

        let xtb_zpe = 0.024943563139;
        assert!((result.thermo.zpe - xtb_zpe).abs() / xtb_zpe < 0.02, "{}", result.thermo.zpe);
    }

    /// A reaction-path request with no gradient at all (rather than a small
    /// one) must also fall back cleanly instead of panicking or erroring.
    #[test]
    fn reaction_path_without_a_gradient_falls_back_to_vibrot() {
        let result = analyze(h2o2_raw(), 1, EckartMode::ReactionPath, 298.15, 100.0).unwrap();
        assert_eq!(result.eckart_used, EckartMode::VibRot);
        assert!(result.reaction_path_fallback_note.is_some());
    }

    /// A genuinely large gradient (not near a stationary point) must let
    /// the reaction-path projection actually run rather than falling back.
    #[test]
    fn reaction_path_with_a_real_gradient_is_not_downgraded() {
        let mut raw = h2o2_raw();
        // A synthetic, clearly non-zero gradient -- far above the
        // near-stationary-point floor this is meant to distinguish from.
        raw.gradient_bohr = Some(vec![0.05; 12]);
        let result = analyze(raw, 1, EckartMode::ReactionPath, 298.15, 100.0).unwrap();
        assert_eq!(result.eckart_used, EckartMode::ReactionPath);
        assert!(result.reaction_path_fallback_note.is_none());
    }

    /// Reproduces the exact bug caught while building this feature: feeding
    /// the Eckart projector a geometry that does not match the one the
    /// Hessian was computed at (the pre-optimization structure, while
    /// `--ohess` computed the Hessian post-optimization) silently corrupts
    /// several frequencies rather than erroring -- which is precisely why
    /// the real feature always reads `xtbhess.xyz` instead of reusing
    /// whatever geometry was on screen.
    #[test]
    fn a_mismatched_geometry_produces_visibly_wrong_frequencies() {
        let mut raw = h2o2_raw();
        let (_n, _atoms, wrong_coords_angstrom) =
            parse_xyz_angstrom(&fixture_example("h2o2.xyz"));
        raw.coords_bohr = wrong_coords_angstrom
            .into_iter()
            .flatten()
            .map(|c| c * ANGSTROM_TO_BOHR)
            .collect();

        let result = analyze(raw, 1, EckartMode::VibRot, 298.15, 100.0).unwrap();
        let mut computed: Vec<f64> = result
            .positive_indices
            .iter()
            .map(|&i| result.frequencies_cm1[i])
            .collect();
        computed.sort_by(|a, b| a.total_cmp(b));
        let expected = [259.56, 1108.86, 1158.96, 1355.45, 3531.45, 3534.68];
        let worst_rel_err = computed
            .iter()
            .zip(&expected)
            .map(|(c, e)| (c - e).abs() / e)
            .fold(0.0f64, f64::max);
        assert!(
            worst_rel_err > 0.05,
            "expected the mismatched geometry to visibly corrupt at least one \
             frequency (as it did when this bug was first found), but the worst \
             relative error was only {:.2}%",
            worst_rel_err * 100.0
        );
    }

    /// The exact bug the user reported: `xtb --hess` on a real-sized,
    /// lower-symmetry structure (`examples/phenyl-OCH3.xyz`, Cs point group)
    /// does not write `xtbhess.xyz` at all, which the first version of this
    /// feature depended on unconditionally. Runs the real binary end to end
    /// and checks the analysis still succeeds using `g98.out` instead.
    /// Skipped gracefully when no xTB binary is configured.
    #[test]
    fn a_real_run_on_a_larger_asymmetric_molecule_does_not_need_xtbhess_xyz() {
        let configured = super::super::xtb_optimize::default_xtb_executable();
        let Ok(xtb_path) = super::super::xtb_optimize::resolve_xtb_executable(&configured) else {
            eprintln!("skipping: no xTB executable found (set BEAVYR_XTB_PATH to run this test)");
            return;
        };

        let (_n, atoms, coords_angstrom) =
            parse_xyz_angstrom(&fixture_example("phenyl-OCH3.xyz"));
        let pos: Vec<Vec3> = coords_angstrom
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let input_xyz = super::super::xtb_optimize::write_xyz_string(&atoms, &pos);

        let cancel = AtomicBool::new(false);
        let child_slot = Mutex::new(None);
        let workdir = xtb_scratch_dir().join(format!("test_freq_asym_{}", std::process::id()));

        let raw = run_xtb_hess_cancellable(
            &xtb_path, &workdir, &input_xyz, 0, 0, false, &cancel, &child_slot,
        )
        .expect("a plain --hess run on this real molecule should succeed");
        assert_eq!(raw.atoms.len(), 16);
        assert!(
            !workdir.join("xtbhess.xyz").exists(),
            "this test fixture was specifically chosen because xTB does not \
             write xtbhess.xyz for it -- if it now exists, xTB's behaviour \
             changed and this test no longer exercises the bug it was written for"
        );

        let result = analyze(raw, 1, EckartMode::VibRot, 298.15, 100.0).unwrap();
        assert_eq!(result.atoms.len(), 16);
        assert!(result.frequencies_cm1.iter().any(|&f| f.abs() > 100.0));

        let _ = fs::remove_dir_all(&workdir);
    }

    /// Runs the real xTB binary end to end: `--hess`, parse its output, run
    /// the full analysis, and check the result against xTB's own frequencies.
    /// Skipped gracefully when no xTB binary is configured.
    #[test]
    fn a_real_xtb_hess_run_reproduces_its_own_vibspectrum() {
        let configured = super::super::xtb_optimize::default_xtb_executable();
        let Ok(xtb_path) = super::super::xtb_optimize::resolve_xtb_executable(&configured) else {
            eprintln!("skipping: no xTB executable found (set BEAVYR_XTB_PATH to run this test)");
            return;
        };

        let (_n, atoms, coords_angstrom) = parse_xyz_angstrom(&fixture_example("h2o2.xyz"));
        let pos: Vec<Vec3> = coords_angstrom
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();
        let input_xyz = super::super::xtb_optimize::write_xyz_string(&atoms, &pos);

        let cancel = AtomicBool::new(false);
        let child_slot = Mutex::new(None);
        let workdir = xtb_scratch_dir().join(format!("test_freq_{}", std::process::id()));

        let raw = run_xtb_hess_cancellable(
            &xtb_path, &workdir, &input_xyz, 0, 0, false, &cancel, &child_slot,
        )
        .expect("a plain --hess run on a stable closed-shell molecule should succeed");
        assert_eq!(raw.atoms.len(), 4);
        assert_eq!(raw.hessian.dim(), (12, 12));

        let result = analyze(raw, 1, EckartMode::VibRot, 298.15, 100.0).unwrap();
        // `examples/h2o2.xyz` is a hand-placed, not fully relaxed geometry
        // (unlike the `--ohess`-derived fixture the other test in this file
        // uses), so its six low modes are only approximately zero rather
        // than the near-machine-zero of a true stationary point -- one can
        // land just outside the classification tolerance either way. The
        // strict "exactly 6 clean low modes" check belongs on the converged
        // geometry (`analyze_reproduces_xtbs_frequencies_and_zpe`, which
        // does assert it); this real-run test only checks the pipeline
        // itself works end to end against the live binary.
        assert!(
            result.positive_indices.len() >= 5,
            "{:?}",
            result.positive_indices
        );

        let mut computed: Vec<f64> = result
            .positive_indices
            .iter()
            .map(|&i| result.frequencies_cm1[i])
            .collect();
        computed.sort_by(|a, b| a.total_cmp(b));
        // `--hess` at the unoptimized geometry gives different (still real)
        // frequencies than the `--ohess` fixture's post-optimization ones,
        // so this only checks internal sanity: all positive, ascending, and
        // in a chemically plausible O-H/O-O stretch range at the top end.
        assert!(computed.iter().all(|&f| f > 0.0));
        assert!(computed.windows(2).all(|w| w[0] <= w[1]));
        assert!(*computed.last().unwrap() > 2000.0, "{computed:?}");

        let _ = fs::remove_dir_all(&workdir);
    }

    /// Drives the exact same resource/task machinery the real app uses --
    /// `XtbFrequencyTask::start`, then repeatedly `poll_xtb_frequencies` in
    /// a real Bevy `App` update loop -- against a molecule loaded from
    /// `examples/h2o2.xyz`, the same way the panel would from `Molecule`.
    /// Skipped gracefully when no xTB binary is configured.
    #[test]
    fn the_full_bevy_task_pipeline_produces_a_result_for_the_loaded_molecule() {
        let configured = super::super::xtb_optimize::default_xtb_executable();
        let Ok(xtb_path) = super::super::xtb_optimize::resolve_xtb_executable(&configured) else {
            eprintln!("skipping: no xTB executable found (set BEAVYR_XTB_PATH to run this test)");
            return;
        };

        let (_n, atoms, coords_angstrom) = parse_xyz_angstrom(&fixture_example("h2o2.xyz"));
        let pos: Vec<Vec3> = coords_angstrom
            .into_iter()
            .map(|v| Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
            .collect();

        // Outside a full Bevy `App`, nothing initializes the task pool that
        // `TaskPoolPlugin` normally sets up at startup -- do it directly.
        AsyncComputeTaskPool::get_or_init(bevy::tasks::TaskPool::default);

        let mut task = XtbFrequencyTask::default();
        task.start(&xtb_path, &atoms, &pos, 0, 0, 1, false);
        assert!(task.is_running(), "start() should have spawned a task");

        let panel_state = XtbFreqPanelState::default();
        let deadline = Instant::now() + Duration::from_secs(60);
        loop {
            // Mirrors `poll_xtb_frequencies`'s own body exactly, since that
            // function takes Bevy `ResMut`/`Res` params that are awkward to
            // construct standalone -- the logic under test is identical.
            if let Some(t) = task.task.as_mut() {
                if let Some(result) = block_on(future::poll_once(t)) {
                    task.task = None;
                    task.cancel = None;
                    task.child = None;
                    task.started_at = None;
                    let run_dir = task.run_dir.take();
                    let multiplicity = task.multiplicity_used;
                    match result {
                        Ok(raw) => {
                            if let Some(dir) = &run_dir {
                                discard_xtb_run_dir(dir);
                            }
                            match analyze(
                                raw,
                                multiplicity,
                                panel_state.eckart_mode,
                                panel_state.thermo_temp_k,
                                panel_state.thermo_freq_cutoff_cm1,
                            ) {
                                Ok(analyzed) => task.result = Some(analyzed),
                                Err(err) => panic!("analyze failed: {err}"),
                            }
                        }
                        Err(err) => panic!("xtb run failed: {err}"),
                    }
                    break;
                }
            }
            assert!(Instant::now() < deadline, "frequency task did not finish within 60s");
            std::thread::sleep(Duration::from_millis(50));
        }

        let result = task.result.expect("the task must have produced a result");
        assert_eq!(result.atoms, atoms, "the analyzed atoms must match the loaded molecule");
        assert!(!result.frequencies_cm1.is_empty());
        assert!(result.positive_indices.len() >= 5, "{:?}", result.positive_indices);
    }

    #[test]
    fn is_linear_detects_a_straight_triatomic() {
        let coords = vec![
            0.0, 0.0, 0.0,
            0.0, 0.0, 2.0,
            0.0, 0.0, -2.0,
        ];
        assert!(is_linear(&coords, 3));
    }

    #[test]
    fn is_linear_rejects_a_bent_triatomic() {
        let coords = vec![
            0.0, 0.0, 0.0,
            1.5, 0.0, 0.0,
            0.0, 1.5, 0.0,
        ];
        assert!(!is_linear(&coords, 3));
    }

    #[test]
    fn animate_mode_starts_at_the_equilibrium_geometry() {
        let atoms = vec!["H".to_string(), "H".to_string()];
        let coords = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.74)];
        let mut modes = Array2::<f64>::zeros((6, 6));
        // A pure stretch along z for both atoms, opposite signs.
        modes[(2, 0)] = 1.0;
        modes[(5, 0)] = -1.0;
        let frames = animate_mode(&atoms, &coords, &modes, 0, 0.3, 8);
        // theta = 0 -> cos(0) = 1 -> maximum displacement, not equilibrium;
        // frame 0 is the amplitude peak by this construction (matching
        // Smite's own convention, where the animation starts at +amplitude).
        assert_eq!(frames.len(), 8);
        assert_eq!(frames[0].atoms, atoms);
    }

    #[test]
    fn animate_mode_respects_the_requested_amplitude() {
        let atoms = vec!["H".to_string(), "H".to_string()];
        let coords = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.74)];
        let mut modes = Array2::<f64>::zeros((6, 6));
        modes[(2, 0)] = 1.0;
        modes[(5, 0)] = -1.0;
        let frames = animate_mode(&atoms, &coords, &modes, 0, 0.5, 4);
        // Frame 0: theta=0, cos=1, so displacement = amplitude exactly.
        let disp = (frames[0].pos[0] - coords[0]).length();
        assert!((disp - 0.5).abs() < 1.0e-5, "{disp}");
    }

    #[test]
    fn animate_mode_normalizes_by_the_largest_atomic_displacement() {
        // One atom moves twice as far as the other in the raw eigenvector;
        // after normalization the amplitude is set by the *larger* mover,
        // so it (not the smaller one) reaches exactly the requested amplitude.
        let atoms = vec!["H".to_string(), "H".to_string()];
        let coords = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.74)];
        let mut modes = Array2::<f64>::zeros((6, 6));
        modes[(2, 0)] = 2.0;
        modes[(5, 0)] = -1.0;
        let frames = animate_mode(&atoms, &coords, &modes, 0, 0.4, 4);
        let disp0 = (frames[0].pos[0] - coords[0]).length();
        let disp1 = (frames[0].pos[1] - coords[1]).length();
        assert!((disp0 - 0.4).abs() < 1.0e-5, "{disp0}");
        assert!((disp1 - 0.2).abs() < 1.0e-5, "{disp1}");
    }

    #[test]
    fn load_mode_animation_sets_speed_amplitude_and_loops() {
        let atoms = vec!["H".to_string(), "H".to_string()];
        let coords = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.74)];
        let mut modes = Array2::<f64>::zeros((6, 6));
        modes[(2, 0)] = 1.0;
        modes[(5, 0)] = -1.0;
        let mut traj = TrajectoryState::default();
        traj.mode = PlaybackMode::Once; // as an optimization would have left it
        traj.playing = false;

        load_mode_animation(&mut traj, &atoms, &coords, &modes, 0, 0.4, 25.0);

        assert_eq!(traj.mode, PlaybackMode::Loop, "a vibration genuinely repeats");
        assert!(traj.playing);
        assert_eq!(traj.fps, 25.0);
        assert_eq!(traj.current_frame, 0);
        assert!(!traj.frames.is_empty());
    }

    #[test]
    fn format_spectrum_dat_uses_space_separated_two_columns_at_the_requested_resolution() {
        let text = format_spectrum_dat(
            &[500.0],
            &[10.0],
            BroadeningKind::Gaussian,
            5.0,
            0.0,
            10.0,
            0.5,
        );
        let lines: Vec<&str> = text.lines().collect();
        // 0.0..=10.0 at 0.5 resolution: 21 points.
        assert_eq!(lines.len(), 21);
        for line in &lines {
            let cols: Vec<&str> = line.split(' ').collect();
            assert_eq!(cols.len(), 2, "expected exactly two space-separated columns: {line:?}");
            assert!(cols[0].parse::<f64>().is_ok());
            assert!(cols[1].parse::<f64>().is_ok());
        }
        // Points land exactly on the 0.5 cm^-1 grid: 0.0, 0.5, 1.0, ...
        let x0: f64 = lines[0].split(' ').next().unwrap().parse().unwrap();
        let x1: f64 = lines[1].split(' ').next().unwrap().parse().unwrap();
        assert!((x0 - 0.0).abs() < 1e-9);
        assert!((x1 - 0.5).abs() < 1e-9);
    }

    #[test]
    fn format_spectrum_dat_reflects_the_chosen_lineshape() {
        let gaussian = format_spectrum_dat(
            &[500.0], &[10.0], BroadeningKind::Gaussian, 5.0, 495.0, 505.0, 1.0,
        );
        let lorentzian = format_spectrum_dat(
            &[500.0], &[10.0], BroadeningKind::Lorentzian, 5.0, 495.0, 505.0, 1.0,
        );
        assert_ne!(
            gaussian, lorentzian,
            "different line shapes must export different curves"
        );
    }

    #[test]
    fn format_spectrum_sticks_prints_the_raw_predicted_pairs_unresampled() {
        let text = format_spectrum_sticks(&[259.56, 1108.86, 3534.68], &[413.48, 0.00066, 101.16]);
        let lines: Vec<&str> = text.lines().collect();
        // One line per mode, not a resampled grid.
        assert_eq!(lines.len(), 3);
        for line in &lines {
            let cols: Vec<&str> = line.split(' ').collect();
            assert_eq!(cols.len(), 2);
        }
        assert!(lines[0].starts_with("259.5"));
        assert!(lines[2].starts_with("3534.6"));
    }

    #[test]
    fn format_spectrum_sticks_of_no_modes_is_empty() {
        assert_eq!(format_spectrum_sticks(&[], &[]), "");
    }

    #[test]
    fn sticks_lineshape_is_never_used_for_a_continuous_curve() {
        // Sticks are handled by a dedicated code path that never samples a
        // curve at all; `lineshape` itself just needs to stay exhaustive
        // and inert for this variant rather than panic.
        assert_eq!(lineshape(BroadeningKind::Sticks, 500.0, 500.0, 10.0, 5.0), 0.0);
    }

    #[test]
    fn is_sticks_identifies_only_the_sticks_variant() {
        assert!(BroadeningKind::Sticks.is_sticks());
        assert!(!BroadeningKind::Gaussian.is_sticks());
        assert!(!BroadeningKind::Lorentzian.is_sticks());
        assert!(!BroadeningKind::Voigt.is_sticks());
    }

    #[test]
    fn gaussian_peak_matches_its_height_at_the_centre() {
        let y = lineshape(BroadeningKind::Gaussian, 500.0, 500.0, 10.0, 5.0);
        assert!((y - 10.0).abs() < 1.0e-9);
    }

    #[test]
    fn lorentzian_peak_matches_its_height_at_the_centre() {
        let y = lineshape(BroadeningKind::Lorentzian, 500.0, 500.0, 10.0, 5.0);
        assert!((y - 10.0).abs() < 1.0e-9);
    }

    #[test]
    fn spectrum_sampling_places_a_visible_peak_near_each_mode() {
        let freqs = [500.0, 1500.0];
        let intens = [10.0, 20.0];
        let curve = sample_spectrum(&freqs, &intens, BroadeningKind::Gaussian, 10.0, 0.0, 2000.0, 400);
        let (peak_x, peak_y) = curve.iter().cloned().fold((0.0, f64::MIN), |acc, p| {
            if p.1 > acc.1 { p } else { acc }
        });
        assert!((peak_x - 1500.0).abs() < 10.0, "peak at {peak_x}, expected near 1500");
        assert!(peak_y > 15.0);
    }

    #[test]
    fn spectrum_of_no_modes_is_flat_zero() {
        let curve = sample_spectrum(&[], &[], BroadeningKind::Gaussian, 10.0, 0.0, 100.0, 5);
        assert!(curve.iter().all(|&(_, y)| y == 0.0));
    }
}
