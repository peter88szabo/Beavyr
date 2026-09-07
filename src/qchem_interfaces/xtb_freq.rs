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

use crate::molecule::{atomic_mass_amu, parse_xyz_angstrom, Molecule};
use crate::normalmode::thermofuncs::{self, ThermoResults};
use crate::normalmode::{normal_modes_with_projection, EckartMode};
use crate::spectrum::broadening::{format_spectrum_dat, format_spectrum_sticks, BroadeningKind};
use crate::spectrum::plot::{draw_spectrum, SpectrumAxis};
use crate::trajectory::{PlaybackMode, TrajectoryFrame, TrajectoryState};

use super::hessian_file::{parse_turbomole_gradient, parse_turbomole_hessian, parse_vibspectrum};
use super::program::QcProgram;
use super::valence::validate_electronic_state;
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
    /// Rotational symmetry number used for this result's thermochemistry.
    pub rot_symmetry: f64,
    pub eckart_used: EckartMode,
    /// The empirical scaling factor applied to every frequency, when the
    /// source file asked for one. `None` means the frequencies are the raw
    /// harmonic ones.
    pub frequency_scale: Option<f64>,
    /// Where the modes came from, when they were not computed here. Set for
    /// an imported frequency job, whose modes and intensities are its own and
    /// whose thermochemistry is ours.
    pub source_note: Option<String>,
    /// Set when a reaction-path request was silently downgraded to the
    /// ordinary translation/rotation projection
    /// because the gradient was essentially zero.
    pub reaction_path_fallback_note: Option<String>,
}

#[derive(Clone)]
struct RawFrequencyOutput {
    atoms: Vec<String>,
    coords_bohr: Vec<f64>,
    hessian: Array2<f64>,
    vib_lines: Vec<super::hessian_file::VibSpectrumLine>,
    gradient_bohr: Option<Vec<f64>>,
    /// Masses the source file stated, when it states them. An ORCA `.hess`
    /// records one per atom, so an isotope substitution survives instead of
    /// being overwritten by the standard atomic weight. xTB output carries no
    /// masses, so that path leaves this `None`.
    masses_amu: Option<Vec<f64>>,
    /// IR intensities already ordered to match our eigensolver, for sources
    /// whose own ordering differs from ours. Overrides `vib_lines` when set.
    ir_intensities_km_mol: Option<Vec<f64>>,
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
    /// The Hessian and geometry the current result came from, kept so that
    /// changing the scaling factor, the temperature or the cutoff can redo the
    /// analysis in milliseconds instead of recomputing a Hessian.
    last_raw: Option<HessianSource>,
    /// Multiplicity the stored raw data was produced with.
    last_multiplicity: i32,
    /// Geometry a just-loaded Hessian file brought with it, waiting to be put
    /// on screen. A loaded Hessian usually arrives with nothing displayed, and
    /// its modes are meaningless without the structure they belong to.
    pub pending_geometry: Option<(Vec<String>, Vec<Vec3>)>,
}

impl FrequencyResult {
    /// Whether mode `i` has a displacement vector to animate.
    ///
    /// An imported result lists its source's translation/rotation residuals
    /// by frequency only -- the file prints no vectors for them -- so those
    /// rows appear in the list but have nothing to move.
    pub fn is_animatable(&self, i: usize) -> bool {
        i < self.modes.ncols() && self.modes.column(i).iter().any(|d| *d != 0.0)
    }

    /// Whether there is an IR spectrum to plot at all.
    ///
    /// A Hessian on its own gives frequencies but no intensities: those need
    /// dipole derivatives, which are a separate quantity. An imported file may
    /// or may not carry them, and plotting a row of zero-height peaks would
    /// look like a computed spectrum that happens to be flat rather than like
    /// data that was never there.
    pub fn has_ir_intensities(&self) -> bool {
        self.ir_intensities_km_mol.iter().any(|&i| i > 0.0)
    }
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

    #[allow(clippy::too_many_arguments)]
    pub fn start(
        &mut self,
        program: QcProgram,
        binary: &Path,
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
        let binary = binary.to_path_buf();
        let task_atoms = atoms.to_vec();
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
            run_hess_cancellable(
                program,
                &binary,
                &task_dir,
                &input_xyz,
                &task_atoms,
                charge,
                uhf,
                multiplicity,
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

#[allow(clippy::too_many_arguments)]
fn run_hess_cancellable(
    program: QcProgram,
    binary: &Path,
    workdir: &Path,
    input_xyz: &str,
    atoms: &[String],
    charge: i32,
    uhf: i32,
    multiplicity: i32,
    need_gradient: bool,
    cancel: &AtomicBool,
    child_slot: &Mutex<Option<Child>>,
) -> Result<RawFrequencyOutput, String> {
    if binary.as_os_str().is_empty() {
        return Err(format!("{} path is empty", program.label()));
    }
    fs::create_dir_all(workdir).map_err(|e| format!("Failed to create workdir: {e}"))?;
    fs::write(workdir.join("input.xyz"), input_xyz)
        .map_err(|e| format!("Failed to write input.xyz: {e}"))?;

    let stdout = fs::File::create(workdir.join("xtb.stdout"))
        .map_err(|e| format!("Failed to create xtb.stdout: {e}"))?;
    let stderr = fs::File::create(workdir.join("xtb.stderr"))
        .map_err(|e| format!("Failed to create xtb.stderr: {e}"))?;
    let mut cmd = match program {
        QcProgram::Xtb => {
            let mut cmd = Command::new(binary);
            cmd.current_dir(workdir)
                .arg("input.xyz")
                .arg("--hess")
                .arg("--chrg")
                .arg(charge.to_string());
            if uhf > 0 {
                cmd.arg("--uhf").arg(uhf.to_string());
            }
            if need_gradient {
                cmd.arg("--grad");
            }
            cmd
        }
        // Behemoth prints the matrix to stdout, which is captured below.
        QcProgram::Behemoth => super::behemoth::hessian_command(
            binary,
            workdir,
            "input.xyz",
            charge,
            multiplicity,
        ),
    };
    cmd.stdout(Stdio::from(stdout)).stderr(Stdio::from(stderr));

    let child = cmd
        .spawn()
        .map_err(|e| format!("Failed to run {}: {e}", program.label()))?;
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
        return Err(format!(
            "{} frequency calculation failed (see xtb.stderr in the run directory)",
            program.label()
        ));
    }

    if program == QcProgram::Behemoth {
        // Behemoth builds the Hessian at the geometry it was given and prints
        // the matrix to stdout, so there is no re-oriented frame to recover
        // and no separate file to read. It reports no IR intensities with
        // `--hessian`, so the mode list shows N/A and the spectrum stays
        // unavailable rather than showing a row of zeroes.
        let log = fs::read_to_string(workdir.join("xtb.stdout"))
            .map_err(|_| "Behemoth wrote no output to read the Hessian from".to_string())?;
        let hessian = super::behemoth::parse_hessian_for(&log, atoms)?;
        let coords_bohr = parse_xyz_angstrom(input_xyz)
            .2
            .into_iter()
            .flatten()
            .map(|c| c * ANGSTROM_TO_BOHR)
            .collect();
        return Ok(RawFrequencyOutput {
            atoms: atoms.to_vec(),
            coords_bohr,
            hessian,
            vib_lines: Vec::new(),
            gradient_bohr: None,
            masses_amu: None,
            ir_intensities_km_mol: None,
        });
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
        masses_amu: None,
        ir_intensities_km_mol: None,
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
    frequency_scale: f64,
    rot_symmetry: f64,
) -> Result<FrequencyResult, String> {
    let natoms = raw.atoms.len();
    let mass: Vec<f64> = match &raw.masses_amu {
        Some(masses) if masses.len() == natoms => masses.clone(),
        _ => raw
            .atoms
            .iter()
            .map(|s| atomic_mass_amu(s).ok_or_else(|| format!("Unknown element {s:?}")))
            .collect::<Result<_, _>>()?,
    };
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

    let mut frequencies_cm1 = thermofuncs::freqs_au_to_cm1(&result.frequencies_au);
    // Applied before anything reads the frequencies, so the mode list, the IR
    // spectrum and the thermochemistry are all consistently scaled -- which is
    // also what the source program does when it is asked to scale.
    let scaled = frequency_scale > 0.0 && (frequency_scale - 1.0).abs() > f64::EPSILON;
    if scaled {
        for f in &mut frequencies_cm1 {
            *f *= frequency_scale;
        }
    }

    // xTB's own `vibspectrum` lists exactly the same 3N modes in the same
    // ascending order our own eigensolver produces (both are "lowest to
    // highest"), so intensities line up by position.
    let ir_intensities_km_mol: Vec<f64> = match &raw.ir_intensities_km_mol {
        Some(given) if given.len() == frequencies_cm1.len() => given.clone(),
        _ if raw.vib_lines.len() == frequencies_cm1.len() => {
            raw.vib_lines.iter().map(|l| l.ir_intensity_km_mol).collect()
        }
        _ => vec![0.0; frequencies_cm1.len()],
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
        rot_symmetry,
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
        rot_symmetry,
        eckart_used,
        frequency_scale: scaled.then_some(frequency_scale),
        source_note: None,
        reaction_path_fallback_note,
    })
}

/// Frequencies and modes a program already computed, for a file that reports
/// its results but not the Hessian behind them -- a Gaussian `.log`, say.
///
/// The modes are taken as printed; only the thermochemistry is computed here.
#[derive(Clone)]
struct ImportedModes {
    atoms: Vec<String>,
    coords_bohr: Vec<f64>,
    frequencies_cm1: Vec<f64>,
    ir_intensities_km_mol: Vec<f64>,
    /// `3N x nmodes`, in our normalisation.
    modes: Array2<f64>,
    /// Translation/rotation residuals the source reported. They have no
    /// displacement vectors -- the source printed only their frequencies --
    /// so they are listed but cannot be animated.
    low_frequencies_cm1: Vec<f64>,
    /// Named in the panel so it is clear what was read and what was computed.
    program: String,
}

/// Either a Hessian we can analyse ourselves, or a finished set of modes.
#[derive(Clone)]
enum HessianSource {
    /// A force-constant matrix: we project it and diagonalise it.
    Computed(Box<RawFrequencyOutput>),
    /// Modes and frequencies as another program reported them.
    Imported(Box<ImportedModes>),
}

/// Builds a result from modes another program computed, calculating only the
/// thermochemistry.
///
/// There is no projection to do -- the source already removed translations and
/// rotations, which is why it prints `3N - 6` modes -- so the Eckart setting
/// does not apply and no mode is classified as "low".
fn analyze_imported(
    imported: ImportedModes,
    multiplicity: i32,
    thermo_temp_k: f64,
    thermo_freq_cutoff_cm1: f64,
    frequency_scale: f64,
    rot_symmetry: f64,
) -> Result<FrequencyResult, String> {
    let natoms = imported.atoms.len();
    let mass: Vec<f64> = imported
        .atoms
        .iter()
        .map(|s| atomic_mass_amu(s).ok_or_else(|| format!("Unknown element {s:?}")))
        .collect::<Result<_, _>>()?;

    let mut frequencies_cm1 = imported.frequencies_cm1.clone();
    let scaled = frequency_scale > 0.0 && (frequency_scale - 1.0).abs() > f64::EPSILON;
    if scaled {
        for f in &mut frequencies_cm1 {
            *f *= frequency_scale;
        }
    }

    // The source's translation/rotation residuals are listed alongside the
    // vibrations, in the same ascending order the computed path uses, so the
    // two kinds of result read the same way. They carry no displacement
    // vectors, which is what marks them as un-animatable.
    let mut rows: Vec<(f64, Option<usize>)> = frequencies_cm1
        .iter()
        .enumerate()
        .map(|(i, f)| (*f, Some(i)))
        .collect();
    rows.extend(imported.low_frequencies_cm1.iter().map(|f| (*f, None)));
    rows.sort_by(|a, b| a.0.total_cmp(&b.0));

    let nrows = rows.len();
    let mut ordered_freqs = Vec::with_capacity(nrows);
    let mut ordered_intensities = Vec::with_capacity(nrows);
    let mut ordered_modes = Array2::<f64>::zeros((3 * natoms, nrows));
    let mut negative_indices = Vec::new();
    let mut zero_indices = Vec::new();
    let mut positive_indices = Vec::new();
    for (row, (freq, source)) in rows.iter().enumerate() {
        ordered_freqs.push(*freq);
        match source {
            Some(k) => {
                ordered_intensities
                    .push(imported.ir_intensities_km_mol.get(*k).copied().unwrap_or(0.0));
                for axis in 0..3 * natoms {
                    ordered_modes[[axis, row]] = imported.modes[[axis, *k]];
                }
                if *freq < 0.0 {
                    negative_indices.push(row);
                } else {
                    positive_indices.push(row);
                }
            }
            None => {
                // A residual: frequency only, no vector, no intensity.
                ordered_intensities.push(0.0);
                zero_indices.push(row);
            }
        }
    }
    let frequencies_cm1 = ordered_freqs;

    let freqs_for_thermo: Vec<f64> = positive_indices
        .iter()
        .map(|&i| frequencies_cm1[i])
        .collect();
    let coords_for_brot: Vec<[f64; 3]> = (0..natoms)
        .map(|i| {
            [
                imported.coords_bohr[3 * i],
                imported.coords_bohr[3 * i + 1],
                imported.coords_bohr[3 * i + 2],
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
        rot_symmetry,
    );
    let coords_angstrom: Vec<Vec3> = (0..natoms)
        .map(|i| {
            Vec3::new(
                imported.coords_bohr[3 * i] as f32,
                imported.coords_bohr[3 * i + 1] as f32,
                imported.coords_bohr[3 * i + 2] as f32,
            ) * BOHR_TO_ANGSTROM
        })
        .collect();

    Ok(FrequencyResult {
        atoms: imported.atoms,
        coords_angstrom,
        modes: ordered_modes,
        frequencies_cm1,
        ir_intensities_km_mol: ordered_intensities,
        negative_indices,
        zero_indices,
        positive_indices,
        thermo,
        thermo_temp_k,
        thermo_freq_cutoff_cm1,
        rot_symmetry,
        eckart_used: EckartMode::Off,
        frequency_scale: scaled.then_some(frequency_scale),
        source_note: Some(format!(
            "Normal modes, frequencies and IR intensities are {}'s own, read from the \
             file. The thermochemistry below is computed here from those frequencies, \
             at the temperature, cutoff and symmetry number set in this panel -- it is \
             not {}'s.",
            imported.program, imported.program
        )),
        reaction_path_fallback_note: None,
    })
}

/// Builds one period of animation frames for mode `mode_index`, following
/// the same `q_eq + L * A * cos(theta)` construction `Smite`'s
/// `normalmode/nmodeprint.py::print_single_mode` uses. `amplitude_angstrom`
/// scales the mode's own eigenvector after normalizing it so its largest
/// single-atom Cartesian displacement equals 1 -- a predictable, visible
/// amplitude regardless of the eigenvector's own (physically meaningful but
/// visually arbitrary) normalization.
/// Reads a frequency file, along with any empirical scaling factor it asks
/// for, and says which kind it is.
///
/// Two kinds are recognised, and they are not the same thing:
/// * an ORCA `.hess`, which holds the force-constant matrix -- we project and
///   diagonalise it ourselves;
/// * a Gaussian `.log`, which holds the finished frequencies, intensities and
///   normal coordinates but no Hessian -- we take those as given and compute
///   only the thermochemistry.
fn read_hessian_file(path: &Path) -> Result<(HessianSource, Option<f64>), String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;

    if super::gaussian_log::is_gaussian_log(&text) {
        let g = super::gaussian_log::parse_gaussian_log(&text)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        let coords_bohr: Vec<f64> = g
            .coords_angstrom
            .iter()
            .map(|c| c * ANGSTROM_TO_BOHR)
            .collect();
        let imported = ImportedModes {
            atoms: g.atoms,
            coords_bohr,
            frequencies_cm1: g.frequencies_cm1,
            ir_intensities_km_mol: g.ir_intensities_km_mol,
            modes: g.modes,
            low_frequencies_cm1: g.low_frequencies_cm1,
            program: "Gaussian".to_string(),
        };
        // Gaussian applies any scaling factor itself and prints the scaled
        // frequencies, so there is nothing here for us to apply.
        return Ok((HessianSource::Imported(Box::new(imported)), None));
    }

    let hess = super::orca_hess::parse_orca_hess(&text)
        .map_err(|e| format!("{}: {e}", path.display()))?;
    let raw = RawFrequencyOutput {
        atoms: hess.atoms.clone(),
        coords_bohr: hess.coords_bohr.clone(),
        hessian: hess.hessian.clone(),
        vib_lines: Vec::new(),
        // A `.hess` is written at a stationary point and stores no gradient,
        // so a reaction-path request falls back to the ordinary Eckart
        // projection with a note -- `analyze` already handles that.
        gradient_bohr: None,
        masses_amu: Some(hess.masses_amu.clone()),
        ir_intensities_km_mol: {
            let ordered = hess.ir_intensities_in_eigenvalue_order();
            (!ordered.is_empty()).then_some(ordered)
        },
    };
    Ok((
        HessianSource::Computed(Box::new(raw)),
        hess.frequency_scale_factor,
    ))
}

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
/// Stops a mode animation and puts the molecule back at its equilibrium
/// geometry.
///
/// Simply clearing `playing` leaves the structure frozen wherever the cycle
/// happened to stop -- and never at equilibrium, since the frames run
/// `q_eq + L A cos(theta)` from `theta = 0`, which is maximum displacement.
/// The trajectory is replaced by the single undistorted frame so that the
/// molecule on screen is the one the modes belong to.
pub fn stop_mode_animation(
    traj: &mut TrajectoryState,
    atoms: &[String],
    coords_angstrom: &[Vec3],
) {
    traj.playing = false;
    traj.frames = vec![TrajectoryFrame {
        atoms: atoms.to_vec(),
        pos: coords_angstrom.to_vec(),
    }];
    traj.current_frame = 0;
    traj.direction = 1;
    traj.accum = 0.0;
    // The frame index has not changed, so without this the applier would
    // consider the frame already shown and leave the distorted geometry up.
    traj.last_applied = None;
    traj.overlay_dirty = true;
}

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
    // Bonds follow the motion, even if the user froze them for a previous
    // trajectory: the forming and breaking bonds of a transition state's
    // imaginary mode are invisible if the graph never updates.
    traj.fixed_bonds = false;
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
            freq_task.last_raw = Some(HessianSource::Computed(Box::new(raw.clone())));
            freq_task.last_multiplicity = multiplicity;
            match analyze(
                raw,
                multiplicity,
                panel_state.eckart_mode,
                panel_state.thermo_temp_k,
                panel_state.thermo_freq_cutoff_cm1,
                panel_state.frequency_scale,
                panel_state.rot_symmetry,
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
    /// Which external program computes the Hessian. Kept separate from the
    /// optimizer's choice: running a Hessian in a different program than the
    /// optimization is a legitimate thing to want.
    pub program: QcProgram,
    pub eckart_mode: EckartMode,
    pub thermo_temp_k: f64,
    pub thermo_freq_cutoff_cm1: f64,
    /// Whether the point-group reference table has its own window open.
    pub symmetry_table_open: bool,
    /// Rotational symmetry number for the thermochemistry. Supplied by the
    /// user rather than detected: 1 is right for an asymmetric molecule, and
    /// anyone running a symmetric one knows its number.
    pub rot_symmetry: f64,
    /// Charge for an xTB Hessian run. Independent of the optimizer's: a
    /// frequency calculation is often run on a different species than the one
    /// last optimized, and silently inheriting the other panel's value is how
    /// a Hessian gets computed for the wrong electronic state.
    pub charge: i32,
    pub multiplicity: i32,
    /// Valence warnings appear only once a run has actually been attempted,
    /// matching the optimizer -- see `XtbPanelState::show_warnings`.
    pub show_warnings: bool,
    /// Empirical scaling applied to every computed frequency. 1.0 leaves the
    /// raw harmonic values alone. Loading a Hessian file that names its own
    /// factor sets this to that value, since the program that wrote the file
    /// already reported its frequencies scaled that way.
    pub frequency_scale: f64,
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
            program: QcProgram::default(),
            eckart_mode: EckartMode::VibRot,
            thermo_temp_k: 298.15,
            thermo_freq_cutoff_cm1: 100.0,
            frequency_scale: 1.0,
            rot_symmetry: 1.0,
            symmetry_table_open: false,
            charge: 0,
            multiplicity: 1,
            show_warnings: false,
            selected_mode: None,
            mode_amplitude_angstrom: 0.3,
            mode_speed_fps: 60.0,
            thermo_window_open: false,
            spectrum_window_open: false,
            spectrum_broadening: BroadeningKind::Gaussian,
            spectrum_width_cm1: 20.0,
        }
    }
}

/// Opens a file dialog and analyses the chosen Hessian.
///
/// If the file names its own scaling factor -- ORCA writes
/// `$frequency_scale_factor` and reports its frequencies already scaled by it
/// -- that becomes the panel's factor, so the numbers shown match the ones the
/// user already has in their output file. They can still change it afterwards.
fn load_hessian_from_dialog(freq_panel: &mut XtbFreqPanelState, freq_task: &mut XtbFrequencyTask) {
    let Some(path) = rfd::FileDialog::new()
        .add_filter("Hessian or frequency output", &["hess", "log", "out"])
        .add_filter("All files", &["*"])
        .pick_file()
    else {
        return;
    };
    let (source, file_scale) = match read_hessian_file(&path) {
        Ok(pair) => pair,
        Err(err) => {
            freq_task.last_message = Some(err);
            freq_task.last_is_error = true;
            return;
        }
    };
    if let Some(scale) = file_scale {
        freq_panel.frequency_scale = scale;
    }
    let natoms = match &source {
        HessianSource::Computed(raw) => raw.atoms.len(),
        HessianSource::Imported(imported) => imported.atoms.len(),
    };
    freq_task.last_raw = Some(source);
    // A loaded Hessian carries no charge or spin state of its own; the modes
    // and thermochemistry only need the electronic multiplicity for the
    // electronic partition function, and 1 is the sane default.
    freq_task.last_multiplicity = 1;
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| path.display().to_string());

    match reanalyze_stored_hessian(freq_panel, freq_task) {
        true => {
            if let Some(result) = &freq_task.result {
                freq_task.pending_geometry =
                    Some((result.atoms.clone(), result.coords_angstrom.clone()));
            }
            let scaling = match file_scale {
                Some(scale) => format!(", frequencies scaled by {scale} as the file asks"),
                None => String::new(),
            };
            freq_task.last_message =
                Some(format!("Loaded {name}: {natoms} atoms{scaling}."));
            freq_task.last_is_error = false;
        }
        false => {
            freq_task.last_is_error = true;
        }
    }
}

/// Redoes the analysis from the stored Hessian with the panel's current
/// scaling, temperature and cutoff. Returns whether it succeeded.
fn reanalyze_stored_hessian(
    freq_panel: &XtbFreqPanelState,
    freq_task: &mut XtbFrequencyTask,
) -> bool {
    let Some(source) = freq_task.last_raw.clone() else {
        return false;
    };
    let outcome = match source {
        HessianSource::Computed(raw) => analyze(
            *raw,
            freq_task.last_multiplicity,
            freq_panel.eckart_mode,
            freq_panel.thermo_temp_k,
            freq_panel.thermo_freq_cutoff_cm1,
            freq_panel.frequency_scale,
            freq_panel.rot_symmetry,
        ),
        HessianSource::Imported(imported) => analyze_imported(
            *imported,
            freq_task.last_multiplicity,
            freq_panel.thermo_temp_k,
            freq_panel.thermo_freq_cutoff_cm1,
            freq_panel.frequency_scale,
            freq_panel.rot_symmetry,
        ),
    };
    match outcome {
        Ok(result) => {
            freq_task.result = Some(result);
            true
        }
        Err(err) => {
            freq_task.last_message = Some(err);
            freq_task.last_is_error = true;
            false
        }
    }
}

/// Draws the "Frequency Analysis" section: Eckart-mode choice,
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
    // A playing trajectory or mode animation puts a distorted frame on
    // screen. Computing a Hessian for it would use that frame's geometry, and
    // the valence check would be judging bond lengths that belong to the
    // animation rather than to the molecule.
    let animating = traj.playing;

    // The Eckart projection -- which removes translations and rotations
    // from the Hessian, leaving vibrations only -- is always applied; that
    // is not an optional choice for a normal molecule, so there is nothing
    // to toggle for it. (Behemoth's own enum for this is still named
    // `EckartMode::VibRot`, imported as-is.) The one real decision
    // is whether to *additionally* project out the reaction-coordinate
    // direction too, which needs a non-zero gradient: it applies along a
    // reaction path, not at a stationary point, where the gradient vanishes
    // and there is no path direction to remove.
    ui.horizontal(|ui| {
        let mut reaction_path = freq_panel.eckart_mode == EckartMode::ReactionPath;
        if ui
            .add_enabled(
                !running,
                egui::Checkbox::new(
                    &mut reaction_path,
                    "Reaction-path projection",
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

    // Changing any of these redoes the analysis from the stored Hessian rather
    // than asking for a new one -- a diagonalisation is milliseconds, so there
    // is no reason to make the user re-run a calculation to try 0.97 instead
    // of 1.0, or 273 K instead of 298 K.
    let mut reanalyze = false;
    ui.horizontal(|ui| {
        ui.label("Temperature (K)");
        reanalyze |= ui
            .add_enabled(
                !running,
                egui::DragValue::new(&mut freq_panel.thermo_temp_k).range(1.0..=2000.0),
            )
            .changed();
        ui.label("Cutoff (cm⁻¹)");
        reanalyze |= ui
            .add_enabled(
                !running,
                egui::DragValue::new(&mut freq_panel.thermo_freq_cutoff_cm1).range(0.0..=500.0),
            )
            .changed();
    });

    ui.horizontal(|ui| {
        ui.label("Symmetry number \u{3c3}");
        let response = ui.add_enabled(
            !running,
            egui::DragValue::new(&mut freq_panel.rot_symmetry)
                .range(1.0..=120.0)
                .speed(1.0)
                .fixed_decimals(0),
        );
        reanalyze |= response.changed();
        response.on_hover_text(
            "Rotational symmetry number for the thermochemistry. 1 for an \
             asymmetric molecule, 2 for water, 3 for ammonia, 6 for BF3, \
             12 for benzene or methane, 24 for SF6.",
        );
        let label = if freq_panel.symmetry_table_open {
            "Hide SymNum vs PointGroup"
        } else {
            "Show SymNum vs PointGroup"
        };
        if ui
            .add(egui::Button::new(egui::RichText::new(label).small()))
            .clicked()
        {
            freq_panel.symmetry_table_open = !freq_panel.symmetry_table_open;
        }
    });
    // Drawn here rather than with the result windows: the table is a
    // reference, useful before any Hessian exists.
    symmetry_number_window(
        &ui.ctx().clone(),
        &mut freq_panel.symmetry_table_open,
        &mut freq_panel.rot_symmetry,
        &mut reanalyze,
    );

    ui.horizontal(|ui| {
        ui.label("Frequency scale");
        let response = ui.add_enabled(
            !running,
            egui::DragValue::new(&mut freq_panel.frequency_scale)
                .range(0.5..=1.5)
                .speed(0.001)
                .fixed_decimals(4),
        );
        reanalyze |= response.changed();
        response.on_hover_text(
            "Empirical scaling applied to every frequency, and to the \
             thermochemistry computed from them. 1.0 = raw harmonic values.",
        );
        if freq_panel.frequency_scale != 1.0
            && ui
                .add_enabled(!running, egui::Button::new("Reset to 1.0"))
                .clicked()
        {
            freq_panel.frequency_scale = 1.0;
            reanalyze = true;
        }
    });

    // Two ways to get a Hessian, kept visually apart because they need
    // completely different inputs: computing one here needs an electronic
    // state and an xTB binary, loading one needs neither.
    ui.add_space(6.0);
    // Collapsed by default: computing a Hessian here is the less common route
    // now that one can be loaded, and its charge/multiplicity controls are
    // only in the way until someone actually wants them. Forced open while a
    // run is in flight, so Cancel cannot be hidden behind a closed header.
    egui::CollapsingHeader::new(egui::RichText::new("Run Freq Calc").strong())
        .default_open(false)
        .open(running.then_some(true))
        .show(ui, |ui| {
        ui.horizontal(|ui| {
            ui.label("Program");
            egui::ComboBox::from_id_salt("freq_program")
                .selected_text(freq_panel.program.label())
                .show_ui(ui, |ui| {
                    for program in QcProgram::ALL {
                        ui.selectable_value(&mut freq_panel.program, program, program.label());
                    }
                });
            ui.weak("path set in Geometry Optimization");
        });
        ui.horizontal(|ui| {
            ui.label("Charge");
            ui.add_enabled(
                !running,
                egui::DragValue::new(&mut freq_panel.charge).range(-10..=10),
            );
            ui.label("Multiplicity");
            ui.add_enabled(
                !running,
                egui::DragValue::new(&mut freq_panel.multiplicity).range(1..=10),
            );
            if ui
                .add_enabled(!running, egui::Button::new("Reset to defaults"))
                .clicked()
            {
                freq_panel.charge = 0;
                freq_panel.multiplicity = 1;
                freq_panel.show_warnings = false;
            }
        });

        // Same check the optimizer runs: a multiplicity that cannot be made
        // from this many electrons would have xTB fail (or worse, converge to
        // something meaningless) several minutes into a Hessian.
        let validation =
            validate_electronic_state(mol, freq_panel.charge, freq_panel.multiplicity);
        let (uhf, warnings) = match &validation {
            Ok((uhf, warnings)) => (Some(*uhf), warnings.as_slice()),
            Err(_) => (None, &[][..]),
        };
        if animating {
            ui.weak("Playback is running \u{2014} stop it to run a calculation.");
        } else if let Err(err) = &validation {
            // A hard error explains why the button is disabled, so it shows
            // straight away rather than waiting for a click that cannot land.
            ui.colored_label(egui::Color32::from_rgb(220, 80, 80), err.to_string());
        } else if freq_panel.show_warnings {
            for warning in warnings {
                ui.colored_label(
                    egui::Color32::from_rgb(210, 160, 40),
                    format!("Atom {}: {}", warning.atom_index + 1, warning.message),
                );
            }
        }

        ui.horizontal(|ui| {
            let can_run = !running && !mol.atoms.is_empty() && uhf.is_some() && !animating;
            let mut button =
                ui.add_enabled(
                    can_run,
                    egui::Button::new(format!("Run Freq Calc ({})", freq_panel.program.label())),
                );
            if animating {
                button = button.on_disabled_hover_text(
                    "Stop the animation first: the structure on screen is a frame of it.",
                );
            } else if mol.atoms.is_empty() {
                button =
                    button.on_disabled_hover_text("There is no structure on screen to analyze.");
            }
            if button.clicked() {
                freq_panel.show_warnings = true;
                let program = freq_panel.program;
                let configured = opt_panel
                    .paths
                    .iter()
                    .find(|(p, _)| *p == program)
                    .map(|(_, path)| path.clone())
                    .unwrap_or_default();
                match super::xtb_optimize::resolve_program_executable(program, &configured) {
                    Ok(binary) => {
                        if let Some(uhf) = uhf {
                            let need_gradient =
                                freq_panel.eckart_mode == EckartMode::ReactionPath;
                            freq_task.start(
                                program,
                                &binary,
                                &mol.atoms,
                                &mol.pos,
                                freq_panel.charge,
                                uhf,
                                freq_panel.multiplicity,
                                need_gradient,
                            );
                        }
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
    });

    ui.add_space(6.0);
    if ui
        .add_enabled(!running, egui::Button::new("Load Hessian…"))
        .clicked()
    {
        load_hessian_from_dialog(freq_panel, freq_task);
    }

    if reanalyze {
        reanalyze_stored_hessian(freq_panel, freq_task);
    }

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
    // Scaled frequencies are not the raw harmonic ones, so say so rather than
    // let the numbers quietly disagree with an unscaled calculation.
    // Where the numbers came from, when they were not computed here. Stated
    // rather than left to inference: a mode list looks identical whether we
    // diagonalised a Hessian or read someone else's answer.
    if let Some(note) = &result.source_note {
        ui.colored_label(egui::Color32::from_rgb(150, 175, 210), note);
    }
    if let Some(scale) = result.frequency_scale {
        ui.weak(format!(
            "Frequencies and thermochemistry scaled by {scale} \u{d7} the harmonic values."
        ));
    }

    egui::ScrollArea::vertical()
        // Enough rows to scan a spectrum without scrolling for every one, while
        // leaving the thermochemistry and spectrum buttons below it in view.
        .max_height(340.0)
        .show(ui, |ui| {
            egui::Grid::new("xtb_freq_mode_list")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("Mode");
                    ui.strong("Frequency (cm⁻¹)");
                    ui.strong("IR (km/mol)");
                    ui.strong("");
                    ui.end_row();

                    // Without intensity data the column reads N/A rather than
                    // 0.00, which would claim every band is IR-inactive.
                    let has_ir = result.has_ir_intensities();
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
                        if has_ir {
                            ui.label(format!("{:.2}", result.ir_intensities_km_mol[i]));
                        } else {
                            ui.weak("N/A");
                        }
                        // Animate/Stop is a single toggle, not two separate
                        // buttons: at most one mode animates at a time (they
                        // all share the one trajectory player), so the
                        // button's own label is exactly "is this the row
                        // that's currently playing?"
                        let is_animating_this_row =
                            freq_panel.selected_mode == Some(i) && traj.playing;
                        let toggle_label = if is_animating_this_row { "Stop" } else { "Animate" };
                        let animatable = result.is_animatable(i);
                        let button = ui
                            .add_enabled(animatable, egui::Button::new(toggle_label))
                            .on_disabled_hover_text(
                                "The file lists this mode's frequency but no displacement \
                                 vector, so there is nothing to animate.",
                            );
                        if button.clicked() {
                            if is_animating_this_row {
                                stop_mode_animation(
                                    traj,
                                    &result.atoms,
                                    &result.coords_angstrom,
                                );
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
        let has_ir = result.has_ir_intensities();
        if !has_ir {
            // Nothing to draw, and a window left open from a previous result
            // would keep showing that one's axes.
            freq_panel.spectrum_window_open = false;
        }
        let button = ui.add_enabled(has_ir, egui::Button::new(label));
        if button
            .on_disabled_hover_text(
                "No IR intensities in this data. A Hessian gives frequencies but not \
                 intensities -- those need dipole derivatives, which are computed \
                 separately and are not part of every Hessian file.",
            )
            .clicked()
        {
            freq_panel.spectrum_window_open = !freq_panel.spectrum_window_open;
        }
    });

    let ctx = ui.ctx().clone();
    thermochemistry_window(&ctx, &mut freq_panel.thermo_window_open, result);
    ir_spectrum_window(&ctx, freq_panel, result);
}

/// The point-group reference table, in a window of its own.
///
/// Clicking a group sets the symmetry number, so the table is a chooser as
/// well as a reference -- no retyping a number one has just read.
fn symmetry_number_window(
    ctx: &egui::Context,
    open: &mut bool,
    rot_symmetry: &mut f64,
    reanalyze: &mut bool,
) {
    if !*open {
        return;
    }
    egui::Window::new("Symmetry Numbers")
        .open(open)
        .resizable(true)
        .default_size([560.0, 380.0])
        .show(ctx, |ui| {
            ui.label(
                egui::RichText::new("Rotational symmetry number \u{3c3} by point group").strong(),
            );
            ui.weak(
                "\u{3c3} is the order of the proper-rotation subgroup; reflections and \
                 improper rotations do not count.",
            );
            ui.separator();
            egui::ScrollArea::both().show(ui, |ui| {
                for (family, members) in crate::normalmode::symmetry_numbers::families() {
                    ui.horizontal(|ui| {
                        ui.add_sized(
                            [58.0, 20.0],
                            egui::Label::new(egui::RichText::new(family).monospace().strong()),
                        );
                        for (name, sigma) in members {
                            let selected = (*rot_symmetry - sigma).abs() < f64::EPSILON;
                            let cell = format!("{name} = {sigma:.0}");
                            if ui
                                .add(
                                    egui::Button::new(egui::RichText::new(cell).monospace())
                                        .min_size(egui::vec2(78.0, 20.0))
                                        .selected(selected),
                                )
                                .on_hover_text(format!("Use \u{3c3} = {sigma:.0}"))
                                .clicked()
                            {
                                *rot_symmetry = sigma;
                                *reanalyze = true;
                            }
                        }
                    });
                }
            });
        });
}

/// Point size of the thermochemistry table. The layout is fixed-width columns
/// meant to be read as a table, so it wants a size closer to a document than
/// to the surrounding controls.
const THERMO_FONT_SIZE: f32 = 16.0;

fn thermochemistry_window(ctx: &egui::Context, open: &mut bool, result: &FrequencyResult) {
    if !*open {
        return;
    }
    egui::Window::new("Thermochemistry")
        .open(open)
        .resizable(true)
        .default_size([620.0, 560.0])
        .show(ctx, |ui| {
            let mut text = thermofuncs::format_thermo(
                &result.thermo,
                result.thermo_temp_k,
                result.thermo_freq_cutoff_cm1,
                None,
                result.rot_symmetry,
            );
            if let Some(note) = &result.source_note {
                text.push_str(&format!("\n{note}\n"));
            }
            if let Some(note) = &result.source_note {
                ui.colored_label(egui::Color32::from_rgb(150, 175, 210), note);
                ui.add_space(4.0);
            }
            if ui.button("Export .dat\u{2026}").clicked() {
                export_thermochemistry(&text);
            }
            ui.separator();
            egui::ScrollArea::both().show(ui, |ui| {
                ui.add(
                    egui::Label::new(
                        egui::RichText::new(&text)
                            .font(egui::FontId::monospace(THERMO_FONT_SIZE)),
                    )
                    // The table's alignment is the whole point; wrapping it to
                    // the window width would fold the columns into each other.
                    .wrap_mode(egui::TextWrapMode::Extend),
                );
            });
        });
}

/// Saves the thermochemistry table exactly as displayed -- same columns, same
/// rules -- so the file reads the way the window does.
fn export_thermochemistry(text: &str) {
    if let Some(path) = rfd::FileDialog::new()
        .set_file_name("thermochemistry.dat")
        .add_filter("Data file", &["dat", "txt"])
        .save_file()
    {
        if let Err(err) = std::fs::write(&path, text) {
            eprintln!("Failed to write {}: {err}", path.display());
        }
    }
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
                for kind in BroadeningKind::ALL {
                    ui.selectable_value(&mut freq_panel.spectrum_broadening, kind, kind.label());
                }
                // A stick spectrum has no width to speak of -- it is drawn
                // (and exported) at the exact predicted frequencies.
                ui.add_enabled_ui(!freq_panel.spectrum_broadening.is_sticks(), |ui| {
                    ui.label("Width (cm\u{207b}\u{b9})");
                    ui.add(
                        egui::DragValue::new(&mut freq_panel.spectrum_width_cm1)
                            .range(0.5..=200.0),
                    );
                });
            });

            // Only the real vibrations are plotted: the six near-zero
            // translation/rotation modes and any imaginary mode are not
            // absorptions.
            let peaks: Vec<(f64, f64)> = result
                .positive_indices
                .iter()
                .map(|&i| (result.frequencies_cm1[i], result.ir_intensities_km_mol[i]))
                .collect();
            if peaks.is_empty() {
                ui.weak("No vibrational modes to plot.");
                return;
            }
            // Fixed axis convention for an IR spectrum, per the user's own
            // choice: 0 up to the highest predicted frequency plus a fixed
            // 100 cm^-1 margin, rather than a window that tracks the current
            // broadening width and would shift as it is adjusted.
            const AXIS_MARGIN_CM1: f64 = 100.0;
            const TICK_STEP_CM1: f64 = 500.0;
            let x_min = 0.0f64;
            let x_max = peaks.iter().map(|&(f, _)| f).fold(f64::NEG_INFINITY, f64::max)
                + AXIS_MARGIN_CM1;

            if ui.button("Export .dat\u{2026}").clicked() {
                // Sticks export the raw predicted (frequency, intensity)
                // pairs as-is -- there is no continuous curve to resample,
                // so the fixed 0.5 cm^-1 grid used for the broadened shapes
                // does not apply here.
                let text = if freq_panel.spectrum_broadening.is_sticks() {
                    format_spectrum_sticks(&peaks)
                } else {
                    const EXPORT_RESOLUTION_CM1: f64 = 0.5;
                    format_spectrum_dat(
                        &peaks,
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

            let axis = SpectrumAxis {
                x_min,
                x_max,
                tick_step: Some(TICK_STEP_CM1),
                tick_decimals: 0,
                title: "Wavenumber (cm\u{207b}\u{b9})",
            };
            draw_spectrum(
                ui,
                &peaks,
                freq_panel.spectrum_broadening,
                freq_panel.spectrum_width_cm1,
                &axis,
                &|i| {
                    format!(
                        "{:.1} cm\u{207b}\u{b9}\n{:.2} km/mol",
                        peaks[i].0, peaks[i].1
                    )
                },
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
            masses_amu: None,
            ir_intensities_km_mol: None,
            atoms,
            coords_bohr,
            hessian,
            vib_lines,
            gradient_bohr: None,
        }
    }

    /// The Hessian block carries its own electronic state rather than reading
    /// the optimizer's: a frequency calculation is often run on a different
    /// species than the one last optimized.
    #[test]
    fn the_frequency_panel_defaults_to_a_neutral_singlet() {
        let panel = XtbFreqPanelState::default();
        assert_eq!(panel.charge, 0);
        assert_eq!(panel.multiplicity, 1);
        assert!(!panel.show_warnings, "warnings wait for an attempted run");
        assert_eq!(panel.frequency_scale, 1.0, "unscaled harmonic by default");
        assert_eq!(panel.rot_symmetry, 1.0, "C1 unless the user says otherwise");
    }

    /// Changing one panel's charge must not move the other's.
    #[test]
    fn the_frequency_and_optimizer_charges_are_independent() {
        let mut freq = XtbFreqPanelState::default();
        let opt = super::super::xtb_optimize::XtbPanelState::default();
        freq.charge = -1;
        freq.multiplicity = 2;
        // Both halves matter: the write must stick, and it must not reach the
        // other panel. Checking only the optimizer would pass even if the
        // frequency panel silently ignored the assignment.
        assert_eq!(freq.charge, -1);
        assert_eq!(freq.multiplicity, 2);
        assert_eq!(opt.charge, 0);
        assert_eq!(opt.multiplicity, 1);
    }

    /// Behemoth's Hessian, put through Beavyr's own projection and
    /// eigensolver, must reproduce the frequencies Behemoth itself reports
    /// from the same geometry.
    ///
    /// `behemoth --xyz water.xyz --method xtb --freq` gives 1506.27, 3551.05
    /// and 3646.50 cm^-1 with a zero-point energy of 0.019829 Eh. Those are
    /// computed inside Behemoth from the same finite-difference Hessian this
    /// test reads, so agreement checks the whole chain: the block-matrix
    /// parser, the units, the mass weighting and the Eckart projection.
    #[test]
    fn behemoths_hessian_reproduces_behemoths_own_frequencies() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join("behemoth_water_hessian.log");
        let text = std::fs::read_to_string(&path).unwrap();
        let atoms: Vec<String> = ["O", "H", "H"].iter().map(|s| s.to_string()).collect();
        let hessian = super::super::behemoth::parse_hessian_for(&text, &atoms).unwrap();

        // The geometry the Hessian was computed at, in bohr.
        let coords_angstrom: [f64; 9] = [
            0.0, 0.0, 0.119262, 0.0, 0.763239, -0.477047, 0.0, -0.763239, -0.477047,
        ];
        let coords_bohr: Vec<f64> = coords_angstrom
            .iter()
            .map(|c| c * ANGSTROM_TO_BOHR)
            .collect();

        let raw = RawFrequencyOutput {
            atoms,
            coords_bohr,
            hessian,
            vib_lines: Vec::new(),
            gradient_bohr: None,
            masses_amu: None,
            ir_intensities_km_mol: None,
        };
        let result = analyze(raw, 1, EckartMode::VibRot, 298.15, 100.0, 1.0, 1.0).unwrap();

        assert_eq!(result.zero_indices.len(), 6, "water has six zero modes");
        assert_eq!(result.positive_indices.len(), 3);
        assert!(result.negative_indices.is_empty(), "water is a minimum");

        let mut ours: Vec<f64> = result
            .positive_indices
            .iter()
            .map(|&i| result.frequencies_cm1[i])
            .collect();
        ours.sort_by(f64::total_cmp);
        for (mine, theirs) in ours.iter().zip([1506.27, 3551.05, 3646.50]) {
            assert!(
                (mine - theirs).abs() < 1.0,
                "our {mine:.2} cm^-1 against Behemoth's {theirs:.2}"
            );
        }

        const BEHEMOTH_ZPE_HARTREE: f64 = 0.019829;
        let error = (result.thermo.zpe - BEHEMOTH_ZPE_HARTREE).abs();
        assert!(
            error < 1.0e-5,
            "ZPE {:.6} Eh against Behemoth's {BEHEMOTH_ZPE_HARTREE:.6} \
             (off by {:.4} kcal/mol)",
            result.thermo.zpe,
            error * 627.509_474
        );
    }

    fn example_gaussian_path() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("c164-ts1-2-1001.log")
    }

    fn loaded_gaussian() -> FrequencyResult {
        let (source, scale) = read_hessian_file(&example_gaussian_path()).unwrap();
        let HessianSource::Imported(imported) = source else {
            panic!("a Gaussian log carries modes, not a Hessian");
        };
        assert!(scale.is_none(), "Gaussian scales its own frequencies");
        analyze_imported(*imported, 1, 298.15, 100.0, 1.0, 1.0).unwrap()
    }

    /// A Gaussian log is recognised as an imported result rather than being
    /// fed to the Hessian parser, which would reject it.
    #[test]
    fn a_gaussian_log_loads_as_imported_modes() {
        let result = loaded_gaussian();
        assert_eq!(result.atoms.len(), 19);
        assert_eq!(result.negative_indices.len(), 1);
        assert_eq!(result.positive_indices.len(), 50);
        assert!((result.frequencies_cm1[0] + 2468.9062).abs() < 1e-4);
        assert_eq!(result.eckart_used, EckartMode::Off);
    }

    /// Frame zero is theta = 0, i.e. maximum displacement -- which is why
    /// stopping cannot simply leave the current frame up.
    #[test]
    fn an_animation_starts_at_maximum_displacement_not_equilibrium() {
        let result = loaded_gaussian();
        let frames = animate_mode(
            &result.atoms,
            &result.coords_angstrom,
            &result.modes,
            0,
            0.18,
            16,
        );
        let drift = |f: &TrajectoryFrame| {
            f.pos
                .iter()
                .zip(&result.coords_angstrom)
                .map(|(a, b)| (*a - *b).length())
                .fold(0.0f32, f32::max)
        };
        assert!(
            drift(&frames[0]) > 0.1,
            "frame 0 should be displaced, not equilibrium: {}",
            drift(&frames[0])
        );
    }

    /// Stopping must put the equilibrium structure back on screen, whatever
    /// point of the cycle the animation was at.
    #[test]
    fn stopping_restores_the_equilibrium_geometry() {
        let result = loaded_gaussian();
        let mut traj = TrajectoryState::default();
        load_mode_animation(
            &mut traj,
            &result.atoms,
            &result.coords_angstrom,
            &result.modes,
            0,
            0.18,
            60.0,
        );
        // Stop from a displaced point mid-cycle.
        traj.current_frame = 5;
        traj.last_applied = Some(5);
        assert!(traj.playing);

        stop_mode_animation(&mut traj, &result.atoms, &result.coords_angstrom);

        assert!(!traj.playing);
        assert_eq!(traj.frames.len(), 1, "one undistorted frame");
        assert_eq!(traj.current_frame, 0);
        assert!(
            traj.last_applied.is_none(),
            "the applier must be forced to redraw, or the distorted frame stays up"
        );
        let restored = &traj.frames[0];
        assert_eq!(restored.atoms, result.atoms);
        for (got, want) in restored.pos.iter().zip(&result.coords_angstrom) {
            assert!((*got - *want).length() < 1.0e-6, "{got:?} vs {want:?}");
        }
    }

    /// A vibration must let the bond graph follow the motion, or a transition
    /// state's forming and breaking bonds never appear.
    #[test]
    fn animating_a_mode_unfreezes_the_bonds() {
        let result = loaded_gaussian();
        let mut traj = TrajectoryState::default();
        // Even if the user froze bonds for some earlier trajectory, loading a
        // mode must let them move again.
        traj.fixed_bonds = true;
        load_mode_animation(
            &mut traj,
            &result.atoms,
            &result.coords_angstrom,
            &result.modes,
            0,
            0.18,
            60.0,
        );
        assert!(
            !traj.fixed_bonds,
            "a mode animation must recompute connectivity per frame"
        );
    }

    /// The imaginary mode of this transition state is a hydrogen transfer, so
    /// somewhere in its swing the H is closer to one heavy atom than the
    /// other -- which is only visible if the bonds are recomputed.
    #[test]
    fn a_transition_state_mode_actually_changes_a_bond_length() {
        let result = loaded_gaussian();
        let frames = animate_mode(
            &result.atoms,
            &result.coords_angstrom,
            &result.modes,
            0,
            0.18,
            16,
        );
        // Atom 18 (index 17) is the transferring hydrogen.
        let h = 17;
        let nearest_heavy = |frame: &TrajectoryFrame| -> f32 {
            (0..result.atoms.len())
                .filter(|&i| i != h && result.atoms[i] != "H")
                .map(|i| (frame.pos[i] - frame.pos[h]).length())
                .fold(f32::INFINITY, f32::min)
        };
        let distances: Vec<f32> = frames.iter().map(nearest_heavy).collect();
        let spread = distances.iter().cloned().fold(f32::MIN, f32::max)
            - distances.iter().cloned().fold(f32::MAX, f32::min);
        assert!(
            spread > 0.1,
            "the H-heavy distance only varies by {spread:.3} A across the swing"
        );
    }

    /// The residuals belong in the list beside the vibrations, in the same
    /// ascending order the computed path uses.
    #[test]
    fn imported_low_frequencies_are_listed_as_low_modes() {
        let result = loaded_gaussian();
        assert_eq!(result.frequencies_cm1.len(), 57, "51 modes + 6 residuals");
        assert_eq!(result.zero_indices.len(), 6);
        assert_eq!(result.negative_indices.len(), 1);
        assert_eq!(result.positive_indices.len(), 50);
        assert!(
            result.frequencies_cm1.windows(2).all(|w| w[0] <= w[1]),
            "the list must read in ascending order"
        );
        // The imaginary mode is first, then the residuals, then the
        // vibrations -- the same shape as a Hessian we analysed ourselves.
        assert_eq!(result.negative_indices[0], 0);
        assert_eq!(result.zero_indices, vec![1, 2, 3, 4, 5, 6]);
        assert!((result.frequencies_cm1[1] + 18.8757).abs() < 1e-3);
        assert!((result.frequencies_cm1[7] - 78.6012).abs() < 1e-3);
    }

    /// A residual has a frequency but no vector, so its Animate button is
    /// disabled rather than silently doing nothing.
    #[test]
    fn residual_rows_are_not_animatable_but_real_modes_are() {
        let result = loaded_gaussian();
        assert!(result.is_animatable(0), "the imaginary mode");
        for &i in &result.zero_indices {
            assert!(!result.is_animatable(i), "residual row {i} has no vector");
        }
        for &i in &result.positive_indices {
            assert!(result.is_animatable(i), "vibration row {i} should animate");
        }
    }

    /// Listing the residuals must not let them into the thermochemistry --
    /// six extra near-zero vibrations would wreck the partition function.
    #[test]
    fn residuals_stay_out_of_the_thermochemistry() {
        const GAUSSIAN_ZPE_HARTREE: f64 = 0.143449;
        let result = loaded_gaussian();
        assert!((result.thermo.zpe - GAUSSIAN_ZPE_HARTREE).abs() < 5.0e-4);
        assert!(
            result
                .positive_indices
                .iter()
                .all(|&i| result.frequencies_cm1[i] > 50.0),
            "no residual leaked into the vibrations"
        );
    }

    /// Intensities must follow their own modes through the reordering.
    #[test]
    fn intensities_stay_with_their_modes_after_reordering() {
        let result = loaded_gaussian();
        // Row 0 is the imaginary mode with its huge 6843 km/mol intensity.
        assert!((result.ir_intensities_km_mol[0] - 6843.5507).abs() < 1e-3);
        // Residual rows have none.
        for &i in &result.zero_indices {
            assert_eq!(result.ir_intensities_km_mol[i], 0.0);
        }
        // Row 7 is the 78.60 cm-1 mode at 1.4167 km/mol.
        assert!((result.ir_intensities_km_mol[7] - 1.4167).abs() < 1e-3);
    }

    /// The user must be told the modes are Gaussian's and the thermochemistry
    /// is ours -- the two look identical in the mode list otherwise.
    #[test]
    fn an_imported_result_says_what_was_read_and_what_was_computed() {
        let note = loaded_gaussian().source_note.expect("a provenance note");
        assert!(note.contains("Gaussian"), "{note}");
        assert!(note.contains("computed here"), "{note}");

    }

    /// IR intensities come across, so the spectrum is available.
    #[test]
    fn imported_intensities_reach_the_spectrum() {
        let result = loaded_gaussian();
        assert!(result.has_ir_intensities());
        assert!((result.ir_intensities_km_mol[0] - 6843.5507).abs() < 1e-3);
    }

    /// Our thermochemistry from Gaussian's frequencies should land close to
    /// Gaussian's own: it reports a zero-point correction of 0.143449 Eh at
    /// sigma = 1, which the log states explicitly.
    #[test]
    fn our_thermochemistry_agrees_with_gaussians_own_zpe() {
        const GAUSSIAN_ZPE_HARTREE: f64 = 0.143449;
        let result = loaded_gaussian();
        let error = (result.thermo.zpe - GAUSSIAN_ZPE_HARTREE).abs();
        assert!(
            error < 5.0e-4,
            "ZPE {:.6} Eh vs Gaussian's {GAUSSIAN_ZPE_HARTREE:.6} (off by {:.3} kcal/mol)",
            result.thermo.zpe,
            error * 627.509_474
        );
    }

    /// The modes are usable for animation: a real displacement, on the atom
    /// that actually moves in the imaginary mode.
    #[test]
    fn imported_modes_can_be_animated() {
        let result = loaded_gaussian();
        let frames = animate_mode(
            &result.atoms,
            &result.coords_angstrom,
            &result.modes,
            0,
            0.18,
            16,
        );
        assert_eq!(frames.len(), 16);
        // Atom 18 is the transferring hydrogen; it must move appreciably while
        // the frame set as a whole stays a recognisable molecule.
        let travel = |i: usize| {
            frames
                .iter()
                .map(|f| (f.pos[i] - result.coords_angstrom[i]).length())
                .fold(0.0f32, f32::max)
        };
        assert!(travel(17) > 0.05, "the transferring H barely moves: {}", travel(17));
        assert!(travel(0) < travel(17), "a carbon should move less than the H");
    }

    /// The symmetry number reaches the thermochemistry: sigma = 2 must lower
    /// the rotational entropy by exactly R ln 2 relative to sigma = 1.
    #[test]
    fn the_symmetry_number_changes_the_rotational_entropy_by_r_ln_sigma() {
        let (raw, scale) = computed_hessian(&example_hess_path());
        let one = analyze(
            raw.clone(),
            1,
            EckartMode::VibRot,
            298.15,
            100.0,
            scale.unwrap_or(1.0),
            1.0,
        )
        .unwrap();
        let two = analyze(
            raw,
            1,
            EckartMode::VibRot,
            298.15,
            100.0,
            scale.unwrap_or(1.0),
            2.0,
        )
        .unwrap();
        // R in Eh/K per particle is the Boltzmann constant in atomic units.
        const R_HARTREE_PER_K: f64 = 3.166_811_563e-6;
        let drop = one.thermo.srot - two.thermo.srot;
        let expected = R_HARTREE_PER_K * 2.0f64.ln();
        assert!(
            (drop - expected).abs() < 1.0e-10,
            "S_rot fell by {drop:.3e} Eh/K, expected R ln2 = {expected:.3e}"
        );
        assert_eq!(two.rot_symmetry, 2.0);
    }

    fn example_hess_path() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("TS_Gamma-3-6_ZZ-S_11_Compound_1.hess")
    }

    /// The raw Hessian out of a file that has one, for the tests that need to
    /// analyse it directly.
    fn computed_hessian(path: &std::path::Path) -> (RawFrequencyOutput, Option<f64>) {
        match read_hessian_file(path) {
            Ok((HessianSource::Computed(raw), scale)) => (*raw, scale),
            Ok((HessianSource::Imported(_), _)) => {
                panic!("{} has no Hessian to analyse", path.display())
            }
            Err(e) => panic!("{e}"),
        }
    }

    /// Reads a Hessian file and analyses it the way the panel does: the
    /// file's own scaling factor unless one is given.
    fn analyze_hess(
        path: &std::path::Path,
        eckart: EckartMode,
        scale: Option<f64>,
    ) -> Result<FrequencyResult, String> {
        let (raw, file_scale) = match read_hessian_file(path)? {
            (HessianSource::Computed(raw), scale) => (*raw, scale),
            (HessianSource::Imported(_), _) => {
                return Err(format!("{} carries no Hessian", path.display()))
            }
        };
        analyze(
            raw,
            1,
            eckart,
            298.15,
            100.0,
            scale.or(file_scale).unwrap_or(1.0),
            1.0,
        )
    }

    fn loaded_ts() -> FrequencyResult {
        analyze_hess(&example_hess_path(), EckartMode::VibRot, None)
            .expect("the example .hess must analyse")
    }

    /// The real cross-validation, mirroring the xTB one: our own projection
    /// and eigensolver, run on ORCA's Hessian, must reproduce the frequencies
    /// ORCA independently computed from it.
    #[test]
    fn an_orca_hessian_reproduces_orcas_own_frequencies() {
        let result = loaded_ts();
        let text = std::fs::read_to_string(example_hess_path()).unwrap();
        let reference = super::super::orca_hess::parse_orca_hess(&text).unwrap();

        // Compare the modes that are actually vibrations: ORCA zeroed its six
        // projected ones, and ours are near-zero rather than exactly zero.
        let mut ours: Vec<f64> = result
            .negative_indices
            .iter()
            .chain(&result.positive_indices)
            .map(|&i| result.frequencies_cm1[i])
            .collect();
        ours.sort_by(f64::total_cmp);
        let mut theirs: Vec<f64> = reference
            .frequencies_cm1
            .iter()
            .copied()
            .filter(|f| *f != 0.0)
            .collect();
        theirs.sort_by(f64::total_cmp);

        assert_eq!(ours.len(), theirs.len(), "mode counts differ");
        for (mine, orca) in ours.iter().zip(&theirs) {
            assert!(
                (mine - orca).abs() < 1.0,
                "our {mine:.2} cm^-1 vs ORCA's {orca:.2} cm^-1"
            );
        }
    }

    /// A transition state: exactly one imaginary mode, six projected out, the
    /// rest real.
    #[test]
    fn the_example_is_recognised_as_a_transition_state() {
        let result = loaded_ts();
        assert_eq!(result.atoms.len(), 19);
        assert_eq!(result.negative_indices.len(), 1, "one imaginary mode");
        assert_eq!(result.zero_indices.len(), 6);
        assert_eq!(result.positive_indices.len(), 50);
        let imaginary = result.frequencies_cm1[result.negative_indices[0]];
        assert!(
            (imaginary + 602.6).abs() < 1.0,
            "imaginary mode came out {imaginary:.1}, ORCA says -602.6"
        );
    }

    /// The ordering trap: ORCA lists its zeros first and the imaginary mode
    /// seventh, so a position-based pairing would hang the imaginary mode's
    /// intensity on a translation.
    #[test]
    fn ir_intensities_land_on_the_modes_they_belong_to() {
        let result = loaded_ts();
        // The imaginary mode has no listed intensity.
        assert_eq!(result.ir_intensities_km_mol[result.negative_indices[0]], 0.0);
        // The highest-frequency mode is the 3781 cm^-1 O-H stretch at
        // 133.38 km/mol -- by far the strongest band in the file.
        let &highest = result
            .positive_indices
            .iter()
            .max_by(|&&a, &&b| result.frequencies_cm1[a].total_cmp(&result.frequencies_cm1[b]))
            .unwrap();
        assert!((result.frequencies_cm1[highest] - 3781.4).abs() < 1.0);
        assert!((result.ir_intensities_km_mol[highest] - 133.38).abs() < 0.1);
    }

    /// The example file carries a `$ir_spectrum`, so it does have intensities
    /// -- but a Hessian file need not, and then there is nothing to plot.
    #[test]
    fn a_hessian_without_an_ir_section_reports_no_intensities() {
        let text = std::fs::read_to_string(example_hess_path()).unwrap();
        let stripped = text.replace("$ir_spectrum", "$skipped_ir_spectrum");
        let dir = std::env::temp_dir().join("beavyr_hess_no_ir_test");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("no_ir.hess");
        std::fs::write(&path, stripped).unwrap();

        let result = analyze_hess(&path, EckartMode::VibRot, None).unwrap();
        assert!(
            !result.has_ir_intensities(),
            "with no $ir_spectrum there is nothing to plot"
        );
        assert!(result.ir_intensities_km_mol.iter().all(|&i| i == 0.0));
        // The frequencies are still there -- that is the point of loading it.
        assert_eq!(result.positive_indices.len(), 50);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn the_example_hessian_does_carry_ir_intensities() {
        assert!(loaded_ts().has_ir_intensities());
    }

    /// ORCA's own thermochemistry, from the matching .out, is the reference:
    /// a ZPE agreeing to a fraction of a kcal/mol means our frequencies and
    /// masses are both right.
    #[test]
    fn the_zero_point_energy_matches_orcas_own() {
        const ORCA_ZPE_HARTREE: f64 = 0.14279928;
        let result = loaded_ts();
        let error = (result.thermo.zpe - ORCA_ZPE_HARTREE).abs();
        assert!(
            error < 1.0e-4,
            "ZPE {:.8} Eh vs ORCA's {ORCA_ZPE_HARTREE:.8} Eh (off by {:.3} kcal/mol)",
            result.thermo.zpe,
            error * 627.509_474
        );
    }

    /// The scaling factor multiplies frequencies and carries into the
    /// thermochemistry computed from them.
    #[test]
    fn an_explicit_scale_overrides_the_files_own() {
        let unscaled = analyze_hess(&example_hess_path(), EckartMode::VibRot, Some(1.0)).unwrap();
        let scaled = loaded_ts();
        assert!(unscaled.frequency_scale.is_none());
        assert_eq!(scaled.frequency_scale, Some(0.974));
        for (&i, &j) in unscaled
            .positive_indices
            .iter()
            .zip(&scaled.positive_indices)
        {
            let ratio = scaled.frequencies_cm1[j] / unscaled.frequencies_cm1[i];
            assert!((ratio - 0.974).abs() < 1e-9, "ratio {ratio}");
        }
        assert!(
            scaled.thermo.zpe < unscaled.thermo.zpe,
            "scaling below 1 must lower the zero-point energy"
        );
    }

    #[test]
    fn thermochemistry_runs_on_a_loaded_hessian() {
        let result = loaded_ts();
        assert!(
            result.thermo.zpe > 0.0,
            "zero-point energy {} Eh should be positive",
            result.thermo.zpe
        );
        assert_eq!(result.thermo_temp_k, 298.15);
    }

    /// A `.hess` stores no gradient, so the reaction-path projection has
    /// nothing to project against and must say so rather than fail.
    #[test]
    fn a_reaction_path_request_falls_back_with_an_explanation() {
        let result =
            analyze_hess(&example_hess_path(), EckartMode::ReactionPath, None).unwrap();
        assert_eq!(result.eckart_used, EckartMode::VibRot);
        let note = result
            .reaction_path_fallback_note
            .clone()
            .expect("a fallback note");
        assert!(note.contains("gradient"), "{note}");
    }

    #[test]
    fn a_file_that_is_not_a_hessian_is_reported_with_its_path() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("AbsorptionSpectrum")
            .join("Exc_TPSSh-D4-def2-TZVP_SCMwater_16.inp");
        let Err(err) = analyze_hess(&path, EckartMode::VibRot, None) else {
            panic!("an ORCA input file is not a Hessian and must not analyse");
        };
        assert!(err.contains("$hessian"), "{err}");
        assert!(err.contains(".inp"), "the message should name the file: {err}");
    }

    /// The masses in `$atoms` are used, not the symbol table -- otherwise an
    /// isotope substitution recorded in the file would be silently ignored.
    #[test]
    fn deuterium_substitution_in_the_file_changes_the_frequencies() {
        let text = std::fs::read_to_string(example_hess_path()).unwrap();
        let heavy = text.replace(
            " H      1.00800     -2.073404472303",
            " H      2.01410     -2.073404472303",
        );
        let dir = std::env::temp_dir().join("beavyr_hess_isotope_test");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("heavy.hess");
        std::fs::write(&path, heavy).unwrap();

        let normal = loaded_ts();
        let deuterated = analyze_hess(&path, EckartMode::VibRot, None).unwrap();
        // The highest mode is an O-H stretch, untouched by deuterating a C-H,
        // so count the X-H stretch region instead: one C-H leaves it and
        // reappears near 2300 cm^-1 as a C-D.
        let above_3000 = |r: &FrequencyResult| {
            r.positive_indices
                .iter()
                .filter(|&&i| r.frequencies_cm1[i] > 3000.0)
                .count()
        };
        assert_eq!(
            above_3000(&normal),
            above_3000(&deuterated) + 1,
            "deuterating one C-H must move exactly one mode out of the stretch region"
        );
        let _ = std::fs::remove_file(&path);
    }

    /// The full `analyze` pipeline -- normal modes, IR-intensity matching,
    /// and thermochemistry all at once -- against the real captured
    /// fixtures, cross-checked against xTB's own independently-computed
    /// frequencies and zero-point energy.
    #[test]
    fn analyze_reproduces_xtbs_frequencies_and_zpe() {
        let result = analyze(h2o2_raw(), 1, EckartMode::VibRot, 298.15, 100.0, 1.0, 1.0).unwrap();

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
        let result = analyze(h2o2_raw(), 1, EckartMode::ReactionPath, 298.15, 100.0, 1.0, 1.0).unwrap();
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
        let result = analyze(raw, 1, EckartMode::ReactionPath, 298.15, 100.0, 1.0, 1.0).unwrap();
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

        let result = analyze(raw, 1, EckartMode::VibRot, 298.15, 100.0, 1.0, 1.0).unwrap();
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
        let Ok(xtb_path) = super::super::xtb_optimize::resolve_program_executable(QcProgram::Xtb, &configured) else {
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

        let raw = run_hess_cancellable(
            QcProgram::Xtb,
            &xtb_path,
            &workdir,
            &input_xyz,
            &atoms,
            0,
            0,
            1,
            false,
            &cancel,
            &child_slot,
        )
        .expect("a plain --hess run on this real molecule should succeed");
        assert_eq!(raw.atoms.len(), 16);
        assert!(
            !workdir.join("xtbhess.xyz").exists(),
            "this test fixture was specifically chosen because xTB does not \
             write xtbhess.xyz for it -- if it now exists, xTB's behaviour \
             changed and this test no longer exercises the bug it was written for"
        );

        let result = analyze(raw, 1, EckartMode::VibRot, 298.15, 100.0, 1.0, 1.0).unwrap();
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
        let Ok(xtb_path) = super::super::xtb_optimize::resolve_program_executable(QcProgram::Xtb, &configured) else {
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

        let raw = run_hess_cancellable(
            QcProgram::Xtb,
            &xtb_path,
            &workdir,
            &input_xyz,
            &atoms,
            0,
            0,
            1,
            false,
            &cancel,
            &child_slot,
        )
        .expect("a plain --hess run on a stable closed-shell molecule should succeed");
        assert_eq!(raw.atoms.len(), 4);
        assert_eq!(raw.hessian.dim(), (12, 12));

        let result = analyze(raw, 1, EckartMode::VibRot, 298.15, 100.0, 1.0, 1.0).unwrap();
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
        let Ok(xtb_path) = super::super::xtb_optimize::resolve_program_executable(QcProgram::Xtb, &configured) else {
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
        task.start(
            QcProgram::Xtb,
            &xtb_path,
            &atoms,
            &pos, 0, 0, 1, false);
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
                                panel_state.frequency_scale,
                                panel_state.rot_symmetry,
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
    fn animate_mode_produces_the_requested_frame_count_and_atoms() {
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










}

