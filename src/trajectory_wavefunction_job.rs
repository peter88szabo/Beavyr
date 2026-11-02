//! One wavefunction per geometry, with portable relative paths and incremental manifests.
use crate::{
    orbitals::molden::{self, OrbitalData},
    qchem_interfaces::{
        job::{JobType, QcMethod},
        method::MethodConfig,
        program::QcProgram,
    },
    trajectory::TrajectoryFrame,
};
use bevy::prelude::*;
use serde::{Deserialize, Serialize};
use std::{
    fs,
    path::{Component, Path, PathBuf},
    process::{Child, Command, Stdio},
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering},
        Arc, Mutex,
    },
    time::Duration,
};
#[path="wavefunction_resume.rs"]
pub mod resume;
pub const MANIFEST: &str = "wavefunctions.beavyr-wfn";
#[derive(Clone, Serialize, Deserialize)]
pub enum Engine {
    Xtb {
        binary: PathBuf,
        gfn: u8,
        charge: i32,
        multiplicity: i32,
        cores: usize,
    },
    QuantumChemistry {
        binary: PathBuf,
        program: QcProgram,
        method: Option<QcMethod>,
        config: MethodConfig,
        charge: i32,
        multiplicity: i32,
    },
}
impl Engine {
    fn state(&self) -> (i32, i32) {
        match self {
            Self::Xtb {
                charge,
                multiplicity,
                ..
            }
            | Self::QuantumChemistry {
                charge,
                multiplicity,
                ..
            } => (*charge, *multiplicity),
        }
    }
    pub fn label(&self) -> String {
        match self {
            Self::Xtb { gfn, .. } => format!("GFN{gfn}-xTB"),
            Self::QuantumChemistry {
                program,
                method,
                config,
                ..
            } => format!(
                "{} / {} / {} / {}",
                program.label(),
                method.map(|m| m.label()).unwrap_or("unspecified"),
                config.functional,
                config.basis
            ),
        }
    }
}
#[derive(Clone, Serialize, Deserialize)]
pub struct Frame {
    pub source_index: usize,
    pub geometry: TrajectoryFrame,
    pub molden: PathBuf,
    pub spin_resolved: bool,
    #[serde(default)]
    pub checksum: Option<u64>,
}
#[derive(Clone, Serialize, Deserialize)]
pub struct Manifest {
    pub version: u32,
    pub method: String,
    pub charge: i32,
    pub multiplicity: i32,
    pub frames: Vec<Frame>,
    pub requested_frames: usize,
    pub complete: bool,
}
#[derive(Clone)]
pub struct Series {
    pub root: PathBuf,
    pub manifest: Manifest,
}
impl Series {
    pub fn frames(&self) -> Vec<TrajectoryFrame> {
        self.manifest
            .frames
            .iter()
            .map(|f| f.geometry.clone())
            .collect()
    }
    pub fn matches(&self, traj: &crate::trajectory::TrajectoryState) -> bool {
        self.manifest.frames.len() == traj.frames.len()
            && self
                .manifest
                .frames
                .get(traj.current_frame)
                .zip(traj.frames.get(traj.current_frame))
                .is_some_and(|(a, b)| geometry_matches(&a.geometry, b))
    }
}
pub fn geometry_matches(a: &TrajectoryFrame, b: &TrajectoryFrame) -> bool {
    a.atoms == b.atoms
        && a.pos.len() == b.pos.len()
        && a.pos.iter().zip(&b.pos).all(|(x, y)| x.distance(*y) < 1e-5)
}
#[derive(Clone, Default)]
pub struct Control {
    pub cancel: Arc<AtomicBool>,
    pub progress: Arc<AtomicUsize>,
    pub active: Arc<AtomicUsize>,
    pub reused: Arc<AtomicUsize>,
    pub phase: Arc<Mutex<String>>,
    pub child: Arc<Mutex<Option<Child>>>,
}
impl Control {
    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::SeqCst);
        if let Ok(mut child) = self.child.lock() {
            if let Some(c) = child.as_mut() {
                let _ = c.kill();
            }
        }
    }
    fn check(&self) -> Result<(), String> {
        if self.cancel.load(Ordering::SeqCst) {
            Err("Wavefunction calculation cancelled. Completed frames remain on disk.".into())
        } else {
            Ok(())
        }
    }
}
pub fn prepared(
    frames: &[TrajectoryFrame],
    first: usize,
    last: usize,
    stride: usize,
) -> Result<Vec<(usize, TrajectoryFrame)>, String> {
    if first == 0 || last < first || last > frames.len() || stride == 0 {
        return Err("Select a valid frame range.".into());
    }
    let mut selected = vec![];
    for i in (first - 1..last).step_by(stride) {
        let f = &frames[i];
        let xyz: Vec<_> = f
            .pos
            .iter()
            .flat_map(|p| p.to_array().map(f64::from))
            .collect();
        let (atoms, coordinates) = crate::calculation_geometry::real_coordinates(&f.atoms, &xyz)?;
        let pos = coordinates
            .chunks_exact(3)
            .map(|p| Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32))
            .collect();
        selected.push((
            i,
            TrajectoryFrame {
                atoms,
                pos,
                metadata: f.metadata.clone(),
            },
        ));
    }
    if selected.iter().any(|(_, f)| f.atoms != selected[0].1.atoms) {
        return Err(
            "Wavefunction trajectories require the same atom identities and order in every frame."
                .into(),
        );
    }
    Ok(selected)
}
fn run_xtb(
    binary: &Path,
    gfn: u8,
    charge: i32,
    multiplicity: i32,
    cores: usize,
    dir: &Path,
    control: &Control,
) -> Result<(PathBuf, f64), String> {
    let cores = cores.clamp(1, 4).to_string();
    let stdout = fs::File::create(dir.join("program.stdout")).map_err(|e| e.to_string())?;
    let stderr = fs::File::create(dir.join("program.stderr")).map_err(|e| e.to_string())?;
    let mut cmd = Command::new(binary);
    cmd.current_dir(dir)
        .arg("input.xyz")
        .args([
            "--gfn",
            &gfn.to_string(),
            "--chrg",
            &charge.to_string(),
            "--uhf",
            &(multiplicity - 1).to_string(),
            "--molden",
            "--parallel",
            &cores,
        ])
        .env("OMP_NUM_THREADS", &cores)
        .env("OPENBLAS_NUM_THREADS", &cores)
        .env("MKL_NUM_THREADS", &cores)
        .stdout(Stdio::from(stdout))
        .stderr(Stdio::from(stderr));
    control.check()?;
    let child = cmd
        .spawn()
        .map_err(|e| format!("Cannot start xTB at {}: {e}", binary.display()))?;
    *control.child.lock().map_err(|e| e.to_string())? = Some(child);
    let status = loop {
        if control.cancel.load(Ordering::SeqCst) {
            control.cancel();
        }
        let result = control
            .child
            .lock()
            .map_err(|e| e.to_string())?
            .as_mut()
            .ok_or("Missing xTB process")?
            .try_wait()
            .map_err(|e| e.to_string())?;
        if let Some(s) = result {
            break s;
        }
        std::thread::sleep(Duration::from_millis(30));
    };
    control.child.lock().map_err(|e| e.to_string())?.take();
    control.check()?;
    let log = fs::read_to_string(dir.join("program.stdout")).map_err(|e| e.to_string())?;
    if !status.success() {
        return Err(format!(
            "xTB failed ({status}) in {}:\n{}\n{}",
            dir.display(),
            log,
            fs::read_to_string(dir.join("program.stderr")).unwrap_or_default()
        ));
    }
    if log.to_ascii_uppercase().contains("SCC NOT CONVERGED") {
        return Err(format!("xTB SCC did not converge in {}.", dir.display()));
    }
    let energy = crate::qchem_interfaces::xtbrun::parse_xtb_energy(&log)?;
    if !energy.is_finite() {
        return Err("Nonfinite xTB energy.".into());
    }
    Ok((dir.join("molden.input"), energy))
}
pub fn calculate(
    engine: Engine,
    frames: Vec<(usize, TrajectoryFrame)>,
    root: PathBuf,
    control: Control,
) -> Result<Series, String> {
    let _lock=resume::lock(&root)?;
    resume::save(&root,&engine,&frames)?;
    calculate_inner(engine,frames,root,control,vec![])
}
fn calculate_inner(engine: Engine, frames: Vec<(usize,TrajectoryFrame)>, root: PathBuf, control: Control, completed: Vec<Frame>) -> Result<Series,String> {
    let (charge, multiplicity) = engine.state();
    if frames.is_empty() {
        return Err("No frames to calculate.".into());
    }
    let mut manifest = Manifest {
        version: 1,
        method: engine.label(),
        charge,
        multiplicity,
        frames: completed,
        requested_frames: frames.len(),
        complete: false,
    };
    fs::write(
        root.join("settings.ron"),
        ron::ser::to_string(&engine).map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;
    let skip=manifest.frames.len();
    control.progress.store(skip,Ordering::Relaxed);
    for (n, (source_index, mut geometry)) in frames.into_iter().enumerate().skip(skip) {
        control.active.store(n+1,Ordering::Relaxed);
        *control.phase.lock().map_err(|e|e.to_string())?=format!("Computing source frame {}",source_index+1);
        control.check()?;
        let mol = crate::molecule::Molecule {
            atoms: geometry.atoms.clone(),
            pos: geometry.pos.clone(),
            bonds: vec![],
            hydrogen_bonds: vec![],
            topology: None,
        };
        crate::qchem_interfaces::valence::validate_electronic_state(&mol, charge, multiplicity)
            .map_err(|e| e.to_string())?;
        let mut attempt=0;
        let dir=loop {
            let name=if attempt==0 {format!("frame-{:06}",n+1)} else {format!("frame-{:06}-retry-{attempt}",n+1)};
            let path=root.join(name);
            match fs::create_dir(&path) {Ok(())=>break path,Err(e) if e.kind()==std::io::ErrorKind::AlreadyExists=>attempt+=1,Err(e)=>return Err(e.to_string())}
        };
        fs::write(
            dir.join("input.xyz"),
            crate::trajectory_data::xyz(&[geometry.clone()], [0]),
        )
        .map_err(|e| e.to_string())?;
        let (path, energy) = match &engine {
            Engine::Xtb {
                binary, gfn, cores, ..
            } => {
                let (p, e) = run_xtb(binary, *gfn, charge, multiplicity, *cores, &dir, &control)?;
                (p, Some(e))
            }
            Engine::QuantumChemistry {
                binary,
                program,
                method,
                config,
                ..
            } => {
                if !matches!(method, Some(QcMethod::Hf | QcMethod::Dft)) {
                    return Err("Choose HF or DFT for a wavefunction trajectory; correlated densities are not exported by this workflow.".into());
                }
                if !matches!(
                    program,
                    QcProgram::Orca | QcProgram::PySCF | QcProgram::Psi4 | QcProgram::Psi4Py
                ) {
                    return Err("Choose ORCA, PySCF or Psi4 in Quantum Chemistry for this route, or use GFN2-xTB.".into());
                }
                let coords: Vec<_> = geometry
                    .pos
                    .iter()
                    .flat_map(|p| p.to_array().map(f64::from))
                    .collect();
                let mut config = config.clone();
                config.nproc = config.nproc.clamp(1, 4);
                let output = crate::qchem_panel::run::run_cancellable(
                    *program,
                    *method,
                    JobType::SinglePoint,
                    binary,
                    &dir,
                    &geometry.atoms,
                    &coords,
                    charge,
                    multiplicity,
                    &config,
                    Default::default(),
                    true,
                    None,
                    &control.cancel,
                    &control.child,
                )?;
                if output.converged == Some(false) {
                    return Err(format!(
                        "Electronic calculation did not converge in {}",
                        dir.display()
                    ));
                }
                (output.molden_path.ok_or("This calculation did not produce a Molden wavefunction; choose a supported HF/DFT method and enable its Molden export dependencies.")?,output.energetics.final_energy_hartree)
            }
        };
        control.check()?;
        let mut text = fs::read_to_string(&path)
            .map_err(|e| format!("No readable Molden file {}: {e}", path.display()))?;
        if matches!(engine, Engine::Xtb { .. }) {
            text = text.replacen(
                "[Title]",
                "[Title]\nBeavyr source: xTB (primitive-normalized coefficients)",
                1,
            );
        }
        let data = molden::parse_molden(&text).map_err(|e| e.to_string())?;
        check_geometry(&data, &geometry)?;
        let target = dir.join("wavefunction.molden");
        let checksum=Some(resume::checksum(text.as_bytes()));
        fs::write(&target, text).map_err(|e| e.to_string())?;
        if let Some(e) = energy {
            geometry.metadata = crate::trajectory_data::Metadata::energy(
                e,
                false,
                format!(
                    "Source_frame={} E={e:.14} Eh | {}",
                    source_index + 1,
                    geometry.metadata.comment
                ),
            );
        }
        manifest.frames.push(Frame {
            source_index,
            geometry,
            molden: target
                .strip_prefix(&root)
                .map_err(|e| e.to_string())?
                .to_owned(),
            spin_resolved: data.is_open_shell(),
            checksum,
        });
        manifest.complete = manifest.frames.len() == manifest.requested_frames;
        crate::project::atomic_write(
            &root.join(MANIFEST),
            ron::ser::to_string(&manifest)
                .map_err(|e| e.to_string())?
                .as_bytes(),
        )
        .map_err(|e| e.to_string())?;
        control.progress.store(n + 1, Ordering::Relaxed);
    }
    control.check()?;
    manifest.complete=true;
    control.active.store(0,Ordering::Relaxed);
    Ok(Series { root, manifest })
}
pub fn check_geometry(data: &OrbitalData, frame: &TrajectoryFrame) -> Result<(), String> {
    if data.atoms != frame.atoms
        || data.positions.len() != frame.pos.len()
        || data
            .positions
            .iter()
            .zip(&frame.pos)
            .any(|(a, b)| !a.is_finite() || a.distance(*b) > 1e-3)
    {
        return Err("Molden geometry does not match its trajectory frame. Surfaces cannot be assigned to a different geometry.".into());
    }
    Ok(())
}
pub fn read(path: &Path) -> Result<Series, String> {
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let manifest: Manifest = ron::from_str(&text).map_err(|e| e.to_string())?;
    if manifest.version != 1 || manifest.frames.is_empty() {
        return Err("Unsupported or empty wavefunction trajectory.".into());
    }
    for f in &manifest.frames {
        if f.molden
            .components()
            .any(|c| !matches!(c, Component::Normal(_)))
            || f.molden.as_os_str().is_empty()
        {
            return Err("Invalid wavefunction path in manifest.".into());
        }
        if f.geometry.atoms.is_empty()
            || f.geometry.atoms.len() != f.geometry.pos.len()
            || f.geometry.pos.iter().any(|p| !p.is_finite())
            || f.geometry.atoms != manifest.frames[0].geometry.atoms
        {
            return Err("Invalid trajectory geometry in manifest.".into());
        }
    }
    Ok(Series {
        root: path.parent().unwrap_or(Path::new(".")).to_owned(),
        manifest,
    })
}
pub fn load(series: &Series, index: usize) -> Result<(OrbitalData, String), String> {
    let frame = series
        .manifest
        .frames
        .get(index)
        .ok_or("No wavefunction for this frame.")?;
    let text = fs::read_to_string(series.root.join(&frame.molden)).map_err(|e| e.to_string())?;
    if frame.checksum.is_some_and(|sum|sum!=resume::checksum(text.as_bytes())) {return Err(format!("The wavefunction file for frame {} was modified after calculation.",index+1));}
    let data = molden::parse_molden(&text).map_err(|e| e.to_string())?;
    check_geometry(&data, &frame.geometry)?;
    Ok((data, text))
}

#[cfg(test)]
#[path = "trajectory_wavefunction_tests.rs"]
mod tests;
