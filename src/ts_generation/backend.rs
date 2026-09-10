//! Cancellable energy/gradient calls for RDA and constrained NEB images.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;

use anyhow::{bail, Context, Result};

use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
use crate::optimizer::traits::Objective;
use crate::qchem_interfaces::{behemoth, method, program::QcProgram, xtbrun};

#[derive(Default)]
pub struct RunControl {
    pub cancelled: AtomicBool,
    pub evaluations: AtomicUsize,
    pub latest_energy: Mutex<Option<f64>>,
    child: Mutex<Option<Child>>,
}

impl RunControl {
    pub fn cancel(&self) {
        self.cancelled.store(true, Ordering::Relaxed);
        if let Ok(mut slot) = self.child.lock() {
            if let Some(child) = slot.as_mut() {
                // Reaping happens on the worker; cancelling must not block the UI.
                let _ = child.kill();
            }
        }
    }

    pub fn check(&self) -> Result<()> {
        if self.cancelled.load(Ordering::Relaxed) {
            bail!("TS generation cancelled.");
        }
        Ok(())
    }
}

#[derive(Clone)]
pub struct BackendConfig {
    pub program: QcProgram,
    pub binary: PathBuf,
    pub method: method::MethodConfig,
    pub charge: i32,
    pub multiplicity: i32,
}

pub struct EnergyGradient {
    config: BackendConfig,
    atoms: Vec<String>,
    workdir: PathBuf,
    control: Arc<RunControl>,
    cached: Option<(Vec<f64>, f64, Vec<f64>)>,
}

impl EnergyGradient {
    pub fn new(
        config: BackendConfig,
        atoms: Vec<String>,
        workdir: PathBuf,
        control: Arc<RunControl>,
    ) -> Self {
        Self {
            config,
            atoms,
            workdir,
            control,
            cached: None,
        }
    }

    fn evaluate(&mut self, x: &[f64]) -> Result<(f64, Vec<f64>)> {
        self.control.check()?;
        if x.len() != self.atoms.len() * 3 || x.iter().any(|v| !v.is_finite()) {
            bail!("The optimizer requested an invalid geometry.");
        }
        if let Some((previous, energy, gradient)) = &self.cached {
            if previous == x {
                return Ok((*energy, gradient.clone()));
            }
        }
        fs::create_dir_all(&self.workdir)?;
        // Preserve f64 precision at the optimizer boundary; viewer coordinates are f32.
        fs::write(
            self.workdir.join("geometry.xyz"),
            xyz(&self.atoms, x, "TS generation energy/gradient"),
        )?;
        let gradient_path = self.workdir.join("gradient");
        if gradient_path.exists() {
            fs::remove_file(&gradient_path)?;
        }
        let mut command = match self.config.program {
            QcProgram::Xtb => {
                let mut command = Command::new(&self.config.binary);
                command
                    .current_dir(&self.workdir)
                    .arg("geometry.xyz")
                    .arg("--grad")
                    .arg("--chrg")
                    .arg(self.config.charge.to_string())
                    .arg("--uhf")
                    .arg((self.config.multiplicity - 1).to_string())
                    .args(method::xtb_method_args(&self.config.method))
                    .arg("-P")
                    .arg(self.config.method.nproc.max(1).to_string());
                command
            }
            QcProgram::Behemoth => behemoth::gradient_command(
                &self.config.binary,
                &self.workdir,
                "geometry.xyz",
                self.config.charge,
                self.config.multiplicity,
                &self.config.method,
            ),
            QcProgram::Dreiding => bail!("Choose xTB or Behemoth for the reaction energy surface."),
        };
        let stdout_path = self.workdir.join("gradient.stdout");
        let stderr_path = self.workdir.join("gradient.stderr");
        command
            .stdout(Stdio::from(fs::File::create(&stdout_path)?))
            .stderr(Stdio::from(fs::File::create(&stderr_path)?));
        self.control.check()?;
        let child = command
            .spawn()
            .with_context(|| format!("Could not start {}", self.config.program.label()))?;
        *self
            .control
            .child
            .lock()
            .map_err(|_| anyhow::anyhow!("Process state lock failed"))? = Some(child);
        self.control.evaluations.fetch_add(1, Ordering::Relaxed);
        let status = loop {
            let mut slot = self
                .control
                .child
                .lock()
                .map_err(|_| anyhow::anyhow!("Process state lock failed"))?;
            if self.control.cancelled.load(Ordering::Relaxed) {
                if let Some(mut child) = slot.take() {
                    let _ = child.kill();
                    let _ = child.wait();
                }
                bail!("TS generation cancelled.");
            }
            let process = slot
                .as_mut()
                .context("Energy calculation process disappeared")?;
            if let Some(status) = process.try_wait()? {
                slot.take();
                break status;
            }
            drop(slot);
            std::thread::sleep(Duration::from_millis(40));
        };
        self.control.check()?;
        if !status.success() {
            let stderr = fs::read_to_string(&stderr_path).unwrap_or_default();
            let detail = stderr
                .lines()
                .rev()
                .find(|l| !l.trim().is_empty())
                .unwrap_or("see gradient.stdout and gradient.stderr");
            // A reaction path must use one energy surface throughout. Do not fall
            // back to a different xTB method when a single-point calculation fails.
            bail!(
                "{} energy/gradient failed ({status}): {detail}",
                self.config.program.label()
            );
        }
        let stdout = fs::read_to_string(&stdout_path)?;
        let (energy, gradient) = match self.config.program {
            QcProgram::Xtb => (
                xtbrun::parse_xtb_energy(&stdout).map_err(anyhow::Error::msg)?,
                xtbrun::parse_xtb_grad(&gradient_path, self.atoms.len())
                    .map_err(anyhow::Error::msg)?,
            ),
            QcProgram::Behemoth => {
                behemoth::parse_gradient(&stdout, self.atoms.len()).map_err(anyhow::Error::msg)?
            }
            QcProgram::Dreiding => unreachable!(),
        };
        if !energy.is_finite()
            || gradient.len() != x.len()
            || gradient.iter().any(|v| !v.is_finite())
        {
            bail!("The backend returned an invalid energy or gradient.");
        }
        if let Ok(mut last) = self.control.latest_energy.lock() {
            *last = Some(energy);
        }
        self.cached = Some((x.to_vec(), energy, gradient.clone()));
        Ok((energy, gradient))
    }
}

impl Objective for EnergyGradient {
    fn energy(&mut self, x: &[f64]) -> Result<f64> {
        Ok(self.evaluate(x)?.0)
    }

    fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
        self.energy_gradient(x, grad)?;
        Ok(())
    }

    fn energy_gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<f64> {
        if grad.len() != x.len() {
            bail!("Gradient buffer has the wrong size.");
        }
        let (energy, gradient) = self.evaluate(x)?;
        grad.copy_from_slice(&gradient);
        Ok(energy)
    }
}

pub fn xyz(atoms: &[String], q_bohr: &[f64], comment: &str) -> String {
    let mut text = format!("{}\n{comment}\n", atoms.len());
    for (atom, q) in atoms.iter().zip(q_bohr.as_chunks::<3>().0) {
        text.push_str(&format!(
            "{atom:<3} {:16.10} {:16.10} {:16.10}\n",
            q[0] * BOHR_TO_ANGSTROM,
            q[1] * BOHR_TO_ANGSTROM,
            q[2] * BOHR_TO_ANGSTROM
        ));
    }
    text
}

pub fn save_text(path: &Path, text: &str) -> Result<()> {
    fs::write(path, text).with_context(|| format!("Could not save {}", path.display()))
}
