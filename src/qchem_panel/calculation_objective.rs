//! Cancellable energy and gradient evaluations shared by scans and MECP.
use super::*;
use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
use crate::optimizer::traits::Objective;
pub(super) fn configure_method(config: &mut MethodConfig, method: Option<QcMethod>) {
    use crate::qchem_interfaces::method::{BehemothMethod as B, XtbMethod as X};
    match method {
        Some(QcMethod::Hf) => {
            config.behemoth = B::Hf;
            config.functional = "hf".into();
        }
        Some(QcMethod::Dft) => config.behemoth = B::Dft,
        Some(QcMethod::Mp2) => config.behemoth = B::Mp2,
        Some(QcMethod::RiMp2) => config.behemoth = B::RiMp2,
        Some(QcMethod::Tasi) => config.behemoth = B::Tasi,
        Some(QcMethod::Gfn1) => {
            config.behemoth = B::Gfn1Xtb;
            config.xtb = X::Gfn1;
        }
        Some(QcMethod::Gfn2) => config.xtb = X::Gfn2,
        Some(QcMethod::GfnFf) => config.xtb = X::GfnFf,
        _ => {}
    }
}
pub(super) struct Evaluator<'a> {
    pub(super) program: QcProgram,
    pub(super) method: Option<QcMethod>,
    pub(super) binary: &'a Path,
    pub(super) dir: &'a Path,
    pub(super) atoms: &'a [String],
    pub(super) charge: i32,
    pub(super) multiplicity: i32,
    pub(super) config: &'a MethodConfig,
    pub(super) cancel: &'a AtomicBool,
    pub(super) child: &'a Mutex<Option<Child>>,
    pub(super) evaluations: usize,
    pub(super) cached: Option<(Vec<f64>, f64, Option<Vec<f64>>)>,
}
impl Evaluator<'_> {
    fn evaluate(&mut self, x: &[f64], gradient: bool) -> anyhow::Result<(f64, Option<Vec<f64>>)> {
        if self.cancel.load(Ordering::SeqCst) {
            anyhow::bail!("Cancelled.");
        }
        if x.len() != self.atoms.len() * 3 || x.iter().any(|v| !v.is_finite()) {
            anyhow::bail!("Invalid calculation geometry.");
        }
        if let Some((old, e, g)) = &self.cached {
            if old == x && (!gradient || g.is_some()) {
                return Ok((*e, g.clone()));
            }
        }
        self.evaluations += 1;
        let dir = self.dir.join(format!("evaluation_{:06}", self.evaluations));
        let coords: Vec<f64> = x.iter().map(|v| v * BOHR_TO_ANGSTROM).collect();
        let job = if gradient {
            JobType::Gradient
        } else {
            JobType::SinglePoint
        };
        let generated = if self.program == QcProgram::PySCF {
            Some(pyscf_script(
                job,
                self.method.unwrap(),
                self.charge,
                self.multiplicity,
                self.config,
                self.atoms,
                &coords,
            ))
        } else {
            None
        };
        crate::qchem_panel::calculation_record::append(self.dir, &format!(
            "{}: evaluation {} ({})", dir.display(), self.evaluations, job.label()));
        let output = run_cancellable(
            self.program,
            self.method,
            job,
            self.binary,
            &dir,
            self.atoms,
            &coords,
            self.charge,
            self.multiplicity,
            self.config,
            Default::default(),
            false,
            generated.as_deref(),
            self.cancel,
            self.child,
        )
        .map_err(anyhow::Error::msg)?;
        if output.converged == Some(false) || output.log.contains("SCF NOT CONVERGED") {
            anyhow::bail!("{} did not converge at this geometry", self.program.label());
        }
        let energy = output
            .energetics
            .final_energy_hartree
            .filter(|e| e.is_finite())
            .ok_or_else(|| anyhow::anyhow!("{} returned no finite energy", self.program.label()))?;
        let grad = if gradient {
            let path = if matches!(self.program, QcProgram::Behemoth | QcProgram::Xtb) {
                dir.join("energy")
            } else {
                dir.clone()
            };
            Some(
                crate::qchem_interfaces::gradient::read(
                    self.program,
                    &path,
                    self.atoms.len(),
                    &output.log,
                )
                .map_err(anyhow::Error::msg)?
                .1,
            )
        } else {
            None
        };
        crate::qchem_panel::calculation_record::append(self.dir, &format!(
            "Evaluation {}: energy = {:.14} Hartree{}", self.evaluations, energy,
            grad.as_ref().map(|g| format!(", gradient norm = {:.8} Hartree/bohr", g.iter().map(|v| v * v).sum::<f64>().sqrt())).unwrap_or_default()));
        self.cached = Some((x.to_vec(), energy, grad.clone()));
        Ok((energy, grad))
    }
}
impl Objective for Evaluator<'_> {
    fn energy(&mut self, x: &[f64]) -> anyhow::Result<f64> {
        Ok(self.evaluate(x, false)?.0)
    }
    fn gradient(&mut self, x: &[f64], g: &mut [f64]) -> anyhow::Result<()> {
        self.energy_gradient(x, g)?;
        Ok(())
    }
    fn energy_gradient(&mut self, x: &[f64], g: &mut [f64]) -> anyhow::Result<f64> {
        let (e, grad) = self.evaluate(x, true)?;
        g.copy_from_slice(&grad.unwrap());
        Ok(e)
    }
}
pub(super) fn pyscf_script(
    job: JobType,
    method: QcMethod,
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    atoms: &[String],
    coords: &[f64],
) -> String {
    let mut config = config.clone();
    if matches!(method, QcMethod::Hf | QcMethod::Mp2) {
        config.functional = "hf".into();
    }
    let script = crate::qchem_interfaces::pyscf_run::panel_script_text(
        job,
        method,
        charge,
        multiplicity,
        &config,
        Default::default(),
        false,
        atoms,
        coords,
    );
    if method != QcMethod::Mp2 {
        return script;
    }
    let marker = "write_xyz(\"final.xyz\", mol, \"PySCF geometry\")";
    script.replace(marker,&format!("{marker}\nfrom pyscf import mp\ncorrelation = mp.MP2(mf)\ncorrelation.kernel()\nenergy = float(correlation.e_tot)\nif not numpy.isfinite(energy):\n    raise RuntimeError('Non-finite MP2 energy')\n"))
        .replace("g = mf.nuc_grad_method().kernel()","g = correlation.nuc_grad_method().kernel()")
}
