//! A completed QC result is data, never a resumable task or process handle.
use crate::qchem_interfaces::{
    job::JobType,
    orca_run::{Energetics, ScfIteration},
    program::QcProgram,
};
use crate::qchem_panel::{
    run::{QcRunOutput, QcRunTask},
    scan_results::ScanReport,
};
use serde::{Deserialize, Serialize};
use std::{path::PathBuf, sync::Arc, time::Duration};

#[derive(Clone, Serialize, Deserialize)]
pub struct Calculation {
    scan: Option<ScanReport>,
    job: JobType,
    program: QcProgram,
    log: String,
    trajectory_text: Option<String>,
    hessian_path: Option<PathBuf>,
    molden_path: Option<PathBuf>,
    files: Vec<PathBuf>,
    kept_dir: Option<PathBuf>,
    energetics: Energetics,
    scf: Vec<ScfIteration>,
    converged: Option<bool>,
    duration: Option<Duration>,
}
impl Calculation {
    pub fn capture(task: &QcRunTask) -> Option<Self> {
        let r = task.last.as_ref()?;
        Some(Self {
            scan: r.scan.as_deref().cloned(),
            job: r.job,
            program: r.program,
            log: r.log.clone(),
            trajectory_text: r.trajectory_text.clone(),
            hessian_path: r.hessian_path.clone(),
            molden_path: r.molden_path.clone(),
            files: r.saveable.iter().map(|(_, p)| p.clone()).collect(),
            kept_dir: r.kept_dir.clone(),
            energetics: r.energetics.clone(),
            scf: r.scf.clone(),
            converged: r.converged,
            duration: task.last_duration,
        })
    }
    pub fn restore(self) -> QcRunTask {
        let mut task = QcRunTask::default();
        task.last = Some(QcRunOutput {
            scan: self.scan.map(Arc::new),
            job: self.job,
            program: self.program,
            log: self.log,
            trajectory_text: self.trajectory_text,
            hessian_path: self.hessian_path,
            molden_path: self.molden_path,
            saveable: self
                .files
                .into_iter()
                .map(|p| ("Recovered file", p))
                .collect(),
            kept_dir: self.kept_dir,
            energetics: self.energetics,
            scf: self.scf,
            converged: self.converged,
        });
        task.last_duration = self.duration;
        task.last_message = Some(
            "Restored completed calculation; native files remain at their original paths.".into(),
        );
        task
    }
}
