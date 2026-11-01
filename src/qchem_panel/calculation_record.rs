//! On-disk record of exactly the submitted calculation, independent of later UI edits.
use super::run::QcRunOutput;
use crate::qchem_interfaces::{
    job::{ExcitedStateOptions, JobType, QcMethod},
    method::MethodConfig,
    program::QcProgram,
};
use serde::Serialize;
use std::{
    fs,
    io::Write,
    path::Path,
    time::{Duration, SystemTime, UNIX_EPOCH},
};
pub const RECORD: &str = "beavyr_record";
#[derive(Serialize)]
pub struct Request<'a> {
    pub version: &'static str,
    pub program: QcProgram,
    pub method: Option<QcMethod>,
    pub job: JobType,
    pub executable: &'a Path,
    pub atoms: &'a [String],
    pub coordinates_angstrom: &'a [f64],
    pub charge: i32,
    pub multiplicity: i32,
    pub config: &'a MethodConfig,
    pub excited: ExcitedStateOptions,
    pub write_orbitals: bool,
    pub custom_input: Option<&'a str>,
}
fn write(dir: &Path, name: &str, text: impl AsRef<[u8]>) -> Result<(), String> {
    fs::write(dir.join(name), text).map_err(|e| {
        format!(
            "Cannot write calculation record {}: {e}",
            dir.join(name).display()
        )
    })
}
pub fn append(root: &Path, message: &str) {
    // Nested gradient evaluations share the outer calculation's progress log.
    if let Some(record) = root
        .ancestors()
        .map(|p| p.join(RECORD))
        .find(|p| p.is_dir())
    {
        if let Ok(mut file) = fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(record.join("run.log"))
        {
            let _ = writeln!(file, "{message}");
        }
    }
}
pub fn begin(root: &Path, r: &Request<'_>) -> Result<(), String> {
    let dir = root.join(RECORD);
    fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
    let settings = ron::ser::to_string_pretty(r, ron::ser::PrettyConfig::default())
        .map_err(|e| e.to_string())?;
    write(&dir, "settings.ron", settings)?;
    write(
        &dir,
        "started.txt",
        format!(
            "Unix timestamp: {}\nOS: {}\nArchitecture: {}\n",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default()
                .as_secs(),
            std::env::consts::OS,
            std::env::consts::ARCH
        ),
    )?;
    if let Some(module) = r.program.python_module() {
        let launch = crate::qchem_interfaces::python_env::launch::settings(module)?;
        write(
            &dir,
            "python_environment.ron",
            ron::ser::to_string_pretty(&launch, ron::ser::PrettyConfig::default())
                .map_err(|e| e.to_string())?,
        )?;
    }
    write(
        &dir,
        "status.txt",
        format!(
            "Running\nStarted Unix timestamp: {}\n",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default()
                .as_secs()
        ),
    )?;
    let (atoms, coords) =
        crate::calculation_geometry::real_coordinates(r.atoms, r.coordinates_angstrom)?;
    write(
        &dir,
        "initial.xyz",
        format!(
            "{}\nSubmitted geometry (angstrom)\n{}",
            atoms.len(),
            crate::qchem_interfaces::orca_run::xyz_body(&atoms, &coords)
        ),
    )?;
    if let Some(input) = r.custom_input {
        write(&dir, "custom_input.txt", input)?;
    }
    write(&dir,"README.txt","Calculation files copied by Beavyr. Native program files retain their original names and directory hierarchy.\nbeavyr_record/settings.ron records the submitted settings and original coordinates; initial.xyz contains real atoms only.\nA custom input, if present, was passed as written and may specify a different geometry/settings.\nstatus.txt distinguishes successful execution, nonconvergence, failure and cancellation.\nsummary.txt is Beavyr's returned output; native stdout/stderr files remain beside the input files.\ntrajectory.xyz and final.xyz are provided when the calculation produced geometry frames.\nNo restart compatibility beyond the saved native program files is implied.\n")?;
    append(
        root,
        &format!(
            "Starting {} — {}\nCharge: {} · multiplicity: {}\n",
            r.program.label(),
            r.job.label(),
            r.charge,
            r.multiplicity
        ),
    );
    Ok(())
}
pub fn finish(
    root: &Path,
    result: &mut Result<QcRunOutput, String>,
    duration: Duration,
) -> Result<(), String> {
    let dir = root.join(RECORD);
    fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
    match result {
        Ok(output) => {
            output.kept_dir = Some(root.to_owned());
            write(&dir, "summary.txt", &output.log)?;
            if let Some(trajectory) = &output.trajectory_text {
                write(&dir, "trajectory.xyz", trajectory)?;
                if let Some(last) = super::run::last_frame_xyz(trajectory) {
                    write(&dir, "final.xyz", last)?;
                }
            }
            if let Some(scan) = &output.scan {
                write(&dir, "scan.csv", scan.csv())?;
                let xyz = scan.xyz();
                write(&dir, "trajectory.xyz", &xyz)?;
                if let Some(last) = super::run::last_frame_xyz(&xyz) {
                    write(&dir, "final.xyz", last)?;
                }
            }
            let mut files = String::from("Artifact\tPath\n");
            let mut artifacts = output
                .saveable
                .iter()
                .map(|(name, path)| ((*name).to_string(), path.clone()))
                .collect::<Vec<_>>();
            if let Some(path) = &output.hessian_path {
                artifacts.push(("Hessian".into(), path.clone()));
            }
            if let Some(path) = &output.molden_path {
                artifacts.push(("Orbitals".into(), path.clone()));
            }
            for (i, (label, path)) in artifacts.iter().enumerate() {
                if !path.is_file() {
                    return Err(format!(
                        "Calculation artifact is missing: {}",
                        path.display()
                    ));
                }
                let saved = if let Ok(relative) = path.strip_prefix(root) {
                    relative.to_owned()
                } else if path.is_file() {
                    let target = dir.join(format!(
                        "artifact_{i:03}_{}",
                        path.file_name().unwrap_or_default().to_string_lossy()
                    ));
                    fs::copy(path, &target)
                        .map_err(|e| format!("Cannot retain {}: {e}", path.display()))?;
                    target.strip_prefix(root).unwrap().to_owned()
                } else {
                    return Err(format!(
                        "Calculation artifact is missing: {}",
                        path.display()
                    ));
                };
                files.push_str(&format!(
                    "{}\t{}\n",
                    label.replace(['\n', '\t'], " "),
                    saved.display()
                ));
            }
            write(&dir, "artifacts.tsv", files)?;
            let status = if output.scan.as_ref().is_some_and(|s| !s.complete) {
                "Incomplete scan"
            } else if output.converged == Some(false) {
                "Finished without convergence"
            } else {
                "Completed"
            };
            let message = format!(
                "{status}\nElapsed seconds: {:.3}\nConverged: {:?}\nFinal energy (Hartree): {:?}\n",
                duration.as_secs_f64(),
                output.converged,
                output.energetics.final_energy_hartree
            );
            write(&dir, "status.txt", &message)?;
            append(root, &message);
        }
        Err(error) => {
            let status = if error.starts_with("Cancelled.") {
                "Cancelled"
            } else {
                "Failed"
            };
            write(
                &dir,
                "status.txt",
                format!(
                    "{status}\nElapsed seconds: {:.3}\n{error}\n",
                    duration.as_secs_f64()
                ),
            )?;
            write(&dir, "error.txt", error.as_bytes())?;
            append(root, &format!("{status}: {error}"));
        }
    }
    Ok(())
}
pub fn mark_cancelled(root: &Path) {
    let _ = write(
        &root.join(RECORD),
        "status.txt",
        "Cancelled\nResult was cancelled before application to the viewer.\n",
    );
    append(root, "Cancelled before application to the viewer.");
}
