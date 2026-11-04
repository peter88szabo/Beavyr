//! Preserve native backend files, but offer one self-contained Hessian export.
use super::QcRunOutput;
use crate::qchem_interfaces::{
    hessian_io::{self, HessianData},
    program::QcProgram,
};
use std::path::Path;

pub(crate) fn load(output: &QcRunOutput, path: &Path) -> Result<HessianData, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("Cannot read {}: {e}", path.display()))?;
    if hessian_io::is_self_contained(&text) {
        return hessian_io::parse(&text, None);
    }
    let dir = path
        .parent()
        .ok_or("Hessian has no calculation directory")?;
    let read = |p: &Path| {
        std::fs::read_to_string(p).map_err(|e| format!("Cannot read {}: {e}", p.display()))
    };
    let geometry = match output.program {
        // The xTB Hessian may be in a reoriented frame: never use the viewer
        // or optimization trajectory as a substitute for this geometry.
        QcProgram::Xtb => crate::qchem_interfaces::hessian_file::parse_g98_geometry(&read(
            &dir.join("hess/g98.out"),
        )?)?,
        QcProgram::Behemoth => hessian_io::geometry(&read(&dir.join("hess/input.xyz"))?)?,
        QcProgram::Dreiding => hessian_io::geometry(&read(&dir.join("hessian_geometry.xyz"))?)?,
        QcProgram::PySCF | QcProgram::Psi4 | QcProgram::Psi4Py => {
            hessian_io::geometry(&read(&dir.join("final.xyz"))?)?
        }
        QcProgram::SparrowPy => hessian_io::geometry(&read(&dir.join("geometry.xyz"))?)?,
        _ => {
            return Err(
                "Hessian has no matching geometry. Load its equilibrium XYZ alongside it.".into(),
            )
        }
    };
    hessian_io::parse(&text, Some(&geometry))
}

pub(crate) fn prepare(output: &mut QcRunOutput) -> Result<(), String> {
    let entries: Vec<_> = output
        .saveable
        .iter()
        .filter(|(label, _)| label.contains("Hessian"))
        .map(|(_, path)| path.clone())
        .collect();
    for path in entries {
        let data = load(output, &path)?;
        let already_standard = std::fs::read_to_string(&path)
            .map(|text| hessian_io::is_standard(&text))
            .unwrap_or(false);
        let standard = if already_standard {
            path.clone()
        } else {
            path.with_file_name(format!(
                "{}_standard.hess",
                path.file_stem().unwrap_or_default().to_string_lossy()
            ))
        };
        if !already_standard {
            hessian_io::save(&standard, &data)?;
        }
        for (label, file) in &mut output.saveable {
            if *file == path {
                *label = "Hessian (Smite-compatible)";
                *file = standard.clone();
            }
        }
        if output.hessian_path.as_ref() == Some(&path) {
            output.hessian_path = Some(standard);
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "hessian_output_tests.rs"]
mod tests;
