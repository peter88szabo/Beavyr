//! Isolated property jobs launched with the user's saved Python environment.
use super::{
    input::StateModes,
    spectral_job::Controls,
    transition_properties::{self, Dipole, PropertySource},
};
use std::sync::atomic::AtomicBool;
pub fn collect(
    c: &Controls,
    a: &StateModes,
    cancel: &AtomicBool,
) -> Result<(Dipole, Option<String>), String> {
    let settings = &c.properties;
    match settings.source {
        PropertySource::Manual => {
            if settings.ht {
                return Err("Herzberg–Teller spectra need transition-dipole derivatives. Choose a property file, PySCF or Psi4.".into());
            }
            Ok((Dipole::constant(c.raman.dipole, a.selected.len()), None))
        }
        PropertySource::File => {
            let text = std::fs::read_to_string(settings.path.trim())
                .map_err(|e| format!("Cannot read transition properties: {e}"))?;
            Ok((transition_properties::parse(&text, a, settings.ht)?, None))
        }
        PropertySource::PyScf | PropertySource::Psi4 => {
            if c.root < 1
                || c.roots < c.root
                || c.multiplicity < 1
                || !settings.step.is_finite()
                || settings.step <= 0.
                || settings.step > 0.1
            {
                return Err("Invalid root, multiplicity or finite-difference step".into());
            }
            let module = if settings.source == PropertySource::PyScf {
                "pyscf"
            } else {
                "psi4"
            };
            let dir = super::process::directory()?;
            std::fs::write(
                dir.join("transition_driver.py"),
                include_str!("transition_driver.py"),
            )
            .map_err(|e| e.to_string())?;
            let geometry = a
                .atoms
                .iter()
                .enumerate()
                .map(|(i, s)| {
                    format!(
                        "{} {:.14e} {:.14e} {:.14e}\n",
                        s,
                        a.coords[3 * i],
                        a.coords[3 * i + 1],
                        a.coords[3 * i + 2]
                    )
                })
                .collect::<String>();
            std::fs::write(dir.join("geometry.bohr"), geometry).map_err(|e| e.to_string())?;
            let mut cmd = crate::qchem_interfaces::python_env::launch::settings(module)?
                .command([dir.join("transition_driver.py")]);
            cmd.current_dir(&dir);
            for (key, value) in [
                ("BACKEND", module.to_string()),
                ("ROOT", c.root.to_string()),
                ("ROOTS", c.roots.to_string()),
                ("CHARGE", c.charge.to_string()),
                ("MULTIPLICITY", c.multiplicity.to_string()),
                ("BASIS", c.method.basis.clone()),
                ("FUNCTIONAL", c.method.functional.clone()),
                ("TDA", c.tda.to_string()),
                (
                    "HF",
                    (c.qc_method == crate::qchem_interfaces::job::QcMethod::Hf).to_string(),
                ),
                ("HT", settings.ht.to_string()),
                ("STEP", settings.step.to_string()),
                ("CORES", c.method.nproc.clamp(1, 4).to_string()),
                ("MEMORY", c.method.memory_mb.to_string()),
            ] {
                let key = if key == "BACKEND" {
                    "BEAVYR_PROPERTY_BACKEND".to_string()
                } else {
                    format!("BEAVYR_{key}")
                };
                cmd.env(key, value);
            }
            cmd.env("OMP_NUM_THREADS", c.method.nproc.clamp(1, 4).to_string())
                .env("OPENBLAS_NUM_THREADS", "1")
                .env("MKL_NUM_THREADS", c.method.nproc.clamp(1, 4).to_string());
            let output=super::process::run(cmd,&dir,"properties",cancel).map_err(|e|format!("{module} transition properties failed. Check the configured Python environment and the calculation output.\n{e}"))?;
            if !output.contains("BEAVYR_TRANSITION_SUCCESS") {
                return Err(format!(
                    "Transition properties did not finish. Files: {}",
                    dir.display()
                ));
            }
            let text = std::fs::read_to_string(dir.join("transition.properties"))
                .map_err(|e| e.to_string())?;
            Ok((
                transition_properties::parse(&text, a, settings.ht)?,
                Some(dir.display().to_string()),
            ))
        }
    }
}
