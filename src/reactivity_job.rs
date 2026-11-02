//! Three fixed-geometry single points per frame, retaining all wavefunctions and native outputs.
use super::data::{self, Results, Sources, Wave};
use crate::{
    trajectory::TrajectoryFrame,
    trajectory_wavefunctions::job::{self, Control, Engine},
};
use std::{
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
};
pub fn engine(
    qc: &crate::qchem_panel::QcPanelState,
    paths: &crate::qchem_interfaces::xtb_optimize::XtbPanelState,
    charge: i32,
    multiplicity: i32,
    cores: u32,
) -> Result<Engine, String> {
    use crate::qchem_interfaces::{
        job::QcMethod, program::QcProgram, xtb_optimize::resolve_program_executable,
    };
    let path = paths
        .paths
        .iter()
        .find(|(p, _)| *p == qc.program)
        .map(|(_, s)| s.as_str())
        .unwrap_or("");
    let binary = if qc.program.needs_binary_path() {
        resolve_program_executable(qc.program, path)?
    } else {
        PathBuf::new()
    };
    if qc.program == QcProgram::Xtb {
        if !matches!(qc.method, Some(QcMethod::Gfn1 | QcMethod::Gfn2)) {
            return Err("Choose GFN1 or GFN2; a force field has no electron density.".into());
        }
        return Ok(Engine::Xtb {
            binary,
            gfn: if qc.method == Some(QcMethod::Gfn1) {
                1
            } else {
                2
            },
            charge,
            multiplicity,
            cores: cores.clamp(1, 4) as usize,
        });
    }
    if !matches!(
        qc.program,
        QcProgram::Orca | QcProgram::PySCF | QcProgram::Psi4 | QcProgram::Psi4Py
    ) || !matches!(qc.method, Some(QcMethod::Hf | QcMethod::Dft))
    {
        return Err("Select HF or DFT with ORCA, PySCF or Psi4 in Quantum Chemistry, or select GFN1/GFN2-xTB. Correlated energies combined with an SCF Molden density are not a consistent Fukui calculation.".into());
    }
    let mut config = qc.config.clone();
    config.nproc = cores.clamp(1, 4);
    // Preserve the chosen backend/environment; only spin/charge differ between the three jobs.
    Ok(Engine::QuantumChemistry {
        binary,
        program: qc.program,
        method: qc.method,
        config,
        charge,
        multiplicity,
    })
}
fn changed_state(mut engine: Engine, charge: i32, multiplicity: i32) -> Engine {
    match &mut engine {
        Engine::Xtb {
            charge: q,
            multiplicity: m,
            ..
        }
        | Engine::QuantumChemistry {
            charge: q,
            multiplicity: m,
            ..
        } => {
            *q = charge;
            *m = multiplicity;
        }
    }
    engine
}
pub fn finite(
    engine: Engine,
    frames: Vec<(usize, TrajectoryFrame)>,
    charge: i32,
    multiplicities: [i32; 3],
    root: &Path,
    control: Control,
    progress: Arc<AtomicUsize>,
) -> Result<Results, String> {
    // Validate all electronic states before spending time on the first calculation.
    for (_, f) in &frames {
        let mol = crate::molecule::Molecule {
            atoms: f.atoms.clone(),
            pos: f.pos.clone(),
            bonds: vec![],
            hydrogen_bonds: vec![],
            topology: None,
        };
        for (q, m) in [
            (charge, multiplicities[0]),
            (charge - 1, multiplicities[1]),
            (charge + 1, multiplicities[2]),
        ] {
            crate::qchem_interfaces::valence::validate_electronic_state(&mol, q, m)
                .map_err(|e| e.to_string())?;
        }
    }
    let mut results = Results {
        version: 1,
        method: format!("Finite difference · {} · Mulliken", engine.label()),
        frames: vec![],
    };
    for (i, (original, geometry)) in frames.into_iter().enumerate() {
        let mut series = vec![];
        for (j, (label, q, m)) in [
            ("N", charge, multiplicities[0]),
            ("N-plus-1", charge - 1, multiplicities[1]),
            ("N-minus-1", charge + 1, multiplicities[2]),
        ]
        .into_iter()
        .enumerate()
        {
            if control.cancel.load(Ordering::Relaxed) {
                return Err(
                    "Reactivity calculation cancelled; completed native files remain on disk."
                        .into(),
                );
            }
            let dir = root.join(format!("frame-{:06}-{label}", i + 1));
            std::fs::create_dir(&dir).map_err(|e| e.to_string())?;
            series.push(job::calculate(
                changed_state(engine.clone(), q, m),
                vec![(original, geometry.clone())],
                dir,
                control.clone(),
            )?);
            progress.store(i * 3 + j + 1, Ordering::Relaxed);
        }
        let source =
            |i: usize| Wave::File(series[i].root.join(&series[i].manifest.frames[0].molden));
        let energy = |i: usize| {
            series[i].manifest.frames[0]
                .geometry
                .metadata
                .energy_hartree
                .ok_or("Missing state energy")
        };
        let mut geometry = geometry;
        geometry.metadata = series[0].manifest.frames[0].geometry.metadata.clone();
        let frame = data::finite(
            geometry,
            Sources::Finite {
                reference: source(0),
                added: source(1),
                removed: source(2),
            },
            [energy(0)?, energy(1)?, energy(2)?],
        )?;
        results.frames.push(frame);
        save(&results, &root.join("reactivity.beavyr-reactivity"))?;
    }
    if control.cancel.load(Ordering::Relaxed) {
        return Err("Reactivity calculation cancelled.".into());
    }
    Ok(results)
}
pub fn save(results: &Results, path: &Path) -> Result<(), String> {
    results.validate()?;
    let text = ron::ser::to_string(results).map_err(|e| e.to_string())?;
    crate::project::atomic_write(path, text.as_bytes()).map_err(|e| e.to_string())
}
pub fn read(path: &Path) -> Result<Results, String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let result: Results = ron::from_str(&text).map_err(|e| e.to_string())?;
    result.validate()?;
    Ok(result)
}
