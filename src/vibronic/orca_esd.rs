//! Native ORCA ESD route. Uses exactly the configured executable and an input file.
use super::{
    input::StateModes, spectral_job::Output, state::Settings, thermal_spectrum::HarmonicSpectrum,
    CM1_TO_AU, HARTREE_TO_EV,
};
use crate::qchem_interfaces::{
    job::{ExcitedStateOptions, JobType},
    orca_hess, orca_run,
};
use std::{path::Path, sync::atomic::AtomicBool};
pub fn input_text(s: &Settings, a: &StateModes) -> Result<String, String> {
    let c = &s.spectrum;
    if c.root == 0 || c.roots < c.root || c.multiplicity < 1 {
        return Err(
            "Choose a positive excited root, sufficient roots, and valid multiplicity".into(),
        );
    }
    super::thermal_spectrum::validate_grid(&c.grid)?;
    let xyz = orca_run::xyz_body(
        &a.atoms,
        &a.coords
            .iter()
            .map(|v| v * 0.529177210903)
            .collect::<Vec<_>>(),
    );
    let mut method = c.method.clone();
    method.nproc = method.nproc.clamp(1, 4);
    let extra = format!("{} ESD({})", method.orca_extra, c.process.orca());
    let base = orca_run::panel_input_text(
        JobType::SinglePoint,
        c.qc_method,
        c.charge,
        c.multiplicity,
        &method,
        ExcitedStateOptions::default(),
        &extra,
        &xyz,
    );
    let block=format!("%tddft\n NRoots {}\n IRoot {}\n TDA {}\nend\n%esd\n GSHESSIAN \"ground.hess\"\n ESHESSIAN \"excited.hess\"\n USEJ true\n SAMEFREQ false\n COORDSYS CART\n TCUTFREQ 0\n DOHT false\n DELE {:.12}\n TEMP {:.8}\n LINES {}\n LINEW {:.8}\n INLINEW {:.8}\n SPECRES {:.8}\n UNIT WN\nend\n",c.roots,c.root,c.tda,c.gap_ev/(HARTREE_TO_EV*CM1_TO_AU),c.grid.temperature,if c.grid.gaussian_fwhm_cm1==0.{"LORENTZ"}else if c.grid.lorentz_fwhm_cm1==0.{"GAUSS"}else{"VOIGT"},c.grid.lorentz_fwhm_cm1/2.,c.grid.gaussian_fwhm_cm1/(8.*2_f64.ln()).sqrt(),c.grid.spacing_cm1);
    let block = if c.process == super::optical_spectrum::Process::Raman {
        if c.grid.temperature != 0. {
            return Err("Native ORCA resonance Raman requires 0 K".into());
        }
        let lasers = c
            .raman
            .wavelengths()?
            .iter()
            .map(|x| format!("{:.10}", 1e7 / x))
            .collect::<Vec<_>>()
            .join(",");
        block.replace(
            " UNIT WN",
            &format!(
                " LASERE {lasers}\n RORDER {}\n RRSLINEW {:.8}\n RRINTENS false\n UNIT WN",
                c.raman.order, c.raman.line_fwhm
            ),
        )
    } else {
        block
    };
    let block = block.replace(
        " DOHT false",
        if c.properties.ht {
            " DOHT true\n CENTRALDIFF true\n FASTDER false"
        } else {
            " DOHT false"
        },
    );
    Ok(format!("{block}{base}"))
}
pub fn parse_spectrum(text: &str, temperature: f64) -> Result<HarmonicSpectrum, String> {
    let mut profile = vec![];
    for line in text.lines() {
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() < 2 {
            continue;
        }
        if let (Ok(x), Ok(y)) = (
            fields[0].replace(['D', 'd'], "E").parse::<f64>(),
            fields[1].replace(['D', 'd'], "E").parse::<f64>(),
        ) {
            if !x.is_finite() || !y.is_finite() {
                return Err("ORCA spectrum contains non-finite data".into());
            }
            profile.push((x * CM1_TO_AU * HARTREE_TO_EV, y));
        }
    }
    if profile.len() < 2 {
        return Err("ORCA did not produce a readable ESD spectrum".into());
    }
    profile.sort_by(|a, b| a.0.total_cmp(&b.0));
    let spacing = profile[1].0 - profile[0].0;
    Ok(HarmonicSpectrum {
        profile,
        correlation: vec![],
        zero_zero_ev: f64::NAN,
        spacing_ev: spacing,
        temperature,
        integral: f64::NAN,
        notes: vec![
            "Native ORCA ESD intensities; see the retained output for units and rates.".into(),
        ],
    })
}
pub fn run(
    s: &Settings,
    a: StateModes,
    b: StateModes,
    cancel: &AtomicBool,
) -> Result<Output, String> {
    if s.spectrum.orca_path.trim().is_empty() {
        return Err("Set the ORCA executable path before running native ESD".into());
    }
    let input = input_text(s, &a)?;
    let directory = super::process::directory()?;
    for (source, name) in [
        (&s.ground.path, "ground.hess"),
        (&s.excited.path, "excited.hess"),
    ] {
        let text = std::fs::read_to_string(source.trim()).map_err(|e| e.to_string())?;
        orca_hess::parse_orca_hess(&text)
            .map_err(|e| format!("Native ORCA ESD needs an ORCA-format Hessian: {e}"))?;
        std::fs::write(directory.join(name), text).map_err(|e| e.to_string())?;
    }
    std::fs::write(
        directory.join(format!("{}.inp", orca_run::JOB_BASENAME)),
        input,
    )
    .map_err(|e| e.to_string())?;
    let given = Path::new(s.spectrum.orca_path.trim());
    let binary = given.canonicalize().unwrap_or_else(|_| given.to_path_buf());
    let output = super::process::run(
        orca_run::run_command(&binary, &directory),
        &directory,
        orca_run::JOB_BASENAME,
        cancel,
    )?;
    orca_run::check_output(&output, true)
        .map_err(|e| format!("{e}\nRun directory: {}", directory.display()))?;
    if !output.contains("ORCA ESD FINISHED WITHOUT ERROR") {
        return Err(format!(
            "ORCA ended without successful ESD completion. Run directory: {}",
            directory.display()
        ));
    }
    let mut files = if s.spectrum.process == super::optical_spectrum::Process::Raman {
        std::fs::read_dir(&directory)
            .map_err(|e| e.to_string())?
            .filter_map(Result::ok)
            .map(|e| e.path())
            .filter(|p| {
                p.file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .starts_with(&format!("{}.spectrum.", orca_run::JOB_BASENAME))
            })
            .collect::<Vec<_>>()
    } else {
        vec![directory.join(format!("{}.spectrum", orca_run::JOB_BASENAME))]
    };
    files.sort();
    if s.spectrum.process == super::optical_spectrum::Process::Raman
        && files.len() != s.spectrum.raman.wavelengths()?.len()
    {
        return Err(format!("ORCA did not produce one Raman spectrum per requested laser. Its filenames round laser wavenumbers; use distinct, sufficiently separated lasers. Run directory: {}",directory.display()));
    }
    if files.is_empty() {
        return Err(format!(
            "ORCA produced no spectrum files. Run directory: {}",
            directory.display()
        ));
    }
    let mut spectra = vec![];
    for file in files {
        let text = std::fs::read_to_string(&file)
            .map_err(|e| format!("Cannot read {}: {e}", file.display()))?;
        let mut spectrum = parse_spectrum(&text, s.spectrum.grid.temperature)?;
        spectrum.zero_zero_ev = s.spectrum.gap_ev
            + 0.5
                * (b.selected.iter().map(|&i| b.frequencies[i]).sum::<f64>()
                    - a.selected.iter().map(|&i| a.frequencies[i]).sum::<f64>())
                * CM1_TO_AU
                * HARTREE_TO_EV;
        let suffix = if s.spectrum.process == super::optical_spectrum::Process::Raman {
            format!(
                " · laser {} cm⁻¹",
                file.extension().unwrap_or_default().to_string_lossy()
            )
        } else {
            String::new()
        };
        spectra.push((
            format!("ORCA {}{suffix}", s.spectrum.process.label()),
            spectrum,
        ));
    }
    Ok(Output {
        raman_lines: vec![],
        spectra,
        ground: a,
        excited: b,
        workdir: Some(directory.display().to_string()),
        notes: vec![],
    })
}
