//! Spectrum job settings and backend dispatch. ORCA uses native ESD.
use super::{
    harmonic_model,
    input::StateModes,
    state::Settings,
    thermal_spectrum::{HarmonicSpectrum, SpectrumGrid},
    HARTREE_TO_EV,
};
use crate::qchem_interfaces::{job::QcMethod, method::MethodConfig};
use std::sync::atomic::{AtomicBool, AtomicUsize};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Engine {
    Internal,
    Orca,
}
#[derive(Clone, Debug)]
pub struct Controls {
    pub process: super::optical_spectrum::Process,
    pub frequency_factor: bool,
    pub properties: super::transition_properties::PropertySettings,
    pub raman: super::raman_spectrum::RamanSettings,
    pub engine: Engine,
    pub gap_ev: f64,
    pub grid: SpectrumGrid,
    pub orca_path: String,
    pub method: MethodConfig,
    pub qc_method: QcMethod,
    pub charge: i32,
    pub multiplicity: i32,
    pub root: u32,
    pub roots: u32,
    pub tda: bool,
}
impl Default for Controls {
    fn default() -> Self {
        Self {
            process: super::optical_spectrum::Process::Absorption,
            frequency_factor: true,
            raman: Default::default(),
            properties: Default::default(),
            engine: Engine::Internal,
            gap_ev: 3.,
            grid: Default::default(),
            orca_path: crate::qchem_interfaces::config::load_path(
                crate::qchem_interfaces::program::QcProgram::Orca,
            ),
            method: Default::default(),
            qc_method: QcMethod::Dft,
            charge: 0,
            multiplicity: 1,
            root: 1,
            roots: 5,
            tda: true,
        }
    }
}
pub struct Output {
    pub raman_lines: Vec<super::raman_spectrum::RamanLine>,
    pub spectra: Vec<(String, HarmonicSpectrum)>,
    pub ground: StateModes,
    pub excited: StateModes,
    pub workdir: Option<String>,
    pub notes: Vec<String>,
}
pub fn run(
    settings: &Settings,
    cancel: &AtomicBool,
    progress: &AtomicUsize,
) -> Result<Output, String> {
    let control = &settings.spectrum;
    if !control.gap_ev.is_finite() || control.gap_ev <= 0. {
        return Err(
            "The electronic minimum-to-minimum gap must be positive, in eV (without ZPE)".into(),
        );
    }
    let (model, a, b) = harmonic_model::load(
        &settings.ground,
        &settings.excited,
        control.gap_ev / HARTREE_TO_EV,
    )?;
    if model.zero_zero() <= 0. {
        return Err(
            "The 0–0 excitation energy is not positive; check the two electronic energies".into(),
        );
    }
    match control.engine {
        Engine::Internal => {
            let (dipole, workdir) = if control.properties.ht
                || control.process == super::optical_spectrum::Process::Raman
            {
                let (d, p) = super::property_job::collect(control, &a, cancel)?;
                (Some(d), p)
            } else {
                (None, None)
            };
            let (spectra, raman_lines) =
                if control.process == super::optical_spectrum::Process::Raman {
                    super::raman_spectrum::calculate_with_dipole(
                        model,
                        &control.grid,
                        &control.raman,
                        dipole.as_ref().unwrap(),
                        cancel,
                        progress,
                    )?
                } else {
                    (
                        vec![(
                            control.process.label().into(),
                            super::optical_spectrum::calculate_with_dipole(
                                model,
                                control.process,
                                &control.grid,
                                control.frequency_factor,
                                dipole.as_ref(),
                                cancel,
                                progress,
                            )?,
                        )],
                        vec![],
                    )
                };
            Ok(Output {
                spectra,
                raman_lines,
                ground: a,
                excited: b,
                workdir,
                notes: vec![],
            })
        }
        Engine::Orca => super::orca_esd::run(settings, a, b, cancel),
    }
}
impl Output {
    pub fn csv(&self) -> String {
        use std::fmt::Write;
        let mut text = String::new();
        for note in &self.notes {
            let _ = writeln!(text, "# {note}");
        }
        if !self.raman_lines.is_empty() {
            text.push_str("laser_nm,modes,shift_cm-1,scattering_strength\n");
            for l in &self.raman_lines {
                let _ = writeln!(
                    text,
                    "{},\"{:?}\",{:.8},{:.12e}",
                    l.laser_nm, l.modes, l.shift_cm1, l.strength
                );
            }
        }
        for (label, spectrum) in &self.spectra {
            let _ = writeln!(
                text,
                "# {label}\n# Temperature_K: {}\n# 0-0_eV: {}\nenergy_eV,intensity",
                spectrum.temperature, spectrum.zero_zero_ev
            );
            for note in &spectrum.notes {
                let _ = writeln!(text, "# {}", note.replace(['\r', '\n'], " "));
            }
            for (x, y) in &spectrum.profile {
                let _ = writeln!(text, "{x:.10},{y:.12e}");
            }
        }
        text
    }
}
