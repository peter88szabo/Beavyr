//! Analysis jobs and reproducible, copyable exports, independent of the GUI.
use super::{
    input,
    state::{AnalysisKind, Settings},
    *,
};
use std::fmt::Write;
pub enum AnalysisData {
    Spectrum(super::spectral_job::Output),
    FranckCondon {
        activity: FranckCondonAnalysis,
        band: Option<VibronicBand>,
    },
    Duschinsky(Duschinsky),
}
pub struct AnalysisResult {
    pub settings: Settings,
    pub data: AnalysisData,
    pub notices: Vec<String>,
}
#[cfg(test)]
pub fn run(settings: Settings) -> Result<AnalysisResult, String> {
    run_with_control(
        settings,
        &std::sync::atomic::AtomicBool::new(false),
        &std::sync::atomic::AtomicUsize::new(0),
    )
}
pub fn run_with_control(
    settings: Settings,
    cancel: &std::sync::atomic::AtomicBool,
    progress: &std::sync::atomic::AtomicUsize,
) -> Result<AnalysisResult, String> {
    if settings.kind == AnalysisKind::Spectrum {
        let output = super::spectral_job::run(&settings, cancel, progress)?;
        return Ok(AnalysisResult {
            settings,
            data: AnalysisData::Spectrum(output),
            notices: vec![],
        });
    }
    if !settings.cutoff.is_finite() || settings.cutoff < 0. {
        return Err("Frequency cutoff must be finite and nonnegative".into());
    }
    let a = input::load(&settings.ground, settings.cutoff)?;
    let mut notices = vec![];
    if a.imaginary > 0 {
        notices.push(format!("Ground state has {} imaginary modes. These are excluded; the harmonic minimum approximation may be inappropriate.",a.imaginary));
    }
    let data = match settings.kind {
        AnalysisKind::Spectrum => unreachable!(),
        AnalysisKind::FranckCondon => {
            let gradient = input::gradient(&settings.gradient, &a)?;
            let activity =
                franck_condon_activity(&a.modes, &a.frequencies, &gradient, &a.masses, &a.selected)
                    .map_err(|e| e.to_string())?;
            if activity
                .modes
                .iter()
                .any(|m| !m.huang_rhys.is_finite() || !m.resonance_raman.is_finite())
            {
                return Err(
                    "Mode activity overflow: check the gradient and very low frequencies".into(),
                );
            }
            if activity.unaccounted_fraction > 0.01 {
                notices.push(format!("The selected vibrations do not account for {:.2}% of the gradient's mass-weighted squared norm. Check the state, excluded modes and geometry.",100.*activity.unaccounted_fraction));
            }
            let band = if settings.make_band {
                match absorption_band(&activity, &settings) {
                    Ok(b) => Some(b),
                    Err(e) => {
                        notices.push(format!("Franck–Condon mode activities are available, but the band could not be generated: {e}"));
                        None
                    }
                }
            } else {
                None
            };
            AnalysisData::FranckCondon { activity, band }
        }
        AnalysisKind::Duschinsky => {
            let b = input::load(&settings.excited, settings.cutoff)?;
            input::same_atoms(&a.atoms, &b.atoms)?;
            if a.masses
                .iter()
                .zip(&b.masses)
                .any(|(a, b)| (a - b).abs() > 1e-5 * a.abs().max(1.))
            {
                return Err("The two Hessians use different isotope masses. Use identical masses for both states.".into());
            }
            if b.imaginary > 0 {
                notices.push(format!(
                    "Excited state has {} imaginary modes; these are excluded.",
                    b.imaginary
                ));
            }
            AnalysisData::Duschinsky(
                duschinsky(
                    &a.modes,
                    &a.frequencies,
                    &a.coords,
                    &b.modes,
                    &b.frequencies,
                    &b.coords,
                    &a.masses,
                    &a.selected,
                    &b.selected,
                )
                .map_err(|e| e.to_string())?,
            )
        }
    };
    Ok(AnalysisResult {
        settings,
        data,
        notices,
    })
}
fn absorption_band(
    activity: &FranckCondonAnalysis,
    settings: &Settings,
) -> Result<VibronicBand, String> {
    if settings.quanta == 0 || settings.quanta > 1000 {
        return Err("Use 1–1000 quanta per mode".into());
    }
    // The core normalises truncated progressions. Reject an appreciable lost
    // tail instead of making that truncation invisible to the user.
    let mut retained = 1.;
    for m in &activity.modes {
        let mut w = (-m.huang_rhys).exp();
        let mut sum = w;
        for n in 1..=settings.quanta {
            w *= m.huang_rhys / n as f64;
            sum += w;
        }
        retained *= sum.min(1.);
    }
    if retained < 0.999 {
        return Err(format!("The quanta limit retains only {:.4}% of the joint population. Increase maximum quanta per mode (at least 99.9% is required).",100.*retained));
    }
    if settings.energy_ev / HARTREE_TO_EV <= activity.reorganisation_cm1 * CM1_TO_AU {
        return Err("Vertical energy must exceed the reorganisation energy for a positive absorption origin".into());
    }
    vibronic_band(
        activity,
        settings.energy_ev / HARTREE_TO_EV,
        BandKind::Absorption,
        settings.resolution,
        settings.quanta,
    )
    .map_err(|e| e.to_string())
}
impl AnalysisResult {
    pub fn csv(&self) -> String {
        let mut out = String::new();
        // Metadata are comments; line breaks are escaped so filenames cannot
        // alter the exported table structure.
        let clean = |s: &str| s.replace('\r', "\\r").replace('\n', "\\n");
        let _ = writeln!(
            out,
            "# Ground Hessian: {}\n# Ground XYZ: {}\n# Frequency cutoff (cm^-1): {}",
            clean(&self.settings.ground.path),
            clean(&self.settings.ground.xyz),
            self.settings.cutoff
        );
        for note in &self.notices {
            let _ = writeln!(out, "# Warning: {note}");
        }
        match &self.data {
            AnalysisData::Spectrum(output) => {
                let c = &self.settings.spectrum;
                let _ = writeln!(out, "# Excited Hessian: {}\n# Excited XYZ: {}\n# Process: {}\n# Engine: {:?}\n# Electronic minimum gap_eV: {:.10}\n# Broadening FWHM_cm-1: Lorentzian={}, Gaussian={}\n# Grid spacing_cm-1: {}\n# Herzberg-Teller: {}\n# Property source: {:?}\n# Property file: {}\n# Central-difference step_bohr: {}\n# Charge: {}\n# Multiplicity: {}\n# Root: {} of {}\n# TDA: {}\n# Method: {:?}\n# Functional: {}\n# Basis: {}", clean(&self.settings.excited.path), clean(&self.settings.excited.xyz), c.process.label(),c.engine,c.gap_ev,c.grid.lorentz_fwhm_cm1,c.grid.gaussian_fwhm_cm1,c.grid.spacing_cm1,c.properties.ht,c.properties.source,clean(&c.properties.path),c.properties.step,c.charge,c.multiplicity,c.root,c.roots,c.tda,c.qc_method,clean(&c.method.functional),clean(&c.method.basis));
                if c.process == super::optical_spectrum::Process::Raman {
                    let _ = writeln!(out, "# Energy axis is Raman shift (eV), not incident photon energy.\n# Laser wavelengths_nm: {}\n# Raman order: {}\n# Raman linewidth_FWHM_cm-1: {}\n# Scattered-photon frequency fourth-power factor: {}",clean(&c.raman.lasers),c.raman.order,c.raman.line_fwhm,c.raman.photon_factor);
                } else {
                    let _ = writeln!(out, "# Optical photon-energy prefactor: {}\n# Relative spectrum; no absolute lifetime inferred",c.frequency_factor);
                }
                out.push_str(&output.csv());
            }
            AnalysisData::FranckCondon { activity, band } => {
                let _=writeln!(out,"# Excited gradient: {}\n# Model: displaced harmonic oscillators; Raman short-time approximation, relative intensities\n# Reorganisation (cm^-1): {:.10}\n# Unaccounted gradient fraction: {:.10}",clean(&self.settings.gradient),activity.reorganisation_cm1,activity.unaccounted_fraction);
                out.push_str("mode,frequency_cm^-1,projected_gradient_au,dimensionless_displacement,huang_rhys,reorganisation_cm^-1,relative_resonance_raman\n");
                let relative = activity.normalised_resonance_raman();
                for (m, (_, rr)) in activity.modes.iter().zip(relative) {
                    let _ = writeln!(
                        out,
                        "{},{:.8},{:.12e},{:.10},{:.10},{:.10},{:.10}",
                        m.index + 1,
                        m.frequency_cm1,
                        m.gradient_au,
                        m.displacement,
                        m.huang_rhys,
                        m.reorganisation_cm1,
                        rr
                    );
                }
                if let Some(band) = band {
                    let _=writeln!(out,"\n# Zero-temperature absorption; no Duschinsky mixing or Herzberg-Teller term\n# Vertical energy (eV): {}\n# Grid spacing (cm^-1): {}\n# Maximum quanta per mode: {}\n# 0-0 energy (eV): {:.10}\nenergy_eV,bin_probability",self.settings.energy_ev,self.settings.resolution,self.settings.quanta,band.zero_zero_energy*HARTREE_TO_EV);
                    for (e, i) in band.profile_ev() {
                        let _ = writeln!(out, "{e:.10},{i:.12e}");
                    }
                }
            }
            AnalysisData::Duschinsky(d) => {
                let _=writeln!(out,"# Excited Hessian: {}\n# Excited XYZ: {}\n# J = L_ground^T M R L_excited; K dimensionless along ground modes\nground_mode,frequency_cm^-1,K,huang_rhys_from_K",clean(&self.settings.excited.path),clean(&self.settings.excited.xyz));
                for (i, &mode) in d.rows.iter().enumerate() {
                    let k = d.displacement[i];
                    let _ = writeln!(
                        out,
                        "{},{:.8},{k:.10},{:.10}",
                        mode + 1,
                        d.row_frequencies_cm1[i],
                        0.5 * k * k
                    );
                }
                out.push_str("\nexcited_mode,frequency_cm^-1,best_ground_mode,ground_frequency_cm^-1,overlap_fraction,participation,shift_cm^-1\n");
                for c in d.correspondence() {
                    let _ = writeln!(
                        out,
                        "{},{:.8},{},{:.8},{:.10},{:.10},{:.8}",
                        c.column_index + 1,
                        c.column_frequency_cm1,
                        c.row_index + 1,
                        c.row_frequency_cm1,
                        c.overlap,
                        c.participation,
                        c.shift_cm1
                    );
                }
                out.push_str("\nground_mode");
                for mode in &d.columns {
                    let _ = write!(out, ",excited_mode_{}", mode + 1);
                }
                out.push('\n');
                for (i, mode) in d.rows.iter().enumerate() {
                    let _ = write!(out, "{}", mode + 1);
                    for j in 0..d.columns.len() {
                        let _ = write!(out, ",{:.12e}", d.rotation[(i, j)]);
                    }
                    out.push('\n');
                }
            }
        }
        out
    }
}
