//! Standalone MRCC MINP adapter. Capability restrictions follow the installed
//! manual sections 6.2/6.4 and native MTEST inputs (including NEVPT2).
use super::{
    job::{ExcitedStateOptions, JobType, QcMethod},
    method::MethodConfig,
};
use std::{path::Path, process::Command};
#[path = "mrcc_output.rs"]
mod output;
pub use output::{parse_energy, parse_gradient, parse_spectrum};

pub const METHODS: &[QcMethod] = &[
    QcMethod::Hf,
    QcMethod::Dft,
    QcMethod::Mp2,
    QcMethod::RiMp2,
    QcMethod::Mp2F12,
    QcMethod::Ccsd,
    QcMethod::CcsdT,
    QcMethod::UccsdT,
    QcMethod::CcsdF12,
    QcMethod::CcsdTF12,
    QcMethod::Ccsdt,
    QcMethod::CcsdtQ,
    QcMethod::Ccsdtq,
    QcMethod::Casscf,
    QcMethod::Nevpt2,
    QcMethod::Drpa,
    QcMethod::Adc2,
    QcMethod::ScsAdc2,
    QcMethod::SosAdc2,
];
pub const BASES: &[&str] = &[
    "STO-3G",
    "def2-SVP",
    "def2-TZVP",
    "def2-TZVPP",
    "def2-QZVPP",
    "cc-pVDZ",
    "cc-pVTZ",
    "cc-pVQZ",
    "aug-cc-pVDZ",
    "aug-cc-pVTZ",
    "cc-pCVDZ",
    "cc-pCVTZ",
    "cc-pCVQZ",
    "cc-pwCVDZ",
    "cc-pwCVTZ",
    "cc-pwCVQZ",
    "aug-cc-pCVTZ",
    "aug-cc-pwCVTZ",
    "cc-pVDZ-F12",
    "cc-pVTZ-F12",
    "cc-pVQZ-F12",
];
pub const FUNCTIONALS: &[&str] = &[
    "PBE", "PBE0", "BLYP", "B3LYP", "BP86", "TPSS", "TPSSh", "M06", "M06-2X", "B2PLYP",
];
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct Options {
    pub correlate_core: bool,
    pub unrestricted: bool,
    pub method: QcMethod,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            correlate_core: false,
            unrestricted: false,
            method: QcMethod::Hf,
        }
    }
}
pub fn adc(m: QcMethod) -> bool {
    matches!(m, QcMethod::Adc2 | QcMethod::ScsAdc2 | QcMethod::SosAdc2)
}
pub fn analytical_gradient(m: QcMethod) -> bool {
    matches!(
        m,
        QcMethod::Hf
            | QcMethod::Dft
            | QcMethod::RiMp2
            | QcMethod::Ccsd
            | QcMethod::Ccsdt
            | QcMethod::Ccsdtq
    )
}
pub fn job_problem(method: Option<QcMethod>, job: JobType) -> Option<&'static str> {
    if !matches!(
        job,
        JobType::SinglePoint | JobType::TdDft | JobType::Gradient | JobType::Optimize
    ) {
        return Some("MRCC currently supports energies, excited states and optimization with analytical gradients in Beavyr. Choose ORCA or another supported engine for this job.");
    }
    let Some(m) = method else { return None };
    if !METHODS.contains(&m) {
        return Some("Choose a method supported by MRCC.");
    }
    match job {
        JobType::TdDft if !adc(m) && !matches!(m,QcMethod::Hf|QcMethod::Dft) => Some("Choose HF, DFT or an ADC(2) variant for MRCC excited states."),
        JobType::SinglePoint if adc(m) => Some("ADC calculates excited states; choose UV-Vis / excited states."),
        JobType::Gradient|JobType::Optimize if !analytical_gradient(m) => Some("This MRCC method has no supported analytical gradient. Choose HF, DFT, RI-MP2, CCSD, CCSDT or CCSDTQ, or run an Energy calculation."),
        _ => None,
    }
}
pub fn core_basis(name: &str) -> bool {
    let n = name.trim().to_ascii_lowercase();
    let n = n.strip_prefix("aug-").unwrap_or(&n);
    ["cc-pcv", "cc-pwcv"].iter().any(|prefix| {
        n.strip_prefix(prefix).is_some_and(|tail| {
            ["dz", "tz", "qz", "5z", "6z", "dz-f12", "tz-f12", "qz-f12"].contains(&tail)
        })
    })
}
pub fn validate(
    job: JobType,
    m: QcMethod,
    charge: i32,
    mult: i32,
    c: &MethodConfig,
    atoms: &[String],
) -> Result<(), String> {
    let real_atoms: Vec<String> = atoms
        .iter()
        .filter(|a| !crate::molecule::is_dummy(a))
        .cloned()
        .collect();
    let atoms = real_atoms.as_slice();
    if let Some(e) = job_problem(Some(m), job) {
        return Err(e.into());
    }
    if c.basis.trim().is_empty()
        || c.basis.chars().any(char::is_whitespace)
        || c.basis.contains(['=', '#'])
    {
        return Err("Choose a single MRCC basis-set name.".into());
    }
    if c.mrcc.correlate_core && m.is_f12() {
        return Err("Core-correlated F12 calculations are not enabled in this MRCC interface. Choose a conventional correlated method with a cc-pCVXZ or cc-pwCVXZ orbital basis for core correlation.".into());
    }
    if c.mrcc.correlate_core && !core_basis(&c.basis) {
        return Err("Core correlation requires a core-correlation basis: choose cc-pCVXZ or cc-pwCVXZ (including augmented variants). Valence-only basis sets are not allowed with Correlate core electrons.".into());
    }
    if m.is_f12() && !c.basis.to_ascii_uppercase().ends_with("-F12") {
        return Err("Choose an F12 orbital basis for this MRCC method.".into());
    }
    if mult < 1 || atoms.is_empty() {
        return Err("MRCC requires atoms and a positive spin multiplicity.".into());
    }
    let total = atoms.iter().try_fold(-i64::from(charge), |n, a| {
        crate::molecule::atomic_number(a)
            .map(|z| n + i64::from(z))
            .ok_or_else(|| format!("Unknown MRCC element: {a}"))
    })?;
    if total < i64::from(mult - 1) || (total - i64::from(mult - 1)) % 2 != 0 {
        return Err("Charge and multiplicity are incompatible with the electron count.".into());
    }
    if job == JobType::TdDft && (mult != 1 || c.mrcc.unrestricted) {
        return Err("The MRCC ADC/TD-DFT spectrum interface currently requires a restricted closed-shell singlet.".into());
    }
    if matches!(job, JobType::Gradient | JobType::Optimize) {
        if m == QcMethod::RiMp2 && (mult != 1 || c.mrcc.unrestricted) {
            return Err(
                "MRCC RI-MP2 analytical gradients require an RHF closed-shell singlet.".into(),
            );
        }
        if m == QcMethod::Dft
            && ![
                "PBE", "PBE0", "BLYP", "B3LYP", "BP86", "TPSS", "TPSSH", "M06", "M06-2X",
            ]
            .contains(&c.functional.to_ascii_uppercase().as_str())
        {
            return Err("MRCC analytical-gradient validation currently covers PBE, PBE0, BLYP, B3LYP, BP86, TPSS, TPSSh, M06 and M06-2X. Choose one of these for optimization.".into());
        }
    }
    if matches!(m, QcMethod::Dft | QcMethod::Drpa)
        && (c.functional.trim().is_empty()
            || c.functional.chars().any(char::is_whitespace)
            || c.functional.contains(['=', '#']))
    {
        return Err("Enter a single MRCC functional name.".into());
    }
    super::job::correlated::validate_active_space(m, &c.active_space, atoms, charge, mult)
}
#[allow(clippy::too_many_arguments)]
pub fn input(
    job: JobType,
    m: QcMethod,
    charge: i32,
    mult: i32,
    c: &MethodConfig,
    excited: ExcitedStateOptions,
    orbitals: bool,
    atoms: &[String],
    xyz: &[f64],
) -> Result<String, String> {
    let (real_atoms, real_xyz) = crate::calculation_geometry::real_coordinates(atoms, xyz)?;
    let atoms = real_atoms.as_slice();
    let xyz = real_xyz.as_slice();
    validate(job, m, charge, mult, c, atoms)?;
    if xyz.len() != 3 * atoms.len() || xyz.iter().any(|v| !v.is_finite()) {
        return Err("Invalid MRCC coordinates.".into());
    }
    let calc = match m {
        QcMethod::Dft if job == JobType::TdDft => {
            if excited.tda {
                "TDA"
            } else {
                "TDDFT"
            }
        }
        QcMethod::Hf if job == JobType::TdDft => {
            if excited.tda {
                "CIS"
            } else {
                "TDHF"
            }
        }
        QcMethod::Dft => "SCF",
        QcMethod::RiMp2 => "DF-MP2",
        QcMethod::Casscf | QcMethod::Nevpt2 => "MCSCF",
        QcMethod::UccsdT => "CCSD(T)",
        _ => m.label(),
    };
    let scf = if matches!(m, QcMethod::Casscf | QcMethod::Nevpt2) {
        "MCSCF"
    } else if mult > 1 || c.mrcc.unrestricted || m == QcMethod::UccsdT {
        "UHF"
    } else {
        "RHF"
    };
    // H and He have no core shell; use the matching valence member for them.
    let basis = if core_basis(&c.basis) && atoms.iter().any(|a| matches!(a.as_str(), "H" | "He")) {
        let mut text = String::from("atomtype");
        let mut seen = std::collections::BTreeSet::new();
        for atom in atoms {
            if seen.insert(atom) {
                let name = if matches!(atom.as_str(), "H" | "He") {
                    c.basis
                        .to_ascii_lowercase()
                        .replace("cc-pwcv", "cc-pv")
                        .replace("cc-pcv", "cc-pv")
                } else {
                    c.basis.clone()
                };
                text.push_str(&format!("\n{atom}:{name}"));
            }
        }
        text
    } else {
        c.basis.clone()
    };
    let mut text=format!("calc={calc}\nbasis={}\ncharge={charge}\nmult={mult}\nscftype={scf}\nmem={}MB\ncore={}\nsymm=off\nunit=angs\nmolden={}\n",basis,c.memory_mb.max(128),if c.mrcc.correlate_core {"corr"}else{"frozen"},if orbitals{"on"}else{"off"});
    if matches!(m, QcMethod::Dft | QcMethod::Drpa) {
        text.push_str(&format!("dft={}\n", c.functional));
    }
    if job == JobType::Gradient {
        text.push_str("dens=2\n");
        if matches!(m, QcMethod::Ccsd | QcMethod::Ccsdt | QcMethod::Ccsdtq) {
            text.push_str("ccprog=mrcc\n");
        }
    }
    if job == JobType::TdDft {
        if excited.states == 0 || excited.states > 200 {
            return Err("Choose 1–200 excited states.".into());
        }
        text.push_str(&format!("nstate={}\n", excited.states + 1));
    }
    if matches!(m, QcMethod::Casscf | QcMethod::Nevpt2) {
        let electrons: i64 = atoms
            .iter()
            .map(|a| i64::from(crate::molecule::atomic_number(a).unwrap()))
            .sum::<i64>()
            - i64::from(charge);
        text.push_str(&format!(
            "docc={}\nmact={}\n",
            (electrons - i64::from(c.active_space.electrons)) / 2,
            c.active_space.orbitals
        ));
        if m == QcMethod::Nevpt2 {
            text.push_str("nevpt=on\n");
        }
    }
    text.push_str(&format!("geom=xyz\n{}\n\n", atoms.len()));
    for (a, r) in atoms.iter().zip(xyz.chunks_exact(3)) {
        text.push_str(&format!("{a} {:.14} {:.14} {:.14}\n", r[0], r[1], r[2]));
    }
    Ok(text)
}
pub fn command(binary: &Path, dir: &Path, nproc: u32) -> Result<Command, String> {
    // Absolute selected executable; companion programs are resolved only in the child.
    let executable = std::fs::canonicalize(binary).map_err(|e| {
        format!(
            "Cannot open the selected MRCC executable {}: {e}",
            binary.display()
        )
    })?;
    let mut paths = vec![executable
        .parent()
        .ok_or("MRCC executable has no parent directory")?
        .to_owned()];
    if let Some(path) = std::env::var_os("PATH") {
        paths.extend(std::env::split_paths(&path));
    }
    let mut cmd = Command::new(executable);
    cmd.current_dir(dir)
        .env(
            "PATH",
            std::env::join_paths(paths).map_err(|e| e.to_string())?,
        )
        .env("OMP_NUM_THREADS", nproc.max(1).to_string())
        .env("MKL_NUM_THREADS", nproc.max(1).to_string());
    Ok(cmd)
}

/// Scientific settings stay controlled by the panel; a manual input may adjust
/// convergence/printing without silently changing the job or bypassing core validation.
pub fn validate_custom(generated: &str, custom: &str) -> Result<(), String> {
    fn sections(s: &str) -> Result<(std::collections::BTreeMap<String, String>, String), String> {
        let mut keys = std::collections::BTreeMap::new();
        let mut geom = String::new();
        let mut geometry = false;
        for line in s.lines() {
            let l = line.split('#').next().unwrap_or("").trim();
            if l.is_empty() {
                continue;
            }
            if geometry {
                geom.push_str(&l.split_whitespace().collect::<Vec<_>>().join(" "));
                geom.push('\n');
                continue;
            }
            if keys.get("basis").is_some_and(|v| v == "atomtype") && !l.contains('=') {
                let (atom, basis) = l
                    .split_once(':')
                    .ok_or("MRCC atomtype basis lines must use element:basis.")?;
                let key = format!("basis:{}", atom.trim().to_ascii_lowercase());
                if keys
                    .insert(key, basis.trim().to_ascii_lowercase())
                    .is_some()
                {
                    return Err("Duplicate MRCC element basis.".into());
                }
                continue;
            }
            let (k, v) = l
                .split_once('=')
                .ok_or("MRCC input must use keyword=value before geom=xyz.")?;
            let key = k.trim().to_ascii_lowercase();
            let value = v.trim().to_ascii_lowercase();
            if keys.insert(key.clone(), value).is_some() {
                return Err(format!("Duplicate MRCC keyword: {key}"));
            }
            geometry = key == "geom";
        }
        Ok((keys, geom))
    }
    let (expected, geometry) = sections(generated)?;
    let (actual, custom_geometry) = sections(custom)?;
    if geometry != custom_geometry {
        return Err(
            "Change the geometry in the viewer, then sync the MRCC input from settings.".into(),
        );
    }
    for (key, value) in &expected {
        if actual.get(key) != Some(value) {
            return Err(format!("Set MRCC {key} through the panel and sync the input. This keeps method, geometry, charge/spin, core and gradient validation consistent."));
        }
    }
    for key in actual.keys() {
        if !expected.contains_key(key)
            && ![
                "scfmaxit",
                "ccmaxit",
                "scftol",
                "cctol",
                "scfiguess",
                "verbosity",
            ]
            .contains(&key.as_str())
        {
            return Err(format!("Extra MRCC keyword '{key}' is not enabled in the input editor. Scientific job controls must be set in the panel."));
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "mrcc_tests.rs"]
mod tests;
