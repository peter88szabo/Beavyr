//! Capabilities and active-space inputs for the calculator's correlated energy jobs.
use super::{JobType, QcMethod};
use crate::qchem_interfaces::program::QcProgram;

#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ActiveSpace {
    pub electrons: u32,
    pub orbitals: u32,
}

pub fn needs_active_space(method: QcMethod) -> bool {
    matches!(method, QcMethod::Casscf | QcMethod::Nevpt2)
}

pub fn energy_only(program: QcProgram, method: QcMethod) -> bool {
    needs_active_space(method)
        || (method == QcMethod::CcsdT && matches!(program, QcProgram::PySCF | QcProgram::Behemoth))
}

pub fn unsupported_job(
    program: QcProgram,
    method: Option<QcMethod>,
    job: JobType,
) -> Option<&'static str> {
    if program == QcProgram::Mrcc { return super::super::mrcc::job_problem(method, job); }
    if method.is_some_and(|m| energy_only(program, m)) && !matches!(job, JobType::SinglePoint | JobType::Scan) {
        return Some("This method is available for Energy calculations in Beavyr.");
    }
    super::unsupported_reason(program, job)
}

pub fn validate_active_space(
    method: QcMethod,
    active: &ActiveSpace,
    atoms: &[String],
    charge: i32,
    multiplicity: i32,
) -> Result<(), String> {
    if !needs_active_space(method) {
        return Ok(());
    }
    if active.electrons == 0 || active.orbitals == 0 {
        return Err("Set the active electrons and active orbitals for CASSCF/NEVPT2.".into());
    }
    let electrons = i64::from(active.electrons);
    let orbitals = i64::from(active.orbitals);
    let unpaired = i64::from(multiplicity) - 1;
    if unpaired < 0
        || electrons < unpaired
        || (electrons - unpaired) % 2 != 0
        || (electrons + unpaired) / 2 > orbitals
    {
        return Err(
            "The active electrons, active orbitals and multiplicity are incompatible.".into(),
        );
    }
    let mut total = -i64::from(charge);
    for atom in atoms.iter().filter(|a| !crate::molecule::is_dummy(a)) {
        total += crate::molecule::atomic_number(atom)
            .map(i64::from)
            .ok_or_else(|| format!("Unknown element {atom} in the active-space calculation."))?;
    }
    if electrons > total || (total - electrons) % 2 != 0 {
        return Err("Active electrons must not exceed the molecule's electron count; inactive electrons must form pairs.".into());
    }
    Ok(())
}

pub fn orca_active_space_block(
    method: QcMethod,
    multiplicity: i32,
    active: &ActiveSpace,
) -> String {
    if !needs_active_space(method) {
        return String::new();
    }
    format!(
        "%casscf\n  nel {}\n  norb {}\n  mult {}\n  nroots 1\n{}end\n",
        active.electrons,
        active.orbitals,
        multiplicity,
        if method == QcMethod::Nevpt2 {
            "  PTMethod SC_NEVPT2\n"
        } else {
            ""
        }
    )
}

/// SCF reference of a CCSD(T) calculation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Reference {
    Rhf,
    Rohf,
    Uhf,
}
impl Reference {
    pub fn label(self) -> &'static str {
        match self {
            Self::Rhf => "RHF",
            Self::Rohf => "ROHF",
            Self::Uhf => "UHF",
        }
    }
    /// The choices that make sense for this spin state.
    pub fn options(multiplicity: i32) -> [Self; 2] {
        if multiplicity > 1 {
            [Self::Uhf, Self::Rohf]
        } else {
            [Self::Rhf, Self::Uhf]
        }
    }
}

/// Programs where the CCSD(T) reference can be chosen.
pub fn reference_choice(program: QcProgram, method: Option<QcMethod>) -> bool {
    method == Some(QcMethod::CcsdT)
        && matches!(program, QcProgram::Orca | QcProgram::Mrcc | QcProgram::Behemoth)
}

/// The reference actually used: the stored choice when it fits the spin
/// state, otherwise RHF for a closed shell and UHF for an open shell (ORCA's
/// default).
pub fn reference(chosen: Option<Reference>, multiplicity: i32) -> Reference {
    let open = multiplicity > 1;
    match chosen {
        Some(Reference::Uhf) => Reference::Uhf,
        Some(Reference::Rohf) if open => Reference::Rohf,
        Some(Reference::Rhf) if !open => Reference::Rhf,
        _ if open => Reference::Uhf,
        _ => Reference::Rhf,
    }
}

/// Programs and methods where the number of correlated electrons can be set.
pub fn correlated_electrons_choice(program: QcProgram, method: Option<QcMethod>) -> bool {
    let Some(m) = method else { return false };
    match program {
        QcProgram::Orca => matches!(
            m,
            QcMethod::Mp2
                | QcMethod::RiMp2
                | QcMethod::Mp2F12
                | QcMethod::CcsdT
                | QcMethod::DlpnoCcsdT
                | QcMethod::CcsdTF12
                | QcMethod::Nevpt2
        ),
        QcProgram::Mrcc => !matches!(m, QcMethod::Hf | QcMethod::Dft | QcMethod::Casscf),
        _ => false,
    }
}

/// Number of electrons of the real atoms at this charge.
pub fn total_electrons(atoms: &[String], charge: i32) -> Result<i64, String> {
    atoms
        .iter()
        .filter(|a| !crate::molecule::is_dummy(a))
        .try_fold(-i64::from(charge), |n, a| {
            crate::molecule::atomic_number(a)
                .map(|z| n + i64::from(z))
                .ok_or_else(|| format!("Unknown element: {a}"))
        })
}

/// Electrons to freeze for `correlated` correlated electrons, or `None` to
/// keep the program's own frozen-core default.
pub fn frozen_electrons(
    correlated: Option<u32>,
    atoms: &[String],
    charge: i32,
    multiplicity: i32,
) -> Result<Option<u32>, String> {
    let Some(n) = correlated else { return Ok(None) };
    let total = total_electrons(atoms, charge)?;
    let n = i64::from(n);
    let unpaired = i64::from(multiplicity.max(1) - 1);
    if n > total || n < unpaired.max(1) || (total - n) % 2 != 0 {
        return Err(format!(
            "Correlated electrons must lie between {} and {total} and leave an even number of frozen core electrons.",
            unpaired.max(1)
        ));
    }
    Ok(Some((total - n) as u32))
}
