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
