//! PySCF energy-only CCSD(T), CASSCF and strongly contracted NEVPT2 scripts.
//! API references: https://pyscf.org/user/cc.html and https://pyscf.org/user/mrpt.html
use super::*;
use crate::qchem_interfaces::job::{correlated, ExcitedStateOptions, QcMethod};

#[allow(clippy::too_many_arguments)]
pub fn panel_script_text(
    job: JobType,
    selected: QcMethod,
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    excited: ExcitedStateOptions,
    wavefunction: bool,
    atoms: &[String],
    coordinates: &[f64],
) -> String {
    if !correlated::energy_only(crate::qchem_interfaces::program::QcProgram::PySCF, selected) {
        return script_text(
            job,
            charge,
            multiplicity,
            config,
            excited,
            wavefunction,
            atoms,
            coordinates,
        );
    }
    if job != JobType::SinglePoint {
        return "raise RuntimeError('This correlated method is available for Energy calculations in Beavyr.')\n".into();
    }
    let mut reference = config.clone();
    reference.functional = "hf".into();
    let mut script = script_text(
        job,
        charge,
        multiplicity,
        &reference,
        excited,
        false,
        atoms,
        coordinates,
    );
    if correlated::needs_active_space(selected) {
        script = script.replace("scf.UHF(molecule)", "scf.ROHF(molecule)");
    }
    script = script.replace(
        "\nBOHR =",
        &format!("\nmol.max_memory = {}\n\nBOHR =", config.memory_mb),
    );
    let mut correction = match selected {
        QcMethod::CcsdT => r#"
from pyscf import cc
correlation = cc.CCSD(mf)
correlation.kernel()
if not correlation.converged:
    raise RuntimeError("PySCF CCSD did not converge; CCSD(T) energy is unavailable.")
energy = float(correlation.e_tot + correlation.ccsd_t())
"#
        .to_owned(),
        QcMethod::Casscf | QcMethod::Nevpt2 => {
            let active = &config.active_space;
            format!(
                r#"
from pyscf import mcscf
active_electrons = {electrons}
active_orbitals = {orbitals}
if active_electrons <= 0 or active_orbitals <= 0:
    raise RuntimeError("Set the active electrons and active orbitals.")
if (active_electrons - mol.spin) % 2 or active_electrons < mol.spin:
    raise RuntimeError("Active electrons are incompatible with the spin multiplicity.")
active_alpha = (active_electrons + mol.spin) // 2
active_beta = (active_electrons - mol.spin) // 2
cas = mcscf.CASSCF(mf, active_orbitals, (active_alpha, active_beta))
cas.fix_spin_(ss=(mol.spin / 2) * (mol.spin / 2 + 1))
cas.kernel()
if not cas.converged:
    raise RuntimeError("PySCF CASSCF did not converge; no correlated energy is accepted.")
energy = float(cas.e_tot)
"#,
                electrons = active.electrons,
                orbitals = active.orbitals
            )
        }
        _ => unreachable!("energy-only PySCF methods handled above"),
    };
    if selected == QcMethod::Nevpt2 {
        correction
            .push_str("\nfrom pyscf import mrpt\nenergy += float(mrpt.NEVPT(cas).kernel())\n");
    }
    if wavefunction {
        correction.push_str("\nfrom pyscf.tools import molden\n");
        if correlated::needs_active_space(selected) {
            correction.push_str("molden.from_mcscf(cas, 'orbitals.molden.input')\n");
        } else {
            correction.push_str("molden.from_scf(mf, 'orbitals.molden.input')\n");
        }
    }
    correction.push_str("\nif not numpy.isfinite(energy):\n    raise RuntimeError('PySCF returned a non-finite correlated energy.')\n");
    // Replace the single final report, so it contains the correlated total, never the HF energy.
    script.replace(
        "\nprint(\"energy %.12f\" % energy)",
        &format!("{correction}\nprint(\"energy %.12f\" % energy)"),
    )
}
