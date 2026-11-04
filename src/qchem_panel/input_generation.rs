//! Input preview shared with the actual calculation runners.
use super::*;
impl QcPanelState {
    /// The name the input is written under, for the save dialog's suggestion.
    pub fn input_file_name(&self) -> &'static str {
        match self.program {
            QcProgram::Mrcc => "MINP",
            QcProgram::Orca => "orca_job.inp",
            QcProgram::Psi4 => "input.dat",
            QcProgram::Psi4Py => "run_psi4.py",
            QcProgram::PySCF => "run_pyscf.py",
            _ => "input.txt",
        }
    }


    /// Whether this program is genuinely driven by an input file.
    ///
    /// Four are: ORCA by its `.inp`, Psi4 by an input file that is itself
    /// Python, and PySCF by a Python script -- which is simply how PySCF is
    /// used. Each of those is a file a chemist would recognise, write by hand
    /// and run outside Beavyr.
    ///
    /// Three are not. xTB and Behemoth take everything on their command lines,
    /// and DREIDING runs in this process. Sparrow belongs with them: it is
    /// driven by settings, through its own command line or through SCINE's
    /// bindings, and the Python script Beavyr generates for it is this
    /// program's scaffolding rather than Sparrow's input format. Offering that
    /// scaffolding as "the input" would invite editing something that is not a
    /// thing Sparrow has.

    /// Scans generate one or more inputs from their own settings.
    pub fn has_input_file(&self) -> bool {
        !matches!(self.job, JobType::Scan | JobType::Mecp) && !crate::qchem_interfaces::job::search::local(self.program, self.job, &self.config.search) && matches!(
            self.program,
            QcProgram::Mrcc | QcProgram::Orca | QcProgram::Psi4 | QcProgram::Psi4Py | QcProgram::PySCF
        )
    }


    /// The input the current settings would write.
    ///
    /// One function for three callers: the preview the editor starts from, the
    /// Save button, and the run itself. They cannot disagree about what would
    /// have been sent, which is the point -- an editor showing something other
    /// than what runs would be worse than no editor.
    pub fn generated_input(&self, atoms: &[String], positions_angstrom: &[f64]) -> String {
        use crate::qchem_interfaces::{orca_run, psi4_run, pyscf_run};

        let method = self.method.unwrap_or(QcMethod::Dft);
        match self.program {
            QcProgram::Mrcc => crate::qchem_interfaces::mrcc::input(self.job, method, self.charge, self.multiplicity, &self.config, self.excited, self.wavefunction, atoms, positions_angstrom).unwrap_or_else(|e| format!("# {e}")),
            QcProgram::Orca => orca_run::panel_input_text(
                self.job,
                method,
                self.charge,
                self.multiplicity,
                &self.config,
                self.excited,
                &self.config.orca_extra,
                &orca_run::xyz_body(atoms, positions_angstrom),
            ),
            QcProgram::Psi4 => psi4_run::input_text(
                self.job,
                method,
                self.charge,
                self.multiplicity,
                &self.config,
                self.excited,
                self.wavefunction,
                atoms,
                positions_angstrom,
            ),
            QcProgram::Psi4Py => psi4_run::script_text(
                self.job,
                method,
                self.charge,
                self.multiplicity,
                &self.config,
                self.excited,
                self.wavefunction,
                atoms,
                positions_angstrom,
            ),
            QcProgram::PySCF => pyscf_run::panel_script_text(
                self.job,
                method,
                self.charge,
                self.multiplicity,
                &self.config,
                self.excited,
                self.wavefunction,
                atoms,
                positions_angstrom,
            ),
            // Everything else is driven by settings rather than a file, so
            // there is nothing to show. See `has_input_file`.
            _ => String::new(),
        }
    }
}
