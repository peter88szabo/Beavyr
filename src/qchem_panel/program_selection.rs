//! Backend changes keep the selected job within the backend's capabilities.
use super::{JobType, QcPanelState, QcProgram};
impl QcPanelState {
    /// Switches program, resetting everything that does not carry across.
    ///
    /// The method, functional and basis all reset, because the catalogues do
    /// not translate between codes: ORCA's `M062X` is an input error to
    /// Behemoth and Behemoth's `m06-2x` is one to ORCA. The charge, the
    /// multiplicity and the resources stay, because they describe the molecule
    /// and the machine rather than the program.
    ///
    /// The job stays too where the new program can do it, and falls back to an
    /// optimisation where available, or a spectrum for standalone sTDA-xTB.
    pub fn set_program(&mut self, program: QcProgram) {
        if program == self.program {
            return;
        }
        self.program = program;
        self.method = crate::qchem_interfaces::job::default_method(program);
        crate::qchem_interfaces::method::adopt_program_defaults(&mut self.config, program);
        if crate::qchem_interfaces::job::correlated::unsupported_job(program, self.method, self.job).is_some() {
            self.job = if program == QcProgram::StdaXtb {
                JobType::TdDft
            } else {
                JobType::Optimize
            };
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn spectrum_only_selection_chooses_a_supported_job_immediately() {
        let mut panel = QcPanelState::default();
        panel.set_program(QcProgram::StdaXtb);
        assert_eq!(panel.job, JobType::TdDft);
        panel.set_program(QcProgram::Xtb);
        assert_eq!(panel.job, JobType::Optimize);
        panel.set_program(QcProgram::StdaXtb);
        panel.set_program(QcProgram::Behemoth);
        assert_eq!(panel.job, JobType::TdDft);
    }
}
