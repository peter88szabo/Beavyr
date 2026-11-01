//! Explicit native-optimizer overrides. Absence preserves the backend default.
use crate::qchem_interfaces::program::QcProgram;
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Criteria {
    pub energy: f64,
    pub max_step: f64,
    pub rms_step: f64,
    pub max_gradient: f64,
    pub rms_gradient: f64,
}
impl Criteria {
    pub fn validate(&self) -> Result<(), String> {
        if [
            self.energy,
            self.max_step,
            self.rms_step,
            self.max_gradient,
            self.rms_gradient,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.)
        {
            Err("All tolerances must be positive and finite.".into())
        } else {
            Ok(())
        }
    }
    /// Common storage uses bohr, Hartree/bohr and Hartree.
    pub fn defaults(program: QcProgram) -> Self {
        let mut result = Self {
            energy: 5e-6,
            max_step: 4e-3,
            rms_step: 2e-3,
            max_gradient: 3e-4,
            rms_gradient: 1e-4,
        };
        match program {
            QcProgram::PySCF => {
                result.energy = 1e-6;
                result.max_gradient = 4.5e-4;
                result.rms_gradient = 3e-4;
                result.max_step = 1.8e-3 / 0.529177210903;
                result.rms_step = 1.2e-3 / 0.529177210903;
            }
            QcProgram::Dreiding => {
                let o = crate::qchem_interfaces::dreiding_run::optimization_options();
                result = Self {
                    energy: o.energy_tol,
                    max_step: o.max_step,
                    rms_step: o.rms_step,
                    max_gradient: o.max_gradient,
                    rms_gradient: o.rms_gradient,
                };
            }
            QcProgram::Psi4 | QcProgram::Psi4Py => {
                result.energy = 1e-6;
                result.max_step = 1.2e-3;
                // RMS criteria are not active in the default QCHEM preset;
                // these are explicit custom thresholds only when saved.
            }
            _ => {}
        }
        result
    }
    pub fn orca(&self) -> String {
        format!(" TolE {:.12e}\n TolMaxD {:.12e}\n TolRMSD {:.12e}\n TolMaxG {:.12e}\n TolRMSG {:.12e}\n",
            self.energy, self.max_step, self.rms_step, self.max_gradient, self.rms_gradient)
    }
    pub fn behemoth(&self, command: &mut std::process::Command) {
        for (flag, value) in [
            ("--opt-energy-tol", self.energy),
            ("--opt-max-step", self.max_step),
            ("--opt-rms-step", self.rms_step),
            ("--opt-max-gradient", self.max_gradient),
            ("--opt-rms-gradient", self.rms_gradient),
        ] {
            command.arg(flag).arg(value.to_string());
        }
    }
}
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct Settings {
    pub native: Vec<(QcProgram, Criteria)>,
    pub xtb_level: Option<String>,
}
impl Settings {
    fn key(program: QcProgram) -> QcProgram {
        if program == QcProgram::Psi4Py {
            QcProgram::Psi4
        } else {
            program
        }
    }
    pub fn get(&self, program: QcProgram) -> Option<&Criteria> {
        self.native
            .iter()
            .find(|(p, _)| *p == Self::key(program))
            .map(|(_, t)| t)
    }
    pub fn set(&mut self, program: QcProgram, value: Option<Criteria>) {
        self.native.retain(|(p, _)| *p != Self::key(program));
        if let Some(value) = value {
            self.native.push((Self::key(program), value));
        }
    }
}
