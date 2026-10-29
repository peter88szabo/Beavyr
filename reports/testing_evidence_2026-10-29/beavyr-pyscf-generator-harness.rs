// Standalone test adapter: compile the unchanged production script generator.
// Minimal configuration types let this run without changing the Beavyr crate.
mod qchem_interfaces {
    pub mod job {
        #[derive(Clone, Copy, PartialEq)]
        pub enum JobType { SinglePoint, Gradient, Hessian, TdDft, Optimize }
        impl JobType {
            pub fn produces_trajectory(self) -> bool { self == Self::Optimize }
            pub fn produces_frequencies(self) -> bool { false }
        }
        pub struct ExcitedStateOptions { pub states: usize, pub tda: bool }
    }
    pub mod method {
        pub struct MethodConfig {
            pub functional: String, pub basis: String, pub nproc: usize,
        }
    }
    pub mod python_env { pub mod launch {
        pub fn command(_: &str, args: [&str; 1]) -> std::process::Command {
            let mut command = std::process::Command::new("python3");
            command.args(args);
            command
        }
    }}
    #[path = "/home/peter/DevCode/Beavyr/src/qchem_interfaces/pyscf_method.rs"]
    pub mod pyscf_method;
    #[path = "/home/peter/DevCode/Beavyr/src/qchem_interfaces/pyscf_run.rs"]
    pub mod pyscf_run;
}
fn main() {
    use qchem_interfaces::{job::{JobType, ExcitedStateOptions}, method::MethodConfig, pyscf_run};
    let root = std::path::PathBuf::from(std::env::args().nth(1).unwrap());
    for (name, job, functional) in [
        ("energy", JobType::SinglePoint, "hf"),
        ("gradient", JobType::Gradient, "hf"),
        ("hessian", JobType::Hessian, "hf"),
        ("dft_gradient", JobType::Gradient, "b3lyp"),
        ("optimization_missing_dependency", JobType::Optimize, "hf"),
    ] {
        let directory = root.join(name);
        std::fs::create_dir_all(&directory).unwrap();
        let script = pyscf_run::script_text(job, 0, 1,
            &MethodConfig { functional: functional.into(), basis: "sto-3g".into(), nproc: 2 },
            ExcitedStateOptions { states: 3, tda: false }, false,
            &["O".into(), "H".into(), "H".into()],
            &[0., 0., 0., 0., 0., 0.96, 0.92, 0., -0.24]);
        std::fs::write(directory.join(pyscf_run::SCRIPT_FILE), script).unwrap();
    }
}
