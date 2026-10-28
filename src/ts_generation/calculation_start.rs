//! Validate and start a reaction-path calculation using real atoms.

use super::*;

impl TsGeneration {
    pub fn input_options(&self) -> Result<RdaOptions> {
        let (reactant, product, parameters) = super::endpoint_preparation::prepare(self)?;
        validate_endpoints(&reactant, &product)?;
        if self.program == QcProgram::Dreiding {
            bail!("Choose xTB or Behemoth for TS generation.");
        }
        let mol = Molecule {
            atoms: reactant.atoms.clone(),
            pos: reactant.positions.clone(),
            bonds: vec![],
            hydrogen_bonds: vec![],
            topology: None,
        };
        valence::validate_electronic_state(&mol, self.charge, self.multiplicity)
            .map_err(|error| anyhow::anyhow!("{error}"))?;
        let issues = method::validate(
            self.program,
            &self.method,
            self.multiplicity,
            &reactant.atoms,
        );
        if let Some(issue) = issues
            .iter()
            .find(|i| i.severity == method::Severity::Block)
        {
            bail!("{}", issue.message);
        }
        if self.method.nproc == 0 {
            bail!("Use at least one calculation thread.");
        }
        parameters.rda_options(reactant.atoms.len())
    }

    pub fn start(&mut self) -> Result<()> {
        if self.is_running() {
            bail!("A TS generation is already running.");
        }
        let options = self.input_options()?;
        let (prepared_reactant, prepared_product, prepared_parameters) =
            super::endpoint_preparation::prepare(self)?;
        let configured = self.path_for(self.program).to_string();
        let configured = &configured;
        let binary =
            resolve_program_executable(self.program, configured).map_err(anyhow::Error::msg)?;
        // Subprocesses change directory; resolve relative executable paths first.
        let binary =
            std::fs::canonicalize(binary).context("Could not resolve the executable path")?;
        let backend = BackendConfig {
            program: self.program,
            binary,
            method: self.method.clone(),
            charge: self.charge,
            multiplicity: self.multiplicity,
        };
        config::save_path(self.program, configured);
        let reactant = prepared_reactant;
        let product = prepared_product;
        self.reactant = Some(reactant.clone());
        self.product = Some(product.clone());
        self.parameters = prepared_parameters;
        let algorithm = self.parameters.algorithm;
        let method = format!(
            "{} / {}",
            self.program.label(),
            method_ui::method_summary(self.program, &self.method)
        );
        let run_dir = create_xtb_run_dir(&std::env::temp_dir().join("beavyr_ts"))
            .map_err(anyhow::Error::msg)?;
        self.control = Arc::new(RunControl::default());
        let control = self.control.clone();
        self.task = Some(AsyncComputeTaskPool::get().spawn(async move {
            let mut objective = EnergyGradient::new(
                backend,
                reactant.atoms.clone(),
                run_dir.join("energy"),
                control.clone(),
            );
            let result = generate(
                &mut objective,
                &reactant,
                &product,
                options,
                algorithm,
                method,
                &run_dir,
            );
            if control.cancelled.load(std::sync::atomic::Ordering::Relaxed) {
                return Err(format!(
                    "TS generation cancelled. Scratch files: {}",
                    run_dir.display()
                ));
            }
            match result {
                Ok(output) => {
                    let _ = std::fs::remove_dir_all(&run_dir);
                    Ok(output)
                }
                Err(error) => Err(format!("{error:#}\nScratch files: {}", run_dir.display())),
            }
        }));
        self.started = Some(Instant::now());
        self.duration = None;
        self.output = None;
        self.energy_plot_open = false;
        self.message = None;
        self.is_error = false;
        Ok(())
    }
}
