//! Writing ORCA's input, running it, and reading what it leaves behind.
//!
//! The two files that matter most are already read elsewhere:
//! [`super::orca_engrad`] takes the Cartesian gradient and
//! [`super::orca_hess`] takes the Hessian, the geometry, the masses and the
//! infrared spectrum. Both were written for files a user had produced
//! themselves and both work unchanged on files Beavyr has just caused to be
//! produced, so nothing here parses either of them a second time.
//!
//! What is here is the input file, the run, and the optimisation trajectory.
//!
//! Verified against ORCA 6.1.1 by running an `Opt Freq` job on water at
//! HF/STO-3G and reading what appeared.

use std::path::{Path, PathBuf};
use std::process::Command;

use super::job::{ExcitedStateOptions, JobType, QcMethod};
use super::method::MethodConfig;
use super::orca_method;

/// The basename every file of a run is built from.
///
/// Deliberately not `input`: xTB is given `input.xyz` as its geometry in the
/// same kind of scratch directory, and ORCA writes its *optimised* geometry to
/// `<basename>.xyz`. Sharing a basename would have one overwrite the other.
pub const JOB_BASENAME: &str = "orca_job";

/// What ORCA prints when a job has failed, whatever its exit status says.
///
/// ORCA can exit zero on a failed run, so the status alone is not enough. This
/// is the string the Smite interface keys on for the same reason.
const ERROR_TERMINATION: &str = "ORCA finished by error termination";

/// What ORCA prints when a job has genuinely finished.
const NORMAL_TERMINATION: &str = "ORCA TERMINATED NORMALLY";

/// What ORCA prints once a geometry optimisation has converged.
const CONVERGED: &str = "THE OPTIMIZATION HAS CONVERGED";

/// Paths to everything a finished run may have written.
#[derive(Debug, Clone)]
pub struct JobFiles {
    pub input: PathBuf,
    pub output: PathBuf,
    /// The optimisation path, written only by a job that optimised.
    pub trajectory: PathBuf,
    /// The final geometry, written only by a job that optimised.
    pub final_geometry: PathBuf,
    pub hessian: PathBuf,
    pub engrad: PathBuf,
}

/// Where each file of a run in `run_dir` lives.
pub fn job_files(run_dir: &Path) -> JobFiles {
    JobFiles {
        input: run_dir.join(format!("{JOB_BASENAME}.inp")),
        output: run_dir.join(format!("{JOB_BASENAME}.out")),
        trajectory: run_dir.join(format!("{JOB_BASENAME}_trj.xyz")),
        final_geometry: run_dir.join(format!("{JOB_BASENAME}.xyz")),
        hessian: run_dir.join(format!("{JOB_BASENAME}.hess")),
        engrad: run_dir.join(format!("{JOB_BASENAME}.engrad")),
    }
}

/// Builds ORCA's input file.
///
/// The layout follows the Smite interface, which is ORCA's ordinary one: the
/// parallel block, the memory block, the keyword line, then the coordinates
/// between `* xyz` and `*`.
///
/// `xyz_body` is one `symbol x y z` line per atom, in angstrom, with no count
/// and no comment line: ORCA's coordinate block is not an XYZ file.
pub fn input_text(
    job: JobType,
    charge: i32,
    multiplicity: i32,
    method: &MethodConfig,
    extra_keywords: &str,
    xyz_body: &str,
) -> String {
    let multiplicity = multiplicity.max(1);
    let nproc = method.nproc.max(1);

    let mut text = String::new();
    text.push_str(&format!("%pal nprocs {nproc} end\n"));
    // ORCA's `%maxcore` is memory *per core*, where the panel asks for a total.
    // Dividing here is what makes the number on screen mean what it says; an
    // undivided figure would let an eight-thread job take eight times what was
    // asked for. A floor keeps a small total from becoming an unusable share.
    let per_core = (method.memory_mb / nproc).max(256);
    text.push_str(&format!("%maxcore {per_core}\n"));

    text.push_str("! ");
    text.push_str(orca_method::functional(&method.functional).map_or(
        method.functional.as_str(),
        |f| f.keyword,
    ));
    // A composite brings its own basis, so passing one alongside would be both
    // redundant and contradictory.
    if !orca_method::is_composite(&method.functional) {
        text.push(' ');
        text.push_str(&method.basis);
    }
    // A double hybrid computes an MP2-like correlation energy, which ORCA
    // refuses to do without a correlation-fitting basis. Tested directly:
    // `! B2PLYP def2-SVP` is an input error and `! B2PLYP def2-SVP def2-SVP/C`
    // runs. Adding it here is what makes those functionals usable at all.
    if let Some(aux) = orca_method::auxiliary_basis(&method.functional, &method.basis) {
        text.push(' ');
        text.push_str(&aux);
    }
    if let Some(dispersion) =
        orca_method::dispersion_keyword(&method.functional, method.orca_dispersion)
    {
        text.push(' ');
        text.push_str(dispersion);
    }
    let job_words = job_keywords(job);
    if !job_words.is_empty() {
        text.push(' ');
        text.push_str(job_words);
    }
    let extra = extra_keywords.trim();
    if !extra.is_empty() {
        text.push(' ');
        text.push_str(extra);
    }
    text.push('\n');

    text.push_str(&format!("* xyz {charge} {multiplicity}\n"));
    text.push_str(xyz_body);
    if !xyz_body.ends_with('\n') {
        text.push('\n');
    }
    text.push_str("*\n");
    text
}

/// The basis sets an explicitly correlated method needs.
///
/// F12 will not run on an ordinary basis: it needs one of the purpose-built
/// `cc-pVnZ-F12` sets, and alongside it a complementary auxiliary basis, the
/// `-CABS` set of the same name. Both are supplied from one choice, because
/// they always travel together and nothing is gained by asking twice.
pub const F12_BASIS_SETS: [(&str, &str); 3] = [
    ("cc-pVDZ-F12", "cc-pVDZ-F12"),
    ("cc-pVTZ-F12", "cc-pVTZ-F12"),
    ("cc-pVQZ-F12", "cc-pVQZ-F12"),
];

/// The default orbital basis for an F12 method.
pub const DEFAULT_F12_BASIS: &str = "cc-pVDZ-F12";

/// The keyword that selects `method`, or the functional where the method is
/// DFT.
///
/// A composite or an ordinary functional travels as its own keyword; everything
/// else is a method name ORCA knows directly.
pub fn method_keyword<'a>(method: QcMethod, functional: &'a str) -> &'a str {
    match method {
        QcMethod::Dft => orca_method::functional(functional).map_or(functional, |f| f.keyword),
        QcMethod::Hf => "HF",
        QcMethod::Mp2 => "MP2",
        QcMethod::RiMp2 => "RI-MP2",
        QcMethod::Mp2F12 => "MP2-F12",
        QcMethod::CcsdT => "CCSD(T)",
        QcMethod::DlpnoCcsdT => "DLPNO-CCSD(T)",
        QcMethod::CcsdTF12 => "CCSD(T)-F12",
        // Not ORCA methods; the panel never offers them for this program.
        QcMethod::Gfn1
        | QcMethod::Gfn2
        | QcMethod::GfnFf
        | QcMethod::Tasi
        | QcMethod::Semiempirical => "HF",
    }
}

/// Every basis-related keyword this method and basis need, in the order they go
/// on the `!` line.
///
/// One function because the rule is not "the basis": a correlated method also
/// needs a correlation-fitting set, and an F12 method needs a complementary
/// auxiliary basis on top of that. Leaving any of them out is an ORCA input
/// error rather than a slower calculation, so they are derived rather than
/// asked for.
pub fn basis_keywords(method: QcMethod, functional: &str, basis: &str) -> Vec<String> {
    let mut keywords = Vec::new();
    if !method.needs_basis() {
        return keywords;
    }
    // A composite functional brings its own basis, so nothing is added.
    if method == QcMethod::Dft && orca_method::is_composite(functional) {
        return keywords;
    }

    let basis = basis.trim();
    if basis.is_empty() {
        return keywords;
    }
    keywords.push(basis.to_string());

    if method.is_f12() {
        // The complementary auxiliary basis, named after the orbital basis.
        keywords.push(format!("{basis}-CABS"));
        // And a correlation-fitting set. The F12 orbital bases have no `/C` of
        // their own, so a matching ordinary one is used.
        keywords.push("cc-pVTZ/C".to_string());
    } else if method.needs_auxiliary_basis() {
        keywords.push(format!("{basis}/C"));
    } else if method == QcMethod::Dft {
        // A double hybrid carries an MP2-like term and needs the same set.
        if let Some(aux) = orca_method::auxiliary_basis(functional, basis) {
            keywords.push(aux);
        }
    }

    keywords
}

/// The keywords that select `job` on the `!` line.
///
/// A single point and a TD-DFT run add nothing: ORCA computes an energy by
/// default, and excited states are asked for by the `%tddft` block alone.
/// `! TDDFT` is an input error, which is worth stating because it is the
/// obvious guess.
pub fn job_keywords(job: JobType) -> &'static str {
    match job {
        JobType::SinglePoint | JobType::TdDft => "",
        // `EnGrad` writes the energy and the Cartesian gradient to `.engrad`,
        // which Beavyr already has a reader for.
        JobType::Gradient => "EnGrad",
        // There is no "matrix only" keyword: `Freq` is what produces the
        // `.hess`, and the difference between this job and Frequencies is what
        // Beavyr then does with the file, not what ORCA was asked for.
        JobType::Hessian => "Freq",
        JobType::Optimize => "Opt",
        JobType::Frequencies => "Freq",
        JobType::OptimizeThenFrequencies => "Opt Freq",
        JobType::TransitionState => "OptTS",
        JobType::Irc => "IRC",
    }
}

/// The `%` blocks a job needs under its keyword line.
///
/// Two jobs need one. A transition-state search wants a real Hessian to start
/// from and again periodically, or it wanders off to whatever saddle is
/// nearest; those settings are the ones the Smite interface already used. An
/// excited-state run is requested entirely by its block.
pub fn job_blocks(job: JobType, excited: ExcitedStateOptions) -> String {
    match job {
        JobType::TransitionState => concat!(
            "%geom\n",
            // An exact Hessian at the first step, and again every five, which
            // is what makes an eigenvector-following search land on the saddle
            // it was aimed at rather than the nearest one downhill.
            "  Calc_Hess true\n",
            "  Recalc_Hess 5\n",
            "  MaxIter 100\n",
            "end\n"
        )
        .to_string(),
        JobType::TdDft => format!(
            "%tddft\n  nroots {}\n  tda {}\nend\n",
            excited.states.max(1),
            if excited.tda { "true" } else { "false" }
        ),
        _ => String::new(),
    }
}

/// The coordinate block ORCA wants: no atom count, no comment line.
pub fn xyz_body(atoms: &[String], positions_angstrom: &[f64]) -> String {
    let mut body = String::new();
    for (index, symbol) in atoms.iter().enumerate() {
        let base = index * 3;
        body.push_str(&format!(
            "{symbol:<3} {:>18.10} {:>18.10} {:>18.10}\n",
            positions_angstrom[base],
            positions_angstrom[base + 1],
            positions_angstrom[base + 2]
        ));
    }
    body
}

/// Builds the input for any job the general panel offers.
///
/// The richer sibling of [`input_text`], which the two single-purpose panels
/// still use. This one knows about methods above the functional, about the
/// auxiliary bases they need, and about the blocks a saddle search or an
/// excited-state run adds underneath.
#[allow(clippy::too_many_arguments)]
pub fn panel_input_text(
    job: JobType,
    method: QcMethod,
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    excited: ExcitedStateOptions,
    extra_keywords: &str,
    xyz_body: &str,
) -> String {
    let multiplicity = multiplicity.max(1);
    let nproc = config.nproc.max(1);

    let mut text = String::new();
    text.push_str(&format!("%pal nprocs {nproc} end\n"));
    // `%maxcore` is memory per core where the panel asks for a total.
    let per_core = (config.memory_mb / nproc).max(256);
    text.push_str(&format!("%maxcore {per_core}\n"));

    text.push('!');
    text.push(' ');
    text.push_str(method_keyword(method, &config.functional));
    for keyword in basis_keywords(method, &config.functional, &config.basis) {
        text.push(' ');
        text.push_str(&keyword);
    }
    // Only a density functional takes one: a wavefunction method has no
    // dispersion term to add, and a composite brings its own.
    if method.needs_functional() {
        if let Some(dispersion) =
            orca_method::dispersion_keyword(&config.functional, config.orca_dispersion)
        {
            text.push(' ');
            text.push_str(dispersion);
        }
    }
    let job_words = job_keywords(job);
    if !job_words.is_empty() {
        text.push(' ');
        text.push_str(job_words);
    }
    let extra = extra_keywords.trim();
    if !extra.is_empty() {
        text.push(' ');
        text.push_str(extra);
    }
    text.push('\n');

    text.push_str(&job_blocks(job, excited));

    text.push_str(&format!("* xyz {charge} {multiplicity}\n"));
    text.push_str(xyz_body);
    if !xyz_body.ends_with('\n') {
        text.push('\n');
    }
    text.push_str("*\n");
    text
}

/// One self-consistent field iteration, for the run summary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScfIteration {
    pub iteration: u32,
    pub energy_hartree: f64,
    /// Change from the previous iteration, in Hartree.
    pub delta_e: f64,
}

/// Reads ORCA's self-consistent field history out of its output.
///
/// The rows look like this, under a banner naming the algorithm:
///
/// ```text
///                ***  Starting incremental Fock matrix formation  ***
///   0    -75.3062011796   0.000000000000 0.03317880  0.00427796  0.1149037 0.7000
///   1    -75.3213834576  -0.015182277969 0.02182200  0.00272513  0.0674985 0.7000
/// ```
///
/// Rows are recognised by shape -- an iteration number then at least two
/// numbers, the first of which is a large negative energy -- rather than by
/// counting lines after the banner, so surrounding text can change without this
/// silently reading the wrong rows. That is the same rule Behemoth's
/// convergence-table reader follows.
pub fn parse_scf_iterations(output_text: &str) -> Vec<ScfIteration> {
    let mut iterations = Vec::new();
    for line in output_text.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            continue;
        }
        let Ok(iteration) = fields[0].parse::<u32>() else {
            continue;
        };
        let (Ok(energy), Ok(delta)) = (fields[1].parse::<f64>(), fields[2].parse::<f64>()) else {
            continue;
        };
        // A total electronic energy is large and negative; this is what keeps
        // ordinary numbered output out of the table.
        if energy > -0.5 || !energy.is_finite() || !delta.is_finite() {
            continue;
        }
        iterations.push(ScfIteration { iteration, energy_hartree: energy, delta_e: delta });
    }
    iterations
}

/// The final single-point energy ORCA printed, in Hartree.
///
/// Taken as the last one, because a job that optimised prints one per cycle and
/// the last is the answer.
pub fn parse_final_energy(output_text: &str) -> Option<f64> {
    output_text
        .lines()
        .filter(|line| line.contains("FINAL SINGLE POINT ENERGY"))
        .next_back()
        .and_then(|line| {
            line.split_whitespace()
                .filter_map(|token| token.parse::<f64>().ok())
                .next_back()
        })
}

/// The energy components worth showing under a finished run.
///
/// Each is `None` when the job did not produce it: a Hartree-Fock single point
/// has no correlation energy, and only a job that optimised has a final
/// gradient worth quoting.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Energetics {
    /// ORCA's own `FINAL SINGLE POINT ENERGY`.
    ///
    /// Read [`Energetics::final_includes_excitation`] before showing this as
    /// "the energy": after an excited-state run it is not the ground state.
    pub final_energy_hartree: Option<f64>,
    /// The converged self-consistent field total energy: the ground state.
    pub scf_energy_hartree: Option<f64>,
    pub correlation_energy_hartree: Option<f64>,
    pub dispersion_correction_hartree: Option<f64>,
    /// Set when ORCA's final figure is an excited state rather than the ground
    /// state.
    ///
    /// After a TD-DFT run ORCA adds the first root's excitation energy to the
    /// self-consistent field energy and prints the sum as
    /// `FINAL SINGLE POINT ENERGY`. Verified on a real run: the output states
    /// `E(SCF) = -76.321074728`, `DE(CIS) = 0.279708561 (Root 1)` and
    /// `E(tot) = -76.041366167`, and it is that last figure ORCA calls final.
    ///
    /// Showing it as the molecule's energy would be wrong by the whole
    /// excitation energy, which is why this flag exists rather than the panel
    /// simply printing whatever came last.
    pub final_includes_excitation: bool,
}

impl Energetics {
    /// Whether there is anything at all to show.
    pub fn has_content(&self) -> bool {
        self.final_energy_hartree.is_some()
            || self.scf_energy_hartree.is_some()
            || self.correlation_energy_hartree.is_some()
            || self.dispersion_correction_hartree.is_some()
    }

    /// The ground-state energy, whatever the job was.
    ///
    /// The self-consistent field energy when the final figure carries an
    /// excitation, and the final figure otherwise -- which is the correlated
    /// total for a method that has one.
    pub fn ground_state_hartree(&self) -> Option<f64> {
        if self.final_includes_excitation {
            self.scf_energy_hartree
        } else {
            self.final_energy_hartree.or(self.scf_energy_hartree)
        }
    }
}

/// Reads the energy breakdown out of an ORCA output.
pub fn parse_energetics(output_text: &str) -> Energetics {
    /// The **first** number after `needle`, on the last line carrying it.
    ///
    /// First, not last, and that is the whole point: ORCA prints a total energy
    /// in Hartree and then repeats it in electronvolts on the same line, as in
    /// `Total Energy       :   -76.32107472787662 Eh   -2076.80203 eV`. Taking
    /// the last number reads the electronvolt figure and reports it as Hartree,
    /// which is wrong by a factor of twenty-seven.
    fn first_number_after(text: &str, needle: &str) -> Option<f64> {
        text.lines()
            .filter(|line| line.contains(needle))
            .next_back()
            .and_then(|line| {
                let tail = line.split(needle).nth(1)?;
                tail.split_whitespace()
                    .find_map(|token| token.parse::<f64>().ok())
            })
    }

    Energetics {
        final_energy_hartree: parse_final_energy(output_text),
        scf_energy_hartree: first_number_after(output_text, "Total Energy       :"),
        correlation_energy_hartree: first_number_after(output_text, "E(CORR)")
            .or_else(|| first_number_after(output_text, "Total correlation energy")),
        dispersion_correction_hartree: first_number_after(output_text, "Dispersion correction"),
        // ORCA's own summary of an excited-state run names the root it added.
        // Either spelling appears depending on the response method used.
        final_includes_excitation: output_text.contains("DE(CIS)")
            || output_text.contains("DE(TD-DFT)"),
    }
}

/// The Molden file ORCA's converter writes, given this run's basename.
pub const MOLDEN_FILE: &str = "orca_job.molden.input";

/// Turns a finished run's binary orbital file into Molden format.
///
/// ORCA writes no Molden file itself. Its orbitals live in a binary `.gbw`, and
/// `orca_2mkl <basename> -molden` is the converter that reads one -- a separate
/// program shipped beside ORCA, which is why this is a second process rather
/// than a keyword on the first. The Smite interface does exactly this.
///
/// `binary` is the ORCA executable; the converter sits beside it under a name
/// derived from its own, so a user who pointed at ORCA has already pointed at
/// this.
pub fn convert_to_molden(binary: &Path, run_dir: &Path) -> Result<PathBuf, String> {
    let converter = {
        let name = binary
            .file_name()
            .map(|n| format!("{}_2mkl", n.to_string_lossy()))
            .unwrap_or_else(|| "orca_2mkl".to_string());
        binary.with_file_name(name)
    };
    if !converter.is_file() {
        return Err(format!(
            "ORCA's Molden converter was not found at {}. It ships beside the ORCA              executable.",
            converter.display()
        ));
    }

    let output = Command::new(&converter)
        .current_dir(run_dir)
        .arg(JOB_BASENAME)
        .arg("-molden")
        .output()
        .map_err(|e| format!("Could not run ORCA's Molden converter: {e}"))?;

    let written = run_dir.join(MOLDEN_FILE);
    if !output.status.success() || !written.is_file() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!(
            "ORCA's Molden converter wrote nothing. {}",
            stderr.lines().next_back().unwrap_or("").trim()
        ));
    }
    Ok(written)
}

/// The command that runs a prepared input.
pub fn run_command(binary: &Path, run_dir: &Path) -> Command {
    let mut command = Command::new(binary);
    // ORCA is given the input's name relative to its working directory, and
    // writes every one of its files beside it.
    command
        .current_dir(run_dir)
        .arg(format!("{JOB_BASENAME}.inp"));
    command
}

/// Whether a finished run actually succeeded.
///
/// The exit status is necessary but not sufficient: ORCA can exit zero having
/// failed, printing its error into the output instead. Both are checked, and
/// a run that never reached its normal-termination banner is a failure too --
/// which is what a job killed by the operating system looks like.
pub fn check_output(output_text: &str, exit_ok: bool) -> Result<(), String> {
    if output_text.contains(ERROR_TERMINATION) {
        let detail = output_text
            .lines()
            .rev()
            .take(40)
            .find(|line| {
                let trimmed = line.trim();
                !trimmed.is_empty() && !trimmed.starts_with("ORCA finished by error")
            })
            .unwrap_or("see the run directory")
            .trim()
            .to_string();
        return Err(format!("ORCA stopped with an error: {detail}"));
    }
    if !output_text.contains(NORMAL_TERMINATION) {
        return Err(if exit_ok {
            "ORCA stopped before finishing and printed no error. See orca_job.out in the run \
             directory."
                .to_string()
        } else {
            "ORCA failed. See orca_job.out in the run directory.".to_string()
        });
    }
    Ok(())
}

/// Whether a geometry optimisation reached ORCA's convergence criteria.
///
/// An unconverged final geometry is still worth showing -- it is where the
/// optimiser got to -- so this is reported alongside the result rather than
/// turned into a failure.
pub fn converged(output_text: &str) -> bool {
    output_text.contains(CONVERGED)
}

/// Rewrites ORCA's optimisation trajectory into the form Beavyr reads.
///
/// ORCA writes each frame's comment as
/// `Coordinates from ORCA-job orca_job E -74.963023139028`, confirmed on a
/// real run. Beavyr's trajectory player, energy plot and summary window all
/// read `cycle <n> E = <value> Eh`, which is Behemoth's format and the one
/// `dreiding_run` also writes into. Translating here means every consumer
/// downstream stays untouched, and there is one trajectory format in Beavyr
/// rather than one per backend.
///
/// A frame whose comment carries no energy keeps its place in the trajectory
/// and is given no energy, so a partial file from a cancelled run still plays.
pub fn trajectory_to_beavyr(orca_trajectory: &str) -> String {
    let mut out = String::new();
    let mut lines = orca_trajectory.lines().peekable();
    let mut cycle = 0usize;

    while let Some(count_line) = lines.next() {
        let trimmed = count_line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let Ok(natoms) = trimmed.parse::<usize>() else {
            continue;
        };
        let comment = lines.next().unwrap_or("");

        // The frame is assembled before any of it is kept, so a file cut short
        // by a cancelled run contributes its whole frames and nothing else. A
        // frame written header-first would leave a final one claiming more
        // atoms than it carries, which no reader could use.
        let mut frame = String::new();
        let mut complete = true;
        for _ in 0..natoms {
            match lines.next() {
                Some(atom) => {
                    frame.push_str(atom.trim_end());
                    frame.push('\n');
                }
                None => {
                    complete = false;
                    break;
                }
            }
        }
        if !complete {
            return out;
        }

        cycle += 1;
        out.push_str(&format!("{natoms}\n"));
        match energy_from_comment(comment) {
            Some(energy) => out.push_str(&format!("cycle {cycle} E = {energy:.10} Eh\n")),
            None => out.push_str(&format!("cycle {cycle}\n")),
        }
        out.push_str(&frame);
    }
    out
}

/// The energy out of an ORCA trajectory comment line.
///
/// Taken as the last token that parses as a number, so the surrounding wording
/// can change without this breaking. The job name sits between the words and
/// the number and never parses, which is what makes "last number" safe here.
fn energy_from_comment(comment: &str) -> Option<f64> {
    comment
        .split_whitespace()
        .filter_map(|token| token.parse::<f64>().ok())
        .next_back()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn method() -> MethodConfig {
        MethodConfig {
            functional: "B3LYP".to_string(),
            basis: "def2-SVP".to_string(),
            memory_mb: 4096,
            nproc: 4,
            ..Default::default()
        }
    }

    fn water_body() -> String {
        xyz_body(
            &["O".to_string(), "H".to_string(), "H".to_string()],
            &[0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692],
        )
    }

    #[test]
    fn the_input_carries_the_method_the_basis_and_the_job() {
        let text = input_text(JobType::OptimizeThenFrequencies, 0, 1, &method(), "", &water_body());
        assert!(text.contains("! B3LYP def2-SVP Opt Freq"), "{text}");
        assert!(text.contains("* xyz 0 1"), "{text}");
        assert!(text.trim_end().ends_with('*'), "{text}");
    }

    /// The panel asks for a total; ORCA's `%maxcore` is per core. Passing the
    /// total unchanged would let a four-thread job take four times what the
    /// user asked for.
    #[test]
    fn memory_is_divided_between_the_cores() {
        let text = input_text(JobType::Optimize, 0, 1, &method(), "", &water_body());
        assert!(text.contains("%maxcore 1024"), "4096 MB over 4 cores: {text}");
        assert!(text.contains("%pal nprocs 4 end"), "{text}");
    }

    /// A tiny total must not become an unusable per-core share.
    #[test]
    fn a_small_memory_total_keeps_a_workable_floor() {
        let mut small = method();
        small.memory_mb = 64;
        small.nproc = 8;
        let text = input_text(JobType::Optimize, 0, 1, &small, "", &water_body());
        assert!(text.contains("%maxcore 256"), "{text}");
    }

    /// A composite defines its own basis, so none is written beside it.
    #[test]
    fn a_composite_carries_no_basis() {
        let mut composite = method();
        composite.functional = "r2SCAN-3c".to_string();
        let text = input_text(JobType::Optimize, 0, 1, &composite, "", &water_body());
        assert!(text.contains("! r2SCAN-3c Opt"), "{text}");
        assert!(!text.contains("def2-SVP"), "{text}");
    }

    /// The catalogue's spelling is what reaches ORCA, not the label. Sending
    /// the readable `M06-2X` would be an input error.
    #[test]
    fn the_functional_is_written_in_orcas_own_spelling() {
        let mut minnesota = method();
        minnesota.functional = "M062X".to_string();
        let text = input_text(JobType::Gradient, 0, 1, &minnesota, "", &water_body());
        assert!(text.contains("! M062X def2-SVP EnGrad"), "{text}");
        assert!(!text.contains("M06-2X"), "{text}");
    }

    /// Extra keywords are appended to the `!` line, so anything the dropdowns
    /// do not cover can still be asked for.
    #[test]
    fn extra_keywords_reach_the_keyword_line() {
        let text = input_text(JobType::Optimize, -1, 3, &method(), "TightSCF D4", &water_body());
        assert!(text.contains("! B3LYP def2-SVP Opt TightSCF D4"), "{text}");
        assert!(text.contains("* xyz -1 3"), "{text}");
    }

    /// ORCA can exit zero having failed. Trusting the status alone would treat
    /// an error message as a finished calculation.
    #[test]
    fn an_error_termination_is_a_failure_even_on_a_zero_exit() {
        let output = concat!(
            "   Some ordinary output\n",
            "   SCF NOT CONVERGED AFTER 125 CYCLES\n",
            "ORCA finished by error termination in SCF\n"
        );
        let err = check_output(output, true).expect_err("this run failed");
        assert!(err.contains("ORCA stopped with an error"), "{err}");
        assert!(err.contains("NOT CONVERGED"), "the reason is carried: {err}");
    }

    /// Output that never reached the normal-termination banner is a failure
    /// too, which is what a killed job looks like.
    #[test]
    fn output_without_a_normal_termination_is_a_failure() {
        assert!(check_output("started, then nothing\n", true).is_err());
        assert!(check_output("", false).is_err());
    }

    #[test]
    fn a_normal_termination_is_a_success() {
        let output = "FINAL SINGLE POINT ENERGY -74.9\n\n****ORCA TERMINATED NORMALLY****\n";
        assert!(check_output(output, true).is_ok());
        assert!(!converged(output), "no optimisation was asked for");
        let optimised = format!("{CONVERGED}\n{output}");
        assert!(converged(&optimised));
    }

    /// Real ORCA trajectory text, verbatim from an `Opt Freq` run on water, is
    /// rewritten into the one format every consumer in Beavyr already reads.
    #[test]
    fn a_real_trajectory_is_rewritten_into_beavyrs_format() {
        let orca = concat!(
            "3\n",
            "Coordinates from ORCA-job orca_job E -74.963023139028\n",
            "  O       0.000000      0.000000      0.117300\n",
            "  H       0.000000      0.757200     -0.469200\n",
            "  H       0.000000     -0.757200     -0.469200\n",
            "3\n",
            "Coordinates from ORCA-job orca_job E -74.965725493874\n",
            "  O       0.000000      0.000000      0.120000\n",
            "  H       0.000000      0.760000     -0.470000\n",
            "  H       0.000000     -0.760000     -0.470000\n"
        );
        let beavyr = trajectory_to_beavyr(orca);
        assert!(beavyr.contains("cycle 1 E = -74.9630231390 Eh"), "{beavyr}");
        assert!(beavyr.contains("cycle 2 E = -74.9657254939 Eh"), "{beavyr}");
        // The same energies Beavyr's own trajectory reader takes out.
        let energies = super::super::behemoth::parse_trajectory_energies(&beavyr);
        assert_eq!(energies.len(), 2);
        assert!((energies[0] + 74.963_023_139_0).abs() < 1.0e-9);
        assert_eq!(beavyr.lines().filter(|l| l.trim() == "3").count(), 2);
    }

    /// The converted trajectory has to satisfy the parser the *panel* uses to
    /// build frames, not merely the one that reads energies out of it.
    ///
    /// This gap is why an ORCA optimisation could finish and leave the
    /// trajectory player empty: the energy reader was tested, the frame reader
    /// was not, and only the second one loads anything on screen.
    #[test]
    fn the_converted_trajectory_loads_as_frames() {
        let orca = concat!(
            "3\n",
            "Coordinates from ORCA-job orca_job E -74.963023139028\n",
            "  O       0.000000      0.000000      0.117300\n",
            "  H       0.000000      0.757200     -0.469200\n",
            "  H       0.000000     -0.757200     -0.469200\n",
            "3\n",
            "Coordinates from ORCA-job orca_job E -74.965725493874\n",
            "  O       0.000000      0.000000      0.120000\n",
            "  H       0.000000      0.760000     -0.470000\n",
            "  H       0.000000     -0.760000     -0.470000\n"
        );
        let beavyr = trajectory_to_beavyr(orca);
        let frames = crate::trajectory::parse_multi_xyz(&beavyr)
            .unwrap_or_else(|e| panic!("the player must be able to load this: {e}\n{beavyr}"));
        assert_eq!(frames.len(), 2, "{beavyr}");
        assert_eq!(frames[0].atoms.len(), 3);
    }

    /// A trajectory cut short by a cancelled run still yields the frames that
    /// were complete, rather than nothing at all.
    #[test]
    fn a_truncated_trajectory_keeps_its_whole_frames() {
        let orca = concat!(
            "3\n",
            "Coordinates from ORCA-job orca_job E -74.963023139028\n",
            "  O       0.000000      0.000000      0.117300\n",
            "  H       0.000000      0.757200     -0.469200\n",
            "  H       0.000000     -0.757200     -0.469200\n",
            "3\n",
            "Coordinates from ORCA-job orca_job E -74.965725493874\n",
            "  O       0.000000      0.000000      0.120000\n"
        );
        let beavyr = trajectory_to_beavyr(orca);
        assert_eq!(
            super::super::behemoth::parse_trajectory_energies(&beavyr).len(),
            1,
            "only the complete frame survives: {beavyr}"
        );
    }

    /// Every file of a run is named from one basename, and none of them
    /// collides with the `input.xyz` the other backends write in the same
    /// kind of directory.
    #[test]
    fn the_job_files_do_not_collide_with_other_backends() {
        let files = job_files(Path::new("/tmp/run"));
        assert!(files.trajectory.ends_with("orca_job_trj.xyz"));
        assert!(files.final_geometry.ends_with("orca_job.xyz"));
        assert_ne!(files.final_geometry, Path::new("/tmp/run/input.xyz"));
        assert!(files.hessian.ends_with("orca_job.hess"));
    }

    /// Writes a real input, runs the real ORCA on it, and reads back what it
    /// wrote. Ignored by default because it needs ORCA installed and takes a
    /// few seconds; run it with
    /// `cargo test -- --ignored orca_really_runs`, or point it somewhere else
    /// with `BEAVYR_ORCA_PATH`.
    ///
    /// The other tests here assert on text. This is the one that would notice
    /// if ORCA changed a file name or a banner.
    #[test]
    #[ignore = "needs ORCA installed"]
    fn a_generated_input_really_runs() {
        use std::fs;

        let binary = std::env::var("BEAVYR_ORCA_PATH").unwrap_or_else(|_| {
            "/home/peter/Programs/orca_6_1_1_linux_x86-64_shared_openmpi418/orca".to_string()
        });
        if !Path::new(&binary).is_file() {
            eprintln!("no ORCA at {binary}; skipping");
            return;
        }

        let dir = std::env::temp_dir().join(format!("beavyr_orca_e2e_{}", std::process::id()));
        fs::create_dir_all(&dir).expect("scratch directory");

        let mut cheap = method();
        cheap.functional = "HF".to_string();
        cheap.basis = "STO-3G".to_string();
        cheap.nproc = 1;
        let files = job_files(&dir);
        fs::write(
            &files.input,
            input_text(JobType::OptimizeThenFrequencies, 0, 1, &cheap, "", &water_body()),
        )
        .expect("write the input");

        // Redirected exactly as the panels do it, into ORCA's own `.out`.
        // Capturing standard output and writing that file by hand instead is
        // what hid a real bug: the app redirected somewhere else and then read
        // an empty `.out`, so every ORCA run failed while this test passed.
        let handle = fs::File::create(&files.output).expect("create the output file");
        let status = run_command(Path::new(&binary), &dir)
            .stdout(std::process::Stdio::from(handle))
            .status()
            .expect("run ORCA");
        let printed = fs::read_to_string(&files.output).expect("ORCA wrote its output");
        check_output(&printed, status.success()).expect("the job should succeed");
        assert!(converged(&printed), "a water optimisation should converge");

        // The trajectory, translated into the one format Beavyr reads.
        let orca_trajectory = fs::read_to_string(&files.trajectory).expect("a trajectory");
        let beavyr = trajectory_to_beavyr(&orca_trajectory);
        let energies = super::super::behemoth::parse_trajectory_energies(&beavyr);
        assert!(energies.len() >= 2, "expected several steps: {energies:?}");
        assert!(
            energies.last().unwrap() <= energies.first().unwrap(),
            "an optimisation should not end higher than it began: {energies:?}"
        );

        // The Hessian, through the reader Beavyr already had.
        let hess_text = fs::read_to_string(&files.hessian).expect("a .hess");
        let hess = super::super::orca_hess::parse_orca_hess(&hess_text)
            .expect("the existing ORCA reader takes it");
        assert_eq!(hess.natoms, 3);
        assert_eq!(hess.hessian.dim(), (9, 9));
        assert_eq!(hess.masses_amu.len(), 3);
        assert!(
            !hess.ir_spectrum.is_empty(),
            "a Freq job reports an infrared spectrum"
        );

        let _ = fs::remove_dir_all(&dir);
    }

    /// A double hybrid, run for real, to prove the auxiliary basis is actually
    /// what makes it work. Without the `/C` set ORCA rejects the input
    /// outright, so this fails loudly if the wiring is ever removed.
    ///
    /// Ignored by default; run with `cargo test -- --ignored double_hybrid`.
    #[test]
    #[ignore = "needs ORCA installed"]
    fn a_double_hybrid_really_runs_with_its_auxiliary_basis() {
        use std::fs;

        let binary = std::env::var("BEAVYR_ORCA_PATH").unwrap_or_else(|_| {
            "/home/peter/Programs/orca_6_1_1_linux_x86-64_shared_openmpi418/orca".to_string()
        });
        if !Path::new(&binary).is_file() {
            eprintln!("no ORCA at {binary}; skipping");
            return;
        }

        let mut double = method();
        double.functional = "B2PLYP".to_string();
        double.basis = "def2-SVP".to_string();
        double.nproc = 1;

        let text = input_text(JobType::Gradient, 0, 1, &double, "", &water_body());
        assert!(
            text.contains("! B2PLYP def2-SVP def2-SVP/C EnGrad"),
            "the auxiliary basis must be on the keyword line: {text}"
        );

        let dir = std::env::temp_dir().join(format!("beavyr_orca_dh_{}", std::process::id()));
        fs::create_dir_all(&dir).expect("scratch directory");
        let files = job_files(&dir);
        fs::write(&files.input, &text).expect("write the input");

        let handle = fs::File::create(&files.output).expect("create the output file");
        let status = run_command(Path::new(&binary), &dir)
            .stdout(std::process::Stdio::from(handle))
            .status()
            .expect("run ORCA");
        let printed = fs::read_to_string(&files.output).expect("ORCA wrote its output");
        check_output(&printed, status.success())
            .expect("a double hybrid with its auxiliary basis should run");

        // And the gradient it was asked for, through the reader Beavyr had.
        let engrad = fs::read_to_string(&files.engrad).expect("an .engrad");
        let gradient = super::super::orca_engrad::parse_orca_engrad(&engrad)
            .expect("the existing reader takes it");
        assert_eq!(gradient.gradient_bohr.len(), 9);
        assert!(gradient.energy_hartree < 0.0);

        let _ = fs::remove_dir_all(&dir);
    }

    /// Real ORCA output from a B3LYP/def2-SVP TD-DFT run on water, captured
    /// from the installed 6.1.1.
    const TDDFT_OUTPUT: &str = include_str!("../../tests/fixtures/orca_water_tddft.out");

    /// The trap this test exists for: after an excited-state run ORCA's
    /// `FINAL SINGLE POINT ENERGY` is the *first excited state*, being the
    /// self-consistent field energy plus that root's excitation energy. Showing
    /// it as the molecule's energy would be wrong by the whole excitation.
    #[test]
    fn a_tddft_final_energy_is_not_the_ground_state() {
        let energetics = parse_energetics(TDDFT_OUTPUT);
        assert!(
            energetics.final_includes_excitation,
            "the run added a root to its final figure"
        );

        let scf = energetics.scf_energy_hartree.expect("an SCF energy");
        let final_energy = energetics.final_energy_hartree.expect("a final energy");
        assert!((scf + 76.321_074_727_876_62).abs() < 1.0e-9, "{scf}");
        assert!((final_energy + 76.041_366_166_992).abs() < 1.0e-9, "{final_energy}");

        // The difference is the first root's excitation energy, as the output
        // itself states.
        let excitation = final_energy - scf;
        assert!(
            (excitation - 0.279_708_561).abs() < 1.0e-6,
            "the gap should be root 1: {excitation}"
        );

        // And the ground state is the SCF energy, not the final figure.
        let ground = energetics.ground_state_hartree().expect("a ground state");
        assert!((ground - scf).abs() < 1.0e-12, "{ground} vs {scf}");
    }

    /// ORCA repeats a total energy in electronvolts on the same line, so the
    /// Hartree figure is the first number after the label, never the last.
    /// Reading the last one is wrong by a factor of twenty-seven and looks
    /// entirely plausible.
    #[test]
    fn an_energy_is_read_in_hartree_not_electronvolts() {
        let line = "Total Energy       :        -76.32107472787662 Eh           -2076.80203 eV\n";
        let energetics = parse_energetics(line);
        let scf = energetics.scf_energy_hartree.expect("an SCF energy");
        assert!((scf + 76.321_074_727_876_62).abs() < 1.0e-9, "got {scf}");
        assert!(scf > -100.0, "that is the electronvolt figure: {scf}");
    }

    /// For an ordinary job the final figure *is* the answer, including the
    /// correlated total for a method that has one.
    #[test]
    fn an_ordinary_final_energy_is_the_ground_state() {
        let output = concat!(
            "Total Energy       :        -76.32107472787662 Eh\n",
            "FINAL SINGLE POINT ENERGY       -76.500000000000\n",
            "****ORCA TERMINATED NORMALLY****\n"
        );
        let energetics = parse_energetics(output);
        assert!(!energetics.final_includes_excitation);
        let ground = energetics.ground_state_hartree().expect("a ground state");
        assert!((ground + 76.5).abs() < 1.0e-9, "{ground}");
    }

    /// The self-consistent field history, read from the real output.
    #[test]
    fn the_scf_history_is_read_from_real_output() {
        let iterations = parse_scf_iterations(TDDFT_OUTPUT);
        assert!(!iterations.is_empty(), "no iterations found");
        for step in &iterations {
            assert!(
                step.energy_hartree < -1.0,
                "an electronic energy should be large and negative: {step:?}"
            );
        }
        // Iteration numbers rise.
        for pair in iterations.windows(2) {
            assert!(pair[1].iteration >= pair[0].iteration, "{pair:?}");
        }
    }

    /// Ordinary numbered output must not be mistaken for an iteration row.
    /// Without the energy test, any line starting with a number would land in
    /// the table.
    #[test]
    fn ordinary_numbered_lines_are_not_iterations() {
        let noise = concat!(
            "  1   2   3\n",
            "  1    0.5    0.25\n",
            "Cartesian coordinates:\n",
            "  1   8   0.000000   0.000000   0.117300\n"
        );
        assert!(parse_scf_iterations(noise).is_empty(), "{noise}");
    }

    /// Excited states asked for by the block alone. `! TDDFT` is an ORCA input
    /// error, which is the obvious guess and the reason this is tested.
    #[test]
    fn a_tddft_job_is_requested_by_its_block_not_a_keyword() {
        let excited = ExcitedStateOptions { states: 12, tda: true };
        let text = panel_input_text(
            JobType::TdDft,
            QcMethod::Dft,
            0,
            1,
            &method(),
            excited,
            "",
            &water_body(),
        );
        assert!(!text.contains("TDDFT "), "TDDFT is not a keyword: {text}");
        assert!(text.contains("%tddft"), "{text}");
        assert!(text.contains("nroots 12"), "{text}");
        assert!(text.contains("tda true"), "{text}");
        assert!(text.contains("! B3LYP def2-SVP"), "{text}");
    }

    /// Turning the Tamm-Dancoff approximation off gives the full response.
    #[test]
    fn the_tamm_dancoff_switch_reaches_the_block() {
        let text = panel_input_text(
            JobType::TdDft,
            QcMethod::Dft,
            0,
            1,
            &method(),
            ExcitedStateOptions { states: 5, tda: false },
            "",
            &water_body(),
        );
        assert!(text.contains("tda false"), "{text}");
    }

    /// A saddle search needs a real Hessian to start from and again as it goes,
    /// or it converges to whatever saddle happens to be nearest.
    #[test]
    fn a_transition_state_search_asks_for_hessians() {
        let text = panel_input_text(
            JobType::TransitionState,
            QcMethod::Dft,
            0,
            1,
            &method(),
            ExcitedStateOptions::default(),
            "",
            &water_body(),
        );
        assert!(text.contains(" OptTS"), "{text}");
        assert!(text.contains("%geom"), "{text}");
        assert!(text.contains("Calc_Hess true"), "{text}");
        assert!(text.contains("Recalc_Hess 5"), "{text}");
    }

    /// Each job puts the right words on the keyword line, and a single point
    /// puts none: ORCA computes an energy by default.
    #[test]
    fn each_job_asks_for_itself() {
        assert_eq!(job_keywords(JobType::SinglePoint), "");
        assert_eq!(job_keywords(JobType::TdDft), "");
        assert_eq!(job_keywords(JobType::Optimize), "Opt");
        assert_eq!(job_keywords(JobType::Frequencies), "Freq");
        assert_eq!(job_keywords(JobType::OptimizeThenFrequencies), "Opt Freq");
        assert_eq!(job_keywords(JobType::TransitionState), "OptTS");
        assert_eq!(job_keywords(JobType::Irc), "IRC");
    }

    /// A correlated method gets its fitting basis, or ORCA refuses the input.
    #[test]
    fn a_correlated_method_gets_its_fitting_basis() {
        for method in [QcMethod::Mp2, QcMethod::RiMp2, QcMethod::CcsdT, QcMethod::DlpnoCcsdT] {
            let keywords = basis_keywords(method, "B3LYP", "def2-TZVP");
            assert_eq!(
                keywords,
                vec!["def2-TZVP".to_string(), "def2-TZVP/C".to_string()],
                "{method:?}"
            );
        }
    }

    /// An F12 method needs more: its own orbital basis, a complementary
    /// auxiliary basis named after it, and a correlation-fitting set.
    #[test]
    fn an_f12_method_gets_its_cabs_as_well() {
        let keywords = basis_keywords(QcMethod::Mp2F12, "B3LYP", "cc-pVDZ-F12");
        assert_eq!(
            keywords,
            vec![
                "cc-pVDZ-F12".to_string(),
                "cc-pVDZ-F12-CABS".to_string(),
                "cc-pVTZ/C".to_string()
            ]
        );
        assert!(QcMethod::Mp2F12.is_f12());
    }

    /// An uncorrelated method gets the basis and nothing else.
    #[test]
    fn an_uncorrelated_method_gets_only_its_basis() {
        assert_eq!(
            basis_keywords(QcMethod::Hf, "B3LYP", "def2-SVP"),
            vec!["def2-SVP".to_string()]
        );
        assert_eq!(
            basis_keywords(QcMethod::Dft, "B3LYP", "def2-SVP"),
            vec!["def2-SVP".to_string()]
        );
        // A composite brings its own, so nothing at all.
        assert!(basis_keywords(QcMethod::Dft, "r2SCAN-3c", "def2-SVP").is_empty());
        // And a semi-empirical method takes none.
        assert!(basis_keywords(QcMethod::Gfn2, "B3LYP", "def2-SVP").is_empty());
    }

    /// A DFT double hybrid still gets its auxiliary set through this path.
    #[test]
    fn a_double_hybrid_keeps_its_auxiliary_basis_here_too() {
        let keywords = basis_keywords(QcMethod::Dft, "B2PLYP", "def2-SVP");
        assert_eq!(
            keywords,
            vec!["def2-SVP".to_string(), "def2-SVP/C".to_string()]
        );
    }

    /// The method keyword is the method's, not the functional's, except for
    /// DFT where the functional is the method.
    #[test]
    fn the_method_keyword_is_the_right_one() {
        assert_eq!(method_keyword(QcMethod::Hf, "B3LYP"), "HF");
        assert_eq!(method_keyword(QcMethod::RiMp2, "B3LYP"), "RI-MP2");
        assert_eq!(method_keyword(QcMethod::DlpnoCcsdT, "B3LYP"), "DLPNO-CCSD(T)");
        assert_eq!(method_keyword(QcMethod::Mp2F12, "B3LYP"), "MP2-F12");
        // DFT travels as its functional, in ORCA's own spelling.
        assert_eq!(method_keyword(QcMethod::Dft, "M062X"), "M062X");
        assert_eq!(method_keyword(QcMethod::Dft, "r2SCAN-3c"), "r2SCAN-3c");
    }

    /// A whole correlated input, end to end.
    #[test]
    fn a_correlated_job_builds_a_complete_input() {
        let mut config = method();
        config.basis = "def2-TZVP".to_string();
        let text = panel_input_text(
            JobType::SinglePoint,
            QcMethod::DlpnoCcsdT,
            -1,
            2,
            &config,
            ExcitedStateOptions::default(),
            "TightSCF",
            &water_body(),
        );
        assert!(
            text.contains("! DLPNO-CCSD(T) def2-TZVP def2-TZVP/C TightSCF"),
            "{text}"
        );
        assert!(text.contains("* xyz -1 2"), "{text}");
        // A single point adds no job keyword and no block.
        assert!(!text.contains("%geom"), "{text}");
        assert!(!text.contains("%tddft"), "{text}");
    }

    /// The converter call, exactly. ORCA writes no Molden file itself; this
    /// separate program, shipped beside the ORCA executable, reads the binary
    /// `.gbw` and writes `<job>.molden.input` next to it.
    ///
    /// Pinned because all three parts matter and none is guessable: the
    /// converter's name is derived from ORCA's own, the argument is the job
    /// basename with no extension, and the flag is `-molden` with one dash.
    #[test]
    fn the_molden_converter_is_called_the_one_way_that_works() {
        let orca = std::env::temp_dir().join("beavyr_fake_orca_dir").join("orca");
        let expected = orca.with_file_name("orca_2mkl");

        // The converter is looked for beside ORCA, under ORCA's own name plus
        // the suffix. A user who pointed at ORCA has already pointed at this.
        let err = convert_to_molden(&orca, Path::new("/tmp")).expect_err("no such converter");
        assert!(err.contains(&expected.display().to_string()), "{err}");
        assert!(err.contains("beside the ORCA"), "{err}");

        // And the file it produces is named from the job basename.
        assert_eq!(MOLDEN_FILE, format!("{JOB_BASENAME}.molden.input"));
    }

    #[test]
    fn the_command_runs_the_input_in_its_own_directory() {
        let command = run_command(Path::new("/bin/orca"), Path::new("/tmp/run"));
        let args: Vec<String> = command
            .get_args()
            .map(|a| a.to_string_lossy().to_string())
            .collect();
        assert_eq!(args, vec!["orca_job.inp"]);
        assert_eq!(command.get_current_dir(), Some(Path::new("/tmp/run")));
    }
}
