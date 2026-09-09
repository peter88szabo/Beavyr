#![allow(dead_code)]

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_92;

#[derive(Clone, Debug)]
pub struct QcInput {
    pub path: PathBuf,
    pub functional: String,
    pub charge: i32,
    pub multiplicity: i32,
    pub additional: String,
    pub wfu: bool,
}

/// Reads the total energy, in Hartree, from xTB's standard output.
///
/// xTB prints it twice and boxes both, so the number is never the last thing on its line:
///
/// ```text
///          :: total energy              -0.327368812180 Eh    ::
///           | TOTAL ENERGY               -0.327368812180 Eh   |
/// ```
///
/// Taking the last whitespace token therefore yields `::` or `|`. The number is the first token
/// on the line that parses as one, which holds for both layouts and survives whatever decoration
/// a future version puts around them.
pub fn parse_xtb_energy(output: &str) -> Result<f64, String> {
    for line in output.lines() {
        if !line.to_ascii_uppercase().contains("TOTAL ENERGY") {
            continue;
        }
        if let Some(energy) = line
            .split_whitespace()
            .filter_map(|token| token.parse::<f64>().ok())
            .next()
        {
            return Ok(energy);
        }
    }
    Err("Total energy not found in XTB output".to_string())
}

pub fn parse_xtb_grad(gradfile: &Path, natom: usize) -> Result<Vec<f64>, String> {
    let data = fs::read_to_string(gradfile)
        .map_err(|e| format!("Failed to read gradient file: {e}"))?;
    let lines: Vec<&str> = data.lines().collect();
    if lines.len() < 2 * natom + 2 {
        return Err("Gradient file is shorter than expected".to_string());
    }

    let mut grad = Vec::with_capacity(3 * natom);
    for line in lines.iter().skip(natom + 2).take(natom) {
        let mut it = line.split_whitespace();
        let x = it
            .next()
            .ok_or("Missing gradient x component")?
            .parse::<f64>()
            .map_err(|e| format!("Failed to parse gradient x: {e}"))?;
        let y = it
            .next()
            .ok_or("Missing gradient y component")?
            .parse::<f64>()
            .map_err(|e| format!("Failed to parse gradient y: {e}"))?;
        let z = it
            .next()
            .ok_or("Missing gradient z component")?
            .parse::<f64>()
            .map_err(|e| format!("Failed to parse gradient z: {e}"))?;
        grad.extend_from_slice(&[x, y, z]);
    }

    Ok(grad)
}

pub fn print_structure(atoms: &[String], q_bohr: &[f64], filename: &Path) -> Result<(), String> {
    if q_bohr.len() != atoms.len() * 3 {
        return Err("Coordinate length does not match atom count".to_string());
    }
    let mut file = fs::File::create(filename)
        .map_err(|e| format!("Failed to create XYZ file: {e}"))?;
    writeln!(file, "{}", atoms.len()).map_err(|e| e.to_string())?;
    writeln!(
        file,
        "Temporary structure for XTB energy/gradient"
    )
    .map_err(|e| e.to_string())?;

    for (i, atom) in atoms.iter().enumerate() {
        let jx = 3 * i;
        let jy = 3 * i + 1;
        let jz = 3 * i + 2;
        writeln!(
            file,
            "{:>3} {:>15.5} {:>15.5} {:>15.5}",
            atom,
            q_bohr[jx] * BOHR_TO_ANGSTROM,
            q_bohr[jy] * BOHR_TO_ANGSTROM,
            q_bohr[jz] * BOHR_TO_ANGSTROM
        )
        .map_err(|e| e.to_string())?;
    }

    Ok(())
}

fn build_command_args(
    inputfile: &str,
    charge: i32,
    multiplicity: i32,
    arg: &str,
    additional: &str,
    restart: bool,
) -> Vec<String> {
    let mut args = Vec::new();
    args.push(inputfile.to_string());
    args.push("-c".to_string());
    args.push(charge.to_string());
    if multiplicity > 1 {
        args.push("-u".to_string());
        args.push((multiplicity - 1).to_string());
    }
    if !arg.is_empty() {
        args.push(arg.to_string());
    }
    if !additional.trim().is_empty() {
        args.extend(additional.split_whitespace().map(|s| s.to_string()));
    }
    if restart {
        args.extend([
            "--etemp",
            "1000.0",
            "--gfnff",
            "--acc",
            "200",
        ]
        .iter()
        .map(|s| s.to_string()));
    }
    args
}

/// Runs xTB once, in `workdir`.
///
/// The directory is a parameter rather than the process's own, and has to be: xTB is given a fixed
/// input file name and writes its gradient to a fixed name beside it, so two calls sharing a
/// directory overwrite each other's files. Passing one in means a caller can give every
/// calculation its own -- which is what makes it safe to drive xTB in a loop, or from several
/// threads at once, as re-ranking a set of conformers does.
pub fn call_xtb_in(
    workdir: &Path,
    q_bohr: &[f64],
    atoms: &[String],
    qcinput: &QcInput,
    arg: &str,
) -> Result<(f64, Vec<f64>), String> {
    if qcinput.path.as_os_str().is_empty() {
        return Err("XTB path is empty".to_string());
    }
    fs::create_dir_all(workdir)
        .map_err(|e| format!("Failed to create {}: {e}", workdir.display()))?;

    // Named relative to `workdir`, which is also the process's working directory, so xTB reads
    // and writes inside it.
    let inputfile = "xtb_geom_file_for_abinitioMD.xyz";
    let gradfile = workdir.join("gradient");

    print_structure(atoms, q_bohr, &workdir.join(inputfile))?;

    let base_args = build_command_args(
        inputfile,
        qcinput.charge,
        qcinput.multiplicity,
        arg,
        &qcinput.additional,
        false,
    );
    let output = Command::new(&qcinput.path)
        .current_dir(workdir)
        .args(&base_args)
        .output()
        .map_err(|e| format!("Failed to run XTB: {e}"))?;

    let mut stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let mut ok = output.status.success();
    if !ok {
        let restart_args = build_command_args(
            inputfile,
            qcinput.charge,
            qcinput.multiplicity,
            arg,
            &qcinput.additional,
            true,
        );
        let restart = Command::new(&qcinput.path)
            .current_dir(workdir)
            .args(&restart_args)
            .output()
            .map_err(|e| format!("Failed to run XTB restart: {e}"))?;
        stdout = String::from_utf8_lossy(&restart.stdout).to_string();
        ok = restart.status.success();
    }

    if !ok {
        return Err("XTB command failed".to_string());
    }

    let energy = parse_xtb_energy(&stdout)?;
    let grad = if arg == "--grad" {
        parse_xtb_grad(&gradfile, atoms.len())?
    } else {
        vec![0.0; atoms.len() * 3]
    };

    Ok((energy, grad))
}

/// Runs xTB in the process's own working directory.
///
/// Kept for the imported code that expects it. Prefer [`call_xtb_in`]: this one cannot be used
/// twice concurrently, and leaves its scratch files wherever Beavyr happens to be running.
pub fn call_xtb(
    q_bohr: &[f64],
    atoms: &[String],
    qcinput: &QcInput,
    arg: &str,
) -> Result<(f64, Vec<f64>), String> {
    call_xtb_in(Path::new("."), q_bohr, atoms, qcinput, arg)
}

pub fn xtb_force(q_bohr: &[f64], atoms: &[String], qcinput: &QcInput) -> Result<Vec<f64>, String> {
    let (_energy, grad) = call_xtb(q_bohr, atoms, qcinput, "--grad")?;
    Ok(grad.into_iter().map(|g| -g).collect())
}

pub fn xtb_energy(
    file_wf: &Path,
    q_bohr: &[f64],
    atoms: &[String],
    qcinput: &QcInput,
) -> Result<f64, String> {
    let arg = if qcinput.wfu { "--molden" } else { "--dipole" };
    let (energy, _grad) = call_xtb(q_bohr, atoms, qcinput, arg)?;

    if qcinput.wfu {
        fs::copy("molden.input", file_wf)
            .map_err(|e| format!("Failed to copy molden input: {e}"))?;
    }

    Ok(energy)
}

pub fn xtb_hessian(
    q_bohr: &mut [f64],
    atoms: &[String],
    qcinput: &QcInput,
) -> Result<Vec<Vec<f64>>, String> {
    let dx = 0.002_f64;
    let ndim = q_bohr.len();
    let mut hess = vec![vec![0.0_f64; ndim]; ndim];

    for i in 0..ndim {
        q_bohr[i] += dx;
        let grad_p = xtb_force(q_bohr, atoms, qcinput)?;

        q_bohr[i] -= 2.0 * dx;
        let grad_m = xtb_force(q_bohr, atoms, qcinput)?;

        for j in 0..ndim {
            hess[i][j] = 0.5 * (grad_p[j] - grad_m[j]) / dx;
        }

        q_bohr[i] += dx;
    }

    Ok(hess)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Both layouts xTB prints, taken verbatim from its output. The number is never last on the
    /// line, so a parser that takes the last token reads `::` or `|` and fails.
    #[test]
    fn the_total_energy_is_read_from_either_layout() {
        let boxed = concat!(
            "         :::::::::::::::::::::::::::::::::::::::::::::::::::::\n",
            "         :: total energy              -0.327368812180 Eh    ::\n",
            "         :: gradient norm              0.017706027413 Eh/a0 ::\n"
        );
        assert_eq!(parse_xtb_energy(boxed).unwrap(), -0.327_368_812_180);

        let summary = concat!(
            "           -------------------------------------------------\n",
            "          | TOTAL ENERGY               -0.327368812180 Eh   |\n",
            "          | GRADIENT NORM               0.017706027413 Eh/a |\n"
        );
        assert_eq!(parse_xtb_energy(summary).unwrap(), -0.327_368_812_180);
    }

    /// Output with no energy must be an error, not a zero -- which would quietly rank a conformer
    /// far below every other and look like a discovery.
    #[test]
    fn output_without_an_energy_is_an_error() {
        assert!(parse_xtb_energy("normal termination of xtb\n").is_err());
        assert!(parse_xtb_energy("| TOTAL ENERGY  not available |\n").is_err());
        assert!(parse_xtb_energy("").is_err());
    }
}
