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

pub fn parse_xtb_energy(output: &str) -> Result<f64, String> {
    for line in output.lines() {
        if line.contains("TOTAL ENERGY") {
            if let Some(tok) = line.split_whitespace().last() {
                return tok
                    .parse::<f64>()
                    .map_err(|e| format!("Failed to parse energy token '{tok}': {e}"));
            }
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

pub fn call_xtb(
    q_bohr: &[f64],
    atoms: &[String],
    qcinput: &QcInput,
    arg: &str,
) -> Result<(f64, Vec<f64>), String> {
    if qcinput.path.as_os_str().is_empty() {
        return Err("XTB path is empty".to_string());
    }

    let inputfile = "xtb_geom_file_for_abinitioMD.xyz";
    let gradfile = Path::new("gradient");

    print_structure(atoms, q_bohr, Path::new(inputfile))?;

    let base_args = build_command_args(
        inputfile,
        qcinput.charge,
        qcinput.multiplicity,
        arg,
        &qcinput.additional,
        false,
    );
    let output = Command::new(&qcinput.path)
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
        parse_xtb_grad(gradfile, atoms.len())?
    } else {
        vec![0.0; atoms.len() * 3]
    };

    Ok((energy, grad))
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
