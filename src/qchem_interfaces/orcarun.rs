#![allow(dead_code)]

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_92;

#[derive(Clone, Debug)]
pub struct QcInput {
    pub path: PathBuf,
    pub nproc: usize,
    pub functional: String,
    pub basis: String,
    pub charge: i32,
    pub multiplicity: i32,
    pub additional: String,
    pub wfu: bool,
    pub what: String, // for optimizer ("TS" or "minimum")
}

pub fn parse_xyz(xyz: &str) -> Result<(usize, Vec<String>, Vec<f64>), String> {
    let mut atoms = Vec::new();
    let mut q = Vec::new();
    for line in xyz.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err("Invalid XYZ line".to_string());
        }
        atoms.push(parts[0].to_string());
        q.push(parts[1].parse::<f64>().map_err(|e| e.to_string())?);
        q.push(parts[2].parse::<f64>().map_err(|e| e.to_string())?);
        q.push(parts[3].parse::<f64>().map_err(|e| e.to_string())?);
    }
    Ok((atoms.len(), atoms, q))
}

pub fn make_xyz(atoms: &[String], q_bohr: &[f64]) -> Result<String, String> {
    if q_bohr.len() != atoms.len() * 3 {
        return Err("Coordinate length does not match atom count".to_string());
    }
    let mut xyz = String::new();
    for (i, atom) in atoms.iter().enumerate() {
        let jx = 3 * i;
        let jy = 3 * i + 1;
        let jz = 3 * i + 2;
        xyz.push_str(&format!(
            "{atom}      {:.8}      {:.8}      {:.8}\n",
            q_bohr[jx] * BOHR_TO_ANGSTROM,
            q_bohr[jy] * BOHR_TO_ANGSTROM,
            q_bohr[jz] * BOHR_TO_ANGSTROM
        ));
    }
    Ok(xyz)
}

fn generate_input(
    inputfile: &Path,
    nproc: usize,
    xyz: &str,
    method: &str,
    basis: &str,
    charge: i32,
    multiplicity: i32,
    additional: &str,
) -> Result<(), String> {
    let mut f = fs::File::create(inputfile.with_extension("inp"))
        .map_err(|e| format!("Failed to write ORCA input: {e}"))?;
    writeln!(f, "%pal nprocs {nproc} end").map_err(|e| e.to_string())?;
    writeln!(f, "! {method} {basis} engrad {additional}").map_err(|e| e.to_string())?;
    writeln!(f, "* xyz {charge} {multiplicity}").map_err(|e| e.to_string())?;
    write!(f, "{xyz}").map_err(|e| e.to_string())?;
    writeln!(f, "*").map_err(|e| e.to_string())?;
    Ok(())
}

fn generate_optimizer_input(
    what: &str,
    inputfile: &Path,
    nproc: usize,
    xyz: &str,
    method: &str,
    basis: &str,
    charge: i32,
    multiplicity: i32,
    additional: &str,
) -> Result<(), String> {
    let mut f = fs::File::create(inputfile.with_extension("inp"))
        .map_err(|e| format!("Failed to write ORCA input: {e}"))?;
    writeln!(f, "%pal nprocs {nproc} end").map_err(|e| e.to_string())?;
    writeln!(f, "! {method} {basis} opt freq {additional}").map_err(|e| e.to_string())?;
    if what == "TS" {
        writeln!(f, "%geom").map_err(|e| e.to_string())?;
        writeln!(f, "TS_search EF").map_err(|e| e.to_string())?;
        writeln!(f, "Calc_Hess true").map_err(|e| e.to_string())?;
        writeln!(f, "Recalc_Hess 5").map_err(|e| e.to_string())?;
        writeln!(f, "MaxIter 100").map_err(|e| e.to_string())?;
        writeln!(f, "ProjectTR true").map_err(|e| e.to_string())?;
        writeln!(f, "end").map_err(|e| e.to_string())?;
    }
    writeln!(f, "* xyz {charge} {multiplicity}").map_err(|e| e.to_string())?;
    write!(f, "{xyz}").map_err(|e| e.to_string())?;
    writeln!(f, "*").map_err(|e| e.to_string())?;
    Ok(())
}

fn call_orca(inputfile: &Path, orcapath: &Path) -> Result<(), String> {
    if orcapath.as_os_str().is_empty() {
        return Err("Cannot determine ORCA path".to_string());
    }
    let output = Command::new(orcapath)
        .arg(inputfile.with_extension("inp"))
        .output()
        .map_err(|e| format!("Failed to run ORCA: {e}"))?;
    if !output.status.success() {
        return Err("ORCA command failed".to_string());
    }
    fs::write(inputfile.with_extension("out"), output.stdout)
        .map_err(|e| format!("Failed to write ORCA output: {e}"))?;
    Ok(())
}

fn read_energy_and_grad(inputfile: &Path) -> Result<(f64, Vec<f64>), String> {
    let data = fs::read_to_string(inputfile.with_extension("engrad"))
        .map_err(|e| format!("Failed to read engrad: {e}"))?;
    let mut lines = data.lines();
    for _ in 0..3 {
        lines.next();
    }
    let natoms: usize = lines
        .next()
        .ok_or("Missing natoms")?
        .trim()
        .parse()
        .map_err(|e| e.to_string())?;
    for _ in 0..3 {
        lines.next();
    }
    let energy: f64 = lines
        .next()
        .ok_or("Missing energy")?
        .trim()
        .parse()
        .map_err(|e| e.to_string())?;
    for _ in 0..3 {
        lines.next();
    }
    let mut grad = Vec::with_capacity(natoms * 3);
    for _ in 0..natoms * 3 {
        let val: f64 = lines
            .next()
            .ok_or("Missing gradient")?
            .trim()
            .parse()
            .map_err(|e| e.to_string())?;
        grad.push(val);
    }
    Ok((energy, grad))
}

fn read_energy_and_optimized_geoms(inputfile: &Path) -> Result<(bool, f64, Vec<f64>), String> {
    let out = fs::read_to_string(inputfile.with_extension("out"))
        .map_err(|e| format!("Failed to read output: {e}"))?;
    let converged = out.contains("THE OPTIMIZATION HAS CONVERGED");

    let xyz_data = fs::read_to_string(inputfile.with_extension("xyz"))
        .map_err(|e| format!("Failed to read optimized xyz: {e}"))?;
    let mut lines = xyz_data.lines();
    let natoms: usize = lines
        .next()
        .ok_or("Missing natoms")?
        .trim()
        .parse()
        .map_err(|e| e.to_string())?;
    let comment = lines.next().unwrap_or("");
    let energy = comment
        .split_whitespace()
        .find_map(|t| t.parse::<f64>().ok())
        .unwrap_or(0.0);

    let mut coords = Vec::with_capacity(natoms * 3);
    for _ in 0..natoms {
        if let Some(line) = lines.next() {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                coords.push(parts[1].parse::<f64>().map_err(|e| e.to_string())?);
                coords.push(parts[2].parse::<f64>().map_err(|e| e.to_string())?);
                coords.push(parts[3].parse::<f64>().map_err(|e| e.to_string())?);
            }
        }
    }
    Ok((converged, energy, coords))
}

pub fn run_orca(
    natoms: usize,
    xyz: &str,
    path: &Path,
    nproc: usize,
    method: &str,
    basis: &str,
    charge: i32,
    multiplicity: i32,
    additional: &str,
) -> Result<(f64, Vec<f64>), String> {
    let fname = "orca_file_for_abinitioMD";
    let directory = Path::new("orca_tmp");
    fs::create_dir_all(directory).map_err(|e| format!("Failed to create orca_tmp: {e}"))?;
    let inputfile = directory.join(fname);

    generate_input(
        &inputfile,
        nproc,
        xyz,
        method,
        basis,
        charge,
        multiplicity,
        additional,
    )?;
    call_orca(&inputfile, path)?;
    read_energy_and_grad(&inputfile)
}

pub fn run_orca_optimizer(
    what: &str,
    natoms: usize,
    xyz: &str,
    path: &Path,
    nproc: usize,
    method: &str,
    basis: &str,
    charge: i32,
    multiplicity: i32,
    additional: &str,
) -> Result<(bool, f64, Vec<f64>), String> {
    let fname = "orca_file_for_geomopt";
    let directory = Path::new("orca_tmp");
    fs::create_dir_all(directory).map_err(|e| format!("Failed to create orca_tmp: {e}"))?;
    let inputfile = directory.join(fname);

    let _ = natoms;
    generate_optimizer_input(
        what,
        &inputfile,
        nproc,
        xyz,
        method,
        basis,
        charge,
        multiplicity,
        additional,
    )?;
    call_orca(&inputfile, path)?;
    read_energy_and_optimized_geoms(&inputfile)
}

fn read_orca_hessian(inputfile: &Path) -> Result<Vec<Vec<f64>>, String> {
    let data = fs::read_to_string(inputfile.with_extension("hess"))
        .map_err(|e| format!("Failed to read hessian: {e}"))?;
    let mut lines = data.lines();
    while let Some(line) = lines.next() {
        if line.trim() == "$hessian" {
            break;
        }
    }
    let n: usize = lines
        .next()
        .ok_or("Missing hessian dimension")?
        .trim()
        .parse()
        .map_err(|e| e.to_string())?;
    let nblocks = (n + 4) / 5;
    let mut hess = vec![vec![0.0_f64; n]; n];
    for b in 0..nblocks {
        lines.next();
        let ncols = std::cmp::min(5, n - b * 5);
        for i in 0..n {
            let line = lines.next().ok_or("Short hessian block")?;
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 1 + ncols {
                return Err("Short hessian line".to_string());
            }
            for c in 0..ncols {
                hess[i][b * 5 + c] = parts[1 + c].parse::<f64>().map_err(|e| e.to_string())?;
            }
        }
    }
    Ok(hess)
}

pub fn get_orca_hessian(
    natoms: usize,
    xyz: &str,
    path: &Path,
    nproc: usize,
    method: &str,
    basis: &str,
    charge: i32,
    multiplicity: i32,
    additional: &str,
) -> Result<Vec<Vec<f64>>, String> {
    let fname = "orca_file_for_abinitioMD";
    let directory = Path::new("orca_tmp");
    fs::create_dir_all(directory).map_err(|e| format!("Failed to create orca_tmp: {e}"))?;
    let inputfile = directory.join(fname);

    let _ = natoms;
    let mut f = fs::File::create(inputfile.with_extension("inp"))
        .map_err(|e| format!("Failed to write ORCA input: {e}"))?;
    writeln!(f, "%pal nprocs {nproc} end").map_err(|e| e.to_string())?;
    writeln!(f, "! {method} {basis} numfreq {additional}").map_err(|e| e.to_string())?;
    writeln!(f, "* xyz {charge} {multiplicity}").map_err(|e| e.to_string())?;
    write!(f, "{xyz}").map_err(|e| e.to_string())?;
    writeln!(f, "*").map_err(|e| e.to_string())?;

    call_orca(&inputfile, path)?;
    read_orca_hessian(&inputfile)
}

pub fn orca_hessian(natoms: usize, xyz: &str, qcinput: &QcInput) -> Result<Vec<Vec<f64>>, String> {
    get_orca_hessian(
        natoms,
        xyz,
        &qcinput.path,
        qcinput.nproc,
        &qcinput.functional,
        &qcinput.basis,
        qcinput.charge,
        qcinput.multiplicity,
        &qcinput.additional,
    )
}

pub fn orca_geom_opt(
    natoms: usize,
    xyz: &str,
    qcinput: &QcInput,
) -> Result<(bool, f64, Vec<f64>), String> {
    run_orca_optimizer(
        &qcinput.what,
        natoms,
        xyz,
        &qcinput.path,
        qcinput.nproc,
        &qcinput.functional,
        &qcinput.basis,
        qcinput.charge,
        qcinput.multiplicity,
        &qcinput.additional,
    )
}

pub fn orca_force(
    q_bohr: &[f64],
    atoms: &[String],
    qcinput: &QcInput,
) -> Result<Vec<f64>, String> {
    let xyz = make_xyz(atoms, q_bohr)?;
    let natoms = atoms.len();
    let (_ene, grad) = run_orca(
        natoms,
        &xyz,
        &qcinput.path,
        qcinput.nproc,
        &qcinput.functional,
        &qcinput.basis,
        qcinput.charge,
        qcinput.multiplicity,
        &qcinput.additional,
    )?;
    Ok(grad.into_iter().map(|g| -g).collect())
}

pub fn orca_energy(
    filename: &Path,
    q_bohr: &[f64],
    atoms: &[String],
    qcinput: &QcInput,
) -> Result<f64, String> {
    let xyz = make_xyz(atoms, q_bohr)?;
    let natoms = atoms.len();
    let directory = Path::new("orca_tmp");
    let inputfile = directory.join("orca_file_for_abinitioMD");
    let engrad_path = inputfile.with_extension("engrad");

    if engrad_path.exists() {
        let (ene, _grad) = read_energy_and_grad(&inputfile)?;
        return Ok(ene);
    }

    let (ene, _grad) = run_orca(
        natoms,
        &xyz,
        &qcinput.path,
        qcinput.nproc,
        &qcinput.functional,
        &qcinput.basis,
        qcinput.charge,
        qcinput.multiplicity,
        &qcinput.additional,
    )?;

    if qcinput.wfu {
        let gbwfile = inputfile.clone();
        let mut cmd = Command::new(format!(
            "{}_2mkl",
            qcinput.path.to_string_lossy()
        ));
        cmd.arg(gbwfile)
            .arg("-molden");
        let output = cmd.output().map_err(|e| format!("Failed to run orca_2mkl: {e}"))?;
        if !output.status.success() {
            return Err("orca_2mkl failed".to_string());
        }

        let moldenfile = inputfile.with_extension("molden.input");
        fs::copy(moldenfile, filename)
            .map_err(|e| format!("Failed to copy molden file: {e}"))?;
    }

    Ok(ene)
}
