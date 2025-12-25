#![allow(dead_code)]

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_92;

#[derive(Clone, Debug)]
pub struct QcInput {
    pub path: PathBuf, // python executable path (empty => python3)
    pub functional: String,
    pub basis: String,
    pub charge: i32,
    pub multiplicity: i32,
    pub wfu: bool,
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

fn write_params(path: &Path, qcinput: &QcInput) -> Result<(), String> {
    let mut f = fs::File::create(path).map_err(|e| format!("Failed to write params: {e}"))?;
    writeln!(f, "functional={}", qcinput.functional).map_err(|e| e.to_string())?;
    writeln!(f, "basis={}", qcinput.basis).map_err(|e| e.to_string())?;
    writeln!(f, "charge={}", qcinput.charge).map_err(|e| e.to_string())?;
    writeln!(f, "multiplicity={}", qcinput.multiplicity).map_err(|e| e.to_string())?;
    writeln!(f, "wfu={}", if qcinput.wfu { 1 } else { 0 }).map_err(|e| e.to_string())?;
    Ok(())
}

fn run_pyscf_python(
    mode: &str,
    xyz: &str,
    qcinput: &QcInput,
    wfn_path: Option<&Path>,
) -> Result<String, String> {
    let workdir = Path::new("pyscf_tmp");
    fs::create_dir_all(workdir).map_err(|e| format!("Failed to create pyscf_tmp: {e}"))?;

    let xyz_path = workdir.join("input.xyz");
    fs::write(&xyz_path, xyz).map_err(|e| format!("Failed to write XYZ: {e}"))?;

    let params_path = workdir.join("params.txt");
    write_params(&params_path, qcinput)?;

    let out_path = workdir.join("result.txt");
    let wfn_file = wfn_path
        .map(|p| p.to_string_lossy().to_string())
        .unwrap_or_else(|| String::new());

    let script = r#"
import sys
import numpy as np
from pyscf import gto, dft
from pyscf.tools import molden

params_path, xyz_path, mode, out_path, wfn_file = sys.argv[1:6]

params = {}
with open(params_path, "r") as f:
    for line in f:
        if "=" in line:
            k, v = line.strip().split("=", 1)
            params[k.strip()] = v.strip()

functional = params.get("functional", "")
basis = params.get("basis", "")
charge = int(params.get("charge", "0"))
multiplicity = int(params.get("multiplicity", "1"))
wfu = params.get("wfu", "0") == "1"

with open(xyz_path, "r") as f:
    xyz = f.read()

mol = gto.M(
    atom=xyz,
    basis=basis,
    charge=charge,
    spin=multiplicity - 1,
    verbose=0,
)

if multiplicity != 1:
    mf = dft.UKS(mol).run(xc=functional)
else:
    mf = dft.RKS(mol).run(xc=functional)

with open(out_path, "w") as f:
    if mode == "energy":
        e = mf.kernel()
        f.write(str(e) + "\n")
        if wfu and wfn_file:
            with open(wfn_file, "w") as wf:
                molden.header(mol, wf)
                molden.orbital_coeff(mol, wf, mf.mo_coeff, ene=mf.mo_energy, occ=mf.mo_occ)
    elif mode == "force":
        g = mf.nuc_grad_method()
        grad = sum(g.kernel().tolist(), [])
        force = [-x for x in grad]
        f.write(" ".join([str(x) for x in force]) + "\n")
    elif mode == "hessian":
        h = mf.Hessian().kernel()
        hess = h.transpose(0, 2, 1, 3).reshape(mol.natm * 3, mol.natm * 3)
        flat = hess.flatten(order="C")
        f.write(str(hess.shape[0]) + "\n")
        f.write(" ".join([str(x) for x in flat]) + "\n")
"#;

    let python = if qcinput.path.as_os_str().is_empty() {
        PathBuf::from("python3")
    } else {
        qcinput.path.clone()
    };

    let status = Command::new(python)
        .arg("-c")
        .arg(script)
        .arg(params_path)
        .arg(xyz_path)
        .arg(mode)
        .arg(&out_path)
        .arg(wfn_file)
        .status()
        .map_err(|e| format!("Failed to run pyscf python: {e}"))?;

    if !status.success() {
        return Err("PySCF Python invocation failed".to_string());
    }

    let out = fs::read_to_string(&out_path)
        .map_err(|e| format!("Failed to read pyscf output: {e}"))?;
    Ok(out)
}

pub fn pyscf_force(
    q_bohr: &[f64],
    atoms: &[String],
    qcinput: &QcInput,
) -> Result<Vec<f64>, String> {
    let xyz = make_xyz(atoms, q_bohr)?;
    let out = run_pyscf_python("force", &xyz, qcinput, None)?;
    let mut force = Vec::new();
    for tok in out.split_whitespace() {
        force.push(tok.parse::<f64>().map_err(|e| format!("Bad force value: {e}"))?);
    }
    Ok(force)
}

pub fn pyscf_energy(
    filename: &Path,
    q_bohr: &[f64],
    atoms: &[String],
    qcinput: &QcInput,
) -> Result<f64, String> {
    let xyz = make_xyz(atoms, q_bohr)?;
    let out = run_pyscf_python("energy", &xyz, qcinput, Some(filename))?;
    let val = out
        .split_whitespace()
        .next()
        .ok_or("PySCF output missing energy")?
        .parse::<f64>()
        .map_err(|e| format!("Failed to parse energy: {e}"))?;
    Ok(val)
}

pub fn pyscf_hessian(
    natoms: usize,
    xyz_angstrom: &str,
    qcinput: &QcInput,
) -> Result<Vec<Vec<f64>>, String> {
    let out = run_pyscf_python("hessian", xyz_angstrom, qcinput, None)?;
    let mut it = out.split_whitespace();
    let n = it
        .next()
        .ok_or("PySCF output missing dimension")?
        .parse::<usize>()
        .map_err(|e| format!("Failed to parse dimension: {e}"))?;
    if n != natoms * 3 {
        return Err("Hessian dimension does not match natoms".to_string());
    }
    let mut vals = Vec::new();
    for tok in it {
        vals.push(tok.parse::<f64>().map_err(|e| format!("Bad hessian value: {e}"))?);
    }
    if vals.len() != n * n {
        return Err("Hessian size mismatch".to_string());
    }
    let mut hess = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            hess[i][j] = vals[i * n + j];
        }
    }
    Ok(hess)
}
