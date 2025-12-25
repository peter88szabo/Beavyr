#![allow(dead_code)]

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_92;

#[derive(Clone, Debug)]
pub struct QcInput {
    pub path: PathBuf, // python executable path (empty => python3)
    pub nproc: usize,
    pub functional: String,
    pub basis: String,
    pub charge: i32,
    pub multiplicity: i32,
    pub additional: String,
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
    writeln!(f, "nproc={}", qcinput.nproc).map_err(|e| e.to_string())?;
    writeln!(f, "additional={}", qcinput.additional).map_err(|e| e.to_string())?;
    writeln!(f, "wfu={}", if qcinput.wfu { 1 } else { 0 }).map_err(|e| e.to_string())?;
    Ok(())
}

fn run_psi4_python(
    mode: &str,
    xyz: &str,
    qcinput: &QcInput,
    wfn_path: Option<&Path>,
) -> Result<String, String> {
    let workdir = Path::new("psi4_tmp");
    fs::create_dir_all(workdir).map_err(|e| format!("Failed to create psi4_tmp: {e}"))?;

    let xyz_path = workdir.join("input.xyz");
    fs::write(&xyz_path, xyz).map_err(|e| format!("Failed to write XYZ: {e}"))?;

    let params_path = workdir.join("params.txt");
    write_params(&params_path, qcinput)?;

    let out_path = workdir.join("result.txt");
    let wfn_base = wfn_path
        .map(|p| p.to_string_lossy().to_string())
        .unwrap_or_else(|| String::new());

    let script = r#"
import sys
import re
import numpy as np
import psi4

params_path, xyz_path, mode, out_path, wfn_base = sys.argv[1:6]

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
nproc = int(params.get("nproc", "1"))
additional = params.get("additional", "")
wfu = params.get("wfu", "0") == "1"

with open(xyz_path, "r") as f:
    xyz = f.read()

psi4.core.set_output_file("/dev/null", False)
psi4.set_options({"PRINT": 0})
psi4.set_num_threads(nproc)

mol = psi4.geometry(f"""{charge} {multiplicity}\n{xyz}""")

fct = functional.lower()
is_mgga = any(tag in fct for tag in ["scan", "r2scan", "m06-l", "tpss", "b97m"])

add = additional.lower() if additional else ""
tight = ("tightscf" in add)
soscf_requested = ("soscf" in add)
soscf_allowed = (soscf_requested and not is_mgga)

maxiter = 200
m = re.search(r"maxiter\s*=\s*(\d+)", add)
if m:
    maxiter = int(m.group(1))

def _run_psi4(multiref=False):
    opts = {
        "basis": basis,
        "reference": "uks" if multiplicity > 1 else "rks",
        "scf_type": "df",
        "maxiter": maxiter,
        "guess": "read",
        "diis": True,
        "soscf": soscf_allowed,
        "e_convergence": 1e-7 if tight else 1e-6,
        "d_convergence": 1e-7 if tight else 1e-6,
    }
    if multiref:
        opts["scf_type"] = "pk"
        opts["maxiter"] = max(maxiter * 2, 400)
        opts["soscf"] = (not is_mgga)
        opts["guess"] = "sad"
        opts["diis_start"] = 5
        opts["damping_percentage"] = 50
        opts["level_shift"] = 0.5
    psi4.set_options(opts)
    if mode == "energy":
        return psi4.energy(functional, return_wfn=True)
    return psi4.gradient(functional)

try:
    if mode == "energy":
        ene, wfn = _run_psi4(multiref=False)
    else:
        grad = _run_psi4(multiref=False)
except Exception:
    if mode == "energy":
        ene, wfn = _run_psi4(multiref=True)
    else:
        grad = _run_psi4(multiref=True)

with open(out_path, "w") as f:
    if mode == "energy":
        f.write(str(ene) + "\n")
        if wfu and wfn_base:
            psi4.molden(wfn, wfn_base + ".molden")
    else:
        grad_matrix = np.array(grad)
        force = -grad_matrix.flatten(order="C")
        f.write(" ".join([str(x) for x in force]) + "\n")
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
        .arg(wfn_base)
        .status()
        .map_err(|e| format!("Failed to run psi4 python: {e}"))?;

    if !status.success() {
        return Err("Psi4 Python invocation failed".to_string());
    }

    let out = fs::read_to_string(&out_path)
        .map_err(|e| format!("Failed to read psi4 output: {e}"))?;
    Ok(out)
}

pub fn psi4_energy(
    filename_base: &Path,
    q_bohr: &[f64],
    atoms: &[String],
    qcinput: &QcInput,
) -> Result<f64, String> {
    let xyz = make_xyz(atoms, q_bohr)?;
    let out = run_psi4_python("energy", &xyz, qcinput, Some(filename_base))?;
    let val = out
        .split_whitespace()
        .next()
        .ok_or("Psi4 output missing energy")?
        .parse::<f64>()
        .map_err(|e| format!("Failed to parse energy: {e}"))?;
    Ok(val)
}

pub fn psi4_force(
    q_bohr: &[f64],
    atoms: &[String],
    qcinput: &QcInput,
) -> Result<Vec<f64>, String> {
    let xyz = make_xyz(atoms, q_bohr)?;
    let out = run_psi4_python("force", &xyz, qcinput, None)?;
    let mut force = Vec::new();
    for tok in out.split_whitespace() {
        force.push(tok.parse::<f64>().map_err(|e| format!("Bad force value: {e}"))?);
    }
    Ok(force)
}

pub fn psi4_hessian(
    q_bohr: &mut [f64],
    atoms: &[String],
    qcinput: &QcInput,
) -> Result<Vec<Vec<f64>>, String> {
    let dx = 0.002_f64;
    let ndim = q_bohr.len();
    let mut hess = vec![vec![0.0_f64; ndim]; ndim];

    for i in 0..ndim {
        q_bohr[i] += dx;
        let grad_p = psi4_force(q_bohr, atoms, qcinput)?;

        q_bohr[i] -= 2.0 * dx;
        let grad_m = psi4_force(q_bohr, atoms, qcinput)?;

        for j in 0..ndim {
            hess[i][j] = 0.5 * (grad_p[j] - grad_m[j]) / dx;
        }

        q_bohr[i] += dx;
    }

    Ok(hess)
}
