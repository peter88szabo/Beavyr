#!/usr/bin/env python3
"""Offline generator: curated CCD monomers -> xTB minima -> static Rust data.

Generation dependencies only: gemmi==0.7.5, rdkit==2025.9.6, numpy.
The application never imports Python or reads these source/audit files.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import urllib.request
from datetime import datetime, timezone

import gemmi
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/amino_acids"
RESIDUES = [
    ("GLY", "Glycine", "G"), ("ALA", "Alanine", "A"),
    ("VAL", "Valine", "V"), ("LEU", "Leucine", "L"),
    ("ILE", "Isoleucine", "I"), ("SER", "Serine", "S"),
    ("THR", "Threonine", "T"), ("CYS", "Cysteine", "C"),
    ("MET", "Methionine", "M"), ("PRO", "Proline", "P"),
    ("ASP", "AsparticAcid", "D"), ("GLU", "GlutamicAcid", "E"),
    ("ASN", "Asparagine", "N"), ("GLN", "Glutamine", "Q"),
    ("HIS", "Histidine", "H"), ("LYS", "Lysine", "K"),
    ("ARG", "Arginine", "R"), ("PHE", "Phenylalanine", "F"),
    ("TYR", "Tyrosine", "Y"), ("TRP", "Tryptophan", "W"),
    ("SEC", "Selenocysteine", "U"), ("PYL", "Pyrrolysine", "O"),
]
ACIDIC = {"ASP": "OD2", "GLU": "OE2", "CYS": "SG", "TYR": "OH", "SEC": "SE"}
BASIC = {"LYS": ("NZ", "HZ3"), "ARG": ("NH2", "HH22")}
METHOD_ARGS = ["--gfn", "2", "--alpb", "water", "--uhf", "0", "--ohess", "verytight", "--acc", "0.1", "--parallel", "4"]
PROTOCOL_VERSION = 1


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def source(code: str, *, offline: bool = False) -> dict:
    path = DATA / "ccd" / f"{code}.cif"
    url = f"https://files.rcsb.org/ligands/download/{code}.cif"
    if not path.exists():
        if offline:
            raise FileNotFoundError(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(urllib.request.urlopen(url, timeout=30).read())
    block = gemmi.cif.read_file(str(path)).sole_block()
    fields = ["atom_id", "type_symbol", "charge", "pdbx_stereo_config",
              "pdbx_model_Cartn_x_ideal", "pdbx_model_Cartn_y_ideal", "pdbx_model_Cartn_z_ideal"]
    atoms = [dict(name=r[0], element=r[1].capitalize(), charge=int(r[2]), stereo=r[3],
                  xyz=[float(x) for x in list(r)[4:7]]) for r in block.find("_chem_comp_atom.", fields)]
    bonds = [dict(a=r[0], b=r[1], order=r[2], aromatic=r[3] == "Y") for r in
             block.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order", "pdbx_aromatic_flag"])]
    return dict(code=code, atoms=atoms, bonds=bonds, source_url=url, source_sha256=digest(path.read_bytes()))


def remove_hydrogen(model: dict, name: str) -> None:
    atom = next(a for a in model["atoms"] if a["name"] == name)
    assert atom["element"] == "H", name
    model["atoms"].remove(atom)
    model["bonds"] = [b for b in model["bonds"] if name not in (b["a"], b["b"])]


def change_charge(model: dict, name: str, charge: int) -> None:
    next(a for a in model["atoms"] if a["name"] == name)["charge"] = charge


def neighbours(model: dict, name: str) -> list[str]:
    return [b["b"] if b["a"] == name else b["a"] for b in model["bonds"] if name in (b["a"], b["b"])]


def variant(base: dict, termini: str, side_chain: str) -> dict:
    model = deepcopy(base)
    code = model["code"]
    if code in BASIC and side_chain == "Neutral":
        heavy, hydrogen = BASIC[code]
        remove_hydrogen(model, hydrogen)
        change_charge(model, heavy, 0)
    if code == "HIS" and side_chain != "Protonated":
        remove_hydrogen(model, "HE2" if side_chain == "HistidineDelta" else "HD1")
        change_charge(model, "ND1", 0)
    if side_chain == "Deprotonated":
        heavy = ACIDIC[code]
        names = {a["name"]: a for a in model["atoms"]}
        hydrogens = [n for n in neighbours(model, heavy) if names[n]["element"] == "H"]
        assert len(hydrogens) == 1, (code, heavy, hydrogens)
        remove_hydrogen(model, hydrogens[0])
        change_charge(model, heavy, -1)
    if termini == "Zwitterionic":
        remove_hydrogen(model, "HXT")
        change_charge(model, "OXT", -1)
        change_charge(model, "N", 1)
        positions = {a["name"]: np.array(a["xyz"]) for a in model["atoms"]}
        n = positions["N"]
        vectors = [positions[name] - n for name in neighbours(model, "N")]
        direction = -sum(v / np.linalg.norm(v) for v in vectors)
        direction /= np.linalg.norm(direction)
        name = "H2" if code == "PRO" else "H3"
        assert name not in positions
        model["atoms"].append(dict(name=name, element="H", charge=0, stereo="N", xyz=(n + direction * 1.03).tolist()))
        model["bonds"].append(dict(a="N", b=name, order="SING", aromatic=False))
    model.update(termini=termini, side_chain=side_chain,
                 id=f"{code.lower()}_{termini.lower()}_{side_chain.lower()}")
    model["charge"] = sum(a["charge"] for a in model["atoms"])
    return model


def rdkit_molecule(model: dict, positions=None):
    rw = Chem.RWMol()
    indices = {}
    for item in model["atoms"]:
        atom = Chem.Atom(item["element"])
        atom.SetFormalCharge(item["charge"])
        atom.SetNoImplicit(True)
        indices[item["name"]] = rw.AddAtom(atom)
    kinds = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE}
    for b in model["bonds"]:
        rw.AddBond(indices[b["a"]], indices[b["b"]], Chem.BondType.AROMATIC if b["aromatic"] else kinds[b["order"]])
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    conf = Chem.Conformer(len(model["atoms"]))
    conf.Set3D(True)
    for i, atom in enumerate(model["atoms"]):
        conf.SetAtomPosition(i, atom["xyz"] if positions is None else positions[i])
    mol.AddConformer(conf)
    Chem.AssignAtomChiralTagsFromStructure(mol, replaceExistingTags=True)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    for i, atom in enumerate(model["atoms"]):
        if atom["stereo"] in ("R", "S"):
            actual = mol.GetAtomWithIdx(i).GetProp("_CIPCode") if mol.GetAtomWithIdx(i).HasProp("_CIPCode") else None
            if actual != atom["stereo"]:
                raise ValueError(f"{model['id']} {atom['name']}: expected {atom['stereo']}, found {actual}")
    return mol


def xyz_text(model: dict, positions, comment: str) -> str:
    return f"{len(model['atoms'])}\n{comment}\n" + "".join(
        f"{atom['element']:<2} {p[0]: .12f} {p[1]: .12f} {p[2]: .12f}\n" for atom, p in zip(model["atoms"], positions))


def check_geometry(model: dict, positions) -> None:
    rdkit_molecule(model, positions)
    xyz = np.array(positions)
    if xyz.shape != (len(model["atoms"]), 3) or not np.isfinite(xyz).all():
        raise ValueError("Invalid coordinate array")
    names = {a["name"]: i for i, a in enumerate(model["atoms"])}
    table = Chem.GetPeriodicTable()
    radii = [table.GetRcovalent(a["element"]) for a in model["atoms"]]
    bonded = set()
    for b in model["bonds"]:
        i, j = names[b["a"]], names[b["b"]]
        bonded.add(tuple(sorted((i, j))))
        distance = float(np.linalg.norm(xyz[i] - xyz[j]))
        limit = radii[i] + radii[j]
        if not 0.65 * limit < distance < 1.35 * limit:
            raise ValueError(f"Bond changed: {b['a']}-{b['b']} = {distance:.3f} A")
    heavy = [i for i, a in enumerate(model["atoms"]) if a["element"] != "H"]
    for i, atom in enumerate(model["atoms"]):
        if atom["element"] == "H":
            nearest = min(heavy, key=lambda j: np.linalg.norm(xyz[i] - xyz[j]))
            if tuple(sorted((i, nearest))) not in bonded:
                raise ValueError(f"Proton transfer: {atom['name']} is now closest to {model['atoms'][nearest]['name']}")
        for j in range(i):
            if (j, i) not in bonded and np.linalg.norm(xyz[i] - xyz[j]) < 0.95 * (radii[i] + radii[j]):
                raise ValueError(f"Unexpected close contact: {atom['name']}-{model['atoms'][j]['name']}")


def optimize(model: dict, binary: str, work: Path, attempts: int) -> dict:
    molecule = rdkit_molecule(model)
    failures = []
    for attempt in range(attempts):
        directory = work / model["id"] / str(attempt)
        directory.mkdir(parents=True, exist_ok=True)
        if attempt == 0:
            positions = [a["xyz"] for a in model["atoms"]]
        else:
            candidate = Chem.Mol(molecule)
            options = AllChem.ETKDGv3()
            options.randomSeed = 20261003 + attempt
            options.numThreads = 1
            options.useRandomCoords = attempt > 5
            if AllChem.EmbedMolecule(candidate, options) != 0:
                failures.append(f"{attempt}: conformer embedding failed")
                continue
            positions = candidate.GetConformer().GetPositions()
        input_text = xyz_text(model, positions, model["id"])
        (directory / "input.xyz").write_text(input_text)
        command = [binary, "input.xyz", *METHOD_ARGS, "--chrg", str(model["charge"])]
        env = dict(os.environ, OMP_NUM_THREADS="4", MKL_NUM_THREADS="4", OPENBLAS_NUM_THREADS="1", OMP_STACKSIZE="256M")
        with (directory / "xtb.log").open("w") as log:
            process = subprocess.run(command, cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=300)
        log_text = (directory / "xtb.log").read_text()
        try:
            if process.returncode != 0 or "GEOMETRY OPTIMIZATION CONVERGED" not in log_text or "normal termination of xtb" not in log_text:
                raise ValueError("xTB did not converge and terminate normally")
            lines = (directory / "xtbopt.xyz").read_text().splitlines()
            assert int(lines[0]) == len(model["atoms"])
            assert [line.split()[0] for line in lines[2:]] == [a["element"] for a in model["atoms"]]
            positions = [[float(x) for x in line.split()[1:4]] for line in lines[2:]]
            check_geometry(model, positions)
            frequencies, vibrations, rigid = [], [], []
            for line in (directory / "vibspectrum").read_text().splitlines():
                fields = line.split()
                if fields and fields[0].isdigit():
                    frequency = float(fields[2] if fields[1] == "a" else fields[1])
                    frequencies.append(frequency)
                    (vibrations if fields[1] == "a" else rigid).append(frequency)
            assert len(frequencies) == 3 * len(model["atoms"])
            assert len(rigid) == 6 and max(abs(v) for v in rigid) < 1.0
            assert len(vibrations) == 3 * len(model["atoms"]) - 6
            if min(vibrations) <= 0:
                raise ValueError(f"Not a verified minimum: lowest vibrational frequency {min(vibrations)} cm^-1")
            energy = float(re.search(r"energy:\s*([-\d.]+)", lines[1])[1])
            gradient = float(re.search(r"gnorm:\s*([-\d.]+)", lines[1])[1])
            if not np.isfinite([energy, gradient]).all() or gradient > 1e-4:
                raise ValueError(f"Excessive final gradient: {gradient}")
        except (ValueError, AssertionError, FileNotFoundError) as error:
            failures.append(f"{attempt}: {error}")
            print(f"  retry {model['id']}: {error}", flush=True)
            continue
        return dict(model, positions=positions, energy_hartree=energy, gradient_norm=gradient,
                    frequencies_cm1=frequencies, lowest_frequency_cm1=min(vibrations), attempt=attempt,
                    rejected_attempts=failures, input_xyz=input_text,
                    optimized_xyz=(directory / "xtbopt.xyz").read_text(),
                    log_sha256=digest(log_text.encode()), log_path=str(directory / "xtb.log"))
    raise RuntimeError(f"No valid minimum for {model['id']}: {failures}")


def rust_string(value) -> str:
    return json.dumps(value, ensure_ascii=True)


def formula(model) -> str:
    counts = Counter(a["element"] for a in model["atoms"])
    elements = [e for e in ["C", "H"] if e in counts] + sorted(set(counts) - {"C", "H"})
    return "".join(e + (str(counts[e]) if counts[e] > 1 else "") for e in elements)


def emit(results: list[dict], run_info: dict) -> None:
    """Publish only after every requested form has passed all checks."""
    lines = ["""//! Embedded free-amino-acid monomers, generated by scripts/generate_amino_acids.py.
//! GFN2-xTB / ALPB(water), verytight optimization, accuracy 0.1, closed-shell.
//! All entries are checked local minima, not claimed global conformer minima.
//! These are complete monomers: peptide assembly must handle terminal atoms and
//! formal charges explicitly. No runtime files, parsers, xTB, or allocations.
// Static templates shared by the peptide builder and its validation tests.
#![allow(dead_code)]

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AminoAcid {"""]
    lines += [f"    {name}," for _, name, _ in RESIDUES]
    lines.append("}\nimpl AminoAcid {")
    lines.append("    pub const ALL: [Self; 22] = [" + ", ".join("Self::" + name for _, name, _ in RESIDUES) + "];")
    for method, return_type, values in [
        ("name", "&'static str", [rust_string(re.sub(r"([a-z])([A-Z])", r"\1 \2", name)) for _, name, _ in RESIDUES]),
        ("three_letter", "&'static str", [rust_string(code) for code, _, _ in RESIDUES]),
        ("one_letter", "char", [repr(letter) for _, _, letter in RESIDUES]),
    ]:
        lines.append(f"    pub const fn {method}(self) -> {return_type} {{ match self {{")
        lines += [f"        Self::{name} => {value}," for (_, name, _), value in zip(RESIDUES, values)]
        lines.append("    }}")
    lines.append("""
    pub fn from_one_letter(code: char) -> Option<Self> {
        Self::ALL.into_iter().find(|a| a.one_letter() == code.to_ascii_uppercase())
    }
    pub fn from_three_letter(code: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|a| a.three_letter().eq_ignore_ascii_case(code))
    }
    pub fn templates(self) -> &'static [AminoAcidTemplate] { match self {""")
    offset = 0
    for code, name, _ in RESIDUES:
        count = sum(t["code"] == code for t in results)
        lines.append(f"        Self::{name} => &TEMPLATES[{offset}..{offset + count}],")
        offset += count
    lines.append("""    }}
    pub fn template(self, termini: TerminalState, side_chain: SideChainState) -> Option<&'static AminoAcidTemplate> {
        self.templates().iter().find(|t| t.termini == termini && t.side_chain == side_chain)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TerminalState {
    /// Uncharged amino and carboxylic-acid termini; proline has a secondary amine.
    Neutral,
    /// Terminal nitrogen +1 and terminal carboxylate -1. Side chains can add charge.
    Zwitterionic,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SideChainState {
    Neutral,
    Protonated,
    Deprotonated,
    /// Neutral histidine with the ring proton on ND1.
    HistidineDelta,
    /// Neutral histidine with the ring proton on NE2.
    HistidineEpsilon,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondOrder { Single, Double, Triple, Aromatic }
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Configuration { R, S }
#[derive(Debug, Clone, Copy)]
pub struct Atom {
    pub name: &'static str,
    pub element: &'static str,
    pub position_angstrom: [f64; 3],
    pub formal_charge: i8,
}
#[derive(Debug, Clone, Copy)]
pub struct Bond {
    pub atoms: [usize; 2],
    pub order: BondOrder,
}
#[derive(Debug, Clone, Copy)]
pub struct Backbone {
    pub nitrogen: usize,
    pub alpha_carbon: usize,
    pub carbonyl_carbon: usize,
    pub carbonyl_oxygen: usize,
    pub terminal_oxygen: usize,
}
#[derive(Debug, Clone, Copy)]
pub struct Stereocenter {
    pub atom: usize,
    pub configuration: Configuration,
    pub neighbours: [usize; 4],
    /// Sign of the oriented tetrahedron in the independent CCD seed geometry.
    pub reference_handedness: i8,
}
#[derive(Debug, Clone, Copy)]
pub struct Optimization {
    pub energy_hartree: f64,
    pub gradient_norm_hartree_per_bohr: f64,
    pub lowest_vibrational_wavenumber_cm1: f64,
    pub vibrational_mode_count: usize,
    pub ccd_sha256: &'static str,
    pub xtb_log_sha256: &'static str,
}
#[derive(Debug)]
pub struct AminoAcidTemplate {
    pub id: &'static str,
    pub amino_acid: AminoAcid,
    pub termini: TerminalState,
    pub side_chain: SideChainState,
    pub formula: &'static str,
    pub charge: i8,
    pub multiplicity: u8,
    pub atoms: &'static [Atom],
    pub bonds: &'static [Bond],
    pub backbone: Backbone,
    pub stereocenters: &'static [Stereocenter],
    pub optimization: Optimization,
}
impl AminoAcidTemplate {
    pub fn atom_index(&self, name: &str) -> Option<usize> {
        self.atoms.iter().position(|a| a.name == name)
    }
}
pub const METHOD: &str = "GFN2-xTB / ALPB(water); verytight; accuracy 0.1; unconstrained";
""")
    lines.append(f"pub const XTB_VERSION: &str = {rust_string(run_info['xtb_version'])};")
    lines.append(f"pub static TEMPLATES: [AminoAcidTemplate; {len(results)}] = [")
    residue_names = {code: name for code, name, _ in RESIDUES}
    for t in results:
        names = {a["name"]: i for i, a in enumerate(t["atoms"])}
        lines += ["    AminoAcidTemplate {", f"        id: {rust_string(t['id'])},",
                  f"        amino_acid: AminoAcid::{residue_names[t['code']]},",
                  f"        termini: TerminalState::{t['termini']},",
                  f"        side_chain: SideChainState::{t['side_chain']},",
                  f"        formula: {rust_string(formula(t))}, charge: {t['charge']}, multiplicity: 1,",
                  "        atoms: &["]
        for atom, p in zip(t["atoms"], t["positions"]):
            lines.append(f"            Atom {{ name: {rust_string(atom['name'])}, element: {rust_string(atom['element'])}, position_angstrom: [{p[0]:.12f}, {p[1]:.12f}, {p[2]:.12f}], formal_charge: {atom['charge']} }},")
        lines += ["        ],", "        bonds: &["]
        for b in t["bonds"]:
            kind = "Aromatic" if b["aromatic"] else {"SING": "Single", "DOUB": "Double", "TRIP": "Triple"}[b["order"]]
            lines.append(f"            Bond {{ atoms: [{names[b['a']]}, {names[b['b']]}], order: BondOrder::{kind} }},")
        lines += ["        ],", "        backbone: Backbone { " + ", ".join(f"{field}: {names[label]}" for field, label in [
            ("nitrogen", "N"), ("alpha_carbon", "CA"), ("carbonyl_carbon", "C"), ("carbonyl_oxygen", "O"), ("terminal_oxygen", "OXT")]) + " },", "        stereocenters: &["]
        for atom in t["atoms"]:
            if atom["stereo"] not in ("R", "S"): continue
            adjacent = sorted(names[n] for n in neighbours(t, atom["name"]))
            assert len(adjacent) == 4
            xyz = np.array([t["atoms"][i]["xyz"] for i in adjacent])
            sign = int(np.sign(np.linalg.det(xyz[:3] - xyz[3])))
            assert sign in (-1, 1)
            lines.append(f"            Stereocenter {{ atom: {names[atom['name']]}, configuration: Configuration::{atom['stereo']}, neighbours: {adjacent}, reference_handedness: {sign} }},")
        lines += ["        ],", "        optimization: Optimization {",
                  f"            energy_hartree: {t['energy_hartree']:.12f}, gradient_norm_hartree_per_bohr: {t['gradient_norm']:.12f},",
                  f"            lowest_vibrational_wavenumber_cm1: {t['lowest_frequency_cm1']:.2f}, vibrational_mode_count: {3*len(t['atoms'])-6},",
                  f"            ccd_sha256: {rust_string(t['source_sha256'])}, xtb_log_sha256: {rust_string(t['log_sha256'])},",
                  "        },", "    },"]
    lines += ["];", '#[cfg(test)]\n#[path = "amino_acid_tests.rs"]\nmod tests;']
    target = ROOT / "src/molecule_builder/amino_acids.rs"
    temporary = target.with_suffix(".rs.tmp")
    temporary.write_text("\n".join(lines) + "\n")
    subprocess.run(["rustfmt", "--edition", "2021", "--config", "skip_children=true", str(temporary)], check=True)
    temporary.replace(target)
    audit = []
    for t in results:
        entry = dict(t)
        log_path = Path(entry.pop("log_path"))
        log_directory = DATA / "logs"
        log_directory.mkdir(exist_ok=True)
        with (log_directory / (t["id"] + ".log.gz")).open("wb") as output:
            with gzip.GzipFile(filename="", fileobj=output, mode="wb", mtime=0) as compressed:
                compressed.write(log_path.read_bytes())
        geometries = DATA / "optimized"
        geometries.mkdir(exist_ok=True)
        (geometries / (t["id"] + ".xyz")).write_text(t["optimized_xyz"])
        entry["optimized_xyz_sha256"] = digest(t["optimized_xyz"].encode())
        spectra = DATA / "spectra"
        spectra.mkdir(exist_ok=True)
        spectrum = (log_path.parent / "vibspectrum").read_bytes()
        (spectra / (t["id"] + ".vibspectrum")).write_bytes(spectrum)
        entry["vibspectrum_sha256"] = digest(spectrum)
        audit.append(entry)
    manifest = dict(run_info, generated_utc=datetime.now(timezone.utc).isoformat(),
                    method_args=METHOD_ARGS, protocol_version=PROTOCOL_VERSION,
                    rdkit_version=rdBase.rdkitVersion, gemmi_version=gemmi.__version__,
                    numpy_version=np.__version__, templates=audit)
    (DATA / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Published {target.relative_to(ROOT)} and reviewable calculation records", flush=True)


def verify(models: list[dict]) -> None:
    """Audit the published records without network access, xTB, or file writes."""
    manifest = json.loads((DATA / "manifest.json").read_text())
    assert manifest["method_args"] == METHOD_ARGS
    assert manifest["protocol_version"] == PROTOCOL_VERSION
    entries = manifest["templates"]
    assert [t["id"] for t in entries] == [m["id"] for m in models]
    for model, entry in zip(models, entries):
        for key, value in model.items():
            assert entry[key] == value, (model["id"], key)
        identifier = entry["id"]
        log = gzip.decompress((DATA / "logs" / (identifier + ".log.gz")).read_bytes())
        assert digest(log) == entry["log_sha256"]
        assert b"GEOMETRY OPTIMIZATION CONVERGED" in log
        assert b"normal termination of xtb" in log
        xyz = (DATA / "optimized" / (identifier + ".xyz")).read_bytes()
        assert digest(xyz) == entry["optimized_xyz_sha256"]
        assert xyz.decode() == entry["optimized_xyz"]
        lines = xyz.decode().splitlines()
        assert int(lines[0]) == len(model["atoms"])
        assert [line.split()[0] for line in lines[2:]] == [a["element"] for a in model["atoms"]]
        positions = [[float(v) for v in line.split()[1:4]] for line in lines[2:]]
        assert positions == entry["positions"]
        check_geometry(model, positions)
        energy = float(re.search(r"energy:\s*([-\d.]+)", lines[1])[1])
        gradient = float(re.search(r"gnorm:\s*([-\d.]+)", lines[1])[1])
        assert np.isfinite([energy, gradient]).all() and 0 <= gradient <= 1e-4
        assert energy == entry["energy_hartree"] and gradient == entry["gradient_norm"]
        spectrum = (DATA / "spectra" / (identifier + ".vibspectrum")).read_bytes()
        assert digest(spectrum) == entry["vibspectrum_sha256"]
        frequencies, vibrations, rigid = [], [], []
        for line in spectrum.decode().splitlines():
            fields = line.split()
            if fields and fields[0].isdigit():
                frequency = float(fields[2] if fields[1] == "a" else fields[1])
                frequencies.append(frequency)
                (vibrations if fields[1] == "a" else rigid).append(frequency)
        assert frequencies == entry["frequencies_cm1"]
        assert len(rigid) == 6 and max(abs(v) for v in rigid) < 1.0
        assert len(vibrations) == 3 * len(model["atoms"]) - 6
        assert min(vibrations) == entry["lowest_frequency_cm1"] > 0
    print(f"Verified all {len(entries)} published templates, source graphs, geometries, logs, and spectra")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xtb")
    parser.add_argument("--work", type=Path, default=Path("/tmp/beavyr-amino-acids/runs"))
    parser.add_argument("--attempts", type=int, default=12)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    if not args.xtb and not (args.prepare_only or args.verify_only):
        parser.error("--xtb is required for optimization")
    if hasattr(os, "sched_getaffinity"):
        os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:4])
    models = []
    for code, _, _ in RESIDUES:
        base = source(code, offline=args.verify_only)
        states = ["HistidineEpsilon", "HistidineDelta", "Protonated"] if code == "HIS" else ["Neutral"]
        if code in BASIC: states.append("Protonated")
        if code in ACIDIC: states.append("Deprotonated")
        for termini in ["Neutral", "Zwitterionic"]:
            for side_chain in states:
                model = variant(base, termini, side_chain)
                rdkit_molecule(model)
                models.append(model)
    if args.verify_only:
        verify(models)
        return
    args.work.mkdir(parents=True, exist_ok=True)
    (args.work / "models.json").write_text(json.dumps(models, indent=2) + "\n")
    print(f"Prepared {len(models)} chemically checked templates across {len(RESIDUES)} amino acids", flush=True)
    if args.prepare_only: return
    version_output = subprocess.run([args.xtb, "--version"], capture_output=True, text=True, check=True)
    version = re.search(r"xtb version\s+([^\n]+)", version_output.stdout)[1].strip()
    binary_hash = digest(Path(args.xtb).read_bytes())
    run_info = dict(xtb_version=version, xtb_binary_sha256=binary_hash)
    results, errors = [], []
    for i, model in enumerate(models):
        cache = args.work / (model["id"] + ".json")
        try:
            fingerprint = digest(json.dumps([PROTOCOL_VERSION, model, METHOD_ARGS, run_info], sort_keys=True).encode())
            result = json.loads(cache.read_text()) if cache.exists() else None
            if result is None or result.get("fingerprint") != fingerprint:
                result = optimize(model, args.xtb, args.work, args.attempts)
                result["fingerprint"] = fingerprint
            check_geometry(model, result["positions"])
            cache.write_text(json.dumps(result, indent=2) + "\n")
            results.append(result)
            print(f"[{i+1}/{len(models)}] {model['id']}: q={model['charge']:+d}, E={result['energy_hartree']:.10f}, min frequency={result['lowest_frequency_cm1']:.2f}", flush=True)
        except Exception as error:
            errors.append(str(error))
            print(f"FAILED {model['id']}: {error}", flush=True)
    if errors:
        (args.work / "failures.json").write_text(json.dumps(errors, indent=2))
        raise RuntimeError(f"{len(errors)} templates failed; no Rust library published")
    (args.work / "validated.json").write_text(json.dumps(results, indent=2) + "\n")
    emit(results, run_info)


if __name__ == "__main__":
    main()
