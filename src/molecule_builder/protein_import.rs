//! PDB/mmCIF coordinates via pdbtbx; Beavyr owns residue chemistry and explicit links.
use super::{
    amino_acids::{AminoAcid, BondOrder, SideChainState, TerminalState},
    peptide::Residue,
    peptide_edit,
    topology::{AtomInfo, ResidueInfo, ResidueOrigin, Topology},
};
use crate::molecule::{covalent_radius_angstrom, Molecule};
use bevy::prelude::Vec3;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufReader, Cursor};

pub struct Imported {
    pub molecule: Molecule,
    pub warnings: Vec<String>,
}
pub fn is_protein_path(path: &std::path::Path) -> bool {
    path.extension().and_then(|s| s.to_str()).is_some_and(|s| {
        matches!(
            s.to_ascii_lowercase().as_str(),
            "pdb" | "ent" | "cif" | "mmcif"
        )
    })
}
fn field(s: &str, a: usize, b: usize) -> &str {
    s.get(a..b.min(s.len())).unwrap_or("").trim()
}
fn empty(s: &str) -> bool {
    matches!(s, "" | "." | "?")
}
fn identity(s: &str) -> bool {
    empty(s) || matches!(s, "1_555" | "1555")
}
fn canonical_name(s: &str) -> String {
    if s.len() > 1 && s.as_bytes()[0].is_ascii_digit() {
        format!("{}{}", &s[1..], &s[..1])
    } else {
        match s {
            "HN" => "H".into(),
            "OT1" => "O".into(),
            "OT2" => "OXT".into(),
            _ => s.into(),
        }
    }
}
fn residue_template(
    name: &str,
    names: &HashSet<String>,
) -> Option<&'static super::amino_acids::AminoAcidTemplate> {
    let (base, state) = match name {
        "HID" | "HSD" => ("HIS", Some(SideChainState::HistidineDelta)),
        "HIE" | "HSE" => ("HIS", Some(SideChainState::HistidineEpsilon)),
        "HIP" | "HSP" => ("HIS", Some(SideChainState::Protonated)),
        "ASH" => ("ASP", Some(SideChainState::Neutral)),
        "GLH" => ("GLU", Some(SideChainState::Neutral)),
        "LYN" => ("LYS", Some(SideChainState::Neutral)),
        "CYM" => ("CYS", Some(SideChainState::Deprotonated)),
        "CYX" => ("CYS", Some(SideChainState::Neutral)),
        x => (x, None),
    };
    let a = AminoAcid::from_three_letter(base)?;
    let state = state.unwrap_or_else(|| {
        if a == AminoAcid::Histidine {
            match (names.contains("HD1"), names.contains("HE2")) {
                (true, true) => SideChainState::Protonated,
                (true, false) => SideChainState::HistidineDelta,
                _ => SideChainState::HistidineEpsilon,
            }
        } else {
            Residue::new(a).side_chain
        }
    });
    a.template(TerminalState::Neutral, state)
}

pub fn parse(text: &str, cif: bool) -> Result<Imported, String> {
    if text.len() > 128 * 1024 * 1024 {
        return Err("Protein file exceeds the 128 MiB import limit.".into());
    }
    let mut options = pdbtbx::ReadOptions::default();
    options
        .set_format(if cif {
            pdbtbx::Format::Mmcif
        } else {
            pdbtbx::Format::Pdb
        })
        .set_level(pdbtbx::StrictnessLevel::Loose)
        .set_only_first_model(true)
        .set_capitalise_chains(false)
        .set_only_atomic_coords(true);
    // The parser uses assertions for some malformed required mmCIF values.
    // A bad file must not terminate the desktop application.
    let parsed =
        std::panic::catch_unwind(|| options.read_raw(BufReader::new(Cursor::new(text.as_bytes()))))
            .map_err(|_| "Malformed protein coordinate records.".to_string())?;
    let (pdb, errors) = parsed.map_err(|e| {
        e.iter()
            .take(3)
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join("\n")
    })?;
    let model = pdb
        .models()
        .next()
        .ok_or("No coordinate model in this file.")?;
    let mut mol = Molecule::empty();
    let mut t = Topology::default();
    let mut warnings = vec![];
    if !errors.is_empty() {
        warnings.push(format!(
            "Parser reported {} nonfatal issue(s).",
            errors.len()
        ));
    }
    let mut serials = HashMap::new();
    let mut alternate_count = 0;
    let mut next_chain = 0;
    let mut breaks = HashSet::new();
    if !cif {
        let mut last = None;
        let mut segment = 0;
        let mut seen = HashMap::new();
        for l in text.lines() {
            if l.starts_with("ENDMDL") {
                break;
            }
            if l.starts_with("ATOM  ") || l.starts_with("HETATM") {
                last = field(l, 6, 11).parse::<usize>().ok();
                let key = (field(l, 21, 22), field(l, 22, 26), field(l, 26, 27));
                if seen.insert(key, segment).is_some_and(|old| old != segment) {
                    return Err("Repeated chain/residue identifier across TER records; use distinct chain IDs before importing.".into());
                }
            } else if l.starts_with("TER") {
                if let Some(last) = last {
                    breaks.insert(last);
                }
                segment += 1;
            }
        }
    }
    for chain in model.chains() {
        next_chain += 1;
        let mut previous_terminated = false;
        for res in chain.residues() {
            if previous_terminated {
                next_chain += 1;
            }
            // Select one coherent alternate per residue, sharing blank-alt atoms.
            let mut choices: Vec<_> = res
                .conformers()
                .filter(|c| c.alternative_location().is_some())
                .collect();
            choices.sort_by(|a, b| {
                let occ = |c: &pdbtbx::Conformer| {
                    let n = c.atoms().count().max(1);
                    c.atoms().map(|a| a.occupancy()).sum::<f64>() / n as f64
                };
                occ(b)
                    .total_cmp(&occ(a))
                    .then_with(|| a.alternative_location().cmp(&b.alternative_location()))
            });
            let choice = choices.first().copied();
            if !choices.is_empty() {
                alternate_count += 1;
            }
            let name = choice
                .map(|c| c.name())
                .or_else(|| res.name())
                .unwrap_or("UNK");
            let alt = choice.and_then(|c| c.alternative_location()).unwrap_or("");
            let mut atoms: BTreeMap<String, &pdbtbx::Atom> = BTreeMap::new();
            let mut conformers: Vec<_> = res
                .conformers()
                .filter(|c| {
                    c.alternative_location().is_none()
                        || (c.name() == name && c.alternative_location() == Some(alt))
                })
                .collect();
            conformers.sort_by_key(|c| c.alternative_location().is_some());
            for c in conformers {
                for a in c.atoms() {
                    atoms.insert(canonical_name(a.name()), a);
                }
            }
            previous_terminated = res.atoms().any(|a| breaks.contains(&a.serial_number()));
            let rid = t.residues.len() as u64 + 1;
            let names = atoms.keys().cloned().collect();
            let template = residue_template(name, &names);
            t.residues.push(ResidueInfo {
                id: rid,
                chain: next_chain,
                template: template.map(|t| t.id).unwrap_or(name).into(),
                origin: Some(ResidueOrigin {
                    chain: chain.id().into(),
                    number: res.serial_number(),
                    insertion: res.insertion_code().unwrap_or("").into(),
                    name: name.into(),
                    alternate: alt.into(),
                }),
            });
            for (aname, a) in atoms {
                let element = a
                    .element()
                    .ok_or_else(|| format!("Missing element for {name} {aname}."))?
                    .to_string();
                let element = crate::molecule::element_symbol(
                    crate::molecule::atomic_number(&element)
                        .ok_or_else(|| format!("Unknown element {element}."))?,
                )
                .ok_or("Unknown element")?
                .to_string();
                let (x, y, z) = a.pos();
                let pos = Vec3::new(x as f32, y as f32, z as f32);
                if !pos.is_finite() {
                    return Err("Non-finite protein coordinates.".into());
                }
                let index = mol.atoms.len();
                if serials.insert(a.serial_number(), index).is_some() {
                    return Err("Duplicate atom serial in the selected model.".into());
                }
                t.atoms.push(AtomInfo {
                    id: index as u64 + 1,
                    residue: Some(rid),
                    name: aname,
                    charge: if a.charge() != 0 {
                        i8::try_from(a.charge()).ok()
                    } else {
                        None
                    },
                });
                mol.atoms.push(element);
                mol.pos.push(pos);
            }
        }
    }
    if mol.atoms.is_empty() {
        return Err("No atoms in the selected model.".into());
    }
    let mut edges = BTreeMap::new();
    let mut lookup = HashMap::new();
    for (i, a) in t.atoms.iter().enumerate() {
        lookup.insert((a.residue.unwrap(), a.name.clone()), i);
    }
    let mut missing = 0;
    let mut unknown = 0;
    for r in &t.residues {
        let indices: Vec<_> = t
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.residue == Some(r.id))
            .map(|(i, _)| i)
            .collect();
        if let Some(template) = peptide_edit::template(&r.template) {
            for a in template
                .atoms
                .iter()
                .filter(|a| a.element != "H" && a.name != "OXT")
            {
                if !lookup.contains_key(&(r.id, a.name.into())) {
                    missing += 1;
                }
            }
            for bond in template.bonds {
                if let (Some(&a), Some(&b)) = (
                    lookup.get(&(r.id, template.atoms[bond.atoms[0]].name.into())),
                    lookup.get(&(r.id, template.atoms[bond.atoms[1]].name.into())),
                ) {
                    edge(
                        &mut edges,
                        a,
                        b,
                        match bond.order {
                            BondOrder::Single => 1,
                            BondOrder::Double => 2,
                            BondOrder::Triple => 3,
                            BondOrder::Aromatic => 4,
                        },
                    );
                }
            }
        } else {
            unknown += 1;
        }
        // Unknown heavy-atom chemistry is distance-inferred only within a residue.
        if peptide_edit::template(&r.template).is_none() {
            for (offset, &a) in indices.iter().enumerate() {
                for &b in &indices[offset + 1..] {
                    if mol.atoms[a] == "H" || mol.atoms[b] == "H" {
                        continue;
                    }
                    let d = mol.pos[a].distance(mol.pos[b]);
                    if d > 0.4
                        && d < 1.2
                            * (covalent_radius_angstrom(&mol.atoms[a])
                                + covalent_radius_angstrom(&mol.atoms[b]))
                    {
                        edge(&mut edges, a, b, 0);
                    }
                }
            }
        }
        // A hydrogen with an unrecognized name gets at most one nearest heavy
        // partner. Never infer H–H bonds or add a second bond to a named H.
        for &h in indices.iter().filter(|&&i| mol.atoms[i] == "H") {
            if edges.keys().any(|&(a, b)| a == h || b == h) {
                continue;
            }
            let nearest = indices
                .iter()
                .copied()
                .filter(|&i| mol.atoms[i] != "H")
                .filter_map(|i| {
                    let d = mol.pos[h].distance(mol.pos[i]);
                    (d > 0.4
                        && d < 1.2
                            * (covalent_radius_angstrom("H")
                                + covalent_radius_angstrom(&mol.atoms[i])))
                    .then_some((i, d))
                })
                .min_by(|a, b| a.1.total_cmp(&b.1));
            if let Some((i, _)) = nearest {
                edge(&mut edges, h, i, 1);
            }
        }
    }
    for pair in t.residues.windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        if a.chain != b.chain {
            continue;
        }
        let (Some(oa), Some(ob)) = (&a.origin, &b.origin) else {
            continue;
        };
        if ob.number != oa.number + 1 && !(ob.number == oa.number && ob.insertion > oa.insertion) {
            continue;
        }
        if let (Some(&c), Some(&n)) = (
            lookup.get(&(a.id, "C".into())),
            lookup.get(&(b.id, "N".into())),
        ) {
            if (1.0..1.8).contains(&mol.pos[c].distance(mol.pos[n])) {
                edge(&mut edges, c, n, 1);
            }
        }
    }
    if cif {
        cif_links(text, &serials, &mut edges, &mut warnings)?;
    } else {
        pdb_links(text, &serials, &t, &mut edges, &mut warnings);
    }
    t.bonds = edges.into_iter().map(|((a, b), o)| (a, b, o)).collect();
    t.elements = mol.atoms.clone();
    mol.topology = Some(t);
    mol.recompute_bonds(1.2, 2.5);
    if alternate_count > 0 {
        warnings.push(format!(
            "Selected the highest-occupancy alternate in {alternate_count} residue(s)."
        ));
    }
    if missing > 0 {
        warnings.push(format!(
            "{missing} missing residue heavy atom(s); incomplete torsions cannot be edited."
        ));
    }
    if unknown > 0 {
        warnings.push(format!("{unknown} nonstandard residue/solvent group(s) retained; their internal bonds may be distance-inferred."));
    }
    warnings.push("Loaded the first model. Protonation is not determined by coordinates; residue defaults are neutral unless named explicitly. Review states before adding H. Total charge stays unknown where unspecified.".into());
    Ok(Imported {
        molecule: mol,
        warnings,
    })
}
fn edge(edges: &mut BTreeMap<(usize, usize), u8>, a: usize, b: usize, o: u8) {
    if a != b {
        let k = (a.min(b), a.max(b));
        edges
            .entry(k)
            .and_modify(|old| {
                if *old == 0 {
                    *old = o;
                }
            })
            .or_insert(o);
    }
}
fn pdb_links(
    text: &str,
    serials: &HashMap<usize, usize>,
    t: &Topology,
    edges: &mut BTreeMap<(usize, usize), u8>,
    warnings: &mut Vec<String>,
) {
    let find = |chain: &str, number: &str, ins: &str, name: &str, alt: &str| -> Option<usize> {
        let number = number.parse::<isize>().ok()?;
        let ids: Vec<_> = t
            .residues
            .iter()
            .filter(|r| {
                r.origin.as_ref().is_some_and(|o| {
                    o.chain == chain
                        && o.number == number
                        && o.insertion == ins
                        && (alt.is_empty() || o.alternate == alt)
                })
            })
            .map(|r| r.id)
            .collect();
        if ids.len() != 1 {
            return None;
        }
        t.atoms
            .iter()
            .position(|a| a.residue == Some(ids[0]) && a.name == canonical_name(name))
    };
    for l in text.lines() {
        if l.starts_with("CONECT") {
            let Some(a) = field(l, 6, 11)
                .parse::<usize>()
                .ok()
                .and_then(|s| serials.get(&s))
                .copied()
            else {
                continue;
            };
            for start in (11..l.len().min(31)).step_by(5) {
                if let Some(&b) = field(l, start, start + 5)
                    .parse::<usize>()
                    .ok()
                    .and_then(|s| serials.get(&s))
                {
                    edge(edges, a, b, 1);
                }
            }
        }
        let pair = if l.starts_with("SSBOND") {
            if !identity(field(l, 59, 65)) || !identity(field(l, 66, 72)) {
                warnings.push(
                    "Skipped a symmetry-related disulfide outside the loaded coordinates.".into(),
                );
                continue;
            }
            Some((
                find(
                    field(l, 15, 16),
                    field(l, 17, 21),
                    field(l, 21, 22),
                    "SG",
                    "",
                ),
                find(
                    field(l, 29, 30),
                    field(l, 31, 35),
                    field(l, 35, 36),
                    "SG",
                    "",
                ),
            ))
        } else if l.starts_with("LINK  ") {
            if !identity(field(l, 59, 65)) || !identity(field(l, 66, 72)) {
                warnings
                    .push("Skipped a symmetry-related LINK outside the loaded coordinates.".into());
                continue;
            }
            Some((
                find(
                    field(l, 21, 22),
                    field(l, 22, 26),
                    field(l, 26, 27),
                    field(l, 12, 16),
                    field(l, 16, 17),
                ),
                find(
                    field(l, 51, 52),
                    field(l, 52, 56),
                    field(l, 56, 57),
                    field(l, 42, 46),
                    field(l, 46, 47),
                ),
            ))
        } else {
            None
        };
        if let Some(pair) = pair {
            if let (Some(a), Some(b)) = pair {
                edge(edges, a, b, 1);
            } else {
                warnings.push("A recorded covalent link could not be resolved in the selected model/alternate.".into());
            }
        }
    }
}

// Small CIF tokenizer only for connection metadata missing from pdbtbx. Coordinate
// parsing and validation remain delegated to the library. Quoted values and
// semicolon text are handled, so unrelated free-text records cannot become tags.
#[derive(Clone)]
struct Token {
    value: String,
    quoted: bool,
}
fn tokens(text: &str) -> Result<Vec<Token>, String> {
    let b = text.as_bytes();
    let mut i = 0;
    let mut out = vec![];
    while i < b.len() {
        if b[i].is_ascii_whitespace() {
            i += 1;
            continue;
        }
        if b[i] == b'#' {
            while i < b.len() && b[i] != b'\n' {
                i += 1;
            }
            continue;
        }
        if b[i] == b';' && (i == 0 || b[i - 1] == b'\n') {
            let start = i + 1;
            i += 1;
            while i < b.len() && !(b[i] == b';' && b[i - 1] == b'\n') {
                i += 1;
            }
            if i == b.len() {
                return Err("Unterminated CIF text field.".into());
            }
            out.push(Token {
                value: text[start..i].into(),
                quoted: true,
            });
            i += 1;
            continue;
        }
        if matches!(b[i], b'\'' | b'"') {
            let quote = b[i];
            i += 1;
            let start = i;
            while i < b.len()
                && !(b[i] == quote && (i + 1 == b.len() || b[i + 1].is_ascii_whitespace()))
            {
                i += 1;
            }
            if i == b.len() {
                return Err("Unterminated CIF quote.".into());
            }
            out.push(Token {
                value: text[start..i].into(),
                quoted: true,
            });
            i += 1;
            continue;
        }
        let start = i;
        while i < b.len() && !b[i].is_ascii_whitespace() {
            i += 1;
        }
        out.push(Token {
            value: text[start..i].into(),
            quoted: false,
        });
    }
    Ok(out)
}
fn tables(text: &str) -> Result<Vec<Vec<HashMap<String, String>>>, String> {
    let ts = tokens(text)?;
    let mut i = 0;
    let mut out = vec![];
    let mut scalar = HashMap::new();
    while i < ts.len() {
        if ts[i].quoted || ts[i].value != "loop_" {
            if !ts[i].quoted && ts[i].value.starts_with('_') && i + 1 < ts.len() {
                let tag = ts[i].value.to_ascii_lowercase();
                if tag.starts_with("_struct_conn.") || tag.starts_with("_atom_site.") {
                    scalar.insert(tag, ts[i + 1].value.clone());
                }
                i += 2;
            } else {
                i += 1;
            }
            continue;
        }
        i += 1;
        let mut tags = vec![];
        while i < ts.len() && !ts[i].quoted && ts[i].value.starts_with('_') {
            tags.push(ts[i].value.to_ascii_lowercase());
            i += 1;
        }
        if tags.is_empty() {
            return Err("CIF loop without column names.".into());
        }
        let mut rows = vec![];
        while i < ts.len() {
            let t = &ts[i];
            if !t.quoted
                && (t.value.starts_with('_')
                    || t.value == "loop_"
                    || t.value == "stop_"
                    || t.value.starts_with("data_")
                    || t.value.starts_with("save_"))
            {
                break;
            }
            if i + tags.len() > ts.len() {
                return Err("Incomplete CIF loop row.".into());
            }
            if ts[i..i + tags.len()]
                .iter()
                .any(|t| !t.quoted && (t.value.starts_with('_') || t.value == "loop_"))
            {
                return Err("Incomplete CIF loop row.".into());
            }
            if tags[0].starts_with("_atom_site.") || tags[0].starts_with("_struct_conn.") {
                rows.push(
                    tags.iter()
                        .cloned()
                        .zip(ts[i..i + tags.len()].iter().map(|t| t.value.clone()))
                        .collect(),
                );
            }
            i += tags.len();
        }
        out.push(rows);
    }
    if !scalar.is_empty() {
        out.push(vec![scalar]);
    }
    Ok(out)
}
fn cif_links(
    text: &str,
    serials: &HashMap<usize, usize>,
    edges: &mut BTreeMap<(usize, usize), u8>,
    warnings: &mut Vec<String>,
) -> Result<(), String> {
    let tables = tables(text)?;
    let rows: Vec<_> = tables.iter().flatten().collect();
    let get = |r: &HashMap<String, String>, k: &str| {
        r.get(k).filter(|v| !empty(v)).cloned().unwrap_or_default()
    };
    let atoms: Vec<_> = rows
        .iter()
        .copied()
        .filter(|r| r.contains_key("_atom_site.id"))
        .collect();
    for r in rows
        .iter()
        .filter(|r| r.contains_key("_struct_conn.conn_type_id"))
    {
        let kind = get(r, "_struct_conn.conn_type_id");
        if !kind.starts_with("covale") && kind != "disulf" {
            continue;
        }
        if !(1..=2).all(|p| identity(&get(r, &format!("_struct_conn.ptnr{p}_symmetry")))) {
            warnings
                .push("Skipped a symmetry-related CIF link outside the loaded coordinates.".into());
            continue;
        }
        let find = |p: usize| -> Option<usize> {
            let mut matches = HashSet::new();
            for a in &atoms {
                let mut specified = 0;
                let mut matched = true;
                for suffix in [
                    "label_asym_id",
                    "label_seq_id",
                    "label_atom_id",
                    "label_comp_id",
                    "auth_asym_id",
                    "auth_seq_id",
                    "auth_atom_id",
                    "auth_comp_id",
                ] {
                    let expected = get(r, &format!("_struct_conn.ptnr{p}_{suffix}"));
                    if !expected.is_empty() {
                        specified += 1;
                        let actual = get(a, &format!("_atom_site.{suffix}"));
                        if actual != expected {
                            matched = false;
                            break;
                        }
                    }
                }
                for (link, site) in [
                    (
                        format!("_struct_conn.pdbx_ptnr{p}_pdb_ins_code"),
                        "_atom_site.pdbx_pdb_ins_code",
                    ),
                    (
                        format!("_struct_conn.pdbx_ptnr{p}_label_alt_id"),
                        "_atom_site.label_alt_id",
                    ),
                ] {
                    let expected = get(r, &link);
                    if !expected.is_empty() && get(a, site) != expected {
                        matched = false;
                    }
                }
                if matched && specified >= 3 {
                    if let Some(&i) = get(a, "_atom_site.id")
                        .parse::<usize>()
                        .ok()
                        .and_then(|s| serials.get(&s))
                    {
                        matches.insert(i);
                    }
                }
            }
            if matches.len() == 1 {
                matches.into_iter().next()
            } else {
                None
            }
        };
        if let (Some(a), Some(b)) = (find(1), find(2)) {
            let order = match get(r, "_struct_conn.pdbx_value_order").as_str() {
                "doub" | "DOUB" => 2,
                "trip" | "TRIP" => 3,
                _ => 1,
            };
            edge(edges, a, b, order);
        } else {
            warnings.push(
                "A CIF covalent link could not be resolved in the selected model/alternate.".into(),
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod connection_tests {
    use super::*;
    #[test]
    fn protein_cif_scalar_connection_resolves_and_symmetry_does_not() {
        let atom_rows="data_test\nloop_\n_atom_site.id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.label_atom_id\n1 A 1 SG\n2 A 5 SG\n#\n";
        let scalar="_struct_conn.conn_type_id disulf\n_struct_conn.ptnr1_label_asym_id A\n_struct_conn.ptnr1_label_seq_id 1\n_struct_conn.ptnr1_label_atom_id 'SG'\n_struct_conn.ptnr2_label_asym_id A\n_struct_conn.ptnr2_label_seq_id 5\n_struct_conn.ptnr2_label_atom_id SG\n";
        let serials = HashMap::from([(1, 0), (2, 1)]);
        let mut edges = BTreeMap::new();
        let mut warnings = vec![];
        cif_links(
            &(atom_rows.to_string() + scalar),
            &serials,
            &mut edges,
            &mut warnings,
        )
        .unwrap();
        assert_eq!(edges.get(&(0, 1)), Some(&1));
        edges.clear();
        cif_links(
            &(atom_rows.to_string() + scalar + "_struct_conn.ptnr2_symmetry 2_555\n"),
            &serials,
            &mut edges,
            &mut warnings,
        )
        .unwrap();
        assert!(edges.is_empty());
        assert_eq!(warnings.len(), 1);
    }
}
