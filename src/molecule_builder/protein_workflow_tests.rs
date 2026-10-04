//! End-to-end protein import and conformational editing regressions.
use super::{
    peptide::{self, parse_sequence, Recipe},
    peptide_conformation::{self as conf, Kind, Preset},
    peptide_edit as edit, peptide_rotamers as rot, protein_import as import,
    topology::Topology,
};
use crate::molecule::Molecule;

fn chain(sequence: &str) -> Molecule {
    peptide::build(
        &Recipe {
            residues: parse_sequence(sequence).unwrap(),
            ..Default::default()
        },
        None,
    )
    .unwrap()
    .molecule
}
fn bonds_preserved(a: &Molecule, b: &Molecule) {
    assert_eq!(a.atoms, b.atoms);
    assert_eq!(a.topology, b.topology);
    for &(i, j, _) in &a.bonds {
        assert!(
            (a.pos[i].distance(a.pos[j]) - b.pos[i].distance(b.pos[j])).abs() < 0.0001,
            "bond {i}-{j} moved"
        );
    }
}
fn angle(m: &Molecule, id: u64, k: Kind, expected: f32) {
    assert!(
        conf::wrap(conf::measured(m, id, k).unwrap() - expected).abs() < 0.02,
        "{id} {k:?}"
    );
}
fn count_links(m: &Molecule, n1: &str, n2: &str) -> usize {
    let t = m.topology.as_ref().unwrap();
    t.bonds
        .iter()
        .filter(|&&(a, b, _)| {
            t.atoms[a].residue != t.atoms[b].residue
                && ((t.atoms[a].name == n1 && t.atoms[b].name == n2)
                    || (t.atoms[b].name == n1 && t.atoms[a].name == n2))
        })
        .count()
}

#[test]
fn protein_import_crambin_pdb_and_cif_agree_including_disulfides() {
    let pdb = import::parse(
        include_str!("../../tests/fixtures/proteins/1CRN.pdb"),
        false,
    )
    .unwrap();
    let cif = import::parse(include_str!("../../tests/fixtures/proteins/1CRN.cif"), true).unwrap();
    for m in [&pdb.molecule, &cif.molecule] {
        let t = m.topology.as_ref().unwrap();
        assert_eq!(t.residues.len(), 46);
        assert_eq!(m.atoms.len(), 327);
        assert_eq!(count_links(m, "SG", "SG"), 3);
        assert_eq!(count_links(m, "C", "N"), 45);
        assert!(t.valid_for(&m.atoms));
        assert_eq!(t.residues[0].origin.as_ref().unwrap().chain, "A");
        assert_eq!(t.residues[45].origin.as_ref().unwrap().number, 46);
        let saved = ron::to_string(t).unwrap();
        assert_eq!(*t, ron::from_str::<Topology>(&saved).unwrap());
    }
    assert_eq!(pdb.molecule.atoms, cif.molecule.atoms);
    assert_eq!(pdb.molecule.pos, cif.molecule.pos);
    assert_eq!(
        pdb.molecule.topology.as_ref().unwrap().bonds,
        cif.molecule.topology.as_ref().unwrap().bonds
    );
}
#[test]
fn protein_import_ubiquitin_keeps_water_and_backbone_connectivity() {
    let loaded = import::parse(
        include_str!("../../tests/fixtures/proteins/1UBQ.pdb"),
        false,
    )
    .unwrap();
    let t = loaded.molecule.topology.as_ref().unwrap();
    assert_eq!(
        t.residues
            .iter()
            .filter(|r| edit::template(&r.template).is_some())
            .count(),
        76
    );
    assert_eq!(
        t.residues.iter().filter(|r| r.template == "HOH").count(),
        58
    );
    assert_eq!(count_links(&loaded.molecule, "C", "N"), 75);
    assert!(loaded.warnings.iter().any(|s| s.contains("nonstandard")));
    let mut m = loaded.molecule;
    let before = m.pos.clone();
    super::hydrogens::complete(&mut m).unwrap();
    assert_eq!(&m.pos[..before.len()], &before);
    let t = m.topology.as_ref().unwrap();
    let n = edit::atom(&m, 2, "N").unwrap();
    assert!(t
        .bonds
        .iter()
        .any(|&(a, b, _)| (a == n && m.atoms[b] == "H") || (b == n && m.atoms[a] == "H")));
}
#[test]
fn protein_backbone_edit_reaches_angles_without_stretching_or_moving_upstream() {
    let m = chain("AAAAA");
    let original = m.pos.clone();
    let changed = conf::edit(&m, 3, Some(-61.0), Some(-42.0), Some(175.0)).unwrap();
    for (k, v) in [(Kind::Phi, -61.0), (Kind::Psi, -42.0), (Kind::Omega, 175.0)] {
        angle(&changed.molecule, 3, k, v);
    }
    bonds_preserved(&m, &changed.molecule);
    assert_eq!(m.pos, original);
    let t = m.topology.as_ref().unwrap();
    for (i, a) in t.atoms.iter().enumerate() {
        if a.residue == Some(1) {
            assert_eq!(m.pos[i], changed.molecule.pos[i]);
        }
    }
    assert!(conf::measured(&m, 1, Kind::Phi).is_none());
    assert!(conf::measured(&m, 5, Kind::Psi).is_none());
    assert!(conf::edit(&m, 3, Some(f32::NAN), None, None).is_err());
}
#[test]
fn protein_backbone_range_presets_and_ring_guards() {
    let m = chain("AAAAA");
    for (p, phi, psi) in [
        (Preset::Helix, -57.0, -47.0),
        (Preset::Strand, -135.0, 135.0),
    ] {
        let b = conf::preset(&m, 1, 5, p).unwrap();
        bonds_preserved(&m, &b.molecule);
        for id in 2..5 {
            angle(&b.molecule, id, Kind::Phi, phi);
            angle(&b.molecule, id, Kind::Psi, psi);
        }
    }
    let b = conf::preset(&m, 2, 3, Preset::TurnI).unwrap();
    angle(&b.molecule, 2, Kind::Phi, -60.0);
    angle(&b.molecule, 2, Kind::Psi, -30.0);
    angle(&b.molecule, 3, Kind::Phi, -90.0);
    angle(&b.molecule, 3, Kind::Psi, 0.0);
    assert!(conf::preset(&m, 1, 5, Preset::TurnI).is_err());
    let pro = chain("APA");
    assert!(conf::edit(&pro, 2, Some(-100.0), None, None).is_err());
    let cross = import::parse(
        include_str!("../../tests/fixtures/proteins/1CRN.pdb"),
        false,
    )
    .unwrap()
    .molecule;
    assert!(conf::edit(&cross, 10, Some(0.0), None, None).is_err());
}
#[test]
fn protein_rotamers_preserve_backbone_bonds_and_reject_stale_choices() {
    for aa in [
        "V", "I", "L", "S", "T", "C", "D", "N", "E", "Q", "H", "F", "Y", "W", "M", "K", "R",
    ] {
        let m = chain(&format!("A{aa}A"));
        let choices = rot::choices(&m, 2).unwrap();
        assert!(!choices.marginal);
        assert!(choices.candidates.len() > 1, "{aa}");
        assert!(choices.candidates[0].current);
        assert!(choices
            .candidates
            .iter()
            .skip(1)
            .all(|c| c.probability.unwrap().is_finite() && c.score.is_finite()));
        assert!(choices.candidates[1..]
            .windows(2)
            .all(|c| c[0].score <= c[1].score));
        for i in 1..choices.candidates.len() {
            let changed = rot::preview(&m, &choices, i).unwrap().molecule;
            bonds_preserved(&m, &changed);
            for (j, a) in m.topology.as_ref().unwrap().atoms.iter().enumerate() {
                if a.residue != Some(2)
                    || matches!(a.name.as_str(), "N" | "CA" | "C" | "O" | "H" | "HA")
                {
                    assert_eq!(m.pos[j], changed.pos[j], "{aa} {}", a.name);
                }
            }
        }
        let mut stale = m.clone();
        stale.pos[0].x += 0.1;
        assert!(rot::preview(&stale, &choices, 1).is_err());
    }
}
#[test]
fn protein_rotamers_terminals_repacking_and_crosslinks() {
    let m = chain("KVLKS");
    assert!(rot::choices(&m, 1).unwrap().marginal);
    let b = rot::repack_neighbors(&m, 3, 8.0).unwrap();
    bonds_preserved(&m, &b.molecule);
    for (i, a) in m.topology.as_ref().unwrap().atoms.iter().enumerate() {
        if a.residue == Some(3) || matches!(a.name.as_str(), "N" | "CA" | "C" | "O") {
            assert_eq!(m.pos[i], b.molecule.pos[i]);
        }
    }
    let cross = import::parse(
        include_str!("../../tests/fixtures/proteins/1CRN.pdb"),
        false,
    )
    .unwrap()
    .molecule;
    assert!(rot::choices(&cross, 3).is_err());
    assert!(rot::choices(&chain("APA"), 2).is_err());
}
#[test]
fn protein_import_export_keeps_all_peptide_hydrogens_monovalent() {
    let source = chain("AKD");
    let text = edit::pdb(&source).unwrap();
    let loaded = import::parse(&text, false).unwrap().molecule;
    assert_eq!(source.atoms.len(), loaded.atoms.len());
    let adj = peptide::adjacency(&loaded);
    for (i, a) in loaded.atoms.iter().enumerate() {
        if a == "H" {
            assert_eq!(
                adj[i].len(),
                1,
                "{}",
                loaded.topology.as_ref().unwrap().atoms[i].name
            );
        }
    }
    assert_eq!(source.bonds.len(), loaded.bonds.len());
}

fn atom_line(
    serial: usize,
    name: &str,
    alt: char,
    res: &str,
    number: isize,
    ins: char,
    x: f32,
    occupancy: f32,
    element: &str,
) -> String {
    format!("ATOM  {serial:>5} {name:^4}{alt}{res:>3} A{number:>4}{ins}   {x:>8.3}{:>8.3}{:>8.3}{occupancy:>6.2}{:>6.2}          {element:>2}  \n",0.0,0.0,20.0)
}
#[test]
fn protein_import_selects_coherent_alternate_and_first_model() {
    let mut text = String::from("MODEL        1\n");
    text += &atom_line(1, "N", ' ', "ALA", 8, 'A', 0.0, 1.0, "N");
    text += &atom_line(2, "CA", 'A', "ALA", 8, 'A', 1.4, 0.3, "C");
    text += &atom_line(3, "CB", 'A', "ALA", 8, 'A', 2.5, 0.3, "C");
    text += &atom_line(4, "CA", 'B', "ALA", 8, 'A', -1.4, 0.7, "C");
    text += &atom_line(5, "CB", 'B', "ALA", 8, 'A', -2.5, 0.7, "C");
    text += "ENDMDL\nMODEL        2\n";
    text += &atom_line(1, "N", ' ', "ALA", 8, 'A', 50.0, 1.0, "N");
    text += "ENDMDL\nEND\n";
    let loaded = import::parse(&text, false).unwrap();
    let m = loaded.molecule;
    assert_eq!(m.atoms.len(), 3);
    let t = m.topology.as_ref().unwrap();
    let o = t.residues[0].origin.as_ref().unwrap();
    assert_eq!(
        (&o.chain, o.number, &o.insertion, &o.alternate),
        (&"A".to_string(), 8, &"A".to_string(), &"B".to_string())
    );
    assert_eq!(m.pos[edit::atom(&m, 1, "CA").unwrap()].x, -1.4);
    let exported = edit::pdb(&m).unwrap();
    let again = import::parse(&exported, false).unwrap().molecule;
    let o = again.topology.as_ref().unwrap().residues[0]
        .origin
        .as_ref()
        .unwrap();
    assert_eq!((o.number, o.insertion.as_str()), (8, "A"));
}
#[test]
fn protein_import_ter_and_number_gaps_block_implicit_peptide_bonds() {
    for middle in ["TER       2      ALA A   1\n", ""] {
        let text = atom_line(1, "C", ' ', "ALA", 1, ' ', 0.0, 1.0, "C")
            + middle
            + &atom_line(3, "N", ' ', "ALA", 2, ' ', 1.34, 1.0, "N")
            + "END\n";
        let m = import::parse(&text, false).unwrap().molecule;
        assert_eq!(count_links(&m, "C", "N"), usize::from(middle.is_empty()));
    }
    let text = atom_line(1, "C", ' ', "ALA", 1, ' ', 0.0, 1.0, "C")
        + &atom_line(2, "N", ' ', "ALA", 3, ' ', 1.34, 1.0, "N")
        + "END\n";
    assert_eq!(
        count_links(&import::parse(&text, false).unwrap().molecule, "C", "N"),
        0
    );
}
#[test]
fn protein_import_explicit_conect_and_unknown_ligands() {
    let text = atom_line(1, "C1", ' ', "LIG", 1, ' ', 0.0, 1.0, "C")
        + &atom_line(2, "N1", ' ', "LIG", 2, ' ', 5.0, 1.0, "N")
        + "CONECT    1    2\nEND\n";
    let m = import::parse(&text, false).unwrap().molecule;
    assert_eq!(m.bonds.len(), 1);
    assert_eq!(m.topology.as_ref().unwrap().residues[0].template, "LIG");
    assert!(import::parse("not a PDB file", false).is_err());
    assert!(import::parse("data_bad\nloop_\n_atom_site.id\n'broken", true).is_err());
}
