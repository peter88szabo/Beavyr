//! TStrail: transition-state conformers with diverse active bond lengths.
//!
//! After Matsutani, Shomura, Tsuji and Maeda, *J. Comput. Chem.* **2026**, *47*, e70480.
//! Anyone publishing results from this tool should cite that paper -- see [`CITATION`].
//!
//! # The idea
//!
//! The conventional way to find a transition state's conformers is to take one transition state,
//! hold its active bond lengths, and sample the rest of the molecule. Every conformer it returns
//! therefore has the *same* active bond length. That is only right if the real transition states
//! share one, and for large, flexible or ionic systems they do not: TStrail's enolate example
//! found 157 transition-state conformers with C-C distances from 1.96 to 2.42 Å.
//!
//! TStrail turns the problem around. Sample the conformers of the **reactant**, where the reacting
//! atoms are still far apart and nothing needs holding, then push each one over the barrier
//! separately. Because each pathway starts from a different reactant shape, each arrives at its
//! own transition state -- with its own active bond length.
//!
//! # What runs here
//!
//! 1. Conformers of the reactant, by the ordinary search: DREIDING, Cartesian, seconds.
//! 2. Any conformer whose bonding differs from the reactant is discarded as broken.
//! 3. From each survivor, an [`super::afir`] pathway with a constant force on the active bonds,
//!    driven by GFN2-xTB -- a force field cannot break a bond, so it cannot do this step.
//! 4. The highest point of each pathway, on the true energy, is a transition-state-like structure.
//! 5. Those are de-duplicated by heavy-atom RMSD and reported with their active bond lengths.
//!
//! # What does not run here
//!
//! Stated rather than glossed, because the difference matters to anyone comparing against the
//! paper:
//!
//! * **No LUP path refinement.** The paper optimises each pathway before taking its maximum, and
//!   notes that the unrefined maximum also works and simply "converges to the TS structure" less
//!   quickly and accurately. What comes out here is therefore a transition-state *guess*.
//! * **No transition-state optimisation or IRC.** Beavyr has both -- `optimizer::transition_state`
//!   and `optimizer::irc` -- but they are not wired to this yet. These structures are the input to
//!   that step, not the end of it.
//! * **No stereochemistry filter.** The paper discards pathways leading to the wrong stereoisomer
//!   by checking dihedrals at the forming centres. Every pathway is kept here.

use std::path::{Path, PathBuf};

use anyhow::{anyhow, Result};

use crate::forcefield::dreiding::objective::{BOHR_TO_ANGSTROM, HARTREE_TO_KCAL};
use crate::molecule::Molecule;
use crate::optimizer::conformer_search::rmsd::{aligned_rmsd, heavy_atom_indices};

use super::afir::{afir_path, ActiveBond, DEFAULT_FORCE_KCAL_PER_ANGSTROM};
use super::refine::{XtbObjective, GFN2_SINGLE_THREAD};
use super::{search, ConformerSettings};

/// The reference to give anyone who publishes results from this.
pub const CITATION: &str = "S. Matsutani, T. Shomura, N. Tsuji, S. Maeda, \"Development of a \
Conformational Sampling Method for Transition State Structures With Diverse Active Bond \
Lengths\", J. Comput. Chem. 2026, 47, e70480. doi:10.1002/jcc.70480";

/// How hard each pathway is pushed, and for how long.
#[derive(Debug, Clone, PartialEq)]
pub struct TrailSettings {
    /// The bonds forming or breaking, which the artificial force acts along.
    pub active: Vec<ActiveBond>,
    /// Artificial force, kcal/mol/Å.
    pub force_kcal_per_angstrom: f64,
    /// Cycles per pathway.
    ///
    /// TStrail stops at 30 rather than converging: the path only has to cross the barrier, and
    /// what happens afterwards is the product.
    pub cycles: usize,
    /// How many reactant conformers to push. Each is an independent pathway.
    pub max_pathways: usize,
    /// Heavy-atom RMSD below which two transition-state guesses are the same one, in Å.
    pub rmsd_threshold: f64,
    pub charge: i32,
    pub multiplicity: i32,
}

impl Default for TrailSettings {
    fn default() -> Self {
        Self {
            active: Vec::new(),
            force_kcal_per_angstrom: DEFAULT_FORCE_KCAL_PER_ANGSTROM,
            cycles: 30,
            max_pathways: 20,
            // The paper clusters reactant conformers at 0.75 Å. Transition-state guesses are
            // tighter together than reactant conformers, so this is smaller.
            rmsd_threshold: 0.35,
            charge: 0,
            multiplicity: 1,
        }
    }
}

/// One transition-state guess.
#[derive(Debug, Clone)]
pub struct TransitionStateGuess {
    /// Flat `[x0, y0, z0, ...]` in Å.
    pub positions_angstrom: Vec<f64>,
    /// kcal/mol above the lowest guess.
    pub relative_energy_kcal: f64,
    /// The active bond lengths, Å, in the order they were given. The spread across guesses is the
    /// thing this method exists to produce.
    pub active_lengths: Vec<f64>,
    /// Which reactant conformer this pathway started from.
    pub from_conformer: usize,
    /// Whether the pathway went over the barrier and came down, or was still climbing when it
    /// stopped. A guess from a pathway that never crossed is a weaker starting point.
    pub crossed: bool,
}

/// Everything a run produced.
#[derive(Debug, Clone)]
pub struct TrailOutcome {
    pub guesses: Vec<TransitionStateGuess>,
    /// Reactant conformers found, before filtering.
    pub reactant_conformers: usize,
    /// Conformers discarded because their bonding differed from the reactant's.
    pub discarded_broken: usize,
    /// Pathways actually run.
    pub pathways: usize,
    /// Pathways that failed, with the reason for the first one.
    pub failed: usize,
    pub first_failure: Option<String>,
    /// xTB calls made.
    pub xtb_calls: usize,
    pub atoms: Vec<String>,
}

impl TrailOutcome {
    /// The range of active bond lengths across the guesses, per active bond.
    ///
    /// This is the number to look at: a spread means the method did what it is for, and a single
    /// value means a fixed-length search would have found the same thing.
    pub fn active_length_range(&self) -> Vec<(f64, f64)> {
        let count = self
            .guesses
            .first()
            .map(|g| g.active_lengths.len())
            .unwrap_or(0);
        (0..count)
            .map(|which| {
                let lengths = self.guesses.iter().map(|g| g.active_lengths[which]);
                lengths.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
                    (lo.min(v), hi.max(v))
                })
            })
            .collect()
    }
}

/// Runs TStrail's first stage: reactant conformers, then a pathway from each.
///
/// Blocking, and slow in wall-clock terms -- one pathway is tens of xTB calls -- so call it on a
/// background task.
pub fn run(
    reactant: &Molecule,
    settings: &TrailSettings,
    conformer_settings: &ConformerSettings,
    binary: &Path,
    scratch: &Path,
) -> Result<TrailOutcome> {
    if settings.active.is_empty() {
        return Err(anyhow!(
            "no active bonds were given; TStrail needs to know which bonds are forming or breaking"
        ));
    }
    for bond in &settings.active {
        if bond.a >= reactant.atoms.len() || bond.b >= reactant.atoms.len() || bond.a == bond.b {
            return Err(anyhow!("an active bond names an atom that is not there"));
        }
    }

    // Stage 1: conformers of the reactant. Unconstrained, so this is the fast Cartesian search.
    let conformers = search(reactant, conformer_settings.clone())
        .map_err(|e| anyhow!("the reactant conformer search failed: {e}"))?;
    let reactant_conformers = conformers.conformers.len();

    // Stage 2: discard anything whose bonding no longer matches the reactant. A conformer search
    // works on torsions and ring puckers, so a broken structure means a generated geometry that
    // collided; pushing it over a barrier would produce nonsense.
    let reference = bond_set(reactant);
    let mut usable = Vec::new();
    let mut discarded_broken = 0usize;
    for (index, conformer) in conformers.conformers.iter().enumerate() {
        let mut candidate = reactant.clone();
        candidate.set_pos(
            conformer
                .positions_angstrom
                .chunks_exact(3)
                .map(|c| bevy::prelude::Vec3::new(c[0] as f32, c[1] as f32, c[2] as f32))
                .collect(),
        );
        candidate.recompute_bonds(1.2, 2.5);
        if bond_set(&candidate) == reference {
            usable.push((index, conformer.positions_angstrom.clone()));
        } else {
            discarded_broken += 1;
        }
        if usable.len() >= settings.max_pathways {
            break;
        }
    }
    if usable.is_empty() {
        return Err(anyhow!(
            "no reactant conformer survived the bonding check; all {reactant_conformers} differed \
             from the input structure"
        ));
    }

    // Stage 3: one AFIR pathway per surviving conformer.
    let mut raw: Vec<TransitionStateGuess> = Vec::new();
    let mut xtb_calls = 0usize;
    let mut failed = 0usize;
    let mut first_failure = None;
    let pathways = usable.len();

    for (slot, (conformer_index, positions)) in usable.iter().enumerate() {
        let mut objective = XtbObjective::new(
            scratch.join(format!("pathway_{slot:03}")),
            &reactant.atoms,
            binary,
            GFN2_SINGLE_THREAD,
            settings.charge,
            settings.multiplicity,
        );
        let start: Vec<f64> = positions.iter().map(|a| a / BOHR_TO_ANGSTROM).collect();
        match afir_path(
            &mut objective,
            start,
            &settings.active,
            settings.force_kcal_per_angstrom,
            settings.cycles,
        ) {
            Ok(path) => {
                let top = path.transition_state_like();
                raw.push(TransitionStateGuess {
                    positions_angstrom: top
                        .coordinates
                        .iter()
                        .map(|b| b * BOHR_TO_ANGSTROM)
                        .collect(),
                    // Absolute for now; made relative once the whole set is known.
                    relative_energy_kcal: top.energy,
                    active_lengths: top.active_lengths.clone(),
                    from_conformer: *conformer_index,
                    crossed: path.crossed,
                });
            }
            Err(problem) => {
                failed += 1;
                if first_failure.is_none() {
                    first_failure = Some(problem.to_string());
                }
            }
        }
        xtb_calls += objective.calls;
    }

    if raw.is_empty() {
        return Err(anyhow!(
            "every pathway failed{}",
            match &first_failure {
                Some(problem) => format!(": {problem}"),
                None => String::new(),
            }
        ));
    }

    // Stage 4: de-duplicate. Two pathways from similar reactant conformers often arrive at the
    // same transition state, and reporting it twice would overstate the diversity found.
    let heavy = heavy_atom_indices(&reactant.atoms);
    raw.sort_by(|a, b| a.relative_energy_kcal.total_cmp(&b.relative_energy_kcal));
    let mut guesses: Vec<TransitionStateGuess> = Vec::new();
    for candidate in raw {
        let duplicate = guesses.iter().any(|kept| {
            same_structure(
                &candidate.positions_angstrom,
                &kept.positions_angstrom,
                &heavy,
                settings.rmsd_threshold,
            )
        });
        if !duplicate {
            guesses.push(candidate);
        }
    }

    let lowest = guesses
        .iter()
        .map(|g| g.relative_energy_kcal)
        .fold(f64::INFINITY, f64::min);
    for guess in guesses.iter_mut() {
        guess.relative_energy_kcal = (guess.relative_energy_kcal - lowest) * HARTREE_TO_KCAL;
    }

    Ok(TrailOutcome {
        guesses,
        reactant_conformers,
        discarded_broken,
        pathways,
        failed,
        first_failure,
        xtb_calls,
        atoms: reactant.atoms.clone(),
    })
}

/// The bonding of a structure, as a set of index pairs, for comparing two conformers.
fn bond_set(mol: &Molecule) -> std::collections::BTreeSet<(usize, usize)> {
    mol.bonds
        .iter()
        .map(|&(i, j, _)| if i < j { (i, j) } else { (j, i) })
        .collect()
}

/// Whether two geometries are the same structure, by heavy-atom RMSD after alignment.
///
/// The same measure the conformer search uses to decide two conformers are one, so a transition
/// state judged distinct here would also have been judged distinct there.
fn same_structure(left: &[f64], right: &[f64], heavy: &[usize], threshold: f64) -> bool {
    // `mirror: false` -- a mirror image is a different transition state when there is any
    // stereochemistry to speak of, and merging the two would hide exactly the pathways a
    // stereoselectivity study is looking for.
    aligned_rmsd(left, right, heavy, false) <= threshold
}

#[cfg(test)]
mod tests {
    use super::*;

    fn molecule(xyz: &str) -> Molecule {
        let mut mol = Molecule::from_xyz(xyz);
        mol.recompute_bonds(1.2, 2.5);
        mol
    }

    const ETHANOL: &str = "9\n\n\
        C -1.20 0.20 0.00\nC 0.00 -0.60 0.00\nO 1.15 0.20 0.00\n\
        H -2.10 -0.42 0.00\nH -1.24 0.84 0.89\nH -1.24 0.84 -0.89\n\
        H 0.04 -1.25 0.88\nH 0.04 -1.25 -0.88\nH 1.90 -0.36 0.00\n";

    /// A bimolecular reactant -- two fragments not bonded to each other -- has no torsion
    /// connecting them, so the conformer search has nothing to sample and refuses.
    ///
    /// This is a real limitation of doing it this way, not a bug: what varies between conformers
    /// of a two-fragment complex is the fragments' relative position and orientation, which is a
    /// rigid-body degree of freedom the torsion-space search has no notion of. TStrail's own
    /// examples start from associated reactant complexes with internal flexibility.
    ///
    /// Pinned so the behaviour is known and the message is checked, rather than discovered by a
    /// user with an SN2 in front of them.
    #[test]
    fn a_two_fragment_reactant_has_nothing_to_sample() {
        // Methane and a hydroxyl radical, 3 Å apart: the atmospheric H-abstraction reactant.
        let mol = molecule(
            "7\n\n\
             C  0.000  0.000  0.000\n\
             H  0.630  0.890  0.000\n\
             H  0.630 -0.890  0.000\n\
             H -0.630  0.890  0.000\n\
             H -0.630 -0.890  0.000\n\
             O  3.000  0.000  0.000\n\
             H  3.970  0.000  0.000\n",
        );
        // Beavyr sees two separate molecules.
        assert_eq!(mol.bonds.len(), 5, "four C-H and one O-H, nothing between them");

        let error = run(
            &mol,
            &TrailSettings {
                // Pull the hydroxyl oxygen onto one of methane's hydrogens.
                active: vec![ActiveBond { a: 5, b: 1 }],
                ..TrailSettings::default()
            },
            &ConformerSettings::default(),
            Path::new("/nonexistent/xtb"),
            Path::new("/tmp"),
        )
        .expect_err("two rigid fragments have no conformers");
        let text = error.to_string();
        assert!(
            text.contains("conformer search failed"),
            "the message should point at the conformer search: {text}"
        );
    }

    /// The whole pipeline against the real xTB: reactant conformers, filtering, one AFIR pathway
    /// per conformer, de-duplication, and the spread of active bond lengths.
    ///
    /// Ignored: it needs the configured xTB and takes tens of seconds.
    ///
    /// ```text
    /// cargo test --release --bins tstrail_against_real_xtb -- --ignored --nocapture
    /// ```
    #[test]
    #[ignore = "needs the configured xTB executable"]
    fn tstrail_against_real_xtb() {
        use crate::conformer::Thoroughness;
        use crate::qchem_interfaces::program::QcProgram;

        let configured = crate::qchem_interfaces::config::load_path(QcProgram::Xtb);
        let Ok(binary) = crate::qchem_interfaces::xtb_optimize::resolve_program_executable(
            QcProgram::Xtb,
            &configured,
        ) else {
            println!("  no xTB configured; nothing to test against");
            return;
        };
        println!("  xTB: {}", binary.display());

        // Ethanol, pulling the hydroxyl oxygen onto a methyl hydrogen: a strained 1,3-shift, and
        // a small enough system to exercise the whole pipeline quickly. The chemistry is beside
        // the point here -- what is being tested is that pathways run, produce maxima, and give a
        // spread of active bond lengths.
        let mol = molecule(ETHANOL);
        println!("  reactant: {} atoms, {} bonds", mol.atoms.len(), mol.bonds.len());

        let scratch = std::env::temp_dir().join(format!("beavyr_tstrail_{}", std::process::id()));
        let outcome = run(
            &mol,
            &TrailSettings {
                // Oxygen is atom index 2, a methyl hydrogen index 3.
                active: vec![ActiveBond { a: 2, b: 3 }],
                cycles: 25,
                max_pathways: 4,
                ..TrailSettings::default()
            },
            &ConformerSettings {
                thoroughness: Thoroughness::Quick,
                ncore: 2,
                ..ConformerSettings::default()
            },
            &binary,
            &scratch,
        );
        let _ = std::fs::remove_dir_all(&scratch);
        let outcome = outcome.expect("the pipeline runs");

        println!(
            "  {} reactant conformers, {} discarded, {} pathways, {} failed, {} xTB calls",
            outcome.reactant_conformers,
            outcome.discarded_broken,
            outcome.pathways,
            outcome.failed,
            outcome.xtb_calls
        );
        for (index, guess) in outcome.guesses.iter().enumerate() {
            println!(
                "    guess {}: dE {:>6.2} kcal/mol, active {:.2} Å, crossed {}",
                index + 1,
                guess.relative_energy_kcal,
                guess.active_lengths[0],
                guess.crossed
            );
        }
        for (which, (lo, hi)) in outcome.active_length_range().iter().enumerate() {
            println!("  active bond {} spans {lo:.2} to {hi:.2} Å", which + 1);
        }

        assert!(!outcome.guesses.is_empty(), "no transition-state guess came back");
        assert!(outcome.xtb_calls > 0, "xTB was never called");
        // The pull must have shortened the active bond from its reactant value.
        let start = (mol.pos[2] - mol.pos[3]).length() as f64;
        for guess in &outcome.guesses {
            assert!(
                guess.active_lengths[0] < start,
                "the active bond went from {start:.2} to {:.2} Å; the force did not pull",
                guess.active_lengths[0]
            );
        }
        // Ordered by energy, since the table shows them ranked.
        for pair in outcome.guesses.windows(2) {
            assert!(pair[0].relative_energy_kcal <= pair[1].relative_energy_kcal + 1e-9);
        }
    }

    #[test]
    fn a_run_without_active_bonds_is_refused() {
        let mol = molecule(ETHANOL);
        let error = run(
            &mol,
            &TrailSettings::default(),
            &ConformerSettings::default(),
            Path::new("/nonexistent/xtb"),
            Path::new("/tmp"),
        )
        .expect_err("no active bonds means nothing to push");
        assert!(error.to_string().contains("active bonds"), "{error}");
    }

    #[test]
    fn an_active_bond_naming_a_missing_atom_is_refused() {
        let mol = molecule(ETHANOL);
        let error = run(
            &mol,
            &TrailSettings {
                active: vec![ActiveBond { a: 0, b: 99 }],
                ..TrailSettings::default()
            },
            &ConformerSettings::default(),
            Path::new("/nonexistent/xtb"),
            Path::new("/tmp"),
        )
        .expect_err("atom 99 is not there");
        assert!(error.to_string().contains("not there"), "{error}");

        // And an atom bonded to itself.
        let error = run(
            &mol,
            &TrailSettings {
                active: vec![ActiveBond { a: 2, b: 2 }],
                ..TrailSettings::default()
            },
            &ConformerSettings::default(),
            Path::new("/nonexistent/xtb"),
            Path::new("/tmp"),
        )
        .expect_err("an atom cannot bond to itself");
        assert!(error.to_string().contains("not there"), "{error}");
    }

    /// The range across guesses is what the method exists to produce, so it must be computed over
    /// the right axis: per active bond, not across them.
    #[test]
    fn the_active_length_range_is_per_bond() {
        let outcome = TrailOutcome {
            guesses: vec![
                TransitionStateGuess {
                    positions_angstrom: vec![0.0; 6],
                    relative_energy_kcal: 0.0,
                    active_lengths: vec![2.10, 3.00],
                    from_conformer: 0,
                    crossed: true,
                },
                TransitionStateGuess {
                    positions_angstrom: vec![0.0; 6],
                    relative_energy_kcal: 1.0,
                    active_lengths: vec![2.42, 2.80],
                    from_conformer: 1,
                    crossed: true,
                },
            ],
            reactant_conformers: 2,
            discarded_broken: 0,
            pathways: 2,
            failed: 0,
            first_failure: None,
            xtb_calls: 0,
            atoms: vec!["C".into(), "C".into()],
        };
        let ranges = outcome.active_length_range();
        assert_eq!(ranges.len(), 2);
        assert!((ranges[0].0 - 2.10).abs() < 1e-12 && (ranges[0].1 - 2.42).abs() < 1e-12);
        assert!((ranges[1].0 - 2.80).abs() < 1e-12 && (ranges[1].1 - 3.00).abs() < 1e-12);
    }

    /// Two identical geometries are one structure; a twisted one is not.
    #[test]
    fn duplicate_detection_uses_aligned_heavy_atom_rmsd() {
        const THRESHOLD: f64 = 0.35;
        let heavy = vec![0usize, 1, 2];
        let a = vec![0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 2.5, 1.0, 0.0];
        // The same structure, translated and rotated: RMSD after alignment is zero.
        let b = vec![5.0, 5.0, 5.0, 5.0, 6.5, 5.0, 4.0, 7.5, 5.0];
        assert!(same_structure(&a, &b, &heavy, THRESHOLD));

        // A genuinely different shape: the third atom swung right round. With only three atoms
        // the internal distances are the whole structure, so this changes two of the three.
        let c = vec![0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 1.5, 2.5, 0.0];
        // Comfortably past the threshold, so the test is about the comparison rather than about
        // a borderline number.
        let apart = aligned_rmsd(&a, &c, &heavy, false);
        assert!(
            apart > THRESHOLD * 1.3,
            "the test structures are {apart:.3} Å apart against a {THRESHOLD} Å threshold -- too \
             close to the boundary to be a meaningful test"
        );
        assert!(!same_structure(&a, &c, &heavy, THRESHOLD));
    }

    #[test]
    fn the_citation_names_the_paper() {
        assert!(CITATION.contains("Matsutani"));
        assert!(CITATION.contains("10.1002/jcc.70480"));
    }
}
