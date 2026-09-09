//! Running the built-in DREIDING force field as if it were one of the external programs.
//!
//! The optimizer and frequency panels are built around a backend that is launched, watched and
//! read back from files. DREIDING is none of those things -- it runs in this process and answers
//! in milliseconds -- but it is offered in the same selectors because from a chemist's point of
//! view it answers the same question.
//!
//! Rather than teach the panels a second set of rules, the results are written in the format
//! [`super::behemoth`] already produces: a multi-frame XYZ whose comment lines read
//! `cycle <n> E = <value> Eh`. Every consumer downstream -- the energy plot, the trajectory
//! player, the summary window -- then works unchanged.

use std::fs;
use std::path::Path;

use crate::forcefield::dreiding::objective::{
    relax_angstrom, Relaxation, HARTREE_TO_KCAL,
};
use crate::forcefield::dreiding::DreidingTopology;
use crate::molecule::Molecule;
use crate::optimizer::minimum_cartesian::{
    CartesianMinimumMethod, CartesianMinimumOptions,
};

use super::behemoth::{OPTIMIZED_GEOMETRY_FILE, TRAJECTORY_FILE};
use super::xtb_optimize::{EnergyHistoryPoint, XtbOptimizationOutput};

/// Convergence for a geometry optimisation, as against a quick cleanup.
///
/// Tighter than `dreiding::objective::cleanup_options`, because this is the user explicitly asking
/// for an optimised structure rather than for a sketch to be tidied. It is still not
/// quantum-chemistry tight: there is no point converging a generic force field's geometry beyond
/// the accuracy of its own parameters.
fn optimization_options() -> CartesianMinimumOptions {
    CartesianMinimumOptions {
        method: CartesianMinimumMethod::Bfgs,
        max_cycles: 1_000,
        // ~0.12 kcal/mol/Å.
        max_gradient: 1.0e-4,
        rms_gradient: 5.0e-5,
        verbosity: 0,
        ..CartesianMinimumOptions::default()
    }
}

/// Reads a structure the way the force field needs it.
///
/// The bond threshold is the standard one rather than whatever the viewport is set to. A generous
/// display setting bonds atoms that are merely close -- in water, the two hydrogens -- and that
/// mistypes every atom downstream, so the force field would confidently return a wrong geometry.
pub(crate) fn molecule_for_force_field(input_xyz: &str) -> Molecule {
    let mut mol = Molecule::from_xyz(input_xyz);
    mol.recompute_bonds(1.2, 2.5);
    mol
}

/// Optimises a geometry with DREIDING, writing its output where the panels expect to find it.
pub fn optimize(workdir: &Path, input_xyz: &str) -> Result<XtbOptimizationOutput, String> {
    let mol = molecule_for_force_field(input_xyz);
    let topology = DreidingTopology::build(&mol).map_err(|e| e.to_string())?;
    let atoms = &mol.atoms;

    let start: Vec<f64> = mol
        .pos
        .iter()
        .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
        .collect();
    let (_relaxed, relaxation) = relax_angstrom(&topology, &start, optimization_options())
        .map_err(|e| format!("DREIDING optimisation failed: {e}"))?;

    let trajectory_text = trajectory_xyz(atoms, &relaxation);
    let energy_history = relaxation
        .trajectory_energies
        .iter()
        .enumerate()
        .map(|(index, energy_kcal)| EnergyHistoryPoint {
            iteration: index as u32 + 1,
            energy_hartree: energy_kcal / HARTREE_TO_KCAL,
        })
        .collect();

    // Written for the same reason the external backends write them: so a run's directory can be
    // inspected, and so the trajectory player and energy plot read from a file as they always do.
    let _ = fs::create_dir_all(workdir);
    let _ = fs::write(workdir.join(TRAJECTORY_FILE), &trajectory_text);
    if let Some(last) = relaxation.trajectory.last() {
        let _ = fs::write(
            workdir.join(OPTIMIZED_GEOMETRY_FILE),
            frame_xyz(atoms, last, "optimised with DREIDING"),
        );
    }

    Ok(XtbOptimizationOutput {
        trajectory_text,
        final_energy_hartree: Some(relaxation.final_energy / HARTREE_TO_KCAL),
        energy_history,
        log: summary_log(&topology, &relaxation),
    })
}

/// One XYZ frame, with `comment` on the second line.
fn frame_xyz(atoms: &[String], positions_angstrom: &[f64], comment: &str) -> String {
    let mut text = format!("{}\n{comment}\n", atoms.len());
    for (symbol, xyz) in atoms.iter().zip(positions_angstrom.chunks_exact(3)) {
        text += &format!(
            "{:<3} {:>14.8} {:>14.8} {:>14.8}\n",
            symbol, xyz[0], xyz[1], xyz[2]
        );
    }
    text
}

/// The whole optimisation as a multi-frame XYZ, in Behemoth's comment-line format so that
/// `behemoth::parse_trajectory_energies` reads the energies straight back out.
fn trajectory_xyz(atoms: &[String], relaxation: &Relaxation) -> String {
    let mut text = String::new();
    for (index, frame) in relaxation.trajectory.iter().enumerate() {
        let energy_hartree = relaxation
            .trajectory_energies
            .get(index)
            .copied()
            .unwrap_or(relaxation.final_energy)
            / HARTREE_TO_KCAL;
        text += &frame_xyz(
            atoms,
            frame,
            &format!("cycle {} E = {:.10} Eh", index + 1, energy_hartree),
        );
    }
    // A minimisation that converged immediately still needs one frame, or there is nothing to
    // show and nothing to plot.
    if text.is_empty() {
        text = frame_xyz(
            atoms,
            &[],
            &format!(
                "cycle 1 E = {:.10} Eh",
                relaxation.final_energy / HARTREE_TO_KCAL
            ),
        );
    }
    text
}

/// The text the summary window shows.
///
/// The energy breakdown is the useful part: a term that has gone wrong -- a mistyped atom, a
/// stretched bond read as the wrong order -- shows up here as an absurd component long before the
/// geometry itself explains anything.
fn summary_log(topology: &DreidingTopology, relaxation: &Relaxation) -> String {
    let mut log = String::new();
    log += "DREIDING force field optimisation\n";
    log += "=================================\n\n";
    log += &format!("cycles                {}\n", relaxation.cycles);
    log += &format!(
        "converged             {}\n",
        if relaxation.converged { "yes" } else { "no" }
    );
    if !relaxation.message.is_empty() {
        log += &format!("optimiser             {}\n", relaxation.message);
    }
    log += &format!(
        "initial energy        {:>14.4} kcal/mol\n",
        relaxation.initial_energy
    );
    log += &format!(
        "final energy          {:>14.4} kcal/mol\n",
        relaxation.final_energy
    );
    log += &format!(
        "change                {:>14.4} kcal/mol\n",
        relaxation.final_energy - relaxation.initial_energy
    );
    log += &format!(
        "largest force         {:>14.4} kcal/mol/A\n\n",
        relaxation.max_force
    );

    let b = &relaxation.breakdown;
    log += "Final energy by term (kcal/mol)\n";
    log += &format!("  bond stretch        {:>14.4}\n", b.bond);
    log += &format!("  angle bend          {:>14.4}\n", b.angle);
    log += &format!("  torsion             {:>14.4}\n", b.torsion);
    log += &format!("  inversion           {:>14.4}\n", b.inversion);
    log += &format!("  van der Waals       {:>14.4}\n", b.vdw);
    log += &format!("  hydrogen bond       {:>14.4}\n", b.hbond);
    log += &format!("  total               {:>14.4}\n\n", b.total());

    // Bond orders are perceived from geometry, so a badly drawn structure can be mistyped. Listing
    // the types makes a surprising result diagnosable instead of mysterious.
    log += "Perceived DREIDING atom types\n";
    for (index, atom_type) in topology.atom_types().iter().enumerate() {
        log += &format!("  {:>4}  {}\n", index + 1, atom_type);
    }
    if let Some(warning) = topology.radical_warning() {
        log += "\nOpen shell\n";
        log += &format!("  {warning}\n");
    }

    log += "\nMayo, Olafson & Goddard, J. Phys. Chem. 1990, 94, 8897.\n";
    log
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: &str = "3\n\nO 0.000 0.000 0.000\nH 0.980 0.000 0.000\nH -0.245 0.949 0.000\n";

    #[test]
    fn optimising_water_writes_a_readable_trajectory() {
        let dir = std::env::temp_dir().join(format!("beavyr_dreiding_test_{}", std::process::id()));
        let output = optimize(&dir, WATER).expect("water is optimisable");

        // The energies must be recoverable by the same parser Behemoth's output uses, or the
        // energy plot silently shows nothing.
        let energies = crate::qchem_interfaces::behemoth::parse_trajectory_energies(
            &output.trajectory_text,
        );
        assert!(!energies.is_empty(), "no cycle lines were written");
        assert_eq!(energies.len(), output.energy_history.len());

        // And it must be loadable as a trajectory.
        let frames = crate::trajectory::parse_multi_xyz(&output.trajectory_text)
            .expect("the trajectory must parse");
        assert_eq!(frames.len(), energies.len());
        assert_eq!(frames[0].atoms.len(), 3);

        assert!(output.final_energy_hartree.is_some());
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Optimisation must lower the energy, and the last trajectory energy must be the one
    /// reported as final -- otherwise the plot and the summary disagree.
    #[test]
    fn the_energy_falls_and_the_plot_agrees_with_the_summary() {
        let dir = std::env::temp_dir()
            .join(format!("beavyr_dreiding_test_e_{}", std::process::id()));
        let output = optimize(&dir, WATER).unwrap();

        let history = &output.energy_history;
        assert!(history.first().unwrap().energy_hartree >= history.last().unwrap().energy_hartree);
        let final_energy = output.final_energy_hartree.unwrap();
        assert!(
            (history.last().unwrap().energy_hartree - final_energy).abs() < 1.0e-9,
            "plot ends at {} but the summary says {final_energy}",
            history.last().unwrap().energy_hartree
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn the_summary_reports_the_energy_breakdown_and_the_atom_types() {
        let dir = std::env::temp_dir()
            .join(format!("beavyr_dreiding_test_s_{}", std::process::id()));
        let output = optimize(&dir, WATER).unwrap();

        for expected in [
            "bond stretch",
            "angle bend",
            "van der Waals",
            "Perceived DREIDING atom types",
            "O_3",
        ] {
            assert!(output.log.contains(expected), "missing {expected:?}");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A molecule the force field cannot describe must fail with its own message, naming the
    /// atom, rather than a generic optimiser error.
    #[test]
    fn an_unparameterised_molecule_is_refused_by_name() {
        let dir = std::env::temp_dir()
            .join(format!("beavyr_dreiding_test_r_{}", std::process::id()));
        let error = match optimize(&dir, "2\n\nLi 0.0 0.0 0.0\nF 0.0 0.0 1.564\n") {
            Err(message) => message,
            Ok(_) => panic!("lithium has no DREIDING parameters"),
        };
        assert!(error.contains("Li") || error.contains("DREIDING"), "{error}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The threshold used for the force field's own perception must be the standard one, not the
    /// loose default `Molecule::from_xyz` applies -- at that setting water's two hydrogens are
    /// "bonded" and every atom is mistyped.
    #[test]
    fn connectivity_uses_the_standard_threshold() {
        let mol = molecule_for_force_field(WATER);
        assert_eq!(mol.bonds.len(), 2, "bonds: {:?}", mol.bonds);
        assert!(
            Molecule::from_xyz(WATER).bonds.len() > 2,
            "the loose default should have been the thing we avoided"
        );
    }
}
