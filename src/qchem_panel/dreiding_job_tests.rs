use super::*;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(0);

struct RunDir(PathBuf);

impl RunDir {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "beavyr-dreiding-job-{}-{}",
            std::process::id(),
            NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed),
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
}

impl Drop for RunDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn water() -> (Vec<String>, Vec<f64>) {
    (
        ["O", "H", "H"].map(str::to_owned).to_vec(),
        vec![0.0, 0.0, 0.0, 1.3, 0.0, 0.0, 0.0, 1.3, 0.0],
    )
}

fn run(
    job: JobType,
    path: &Path,
    atoms: &[String],
    positions: &[f64],
) -> Result<QcRunOutput, String> {
    super::run(job, path, atoms, positions, &AtomicBool::new(false))
}

#[test]
fn calculator_preserves_nonconvergence_after_the_bounded_restart() {
    let (atoms, positions) = water();
    for job in [JobType::Optimize, JobType::OptimizeThenFrequencies] {
        let dir = RunDir::new();
        let options = CartesianMinimumOptions {
            max_cycles: 1,
            ..dreiding_run::optimization_options()
        };
        let result = run_with_options(
            job,
            &dir.0,
            &atoms,
            &positions,
            options,
            &AtomicBool::new(false),
        )
        .unwrap();
        assert_eq!(result.converged, Some(false));
        assert!(result.log.contains("converged             no"));
        assert!(result.log.contains("cycles                2"));
        assert!(result.energetics.final_energy_hartree.unwrap().is_finite());
        assert!(
            crate::trajectory::parse_multi_xyz(result.trajectory_text.as_deref().unwrap()).is_ok()
        );
        assert_eq!(
            result.hessian_path.is_some(),
            job == JobType::OptimizeThenFrequencies
        );
    }
}

#[test]
fn calculator_preserves_a_successful_optimization_and_its_log() {
    let dir = RunDir::new();
    let (atoms, positions) = water();
    let result = run(JobType::Optimize, &dir.0, &atoms, &positions).unwrap();
    assert_eq!(result.converged, Some(true));
    assert!(result.log.contains("converged             yes"));
    assert!(result.log.contains("Final energy by term"));
}

#[test]
fn single_point_reports_the_force_field_energy_without_optimizing() {
    use crate::forcefield::dreiding::objective::HARTREE_TO_KCAL;
    use crate::forcefield::dreiding::{DreidingTopology, Workspace};

    let dir = RunDir::new();
    let (atoms, positions) = water();
    let molecule = dreiding_run::molecule_for_force_field("3\n\nO 0 0 0\nH 1.3 0 0\nH 0 1.3 0\n");
    let topology = DreidingTopology::build(&molecule).unwrap();
    let coords: Vec<_> = molecule
        .pos
        .iter()
        .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
        .collect();
    let expected = topology
        .energy_breakdown(&coords, &mut Workspace::new())
        .total()
        / HARTREE_TO_KCAL;
    let result = run(JobType::SinglePoint, &dir.0, &atoms, &positions).unwrap();
    assert!((result.energetics.final_energy_hartree.unwrap() - expected).abs() < 1e-10);
    assert!(
        expected.abs() > 1e-3,
        "use a geometry with a nonzero energy"
    );
    assert!(result.trajectory_text.is_none());
    assert!(result.hessian_path.is_none());
    assert_eq!(result.converged, None);
    assert!(!result.log.is_empty());
}

#[test]
fn saved_gradient_matches_energy_differences_in_hartree_per_bohr() {
    use crate::forcefield::dreiding::objective::BOHR_TO_ANGSTROM;
    use crate::qchem_interfaces::hessian_file::parse_turbomole_gradient;

    let dir = RunDir::new();
    let (atoms, positions) = water();
    let result = run(JobType::Gradient, &dir.0, &atoms, &positions).unwrap();
    let energy = run(JobType::SinglePoint, &dir.0, &atoms, &positions).unwrap();
    assert_eq!(
        result.energetics.final_energy_hartree,
        energy.energetics.final_energy_hartree
    );
    assert_eq!(result.converged, None);
    assert!(result.trajectory_text.is_none());
    assert_eq!(result.kept_dir.as_deref(), Some(dir.0.as_path()));
    let (_, path) = result
        .saveable
        .iter()
        .find(|(label, _)| label.contains("Gradient"))
        .unwrap();
    let text = fs::read_to_string(path).unwrap();
    let gradient = parse_turbomole_gradient(&text, atoms.len()).unwrap();
    // Check the export's coordinate units separately from its gradient units.
    let hydrogen_x: f64 = text
        .lines()
        .nth(3)
        .unwrap()
        .split_whitespace()
        .next()
        .unwrap()
        .parse()
        .unwrap();
    assert!((hydrogen_x - positions[3] / BOHR_TO_ANGSTROM).abs() < 1e-6);
    let step_bohr = 2e-3;
    for (index, actual) in gradient.iter().enumerate() {
        let mut plus = positions.clone();
        let mut minus = positions.clone();
        plus[index] += step_bohr * BOHR_TO_ANGSTROM;
        minus[index] -= step_bohr * BOHR_TO_ANGSTROM;
        let e_plus = run(JobType::SinglePoint, &dir.0, &atoms, &plus)
            .unwrap()
            .energetics
            .final_energy_hartree
            .unwrap();
        let e_minus = run(JobType::SinglePoint, &dir.0, &atoms, &minus)
            .unwrap()
            .energetics
            .final_energy_hartree
            .unwrap();
        let expected = (e_plus - e_minus) / (2.0 * step_bohr);
        assert!(
            (actual - expected).abs() < 2e-5,
            "component {index}: {actual} vs {expected}"
        );
    }
}

#[test]
fn energy_and_gradient_reject_unsupported_chemistry_instead_of_succeeding() {
    let dir = RunDir::new();
    let atoms = ["Li".into(), "F".into()];
    for job in [JobType::SinglePoint, JobType::Gradient] {
        let result = run(job, &dir.0, &atoms, &[0.0, 0.0, 0.0, 0.0, 0.0, 1.564]);
        assert!(
            result.is_err(),
            "{} must actually evaluate the structure",
            job.label()
        );
    }
}

#[test]
fn all_dreiding_jobs_honor_cancellation_before_starting() {
    let (atoms, positions) = water();
    for job in [
        JobType::SinglePoint,
        JobType::Gradient,
        JobType::Optimize,
        JobType::Hessian,
        JobType::Frequencies,
        JobType::OptimizeThenFrequencies,
    ] {
        let dir = RunDir::new();
        let result = super::run(job, &dir.0, &atoms, &positions, &AtomicBool::new(true));
        assert_eq!(
            result.err().as_deref(),
            Some("Cancelled."),
            "{}",
            job.label()
        );
        assert_eq!(fs::read_dir(&dir.0).unwrap().count(), 0);
    }
}
