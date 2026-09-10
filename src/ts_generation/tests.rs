use super::parameters::{atom_indices, ReactionKind};
use super::*;
use crate::optimizer::rda::{RdaConditionalCoordinates, RdaDistanceMode};
use std::path::PathBuf;
use std::sync::atomic::Ordering;

struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        Self(create_xtb_run_dir(&std::env::temp_dir().join("beavyr_ts_tests")).unwrap())
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn endpoints() -> (Endpoint, Endpoint) {
    let endpoint = |distance: f32| Endpoint {
        name: "H2".into(),
        atoms: vec!["H".into(), "H".into()],
        positions: vec![Vec3::ZERO, Vec3::X * distance * BOHR_TO_ANGSTROM as f32],
    };
    (endpoint(2.0), endpoint(4.0))
}

#[test]
fn transfer_and_angle_presets_preserve_atom_order_and_use_zero_based_indices() {
    let mut parameters = Parameters {
        reaction: ReactionKind::Transfer,
        indices: "3 1 2".into(),
        ..Default::default()
    };
    let opts = parameters.rda_options(3).unwrap();
    assert_eq!(opts.reactive_bonds, [[2, 0], [0, 1]]);
    parameters.reaction = ReactionKind::Angle;
    let opts = parameters.rda_options(3).unwrap();
    assert_eq!(opts.reactive_angles, [[2, 0, 1]]);
    assert!(opts.reactive_bonds.is_empty());
}

#[test]
fn custom_coordinates_and_weights_follow_bond_angle_dihedral_order() {
    let parameters = Parameters {
        reaction: ReactionKind::Custom,
        bonds: "1 2; 2 3".into(),
        angles: "1 2 3".into(),
        dihedrals: "1 2 3 4".into(),
        weights: "1, 2, 3, 4".into(),
        active_atoms: "2, 4".into(),
        align_atoms: "1 3".into(),
        ..Default::default()
    };
    let opts = parameters.rda_options(4).unwrap();
    assert_eq!(opts.reactive_bonds, [[0, 1], [1, 2]]);
    assert_eq!(opts.reactive_angles, [[0, 1, 2]]);
    assert_eq!(opts.reactive_dihedrals, [[0, 1, 2, 3]]);
    assert_eq!(
        opts.reactive_coordinate_weights.unwrap(),
        [1.0, 2.0, 3.0, 4.0]
    );
    assert_eq!(opts.active_atoms.unwrap(), [1, 3]);
    assert_eq!(opts.align_atoms.unwrap(), [0, 2]);
}

#[test]
fn invalid_atom_numbers_are_rejected_before_the_optimizer_indexes_arrays() {
    for text in ["0 1", "1 4", "1 1", "-1 2", "a 2"] {
        assert!(atom_indices(text, 3).is_err(), "{text}");
    }
    let mut parameters = Parameters {
        indices: "1 2 3".into(),
        ..Default::default()
    };
    assert!(parameters.rda_options(3).is_err());
    parameters.reaction = ReactionKind::Custom;
    parameters.dihedrals = "1 2 3 5".into();
    assert!(parameters.rda_options(4).is_err());
}

#[test]
fn neb_requires_reactive_coordinates_and_at_least_three_images() {
    let mut parameters = Parameters {
        algorithm: Algorithm::PoorMansNeb,
        reaction: ReactionKind::ActiveAtoms,
        ..Default::default()
    };
    assert!(parameters.rda_options(3).is_err());
    parameters.reaction = ReactionKind::Bond;
    parameters.options.poormans_neb_images = 2;
    assert!(parameters.rda_options(3).is_err());
    parameters.options.poormans_neb_images = 5;
    assert!(parameters.rda_options(3).unwrap().poormans_neb);
}

#[test]
fn bad_weights_and_search_parameters_cannot_produce_nan_distances() {
    let mut parameters = Parameters::default();
    for weights in ["NaN", "-1", "0", "1 2"] {
        parameters.weights = weights.into();
        assert!(parameters.rda_options(3).is_err(), "{weights}");
    }
    parameters.weights.clear();
    parameters.options.distance_tol = f64::NAN;
    assert!(parameters.rda_options(3).is_err());
    parameters.options.distance_tol = 0.05;
    parameters.gamma_values = "0.1 1.0".into();
    assert!(parameters.rda_options(3).is_err());
}

#[test]
fn endpoints_must_match_in_order_and_have_a_real_geometry_change() {
    let (reactant, mut product) = endpoints();
    validate_endpoints(&reactant, &product).unwrap();
    assert!(validate_endpoints(&reactant, &reactant).is_err());
    product.atoms[0] = "He".into();
    assert!(validate_endpoints(&reactant, &product).is_err());
    product.atoms[0] = "H".into();
    product.positions[0].x = f32::NAN;
    assert!(validate_endpoints(&reactant, &product).is_err());
}

#[test]
fn multi_frame_files_require_an_explicit_endpoint_choice() {
    let scratch = Scratch::new();
    let path = scratch.0.join("endpoints.xyz");
    std::fs::write(
        &path,
        "2\nfirst\nH 0 0 0\nH 1 0 0\n2\nlast\nH 0 0 0\nH 2 0 0\n",
    )
    .unwrap();
    let err = Endpoint::from_path(&path).err().unwrap().to_string();
    assert!(err.contains("single-structure"));
}

/// Two wells at bond lengths 2 and 4 bohr with a known saddle at 3 bohr.
struct DoubleWell;
impl Objective for DoubleWell {
    fn energy(&mut self, x: &[f64]) -> Result<f64> {
        let mut grad = vec![0.0; x.len()];
        self.energy_gradient(x, &mut grad)
    }
    fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
        self.energy_gradient(x, grad)?;
        Ok(())
    }
    fn energy_gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<f64> {
        let delta = [x[3] - x[0], x[4] - x[1], x[5] - x[2]];
        let distance = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
        let s = distance - 3.0;
        let energy = (s * s - 1.0).powi(2);
        let derivative = 4.0 * s * (s * s - 1.0);
        for axis in 0..3 {
            grad[axis] = -derivative * delta[axis] / distance;
            grad[axis + 3] = -grad[axis];
        }
        Ok(energy)
    }
}

fn run_model(algorithm: Algorithm, directory: &Path) -> GenerationOutput {
    let (reactant, product) = endpoints();
    let mut parameters = Parameters {
        algorithm,
        ..Default::default()
    };
    parameters.options.align = false;
    parameters.options.conditional_coordinates = RdaConditionalCoordinates::Cartesian;
    parameters.options.distance_mode = RdaDistanceMode::ReactiveInternal;
    parameters.options.poormans_neb_images = 5;
    generate(
        &mut DoubleWell,
        &reactant,
        &product,
        parameters.rda_options(2).unwrap(),
        algorithm,
        "analytic double well".into(),
        directory,
    )
    .unwrap()
}

#[test]
fn rda_generates_and_exports_the_known_saddle_in_angstrom() {
    let scratch = Scratch::new();
    let output = run_model(Algorithm::Rda, &scratch.0);
    let frame = parse_multi_xyz(&output.guess_xyz).unwrap().remove(0);
    let distance = frame.pos[0].distance(frame.pos[1]) as f64;
    assert!(
        (distance - 3.0 * BOHR_TO_ANGSTROM).abs() < 1.0e-4,
        "{distance}"
    );
    assert!((output.energy_hartree - 1.0).abs() < 1.0e-5);
    assert!(output.frames.len() >= 3);
    assert!(output.profile_dat.is_none());
    assert!(!output.search_xyz.is_empty());
}

#[test]
fn neb_xyz_is_the_complete_relaxed_trajectory_with_the_peak_selected() {
    let scratch = Scratch::new();
    let output = run_model(Algorithm::PoorMansNeb, &scratch.0);
    assert_eq!(output.frames.len(), 5);
    assert_eq!(output.guess.points.len(), 5);
    assert_eq!(parse_multi_xyz(&output.guess_xyz).unwrap().len(), 1);
    for (index, frame) in output.frames.iter().enumerate() {
        let distance = frame.pos[0].distance(frame.pos[1]) as f64 / BOHR_TO_ANGSTROM;
        assert!(
            (distance - (2.0 + index as f64 * 0.5)).abs() < 1.0e-3,
            "image {index}: {distance}"
        );
    }
    assert!((output.energy_hartree - 1.0).abs() < 1.0e-5);
    assert!(output.path_xyz.contains("selected=true"));
    assert_eq!(
        output
            .profile_dat
            .as_ref()
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .count(),
        5
    );
    let saved = scratch.0.join("saved-path.xyz");
    backend::save_text(&saved, &output.path_xyz).unwrap();
    assert_eq!(
        parse_multi_xyz(&std::fs::read_to_string(saved).unwrap())
            .unwrap()
            .len(),
        5
    );
}

#[test]
fn showing_a_guess_stops_old_playback_and_keeps_reactive_bonds_dynamic() {
    let mut mol = Molecule::empty();
    let mut traj = TrajectoryState {
        playing: true,
        fixed_bonds: true,
        overlay_enabled: true,
        ..Default::default()
    };
    let (endpoint, _) = endpoints();
    assert!(show_frames(
        vec![TrajectoryFrame {
            atoms: endpoint.atoms,
            pos: endpoint.positions
        }],
        &mut mol,
        &mut traj
    ));
    assert!(!traj.playing);
    assert!(!traj.fixed_bonds);
    assert!(!traj.overlay_enabled);
    assert_eq!(mol.atoms.len(), 2);
    assert_eq!(traj.frames.len(), 1);
}

#[cfg(unix)]
fn fake_backend(scratch: &Scratch, body: &str) -> BackendConfig {
    use std::os::unix::fs::PermissionsExt;
    let binary = scratch.0.join("fake-xtb");
    std::fs::write(&binary, format!("#!/bin/sh\n{body}\n")).unwrap();
    std::fs::set_permissions(&binary, std::fs::Permissions::from_mode(0o755)).unwrap();
    BackendConfig {
        program: QcProgram::Xtb,
        binary,
        method: method::MethodConfig::default(),
        charge: 0,
        multiplicity: 1,
    }
}

#[cfg(unix)]
#[test]
fn failed_gradient_does_not_switch_to_a_different_energy_surface() {
    let scratch = Scratch::new();
    let config = fake_backend(
        &scratch,
        "echo called >> calls.txt\necho 'SCF failed' >&2\nexit 1",
    );
    let control = Arc::new(RunControl::default());
    let mut objective = EnergyGradient::new(
        config,
        vec!["H".into(), "H".into()],
        scratch.0.join("energy"),
        control.clone(),
    );
    assert!(objective
        .energy(&[0.0, 0.0, 0.0, 2.0, 0.0, 0.0])
        .unwrap_err()
        .to_string()
        .contains("SCF failed"));
    assert_eq!(control.evaluations.load(Ordering::Relaxed), 1);
    assert_eq!(
        std::fs::read_to_string(scratch.0.join("energy/calls.txt"))
            .unwrap()
            .lines()
            .count(),
        1
    );
}

#[cfg(unix)]
#[test]
fn cancelling_terminates_the_running_energy_process() {
    let scratch = Scratch::new();
    let config = fake_backend(&scratch, "exec sleep 30");
    let control = Arc::new(RunControl::default());
    let task_control = control.clone();
    let workdir = scratch.0.join("energy");
    let worker = std::thread::spawn(move || {
        let mut objective =
            EnergyGradient::new(config, vec!["H".into(), "H".into()], workdir, task_control);
        objective.energy(&[0.0, 0.0, 0.0, 2.0, 0.0, 0.0])
    });
    let wait = Instant::now();
    while control.evaluations.load(Ordering::Relaxed) == 0
        && wait.elapsed() < Duration::from_secs(3)
    {
        std::thread::sleep(Duration::from_millis(10));
    }
    let cancellation = Instant::now();
    control.cancel();
    assert!(worker
        .join()
        .unwrap()
        .unwrap_err()
        .to_string()
        .contains("cancelled"));
    assert!(cancellation.elapsed() < Duration::from_secs(5));
}

#[test]
fn fully_constrained_bond_is_fitted_before_its_energy_is_evaluated() {
    use crate::optimizer::constrained::{
        optimize_constrained_with_internals, ConstrainedOptimizationOptions, ConstraintTarget,
    };
    use crate::optimizer::internal_coords::InternalCoords;
    let result = optimize_constrained_with_internals(
        &mut DoubleWell,
        vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0],
        InternalCoords {
            nat: 2,
            bonds: vec![[0, 1]],
            angles: vec![],
            dihedrals: vec![],
        },
        vec![ConstraintTarget::bond_bohr([0, 1], 3.0)],
        ConstrainedOptimizationOptions::default(),
    )
    .unwrap();
    assert!(result.converged);
    assert_eq!(result.cycles, 0);
    assert!((result.final_values[0] - 3.0 * BOHR_TO_ANGSTROM).abs() < 1.0e-6);
    assert!((result.energy - 1.0).abs() < 1.0e-8);
}

#[test]
fn neb_relaxes_the_unconstrained_bond_in_each_interior_image() {
    struct Surface;
    impl Objective for Surface {
        fn energy(&mut self, x: &[f64]) -> Result<f64> {
            let distance = |offset: usize| {
                (0..3)
                    .map(|axis| (x[offset + axis] - x[axis]).powi(2))
                    .sum::<f64>()
                    .sqrt()
            };
            Ok(((distance(3) - 3.0).powi(2) - 1.0).powi(2) + (distance(6) - 2.5).powi(2))
        }
        fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
            let mut trial = x.to_vec();
            for i in 0..x.len() {
                trial[i] = x[i] + 1.0e-5;
                let plus = self.energy(&trial)?;
                trial[i] = x[i] - 1.0e-5;
                let minus = self.energy(&trial)?;
                trial[i] = x[i];
                grad[i] = (plus - minus) / 2.0e-5;
            }
            Ok(())
        }
    }
    let (mut reactant, mut product) = endpoints();
    for endpoint in [&mut reactant, &mut product] {
        endpoint.atoms[0] = "C".into();
        endpoint.atoms.push("H".into());
        endpoint
            .positions
            .push(Vec3::Y * 2.0 * BOHR_TO_ANGSTROM as f32);
    }
    let parameters = Parameters {
        algorithm: Algorithm::PoorMansNeb,
        options: RdaOptions {
            align: false,
            poormans_neb_images: 5,
            ..Parameters::default().options
        },
        ..Default::default()
    };
    let scratch = Scratch::new();
    let output = generate(
        &mut Surface,
        &reactant,
        &product,
        parameters.rda_options(3).unwrap(),
        Algorithm::PoorMansNeb,
        "double well and spectator bond".into(),
        &scratch.0,
    )
    .unwrap();
    for (index, frame) in output.frames.iter().enumerate().take(4).skip(1) {
        let reactive = frame.pos[0].distance(frame.pos[1]) as f64 / BOHR_TO_ANGSTROM;
        let relaxed = frame.pos[0].distance(frame.pos[2]) as f64 / BOHR_TO_ANGSTROM;
        assert!(
            (reactive - (2.0 + index as f64 * 0.5)).abs() < 1.0e-3,
            "image {index}: reactive {reactive}"
        );
        assert!(
            (relaxed - 2.5).abs() < 1.0e-3,
            "image {index}: spectator {relaxed}"
        );
    }
}

#[cfg(unix)]
#[test]
fn both_backends_preserve_units_method_and_cache_a_repeated_geometry() {
    for program in [QcProgram::Xtb, QcProgram::Behemoth] {
        let scratch = Scratch::new();
        let output = if program == QcProgram::Xtb {
            "cat > gradient <<'GRAD'\n$grad\ncycle = 1\n0 0 0 H\n2 0 0 H\n0.125 0 0\n-0.125 0 0\n$end\nGRAD\necho '| TOTAL ENERGY -1.25 Eh |'"
        } else {
            "echo 'Analytical TASI nuclear gradient\nUnits: Hartree/bohr\nEnergy: -1.25 Hartree\n1 0.125 0 0\n2 -0.125 0 0'"
        };
        let mut config = fake_backend(
            &scratch,
            &format!("echo called >> calls.txt\nprintf '%s\\n' \"$@\" > args.txt\n{output}"),
        );
        config.program = program;
        config.method.behemoth = method::BehemothMethod::Tasi;
        let control = Arc::new(RunControl::default());
        let workdir = scratch.0.join("energy");
        let mut objective = EnergyGradient::new(
            config,
            vec!["H".into(), "H".into()],
            workdir.clone(),
            control.clone(),
        );
        let x = [0.0, 0.0, 0.0, 2.123456789, 0.0, 0.0];
        assert_eq!(objective.energy(&x).unwrap(), -1.25);
        let mut gradient = [0.0; 6];
        objective.gradient(&x, &mut gradient).unwrap();
        assert_eq!(gradient, [0.125, 0.0, 0.0, -0.125, 0.0, 0.0]);
        assert_eq!(control.evaluations.load(Ordering::Relaxed), 1);
        let geometry = std::fs::read_to_string(workdir.join("geometry.xyz")).unwrap();
        let coordinate = geometry
            .lines()
            .nth(3)
            .unwrap()
            .split_whitespace()
            .nth(1)
            .unwrap()
            .parse::<f64>()
            .unwrap();
        assert!((coordinate - x[3] * BOHR_TO_ANGSTROM).abs() < 1.0e-10);
        let args = std::fs::read_to_string(workdir.join("args.txt")).unwrap();
        assert!(
            args.contains(if program == QcProgram::Xtb {
                "--gfn\n2"
            } else {
                "--method\ntasi"
            }),
            "{args}"
        );
        let mut changed = x;
        changed[3] += 0.01;
        objective.energy(&changed).unwrap();
        assert_eq!(control.evaluations.load(Ordering::Relaxed), 2);
    }
}

#[test]
fn both_panel_modes_render_without_changing_the_viewer() {
    use bevy_egui::egui;
    let ctx = egui::Context::default();
    let mut state = TsGeneration::default();
    let mut mol = Molecule::empty();
    let mut traj = TrajectoryState::default();
    for algorithm in [Algorithm::Rda, Algorithm::PoorMansNeb] {
        state.parameters.algorithm = algorithm;
        let (reactant, product) = endpoints();
        state.reactant = Some(reactant);
        state.product = Some(product);
        let scratch = Scratch::new();
        state.output = Some(run_model(algorithm, &scratch.0));
        let mut output = ctx.run_ui(egui::RawInput::default(), |ui| {
            assert!(!super::ui::panel(ui, &mut state, &mut mol, &mut traj));
        });
        assert!(!output.shapes.is_empty());
        output.textures_delta.clear();
        assert!(mol.atoms.is_empty());
        assert!(traj.frames.is_empty());
    }
}

#[test]
#[ignore = "requires an installed Behemoth executable; set BEAVYR_TEST_BEHEMOTH"]
fn real_behemoth_runs_both_generators_and_exports_xyz() {
    let binary =
        std::env::var("BEAVYR_TEST_BEHEMOTH").expect("set BEAVYR_TEST_BEHEMOTH to the executable");
    let reactant = Endpoint {
        name: "water".into(),
        atoms: vec!["O".into(), "H".into(), "H".into()],
        positions: vec![
            Vec3::ZERO,
            Vec3::new(0.96, 0.0, 0.0),
            Vec3::new(-0.24, 0.93, 0.0),
        ],
    };
    let mut product = reactant.clone();
    product.positions[1].x = 1.16;
    // These deformed water geometries exercise the external gradient pipeline;
    // they are not intended to represent a chemical transition state.
    for algorithm in [Algorithm::Rda, Algorithm::PoorMansNeb] {
        let scratch = Scratch::new();
        let mut parameters = Parameters {
            algorithm,
            gamma_values: "0.25 0.5 0.75".into(),
            ..Default::default()
        };
        parameters.options.max_conditional_steps = 3;
        parameters.options.max_bracket_attempts = 2;
        parameters.options.poormans_neb_images = 3;
        parameters.options.poormans_neb_maxiter = 5;
        let mut objective = EnergyGradient::new(
            BackendConfig {
                program: QcProgram::Behemoth,
                binary: std::fs::canonicalize(&binary).unwrap(),
                method: method::MethodConfig {
                    behemoth: method::BehemothMethod::Tasi,
                    ..Default::default()
                },
                charge: 0,
                multiplicity: 1,
            },
            reactant.atoms.clone(),
            scratch.0.join("energy"),
            Arc::new(RunControl::default()),
        );
        let output = generate(
            &mut objective,
            &reactant,
            &product,
            parameters.rda_options(3).unwrap(),
            algorithm,
            "Behemoth / TASI".into(),
            &scratch.0,
        )
        .unwrap();
        assert!(output.energy_hartree.is_finite());
        assert_eq!(output.guess.q.len(), 9);
        assert_eq!(parse_multi_xyz(&output.guess_xyz).unwrap().len(), 1);
        assert!(parse_multi_xyz(&output.path_xyz).unwrap().len() >= 3);
    }
}

#[test]
fn captured_endpoints_stay_independent_while_the_viewer_is_edited() {
    let (first, second) = endpoints();
    let mut mol = Molecule {
        atoms: first.atoms.clone(),
        pos: first.positions.clone(),
        bonds: vec![],
        hydrogen_bonds: vec![],
    };
    let a = Endpoint::capture(&mol, "A").unwrap();
    // Drawing or a completed geometry optimization replaces the displayed coordinates.
    mol.pos.clone_from(&second.positions);
    let b = Endpoint::capture(&mol, "B").unwrap();
    mol.pos[1] += Vec3::Y;
    assert_eq!(a.positions, first.positions);
    assert_eq!(b.positions, second.positions);
    let mut traj = TrajectoryState::default();
    for stored in [&a, &b] {
        assert!(show_frames(
            vec![TrajectoryFrame {
                atoms: stored.atoms.clone(),
                pos: stored.positions.clone(),
            }],
            &mut mol,
            &mut traj
        ));
        assert_eq!(mol.pos, stored.positions);
        mol.pos[1] += Vec3::Z;
    }
    assert_eq!(a.positions, first.positions);
    assert_eq!(b.positions, second.positions);
}

#[test]
fn neb_energy_plot_renders_from_a_completed_path_without_changing_the_viewer() {
    use bevy_egui::egui;
    let scratch = Scratch::new();
    let mut state = TsGeneration::default();
    state.output = Some(run_model(Algorithm::PoorMansNeb, &scratch.0));
    state.energy_plot_open = true;
    let before = state.output.as_ref().unwrap().path_xyz.clone();
    let ctx = egui::Context::default();
    let mut rect = None;
    for _ in 0..2 {
        let mut output = ctx.run_ui(egui::RawInput::default(), |ui| {
            rect = super::plot::window(ui.ctx(), &mut state);
        });
        assert!(!output.shapes.is_empty());
        output.textures_delta.clear();
    }
    assert!(rect.unwrap().is_finite());
    assert!(state.energy_plot_open);
    assert_eq!(state.output.as_ref().unwrap().path_xyz, before);
    assert!(!state.is_running());

    state.output = Some(run_model(Algorithm::Rda, &scratch.0));
    assert!(super::plot::window(&ctx, &mut state).is_none());
    assert!(!state.energy_plot_open);
}
