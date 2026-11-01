use super::input::{aligned_gradient, HessianSource, StateModes};
use super::*;
use ndarray::Array2;
fn reference() -> StateModes {
    StateModes {
        atoms: vec!["O".into(), "H".into(), "H".into()],
        masses: vec![16., 1., 1.],
        coords: vec![0., 0., 0., 1.5, 1., 0., -1.5, 1., 0.],
        modes: Array2::zeros((9, 3)),
        frequencies: vec![1000.; 3],
        selected: vec![0, 1, 2],
        imaginary: 0,
    }
}
#[test]
fn exact_half_turn_alignment_preserves_geometry_and_gradient() {
    let r = reference();
    let gradient = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
    let grad = crate::qchem_interfaces::orca_engrad::OrcaGradient {
        atoms: r.atoms.clone(),
        coords_bohr: r
            .coords
            .chunks_exact(3)
            .flat_map(|p| [-p[0] + 3., -p[1] + 2., p[2] - 5.])
            .collect(),
        gradient_bohr: gradient
            .chunks_exact(3)
            .flat_map(|p| [-p[0], -p[1], p[2]])
            .collect(),
        energy_hartree: 0.,
    };
    let aligned = aligned_gradient(&grad, &r).unwrap();
    for (actual, expected) in aligned.iter().zip(gradient) {
        assert!((actual - expected).abs() < 1e-10, "{actual} != {expected}");
    }
}
#[test]
fn a_different_gradient_geometry_is_rejected() {
    let r = reference();
    let mut coords = r.coords.clone();
    coords[3] += 0.02;
    let grad = crate::qchem_interfaces::orca_engrad::OrcaGradient {
        atoms: r.atoms.clone(),
        coords_bohr: coords,
        gradient_bohr: vec![0.; 9],
        energy_hartree: 0.,
    };
    assert!(aligned_gradient(&grad, &r)
        .unwrap_err()
        .contains("not calculated at"));
}
#[test]
fn turbomole_gradient_retains_geometry_and_fortran_exponents() {
    let g = input::turbomole_gradient(
        "$grad\ncycle = 1\n0 0 0 o\n1.5 1 0 h\n-1.5 1 0 h\n1.0D-3 2.0d-3 0\n0 0 0\n0 0 0\n$end\n",
    )
    .unwrap();
    assert_eq!(g.atoms, vec!["o", "h", "h"]);
    assert_eq!(g.gradient_bohr[0], 0.001);
    assert_eq!(g.coords_bohr[3], 1.5);
}
#[test]
fn malformed_turbomole_gradient_is_rejected() {
    assert!(input::turbomole_gradient("$grad\n0 0 0 h\nNaN 0 0\n$end").is_err());
    assert!(input::turbomole_gradient("$grad\n0 0 h\n0 0 0\n$end").is_err());
}
#[test]
fn matrix_inputs_are_checked_before_indexing() {
    let r = reference();
    assert!(duschinsky(
        &r.modes,
        &r.frequencies,
        &r.coords,
        &r.modes,
        &r.frequencies,
        &r.coords,
        &r.masses,
        &[99],
        &[0]
    )
    .is_err());
    assert!(duschinsky(
        &r.modes,
        &r.frequencies,
        &r.coords,
        &r.modes,
        &r.frequencies,
        &r.coords,
        &r.masses,
        &[],
        &[0]
    )
    .is_err());
}
#[test]
fn band_rejects_unbounded_grid_before_allocation() {
    let a = FranckCondonAnalysis {
        modes: vec![ModeActivity {
            index: 0,
            frequency_cm1: 1000.,
            gradient_au: 0.,
            displacement: 1.,
            huang_rhys: 0.5,
            reorganisation_cm1: 500.,
            resonance_raman: 1.,
        }],
        reorganisation_cm1: 500.,
        unaccounted_fraction: 0.,
    };
    assert!(vibronic_band(&a, 0.2, BandKind::Absorption, 1e-12, 10)
        .unwrap_err()
        .to_string()
        .contains("grid"));
    assert!(vibronic_band(&a, 0.2, BandKind::Absorption, f64::NAN, 10).is_err());
}
#[test]
fn new_analysis_window_opens_without_a_uvvis_spectrum() {
    use bevy_egui::egui;
    let ctx = egui::Context::default();
    let mut state = VibronicState::default();
    let uvvis = crate::uvvis::UvVisState::default();
    ui::open(&ctx);
    let mut rect = None;
    let mut output = ctx.run_ui(egui::RawInput::default(), |_| {
        rect = ui::window(&ctx, &mut state, &uvvis);
    });
    output.textures_delta.clear();
    assert!(state.open);
    assert!(rect.is_some());
}
#[test]
fn bare_hessian_requires_explicit_matching_geometry() {
    let path =
        std::env::temp_dir().join(format!("beavyr-vibronic-bare-{}.hess", std::process::id()));
    std::fs::write(&path, "$hessian\n0\n$end\n").unwrap();
    let result = input::load(
        &HessianSource {
            path: path.display().to_string(),
            xyz: String::new(),
        },
        0.,
    );
    std::fs::remove_file(path).unwrap();
    assert!(result.unwrap_err().contains("matching XYZ"));
}
#[test]
fn provided_orca_files_load_and_export_finite_activities() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("examples");
    let settings = state::Settings {
        ground: HessianSource {
            path: root
                .join("Exc_root_OnlyFreq_from_optimized.hess")
                .display()
                .to_string(),
            xyz: String::new(),
        },
        gradient: root
            .join("Exc_root_OnlyFreq_from_optimized.engrad")
            .display()
            .to_string(),
        make_band: false,
        ..Default::default()
    };
    let result = analysis::run(settings).unwrap();
    let csv = result.csv();
    assert!(csv.contains("relative_resonance_raman"));
    assert!(!csv.contains("NaN"));
    let analysis::AnalysisData::FranckCondon { activity, .. } = result.data else {
        panic!("wrong analysis")
    };
    assert!(!activity.modes.is_empty());
    assert!(activity.reorganisation_cm1.is_finite());
}
#[test]
fn real_hessian_against_itself_has_zero_displacement_and_identity_rotation() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples/TS_Gamma-3-6_ZZ-S_11_Compound_1.hess");
    let source = HessianSource {
        path: path.display().to_string(),
        xyz: String::new(),
    };
    let settings = state::Settings {
        kind: state::AnalysisKind::Duschinsky,
        ground: source.clone(),
        excited: source,
        ..Default::default()
    };
    let result = analysis::run(settings).unwrap();
    assert!(result.csv().contains("excited_mode_"));
    let analysis::AnalysisData::Duschinsky(d) = result.data else {
        panic!("wrong analysis")
    };
    assert!(d.displacement.iter().all(|k| k.abs() < 1e-8));
    for ((i, j), v) in d.rotation.indexed_iter() {
        assert!((v - if i == j { 1. } else { 0. }).abs() < 1e-8);
    }
}

#[test]
fn analysis_window_retains_user_resize_in_both_dimensions() {
    use bevy_egui::egui;
    let ctx = egui::Context::default();
    let mut state = VibronicState::default();
    state.open = true;
    fn frame(
        ctx: &egui::Context,
        state: &mut VibronicState,
        events: Vec<egui::Event>,
    ) -> egui::Rect {
        let mut rect = None;
        let mut output = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1500., 1200.),
                )),
                events,
                ..Default::default()
            },
            |_| {
                rect = super::ui::window(ctx, state, &crate::uvvis::UvVisState::default());
            },
        );
        output.textures_delta.clear();
        rect.unwrap()
    }
    frame(&ctx, &mut state, vec![]);
    let initial = frame(&ctx, &mut state, vec![]);
    let corner = initial.right_bottom() - egui::vec2(2., 2.);
    frame(&ctx, &mut state, vec![egui::Event::PointerMoved(corner)]);
    frame(
        &ctx,
        &mut state,
        vec![egui::Event::PointerButton {
            pos: corner,
            button: egui::PointerButton::Primary,
            pressed: true,
            modifiers: egui::Modifiers::NONE,
        }],
    );
    let target = corner + egui::vec2(140., 100.);
    frame(&ctx, &mut state, vec![egui::Event::PointerMoved(target)]);
    frame(
        &ctx,
        &mut state,
        vec![egui::Event::PointerButton {
            pos: target,
            button: egui::PointerButton::Primary,
            pressed: false,
            modifiers: egui::Modifiers::NONE,
        }],
    );
    let enlarged = frame(&ctx, &mut state, vec![]);
    assert!(
        enlarged.width() > initial.width() + 100.,
        "{initial:?} -> {enlarged:?}"
    );
    assert!(
        enlarged.height() > initial.height() + 60.,
        "{initial:?} -> {enlarged:?}"
    );
}

#[test]
fn mismatched_atom_order_is_rejected() {
    assert!(input::same_atoms(&["O".into(), "H".into()], &["H".into(), "O".into()]).is_err());
    assert!(input::same_atoms(&["O".into(), "H".into()], &["o".into(), "h".into()]).is_ok());
}
#[test]
fn bare_hessian_preserves_all_atoms_and_uses_xyz_angstroms() {
    let directory = std::env::temp_dir().join(format!("beavyr-vibronic-h2-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let mut h = Array2::<f64>::zeros((6, 6));
    h[(2, 2)] = 0.4;
    h[(5, 5)] = 0.4;
    h[(2, 5)] = -0.4;
    h[(5, 2)] = -0.4;
    std::fs::write(
        directory.join("hessian"),
        format!(
            "$hessian\n{}\n$end\n",
            h.iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(" ")
        ),
    )
    .unwrap();
    std::fs::write(
        directory.join("equilibrium.xyz"),
        "2\nH2\nH 0 0 0\nH 0 0 0.74\n",
    )
    .unwrap();
    let source = HessianSource {
        path: directory.join("hessian").display().to_string(),
        xyz: directory.join("equilibrium.xyz").display().to_string(),
    };
    let modes = input::load(&source, 0.).unwrap();
    assert_eq!(modes.selected.len(), 1);
    assert!((modes.coords[5] - 0.74 / 0.529177210903).abs() < 1e-12);
    let omega = (0.8 / (modes.masses[0] * AMU_TO_ELECTRON_MASS)).sqrt();
    assert!((modes.frequencies[modes.selected[0]] * CM1_TO_AU - omega).abs() < 1e-10);
    let gradient = directory.join("gradient");
    std::fs::write(
        &gradient,
        format!(
            "$grad\n0 0 0 h\n0 0 {} h\n0 0 -0.1\n0 0 0.1\n$end\n",
            0.74 / 0.529177210903
        ),
    )
    .unwrap();
    let mut settings = state::Settings {
        ground: source.clone(),
        gradient: gradient.display().to_string(),
        quanta: 1,
        ..Default::default()
    };
    let truncated = analysis::run(settings.clone()).unwrap();
    assert!(truncated.notices.iter().any(|n| n.contains("retains only")));
    assert!(matches!(
        truncated.data,
        analysis::AnalysisData::FranckCondon { band: None, .. }
    ));
    settings.quanta = 12;
    let complete = analysis::run(settings).unwrap();
    assert!(complete.notices.is_empty());
    assert!(complete.csv().contains("energy_eV,bin_probability"));
    let analysis::AnalysisData::FranckCondon {
        band: Some(band), ..
    } = complete.data
    else {
        panic!("absorption band missing")
    };
    assert!((band.profile.iter().map(|p| p.1).sum::<f64>() - 1.).abs() < 1e-12);
    std::fs::remove_file(gradient).unwrap();
    std::fs::write(
        directory.join("equilibrium.xyz"),
        "2\nNo silent dummy removal from Hessian\nX 0 0 0\nH 0 0 0.74\n",
    )
    .unwrap();
    assert!(input::load(&source, 0.).is_err());
    for name in ["hessian", "equilibrium.xyz"] {
        std::fs::remove_file(directory.join(name)).unwrap();
    }
    std::fs::remove_dir(directory).unwrap();
}
