use super::*;

#[test]
fn tolerance_drafts_save_and_reset_only_the_selected_optimizer() {
    let mut panel = QcPanelState::default();
    for job in [
        JobType::Optimize,
        JobType::OptimizeThenFrequencies,
        JobType::TransitionState,
        JobType::Mecp,
        JobType::Scan,
    ] {
        panel.job = job;
        panel.config.search.driver = job::search::Driver::Beavyr;
        let source = panel.config.clone();
        let target = target(&panel).unwrap();
        let mut edit = Editor::new(&panel, target, false);
        edit.values.energy = 2e-8;
        edit.values.max_gradient = 4e-6;
        assert_eq!(panel.config, source, "Drafts must not change live settings");
        edit.save(&mut panel.config).unwrap();
        let reread = Editor::new(&panel, target, false);
        assert_eq!(reread.values.energy, 2e-8);
        assert_eq!(reread.values.max_gradient, 4e-6);
        let defaults = Editor::new(&panel, target, true);
        defaults.save(&mut panel.config).unwrap();
        assert_eq!(Editor::new(&panel, target, false).values, defaults.values);
        assert_eq!(panel.config.mecp.multiplicity_b, source.mecp.multiplicity_b);
        assert_eq!(panel.config.search.driver, source.search.driver);
    }
}

#[test]
fn tolerance_native_defaults_and_programs_are_independent() {
    let mut panel = QcPanelState::default();
    for program in [
        QcProgram::Orca,
        QcProgram::PySCF,
        QcProgram::Psi4,
        QcProgram::Dreiding,
        QcProgram::Behemoth,
        QcProgram::SparrowPy,
    ] {
        panel.program = program;
        let mut edit = Editor::new(&panel, Target::Native, false);
        assert!(edit.native_defaults);
        edit.native_defaults = false;
        edit.values.rms_gradient = 8e-7;
        edit.save(&mut panel.config).unwrap();
        assert_eq!(
            panel.config.optimization.get(program).unwrap().rms_gradient,
            8e-7
        );
        let reset = Editor::new(&panel, Target::Native, true);
        reset.save(&mut panel.config).unwrap();
        assert!(panel.config.optimization.get(program).is_none());
    }
    let mut edit = Editor::new(&panel, Target::Native, false);
    edit.values.energy = f64::NAN;
    let original = panel.config.clone();
    assert!(edit.save(&mut panel.config).is_err());
    assert_eq!(panel.config, original);
}

#[test]
fn tolerance_mecp_defaults_keep_original_units_and_settings_serialize() {
    let options = job::mecp::Options::default();
    let original = crate::optimizer::spin_crossing::SpinCrossingOptions::default();
    assert_eq!(options.max_step, original.max_step);
    assert_eq!(options.rms_gradient, original.rms_gradient);
    assert_eq!(original.initial_inverse_hessian, 0.7);
    let legacy =
        ron::from_str::<job::mecp::Options>("(multiplicity_b:5,gap_tolerance:0.0001)").unwrap();
    assert_eq!(legacy.multiplicity_b, 5);
    assert_eq!(legacy.max_gradient, original.max_gradient);
    let mut config = MethodConfig::default();
    config.mecp.max_gradient = 9e-6;
    let serialized = ron::to_string(&config).unwrap();
    assert_eq!(ron::from_str::<MethodConfig>(&serialized).unwrap(), config);
}

#[test]
fn tolerance_native_input_keywords_and_displacement_units_reach_backends() {
    use crate::qchem_interfaces::{orca_run, psi4_run, pyscf_run};
    let mut config = MethodConfig::default();
    let t = Criteria {
        energy: 2e-7,
        max_step: 0.003,
        rms_step: 0.001,
        max_gradient: 4e-5,
        rms_gradient: 3e-5,
    };
    for p in [QcProgram::Orca, QcProgram::PySCF, QcProgram::Psi4] {
        config.optimization.set(p, Some(t.clone()));
    }
    let orca = orca_run::panel_input_text(
        JobType::Optimize,
        job::QcMethod::Hf,
        0,
        1,
        &config,
        Default::default(),
        "",
        "H 0 0 0\nH 0 0 0.74\n",
    );
    for keyword in ["TolE", "TolMaxD", "TolRMSD", "TolMaxG", "TolRMSG"] {
        assert!(orca.contains(keyword));
    }
    let atoms = vec!["H".into(), "H".into()];
    let x = [0., 0., 0., 0., 0., 0.74];
    let python = pyscf_run::script_text(
        JobType::Optimize,
        0,
        1,
        &config,
        Default::default(),
        false,
        &atoms,
        &x,
    );
    assert!(python.contains(&format!(
        "'convergence_dmax': {:.12e}",
        t.max_step * 0.529177210903
    )));
    assert!(python.contains("kernel(mf, callback=record, **geometry_tolerances)"));
    let psi4 = psi4_run::body_text(
        JobType::Optimize,
        job::QcMethod::Hf,
        0,
        1,
        &config,
        Default::default(),
        false,
        &atoms,
        &x,
    );
    assert!(psi4.contains("'rms_force_g_convergence': 3.000000000000e-5"));
    let mut command = std::process::Command::new("behemoth");
    t.behemoth(&mut command);
    let args: Vec<_> = command.get_args().map(|a| a.to_str().unwrap()).collect();
    assert_eq!(
        args,
        [
            "--opt-energy-tol",
            "0.0000002",
            "--opt-max-step",
            "0.003",
            "--opt-rms-step",
            "0.001",
            "--opt-max-gradient",
            "0.00004",
            "--opt-rms-gradient",
            "0.00003"
        ]
    );
}

fn frame(
    ctx: &egui::Context,
    panel: &mut QcPanelState,
    events: Vec<egui::Event>,
) -> egui::FullOutput {
    ctx.run_ui(
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1100., 900.),
            )),
            events,
            ..Default::default()
        },
        |_| {
            window(ctx, panel, false);
        },
    )
}
#[test]
fn tolerance_window_buttons_save_reset_and_cancel_drafts() {
    let ctx = egui::Context::default();
    let mut panel = QcPanelState::default();
    panel.job = JobType::Mecp;
    fn click(ctx: &egui::Context, panel: &mut QcPanelState, label: &str) {
        let mut first = frame(ctx, panel, vec![]);
        first.textures_delta.clear();
        let mut out = frame(ctx, panel, vec![]);
        let pos = out
            .shapes
            .iter()
            .find_map(|s| match &s.shape {
                egui::Shape::Text(t) if t.galley.text() == label => {
                    Some(t.pos + t.galley.size() / 2.)
                }
                _ => None,
            })
            .unwrap_or_else(|| panic!("Missing button {label}"));
        out.textures_delta.clear();
        for pressed in [true, false] {
            let mut out = frame(
                ctx,
                panel,
                vec![
                    egui::Event::PointerMoved(pos),
                    egui::Event::PointerButton {
                        pos,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::NONE,
                    },
                ],
            );
            out.textures_delta.clear();
        }
    }
    let mut draft = Editor::new(&panel, Target::Mecp, false);
    draft.values.max_gradient = 7e-6;
    ctx.data_mut(|d| d.insert_temp(id(), draft.clone()));
    click(&ctx, &mut panel, "Save");
    assert_eq!(panel.config.mecp.max_gradient, 7e-6);
    draft.values.max_gradient = 9e-6;
    ctx.data_mut(|d| d.insert_temp(id(), draft.clone()));
    click(&ctx, &mut panel, "Cancel");
    assert_eq!(panel.config.mecp.max_gradient, 7e-6);
    ctx.data_mut(|d| d.insert_temp(id(), draft));
    click(&ctx, &mut panel, "Reset defaults");
    assert_eq!(
        panel.config.mecp.max_gradient, 7e-6,
        "Reset edits the draft until saved"
    );
    click(&ctx, &mut panel, "Save");
    assert_eq!(
        panel.config.mecp.max_gradient,
        job::mecp::Options::default().max_gradient
    );
}
