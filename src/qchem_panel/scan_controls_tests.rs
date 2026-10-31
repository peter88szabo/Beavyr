use super::*;
#[test]
fn scan_picking_preserves_order_and_ignores_dummies_and_duplicates() {
    let mut app = App::new();
    app.add_message::<crate::events::ViewportClicked>();
    let mut p = QcPanelState::default();
    p.open = true;
    p.job = JobType::Scan;
    p.config.scan.picking = true;
    p.config.scan.coordinate = Coordinate::Angle;
    app.insert_resource(p)
        .init_resource::<QcRunTask>()
        .insert_resource(crate::molecule::Molecule::from_xyz(
            "4\n\nX 0 0 0\nO 0 0 1\nH 1 0 1\nH 0 1 1\n",
        ))
        .add_systems(Update, collect_picks);
    for hit in [Some(0), Some(3), Some(3), None, Some(1), Some(2)] {
        app.world_mut()
            .write_message(crate::events::ViewportClicked { hit });
    }
    app.update();
    let s = &app.world().resource::<QcPanelState>().config.scan;
    assert_eq!(s.atoms, "4 2 3");
    assert!(!s.picking);
}
#[test]
fn scan_settings_offer_native_driver_only_for_orca() {
    let ctx = egui::Context::default();
    let mol = crate::molecule::Molecule::from_xyz("2\n\nH 0 0 0\nH 1 0 0\n");
    let task = QcRunTask::default();
    for program in [QcProgram::Orca, QcProgram::Behemoth, QcProgram::Xtb] {
        let mut panel = QcPanelState::default();
        panel.program = program;
        panel.job = JobType::Scan;
        let mut output = None;
        for _ in 0..2 {
            let mut frame = ctx.run_ui(Default::default(), |ui| {
                settings(ui, &mut panel, &mol, &task)
            });
            frame.textures_delta.clear();
            output = Some(frame);
        }
        let labels: Vec<_> = output
            .unwrap()
            .shapes
            .into_iter()
            .filter_map(|s| match s.shape {
                egui::Shape::Text(t) => Some(t.galley.text().to_owned()),
                _ => None,
            })
            .collect();
        assert_eq!(
            labels.iter().any(|t| t == "ORCA built-in"),
            program == QcProgram::Orca
        );
        assert_eq!(
            labels.iter().any(|t| t == "Behemoth built-in"),
            program == QcProgram::Behemoth
        );
        assert!(labels.iter().any(|t| t == "Scan settings"));
        assert!(labels.iter().any(|t| t == "Pick atoms"));
    }
}
