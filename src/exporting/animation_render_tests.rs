//! Opt-in GPU check of the real capture/encoding pipeline, without a desktop window.
use super::*;

#[test]
#[ignore = "requires a working GPU; run explicitly with four-core affinity"]
fn offscreen_rotation_exports_gif_without_changing_editor() {
    let mut app = App::new();
    app.add_plugins(
        DefaultPlugins
            .set(WindowPlugin {
                primary_window: None,
                exit_condition: bevy::window::ExitCondition::DontExit,
                ..default()
            })
            .set(bevy::render::RenderPlugin {
                synchronous_pipeline_compilation: true,
                ..default()
            })
            .disable::<bevy::winit::WinitPlugin>(),
    );
    let mut settings = MolSettings::default();
    settings.bg_color = Color::srgb(0.05, 0.08, 0.12);
    let mut molecule = Molecule::from_xyz("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n");
    molecule.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
    let original_positions = molecule.pos.clone();
    let mut orbit = OrbitCamera::default();
    orbit.radius = 5.;
    orbit.target = Vec3::ZERO;
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/animation-render-test.gif");
    app.insert_resource(molecule.clone())
        .insert_resource(orbit.clone())
        .insert_resource(State {
            width: 160,
            height: 120,
            request: Some(Request {
                frames: vec![TrajectoryFrame {
                    metadata: Default::default(),
                    atoms: molecule.atoms.clone(),
                    pos: molecule.pos.clone(),
                }],
                count: 4,
                rotation: true,
                template: molecule,
                settings,
                hidden: Default::default(),
                orbit,
                projection: Projection::default(),
                path: path.clone(),
                fixed_bonds: false,
            }),
            ..default()
        })
        .add_systems(Update, tick);
    app.finish();
    app.cleanup();
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(90);
    loop {
        app.update();
        let state = app.world().resource::<State>();
        assert!(state.error.is_none(), "{:?}", state.error);
        if state.job.is_none() && state.request.is_none() {
            assert!(state
                .message
                .as_deref()
                .is_some_and(|m| m.starts_with("Saved ")));
            break;
        }
        assert!(
            std::time::Instant::now() < deadline,
            "GPU capture timed out"
        );
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
    assert_eq!(app.world().resource::<Molecule>().pos, original_positions);
    assert_eq!(app.world().resource::<OrbitCamera>().radius, 5.);
    let bytes = std::fs::read(path).unwrap();
    assert_eq!(&bytes[..6], b"GIF89a");
    assert_eq!(bytes.last(), Some(&0x3b));
}
