//! Opt-in GPU check of the real capture/encoding pipeline, without a desktop window.
use super::*;

#[test]
#[ignore = "requires a working GPU; run explicitly with four-core affinity"]
fn offscreen_rotation_exports_gif_without_changing_editor() { run_export(false); }
#[test]
#[ignore="requires a working GPU; run explicitly with four-core affinity"]
fn offscreen_rich_export_includes_surface_halos_and_labels() { run_export(true); }
fn run_export(rich:bool) {
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
    let mut orbital=crate::orbitals::OrbitalState::default();
    if rich { let text=include_str!("../../tests/fixtures/behemoth_ch2o_canonicalMO.molden");let data=crate::orbitals::molden::parse_molden(text).unwrap();molecule.atoms=data.atoms.clone();molecule.pos=data.positions.clone();orbital.adopt(data,None,None,Some(text.into()));orbital.resolution=1;orbital.isovalue=0.07; }
    molecule.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
    let original_positions = molecule.pos.clone();
    let mut orbit = OrbitCamera::default();
    orbit.radius = 5.;
    orbit.target = Vec3::ZERO;
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(if rich{"target/animation-rich-render-test.gif"}else{"target/animation-render-test.gif"});
    let frame=TrajectoryFrame{atoms:molecule.atoms.clone(),pos:molecule.pos.clone(),metadata:crate::trajectory_data::Metadata::energy(-113.,false,"test".into())};
    let mut polar=crate::polarizability::State::default();
    let mut measurements=crate::measurements::Measurements::default();measurements.new_pair(0,1,molecule.pos[0].distance(molecule.pos[1]));
    let mut labels=crate::saved_views::Library::default();labels.show_labels=true;labels.atoms=molecule.atoms.clone();labels.labels.push(crate::saved_views::AtomLabel{atom:0,text:"Probe".into()});
    if rich {polar.results=Some(Arc::new(crate::polarizability::Results{frames:vec![frame.clone()],values:vec![crate::polarizability::data::compute(orbital.data.as_ref().unwrap()).unwrap()],method:"test".into()}));polar.labels=true;}
    let extras=if rich{extras::Extras::snapshot(extras::Options{surface:true,halos:true,measurements:true,labels:true,caption:true,..Default::default()},&orbital,&crate::trajectory_wavefunctions::State::default(),&polar,&crate::reactivity::State::default(),&measurements,&labels,&[frame],true).unwrap()}else{extras::Extras::default()};
    app.insert_resource(molecule.clone())
        .insert_resource(orbit.clone())
        .insert_resource(State {
            width: if rich{480}else{160},
            height: if rich{360}else{120},
            request: Some(Request {
                extras,
                source_indices: vec![0],
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
