// src/main.rs
use bevy::prelude::*;
use bevy_egui::{EguiGlobalSettings, EguiPlugin, EguiPrimaryContextPass};

mod bond_order;
mod camera;
mod cli;
mod color_schemes;
mod custom_schemes;
mod conformer;
mod diagnostics;
mod events;
mod export_image;
mod forcefield;
mod hbonds;
mod measurements;
mod molecule;
mod molecule_builder;
mod normalmode;
mod orientation;
mod numerics;
mod optimizer;
mod orbitals;
mod picking;
mod project;
mod qchem_interfaces;
mod qchem_panel;
mod recent_dir;
mod rmsd;
mod scene;
mod settings;
mod spectrum;
mod trajectory;
mod structure_history;
mod ts_generation;
mod ui;
mod ui_background_menu;
mod ui_icons;
mod ui_layout;
mod ui_measurements;
mod ui_structure;
mod ui_style;
mod ui_view_panels;
mod uvvis;
mod vibronic;

#[derive(SystemSet, Debug, Hash, PartialEq, Eq, Clone)]
enum UpdateSet {
    Core,
    Extras,
}

use diagnostics::DiagnosticsCache;
use export_image::{ExportCounter, ExportSelectArea, ExportSettings};
use hbonds::{draw_hydrogen_bonds_dashed, HbondGizmos};
use scene::BondLineGizmos;
use scene::{
    center_camera_on_startup, draw_bond_lines, react_to_molecule_changed_mark_dirty,
    rebuild_if_dirty, setup, sync_axis_camera_to_main, sync_coordinates_if_dirty,
    update_axis_viewport_on_resize, update_light_positions, update_lighting_if_dirty,
    update_materials_if_dirty, update_meshes_if_dirty, GeometryCache,
};
use settings::MolSettings;
use trajectory::{TrajectoryOverlayAssets, TrajectoryState};
use ui::XyzBuffer;

use molecule_builder::builder_ui::{
    configure_builder_gizmos, draw_builder_highlights, BuilderGizmos, BuilderHighlightGizmos,
    EditorRotateState, ZMatrixBuilderState,
};

fn main() {
    // `--help` is answered before Bevy is built. Handling it in a `Startup`
    // system, as this first did, means the window is already on screen by the
    // time the text prints -- the user asked a question and got a flash of an
    // application.
    let cmd = cli::parse_args(std::env::args().skip(1));
    if cmd.help {
        println!("{}", cli::help_text());
        return;
    }

    App::new()
        .add_plugins(
            DefaultPlugins
                .set(ImagePlugin::default_nearest())
                .set(bevy::window::WindowPlugin {
                    primary_window: Some(Window {
                        title: format!("Beavyr {}", env!("CARGO_PKG_VERSION")),
                        ..default()
                    }),
                    ..default()
                })
                .set(bevy::log::LogPlugin {
                    // Bevy's own defaults, plus one addition.
                    //
                    // `sctk-adwaita` draws the window decorations on Wayland
                    // and warns about every button in the title bar it does
                    // not recognise. There is nothing to fix on our side and
                    // nothing a user can do about it, so it is silenced rather
                    // than left to scroll past on every launch. Its errors
                    // still show, in case a decoration ever genuinely fails.
                    filter: format!("{}sctk_adwaita=error", bevy::log::DEFAULT_FILTER),
                    ..default()
                }),
        )
        .add_plugins(EguiPlugin::default())
        .insert_resource(EguiGlobalSettings {
            auto_create_primary_context: false,
            // Let egui consume pointer and keyboard input before anything else
            // sees it.  `EguiWantsInput` is written in PostUpdate, so a guard
            // read in Update is always a frame stale; this clears the Bevy
            // message queues in PreUpdate instead, which is the only way the
            // viewport can be reliably kept out of input meant for a panel.
            enable_absorb_bevy_input_system: true,
            ..default()
        })
        // picker plugin (now active)
        .add_plugins(picking::PickerPlugin)
        .add_plugins(orbitals::OrbitalPlugin)
        // messages
        .add_message::<events::AtomPicked>()
        .add_message::<events::ViewportClicked>()
        .add_message::<events::MoleculeChanged>()
        // resources
        .insert_resource(MolSettings::default())
        .insert_resource(camera::OrbitCamera::default())
        .insert_resource(molecule::Molecule::empty())
        .insert_resource(XyzBuffer {
            text: String::new(),
            last_dir: None,
            current_file: None,
            warning: None,
            live_xyz_atoms: Vec::new(),
            live_xyz_pos: Vec::new(),
            live_xyz_lines: Vec::new(),
            label_atoms: Vec::new(),
            type_label_galleys: Vec::new(),
            index_label_galleys: Vec::new(),
            type_index_label_galleys: Vec::new(),
            element_filter_atoms: Vec::new(),
            element_filter_extras: Vec::new(),
        })
        // export resources
        .init_resource::<ui::UiPanelRegions>()
        .insert_resource(ExportSettings::default())
        .insert_resource(ExportCounter::default())
        .insert_resource(ExportSelectArea::default())
        // measurements
        .init_resource::<conformer::ConformerRun>()
        .init_resource::<ts_generation::TsGeneration>()
        .init_resource::<structure_history::StructureHistory>()
        .init_resource::<measurements::Measurements>()
        .init_resource::<DiagnosticsCache>()
        .init_gizmo_group::<measurements::MeasurementGizmos>()
        .init_gizmo_group::<measurements::MeasurementHighlightGizmos>()
        // hbonds
        .init_gizmo_group::<HbondGizmos>()
        // bond lines
        .init_gizmo_group::<BondLineGizmos>()
        .init_resource::<GeometryCache>()
        // molecule editor
        .init_resource::<EditorRotateState>()
        .init_resource::<ZMatrixBuilderState>()
        .init_gizmo_group::<BuilderGizmos>()
        .init_gizmo_group::<molecule_builder::cleanup::CleanupGizmos>()
        .init_gizmo_group::<BuilderHighlightGizmos>()
        // trajectory
        .init_resource::<TrajectoryState>()
        .init_resource::<TrajectoryOverlayAssets>()
        // xtb geometry optimization
        .init_resource::<qchem_interfaces::xtb_optimize::XtbOptimizationTask>()
        .init_resource::<qchem_interfaces::xtb_optimize::XtbPanelState>()
        .init_resource::<qchem_interfaces::xtb_freq::XtbFrequencyTask>()
        .init_resource::<qchem_interfaces::xtb_freq::XtbFreqPanelState>()
        // structure comparison (RMSD / Kabsch)
        .init_resource::<rmsd::RmsdState>()
        .init_resource::<qchem_panel::QcPanelState>()
        .init_resource::<qchem_panel::run::QcRunTask>()
        .init_resource::<uvvis::UvVisState>()
        .init_resource::<uvvis::run::SpectrumTask>()
        .init_resource::<ui_layout::UiLayout>()
        .init_resource::<project::ProjectState>()
        .add_systems(Update, project::process.before(UpdateSet::Core))
        // scene
        // Files named on the command line load before the camera is centred,
        // so the view frames what was loaded rather than an empty scene.
        .init_resource::<cli::StartupLoadReport>()
        .add_systems(
            Startup,
            (setup, cli::load_command_line_files, center_camera_on_startup).chain(),
        )
        // egui pass
        .add_systems(EguiPrimaryContextPass, ui::ui_panel)
        .add_systems(
            EguiPrimaryContextPass,
            ui_measurements::distance_labels_overlay,
        )
        // update
        .configure_sets(Update, UpdateSet::Core.before(UpdateSet::Extras))
        .add_systems(
            Update,
            (
                camera::orbit_camera_system,
                export_image::update_export_select_area,
                react_to_molecule_changed_mark_dirty,
                measurements::clear_measurements_on_structure_load,
                diagnostics::refresh_diagnostics_cache,
                sync_coordinates_if_dirty,
                update_lighting_if_dirty,
                update_materials_if_dirty,
                update_meshes_if_dirty,
                molecule_builder::builder_ui::sync_builder_on_structure_load,
                molecule_builder::selection::sync,
                rebuild_if_dirty,
                update_light_positions,
                update_axis_viewport_on_resize,
                sync_axis_camera_to_main,
            )
                .chain()
                .in_set(UpdateSet::Core),
        )
        .add_systems(
            Update,
            (
                // measurements
                measurements::configure_measurement_gizmos,
                measurements::handle_measurement_picking, // keep if you haven't migrated to events
                measurements::update_measurement_distances,
                measurements::draw_measurement_lines,
                // hbonds
                draw_hydrogen_bonds_dashed,
                hbonds::configure_hbond_gizmos,
                scene::configure_bond_line_gizmos,
                draw_bond_lines,
                // builder (no more local click handler)
                configure_builder_gizmos,
                molecule_builder::builder_ui::handle_builder_atom_picked,
                molecule_builder::builder_ui::handle_viewport_click,
                draw_builder_highlights,
                // trajectory
                trajectory::advance_trajectory,
                trajectory::rebuild_trajectory_overlay,
                // xtb geometry optimization
                qchem_interfaces::xtb_optimize::poll_xtb_optimization,
                qchem_interfaces::xtb_freq::poll_xtb_frequencies,
                // Nested because the outer tuple is at Bevy's 20-element
                // limit for a system set; nesting costs nothing and keeps
                // the ordering.
                (
                    molecule_builder::cleanup::configure,
                    molecule_builder::cleanup::poll,
                    molecule_builder::cleanup::draw,
                    qchem_panel::run::poll_qc_run,
                    uvvis::run::poll_spectrum_run,
                    conformer::poll_conformer_search,
                    ts_generation::poll_generation,
                    export_image::process_pending_canvas_captures,
                    export_image::cleanup_export_captures,
                ),
            )
                .in_set(UpdateSet::Extras),
        )
        .run();
}
