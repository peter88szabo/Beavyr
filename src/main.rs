// src/main.rs
use bevy::prelude::*;
use bevy_egui::{EguiPlugin, EguiPrimaryContextPass};

mod camera;
mod color_schemes;
mod molecule;
mod scene;
mod settings;
mod ui;
mod export_image;
mod hbonds;
mod measurements;
mod ui_measurements;
mod molecule_builder;
mod events;
mod picking;
mod trajectory;

#[derive(SystemSet, Debug, Hash, PartialEq, Eq, Clone)]
enum UpdateSet {
    Core,
    Extras,
}

use scene::{
    center_camera_on_startup,
    draw_bond_lines,
    react_to_molecule_changed_mark_dirty,
    rebuild_if_dirty,
    setup,
    sync_axis_camera_to_main,
    update_axis_viewport_on_resize,
    update_light_positions,
};
use settings::MolSettings;
use ui::XyzBuffer;
use export_image::{ExportCounter, ExportSelectArea, ExportSettings};
use hbonds::{HbondGizmos, draw_hydrogen_bonds_dashed};
use trajectory::TrajectoryState;
use scene::BondLineGizmos;

use molecule_builder::builder_ui::{
    builder_ui_panel, configure_builder_gizmos, draw_builder_highlights,
    EditorRotateState, BuilderGizmos, ZMatrixBuilderState
};

fn main() {
    App::new()
        .add_plugins(
            DefaultPlugins
                .set(ImagePlugin::default_nearest())
                .disable::<bevy::audio::AudioPlugin>(),
        )
        .add_plugins(EguiPlugin::default())
        // picker plugin (now active)
        .add_plugins(picking::PickerPlugin)
        // messages
        .add_message::<events::AtomPicked>()
        .add_message::<events::MoleculeChanged>()
        // resources
        .insert_resource(MolSettings::default())
        .insert_resource(camera::OrbitCamera::default())
        .insert_resource(molecule::Molecule::from_xyz(scene::DEFAULT_WATER))
        .insert_resource(XyzBuffer {
            text: scene::DEFAULT_WATER.trim().to_string(),
            last_dir: None,
            warning: None,
        })
        // export resources
        .insert_resource(ExportSettings::default())
        .insert_resource(ExportCounter::default())
        .insert_resource(ExportSelectArea::default())
        // measurements
        .init_resource::<measurements::Measurements>()
        .init_gizmo_group::<measurements::MeasurementGizmos>()
        // hbonds
        .init_gizmo_group::<HbondGizmos>()
        // bond lines
        .init_gizmo_group::<BondLineGizmos>()
        // molecule editor
        .init_resource::<EditorRotateState>()
        .init_resource::<ZMatrixBuilderState>()
        .init_gizmo_group::<BuilderGizmos>()
        // trajectory
        .init_resource::<TrajectoryState>()
        // scene
        .add_systems(Startup, (setup, center_camera_on_startup))
        // egui pass
        .add_systems(EguiPrimaryContextPass, ui::ui_panel)
        .add_systems(EguiPrimaryContextPass, ui_measurements::distance_labels_overlay)
        .add_systems(EguiPrimaryContextPass, builder_ui_panel)
        // update
        .configure_sets(Update, UpdateSet::Core.before(UpdateSet::Extras))
        .add_systems(
            Update,
            (
                camera::orbit_camera_system,
                export_image::update_export_select_area,
                react_to_molecule_changed_mark_dirty,
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
                hbonds::draw_hydrogen_bonds_dashed,
                scene::configure_bond_line_gizmos,
                draw_bond_lines,

                // builder (no more local click handler)
                configure_builder_gizmos,
                molecule_builder::builder_ui::handle_builder_atom_picked,
                draw_builder_highlights,

                // trajectory
                trajectory::advance_trajectory,
                trajectory::rebuild_trajectory_overlay,
                export_image::process_pending_canvas_captures,
                export_image::cleanup_export_captures,
            )
                .chain()
                .in_set(UpdateSet::Extras),
        )
        .run();
}
