// src/main.rs
use bevy::prelude::*;
use bevy_egui::{EguiPlugin, EguiPrimaryContextPass};

mod camera;
mod color_schemes;
mod molecule;
mod scene;
mod settings;
mod ui;
mod export;
mod export_ui;
mod hbonds;
mod measurements;
mod ui_measurements;
mod molecule_builder;
mod events;
mod picking;

use camera::orbit_camera_system;
use scene::{
    center_camera_on_startup,
    react_to_molecule_changed_mark_dirty,
    rebuild_if_dirty,
    setup,
    sync_axis_camera_to_main,
    update_axis_viewport_on_resize,
};
use settings::MolSettings;
use ui::XyzBuffer;
use hbonds::{HbondGizmos, draw_hydrogen_bonds_dashed};

use molecule_builder::builder_ui::{
    builder_ui_panel, configure_builder_gizmos, draw_builder_highlights,
    EditorRotateState, BuilderGizmos
};

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(ImagePlugin::default_nearest()))
        .add_plugins(EguiPlugin::default())
        .add_plugins(export::ExportPlugin)
        // picker plugin (now active)
        .add_plugins(picking::PickerPlugin)
        // events
        .add_event::<events::AtomPicked>()
        .add_event::<events::MoleculeChanged>()
        // resources
        .insert_resource(MolSettings::default())
        .insert_resource(camera::OrbitCamera::default())
        .insert_resource(molecule::Molecule::from_xyz(scene::DEFAULT_WATER))
        .insert_resource(XyzBuffer {
            text: scene::DEFAULT_WATER.trim().to_string(),
            last_dir: None,
        })
        // export UI resources
        .init_resource::<export::ExportUiState>()
        .init_resource::<export::ExportRequestQueue>()
        // measurements
        .init_resource::<measurements::Measurements>()
        .init_gizmo_group::<measurements::MeasurementGizmos>()
        // hbonds
        .init_gizmo_group::<HbondGizmos>()
        // molecule editor
        .init_resource::<EditorRotateState>()
        .init_gizmo_group::<BuilderGizmos>()
        // scene
        .add_systems(Startup, (setup, center_camera_on_startup))
        // egui pass
        .add_systems(EguiPrimaryContextPass, ui::ui_panel)
        .add_systems(EguiPrimaryContextPass, ui_measurements::distance_labels_overlay)
        .add_systems(EguiPrimaryContextPass, builder_ui_panel)
        // update
        .add_systems(
            Update,
            (
                camera::orbit_camera_system,
                react_to_molecule_changed_mark_dirty,
                rebuild_if_dirty,
                update_axis_viewport_on_resize,
                sync_axis_camera_to_main,

                // measurements
                measurements::configure_measurement_gizmos,
                measurements::handle_measurement_picking, // keep if you haven't migrated to events
                measurements::update_measurement_distances,
                measurements::draw_measurement_lines,

                // hbonds
                draw_hydrogen_bonds_dashed,
                hbonds::configure_hbond_gizmos,
                hbonds::draw_hydrogen_bonds_dashed,

                // builder (no more local click handler)
                configure_builder_gizmos,
                molecule_builder::builder_ui::handle_builder_atom_picked,
                draw_builder_highlights,
            ).chain(),
        )
        .run();
}
