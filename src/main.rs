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
// --- added ---
mod measurements;
mod ui_measurements;
// --- end added ---

use camera::orbit_camera_system;
use scene::{rebuild_if_dirty, setup, sync_axis_camera_to_main, update_axis_viewport_on_resize};
use settings::MolSettings;
use ui::XyzBuffer;
use hbonds::{HbondGizmos, draw_hydrogen_bonds_dashed};

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(ImagePlugin::default_nearest()))
        .add_plugins(EguiPlugin::default())
        .add_plugins(export::ExportPlugin) 
        // resources
        .insert_resource(MolSettings::default())
        .insert_resource(camera::OrbitCamera::default())
        .insert_resource(molecule::Molecule::from_xyz(scene::DEFAULT_WATER))
        .insert_resource(XyzBuffer {
            text: scene::DEFAULT_WATER.trim().to_string(),
        })
        // ensure export UI resources exist before any egui systems run
        .init_resource::<export::ExportUiState>()
        .init_resource::<export::ExportRequestQueue>()
        // --- added: measurements resource & gizmo group ---
        .init_resource::<measurements::Measurements>()
        .init_gizmo_group::<measurements::MeasurementGizmos>()
        // --- end added ---
        // gizmos for dashed hydrogen bonds
        .init_gizmo_group::<HbondGizmos>()
        // scene
        .add_systems(Startup, setup)
        // egui must run in this pass on bevy_egui 0.36
        .add_systems(EguiPrimaryContextPass, ui::ui_panel)
        .add_systems(EguiPrimaryContextPass, ui_measurements::distance_labels_overlay)//to display
                                                                                      //measured
                                                                                      //distances
                                                                                      //on 3D
                                                                                      //structure
        .add_systems(
            Update,
            (
                orbit_camera_system,
                rebuild_if_dirty,
                update_axis_viewport_on_resize,
                sync_axis_camera_to_main,

                // --- added: measurements systems (order matters; before drawing gizmos) ---
                measurements::configure_measurement_gizmos,
                measurements::handle_measurement_picking,
                measurements::update_measurement_distances,
                measurements::draw_measurement_lines,
                // --- end added ---

                draw_hydrogen_bonds_dashed,
                hbonds::configure_hbond_gizmos,
                hbonds::draw_hydrogen_bonds_dashed,
            ).chain(),
        )
        .run();
}

