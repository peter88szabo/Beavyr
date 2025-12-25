use bevy::prelude::*;
use std::collections::HashMap;

use crate::color_schemes::color_scheme_map;

/// Color schemes exactly as in the single-file version
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ColorScheme {
    Custom,
    CPK,
    Jmol,
    VMD,
    Molden,
    Molden0,
}

#[allow(dead_code)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ShadingPreset {
    MoldenClassic,
    NeutralPbr,
    Glossy,
    Matte,
}

/// Lighting modes (old single-point or a classic three-point rig)
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LightingMode {
    SinglePoint,
    ThreePoint,
}

/// Bond color modes: a single uniform color or split by the colors of the two atoms.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BondColorMode {
    Uniform,
    AtomSplit,
}

#[derive(Resource, Clone)]
pub struct MolSettings {
    // geometry
    pub atom_scale: f32,        // multiplier on covalent radius for spheres
    pub atom_resolution: u32,   // atom resolution 
    pub bond_thresh_scale: f32, // neighbor cutoff factor
    pub hbond_cutoff: f32,      // neighbor cutoff factor
    pub bond_radius_pct: f32,   // fraction (0.05..1.0) of the smaller covalent radius

    // appearance
    pub scheme: ColorScheme,
    pub element_colors: HashMap<String, Color>,
    pub bg_color: Color,

    // material (shared by atoms and bonds)
    pub metallic: f32,
    pub roughness: f32,   // perceptual_roughness in Bevy
    pub reflectance: f32, // specular reflectance
    pub use_emissive: bool,
    pub emissive_color: Color,
    pub emissive_strength: f32,

    // lighting (ambient + key/fill/rim)
    pub lighting_mode: LightingMode,
    pub use_ambient: bool,
    pub ambient_color: Color,
    pub ambient_brightness: f32, // AmbientLight brightness

    // Key (this reuses your original sliders)
    pub light_intensity: f32,
    pub light_distance: f32,

    // Fill / Rim (only used in ThreePoint)
    pub fill_intensity: f32,
    pub fill_distance: f32,
    pub rim_intensity: f32,
    pub rim_distance: f32,

    // --- new: bond coloring options ---
    pub bond_color_mode: BondColorMode,
    pub uniform_bond_color: Color,

    // --- new: hydrogen-bond (gizmo) options ---
    pub show_hbonds: bool,
    pub hbond_color: Color,
    pub hbond_thickness: f32,
    pub hbond_gap_scale: f32,
    pub hbond_line_scale: f32,

    // dirty flag
    pub dirty: bool,
}

impl Default for MolSettings {
    fn default() -> Self {
        MolSettings {
            // requested defaults
            atom_scale: 0.6,
            atom_resolution: 7,
            bond_thresh_scale: 1.2,
            hbond_cutoff: 3.0,
            bond_radius_pct: 0.40, // 50% of smaller covalent radius

            scheme: ColorScheme::Jmol,
            element_colors: color_scheme_map(ColorScheme::Jmol),
            bg_color: Color::srgb(0.0, 0.0, 0.0),

            metallic: 0.05,
            roughness: 0.25,
            reflectance: 0.2,
            use_emissive: false,
            emissive_color: Color::srgb(1.0, 1.0, 1.0),
            emissive_strength: 0.0,

            // lighting mode + ambient
            lighting_mode: LightingMode::ThreePoint,
            use_ambient: false,
            ambient_color: Color::srgb(1.0, 1.0, 1.0),
            ambient_brightness: 0.0, // raise to “wash” the scene

            // key (original)
            light_intensity: 6_200_000.0,
            light_distance: 10.0,

            // three-point defaults (gentle)
            fill_intensity: 4_000_000.0,
            fill_distance: 10.0,
            rim_intensity: 4_000_000.0,
            rim_distance: 10.0,

            // --- new defaults ---
            bond_color_mode: BondColorMode::AtomSplit,
            uniform_bond_color: Color::srgb(0.85, 0.85, 0.88),

            // pleasant blue for H-bonds
            show_hbonds: true,
            hbond_color: Color::srgb(0.20, 0.60, 1.00),
            hbond_thickness: 30.0,
            hbond_gap_scale: 1.0,
            hbond_line_scale: 1.0,

            dirty: true,
        }
    }
}
