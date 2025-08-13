use bevy::prelude::*;
use std::collections::HashMap;

use crate::settings::ColorScheme;

/// Handy: get a full per-element map (used to seed `element_colors`)
pub fn color_scheme_map(s: ColorScheme) -> HashMap<String, Color> {
    use ColorScheme::*;
    let mut m = HashMap::new();
    match s {
        // Simplified sets; tune as you like
        CPK => {
            m.insert("H".into(), Color::srgb(1.0, 1.0, 1.0));
            m.insert("C".into(), Color::srgb(0.2, 0.2, 0.2));
            m.insert("N".into(), Color::srgb(0.2, 0.2, 1.0));
            m.insert("O".into(), Color::srgb(1.0, 0.0, 0.0));
            m.insert("F".into(), Color::srgb(0.0, 1.0, 0.0));
            m.insert("S".into(), Color::srgb(1.0, 1.0, 0.0));
            m.insert("Cl".into(), Color::srgb(0.0, 0.7, 0.0));
        }
        Jmol => {
            m.insert("H".into(), Color::srgb(1.0, 1.0, 1.0));
            m.insert("C".into(), Color::srgb(0.56, 0.56, 0.56));
            m.insert("N".into(), Color::srgb(0.19, 0.31, 0.97));
            m.insert("O".into(), Color::srgb(1.0, 0.05, 0.05));
            m.insert("F".into(), Color::srgb(0.56, 0.88, 0.31));
            m.insert("S".into(), Color::srgb(1.0, 0.78, 0.18));
            m.insert("Cl".into(), Color::srgb(0.12, 0.94, 0.12));
        }
        VMD => {
            m.insert("H".into(), Color::srgb(0.85, 0.85, 0.85));
            m.insert("C".into(), Color::srgb(0.30, 0.30, 0.30));
            m.insert("N".into(), Color::srgb(0.10, 0.30, 0.90));
            m.insert("O".into(), Color::srgb(0.90, 0.10, 0.10));
            m.insert("F".into(), Color::srgb(0.15, 0.85, 0.15));
            m.insert("S".into(), Color::srgb(1.0, 0.90, 0.10));
            m.insert("Cl".into(), Color::srgb(0.05, 0.70, 0.05));
        }
        Molden | Molden0 => {
            // Two similar variants; Molden0 slightly softer O red + brighter H
            let (o_r, h_w) = if matches!(s, Molden0) { (0.90, 0.95) } else { (1.0, 0.90) };
            m.insert("H".into(), Color::srgb(h_w, h_w, h_w));
            m.insert("C".into(), Color::srgb(0.20, 0.20, 0.20));
            m.insert("N".into(), Color::srgb(0.15, 0.25, 0.95));
            m.insert("O".into(), Color::srgb(o_r, 0.05, 0.05));
            m.insert("F".into(), Color::srgb(0.10, 0.85, 0.10));
            m.insert("S".into(), Color::srgb(0.95, 0.90, 0.10));
            m.insert("Cl".into(), Color::srgb(0.00, 0.70, 0.00));
        }
        Custom => {
            // Neutral start
            m.insert("H".into(), Color::srgb(1.0, 1.0, 1.0));
            m.insert("C".into(), Color::srgb(0.2, 0.2, 0.2));
            m.insert("N".into(), Color::srgb(0.2, 0.2, 1.0));
            m.insert("O".into(), Color::srgb(1.0, 0.0, 0.0));
            m.insert("F".into(), Color::srgb(0.0, 1.0, 0.0));
            m.insert("S".into(), Color::srgb(1.0, 1.0, 0.0));
            m.insert("Cl".into(), Color::srgb(0.0, 0.7, 0.0));
        }
    }
    m
}

/// A convenience single-atom color query (mirrors the maps above)
pub fn color_for(atom: &str, scheme: ColorScheme) -> Color {
    match scheme {
        ColorScheme::Custom => default_custom(atom),
        ColorScheme::CPK    => cpk(atom),
        ColorScheme::Jmol   => jmol(atom),
        ColorScheme::VMD    => vmd(atom),
        ColorScheme::Molden => molden(atom),
        ColorScheme::Molden0=> molden0(atom),
    }
}

fn default_custom(a: &str) -> Color {
    match a {
        "H"  => Color::srgb(1.0, 1.0, 1.0),
        "C"  => Color::srgb(0.2, 0.2, 0.2),
        "N"  => Color::srgb(0.2, 0.2, 1.0),
        "O"  => Color::srgb(1.0, 0.0, 0.0),
        "F"  => Color::srgb(0.0, 1.0, 0.0),
        "S"  => Color::srgb(1.0, 1.0, 0.0),
        "Cl" => Color::srgb(0.0, 0.7, 0.0),
        _    => Color::srgb(0.7, 0.7, 0.7),
    }
}
fn cpk(a: &str) -> Color { default_custom(a) }
fn jmol(a: &str) -> Color { default_custom(a) }
fn vmd(a: &str)  -> Color { default_custom(a) }

fn molden(a: &str) -> Color {
    match a {
        "H"  => Color::srgb(0.95, 0.95, 0.95),
        "C"  => Color::srgb(0.12, 0.12, 0.12),
        "N"  => Color::srgb(0.2, 0.35, 1.0),
        "O"  => Color::srgb(1.0, 0.2, 0.2),
        "F"  => Color::srgb(0.3, 1.0, 0.3),
        "S"  => Color::srgb(1.0, 1.0, 0.2),
        "Cl" => Color::srgb(0.0, 0.8, 0.0),
        _    => Color::srgb(0.8, 0.8, 0.85),
    }
}
fn molden0(a: &str) -> Color {
    match a {
        "H"  => Color::srgb(0.98, 0.98, 0.98),
        "C"  => Color::srgb(0.15, 0.15, 0.15),
        "N"  => Color::srgb(0.25, 0.45, 1.0),
        "O"  => Color::srgb(1.0, 0.25, 0.25),
        "F"  => Color::srgb(0.35, 1.0, 0.35),
        "S"  => Color::srgb(1.0, 1.0, 0.25),
        "Cl" => Color::srgb(0.0, 0.85, 0.0),
        _    => Color::srgb(0.86, 0.86, 0.9),
    }
}

