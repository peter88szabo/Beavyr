use bevy::prelude::*;
use std::collections::HashMap;

use crate::settings::ColorScheme;

pub const ELEMENT_SYMBOLS: &[&str] = &[
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

const HALO_F: Color = Color::srgb(0.22, 1.0, 0.08); // neon green
const HALO_CL: Color = Color::srgb(0.70, 1.0, 0.05); // yellowish green
const HALO_BR: Color = Color::srgb(0.65, 0.20, 0.05); // rusty red-brown
const HALO_I: Color = Color::srgb(0.80, 0.20, 0.80); // purple-pink

fn base_map() -> HashMap<String, Color> {
    let mut m = HashMap::new();
    for &sym in ELEMENT_SYMBOLS {
        m.insert(sym.to_string(), Color::srgb(0.7, 0.7, 0.7));
    }
    m
}

/// Handy: get a full per-element map (used to seed `element_colors`)
pub fn color_scheme_map(s: ColorScheme) -> HashMap<String, Color> {
    use ColorScheme::*;
    let mut m = base_map();
    match s {
        // Simplified sets; tune as you like
        CPK => {
            m.insert("H".into(), Color::srgb(1.0, 1.0, 1.0));
            m.insert("C".into(), Color::srgb(0.2, 0.2, 0.2));
            m.insert("N".into(), Color::srgb(0.2, 0.2, 1.0));
            m.insert("O".into(), Color::srgb(1.0, 0.0, 0.0));
            m.insert("F".into(), HALO_F);
            m.insert("S".into(), Color::srgb(1.0, 1.0, 0.0));
            m.insert("Cl".into(), HALO_CL);
            m.insert("Br".into(), HALO_BR);
            m.insert("I".into(), HALO_I);
        }
        Jmol => {
            m.insert("H".into(), Color::srgb(1.0, 1.0, 1.0));
            m.insert("C".into(), Color::srgb(0.56, 0.56, 0.56));
            m.insert("N".into(), Color::srgb(0.19, 0.31, 0.97));
            m.insert("O".into(), Color::srgb(1.0, 0.05, 0.05));
            m.insert("F".into(), HALO_F);
            m.insert("S".into(), Color::srgb(1.0, 0.78, 0.18));
            m.insert("Cl".into(), HALO_CL);
            m.insert("Br".into(), HALO_BR);
            m.insert("I".into(), HALO_I);
        }
        VMD => {
            m.insert("H".into(), Color::srgb(0.85, 0.85, 0.85));
            m.insert("C".into(), Color::srgb(0.30, 0.30, 0.30));
            m.insert("N".into(), Color::srgb(0.10, 0.30, 0.90));
            m.insert("O".into(), Color::srgb(0.90, 0.10, 0.10));
            m.insert("F".into(), HALO_F);
            m.insert("S".into(), Color::srgb(1.0, 0.90, 0.10));
            m.insert("Cl".into(), HALO_CL);
            m.insert("Br".into(), HALO_BR);
            m.insert("I".into(), HALO_I);
        }
        Molden | Molden0 => {
            // Two similar variants; Molden0 slightly softer O red + brighter H
            let (o_r, h_w) = if matches!(s, Molden0) { (0.90, 0.95) } else { (1.0, 0.90) };
            m.insert("H".into(), Color::srgb(h_w, h_w, h_w));
            m.insert("C".into(), Color::srgb(0.20, 0.20, 0.20));
            m.insert("N".into(), Color::srgb(0.15, 0.25, 0.95));
            m.insert("O".into(), Color::srgb(o_r, 0.05, 0.05));
            m.insert("F".into(), HALO_F);
            m.insert("S".into(), Color::srgb(0.95, 0.90, 0.10));
            m.insert("Cl".into(), HALO_CL);
            m.insert("Br".into(), HALO_BR);
            m.insert("I".into(), HALO_I);
        }
        Custom => {
            // Neutral start
            m.insert("H".into(), Color::srgb(1.0, 1.0, 1.0));
            m.insert("C".into(), Color::srgb(0.2, 0.2, 0.2));
            m.insert("N".into(), Color::srgb(0.2, 0.2, 1.0));
            m.insert("O".into(), Color::srgb(1.0, 0.0, 0.0));
            m.insert("F".into(), HALO_F);
            m.insert("S".into(), Color::srgb(1.0, 1.0, 0.0));
            m.insert("Cl".into(), HALO_CL);
            m.insert("Br".into(), HALO_BR);
            m.insert("I".into(), HALO_I);
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
        "F"  => HALO_F,
        "S"  => Color::srgb(1.0, 1.0, 0.0),
        "Cl" => HALO_CL,
        "Br" => HALO_BR,
        "I"  => HALO_I,
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
        "F"  => HALO_F,
        "S"  => Color::srgb(1.0, 1.0, 0.2),
        "Cl" => HALO_CL,
        "Br" => HALO_BR,
        "I"  => HALO_I,
        _    => Color::srgb(0.8, 0.8, 0.85),
    }
}
fn molden0(a: &str) -> Color {
    match a {
        "H"  => Color::srgb(0.98, 0.98, 0.98),
        "C"  => Color::srgb(0.15, 0.15, 0.15),
        "N"  => Color::srgb(0.25, 0.45, 1.0),
        "O"  => Color::srgb(1.0, 0.25, 0.25),
        "F"  => HALO_F,
        "S"  => Color::srgb(1.0, 1.0, 0.25),
        "Cl" => HALO_CL,
        "Br" => HALO_BR,
        "I"  => HALO_I,
        _    => Color::srgb(0.86, 0.86, 0.9),
    }
}
