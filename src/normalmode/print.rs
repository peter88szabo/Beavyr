#![allow(dead_code)]
use crate::normalmode::{EckartMode, HessianWeight};
use crate::normalmode::print_matrix_in_ao::print_ao_matrix;
use ndarray::Array2;

pub fn coord_labels(elements: &[String]) -> Vec<String> {
    let mut labels = Vec::with_capacity(elements.len() * 3);
    for (i, element) in elements.iter().enumerate() {
        let idx = i + 1;
        labels.push(format!("{}{}x", element, idx));
        labels.push(format!("{}{}y", element, idx));
        labels.push(format!("{}{}z", element, idx));
    }
    labels
}

pub fn print_hessian_matrix(title: &str, hess: &Array2<f64>, labels: &[String]) {
    let cols = std::env::var("HESS_COLS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|v| *v > 0)
        .unwrap_or(8);
    let val_width = 14;
    let val_prec = 8;
    print_ao_matrix(title, hess, labels, cols, val_width, val_prec)
        .expect("failed to print Hessian matrix");
}

pub fn format_hessian_title(base: &str, eckart: EckartMode, weight: HessianWeight) -> String {
    let eckart_label = match eckart {
        EckartMode::Off => "none",
        EckartMode::VibRot => "VibRot",
        EckartMode::ReactionPath => "ReactionPath",
    };
    let weight_label = match weight {
        HessianWeight::Cartesian => "Cartesian",
        HessianWeight::MassWeighted => "Mass-weighted",
    };
    format!("{base} (Eckart={eckart_label}, {weight_label})")
}
