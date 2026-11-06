//! Geometry pairing, mode overview, histograms, and the clash rule.
use super::{
    analysis::{self, HistogramView, Quantity},
    input::{self, BOHR},
    overlap::OverlapFilter,
    sampling::generate,
    state::*,
};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}
fn water_partner() -> Partner {
    let mut p = Partner::default();
    let text = std::fs::read_to_string(fixtures().join("pyscf_pbe_water_hessian.xyz")).unwrap();
    p.geometry = Some(input::geometry(&text).unwrap());
    p.geometry_name = Some("water.xyz".into());
    let f = input::load_partner_hessian(&mut p, &fixtures().join("pyscf_pbe_water_hessian.txt")).unwrap();
    p.fragment = Some(f);
    p
}
fn run(s: Settings) -> super::output::SamplingOutput {
    generate(s, 5, &AtomicBool::new(false), &AtomicUsize::new(0)).unwrap()
}

#[test]
fn an_equilibrium_xyz_is_used_once_and_recorded() {
    let mut p = water_partner();
    assert_eq!(p.coordinates_source.as_deref(), Some("coordinates from water.xyz"));
    assert!(p.geometry.is_none() && p.geometry_name.is_none());
    // A second coordinate-free matrix must not silently reuse the old XYZ.
    let err = input::load_partner_hessian(&mut p, &fixtures().join("pyscf_pbe_water_hessian.txt"));
    assert!(err.is_err());
}

#[test]
fn mode_overview_bounds_the_sampled_displacements() {
    let mut p = water_partner();
    p.excitation.population = Population::Quantum;
    p.excitation.quantum = 2;
    let f = p.fragment.clone().unwrap();
    let rows = analysis::mode_overview(&f, &p);
    assert_eq!(rows.len(), 3);
    for (r, &w) in rows.iter().zip(&f.frequencies) {
        assert!((r.energy - 2.5 * w).abs() < 1e-12 && r.is_maximum);
        assert!(r.max_displacement > 0.05 && r.max_displacement < 0.6, "{r:?}");
    }
    // Uniform samples never exceed the sum of the per-mode turning points.
    let bound: f64 = rows.iter().map(|r| r.max_displacement).sum();
    let mut s = Settings::default();
    s.partners[0] = p.clone();
    s.points = 50;
    s.diversity.enabled = false;
    let out = run(s);
    let eq: Vec<f64> = f.equilibrium.iter().map(|v| v * BOHR).collect();
    let mut largest = 0f64;
    for x in &out.samples {
        // Remove rotation/translation-free comparison: water samples stay
        // centred and unrotated in the unimolecular workflow.
        for i in 0..3 {
            let d: f64 = (0..3).map(|k| (x.coordinates[3 * i + k] - eq[3 * i + k]).powi(2)).sum::<f64>().sqrt();
            largest = largest.max(d);
        }
    }
    assert!(largest <= bound + 1e-6, "{largest} > {bound}");
    // Smite thermal reports a mean, not a maximum.
    p.excitation.population = Population::Thermal;
    p.excitation.temperature_rule = TemperatureRule::SmiteThermal;
    assert!(analysis::mode_overview(&f, &p).iter().all(|r| !r.is_maximum));
}

#[test]
fn sample_energies_and_histograms_are_consistent() {
    let p = water_partner();
    let f = p.fragment.clone().unwrap();
    let zpe: f64 = f.frequencies.iter().map(|w| 0.5 * w).sum();
    let mut s = Settings::default();
    s.partners[0] = p;
    s.points = 60;
    let out = run(s);
    assert!(out.samples.iter().all(|x| x.energy >= 0. && x.energy <= zpe + 1e-12));
    assert!(out.samples.iter().any(|x| x.energy > 0.1 * zpe));
    let mut view = HistogramView::default();
    view.atoms = [1, 2, 3];
    let (d, unit) = analysis::values(&out, &view).unwrap();
    assert_eq!(unit, "Å");
    assert!(d.iter().all(|r| (0.8..1.2).contains(r)));
    view.quantity = Quantity::Angle;
    view.atoms = [2, 1, 3]; // H–O–H
    let (a, _) = analysis::values(&out, &view).unwrap();
    assert!(a.iter().all(|t| (80.0..130.).contains(t)));
    view.quantity = Quantity::Energy;
    let (e, unit) = analysis::values(&out, &view).unwrap();
    assert_eq!(unit, "kcal/mol");
    assert!((e[0] - out.samples[0].energy * 627.51).abs() < 1e-3 * e[0].max(1e-9));
    view.quantity = Quantity::Distance;
    view.atoms = [1, 1, 3];
    assert!(analysis::values(&out, &view).is_err());
    view.atoms = [1, 9, 3];
    assert!(analysis::values(&out, &view).is_err());
    let (lo, hi, counts) = analysis::histogram(&[0., 1., 1., 2.], 2);
    assert_eq!((lo, hi, counts), (0., 2., vec![1, 3]));
}

#[test]
fn covalent_rule_catches_heavy_atom_clashes() {
    let atoms = vec!["I".to_string(), "I".to_string()];
    let x = [0., 0., 0., 1.2, 0., 0.];
    let mut filter = OverlapFilter::default();
    assert!(filter.covalent);
    assert!(!filter.accepts(&atoms, &x));
    filter.covalent = false;
    assert!(filter.accepts(&atoms, &x), "fixed Pebble cutoff alone accepts 1.2 Å I–I");
    // A normal C–H bond is never rejected.
    let mut filter = OverlapFilter::default();
    filter.covalent_factor = 0.5;
    assert!(filter.accepts(&["C".into(), "H".into()], &[0., 0., 0., 1.09, 0., 0.]));
    filter.covalent_factor = 2.;
    assert!(filter.validate().is_err());
}
