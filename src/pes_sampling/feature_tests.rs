//! Quasi-random draws, diverse subsets, reaction-path and curvilinear sampling.
use super::{
    curvilinear, diverse,
    input::{self, Fragment, BOHR},
    path,
    quasi::{scrambled_radical_inverse, Draws, Unit},
    sampling::{self, generate},
    state::*,
};
use crate::qchem_interfaces::hessian_io::{self, HessianData};
use ndarray::Array2;
use std::sync::atomic::{AtomicBool, AtomicUsize};

fn water() -> HessianData {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let xyz = hessian_io::geometry(
        &std::fs::read_to_string(root.join("pyscf_pbe_water_hessian.xyz")).unwrap(),
    )
    .unwrap();
    hessian_io::load(&root.join("pyscf_pbe_water_hessian.txt"), Some(&xyz)).unwrap()
}
fn water_fragment() -> Fragment {
    let h = water();
    input::from_hessian(h.atoms, h.masses_amu, h.coords_bohr, h.hessian, "water".into()).unwrap()
}
fn run(s: Settings) -> super::output::SamplingOutput {
    generate(s, 11, &AtomicBool::new(false), &AtomicUsize::new(0)).unwrap()
}
fn distance(x: &[f64], a: usize, b: usize) -> f64 {
    (0..3).map(|k| (x[3 * a + k] - x[3 * b + k]).powi(2)).sum::<f64>().sqrt()
}
fn chi2(points: &[[f64; 2]], bins: usize) -> f64 {
    let mut count = vec![0f64; bins * bins];
    for p in points {
        let (i, j) = ((p[0] * bins as f64) as usize, (p[1] * bins as f64) as usize);
        count[i * bins + j] += 1.;
    }
    let expected = points.len() as f64 / (bins * bins) as f64;
    count.iter().map(|c| (c - expected).powi(2) / expected).sum()
}
fn pairs(kind: RandomNumbers, seed: u64, dims: (usize, usize), n: usize) -> Vec<[f64; 2]> {
    let mut d = Draws::new(kind, seed);
    (0..n)
        .map(|_| {
            d.next_point();
            let v: Vec<f64> = (0..=dims.1).map(|_| d.unit()).collect();
            [v[dims.0], v[dims.1]]
        })
        .collect()
}

#[test]
fn quasi_random_is_default_repeatable_and_more_even_than_pseudo_random() {
    assert_eq!(Settings::default().random, RandomNumbers::QuasiRandom);
    // Identity permutation reproduces the plain van der Corput sequence.
    let identity: Vec<u32> = (0..2).collect();
    assert_eq!(scrambled_radical_inverse(6, 2, &identity), 0.375);
    let a = pairs(RandomNumbers::QuasiRandom, 5, (0, 1), 64);
    assert_eq!(a, pairs(RandomNumbers::QuasiRandom, 5, (0, 1), 64));
    assert_ne!(a, pairs(RandomNumbers::QuasiRandom, 6, (0, 1), 64));
    // Leading dimensions: 2-D projections are far more even (chi² for 64
    // bins is ~63 ± 11 for independent draws).
    for dims in [(0, 1), (2, 3), (4, 5)] {
        let q = pairs(RandomNumbers::QuasiRandom, 3, dims, 4096);
        let p = pairs(RandomNumbers::PseudoRandom, 3, dims, 4096);
        assert!(q.iter().flatten().all(|v| (0.0..1.0).contains(v)));
        let (cq, cp) = (chi2(&q, 8), chi2(&p, 8));
        assert!(cq < 0.5 * cp, "dims {dims:?}: quasi {cq:.1} pseudo {cp:.1}");
    }
    // High dimensions (bases ~180-660): each coordinate on its own is still
    // stratified, although pairs are only mildly better than independent
    // draws until N exceeds the product of their bases.
    for k in [40, 120] {
        let q = pairs(RandomNumbers::QuasiRandom, 3, (k, k), 4096);
        let p = pairs(RandomNumbers::PseudoRandom, 3, (k, k), 4096);
        let one = |v: &[[f64; 2]]| {
            let mut c = [0f64; 64];
            v.iter().for_each(|x| c[(x[0] * 64.) as usize] += 1.);
            c.iter().map(|n| (n - 64.).powi(2) / 64.).sum::<f64>()
        };
        assert!(one(&q) < 0.25 * one(&p), "dim {k}: {} vs {}", one(&q), one(&p));
    }
}

#[test]
fn farthest_point_keeps_extremes_and_drops_duplicates() {
    // Diatomic candidates: 40 copies at 1.5 Å plus a 1.0–2.0 Å ladder.
    let mut candidates = vec![vec![0., 0., 0., 1.5, 0., 0.]; 40];
    candidates.extend((0..=10).map(|i| vec![0., 0., 0., 1. + 0.1 * i as f64, 0., 0.]));
    let chosen = diverse::farthest_point(&candidates, 6, &AtomicBool::new(false), &AtomicUsize::new(0)).unwrap();
    let mut r: Vec<f64> = chosen.iter().map(|&i| candidates[i][3]).collect();
    r.sort_by(f64::total_cmp);
    r.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    assert_eq!(r.len(), 6, "a duplicate was selected: {r:?}");
    assert!(r.contains(&1.0) && r.contains(&2.0), "{r:?}");
}

#[test]
fn diverse_subset_is_default_keeps_counts_and_widens_coverage() {
    assert!(Settings::default().diversity.enabled);
    let mut s = Settings::default();
    s.partners[0].source = Source::Diatomic;
    s.partners[0].excitation.population = Population::Quantum;
    s.partners[0].excitation.quantum = 3;
    s.points = 30;
    s.random = RandomNumbers::PseudoRandom;
    // Smallest gap between sorted bond lengths: near-duplicates are removed.
    let gap = |s: Settings| {
        let out = run(s);
        assert_eq!(out.samples.len(), 30);
        let mut r: Vec<f64> = out.samples.iter().map(|x| distance(&x.coordinates, 0, 1)).collect();
        r.sort_by(f64::total_cmp);
        (r.windows(2).map(|w| w[1] - w[0]).fold(f64::INFINITY, f64::min), out)
    };
    let (diverse_gap, out) = gap(s.clone());
    assert!(out.sample_comment(0).contains("diverse_subset=farthest_point_inverse_distances:oversampling=3"));
    s.diversity.enabled = false;
    let (plain_gap, _) = gap(s.clone());
    assert!(diverse_gap > 3. * plain_gap, "{diverse_gap} vs {plain_gap}");
    // Bimolecular selection is per distance: the grid keeps its exact layout.
    s.diversity.enabled = true;
    s.workflow = Workflow::Bimolecular;
    s.partners[1].source = Source::Atom;
    s.partners[1].symbols[0] = "Ar".into();
    s.distances = 3;
    s.samples_per_distance = 4;
    let out = run(s.clone());
    let index: Vec<_> = out.samples.iter().map(|x| (x.distance_index.unwrap(), x.sample_at_distance.unwrap())).collect();
    assert_eq!(index, (1..=3).flat_map(|d| (1..=4).map(move |i| (d, i))).collect::<Vec<_>>());
    assert_eq!(sampling::work(&s), 3 * (12 + 4));
    // Oversampling outside 1–100 is rejected.
    s.diversity.oversampling = 0.5;
    assert!(sampling::validate(&s).is_err());
}

#[test]
fn internal_coordinates_reproduce_first_order_bond_and_angle_changes() {
    let f = water_fragment();
    let internal = curvilinear::prepare(&f).unwrap();
    let x0 = &f.equilibrium;
    // Large-amplitude bend (lowest mode) plus a stretch.
    let q = [60., 0., 15.];
    let dx: Vec<f64> = (0..9)
        .map(|i| (0..3).map(|j| f.modes[(i, j)] * q[j]).sum())
        .collect();
    let curved = curvilinear::displace(x0, &dx, &internal).unwrap();
    let straight: Vec<f64> = x0.iter().zip(&dx).map(|(a, b)| a + b).collect();
    for h in [1, 2] {
        let e: Vec<f64> = (0..3).map(|k| (x0[3 * h + k] - x0[k]) / distance(x0, 0, h)).collect();
        let first_order = distance(x0, 0, h) + (0..3).map(|k| e[k] * (dx[3 * h + k] - dx[k])).sum::<f64>();
        // Water has exactly 3 internals = 3N-6, so the target is met exactly.
        assert!((distance(&curved, 0, h) - first_order).abs() < 1e-8);
        // A straight-line bend stretches the bond at second order.
        assert!((distance(&straight, 0, h) - first_order).abs() > 1e-4);
    }
    // Linear and disconnected structures are refused with a clear message.
    let mut linear = f.clone();
    linear.equilibrium = vec![0., 0., -2.2, 0., 0., 0., 0., 0., 2.2];
    assert!(curvilinear::prepare(&linear).unwrap_err().contains("nonlinear"));
    let mut apart = f.clone();
    apart.equilibrium[6..9].copy_from_slice(&[30., 0., 0.]);
    assert!(curvilinear::prepare(&apart).unwrap_err().contains("Cartesian displacement"));
}

#[test]
fn internal_coordinate_option_runs_end_to_end_and_is_recorded() {
    assert_eq!(Settings::default().coordinates, Coordinates::Cartesian);
    let mut s = Settings::default();
    s.partners[0].fragment = Some(water_fragment());
    s.partners[0].excitation.population = Population::Thermal;
    s.partners[0].excitation.temperature = 2000.;
    s.coordinates = Coordinates::Internal;
    s.points = 40;
    let out = run(s);
    assert_eq!(out.samples.len(), 40);
    assert!(out.sample_comment(0).contains("coordinates=internal_curvilinear"));
    assert!(out.sample_comment(0).contains("random=scrambled_halton"));
    for x in &out.samples {
        for h in [1, 2] {
            let r = distance(&x.coordinates, 0, h);
            // n = 2 of the O–H stretch at 2000 K: turning points ~0.7 and ~1.25 Å.
            assert!((0.6..1.4).contains(&r), "O-H {r}");
        }
    }
}

fn rotation(a: f64, b: f64) -> [[f64; 3]; 3] {
    let (sa, ca) = a.sin_cos();
    let (sb, cb) = b.sin_cos();
    [[ca, -sa * cb, sa * sb], [sa, ca * cb, -ca * sb], [0., sb, cb]]
}
fn rotate(x: &[f64], r: [[f64; 3]; 3]) -> Vec<f64> {
    x.chunks_exact(3)
        .flat_map(|p| (0..3).map(move |i| (0..3).map(|k| r[i][k] * p[k]).sum::<f64>()))
        .collect()
}
/// Water with one O–H stretched step by step: a model reaction path.
fn water_path(rotated: bool) -> Vec<(String, HessianData)> {
    let h = water();
    (0..5)
        .map(|k| {
            let mut d = h.clone();
            let x = &mut d.coords_bohr;
            let e: Vec<f64> = (0..3).map(|i| (x[3 + i] - x[i]) / distance(x, 0, 1)).collect();
            for i in 0..3 {
                x[3 + i] += e[i] * 0.05 * k as f64 / BOHR;
            }
            if rotated {
                let r = rotation(0.7 * k as f64, 0.3 + 0.4 * k as f64);
                d.coords_bohr = rotate(&d.coords_bohr, r);
                let big = Array2::from_shape_fn((9, 9), |(i, j)| if i / 3 == j / 3 { r[i % 3][j % 3] } else { 0. });
                d.hessian = big.dot(&d.hessian).dot(&big.t());
            }
            (format!("irc_{}.hess", k + 1), d)
        })
        .collect()
}

#[test]
fn superposition_and_tangents_are_frame_independent() {
    let x = water().coords_bohr;
    let m = water().masses_amu;
    let mut centered = x.clone();
    input::center(&mut centered, &m);
    let turned = rotate(&centered, rotation(1.1, -0.4));
    let back = path::superpose(&turned, &centered, &m).unwrap();
    assert!(back.iter().zip(&centered).all(|(a, b)| (a - b).abs() < 1e-10));
    let plain = path::tangents(&water_path(false).into_iter().map(|p| p.1.coords_bohr).collect::<Vec<_>>(), &m).unwrap();
    let turned = path::tangents(&water_path(true).into_iter().map(|p| p.1.coords_bohr).collect::<Vec<_>>(), &m).unwrap();
    for (a, b) in plain.iter().zip(&turned) {
        let n = |t: &[f64]| t.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!((n(a) - n(b)).abs() < 1e-10);
    }
    assert!(path::tangents(&[x.clone()], &m).is_err());
    assert!(path::tangents(&[x.clone(), x], &m).unwrap_err().contains("same geometry"));
    let mut names = vec!["irc_10.hess", "irc_2.hess", "IRC_1.hess"];
    names.sort_by_key(|n| path::natural_key(n));
    assert_eq!(names, ["IRC_1.hess", "irc_2.hess", "irc_10.hess"]);
}

#[test]
fn reaction_path_modes_are_orthogonal_to_the_path_and_sampled_per_point() {
    let points = path::build(water_path(true)).unwrap();
    let m = water().masses_amu;
    let geometries: Vec<Vec<f64>> = water_path(true).into_iter().map(|p| p.1.coords_bohr).collect();
    let tangents = path::tangents(&geometries, &m).unwrap();
    for (f, t) in points.iter().zip(&tangents) {
        // 3N-6 = 3 vibrations minus the path direction.
        assert_eq!(f.frequencies.len() + f.imaginary_cm1.len(), 2);
        for j in 0..f.frequencies.len() {
            let overlap: f64 = (0..9).map(|i| m[i / 3] * f.modes[(i, j)] * t[i]).sum();
            let norm: f64 = (0..9).map(|i| m[i / 3] * f.modes[(i, j)].powi(2)).sum::<f64>().sqrt()
                * (0..9).map(|i| m[i / 3] * t[i] * t[i]).sum::<f64>().sqrt();
            assert!(overlap.abs() / norm < 1e-5, "mode {j} overlaps the path: {}", overlap / norm);
        }
    }
    // Without the path projection the stretch mode lies largely along the path.
    let h = &water_path(true)[2].1;
    let free = input::from_hessian(h.atoms.clone(), m.clone(), h.coords_bohr.clone(), h.hessian.clone(), "free".into()).unwrap();
    let t = &tangents[2];
    let largest = (0..free.frequencies.len()).map(|j| {
        let o: f64 = (0..9).map(|i| m[i / 3] * free.modes[(i, j)] * t[i]).sum();
        let n: f64 = (0..9).map(|i| m[i / 3] * free.modes[(i, j)].powi(2)).sum::<f64>().sqrt()
            * (0..9).map(|i| m[i / 3] * t[i] * t[i]).sum::<f64>().sqrt();
        (o / n).abs()
    }).fold(0., f64::max);
    assert!(largest > 0.5, "{largest}");
    let mut s = Settings::default();
    s.workflow = Workflow::ReactionPath;
    assert!(sampling::total(&s).is_err());
    s.path.points = points;
    s.path.samples_per_point = 3;
    let out = run(s);
    assert_eq!(out.samples.len(), 15);
    assert_eq!(out.samples.iter().map(|x| x.path_point.unwrap()).collect::<Vec<_>>(), (1..=5).flat_map(|k| [k; 3]).collect::<Vec<_>>());
    assert!(out.sample_comment(4).contains("path_point=2 path_points=5 path_file=irc_2.hess"));
    assert!(out.sample_comment(4).contains("path_projection=reaction_path_hamiltonian"));
}

#[test]
fn internal_coordinates_converge_for_a_real_19_atom_transition_state() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples/TS_Gamma-3-6_ZZ-S_11_Compound_1.hess");
    let f = input::load_hessian(&path, None).unwrap();
    let mut s = Settings::default();
    s.partners[0].fragment = Some(f.clone());
    s.partners[0].excitation.population = Population::Thermal;
    s.partners[0].excitation.temperature = 500.;
    s.coordinates = Coordinates::Internal;
    s.points = 20;
    let out = run(s.clone());
    assert_eq!(out.samples.len(), 20);
    // Bonds detected at equilibrium stay within a factor of their length.
    let x0: Vec<f64> = f.equilibrium.iter().map(|v| v * BOHR).collect();
    let bonded: Vec<(usize, usize)> = (0..19)
        .flat_map(|i| (0..i).map(move |j| (i, j)))
        .filter(|&(i, j)| distance(&x0, i, j) < 1.6)
        .collect();
    for x in &out.samples {
        for &(i, j) in &bonded {
            let ratio = distance(&x.coordinates, i, j) / distance(&x0, i, j);
            assert!((0.7..1.4).contains(&ratio), "bond {i}-{j} ratio {ratio}");
        }
    }
}
