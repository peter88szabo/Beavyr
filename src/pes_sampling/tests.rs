use super::{input::*, output::*, sampling::*, state::*};
use crate::optimizer::prng::StdRng;
use std::sync::atomic::{AtomicBool, AtomicUsize};

#[test]
fn pes_window_retains_user_resize_in_both_dimensions() {
    use bevy_egui::egui;
    let ctx = egui::Context::default();
    let mut state = PesSampling::default();
    state.open = true;
    fn frame(ctx: &egui::Context, state: &mut PesSampling, events: Vec<egui::Event>) -> egui::Rect {
        let mut rect = None;
        let mut output = ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1500., 1200.),
                )),
                events,
                ..Default::default()
            },
            |_| {
                rect = super::ui::window(ctx, state, &crate::molecule::Molecule::empty());
            },
        );
        output.textures_delta.clear();
        rect.unwrap()
    }
    frame(&ctx, &mut state, vec![]);
    let initial = frame(&ctx, &mut state, vec![]);
    let corner = initial.right_bottom() - egui::vec2(2., 2.);
    frame(&ctx, &mut state, vec![egui::Event::PointerMoved(corner)]);
    frame(
        &ctx,
        &mut state,
        vec![egui::Event::PointerButton {
            pos: corner,
            button: egui::PointerButton::Primary,
            pressed: true,
            modifiers: egui::Modifiers::NONE,
        }],
    );
    let target = corner + egui::vec2(140., 100.);
    frame(&ctx, &mut state, vec![egui::Event::PointerMoved(target)]);
    frame(
        &ctx,
        &mut state,
        vec![egui::Event::PointerButton {
            pos: target,
            button: egui::PointerButton::Primary,
            pressed: false,
            modifiers: egui::Modifiers::NONE,
        }],
    );
    let enlarged = frame(&ctx, &mut state, vec![]);
    assert!(
        enlarged.width() > initial.width() + 100.,
        "{initial:?} -> {enlarged:?}"
    );
    assert!(
        enlarged.height() > initial.height() + 60.,
        "{initial:?} -> {enlarged:?}"
    );
}

#[test]
fn pes_pebble_overlap_thresholds_and_pair_order() {
    let mut filter = super::overlap::OverlapFilter::default();
    assert!(filter.enabled);
    for (a, b, cutoff) in [
        ("H", "H", 0.4),
        ("H", "C", 0.5),
        ("O", "H", 0.5),
        ("C", "O", 0.8),
        ("Fe", "Br", 0.8),
    ] {
        let atoms = vec![a.into(), b.into()];
        assert!(!filter.accepts(&atoms, &[0., 0., 0., cutoff - 1e-6, 0., 0.]));
        assert!(filter.accepts(&atoms, &[0., 0., 0., cutoff, 0., 0.]));
        assert!(filter.accepts(&atoms, &[0., 0., 0., cutoff + 1e-6, 0., 0.]));
    }
    filter.enabled = false;
    assert!(filter.accepts(&vec!["H".into(), "O".into()], &[0.; 6]));
    filter.hydrogen_other = f64::NAN;
    assert!(filter.validate().is_err());
}

#[test]
fn pes_default_overlap_filter_redraws_unimolecular_samples() {
    let mut s = Settings::default();
    s.partners[0] = diatomic();
    s.partners[0].excitation.population = Population::Quantum;
    s.partners[0].excitation.quantum = 8;
    s.points = 100;
    let filtered = run(s.clone());
    assert_eq!(filtered.samples.len(), 100);
    assert!(filtered.rejected > 0);
    assert!(filtered
        .samples
        .iter()
        .all(|p| s.overlap.accepts(&filtered.atoms, &p.coordinates)));
    s.overlap.enabled = false;
    let unfiltered = run(s);
    let defaults = super::overlap::OverlapFilter::default();
    assert_eq!(unfiltered.rejected, 0);
    assert!(unfiltered
        .samples
        .iter()
        .any(|p| !defaults.accepts(&unfiltered.atoms, &p.coordinates)));
}

#[test]
fn pes_default_overlap_filter_checks_between_collision_partners() {
    let mut s = Settings::default();
    s.workflow = Workflow::Bimolecular;
    s.distances = 1;
    s.samples_per_distance = 3;
    for p in &mut s.partners {
        p.source = Source::Atom;
        p.symbols[0] = "C".into();
    }
    s.distance_start = 0.79;
    assert!(
        generate(s.clone(), 17, &AtomicBool::new(false), &AtomicUsize::new(0))
            .unwrap_err()
            .contains("1000 attempts")
    );
    s.partners[0].symbols[0] = "H".into();
    // The same separation is accepted for H-C, whose cutoff is 0.5 Å.
    assert_eq!(run(s.clone()).samples.len(), 3);
    s.distance_start = 0.49;
    assert!(
        generate(s.clone(), 17, &AtomicBool::new(false), &AtomicUsize::new(0))
            .unwrap_err()
            .contains("1000 attempts")
    );
    s.overlap.enabled = false;
    assert_eq!(run(s).samples.len(), 3);
}

fn diatomic() -> Partner {
    Partner {
        source: Source::Diatomic,
        ..Default::default()
    }
}
fn run(s: Settings) -> SamplingOutput {
    generate(s, 17, &AtomicBool::new(false), &AtomicUsize::new(0)).unwrap()
}
fn com(x: &[f64], m: &[f64]) -> [f64; 3] {
    let mut c = [0.; 3];
    let total: f64 = m.iter().sum();
    for (p, w) in x.chunks_exact(3).zip(m) {
        for k in 0..3 {
            c[k] += p[k] * w / total;
        }
    }
    c
}
#[test]
fn pes_uniform_phase_has_arcsine_displacements_and_conserves_mode_energy() {
    let w: f64 = 0.01;
    let e = 1.5 * w;
    let amplitude = (2. * e).sqrt() / w;
    let mut rng = StdRng::seed_from_u64(2026);
    let mut mean = 0.;
    let mut square = 0.;
    for _ in 0..100000 {
        let (q, v) = harmonic_state(w, e, std::f64::consts::TAU * uniform(&mut rng));
        assert!(q.abs() <= amplitude);
        assert!((0.5 * (v * v + w * w * q * q) - e).abs() < 1e-14);
        mean += q / amplitude;
        square += (q / amplitude).powi(2);
    }
    assert!(mean.abs() / 100000. < 0.005);
    assert!((square / 100000. - 0.5).abs() < 0.005);
}
#[test]
fn pes_temperature_rules_have_distinct_correct_limits() {
    let w = 0.001;
    let t = 1000.;
    let n = tail_quantum(w, t, 99.9);
    let ratio = (-w / (KB * t)).exp();
    assert!(1. - ratio.powf(n + 1.) >= 0.999 - 1e-14);
    assert!(n == 0. || 1. - ratio.powf(n) < 0.999);
    let mut e = Excitation {
        population: Population::Thermal,
        temperature: t,
        ..Default::default()
    };
    let mut rng = StdRng::seed_from_u64(17);
    assert_eq!(energy(&e, w, &mut rng), w * (n + 0.5));
    e.temperature_rule = TemperatureRule::ZpePlusKt;
    assert_eq!(energy(&e, w, &mut rng), 0.5 * w + KB * t);
    e.temperature_rule = TemperatureRule::SmiteThermal;
    let mut sum = 0.;
    let mut over = 0;
    for _ in 0..100000 {
        let quantum = energy(&e, w, &mut rng) / w - 0.5;
        sum += quantum;
        if quantum > n + 0.1 {
            over += 1;
        }
    }
    assert!((sum / 100000. - ratio / (1. - ratio)).abs() < 0.05);
    assert!(over > 20);
    assert_eq!(tail_quantum(w, 0., 99.9), 0.);
}
#[test]
fn pes_euler_rotations_are_isotropic_and_preserve_geometry() {
    let mut rng = StdRng::seed_from_u64(18);
    let mut sums = [0.; 3];
    let mut square = [0.; 3];
    for _ in 0..50000 {
        let mut x = [0., 0., 1.];
        rotate(
            &mut x,
            [uniform(&mut rng), uniform(&mut rng), uniform(&mut rng)],
        );
        assert!((x.iter().map(|v| v * v).sum::<f64>() - 1.).abs() < 1e-12);
        for k in 0..3 {
            sums[k] += x[k];
            square[k] += x[k] * x[k];
        }
    }
    for k in 0..3 {
        assert!(sums[k].abs() / 50000. < 0.01);
        assert!((square[k] / 50000. - 1. / 3.).abs() < 0.01);
    }
}
#[test]
fn pes_bimolecular_count_order_exact_com_and_independent_rotations() {
    let mut s = Settings::default();
    s.workflow = Workflow::Bimolecular;
    s.partners = [diatomic(), diatomic()];
    s.partners[1].symbols = ["H".into(), "Cl".into()];
    s.partners[1].frequency = 3000.;
    s.partners[1].equilibrium = 1.27;
    s.distance_start = 6.;
    s.distance_end = 2.;
    s.distances = 5;
    s.samples_per_distance = 13;
    let out = run(s.clone());
    assert_eq!(out.samples.len(), 65);
    let ma = masses(&out.atoms[..2]).unwrap();
    let mb = masses(&out.atoms[2..]).unwrap();
    let mut axis_correlation = 0.;
    for (i, f) in out.samples.iter().enumerate() {
        let a = com(&f.coordinates[..6], &ma);
        let b = com(&f.coordinates[6..], &mb);
        assert!(a[0].abs() < 1e-12 && b[1].abs() < 1e-12);
        let d = 6. - (i / 13) as f64;
        assert!((b[2] - a[2] - d).abs() < 1e-12);
        assert_eq!(f.distance, Some(d));
        let mut all = ma.clone();
        all.extend(&mb);
        assert!(com(&f.coordinates, &all).iter().all(|v| v.abs() < 1e-12));
        let x = &f.coordinates;
        let ua = [x[0] - x[3], x[1] - x[4], x[2] - x[5]];
        let ub = [x[6] - x[9], x[7] - x[10], x[8] - x[11]];
        let dot: f64 = (0..3).map(|k| ua[k] * ub[k]).sum();
        axis_correlation += dot
            / (ua.iter().map(|v| v * v).sum::<f64>() * ub.iter().map(|v| v * v).sum::<f64>())
                .sqrt();
    }
    assert!(axis_correlation.abs() / 65. < 0.4);
    let again = run(s);
    assert_eq!(out.samples[12].coordinates, again.samples[12].coordinates);
    let mut xyz = vec![];
    out.write_xyz(&mut xyz).unwrap();
    let text = String::from_utf8(xyz).unwrap();
    let frames = crate::trajectory::parse_multi_xyz(&text).unwrap();
    assert_eq!(frames.len(), 65);
    assert!(text.contains("COM_distance_A=2.0000000000 distance_index=5"));
}
#[test]
fn pes_diatomic_modes_mass_normalization_amplitude_and_morse_bounds() {
    let mut p = diatomic();
    p.symbols = ["H".into(), "Cl".into()];
    let f = prepare(&p).unwrap();
    let norm: f64 = (0..6)
        .map(|i| f.modes[(i, 0)].powi(2) * f.masses[i / 3] * AMU)
        .sum();
    assert!((norm - 1.).abs() < 1e-12);
    let mut rng = StdRng::seed_from_u64(1);
    let x = sample_fragment(&f, &p, &mut rng).unwrap();
    assert!(com(&x, &f.masses).iter().all(|v| v.abs() < 1e-12));
    let re = 1.4;
    let mu: f64 = 900.;
    let w = 0.02;
    let de = 0.2;
    let en = 0.05;
    let beta = w * (mu / (2. * de)).sqrt();
    for i in 0..1000 {
        let r =
            morse_distance(re, mu, w, de, en, i as f64 / 1000. * std::f64::consts::TAU).unwrap();
        let potential = de * (1. - (-beta * (r - re)).exp()).powi(2);
        assert!(potential <= en + 1e-12);
    }
    assert!(morse_distance(re, mu, w, de, de, 0.).is_err());
    p.morse = true;
    p.excitation.population = Population::Quantum;
    p.excitation.quantum = 100000;
    let mut s = Settings::default();
    s.partners[0] = p;
    assert!(validate(&s).is_err());
}
#[test]
fn pes_real_example_hessian_retains_ts_direction_separately_from_stable_modes() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples/TS_Gamma-3-6_ZZ-S_11_Compound_1.hess");
    let f = load_hessian(&path, None).unwrap();
    assert_eq!(f.atoms.len(), 19);
    assert_eq!(f.imaginary_cm1.len(), 1);
    assert_eq!(f.frequencies.len(), 50);
    for j in 0..f.modes.ncols() {
        let norm: f64 = (0..f.modes.nrows())
            .map(|i| f.modes[(i, j)].powi(2) * f.masses[i / 3] * AMU)
            .sum();
        assert!((norm - 1.).abs() < 1e-8);
    }
    let mut s = Settings::default();
    s.points = 20;
    s.partners[0].fragment = Some(f);
    let o = run(s);
    assert_eq!(o.samples.len(), 20);
    assert!(!o.warnings.is_empty());
}
#[test]
fn pes_cancellation_bad_input_and_impossible_filter_are_reported() {
    let mut s = Settings::default();
    s.partners[0] = diatomic();
    assert!(
        generate(s.clone(), 1, &AtomicBool::new(true), &AtomicUsize::new(0))
            .err()
            .unwrap()
            .starts_with("Cancelled")
    );
    s.overlap.hydrogen_hydrogen = 100.;
    assert!(
        generate(s.clone(), 1, &AtomicBool::new(false), &AtomicUsize::new(0))
            .err()
            .unwrap()
            .contains("1000 attempts")
    );
    s.overlap.enabled = false;
    s.points = 0;
    assert!(validate(&s).is_err());
    assert!(geometry("2\nx\nH 0 0 0\nH NaN 0 1\n").is_err());
    let (a, x) = geometry("3\nx\nX 0 0 0\nH 0 0 1\nH 0 0 2\n").unwrap();
    assert_eq!(a.len(), 2);
    assert_eq!(x.len(), 6);
}
#[test]
fn pes_result_loads_automatically_stopped_and_can_be_saved() {
    let mut s = Settings::default();
    s.points = 7;
    s.partners[0] = diatomic();
    let o = run(s);
    let mut state = PesSampling::default();
    state.output = Some(o.clone());
    state.load_requested = true;
    let mut mol = crate::molecule::Molecule::empty();
    let mut traj = crate::trajectory::TrajectoryState::default();
    traj.playing = true;
    assert!(super::ui::load_if_requested(
        &mut state, &mut mol, &mut traj
    ));
    assert_eq!(traj.frames.len(), 7);
    assert!(!traj.playing);
    assert_eq!(traj.mode, crate::trajectory::PlaybackMode::Once);
    assert_eq!(mol.atoms, vec!["H", "H"]);
    assert!(!state.load_requested);
    let path = std::env::temp_dir().join(format!("beavyr-pes-export-{}.xyz", std::process::id()));
    o.save(&path).unwrap();
    assert_eq!(
        crate::trajectory::parse_multi_xyz(&std::fs::read_to_string(path).unwrap())
            .unwrap()
            .len(),
        7
    );
}
#[test]
fn pes_ui_labels_and_old_layout_serialization_stay_compatible() {
    let ctx = bevy_egui::egui::Context::default();
    let mut state = PesSampling::default();
    state.open = true;
    state.settings.partners[0] = diatomic();
    state.settings.partners[0].excitation.population = Population::Thermal;
    let mut labels = Vec::new();
    for _ in 0..2 {
        let mut frame = ctx.run_ui(Default::default(), |_| {
            super::ui::window(&ctx, &mut state, &crate::molecule::Molecule::empty());
        });
        frame.textures_delta.clear();
        labels = frame
            .shapes
            .into_iter()
            .filter_map(|s| {
                if let bevy_egui::egui::Shape::Text(t) = s.shape {
                    Some(t.galley.text().to_owned())
                } else {
                    None
                }
            })
            .collect();
    }
    for expected in [
        "PES Sampling",
        "Unimolecular",
        "Bimolecular",
        "Harmonic",
        "Morse",
        "Boltzmann-tail cutoff",
        "Retained population",
    ] {
        assert!(
            labels.iter().any(|s| s == expected),
            "Missing control: {expected}"
        );
    }
    let layout:crate::ui_layout::UiLayout=ron::from_str("(windowed:true,open:(false,false,false,false,false,false,false,false,false,false,false,false,false,false,false),builder_open:false)").unwrap();
    assert_eq!(layout.open.len(), 15);
    assert_eq!(
        state.settings.partners[0].excitation.temperature_rule,
        TemperatureRule::BoltzmannTail
    );
}
#[test]
fn pes_matches_mature_smite_fixed_random_draws() {
    // Values obtained by running the original Smite functions with patched
    // random.uniform (see /tmp/beavyr-pes-smite-reference.py).
    let mu: f64 = 1837. * 3674. / (1837. + 3674.);
    let (q, v) = harmonic_state(0.014, 0.014 * 2.5, std::f64::consts::PI / 3.);
    let r = 2. + q / mu.sqrt();
    assert!((r * 2. / 3. - 1.513340908211463).abs() < 1e-13);
    assert!((mu.sqrt() * v + 8.018416302487667).abs() < 1e-12);
    let r = morse_distance(1.4, 900., 0.02, 0.2, 0.05, std::f64::consts::TAU * 0.37).unwrap();
    assert!((r - 2.03083079335721).abs() < 1e-13);
    let mut x = [1., 2., 3., -2., 1., 0.5];
    rotate(&mut x, [0.13, 0.87, 0.53]);
    let expected = [
        2.6422509100666933,
        -0.3299083362896216,
        2.628625233443966,
        1.7631871254540965,
        0.9901195568603351,
        -1.0774202632935825,
    ];
    for (a, b) in x.iter().zip(expected) {
        assert!((a - b).abs() < 1e-13);
    }
    assert_eq!(thermal_quantum(0.001, 1000., 0.873), 6.);
}
#[test]
fn pes_linear_hessian_retains_one_mode_and_respects_isotope_mass() {
    let p = diatomic();
    let f = prepare(&p).unwrap();
    let mu = f.masses[0] * f.masses[1] / (f.masses[0] + f.masses[1]) * AMU;
    let k = mu * f.frequencies[0].powi(2);
    let mut h = ndarray::Array2::zeros((6, 6));
    h[(0, 0)] = k;
    h[(3, 3)] = k;
    h[(0, 3)] = -k;
    h[(3, 0)] = -k;
    let light = from_hessian(
        f.atoms.clone(),
        f.masses.clone(),
        f.equilibrium.clone(),
        h.clone(),
        "test".into(),
    )
    .unwrap();
    assert_eq!(light.frequencies.len(), 1);
    assert!((light.frequencies[0] - f.frequencies[0]).abs() < 1e-10);
    let heavy = from_hessian(
        f.atoms,
        vec![2. * f.masses[0], 2. * f.masses[1]],
        f.equilibrium,
        h,
        "isotope".into(),
    )
    .unwrap();
    assert!((heavy.frequencies[0] / light.frequencies[0] - 1. / 2_f64.sqrt()).abs() < 1e-10);
}
#[test]
fn pes_uniform_displacement_covers_the_well_evenly_and_is_default() {
    assert_eq!(Settings::default().displacement, Displacement::Uniform);
    let (w, e): (f64, f64) = (0.01, 0.005);
    let a = (2. * e).sqrt() / w;
    let mut rng = StdRng::seed_from_u64(3);
    let mut bins = [0usize; 10];
    let mut classical_centre = 0;
    for _ in 0..20000 {
        let q = harmonic_displacement(w, e, Displacement::Uniform, uniform(&mut rng));
        assert!(q.abs() <= a);
        bins[(((q / a + 1.) / 2. * 10.) as usize).min(9)] += 1;
        let c = harmonic_displacement(w, e, Displacement::ClassicalPhase, uniform(&mut rng));
        if c.abs() < 0.2 * a {
            classical_centre += 1;
        }
    }
    // Uniform: every tenth of [-A, A] gets ~2000; classical arcsine gives the
    // central fifth only 2/π·asin(0.2) ≈ 12.8 %, versus 20 % here.
    assert!(bins.iter().all(|n| (1800..2200).contains(n)), "{bins:?}");
    assert!((2300..2850).contains(&classical_centre), "{classical_centre}");
    // Morse: uniform between the exact classical turning points.
    let (re, mu, de, en) = (1.4, 900., 0.2, 0.05);
    let inner = morse_displacement(re, mu, w, de, en, Displacement::Uniform, 0.).unwrap();
    let outer = morse_displacement(re, mu, w, de, en, Displacement::Uniform, 1.).unwrap();
    let beta = w * (mu / (2. * de)).sqrt();
    for r in [inner, outer] {
        assert!((de * (1. - (-beta * (r - re)).exp()).powi(2) - en).abs() < 1e-12);
    }
    assert!(inner < re && outer > re);
    let mid = morse_displacement(re, mu, w, de, en, Displacement::Uniform, 0.5).unwrap();
    assert!((mid - 0.5 * (inner + outer)).abs() < 1e-12);
}
#[test]
fn pes_mode_limits_override_global_values() {
    let mut p = diatomic();
    p.excitation.population = Population::Quantum;
    p.excitation.quantum = 5;
    let f = prepare(&p).unwrap();
    let mut a = StdRng::seed_from_u64(100);
    let excited = sample_fragment(&f, &p, &mut a).unwrap();
    p.use_overrides = true;
    p.overrides = vec![Excitation::default()];
    let mut b = StdRng::seed_from_u64(100);
    let ground = sample_fragment(&f, &p, &mut b).unwrap();
    let de = excited[0] - excited[3] - (f.equilibrium[0] - f.equilibrium[3]);
    let dg = ground[0] - ground[3] - (f.equilibrium[0] - f.equilibrium[3]);
    assert!((de / dg - 11_f64.sqrt()).abs() < 1e-10);
}

#[test]
fn hessian_export_reloads_for_both_pes_partners_and_frequency_analysis() {
    use crate::qchem_interfaces::hessian_io;
    use super::{input, state::{PesSampling, Workflow}};
    use crate::qchem_interfaces::xtb_freq::{load_hessian_from_path, XtbFreqPanelState, XtbFrequencyTask};
    let mut h = ndarray::Array2::zeros((6,6));
    h[(0,0)] = 0.3; h[(3,3)] = 0.3;
    h[(0,3)] = -0.3; h[(3,0)] = -0.3 + 2e-6;
    let data = hessian_io::from_geometry(vec!["H".into(), "H".into()], vec![0.,0.,0.,1.4,0.,0.], h).unwrap();
    let path = std::env::temp_dir().join(format!("beavyr-hessian-roundtrip-{}.hess", std::process::id()));
    hessian_io::save(&path, &data).unwrap();
    let mut state = PesSampling::default();
    state.settings.workflow = Workflow::Bimolecular;
    for partner in &mut state.settings.partners {
        partner.fragment = Some(input::load_hessian(&path, None).unwrap());
        assert_eq!(input::prepare(partner).unwrap().frequencies.len(), 1);
    }
    let fragment = input::load_hessian(&path, None).unwrap();
    let mut panel = XtbFreqPanelState::default();
    let mut task = XtbFrequencyTask::default();
    assert!(load_hessian_from_path(&path, &mut panel, &mut task), "{:?}", task.last_message);
    let frequency = fragment.frequencies[0] / input::CM_TO_AU;
    let result = task.result.unwrap();
    assert!(result.frequencies_cm1.iter().any(|f| (f-frequency).abs() < 1e-6));
    std::fs::remove_file(path).unwrap();
}
