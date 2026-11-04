use super::*;
use crate::optimizer::prng::StdRng;
use crate::pes_sampling::{sampling, state::*};
use ndarray::Array2;

fn ethane() -> Fragment {
    let atoms = ["C", "C", "H", "H", "H", "H", "H", "H"]
        .map(String::from)
        .to_vec();
    let x = [
        0., 0., 0., 1.54, 0., 0., -0.36, 1., 0., -0.36, -0.5, 0.866, -0.36, -0.5, -0.866, 1.90, 1.,
        0., 1.90, -0.5, 0.866, 1.90, -0.5, -0.866,
    ];
    Fragment {
        masses: super::super::input::masses(&atoms).unwrap(),
        atoms,
        equilibrium: x.iter().map(|v| v / BOHR).collect(),
        frequencies: vec![0.01],
        modes: Array2::zeros((24, 1)),
        imaginary_cm1: vec![],
        imaginary_modes: Array2::zeros((24, 0)),
        name: "ethane".into(),
        diatomic: None,
    }
}
fn rotor() -> Rotor {
    Rotor {
        atoms: [3, 1, 2, 6],
        side: FragmentSide::Second,
        ..Default::default()
    }
}
fn distance(x: &[f64], i: usize, j: usize) -> f64 {
    let d = sub(point(x, i), point(x, j));
    dot(d, d).sqrt()
}

#[test]
fn rotor_fragment_rotation_preserves_bonds_axis_and_spectator() {
    let f = ethane();
    let r = prepare(&f, &rotor()).unwrap();
    assert_eq!(r.moving, vec![5, 6, 7]);
    let mut x = f.equilibrium.clone();
    rotate(&mut x, &r, std::f64::consts::FRAC_PI_2).unwrap();
    assert_eq!(&x[..15], &f.equilibrium[..15]);
    for i in [1, 5, 6, 7] {
        for j in [1, 5, 6, 7] {
            assert!((distance(&x, i, j) - distance(&f.equilibrium, i, j)).abs() < 1e-12);
        }
    }
    assert!(x[17] > 1.8); // y=1 Å rotates into +z about +x.
    let mut manual = rotor();
    manual.side = FragmentSide::Explicit;
    manual.atom_list = "2,6-8".into();
    assert_eq!(prepare(&f, &manual).unwrap().moving, r.moving);
    manual.atom_list = "6,7".into();
    assert!(prepare(&f, &manual)
        .unwrap_err()
        .contains("complete fragment"));
}

#[test]
fn rotor_ring_invalid_atoms_duplicate_axes_and_modes_are_rejected() {
    let mut f = ethane();
    let r = rotor();
    let mut invalid = r.clone();
    invalid.atoms[2] = invalid.atoms[1];
    assert!(prepare(&f, &invalid).is_err());
    assert!(prepare_all(&f, &[r.clone(), r.clone()])
        .unwrap_err()
        .contains("more than once"));
    // Add a second heavy-atom path around the chosen bond.
    f.atoms[2] = "C".into();
    f.equilibrium[6..9].copy_from_slice(&[0.77 / BOHR, 1. / BOHR, 0.]);
    let mut ring = r;
    ring.dihedral = false;
    assert!(prepare(&f, &ring).unwrap_err().contains("ring"));
}

#[test]
fn rotor_angles_are_uniform_and_assigned_mode_can_be_replaced() {
    let f = ethane();
    let mut partner = Partner::default();
    partner.fragment = Some(f.clone());
    partner.rotors = vec![rotor()];
    let mut rng = StdRng::seed_from_u64(17);
    let mut bins = [0usize; 8];
    let mut sin_sum = 0.;
    let mut cos_sum = 0.;
    for _ in 0..4000 {
        let x = sampling::sample_fragment(&f, &partner, &mut rng).unwrap();
        let v = sub(point(&x, 5), point(&x, 1));
        let angle = v[2].atan2(v[1]).rem_euclid(std::f64::consts::TAU);
        bins[(angle / std::f64::consts::TAU * 8.).floor() as usize] += 1;
        sin_sum += angle.sin();
        cos_sum += angle.cos();
    }
    assert!(bins.iter().all(|n| (400..600).contains(n)), "{bins:?}");
    assert!(sin_sum.abs() / 4000. < 0.05 && cos_sum.abs() / 4000. < 0.05);
    let mut f = f;
    f.modes[(0, 0)] = 0.01;
    f.modes[(3, 0)] = -0.01;
    partner.rotors[0].mode = 1;
    partner.rotors[0].treatment = ModeTreatment::Replace;
    let mut rng = StdRng::seed_from_u64(19);
    let x = sampling::sample_fragment(&f, &partner, &mut rng).unwrap();
    assert!((distance(&x, 0, 1) - distance(&f.equilibrium, 0, 1)).abs() < 1e-12);
    partner.rotors[0].treatment = ModeTreatment::ZpeOnly;
    let x = sampling::sample_fragment(&f, &partner, &mut rng).unwrap();
    assert!((distance(&x, 0, 1) - distance(&f.equilibrium, 0, 1)).abs() > 1e-6);
}

#[test]
fn rotors_work_for_both_collision_partners_and_keep_the_com_distance() {
    use std::sync::atomic::{AtomicBool, AtomicUsize};
    let mut settings = Settings::default();
    settings.workflow = Workflow::Bimolecular;
    settings.distances = 2;
    settings.samples_per_distance = 8;
    settings.distance_start = 8.;
    settings.distance_end = 10.;
    for p in &mut settings.partners {
        p.fragment = Some(ethane());
        p.rotors = vec![rotor()];
    }
    let out =
        sampling::generate(settings, 42, &AtomicBool::new(false), &AtomicUsize::new(0)).unwrap();
    assert_eq!(out.samples.len(), 16);
    for sample in &out.samples {
        let masses = ethane().masses;
        let com = |offset: usize| {
            let mut c = [0.; 3];
            let m = masses.iter().sum::<f64>();
            for (i, mass) in masses.iter().enumerate() {
                for k in 0..3 {
                    c[k] += mass / m * sample.coordinates[3 * (offset + i) + k];
                }
            }
            c
        };
        let d = sub(com(0), com(8));
        assert!((dot(d, d).sqrt() - sample.distance.unwrap()).abs() < 1e-10);
    }
    let mut xyz = vec![];
    out.write_xyz(&mut xyz).unwrap();
    assert_eq!(
        String::from_utf8(xyz)
            .unwrap()
            .lines()
            .filter(|l| *l == "16")
            .count(),
        16
    );
}
