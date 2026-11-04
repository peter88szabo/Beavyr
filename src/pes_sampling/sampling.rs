//! Smite's harmonic QCT phases, thermal quanta and isotropic Euler rotations.
use super::{
    input::{self, Fragment, BOHR},
    output::{Sample, SamplingOutput},
    state::*,
};
use super::quasi::{Draws, Unit};
use crate::optimizer::prng::StdRng;
use std::{
    f64::consts::TAU,
    sync::atomic::{AtomicBool, AtomicUsize, Ordering},
};
pub const KB: f64 = 8.3144598 / 1000. / 2625.5;
#[cfg(test)]
pub fn uniform(rng: &mut StdRng) -> f64 {
    rng.gen_range(0.0..1.0)
}
pub fn thermal_quantum(omega: f64, temperature: f64, u: f64) -> f64 {
    (-KB * temperature * (-u).ln_1p() / omega).floor()
}
pub fn harmonic_state(omega: f64, energy: f64, phase: f64) -> (f64, f64) {
    let amplitude = (2. * energy).sqrt() / omega;
    (amplitude * phase.cos(), -amplitude * omega * phase.sin())
}
fn excitation_valid(e: &Excitation) -> Result<(), String> {
    if !e.temperature.is_finite()
        || e.temperature < 0.
        || !e.energy_kcal.is_finite()
        || e.energy_kcal < 0.
        || !e.retained_population.is_finite()
        || e.retained_population <= 0.
        || e.retained_population >= 100.
    {
        Err("Temperature and energy must be finite and nonnegative; retained population must lie strictly between 0 and 100%".into())
    } else {
        Ok(())
    }
}
pub fn tail_quantum(w: f64, t: f64, population: f64) -> f64 {
    if t == 0. {
        return 0.;
    }
    ((-(1. - population / 100.).ln() * KB * t / w).ceil() - 1.).max(0.)
}
pub fn energy(e: &Excitation, w: f64, rng: &mut impl Unit) -> f64 {
    match e.population {
        Population::Zpe => 0.5 * w,
        Population::Quantum => w * (e.quantum as f64 + 0.5),
        Population::Thermal => match e.temperature_rule {
            TemperatureRule::BoltzmannTail => {
                w * (tail_quantum(w, e.temperature, e.retained_population) + 0.5)
            }
            TemperatureRule::ZpePlusKt => 0.5 * w + KB * e.temperature,
            TemperatureRule::SmiteThermal => {
                w * (thermal_quantum(w, e.temperature, rng.unit()) + 0.5)
            }
        },
        Population::Energy => e.energy_kcal / 627.51,
    }
}
pub fn rotate(x: &mut [f64], u: [f64; 3]) {
    // Identical Euler matrix to Smite utils/euler.py, cos(theta) uniform [-1,1].
    let (sp, cp) = (TAU * u[0]).sin_cos();
    let ct = 2. * u[1] - 1.;
    let st = (1. - ct * ct).max(0.).sqrt();
    let (sc, cc) = (TAU * u[2]).sin_cos();
    let r = [
        [ct * cp * cc - sp * sc, -ct * cp * sc - sp * cc, st * cp],
        [ct * sp * cc + cp * sc, -ct * sp * sc + cp * cc, st * sp],
        [-st * cc, st * sc, ct],
    ];
    for p in x.chunks_exact_mut(3) {
        let a = [p[0], p[1], p[2]];
        for i in 0..3 {
            p[i] = (0..3).map(|k| r[i][k] * a[k]).sum();
        }
    }
}
fn morse_energy(w: f64, de: f64, n: f64) -> f64 {
    let eh = w * (n + 0.5);
    eh - eh * eh / (4. * de)
}
fn morse_bound(w: f64, de: f64, n: f64) -> bool {
    n >= 0. && w * (n + 0.5) < 2. * de && morse_energy(w, de, n) < de
}
fn morse_level(e: &Excitation, w: f64, de: f64, rng: &mut impl Unit) -> Result<f64, String> {
    match e.population {
        Population::Energy => {
            let value = e.energy_kcal / 627.51;
            if value >= de {
                return Err("Morse energy must be below the dissociation energy".into());
            }
            Ok(value)
        }
        Population::Thermal if e.temperature_rule == TemperatureRule::ZpePlusKt => {
            let value = morse_energy(w, de, 0.) + KB * e.temperature;
            if value >= de {
                return Err(
                    "ZPE + kB × Tmax reaches the Morse dissociation energy. Reduce Tmax.".into(),
                );
            }
            Ok(value)
        }
        Population::Thermal if e.temperature > 0. => {
            let limit = (2. * de / w - 0.5).ceil().max(0.) as usize;
            if limit > 100000 {
                return Err("Morse parameters imply too many bound states; check the units".into());
            }
            let e0 = morse_energy(w, de, 0.);
            let mut weights = vec![];
            for n in 0..limit {
                if morse_bound(w, de, n as f64) {
                    weights.push((
                        morse_energy(w, de, n as f64),
                        (-(morse_energy(w, de, n as f64) - e0) / (KB * e.temperature)).exp(),
                    ));
                }
            }
            let total: f64 = weights.iter().map(|(_, p)| p).sum();
            if total <= 0. {
                return Err("The Morse oscillator has no bound vibrational ground state".into());
            }
            let target = if e.temperature_rule == TemperatureRule::BoltzmannTail {
                e.retained_population / 100. * total
            } else {
                rng.unit() * total
            };
            let mut c = 0.;
            for (en, p) in &weights {
                c += p;
                if target < c {
                    return Ok(*en);
                }
            }
            Ok(weights.last().unwrap().0)
        }
        _ => {
            let n = if e.population == Population::Quantum {
                e.quantum as f64
            } else {
                0.
            };
            if !morse_bound(w, de, n) {
                return Err("Requested Morse vibrational level is not bound".into());
            }
            Ok(morse_energy(w, de, n))
        }
    }
}
pub fn morse_distance(
    re: f64,
    mu: f64,
    w: f64,
    de: f64,
    en: f64,
    phase: f64,
) -> Result<f64, String> {
    // Smite prmconst/rprsetup at J=0. beta is fixed by the small-amplitude frequency.
    if en < 0. || en >= de {
        return Err("Morse sampling requires 0 ≤ energy < De".into());
    }
    let beta = w * (mu / (2. * de)).sqrt();
    if en == 0. {
        return Ok(re);
    }
    let a = en - de;
    let b = 2. * de;
    let c = -de;
    let xi = -2. * a / (b + phase.sin() * (b * b - 4. * a * c).sqrt());
    let r = re - xi.ln() / beta;
    if !r.is_finite() || r <= 0. {
        return Err("Morse parameters produced an invalid bond distance".into());
    }
    Ok(r)
}
/// Harmonic normal-coordinate displacement at mode energy `energy`, drawn
/// within the classical turning points ±A, A = sqrt(2E)/ω.
pub fn harmonic_displacement(omega: f64, energy: f64, d: Displacement, u: f64) -> f64 {
    match d {
        Displacement::Uniform => (2. * energy).sqrt() / omega * (2. * u - 1.),
        Displacement::ClassicalPhase => harmonic_state(omega, energy, TAU * u).0,
    }
}
/// Morse bond length at energy `en`: uniform between the inner and outer
/// turning points, or the exact classical time-weighted distribution.
pub fn morse_displacement(
    re: f64,
    mu: f64,
    w: f64,
    de: f64,
    en: f64,
    d: Displacement,
    u: f64,
) -> Result<f64, String> {
    match d {
        Displacement::Uniform => {
            let half_pi = std::f64::consts::FRAC_PI_2;
            let inner = morse_distance(re, mu, w, de, en, -half_pi)?;
            let outer = morse_distance(re, mu, w, de, en, half_pi)?;
            Ok(inner + (outer - inner) * u)
        }
        Displacement::ClassicalPhase => morse_distance(re, mu, w, de, en, TAU * u),
    }
}
/// Per-fragment data prepared once per run.
pub struct Prepared {
    pub rotors: Vec<super::rotors::Prepared>,
    pub internal: Option<super::curvilinear::Internal>,
}

pub fn prepare_fragment(f: &Fragment, p: &Partner, coordinates: Coordinates) -> Result<Prepared, String> {
    if p.rigid {
        return Ok(Prepared { rotors: vec![], internal: None });
    }
    let rotors = super::rotors::prepare_all(f, &p.rotors)?;
    let internal = if coordinates == Coordinates::Internal
        && f.diatomic.is_none()
        && f.atoms.len() >= 3
    {
        Some(super::curvilinear::prepare(f)?)
    } else {
        None
    };
    Ok(Prepared { rotors, internal })
}

#[cfg(test)]
pub fn sample_fragment(f: &Fragment, p: &Partner, rng: &mut StdRng) -> Result<Vec<f64>, String> {
    let prepared = prepare_fragment(f, p, Coordinates::Cartesian)?;
    sample_prepared_fragment(f, p, &prepared, Displacement::default(), rng)
}

pub fn sample_prepared_fragment(
    f: &Fragment,
    p: &Partner,
    prepared: &Prepared,
    displacement: Displacement,
    rng: &mut impl Unit,
) -> Result<Vec<f64>, String> {
    let mut x = f.equilibrium.clone();
    if !p.rigid {
        if let Some((re, mu, Some(de))) = f.diatomic {
            let e = if p.use_overrides {
                p.overrides.first().unwrap_or(&p.excitation)
            } else {
                &p.excitation
            };
            let en = morse_level(e, f.frequencies[0], de, rng)?;
            let r = morse_displacement(re, mu, f.frequencies[0], de, en, displacement, rng.unit())?;
            x = vec![r, 0., 0., 0., 0., 0.];
        } else {
            let mut dx = vec![0.; x.len()];
            for (j, w) in f.frequencies.iter().enumerate() {
                let e = if p.use_overrides {
                    p.overrides.get(j).unwrap_or(&p.excitation)
                } else {
                    &p.excitation
                };
                let en = match prepared
                    .rotors
                    .iter()
                    .find(|r| r.mode == Some(j))
                    .map(|r| r.treatment)
                {
                    Some(super::rotors::ModeTreatment::ZpeOnly) => 0.5 * w,
                    Some(super::rotors::ModeTreatment::Replace) => 0.,
                    None => energy(e, *w, rng),
                };
                let q = harmonic_displacement(*w, en, displacement, rng.unit());
                for (i, v) in dx.iter_mut().enumerate() {
                    *v += f.modes[(i, j)] * q;
                }
            }
            if p.ts_mode.enabled && f.imaginary_cm1.len() == 1 {
                super::ts_mode::apply(&mut dx, f, &p.ts_mode, rng.unit())?;
            }
            match &prepared.internal {
                Some(internal) => x = super::curvilinear::displace(&f.equilibrium, &dx, internal)?,
                None => x.iter_mut().zip(&dx).for_each(|(v, d)| *v += d),
            }
            if f.diatomic.is_some() && x[0] - x[3] <= 0. {
                return Err("Sampled harmonic diatomic bond length is nonpositive. Reduce its excitation or use Morse sampling.".into());
            }
        }
        for rotor in &prepared.rotors {
            super::rotors::rotate(&mut x, rotor, TAU * rng.unit())?;
        }
    }
    input::center(&mut x, &f.masses);
    if x.iter().any(|v| !v.is_finite()) {
        return Err(
            "Sampling produced non-finite coordinates; check frequencies and excitation".into(),
        );
    }
    Ok(x)
}

/// One group of samples sharing a COM distance or path point. Diverse-subset
/// selection works within a block, so grids keep their exact counts.
#[derive(Clone, Debug)]
pub struct Block {
    pub unit: usize,
    pub distance: Option<f64>,
    pub distance_index: Option<usize>,
    pub path_point: Option<usize>,
    pub count: usize,
}

pub fn blocks(s: &Settings) -> Vec<Block> {
    let plain = |unit, count| Block { unit, distance: None, distance_index: None, path_point: None, count };
    match s.workflow {
        Workflow::Unimolecular => vec![plain(0, s.points)],
        Workflow::Bimolecular => (0..s.distances)
            .map(|i| Block {
                distance: Some(if s.distances == 1 {
                    s.distance_start
                } else {
                    s.distance_start
                        + (s.distance_end - s.distance_start) * i as f64 / (s.distances - 1) as f64
                }),
                distance_index: Some(i + 1),
                ..plain(0, s.samples_per_distance)
            })
            .collect(),
        Workflow::ReactionPath => (0..s.path.points.len())
            .map(|k| Block { path_point: Some(k + 1), ..plain(k, s.path.samples_per_point) })
            .collect(),
    }
}

pub fn total(s: &Settings) -> Result<usize, String> {
    if s.workflow == Workflow::ReactionPath && s.path.points.is_empty() {
        return Err("Load the Hessians along the reaction path".into());
    }
    let n = blocks(s)
        .iter()
        .try_fold(0usize, |a, b| a.checked_add(b.count))
        .ok_or("Too many samples")?;
    if n == 0 {
        return Err("Sample, distance and path-point counts must be at least one".into());
    }
    Ok(n)
}

/// Progress units: candidate generation plus farthest-point selection steps.
pub fn work(s: &Settings) -> usize {
    blocks(s)
        .iter()
        .map(|b| {
            let c = s.diversity.candidates(b.count);
            c + if c > b.count { b.count } else { 0 }
        })
        .sum()
}

fn path_partner(s: &Settings) -> Partner {
    let mut p = Partner { excitation: s.path.excitation.clone(), ..Default::default() };
    p.ts_mode.enabled = false;
    p
}

/// Fragments, the partner settings that drive them, and prepared data.
pub type Units = Vec<(Fragment, Partner, Prepared)>;

pub fn validate(s: &Settings) -> Result<Units, String> {
    let count = total(s)?;
    if s.workflow == Workflow::Bimolecular
        && (![s.distance_start, s.distance_end]
            .iter()
            .all(|v| v.is_finite() && *v > 0.))
    {
        return Err("COM distances must be positive and finite (Å)".into());
    }
    s.overlap.validate()?;
    s.diversity.validate()?;
    let mut units = vec![];
    if s.workflow == Workflow::ReactionPath {
        excitation_valid(&s.path.excitation)?;
        let p = path_partner(s);
        for (k, f) in s.path.points.iter().enumerate() {
            let prepared = prepare_fragment(f, &p, s.coordinates)
                .map_err(|e| format!("Path point {}: {e}", k + 1))?;
            units.push((f.clone(), p.clone(), prepared));
        }
    } else {
        let n = if s.workflow == Workflow::Bimolecular { 2 } else { 1 };
        for p in &s.partners[..n] {
            excitation_valid(&p.excitation)?;
            for e in &p.overrides {
                excitation_valid(e)?;
            }
            let f = input::prepare(p)?;
            if !p.rigid {
                super::ts_mode::validate(&f, &p.ts_mode)?;
            }
            if let Some((_, _, Some(de))) = f.diatomic {
                let e = if p.use_overrides {
                    p.overrides.first().unwrap_or(&p.excitation)
                } else {
                    &p.excitation
                };
                if !p.rigid {
                    morse_level(e, f.frequencies[0], de, &mut StdRng::seed_from_u64(0))?;
                }
            }
            let prepared = prepare_fragment(&f, p, s.coordinates)?;
            units.push((f, p.clone(), prepared));
        }
    }
    let nat: usize = if s.workflow == Workflow::Bimolecular {
        units.iter().map(|u| u.0.atoms.len()).sum()
    } else {
        units[0].0.atoms.len()
    };
    if count.checked_mul(nat).is_none_or(|n| n > 5_000_000) {
        return Err("This dataset exceeds the automatic trajectory preview limit (5 million atom-frames). Reduce the number of samples.".into());
    }
    if s.diversity.enabled {
        let d = super::diverse::descriptor_len(nat);
        let mut selection = 0f64;
        for b in blocks(s) {
            let c = s.diversity.candidates(b.count);
            if c.saturating_mul(nat) > 20_000_000 || c.saturating_mul(d) > super::diverse::MAX_DESCRIPTOR_VALUES {
                return Err("Too many candidate geometries for the diverse subset. Reduce the samples per block or the oversampling, or disable the diverse subset.".into());
            }
            selection += b.count as f64 * c as f64 * d as f64;
        }
        if selection > super::diverse::MAX_WORK {
            return Err("Diverse-subset selection would take too long for this many samples. Reduce the samples per block or the oversampling, or disable the diverse subset.".into());
        }
    }
    Ok(units)
}

pub fn generate(
    s: Settings,
    seed: u64,
    cancel: &AtomicBool,
    progress: &AtomicUsize,
) -> Result<SamplingOutput, String> {
    let units = validate(&s)?;
    let bimolecular = s.workflow == Workflow::Bimolecular;
    let atoms: Vec<String> = if bimolecular {
        units.iter().flat_map(|u| u.0.atoms.clone()).collect()
    } else {
        units[0].0.atoms.clone()
    };
    let mut draws = Draws::new(s.random, seed);
    let mut samples = Vec::with_capacity(total(&s)?);
    let mut rejected = 0;
    let mut warnings = vec![];
    if s.workflow == Workflow::ReactionPath {
        let bent: Vec<String> = units
            .iter()
            .enumerate()
            .filter(|(_, u)| !u.0.imaginary_cm1.is_empty())
            .map(|(k, u)| format!("{} ({:?})", k + 1, u.0.imaginary_cm1))
            .collect();
        if !bent.is_empty() {
            warnings.push(format!("Path points with negative curvature orthogonal to the path; those modes are excluded (cm^-1): {}", bent.join(", ")));
        }
    } else {
        for (i, (f, p, prepared)) in units.iter().enumerate() {
            if f.imaginary_cm1.len() == 1 && p.ts_mode.enabled && !p.rigid {
                warnings.push(format!("Partner {}: TS mode {:.1} cm^-1 sampled with bounded uniform displacement, {}, maximum atom displacement {} Å (no thermal/ZPE assignment).", i+1, f.imaginary_cm1[0], p.ts_mode.side.label(), p.ts_mode.max_displacement_angstrom));
            } else if !f.imaginary_cm1.is_empty() {
                warnings.push(format!(
                    "Partner {}: excluded imaginary modes (cm^-1): {:?}",
                    i + 1,
                    f.imaginary_cm1
                ));
            }
            if !prepared.rotors.is_empty() {
                warnings.push(format!("Partner {}: {} free internal rotor(s), uniform 0–360°; bond detection from equilibrium geometry.", i+1, prepared.rotors.len()));
            }
        }
    }
    let cancelled = || Err("Cancelled. Previous results were retained.".to_string());
    for block in blocks(&s) {
        let candidates = s.diversity.candidates(block.count);
        let mut geometries = Vec::with_capacity(candidates);
        for index in 0..candidates {
            let mut accepted = None;
            for _ in 0..1000 {
                if cancel.load(Ordering::Relaxed) {
                    return cancelled();
                }
                draws.next_point();
                let (f, p, prepared) = &units[block.unit];
                let mut x = sample_prepared_fragment(f, p, prepared, s.displacement, &mut draws)?;
                if let Some(d) = block.distance {
                    rotate(&mut x, [draws.unit(), draws.unit(), draws.unit()]);
                    let (g, q, prepared_b) = &units[1];
                    let mut b = sample_prepared_fragment(g, q, prepared_b, s.displacement, &mut draws)?;
                    rotate(&mut b, [draws.unit(), draws.unit(), draws.unit()]);
                    let ma: f64 = f.masses.iter().sum();
                    let mb: f64 = g.masses.iter().sum();
                    for a in x.chunks_exact_mut(3) {
                        a[2] -= d / BOHR * mb / (ma + mb);
                    }
                    for a in b.chunks_exact_mut(3) {
                        a[2] += d / BOHR * ma / (ma + mb);
                    }
                    x.extend(b);
                }
                for v in &mut x {
                    *v *= BOHR;
                }
                if s.overlap.accepts(&atoms, &x) {
                    accepted = Some(x);
                    break;
                }
                rejected += 1;
            }
            let x = accepted.ok_or_else(|| format!("Could not generate sample {} after 1000 attempts. The overlap filter may be incompatible with this distance or excitation.", samples.len() + index + 1))?;
            geometries.push(x);
            progress.fetch_add(1, Ordering::Relaxed);
        }
        let order = if candidates > block.count {
            super::diverse::farthest_point(&geometries, block.count, cancel, progress)?
        } else {
            (0..block.count).collect()
        };
        for (j, i) in order.into_iter().enumerate() {
            samples.push(Sample {
                coordinates: std::mem::take(&mut geometries[i]),
                distance: block.distance,
                distance_index: block.distance_index,
                sample_at_distance: block.distance.map(|_| j + 1),
                path_point: block.path_point,
            });
        }
    }
    if cancel.load(Ordering::Relaxed) {
        return cancelled();
    }
    let mut notes = vec![
        s.random.comment().to_string(),
        s.coordinates.comment().to_string(),
    ];
    if s.diversity.enabled {
        notes.push(format!(
            "diverse_subset=farthest_point_inverse_distances:oversampling={}",
            s.diversity.oversampling
        ));
    }
    if s.workflow == Workflow::ReactionPath {
        notes.push("path_projection=reaction_path_hamiltonian".into());
    }
    for (partner, (f, p, prepared)) in units.iter().enumerate() {
        if p.rigid || s.workflow == Workflow::ReactionPath {
            continue;
        }
        if f.imaginary_cm1.len() == 1 && p.ts_mode.enabled {
            notes.push(format!(
                "partner{}_TS=uniform_displacement:{:?}:max_atom_A={:.8}",
                partner + 1,
                p.ts_mode.side,
                p.ts_mode.max_displacement_angstrom
            ));
        }
        for (index, r) in prepared.rotors.iter().enumerate() {
            notes.push(format!(
                "partner{}_rotor{}=axis_{}-{}:atoms_{}:mode_{}:{:?}",
                partner + 1,
                index + 1,
                r.axis[0] + 1,
                r.axis[1] + 1,
                r.moving
                    .iter()
                    .map(|i| (i + 1).to_string())
                    .collect::<Vec<_>>()
                    .join(","),
                r.mode.map_or(0, |m| m + 1),
                r.treatment
            ));
        }
    }
    Ok(SamplingOutput {
        atoms,
        samples,
        seed,
        rejected,
        warnings,
        partner_sizes: if bimolecular {
            units.iter().map(|u| u.0.atoms.len()).collect()
        } else {
            vec![units[0].0.atoms.len()]
        },
        sampling_note: notes.join(" "),
        displacement: s.displacement,
        path_names: if s.workflow == Workflow::ReactionPath {
            s.path.points.iter().map(|f| f.name.clone()).collect()
        } else {
            vec![]
        },
    })
}
