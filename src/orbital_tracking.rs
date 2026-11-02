//! Adjacent-frame tracking of one selected real MO by normalized cross-AO overlap.
//! O_ij = C_i(A)^T S(A,B) C_j(B). No energy-order or occupation assumption.
//! References: PySCF gto.intor_cross and its two-molecule overlap examples.
use crate::ui_style::ResponseExt;
use crate::{
    orbitals::{
        molden::{OrbitalData, Spin},
        overlap,
    },
    trajectory_wavefunctions::job::{self, Control, Series},
};
use bevy::{
    prelude::*,
    tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task},
};
use bevy_egui::egui;
use std::{
    path::PathBuf,
    sync::{atomic::Ordering, Arc},
};
#[derive(Clone, Copy, Debug)]
pub struct Choice {
    pub index: usize,
    pub phase: f64,
    pub overlap: f64,
}
#[derive(Clone)]
pub struct Track {
    pub root: PathBuf,
    pub spin: Spin,
    pub anchor: usize,
    pub original: usize,
    pub choices: Vec<Option<Choice>>,
    pub issues: Vec<String>,
}
pub fn norm(matrix: &overlap::Overlap, c: &[f64]) -> Result<f64, String> {
    if c.len() != matrix.n {
        return Err("Invalid orbital coefficient count".into());
    }
    let n = c
        .iter()
        .enumerate()
        .map(|(i, a)| {
            a * matrix.matrix[i * matrix.n..(i + 1) * matrix.n]
                .iter()
                .zip(c)
                .map(|(s, b)| s * b)
                .sum::<f64>()
        })
        .sum::<f64>();
    if !n.is_finite() || (n - 1.).abs() > 0.03 {
        return Err(format!(
            "Orbital norm {n:.5}: cannot track this Molden convention reliably."
        ));
    }
    Ok(n)
}
fn cancelled(c: &Control) -> Result<(), String> {
    if c.cancel.load(Ordering::Relaxed) {
        Err("Orbital tracking cancelled.".into())
    } else {
        Ok(())
    }
}
pub fn projection(
    a: &OrbitalData,
    b: &OrbitalData,
    c: &[f64],
    control: &Control,
) -> Result<Vec<f64>, String> {
    if a.atoms != b.atoms
        || c.len() != a.layout.n_basis
        || a.layout.n_basis > 1000
        || b.layout.n_basis > 1000
    {
        return Err("Orbital tracking requires identical atom order and at most 1000 basis functions per frame.".into());
    }
    let expand = |d: &OrbitalData| {
        d.shells
            .iter()
            .enumerate()
            .flat_map(|(i, s)| (0..s.n_functions()).map(move |slot| (i, overlap::terms(s, slot))))
            .collect::<Vec<_>>()
    };
    let left = expand(a);
    let right = expand(b);
    let mut out = vec![0.; b.layout.n_basis];
    for (mu, (si, ta)) in left.iter().enumerate() {
        cancelled(control)?;
        if c[mu].abs() < 1e-16 {
            continue;
        }
        let sa = &a.shells[*si];
        let ca = a.positions[sa.center].as_dvec3() / crate::orbitals::grid::BOHR;
        for (nu, (sj, tb)) in right.iter().enumerate() {
            let sb = &b.shells[*sj];
            let cb = b.positions[sb.center].as_dvec3() / crate::orbitals::grid::BOHR;
            let mut value = 0.;
            for (&alpha, &ac) in sa.exponents.iter().zip(&sa.coefficients) {
                for (&beta, &bc) in sb.exponents.iter().zip(&sb.coefficients) {
                    for &(pa, na) in ta {
                        for &(pb, nb) in tb {
                            value += ac
                                * bc
                                * na
                                * nb
                                * (0..3)
                                    .map(|k| {
                                        overlap::axis_overlap(
                                            pa[k], pb[k], ca[k], cb[k], alpha, beta,
                                        )
                                    })
                                    .product::<f64>();
                        }
                    }
                }
            }
            out[nu] += c[mu] * value;
        }
    }
    Ok(out)
}
pub fn next(
    a: &OrbitalData,
    b: &OrbitalData,
    spin: Spin,
    previous: Choice,
    min: f64,
    gap: f64,
    control: &Control,
) -> Result<Choice, String> {
    let old = a
        .orbitals(spin)
        .get(previous.index)
        .ok_or("Selected orbital is missing")?;
    let p = projection(a, b, &old.coefficients, control)?;
    let na = norm(&overlap::Overlap::build(a)?, &old.coefficients)?;
    let sb = overlap::Overlap::build(b)?;
    let mut scores = vec![];
    for (i, o) in b.orbitals(spin).iter().enumerate() {
        cancelled(control)?;
        let nb = norm(&sb, &o.coefficients)?;
        let value = o
            .coefficients
            .iter()
            .zip(&p)
            .map(|(c, s)| c * s)
            .sum::<f64>()
            / (na * nb).sqrt();
        if !value.is_finite() || value.abs() > 1.001 {
            return Err("Invalid normalized cross-orbital overlap".into());
        }
        scores.push((i, value));
    }
    scores.sort_by(|a, b| b.1.abs().total_cmp(&a.1.abs()));
    let &(index, value) = scores
        .first()
        .ok_or("This frame has no orbitals for the selected spin")?;
    let runner = scores.get(1).map_or(0., |v| v.1.abs());
    if value.abs() < min || value.abs() - runner < gap {
        return Err(format!("Ambiguous orbital match: best |overlap| {:.3}, next {:.3}. Use closer-spaced frames or choose a new anchor; near-degenerate orbitals may mix.",value.abs(),runner));
    }
    Ok(Choice {
        index,
        phase: previous.phase * if value < 0. { -1. } else { 1. },
        overlap: value.abs(),
    })
}
pub fn calculate(
    series: &Series,
    anchor: usize,
    index: usize,
    spin: Spin,
    min: f64,
    gap: f64,
    control: &Control,
) -> Result<Track, String> {
    let (origin, _) = job::load(series, anchor)?;
    let orbital = origin
        .orbitals(spin)
        .get(index)
        .ok_or("Select an available orbital first")?;
    norm(&overlap::Overlap::build(&origin)?, &orbital.coefficients)?;
    let mut track = Track {
        root: series.root.clone(),
        spin,
        anchor,
        original: index,
        choices: vec![None; series.manifest.frames.len()],
        issues: vec![],
    };
    let first = Choice {
        index,
        phase: 1.,
        overlap: 1.,
    };
    track.choices[anchor] = Some(first);
    control.progress.store(1, Ordering::Relaxed);
    for direction in [-1isize, 1] {
        let mut old = origin.clone();
        let mut choice = first;
        let mut i = anchor as isize + direction;
        while i >= 0 && (i as usize) < track.choices.len() {
            cancelled(control)?;
            let n = i as usize;
            let (current, _) = job::load(series, n)?;
            match next(&old, &current, spin, choice, min, gap, control) {
                Ok(c) => {
                    track.choices[n] = Some(c);
                    choice = c;
                    old = current;
                    control.progress.fetch_add(1, Ordering::Relaxed);
                }
                Err(e) => {
                    track.issues.push(format!("Frame {}: {e}", n + 1));
                    break;
                }
            }
            i += direction;
        }
    }
    cancelled(control)?;
    Ok(track)
}
pub fn apply(data: &mut OrbitalData, track: &Track, frame: usize) -> Result<Choice, String> {
    let choice=track.choices.get(frame).copied().flatten().ok_or("This frame has no unambiguous tracked orbital. Playback stopped; choose a new anchor or switch to orbital index.")?;
    let list = if track.spin == Spin::Alpha {
        &mut data.alpha
    } else {
        &mut data.beta
    };
    let orbital = list
        .get_mut(choice.index)
        .ok_or("Tracked orbital is missing")?;
    for c in &mut orbital.coefficients {
        *c *= choice.phase;
    }
    Ok(choice)
}
struct Running {
    task: Task<Result<Track, String>>,
    control: Control,
}
impl Drop for Running {
    fn drop(&mut self) {
        self.control.cancel();
    }
}
pub struct State {
    pub enabled: bool,
    pub track: Option<Arc<Track>>,
    running: Option<Running>,
    min: f64,
    gap: f64,
    pub message: String,
}
impl Default for State {
    fn default() -> Self {
        Self {
            enabled: false,
            track: None,
            running: None,
            min: 0.5,
            gap: 0.1,
            message: String::new(),
        }
    }
}
impl State {
    pub fn poll(&mut self, series: Option<&Series>) {
        if let Some(run) = &mut self.running {
            if let Some(result) = block_on(future::poll_once(&mut run.task)) {
                self.running = None;
                match result {
                    Ok(t)
                        if series.is_some_and(|s| {
                            s.root == t.root && s.manifest.frames.len() == t.choices.len()
                        }) =>
                    {
                        self.message = if t.issues.is_empty() {
                            "Orbital character tracked through all frames.".into()
                        } else {
                            t.issues.join("\n")
                        };
                        self.track = Some(Arc::new(t));
                        self.enabled = true;
                    }
                    Ok(_) => {
                        self.message="The wavefunction series changed; the old tracking result was discarded.".into();
                    }
                    Err(e) => self.message = e,
                }
            }
        }
        if self.track.as_ref().is_some_and(|t| {
            !series.is_some_and(|s| s.root == t.root && s.manifest.frames.len() == t.choices.len())
        }) {
            self.track = None;
            self.enabled = false;
        }
    }
    pub fn ui(
        &mut self,
        ui: &mut egui::Ui,
        series: &Arc<Series>,
        orbital: &crate::orbitals::OrbitalState,
        frame: usize,
    ) {
        ui.horizontal_wrapped(|ui| {
            ui.add_enabled(
                self.track.is_some(),
                egui::Checkbox::new(&mut self.enabled, "Follow orbital character"),
            );
            let valid = orbital.quantity == crate::orbitals::density::Quantity::Orbital
                && orbital.selected.is_some()
                && orbital.data.as_ref().is_some_and(|d| {
                    series
                        .manifest
                        .frames
                        .get(frame)
                        .is_some_and(|f| job::check_geometry(d, &f.geometry).is_ok())
                });
            if ui
                .add_enabled(
                    valid && self.running.is_none(),
                    egui::Button::new("Track selected orbital"),
                )
                .clicked_once()
            {
                let s = series.clone();
                let i = orbital.selected.unwrap();
                let spin = orbital.spin;
                let min = self.min;
                let gap = self.gap;
                let control = Control::default();
                let c = control.clone();
                self.enabled = false;
                self.message.clear();
                self.running = Some(Running {
                    task: AsyncComputeTaskPool::get()
                        .spawn(async move { calculate(&s, frame, i, spin, min, gap, &c) }),
                    control,
                });
            }
        });
        egui::CollapsingHeader::new("Tracking settings").default_open(false).show(ui,|ui|{ui.horizontal(|ui|{ui.label("Minimum |overlap|");ui.add(egui::DragValue::new(&mut self.min).range(0.1..=0.99).speed(0.01));ui.label("Separation from next match");ui.add(egui::DragValue::new(&mut self.gap).range(0.01..=0.5).speed(0.01));});ui.small("Matches adjacent frames in their stored orientation; large rotations and degenerate mixtures can be ambiguous. Raw wavefunctions are unchanged.");});
        if let Some(run) = &self.running {
            ui.label(format!(
                "Tracking: {} / {} frames",
                run.control.progress.load(Ordering::Relaxed),
                series.manifest.frames.len()
            ));
            if ui.button("Cancel tracking").clicked_once() {
                run.control.cancel();
            }
            ui.ctx().request_repaint();
        }
        if let Some(t) = &self.track {
            if let Some(c) = t.choices.get(frame).copied().flatten() {
                ui.label(format!(
                    "Anchor: frame {}, orbital {} → orbital {} · |overlap| {:.3}",
                    t.anchor + 1,
                    t.original + 1,
                    c.index + 1,
                    c.overlap
                ));
            }
        }
        if !self.message.is_empty() {
            ui.label(&self.message);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn molecule() -> OrbitalData {
        let overlap = (-0.5_f64).exp();
        let a = 1. / (2. * (1. + overlap)).sqrt();
        let b = 1. / (2. * (1. - overlap)).sqrt();
        let text=format!("[Molden Format]\n[Atoms] AU\nH 1 1 0 0 0\nH 2 1 1 0 0\n[GTO]\n1 0\ns 1 1.0\n1.0 1.0\n\n2 0\ns 1 1.0\n1.0 1.0\n\n[MO]\nSym= 1\nEne= -0.5\nSpin= Alpha\nOccup= 2.0\n1 {a:.15}\n2 {a:.15}\nSym= 2\nEne= 0.2\nSpin= Alpha\nOccup= 0.0\n1 {b:.15}\n2 {:.15}\n",-b);
        crate::orbitals::molden::parse_molden(&text).unwrap()
    }
    #[test]
    fn tracks_reordered_and_sign_flipped_orbitals() {
        let a = molecule();
        let mut b = a.clone();
        b.alpha.swap(0, 1);
        for c in &mut b.alpha[1].coefficients {
            *c *= -1.;
        }
        let c = next(
            &a,
            &b,
            Spin::Alpha,
            Choice {
                index: 0,
                phase: 1.,
                overlap: 1.,
            },
            0.5,
            0.1,
            &Control::default(),
        )
        .unwrap();
        assert_eq!(c.index, 1);
        assert_eq!(c.phase, -1.);
        assert!((c.overlap - 1.).abs() < 1e-7);
        let back = next(&b, &a, Spin::Alpha, c, 0.5, 0.1, &Control::default()).unwrap();
        assert_eq!(back.index, 0);
        assert_eq!(back.phase, 1.);
        let track = Track {
            root: PathBuf::new(),
            spin: Spin::Alpha,
            anchor: 0,
            original: 0,
            choices: vec![Some(c), None],
            issues: vec![],
        };
        apply(&mut b, &track, 0).unwrap();
        assert_eq!(b.alpha[1].coefficients, a.alpha[0].coefficients);
        assert!(apply(&mut b, &track, 1).is_err());
    }
    #[test]
    fn rejects_degenerate_mixtures_and_wrong_atom_order() {
        let a = molecule();
        let mut b = a.clone();
        for i in 0..2 {
            b.alpha[0].coefficients[i] =
                (a.alpha[0].coefficients[i] + a.alpha[1].coefficients[i]) / 2_f64.sqrt();
            b.alpha[1].coefficients[i] =
                (a.alpha[0].coefficients[i] - a.alpha[1].coefficients[i]) / 2_f64.sqrt();
        }
        let start = Choice {
            index: 0,
            phase: 1.,
            overlap: 1.,
        };
        assert!(
            next(&a, &b, Spin::Alpha, start, 0.5, 0.1, &Control::default())
                .err()
                .unwrap()
                .contains("Ambiguous")
        );
        b.atoms[0] = "He".into();
        assert!(next(&a, &b, Spin::Alpha, start, 0.5, 0.1, &Control::default()).is_err());
    }
    #[test]
    fn cancellation_stops_overlap_work() {
        let a = molecule();
        let control = Control::default();
        control.cancel();
        assert!(projection(&a, &a, &a.alpha[0].coefficients, &control).is_err());
    }
}

#[cfg(test)]
mod independent_reference {
    use super::*;
    #[test]
    fn cross_overlap_matches_pyscf_for_displaced_spd_basis() {
        let a = crate::orbitals::molden::parse_molden(include_str!(
            "../tests/fixtures/tracking_water_a.molden"
        ))
        .unwrap();
        let b = crate::orbitals::molden::parse_molden(include_str!(
            "../tests/fixtures/tracking_water_b.molden"
        ))
        .unwrap();
        let expected: Vec<f64> = include_str!("../tests/fixtures/tracking_water_overlap.txt")
            .lines()
            .map(|s| s.parse().unwrap())
            .collect();
        let p = projection(&a, &b, &a.alpha[4].coefficients, &Control::default()).unwrap();
        for (o, e) in b.alpha.iter().zip(expected) {
            let actual = o
                .coefficients
                .iter()
                .zip(&p)
                .map(|(c, s)| c * s)
                .sum::<f64>();
            assert!((actual - e).abs() < 2e-6, "overlap: {actual} vs PySCF {e}");
        }
        let choice = next(
            &a,
            &b,
            Spin::Alpha,
            Choice {
                index: 4,
                phase: 1.,
                overlap: 1.,
            },
            0.5,
            0.1,
            &Control::default(),
        )
        .unwrap();
        assert_eq!(choice.index, 4);
        assert!(choice.overlap > 0.99);
    }
}
