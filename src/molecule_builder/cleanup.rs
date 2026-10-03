//! Cancellable, constrained cleanup on a copy. Only Accept changes the molecule.
use std::collections::BTreeSet;
use std::sync::{
    atomic::{AtomicBool, Ordering},
    Arc, Mutex,
};

use anyhow::{bail, Result};
use bevy::prelude::*;
use bevy_egui::egui;

use super::{builder_ui::ZMatrixBuilderState, selection::Selection};
use crate::forcefield::dreiding::{
    objective::{cleanup_options, DreidingObjective, BOHR_TO_ANGSTROM, HARTREE_TO_KCAL},
    DreidingTopology,
};
use crate::molecule::{is_dummy, Molecule};
use crate::optimizer::{minimum_cartesian::minimize_cartesian, traits::Objective};

#[derive(Clone, Copy, Default, PartialEq, Eq)]
pub enum Scope {
    #[default]
    Whole,
    Selected,
}

#[derive(Clone)]
pub struct Preview {
    atoms: Vec<String>,
    original: Vec<Vec3>,
    pub positions: Vec<Vec3>,
    pub mobile: BTreeSet<usize>,
    pub report: String,
}

impl Preview {
    pub fn matches(&self, mol: &Molecule) -> bool {
        self.atoms == mol.atoms && self.original == mol.pos
    }
}

struct Job {
    mobile: BTreeSet<usize>,
    cancel: Arc<AtomicBool>,
    result: Arc<Mutex<Option<std::result::Result<Preview, String>>>>,
}
impl Drop for Job {
    fn drop(&mut self) {
        self.cancel.store(true, Ordering::Relaxed);
    }
}

#[derive(Clone, Default)]
pub struct Cleanup {
    pub open: bool,
    pub rect: Option<egui::Rect>,
    pub scope: Scope,
    pub preview: Option<Preview>,
    job: Option<Arc<Job>>,
    error: Option<String>,
}

impl Cleanup {
    pub fn cancel(&mut self) {
        if let Some(job) = self.job.take() {
            job.cancel.store(true, Ordering::Relaxed);
        }
        self.preview = None;
    }

    fn start(&mut self, mol: &Molecule, selection: &Selection) {
        self.cancel();
        self.error = None;
        let mol = mol.clone();
        let selection = selection.clone();
        let scope = self.scope;
        let cancel = Arc::new(AtomicBool::new(false));
        let result = Arc::new(Mutex::new(None));
        self.job = Some(Arc::new(Job {
            mobile: movable(&mol, &selection, scope),
            cancel: cancel.clone(),
            result: result.clone(),
        }));
        bevy::tasks::AsyncComputeTaskPool::get()
            .spawn(async move {
                let output = relax(&mol, &selection, scope, &cancel).map_err(|e| e.to_string());
                if let Ok(mut slot) = result.lock() {
                    *slot = Some(output);
                }
            })
            .detach();
    }

    pub fn poll(&mut self, mol: &Molecule, selection: &Selection) {
        if self.job.is_some() || self.preview.is_some() {
            let current = movable(mol, selection, self.scope);
            if self.job.as_ref().is_some_and(|j| j.mobile != current)
                || self.preview.as_ref().is_some_and(|p| p.mobile != current)
            {
                self.cancel();
                self.error =
                    Some("Selection or frozen atoms changed. Preview cleanup again.".into());
                return;
            }
        }
        let output = self
            .job
            .as_ref()
            .and_then(|job| job.result.lock().ok()?.take());
        if let Some(output) = output {
            self.job = None;
            match output {
                Ok(preview) if preview.matches(mol) => self.preview = Some(preview),
                Ok(_) => {
                    self.error = Some(
                        "The structure changed. Run cleanup again on the current geometry.".into(),
                    )
                }
                Err(error) => self.error = Some(error),
            }
        }
        if self.preview.as_ref().is_some_and(|p| !p.matches(mol)) {
            self.preview = None;
            self.error = Some("The structure changed; the old preview was discarded.".into());
        }
    }

    pub fn show(
        &mut self,
        ctx: &egui::Context,
        mol: &Molecule,
        selection: &Selection,
        available: bool,
    ) -> bool {
        self.rect = None;
        if !self.open {
            return false;
        }
        self.poll(mol, selection);
        let mut open = true;
        let mut accept = false;
        let mut discard = false;
        let response = egui::Window::new("Quick geometry cleanup")
            .id(egui::Id::new("quick_geometry_cleanup"))
            .open(&mut open).default_width(370.0).resizable(true)
            .show(ctx, |ui| {
                ui.strong("DREIDING · built-in force field");
                ui.small("For preparing a starting geometry. Metal coordination, unusual bonding and electronic states need a suitable quantum method.");
                ui.separator();
                ui.add_enabled_ui(self.job.is_none() && self.preview.is_none(), |ui| {
                    ui.horizontal(|ui| {
                        ui.selectable_value(&mut self.scope, Scope::Whole, "Whole structure");
                        ui.selectable_value(&mut self.scope, Scope::Selected, "Selected atoms");
                    });
                });
                let mobile = movable(mol, selection, self.scope);
                if let Some(preview) = &self.preview {
                    ui.label(format!("{} moving · {} fixed in this preview", preview.mobile.len(), mol.atoms.len() - preview.mobile.len()));
                    ui.label(&preview.report);
                    ui.small("Cyan: proposed geometry and movement. Current coordinates remain unchanged until Accept.");
                    ui.horizontal(|ui| {
                        accept = ui.add_enabled(available, egui::Button::new("Accept cleanup")).clicked();
                        discard = ui.button("Discard").clicked();
                    });
                } else if self.job.is_some() {
                    ui.horizontal(|ui| { ui.spinner(); ui.label("Relaxing geometry…"); });
                    if ui.button("Cancel cleanup").clicked() { self.cancel(); }
                    ctx.request_repaint_after(std::time::Duration::from_millis(50));
                } else {
                    ui.label(format!("{} to move · {} fixed", mobile.len(), mol.atoms.len() - mobile.len()));
                    if self.scope == Scope::Selected { ui.small("Unselected atoms stay fixed. Frozen atoms always stay fixed."); }
                    if ui.add_enabled(available && !mobile.is_empty(), egui::Button::new("Preview cleanup")).clicked() { self.start(mol, selection); }
                    if mobile.is_empty() { ui.weak("Select movable atoms or unfreeze them to continue."); }
                }
                if !available { ui.weak("Pause playback and apply or cancel fragment adjustments first."); }
                if let Some(error) = &self.error { ui.colored_label(egui::Color32::LIGHT_RED, error); }
            });
        self.rect = response.map(|r| r.response.rect);
        if !open || discard {
            self.cancel();
            self.open = false;
        }
        accept
    }

    pub fn accepted(
        &mut self,
        mol: &Molecule,
        selection: &Selection,
    ) -> std::result::Result<Preview, String> {
        let preview = self.preview.as_ref().ok_or("Preview cleanup first.")?;
        if !preview.matches(mol) || preview.mobile != movable(mol, selection, self.scope) {
            return Err("The structure or selection changed. Preview cleanup again.".into());
        }
        self.open = false;
        Ok(self.preview.take().unwrap())
    }
}

fn movable(mol: &Molecule, selection: &Selection, scope: Scope) -> BTreeSet<usize> {
    (0..mol.atoms.len())
        .filter(|i| {
            !is_dummy(&mol.atoms[*i])
                && !selection.frozen.contains(i)
                && (scope == Scope::Whole || selection.atoms.contains(i))
        })
        .collect()
}

/// Optimise only free Cartesian coordinates, evaluating energy against the entire
/// structure. This keeps fixed atoms exact, including during line searches.
struct FreeCoordinates<'a> {
    inner: DreidingObjective<'a>,
    full: Vec<f64>,
    indices: Vec<usize>,
    gradient: Vec<f64>,
    cancel: &'a AtomicBool,
}
impl FreeCoordinates<'_> {
    fn update(&mut self, x: &[f64]) -> Result<()> {
        if self.cancel.load(Ordering::Relaxed) {
            bail!("Cleanup cancelled.");
        }
        for (&i, &value) in self.indices.iter().zip(x) {
            self.full[i] = value;
        }
        Ok(())
    }
}
impl Objective for FreeCoordinates<'_> {
    fn energy(&mut self, x: &[f64]) -> Result<f64> {
        self.update(x)?;
        self.inner.energy(&self.full)
    }
    fn gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<()> {
        self.energy_gradient(x, grad).map(|_| ())
    }
    fn energy_gradient(&mut self, x: &[f64], grad: &mut [f64]) -> Result<f64> {
        self.update(x)?;
        let e = self.inner.energy_gradient(&self.full, &mut self.gradient)?;
        for (g, &i) in grad.iter_mut().zip(&self.indices) {
            *g = self.gradient[i];
        }
        Ok(e)
    }
}

pub fn relax(
    mol: &Molecule,
    selection: &Selection,
    scope: Scope,
    cancel: &AtomicBool,
) -> Result<Preview> {
    if mol.pos.len() != mol.atoms.len() || mol.pos.iter().any(|p| !p.is_finite()) {
        bail!("Invalid coordinates.");
    }
    let mobile = movable(mol, selection, scope);
    if mobile.is_empty() {
        bail!("No movable atoms. Select atoms or unfreeze them first.");
    }
    // Dummies are excluded from the force field but retained, fixed, in the editor.
    let real: Vec<_> = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, s)| !is_dummy(s))
        .map(|(i, _)| i)
        .collect();
    let mut perceived = mol.clone();
    perceived.atoms = real.iter().map(|&i| mol.atoms[i].clone()).collect();
    perceived.pos = real.iter().map(|&i| mol.pos[i]).collect();
    perceived.recompute_bonds(1.2, 2.5);
    let topology =
        DreidingTopology::build(&perceived).map_err(|e| anyhow::anyhow!(e.to_string()))?;
    let full: Vec<_> = perceived
        .pos
        .iter()
        .flat_map(|p| [p.x as f64, p.y as f64, p.z as f64])
        .map(|x| x / BOHR_TO_ANGSTROM)
        .collect();
    let indices: Vec<_> = real
        .iter()
        .enumerate()
        .filter(|(_, original)| mobile.contains(original))
        .flat_map(|(i, _)| [3 * i, 3 * i + 1, 3 * i + 2])
        .collect();
    let x: Vec<_> = indices.iter().map(|&i| full[i]).collect();
    let mut objective = FreeCoordinates {
        inner: DreidingObjective::new(&topology),
        gradient: vec![0.0; full.len()],
        full,
        indices,
        cancel,
    };
    let initial = objective.energy(&x)?;
    if !initial.is_finite() {
        bail!("Initial force-field energy is not finite. Check overlapping atoms.");
    }
    let result = minimize_cartesian(x, &mut objective, cleanup_options())?;
    if cancel.load(Ordering::Relaxed) {
        bail!("Cleanup cancelled.");
    }
    if !result.energy.is_finite()
        || result.x.iter().any(|v| !v.is_finite())
        || result.gradient.iter().any(|v| !v.is_finite())
        || result.energy > initial + 1e-9
    {
        bail!("Cleanup did not produce a finite, lower-energy geometry. Check atom overlaps and bonding.");
    }
    objective.update(&result.x)?;
    let mut positions = mol.pos.clone();
    for (i, &original) in real.iter().enumerate() {
        if mobile.contains(&original) {
            positions[original] = Vec3::new(
                objective.full[3 * i] as f32,
                objective.full[3 * i + 1] as f32,
                objective.full[3 * i + 2] as f32,
            ) * BOHR_TO_ANGSTROM as f32;
        }
    }
    if positions.iter().any(|p| !p.is_finite()) {
        bail!("Cleanup coordinates exceed the display range.");
    }
    let report = format!(
        "{} cycles · ΔE {:+.2} kcal/mol\n{}{}",
        result.cycles,
        (result.energy - initial) * HARTREE_TO_KCAL,
        if result.converged {
            "Converged."
        } else {
            "Stopped before convergence; this is an intermediate geometry."
        },
        topology
            .radical_warning()
            .map(|w| format!("\n{w}"))
            .unwrap_or_default()
    );
    Ok(Preview {
        atoms: mol.atoms.clone(),
        original: mol.pos.clone(),
        positions,
        mobile,
        report,
    })
}

pub fn poll(mut builder: ResMut<ZMatrixBuilderState>, mol: Res<Molecule>) {
    let builder = &mut *builder;
    builder.cleanup.poll(&mol, &builder.selection);
}

#[derive(Default, Reflect, GizmoConfigGroup)]
pub struct CleanupGizmos;

pub fn configure(mut store: ResMut<GizmoConfigStore>) {
    let (cfg, _) = store.config_mut::<CleanupGizmos>();
    cfg.enabled = true;
    cfg.line.width = 2.0;
    cfg.line.perspective = false;
    cfg.depth_bias = -1.0;
}

/// Ghost geometry never modifies coordinates, picking, or calculation input.
pub fn draw(
    settings: Res<crate::settings::MolSettings>,
    mut gizmos: Gizmos<CleanupGizmos>,
    builder: Res<ZMatrixBuilderState>,
    mol: Res<Molecule>,
) {
    let frozen = Color::srgb(1.0, 0.70, 0.23);
    for &i in &builder.selection.frozen {
        if !builder.selection.hidden.contains(&i)
            && crate::scene::element_visible(&mol.atoms[i], &settings)
        {
            if let Some(&p) = mol.pos.get(i) {
                gizmos.cube(
                    Transform::from_translation(p).with_scale(Vec3::splat(
                        crate::scene::atom_display_radius(&mol.atoms[i], &settings).max(0.15) * 0.4,
                    )),
                    frozen,
                );
            }
        }
    }
    if let Some(preview) = builder
        .cleanup
        .preview
        .as_ref()
        .filter(|p| builder.cleanup.open && p.matches(&mol))
    {
        let color = Color::srgb(0.15, 0.95, 0.85);
        for &(a, b, _) in &mol.bonds {
            if a < preview.positions.len()
                && b < preview.positions.len()
                && !builder.selection.hidden.contains(&a)
                && !builder.selection.hidden.contains(&b)
                && crate::scene::element_visible(&mol.atoms[a], &settings)
                && crate::scene::element_visible(&mol.atoms[b], &settings)
            {
                gizmos.line(preview.positions[a], preview.positions[b], color);
            }
        }
        for &i in &preview.mobile {
            if !builder.selection.hidden.contains(&i)
                && crate::scene::element_visible(&mol.atoms[i], &settings)
            {
                gizmos.sphere(preview.positions[i], 0.18, color);
                gizmos.line(mol.pos[i], preview.positions[i], color);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn water() -> Molecule {
        Molecule::from_xyz("3\n\nO 0 0 0\nH 1.1 0 0\nH -0.2 1.05 0\n")
    }
    #[test]
    fn selected_cleanup_keeps_every_other_atom_exactly_fixed() {
        let mol = water();
        let mut selection = Selection::default();
        selection.atoms.insert(1);
        let p = relax(&mol, &selection, Scope::Selected, &AtomicBool::new(false)).unwrap();
        assert_eq!(p.positions[0], mol.pos[0]);
        assert_eq!(p.positions[2], mol.pos[2]);
        assert!(p.positions[1].distance(mol.pos[1]) > 0.01);
        assert!(p.matches(&mol));
    }
    #[test]
    fn whole_cleanup_honours_frozen_atoms_and_dummies() {
        let mut mol = water();
        mol.atoms.push("X".into());
        mol.pos.push(Vec3::splat(5.0));
        let mut selection = Selection::default();
        selection.frozen = [0, 2].into();
        let p = relax(&mol, &selection, Scope::Whole, &AtomicBool::new(false)).unwrap();
        assert_eq!(p.positions[0], mol.pos[0]);
        assert_eq!(p.positions[2..], mol.pos[2..]);
        assert_ne!(p.positions[1], mol.pos[1]);
    }
    #[test]
    fn invalid_or_cancelled_cleanup_cannot_be_accepted() {
        let mol = water();
        assert!(relax(
            &mol,
            &Selection::default(),
            Scope::Selected,
            &AtomicBool::new(false)
        )
        .is_err());
        assert!(relax(
            &mol,
            &Selection::default(),
            Scope::Whole,
            &AtomicBool::new(true)
        )
        .is_err());
        let unsupported = Molecule::from_xyz("2\n\nLi 0 0 0\nF 1.5 0 0\n");
        assert!(relax(
            &unsupported,
            &Selection::default(),
            Scope::Whole,
            &AtomicBool::new(false)
        )
        .is_err());
        let mut state = Cleanup::default();
        state.preview = Some(
            relax(
                &mol,
                &Selection::default(),
                Scope::Whole,
                &AtomicBool::new(false),
            )
            .unwrap(),
        );
        let mut changed = mol.clone();
        changed.pos[0].x += 0.01;
        assert!(state.accepted(&changed, &Selection::default()).is_err());
        state.cancel();
        assert!(state.accepted(&mol, &Selection::default()).is_err());
    }

    #[test]
    fn freezing_an_atom_after_preview_invalidates_that_preview() {
        let mol = water();
        let mut selection = Selection::default();
        let mut state = Cleanup::default();
        state.preview =
            Some(relax(&mol, &selection, Scope::Whole, &AtomicBool::new(false)).unwrap());
        selection.frozen.insert(1);
        assert!(state.accepted(&mol, &selection).is_err());
        state.poll(&mol, &selection);
        assert!(state.preview.is_none());
        assert!(state.error.as_deref().unwrap().contains("frozen"));
    }
}
