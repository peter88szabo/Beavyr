//! Named camera views, notes and atom labels, stored in the existing project archive.
use crate::{camera::OrbitCamera, molecule::Molecule, ui_style::ResponseExt};
use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts, EguiPrimaryContextPass};
use serde::{Deserialize, Serialize};
const OPEN: &str = "open_saved_views";
#[derive(Clone, Serialize, Deserialize)]
pub struct AtomLabel {
    pub atom: usize,
    pub text: String,
}
#[derive(Clone, Serialize, Deserialize)]
pub struct View {
    pub name: String,
    pub camera: OrbitCamera,
    pub atoms: Vec<String>,
    pub labels: Vec<AtomLabel>,
    pub note: String,
}
#[derive(Resource, Clone, Default, Serialize, Deserialize)]
pub struct Library {
    pub views: Vec<View>,
    pub atoms: Vec<String>,
    pub labels: Vec<AtomLabel>,
    pub note: String,
    pub show_labels: bool,
}
impl Library {
    pub fn validate(&self) -> Result<(), String> {
        for view in &self.views {
            let c = &view.camera;
            if !c.radius.is_finite()
                || c.radius <= 0.
                || !c.target.is_finite()
                || !c.orientation.is_finite()
                || c.orientation.length_squared() < 1e-8
            {
                return Err(format!("Invalid camera in saved view {}", view.name));
            }
            if view.labels.iter().any(|l| l.atom >= view.atoms.len()) {
                return Err("Saved view has an invalid atom-label index.".into());
            }
        }
        if self.labels.iter().any(|l| l.atom >= self.atoms.len()) {
            return Err("Invalid atom-label index.".into());
        }
        Ok(())
    }
    fn capture(
        &mut self,
        name: &str,
        camera: &OrbitCamera,
        atoms: &[String],
    ) -> Result<(), String> {
        if name.trim().is_empty() {
            return Err("Give the view a name.".into());
        }
        let mut camera = camera.clone();
        camera.rotating = false;
        camera.panning = false;
        self.views.push(View {
            name: name.trim().into(),
            camera,
            atoms: atoms.to_vec(),
            labels: if self.atoms == atoms {
                self.labels.clone()
            } else {
                vec![]
            },
            note: self.note.clone(),
        });
        Ok(())
    }
    fn restore(&mut self, index: usize, camera: &mut OrbitCamera) -> Result<(), String> {
        let view = self
            .views
            .get(index)
            .ok_or("Saved view no longer exists.")?
            .clone();
        *camera = view.camera;
        camera.orientation = camera.orientation.normalize();
        camera.rotating = false;
        camera.panning = false;
        self.atoms = view.atoms;
        self.labels = view.labels;
        self.note = view.note;
        self.show_labels = true;
        Ok(())
    }
}
#[derive(Resource)]
struct Panel {
    open: bool,
    name: String,
    atom: usize,
    label: String,
    undo: Option<(Library, OrbitCamera)>,
}
impl Default for Panel {
    fn default() -> Self {
        Self {
            open: false,
            name: String::new(),
            atom: 1,
            label: String::new(),
            undo: None,
        }
    }
}
pub struct Plugin;
impl bevy::prelude::Plugin for Plugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<Library>()
            .init_resource::<Panel>()
            .add_systems(EguiPrimaryContextPass, ui.after(crate::ui::ui_panel));
    }
}
pub fn button(ui: &mut egui::Ui) {
    if ui.button("Saved views / notes…").clicked_once() {
        ui.ctx()
            .data_mut(|d| d.insert_temp(egui::Id::new(OPEN), true));
    }
}
fn ui(
    mut contexts: EguiContexts,
    mut library: ResMut<Library>,
    mut panel: ResMut<Panel>,
    mut orbit: ResMut<OrbitCamera>,
    mol: Res<Molecule>,
    camera: Query<(&Camera, &GlobalTransform), With<crate::scene::MainCamera>>,
    windows: Query<&Window>,
    mut regions: ResMut<crate::ui::UiPanelRegions>,
) {
    let Ok(ctx) = contexts.ctx_mut() else {
        return;
    };
    if ctx
        .data_mut(|d| d.remove_temp::<bool>(egui::Id::new(OPEN)))
        .unwrap_or(false)
    {
        panel.open = true;
    }
    if panel.open {
        let mut open = true;
        let mut action = None;
        let mut remove_label = None;
        let window = egui::Window::new("Saved views and annotations")
            .open(&mut open)
            .resizable(true)
            .default_size([460., 470.])
            .show(ctx, |ui| {
                crate::ui_style::modern(ui);
                ui.horizontal(|ui| {
                    ui.add(
                        egui::TextEdit::singleline(&mut panel.name)
                            .hint_text("View name")
                            .desired_width(220.),
                    );
                    if ui.button("Save view").clicked_once() {
                        action = Some((0, 0));
                    }
                });
                if ui
                    .add_enabled(
                        panel.undo.is_some(),
                        egui::Button::new("Undo view / annotation change"),
                    )
                    .clicked_once()
                {
                    action = Some((3, 0));
                }
                egui::ScrollArea::vertical()
                    .id_salt("saved_view_list")
                    .max_height(160.)
                    .show(ui, |ui| {
                        for (i, view) in library.views.iter().enumerate() {
                            ui.horizontal(|ui| {
                                if ui.button(&view.name).clicked_once() {
                                    action = Some((1, i));
                                }
                                if ui.small_button("Delete").clicked_once() {
                                    action = Some((2, i));
                                }
                            });
                        }
                    });
                ui.separator();
                ui.label("Notes");
                ui.add(
                    egui::TextEdit::multiline(&mut library.note)
                        .desired_rows(4)
                        .desired_width(f32::INFINITY),
                );
                ui.horizontal_wrapped(|ui| {
                    ui.checkbox(&mut library.show_labels, "Show atom labels");
                    ui.label("Atom");
                    ui.add(egui::DragValue::new(&mut panel.atom).range(1..=mol.atoms.len().max(1)));
                });
                ui.horizontal(|ui| {
                    ui.add(egui::TextEdit::singleline(&mut panel.label).hint_text("Custom label"));
                    if ui.button("Set label").clicked_once() {
                        action = Some((4, 0));
                    }
                });
                if library.atoms != mol.atoms && !library.labels.is_empty() {
                    ui.label("Labels belong to a different atom list and are hidden.");
                }
                egui::ScrollArea::vertical()
                    .id_salt("saved_atom_labels")
                    .max_height(150.)
                    .show(ui, |ui| {
                        for (i, label) in library.labels.iter().enumerate() {
                            ui.horizontal(|ui| {
                                ui.label(format!("{}: {}", label.atom + 1, label.text));
                                if ui.small_button("Remove").clicked_once() {
                                    remove_label = Some(i);
                                }
                            });
                        }
                    });
                if ui.button("Copy notes").clicked_once() {
                    ctx.copy_text(library.note.clone());
                }
                ui.small("Views, labels and notes are saved with the Beavyr project.");
            });
        panel.open = open;
        if let Some(w) = window {
            let r = w.response.rect;
            regions
                .windows
                .push(Rect::new(r.min.x, r.min.y, r.max.x, r.max.y));
        }
        if let Some((kind, index)) = action {
            if kind == 3 {
                if let Some((old, camera)) = panel.undo.take() {
                    panel.undo = Some((library.clone(), orbit.clone()));
                    *library = old;
                    *orbit = camera;
                }
            } else {
                panel.undo = Some((library.clone(), orbit.clone()));
                let result = match kind {
                    0 => library.capture(&panel.name, &orbit, &mol.atoms),
                    1 => library.restore(index, &mut orbit),
                    2 => {
                        library.views.remove(index);
                        Ok(())
                    }
                    4 => {
                        if panel.atom == 0
                            || panel.atom > mol.atoms.len()
                            || panel.label.trim().is_empty()
                        {
                            Err("Choose an atom and enter a label.".into())
                        } else {
                            if library.atoms != mol.atoms {
                                library.labels.clear();
                                library.atoms = mol.atoms.clone();
                            }
                            library.labels.retain(|l| l.atom != panel.atom - 1);
                            library.labels.push(AtomLabel {
                                atom: panel.atom - 1,
                                text: panel.label.trim().into(),
                            });
                            library.show_labels = true;
                            Ok(())
                        }
                    }
                    _ => Ok(()),
                };
                if let Err(error) = result {
                    crate::runtime_warning::push(ctx, "Saved views", &error);
                }
            }
        }
        if let Some(i) = remove_label {
            panel.undo = Some((library.clone(), orbit.clone()));
            library.labels.remove(i);
        }
    }
    if library.show_labels && library.atoms == mol.atoms {
        if let (Ok((camera, transform)), Ok(window)) = (camera.single(), windows.single()) {
            let scale = window.scale_factor() as f32;
            let painter = ctx.layer_painter(egui::LayerId::new(
                egui::Order::Background,
                egui::Id::new("custom_atom_labels"),
            ));
            for label in &library.labels {
                if let Some(&position) = mol.pos.get(label.atom) {
                    if let Ok(p) = camera.world_to_viewport(transform, position) {
                        let p = egui::pos2(p.x / scale, p.y / scale);
                        if !regions.covers(Vec2::new(p.x, p.y)) {
                            painter.text(
                                p + egui::vec2(1., 1.),
                                egui::Align2::LEFT_BOTTOM,
                                &label.text,
                                egui::FontId::proportional(15.),
                                egui::Color32::BLACK,
                            );
                            painter.text(
                                p,
                                egui::Align2::LEFT_BOTTOM,
                                &label.text,
                                egui::FontId::proportional(15.),
                                egui::Color32::WHITE,
                            );
                        }
                    }
                }
            }
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn views_restore_camera_notes_and_labels_and_serialize() {
        let mut l = Library::default();
        let mut c = OrbitCamera::default();
        c.radius = 7.;
        c.target = Vec3::new(1., 2., 3.);
        l.atoms = vec!["Fe".into()];
        l.labels.push(AtomLabel {
            atom: 0,
            text: "Metal centre".into(),
        });
        l.note = "Ferrocene view".into();
        l.capture("Front", &c, &["Fe".into()]).unwrap();
        c.radius = 90.;
        l.note.clear();
        l.restore(0, &mut c).unwrap();
        assert_eq!(c.radius, 7.);
        assert_eq!(l.note, "Ferrocene view");
        let restored: Library = ron::from_str(&ron::to_string(&l).unwrap()).unwrap();
        assert_eq!(restored.views[0].labels[0].text, "Metal centre");
        assert!(restored.validate().is_ok());
    }
    #[test]
    fn bad_camera_or_annotation_is_rejected() {
        let mut l = Library::default();
        l.capture("View", &OrbitCamera::default(), &["C".into()])
            .unwrap();
        l.views[0].camera.radius = f32::NAN;
        assert!(l.validate().is_err());
        assert!(Library::default()
            .capture(" ", &OrbitCamera::default(), &[])
            .is_err());
    }
}
