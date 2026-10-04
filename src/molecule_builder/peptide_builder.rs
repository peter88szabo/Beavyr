//! Independent peptide window. It never changes a structure until Insert.
use bevy::prelude::*;
use bevy_egui::egui;

use super::amino_acids::{AminoAcid, SideChainState, TerminalState};
use super::builder_ui::{self, BuilderChange, EditorRotateState, ZMatrixBuilderState};
use super::peptide::{self, Attachment, Build, End, Recipe, Residue, Shape};
use crate::{molecule::Molecule, settings::MolSettings, trajectory::TrajectoryState};

const ACCENT: egui::Color32 = egui::Color32::from_rgb(104, 174, 255);

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    New,
    Extend,
    Edit,
}

#[derive(Clone)]
struct Target {
    atoms: Vec<String>,
    positions: Vec<Vec3>,
    attachment: Attachment,
}

#[derive(Clone)]
struct Preview {
    build: Build,
    recipe: Recipe,
    mode: Mode,
    target: Option<(usize, End)>,
    atoms: Vec<String>,
    positions: Vec<Vec3>,
    topology: Option<super::topology::Topology>,
    fit: bool,
}

#[derive(Clone)]
pub struct State {
    pub open: bool,
    edit_residue: Option<u64>,
    link_residue: Option<u64>,
    replacement: Residue,
    fit: bool,
    backbone: BackboneControls,
    rotamers: Option<RotamerCache>,
    neighbor_radius: f32,
    recipe: Recipe,
    mode: Mode,
    end: End,
    target: Option<Target>,
    search: String,
    sequence_input: String,
    selected_row: Option<usize>,
    preview: Option<Preview>,
    message: String,
    error: Option<String>,
    yaw: f32,
    pitch: f32,
    zoom: f32,
    show_h: bool,
}

impl Default for State {
    fn default() -> Self {
        Self {
            open: false,
            edit_residue: None,
            link_residue: None,
            replacement: Residue::new(AminoAcid::Alanine),
            fit: true,
            backbone: BackboneControls::default(),
            rotamers: None,
            neighbor_radius: 6.0,
            recipe: Recipe::default(),
            mode: Mode::New,
            end: End::C,
            target: None,
            search: String::new(),
            sequence_input: String::new(),
            selected_row: None,
            preview: None,
            message: String::new(),
            error: None,
            yaw: 0.3,
            pitch: -0.3,
            zoom: 1.0,
            show_h: false,
        }
    }
}

fn state_name(state: SideChainState) -> &'static str {
    match state {
        SideChainState::Neutral => "Neutral",
        SideChainState::Protonated => "+1 · protonated",
        SideChainState::Deprotonated => "−1 · deprotonated",
        SideChainState::HistidineDelta => "His δ · ND1-H",
        SideChainState::HistidineEpsilon => "His ε · NE2-H",
    }
}

fn move_button(ui: &mut egui::Ui, enabled: bool, up: bool) -> bool {
    let label = if up {
        "Move residue toward N terminus"
    } else {
        "Move residue toward C terminus"
    };
    let response = ui
        .add_enabled(
            enabled,
            egui::Button::new("").min_size(egui::vec2(24.0, 20.0)),
        )
        .on_hover_text(label);
    response.widget_info(|| egui::WidgetInfo::labeled(egui::WidgetType::Button, enabled, label));
    let center = response.rect.center();
    let sign = if up { -1.0 } else { 1.0 };
    let color = if enabled {
        ui.visuals().text_color()
    } else {
        ui.visuals().weak_text_color()
    };
    ui.painter().add(egui::Shape::line(
        vec![
            center + egui::vec2(-4.0, -2.0 * sign),
            center + egui::vec2(0.0, 2.0 * sign),
            center + egui::vec2(4.0, -2.0 * sign),
        ],
        egui::Stroke::new(1.5, color),
    ));
    response.clicked()
}

impl State {
    fn target_key(&self) -> Option<(usize, End)> {
        (self.mode == Mode::Extend)
            .then(|| {
                self.target
                    .as_ref()
                    .map(|t| (t.attachment.atom, t.attachment.end))
            })
            .flatten()
    }

    fn preview_current(&self, mol: &Molecule) -> bool {
        self.preview.as_ref().is_some_and(|p| {
            p.recipe == self.recipe
                && p.mode == self.mode
                && p.target == self.target_key()
                && p.atoms == mol.atoms
                && p.positions == mol.pos
                && p.topology == mol.topology
                && p.fit == self.fit
        })
    }

    fn make_preview(&mut self, mol: &Molecule) -> Result<(), String> {
        let host = if self.mode == Mode::Extend {
            let target = self
                .target
                .as_ref()
                .ok_or("Select a terminal atom, then click Use selected atom.")?;
            if target.atoms != mol.atoms
                || target.positions != mol.pos
                || target.attachment.end != self.end
            {
                return Err("The attachment site changed. Use selected atom again.".into());
            }
            Some((mol, &target.attachment))
        } else {
            None
        };
        let mut build = peptide::build(&self.recipe, host)?;
        if self.fit {
            // Fit only newly built residues; never move the retained host.
            let ids: Vec<_> = build.molecule.topology.as_ref().unwrap().residues.iter()
                .filter(|r| build.moving.iter().any(|&i| build.molecule.topology.as_ref().unwrap().atoms[i].residue == Some(r.id)))
                .map(|r| r.id).collect();
            for id in ids { super::peptide_edit::pack_sidechains(&mut build.molecule,Some(id))?; }
            build.clashes = peptide::clashes(&build.molecule,build.moving.first().copied().unwrap_or(0));
        }
        if self.mode == Mode::New {
            build = peptide::beside(build, mol);
        }
        self.preview = Some(Preview {
            build,
            recipe: self.recipe.clone(),
            mode: self.mode,
            target: self.target_key(),
            atoms: mol.atoms.clone(),
            positions: mol.pos.clone(),
            topology: mol.topology.clone(),
            fit: self.fit,
        });
        self.error = None;
        self.message.clear();
        Ok(())
    }

    fn add(&mut self, amino: AminoAcid) {
        if self.recipe.residues.len() >= peptide::MAX_RESIDUES {
            return;
        }
        let index = if self.end == End::C {
            self.recipe.residues.len()
        } else {
            0
        };
        self.recipe.residues.insert(index, Residue::new(amino));
        self.selected_row = Some(index);
        self.error = None;
    }
}

/// Render independently of the Molecular Editor's visibility; return its rect
/// so viewport clicks on this floating window never pick atoms behind it.
pub fn window(
    ctx: &egui::Context,
    zmat: &mut ZMatrixBuilderState,
    editor: &mut EditorRotateState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
    trajectory: &mut TrajectoryState,
) -> (Option<egui::Rect>, BuilderChange) {
    if !zmat.peptide.open {
        return (None, BuilderChange::None);
    }
    let mut state = std::mem::take(&mut zmat.peptide);
    let mut open = state.open;
    let mut insert = false;
    let ready = !trajectory.playing
        && !zmat.confirm_clear
        && !zmat.add_atom_active
        && !editor.active
        && [editor as &EditorRotateState, &zmat.placed_editor]
            .iter()
            .all(|e| e.base_pos.as_ref().is_none_or(|p| p == &mol.pos));
    if state.target.as_ref().is_some_and(|t| {
        t.atoms != mol.atoms || t.positions != mol.pos || t.attachment.end != state.end
    }) {
        state.target = None;
    }
    let response = egui::Window::new("Peptide Builder")
        .id(egui::Id::new("peptide_builder_window"))
        .open(&mut open).default_size([800.0, 650.0]).min_width(640.0)
         .resizable(true).collapsible(false)
        .max_height((ctx.content_rect().height() - 24.0).max(300.0)).vscroll(true)
        .show(ctx, |ui| {
            ui.spacing_mut().item_spacing = egui::vec2(8.0, 7.0);
            ui.spacing_mut().button_padding = egui::vec2(10.0, 6.0);
            let visuals = &mut ui.style_mut().visuals;
            visuals.selection.bg_fill = egui::Color32::from_rgb(35, 75, 117);
            visuals.selection.stroke = egui::Stroke::new(1.0, ACCENT);
            for widget in [&mut visuals.widgets.inactive, &mut visuals.widgets.hovered, &mut visuals.widgets.active] {
                widget.corner_radius = egui::CornerRadius::same(6);
            }
            ui.horizontal(|ui| {
                ui.selectable_value(&mut state.mode, Mode::New, "New peptide");
                ui.selectable_value(&mut state.mode, Mode::Extend, "Extend selected terminus");
                ui.selectable_value(&mut state.mode, Mode::Edit, "Edit residues");
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let count=if state.mode==Mode::Edit {mol.topology.as_ref().map(|t|t.peptide_residue_count()).unwrap_or(0)} else {state.recipe.residues.len()};
                    ui.weak(format!("{count} residues"));
                });
            });
            if state.mode == Mode::Edit {
                edit_controls(ui, &mut state, mol, zmat.selected_index, &zmat.selection.frozen, ready);
                let current = state.preview_current(mol);
                insert = ui.add_enabled(ready && current, egui::Button::new("Apply preview")).clicked();
                if let Some(error) = &state.error { ui.colored_label(egui::Color32::LIGHT_RED,error); }
                if !state.message.is_empty() { ui.label(&state.message); }
                if let Some(preview) = &state.preview {
                    if !current { ui.colored_label(egui::Color32::YELLOW,"Preview is out of date — preview the edit again."); }
                    ui.label(&preview.build.sequence);
                    if let Some(q) = preview.build.charge { ui.small(format!("Formal charge {q:+}")); }
                    if preview.build.clashes > 0 { ui.colored_label(egui::Color32::YELLOW,format!("{} close contacts",preview.build.clashes)); }
                    draw_preview(ui,&preview.build.molecule,&preview.build.moving,&mut state.yaw,&mut state.pitch,&mut state.zoom,state.show_h);
                }
                return;
            }
            ui.horizontal(|ui| {
                ui.label(if state.mode == Mode::New { "Next residue at" } else { "Extend host at" });
                if ui.selectable_value(&mut state.end, End::C, "C terminus").changed()
                    | ui.selectable_value(&mut state.end, End::N, "N terminus").changed() {
                    state.target = None;
                }
            });
            if state.mode == Mode::Extend {
                ui.horizontal_wrapped(|ui| {
                    let selected = zmat.selected_index;
                    let label = selected.filter(|&i| i < mol.atoms.len())
                        .map(|i| format!("Use selected atom · {}{}", mol.atoms[i], i + 1))
                        .unwrap_or_else(|| "Use selected atom".into());
                    if ui.add_enabled(selected.is_some() && ready, egui::Button::new(label)).clicked() {
                        match peptide::attachment(mol, selected.unwrap(), state.end) {
                            Ok(attachment) => {
                                state.target = Some(Target { atoms: mol.atoms.clone(), positions: mol.pos.clone(), attachment });
                                state.error = None;
                            }
                            Err(error) => state.error = Some(error),
                        }
                    }
                    if let Some(target) = &state.target {
                        ui.colored_label(ACCENT, format!("{}{} connects to {} of new sequence",
                            mol.atoms[target.attachment.atom], target.attachment.atom + 1,
                            if state.end == End::C { "N" } else { "COOH" }));
                    } else {
                        ui.weak(if state.end == End::C { "Pick a free COOH/COO- carbon." } else { "Pick a free amino nitrogen." });
                    }
                });
            } else if !mol.atoms.is_empty() {
                ui.small("Adds a separate peptide beside the displayed structure.");
            }
            ui.separator();
            let body_height = (ui.available_height() - 415.0).max(165.0);
            ui.columns(2, |columns| {
                let ui = &mut columns[0];
                ui.add(egui::TextEdit::singleline(&mut state.search).hint_text("Find amino acid").desired_width(f32::INFINITY));
                egui::ScrollArea::vertical().id_salt("peptide_palette").max_height(body_height).auto_shrink([false, false]).show(ui, |ui| {
                    let query = state.search.trim().to_ascii_lowercase();
                    let amino: Vec<_> = AminoAcid::ALL.into_iter().filter(|a|
                        query.is_empty() || a.name().to_ascii_lowercase().contains(&query)
                            || a.three_letter().to_ascii_lowercase().contains(&query)
                            || a.one_letter().to_ascii_lowercase().to_string() == query).collect();
                    let tile_width = ((ui.available_width() - 18.0) / 4.0).floor();
                    egui::Grid::new("peptide_tiles").spacing([6.0, 6.0]).show(ui, |ui| {
                        for (i, acid) in amino.iter().enumerate() {
                            let title = format!("{}  {}", acid.three_letter(), acid.one_letter());
                            if ui.add_enabled(state.recipe.residues.len() < peptide::MAX_RESIDUES,
                                egui::Button::new(egui::RichText::new(title).strong()).min_size(egui::vec2(tile_width, 28.0)))
                                .on_hover_text(acid.name()).clicked() { state.add(*acid); }
                            if i % 4 == 3 { ui.end_row(); }
                        }
                    });
                });
                let ui = &mut columns[1];
                ui.horizontal(|ui| {
                    ui.strong("N to C sequence");
                    if ui.small_button("Clear").clicked() {
                        state.recipe.residues.clear();
                        state.selected_row = None;
                    }
                });
                if state.mode == Mode::Extend {
                    ui.small(if state.end == End::C { "Host — [new sequence]" } else { "[new sequence] — Host" });
                }
                egui::ScrollArea::vertical().id_salt("peptide_sequence_rows").max_height(body_height - 36.0)
                    .auto_shrink([false, false]).show(ui, |ui| {
                    if state.recipe.residues.is_empty() { ui.weak("Choose the first amino acid from the tiles."); }
                    let mut movement = None;
                    let mut remove = None;
                    for (i, residue) in state.recipe.residues.iter().enumerate() {
                        ui.horizontal(|ui| {
                            if ui.selectable_label(state.selected_row == Some(i), format!("{:>2}  {}", i + 1, residue.amino.three_letter())).clicked() {
                                state.selected_row = Some(i);
                            }
                            ui.small(state_name(residue.side_chain));
                            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                if ui.small_button("×").on_hover_text("Remove residue").clicked() { remove = Some(i); }
                                if move_button(ui, i + 1 < state.recipe.residues.len(), false) { movement = Some((i, i + 1)); }
                                if move_button(ui, i > 0, true) { movement = Some((i, i.saturating_sub(1))); }
                            });
                        });
                    }
                    if let Some(i) = remove {
                        state.recipe.residues.remove(i);
                        state.selected_row = None;
                    } else if let Some((a, b)) = movement {
                        state.recipe.residues.swap(a, b);
                        state.selected_row = Some(b);
                    }
                });
            });
            if let Some(residue) = state.selected_row.and_then(|i| state.recipe.residues.get_mut(i))
                .filter(|r| r.amino.templates().iter().filter(|t| t.termini == TerminalState::Neutral).count() > 1) {
                ui.horizontal(|ui| {
                    ui.label(format!("{} side chain", residue.amino.three_letter()));
                    egui::ComboBox::from_id_salt("peptide_residue_state")
                        .selected_text(state_name(residue.side_chain)).show_ui(ui, |ui| {
                        for template in residue.amino.templates().iter().filter(|t| t.termini == TerminalState::Neutral) {
                            ui.selectable_value(&mut residue.side_chain, template.side_chain, state_name(template.side_chain));
                        }
                    });
                });
            }
            egui::CollapsingHeader::new("Paste sequence").id_salt("peptide_paste").show(ui, |ui| {
                ui.add(egui::TextEdit::multiline(&mut state.sequence_input).desired_rows(2)
                    .hint_text("AGSK… or Ala-Gly-Ser-Lys (N to C); FASTA accepted").desired_width(f32::INFINITY));
                if ui.button("Replace sequence list").clicked() {
                    match peptide::parse_sequence(&state.sequence_input) {
                        Ok(residues) => { state.recipe.residues = residues; state.selected_row = None; state.error = None; }
                        Err(error) => state.error = Some(error),
                    }
                }
                ui.small("A bare ALA means A–L–A. Separate three-letter names with hyphens.");
            });
            ui.separator();
            ui.horizontal_wrapped(|ui| {
                ui.label("Backbone");
                ui.selectable_value(&mut state.recipe.shape, Shape::Extended, "Extended");
                ui.selectable_value(&mut state.recipe.shape, Shape::Alpha, "α-helix");
                ui.selectable_value(&mut state.recipe.shape, Shape::Custom, "Custom φ/ψ");
                if state.recipe.shape == Shape::Custom {
                    ui.add(egui::DragValue::new(&mut state.recipe.phi).prefix("φ ").suffix("°").range(-180.0..=180.0));
                    ui.add(egui::DragValue::new(&mut state.recipe.psi).prefix("ψ ").suffix("°").range(-180.0..=180.0));
                }
            });
            ui.horizontal_wrapped(|ui| {
                ui.label("New free ends");
                ui.selectable_value(&mut state.recipe.termini, TerminalState::Neutral, "NH₂ / COOH");
                ui.selectable_value(&mut state.recipe.termini, TerminalState::Zwitterionic, "NH₃+ / COO-");
                ui.weak("Trans peptide bonds");
                ui.checkbox(&mut state.fit,"Fit side chains").on_hover_text("Sample side-chain torsions to reduce local overlaps; preserve rings and the backbone.");
            });
            if state.recipe.residues.iter().any(|r| r.amino == AminoAcid::Proline) {
                ui.small("Proline keeps its ring geometry; its φ is constrained. Its free amine is NH / NH₂+.");
            }
            let current = state.preview_current(mol);
            ui.horizontal(|ui| {
                if ui.add_enabled(ready && !state.recipe.residues.is_empty(), egui::Button::new("Preview peptide")).clicked() {
                    if let Err(error) = state.make_preview(mol) { state.error = Some(error); }
                }
                insert = ui.add_enabled(ready && current, egui::Button::new(egui::RichText::new("Insert peptide").strong()).fill(egui::Color32::from_rgb(35, 75, 117)))
                    .on_hover_text("One undoable edit; the new atoms are ready for the existing Fragment Editor.").clicked();
                ui.checkbox(&mut state.show_h, "Show H");
            });
            if !ready { ui.weak("Pause playback and finish any active fragment adjustment or atom picking first."); }
            if let Some(error) = &state.error { ui.colored_label(egui::Color32::LIGHT_RED, error); }
            if !state.message.is_empty() { ui.colored_label(ACCENT, &state.message); }
            let is_current = state.preview_current(mol);
            if let Some(preview) = &state.preview {
                if !is_current { ui.colored_label(egui::Color32::YELLOW, "Preview is out of date — preview again before inserting."); }
                let build = &preview.build;
                ui.horizontal_wrapped(|ui| {
                    ui.small(format!("{} new atoms", build.moving.len()));
                    if let Some(charge) = build.charge { ui.small(format!("Charge {charge:+}")); }
                    if build.removed_host_atoms > 0 { ui.small(format!("{} leaving atoms removed from host", build.removed_host_atoms)); }
                    if build.clashes > 0 {
                        ui.colored_label(egui::Color32::YELLOW, format!("{} close contacts — adjust or relax after insertion", build.clashes));
                    }
                });
                draw_preview(ui, &build.molecule, &build.moving, &mut state.yaw, &mut state.pitch, &mut state.zoom, state.show_h);
            } else {
                let height = (ui.available_height() - 30.0).clamp(80.0, 140.0);
                let (rect, _) = ui.allocate_exact_size(egui::vec2(ui.available_width(), height), egui::Sense::hover());
                ui.painter().rect_filled(rect, 8, ui.visuals().extreme_bg_color);
                ui.painter().text(rect.center(), egui::Align2::CENTER_CENTER, "Build a sequence, then preview it here", egui::FontId::proportional(14.0), ui.visuals().weak_text_color());
            }
            ui.small("Starting geometry from optimized monomers. The assembled chain still needs relaxation.");
        });
    let mut change = BuilderChange::None;
    if insert && state.preview_current(mol) {
        let build = state.preview.as_ref().unwrap().build.clone();
        match builder_ui::insert_peptide(build, editor, zmat, mol, settings, trajectory) {
            Ok(result) => {
                // A separate new chain may lie outside the old camera frame.
                // Terminal extension keeps the user's current view.
                change = if state.mode == Mode::New {
                    BuilderChange::Replaced
                } else {
                    result
                };
                state.preview = None;
                state.target = None;
                state.message = if state.mode == Mode::Edit {
                    "Applied. Undo is available in Molecular Editor."
                } else {
                    "Inserted. Undo is available in Molecular Editor; Fragment Editor controls the new atoms."
                }.into();
                state.error = None;
            }
            Err(error) => state.error = Some(error),
        }
    }
    state.open = open;
    zmat.peptide = state;
    (response.map(|r| r.response.rect), change)
}

fn draw_preview(
    ui: &mut egui::Ui,
    mol: &Molecule,
    moving: &[usize],
    yaw: &mut f32,
    pitch: &mut f32,
    zoom: &mut f32,
    show_h: bool,
) {
    let (rect, response) = ui.allocate_exact_size(
        egui::vec2(
            ui.available_width(),
            (ui.available_height() - 30.0).clamp(80.0, 155.0),
        ),
        egui::Sense::drag(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 8, ui.visuals().extreme_bg_color);
    if response.dragged() {
        let delta = ui.input(|i| i.pointer.delta());
        *yaw += delta.x * 0.01;
        *pitch += delta.y * 0.01;
    }
    if response.hovered() {
        *zoom = (*zoom * ui.input(|i| (i.smooth_scroll_delta.y * 0.003).exp())).clamp(0.2, 8.0);
    }
    let visible: Vec<_> = (0..mol.atoms.len())
        .filter(|&i| show_h || mol.atoms[i] != "H")
        .collect();
    if visible.is_empty() {
        return;
    }
    let center = visible.iter().map(|&i| mol.pos[i]).sum::<Vec3>() / visible.len() as f32;
    let rotation = Quat::from_rotation_x(*pitch) * Quat::from_rotation_y(*yaw);
    let positions: Vec<_> = mol.pos.iter().map(|p| rotation * (*p - center)).collect();
    let mut span = Vec2::splat(0.1);
    for &i in &visible {
        span = span.max(positions[i].truncate().abs());
    }
    let scale = ((rect.width() - 28.0) / (2.0 * span.x))
        .min((rect.height() - 28.0) / (2.0 * span.y))
        * *zoom;
    let project = |i: usize| rect.center() + egui::vec2(positions[i].x, -positions[i].y) * scale;
    for &(a, b, _) in &mol.bonds {
        if !show_h && (mol.atoms[a] == "H" || mol.atoms[b] == "H") {
            continue;
        }
        painter.line_segment(
            [project(a), project(b)],
            egui::Stroke::new(1.5, egui::Color32::from_gray(120)),
        );
    }
    let mut order = visible;
    order.sort_by(|&a, &b| positions[a].z.total_cmp(&positions[b].z));
    for i in order {
        let color = match mol.atoms[i].as_str() {
            "N" => egui::Color32::from_rgb(95, 155, 255),
            "O" => egui::Color32::from_rgb(245, 102, 113),
            "S" | "Se" => egui::Color32::from_rgb(239, 200, 80),
            "H" => egui::Color32::from_gray(225),
            _ => egui::Color32::from_rgb(151, 176, 185),
        };
        // Build::moving is sorted; keep redraw work linearithmic for long chains.
        let color = if moving.binary_search(&i).is_ok() {
            color
        } else {
            color.gamma_multiply(0.4)
        };
        painter.circle_filled(
            project(i),
            if mol.atoms[i] == "H" { 2.0 } else { 3.8 },
            color,
        );
    }
    painter.text(
        rect.left_bottom() + egui::vec2(8.0, -6.0),
        egui::Align2::LEFT_BOTTOM,
        "Drag to rotate · scroll to zoom",
        egui::FontId::proportional(11.0),
        ui.visuals().weak_text_color(),
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn floating_window_renders_and_previews_without_the_molecular_editor() {
        let ctx = egui::Context::default();
        let mut zmat = ZMatrixBuilderState::default();
        zmat.peptide.open = true;
        zmat.peptide.recipe.residues = peptide::parse_sequence("AGS").unwrap();
        zmat.peptide.selected_row = Some(2);
        let mut mol = Molecule::empty();
        zmat.peptide.make_preview(&mol).unwrap();
        let mut editor = EditorRotateState::default();
        let mut settings = MolSettings::default();
        let mut trajectory = TrajectoryState::default();
        let mut rect = None;
        for _ in 0..3 {
            let mut output = ctx.run_ui(
                egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(1280.0, 720.0),
                    )),
                    ..Default::default()
                },
                |ui| {
                    let (area, change) = window(
                        ui.ctx(),
                        &mut zmat,
                        &mut editor,
                        &mut mol,
                        &mut settings,
                        &mut trajectory,
                    );
                    rect = area;
                    assert!(!change.happened());
                },
            );
            output.textures_delta.clear();
        }
        assert!(zmat.peptide.open && zmat.peptide.preview_current(&mol));
        assert!(mol.atoms.is_empty(), "preview must not mutate the host");
        assert!(rect.unwrap().bottom() <= 720.0);
    }

    #[test]
    fn n_growth_prepends_and_preview_invalidates_on_recipe_or_host_edits() {
        let mut state = State::default();
        state.add(AminoAcid::Alanine);
        state.end = End::N;
        state.add(AminoAcid::Glycine);
        assert_eq!(state.recipe.sequence(), "GA");
        let mut mol = Molecule::empty();
        state.make_preview(&mol).unwrap();
        assert!(state.preview_current(&mol));
        state.recipe.residues[0].amino = AminoAcid::Serine;
        assert!(!state.preview_current(&mol));
        state.make_preview(&mol).unwrap();
        mol.atoms.push("H".into());
        mol.pos.push(Vec3::ZERO);
        assert!(!state.preview_current(&mol));
    }
}

fn edit_controls(
    ui: &mut egui::Ui,
    state: &mut State,
    mol: &Molecule,
    selected: Option<usize>,
    frozen: &std::collections::BTreeSet<usize>,
    ready: bool,
) {
    use super::peptide_edit as edit;
    let before = (
        state.edit_residue,
        state.link_residue,
        state.replacement.clone(),
        state.fit,
    );
    let Ok(t) = edit::metadata(mol) else {
        ui.weak("Build a peptide or open a PDB, mmCIF or Beavyr project to edit residues.");
        return;
    };
    let residues: Vec<_> = t
        .residues
        .iter()
        .filter(|r| edit::template(&r.template).is_some())
        .collect();
    if !residues.iter().any(|r| Some(r.id) == state.edit_residue) {
        state.edit_residue = residues.first().map(|r| r.id);
    }
    let label = |id: Option<u64>| {
        residues
            .iter()
            .find(|r| Some(r.id) == id)
            .map(|r| residue_label(r))
            .unwrap_or_else(|| "Select residue".into())
    };
    ui.horizontal(|ui| {
        egui::ComboBox::from_id_salt("edit_residue")
            .selected_text(label(state.edit_residue))
            .show_ui(ui, |ui| {
                for r in &residues {
                    ui.selectable_value(&mut state.edit_residue, Some(r.id), label(Some(r.id)));
                }
            });
        if ui
            .add_enabled(selected.is_some(), egui::Button::new("From selected atom"))
            .clicked()
        {
            state.edit_residue = selected
                .and_then(|i| t.atoms.get(i))
                .and_then(|a| a.residue);
        }
        if let Some(q) = t.charge() {
            ui.small(format!("Charge {q:+}"));
        }
        if ui.button("Export PDB…").clicked() {
            match edit::pdb(mol) {
                Ok(text) => {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("PDB", &["pdb"])
                        .set_file_name("peptide.pdb")
                        .save_file()
                    {
                        match std::fs::write(path, text) {
                            Ok(()) => {
                                state.message = "Exported residue names and connectivity.".into()
                            }
                            Err(e) => state.error = Some(e.to_string()),
                        }
                    }
                }
                Err(e) => state.error = Some(e),
            }
        }
    });

    let mut result = None;
    ui.add_enabled_ui(ready, |ui| {
        ui.horizontal_wrapped(|ui| {
            ui.label("Replace with");
            let before = state.replacement.amino;
            egui::ComboBox::from_id_salt("replacement_acid")
                .selected_text(before.three_letter())
                .show_ui(ui, |ui| {
                    for a in AminoAcid::ALL {
                        ui.selectable_value(&mut state.replacement.amino, a, a.three_letter());
                    }
                });
            if before != state.replacement.amino {
                state.replacement = Residue::new(state.replacement.amino);
            }
            egui::ComboBox::from_id_salt("replacement_state")
                .selected_text(state_name(state.replacement.side_chain))
                .show_ui(ui, |ui| {
                    for t in state
                        .replacement
                        .amino
                        .templates()
                        .iter()
                        .filter(|t| t.termini == TerminalState::Neutral)
                    {
                        ui.selectable_value(
                            &mut state.replacement.side_chain,
                            t.side_chain,
                            state_name(t.side_chain),
                        );
                    }
                });
            if ui.button("Preview replacement").clicked() {
                result = state
                    .edit_residue
                    .map(|r| edit::replace(mol, r, &state.replacement, state.fit));
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.checkbox(&mut state.fit, "Fit side chains");
            if ui.button("Preview side-chain fit").clicked() {
                result = Some(edit::pack(mol, state.edit_residue));
            }
            if ui.button("Preview ACE cap (N)").clicked() {
                result = state.edit_residue.map(|r| edit::cap(mol, r, End::N));
            }
            if ui.button("Preview NME cap (C)").clicked() {
                result = state.edit_residue.map(|r| edit::cap(mol, r, End::C));
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.label("Disulfide partner");
            egui::ComboBox::from_id_salt("disulfide_partner")
                .selected_text(label(state.link_residue))
                .show_ui(ui, |ui| {
                    for r in &residues {
                        if edit::template(&r.template).unwrap().amino_acid == AminoAcid::Cysteine {
                            ui.selectable_value(
                                &mut state.link_residue,
                                Some(r.id),
                                label(Some(r.id)),
                            );
                        }
                    }
                });
            if ui
                .button("Preview S–S link")
                .on_hover_text(
                    "Requires two positioned cysteine sulfurs; checks distance and bond angles.",
                )
                .clicked()
            {
                result = Some(match (state.edit_residue, state.link_residue) {
                    (Some(a), Some(b)) => edit::disulfide(mol, a, b),
                    _ => Err("Select two cysteine residues.".into()),
                });
            }
        });
    });
    if let Some(operation) = conformation_controls(ui, state, mol, ready) {
        result = Some(operation);
    }
    if before
        != (
            state.edit_residue,
            state.link_residue,
            state.replacement.clone(),
            state.fit,
        )
    {
        state.preview = None;
    }
    if let Some(result) = result {
        let result = result.and_then(|build| {
            if frozen.iter().any(|&i| {
                build
                    .host_mapping
                    .get(i)
                    .copied()
                    .flatten()
                    .is_none_or(|j| build.molecule.pos[j] != mol.pos[i])
            }) {
                Err("This edit would move or remove frozen atoms. Unfreeze them first.".into())
            } else {
                Ok(build)
            }
        });
        match result {
            Ok(build) => {
                state.error = None;
                state.message.clear();
                state.preview = Some(Preview {
                    build,
                    recipe: state.recipe.clone(),
                    mode: state.mode,
                    target: state.target_key(),
                    atoms: mol.atoms.clone(),
                    positions: mol.pos.clone(),
                    topology: mol.topology.clone(),
                    fit: state.fit,
                });
            }
            Err(e) => {
                state.preview = None;
                state.error = Some(e);
            }
        }
    }
}

#[derive(Clone)]
struct BackboneControls {
    loaded: Option<u64>,
    positions: Vec<Vec3>,
    topology: Option<super::topology::Topology>,
    angles: [f32; 3],
    available: [bool; 3],
    end: Option<u64>,
    preset: super::peptide_conformation::Preset,
    points: Vec<(u64, f32, f32)>,
}
impl Default for BackboneControls {
    fn default() -> Self {
        Self {
            loaded: None,
            positions: vec![],
            topology: None,
            angles: [0.0, 0.0, 180.0],
            available: [false; 3],
            end: None,
            preset: super::peptide_conformation::Preset::Helix,
            points: vec![],
        }
    }
}
#[derive(Clone)]
struct RotamerCache {
    choices: super::peptide_rotamers::Choices,
    source: Vec<Vec3>,
    topology: Option<super::topology::Topology>,
    index: usize,
}
fn residue_label(r: &super::topology::ResidueInfo) -> String {
    let name = super::peptide_edit::template(&r.template)
        .map(|t| t.amino_acid.three_letter())
        .unwrap_or(&r.template);
    if let Some(o) = &r.origin {
        format!("{} · {}{} · {name}", o.chain, o.number, o.insertion)
    } else {
        format!("Chain {} · {} · {name}", r.chain, r.id)
    }
}
fn conformation_controls(
    ui: &mut egui::Ui,
    state: &mut State,
    mol: &Molecule,
    ready: bool,
) -> Option<Result<Build, String>> {
    use super::peptide_conformation::{self as conf, Kind, Preset};
    use super::peptide_rotamers as rot;
    let id = state.edit_residue?;
    let t = super::peptide_edit::metadata(mol).ok()?;
    if state.backbone.loaded != Some(id)
        || state.backbone.positions != mol.pos
        || state.backbone.topology != mol.topology
    {
        state.preview = None;
        state.rotamers = None;
        let mut controls = BackboneControls::default();
        controls.loaded = Some(id);
        controls.end = Some(id);
        controls.positions = mol.pos.clone();
        controls.topology = mol.topology.clone();
        for (i, k) in [Kind::Phi, Kind::Psi, Kind::Omega].into_iter().enumerate() {
            if let Some(v) = conf::measured(mol, id, k) {
                controls.angles[i] = v;
                controls.available[i] = true;
            }
        }
        for r in &t.residues {
            if let (Some(phi), Some(psi)) = (
                conf::measured(mol, r.id, Kind::Phi),
                conf::measured(mol, r.id, Kind::Psi),
            ) {
                controls.points.push((r.id, phi, psi));
            }
        }
        state.backbone = controls;
    }
    let mut result = None;
    egui::CollapsingHeader::new("Backbone · φ / ψ / ω").default_open(true).show(ui,|ui| {
        let controls=&mut state.backbone;let before=(controls.angles,controls.end,controls.preset);
        let mut preview=false;let mut preset=false;
        ui.add_enabled_ui(ready,|ui| {
            ui.horizontal(|ui| {
                ui.vertical(|ui| {
                    for (i,label) in ["φ","ψ","ω (incoming)"].into_iter().enumerate() {
                        ui.horizontal(|ui| {ui.label(label);if controls.available[i] {ui.add(egui::DragValue::new(&mut controls.angles[i]).range(-180.0..=180.0).suffix("°").speed(0.5));}else{ui.weak("Undefined at this end/break");}});
                    }
                    preview=ui.button("Preview backbone").clicked();
                    ui.horizontal(|ui| {
                        egui::ComboBox::from_id_salt("backbone_preset").selected_text(match controls.preset{Preset::Helix=>"α-helix",Preset::Strand=>"β-strand",Preset::TurnI=>"Type-I turn"}).show_ui(ui,|ui|{for(p,label)in[(Preset::Helix,"α-helix"),(Preset::Strand,"β-strand"),(Preset::TurnI,"Type-I turn")]{ui.selectable_value(&mut controls.preset,p,label);}});
                        ui.label("through");
                        let chain=t.residues.iter().find(|r|r.id==id).map(|r|r.chain);
                        let label=t.residues.iter().find(|r|Some(r.id)==controls.end).map(residue_label).unwrap_or_default();
                        egui::ComboBox::from_id_salt("backbone_range_end").selected_text(label).show_ui(ui,|ui|{for r in t.residues.iter().filter(|r|Some(r.chain)==chain && super::peptide_edit::template(&r.template).is_some()){ui.selectable_value(&mut controls.end,Some(r.id),residue_label(r));}});
                    });
                    preset=ui.button("Preview range preset").on_hover_text("A type-I turn sets the two central residues. Rings and cross-links are protected.").clicked();
                });
                ramachandran(ui,controls,id);
            });
        });
        if before!=(controls.angles,controls.end,controls.preset){state.preview=None;}
        if preview {let values:Vec<_>=controls.angles.iter().zip(controls.available).map(|(&v,ok)|ok.then_some(v)).collect();result=Some(conf::edit(mol,id,values[0],values[1],values[2]));}
        if preset {result=Some(conf::preset(mol,id,controls.end.unwrap_or(id),controls.preset));}
    });
    egui::CollapsingHeader::new("Side-chain rotamers").default_open(true).show(ui,|ui| {
        ui.add_enabled_ui(ready,|ui| {
            ui.horizontal_wrapped(|ui| {
                if ui.button("Find rotamers").clicked(){match rot::choices(mol,id){Ok(choices)=>{state.rotamers=Some(RotamerCache{choices,source:mol.pos.clone(),topology:mol.topology.clone(),index:0});state.error=None;},Err(e)=>{state.rotamers=None;result=Some(Err(e));}}}
                if ui.add(egui::DragValue::new(&mut state.neighbor_radius).range(2.0..=12.0).suffix(" Å").prefix("Nearby ")).changed() {state.preview=None;}
                if ui.button("Preview nearby repack").on_hover_text("Keep this residue and all backbones fixed; choose lower-clash rotamers for nearby side chains.").clicked(){result=Some(rot::repack_neighbors(mol,id,state.neighbor_radius));}
            });
            if state.rotamers.as_ref().is_some_and(|c|c.source!=mol.pos||c.topology!=mol.topology||c.choices.residue!=id){state.rotamers=None;}
            if let Some(cache)=&mut state.rotamers {
                if cache.choices.marginal {ui.small("Terminal rotamers: averaged over undefined backbone angles.");}
                ui.horizontal_wrapped(|ui|{
                    let mut selected=cache.index;
                    if ui.add_enabled(selected>0,egui::Button::new("←")).on_hover_text("Previous rotamer").clicked(){selected-=1;}
                    egui::ComboBox::from_id_salt("rotamer_choice").selected_text(if selected==0{"Current".into()}else{format!("Rotamer {selected}")}).show_ui(ui,|ui| {
                        for(i,c)in cache.choices.candidates.iter().enumerate(){let label=if c.current {format!("Current · overlap {:.3}",c.score)}else{format!("{i} · overlap {:.3} · prior {:.1}%",c.score,c.probability.unwrap_or(0.0)*100.0)};ui.selectable_value(&mut selected,i,label);}
                    });
                    if ui.add_enabled(selected+1<cache.choices.candidates.len(),egui::Button::new("→")).on_hover_text("Next rotamer").clicked(){selected+=1;}
                    if selected!=cache.index {cache.index=selected;result=Some(rot::preview(mol,&cache.choices,selected));}
                    if ui.button("Preview rotamer").clicked(){result=Some(rot::preview(mol,&cache.choices,cache.index));}
                    let c=&cache.choices.candidates[cache.index];ui.small(c.chi.iter().enumerate().map(|(i,v)|format!("χ{} {v:.0}°",i+1)).collect::<Vec<_>>().join(" · "));
                });
                ui.small("Ranked by local overlap, then library probability. Overlap is a steric score, not energy.");
            }
        });
    });
    result
}
fn ramachandran(ui: &mut egui::Ui, controls: &mut BackboneControls, id: u64) {
    let (rect, response) =
        ui.allocate_exact_size(egui::vec2(156.0, 156.0), egui::Sense::click_and_drag());
    let area = rect.shrink(13.0);
    let p = ui.painter();
    p.rect_filled(rect, 6, ui.visuals().extreme_bg_color);
    let project = |phi: f32, psi: f32| {
        egui::pos2(
            area.left() + (phi + 180.0) / 360.0 * area.width(),
            area.bottom() - (psi + 180.0) / 360.0 * area.height(),
        )
    };
    p.line_segment(
        [project(0.0, -180.0), project(0.0, 180.0)],
        egui::Stroke::new(1.0, ui.visuals().weak_text_color()),
    );
    p.line_segment(
        [project(-180.0, 0.0), project(180.0, 0.0)],
        egui::Stroke::new(1.0, ui.visuals().weak_text_color()),
    );
    for &(r, phi, psi) in &controls.points {
        p.circle_filled(
            project(phi, psi),
            if r == id { 3.5 } else { 2.0 },
            if r == id {
                ACCENT
            } else {
                ui.visuals().weak_text_color()
            },
        );
    }
    if controls.available[0] && controls.available[1] {
        let at = project(controls.angles[0], controls.angles[1]);
        p.circle_stroke(at, 5.0, egui::Stroke::new(1.5, egui::Color32::YELLOW));
    }
    p.text(
        egui::pos2(rect.right() - 5.0, area.center().y),
        egui::Align2::RIGHT_BOTTOM,
        "φ",
        egui::FontId::proportional(11.0),
        ui.visuals().text_color(),
    );
    p.text(
        egui::pos2(area.center().x, rect.top() + 2.0),
        egui::Align2::LEFT_TOP,
        "ψ",
        egui::FontId::proportional(11.0),
        ui.visuals().text_color(),
    );
    for (phi, psi, label) in [(-57.0, -47.0, "α"), (-135.0, 135.0, "β")] {
        p.text(
            project(phi, psi),
            egui::Align2::LEFT_BOTTOM,
            label,
            egui::FontId::proportional(11.0),
            ui.visuals().weak_text_color(),
        );
    }
    if controls.available[0] && controls.available[1] && (response.clicked() || response.dragged())
    {
        if let Some(at) = response.interact_pointer_pos() {
            controls.angles[0] =
                ((at.x - area.left()) / area.width() * 360.0 - 180.0).clamp(-180.0, 180.0);
            controls.angles[1] =
                ((area.bottom() - at.y) / area.height() * 360.0 - 180.0).clamp(-180.0, 180.0);
        }
    }
    response.on_hover_text("Ramachandran plot: φ and ψ from −180° to +180°. Click to set the preview target. α/β are reference positions, not probability contours or a validation map.");
}

#[cfg(test)]
mod protein_panel_tests {
    use super::*;
    #[test]
    fn protein_edit_panel_renders_backbone_and_rotamers_without_changing_structure() {
        let ctx = egui::Context::default();
        let mut mol = peptide::build(
            &Recipe {
                residues: peptide::parse_sequence("AKLA").unwrap(),
                ..Default::default()
            },
            None,
        )
        .unwrap()
        .molecule;
        let original = mol.clone();
        let mut zmat = ZMatrixBuilderState::default();
        zmat.peptide.open = true;
        zmat.peptide.mode = Mode::Edit;
        zmat.peptide.edit_residue = Some(2);
        let mut editor = EditorRotateState::default();
        let mut settings = MolSettings::default();
        let mut traj = TrajectoryState::default();
        let mut rect = None;
        for frame in 0..4 {
            let mut output = ctx.run_ui(
                egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(1280.0, 720.0),
                    )),
                    ..Default::default()
                },
                |ui| {
                    let (r, change) = window(
                        ui.ctx(),
                        &mut zmat,
                        &mut editor,
                        &mut mol,
                        &mut settings,
                        &mut traj,
                    );
                    rect = r;
                    assert!(!change.happened());
                },
            );
            output.textures_delta.clear();
            if frame == 1 {
                zmat.peptide.rotamers = Some(RotamerCache {
                    choices: super::super::peptide_rotamers::choices(&mol, 2).unwrap(),
                    source: mol.pos.clone(),
                    topology: mol.topology.clone(),
                    index: 1,
                });
            }
        }
        assert_eq!(zmat.peptide.backbone.available, [true; 3]);
        assert!(zmat.peptide.rotamers.is_some());
        assert_eq!(mol.pos, original.pos);
        assert_eq!(mol.topology, original.topology);
        assert!(rect.unwrap().bottom() <= 720.0);
        assert!(rect.unwrap().right() <= 1280.0);
    }
}
