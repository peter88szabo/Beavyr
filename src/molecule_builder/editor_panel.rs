//! The molecule/cluster editor, using the existing placement and fragment tools.
//! Layout preferences live in egui, never in a molecular undo snapshot.

use super::*;
use crate::molecule::{atomic_number, composition, element_symbol, is_dummy};
use crate::structure_history::StructureHistory;
use crate::trajectory::TrajectoryState;

const ACCENT: egui::Color32 = egui::Color32::from_rgb(104, 174, 255);
const COMMON_ATOMS: &[&str] = &["H", "C", "N", "O", "S", "F", "Cl"];
const COMMON_FRAGMENTS: &[&str] = &["-CH3", "-OH", "-NH2", "-Phenyl", "-OCH3", "-COOH"];

#[derive(Clone, Default)]
struct PanelState {
    docked: bool,
    close_requested: bool,
    search: String,
    all_fragments: bool,
    my_fragments: bool,
    save_fragment: super::super::custom_fragments::SaveForm,
    status: String,
}

fn panel_state(ctx: &egui::Context) -> PanelState {
    ctx.data_mut(|data| {
        data.get_temp_mut_or_default::<PanelState>(egui::Id::new("molecular_editor_layout"))
            .clone()
    })
}

fn store_panel_state(ctx: &egui::Context, state: PanelState) {
    ctx.data_mut(|data| data.insert_temp(egui::Id::new("molecular_editor_layout"), state));
}

pub(crate) fn restore_editor_docking(ctx: &egui::Context, docked: bool) {
    let mut state = panel_state(ctx);
    state.docked = docked;
    store_panel_state(ctx, state);
}

pub(crate) fn editor_docked(ctx: &egui::Context) -> bool {
    panel_state(ctx).docked
}

pub(crate) fn take_close_request(ctx: &egui::Context) -> bool {
    let mut state = panel_state(ctx);
    let close = std::mem::take(&mut state.close_requested);
    store_panel_state(ctx, state);
    close
}

fn style_panel(ui: &mut egui::Ui) {
    let style = ui.style_mut();
    style.spacing.item_spacing = egui::vec2(8.0, 8.0);
    style.spacing.button_padding = egui::vec2(10.0, 6.0);
    style.spacing.interact_size.y = 30.0;
    style.visuals.selection.bg_fill = egui::Color32::from_rgb(35, 75, 117);
    style.visuals.selection.stroke = egui::Stroke::new(1.0, ACCENT);
    for widget in [
        &mut style.visuals.widgets.inactive,
        &mut style.visuals.widgets.hovered,
        &mut style.visuals.widgets.active,
    ] {
        widget.corner_radius = egui::CornerRadius::same(6);
    }
}

fn card(ui: &egui::Ui) -> egui::Frame {
    egui::Frame::new()
        .fill(ui.visuals().faint_bg_color)
        .stroke(ui.visuals().widgets.noninteractive.bg_stroke)
        .corner_radius(8)
        .inner_margin(10)
}

#[derive(Clone)]
enum Action {
    Peptide(Box<super::super::peptide::Build>),
    Place,
    Add,
    Replace,
    Delete(usize),
    ChangeElement(usize, String),
    Hydrogens,
    Cleanup,
    Undo,
    Redo,
    References,
}

fn pending_fragment_edit(
    state: &EditorRotateState,
    zmat: &ZMatrixBuilderState,
    mol: &Molecule,
) -> bool {
    [state, &zmat.placed_editor].iter().any(|editor| {
        editor
            .base_pos
            .as_ref()
            .is_some_and(|base| base != &mol.pos)
    })
}

/// Only the main editor receives these shortcuts. A text field, popup or the
/// independent Fragment Editor keeps ownership of its keyboard interaction.
fn shortcut_action(ui: &egui::Ui, available: bool, zmat: &ZMatrixBuilderState) -> Option<Action> {
    if !available
        || ui.ctx().egui_wants_keyboard_input()
        || ui.ctx().any_popup_open()
        || !ui.rect_contains_pointer(ui.max_rect())
    {
        return None;
    }
    ui.input_mut(|input| {
        if input.consume_key(
            egui::Modifiers::COMMAND.plus(egui::Modifiers::SHIFT),
            egui::Key::Z,
        ) {
            Some(Action::Redo)
        } else if input.consume_key(egui::Modifiers::COMMAND, egui::Key::Z) {
            Some(Action::Undo)
        } else if input.consume_key(egui::Modifiers::NONE, egui::Key::Delete) {
            zmat.selected_index.map(Action::Delete)
        } else {
            None
        }
    })
}

pub(super) fn contents(
    ui: &mut egui::Ui,
    editor: &mut EditorRotateState,
    zmat: &mut ZMatrixBuilderState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
    history: &mut StructureHistory,
    traj: &mut TrajectoryState,
) -> BuilderChange {
    let mut panel = panel_state(ui.ctx());
    if zmat.selection.sync(mol) {
        settings.geometry_dirty = true;
    }
    let original_style = ui.style().clone();
    let pending = pending_fragment_edit(editor, zmat, mol);
    let available = !pending && !traj.playing && !zmat.confirm_clear;
    let mut action = shortcut_action(
        ui,
        available && !editor.active && !zmat.add_atom_active,
        zmat,
    );

    ui.scope(|ui| {
        style_panel(ui);
        header(ui, &mut panel, zmat, mol, history, available, &mut action);
        selection(ui, zmat, mol, settings, available, &mut action);

        // Reserve the footer before laying out the scroll area: its controls
        // never disappear behind a long library or the placed-fragment form.
        let remaining = ui.available_rect_before_wrap();
        let footer_height = 116.0;
        let body_bottom = (remaining.bottom() - footer_height - 8.0).max(remaining.top() + 80.0);
        let body_rect =
            egui::Rect::from_min_max(remaining.min, egui::pos2(remaining.right(), body_bottom));
        let footer_rect = egui::Rect::from_min_max(
            egui::pos2(remaining.left(), body_bottom + 8.0),
            egui::pos2(
                remaining.right(),
                (body_bottom + 8.0 + footer_height).max(remaining.bottom()),
            ),
        );

        let mut body = ui.new_child(
            egui::UiBuilder::new()
                .id_salt("modern_editor_body")
                .max_rect(body_rect),
        );
        egui::ScrollArea::vertical()
            .id_salt("modern_editor_scroll")
            .auto_shrink([false, false])
            .max_height(body_rect.height())
            .show(&mut body, |ui| {
                if pending {
                    ui.label("Apply or cancel the fragment adjustment before editing atoms.");
                } else if traj.playing {
                    ui.label("Pause trajectory playback to edit this structure.");
                }
                ui.add_enabled_ui(available && !zmat.add_atom_active, |ui| {
                    ui.horizontal_wrapped(|ui| {
                        ui.selectable_value(&mut panel.my_fragments, false, "Built-in");
                        ui.selectable_value(&mut panel.my_fragments, true, "My fragments");
                        if ui.button("Peptide Builder…").clicked() {
                            zmat.peptide.open = true;
                        }
                    });
                    if panel.my_fragments {
                        ui.add(
                            egui::TextEdit::singleline(&mut panel.search)
                                .hint_text("Search my fragments")
                                .desired_width(f32::INFINITY),
                        );
                        let moving = zmat
                            .placed_editor
                            .moving
                            .as_deref()
                            .or(editor.moving.as_deref());
                        super::super::custom_fragments::show(
                            ui,
                            &mut zmat.custom_fragments,
                            &mut panel.save_fragment,
                            &mut zmat.frag_name,
                            mol,
                            zmat.selected_index,
                            moving,
                            &panel.search,
                        );
                    } else {
                        palette(ui, &mut panel, zmat);
                    }
                    placement_controls(ui, zmat, mol, &mut action);
                });
                reference_controls(ui, zmat, editor, mol, available, &mut action);

                if let Some(error) = &zmat.last_error {
                    ui.colored_label(egui::Color32::LIGHT_RED, error);
                }
                if let Some(report) = &zmat.cleanup_report {
                    match report {
                        Ok(message) => {
                            ui.small(message);
                        }
                        Err(message) => {
                            ui.colored_label(egui::Color32::LIGHT_RED, message);
                        }
                    }
                }
                crate::molecule_builder::hydrogens::report(ui, &zmat.hydrogen_state);

                // This is the existing form, with its original style restored.
                // Its sliders, state, Apply/Cancel/Reset and geometry code are
                // deliberately shared, not copied or restyled.
                ui.scope(|ui| {
                    ui.set_style(original_style.clone());
                    placed_fragment_editor(ui, zmat, mol, settings);
                });
            });

        let mut footer_ui = ui.new_child(
            egui::UiBuilder::new()
                .id_salt("modern_editor_footer")
                .max_rect(footer_rect),
        );
        footer(
            &mut footer_ui,
            &panel,
            editor,
            zmat,
            mol,
            available,
            &mut action,
        );
        ui.allocate_rect(remaining, egui::Sense::hover());
    });

    if zmat.cleanup.show(ui.ctx(), mol, &zmat.selection, available) {
        action = Some(Action::Cleanup);
    }
    let mut change = BuilderChange::None;
    if let Some(action) = action {
        match perform(action, editor, zmat, mol, settings, traj) {
            Ok(Some((result, message))) => {
                change = result;
                panel.status = message;
            }
            Ok(None) => {}
            Err(error) => zmat.last_error = Some(error),
        }
    }
    store_panel_state(ui.ctx(), panel);
    change
}

fn header(
    ui: &mut egui::Ui,
    panel: &mut PanelState,
    zmat: &mut ZMatrixBuilderState,
    mol: &Molecule,
    history: &mut StructureHistory,
    available: bool,
    action: &mut Option<Action>,
) {
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new("Molecular Editor").size(21.0).strong());
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if ui.small_button("×").on_hover_text("Close editor").clicked() {
                panel.close_requested = true;
            }
            if ui
                .small_button(if panel.docked { "Float" } else { "Dock" })
                .clicked()
            {
                panel.docked = !panel.docked;
            }
            ui.menu_button("•••", |ui| {
                if ui
                    .add_enabled(
                        available && !mol.atoms.is_empty(),
                        egui::Button::new("Clear Display…"),
                    )
                    .clicked()
                {
                    zmat.confirm_clear = true;
                    ui.close();
                }
            });
        });
    });
    let comp = composition(&mol.atoms);
    let summary = if comp.natoms == 0 {
        "Start with an atom or a fragment".to_owned()
    } else {
        format!(
            "{}  ·  {} atoms",
            comp.formula.replace(" : ", ""),
            comp.natoms
        )
    };
    let response = ui.label(egui::RichText::new(summary).color(ACCENT));
    if let Some(mass) = comp.mass_amu {
        response.on_hover_text(format!("Molecular weight: {mass:.3} g/mol"));
    } else if !comp.unknown.is_empty() {
        response.on_hover_text(format!("Unknown elements: {}", comp.unknown.join(", ")));
    }
    ui.horizontal_wrapped(|ui| {
        if ui
            .add_enabled(
                available && zmat.history.next_undo().is_some(),
                egui::Button::new("Undo"),
            )
            .on_hover_text(
                zmat.history
                    .next_undo()
                    .map(|s| format!("Undo: {s}\nCtrl/Cmd+Z over this panel"))
                    .unwrap_or_else(|| "Nothing to undo".into()),
            )
            .clicked()
        {
            *action = Some(Action::Undo);
        }
        if ui
            .add_enabled(
                available && zmat.history.next_redo().is_some(),
                egui::Button::new("Redo"),
            )
            .on_hover_text(
                zmat.history
                    .next_redo()
                    .map(|s| format!("Redo: {s}"))
                    .unwrap_or_else(|| "Nothing to redo".into()),
            )
            .clicked()
        {
            *action = Some(Action::Redo);
        }
        ui.add_enabled_ui(available && !mol.atoms.is_empty(), |ui| {
            ui.menu_button("Save to history", |ui| {
                ui.label("Keep this version across sessions");
                ui.add(
                    egui::TextEdit::singleline(&mut history.save_name)
                        .hint_text("Name (optional)")
                        .desired_width(210.0),
                );
                if ui.button("Save version").clicked() {
                    history.save_current(mol);
                    panel.status = "Saved to history".into();
                    ui.close();
                }
            });
        });
    });
    if history.is_error {
        if let Some(message) = &history.message {
            ui.colored_label(egui::Color32::LIGHT_RED, message);
        }
    }
}

fn selection(
    ui: &mut egui::Ui,
    zmat: &mut ZMatrixBuilderState,
    mol: &Molecule,
    settings: &mut MolSettings,
    available: bool,
    action: &mut Option<Action>,
) {
    card(ui).show(ui, |ui| {
        ui.set_min_width((ui.available_width() - 1.0).max(0.0));
        if let Some(index) = zmat.selected_index.filter(|&i| i < mol.atoms.len()) {
            ui.horizontal(|ui| {
                ui.strong(format!("Active atom  {}{}", mol.atoms[index], index + 1));
                ui.menu_button("Atom…", |ui| {
                    if ui.button("Deselect").clicked() {
                        set_selected_atom(zmat, &mol.atoms, None);
                        zmat.selection.atoms.clear();
                        ui.close();
                    }
                    ui.add_enabled_ui(available, |ui| {
                        ui.menu_button("Change element…", |ui| {
                            if let Some(symbol) = element_picker(ui, "change_selected_element") {
                                *action = Some(Action::ChangeElement(index, symbol));
                                ui.close();
                            }
                        });
                        if ui.button("Delete atom").clicked() {
                            *action = Some(Action::Delete(index));
                            ui.close();
                        }
                    });
                });
            });
        } else if mol.atoms.is_empty() {
            ui.strong("Build a molecule or cluster");
            ui.small("Choose below, then Place to start.");
        } else {
            ui.strong("Select atoms in the 3D view").on_hover_text(
                "Click an atom to attach a fragment. Shift-click toggles multiple selection.",
            );
        }
        if !mol.atoms.is_empty() {
            zmat.selection
                .controls(ui, mol, settings, zmat.placed_editor.moving.as_deref());
        }
    });
}

fn choose_atom(zmat: &mut ZMatrixBuilderState, symbol: &str) {
    zmat.new_symbol = symbol.to_owned();
    zmat.frag_name = format!("Atom: {symbol}");
    zmat.last_error = None;
}

fn element_picker(ui: &mut egui::Ui, id: &str) -> Option<String> {
    let mut result = None;
    ui.push_id(id, |ui| {
        ui.label("Choose an element");
        egui::ScrollArea::vertical()
            .max_height(300.0)
            .show(ui, |ui| {
                egui::Grid::new("elements")
                    .spacing(egui::vec2(4.0, 4.0))
                    .show(ui, |ui| {
                        for (index, symbol) in (1..=118)
                            .filter_map(element_symbol)
                            .chain(std::iter::once("X"))
                            .enumerate()
                        {
                            if ui
                                .add_sized([35.0, 30.0], egui::Button::new(symbol))
                                .on_hover_text(if symbol == "X" {
                                    "Dummy / attachment centre"
                                } else {
                                    symbol
                                })
                                .clicked()
                            {
                                result = Some(symbol.to_owned());
                            }
                            if index % 8 == 7 {
                                ui.end_row();
                            }
                        }
                    });
            });
    });
    result
}

fn fragment_description(name: &str) -> &'static str {
    match name {
        "-CH3" => "Methyl",
        "-OH" => "Hydroxyl",
        "-NH2" => "Amino",
        "-Phenyl" => "Phenyl / aromatic ring",
        "-Pyrrole" => "Pyrrole ring",
        "-OCH3" => "Methoxy",
        "-COOH" => "Carboxyl",
        "-NO2" => "Nitro",
        "-OOH" => "Hydroperoxy",
        "-CCH" => "Ethynyl / alkyne",
        "-CH=CH2" => "Vinyl / alkene",
        "-CH=O" => "Formyl / aldehyde",
        "-CycloPentane" => "Cyclopentane ring",
        "-CycloHexane" => "Cyclohexane ring",
        "-SH" => "Thiol",
        "5-Ring centered" => "Cyclopentadienyl with dummy centre",
        "6-Ring centered" => "Benzene with dummy centre",
        _ => "Fragment",
    }
}

fn fragment_matches(name: &str, query: &str) -> bool {
    let query = query.trim().to_lowercase();
    name.to_lowercase().contains(&query)
        || fragment_description(name).to_lowercase().contains(&query)
}

fn palette(ui: &mut egui::Ui, panel: &mut PanelState, zmat: &mut ZMatrixBuilderState) {
    ui.horizontal_wrapped(|ui| {
        for symbol in COMMON_ATOMS {
            if ui
                .add_sized(
                    [33.0, 32.0],
                    egui::Button::new(*symbol)
                        .selected(zmat.frag_name == format!("Atom: {symbol}")),
                )
                .clicked()
            {
                choose_atom(zmat, symbol);
            }
        }
        ui.menu_button("More…", |ui| {
            if let Some(symbol) = element_picker(ui, "insert_element") {
                choose_atom(zmat, &symbol);
                ui.close();
            }
        });
    });
    ui.add(
        egui::TextEdit::singleline(&mut panel.search)
            .hint_text("Search fragments — name or formula")
            .desired_width(f32::INFINITY),
    );

    let groups: Vec<_> = tile_fragments()
        .into_iter()
        .filter(|fragment| {
            if !panel.search.trim().is_empty() {
                fragment_matches(fragment.name, &panel.search)
            } else {
                panel.all_fragments || COMMON_FRAGMENTS.contains(&fragment.name)
            }
        })
        .collect();
    // Three equal columns stay legible at the panel's smallest supported width.
    let tile_width = ((ui.available_width() - 16.0) / 3.0).max(70.0);
    egui::Grid::new("modern_fragment_palette")
        .spacing(egui::vec2(8.0, 6.0))
        .show(ui, |ui| {
            for (index, fragment) in groups.iter().enumerate() {
                let name = fragment.name.trim_start_matches('-');
                let label = match name {
                    "CycloPentane" => "5-ring",
                    "CycloHexane" => "6-ring",
                    "5-Ring centered" => "5-ring + centre",
                    "6-Ring centered" => "6-ring + centre",
                    _ => name,
                };
                if ui
                    .add_sized(
                        [tile_width, 36.0],
                        egui::Button::new(label).selected(zmat.frag_name == fragment.name),
                    )
                    .on_hover_text(fragment_description(fragment.name))
                    .clicked()
                {
                    zmat.frag_name = fragment.name.to_owned();
                    zmat.last_error = None;
                }
                if index % 3 == 2 {
                    ui.end_row();
                }
            }
        });
    if groups.is_empty() {
        ui.small("No matching fragments. Try a formula such as CH3 or a name such as methyl.");
    }
    if panel.search.is_empty() {
        ui.checkbox(&mut panel.all_fragments, "Show all fragments");
    } else if ui.small_button("Clear search").clicked() {
        panel.search.clear();
    }
}

fn placement_controls(
    ui: &mut egui::Ui,
    zmat: &mut ZMatrixBuilderState,
    mol: &Molecule,
    action: &mut Option<Action>,
) {
    card(ui).show(ui, |ui| {
        ui.set_min_width((ui.available_width() - 1.0).max(0.0));
        ui.horizontal_wrapped(|ui| {
            ui.strong(format!("Chosen: {}", zmat.frag_name));
            ui.label("Bond");
            egui::ComboBox::from_id_salt("modern_bond_order")
                .width(75.0)
                .selected_text(zmat.frag_bond_order.label())
                .show_ui(ui, |ui| {
                    for order in [BondOrder::Single, BondOrder::Double, BondOrder::Triple] {
                        if ui
                            .selectable_value(&mut zmat.frag_bond_order, order, order.label())
                            .changed()
                        {
                            zmat.frag_angle_deg = order.default_angle();
                            zmat.frag_dihedral_deg = order.default_dihedral();
                        }
                    }
                });
        });
        let selection = zmat.selected_index.filter(|&i| i < mol.atoms.len());
        let chosen = find_selected_fragment(zmat).is_some();
        let primary = if mol.atoms.is_empty() {
            format!("Place {}", zmat.frag_name)
        } else if let Some(i) = selection {
            format!("Add to {}{}", mol.atoms[i], i + 1)
        } else {
            "Add to selection".to_owned()
        };
        ui.horizontal_wrapped(|ui| {
            let button =
                egui::Button::new(egui::RichText::new(primary).color(egui::Color32::WHITE))
                    .fill(egui::Color32::from_rgb(38, 96, 159))
                    .min_size(egui::vec2(150.0, 34.0));
            if ui
                .add_enabled(
                    chosen && (mol.atoms.is_empty() || selection.is_some()),
                    button,
                )
                .on_disabled_hover_text("Select an atom in the 3D view first.")
                .clicked()
            {
                *action = Some(if mol.atoms.is_empty() {
                    Action::Place
                } else {
                    Action::Add
                });
            }
            if !mol.atoms.is_empty() {
                let label = selection
                    .map(|i| format!("Replace {}{}", mol.atoms[i], i + 1))
                    .unwrap_or_else(|| "Replace selection".into());
                if ui
                    .add_enabled(chosen && selection.is_some(), egui::Button::new(label))
                    .on_disabled_hover_text("Select the atom you want to replace.")
                    .clicked()
                {
                    *action = Some(Action::Replace);
                }
            }
        });
    });
}

fn reference_controls(
    ui: &mut egui::Ui,
    zmat: &mut ZMatrixBuilderState,
    editor: &EditorRotateState,
    mol: &Molecule,
    available: bool,
    action: &mut Option<Action>,
) {
    ui.add_enabled_ui(available, |ui| {
        egui::CollapsingHeader::new("Placement settings").id_salt("modern_placement_settings").show(ui, |ui| {
            ui.small("Attachment direction is chosen automatically. Set orientation or pick exact atom references here.");
            ui.horizontal_wrapped(|ui| {
                ui.label("Angle");
                ui.add(egui::DragValue::new(&mut zmat.frag_angle_deg).speed(0.1).suffix("°"));
                ui.label("Dihedral");
                ui.add(egui::DragValue::new(&mut zmat.frag_dihedral_deg).speed(0.1).suffix("°"));
            });
            if zmat.frag_name.starts_with("Atom: ") && !mol.atoms.is_empty() && !zmat.add_atom_active {
                if ui.button("Pick references…").clicked() {
                    zmat.new_angle_deg = zmat.frag_angle_deg;
                    zmat.new_dihedral_deg = zmat.frag_dihedral_deg;
                    zmat.new_bond_order = zmat.frag_bond_order;
                    zmat.add_atom_active = true;
                    zmat.add_atom_auto = false;
                    zmat.add_atom_picks.clear();
                    zmat.add_atom_status = None;
                }
            }
        });
    });
    if zmat.add_atom_active {
        if let Some(prompt) = next_click_meaning(zmat, editor) {
            ui.strong(prompt);
        }
        if let Some(message) = &zmat.add_atom_status {
            ui.label(message);
        }
        if ui.button("Cancel picking").clicked() {
            zmat.add_atom_active = false;
            zmat.add_atom_picks.clear();
            zmat.add_atom_status = None;
        }
        let needed = if zmat.add_atom_auto {
            1
        } else {
            add_atom_picks_required(mol.atoms.len())
        };
        if available && zmat.add_atom_active && needed > 0 && zmat.add_atom_picks.len() >= needed {
            *action = Some(Action::References);
        }
    }
}

fn footer(
    ui: &mut egui::Ui,
    panel: &PanelState,
    editor: &mut EditorRotateState,
    zmat: &mut ZMatrixBuilderState,
    mol: &Molecule,
    available: bool,
    action: &mut Option<Action>,
) {
    ui.separator();
    ui.horizontal(|ui| {
        if ui
            .add_enabled(
                available && !mol.atoms.is_empty(),
                egui::Button::new("Add missing H"),
            )
            .on_hover_text(
                "Restore peptide H atoms and fill supported carbon, nitrogen and oxygen sites. Existing atoms stay fixed.",
            )
            .clicked()
        {
            *action = Some(Action::Hydrogens);
        }
        if ui
            .add_enabled(
                available && !mol.atoms.is_empty(),
                egui::Button::new("Clean up geometry"),
            )
            .on_hover_text(
                "Built-in DREIDING cleanup. Produces a starting geometry for optimisation.",
            )
            .clicked()
        {
            zmat.cleanup.open = true;
        }
    });
    ui.horizontal(|ui| {
        if ui.button("Z-matrix ↗").clicked() {
            zmat.zmat_window_open = !zmat.zmat_window_open;
        }
        if ui.button("Fragment Editor ↗").clicked() {
            if editor.window_open {
                editor.window_open = false;
            } else {
                // Hand the known target to the existing window. Only one
                // editor owns its live preview, so the inline form steps aside.
                if zmat.placed_editor.moving.is_some() {
                    *editor = std::mem::take(&mut zmat.placed_editor);
                }
                editor.window_open = true;
            }
        }
    });
    let status = if zmat.add_atom_active {
        "Picking atom references — see Placement settings"
    } else if editor.active {
        "Fragment Editor is picking atoms in the 3D view"
    } else {
        &panel.status
    };
    ui.add(egui::Label::new(egui::RichText::new(status).small().weak()).truncate())
        .on_hover_text(status);
}

/// Existing geometry operations run within a single transaction boundary,
/// including deletion and hydrogen completion. Failed
/// operations do not add undo entries or destroy the last good geometry.
fn perform(
    action: Action,
    editor: &mut EditorRotateState,
    zmat: &mut ZMatrixBuilderState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
    traj: &mut TrajectoryState,
) -> Result<Option<(BuilderChange, String)>, String> {
    if pending_fragment_edit(editor, zmat, mol) {
        return Err("Apply or cancel the fragment adjustment first.".into());
    }
    if traj.playing {
        return Err("Pause trajectory playback before editing.".into());
    }
    if matches!(action, Action::Undo | Action::Redo) {
        let now = current_step(zmat, mol);
        let undo = matches!(action, Action::Undo);
        let name = if undo {
            zmat.history.next_undo()
        } else {
            zmat.history.next_redo()
        }
        .map(str::to_owned);
        let step = if undo {
            zmat.history.undo(now)
        } else {
            zmat.history.redo(now)
        };
        if let Some(step) = step {
            restore_step(step, zmat, mol, settings);
            zmat.selection = Default::default();
            zmat.selection.sync(mol);
            zmat.cleanup.cancel();
            reset_fragment_targets(editor, zmat);
            zmat.hydrogen_state = Default::default();
            zmat.cleanup_report = None;
            finish_edit(zmat, traj);
            return Ok(Some((
                BuilderChange::Edited,
                format!(
                    "{}: {}",
                    if undo { "Undid" } else { "Redid" },
                    name.unwrap_or_default()
                ),
            )));
        }
        return Ok(None);
    }
    let before = current_step(zmat, mol);
    let empty_before = mol.atoms.is_empty();
    let outcome = perform_mutation(&action, editor, zmat, mol, settings);
    let message = match outcome {
        Ok(message) => message,
        Err(error) => {
            // Some older operations reconstruct before reporting an error.
            // Preserve the exact displayed Cartesian snapshot in that case.
            if mol.atoms != before.atoms || mol.pos != before.pos || mol.topology != before.topology {
                restore_step(before, zmat, mol, settings);
                reset_fragment_targets(editor, zmat);
            }
            return Err(error);
        }
    };
    if before.atoms == mol.atoms && before.pos == mol.pos && before.topology == mol.topology {
        return Ok(None);
    }
    zmat.history
        .record(&message, &before.atoms, &before.pos, &before.zmat);
    zmat.history.set_last_topology(before.topology.clone());
    if !matches!(action, Action::Cleanup) {
        zmat.cleanup_report = None;
    }
    if !matches!(action, Action::Add | Action::Replace | Action::Peptide(_)) {
        reset_fragment_targets(editor, zmat);
    }
    if !matches!(action, Action::Peptide(_))
        && (!mol.atoms.starts_with(&before.atoms) || matches!(action, Action::Replace)) {
        zmat.selection = Default::default();
        zmat.selection.sync(mol);
    }
    settings.geometry_dirty = true;
    settings.bond_topology_dirty = true;
    zmat.hydrogen_state.invalidate_if_changed(mol);
    zmat.edit_refresh = true;
    zmat.last_error = None;
    finish_edit(zmat, traj);
    Ok(Some((
        if empty_before {
            BuilderChange::Replaced
        } else {
            BuilderChange::Edited
        },
        message,
    )))
}

fn finish_edit(zmat: &mut ZMatrixBuilderState, traj: &mut TrajectoryState) {
    zmat.skip_next_sync = true;
    zmat.add_atom_active = false;
    zmat.add_atom_auto = false;
    zmat.add_atom_picks.clear();
    zmat.add_atom_status = None;
    zmat.edit_preview.clear();
    traj.clear_for_structure();
}

fn reset_fragment_targets(editor: &mut EditorRotateState, zmat: &mut ZMatrixBuilderState) {
    let open = editor.window_open;
    *editor = EditorRotateState::default();
    editor.window_open = open;
    zmat.placed_editor = EditorRotateState::default();
    zmat.placed_atom_count = 0;
}

fn perform_mutation(
    action: &Action,
    editor: &mut EditorRotateState,
    zmat: &mut ZMatrixBuilderState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
) -> Result<String, String> {
    match action {
        Action::Peptide(build) => {
            if zmat.selection.frozen.iter().any(|&i| {
                build.host_mapping.get(i).copied().flatten()
                    .is_none_or(|j| build.molecule.pos.get(j) != mol.pos.get(i))
            }) {
                return Err("This edit would move or remove frozen atoms. Unfreeze them first.".into());
            }
            *mol = build.molecule.clone();
            if let Some(t) = &mol.topology { t.display_defaults(settings); }
            mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
            zmat.zmat = zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos);
            let frozen = zmat.selection.frozen.iter().filter_map(|&i|
                build.host_mapping.get(i).copied().flatten()).collect();
            let hidden = zmat.selection.hidden.iter().filter_map(|&i|
                build.host_mapping.get(i).copied().flatten()).collect();
            zmat.selection = Default::default();
            zmat.selection.sync(mol);
            zmat.selection.frozen = frozen;
            zmat.selection.hidden = hidden;
            zmat.cleanup.cancel();
            set_selected_atom(zmat, &mol.atoms, None);
            let open = editor.window_open;
            *editor = EditorRotateState::default();
            editor.window_open = open;
            zmat.placed_editor = EditorRotateState {
                picks: build.picks.to_vec(),
                chosen_side: Some(RotateSide::B),
                moving: Some(build.moving.clone()),
                session_origin: Some(mol.pos.clone()),
                ..Default::default()
            };
            zmat.placed_atom_count = mol.atoms.len();
            Ok(format!("Peptide: {}", build.sequence))
        }
        Action::Place => {
            if !mol.atoms.is_empty() {
                return Err("Select an atom to attach to.".into());
            }
            let fragment =
                find_selected_fragment(zmat).ok_or("Choose an atom or fragment first.")?;
            zmat.zmat.clear();
            add_fragment_to_zmat(
                &mut zmat.zmat,
                &fragment.xyz,
                None,
                0.0,
                zmat.frag_angle_deg,
                zmat.frag_dihedral_deg,
            );
            mol.atoms = zmat.zmat.iter().map(|atom| atom.symbol.clone()).collect();
            mol.pos = zmat2xyz::zmat_to_xyz(&zmat.zmat);
            mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
            Ok(format!("Placed {}", zmat.frag_name))
        }
        Action::Add | Action::Replace => {
            let i = zmat
                .selected_index
                .filter(|&i| i < mol.atoms.len())
                .ok_or(NOTHING_SELECTED)?;
            let target = format!("{}{}", mol.atoms[i], i + 1);
            // The independent Fragment Editor may have moved the Cartesian
            // structure since its Z-matrix was last displayed.
            zmat.zmat = zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos);
            zmat.last_error = None;
            if matches!(action, Action::Add) {
                zmat.frag_mode = FragmentInsertMode::Connect;
                commit_fragment_connect(zmat, mol, settings);
            } else {
                zmat.frag_mode = FragmentInsertMode::Replace;
                commit_fragment_replace(zmat, mol, settings, editor);
            }
            if let Some(error) = &zmat.last_error {
                return Err(error.clone());
            }
            // Keep the auto-targeted placed_editor intact. An old independent
            // target must not later restore coordinates from before insertion.
            let open = editor.window_open;
            *editor = EditorRotateState::default();
            editor.window_open = open;
            Ok(if matches!(action, Action::Add) {
                format!("Added {} to {target}", zmat.frag_name)
            } else {
                format!("Replaced {target} with {}", zmat.frag_name)
            })
        }
        Action::Delete(index) => {
            let symbol = mol
                .atoms
                .get(*index)
                .ok_or("Select an atom first.")?
                .clone();
            let open = editor.window_open;
            if !delete_atom_cartesian(mol, settings, zmat, editor, *index) {
                return Err("Could not delete this atom.".into());
            }
            editor.window_open = open;
            Ok(format!("Deleted {symbol}{}", index + 1))
        }
        Action::ChangeElement(index, symbol) => {
            if atomic_number(symbol).is_none() && !is_dummy(symbol) {
                return Err("Choose a recognised element.".into());
            }
            let atom = mol.atoms.get_mut(*index).ok_or("Select an atom first.")?;
            let old = std::mem::replace(atom, symbol.clone());
            mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
            zmat.zmat = zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos);
            set_selected_atom(zmat, &mol.atoms, Some(*index));
            Ok(format!("Changed {old}{} to {symbol}", index + 1))
        }
        Action::Hydrogens => {
            let before_count = mol.atoms.len();
            zmat.hydrogen_state.add(mol, settings);
            zmat.zmat = crate::molecule_builder::hydrogens::builder_zmat(mol, before_count);
            set_selected_atom(zmat, &mol.atoms, None);
            Ok(format!(
                "Added {} hydrogens",
                mol.atoms.len().saturating_sub(before_count)
            ))
        }
        Action::Cleanup => {
            let preview = zmat.cleanup.accepted(mol, &zmat.selection)?;
            mol.set_pos(preview.positions);
            mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
            let report = preview.report;
            zmat.zmat = zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos);
            zmat.cleanup_report = Some(Ok(report));
            set_selected_atom(zmat, &mol.atoms, None);
            Ok("Cleaned up geometry".into())
        }
        Action::References => {
            if zmat.add_atom_auto {
                commit_add_atom_auto(zmat, mol, settings);
            } else {
                commit_add_atom_with_refs(zmat, mol, settings);
            }
            if let Some(error) = &zmat.last_error {
                return Err(error.clone());
            }
            Ok(format!("Added {} using picked references", zmat.new_symbol))
        }
        Action::Undo | Action::Redo => unreachable!("handled before mutation"),
    }
}

pub(crate) fn insert_peptide(
    build: super::super::peptide::Build,
    editor: &mut EditorRotateState,
    zmat: &mut ZMatrixBuilderState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
    traj: &mut TrajectoryState,
) -> Result<BuilderChange, String> {
    if zmat.confirm_clear || zmat.add_atom_active || editor.active {
        return Err("Finish atom picking or the pending clear action first.".into());
    }
    perform(Action::Peptide(Box::new(build)), editor, zmat, mol, settings, traj)
        .map(|result| result.map_or(BuilderChange::None, |(change, _)| change))
}

#[cfg(test)]
#[path = "editor_panel_tests.rs"]
mod tests;
