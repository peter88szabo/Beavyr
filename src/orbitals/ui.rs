//! egui panel for the molecular orbital viewer.

use crate::ui_style::ResponseExt;
use bevy::prelude::*;
use bevy_egui::egui;

use super::density::Quantity;
use super::molden::{parse_molden, OrbitalData, Spin};
use super::{OrbitalState, SurfaceStyle};
use crate::events::MoleculeChanged;
use crate::molecule::Molecule;
use crate::settings::{MolSettings, RepresentationMode};

fn to_egui(c: Color) -> egui::Color32 {
    let s = c.to_srgba();
    egui::Color32::from_rgb(
        (s.red * 255.0) as u8,
        (s.green * 255.0) as u8,
        (s.blue * 255.0) as u8,
    )
}

fn from_egui(c: egui::Color32) -> Color {
    Color::srgb(
        c.r() as f32 / 255.0,
        c.g() as f32 / 255.0,
        c.b() as f32 / 255.0,
    )
}

/// Replace the displayed molecule with the geometry from the Molden file.
///
/// Emitting `parse_xyz` rather than `set_pos` is deliberate: it is the event
/// that means "a different structure", so the existing systems recentre the
/// camera and clear stale measurements without any new code here.
fn adopt_geometry(
    data: &OrbitalData,
    mol: &mut Molecule,
    ev_changed: &mut MessageWriter<MoleculeChanged>,
) {
    mol.atoms.clone_from(&data.atoms);
    mol.pos.clone_from(&data.positions);
    ev_changed.write(MoleculeChanged::parse_xyz(true));
}

/// Writes the Molden text behind a computed orbital set to a file the user
/// picks.
///
/// The text is kept in memory because the run directory it was written in is
/// discarded as soon as the run succeeds, so this is the only copy -- which is
/// exactly why the button exists.
fn save_molden(state: &mut OrbitalState) {
    let Some(text) = state.molden_text.clone() else {
        return;
    };
    let dialog = crate::recent_dir::save("canonicalMO.molden").add_filter("Molden", &["molden"]);
    let Some(path) = crate::recent_dir::save_file(dialog) else {
        return;
    };
    state.last_dir = path.parent().map(|p| p.to_path_buf());
    match std::fs::write(&path, text) {
        Ok(()) => {
            state.warning = None;
            // The file now exists, so it is named the way a loaded one is.
            state.file_name = path.file_name().map(|n| n.to_string_lossy().into_owned());
        }
        Err(err) => {
            state.warning = Some(format!("Could not write {}: {err}", path.display()));
        }
    }
}

/// Read a Molden file from a path we already have.
///
/// The dialog button and a filename given on the command line both come
/// through here, so an orbital file opens the same way however it arrived.
/// Returns whether it worked; on failure `state.warning` says why.
pub(crate) fn load_molden_from_path(
    path: &std::path::Path,
    state: &mut OrbitalState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
    ev_changed: &mut MessageWriter<MoleculeChanged>,
) -> bool {
    state.last_dir = path.parent().map(|p| p.to_path_buf());
    state.warning = None;
    match std::fs::read_to_string(path) {
        Ok(text) => match parse_molden(&text) {
            Ok(data) => {
                adopt_geometry(&data, mol, ev_changed);
                // Sticks keep the structure readable without filling the
                // volume the isosurface needs to occupy; solid spheres hide
                // the lobes they sit inside.
                settings.representation = RepresentationMode::SticksRounded;
                settings.geometry_dirty = true;
                let name = path.file_name().map(|n| n.to_string_lossy().into_owned());
                // A file the user opened needs no provenance line: its name is
                // the provenance, and the file is already on disk so there is
                // nothing to save.
                state.adopt(data, name, None, None);
                true
            }
            Err(error) => {
                state.warning = Some(format!("Could not read Molden file: {error}"));
                false
            }
        },
        Err(error) => {
            state.warning = Some(format!("Could not open file: {error}"));
            false
        }
    }
}

pub fn orbital_panel(
    ui: &mut egui::Ui,
    windowed: bool,
    state: &mut OrbitalState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
    ev_changed: &mut MessageWriter<MoleculeChanged>,
) -> bool {
    ui.scope(|ui| {
        crate::ui_style::modern(ui);
        panel_body(ui, windowed, state, mol, settings, ev_changed)
    })
    .inner
}

fn panel_body(
    ui: &mut egui::Ui,
    windowed: bool,
    state: &mut OrbitalState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
    ev_changed: &mut MessageWriter<MoleculeChanged>,
) -> bool {
    let mut loaded = false;
    ui.horizontal_wrapped(|ui| {
        if ui
            .add(crate::ui_style::primary("Load Molden…"))
            .clicked_once()
        {
            let dialog =
                crate::recent_dir::open().add_filter("Molden", &["molden", "molden.input", "inp"]);
            if let Some(path) = crate::recent_dir::pick_file(dialog) {
                loaded = load_molden_from_path(&path, state, mol, settings, ev_changed);
            }
        }
        if state.data.is_some() && ui.button("Clear orbitals").clicked_once() {
            state.clear();
        }
        if state.is_busy() {
            ui.spinner().on_hover_text("Sampling surface…");
        }
        if let (Some(data), Some(name)) = (&state.data, &state.file_name) {
            ui.add(egui::Label::new(egui::RichText::new(name).small().weak()).truncate())
                .on_hover_text(format!(
                    "{name}\n{} basis functions · {}",
                    data.layout.n_basis,
                    if data.is_open_shell() {
                        "open shell"
                    } else {
                        "closed shell"
                    },
                ));
        }
    });

    if let Some(warning) = state.warning.clone() {
        ui.add_space(4.0);
        ui.colored_label(egui::Color32::from_rgb(220, 80, 80), warning);
    }

    let Some(data) = state.data.clone() else {
        ui.add_space(4.0);
        ui.weak("No orbital file loaded.");
        return loaded;
    };

    // A computed set has no file behind it, so it says where it came from and
    // offers to become one. Without this the panel would look as though a file
    // had been loaded, with no way to tell which run these orbitals belong to
    // or to keep them once the scratch directory is gone.
    if let Some(provenance) = state.provenance.clone() {
        ui.horizontal_wrapped(|ui| {
            ui.colored_label(egui::Color32::from_rgb(120, 170, 220), provenance)
                .on_hover_text(format!(
                    "{} basis functions · {}",
                    data.layout.n_basis,
                    if data.is_open_shell() {
                        "open shell"
                    } else {
                        "closed shell"
                    }
                ));
            if state.molden_text.is_some() && ui.button("Save .molden\u{2026}").clicked_once() {
                save_molden(state);
            }
        });
    }
    ui.horizontal_wrapped(|ui| {
        let selected_fill = ui.visuals().selection.bg_fill;
        let dim = ui.visuals().widgets.inactive.bg_fill;
        for quantity in [Quantity::Orbital, Quantity::Density, Quantity::SpinDensity] {
            let available = quantity.is_available(&data);
            let fill = if state.quantity == quantity {
                selected_fill
            } else {
                dim
            };
            let button = ui.add_enabled(available, egui::Button::new(quantity.label()).fill(fill));
            if !available {
                button.on_hover_text("Spin density needs an open-shell (alpha and beta) file.");
            } else if button.clicked_once() && state.quantity != quantity {
                state.quantity = quantity;
                state.isovalue = quantity.default_isovalue();
                state.mark_surface_dirty();
            }
        }
        // Spin chooses an orbital set; densities combine both sets.
        if data.is_open_shell() && matches!(state.quantity, Quantity::Orbital) {
            for (spin, label) in [(Spin::Alpha, "Alpha"), (Spin::Beta, "Beta")] {
                if ui
                    .selectable_label(state.spin == spin, label)
                    .on_hover_text("Choose the spin set to browse")
                    .clicked_once()
                    && state.spin != spin
                {
                    state.spin = spin;
                    state.selected = data.homo_index(spin);
                }
            }
        }
    });
    let showing_orbital = matches!(state.quantity, Quantity::Orbital);

    let orbitals = data.orbitals(state.spin);
    let homo = data.homo_index(state.spin);

    let mut surface_changed = false;
    crate::ui_style::group(ui, "", |ui| {
        // Useful isovalues span very different magnitudes for an orbital, a total
        // density and a spin density, so the slider bounds follow the quantity.
        let (isovalue_low, isovalue_high) = state.quantity.isovalue_range();
        surface_changed |= ui
            .add(
                egui::Slider::new(&mut state.isovalue, isovalue_low..=isovalue_high)
                    .logarithmic(true)
                    .text("Isovalue"),
            )
            .on_hover_text(if showing_orbital {
                "Lower encloses more of the orbital; higher shrinks to the core."
            } else {
                "Lower encloses more of the distribution; higher shrinks to the densest region."
            })
            .changed();

        ui.horizontal_wrapped(|ui| {
            surface_changed |= ui
                .checkbox(&mut state.show_surface, "Show surface")
                .changed();
            surface_changed |= ui
                .add(egui::Slider::new(&mut state.opacity, 0.1..=1.0).text("Opacity"))
                .changed();
        });
        if state.field_peak().is_some_and(|peak| state.isovalue > peak) {
            ui.colored_label(
                egui::Color32::from_rgb(220, 160, 60),
                "Isovalue exceeds the peak; nothing will be drawn.",
            );
        }
    });
    egui::CollapsingHeader::new(egui::RichText::new("Surface detail and colors").strong())
        .default_open(false)
        .show(ui, |ui| {
            // Keep advanced settings from pushing the orbital list out of the
            // floating window when expanded.
            egui::ScrollArea::vertical()
                .id_salt("surface_detail")
                .max_height(if windowed { 140.0 } else { 360.0 })
                .show(ui, |ui| {
                    ui.add_space(4.0);
                    // Resolution changes the sampling density, so unlike the controls above it
                    // forces a resample rather than a remesh.  Showing the grid it produces
                    // makes that cost visible before it is paid.
                    ui.add(
                        egui::Slider::new(
                            &mut state.resolution,
                            crate::orbitals::MIN_RESOLUTION..=crate::orbitals::MAX_RESOLUTION,
                        )
                        .text("Surface resolution"),
                    );
                    let planned = state.planned_grid(&data);
                    ui.label(
                        egui::RichText::new(format!(
                            "{:.2} Å grid, {}×{}×{} = {} points",
                            planned.step,
                            planned.dims[0],
                            planned.dims[1],
                            planned.dims[2],
                            planned.point_count()
                        ))
                        .small()
                        .weak(),
                    );
                    ui.label(
                        egui::RichText::new(
                            "Higher resolution gives a smoother surface but is slower to compute.",
                        )
                        .small()
                        .weak(),
                    );

                    ui.add_space(6.0);
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Style:");
                        let selected = ui.visuals().selection.bg_fill;
                        let dim = ui.visuals().widgets.inactive.bg_fill;
                        for (style, label) in
                            [(SurfaceStyle::Solid, "Solid"), (SurfaceStyle::Mesh, "Mesh")]
                        {
                            let fill = if state.style == style { selected } else { dim };
                            if ui.add(egui::Button::new(label).fill(fill)).clicked_once()
                                && state.style != style
                            {
                                state.style = style;
                                surface_changed = true;
                            }
                        }
                    });

                    // Mesh controls only exist while the mesh style is active.
                    if matches!(state.style, SurfaceStyle::Mesh) {
                        surface_changed |= ui
                            .add(
                                egui::Slider::new(
                                    &mut state.mesh_density,
                                    crate::orbitals::MIN_MESH_DENSITY
                                        ..=crate::orbitals::MAX_MESH_DENSITY,
                                )
                                .text("Mesh grid"),
                            )
                            .changed();
                        ui.label(
                            egui::RichText::new(format!(
                                "Higher is denser: contour planes every {:.2} Å.",
                                state.mesh_plane_spacing()
                            ))
                            .small()
                            .weak(),
                        );
                        // Width is a gizmo setting, so it needs no rebuild.
                        ui.add(
                            egui::Slider::new(&mut state.mesh_width, 0.5..=8.0)
                                .text("Mesh line width"),
                        );
                    }

                    ui.horizontal_wrapped(|ui| {
                        ui.label("Lobe colours:");
                        let mut positive = to_egui(state.positive_color);
                        if ui.color_edit_button_srgba(&mut positive).changed() {
                            state.positive_color = from_egui(positive);
                            surface_changed = true;
                        }
                        let mut negative = to_egui(state.negative_color);
                        if ui.color_edit_button_srgba(&mut negative).changed() {
                            state.negative_color = from_egui(negative);
                            surface_changed = true;
                        }
                    });
                    if let Some(peak) = state.field_peak() {
                        let unit = match state.style {
                            SurfaceStyle::Solid => "triangles",
                            SurfaceStyle::Mesh => "segments",
                        };
                        ui.weak(format!(
                            "Peak value {peak:.4}   surface {} {unit}",
                            state.triangle_count
                        ));
                        if state.isovalue > peak {
                            ui.colored_label(
                                egui::Color32::from_rgb(220, 160, 60),
                                "Isovalue exceeds the peak; nothing will be drawn.",
                            );
                        }
                    }
                });
        });
    ui.add_space(2.0);
    // The orbital list is only meaningful when plotting one orbital.
    if showing_orbital {
        crate::ui_style::group(ui, "", |ui| {
            ui.spacing_mut().item_spacing.y = 2.0;
            ui.spacing_mut().button_padding.y = 2.0;
            ui.spacing_mut().interact_size.y = 24.0;
            // A single header carries the column meanings so the rows themselves need
            // no repeated labels.  Header and rows share one format spec, so they
            // cannot drift out of alignment.
            ui.style_mut().override_text_style = Some(egui::TextStyle::Monospace);
            ui.weak(format!("{:>4} {:<5}{:>10}{:>7}", "#", "sym", "E/Eh", "occ"));
            ui.style_mut().override_text_style = None;

            let row_height = (ui.text_style_height(&egui::TextStyle::Monospace) + 4.0)
                .max(ui.spacing().interact_size.y);
            let list_height = if windowed {
                ui.available_height().max(row_height)
            } else {
                360.0
            };
            egui::ScrollArea::vertical()
                .id_salt("surface_orbitals")
                .max_height(list_height)
                .auto_shrink([false; 2])
                // Virtualised: 546 orbitals must not cost 546 laid-out widgets a frame.
                .show_rows(ui, row_height, orbitals.len(), |ui, range| {
                    ui.style_mut().override_text_style = Some(egui::TextStyle::Monospace);
                    for index in range {
                        let orbital = &orbitals[index];
                        let marker = match homo {
                            Some(h) if index == h => "  HOMO",
                            Some(h) if index == h + 1 => "  LUMO",
                            _ => "",
                        };
                        let label = format!(
                            "{:>4} {:<5}{:>10.4}{:>7.2}{}",
                            index + 1,
                            orbital.symmetry,
                            orbital.energy,
                            orbital.occupation,
                            marker
                        );
                        if ui
                            .selectable_label(state.selected == Some(index), label)
                            .clicked_once()
                        {
                            state.selected = Some(index);
                        }
                    }
                });
        });
    }

    if surface_changed {
        // Only the mesh changes; the sampled field is still valid.
        state.mark_surface_dirty();
    }
    loaded
}
