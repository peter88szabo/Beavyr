use super::*;
#[allow(clippy::too_many_arguments)]
pub fn window(
    ctx: &egui::Context,
    state: &mut State,
    traj: &mut TrajectoryState,
    mol: &mut crate::molecule::Molecule,
    settings: &mut crate::settings::MolSettings,
    events: &mut MessageWriter<crate::events::MoleculeChanged>,
    builder: &crate::molecule_builder::builder_ui::ZMatrixBuilderState,
    qc: &crate::qchem_panel::QcPanelState,
    layout: &mut crate::ui_layout::UiLayout,
    regions: &mut crate::ui::UiPanelRegions,
) {
    if !state.plot_open {
        return;
    }
    let Some(result) = state.results.clone() else {
        return;
    };
    let mut open = true;
    let window = egui::Window::new("Reactivity plot")
        .open(&mut open)
        .resizable(true)
        .default_size([760., 520.])
        .show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                egui::ComboBox::from_id_salt("reactivity_plot_quantity")
                    .selected_text(state.descriptor.label())
                    .show_ui(ui, |ui| {
                        for d in data::Descriptor::ALL {
                            ui.selectable_value(&mut state.descriptor, d, d.label());
                        }
                    });
                if state.descriptor.local() {
                    ui.label("Atoms");
                    ui.text_edit_singleline(&mut state.group);
                    if ui.button("Selected atoms").clicked_once() {
                        state.group = builder
                            .selection
                            .atoms
                            .iter()
                            .map(|i| (i + 1).to_string())
                            .collect::<Vec<_>>()
                            .join(",");
                    }
                    if ui.button("All").clicked_once() {
                        state.group = format!("1-{}", result.frames[0].geometry.atoms.len());
                    }
                }
            });
            use crate::polarizability::data::Axis;
            ui.horizontal_wrapped(|ui| {
                for axis in [
                    Axis::Frame,
                    Axis::Time,
                    Axis::Distance,
                    Axis::Angle,
                    Axis::Dihedral,
                ] {
                    ui.selectable_value(&mut state.axis, axis, axis.label());
                }
            });
            if state.axis == Axis::Time {
                ui.horizontal(|ui| {
                    ui.label("Time per sampled frame (fs)");
                    ui.add(egui::DragValue::new(&mut state.dt).range(0.000001..=1e12));
                    ui.label("First time (fs)");
                    ui.add(egui::DragValue::new(&mut state.origin));
                });
            }
            if state.axis.count() > 0 {
                ui.horizontal(|ui| {
                    ui.label("Coordinate atoms");
                    for a in state.coordinate.iter_mut().take(state.axis.count()) {
                        ui.add(
                            egui::DragValue::new(a)
                                .range(1..=result.frames[0].geometry.atoms.len()),
                        );
                    }
                    if ui.button("Use scan coordinate").clicked_once() {
                        let scan = &qc.config.scan;
                        use crate::qchem_interfaces::job::scan::{self, Coordinate};
                        match scan::atom_indices(
                            &scan.atoms,
                            scan.coordinate.count(),
                            &result.frames[0].geometry.atoms,
                        ) {
                            Ok(indices) => {
                                state.axis = match scan.coordinate {
                                    Coordinate::Bond => Axis::Distance,
                                    Coordinate::Angle => Axis::Angle,
                                    Coordinate::Dihedral => Axis::Dihedral,
                                };
                                for (a, i) in state.coordinate.iter_mut().zip(indices) {
                                    *a = i + 1;
                                }
                            }
                            Err(e) => state.error = Some(e),
                        }
                    }
                });
            }
            let group = if state.descriptor.local() {
                match crate::polarizability::data::atom_list(
                    &state.group,
                    result.frames[0].geometry.atoms.len(),
                ) {
                    Ok(v) => v,
                    Err(e) => {
                        ui.label(e);
                        return;
                    }
                }
            } else {
                vec![]
            };
            let points: Vec<_> = result
                .frames
                .iter()
                .enumerate()
                .map(|(i, f)| {
                    Some((
                        state.axis.value(
                            &f.geometry,
                            i,
                            &state.coordinate,
                            state.dt,
                            state.origin,
                        )?,
                        state.descriptor.value(f, &group)?,
                    ))
                })
                .collect();
            ui.horizontal(|ui| {
                let csv = || {
                    let mut s = format!(
                        "# {}\n# {} ; atoms: {}\nframe,{},{} ({})\n",
                        result.method,
                        state.descriptor.label(),
                        state.group,
                        state.axis.label(),
                        state.descriptor.label(),
                        state.descriptor.unit()
                    );
                    for (i, p) in points.iter().enumerate() {
                        if let Some((x, y)) = p {
                            s.push_str(&format!("{},{x:.10},{y:.10}\n", i + 1));
                        }
                    }
                    s
                };
                if ui.button("Copy data").clicked_once() {
                    ctx.copy_text(csv());
                }
                if ui.button("Save CSV…").clicked_once() {
                    if let Some(p) = rfd::FileDialog::new()
                        .set_file_name("reactivity.csv")
                        .save_file()
                    {
                        if let Err(e) = std::fs::write(p, csv()) {
                            state.error = Some(e.to_string());
                        }
                    }
                }
            });
            let matches = result.frames.len() == traj.frames.len()
                && result.frames.iter().zip(&traj.frames).all(|(a, b)| {
                    crate::trajectory_wavefunctions::job::geometry_matches(&a.geometry, b)
                });
            if !matches {
                ui.label("Load the matching trajectory to select frames from this plot.");
            }
            let title = format!(
                "{}{}",
                if state.descriptor.local() {
                    "Atomic sum: "
                } else {
                    ""
                },
                state.descriptor.label()
            );
            if let Some(i) = crate::polarizability::plot::draw_named(
                ui,
                &points,
                if matches {
                    Some(traj.current_frame)
                } else {
                    None
                },
                state.axis.label(),
                state.descriptor.unit(),
                &title,
            ) {
                if matches {
                    traj.playing = false;
                    traj.current_frame = i;
                    crate::trajectory::apply_current_frame(
                        traj, mol, settings, events, false, None,
                    );
                    layout.open[crate::ui_layout::Tab::Trajectory.index()] = true;
                }
            }
        });
    state.plot_open = open;
    if let Some(w) = window {
        super::register(regions, w.response.rect);
    }
}
