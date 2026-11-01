//! Frame-indexed trajectory analysis. Plot clicks preserve the camera and select the original frame.
use crate::{
    trajectory::{TrajectoryFrame, TrajectoryState},
    ui_style::ResponseExt,
};
use bevy::{
    prelude::*,
    tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task},
};
use bevy_egui::{egui, EguiContexts, EguiPrimaryContextPass};
#[path = "trajectory_plot.rs"]
mod plot;
const OPEN: &str = "open_trajectory_analysis";
const ENERGY_SELECTION: &str = "trajectory_energy_selection";
#[derive(Clone, Default)]
struct EnergySelection {
    energies: Vec<f64>,
    index: usize,
}
pub fn energy_point_clicked(ctx: &egui::Context, energies: Vec<f64>, index: usize) {
    ctx.data_mut(|d| {
        d.insert_temp(
            egui::Id::new(ENERGY_SELECTION),
            EnergySelection { energies, index },
        )
    });
}
fn mapped_energy_frame(frames: &[TrajectoryFrame], request: &EnergySelection) -> Option<usize> {
    let available: Vec<_> = frames
        .iter()
        .enumerate()
        .filter_map(|(i, f)| {
            if f.metadata.energy_is_relative {
                None
            } else {
                f.metadata.energy_hartree.map(|e| (i, e))
            }
        })
        .collect();
    if available.len() != request.energies.len()
        || available
            .iter()
            .zip(&request.energies)
            .any(|((_, a), b)| (a - b).abs() > 1e-7)
    {
        return None;
    }
    available.get(request.index).map(|&(index, _)| index)
}
pub struct Plugin;
impl bevy::prelude::Plugin for Plugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<State>()
            .add_systems(EguiPrimaryContextPass, ui.after(crate::ui::ui_panel));
    }
}
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Quantity {
    Energy,
    Distance,
    Angle,
    Dihedral,
    Rmsd,
}
impl Quantity {
    fn label(self) -> &'static str {
        match self {
            Self::Energy => "Energy",
            Self::Distance => "Distance",
            Self::Angle => "Angle",
            Self::Dihedral => "Dihedral",
            Self::Rmsd => "RMSD (aligned)",
        }
    }
    fn atoms(self) -> usize {
        match self {
            Self::Distance => 2,
            Self::Angle => 3,
            Self::Dihedral => 4,
            _ => 0,
        }
    }
}
#[derive(Clone, Copy, PartialEq, Eq)]
enum Unit {
    Kcal,
    Kj,
    Hartree,
}
impl Unit {
    fn label(self) -> &'static str {
        match self {
            Self::Kcal => "kcal/mol",
            Self::Kj => "kJ/mol",
            Self::Hartree => "Hartree",
        }
    }
    fn factor(self) -> f64 {
        match self {
            Self::Kcal => 627.5094740631,
            Self::Kj => 2625.499639479163,
            Self::Hartree => 1.,
        }
    }
}
#[derive(Clone, Copy, PartialEq, Eq)]
enum Reference {
    First,
    Minimum,
    Last,
}
struct Series {
    values: Vec<Option<f64>>,
    fingerprint: u64,
    quantity: Quantity,
    atoms: [usize; 4],
    reference_frame: usize,
}
#[derive(Resource)]
struct State {
    open: bool,
    quantity: Quantity,
    atoms: [usize; 4],
    reference_frame: usize,
    unit: Unit,
    reference: Reference,
    histogram: bool,
    bins: usize,
    data: Option<Series>,
    task: Option<Task<Result<Series, String>>>,
}
impl Default for State {
    fn default() -> Self {
        Self {
            open: false,
            quantity: Quantity::Energy,
            atoms: [1, 2, 3, 4],
            reference_frame: 1,
            unit: Unit::Kcal,
            reference: Reference::First,
            histogram: false,
            bins: 30,
            data: None,
            task: None,
        }
    }
}
pub fn button(ui: &mut egui::Ui) {
    if ui.button("Trajectory analysis…").clicked_once() {
        ui.ctx()
            .data_mut(|d| d.insert_temp(egui::Id::new(OPEN), true));
    }
}
fn compute(
    frames: &[TrajectoryFrame],
    quantity: Quantity,
    atoms: [usize; 4],
    reference_frame: usize,
) -> Result<Series, String> {
    if frames.is_empty() {
        return Err("Load a trajectory first.".into());
    }
    if quantity == Quantity::Energy {
        let mut kinds = frames
            .iter()
            .filter(|f| f.metadata.energy_hartree.is_some())
            .map(|f| f.metadata.energy_is_relative);
        if let Some(first) = kinds.next() {
            if kinds.any(|k| k != first) {
                return Err(
                    "These frames mix absolute and relative energies. Analyse them separately."
                        .into(),
                );
            }
        }
    }
    let indices: Vec<_> = atoms[..quantity.atoms()]
        .iter()
        .map(|i| i.checked_sub(1).ok_or("Atom numbers start at 1."))
        .collect::<Result<_, _>>()?;
    let reference = frames
        .get(if quantity == Quantity::Rmsd {
            reference_frame.saturating_sub(1)
        } else {
            0
        })
        .ok_or("Choose a valid RMSD reference frame.")?;
    let values = frames
        .iter()
        .map(|frame| {
            if frame.atoms.len() != frame.pos.len() || frame.pos.iter().any(|p| !p.is_finite()) {
                return None;
            }
            if quantity == Quantity::Energy {
                return frame.metadata.energy_hartree.filter(|v| v.is_finite());
            }
            if quantity == Quantity::Rmsd {
                if frame.atoms != reference.atoms {
                    return None;
                }
                return crate::molecular_alignment::superpose(&frame.pos, &reference.pos)
                    .ok()
                    .and_then(|s| {
                        crate::rmsd::kabsch::rmsd(&s.apply_all(&frame.pos), &reference.pos).ok()
                    });
            }
            let p: Vec<Vec3> = indices
                .iter()
                .map(|&i| frame.pos.get(i).copied())
                .collect::<Option<_>>()?;
            match quantity {
                Quantity::Distance => Some(p[0].distance(p[1]) as f64),
                Quantity::Angle => {
                    let a = (p[0] - p[1]).try_normalize()?;
                    let b = (p[2] - p[1]).try_normalize()?;
                    Some(a.dot(b).clamp(-1., 1.).acos().to_degrees() as f64)
                }
                Quantity::Dihedral => {
                    let axis = (p[2] - p[1]).try_normalize()?;
                    let a = p[0] - p[1];
                    let b = p[3] - p[2];
                    let a = (a - axis * a.dot(axis)).try_normalize()?;
                    let b = (b - axis * b.dot(axis)).try_normalize()?;
                    Some(axis.cross(a).dot(b).atan2(a.dot(b)).to_degrees() as f64)
                }
                _ => None,
            }
        })
        .collect::<Vec<_>>();
    if values.iter().all(Option::is_none) {
        return Err(if quantity==Quantity::Energy{"No frame energies were found. Energies need explicit units in XYZ comments (for example E = -76.0 Eh), or a supported calculation result."}else{"No valid measurements could be computed. Check atom numbers and geometry."}.into());
    }
    Ok(Series {
        values,
        fingerprint: crate::trajectory_data::fingerprint(frames),
        quantity,
        atoms,
        reference_frame,
    })
}
fn displayed(data: &Series, unit: Unit, reference: Reference) -> (Vec<Option<f64>>, String) {
    if data.quantity != Quantity::Energy {
        return (
            data.values.clone(),
            if matches!(data.quantity, Quantity::Angle | Quantity::Dihedral) {
                "degrees"
            } else {
                "Å"
            }
            .into(),
        );
    }
    let available: Vec<_> = data.values.iter().flatten().copied().collect();
    let zero = match reference {
        Reference::First => available[0],
        Reference::Last => *available.last().unwrap(),
        Reference::Minimum => available.into_iter().fold(f64::INFINITY, f64::min),
    };
    (
        data.values
            .iter()
            .map(|v| v.map(|v| (v - zero) * unit.factor()))
            .collect(),
        unit.label().into(),
    )
}
fn csv(values: &[Option<f64>], quantity: Quantity, unit: &str) -> String {
    use std::fmt::Write;
    let mut text = format!("frame,{} ({})\n", quantity.label(), unit);
    for (i, v) in values.iter().enumerate() {
        writeln!(
            text,
            "{},{}",
            i + 1,
            v.map(|v| format!("{v:.10}")).unwrap_or_default()
        )
        .unwrap();
    }
    text
}
fn ui(
    mut contexts: EguiContexts,
    mut state: ResMut<State>,
    mut traj: ResMut<TrajectoryState>,
    mut mol: ResMut<crate::molecule::Molecule>,
    mut settings: ResMut<crate::settings::MolSettings>,
    mut events: MessageWriter<crate::events::MoleculeChanged>,
    mut regions: ResMut<crate::ui::UiPanelRegions>,
) {
    let Ok(ctx) = contexts.ctx_mut() else {
        return;
    };
    if let Some(request) =
        ctx.data_mut(|d| d.remove_temp::<EnergySelection>(egui::Id::new(ENERGY_SELECTION)))
    {
        if let Some(index) = mapped_energy_frame(&traj.frames, &request) {
            traj.current_frame = index;
            traj.playing = false;
            crate::trajectory::apply_current_frame(
                &mut traj,
                &mut mol,
                &mut settings,
                &mut events,
                false,
                None,
            );
        } else {
            crate::runtime_warning::push(ctx,"Energy plot","This energy history does not match the loaded trajectory. Load its trajectory first, or use Trajectory analysis for the current frames.");
        }
    }
    if ctx
        .data_mut(|d| d.remove_temp::<bool>(egui::Id::new(OPEN)))
        .unwrap_or(false)
    {
        state.open = true;
    }
    if let Some(task) = &mut state.task {
        if let Some(result) = block_on(future::poll_once(task)) {
            state.task = None;
            match result {
                Ok(s) => state.data = Some(s),
                Err(e) => crate::runtime_warning::push(ctx, "Trajectory analysis", &e),
            }
        } else {
            ctx.request_repaint();
        }
    }
    if !state.open {
        return;
    }
    let mut open = true;
    let mut selected = None;
    let window = egui::Window::new("Trajectory analysis")
        .open(&mut open)
        .resizable(true)
        .default_size([720., 520.])
        .show(ctx, |ui| {
            crate::ui_style::modern(ui);
            ui.horizontal_wrapped(|ui| {
                egui::ComboBox::from_id_salt("trajectory_quantity")
                    .selected_text(state.quantity.label())
                    .show_ui(ui, |ui| {
                        for q in [
                            Quantity::Energy,
                            Quantity::Distance,
                            Quantity::Angle,
                            Quantity::Dihedral,
                            Quantity::Rmsd,
                        ] {
                            ui.selectable_value(&mut state.quantity, q, q.label());
                        }
                    });
                for i in 0..state.quantity.atoms() {
                    ui.label(["A", "B", "C", "D"][i]);
                    ui.add(
                        egui::DragValue::new(&mut state.atoms[i]).range(1..=mol.atoms.len().max(1)),
                    );
                }
                if state.quantity == Quantity::Rmsd {
                    ui.label("Reference frame");
                    ui.add(
                        egui::DragValue::new(&mut state.reference_frame)
                            .range(1..=traj.frames.len().max(1)),
                    );
                }
                if ui
                    .add_enabled(
                        state.task.is_none() && !traj.frames.is_empty(),
                        egui::Button::new("Plot / refresh"),
                    )
                    .clicked_once()
                {
                    let frames = traj.frames.clone();
                    let q = state.quantity;
                    let a = state.atoms;
                    let r = state.reference_frame;
                    state.task = Some(
                        AsyncComputeTaskPool::get().spawn(async move { compute(&frames, q, a, r) }),
                    );
                }
            });
            if state.task.is_some() {
                ui.spinner();
            }
            if state
                .data
                .as_ref()
                .is_some_and(|d| d.quantity == Quantity::Energy)
            {
                ui.horizontal_wrapped(|ui| {
                    for u in [Unit::Kcal, Unit::Kj, Unit::Hartree] {
                        ui.selectable_value(&mut state.unit, u, u.label());
                    }
                    ui.add_space(16.);
                    ui.label("Reference");
                    for (r, l) in [
                        (Reference::First, "First"),
                        (Reference::Minimum, "Minimum"),
                        (Reference::Last, "Last"),
                    ] {
                        ui.selectable_value(&mut state.reference, r, l);
                    }
                });
            }
            ui.horizontal_wrapped(|ui| {
                ui.checkbox(&mut state.histogram, "Distribution");
                if state.histogram {
                    ui.label("Bins");
                    ui.add(egui::DragValue::new(&mut state.bins).range(2..=200));
                }
            });
            if let Some(data) = &state.data {
                let (values, unit) = displayed(data, state.unit, state.reference);
                let title = if data.quantity.atoms() > 0 {
                    format!(
                        "{} {:?} ({unit})",
                        data.quantity.label(),
                        &data.atoms[..data.quantity.atoms()]
                    )
                } else if data.quantity == Quantity::Rmsd {
                    format!("RMSD to frame {} ({unit})", data.reference_frame)
                } else {
                    format!("Relative energy ({unit})")
                };
                ui.label(&title);
                ui.horizontal(|ui| {
                    if ui.button("Copy CSV").clicked_once() {
                        ctx.copy_text(csv(&values, data.quantity, &unit));
                    }
                    if ui.button("Save CSV…").clicked_once() {
                        if let Some(path) = rfd::FileDialog::new()
                            .set_file_name("trajectory-analysis.csv")
                            .save_file()
                        {
                            if let Err(e) = crate::project::atomic_write(
                                &path,
                                csv(&values, data.quantity, &unit).as_bytes(),
                            ) {
                                crate::runtime_warning::push(
                                    ctx,
                                    "Save trajectory analysis",
                                    &e.to_string(),
                                );
                            }
                        }
                    }
                });
                if state.histogram {
                    plot::histogram(ui, &values, state.bins, &unit);
                } else {
                    selected = plot::series(ui, &values, traj.current_frame, &unit);
                }
                let missing = values.iter().filter(|v| v.is_none()).count();
                if missing > 0 {
                    ui.small(format!(
                        "{missing} frames have no valid value; their frame numbers are retained."
                    ));
                }
            }
        });
    state.open = open;
    if let Some(window) = window {
        let r = window.response.rect;
        regions
            .windows
            .push(Rect::new(r.min.x, r.min.y, r.max.x, r.max.y));
    }
    if let Some(frame) = selected {
        if state
            .data
            .as_ref()
            .is_some_and(|d| d.fingerprint == crate::trajectory_data::fingerprint(&traj.frames))
        {
            traj.current_frame = frame;
            traj.playing = false;
            crate::trajectory::apply_current_frame(
                &mut traj,
                &mut mol,
                &mut settings,
                &mut events,
                false,
                None,
            );
        } else {
            crate::runtime_warning::push(
                ctx,
                "Trajectory analysis",
                "The trajectory changed. Refresh the plot before selecting a frame.",
            );
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn existing_energy_plots_map_only_to_the_matching_trajectory() {
        let f = frames();
        let request = EnergySelection {
            energies: vec![-2., -1.],
            index: 1,
        };
        assert_eq!(mapped_energy_frame(&f, &request), Some(2));
        assert_eq!(
            mapped_energy_frame(
                &f,
                &EnergySelection {
                    energies: vec![-3., -1.],
                    index: 1
                }
            ),
            None
        );
    }
    fn frames() -> Vec<TrajectoryFrame> {
        crate::trajectory::parse_multi_xyz("2\nE = -2 Eh\nH 0 0 0\nH 1 0 0\n2\nmissing\nH 0 0 0\nH 2 0 0\n2\nE = -1 Eh\nH 0 0 0\nH 3 0 0\n").unwrap()
    }
    #[test]
    fn energy_gaps_keep_frame_indices_and_reference_conversion() {
        let s = compute(&frames(), Quantity::Energy, [1, 2, 3, 4], 1).unwrap();
        assert_eq!(s.values, vec![Some(-2.), None, Some(-1.)]);
        let (v, _) = displayed(&s, Unit::Kcal, Reference::First);
        assert_eq!(v[0], Some(0.));
        assert!((v[2].unwrap() - 627.5094740631).abs() < 1e-9);
        assert!(csv(&v, Quantity::Energy, "kcal/mol").contains("2,\n3,"));
    }
    #[test]
    fn distance_measurements_and_aligned_rmsd_are_valid() {
        let f = frames();
        let d = compute(&f, Quantity::Distance, [1, 2, 3, 4], 1).unwrap();
        assert_eq!(d.values, vec![Some(1.), Some(2.), Some(3.)]);
        let r = compute(&f, Quantity::Rmsd, [1, 2, 3, 4], 1).unwrap();
        assert!(r.values[0].unwrap() < 1e-8);
        assert!((r.values[1].unwrap() - 0.5).abs() < 1e-6);
    }
}
