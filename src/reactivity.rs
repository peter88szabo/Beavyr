//! Surface Tools reactivity analysis: fixed-geometry differences and frontier approximations.
use crate::{
    trajectory::{TrajectoryFrame, TrajectoryState},
    ui_style::ResponseExt,
};
use bevy::{
    prelude::*,
    tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task},
};
use bevy_egui::{egui, EguiContexts, EguiPrimaryContextPass};
use std::{
    path::PathBuf,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
};
#[path = "reactivity_data.rs"]
pub mod data;
#[path = "reactivity_job.rs"]
mod job;
#[path = "reactivity_plot.rs"]
mod plotting;
#[path = "reactivity_render.rs"]
mod render;
struct Running {
    task: Task<Result<data::Results, String>>,
    control: crate::trajectory_wavefunctions::job::Control,
    progress: Arc<AtomicUsize>,
    total: usize,
}
impl Drop for Running {
    fn drop(&mut self) {
        self.control.cancel();
    }
}
#[derive(Resource)]
pub struct State {
    pub results: Option<Arc<data::Results>>,
    pub descriptor: data::Descriptor,
    pub halos: bool,
    pub labels: bool,
    pub surface: bool,
    pub opacity: f32,
    pub colors: crate::property_legend::Scale,
    pub scale: f32,
    pub isovalue: f32,
    pub resolution: u32,
    pub positive: Color,
    pub negative: Color,
    pub ready: Option<usize>,
    open: bool,
    plot_open: bool,
    approx: bool,
    use_trajectory: bool,
    spin: crate::orbitals::molden::Spin,
    charge: i32,
    mult: [i32; 3],
    cores: u32,
    running: Option<Running>,
    error: Option<String>,
    directory: Option<PathBuf>,
    pub group: String,
    pub axis: crate::polarizability::data::Axis,
    pub coordinate: [usize; 4],
    pub dt: f64,
    pub origin: f64,
}
impl Default for State {
    fn default() -> Self {
        Self {
            results: None,
            descriptor: data::Descriptor::Dual,
            halos: false,
            labels: false,
            surface: false,
            opacity: 0.35,
            colors: Default::default(),
            scale: 1.5,
            isovalue: 0.005,
            resolution: 7,
            positive: Color::srgb(0.2, 0.5, 1.),
            negative: Color::srgb(0.95, 0.25, 0.2),
            ready: None,
            open: false,
            plot_open: false,
            approx: true,
            use_trajectory: false,
            spin: crate::orbitals::molden::Spin::Alpha,
            charge: 0,
            mult: [1, 2, 2],
            cores: 4,
            running: None,
            error: None,
            directory: None,
            group: "1".into(),
            axis: crate::polarizability::data::Axis::Frame,
            coordinate: [1, 2, 3, 4],
            dt: 1.,
            origin: 0.,
        }
    }
}
impl State {
    pub fn is_running(&self) -> bool {
        self.running.is_some()
    }
    pub fn current_index(
        &self,
        mol: &crate::molecule::Molecule,
        traj: &TrajectoryState,
    ) -> Option<usize> {
        let r = self.results.as_ref()?;
        let i = if r.frames.len() == 1 {
            0
        } else {
            traj.current_frame
        };
        let f = &r.frames.get(i)?.geometry;
        if f.atoms != mol.atoms
            || f.pos.len() != mol.pos.len()
            || f.pos
                .iter()
                .zip(&mol.pos)
                .any(|(a, b)| a.distance(*b) > 1e-5)
        {
            return None;
        }
        Some(i)
    }
    pub fn hold_playback(&self, mol: &crate::molecule::Molecule, traj: &TrajectoryState) -> bool {
        self.surface
            && self.descriptor.local()
            && self
                .current_index(mol, traj)
                .is_some_and(|i| self.ready != Some(i))
    }
}
pub struct Plugin;
impl bevy::prelude::Plugin for Plugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<State>()
            .add_plugins(render::Plugin)
            .add_systems(EguiPrimaryContextPass, ui.after(crate::ui::ui_panel));
    }
}
pub fn button(ui: &mut egui::Ui) {
    if ui.button("Reactivity / Fukui…").clicked_once() {
        ui.ctx()
            .data_mut(|d| d.insert_temp(egui::Id::new("open_reactivity"), true));
    }
}
fn register(regions: &mut crate::ui::UiPanelRegions, r: egui::Rect) {
    regions
        .windows
        .push(Rect::new(r.min.x, r.min.y, r.max.x, r.max.y));
}
#[allow(clippy::too_many_arguments)]
fn ui(
    mut contexts: EguiContexts,
    mut state: ResMut<State>,
    wave: Res<crate::trajectory_wavefunctions::State>,
    orbital: Res<crate::orbitals::OrbitalState>,
    mut traj: ResMut<TrajectoryState>,
    mut mol: ResMut<crate::molecule::Molecule>,
    mut settings: ResMut<crate::settings::MolSettings>,
    mut events: MessageWriter<crate::events::MoleculeChanged>,
    qc: Res<crate::qchem_panel::QcPanelState>,
    paths: Res<crate::qchem_interfaces::xtb_optimize::XtbPanelState>,
    builder: Res<crate::molecule_builder::builder_ui::ZMatrixBuilderState>,
    mut layout: ResMut<crate::ui_layout::UiLayout>,
    mut history: ResMut<crate::structure_history::StructureHistory>,
    mut regions: ResMut<crate::ui::UiPanelRegions>,
) {
    let Ok(ctx) = contexts.ctx_mut() else {
        return;
    };
    let state = &mut *state;
    if ctx
        .data_mut(|d| d.remove_temp::<bool>(egui::Id::new("open_reactivity")))
        .unwrap_or(false)
    {
        state.open = true;
        state.charge = qc.charge;
        state.mult[0] = qc.multiplicity;
    }
    if let Some(run) = &mut state.running {
        if let Some(result) = block_on(future::poll_once(&mut run.task)) {
            let cancelled = run.control.cancel.load(Ordering::Relaxed);
            state.running = None;
            let result = if cancelled {
                Err("Reactivity analysis cancelled.".into())
            } else {
                result
            };
            match result {
                Ok(r) => {
                    state.results = Some(Arc::new(r));
                    state.halos = true;
                    state.ready = None;
                }
                Err(e) => state.error = Some(e),
            }
        } else {
            ctx.request_repaint();
        }
    }
    if let Some(e) = state.error.take() {
        crate::runtime_warning::push(ctx, "Reactivity analysis", &e);
    }
    let mut load = false;
    if let Some(r)=state.results.clone() {
        let descriptor=state.descriptor;
        state.colors.update((Arc::as_ptr(&r) as usize,descriptor as usize),true,||r.frames.iter().filter_map(|f|descriptor.atoms(f)).flatten().collect());
    }
    if state.halos && state.descriptor.local() && state.current_index(&mol,&traj).is_some() {state.colors.overlay(ctx,"reactivity_legend",state.descriptor.label(),1.,state.descriptor.unit(),[state.positive,state.negative],90.);}
    if state.open {
        let mut open = true;
        let window=egui::Window::new("Reactivity / Fukui analysis").open(&mut open).resizable(true).default_size([650.,590.]).show(ctx,|ui|{
            crate::ui_style::modern(ui);
            ui.add_enabled_ui(!state.is_running(),|ui|{
                ui.horizontal_wrapped(|ui|{ui.selectable_value(&mut state.approx,true,"HOMO/LUMO approximation");ui.selectable_value(&mut state.approx,false,"Finite difference (N, N+1, N−1)");});
                ui.horizontal(|ui|{ui.selectable_value(&mut state.use_trajectory,false,"Current structure");ui.selectable_value(&mut state.use_trajectory,true,"Trajectory");});
                if state.approx{
                    ui.horizontal(|ui|{ui.label("Spin channel");ui.selectable_value(&mut state.spin,crate::orbitals::molden::Spin::Alpha,"Alpha / restricted");ui.selectable_value(&mut state.spin,crate::orbitals::molden::Spin::Beta,"Beta");});
                    egui::CollapsingHeader::new("Approximation details").default_open(false).show(ui,|ui|{ui.small("f+ ≈ |LUMO|²; f− ≈ |HOMO|². Frontier indices are selected separately in each frame. Degenerate orbitals and open shells require care.");});
                    crate::trajectory_wavefunctions::button(ui);
                }else{
                    ui.label(format!("Quantum Chemistry setup: {} · {}",qc.program.label(),qc.method.map(|m|m.label()).unwrap_or("Choose a method")));
                    ui.horizontal(|ui|{ui.label("Reference charge");ui.add(egui::DragValue::new(&mut state.charge).range(-100..=100));ui.label("Cores");ui.add(egui::DragValue::new(&mut state.cores).range(1..=4));});
                    ui.horizontal_wrapped(|ui|{for(i,label)in["N multiplicity","N+1 multiplicity","N−1 multiplicity"].iter().enumerate(){ui.label(*label);ui.add(egui::DragValue::new(&mut state.mult[i]).range(1..=100));}});
                    ui.small(format!("Charges: {}, {}, {}. All three states use identical geometry, method and basis.",state.charge,state.charge-1,state.charge+1));
                }
                if ui.button("Compute reactivity…").clicked_once(){let result=(||->Result<(),String>{
                    let control=crate::trajectory_wavefunctions::job::Control::default();let progress=Arc::new(AtomicUsize::new(0));let c=control.clone();let p=progress.clone();
                    let(task,total)=if state.approx{
                        let spin=state.spin;
                        let inputs:Vec<(TrajectoryFrame,data::Wave)>=if state.use_trajectory{let series=wave.series.as_ref().ok_or("Compute/load trajectory wavefunctions in Surface Tools first.")?;series.manifest.frames.iter().map(|f|(f.geometry.clone(),data::Wave::File(series.root.join(&f.molden)))).collect()}
                        else{let d=orbital.data.as_ref().ok_or("Load or compute a wavefunction in Surface Tools first.")?;vec![(TrajectoryFrame{atoms:d.atoms.clone(),pos:d.positions.clone(),metadata:Default::default()},data::Wave::Embedded(orbital.molden_text.clone().ok_or("No original Molden data")?))]};
                        let total=inputs.len();let task=AsyncComputeTaskPool::get().spawn(async move{let mut frames=vec![];for(i,(f,w))in inputs.into_iter().enumerate(){if c.cancel.load(Ordering::Relaxed){return Err("Reactivity analysis cancelled.".into());}frames.push(data::frontier(f,w,spin)?);p.store(i+1,Ordering::Relaxed);}Ok(data::Results{version:1,method:format!("HOMO/LUMO approximation · {spin:?} · Mulliken"),frames})});(task,total)
                    }else{
                        let input=if state.use_trajectory{traj.frames.clone()}else{vec![TrajectoryFrame{atoms:mol.atoms.clone(),pos:mol.pos.clone(),metadata:Default::default()}]};
                        let frames=crate::trajectory_wavefunctions::job::prepared(&input,1,input.len(),1)?;let total=frames.len()*3;let engine=job::engine(&qc,&paths,state.charge,state.mult[0],state.cores)?;
                        let Some(parent)=rfd::FileDialog::new().set_title("Choose a folder for a new reactivity calculation").pick_folder()else{return Ok(());};let root=crate::qchem_interfaces::xtb_optimize::create_xtb_run_dir(&parent)?;state.directory=Some(root.clone());let charge=state.charge;let mult=state.mult;
                        let task=AsyncComputeTaskPool::get().spawn(async move{job::finite(engine,frames,charge,mult,&root,c,p)});(task,total)
                    };state.running=Some(Running{task,control,progress,total});Ok(())})();if let Err(e)=result{state.error=Some(e);}}
                ui.horizontal(|ui|{
                    if ui.button("Load analysis…").clicked_once(){if let Some(path)=rfd::FileDialog::new().add_filter("Reactivity",&["beavyr-reactivity"]).pick_file(){match job::read(&path){Ok(r)=>{state.results=Some(Arc::new(r));state.ready=None;},Err(e)=>state.error=Some(e)}}}
                    if let Some(r)=&state.results{if ui.button("Save analysis…").clicked_once(){if let Some(path)=rfd::FileDialog::new().set_file_name("reactivity.beavyr-reactivity").save_file(){if let Err(e)=job::save(r,&path){state.error=Some(e);}}}}
                });
            });
            if let Some(run)=&state.running{ui.add(egui::ProgressBar::new(run.progress.load(Ordering::Relaxed)as f32/run.total.max(1)as f32).text(format!("{} / {}",run.progress.load(Ordering::Relaxed),run.total)));if ui.button("Cancel").clicked_once(){run.control.cancel();}}
            if let Some(result)=state.results.clone(){
                ui.separator();ui.label(&result.method);ui.horizontal(|ui|{if ui.button("Load as Trajectory").clicked_once(){load=true;}if ui.button("Plot atom / group…").clicked_once(){state.plot_open=true;}});
                egui::ComboBox::from_id_salt("reactivity_descriptor").selected_text(state.descriptor.label()).show_ui(ui,|ui|{for d in data::Descriptor::ALL{ui.selectable_value(&mut state.descriptor,d,d.label());}});
                ui.add_enabled_ui(state.descriptor.local(),|ui|{
                    ui.horizontal(|ui|{ui.checkbox(&mut state.halos,"Show atomic halos");ui.checkbox(&mut state.labels,"Show values");ui.checkbox(&mut state.surface,"Show surface");});
                    ui.horizontal(|ui|{ui.label("Positive");let mut c=state.positive.to_srgba().to_f32_array();if ui.color_edit_button_rgba_unmultiplied(&mut c).changed(){state.positive=Color::srgba(c[0],c[1],c[2],c[3]);}ui.label("Negative");let mut c=state.negative.to_srgba().to_f32_array();if ui.color_edit_button_rgba_unmultiplied(&mut c).changed(){state.negative=Color::srgba(c[0],c[1],c[2],c[3]);}});
                    crate::surface_controls::rows(ui,"reactivity_controls",|ui|{
                        crate::surface_controls::slider(ui,"Halo radius scale",egui::Slider::new(&mut state.scale,0.1..=5.));
                        crate::surface_controls::slider(ui,"Opacity",egui::Slider::new(&mut state.opacity,0.02..=1.));
                        if state.surface {
                            crate::surface_controls::slider(ui,"Isovalue",egui::Slider::new(&mut state.isovalue,0.00001..=0.5).logarithmic(true));
                            crate::surface_controls::slider(ui,"Surface resolution",egui::Slider::new(&mut state.resolution,1..=10));
                        }
                    });
                    state.colors.controls(ui,1.,state.descriptor.unit());
                });
                if let Some(i)=state.current_index(&mol,&traj){let f=&result.frames[i];let g=&f.global;ui.label(format!("I {:.4} eV · A {:.4} eV · μ {:.4} eV · η {:.4} eV",g.ionization,g.affinity,g.mu,g.eta));if let(Some(s),Some(w))=(g.softness,g.electrophilicity){ui.label(format!("S {s:.5} eV⁻¹ · ω {w:.5} eV"));}else{ui.label("Nonpositive gap: softness and electrophilicity unavailable.");}
                    if let Some(values)=state.descriptor.atoms(f){egui::ScrollArea::vertical().max_height(160.).show(ui,|ui|{egui::Grid::new("fukui_atoms").striped(true).show(ui,|ui|{for (j,v)in values.iter().enumerate(){ui.label(format!("{} {}",j+1,f.geometry.atoms[j]));ui.label(format!("{v:+.6} {}",state.descriptor.unit()));ui.end_row();}});});}
                }else{ui.label("Load the matching geometry to display its results.");}
                egui::CollapsingHeader::new("Method details").default_open(false).show(ui,|ui|{ui.small("Atomic values use Mulliken partitioning and may be negative. η=(I−A)/2, S=1/(2η). Local philicity is ω f. Saved analyses reference their wavefunction files; keep the calculation folder.");});
            }
            if let Some(dir)=&state.directory{if ui.button("Copy calculation folder").clicked_once(){ctx.copy_text(dir.display().to_string());}}
        });
        state.open = open;
        if let Some(w) = window {
            register(&mut regions, w.response.rect);
        }
    }
    if load {
        if let Some(result) = &state.results {
            let before = history.before(&mol, &traj);
            if crate::ts_generation::show_frames(
                result.frames.iter().map(|f| f.geometry.clone()).collect(),
                &mut mol,
                &mut traj,
            ) {
                history.replaced(before, &mol, &traj, "Reactivity trajectory");
                mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                events.write(crate::events::MoleculeChanged::parse_xyz(false));
                layout.open[crate::ui_layout::Tab::Trajectory.index()] = true;
            }
        }
    }
    plotting::window(
        ctx,
        state,
        &mut traj,
        &mut mol,
        &mut settings,
        &mut events,
        &builder,
        &qc,
        &mut layout,
        &mut regions,
    );
}

#[cfg(test)]
#[path = "reactivity_tests.rs"]
mod tests;
