//! Surface Tools controls, disk-backed trajectory CPA analysis and atom/group plots.
use crate::{
    trajectory::{TrajectoryFrame, TrajectoryState},
    ui_style::ResponseExt,
};
use bevy::{
    prelude::*,
    tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task},
};
use bevy_egui::{egui, EguiContexts, EguiPrimaryContextPass};
use std::sync::{
    atomic::{AtomicBool, AtomicUsize, Ordering},
    Arc,
};
#[path = "polarizability_data.rs"]
pub mod data;
#[path = "polarizability_plot.rs"]
pub(crate) mod plot;
#[path = "polarizability_render.rs"]
mod render;
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Results {
    pub frames: Vec<TrajectoryFrame>,
    pub values: Vec<data::Values>,
    pub method: String,
}
struct Running {
    task: Task<Result<Results, String>>,
    cancel: Arc<AtomicBool>,
    progress: Arc<AtomicUsize>,
    total: usize,
}
impl Drop for Running {
    fn drop(&mut self) {
        self.cancel.store(true, Ordering::Relaxed);
    }
}
#[derive(Resource)]
pub struct State {
    pub results: Option<Arc<Results>>,
    pub show: bool,
    pub labels: bool,
    pub angstrom: bool,
    pub scale: f32,
    pub opacity: f32,
    pub colors: crate::property_legend::Scale,
    open: bool,
    plot_open: bool,
    task: Option<Running>,
    group: String,
    axis: data::Axis,
    atoms: [usize; 4],
    dt: f64,
    origin: f64,
    error: Option<String>,
}
impl Default for State {
    fn default() -> Self {
        Self {
            results: None,
            show: false,
            labels: false,
            angstrom: true,
            scale: 1.,
            opacity: 0.2,
            colors: Default::default(),
            open: false,
            plot_open: false,
            task: None,
            group: "1".into(),
            axis: data::Axis::Frame,
            atoms: [1, 2, 3, 4],
            dt: 1.,
            origin: 0.,
            error: None,
        }
    }
}
impl State {
    pub fn factor(&self) -> f64 {
        if self.angstrom {
            data::BOHR3_TO_ANGSTROM3
        } else {
            1.
        }
    }
    pub fn unit(&self) -> &'static str {
        if self.angstrom {
            "Å³"
        } else {
            "bohr³"
        }
    }
    pub fn current(
        &self,
        mol: &crate::molecule::Molecule,
        traj: &TrajectoryState,
    ) -> Option<&data::Values> {
        let r = self.results.as_ref()?;
        let index = if r.frames.len() == 1 {
            0
        } else {
            traj.current_frame
        };
        let frame = r.frames.get(index)?;
        if frame.atoms != mol.atoms
            || frame.pos.len() != mol.pos.len()
            || frame
                .pos
                .iter()
                .zip(&mol.pos)
                .any(|(a, b)| a.distance(*b) > 1e-5)
        {
            return None;
        }
        r.values.get(index)
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
    if ui.button("Atomic polarizabilities…").clicked_once() {
        ui.ctx()
            .data_mut(|d| d.insert_temp(egui::Id::new("open_polarizability"), true));
    }
}
fn csv(
    result: &Results,
    points: &[Option<(f64, f64)>],
    group: &[usize],
    axis: data::Axis,
    unit: &str,
) -> String {
    let mut text = format!(
        "# Stoehr CPA / TS; {}\n# atoms: {}\nframe,{},sum_alpha_{}",
        result.method,
        group
            .iter()
            .map(|i| (i + 1).to_string())
            .collect::<Vec<_>>()
            .join(" "),
        axis.label(),
        unit
    );
    for i in group {
        text.push_str(&format!(",atom_{}_alpha_{}", i + 1, unit));
    }
    text.push('\n');
    let factor = if unit == "Å³" {
        data::BOHR3_TO_ANGSTROM3
    } else {
        1.
    };
    for (i, p) in points.iter().enumerate() {
        if let Some((x, y)) = p {
            text.push_str(&format!("{},{x:.10},{y:.10}", i + 1));
            for &j in group {
                text.push_str(&format!(",{:.10}", result.values[i].alpha[j] * factor));
            }
            text.push('\n');
        }
    }
    text
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
    builder: Res<crate::molecule_builder::builder_ui::ZMatrixBuilderState>,
    mut layout: ResMut<crate::ui_layout::UiLayout>,
    qc: Res<crate::qchem_panel::QcPanelState>,
    mut regions: ResMut<crate::ui::UiPanelRegions>,
) {
    let Ok(ctx) = contexts.ctx_mut() else {
        return;
    };
    let state = &mut *state;
    if ctx
        .data_mut(|d| d.remove_temp::<bool>(egui::Id::new("open_polarizability")))
        .unwrap_or(false)
    {
        state.open = true;
    }
    if let Some(run) = &mut state.task {
        if let Some(result) = block_on(future::poll_once(&mut run.task)) {
            let cancelled = run.cancel.load(Ordering::Relaxed);
            state.task = None;
            let result = if cancelled {
                Err("Polarizability calculation cancelled.".into())
            } else {
                result
            };
            match result {
                Ok(r) => {
                    state.results = Some(Arc::new(r));
                    state.show = true;
                }
                Err(e) => state.error = Some(e),
            }
        } else {
            ctx.request_repaint();
        }
    }
    if let Some(error) = state.error.take() {
        crate::runtime_warning::push(ctx, "Atomic polarizabilities", &error);
    }
    if let Some(r)=state.results.clone() {
        state.colors.update((Arc::as_ptr(&r) as usize,0),false,||r.values.iter().flat_map(|v|v.alpha.iter().copied()).collect());
    }
    if state.show && state.current(&mol,&traj).is_some() { state.colors.overlay(ctx,"polar_legend","Atomic α",state.factor(),state.unit(),[Color::srgb(0.25,0.75,1.);2],0.); }
    if state.open {
        let mut open = true;
        let window=egui::Window::new("Atomic polarizabilities").open(&mut open).resizable(true).default_size([530.,420.]).show(ctx,|ui|{
            crate::ui_style::modern(ui);
            ui.horizontal_wrapped(|ui|{
                if ui.add_enabled(state.task.is_none()&&orbital.data.is_some(),egui::Button::new("Compute current wavefunction")).clicked_once(){
                    let d=orbital.data.as_ref().unwrap();match data::compute(d){Ok(v)=>{state.results=Some(Arc::new(Results{frames:vec![TrajectoryFrame{atoms:d.atoms.clone(),pos:d.positions.clone(),metadata:Default::default()}],values:vec![v],method:orbital.provenance.clone().unwrap_or_else(||"Imported Molden".into())}));state.show=true;},Err(e)=>state.error=Some(e)}
                }
                if ui.add_enabled(state.task.is_none()&&wave.series.is_some(),egui::Button::new("Compute trajectory")).clicked_once(){
                    let series=wave.series.clone().unwrap();let total=series.manifest.frames.len();let cancel=Arc::new(AtomicBool::new(false));let progress=Arc::new(AtomicUsize::new(0));let c=cancel.clone();let p=progress.clone();
                    let task=AsyncComputeTaskPool::get().spawn(async move{let mut values=vec![];for i in 0..total{if c.load(Ordering::Relaxed){return Err("Polarizability calculation cancelled.".into());}let (data,_)=crate::trajectory_wavefunctions::job::load(&series,i)?;values.push(data::compute(&data)?);p.store(i+1,Ordering::Relaxed);}let result=Results{frames:series.frames(),values,method:series.manifest.method.clone()};Ok(result)});
                    state.task=Some(Running{task,cancel,progress,total});
                }
            });
            crate::trajectory_wavefunctions::button(ui);
            if let Some(run)=&state.task{ui.add(egui::ProgressBar::new(run.progress.load(Ordering::Relaxed)as f32/run.total as f32));if ui.button("Cancel").clicked_once(){run.cancel.store(true,Ordering::Relaxed);}}
            ui.horizontal(|ui|{ui.checkbox(&mut state.show,"Show atomic halos");ui.checkbox(&mut state.labels,"Values");});
            ui.horizontal(|ui|{ui.selectable_value(&mut state.angstrom,true,"Å³");ui.selectable_value(&mut state.angstrom,false,"bohr³");});
            crate::surface_controls::rows(ui,"polar_controls",|ui|{
                crate::surface_controls::slider(ui,"Radius scale",egui::Slider::new(&mut state.scale,0.1..=3.));
                crate::surface_controls::slider(ui,"Opacity",egui::Slider::new(&mut state.opacity,0.02..=0.8));
            });
            let factor=state.factor();let unit=state.unit();state.colors.controls(ui,factor,unit);
            if let Some(result)=state.results.clone(){
                ui.label(format!("{} · {} geometries",result.method,result.frames.len()));
                if ui.button("Plot atom / group…").clicked_once(){state.plot_open=true;}
                if let Some(values)=state.current(&mol,&traj){
                    let factor=state.factor();let unit=state.unit();
                    egui::ScrollArea::vertical().max_height(180.).show(ui,|ui|{egui::Grid::new("polar_table").striped(true).show(ui,|ui|{ui.label("Atom");ui.label(format!("α ({unit})"));ui.label("CPA population");ui.end_row();for(i,v)in values.alpha.iter().enumerate(){ui.label(format!("{} {}",i+1,mol.atoms[i]));ui.label(format!("{:.5}",v*factor));ui.label(format!("{:.5}",values.population[i]));ui.end_row();}});});
                }else{ui.label("These results do not match the current geometry.");}
            }
            egui::CollapsingHeader::new("Method details").default_open(false).show(ui,|ui|{ui.small("Stöhr CPA / TS: αA = αfree × hA/ZA. The sum is an atomic descriptor, not molecular response polarizability. Values depend on the basis; xTB/ECP wavefunctions omit core electrons.");});
        });
        state.open = open;
        if let Some(w) = window {
            register(&mut regions, w.response.rect);
        }
    }
    if state.plot_open {
        if let Some(result) = state.results.clone() {
            let mut open = true;
            let window=egui::Window::new("Atomic polarizability plot").open(&mut open).resizable(true).default_size([740.,520.]).show(ctx,|ui|{
            ui.horizontal(|ui|{ui.label("Atoms");ui.text_edit_singleline(&mut state.group);if ui.button("Use selected atoms").clicked_once(){state.group=builder.selection.atoms.iter().map(|i|(i+1).to_string()).collect::<Vec<_>>().join(",");}if ui.button("All").clicked_once(){state.group=format!("1-{}",result.frames[0].atoms.len());}});
            ui.horizontal_wrapped(|ui|{for axis in[data::Axis::Frame,data::Axis::Time,data::Axis::Distance,data::Axis::Angle,data::Axis::Dihedral]{ui.selectable_value(&mut state.axis,axis,axis.label());}ui.add_space(20.);ui.selectable_value(&mut state.angstrom,true,"Å³");ui.selectable_value(&mut state.angstrom,false,"bohr³");});
            if state.axis==data::Axis::Time{ui.horizontal(|ui|{ui.label("Time per sampled frame (fs)");ui.add(egui::DragValue::new(&mut state.dt).range(0.000001..=1e12));ui.label("First time (fs)");ui.add(egui::DragValue::new(&mut state.origin));});}
            if state.axis.count()>0{ui.horizontal(|ui|{ui.label("Coordinate atoms");for a in state.atoms.iter_mut().take(state.axis.count()){ui.add(egui::DragValue::new(a).range(1..=result.frames[0].atoms.len()));}if ui.button("Use scan coordinate").clicked_once(){use crate::qchem_interfaces::job::scan::Coordinate;let scan=&qc.config.scan;state.axis=match scan.coordinate{Coordinate::Bond=>data::Axis::Distance,Coordinate::Angle=>data::Axis::Angle,Coordinate::Dihedral=>data::Axis::Dihedral};match crate::qchem_interfaces::job::scan::atom_indices(&scan.atoms,scan.coordinate.count(),&result.frames[0].atoms){Ok(indices)=>{for (a,i) in state.atoms.iter_mut().zip(indices){*a=i+1;}},Err(e)=>state.error=Some(e)}}});}
            let group=match data::atom_list(&state.group,result.frames[0].atoms.len()){Ok(g)=>g,Err(e)=>{ui.label(e);return;}};
            let points:Vec<_>=result.frames.iter().enumerate().map(|(i,f)|{let x=state.axis.value(f,i,&state.atoms,state.dt,state.origin)?;let y=group.iter().map(|&j|result.values[i].alpha[j]).sum::<f64>()*state.factor();Some((x,y))}).collect();
            ui.horizontal(|ui|{let text=||csv(&result,&points,&group,state.axis,state.unit());if ui.button("Copy data").clicked_once(){ctx.copy_text(text());}if ui.button("Save CSV…").clicked_once(){if let Some(path)=rfd::FileDialog::new().set_file_name("atomic-polarizability.csv").save_file(){if let Err(e)=std::fs::write(path,text()){state.error=Some(e.to_string());}}}});
            let matches=result.frames.len()==traj.frames.len()&&result.frames.iter().zip(&traj.frames).all(|(a,b)|crate::trajectory_wavefunctions::job::geometry_matches(a,b));
            if !matches{ui.label("Load the matching wavefunction trajectory to select frames from this plot.");}
            if let Some(i)=plot::draw(ui,&points,if matches{Some(traj.current_frame)}else{None},state.axis.label(),state.unit()){
                if matches{traj.playing=false;traj.current_frame=i;crate::trajectory::apply_current_frame(&mut traj,&mut mol,&mut settings,&mut events,false,None);layout.open[crate::ui_layout::Tab::Trajectory.index()]=true;}
            }
        });
            state.plot_open = open;
            if let Some(w) = window {
                register(&mut regions, w.response.rect);
            }
        }
    }
}
fn register(regions: &mut crate::ui::UiPanelRegions, r: egui::Rect) {
    regions
        .windows
        .push(Rect::new(r.min.x, r.min.y, r.max.x, r.max.y));
}

impl Results {
    pub fn validate(&self) -> Result<(), String> {
        if self.frames.is_empty() || self.frames.len() != self.values.len() {
            return Err("Invalid polarizability frame count".into());
        }
        for (f, v) in self.frames.iter().zip(&self.values) {
            let n = f.atoms.len();
            if n == 0
                || f.pos.len() != n
                || f.pos.iter().any(|p| !p.is_finite())
                || f.atoms != self.frames[0].atoms
                || v.alpha.len() != n
                || v.population.len() != n
                || v.c6.len() != n
                || v.alpha
                    .iter()
                    .chain(&v.population)
                    .chain(&v.c6)
                    .any(|x| !x.is_finite() || *x < 0.)
            {
                return Err("Invalid atomic polarizability data".into());
            }
        }
        Ok(())
    }
}
