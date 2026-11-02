//! Surface Tools entry point for disk-backed wavefunctions along a trajectory.
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
    sync::{atomic::Ordering, Arc},
};
#[path = "trajectory_wavefunction_job.rs"]
pub mod job;
pub struct Plugin;
impl bevy::prelude::Plugin for Plugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<State>()
            .add_plugins(crate::orbitals::trajectory_display::Plugin)
            .add_systems(EguiPrimaryContextPass, ui.after(crate::ui::ui_panel));
    }
}
struct Running {
    task: Task<Result<job::Series, String>>,
    control: job::Control,
    total: usize,
    source: u64,
    started: std::time::Instant,
}
impl Drop for Running {
    fn drop(&mut self) {
        self.control.cancel();
    }
}
#[derive(Resource)]
pub struct State {
    pub series: Option<Arc<job::Series>>,
    pub follow: bool,
    pub ready: Option<usize>,
    pub tracking: crate::orbital_tracking::State,
    pub error: Option<String>,
    open: bool,
    running: Option<Running>,
    engine: u8,
    charge: i32,
    multiplicity: i32,
    cores: usize,
    first: usize,
    last: usize,
    stride: usize,
    use_trajectory: bool,
    path: String,
    message: String,
    pub directory: Option<PathBuf>,
}
impl Default for State {
    fn default() -> Self {
        Self {
            series: None,
            follow: false,
            ready: None,
            tracking: Default::default(),
            error: None,
            open: false,
            running: None,
            engine: 0,
            charge: 0,
            multiplicity: 1,
            cores: 4,
            first: 1,
            last: 1,
            stride: 1,
            use_trajectory: true,
            path: String::new(),
            message: String::new(),
            directory: None,
        }
    }
}
impl State {
    pub fn is_running(&self) -> bool {
        self.running.is_some()
    }
    pub fn hold_playback(&self, traj: &TrajectoryState) -> bool {
        self.follow
            && self.series.as_ref().is_some_and(|s| s.matches(traj))
            && self.ready != Some(traj.current_frame)
    }
}
pub fn button(ui: &mut egui::Ui) {
    if ui.button("Trajectory wavefunctions…").clicked_once() {
        ui.ctx()
            .data_mut(|d| d.insert_temp(egui::Id::new("open_trajectory_wavefunctions"), true));
    }
}
fn current_frames(
    state: &State,
    traj: &TrajectoryState,
    mol: &crate::molecule::Molecule,
) -> Vec<TrajectoryFrame> {
    if state.use_trajectory {
        traj.frames.clone()
    } else {
        vec![TrajectoryFrame {
            atoms: mol.atoms.clone(),
            pos: mol.pos.clone(),
            metadata: Default::default(),
        }]
    }
}
#[allow(clippy::too_many_arguments)]
fn ui(
    mut contexts: EguiContexts,
    mut state: ResMut<State>,
    mut traj: ResMut<TrajectoryState>,
    mut mol: ResMut<crate::molecule::Molecule>,
    settings: Res<crate::settings::MolSettings>,
    mut events: MessageWriter<crate::events::MoleculeChanged>,
    qc: Res<crate::qchem_panel::QcPanelState>,
    orbital: Res<crate::orbitals::OrbitalState>,
    paths: Res<crate::qchem_interfaces::xtb_optimize::XtbPanelState>,
    mut layout: ResMut<crate::ui_layout::UiLayout>,
    mut history: ResMut<crate::structure_history::StructureHistory>,
    mut regions: ResMut<crate::ui::UiPanelRegions>,
) {
    let Ok(ctx) = contexts.ctx_mut() else {
        return;
    };
    let current_series=state.series.clone();
    let tracking_was_enabled=state.tracking.enabled;
    state.tracking.poll(current_series.as_deref());
    if state.tracking.enabled && !tracking_was_enabled {state.follow=true;state.ready=None;}
    let mut load = false;
    if ctx
        .data_mut(|d| d.remove_temp::<bool>(egui::Id::new("open_trajectory_wavefunctions")))
        .unwrap_or(false)
    {
        state.open = true;
        state.use_trajectory = !traj.frames.is_empty();
        state.last = traj.frames.len().max(1);
        state.charge = qc.charge;
        state.multiplicity = qc.multiplicity;
        if state.path.is_empty() {
            state.path = paths
                .paths
                .iter()
                .find(|(p, _)| *p == crate::qchem_interfaces::program::QcProgram::Xtb)
                .map(|(_, s)| s.clone())
                .unwrap_or_default();
        }
    }
    if let Some(run) = &mut state.running {
        if let Some(result) = block_on(future::poll_once(&mut run.task)) {
            let original = run.source;
            state.running = None;
            match result {
                Ok(series) => {
                    let frames = current_frames(&state, &traj, &mol);
                    load = crate::trajectory_data::fingerprint(&frames) == original;
                    state.message = if load {
                        "Wavefunctions completed.".into()
                    } else {
                        "Wavefunctions completed. The viewer changed during calculation; use Load as Trajectory to display the result.".into()
                    };
                    state.series = Some(Arc::new(series));
                    state.ready = None;
                }
                Err(e) => state.error = Some(e),
            }
        } else {
            ctx.request_repaint();
        }
    }
    if let Some(error) = state.error.take() {
        state.message = error.clone();
        crate::runtime_warning::push(ctx, "Trajectory wavefunctions", &error);
    }
    if state.open {
        let mut open = true;
        let window=egui::Window::new("Trajectory wavefunctions").open(&mut open).resizable(true).default_size([580.,470.]).show(ctx,|ui|{
            crate::ui_style::modern(ui);
            ui.add_enabled_ui(!state.is_running(),|ui|{
                ui.horizontal_wrapped(|ui|{
                    if ui.selectable_value(&mut state.use_trajectory,false,"Current structure").changed(){state.first=1;state.last=1;}
                    if ui.add_enabled(!traj.frames.is_empty(),egui::Button::selectable(state.use_trajectory,"Loaded trajectory")).clicked_once(){state.use_trajectory=true;state.first=1;state.last=traj.frames.len();}
                });
                ui.horizontal_wrapped(|ui|{ui.selectable_value(&mut state.engine,0,"GFN2-xTB");ui.selectable_value(&mut state.engine,1,"GFN1-xTB");ui.selectable_value(&mut state.engine,2,"Quantum Chemistry setup");});
                if state.engine<2{
                    ui.horizontal(|ui|{ui.label("xTB executable");ui.text_edit_singleline(&mut state.path);if ui.button("Browse…").clicked_once(){if let Some(p)=crate::recent_dir::pick_file(crate::recent_dir::open()){state.path=p.display().to_string();}}});
                }else{ui.label(format!("{} · {} · {}",qc.program.label(),qc.method.map(|m|m.label()).unwrap_or("Choose a method in Quantum Chemistry"),qc.config.basis));}
                ui.horizontal_wrapped(|ui|{
                    ui.label("Charge");ui.add(egui::DragValue::new(&mut state.charge).range(-100..=100));
                    ui.label("Multiplicity");ui.add(egui::DragValue::new(&mut state.multiplicity).range(1..=100));
                    ui.label("Cores");ui.add(egui::DragValue::new(&mut state.cores).range(1..=4));
                });
                if state.use_trajectory{ui.horizontal(|ui|{ui.label("Frames");ui.add(egui::DragValue::new(&mut state.first).range(1..=traj.frames.len().max(1)));ui.label("to");ui.add(egui::DragValue::new(&mut state.last).range(1..=traj.frames.len().max(1)));ui.label("Stride");ui.add(egui::DragValue::new(&mut state.stride).range(1..=100000));});}
                if ui.button("Compute wavefunctions…").clicked_once(){
                    let result=(||->Result<(),String>{
                        use crate::qchem_interfaces::{program::QcProgram,xtb_optimize::resolve_program_executable};
                        let frames=current_frames(&state,&traj,&mol);let source=crate::trajectory_data::fingerprint(&frames);
                        let selected=job::prepared(&frames,state.first,state.last,state.stride)?;
                        let engine=if state.engine<2{
                            job::Engine::Xtb{binary:resolve_program_executable(QcProgram::Xtb,&state.path)?,gfn:if state.engine==0{2}else{1},charge:state.charge,multiplicity:state.multiplicity,cores:state.cores}
                        }else{
                            if !matches!(qc.program,QcProgram::Orca|QcProgram::PySCF|QcProgram::Psi4|QcProgram::Psi4Py){return Err("This route supports ORCA, PySCF and Psi4. Select one in Quantum Chemistry, or choose GFN2-xTB here.".into());}
                            let path=paths.paths.iter().find(|(p,_)|*p==qc.program).map(|(_,s)|s.as_str()).unwrap_or("");
                            let binary=if qc.program.needs_binary_path(){resolve_program_executable(qc.program,path)?}else{PathBuf::new()};
                            let mut config=qc.config.clone();config.nproc=state.cores as u32;
                            job::Engine::QuantumChemistry{binary,program:qc.program,method:qc.method,config,charge:state.charge,multiplicity:state.multiplicity}
                        };
                        let Some(parent)=rfd::FileDialog::new().set_title("Choose a folder for a new wavefunction trajectory").pick_folder()else{return Ok(());};
                        let root=crate::qchem_interfaces::xtb_optimize::create_xtb_run_dir(&parent)?;state.directory=Some(root.clone());
                        let control=job::Control::default();let worker_control=control.clone();let total=selected.len();
                        let task=AsyncComputeTaskPool::get().spawn(async move{job::calculate(engine,selected,root,worker_control)});
                        state.running=Some(Running{task,control,total,source,started:std::time::Instant::now()});state.message.clear();Ok(())
                    })();if let Err(e)=result{state.error=Some(e);}
                }
                if ui.button("Resume calculation…").on_hover_text("Uses the saved method, charge, spin and frame plan. Verifies completed files; retries unfinished frames in new folders.").clicked_once(){
                    if let Some(root)=rfd::FileDialog::new().set_title("Choose the calculation folder containing wavefunction-plan.ron").pick_folder(){
                        match job::resume::read_plan(&root){Ok(plan)=>{
                            let total=plan.frames.len();let control=job::Control::default();let worker=control.clone();state.directory=Some(root.clone());
                            let task=AsyncComputeTaskPool::get().spawn(async move{job::resume::run(root,worker)});
                            state.running=Some(Running{task,control,total,source:0,started:std::time::Instant::now()});state.message.clear();
                        },Err(e)=>state.error=Some(e)}
                    }
                }
                if ui.button("Load wavefunction trajectory…").clicked_once(){
                    if let Some(path)=crate::recent_dir::pick_file(crate::recent_dir::open().add_filter("Wavefunction trajectory",&["beavyr-wfn"])){
                        match job::read(&path){Ok(series)=>{state.series=Some(Arc::new(series));state.ready=None;load=true;},Err(e)=>state.error=Some(e)}
                    }
                }
            });
            if let Some(run)=&state.running{
                let done=run.control.progress.load(Ordering::Relaxed);let reused=run.control.reused.load(Ordering::Relaxed);
                let elapsed=run.started.elapsed().as_secs_f64();let new_done=done.saturating_sub(reused);
                ui.add(egui::ProgressBar::new(done as f32/run.total as f32).text(format!("{done} / {} completed · {} remaining",run.total,run.total.saturating_sub(done))));
                if let Ok(phase)=run.control.phase.lock(){ui.label(phase.as_str());}
                let eta=if new_done>0{job::resume::duration(elapsed/new_done as f64*run.total.saturating_sub(done)as f64)}else{"Estimating…".into()};
                ui.label(format!("Elapsed {} · Remaining ≈ {eta}",job::resume::duration(elapsed)));
                if reused>0{ui.small(format!("{reused} completed frames verified and reused"));}
                if ui.button("Cancel").clicked_once(){run.control.cancel();}
            }
            if let Some(series)=state.series.clone(){
                ui.separator();ui.label(format!("{} · {} frames",series.manifest.method,series.manifest.frames.len()));
                if !series.manifest.complete{ui.label("Partial calculation: only completed frames are available.");}
                ui.horizontal(|ui|{
                    if ui.button("Load as Trajectory").clicked_once(){load=true;}
                    if ui.checkbox(&mut state.follow,"Follow trajectory in Surface Tools").changed(){state.ready=None;if !state.follow {state.tracking.enabled=false;}}
                    if ui.button("Copy manifest path").clicked_once(){ctx.copy_text(series.root.join(job::MANIFEST).display().to_string());}
                });
                if series.manifest.frames.iter().any(|f|!f.spin_resolved){ui.small("Some frames lack separate alpha/beta orbitals; spin density is unavailable for those files.");}
                ui.add_enabled_ui(series.matches(&traj),|ui|{state.tracking.ui(ui,&series,&orbital,traj.current_frame);});
                if state.tracking.enabled {state.follow=true;}
            }
            if let Some(dir)=&state.directory{ui.horizontal(|ui|{ui.label("Calculation files");if ui.button("Copy folder path").clicked_once(){ctx.copy_text(dir.display().to_string());}});}
            if !state.message.is_empty(){ui.label(&state.message);}
        });
        state.open = open;
        if let Some(w) = window {
            let r = w.response.rect;
            regions
                .windows
                .push(Rect::new(r.min.x, r.min.y, r.max.x, r.max.y));
        }
    }
    if load {
        if let Some(series) = state.series.clone() {
            let before = history.before(&mol, &traj);
            if crate::ts_generation::show_frames(series.frames(), &mut mol, &mut traj) {
                history.replaced(before, &mol, &traj, "Wavefunction trajectory");
                mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                events.write(crate::events::MoleculeChanged::parse_xyz(false));
                layout.open[crate::ui_layout::Tab::Trajectory.index()] = true;
                state.follow = true;
                state.ready = None;
            }
        }
    }
}
