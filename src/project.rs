//! Versioned, self-contained projects. Parse and validate before replacing a session.
use crate::ui_style::ResponseExt;
use crate::{
    camera::OrbitCamera,
    molecule::Molecule,
    orbitals::{OrbitalState, SurfaceStyle},
    qchem_interfaces::xtb_freq::{FrequencyArchive, XtbFreqPanelState, XtbFrequencyTask},
    qchem_panel::QcPanelState,
    settings::MolSettings,
    trajectory::TrajectoryState,
    ui::XyzBuffer,
    ui_layout::UiLayout,
    uvvis::UvVisState,
};
use anyhow::{bail, Context, Result};
use bevy::prelude::*;
use bevy_egui::egui;
use serde::{Deserialize, Serialize};
use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};

const VERSION: u32 = 1;
static SERIAL: AtomicU64 = AtomicU64::new(0);

#[derive(Clone, Serialize, Deserialize)]
struct Workspace {
    memory: egui::Memory,
    editor_docked: bool,
    qc_presentation: (bool, bool),
    input_editor: (
        bool,
        Option<crate::qchem_interfaces::program::QcProgram>,
        String,
    ),
}

impl Workspace {
    fn capture(ctx: &egui::Context) -> Self {
        Self {
            memory: ctx.memory(Clone::clone),
            editor_docked: crate::molecule_builder::builder_ui::editor_docked(ctx),
            qc_presentation: crate::qchem_panel::ui::project_presentation(ctx),
            input_editor: crate::qchem_panel::ui::project_input_state(ctx),
        }
    }
    fn apply(self, ctx: &egui::Context) {
        // This runs inside an egui pass. Keep its live input/focus bookkeeping;
        // deserialization deliberately omits that state.
        ctx.memory_mut(|m| {
            *m.areas_mut() = self.memory.areas().clone();
            m.options = self.memory.options;
            m.data = self.memory.data;
            m.to_global = self.memory.to_global;
        });
        crate::molecule_builder::builder_ui::restore_editor_docking(ctx, self.editor_docked);
        crate::qchem_panel::ui::restore_project_presentation(ctx, self.qc_presentation);
        crate::qchem_panel::ui::restore_project_input_state(ctx, self.input_editor);
    }
}

enum Request {
    Save(PathBuf, Workspace),
    Open(PathBuf, Workspace),
    Restore(Workspace),
}

#[derive(Resource, Default)]
pub struct ProjectState {
    request: Option<Request>,
    restore_ui: Option<Workspace>,
    previous: Option<Snapshot>,
    pub path: Option<PathBuf>,
    message: Option<String>,
    error: bool,
}

pub fn restore_ui(ctx: &egui::Context, state: &mut ProjectState) {
    if let Some(workspace) = state.restore_ui.take() {
        workspace.apply(ctx);
    }
}

/// Compact entry point next to the Structure file controls.
pub fn menu(ui: &mut egui::Ui, state: &mut ProjectState) {
    ui.menu_button("Project", |ui| {
        if ui.button("Open project…").clicked_once() {
            if let Some(path) = crate::recent_dir::pick_file(
                crate::recent_dir::open().add_filter("Beavyr project", &["beavyr"]),
            ) {
                state.request = Some(Request::Open(path, Workspace::capture(ui.ctx())));
            }
            ui.close();
        }
        if ui.button("Save project…").clicked_once() {
            let name = state
                .path
                .as_ref()
                .and_then(|p| p.file_name())
                .and_then(|s| s.to_str())
                .unwrap_or("project.beavyr");
            if let Some(path) = crate::recent_dir::save_file(
                crate::recent_dir::save(name).add_filter("Beavyr project", &["beavyr"]),
            ) {
                state.request = Some(Request::Save(path, Workspace::capture(ui.ctx())));
            }
            ui.close();
        }
        if ui
            .add_enabled(
                state.previous.is_some(),
                egui::Button::new("Restore previous session"),
            )
            .on_hover_text("Undo the last project open, including its data and settings")
            .clicked_once()
        {
            state.request = Some(Request::Restore(Workspace::capture(ui.ctx())));
            ui.close();
        }
        if let Some(path) = &state.path {
            ui.weak(path.display().to_string());
        }
        if let Some(message) = &state.message {
            if state.error {
                ui.colored_label(egui::Color32::LIGHT_RED, message);
            } else {
                ui.weak(message);
            }
        }
    });
    if state.error {
        ui.colored_label(egui::Color32::LIGHT_RED, "Project error")
            .on_hover_text(state.message.as_deref().unwrap_or_default());
    }
}

#[derive(Clone, Serialize, Deserialize)]
struct Orbitals {
    molden: String,
    file_name: Option<String>,
    provenance: Option<String>,
    quantity: crate::orbitals::density::Quantity,
    spin: crate::orbitals::molden::Spin,
    selected: Option<usize>,
    browser: crate::orbitals::browser::Browser,
    isovalue: f32,
    resolution: u32,
    style: SurfaceStyle,
    mesh_density: u32,
    mesh_width: f32,
    visible: bool,
    opacity: f32,
    positive: Color,
    negative: Color,
}

impl Orbitals {
    fn capture(s: &OrbitalState) -> Result<Option<Self>> {
        if s.data.is_none() {
            return Ok(None);
        }
        Ok(Some(Self {
            molden: s.molden_text.clone().context(
                "The loaded orbitals have no source text; reload their Molden file before saving",
            )?,
            file_name: s.file_name.clone(),
            provenance: s.provenance.clone(),
            quantity: s.quantity,
            spin: s.spin,
            selected: s.selected,
            browser: s.browser.clone(),
            isovalue: s.isovalue,
            resolution: s.resolution,
            style: s.style,
            mesh_density: s.mesh_density,
            mesh_width: s.mesh_width,
            visible: s.show_surface,
            opacity: s.opacity,
            positive: s.positive_color,
            negative: s.negative_color,
        }))
    }
    fn prepare(&self) -> Result<OrbitalState> {
        let data =
            crate::orbitals::molden::parse_molden(&self.molden).map_err(anyhow::Error::msg)?;
        if self
            .selected
            .is_some_and(|i| i >= data.orbitals(self.spin).len())
            || !self.isovalue.is_finite()
            || self.isovalue <= 0.0
            || !(crate::orbitals::MIN_RESOLUTION..=crate::orbitals::MAX_RESOLUTION)
                .contains(&self.resolution)
            || !self.opacity.is_finite()
            || !(0.0..=1.0).contains(&self.opacity)
            || !self.mesh_width.is_finite()
            || self.mesh_width <= 0.0
        {
            bail!("Invalid saved orbital selection or surface settings");
        }
        let mut s = OrbitalState::default();
        s.adopt(
            data,
            self.file_name.clone(),
            self.provenance.clone(),
            Some(self.molden.clone()),
        );
        s.quantity = self.quantity;
        s.spin = self.spin;
        s.selected = self.selected;
        s.browser = self.browser.clone();
        s.browser.reveal = true;
        s.isovalue = self.isovalue;
        s.resolution = self.resolution;
        s.style = self.style;
        s.mesh_density = self.mesh_density.clamp(
            crate::orbitals::MIN_MESH_DENSITY,
            crate::orbitals::MAX_MESH_DENSITY,
        );
        s.mesh_width = self.mesh_width;
        s.show_surface = self.visible;
        s.opacity = self.opacity;
        s.positive_color = self.positive;
        s.negative_color = self.negative;
        Ok(s)
    }
}

#[derive(Clone, Serialize, Deserialize)]
struct Snapshot {
    #[serde(default)]
    topology: Option<crate::molecule_builder::topology::Topology>,
    version: u32,
    application: String,
    atoms: Vec<String>,
    positions: Vec<Vec3>,
    xyz_draft: String,
    settings: MolSettings,
    camera: OrbitCamera,
    trajectory: TrajectoryState,
    calculation: QcPanelState,
    optimization: crate::qchem_interfaces::xtb_optimize::XtbPanelState,
    frequencies: FrequencyArchive,
    frequency_settings: XtbFreqPanelState,
    spectrum: UvVisState,
    orbitals: Option<Orbitals>,
    linked_orbitals: bool,
    layout: UiLayout,
    workspace: Workspace,
}

impl Snapshot {
    fn capture(world: &World, workspace: Workspace) -> Result<Self> {
        let mol = world.resource::<Molecule>();
        let orbitals = world.resource::<OrbitalState>();
        let spectrum = world.resource::<UvVisState>();
        Ok(Self {
            version: VERSION,
            application: env!("CARGO_PKG_VERSION").into(),
            topology: mol.topology.clone(),
            atoms: mol.atoms.clone(),
            positions: mol.pos.clone(),
            xyz_draft: world.resource::<XyzBuffer>().text.clone(),
            settings: world.resource::<MolSettings>().clone(),
            camera: world.resource::<OrbitCamera>().clone(),
            trajectory: world.resource::<TrajectoryState>().clone(),
            calculation: world.resource::<QcPanelState>().clone(),
            optimization: world
                .resource::<crate::qchem_interfaces::xtb_optimize::XtbPanelState>()
                .clone(),
            frequencies: world.resource::<XtbFrequencyTask>().archive(),
            frequency_settings: world.resource::<XtbFreqPanelState>().clone(),
            spectrum: spectrum.clone(),
            orbitals: Orbitals::capture(orbitals)?,
            linked_orbitals: spectrum.orbital_link.matches(orbitals),
            layout: world.resource::<UiLayout>().clone(),
            workspace,
        })
    }

    fn validate(&self) -> Result<OrbitalState> {
        if self.version != VERSION {
            bail!(
                "Unsupported project version {} (this build reads {})",
                self.version,
                VERSION
            );
        }
        validate_geometry(&self.atoms, &self.positions)?;
        if self.topology.as_ref().is_some_and(|t| !t.valid_for(&self.atoms)) { bail!("Invalid saved residue topology"); }
        for frame in &self.trajectory.frames {
            validate_geometry(&frame.atoms, &frame.pos)?;
        }
        if !self.trajectory.frames.is_empty()
            && self.trajectory.current_frame >= self.trajectory.frames.len()
        {
            bail!("Invalid trajectory frame");
        }
        if !self.camera.radius.is_finite()
            || self.camera.radius <= 0.0
            || !self.camera.target.is_finite()
            || !self.camera.orientation.is_finite()
            || self.camera.orientation.length_squared() < 1e-8
        {
            bail!("Invalid saved camera");
        }
        for value in [
            self.spectrum.width_nm,
            self.spectrum.width_ev,
            self.frequency_settings.thermo_temp_k,
            self.frequency_settings.frequency_scale,
            self.frequency_settings.rot_symmetry,
            self.frequency_settings.spectrum_width_cm1,
        ] {
            if !value.is_finite() || value <= 0.0 {
                bail!("Invalid saved analysis settings");
            }
        }
        self.frequencies.validate().map_err(anyhow::Error::msg)?;
        self.orbitals
            .as_ref()
            .map(Orbitals::prepare)
            .transpose()
            .map(Option::unwrap_or_default)
    }

    fn apply(mut self, world: &mut World, orbitals: OrbitalState) {
        self.camera.orientation = self.camera.orientation.normalize();
        self.camera.rotating = false;
        self.camera.panning = false;
        self.trajectory.playing = false;
        self.trajectory.accum = 0.0;
        self.trajectory.overlay_dirty = true;
        self.trajectory.last_applied = Some(self.trajectory.current_frame);
        self.calculation.config.nproc = self.calculation.config.nproc.clamp(1, 4);
        self.calculation.summary_open = false;
        self.calculation.summary_confirm_close = false;
        self.optimization.method.nproc = self.optimization.method.nproc.clamp(1, 4);
        // Executable paths belong to this installation, not to a portable project.
        self.optimization.paths.clone_from(
            &world
                .resource::<crate::qchem_interfaces::xtb_optimize::XtbPanelState>()
                .paths,
        );
        self.optimization.summary_open = false;
        self.optimization.energy_plot_open = false;
        self.frequency_settings.method.nproc = self.frequency_settings.method.nproc.clamp(1, 4);
        self.frequency_settings.selected_mode = None;
        self.spectrum.reference.nproc = self.spectrum.reference.nproc.clamp(1, 4);
        self.spectrum
            .run_config
            .xtb4stda_path
            .clone_from(&world.resource::<UvVisState>().run_config.xtb4stda_path);
        self.spectrum.orbital_link = Default::default();
        if self.linked_orbitals {
            self.spectrum.orbital_link.bind(&orbitals);
        }
        self.settings.geometry_dirty = true;
        self.settings.bond_topology_dirty = true;
        self.settings.lighting_dirty = true;
        let mut mol = Molecule {
            atoms: self.atoms,
            pos: self.positions,
            bonds: vec![],
            hydrogen_bonds: vec![],
            topology: self.topology,
        };
        mol.recompute_bonds(self.settings.bond_thresh_scale, self.settings.hbond_cutoff);
        world.insert_resource(mol);
        world.insert_resource(XyzBuffer {
            text: self.xyz_draft,
            ..Default::default()
        });
        world.insert_resource(self.settings);
        world.insert_resource(self.camera);
        world.insert_resource(self.trajectory);
        world.insert_resource(self.calculation);
        world.insert_resource(self.optimization);
        world.insert_resource(self.frequency_settings);
        world
            .resource_mut::<XtbFrequencyTask>()
            .restore_archive(self.frequencies);
        world.insert_resource(self.spectrum);
        world.insert_resource(orbitals);
        world.insert_resource(self.layout);
        world.insert_resource(crate::molecule_builder::builder_ui::ZMatrixBuilderState::default());
        world.insert_resource(crate::molecule_builder::builder_ui::EditorRotateState::default());
        world
            .resource_mut::<crate::measurements::Measurements>()
            .reset_all();
        world.resource_mut::<ProjectState>().restore_ui = Some(self.workspace);
        world.write_message(crate::events::MoleculeChanged::parse_xyz(false));
    }
}

fn validate_geometry(atoms: &[String], positions: &[Vec3]) -> Result<()> {
    if atoms.len() != positions.len()
        || positions.iter().any(|p| !p.is_finite())
        || atoms
            .iter()
            .any(|a| a.is_empty() || a.chars().any(char::is_whitespace))
    {
        bail!("Invalid structure coordinates in project");
    }
    Ok(())
}

/// Write beside the destination, then replace it only after the complete file is synced.
pub(crate) fn atomic_write(path: &Path, bytes: &[u8]) -> Result<()> {
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let temporary = parent.join(format!(
        ".beavyr-{}-{}.tmp",
        std::process::id(),
        SERIAL.fetch_add(1, Ordering::Relaxed)
    ));
    let result = (|| {
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)?;
        file.write_all(bytes)?;
        file.sync_all()?;
        fs::rename(&temporary, path)?;
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(temporary);
    }
    result
}

fn busy(world: &World) -> bool {
    world
        .resource::<crate::qchem_panel::run::QcRunTask>()
        .is_running()
        || world
            .resource::<crate::qchem_interfaces::xtb_optimize::XtbOptimizationTask>()
            .is_running()
        || world.resource::<XtbFrequencyTask>().is_running()
        || world
            .resource::<crate::uvvis::run::SpectrumTask>()
            .is_running()
        || world
            .resource::<crate::conformer::ConformerRun>()
            .is_running()
        || world
            .resource::<crate::ts_generation::TsGeneration>()
            .is_running()
}

/// Runs before scene refresh, outside the egui drawing pass.
pub fn process(world: &mut World) {
    let Some(request) = world.resource_mut::<ProjectState>().request.take() else {
        return;
    };
    let result: Result<String> = (|| match request {
        Request::Save(path, workspace) => {
            let snapshot = Snapshot::capture(world, workspace)?;
            snapshot.validate()?;
            let text = ron::ser::to_string(&snapshot)?;
            atomic_write(&path, text.as_bytes())
                .with_context(|| format!("Cannot save {}", path.display()))?;
            world.resource_mut::<ProjectState>().path = Some(path);
            Ok("Project saved".into())
        }
        Request::Open(path, workspace) => {
            if busy(world) {
                bail!("Stop running calculations before opening a project");
            }
            let text = fs::read_to_string(&path)
                .with_context(|| format!("Cannot read {}", path.display()))?;
            let snapshot: Snapshot = ron::from_str(&text).context("Invalid Beavyr project")?;
            let orbitals = snapshot.validate()?;
            let previous = Snapshot::capture(world, workspace)?;
            snapshot.apply(world, orbitals);
            let mut state = world.resource_mut::<ProjectState>();
            state.previous = Some(previous);
            state.path = Some(path);
            Ok("Project opened; Restore previous session is available".into())
        }
        Request::Restore(workspace) => {
            if busy(world) {
                bail!("Stop running calculations before restoring a session");
            }
            let Some(snapshot) = world.resource::<ProjectState>().previous.as_ref() else {
                bail!("No previous session");
            };
            let orbitals = snapshot.validate()?;
            let previous = Snapshot::capture(world, workspace)?;
            let snapshot = world
                .resource_mut::<ProjectState>()
                .previous
                .replace(previous)
                .unwrap();
            snapshot.apply(world, orbitals);
            world.resource_mut::<ProjectState>().path = None;
            Ok("Previous session restored".into())
        }
    })();
    let mut state = world.resource_mut::<ProjectState>();
    state.error = result.is_err();
    state.message = Some(match result {
        Ok(message) => message,
        Err(error) => format!("{error:#}"),
    });
}

#[cfg(test)]
mod tests;
