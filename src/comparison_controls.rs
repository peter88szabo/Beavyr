//! Additional controls for the two private comparison views.
use super::{analysis, State};
use crate::ui_style::ResponseExt;
use bevy::{
    prelude::*,
    tasks::{block_on, futures_lite::future, AsyncComputeTaskPool, Task},
};
use bevy_egui::egui;
use std::time::Instant;
type MappingResult = Result<(u64, u64, Vec<usize>), String>;
pub struct Controls {
    pub colors: [Option<Vec<Color>>; 2],
    pub mapping: Vec<usize>,
    text: String,
    mapping_open: bool,
    task: Option<Task<MappingResult>>,
    linked: bool,
    frame: usize,
    playing: bool,
    fps: f32,
    clock: Instant,
    fraction: f32,
    colored: bool,
    maximum: f32,
    rmsd: Option<f64>,
    dirty: bool,
    keys: Option<(u64, u64)>,
}
impl Default for Controls {
    fn default() -> Self {
        Self {
            colors: [None, None],
            mapping: vec![],
            text: String::new(),
            mapping_open: false,
            task: None,
            linked: true,
            frame: 0,
            playing: false,
            fps: 20.,
            clock: Instant::now(),
            fraction: 0.,
            colored: false,
            maximum: 1.,
            rmsd: None,
            dirty: true,
            keys: None,
        }
    }
}
pub fn poll(ctx: &egui::Context, state: &mut State) {
    if let Some(task) = &mut state.extra.task {
        if let Some(result) = block_on(future::poll_once(task)) {
            state.extra.task = None;
            match result {
                Ok((a, b, map)) => {
                    if state.sides[0].molecule.as_ref().map(analysis::key) == Some(a)
                        && state.sides[1].molecule.as_ref().map(analysis::key) == Some(b)
                    {
                        set_mapping(state, map);
                        state.message="Geometric atom mapping found; review the correspondence, especially for equivalent atoms.".into();
                    } else {
                        state.message =
                            "Structures changed while mapping; run automatic mapping again.".into();
                    }
                }
                Err(e) => crate::runtime_warning::push(ctx, "Atom mapping", &e),
            }
        } else {
            ctx.request_repaint();
        }
    }
    if !state.open {
        state.extra.playing = false;
        state.extra.clock = Instant::now();
    }
}
fn set_mapping(state: &mut State, map: Vec<usize>) {
    state.extra.text = map
        .iter()
        .map(|i| (i + 1).to_string())
        .collect::<Vec<_>>()
        .join(" ");
    state.extra.mapping = map;
    state.extra.dirty = true;
}
pub fn refresh(state: &mut State) {
    let (Some(a), Some(b)) = (&state.sides[0].molecule, &state.sides[1].molecule) else {
        return;
    };
    let keys = (analysis::key(a), analysis::key(b));
    if state.extra.keys == Some(keys) && !state.extra.dirty {
        return;
    }
    if analysis::validate(a, b, &state.extra.mapping).is_err() {
        state.extra.mapping = if a.atoms == b.atoms {
            (0..a.atoms.len()).collect()
        } else {
            vec![]
        };
        state.extra.text = state
            .extra
            .mapping
            .iter()
            .map(|i| (i + 1).to_string())
            .collect::<Vec<_>>()
            .join(" ");
    }
    state.extra.keys = Some(keys);
    state.extra.dirty = false;
    state.extra.colors = [None, None];
    state.extra.rmsd = None;
    if let Ok(fit) = analysis::aligned(a, b, &state.extra.mapping) {
        state.extra.rmsd = Some(fit.rmsd);
        if state.extra.colored {
            let colors: Vec<_> = fit
                .distances
                .iter()
                .map(|&d| analysis::color(d, state.extra.maximum))
                .collect();
            let mut other = vec![Color::WHITE; colors.len()];
            for (i, &j) in state.extra.mapping.iter().enumerate() {
                other[j] = colors[i];
            }
            state.extra.colors = [Some(colors), Some(other)];
        }
    }
    for side in &mut state.sides {
        side.dirty = true;
    }
}
pub fn linked(state: &State) -> bool {
    state.extra.linked
}
pub fn controls(ui: &mut egui::Ui, state: &mut State, actions: &mut Vec<(usize, u8)>) {
    ui.horizontal_wrapped(|ui| {
        if ui.button("Atom mapping…").clicked_once() {
            state.extra.mapping_open = true;
        }
        if ui
            .checkbox(&mut state.extra.colored, "Displacement colours")
            .changed()
        {
            state.extra.dirty = true;
        }
        if state.extra.colored {
            ui.colored_label(egui::Color32::from_rgb(0, 64, 255), "0 Å");
            if ui
                .add(
                    egui::DragValue::new(&mut state.extra.maximum)
                        .range(0.001..=100.)
                        .speed(0.01)
                        .suffix(" Å"),
                )
                .changed()
            {
                state.extra.dirty = true;
            }
            ui.colored_label(egui::Color32::from_rgb(255, 64, 0), "maximum");
        }
        if let Some(rmsd) = state.extra.rmsd {
            ui.label(format!("Aligned RMSD {rmsd:.4} Å"));
        }
    });
    let count = analysis::linked_count(state.sides[0].frames.len(), state.sides[1].frames.len());
    state.extra.frame = state.extra.frame.min(count - 1);
    let mut changed = false;
    ui.horizontal_wrapped(|ui| {
        changed |= ui
            .checkbox(&mut state.extra.linked, "Link frames")
            .changed();
        ui.add_enabled_ui(state.extra.linked && count > 1, |ui| {
            if ui
                .button(if state.extra.playing { "⏸" } else { "▶" })
                .clicked_once()
            {
                if !state.extra.playing && state.extra.frame + 1 == count {
                    state.extra.frame = 0;
                    changed = true;
                }
                state.extra.playing = !state.extra.playing;
                state.extra.fraction = 0.;
                state.extra.clock = Instant::now();
            }
            let mut frame = state.extra.frame + 1;
            if ui
                .add(egui::Slider::new(&mut frame, 1..=count).text("Frame"))
                .changed()
            {
                state.extra.frame = frame - 1;
                changed = true;
            }
            ui.add(
                egui::DragValue::new(&mut state.extra.fps)
                    .range(1.0..=100.)
                    .suffix(" fps"),
            );
        });
    });
    if !state.extra.linked {
        state.extra.playing = false;
    }
    let elapsed = state.extra.clock.elapsed().as_secs_f32().min(0.25);
    state.extra.clock = Instant::now();
    if state.extra.playing {
        state.extra.fraction += elapsed * state.extra.fps;
        let steps = state.extra.fraction as usize;
        state.extra.fraction -= steps as f32;
        if steps > 0 {
            state.extra.frame = (state.extra.frame + steps).min(count - 1);
            changed = true;
        }
        if state.extra.frame + 1 == count {
            state.extra.playing = false;
        }
        ui.ctx().request_repaint();
    }
    if state.extra.linked {
        for (i, side) in state.sides.iter_mut().enumerate() {
            if side.frames.len() > 1
                && (changed || side.frame != state.extra.frame.min(side.frames.len() - 1))
            {
                side.frame = state.extra.frame.min(side.frames.len() - 1);
                actions.push((i, 4));
            }
        }
    }
}
pub fn mapping_window(ctx: &egui::Context, state: &mut State) -> Option<Rect> {
    if !state.extra.mapping_open {
        return None;
    }
    let mut open = true;
    let window=egui::Window::new("Comparison atom mapping").open(&mut open).resizable(true).default_size([580.,330.]).show(ctx,|ui|{
        ui.label("B atom indices corresponding to A atoms 1, 2, 3, …");
        ui.add(egui::TextEdit::multiline(&mut state.extra.text).desired_rows(4));
        let mut mode=None;
        ui.horizontal(|ui|{
            if ui.button("Apply mapping").clicked_once(){mode=Some(0);}
            if ui.button("Use atom order").clicked_once(){mode=Some(1);}
            if ui.add_enabled(state.extra.task.is_none(),egui::Button::new("Automatic mapping")).clicked_once(){mode=Some(2);}
            if ui.button("Copy mapping").clicked_once(){ctx.copy_text(state.extra.text.clone());}
        });
        if let Some(mode)=mode {
            let result=(||->Result<(),String>{
                let a=state.sides[0].molecule.as_ref().ok_or("Load structure A.")?;
                let b=state.sides[1].molecule.as_ref().ok_or("Load structure B.")?;
                if mode==2{
                    let a=a.clone();let b=b.clone();state.extra.task=Some(AsyncComputeTaskPool::get().spawn(async move{let map=analysis::automatic(&a,&b)?;Ok((analysis::key(&a),analysis::key(&b),map))}));
                }else{
                    let map=if mode==0{analysis::parse(a,b,&state.extra.text)?}else{let map=(0..a.atoms.len()).collect::<Vec<_>>();analysis::validate(a,b,&map)?;map};
                    set_mapping(state,map);
                }Ok(())
            })();
            if let Err(e)=result{crate::runtime_warning::push(ctx,"Atom mapping",&e);}
        }
        if state.extra.task.is_some(){ui.spinner();ctx.request_repaint();}
        ui.small("Automatic correspondence is based on elements and geometry. Equivalent atoms can have multiple valid mappings; use manual indices when atom identity matters.");
    });
    state.extra.mapping_open = open;
    window.map(|w| {
        let r = w.response.rect;
        Rect::new(r.min.x, r.min.y, r.max.x, r.max.y)
    })
}
