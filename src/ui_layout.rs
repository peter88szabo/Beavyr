//! Layout mode for the control panels.
//!
//! General tools can use either layout at runtime:
//!
//! * `Classic` -- the original always-visible right panel with every section
//!   stacked as collapsing headers.
//! * `Windowed` -- a permanent, narrow icon rail down the left edge, where
//!   each icon toggles that tool's own floating window. Several tools can be
//!   open and arranged at once, the viewport keeps both edges except for the
//!   rail itself.
//!
//! The molecule editor and Quantum Chemistry calculator have their own
//! launchers on the right, float by default, and offer explicit docking.
//!
//! Nothing here draws chemistry -- it only decides *where* the existing
//! panels go, so both layouts call exactly the same section bodies.

use crate::ui_icons::{Icon, IconButton};
use crate::ui_style::ResponseExt;
use bevy::prelude::*;
use bevy_egui::egui;

/// Which section the drawer is showing. One tool per tab -- grouping
/// several tools behind one icon just recreates the scrolling accordion the
/// rail exists to replace.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tab {
    Structure,
    RecentStructures,
    Representation,
    Optimize,
    TsGeneration,
    Conformers,
    Vibrations,
    UvVis,
    Trajectory,
    Compare,
    Measure,
    Surface,
    Diagnostics,
    Appearance,
    Export,
}

impl Tab {
    pub const COUNT: usize = 15;

    /// Rail order, matching the workflow order the panels were arranged in.
    pub const ALL: [Tab; Tab::COUNT] = [
        Tab::Structure,
        Tab::RecentStructures,
        Tab::Representation,
        Tab::Optimize,
        Tab::TsGeneration,
        Tab::Conformers,
        Tab::Vibrations,
        Tab::UvVis,
        Tab::Trajectory,
        Tab::Compare,
        Tab::Measure,
        Tab::Surface,
        Tab::Diagnostics,
        Tab::Appearance,
        Tab::Export,
    ];

    /// A vector icon for the rail, independent of the installed fonts.
    pub fn icon(self) -> Icon {
        match self {
            Tab::Structure => Icon::Structure,
            Tab::RecentStructures => Icon::History,
            Tab::Representation => Icon::Representation,
            Tab::Optimize => Icon::Optimize,
            Tab::TsGeneration => Icon::TransitionState,
            Tab::Conformers => Icon::Conformers,
            Tab::Vibrations => Icon::Vibrations,
            Tab::UvVis => Icon::Spectrum,
            Tab::Trajectory => Icon::Trajectory,
            Tab::Compare => Icon::Compare,
            Tab::Measure => Icon::Measure,
            Tab::Surface => Icon::Surface,
            Tab::Diagnostics => Icon::Diagnostics,
            Tab::Appearance => Icon::Appearance,
            Tab::Export => Icon::Export,
        }
    }

    /// Position in `ALL`, which is also this tab's index into the open-window
    /// flags -- keeping one source of truth for that pairing.
    pub fn index(self) -> usize {
        Tab::ALL
            .iter()
            .position(|&t| t == self)
            .expect("every Tab variant is listed in Tab::ALL")
    }

    /// A sensible opening size per tool -- the wide, table-heavy ones need
    /// more room than the short control lists.
    pub fn default_size(self) -> [f32; 2] {
        match self {
            Tab::Structure => [500.0, 620.0],
            Tab::RecentStructures => [460.0, 620.0],
            Tab::Representation => [460.0, 520.0],
            Tab::Optimize => [430.0, 420.0],
            Tab::TsGeneration => [580.0, 720.0],
            // Wide and tall: the conformer table earns the room, one row per
            // conformer with energies and populations alongside.
            Tab::Conformers => [520.0, 620.0],
            // Tall: two run blocks, the settings above them and a mode list
            // that is worth seeing more than a handful of rows of.
            Tab::Vibrations => [540.0, 660.0],
            Tab::UvVis => [560.0, 560.0],
            Tab::Trajectory => [500.0, 620.0],
            Tab::Compare => [480.0, 540.0],
            Tab::Measure => [500.0, 640.0],
            Tab::Surface => [520.0, 660.0],
            Tab::Diagnostics => [500.0, 520.0],
            Tab::Appearance => [520.0, 640.0],
            Tab::Export => [460.0, 520.0],
        }
    }

    pub fn title(self) -> &'static str {
        match self {
            Tab::Structure => "Structure (XYZ)",
            Tab::RecentStructures => "Recent Structures",
            Tab::Representation => "Representation",
            Tab::Optimize => "Geometry Optimization",
            Tab::TsGeneration => "TS Generation",
            Tab::Conformers => "Conformer Search",
            Tab::Vibrations => "Vibrations",
            Tab::UvVis => "UV-Vis (TD-DFT)",
            Tab::Trajectory => "Trajectory",
            Tab::Compare => "Structure Comparison",
            Tab::Measure => "Measurements",
            Tab::Surface => "Surface Tools",
            Tab::Diagnostics => "Diagnostics",
            Tab::Appearance => "Appearance",
            Tab::Export => "Export Image",
        }
    }
}

#[derive(Resource)]
pub struct UiLayout {
    /// `false` presents the general tools in a side panel.
    pub windowed: bool,
    /// One flag per entry of `Tab::ALL`, in the same order.
    pub open: [bool; Tab::COUNT],
    /// Whether the molecule editor is open. It floats by default in either
    /// layout, and the user can dock it from its header.
    pub builder_open: bool,
}

impl Default for UiLayout {
    fn default() -> Self {
        Self {
            windowed: true,
            open: [false; Tab::COUNT],
            builder_open: false,
        }
    }
}

/// Draws one control section in whichever layout is active, so both layouts
/// share a single copy of every section body.
///
/// In the windowed layout each section is its own floating `Window`, so any
/// number of tools can be open at once and arranged freely; in the classic
/// layout it stays a collapsing header inside the right panel.
///
/// Returns whether the body was actually shown, which a couple of callers
/// need in order to skip expensive work while their section is hidden.
pub fn section<R>(
    ui: &mut egui::Ui,
    windowed: bool,
    open: &mut bool,
    rects: &mut Vec<egui::Rect>,
    default_size: [f32; 2],
    title: &str,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> bool {
    section_impl(
        ui,
        windowed,
        open,
        rects,
        default_size,
        title,
        true,
        add_contents,
    )
}

/// A section whose body owns its scrolling, so toolbars stay visible while
/// an editor uses the remaining window height.
pub fn editor_section<R>(
    ui: &mut egui::Ui,
    windowed: bool,
    open: &mut bool,
    rects: &mut Vec<egui::Rect>,
    default_size: [f32; 2],
    title: &str,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> bool {
    section_impl(
        ui,
        windowed,
        open,
        rects,
        default_size,
        title,
        false,
        add_contents,
    )
}

fn section_impl<R>(
    ui: &mut egui::Ui,
    windowed: bool,
    open: &mut bool,
    rects: &mut Vec<egui::Rect>,
    default_size: [f32; 2],
    title: &str,
    scroll: bool,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> bool {
    if windowed {
        if !*open {
            return false;
        }
        let ctx = ui.ctx().clone();
        let mut still_open = *open;
        let response = egui::Window::new(title)
            .frame(crate::ui_style::window_frame(&ctx))
            .open(&mut still_open)
            .resizable(true)
            .vscroll(scroll)
            .default_size(default_size)
            .max_height((ctx.viewport_rect().height() - 64.0).max(240.0))
            .show(&ctx, |ui| {
                crate::ui_style::modern(ui);
                add_contents(ui);
            });
        *open = still_open;
        // The rect is recorded so the caller can tell the 3D picker this area
        // is covered -- otherwise a click on a window would also select
        // whatever atom happens to sit behind it.
        if let Some(r) = &response {
            rects.push(r.response.rect);
        }
        response.is_some()
    } else {
        ui.add_space(8.0);
        let response = ui.collapsing(title, |ui| {
            crate::ui_style::modern(ui);
            add_contents(ui)
        });
        response.openness > 0.0
    }
}

/// The always-visible vertical rail of tool icons. Each button toggles that
/// tool's window, and stays highlighted while the window is open.
pub fn icon_rail(ui: &mut egui::Ui, open: &mut [bool; Tab::COUNT]) {
    ui.vertical_centered(|ui| {
        ui.add_space(6.0);
        for (index, tab) in Tab::ALL.iter().enumerate() {
            // Geometry optimisation is run from the Quantum Chemistry panel
            // beside the editor, which does the same job for every program.
            // A second button for it only offered two ways into one thing.
            if *tab == Tab::Optimize {
                open[index] = false;
                continue;
            }
            let is_open = open[index];
            let button = IconButton::new(tab.icon())
                .label(tab.title())
                .selected(is_open);
            if ui.add(button).on_hover_text(tab.title()).clicked_once() {
                open[index] = !is_open;
            }
            ui.add_space(3.0);
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_click_keeps_a_new_tool_window_open_across_layout_passes() {
        let ctx = egui::Context::default();
        let mut open = [false; Tab::COUNT];
        let mut frame = |events| {
            let mut output = ctx.run_ui(
                egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(1280.0, 720.0),
                    )),
                    events,
                    ..Default::default()
                },
                |ui| {
                    egui::Panel::left("test_tool_rail")
                        .exact_size(50.0)
                        .show(ui, |ui| {
                            icon_rail(ui, &mut open);
                        });
                    if open[Tab::Structure.index()] {
                        egui::Window::new("New tool window")
                            .default_pos([200.0, 50.0])
                            .show(&ctx, |ui| {
                                ui.label("Tool controls");
                            });
                    }
                },
            );
            output.textures_delta.clear();
            output
        };
        frame(vec![]);
        frame(vec![]);
        let pos = egui::pos2(25.0, 25.0);
        frame(vec![
            egui::Event::PointerMoved(pos),
            egui::Event::PointerButton {
                pos,
                button: egui::PointerButton::Primary,
                pressed: true,
                modifiers: egui::Modifiers::NONE,
            },
        ]);
        frame(vec![egui::Event::PointerButton {
            pos,
            button: egui::PointerButton::Primary,
            pressed: false,
            modifiers: egui::Modifiers::NONE,
        }]);
        frame(vec![]);
        let output = frame(vec![]);
        assert!(
            open[Tab::Structure.index()],
            "the layout retry must not toggle the window closed"
        );
        assert!(output.shapes.iter().any(|shape| matches!(&shape.shape,
            egui::Shape::Text(text) if text.galley.text() == "Tool controls")));
    }

    #[test]
    fn every_tab_has_an_icon_and_a_title() {
        for tab in Tab::ALL {
            assert!(!tab.icon().label().is_empty());
            assert!(!tab.title().is_empty());
        }
    }

    #[test]
    fn tabs_are_all_distinct() {
        let titles: Vec<&str> = Tab::ALL.iter().map(|t| t.title()).collect();
        let mut sorted = titles.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(
            sorted.len(),
            titles.len(),
            "duplicate tab titles: {titles:?}"
        );
    }

    #[test]
    fn every_variant_appears_exactly_once_in_all() {
        // `index()` and the rail both rely on this.
        for tab in Tab::ALL {
            assert_eq!(Tab::ALL.iter().filter(|&&t| t == tab).count(), 1, "{tab:?}");
        }
        assert_eq!(Tab::ALL.len(), Tab::COUNT);
    }

    #[test]
    fn index_matches_position_in_all() {
        for (i, tab) in Tab::ALL.iter().enumerate() {
            assert_eq!(tab.index(), i);
        }
    }

    #[test]
    fn the_open_flags_line_up_with_the_rail() {
        // `icon_rail` indexes the flag array by position in `Tab::ALL`, so a
        // mismatch here would toggle the wrong tool's window.
        let layout = UiLayout::default();
        assert_eq!(layout.open.len(), Tab::ALL.len());
    }

    #[test]
    fn nothing_is_open_on_first_launch() {
        // The rail is the entry point; opening ten windows unprompted would
        // bury the viewport.
        let layout = UiLayout::default();
        assert!(layout.open.iter().all(|&o| !o));
    }
}
