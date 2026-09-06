//! Layout mode for the control panels.
//!
//! Two layouts coexist and are switchable at runtime, so the newer tabbed
//! drawer can be abandoned instantly without a rebuild if it turns out to be
//! worse in practice:
//!
//! * `Classic` -- the original always-visible right panel with every section
//!   stacked as collapsing headers, plus the molecule editor as a permanent
//!   left panel.
//! * `Windowed` -- a permanent, narrow icon rail down the left edge, where
//!   each icon toggles that tool's own floating window. Several tools can be
//!   open and arranged at once, the viewport keeps both edges except for the
//!   rail itself, and the molecule editor is reached from a small round
//!   button on the right.
//!
//! Nothing here draws chemistry -- it only decides *where* the existing
//! panels go, so both layouts call exactly the same section bodies.

use bevy::prelude::*;
use bevy_egui::egui;

/// Which section the drawer is showing. One tool per tab -- grouping
/// several tools behind one icon just recreates the scrolling accordion the
/// rail exists to replace.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tab {
    Structure,
    Representation,
    Optimize,
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
    pub const COUNT: usize = 12;

    /// Rail order, matching the workflow order the panels were arranged in.
    pub const ALL: [Tab; Tab::COUNT] = [
        Tab::Structure,
        Tab::Representation,
        Tab::Optimize,
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

    /// A single glyph for the rail. These come from the emoji set egui
    /// bundles by default (the same family its own demo app uses), so they
    /// render without shipping an extra icon font.
    pub fn icon(self) -> &'static str {
        match self {
            Tab::Structure => "🖹",
            Tab::Representation => "⬢",
            Tab::Optimize => "🔧",
            Tab::Vibrations => "📈",
            Tab::UvVis => "🌈",
            Tab::Trajectory => "▶",
            Tab::Compare => "⚖",
            Tab::Measure => "📏",
            Tab::Surface => "◉",
            Tab::Diagnostics => "⚠",
            Tab::Appearance => "🎨",
            Tab::Export => "💾",
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
            Tab::Structure => [420.0, 520.0],
            Tab::Representation => [420.0, 480.0],
            Tab::Optimize => [430.0, 420.0],
            Tab::Vibrations => [460.0, 520.0],
            Tab::UvVis => [560.0, 560.0],
            Tab::Trajectory => [420.0, 520.0],
            Tab::Compare => [380.0, 360.0],
            Tab::Measure => [400.0, 460.0],
            Tab::Surface => [360.0, 320.0],
            Tab::Diagnostics => [420.0, 380.0],
            Tab::Appearance => [420.0, 560.0],
            Tab::Export => [400.0, 420.0],
        }
    }

    pub fn title(self) -> &'static str {
        match self {
            Tab::Structure => "Structure (XYZ)",
            Tab::Representation => "Representation",
            Tab::Optimize => "Geometry Optimization",
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
    /// `false` restores the original layout exactly.
    pub windowed: bool,
    /// One flag per entry of `Tab::ALL`, in the same order.
    pub open: [bool; Tab::COUNT],
    /// Whether the molecule editor's floating window is open (windowed
    /// layout only; the classic layout keeps it as a docked left panel).
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
    if windowed {
        if !*open {
            return false;
        }
        let ctx = ui.ctx().clone();
        let mut still_open = *open;
        let response = egui::Window::new(title)
            .open(&mut still_open)
            .resizable(true)
            .vscroll(true)
            .default_size(default_size)
            .show(&ctx, |ui| {
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
        let response = ui.collapsing(title, add_contents);
        response.openness > 0.0
    }
}

/// The always-visible vertical rail of tool icons. Each button toggles that
/// tool's window, and stays highlighted while the window is open.
pub fn icon_rail(ui: &mut egui::Ui, open: &mut [bool; Tab::COUNT]) {
    ui.vertical_centered(|ui| {
        ui.add_space(6.0);
        for (index, tab) in Tab::ALL.iter().enumerate() {
            let is_open = open[index];
            let button = egui::Button::new(egui::RichText::new(tab.icon()).size(19.0))
                .min_size(egui::vec2(34.0, 34.0))
                .selected(is_open);
            if ui.add(button).on_hover_text(tab.title()).clicked() {
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
    fn every_tab_has_an_icon_and_a_title() {
        for tab in Tab::ALL {
            assert!(!tab.icon().is_empty());
            assert!(!tab.title().is_empty());
        }
    }

    #[test]
    fn tabs_are_all_distinct() {
        let titles: Vec<&str> = Tab::ALL.iter().map(|t| t.title()).collect();
        let mut sorted = titles.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), titles.len(), "duplicate tab titles: {titles:?}");
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
