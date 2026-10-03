//! Navigation uses original file indices even when the orbital list is filtered.
use super::{molden::MolecularOrbital, OrbitalState};
use crate::ui_style::ResponseExt;
use bevy_egui::egui;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Occupancy {
    #[default]
    All,
    Occupied,
    Virtual,
}

#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct Browser {
    pub jump: String,
    pub occupancy: Occupancy,
    pub energy_range: bool,
    pub minimum: String,
    pub maximum: String,
    pub error: Option<String>,
    pub reveal: bool,
}

impl Browser {
    pub fn visible(&self, orbitals: &[MolecularOrbital]) -> Result<Vec<usize>, String> {
        let bound = |text: &str| -> Result<Option<f64>, String> {
            if text.trim().is_empty() {
                return Ok(None);
            }
            text.trim()
                .parse::<f64>()
                .ok()
                .filter(|v| v.is_finite())
                .map(Some)
                .ok_or_else(|| "Enter a finite energy in Eh.".to_string())
        };
        let (minimum, maximum) = if self.energy_range {
            (bound(&self.minimum)?, bound(&self.maximum)?)
        } else {
            (None, None)
        };
        if minimum.zip(maximum).is_some_and(|(a, b)| a > b) {
            return Err("Minimum energy must not exceed maximum energy.".into());
        }
        Ok(orbitals
            .iter()
            .enumerate()
            .filter_map(|(index, orbital)| {
                let occupied = orbital.occupation > 1.0e-6;
                let keep = match self.occupancy {
                    Occupancy::All => true,
                    Occupancy::Occupied => occupied,
                    Occupancy::Virtual => !occupied,
                };
                (keep
                    && minimum.is_none_or(|v| orbital.energy >= v)
                    && maximum.is_none_or(|v| orbital.energy <= v))
                .then_some(index)
            })
            .collect())
    }
    pub fn jump_index(&self, count: usize) -> Result<usize, String> {
        self.jump
            .trim()
            .parse::<usize>()
            .ok()
            .and_then(|n| n.checked_sub(1))
            .filter(|&index| index < count)
            .ok_or_else(|| format!("Enter an orbital number from 1 to {count}."))
    }
    pub fn reset_filters(&mut self) {
        self.occupancy = Occupancy::All;
        self.energy_range = false;
        self.minimum.clear();
        self.maximum.clear();
        self.error = None;
    }
}

pub fn step(current: Option<usize>, visible: &[usize], delta: isize) -> Option<usize> {
    if visible.is_empty() {
        return None;
    }
    let position = current.and_then(|index| visible.iter().position(|&i| i == index));
    let next = match position {
        Some(p) => (p as isize + delta).clamp(0, visible.len() as isize - 1) as usize,
        None if delta < 0 => visible.len() - 1,
        None => 0,
    };
    Some(visible[next])
}

pub fn select(state: &mut OrbitalState, index: usize) {
    state.selected = Some(index);
    state.browser.jump = (index + 1).to_string();
    state.browser.error = None;
    state.browser.reveal = true;
}

pub fn keyboard_id(ui: &egui::Ui) -> egui::Id {
    ui.make_persistent_id("orbital_navigation")
}

pub fn controls(
    ui: &mut egui::Ui,
    state: &mut OrbitalState,
    orbitals: &[MolecularOrbital],
) -> Vec<usize> {
    let mut jump = false;
    let mut direction = 0;
    let index_id = ui.make_persistent_id("orbital_index");
    ui.horizontal(|ui| {
        ui.label("#");
        let response = ui.add(
            egui::TextEdit::singleline(&mut state.browser.jump)
                .id(index_id)
                .desired_width(54.0)
                .hint_text("Index"),
        );
        jump = response.lost_focus() && ui.input(|input| input.key_pressed(egui::Key::Enter));
        if !response.has_focus() && !jump {
            state.browser.jump = state
                .selected
                .map(|index| (index + 1).to_string())
                .unwrap_or_default();
        }
        response.on_hover_text(
            "Type an orbital number and press Enter. Use ↑/↓ over the list to browse.",
        );
        if ui
            .small_button("◀")
            .on_hover_text("Previous matching orbital")
            .clicked_once()
        {
            direction = -1;
        }
        if ui
            .small_button("▶")
            .on_hover_text("Next matching orbital")
            .clicked_once()
        {
            direction = 1;
        }
        let filtering = state.browser.occupancy != Occupancy::All || state.browser.energy_range;
        ui.menu_button(if filtering { "Filter •" } else { "Filter" }, |ui| {
            for (choice, label) in [
                (Occupancy::All, "All orbitals"),
                (Occupancy::Occupied, "Occupied"),
                (Occupancy::Virtual, "Virtual"),
            ] {
                ui.radio_value(&mut state.browser.occupancy, choice, label);
            }
            ui.separator();
            ui.checkbox(&mut state.browser.energy_range, "Energy range (Eh)");
            ui.add_enabled_ui(state.browser.energy_range, |ui| {
                ui.horizontal(|ui| {
                    ui.add(
                        egui::TextEdit::singleline(&mut state.browser.minimum)
                            .desired_width(78.0)
                            .hint_text("Minimum"),
                    );
                    ui.add(
                        egui::TextEdit::singleline(&mut state.browser.maximum)
                            .desired_width(78.0)
                            .hint_text("Maximum"),
                    );
                });
            });
            if ui.button("Clear filters").clicked_once() {
                state.browser.reset_filters();
            }
        });
    });
    if jump {
        match state.browser.jump_index(orbitals.len()) {
            Ok(index) => {
                state.browser.reset_filters();
                select(state, index);
            }
            Err(error) => state.browser.error = Some(error),
        }
    }
    let visible = match state.browser.visible(orbitals) {
        Ok(indices) => indices,
        Err(error) => {
            ui.colored_label(egui::Color32::from_rgb(230, 160, 75), error);
            Vec::new()
        }
    };
    let keyboard_id = keyboard_id(ui);
    ui.interact(
        ui.available_rect_before_wrap(),
        keyboard_id,
        egui::Sense::focusable_noninteractive(),
    );
    let focused = ui.memory(|memory| memory.has_focus(keyboard_id));
    if (focused || ui.rect_contains_pointer(ui.available_rect_before_wrap()))
        && (focused || !ui.ctx().egui_wants_keyboard_input())
    {
        ui.input_mut(|input| {
            if input.consume_key(egui::Modifiers::NONE, egui::Key::ArrowUp) {
                direction = -1;
            }
            if input.consume_key(egui::Modifiers::NONE, egui::Key::ArrowDown) {
                direction = 1;
            }
            if input.consume_key(egui::Modifiers::NONE, egui::Key::PageUp) {
                direction = -10;
            }
            if input.consume_key(egui::Modifiers::NONE, egui::Key::PageDown) {
                direction = 10;
            }
        });
    }
    if direction != 0 {
        if let Some(index) = step(state.selected, &visible, direction) {
            select(state, index);
            ui.memory_mut(|memory| memory.request_focus(keyboard_id));
        }
    }
    if let Some(error) = &state.browser.error {
        ui.colored_label(egui::Color32::from_rgb(230, 160, 75), error);
    }
    if visible.is_empty() {
        ui.weak("No matching orbitals.");
    }
    visible
}

#[cfg(test)]
mod tests {
    #[test]
    fn keyboard_navigation_and_enter_jump_use_real_egui_focus() {
        use super::*;
        let ctx = egui::Context::default();
        let orbitals = (0..6)
            .map(|i| MolecularOrbital {
                symmetry: String::new(),
                energy: i as f64,
                occupation: if i < 3 { 2.0 } else { 0.0 },
                coefficients: vec![],
            })
            .collect::<Vec<_>>();
        let mut state = OrbitalState::default();
        state.selected = Some(3);
        state.browser.occupancy = Occupancy::Virtual;
        let frame = |state: &mut OrbitalState, events: Vec<egui::Event>, focus: Option<&str>| {
            let mut output = ctx.run_ui(
                egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(500.0, 500.0),
                    )),
                    events,
                    ..Default::default()
                },
                |ui| {
                    controls(ui, state, &orbitals);
                    if let Some(name) = focus {
                        ui.memory_mut(|m| m.request_focus(ui.make_persistent_id(name)));
                    }
                },
            );
            output.textures_delta.clear();
        };
        let key = |key| {
            vec![egui::Event::Key {
                key,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: egui::Modifiers::NONE,
            }]
        };
        frame(&mut state, vec![], Some("orbital_navigation"));
        frame(&mut state, key(egui::Key::ArrowDown), None);
        assert_eq!(state.selected, Some(4));
        frame(&mut state, vec![], Some("orbital_index"));
        state.browser.jump = "2".into();
        frame(&mut state, key(egui::Key::Enter), None);
        assert_eq!(state.selected, Some(1));
        assert_eq!(
            state.browser.occupancy,
            Occupancy::All,
            "explicit jumps reveal hidden orbitals"
        );
    }
    use super::*;
    fn orbitals() -> Vec<MolecularOrbital> {
        [(-2.0, 2.0), (-0.5, 1.0), (0.1, 0.0), (0.3, 0.0)]
            .into_iter()
            .map(|(energy, occupation)| MolecularOrbital {
                symmetry: "a".into(),
                energy,
                occupation,
                coefficients: vec![],
            })
            .collect()
    }
    #[test]
    fn filtering_and_navigation_preserve_file_indices() {
        let browser = Browser {
            occupancy: Occupancy::Virtual,
            ..Default::default()
        };
        let visible = browser.visible(&orbitals()).unwrap();
        assert_eq!(visible, [2, 3]);
        assert_eq!(step(Some(0), &visible, 1), Some(2));
        assert_eq!(step(Some(2), &visible, 1), Some(3));
        assert_eq!(step(Some(3), &visible, 10), Some(3));
        assert_eq!(step(Some(3), &visible, -1), Some(2));
        assert_eq!(step(Some(3), &[], -1), None);
    }
    #[test]
    fn energy_bounds_are_inclusive_and_validate_input() {
        let mut browser = Browser {
            energy_range: true,
            minimum: "-0.5".into(),
            maximum: "0.1".into(),
            ..Default::default()
        };
        assert_eq!(browser.visible(&orbitals()).unwrap(), [1, 2]);
        browser.minimum = "0.2".into();
        assert!(browser.visible(&orbitals()).is_err());
        browser.minimum = "NaN".into();
        assert!(browser.visible(&orbitals()).is_err());
    }
    #[test]
    fn jumping_uses_one_based_numbers_and_rejects_invalid_input() {
        let mut browser = Browser::default();
        for (text, expected) in [("1", 0), ("4", 3), (" 2 ", 1)] {
            browser.jump = text.into();
            assert_eq!(browser.jump_index(4), Ok(expected));
        }
        for text in ["0", "5", "-1", "2.5", "", "abc"] {
            browser.jump = text.into();
            assert!(browser.jump_index(4).is_err());
        }
    }
}
