//! Surface sampling settings. Larger levels use finer grids, subject to the existing point budget.
use super::OrbitalState;
pub const MIN_RESOLUTION: u32 = 1;
pub const MAX_RESOLUTION: u32 = 10;
pub(super) const DEFAULT_RESOLUTION: u32 = 7;
const SPACING: [f32; 11] = [
    0.0, 0.40, 0.30, 0.25, 0.20, 0.15, 0.12, 0.09, 0.07, 0.05, 0.04,
];
impl OrbitalState {
    pub fn grid_spacing(&self) -> f32 {
        SPACING[self.resolution.clamp(MIN_RESOLUTION, MAX_RESOLUTION) as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn extra_resolution_levels_are_finer_and_keep_the_existing_levels() {
        let mut state = OrbitalState::default();
        for (level, expected) in [(6, 0.12), (8, 0.07), (9, 0.05), (10, 0.04), (11, 0.04)] {
            state.resolution = level;
            assert_eq!(state.grid_spacing(), expected);
        }
    }
}
