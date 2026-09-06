//! Line shapes, curve sampling and `.dat` export for simulated spectra.
//!
//! Everything here is unit-agnostic: `x` is whatever the caller's axis is
//! (cm⁻¹ for IR, nm or eV for UV-Vis) and `width` is in those same units. A
//! Gaussian in nm is not a Gaussian in eV, so the caller is responsible for
//! converting its peaks into the axis it wants *before* broadening them.

/// How a stick spectrum is turned into a continuous curve -- or not, in the
/// `Sticks` case.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BroadeningKind {
    Gaussian,
    Lorentzian,
    Voigt,
    /// No broadening at all: a vertical line at each predicted position,
    /// height equal to its raw intensity -- the theoretical stick spectrum, as
    /// opposed to a simulated experimental one.
    Sticks,
}

impl BroadeningKind {
    pub fn is_sticks(self) -> bool {
        matches!(self, BroadeningKind::Sticks)
    }

    pub fn label(self) -> &'static str {
        match self {
            BroadeningKind::Gaussian => "Gaussian",
            BroadeningKind::Lorentzian => "Lorentzian",
            BroadeningKind::Voigt => "Voigt",
            BroadeningKind::Sticks => "Sticks",
        }
    }

    /// Every variant, in the order the selectors offer them.
    pub const ALL: [BroadeningKind; 4] = [
        BroadeningKind::Gaussian,
        BroadeningKind::Lorentzian,
        BroadeningKind::Voigt,
        BroadeningKind::Sticks,
    ];
}

/// A single line-shape value at `x`, centred on `x0`, with intensity `height`
/// and characteristic `width` (interpreted per shape: Gaussian uses it as
/// sigma, Lorentzian as the half-width-at-half-maximum gamma, Voigt as both
/// via the pseudo-Voigt mixing below).
pub fn lineshape(kind: BroadeningKind, x: f64, x0: f64, height: f64, width: f64) -> f64 {
    let width = width.max(1.0e-6);
    match kind {
        BroadeningKind::Gaussian => {
            let z = (x - x0) / width;
            height * (-0.5 * z * z).exp()
        }
        BroadeningKind::Lorentzian => {
            let z = (x - x0) / width;
            height / (1.0 + z * z)
        }
        BroadeningKind::Voigt => {
            // Pseudo-Voigt: a linear mix of Gaussian and Lorentzian with the
            // same width, which is what most spectroscopy software calls a
            // practical "Voigt" profile without the cost of the true Voigt
            // convolution integral.
            0.5 * lineshape(BroadeningKind::Gaussian, x, x0, height, width)
                + 0.5 * lineshape(BroadeningKind::Lorentzian, x, x0, height, width)
        }
        // Sticks have no continuous curve -- callers branch on
        // `BroadeningKind::is_sticks` before ever reaching this function.
        BroadeningKind::Sticks => 0.0,
    }
}

/// The broadened curve's height at a single `x`, summed over every peak.
pub fn intensity_at(peaks: &[(f64, f64)], kind: BroadeningKind, width: f64, x: f64) -> f64 {
    peaks
        .iter()
        .map(|&(x0, height)| lineshape(kind, x, x0, height, width))
        .sum()
}

/// Samples the broadened spectrum at `n_points` evenly spaced points across
/// `[x_min, x_max]`, summing every peak's line shape.
pub fn sample_spectrum(
    peaks: &[(f64, f64)],
    kind: BroadeningKind,
    width: f64,
    x_min: f64,
    x_max: f64,
    n_points: usize,
) -> Vec<(f64, f64)> {
    let n_points = n_points.max(2);
    (0..n_points)
        .map(|i| {
            let x = x_min + (x_max - x_min) * i as f64 / (n_points - 1) as f64;
            (x, intensity_at(peaks, kind, width, x))
        })
        .collect()
}

/// Formats a broadened spectrum as two space-separated columns (position,
/// intensity), one point per line, sampled at a fixed resolution rather than
/// the display's own point count -- so the exported `.dat` does not depend on
/// how large the plot window happened to be when it was exported.
pub fn format_spectrum_dat(
    peaks: &[(f64, f64)],
    kind: BroadeningKind,
    width: f64,
    x_min: f64,
    x_max: f64,
    resolution: f64,
) -> String {
    let resolution = resolution.max(1.0e-9);
    let n_points = (((x_max - x_min) / resolution).round() as usize).max(1) + 1;
    let curve = sample_spectrum(peaks, kind, width, x_min, x_max, n_points);
    let mut out = String::new();
    for (x, y) in curve {
        out.push_str(&format!("{x:.4} {y:.8}\n"));
    }
    out
}

/// Exports the theoretical stick spectrum as-is: one line per predicted peak,
/// the same pairs the list and the stick plot show -- no resampling grid,
/// since there is no continuous curve to resample.
pub fn format_spectrum_sticks(peaks: &[(f64, f64)]) -> String {
    let mut out = String::new();
    for &(x, y) in peaks {
        out.push_str(&format!("{x:.4} {y:.8}\n"));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_spectrum_dat_uses_space_separated_two_columns_at_the_requested_resolution() {
        let text = format_spectrum_dat(
            &[(500.0, 10.0)],
            BroadeningKind::Gaussian,
            5.0,
            0.0,
            10.0,
            0.5,
        );
        let lines: Vec<&str> = text.lines().collect();
        // 0.0..=10.0 at 0.5 resolution: 21 points.
        assert_eq!(lines.len(), 21);
        for line in &lines {
            let cols: Vec<&str> = line.split(' ').collect();
            assert_eq!(cols.len(), 2, "expected exactly two space-separated columns: {line:?}");
            assert!(cols[0].parse::<f64>().is_ok());
            assert!(cols[1].parse::<f64>().is_ok());
        }
        // Points land exactly on the 0.5 grid: 0.0, 0.5, 1.0, ...
        let x0: f64 = lines[0].split(' ').next().unwrap().parse().unwrap();
        let x1: f64 = lines[1].split(' ').next().unwrap().parse().unwrap();
        assert!((x0 - 0.0).abs() < 1e-9);
        assert!((x1 - 0.5).abs() < 1e-9);
    }

    /// The UV-Vis export is specified at 0.1 nm, a finer grid than the IR
    /// export ever asks for -- fine enough that rounding the point count could
    /// drift the last point off the requested end of the range.
    #[test]
    fn a_tenth_of_a_nanometre_grid_lands_exactly_on_its_endpoints() {
        let text = format_spectrum_dat(
            &[(400.0, 1.0)],
            BroadeningKind::Gaussian,
            20.0,
            200.0,
            700.0,
            0.1,
        );
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 5001, "200..=700 at 0.1 nm");
        let first: f64 = lines[0].split(' ').next().unwrap().parse().unwrap();
        let last: f64 = lines[5000].split(' ').next().unwrap().parse().unwrap();
        assert!((first - 200.0).abs() < 1e-6);
        assert!((last - 700.0).abs() < 1e-6);
    }

    #[test]
    fn format_spectrum_dat_reflects_the_chosen_lineshape() {
        let gaussian = format_spectrum_dat(
            &[(500.0, 10.0)], BroadeningKind::Gaussian, 5.0, 495.0, 505.0, 1.0,
        );
        let lorentzian = format_spectrum_dat(
            &[(500.0, 10.0)], BroadeningKind::Lorentzian, 5.0, 495.0, 505.0, 1.0,
        );
        assert_ne!(
            gaussian, lorentzian,
            "different line shapes must export different curves"
        );
    }

    #[test]
    fn format_spectrum_sticks_prints_the_raw_predicted_pairs_unresampled() {
        let text = format_spectrum_sticks(&[
            (259.56, 413.48),
            (1108.86, 0.00066),
            (3534.68, 101.16),
        ]);
        let lines: Vec<&str> = text.lines().collect();
        // One line per peak, not a resampled grid.
        assert_eq!(lines.len(), 3);
        for line in &lines {
            let cols: Vec<&str> = line.split(' ').collect();
            assert_eq!(cols.len(), 2);
        }
        assert!(lines[0].starts_with("259.5"));
        assert!(lines[2].starts_with("3534.6"));
    }

    #[test]
    fn format_spectrum_sticks_of_no_peaks_is_empty() {
        assert_eq!(format_spectrum_sticks(&[]), "");
    }

    #[test]
    fn sticks_lineshape_is_never_used_for_a_continuous_curve() {
        // Sticks are handled by a dedicated code path that never samples a
        // curve at all; `lineshape` itself just needs to stay exhaustive and
        // inert for this variant rather than panic.
        assert_eq!(lineshape(BroadeningKind::Sticks, 500.0, 500.0, 10.0, 5.0), 0.0);
    }

    #[test]
    fn is_sticks_identifies_only_the_sticks_variant() {
        assert!(BroadeningKind::Sticks.is_sticks());
        assert!(!BroadeningKind::Gaussian.is_sticks());
        assert!(!BroadeningKind::Lorentzian.is_sticks());
        assert!(!BroadeningKind::Voigt.is_sticks());
    }

    #[test]
    fn every_variant_is_offered_exactly_once() {
        // The selectors iterate `ALL`; a missing variant would be unreachable
        // from the UI.
        assert_eq!(BroadeningKind::ALL.len(), 4);
        for kind in BroadeningKind::ALL {
            assert_eq!(BroadeningKind::ALL.iter().filter(|&&k| k == kind).count(), 1);
        }
    }

    #[test]
    fn gaussian_peak_matches_its_height_at_the_centre() {
        let y = lineshape(BroadeningKind::Gaussian, 500.0, 500.0, 10.0, 5.0);
        assert!((y - 10.0).abs() < 1.0e-9);
    }

    #[test]
    fn lorentzian_peak_matches_its_height_at_the_centre() {
        let y = lineshape(BroadeningKind::Lorentzian, 500.0, 500.0, 10.0, 5.0);
        assert!((y - 10.0).abs() < 1.0e-9);
    }

    #[test]
    fn spectrum_sampling_places_a_visible_peak_near_each_position() {
        let peaks = [(500.0, 10.0), (1500.0, 20.0)];
        let curve = sample_spectrum(&peaks, BroadeningKind::Gaussian, 10.0, 0.0, 2000.0, 400);
        let (peak_x, peak_y) = curve.iter().cloned().fold((0.0, f64::MIN), |acc, p| {
            if p.1 > acc.1 { p } else { acc }
        });
        assert!((peak_x - 1500.0).abs() < 10.0, "peak at {peak_x}, expected near 1500");
        assert!(peak_y > 15.0);
    }

    #[test]
    fn spectrum_of_no_peaks_is_flat_zero() {
        let curve = sample_spectrum(&[], BroadeningKind::Gaussian, 10.0, 0.0, 100.0, 5);
        assert!(curve.iter().all(|&(_, y)| y == 0.0));
    }

    #[test]
    fn intensity_at_agrees_with_the_sampled_curve() {
        // The plot uses `intensity_at` for its peak markers and
        // `sample_spectrum` for the curve; if they disagreed the markers
        // would float off the line they mark.
        let peaks = [(400.0, 1.0), (410.0, 2.0)];
        let curve = sample_spectrum(&peaks, BroadeningKind::Voigt, 15.0, 300.0, 500.0, 201);
        for &(x, y) in &curve {
            let direct = intensity_at(&peaks, BroadeningKind::Voigt, 15.0, x);
            assert!((y - direct).abs() < 1.0e-12);
        }
    }
}
