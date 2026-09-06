//! Program-agnostic spectrum simulation and plotting.
//!
//! Both the IR spectrum (from an xTB Hessian) and the UV-Vis absorption
//! spectrum (from an ORCA TD-DFT output) are the same picture: a list of
//! `(position, intensity)` peaks, optionally broadened by a line shape, drawn
//! on a linear axis. Only the units and the wording differ, so the maths and
//! the painter live here once instead of once per spectroscopy.
//!
//! A "peak" is deliberately just `(f64, f64)`. Wavenumbers with km/mol
//! intensities and wavelengths with oscillator strengths are the same shape of
//! data, and pretending otherwise would mean two copies of every line shape.

pub mod broadening;
pub mod plot;
