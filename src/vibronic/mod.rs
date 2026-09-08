//! Vibronic analysis: what an excited state does to the ground state's
//! vibrations.
//!
//! Two analyses, sharing one idea -- express an excited-state quantity in the
//! ground state's normal coordinates.
//!
//! * [`franck_condon_activity`] projects an excited-state *gradient*, taken at
//!   the ground-state geometry, onto the ground-state modes. That gives the
//!   displacement each mode would undergo, its Huang-Rhys factor, and its
//!   resonance-Raman intensity in the short-time approximation. It needs no
//!   excited-state optimisation at all, which matters: an excited-state
//!   minimum is often hard to converge, and this analysis is untouched by
//!   that.
//!
//! * [`duschinsky`] compares two full mode sets, saying how the modes of one
//!   state are built from the modes of the other. That is the tool for asking
//!   which ground-state mode an excited-state band descends from, especially
//!   where a near-degenerate pair makes the answer ambiguous by frequency
//!   alone.
//!
//! Conventions throughout: atomic units internally, masses in electron
//! masses, normal-mode columns Cartesian and normalised so that
//! `sum_i m_i L_ik^2 = 1`, which is what `normalmode` produces.

use ndarray::Array2;

use crate::rmsd::kabsch;
use bevy::math::Vec3;

/// Hartree per cm^-1: converts a wavenumber to an atomic-unit frequency.
pub const CM1_TO_AU: f64 = 1.0 / 219_474.631_363_20;
/// Electron masses per amu.
pub const AMU_TO_AU: f64 = 1822.888_486_209;

/// What one ground-state mode does when the molecule is put on an excited
/// state's surface at the Franck-Condon geometry.
#[derive(Debug, Clone, PartialEq)]
pub struct ModeActivity {
    /// Index into the ground state's mode list.
    pub index: usize,
    pub frequency_cm1: f64,
    /// The excited-state gradient projected onto this mode, in atomic units.
    pub gradient_au: f64,
    /// Displacement of the excited-state minimum along this mode, in
    /// dimensionless normal coordinates: `-g_k / w_k^{3/2}`.
    ///
    /// Signed, because the direction matters when comparing two states.
    pub displacement: f64,
    /// Huang-Rhys factor `S_k = delta_k^2 / 2`: the mean number of quanta of
    /// this mode excited in a vertical transition.
    pub huang_rhys: f64,
    /// This mode's share of the reorganisation energy, `S_k * w_k`, in cm^-1.
    pub reorganisation_cm1: f64,
    /// Resonance-Raman intensity in the short-time approximation, before
    /// normalisation: proportional to `w_k^2 * delta_k^2`.
    pub rr_intensity: f64,
}

/// The whole Franck-Condon picture for one excited state.
#[derive(Debug, Clone, PartialEq)]
pub struct FrankCondonAnalysis {
    pub modes: Vec<ModeActivity>,
    /// Total reorganisation energy in cm^-1, the sum over modes.
    pub reorganisation_cm1: f64,
    /// The part of the gradient that no vibration accounts for, as a fraction
    /// of its norm.
    ///
    /// It should be small. A large value means the gradient points somewhere
    /// the mode list cannot represent -- translations and rotations that were
    /// projected out of the Hessian, or a mismatch between the gradient's
    /// geometry and the Hessian's.
    pub unaccounted_fraction: f64,
}

impl FrankCondonAnalysis {
    /// Modes ordered by resonance-Raman intensity, strongest first: the
    /// spectrum this state should show.
    pub fn by_rr_intensity(&self) -> Vec<&ModeActivity> {
        let mut v: Vec<&ModeActivity> = self.modes.iter().collect();
        v.sort_by(|a, b| b.rr_intensity.total_cmp(&a.rr_intensity));
        v
    }

    /// Intensities rescaled so the strongest is 1.0, which is the only form
    /// comparable with a measured spectrum.
    pub fn normalised_rr(&self) -> Vec<(f64, f64)> {
        let peak = self
            .modes
            .iter()
            .map(|m| m.rr_intensity)
            .fold(0.0_f64, f64::max);
        if peak <= 0.0 {
            return self.modes.iter().map(|m| (m.frequency_cm1, 0.0)).collect();
        }
        self.modes
            .iter()
            .map(|m| (m.frequency_cm1, m.rr_intensity / peak))
            .collect()
    }
}

/// Projects an excited-state gradient onto the ground-state normal modes.
///
/// `modes` are the ground state's, Cartesian columns normalised so that
/// `sum_i m_i L_ik^2 = 1`; `gradient_cart` is the excited state's Cartesian
/// gradient in Eh/bohr **at the ground-state geometry**; `which` selects the
/// genuine vibrations, since the projected translations and rotations carry no
/// meaning here.
///
/// The mass weighting cancels. In mass-weighted coordinates the projection is
/// `sum_i (sqrt(m_i) L_ik)(g_i / sqrt(m_i))`, so it reduces to a plain dot
/// product of the stored Cartesian column with the raw gradient.
///
/// Every quantity follows from the vertical-gradient (independent-mode,
/// displaced-harmonic-oscillator) picture: the excited-state surface is taken
/// to have the ground state's curvature, displaced by `-g_k / w_k^2`.
pub fn franck_condon_activity(
    modes: &Array2<f64>,
    frequencies_cm1: &[f64],
    gradient_cart: &[f64],
    masses_amu: &[f64],
    which: &[usize],
) -> Result<FrankCondonAnalysis, String> {
    let (ndim, ncol) = modes.dim();
    if gradient_cart.len() != ndim {
        return Err(format!(
            "the gradient has {} components but the modes have {ndim}",
            gradient_cart.len()
        ));
    }
    if frequencies_cm1.len() != ncol {
        return Err(format!(
            "{} frequencies for {ncol} modes",
            frequencies_cm1.len()
        ));
    }

    let mut out = Vec::with_capacity(which.len());
    let mut recovered = 0.0_f64;
    for &k in which {
        if k >= ncol {
            return Err(format!("mode index {k} is out of range"));
        }
        let omega_cm1 = frequencies_cm1[k];
        // An imaginary mode has no displaced minimum and no Franck-Condon
        // progression; reporting one would be inventing a number.
        if omega_cm1 <= 0.0 {
            continue;
        }
        let g_k: f64 = (0..ndim).map(|i| modes[(i, k)] * gradient_cart[i]).sum();
        recovered += g_k * g_k;

        let omega_au = omega_cm1 * CM1_TO_AU;
        let displacement = -g_k / omega_au.powf(1.5);
        let huang_rhys = 0.5 * displacement * displacement;
        out.push(ModeActivity {
            index: k,
            frequency_cm1: omega_cm1,
            gradient_au: g_k,
            displacement,
            huang_rhys,
            reorganisation_cm1: huang_rhys * omega_cm1,
            // Short-time approximation: the intensity of a fundamental goes as
            // the square of the dimensionless gradient, w_k^2 delta_k^2.
            rr_intensity: omega_cm1 * omega_cm1 * displacement * displacement,
        });
    }

    // The comparison must be in the metric the modes are orthonormal in,
    // which is the mass-weighted one: sum_k g_k^2 recovers |M^-1/2 g|^2, not
    // |g|^2. Comparing against the Cartesian norm reports every gradient as
    // entirely unaccounted for, whatever the projection actually did.
    let total: f64 = (0..ndim)
        .map(|i| {
            let m = masses_amu[i / 3] * AMU_TO_AU;
            gradient_cart[i] * gradient_cart[i] / m
        })
        .sum();
    let unaccounted_fraction = if total > 0.0 {
        (1.0 - (recovered / total)).max(0.0).sqrt()
    } else {
        0.0
    };

    let reorganisation_cm1 = out.iter().map(|m| m.reorganisation_cm1).sum();
    Ok(FrankCondonAnalysis {
        modes: out,
        reorganisation_cm1,
        unaccounted_fraction,
    })
}

/// How one state's normal modes are built out of another's.
#[derive(Debug, Clone)]
pub struct Duschinsky {
    /// `J[i][j]`: the amount of ground-state mode `i` in excited-state mode
    /// `j`. Rows are the ground state, columns the excited state.
    pub j: Array2<f64>,
    /// Ground-state mode indices labelling the rows.
    pub rows: Vec<usize>,
    /// Excited-state mode indices labelling the columns.
    pub cols: Vec<usize>,
    pub row_frequencies_cm1: Vec<f64>,
    pub col_frequencies_cm1: Vec<f64>,
}

/// One excited-state mode and the ground-state mode it most resembles.
#[derive(Debug, Clone, PartialEq)]
pub struct ModeCorrespondence {
    pub excited_index: usize,
    pub excited_frequency_cm1: f64,
    pub ground_index: usize,
    pub ground_frequency_cm1: f64,
    /// `J_ij^2` for the best partner: the fraction of this excited-state mode
    /// that is the ground-state one. 1.0 is a pure correspondence.
    pub overlap: f64,
    /// How many ground-state modes it takes to make up this excited-state
    /// mode, as an inverse participation ratio over `J^2`.
    ///
    /// 1.0 means the mode is preserved intact; 3.0 means it is a mixture of
    /// about three. This is the quantitative measure of Duschinsky mixing.
    pub participation: f64,
    pub shift_cm1: f64,
}

impl Duschinsky {
    /// For each excited-state mode, the ground-state mode it most resembles.
    pub fn correspondence(&self) -> Vec<ModeCorrespondence> {
        let (nrow, ncol) = self.j.dim();
        let mut out = Vec::with_capacity(ncol);
        for c in 0..ncol {
            let mut best = 0usize;
            let mut best_sq = -1.0_f64;
            let mut sum_sq = 0.0_f64;
            let mut sum_quart = 0.0_f64;
            for r in 0..nrow {
                let sq = self.j[(r, c)] * self.j[(r, c)];
                sum_sq += sq;
                sum_quart += sq * sq;
                if sq > best_sq {
                    best_sq = sq;
                    best = r;
                }
            }
            // Normalised so a row that is spread over the block still reports
            // a meaningful fraction even when the block is a subset of modes.
            let overlap = if sum_sq > 0.0 { best_sq / sum_sq } else { 0.0 };
            let participation = if sum_quart > 0.0 {
                (sum_sq * sum_sq) / sum_quart
            } else {
                0.0
            };
            out.push(ModeCorrespondence {
                excited_index: self.cols[c],
                excited_frequency_cm1: self.col_frequencies_cm1[c],
                ground_index: self.rows[best],
                ground_frequency_cm1: self.row_frequencies_cm1[best],
                overlap,
                participation,
                shift_cm1: self.col_frequencies_cm1[c] - self.row_frequencies_cm1[best],
            });
        }
        out
    }

    /// The mean participation over all columns: one number for how much this
    /// state scrambles the ground state's modes.
    ///
    /// 1.0 is no mixing at all. Larger means the excited state's vibrations
    /// are combinations of several ground-state ones, which is the situation
    /// that opens intramolecular energy-redistribution channels.
    pub fn mean_participation(&self) -> f64 {
        let c = self.correspondence();
        if c.is_empty() {
            return 0.0;
        }
        c.iter().map(|m| m.participation).sum::<f64>() / c.len() as f64
    }
}

/// The Duschinsky matrix between two mode sets, `J = L_A^t M L_B`.
///
/// Both mode sets are Cartesian columns normalised against the mass metric, so
/// the overlap that makes them orthonormal is `sum_i m_i L^A_ik L^B_il`, and
/// the mass matrix has to appear explicitly.
///
/// The two geometries are Eckart-aligned first: a normal mode is defined in
/// the molecule's own frame, and two optimisations have no reason to leave
/// their structures in the same one. Without alignment the rotation between
/// frames leaks into `J` and looks like mode mixing.
#[allow(clippy::too_many_arguments)]
pub fn duschinsky(
    ground_modes: &Array2<f64>,
    ground_freqs_cm1: &[f64],
    ground_coords: &[Vec3],
    excited_modes: &Array2<f64>,
    excited_freqs_cm1: &[f64],
    excited_coords: &[Vec3],
    masses_amu: &[f64],
    ground_which: &[usize],
    excited_which: &[usize],
) -> Result<Duschinsky, String> {
    let natoms = masses_amu.len();
    let ndim = 3 * natoms;
    if ground_modes.dim().0 != ndim || excited_modes.dim().0 != ndim {
        return Err("the two mode sets describe different numbers of atoms".to_string());
    }
    if ground_coords.len() != natoms || excited_coords.len() != natoms {
        return Err("coordinates do not match the mass list".to_string());
    }

    // Bring the excited-state frame onto the ground-state one, and rotate its
    // modes by the same rotation so they stay in step with their geometry.
    let superposition = kabsch::superpose(excited_coords, ground_coords)?;
    // Only the rotation applies to a mode: a normal mode is a displacement,
    // so the translation part of the superposition must not be added to it.
    let r = &superposition.rotation;
    let rotate = |v: Vec3| {
        let (x, y, z) = (v.x as f64, v.y as f64, v.z as f64);
        Vec3::new(
            (r[0][0] * x + r[0][1] * y + r[0][2] * z) as f32,
            (r[1][0] * x + r[1][1] * y + r[1][2] * z) as f32,
            (r[2][0] * x + r[2][1] * y + r[2][2] * z) as f32,
        )
    };

    let mut j = Array2::<f64>::zeros((ground_which.len(), excited_which.len()));
    for (r, &a) in ground_which.iter().enumerate() {
        for (c, &b) in excited_which.iter().enumerate() {
            let mut acc = 0.0;
            for atom in 0..natoms {
                let m = masses_amu[atom] * AMU_TO_AU;
                let ga = Vec3::new(
                    ground_modes[(3 * atom, a)] as f32,
                    ground_modes[(3 * atom + 1, a)] as f32,
                    ground_modes[(3 * atom + 2, a)] as f32,
                );
                let eb = rotate(Vec3::new(
                    excited_modes[(3 * atom, b)] as f32,
                    excited_modes[(3 * atom + 1, b)] as f32,
                    excited_modes[(3 * atom + 2, b)] as f32,
                ));
                acc += m * (ga.x * eb.x + ga.y * eb.y + ga.z * eb.z) as f64;
            }
            j[(r, c)] = acc;
        }
    }

    Ok(Duschinsky {
        j,
        rows: ground_which.to_vec(),
        cols: excited_which.to_vec(),
        row_frequencies_cm1: ground_which.iter().map(|&i| ground_freqs_cm1[i]).collect(),
        col_frequencies_cm1: excited_which
            .iter()
            .map(|&i| excited_freqs_cm1[i])
            .collect(),
    })
}

#[cfg(test)]
mod feacac3 {
    use super::*;
    use crate::normalmode::{normal_modes_with_projection, EckartMode};
    use crate::qchem_interfaces::{orca_engrad, orca_hess};
    use crate::vibronic::*;
    use std::path::PathBuf;

    fn dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/FeAcac3")
    }

    struct State { freqs: Vec<f64>, modes: ndarray::Array2<f64>, pos: Vec<f64>,
                   masses: Vec<f64>, vibs: Vec<usize>, imag: usize }

    fn load(rel: &str) -> State {
        let text = std::fs::read_to_string(dir().join(rel)).unwrap();
        let h = orca_hess::parse_orca_hess(&text).unwrap();
        let r = normal_modes_with_projection(&h.masses_amu, &h.hessian, false,
                    Some(&h.coords_bohr), EckartMode::VibRot, None).unwrap();
        let f: Vec<f64> = r.frequencies_au.iter().map(|w| w * 219474.63136320).collect();
        let mut vibs = r.positive_indices.clone();
        vibs.sort();
        State { freqs: f, modes: r.modes, pos: h.coords_bohr,
                masses: h.masses_amu, vibs, imag: r.negative_indices.len() }
    }

    fn grad(rel: &str) -> Vec<f64> {
        let t = std::fs::read_to_string(dir().join(rel)).unwrap();
        orca_engrad::parse_orca_engrad(&t).unwrap().gradient_bohr
    }



    fn window(st: &State, lo: f64, hi: f64) -> Vec<usize> {
        let mut v: Vec<usize> = st.vibs.iter().copied()
            .filter(|&i| st.freqs[i] > lo && st.freqs[i] < hi).collect();
        v.sort_by(|&a,&b| st.freqs[a].total_cmp(&st.freqs[b]));
        v
    }

    fn as_vec3(c: &[f64]) -> Vec<bevy::math::Vec3> {
        c.chunks(3).map(|p| bevy::math::Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32)).collect()
    }

    #[test]
    #[ignore]
    fn duschinsky_analysis() {
        let gs = load("GroundState_Opt/TPSSh-def2-TZVP_SCMwater_NoSym.hess");
        for (name, rel) in [
            ("Root 15  (424 nm, exp 400)", "Excited_State_Optimizations/Root_15_424nm_Opt_Freq/Exc_root_OnlyFreq_from_optimized.hess"),
            ("Root 61  (274 nm, exp 270)", "Excited_State_Optimizations/Root_61_274nm_Opt_Freq/Exc_root_opt_continue.hess")] {
            let es = load(rel);
            println!("\n============ {name} ============");
            println!("excited state: {} real modes, {} imaginary", es.vibs.len(), es.imag);

            // Mid-frequency block only: the low region is contaminated by the
            // unconverged soft torsions.
            let gw = window(&gs, 700.0, 1700.0);
            let ew = window(&es, 700.0, 1700.0);
            let d = duschinsky(&gs.modes, &gs.freqs, &as_vec3(&gs.pos),
                               &es.modes, &es.freqs, &as_vec3(&es.pos),
                               &gs.masses, &gw, &ew).unwrap();
            let corr = d.correspondence();
            println!("block: {} GS x {} ES modes in 700-1700 cm-1", gw.len(), ew.len());
            println!("mean participation (1.0 = no mixing): {:.2}", d.mean_participation());
            let pure = corr.iter().filter(|c| c.overlap > 0.7).count();
            println!("modes preserved intact (overlap>0.7): {} of {}", pure, corr.len());

            println!("--- excited-state modes in the 900-1100 window (the 941 band region) ---");
            println!("{:>9} {:>10} {:>9} {:>8} {:>7} {:>7}", "ES freq", "<- GS freq", "shift", "overlap", "partic", "GS idx");
            for c in corr.iter().filter(|c| c.excited_frequency_cm1 > 900.0 && c.excited_frequency_cm1 < 1100.0) {
                println!("{:>9.1} {:>10.1} {:>9.1} {:>8.2} {:>7.1} {:>7}",
                    c.excited_frequency_cm1, c.ground_frequency_cm1, c.shift_cm1,
                    c.overlap, c.participation, c.ground_index);
            }
        }
    }


    #[test]
    #[ignore]
    fn biggest_shifts() {
        let gs = load("GroundState_Opt/TPSSh-def2-TZVP_SCMwater_NoSym.hess");
        for (name, rel) in [
            ("Root 15", "Excited_State_Optimizations/Root_15_424nm_Opt_Freq/Exc_root_OnlyFreq_from_optimized.hess"),
            ("Root 61", "Excited_State_Optimizations/Root_61_274nm_Opt_Freq/Exc_root_opt_continue.hess")] {
            let es = load(rel);
            let gw = window(&gs, 700.0, 1700.0);
            let ew = window(&es, 700.0, 1700.0);
            let d = duschinsky(&gs.modes, &gs.freqs, &as_vec3(&gs.pos),
                               &es.modes, &es.freqs, &as_vec3(&es.pos),
                               &gs.masses, &gw, &ew).unwrap();
            let mut c = d.correspondence();
            c.sort_by(|a,b| b.shift_cm1.total_cmp(&a.shift_cm1));
            println!("\n{name}: largest BLUE shifts (ES higher than its GS parent)");
            for x in c.iter().take(5) {
                println!("   GS {:>7.1} -> ES {:>7.1}   {:+7.1} cm-1   overlap {:.2}",
                    x.ground_frequency_cm1, x.excited_frequency_cm1, x.shift_cm1, x.overlap);
            }
            // what becomes of the nu8 sextet specifically (GS 948-952)
            let nu8: Vec<&ModeCorrespondence> = c.iter()
                .filter(|x| x.ground_frequency_cm1 > 947.0 && x.ground_frequency_cm1 < 953.0).collect();
            let lo = nu8.iter().map(|x| x.excited_frequency_cm1).fold(f64::MAX, f64::min);
            let hi = nu8.iter().map(|x| x.excited_frequency_cm1).fold(0.0f64, f64::max);
            println!("   nu8 descendants: {} modes, {:.1} to {:.1} cm-1 (GS sextet spans 948.1-951.5)", nu8.len(), lo, hi);
        }
    }


    #[test]
    #[ignore]
    fn write_full_analysis() {
        use std::fmt::Write as _;
        let gs = load("GroundState_Opt/TPSSh-def2-TZVP_SCMwater_NoSym.hess");
        let mut o = String::new();
        writeln!(o, "Fe(acac)3 -- vibronic analysis").unwrap();
        writeln!(o, "TPSSh / def2-TZVP(Fe,O) def2-SVP(C,H) / SMD(water) / TDA, ORCA 6").unwrap();
        writeln!(o, "Ground state: {} real modes, {} imaginary\n", gs.vibs.len(), gs.imag).unwrap();

        for (name, gpath, hpath, ev, nm) in [
          ("ROOT 15", "Excited_Engrad/Root_15/Exc_root_engrad.engrad",
           "Excited_State_Optimizations/Root_15_424nm_Opt_Freq/Exc_root_OnlyFreq_from_optimized.hess", 2.918, 424.8),
          ("ROOT 61", "Excited_Engrad/Root_61/Exc_root_engrad.engrad",
           "Excited_State_Optimizations/Root_61_274nm_Opt_Freq/Exc_root_opt_continue.hess", 4.475, 277.0)] {

            let g = grad(gpath);
            let a = franck_condon_activity(&gs.modes, &gs.freqs, &g, &gs.masses, &gs.vibs).unwrap();
            let peak = a.modes.iter().map(|m| m.rr_intensity).fold(0.0f64,f64::max);
            writeln!(o, "\n{:=<100}", "").unwrap();
            writeln!(o, "{name}   E = {ev} eV = {nm} nm").unwrap();
            writeln!(o, "{:=<100}\n", "").unwrap();
            writeln!(o, "A. FRANCK-CONDON ACTIVITY  (excited-state gradient projected on ground-state modes)").unwrap();
            writeln!(o, "   total reorganisation energy = {:.1} cm-1 = {:.4} eV", a.reorganisation_cm1, a.reorganisation_cm1/8065.544).unwrap();
            writeln!(o, "   gradient unaccounted for by the vibrations = {:.3} %\n", 100.0*a.unaccounted_fraction).unwrap();
            writeln!(o, "{:>5} {:>10} {:>14} {:>12} {:>12} {:>12} {:>10}",
                "mode","freq/cm-1","grad/au","displ","S_k","lambda/cm-1","RR(rel)").unwrap();
            let mut byfreq = a.modes.clone();
            byfreq.sort_by(|x,y| x.frequency_cm1.total_cmp(&y.frequency_cm1));
            for m in &byfreq {
                writeln!(o, "{:>5} {:>10.2} {:>14.6e} {:>12.4} {:>12.5} {:>12.2} {:>10.4}",
                  m.index, m.frequency_cm1, m.gradient_au, m.displacement, m.huang_rhys,
                  m.reorganisation_cm1, m.rr_intensity/peak).unwrap();
            }

            let es = load(hpath);
            let gw = window(&gs, 200.0, 1800.0);
            let ew = window(&es, 200.0, 1800.0);
            let d = duschinsky(&gs.modes, &gs.freqs, &as_vec3(&gs.pos),
                               &es.modes, &es.freqs, &as_vec3(&es.pos),
                               &gs.masses, &gw, &ew).unwrap();
            let corr = d.correspondence();
            writeln!(o, "\nB. DUSCHINSKY MODE CORRESPONDENCE  (200-1800 cm-1 block)").unwrap();
            writeln!(o, "   excited state: {} real modes, {} imaginary", es.vibs.len(), es.imag).unwrap();
            writeln!(o, "   block: {} ground-state x {} excited-state modes", gw.len(), ew.len()).unwrap();
            writeln!(o, "   mean participation = {:.3}   (1.0 = modes preserved intact)", d.mean_participation()).unwrap();
            writeln!(o, "   modes with overlap > 0.7 : {} of {}\n",
                corr.iter().filter(|c| c.overlap>0.7).count(), corr.len()).unwrap();
            writeln!(o, "{:>10} {:>12} {:>10} {:>10} {:>14} {:>8}",
                "ES/cm-1","<-GS/cm-1","shift","overlap","participation","GSidx").unwrap();
            for c in &corr {
                writeln!(o, "{:>10.2} {:>12.2} {:>10.2} {:>10.4} {:>14.2} {:>8}",
                  c.excited_frequency_cm1, c.ground_frequency_cm1, c.shift_cm1,
                  c.overlap, c.participation, c.ground_index).unwrap();
            }
        }
        let path = dir().join("vibronic_analysis.txt");
        std::fs::write(&path, o).unwrap();
        println!("wrote {}", path.display());
    }


    #[test]
    #[ignore]
    fn alt_mode_character() {
        let gs = load("GroundState_Opt/TPSSh-def2-TZVP_SCMwater_NoSym.hess");
        let text = std::fs::read_to_string(dir().join("GroundState_Opt/TPSSh-def2-TZVP_SCMwater_NoSym.hess")).unwrap();
        let h = crate::qchem_interfaces::orca_hess::parse_orca_hess(&text).unwrap();
        let n = h.atoms.len();
        for target in [1031.4, 1040.2, 1026.3, 1295.5, 1035.9, 661.1, 277.6] {
            let k = gs.vibs.iter().copied()
                .min_by(|&a,&b| (gs.freqs[a]-target).abs().total_cmp(&(gs.freqs[b]-target).abs())).unwrap();
            let mut share: std::collections::HashMap<String,f64> = Default::default();
            let mut tot=0.0;
            for a in 0..n {
                let m=h.masses_amu[a];
                let d: f64=(0..3).map(|c| gs.modes[(3*a+c,k)].powi(2)).sum();
                *share.entry(h.atoms[a].clone()).or_default()+=m*d; tot+=m*d;
            }
            let mut v: Vec<(String,f64)>=share.into_iter().map(|(e,x)|(e,100.0*x/tot)).collect();
            v.sort_by(|a,b| b.1.total_cmp(&a.1));
            let sv: Vec<String>=v.iter().filter(|(_,pc)| *pc>3.0).map(|(e,pc)| format!("{e} {pc:.0}%")).collect();
            println!("{:>8.1} : {}", gs.freqs[k], sv.join(", "));
        }
        // how many GS modes lie in each experimental window
        for (lab,lo,hi) in [("780-1010",780.0,1010.0),("1010-1060",1010.0,1060.0),("1400-1520",1400.0,1520.0),("1500-1620",1500.0,1620.0)] {
            let c=gs.vibs.iter().filter(|&&i| gs.freqs[i]>lo && gs.freqs[i]<hi).count();
            println!("window {lab}: {c} modes");
        }
    }


    #[test]
    #[ignore]
    fn cluster_fine_structure() {
        let gs = load("GroundState_Opt/TPSSh-def2-TZVP_SCMwater_NoSym.hess");
        let mut f: Vec<f64> = gs.vibs.iter().map(|&i| gs.freqs[i]).collect();
        f.sort_by(|a,b| a.total_cmp(b));
        // group into clusters separated by >5 cm-1
        let mut groups: Vec<Vec<f64>> = vec![];
        for v in f {
            if let Some(g) = groups.last_mut() {
                if v - g[g.len()-1] < 5.0 { g.push(v); continue; }
            }
            groups.push(vec![v]);
        }
        println!("GROUND-STATE CLUSTERS in the FSRS window (200-1700 cm-1)");
        for g in groups.iter().filter(|g| g[0]>200.0 && g[0]<1700.0 && g.len()>1) {
            let span = g[g.len()-1]-g[0];
            if span < 0.05 && g.len() < 3 { continue; }
            let vals: Vec<String> = g.iter().map(|v| format!("{v:.1}")).collect();
            println!("   n={:2}  span {:5.1}   {}", g.len(), span, vals.join(" "));
        }
    }


    #[test]
    #[ignore]
    fn symmetry_breaking_evidence() {
        use crate::qchem_interfaces::orca_hess;
        let paths = [
            ("Ground state", "GroundState_Opt/TPSSh-def2-TZVP_SCMwater_NoSym.hess"),
            ("Root 15", "Excited_State_Optimizations/Root_15_424nm_Opt_Freq/Exc_root_OnlyFreq_from_optimized.hess"),
            ("Root 61", "Excited_State_Optimizations/Root_61_274nm_Opt_Freq/Exc_root_opt_continue.hess")];

        println!("=== TEST 1: Fe-O bond lengths (bohr) -- is the FeO6 core still symmetric? ===");
        for (name, rel) in paths {
            let t = std::fs::read_to_string(dir().join(rel)).unwrap();
            let h = orca_hess::parse_orca_hess(&t).unwrap();
            let fe = (0..h.atoms.len()).find(|&i| h.atoms[i]=="Fe").unwrap();
            let mut d: Vec<f64> = (0..h.atoms.len()).filter(|&i| h.atoms[i]=="O")
                .map(|i| (0..3).map(|k| (h.coords_bohr[3*i+k]-h.coords_bohr[3*fe+k]).powi(2)).sum::<f64>().sqrt())
                .collect();
            d.sort_by(|a,b| a.total_cmp(b));
            let mean: f64 = d.iter().sum::<f64>()/d.len() as f64;
            let spread = d[d.len()-1]-d[0];
            let vals: Vec<String> = d.iter().map(|v| format!("{v:.4}")).collect();
            println!("  {:<13} {}   mean {:.4}  spread {:.4} bohr ({:.3} A)",
                     name, vals.join(" "), mean, spread, spread*0.529177);
        }

        println!("\n=== TEST 2: does every near-degenerate cluster widen, or only nu8? ===");
        let gs = load("GroundState_Opt/TPSSh-def2-TZVP_SCMwater_NoSym.hess");
        for (name, rel) in [("Root 15", paths[1].1), ("Root 61", paths[2].1)] {
            let es = load(rel);
            let gw = window(&gs, 200.0, 1700.0);
            let ew = window(&es, 200.0, 1700.0);
            let d = duschinsky(&gs.modes, &gs.freqs, &as_vec3(&gs.pos),
                               &es.modes, &es.freqs, &as_vec3(&es.pos),
                               &gs.masses, &gw, &ew).unwrap();
            let corr = d.correspondence();
            println!("  --- {name}");
            // clusters of GS modes within 5 cm-1
            let mut f: Vec<f64> = gw.iter().map(|&i| gs.freqs[i]).collect();
            f.sort_by(|a,b| a.total_cmp(b));
            let mut groups: Vec<Vec<f64>> = vec![];
            for v in f { if let Some(g)=groups.last_mut() { if v-g[g.len()-1]<5.0 { g.push(v); continue; } } groups.push(vec![v]); }
            println!("     {:>16} {:>8} {:>8} {:>7}", "GS cluster", "GS span", "ES span", "factor");
            for g in groups.iter().filter(|g| g.len()>=3) {
                let (lo,hi)=(g[0],g[g.len()-1]);
                let ds: Vec<f64> = corr.iter()
                    .filter(|c| c.ground_frequency_cm1>=lo-0.05 && c.ground_frequency_cm1<=hi+0.05)
                    .map(|c| c.excited_frequency_cm1).collect();
                if ds.len()<2 { continue; }
                let es_span = ds.iter().cloned().fold(f64::MIN,f64::max) - ds.iter().cloned().fold(f64::MAX,f64::min);
                let gs_span = (hi-lo).max(0.1);
                println!("     {:>7.1}-{:<7.1} {:>8.1} {:>8.1} {:>7.1}x", lo, hi, hi-lo, es_span, es_span/gs_span);
            }
        }
    }

    #[test]
    #[ignore]
    fn mode_character() {
        let gs = load("GroundState_Opt/TPSSh-def2-TZVP_SCMwater_NoSym.hess");
        let text = std::fs::read_to_string(dir().join("GroundState_Opt/TPSSh-def2-TZVP_SCMwater_NoSym.hess")).unwrap();
        let h = crate::qchem_interfaces::orca_hess::parse_orca_hess(&text).unwrap();
        let n = h.atoms.len();
        for target in [1559.5, 435.0, 200.6, 948.1, 811.1, 1443.9] {
            let k = gs.vibs.iter().copied()
                .min_by(|&a,&b| (gs.freqs[a]-target).abs().total_cmp(&(gs.freqs[b]-target).abs())).unwrap();
            // per-element share of the mass-weighted amplitude
            let mut share: std::collections::HashMap<String,f64> = Default::default();
            let mut tot=0.0;
            for a in 0..n {
                let m = h.masses_amu[a];
                let d: f64 = (0..3).map(|c| gs.modes[(3*a+c,k)].powi(2)).sum();
                *share.entry(h.atoms[a].clone()).or_default() += m*d; tot += m*d;
            }
            let mut v: Vec<(String,f64)> = share.into_iter().map(|(e,x)| (e, 100.0*x/tot)).collect();
            v.sort_by(|a,b| b.1.total_cmp(&a.1));
            let s: Vec<String> = v.iter().filter(|(_,p)| *p>3.0).map(|(e,p)| format!("{e} {p:.0}%")).collect();
            println!("{:>8.1} cm-1 : {}", gs.freqs[k], s.join(", "));
        }
    }

    #[test]
    #[ignore]
    fn fc_analysis() {
        let gs = load("GroundState_Opt/TPSSh-def2-TZVP_SCMwater_NoSym.hess");
        println!("GROUND STATE: {} modes, {} imaginary", gs.vibs.len(), gs.imag);
        let mut sorted: Vec<f64> = gs.vibs.iter().map(|&i| gs.freqs[i]).collect();
        sorted.sort_by(|a,b| a.total_cmp(b));
        println!("  lowest 6 : {:?}", sorted.iter().take(6).map(|v| format!("{v:.1}")).collect::<Vec<_>>());
        println!("  near 800-1000: {:?}", sorted.iter().filter(|v| **v>780.0 && **v<1010.0)
                 .map(|v| format!("{v:.1}")).collect::<Vec<_>>());
        println!("  near 1050-1120: {:?}", sorted.iter().filter(|v| **v>1050.0 && **v<1120.0)
                 .map(|v| format!("{v:.1}")).collect::<Vec<_>>());
        println!("  near 1400-1520: {:?}", sorted.iter().filter(|v| **v>1400.0 && **v<1520.0)
                 .map(|v| format!("{v:.1}")).collect::<Vec<_>>());

        for (name, gpath) in [("Root 15  (424 nm, exp 400)", "Excited_Engrad/Root_15/Exc_root_engrad.engrad"),
                              ("Root 61  (274 nm, exp 270)", "Excited_Engrad/Root_61/Exc_root_engrad.engrad")] {
            let g = grad(gpath);
            let a = franck_condon_activity(&gs.modes, &gs.freqs, &g, &gs.masses, &gs.vibs).unwrap();
            println!("\n================ {name} ================");
            println!("total reorganisation energy = {:.0} cm-1 ({:.3} eV)",
                     a.reorganisation_cm1, a.reorganisation_cm1/8065.544);
            println!("gradient not accounted for by vibrations: {:.2}%", 100.0*a.unaccounted_fraction);
            println!("{:>8} {:>10} {:>12} {:>10} {:>10}", "freq", "S_k", "displacement", "lambda", "RR(rel)");
            let peak = a.modes.iter().map(|m| m.rr_intensity).fold(0.0f64,f64::max);
            for m in a.by_rr_intensity().iter().take(12) {
                println!("{:>8.1} {:>10.4} {:>12.3} {:>10.1} {:>10.3}",
                    m.frequency_cm1, m.huang_rhys, m.displacement, m.reorganisation_cm1,
                    m.rr_intensity/peak);
            }
            // what happens in the experimental windows
            for (lab, lo, hi) in [("nu8  ~878/896", 780.0, 1010.0), ("nu9  ~1085", 1050.0, 1120.0), ("nu10 ~1456", 1400.0, 1520.0)] {
                let tot: f64 = a.modes.iter().filter(|m| m.frequency_cm1>lo && m.frequency_cm1<hi)
                    .map(|m| m.rr_intensity).sum();
                let s: f64 = a.modes.iter().filter(|m| m.frequency_cm1>lo && m.frequency_cm1<hi)
                    .map(|m| m.huang_rhys).sum();
                println!("  {lab}: summed RR = {:.3} of peak, summed S = {:.4}", tot/peak, s);
            }
        }
    }
}
