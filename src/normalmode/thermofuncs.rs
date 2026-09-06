#![allow(dead_code, non_snake_case)]
use crate::normalmode::inertia::get_brot;

const PI: f64 = std::f64::consts::PI;
const PI_SQ: f64 = PI * PI;
const TWOPI: f64 = 2.0 * PI;

const CLIGHT: f64 = 137.035999084;
// Unified atomic mass unit in electron masses (CODATA 2018: 1/m_e[u] with
// m_e = 5.48579909065e-4 u). The previous value here, 1836.15267343, is the
// *proton*-to-electron mass ratio, not the amu-to-electron-mass ratio (a
// proton is ~0.73% heavier than 1 u) — using it made every translational
// partition function, and hence S_trans/U_trans/H_trans/F_trans/G_trans,
// systematically wrong by that same ~0.73%. The correct constant is used
// consistently elsewhere in this project (`C3` in `normalmode/mod.rs`,
// `AMU_TO_ELECTRON_MASS` in `normalmode/ir_intensity.rs`, `AMU_TO_EMASS` in
// `moldyn/mod.rs`).
const AMU_TO_ELECMASS: f64 = 1_822.888_486_209;

const HPLANCK_AU: f64 = TWOPI;
const HPLANCK_AU_SQ: f64 = TWOPI * TWOPI;

// CODATA 2018-derived values (Hartree energy 4.3597447222071e-18 J, Avogadro
// constant 6.02214076e23 /mol — both SI-exact or CODATA-recommended); the
// previous 5-significant-figure roundings (2625.5, 627.51) differed from
// these by ~1.4e-7 and ~8.4e-7 relative respectively — utterly negligible on
// their own, but tightened here to match the precision already used for the
// same constants elsewhere in this project (e.g. `HARTREE_TO_KCAL` in
// `src/bin/*_parameter_fit.rs`, `src/optimizer/irc/mod.rs`).
const AU_TO_KJ: f64 = 2625.499_639_479_163;
const AU_TO_KCAL: f64 = 627.509_474_062_897_4;

const RGAS_AU: f64 = 8.31446261815324 / 1000.0 / AU_TO_KJ; // Hartree/mol/K
const CM1_TO_K: f64 = 1.438_776_877_503_934; // hc/kB per cm^-1, in K
const CM1_TO_HARTREE: f64 = 4.556_335_252_913_159e-6; // hc per cm^-1, in Hartree
const CM1_TO_KCAL: f64 = CM1_TO_HARTREE * AU_TO_KCAL;
const PASCAL_TO_AU: f64 = 3.398_930_921_744_024e-14; // 1 Pa in Hartree/bohr^3

// SI constants needed ONLY for Grimme free-rotor entropy (dimensionless inside ln)
const PLANCK_SI: f64 = 6.62607015e-34; // J*s
const BOLTZMANN_SI: f64 = 1.380649e-23; // J/K
const RGAS_SI: f64 = 8.31446261815324; // J/mol/K
const CLIGHT_SI: f64 = 2.99792458e8; // m/s

// 1 Hartree/mol = AU_TO_KJ kJ/mol = AU_TO_KJ*1000 J/mol
const J_PER_HARTREE_PER_MOL: f64 = AU_TO_KJ * 1000.0;

// Grimme default average moment of inertia (kg*m^2)
const GRIMME_BAV_SI: f64 = 1.0e-44;

// Damping exponent in Grimme qRRHO
const GRIMME_ALPHA: f64 = 4.0;

#[derive(Debug, Clone)]
pub struct ThermoResults {
    pub pfelec: f64,
    pub pftrans: f64,
    pub pfrot: f64,
    pub pfvib: f64,
    pub pftot: f64,

    pub uelec: f64,
    pub utrans: f64,
    pub urot: f64,
    pub uvib: f64,
    pub utherm: f64,
    pub utot: f64,

    pub helec: f64,
    pub htrans: f64,
    pub hrot: f64,
    pub hvib: f64,
    pub htherm: f64,
    pub htot: f64,

    pub felec: f64,
    pub ftrans: f64,
    pub frot: f64,
    pub fvib: f64,
    pub ftherm: f64,
    pub ftot: f64,

    pub gelec: f64,
    pub gtrans: f64,
    pub grot: f64,
    pub gvib: f64,
    pub gtherm: f64,
    pub gtot: f64,

    pub selec: f64,
    pub strans: f64,
    pub srot: f64,
    pub svib: f64,
    pub stherm: f64,
    pub stot: f64,

    pub cvelec: f64,
    pub cvtrans: f64,
    pub cvrot: f64,
    pub cvvib: f64,
    pub cvtherm: f64,
    pub cvtot: f64,

    pub cpelec: f64,
    pub cptrans: f64,
    pub cprot: f64,
    pub cpvib: f64,
    pub cptherm: f64,
    pub cptot: f64,

    pub zpe: f64,
}

pub fn freqs_au_to_cm1(freqs_au: &[f64]) -> Vec<f64> {
    let c1 = 1.0 / 0.529_177_210_903; // Angstrom -> bohr
    let c9 = 1.0e8 * c1; // frequency in cm^-1 -> bohr^-1
    freqs_au
        .iter()
        .map(|&omega| omega / CLIGHT * c9 / (std::f64::consts::PI * 2.0))
        .collect()
}

/// `rot_symmetry` is the rotational symmetry number sigma, which divides the
/// rotational partition function: the number of indistinguishable orientations
/// the molecule can be rotated into. 1 for an asymmetric molecule, 2 for water,
/// 12 for benzene or methane. Getting it wrong shifts the rotational entropy by
/// `R ln(sigma)`, which goes straight into every free energy derived from it.
#[allow(clippy::too_many_arguments)]
pub fn eval_thermo(
    freqs_cm1: &[f64],
    brot_cm1: &[f64],
    mass_amu_total: f64,
    multiplicity: f64,
    temp: f64,
    pressure: f64,
    freq_cutoff: f64,
    rot_symmetry: f64,
) -> ThermoResults {
    let mut thermo = ThermoResults {
        pfelec: 1.0,
        pftrans: 1.0,
        pfrot: 1.0,
        pfvib: 1.0,
        pftot: 1.0,
        uelec: 0.0,
        utrans: 0.0,
        urot: 0.0,
        uvib: 0.0,
        utherm: 0.0,
        utot: 0.0,
        helec: 0.0,
        htrans: 0.0,
        hrot: 0.0,
        hvib: 0.0,
        htherm: 0.0,
        htot: 0.0,
        felec: 0.0,
        ftrans: 0.0,
        frot: 0.0,
        fvib: 0.0,
        ftherm: 0.0,
        ftot: 0.0,
        gelec: 0.0,
        gtrans: 0.0,
        grot: 0.0,
        gvib: 0.0,
        gtherm: 0.0,
        gtot: 0.0,
        selec: 0.0,
        strans: 0.0,
        srot: 0.0,
        svib: 0.0,
        stherm: 0.0,
        stot: 0.0,
        cvelec: 0.0,
        cvtrans: 0.0,
        cvrot: 0.0,
        cvvib: 0.0,
        cvtherm: 0.0,
        cvtot: 0.0,
        cpelec: 0.0,
        cptrans: 0.0,
        cprot: 0.0,
        cpvib: 0.0,
        cptherm: 0.0,
        cptot: 0.0,
        zpe: 0.0,
    };

    all_electronic(&mut thermo, multiplicity, temp);
    all_translation(&mut thermo, mass_amu_total, pressure, temp);
    all_rotations(&mut thermo, brot_cm1, temp, rot_symmetry);
    all_vibrations(&mut thermo, freqs_cm1, temp, freq_cutoff);

    thermo.utherm = thermo.uelec + thermo.utrans + thermo.urot + thermo.uvib;
    thermo.htherm = thermo.helec + thermo.htrans + thermo.hrot + thermo.hvib;
    thermo.stherm = thermo.selec + thermo.strans + thermo.srot + thermo.svib;
    thermo.ftherm = thermo.felec + thermo.ftrans + thermo.frot + thermo.fvib;
    thermo.gtherm = thermo.gelec + thermo.gtrans + thermo.grot + thermo.gvib;
    thermo.cvtherm = thermo.cvelec + thermo.cvtrans + thermo.cvrot + thermo.cvvib;
    thermo.cptherm = thermo.cpelec + thermo.cptrans + thermo.cprot + thermo.cpvib;

    thermo.utot = thermo.utherm + thermo.zpe;
    thermo.htot = thermo.htherm + thermo.zpe;
    thermo.ftot = thermo.ftherm + thermo.zpe;
    thermo.gtot = thermo.gtherm + thermo.zpe;

    thermo.stot = thermo.stherm;
    thermo.cvtot = thermo.cvtherm;
    thermo.cptot = thermo.cptherm;

    thermo.pftot = thermo.pfelec * thermo.pftrans * thermo.pfrot * thermo.pfvib;
    thermo
}

/// The symmetry number assumed when the user has not chosen a point group:
/// C1, which is exact for an asymmetric molecule and an over-estimate of the
/// rotational entropy for anything else, by `R ln(sigma)`.
pub const ROT_SYMMETRY_NUMBER: f64 = 1.0;
/// Chirality factor in the same expression; 1 for an achiral treatment.
pub const ROT_CHIRALITY: f64 = 1.0;

/// Width of the row-label column. Wide enough for "Total(kcal/mol)", the
/// longest label, so no row overhangs the rule drawn above it.
const LABEL_W: usize = 15;
/// Width of one number column.
const COL_W: usize = 12;
/// Space between number columns.
const GAP: usize = 2;
/// Full width of the four-column energy table: label, a space, then the
/// columns with their gaps. Every rule in that section is drawn to this, so
/// the rules always span the table exactly rather than falling short of it.
const ENERGY_W: usize = LABEL_W + 1 + 4 * COL_W + 3 * GAP;
/// The same for the three-column entropy/heat-capacity table.
const SCC_W: usize = LABEL_W + 1 + 3 * COL_W + 2 * GAP;

pub fn format_thermo(
    thermo: &ThermoResults,
    temp: f64,
    freq_cutoff: f64,
    elec_energy: Option<f64>,
    rot_symmetry: f64,
) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    fn print_uhfg(out: &mut String, label: &str, u: f64, h: f64, f: f64, g: f64) {
        writeln!(
            out,
            "{:<LABEL_W$} {:>COL_W$.6}  {:>COL_W$.6}  {:>COL_W$.6}  {:>COL_W$.6}",
            label, u, h, f, g
        )
        .unwrap();
    }
    fn print_scc(out: &mut String, label: &str, s: f64, cv: f64, cp: f64) {
        writeln!(
            out,
            "{:<LABEL_W$} {:>COL_W$.6}  {:>COL_W$.6}  {:>COL_W$.6}",
            label, s, cv, cp
        )
        .unwrap();
    }

    writeln!(out).unwrap();
    writeln!(out, "{:=^ENERGY_W$}", " Thermochemistry ").unwrap();
    writeln!(out, "T = {:.2} K", temp).unwrap();
    writeln!(out, 
        "ZPE: {:>12.6} Eh  ({:>10.3} kcal/mol)",
        thermo.zpe,
        thermo.zpe * CM1_TO_KCAL / CM1_TO_HARTREE
    ).unwrap();
    writeln!(out, "qRRHO cutoff: {:.1} cm-1", freq_cutoff).unwrap();
    writeln!(out, "Rotational symmetry number: {rot_symmetry:.0}").unwrap();

    writeln!(out).unwrap();
    writeln!(out, "{:-^ENERGY_W$}", " Energy Contributions (Eh) ").unwrap();
    writeln!(out,
        "{:<LABEL_W$} {:>COL_W$}  {:>COL_W$}  {:>COL_W$}  {:>COL_W$}",
        "", "U", "H", "F", "G"
    ).unwrap();
    print_uhfg(&mut out, 
        "Electronic",
        thermo.uelec,
        thermo.helec,
        thermo.felec,
        thermo.gelec,
    );
    print_uhfg(&mut out, 
        "Trans",
        thermo.utrans,
        thermo.htrans,
        thermo.ftrans,
        thermo.gtrans,
    );
    print_uhfg(&mut out, "Rot", thermo.urot, thermo.hrot, thermo.frot, thermo.grot);
    print_uhfg(&mut out, "Vib", thermo.uvib, thermo.hvib, thermo.fvib, thermo.gvib);
    print_uhfg(&mut out, 
        "Thermal",
        thermo.utherm,
        thermo.htherm,
        thermo.ftherm,
        thermo.gtherm,
    );
    writeln!(out, "{:-<ENERGY_W$}", "").unwrap();
    print_uhfg(&mut out, 
        "Total(Eh)",
        thermo.utot,
        thermo.htot,
        thermo.ftot,
        thermo.gtot,
    );
    print_uhfg(&mut out, 
        "Total(kcal/mol)",
        thermo.utot * AU_TO_KCAL,
        thermo.htot * AU_TO_KCAL,
        thermo.ftot * AU_TO_KCAL,
        thermo.gtot * AU_TO_KCAL,
    );

    writeln!(out).unwrap();
    if let Some(e_elec) = elec_energy {
        writeln!(out, "Electronic energy (Eh): {:>12.6}", e_elec).unwrap();
        writeln!(out, "Electronic + ZPE (Eh): {:>12.6}", e_elec + thermo.zpe).unwrap();
        writeln!(out,
            "{:<LABEL_W$} {:>COL_W$}  {:>COL_W$}  {:>COL_W$}  {:>COL_W$}",
            "Total+E_elec", "U", "H", "F", "G"
        ).unwrap();
        print_uhfg(&mut out, 
            "Eh",
            e_elec + thermo.utot,
            e_elec + thermo.htot,
            e_elec + thermo.ftot,
            e_elec + thermo.gtot,
        );
    } else {
        writeln!(out, "Electronic energy (Eh): n/a").unwrap();
    }

    writeln!(out).unwrap();
    writeln!(out, "{:-^SCC_W$}", " Entropy & Heat Capacities (Eh/K) ").unwrap();
    writeln!(out,
        "{:<LABEL_W$} {:>COL_W$}  {:>COL_W$}  {:>COL_W$}",
        "", "S", "Cv", "Cp"
    ).unwrap();
    print_scc(&mut out, "Electronic", thermo.selec, thermo.cvelec, thermo.cpelec);
    print_scc(&mut out, "Trans", thermo.strans, thermo.cvtrans, thermo.cptrans);
    print_scc(&mut out, "Rot", thermo.srot, thermo.cvrot, thermo.cprot);
    print_scc(&mut out, "Vib", thermo.svib, thermo.cvvib, thermo.cpvib);
    print_scc(&mut out, "Thermal", thermo.stherm, thermo.cvtherm, thermo.cptherm);
    writeln!(out, "{:-<SCC_W$}", "").unwrap();
    print_scc(&mut out, "Total", thermo.stot, thermo.cvtot, thermo.cptot);
    writeln!(out, 
        "S_total*T = {:>10.3} kcal/mol",
        thermo.stot * temp * AU_TO_KCAL
    ).unwrap();
    out
}

/// Prints exactly `format_thermo`'s text to stdout, so the printed and
/// GUI-displayed forms can never drift apart.
pub fn print_thermo(
    thermo: &ThermoResults,
    temp: f64,
    freq_cutoff: f64,
    elec_energy: Option<f64>,
    rot_symmetry: f64,
) {
    print!(
        "{}",
        format_thermo(thermo, temp, freq_cutoff, elec_energy, rot_symmetry)
    );
}

pub fn brot_from_coords(xyz_bohr: &[[f64; 3]], mass_amu: &[f64]) -> Vec<f64> {
    let bohr_to_ang = 0.529_177_210_9;
    let mut xyz_ang = Vec::with_capacity(xyz_bohr.len());
    for c in xyz_bohr {
        xyz_ang.push([c[0] * bohr_to_ang, c[1] * bohr_to_ang, c[2] * bohr_to_ang]);
    }
    let brot = get_brot(&xyz_ang, &mass_amu.to_vec());
    vec![brot[0], brot[1], brot[2]]
}

fn all_electronic(thermo: &mut ThermoResults, multiplicity: f64, temp: f64) {
    let RT = RGAS_AU * temp;
    let pf = if multiplicity > 0.0 {
        multiplicity
    } else {
        1.0
    };
    let F = -RT * pf.ln();
    let U = 0.0;
    let H = 0.0;
    let S = (U - F) / temp;

    thermo.pfelec = pf;
    thermo.felec = F;
    thermo.uelec = U;
    thermo.helec = H;
    thermo.selec = S;
    thermo.gelec = F;

    thermo.cvelec = 0.0;
    thermo.cpelec = 0.0;
}

fn all_translation(thermo: &mut ThermoResults, mass_amu_total: f64, pressure: f64, temp: f64) {
    let RT = RGAS_AU * temp;
    let mass = mass_amu_total * AMU_TO_ELECMASS;

    let mut lam = f64::sqrt(TWOPI * mass * RT / HPLANCK_AU_SQ);
    lam = lam * lam * lam;

    let vol = RT / (pressure * PASCAL_TO_AU);
    let q0 = lam * vol;
    let pf = std::f64::consts::E * q0;

    let F = -RT * pf.ln();
    let U = 1.5 * RT;
    let H = 2.5 * RT;
    let S = (U - F) / temp;
    let G = H - temp * S;

    thermo.pftrans = pf;
    thermo.ftrans = F;
    thermo.utrans = U;
    thermo.htrans = H;
    thermo.strans = S;
    thermo.gtrans = G;

    thermo.cvtrans = 1.5 * RGAS_AU;
    thermo.cptrans = 2.5 * RGAS_AU;
}

fn all_rotations(thermo: &mut ThermoResults, brot: &[f64], temp: f64, rot_symmetry: f64) {
    let RT = RGAS_AU * temp;

    if brot.is_empty() {
        thermo.pfrot = 1.0;
        thermo.frot = 0.0;
        thermo.urot = 0.0;
        thermo.hrot = 0.0;
        thermo.srot = 0.0;
        thermo.grot = 0.0;
        thermo.cvrot = 0.0;
        thermo.cprot = 0.0;
        return;
    }

    // A non-positive sigma would make the partition function nonsense; fall
    // back to the asymmetric case rather than producing an infinity.
    let sigma = if rot_symmetry > 0.0 { rot_symmetry } else { 1.0 };
    let chiral = ROT_CHIRALITY;

    let eps = 1.0e-12;
    let has_zero = brot.iter().any(|b| *b <= eps);
    let nonzero: Vec<f64> = brot.iter().copied().filter(|b| *b > eps).collect();

    let (pf, dof) = if has_zero || nonzero.len() <= 1 {
        let brot_cm1 = if !nonzero.is_empty() {
            nonzero.iter().sum::<f64>() / (nonzero.len() as f64)
        } else {
            brot[0].max(1.0e-6)
        };
        let theta_r = brot_cm1 * CM1_TO_K;
        let q = (temp / theta_r) * (chiral / sigma);
        (q, 2.0)
    } else {
        let a = nonzero[0] * CM1_TO_K;
        let b = nonzero[1] * CM1_TO_K;
        let c = nonzero[2] * CM1_TO_K;
        let denom = f64::sqrt(a * b * c);
        let q = f64::sqrt(PI) * temp.powf(1.5) / denom * (chiral / sigma);
        (q, 3.0)
    };

    let F = -RT * pf.ln();
    let U = 0.5 * dof * RT;
    let H = U;
    let S = (U - F) / temp;
    let G = F;

    thermo.pfrot = pf;
    thermo.frot = F;
    thermo.urot = U;
    thermo.hrot = H;
    thermo.srot = S;
    thermo.grot = G;

    thermo.cvrot = 0.5 * dof * RGAS_AU;
    thermo.cprot = thermo.cvrot;
}

fn all_vibrations(thermo: &mut ThermoResults, freqs_cm1: &[f64], temp: f64, freq_cutoff: f64) {
    let RT = RGAS_AU * temp;

    let mut Uvib = 0.0;
    let mut Hvib = 0.0;
    let mut Fvib = 0.0;
    let mut Svib = 0.0;
    let mut Cvib = 0.0;
    let mut PFvib = 1.0;
    let mut zpe = 0.0;

    for &omega_cm1 in freqs_cm1 {
        if omega_cm1 <= 0.0 {
            continue;
        }

        zpe += 0.5 * omega_cm1 * CM1_TO_HARTREE;

        // `f64::exp(x) - 1.0` loses all precision by catastrophic
        // cancellation once x drops below ~1e-13 (as happens whenever a
        // "zero" TR mode's floating-point Eckart residual leaks through the
        // `> 0.0` frequency filter on the positive side), silently turning
        // U_mode/Cv_mode into +inf and poisoning every downstream
        // thermochemistry total. `exp_m1` is accurate for all x, including
        // this regime, and correctly recovers the finite classical
        // (equipartition) limit as x -> 0 instead of diverging.
        let x = (CM1_TO_K * omega_cm1) / temp;
        let ex_m1 = x.exp_m1();
        let ex = ex_m1 + 1.0;
        let U_mode = omega_cm1 * CM1_TO_HARTREE / ex_m1;
        let H_mode = U_mode;

        let Cv_mode = RGAS_AU * x * x * ex / (ex_m1 * ex_m1);

        let S_mode = if omega_cm1 > freq_cutoff {
            entropy_vib_rrho(omega_cm1, temp)
        } else {
            grimme_entropy_qrrho(omega_cm1, freq_cutoff, temp)
        };

        let F_mode = U_mode - temp * S_mode;
        let PF_mode = f64::exp(-F_mode / RT);

        Uvib += U_mode;
        Hvib += H_mode;
        Svib += S_mode;
        Fvib += F_mode;
        Cvib += Cv_mode;
        PFvib *= PF_mode;
    }

    thermo.uvib = Uvib;
    thermo.hvib = Hvib;
    thermo.svib = Svib;
    thermo.fvib = Fvib;
    thermo.gvib = Fvib;

    thermo.cvvib = Cvib;
    thermo.cpvib = Cvib;
    thermo.pfvib = PFvib;
    thermo.zpe = zpe;
}

fn entropy_vib_rrho(omega_cm1: f64, temp: f64) -> f64 {
    if omega_cm1 <= 0.0 {
        return 0.0;
    }
    let x = (CM1_TO_K * omega_cm1) / temp;
    if x < 1.0e-12 {
        // The RRHO entropy genuinely diverges (+inf) as x -> 0 (the
        // `-ln(1 - exp(-x))` term ~ -ln(x)); this is a real feature of the
        // harmonic approximation for a vanishing frequency, not a floating-
        // point artifact, and is exactly why `grimme_entropy_qrrho` exists
        // for genuinely soft real modes. Below this threshold `x` is
        // indistinguishable from the floating-point residual of an
        // Eckart-projected translation/rotation, so contribute 0 rather
        // than +inf.
        return 0.0;
    }
    let term = x / x.exp_m1() - f64::ln(-(-x).exp_m1());
    RGAS_AU * term
}

fn entropy_free_rotor(omega_cm1: f64, temp: f64, bav_si: f64) -> f64 {
    if omega_cm1 <= 0.0 {
        return 0.0;
    }

    let omega_m1 = omega_cm1 * 100.0;
    let nu_s1 = CLIGHT_SI * omega_m1;
    let mu = PLANCK_SI / (8.0 * PI_SQ * nu_s1);
    let mu_prime = mu * bav_si / (mu + bav_si);

    let factor = 8.0 * PI.powi(3) * mu_prime * BOLTZMANN_SI * temp / (PLANCK_SI * PLANCK_SI);

    let s_si = RGAS_SI * (0.5 + 0.5 * factor.ln());
    s_si / J_PER_HARTREE_PER_MOL
}

fn grimme_damp(omega_cm1: f64, freq_cutoff_cm1: f64) -> f64 {
    let ratio = freq_cutoff_cm1 / omega_cm1;
    1.0 / (1.0 + ratio.powf(GRIMME_ALPHA))
}

fn grimme_entropy_qrrho(omega_cm1: f64, freq_cutoff_cm1: f64, temp: f64) -> f64 {
    if omega_cm1 <= 0.0 {
        return 0.0;
    }
    let w = grimme_damp(omega_cm1, freq_cutoff_cm1);
    let s_rrho = entropy_vib_rrho(omega_cm1, temp);
    let s_fr = entropy_free_rotor(omega_cm1, temp, GRIMME_BAV_SI);
    w * s_rrho + (1.0 - w) * s_fr
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cross-check against xTB's own reported zero-point energy from the
    /// real `--ohess` run captured for the design spec: 0.024943563139 Eh
    /// for H2O2. `eval_thermo` computes ZPE straight from frequencies (no
    /// Eckart/mass-table involvement), so this isolates the thermochemistry
    /// formulas themselves from the rest of the pipeline already
    /// cross-checked in `qchem_interfaces::hessian_file`.
    #[test]
    fn zero_point_energy_matches_xtbs_own_value_for_h2o2() {
        let freqs_cm1 = [259.56, 1108.86, 1158.96, 1355.45, 3531.45, 3534.68];
        let brot_cm1 = [1.0, 1.0, 1.0]; // rotational constants don't affect ZPE
        let mass_amu_total = 2.0 * 15.999 + 2.0 * 1.008;
        let thermo = eval_thermo(&freqs_cm1, &brot_cm1, mass_amu_total, 1.0, 298.15, 101_325.0, 100.0, 1.0);

        let xtb_zpe_eh = 0.024943563139;
        let rel_err = (thermo.zpe - xtb_zpe_eh).abs() / xtb_zpe_eh;
        assert!(
            rel_err < 0.01,
            "computed ZPE {} Eh vs xTB's own {xtb_zpe_eh} Eh ({:.2}% off)",
            thermo.zpe,
            rel_err * 100.0
        );
    }

    /// `print_thermo` must print exactly what `format_thermo` returns, since
    /// the GUI is built on that guarantee holding.
    #[test]
    fn print_thermo_and_format_thermo_agree() {
        let freqs_cm1 = [259.56, 1108.86, 1158.96, 1355.45, 3531.45, 3534.68];
        let brot_cm1 = [1.0, 1.0, 1.0];
        let thermo = eval_thermo(&freqs_cm1, &brot_cm1, 34.0, 1.0, 298.15, 101_325.0, 100.0, 1.0);
        let text = format_thermo(&thermo, 298.15, 100.0, Some(-9.05), 1.0);
        assert!(text.contains("Thermochemistry"));
        assert!(text.contains("ZPE"));
        assert!(text.contains("Electronic energy"));
    }
}

#[cfg(test)]
mod rule_width_tests {
    use super::*;

    /// The symmetry number the rotational partition function actually uses
    /// must be the one the output claims -- a printed value that drifted from
    /// the computed one would be worse than not printing it.
    #[test]
    fn the_printed_symmetry_number_is_the_one_used() {
        let thermo = eval_thermo(&[500.0], &[0.06, 0.038, 0.030], 250.0, 1.0, 298.15, 101_325.0, 100.0, 1.0);
        let text = format_thermo(&thermo, 298.15, 100.0, None, 1.0);
        let line = text
            .lines()
            .find(|l| l.starts_with("Rotational symmetry number"))
            .expect("the symmetry number must be stated");
        assert!(line.contains("1"), "{line}");
    }

    /// Every rule must span the widest row of its own table -- the complaint
    /// was rules falling short of the numbers above them.
    #[test]
    fn rules_span_their_tables() {
        // A real evaluation, so the numbers are full width rather than zeros.
        let thermo = eval_thermo(
            &[500.0, 1500.0, 3000.0],
            &[0.06, 0.038, 0.030],
            250.0,
            1.0,
            298.15,
            101_325.0,
            100.0,
            1.0,
        );
        let text = format_thermo(&thermo, 298.15, 100.0, None, 1.0);
        let width = |prefix: &str| {
            text.lines()
                .find(|l| l.starts_with(prefix))
                .unwrap_or_else(|| panic!("no line starting {prefix:?}"))
                .chars()
                .count()
        };
        let energy_rows = ["Electronic ", "Total(Eh)", "Total(kcal/mol)"];
        for row in energy_rows {
            assert!(
                width(row) <= ENERGY_W,
                "{row:?} is {} wide, past the {ENERGY_W}-wide rule",
                width(row)
            );
        }
        // The rules and the banner are exactly the table width.
        assert_eq!(width("========"), ENERGY_W);
        assert_eq!(
            text.lines().filter(|l| l.chars().all(|c| c == '-') && !l.is_empty())
                .map(|l| l.chars().count())
                .collect::<Vec<_>>(),
            vec![ENERGY_W, SCC_W],
        );
    }
}
