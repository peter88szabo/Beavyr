//! DREIDING energy terms and their analytic gradients.
//!
//! Every function accumulates `dE/dx` into `grad` and returns its energy contribution, so the
//! caller sums terms without allocating. Forces are `-grad`; the negation happens once, in
//! [`super::energy_and_forces`].
//!
//! Coordinates are a flat `[x0, y0, z0, x1, ...]` in Å, energies in kcal/mol. Analytic gradients
//! are checked against central finite differences in this module's tests -- a hand-derived
//! gradient that is subtly wrong makes a minimiser converge to the wrong structure without ever
//! failing, so the check is not optional.
//!
//! Equation numbers refer to Mayo, Olafson & Goddard, *J. Phys. Chem.* **1990**, *94*, 8897.

/// A 3-vector view into the flat coordinate array.
type V3 = [f64; 3];

#[inline]
fn get(pos: &[f64], i: usize) -> V3 {
    [pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]]
}

#[inline]
fn add(grad: &mut [f64], i: usize, v: V3) {
    grad[3 * i] += v[0];
    grad[3 * i + 1] += v[1];
    grad[3 * i + 2] += v[2];
}

#[inline]
fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

/// Below this length a vector's direction is numerically meaningless, so the term is skipped
/// rather than dividing by it. Two atoms this close are already a pathological input.
const TINY: f64 = 1.0e-12;

// ---------------------------------------------------------------------------------------------
// Bond stretch -- Eq 4a
// ---------------------------------------------------------------------------------------------

/// Harmonic bond stretch: `E = ½k(R − Rₑ)²`.
///
/// The paper's Morse alternative (Eq 5a, `DREIDING/M`) is deliberately not used. Its own reasoning
/// applies directly here: a sketched starting geometry may be far from equilibrium, and Morse
/// forces vanish at long range where harmonic forces keep pulling, so the harmonic form recovers
/// from a bad guess and Morse does not.
pub fn bond(pos: &[f64], grad: &mut [f64], i: usize, j: usize, k: f64, r_e: f64) -> f64 {
    let d = sub(get(pos, i), get(pos, j));
    let r = norm(d);
    if r < TINY {
        return 0.0;
    }
    let dr = r - r_e;
    // dE/dr = k·dr, and dr/dx_i = d/r
    let f = scale(d, k * dr / r);
    add(grad, i, f);
    add(grad, j, scale(f, -1.0));
    0.5 * k * dr * dr
}

// ---------------------------------------------------------------------------------------------
// Angle bend -- Eq 10a and 10'
// ---------------------------------------------------------------------------------------------

/// Harmonic-cosine angle bend about centre `j`: `E = ½C[cosθ − cosθ⁰]²`, Eq 10a.
///
/// `linear` selects Eq 10' instead, `E = K[1 + cosθ]`, for centres whose equilibrium angle is
/// 180°. The paper prefers the cosine forms over a plain harmonic in θ because the latter does not
/// give zero slope as θ approaches 180°, and because at θ⁰ = 180° the harmonic-cosine constant
/// `C = K/sin²θ⁰` diverges.
pub fn angle(
    pos: &[f64],
    grad: &mut [f64],
    i: usize,
    j: usize,
    k: usize,
    c: f64,
    cos_theta0: f64,
    linear: bool,
) -> f64 {
    let u = sub(get(pos, i), get(pos, j));
    let v = sub(get(pos, k), get(pos, j));
    let (lu, lv) = (norm(u), norm(v));
    if lu < TINY || lv < TINY {
        return 0.0;
    }
    let cos_t = (dot(u, v) / (lu * lv)).clamp(-1.0, 1.0);

    let (energy, de_dcos) = if linear {
        (c * (1.0 + cos_t), c)
    } else {
        let d = cos_t - cos_theta0;
        (0.5 * c * d * d, c * d)
    };

    // d(cosθ)/du = v/(|u||v|) − cosθ·u/|u|²,  symmetrically for v.
    let dcos_du = sub(scale(v, 1.0 / (lu * lv)), scale(u, cos_t / (lu * lu)));
    let dcos_dv = sub(scale(u, 1.0 / (lu * lv)), scale(v, cos_t / (lv * lv)));
    let gi = scale(dcos_du, de_dcos);
    let gk = scale(dcos_dv, de_dcos);
    add(grad, i, gi);
    add(grad, k, gk);
    add(grad, j, scale(add_v(gi, gk), -1.0)); // translation invariance fixes the centre
    energy
}

#[inline]
fn add_v(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

// ---------------------------------------------------------------------------------------------
// Torsion -- Eq 13
// ---------------------------------------------------------------------------------------------

/// Proper torsion `i-j-k-l` about the central bond `j-k`: `E = ½V[1 − cos n(φ − φ⁰)]`, Eq 13.
///
/// `v` must already carry the paper's division by the number of (I,L) pairs sharing the central
/// bond: it states "2/9 for each of the nine possibilities of I and L" for substituted ethane and
/// "45/4 for each of the four possibilities" for ethylene. Without that division ethane's barrier
/// comes out near 18 kcal/mol instead of ~2, against ~2.9 experimentally.
pub fn torsion(
    pos: &[f64],
    grad: &mut [f64],
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    v: f64,
    n: f64,
    phi0: f64,
) -> f64 {
    let b1 = sub(get(pos, j), get(pos, i));
    let b2 = sub(get(pos, k), get(pos, j));
    let b3 = sub(get(pos, l), get(pos, k));

    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let lb2 = norm(b2);
    // Collinear i-j-k or j-k-l leaves the dihedral undefined.
    if dot(n1, n1) < TINY || dot(n2, n2) < TINY || lb2 < TINY {
        return 0.0;
    }
    let b2h = scale(b2, 1.0 / lb2);

    // Signed dihedral from its cosine-like and sine-like parts. atan2 keeps the sign that a
    // cosine alone would lose, which matters because n is odd for sp3 centres.
    let c = dot(n1, n2);
    let m = dot(cross(n1, n2), b2h);
    let phi = m.atan2(c);

    let arg = n * (phi - phi0);
    let energy = 0.5 * v * (1.0 - arg.cos());
    let de_dphi = 0.5 * v * n * arg.sin();

    // dφ/dr, by chain rule through n1, n2 and b̂2 rather than a remembered closed form.
    //
    //   φ = atan2(m, c)  ⇒  dφ = (c·dm − m·dc)/(c² + m²)
    //   c = n1·n2                ⇒  dc = n2·dn1 + n1·dn2
    //   m = (n1×n2)·b̂2           ⇒  dm = (n2×b̂2)·dn1 + (b̂2×n1)·dn2 + (n1×n2)·db̂2
    //
    // then n1 = b1×b2 and n2 = b2×b3 push those onto db1, db2, db3, and finally
    // b1 = r_j − r_i, b2 = r_k − r_j, b3 = r_l − r_k give the four atomic gradients. Deriving it
    // this way is longer but mechanical; the finite-difference test in this module is what
    // established that the closed forms in circulation are easy to misremember.
    let denom = c * c + m * m;
    if denom < TINY {
        return 0.0;
    }
    let coeff_n1 = scale(sub(scale(cross(n2, b2h), c), scale(n2, m)), 1.0 / denom);
    let coeff_n2 = scale(sub(scale(cross(b2h, n1), c), scale(n1, m)), 1.0 / denom);
    let coeff_b2h = scale(cross(n1, n2), c / denom);

    // db̂2/db2 projects out the component along b̂2 and divides by |b2|.
    let along = dot(b2h, coeff_b2h);
    let coeff_b2 = scale(sub(coeff_b2h, scale(b2h, along)), 1.0 / lb2);

    // Gradients with respect to the three bond vectors.
    let g1 = cross(b2, coeff_n1);
    let g2 = add_v(
        add_v(cross(coeff_n1, b1), cross(b3, coeff_n2)),
        coeff_b2,
    );
    let g3 = cross(coeff_n2, b2);

    let dphi_di = scale(g1, -1.0);
    let dphi_dj = sub(g1, g2);
    let dphi_dk = sub(g2, g3);
    let dphi_dl = g3;

    add(grad, i, scale(dphi_di, de_dphi));
    add(grad, j, scale(dphi_dj, de_dphi));
    add(grad, k, scale(dphi_dk, de_dphi));
    add(grad, l, scale(dphi_dl, de_dphi));
    energy
}

// ---------------------------------------------------------------------------------------------
// Inversion -- Eq 28, 29
// ---------------------------------------------------------------------------------------------

/// Inversion (umbrella) term at centre `i` bonded to exactly `j`, `k`, `l`.
///
/// `psi` is the angle between the `i→l` bond and the `j-i-k` plane. For a planar centre
/// (`cos_psi0 == 1`, i.e. ψ⁰ = 0°) the harmonic-cosine constant `K/sin²ψ⁰` diverges, exactly as
/// for a linear angle, so the linear form `E = K[1 − cosψ]` is used instead. Otherwise
/// `E = ½(K/sin²ψ⁰)[cosψ − cosψ⁰]²`.
///
/// The caller supplies the three permutations of (j, k, l) each weighted 1/3, per the paper:
/// "we add together all three possible inversion terms with each weighted by a factor of 1/3".
pub fn inversion(
    pos: &[f64],
    grad: &mut [f64],
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    c: f64,
    cos_psi0: f64,
    planar: bool,
) -> f64 {
    let ri = get(pos, i);
    let u = sub(get(pos, j), ri);
    let v = sub(get(pos, k), ri);
    let w = sub(get(pos, l), ri);

    let nrm = cross(u, v);
    let ln = norm(nrm);
    let lw = norm(w);
    if ln < TINY || lw < TINY {
        return 0.0;
    }

    // sinψ is the projection of the i→l bond onto the plane normal.
    let sin_psi = (dot(nrm, w) / (ln * lw)).clamp(-1.0, 1.0);
    let cos_psi = (1.0 - sin_psi * sin_psi).max(0.0).sqrt();
    if cos_psi < TINY {
        // ψ = ±90°: the i→l bond is along the normal. dcosψ/dsinψ is singular here, and the
        // geometry is already absurd for a three-coordinate centre.
        return 0.0;
    }

    let (energy, de_dcos) = if planar {
        (c * (1.0 - cos_psi), -c)
    } else {
        let d = cos_psi - cos_psi0;
        (0.5 * c * d * d, c * d)
    };
    // cosψ = √(1 − sin²ψ)  ⇒  dcosψ/dsinψ = −sinψ/cosψ
    let de_dsin = de_dcos * (-sin_psi / cos_psi);

    // sinψ = (n·w)/(|n||w|) with n = u × v. Differentiate each factor in turn.
    let n_dot_w = dot(nrm, w);
    let inv = 1.0 / (ln * lw);

    // ∂sinψ/∂w = n/(|n||w|) − (n·w)·w/(|n||w|³)
    let ds_dw = sub(scale(nrm, inv), scale(w, n_dot_w / (ln * lw * lw * lw)));
    // ∂sinψ/∂n = w/(|n||w|) − (n·w)·n/(|n|³|w|)
    let ds_dn = sub(scale(w, inv), scale(nrm, n_dot_w / (ln * ln * ln * lw)));
    // n = u × v  ⇒  ∂/∂u = (∂/∂n) × ... : dn·(v × ds_dn) for u, and (ds_dn × u) for v
    let ds_du = cross(v, ds_dn);
    let ds_dv = cross(ds_dn, u);

    let gj = scale(ds_du, de_dsin);
    let gk = scale(ds_dv, de_dsin);
    let gl = scale(ds_dw, de_dsin);
    add(grad, j, gj);
    add(grad, k, gk);
    add(grad, l, gl);
    add(grad, i, scale(add_v(add_v(gj, gk), gl), -1.0));
    energy
}

// ---------------------------------------------------------------------------------------------
// van der Waals -- Table II
// ---------------------------------------------------------------------------------------------

/// Lennard-Jones 12-6 for one pair, with DREIDING's geometric-mean combination applied by the
/// caller: `R⁰_IJ = √(R⁰_I R⁰_J)`, `D⁰_IJ = √(D⁰_I D⁰_J)`.
///
/// The paper's exponential-6 (X6) alternative is not implemented. 1-2 and 1-3 pairs are excluded
/// by the caller; 1-4 pairs are kept at full strength, DREIDING specifying no 1-4 scaling.
pub fn vdw_pair(pos: &[f64], grad: &mut [f64], i: usize, j: usize, r0: f64, d0: f64) -> f64 {
    let d = sub(get(pos, i), get(pos, j));
    let r2 = dot(d, d);
    if r2 < TINY {
        return 0.0;
    }
    let s2 = r0 * r0 / r2;
    let s6 = s2 * s2 * s2;
    let s12 = s6 * s6;
    let energy = d0 * (s12 - 2.0 * s6);
    // dE/dr² = d0·(−6·s12 + 6·s6)/r², and dr²/dx_i = 2d
    let de_dr2 = d0 * 6.0 * (s6 - s12) / r2;
    let f = scale(d, 2.0 * de_dr2);
    add(grad, i, f);
    add(grad, j, scale(f, -1.0));
    energy
}

// ---------------------------------------------------------------------------------------------
// Hydrogen bond -- Table V
// ---------------------------------------------------------------------------------------------

/// Explicit hydrogen bond between donor `d`, its hydrogen `h`, and acceptor `a`:
/// `E = D_hb[5(R_hb/R)¹² − 6(R_hb/R)¹⁰]cos⁴θ_DHA`, with `R` the donor-acceptor distance.
///
/// `D_hb` depends on the charge convention (Table V). Beavyr has no charge model, so the
/// no-charges value is used -- see [`super::params::D_HB`].
pub fn hbond(
    pos: &[f64],
    grad: &mut [f64],
    d: usize,
    h: usize,
    a: usize,
    d_hb: f64,
    r_hb: f64,
) -> f64 {
    let rda = sub(get(pos, a), get(pos, d));
    let r = norm(rda);
    if r < TINY {
        return 0.0;
    }
    // Angular factor cos⁴θ at the hydrogen, θ = angle D-H-A.
    let u = sub(get(pos, d), get(pos, h));
    let v = sub(get(pos, a), get(pos, h));
    let (lu, lv) = (norm(u), norm(v));
    if lu < TINY || lv < TINY {
        return 0.0;
    }
    let cos_t = (dot(u, v) / (lu * lv)).clamp(-1.0, 1.0);
    // Only the linear-ish half of the range binds; a bent arrangement is not a hydrogen bond.
    if cos_t >= 0.0 {
        return 0.0;
    }
    let ang = cos_t * cos_t * cos_t * cos_t;
    let dang_dcos = 4.0 * cos_t * cos_t * cos_t;

    let s = r_hb / r;
    let s10 = s.powi(10);
    let s12 = s10 * s * s;
    let radial = d_hb * (5.0 * s12 - 6.0 * s10);
    // dE_radial/dR = D·(−60·s12 + 60·s10)/R
    let dradial_dr = d_hb * 60.0 * (s10 - s12) / r;

    // Radial part: acts along the donor-acceptor axis.
    let fr = scale(rda, dradial_dr * ang / r);
    add(grad, a, fr);
    add(grad, d, scale(fr, -1.0));

    // Angular part: same chain rule as `angle`, centred on the hydrogen.
    let de_dcos = radial * dang_dcos;
    let dcos_du = sub(scale(v, 1.0 / (lu * lv)), scale(u, cos_t / (lu * lu)));
    let dcos_dv = sub(scale(u, 1.0 / (lu * lv)), scale(v, cos_t / (lv * lv)));
    let gd = scale(dcos_du, de_dcos);
    let ga = scale(dcos_dv, de_dcos);
    add(grad, d, gd);
    add(grad, a, ga);
    add(grad, h, scale(add_v(gd, ga), -1.0));

    radial * ang
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Central finite differences on every degree of freedom, compared against the analytic
    /// gradient the term reported.
    fn check_gradient<F>(pos: &[f64], mut eval: F, tol: f64, what: &str)
    where
        F: FnMut(&[f64], &mut [f64]) -> f64,
    {
        let n = pos.len();
        let mut analytic = vec![0.0; n];
        eval(pos, &mut analytic);

        const H: f64 = 1.0e-6;
        let mut scratch = vec![0.0; n];
        for c in 0..n {
            let mut plus = pos.to_vec();
            let mut minus = pos.to_vec();
            plus[c] += H;
            minus[c] -= H;
            scratch.iter_mut().for_each(|g| *g = 0.0);
            let e_plus = eval(&plus, &mut scratch);
            scratch.iter_mut().for_each(|g| *g = 0.0);
            let e_minus = eval(&minus, &mut scratch);
            let numeric = (e_plus - e_minus) / (2.0 * H);
            let scale = analytic[c].abs().max(numeric.abs()).max(1.0);
            assert!(
                (analytic[c] - numeric).abs() / scale < tol,
                "{what}: coordinate {c} analytic {} vs numeric {numeric}",
                analytic[c]
            );
        }
    }

    /// A gradient must sum to zero over all atoms, or the term applies a net force to the whole
    /// molecule and an MD run would see it drift for no physical reason.
    fn check_translation_invariance(pos: &[f64], mut eval: impl FnMut(&[f64], &mut [f64]) -> f64, what: &str) {
        let mut g = vec![0.0; pos.len()];
        eval(pos, &mut g);
        for axis in 0..3 {
            let net: f64 = g.iter().skip(axis).step_by(3).sum();
            assert!(net.abs() < 1.0e-8, "{what}: net force on axis {axis} is {net}");
        }
    }

    #[test]
    fn bond_gradient_matches_finite_differences() {
        let pos = [0.1, 0.0, 0.0, 1.6, 0.2, -0.1];
        check_gradient(&pos, |p, g| bond(p, g, 0, 1, 700.0, 1.53), 1.0e-6, "bond");
        check_translation_invariance(&pos, |p, g| bond(p, g, 0, 1, 700.0, 1.53), "bond");
    }

    #[test]
    fn bond_energy_is_zero_at_equilibrium_and_positive_either_side() {
        let mut g = vec![0.0; 6];
        let at_eq = [0.0, 0.0, 0.0, 1.53, 0.0, 0.0];
        assert!(bond(&at_eq, &mut g, 0, 1, 700.0, 1.53).abs() < 1.0e-12);
        g.iter_mut().for_each(|x| *x = 0.0);
        let short = [0.0, 0.0, 0.0, 1.40, 0.0, 0.0];
        assert!(bond(&short, &mut g, 0, 1, 700.0, 1.53) > 0.0);
        g.iter_mut().for_each(|x| *x = 0.0);
        let long = [0.0, 0.0, 0.0, 1.70, 0.0, 0.0];
        assert!(bond(&long, &mut g, 0, 1, 700.0, 1.53) > 0.0);
    }

    #[test]
    fn angle_gradient_matches_finite_differences() {
        // A bent triple, away from both 0 and 180 degrees.
        let pos = [1.0, 0.2, 0.0, 0.0, 0.0, 0.0, -0.3, 0.95, 0.1];
        let cos109 = (109.471_f64).to_radians().cos();
        let c = 100.0 / (109.471_f64).to_radians().sin().powi(2);
        check_gradient(
            &pos,
            |p, g| angle(p, g, 0, 1, 2, c, cos109, false),
            1.0e-6,
            "angle",
        );
        check_translation_invariance(&pos, |p, g| angle(p, g, 0, 1, 2, c, cos109, false), "angle");
    }

    #[test]
    fn linear_angle_gradient_matches_finite_differences() {
        // Near-linear, which is where the Eq 10a form would misbehave and Eq 10' is used instead.
        let pos = [1.2, 0.05, 0.0, 0.0, 0.0, 0.0, -1.1, 0.04, 0.02];
        check_gradient(
            &pos,
            |p, g| angle(p, g, 0, 1, 2, 100.0, -1.0, true),
            1.0e-6,
            "linear angle",
        );
    }

    #[test]
    fn torsion_gradient_matches_finite_differences() {
        // A generic gauche-ish arrangement: no two of the three bond vectors are collinear.
        let pos = [
            1.0, 1.0, 0.2, //
            0.5, 0.0, 0.0, //
            -0.9, 0.1, 0.05, //
            -1.5, 1.0, 0.7,
        ];
        let phi0 = std::f64::consts::PI;
        check_gradient(
            &pos,
            |p, g| torsion(p, g, 0, 1, 2, 3, 2.0 / 9.0, 3.0, phi0),
            1.0e-6,
            "torsion",
        );
        check_translation_invariance(
            &pos,
            |p, g| torsion(p, g, 0, 1, 2, 3, 2.0 / 9.0, 3.0, phi0),
            "torsion",
        );
    }

    /// A torsion term must be periodic in φ with period 2π/n and must not depend on the sign
    /// convention of the dihedral for these symmetric cases.
    #[test]
    fn torsion_energy_is_periodic() {
        let mut g = vec![0.0; 12];
        let energy_at = |phi_deg: f64, g: &mut Vec<f64>| {
            let a = phi_deg.to_radians();
            // Build i-j-k-l with dihedral exactly `phi_deg` about the x axis.
            let pos = [
                0.0,
                1.0,
                0.0, //
                0.0,
                0.0,
                0.0, //
                1.5,
                0.0,
                0.0, //
                1.5,
                a.cos(),
                a.sin(),
            ];
            g.iter_mut().for_each(|x| *x = 0.0);
            torsion(&pos, g, 0, 1, 2, 3, 2.0, 3.0, std::f64::consts::PI)
        };
        for phi in [0.0, 17.0, 45.0, 123.0] {
            let a = energy_at(phi, &mut g);
            let b = energy_at(phi + 120.0, &mut g);
            assert!((a - b).abs() < 1.0e-9, "n=3 torsion not 120°-periodic at {phi}");
        }
    }

    #[test]
    fn inversion_gradient_matches_finite_differences() {
        // A pyramidal centre: out of plane, but well away from the ψ = ±90° singularity.
        let pos = [
            0.0, 0.0, 0.0, // centre
            1.0, 0.0, -0.2, //
            -0.5, 0.9, -0.2, //
            -0.5, -0.9, -0.25,
        ];
        check_gradient(
            &pos,
            |p, g| inversion(p, g, 0, 1, 2, 3, 40.0, 1.0, true),
            1.0e-5,
            "planar inversion",
        );
        check_translation_invariance(
            &pos,
            |p, g| inversion(p, g, 0, 1, 2, 3, 40.0, 1.0, true),
            "planar inversion",
        );
        let cos_psi0 = (54.74_f64).to_radians().cos();
        let c = 40.0 / (54.74_f64).to_radians().sin().powi(2);
        check_gradient(
            &pos,
            |p, g| inversion(p, g, 0, 1, 2, 3, c, cos_psi0, false),
            1.0e-5,
            "pyramidal inversion",
        );
    }

    /// The planar form must be at its minimum when the centre is planar, and rise as it puckers.
    /// This is the term that keeps ethylene and benzene flat, so the sign matters.
    #[test]
    fn planar_inversion_penalises_pyramidalisation() {
        let mut g = vec![0.0; 12];
        let at_height = |h: f64, g: &mut Vec<f64>| {
            let pos = [
                0.0, 0.0, 0.0, //
                1.0, 0.0, 0.0, //
                -0.5, 0.866, 0.0, //
                -0.5, -0.866, h,
            ];
            g.iter_mut().for_each(|x| *x = 0.0);
            inversion(&pos, g, 0, 1, 2, 3, 40.0, 1.0, true)
        };
        let flat = at_height(0.0, &mut g);
        assert!(flat.abs() < 1.0e-12, "planar centre should cost nothing");
        assert!(at_height(0.2, &mut g) > flat);
        assert!(at_height(0.4, &mut g) > at_height(0.2, &mut g));
    }

    #[test]
    fn vdw_gradient_matches_finite_differences() {
        let pos = [0.0, 0.0, 0.0, 3.1, 0.7, -0.4];
        check_gradient(
            &pos,
            |p, g| vdw_pair(p, g, 0, 1, 3.8983, 0.0951),
            1.0e-6,
            "vdw",
        );
        check_translation_invariance(&pos, |p, g| vdw_pair(p, g, 0, 1, 3.8983, 0.0951), "vdw");
    }

    /// The 12-6 well must sit at R⁰ with depth D⁰; a factor-of-2^(1/6) slip here is a classic
    /// Lennard-Jones convention error (σ versus R_min) and would shift every contact distance.
    #[test]
    fn vdw_minimum_is_at_r0_with_depth_d0() {
        let (r0, d0) = (3.8983, 0.0951);
        let mut g = vec![0.0; 6];
        let energy_at = |r: f64, g: &mut Vec<f64>| {
            g.iter_mut().for_each(|x| *x = 0.0);
            vdw_pair(&[0.0, 0.0, 0.0, r, 0.0, 0.0], g, 0, 1, r0, d0)
        };
        let at_min = energy_at(r0, &mut g);
        assert!((at_min + d0).abs() < 1.0e-9, "well depth is {at_min}, want {}", -d0);
        assert!(energy_at(r0 * 0.9, &mut g) > at_min);
        assert!(energy_at(r0 * 1.1, &mut g) > at_min);
        // And the force vanishes at the minimum.
        g.iter_mut().for_each(|x| *x = 0.0);
        vdw_pair(&[0.0, 0.0, 0.0, r0, 0.0, 0.0], &mut g, 0, 1, r0, d0);
        assert!(g[0].abs() < 1.0e-9);
    }

    #[test]
    fn hbond_gradient_matches_finite_differences() {
        // Near-linear O-H···O, the binding arrangement.
        let pos = [
            0.0, 0.0, 0.0, // donor
            0.95, 0.05, 0.0, // hydrogen
            2.80, 0.12, 0.10, // acceptor
        ];
        check_gradient(
            &pos,
            |p, g| hbond(p, g, 0, 1, 2, super::super::params::D_HB, super::super::params::R_HB),
            1.0e-5,
            "hbond",
        );
        check_translation_invariance(
            &pos,
            |p, g| hbond(p, g, 0, 1, 2, super::super::params::D_HB, super::super::params::R_HB),
            "hbond",
        );
    }

    /// The 12-10 well must sit at R_hb with depth D_hb when the arrangement is linear.
    #[test]
    fn hbond_minimum_is_at_r_hb_when_linear() {
        use super::super::params::{D_HB, R_HB};
        let mut g = vec![0.0; 9];
        // Donor at origin, hydrogen on the axis, acceptor at `r` -- exactly linear, cos⁴θ = 1.
        let energy_at = |r: f64, g: &mut Vec<f64>| {
            g.iter_mut().for_each(|x| *x = 0.0);
            hbond(
                &[0.0, 0.0, 0.0, 0.95, 0.0, 0.0, r, 0.0, 0.0],
                g,
                0,
                1,
                2,
                D_HB,
                R_HB,
            )
        };
        let at_min = energy_at(R_HB, g.as_mut());
        assert!(
            (at_min + D_HB).abs() < 1.0e-9,
            "well depth is {at_min}, want {}",
            -D_HB
        );
        assert!(energy_at(R_HB * 0.85, &mut g) > at_min);
        assert!(energy_at(R_HB * 1.20, &mut g) > at_min);
    }

    /// A bent D-H···A is not a hydrogen bond, and must contribute nothing rather than a
    /// spurious attraction.
    #[test]
    fn hbond_vanishes_when_bent_past_ninety_degrees() {
        use super::super::params::{D_HB, R_HB};
        let mut g = vec![0.0; 9];
        // Donor and acceptor on the same side of the hydrogen: θ_DHA < 90°.
        let pos = [0.0, 0.0, 0.0, 0.95, 0.0, 0.0, 0.5, 2.0, 0.0];
        assert_eq!(hbond(&pos, &mut g, 0, 1, 2, D_HB, R_HB), 0.0);
        assert!(g.iter().all(|x| *x == 0.0));
    }
}
