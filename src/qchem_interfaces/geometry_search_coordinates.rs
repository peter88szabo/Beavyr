//! Validate reaction references and map user atom numbers past dummy removal.
use super::*;
use crate::optimizer::transition_state::{
    InternalReferenceKind as K, TransitionStateInitialHessian, TransitionStateReactionMode as R,
    WeightedInternalReference as W,
};
fn indices(text: &str, n: usize, atoms: &[String]) -> Result<Vec<usize>, String> {
    super::super::scan::atom_indices(text, n, atoms)
}
fn entries(text: &str) -> impl Iterator<Item = &str> {
    text.split(';').map(str::trim).filter(|s| !s.is_empty())
}
pub fn prepare_ts(o: &Options, atoms: &[String]) -> Result<TransitionStateOptions, String> {
    let mut ts = o.ts.clone();
    ts.reaction_mode = match o.reaction {
        Reaction::Lowest => R::Lowest,
        Reaction::Bond => {
            let a = indices(&o.reaction_atoms, 2, atoms)?;
            R::Bond([a[0], a[1]])
        }
        Reaction::Angle => {
            let a = indices(&o.reaction_atoms, 3, atoms)?;
            R::Angle([a[0], a[1], a[2]])
        }
        Reaction::Dihedral => {
            let a = indices(&o.reaction_atoms, 4, atoms)?;
            R::Dihedral([a[0], a[1], a[2], a[3]])
        }
        Reaction::Transfer => {
            let a = indices(&o.reaction_atoms, 3, atoms)?;
            if ts.coordinate_mode == TransitionStateCoordinateMode::Cartesian {
                R::Transfer([a[0], a[1], a[2]])
            } else {
                R::Coordinates(vec![W::bond([a[0], a[1]], 1.), W::bond([a[2], a[1]], -1.)])
            }
        }
        Reaction::Mixed => {
            let mut refs = vec![];
            for e in entries(&o.reaction_coordinates) {
                let fields: Vec<_> = e.split(':').collect();
                if !(2..=3).contains(&fields.len()) {
                    return Err("Use kind:atoms:weight, separated by semicolons.".into());
                }
                let (kind, n) = match fields[0].trim().to_lowercase().as_str() {
                    "bond" => (K::Bond, 2),
                    "angle" => (K::Angle, 3),
                    "dihedral" => (K::Dihedral, 4),
                    _ => {
                        return Err("Use bond, angle or dihedral for a reaction coordinate.".into())
                    }
                };
                let weight = if fields.len() == 3 {
                    fields[2]
                        .trim()
                        .parse::<f64>()
                        .map_err(|_| "Invalid reaction-coordinate weight")?
                } else {
                    1.
                };
                if !weight.is_finite() || weight == 0. {
                    return Err("Reaction-coordinate weights must be finite and nonzero.".into());
                }
                refs.push(W {
                    kind,
                    atoms: indices(fields[1], n, atoms)?,
                    weight,
                });
            }
            if refs.is_empty() {
                return Err("Add at least one reaction coordinate.".into());
            }
            R::Coordinates(refs)
        }
        Reaction::Direction => {
            let d = o
                .direction
                .split(|c: char| c.is_whitespace() || c == ',')
                .filter(|s| !s.is_empty())
                .map(|s| s.parse::<f64>().map_err(|_| "Invalid Cartesian direction"))
                .collect::<Result<Vec<_>, _>>()?;
            if d.len() != 3 * atoms.len()
                || d.iter().any(|v| !v.is_finite())
                || d.iter().all(|v| *v == 0.)
            {
                return Err("Give three finite direction components per atom, with at least one nonzero component.".into());
            }
            R::Direction(d)
        }
    };
    ts.extra_bonds.clear();
    for (text, n) in [
        (&o.extra_bonds, 2),
        (&o.extra_angles, 3),
        (&o.extra_dihedrals, 4),
    ] {
        for e in entries(text) {
            let a = indices(e, n, atoms)?;
            for pair in a.windows(2) {
                let mut b = [pair[0], pair[1]];
                b.sort();
                if !ts.extra_bonds.contains(&b) {
                    ts.extra_bonds.push(b);
                }
            }
        }
    }
    if ts.coordinate_mode == TransitionStateCoordinateMode::Cartesian
        && ts.ts_initial_hessian == TransitionStateInitialHessian::Model
    {
        return Err("Cartesian TS searches require a calculated initial Hessian.".into());
    }
    if ts.coordinate_mode == TransitionStateCoordinateMode::Cartesian
        && matches!(
            o.reaction,
            Reaction::Angle | Reaction::Dihedral | Reaction::Mixed
        )
    {
        return Err("Angle, dihedral and mixed references require internal coordinates.".into());
    }
    Ok(ts)
}
pub fn validate(o: &Options, job: JobType, atoms: &[String]) -> Result<(), String> {
    if job == JobType::TransitionState {
        let t = prepare_ts(o, atoms)?;
        let vals = [
            t.energy_tol,
            t.max_step,
            t.rms_step,
            t.max_gradient,
            t.rms_gradient,
            t.max_step_internal,
            t.best_fit_rms_tol,
            t.hessian_step,
            t.trust_radius_min,
            t.g_cutoff,
            t.redundant_penalty,
            t.ts_hessian_eigenvalue_floor,
            t.ts_model_negative_curvature,
        ];
        if t.max_cycles == 0
            || t.best_fit_iters == 0
            || vals.iter().any(|v| !v.is_finite() || *v <= 0.)
        {
            return Err("TS step limits and tolerances must be positive and finite.".into());
        }
        if !t.min_gradient_improvement.is_finite()
            || t.min_gradient_improvement < 0.
            || !t.adaptive_hessian_overlap_min.is_finite()
            || !(0.0..=1.0).contains(&t.adaptive_hessian_overlap_min)
            || !t.adaptive_hessian_trust_fraction.is_finite()
            || !(0.0..=1.0).contains(&t.adaptive_hessian_trust_fraction)
            || !t.adaptive_hessian_rho_min.is_finite()
            || !t.adaptive_hessian_rho_max.is_finite()
            || t.adaptive_hessian_rho_min < 0.
            || t.adaptive_hessian_rho_max <= t.adaptive_hessian_rho_min
            || !t.hessian_recalc_near_convergence_factor.is_finite()
            || t.hessian_recalc_near_convergence_factor < 1.
        {
            return Err("Invalid adaptive TS Hessian thresholds.".into());
        }
        if t.trust_radius
            .into_iter()
            .chain(t.trust_radius_max)
            .any(|v| !v.is_finite() || v < t.trust_radius_min)
        {
            return Err("Trust radii must be finite and at least the minimum trust radius.".into());
        }
        if matches!((t.trust_radius,t.trust_radius_max),(Some(a),Some(b)) if a>b) {
            return Err("The initial trust radius exceeds the maximum.".into());
        }
    } else {
        let m = &o.minimum;
        if m.max_cycles == 0
            || m.best_fit_iters == 0
            || [
                m.energy_tol,
                m.max_step,
                m.rms_step,
                m.max_gradient,
                m.rms_gradient,
                m.max_step_internal,
                m.best_fit_rms_tol,
                m.min_line_search_scale,
                m.g_cutoff,
                m.redundant_penalty,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.)
        {
            return Err(
                "Minimum-search step limits and tolerances must be positive and finite.".into(),
            );
        }
        if m.coordinates != Coordinates::Cartesian
            && matches!(m.algorithm, MinimumAlgorithm::Bb1 | MinimumAlgorithm::Bb2)
        {
            return Err("Two-point gradient methods require Cartesian coordinates.".into());
        }
    }
    Ok(())
}
pub fn remap_dummies(o: &Options, atoms: &[String]) -> Result<Options, String> {
    prepare_ts(o, atoms)?;
    let mut clean = o.clone();
    let remap = |text: &str, n| -> Result<String, String> {
        Ok(indices(text, n, atoms)?
            .iter()
            .map(|&i| {
                atoms[..=i]
                    .iter()
                    .filter(|a| !crate::molecule::is_dummy(a))
                    .count()
                    .to_string()
            })
            .collect::<Vec<_>>()
            .join(","))
    };
    let n = match o.reaction {
        Reaction::Bond => 2,
        Reaction::Angle | Reaction::Transfer => 3,
        Reaction::Dihedral => 4,
        _ => 0,
    };
    if n > 0 {
        clean.reaction_atoms = remap(&o.reaction_atoms, n)?;
    }
    if o.reaction == Reaction::Mixed {
        clean.reaction_coordinates = entries(&o.reaction_coordinates)
            .map(|e| {
                let f: Vec<_> = e.split(':').collect();
                let n = match f[0].trim().to_lowercase().as_str() {
                    "bond" => 2,
                    "angle" => 3,
                    _ => 4,
                };
                Ok(format!(
                    "{}:{}:{}",
                    f[0],
                    remap(f[1], n)?,
                    f.get(2).copied().unwrap_or("1")
                ))
            })
            .collect::<Result<Vec<_>, String>>()?
            .join(";");
    }
    if o.reaction == Reaction::Direction {
        if let R::Direction(d) = prepare_ts(o, atoms)?.reaction_mode {
            clean.direction = d
                .chunks_exact(3)
                .zip(atoms)
                .filter(|(_, a)| !crate::molecule::is_dummy(a))
                .flat_map(|(p, _)| p.iter().map(ToString::to_string))
                .collect::<Vec<_>>()
                .join(" ");
        }
    }
    for (out, text, n) in [
        (&mut clean.extra_bonds, &o.extra_bonds, 2),
        (&mut clean.extra_angles, &o.extra_angles, 3),
        (&mut clean.extra_dihedrals, &o.extra_dihedrals, 4),
    ] {
        *out = entries(text)
            .map(|e| remap(e, n))
            .collect::<Result<Vec<_>, _>>()?
            .join(";");
    }
    Ok(clean)
}
