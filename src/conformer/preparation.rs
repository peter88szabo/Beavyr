//! Keep constraints attached to the same real atoms after placeholders are removed.

use super::{afir::ActiveBond, constraints::Constraint, ConformerRun};
use crate::optimizer::constrained::ConstraintCoordinate;

fn atom(mapping: &[Option<usize>], index: usize) -> Result<usize, String> {
    mapping.get(index).copied().flatten().ok_or_else(|| {
        "A constraint refers to a removed dummy atom or an invalid atom. Select real atoms before starting the search.".into()
    })
}

pub(crate) fn remap_constraints(
    constraints: &[Constraint],
    mapping: &[Option<usize>],
) -> Result<Vec<Constraint>, String> {
    constraints
        .iter()
        .map(|c| {
            let coordinate = match c.coordinate {
                ConstraintCoordinate::Bond([a, b]) => {
                    ConstraintCoordinate::Bond([atom(mapping, a)?, atom(mapping, b)?])
                }
                ConstraintCoordinate::Angle([a, b, c]) => ConstraintCoordinate::Angle([
                    atom(mapping, a)?,
                    atom(mapping, b)?,
                    atom(mapping, c)?,
                ]),
                ConstraintCoordinate::Dihedral([a, b, c, d]) => ConstraintCoordinate::Dihedral([
                    atom(mapping, a)?,
                    atom(mapping, b)?,
                    atom(mapping, c)?,
                    atom(mapping, d)?,
                ]),
            };
            Ok(Constraint {
                coordinate,
                value: c.value,
            })
        })
        .collect()
}

pub(crate) fn remap_active(
    active: &[ActiveBond],
    mapping: &[Option<usize>],
) -> Result<Vec<ActiveBond>, String> {
    if active
        .iter()
        .any(|bond| bond.a >= mapping.len() || bond.b >= mapping.len())
    {
        return Err("an active bond names an atom that is not there".into());
    }
    active
        .iter()
        .map(|bond| {
            Ok(ActiveBond {
                a: atom(mapping, bond.a)?,
                b: atom(mapping, bond.b)?,
            })
        })
        .collect()
}

pub(crate) fn remap_display_settings(run: &mut ConformerRun, mapping: &[Option<usize>]) {
    match remap_constraints(&run.settings.constraints, mapping) {
        Ok(constraints) => run.settings.constraints = constraints,
        Err(error) => {
            run.settings.constraints = run
                .settings
                .constraints
                .iter()
                .filter_map(|c| remap_constraints(std::slice::from_ref(c), mapping).ok())
                .flatten()
                .collect();
            run.error = Some(error);
        }
    }
    match remap_active(&run.trail.active, mapping) {
        Ok(active) => run.trail.active = active,
        Err(error) => {
            run.trail.active = run
                .trail
                .active
                .iter()
                .filter_map(|a| remap_active(std::slice::from_ref(a), mapping).ok())
                .flatten()
                .collect();
            run.trail_error = Some(error);
        }
    }
    let pending = |i: usize| {
        i.checked_sub(1)
            .and_then(|i| mapping.get(i).copied().flatten())
            .map_or(1, |i| i + 1)
    };
    run.pending_active = (pending(run.pending_active.0), pending(run.pending_active.1));
}
