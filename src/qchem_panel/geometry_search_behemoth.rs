//! Invoke Behemoth's own TS optimizer and retain its convergence evidence.
use super::*;
use crate::optimizer::transition_state::{
    TransitionStateCoordinateMode as C, TransitionStateInitialHessian as H,
    TransitionStateModeTracking as T,
};
use std::process::Command;
pub(super) fn command(
    binary: &Path,
    dir: &Path,
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
) -> Result<Command, String> {
    let o = &config.search;
    let t = &o.ts;
    let mut c = Command::new(binary);
    let defaults = crate::optimizer::transition_state::TransitionStateOptions::default();
    if [t.energy_tol, t.max_step, t.rms_step, t.max_gradient, t.rms_gradient]
        != [defaults.energy_tol, defaults.max_step, defaults.rms_step, defaults.max_gradient, defaults.rms_gradient] {
        crate::qchem_interfaces::method::tolerances::Criteria {
            energy:t.energy_tol, max_step:t.max_step, rms_step:t.rms_step,
            max_gradient:t.max_gradient, rms_gradient:t.rms_gradient,
        }.behemoth(&mut c);
    }
    c.current_dir(dir)
        .args([
            "--xyz",
            "geometry.xyz",
            "--method",
            "xtb",
            "--opt-ts",
            "--index-base",
            "1",
        ])
        .arg("--charge")
        .arg(charge.to_string())
        .arg("--multi")
        .arg(multiplicity.to_string())
        .arg("--nproc")
        .arg(config.nproc.max(1).to_string())
        .arg("--memory")
        .arg(config.memory_mb.to_string())
        .arg("--opt-coord")
        .arg(match t.coordinate_mode {
            C::Cartesian => "cartesian",
            C::Internal => "internal",
            C::PulayRedundant => "pulay-redundant-internal",
        })
        .arg("--ts-initial-hessian")
        .arg(if t.ts_initial_hessian == H::Exact {
            "exact"
        } else {
            "model"
        })
        .arg("--ts-mode-tracking")
        .arg(if t.mode_tracking_coordinates == T::Internal {
            "internal"
        } else {
            "cartesian"
        });
    for (flag, value) in [
        ("--ts-max-cycles", t.max_cycles.to_string()),
        (
            "--ts-hessian-recalc-interval",
            t.hessian_recalc_interval.to_string(),
        ),
        (
            "--ts-model-negative-curvature",
            t.ts_model_negative_curvature.to_string(),
        ),
        ("--ts-trust-radius-min", t.trust_radius_min.to_string()),
        ("--ts-max-rejected-steps", t.max_rejected_steps.to_string()),
    ] {
        c.arg(flag).arg(value);
    }
    for (flag, value) in [
        ("--ts-trust-radius", t.trust_radius),
        ("--ts-trust-radius-max", t.trust_radius_max),
    ] {
        if let Some(v) = value {
            c.arg(flag).arg(v.to_string());
        }
    }
    if !t.project_eckart {
        c.arg("--ts-no-eckart");
    }
    let flag = match o.reaction {
        search::Reaction::Lowest => None,
        search::Reaction::Bond => Some("--reaction-bond"),
        search::Reaction::Angle => Some("--reaction-angle"),
        search::Reaction::Dihedral => Some("--reaction-dihedral"),
        search::Reaction::Transfer => Some("--reaction-transfer"),
        search::Reaction::Mixed => None,
        search::Reaction::Direction => {
            return Err("Behemoth's CLI does not accept a Cartesian direction reference.".into())
        }
    };
    if let Some(flag) = flag {
        c.arg(flag).arg(&o.reaction_atoms);
    }
    for (flag, entries) in [
        (
            "--reaction-coordinate",
            if o.reaction == search::Reaction::Mixed {
                &o.reaction_coordinates
            } else {
                ""
            },
        ),
        ("--extra-bond", o.extra_bonds.as_str()),
        ("--extra-angle", o.extra_angles.as_str()),
        ("--extra-dihedral", o.extra_dihedrals.as_str()),
    ] {
        for entry in entries.split(';').map(str::trim).filter(|s| !s.is_empty()) {
            c.arg(flag).arg(entry);
        }
    }
    Ok(c)
}
#[allow(clippy::too_many_arguments)]
pub(super) fn run(
    binary: &Path,
    dir: &Path,
    atoms: &[String],
    coordinates: &[f64],
    charge: i32,
    multiplicity: i32,
    config: &MethodConfig,
    cancel: &AtomicBool,
    child: &Mutex<Option<Child>>,
) -> Result<QcRunOutput, String> {
    fs::write(
        dir.join("geometry.xyz"),
        format!(
            "{}\nTS guess\n{}",
            atoms.len(),
            orca_run::xyz_body(atoms, coordinates)
        ),
    )
    .map_err(|e| e.to_string())?;
    let mut c = command(binary, dir, charge, multiplicity, config)?;
    let out = dir.join("behemoth_ts.out");
    let err = dir.join("behemoth_ts.err");
    let cmd = dir.join("behemoth_ts_command.txt");
    fs::write(&cmd, format!("{c:?}\n")).map_err(|e| e.to_string())?;
    c.stdout(Stdio::from(
        fs::File::create(&out).map_err(|e| e.to_string())?,
    ))
    .stderr(Stdio::from(
        fs::File::create(&err).map_err(|e| e.to_string())?,
    ));
    check(cancel)?;
    let status = spawn_and_wait(c, "Behemoth TS search", cancel, child)?;
    check(cancel)?;
    let log = fs::read_to_string(&out).map_err(|e| e.to_string())?;
    let errors = fs::read_to_string(&err).unwrap_or_default();
    if !status.success() || log.contains("SCF did not converge") {
        return Err(format!(
            "Behemoth TS search failed ({status}).\n{errors}\nOutput: {}",
            out.display()
        ));
    }
    let final_xyz = fs::read_to_string(dir.join("transition_state_optimized.xyz"))
        .map_err(|e| format!("Missing Behemoth TS geometry: {e}"))?;
    let history = fs::read_to_string(dir.join("tsopt_trajectory.xyz")).unwrap_or_default();
    let result = read_result(&log, &final_xyz, &history, atoms)?;
    let mut r = output(
        QcProgram::Behemoth,
        JobType::TransitionState,
        dir,
        atoms,
        result,
    )?;
    r.log = format!("{}\n\n{log}", r.log);
    r.saveable.extend([
        ("Behemoth output", out),
        ("Behemoth errors", err),
        ("Behemoth command", cmd),
    ]);
    Ok(r)
}
fn value<'a>(text: &'a str, key: &str) -> Option<&'a str> {
    text.lines().find_map(|line| {
        let (name, val) = line.trim().split_once(':')?;
        (name.trim() == key).then_some(val.trim())
    })
}
fn geometry(frame: &crate::trajectory::TrajectoryFrame) -> Vec<f64> {
    frame
        .pos
        .iter()
        .flat_map(|p| p.to_array().map(|x| x as f64 / BOHR_TO_ANGSTROM))
        .collect()
}
pub(super) fn read_result(
    log: &str,
    final_xyz: &str,
    history: &str,
    atoms: &[String],
) -> Result<Outcome, String> {
    let summary = log
        .rsplit("Transition-state result")
        .next()
        .ok_or("Missing Behemoth TS result")?;
    let converged = value(summary, "converged")
        .and_then(|s| s.parse::<bool>().ok())
        .ok_or("Missing TS convergence status")?;
    let cycles = value(summary, "cycles")
        .and_then(|s| s.parse().ok())
        .ok_or("Missing TS step count")?;
    let energy = value(summary, "E(final)")
        .and_then(|s| s.split_whitespace().next())
        .and_then(|s| s.parse::<f64>().ok())
        .filter(|e| e.is_finite())
        .ok_or("Missing finite TS energy")?;
    let negative = value(summary, "negative_modes").and_then(|s| s.parse::<usize>().ok());
    let final_frames = crate::trajectory::parse_multi_xyz(final_xyz)?;
    if final_frames.len() != 1 || final_frames[0].atoms != atoms {
        return Err("Behemoth returned a TS geometry with different atoms.".into());
    }
    let mut frames = vec![];
    let mut energies = vec![];
    if !history.trim().is_empty() {
        let parsed = crate::trajectory::parse_multi_xyz(history)?;
        let lines: Vec<_> = history.lines().collect();
        for (i, f) in parsed.iter().enumerate() {
            if f.atoms.len() != atoms.len()
                || !f.atoms.iter().zip(atoms).all(|(a, b)| a == b || a == "X")
            {
                return Err("Unexpected atoms in Behemoth TS history.".into());
            }
            let e = lines
                .get(i * (atoms.len() + 2) + 1)
                .and_then(|s| s.split("Ene =").nth(1))
                .and_then(|s| s.split_whitespace().next())
                .and_then(|s| s.parse::<f64>().ok())
                .filter(|e| e.is_finite())
                .ok_or("Missing TS history energy")?;
            frames.push(geometry(f));
            energies.push(e);
        }
    }
    Ok(Outcome {
        x: geometry(&final_frames[0]),
        energy,
        converged: converged && negative == Some(1),
        cycles,
        message: format!(
            "{}\nNegative modes: {}",
            value(summary, "message").unwrap_or("Native Behemoth TS search finished."),
            negative.map(|n| n.to_string()).unwrap_or_else(|| {
                "not reported; update the Behemoth executable to verify this saddle".into()
            })
        ),
        frames,
        energies,
        hessian: None,
    })
}
