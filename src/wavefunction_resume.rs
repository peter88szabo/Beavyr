//! Immutable input plans and verified restart of completed wavefunction frames.
use super::*;
pub const PLAN: &str = "wavefunction-plan.ron";
#[derive(Serialize, Deserialize)]
pub struct Plan {
    pub version: u32,
    pub engine: Engine,
    pub frames: Vec<(usize, TrajectoryFrame)>,
}
pub fn checksum(bytes: &[u8]) -> u64 {
    // Stable FNV-1a: detects accidental edits, not an authentication mechanism.
    bytes.iter().fold(0xcbf29ce484222325, |h, b| {
        (h ^ *b as u64).wrapping_mul(0x100000001b3)
    })
}
pub fn lock(root: &Path) -> Result<fs::File, String> {
    let lock = fs::File::options()
        .create(true)
        .truncate(false)
        .read(true)
        .write(true)
        .open(root.join("calculation.lock"))
        .map_err(|e| e.to_string())?;
    lock.try_lock().map_err(|_| {
        "This wavefunction calculation is already running in another worker or Beavyr window."
            .to_string()
    })?;
    Ok(lock)
}
pub fn save(
    root: &Path,
    engine: &Engine,
    frames: &[(usize, TrajectoryFrame)],
) -> Result<(), String> {
    use std::io::Write;
    let plan = Plan {
        version: 1,
        engine: engine.clone(),
        frames: frames.to_vec(),
    };
    let text = ron::ser::to_string(&plan).map_err(|e| e.to_string())?;
    let mut f=fs::File::options().write(true).create_new(true).open(root.join(PLAN)).map_err(|e|format!("Cannot create a new calculation plan: {e}. Use Resume for an existing calculation."))?;
    f.write_all(text.as_bytes())
        .and_then(|_| f.sync_all())
        .map_err(|e| e.to_string())?;
    fs::write(
        root.join("plan.checksum"),
        checksum(text.as_bytes()).to_string(),
    )
    .map_err(|e| e.to_string())
}
pub fn read_plan(root: &Path) -> Result<Plan, String> {
    let text=fs::read_to_string(root.join(PLAN)).map_err(|e|format!("Cannot read saved calculation plan: {e}. Calculations created before resume support cannot be resumed."))?;
    let saved = fs::read_to_string(root.join("plan.checksum"))
        .map_err(|e| format!("Cannot verify saved calculation settings: {e}"))?;
    if saved.trim() != checksum(text.as_bytes()).to_string() {
        return Err("The saved calculation plan was modified. Resume was refused to avoid mixing methods or geometries.".into());
    }
    let p: Plan = ron::from_str(&text).map_err(|e| e.to_string())?;
    if p.version != 1 || p.frames.is_empty() || p.frames.windows(2).any(|w| w[0].0 >= w[1].0) {
        return Err("Invalid wavefunction calculation plan.".into());
    }
    let geometry: Vec<_> = p.frames.iter().map(|(_, f)| f.clone()).collect();
    let clean = prepared(&geometry, 1, geometry.len(), 1)?;
    if clean
        .iter()
        .zip(&p.frames)
        .any(|((_, a), (_, b))| !geometry_matches(a, b))
    {
        return Err("Invalid geometry in saved calculation plan.".into());
    }
    Ok(p)
}
pub fn verified(root: &Path, plan: &Plan, control: &Control) -> Result<Vec<Frame>, String> {
    if !root.join(MANIFEST).exists() {
        return Ok(vec![]);
    }
    let series = read(&root.join(MANIFEST))?;
    let m = &series.manifest;
    let (charge, multiplicity) = plan.engine.state();
    if m.method != plan.engine.label()
        || m.charge != charge
        || m.multiplicity != multiplicity
        || m.requested_frames != plan.frames.len()
        || m.frames.len() > plan.frames.len()
    {
        return Err(
            "Saved results do not match the calculation plan. No files were overwritten.".into(),
        );
    }
    for (i, (done, (source, geometry))) in m.frames.iter().zip(&plan.frames).enumerate() {
        control.check()?;
        if done.source_index != *source
            || !geometry_matches(&done.geometry, geometry)
            || !done
                .geometry
                .metadata
                .energy_hartree
                .is_some_and(f64::is_finite)
        {
            return Err(format!(
                "Completed frame {} does not match the saved plan.",
                i + 1
            ));
        }
        let (_, text) = load(&series, i)?;
        if done.checksum != Some(checksum(text.as_bytes())) {
            return Err(format!("Wavefunction file for frame {} changed or has no verification checksum. No files were overwritten.",i+1));
        }
        control.progress.store(i + 1, Ordering::Relaxed);
    }
    Ok(m.frames.clone())
}
pub fn run(root: PathBuf, control: Control) -> Result<Series, String> {
    let _lock = lock(&root)?;
    let plan = read_plan(&root)?;
    *control.phase.lock().map_err(|e| e.to_string())? = "Verifying completed frames".into();
    let completed = verified(&root, &plan, &control)?;
    control.reused.store(completed.len(), Ordering::Relaxed);
    calculate_inner(plan.engine, plan.frames, root, control, completed)
}
pub fn duration(seconds: f64) -> String {
    let s = seconds.max(0.) as u64;
    format!("{}:{:02}:{:02}", s / 3600, (s / 60) % 60, s % 60)
}
