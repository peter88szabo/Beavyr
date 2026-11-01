//! Reversible dataset selection with original XYZ records kept for lossless subset export.
use crate::{pes_sampling::overlap::OverlapFilter, trajectory::TrajectoryFrame};
use std::sync::{
    atomic::{AtomicBool, AtomicUsize, Ordering},
    Arc,
};
#[derive(Clone)]
pub struct Dataset {
    pub source: String,
    pub frames: Arc<Vec<TrajectoryFrame>>,
    pub originals: Arc<Vec<String>>,
}
#[derive(Clone, Debug)]
pub struct Row {
    pub keep: bool,
    pub reason: String,
    pub duplicate_of: Option<usize>,
    pub nearest_rmsd: Option<f64>,
}
impl Default for Row {
    fn default() -> Self {
        Self {
            keep: true,
            reason: String::new(),
            duplicate_of: None,
            nearest_rmsd: None,
        }
    }
}
#[derive(Clone)]
pub struct Options {
    pub overlap: OverlapFilter,
    pub duplicates: bool,
    pub rmsd: f64,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            overlap: Default::default(),
            duplicates: true,
            rmsd: 0.05,
        }
    }
}
impl Dataset {
    pub fn from_xyz(source: String, text: &str) -> Result<Self, String> {
        let frames = crate::trajectory::parse_multi_xyz(text)?;
        // A strict standard XYZ record retains original coordinate precision and comments.
        let lines: Vec<_> = text.lines().collect();
        let mut cursor = 0;
        let mut originals = vec![];
        while cursor < lines.len() {
            if lines[cursor].trim().is_empty() {
                cursor += 1;
                continue;
            }
            let n = lines[cursor]
                .trim()
                .parse::<usize>()
                .map_err(|e| e.to_string())?;
            let end = cursor.checked_add(n + 2).ok_or("XYZ frame too large")?;
            if end > lines.len() {
                return Err("Incomplete XYZ dataset record.".into());
            }
            originals.push(format!("{}\n", lines[cursor..end].join("\n")));
            cursor = end;
        }
        if originals.len() != frames.len() {
            return Err(
                "Dataset requires standard XYZ records with a comment line per frame.".into(),
            );
        }
        Ok(Self {
            source,
            frames: Arc::new(frames),
            originals: Arc::new(originals),
        })
    }
    pub fn from_trajectory(traj: &crate::trajectory::TrajectoryState) -> Result<Self, String> {
        if traj.frames.is_empty() {
            return Err("Load a trajectory or generate PES samples first.".into());
        }
        if let Some(path) = &traj.current_file {
            if let Ok(text) = std::fs::read_to_string(path) {
                if let Ok(data) = Self::from_xyz(path.display().to_string(), &text) {
                    if crate::trajectory_data::fingerprint(&data.frames)
                        == crate::trajectory_data::fingerprint(&traj.frames)
                    {
                        return Ok(data);
                    }
                }
            }
        }
        let originals = (0..traj.frames.len())
            .map(|i| crate::trajectory_data::xyz(&traj.frames, [i]))
            .collect();
        Ok(Self {
            source: traj
                .current_file
                .as_ref()
                .map(|p| p.display().to_string())
                .unwrap_or_else(|| "Viewer trajectory".into()),
            frames: Arc::new(traj.frames.clone()),
            originals: Arc::new(originals),
        })
    }
    pub fn selected_xyz(&self, rows: &[Row]) -> String {
        let mut text = String::new();
        for (i, row) in rows.iter().enumerate() {
            if row.keep {
                let mut lines = self.originals[i].lines();
                text.push_str(lines.next().unwrap());
                text.push('\n');
                text.push_str(lines.next().unwrap_or(""));
                text.push_str(&format!(" source_frame={}\n", i + 1));
                for line in lines {
                    text.push_str(line);
                    text.push('\n');
                }
            }
        }
        text
    }
    pub fn manifest(&self, rows: &[Row]) -> String {
        fn quote(s: &str) -> String {
            format!("\"{}\"", s.replace('"', "\"\""))
        }
        let mut out=String::from("source,source_frame,keep,reason,duplicate_of_frame,nearest_prior_kept_rmsd_A,original_comment\n");
        for (i, r) in rows.iter().enumerate() {
            out.push_str(&format!(
                "{},{},{},{},{},{},{}\n",
                quote(&self.source),
                i + 1,
                r.keep,
                quote(&r.reason),
                r.duplicate_of
                    .map(|i| (i + 1).to_string())
                    .unwrap_or_default(),
                r.nearest_rmsd
                    .map(|v| format!("{v:.10}"))
                    .unwrap_or_default(),
                quote(&self.frames[i].metadata.comment)
            ));
        }
        out
    }
}
fn check(cancel: &AtomicBool) -> Result<(), String> {
    if cancel.load(Ordering::Relaxed) {
        Err("Cancelled; the previous selection is unchanged.".into())
    } else {
        Ok(())
    }
}
fn radii(frame: &TrajectoryFrame) -> Vec<f64> {
    let center = crate::rmsd::kabsch::centroid(&frame.pos);
    frame
        .pos
        .iter()
        .map(|p| (*p - center).length() as f64)
        .collect()
}
fn lower_bound(a: &[f64], b: &[f64]) -> f64 {
    (a.iter().zip(b).map(|(a, b)| (a - b).powi(2)).sum::<f64>() / a.len().max(1) as f64).sqrt()
}
fn distance(a: &TrajectoryFrame, b: &TrajectoryFrame) -> Option<f64> {
    if a.atoms != b.atoms {
        return None;
    }
    let s = crate::molecular_alignment::superpose(&a.pos, &b.pos).ok()?;
    crate::rmsd::kabsch::rmsd(&s.apply_all(&a.pos), &b.pos).ok()
}
pub fn analyze(
    data: &Dataset,
    options: &Options,
    cancel: &AtomicBool,
    progress: &AtomicUsize,
) -> Result<Vec<Row>, String> {
    options.overlap.validate()?;
    if !options.rmsd.is_finite() || options.rmsd < 0. {
        return Err("RMSD threshold must be finite and nonnegative.".into());
    }
    let mut rows = vec![Row::default(); data.frames.len()];
    let mut kept = Vec::<usize>::new();
    let radial: Vec<_> = data.frames.iter().map(radii).collect();
    for (i, f) in data.frames.iter().enumerate() {
        check(cancel)?;
        progress.store(i, Ordering::Relaxed);
        if f.pos.len() != f.atoms.len() || f.pos.is_empty() || f.pos.iter().any(|p| !p.is_finite())
        {
            rows[i].keep = false;
            rows[i].reason = "Invalid geometry".into();
            continue;
        }
        let coords: Vec<f64> = f
            .pos
            .iter()
            .flat_map(|p| p.to_array().map(|x| x as f64))
            .collect();
        if !options.overlap.accepts(&f.atoms, &coords) {
            rows[i].keep = false;
            rows[i].reason = "Atomic overlap".into();
            continue;
        }
        let mut nearest = f64::INFINITY;
        let mut nearest_index = None;
        if options.duplicates {
            for &j in &kept {
                check(cancel)?;
                if f.atoms != data.frames[j].atoms || lower_bound(&radial[i], &radial[j]) > nearest
                {
                    continue;
                }
                if let Some(d) = distance(f, &data.frames[j]) {
                    if d < nearest {
                        nearest = d;
                        nearest_index = Some(j);
                    }
                }
            }
        }
        if nearest.is_finite() {
            rows[i].nearest_rmsd = Some(nearest);
        }
        if options.duplicates && nearest <= options.rmsd {
            rows[i].keep = false;
            rows[i].reason = "RMSD duplicate".into();
            rows[i].duplicate_of = nearest_index;
        } else {
            kept.push(i);
        }
    }
    progress.store(data.frames.len(), Ordering::Relaxed);
    Ok(rows)
}
pub fn diverse(
    data: &Dataset,
    rows: &[Row],
    count: usize,
    cancel: &AtomicBool,
    progress: &AtomicUsize,
) -> Result<Vec<Row>, String> {
    let candidates: Vec<_> = rows
        .iter()
        .enumerate()
        .filter(|(_, r)| r.keep)
        .map(|(i, _)| i)
        .collect();
    if count == 0 || candidates.is_empty() {
        return Err("Choose at least one structure from a nonempty selection.".into());
    }
    let count = count.min(candidates.len());
    let mut selected = vec![candidates[0]];
    let mut nearest = vec![f64::INFINITY; rows.len()];
    while selected.len() < count {
        check(cancel)?;
        let last = *selected.last().unwrap();
        for &i in &candidates {
            check(cancel)?;
            if selected.contains(&i) {
                nearest[i] = f64::NEG_INFINITY;
                continue;
            }
            if let Some(d) = distance(&data.frames[i], &data.frames[last]) {
                nearest[i] = nearest[i].min(d);
            }
        }
        let next = *candidates
            .iter()
            .filter(|i| !selected.contains(i))
            .max_by(|&&a, &&b| nearest[a].total_cmp(&nearest[b]).then_with(|| b.cmp(&a)))
            .unwrap();
        selected.push(next);
        progress.store(selected.len(), Ordering::Relaxed);
    }
    let mut result = rows.to_vec();
    for &i in &candidates {
        if !selected.contains(&i) {
            result[i].keep = false;
            result[i].reason = "Outside diverse subset".into();
        }
    }
    Ok(result)
}
pub fn save(
    data: &Dataset,
    rows: &[Row],
    parent: &std::path::Path,
) -> Result<std::path::PathBuf, String> {
    if !rows.iter().any(|r| r.keep) {
        return Err("No frames are selected.".into());
    }
    let path = crate::qchem_interfaces::xtb_optimize::create_xtb_run_dir(parent)?;
    let result: Result<(), String> = (|| {
        std::fs::write(path.join("selected.xyz"), data.selected_xyz(rows))
            .map_err(|e| e.to_string())?;
        std::fs::write(path.join("selection.csv"), data.manifest(rows))
            .map_err(|e| e.to_string())?;
        Ok(())
    })();
    if let Err(e) = result {
        let _ = std::fs::write(path.join("INCOMPLETE_EXPORT.txt"), &e);
        return Err(format!("{e}\nIncomplete export: {}", path.display()));
    }
    Ok(path)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn duplicates_ignore_translation_and_rotation_but_overlaps_are_rejected() {
        let d=Dataset::from_xyz("test".into(),"2\noriginal\nH 0 0 0\nH 1 0 0\n2\nrotated\nH 2 2 0\nH 2 3 0\n2\nclash\nH 0 0 0\nH .1 0 0\n2\ndiverse\nH 0 0 0\nH 2 0 0\n").unwrap();
        let r = analyze(
            &d,
            &Options::default(),
            &AtomicBool::new(false),
            &AtomicUsize::new(0),
        )
        .unwrap();
        assert_eq!(
            r.iter().map(|r| r.keep).collect::<Vec<_>>(),
            vec![true, false, false, true]
        );
        assert_eq!(r[1].duplicate_of, Some(0));
        assert_eq!(r[2].reason, "Atomic overlap");
        let selected = diverse(&d, &r, 1, &AtomicBool::new(false), &AtomicUsize::new(0)).unwrap();
        assert_eq!(selected.iter().filter(|r| r.keep).count(), 1);
    }
    #[test]
    fn subset_export_preserves_source_precision_and_metadata() {
        let d = Dataset::from_xyz(
            "seeded".into(),
            "1\nPES_sample=1 seed=7\nC 1.1234567890123 0 0\n",
        )
        .unwrap();
        let r = vec![Row::default()];
        let xyz = d.selected_xyz(&r);
        assert!(xyz.contains("1.1234567890123"));
        assert!(xyz.contains("PES_sample=1 seed=7 source_frame=1"));
        assert!(d.manifest(&r).contains("seeded"));
    }
    #[test]
    fn cancellation_does_not_return_a_partial_selection() {
        let d = Dataset::from_xyz("test".into(), "1\nframe\nH 0 0 0\n").unwrap();
        assert!(analyze(
            &d,
            &Options::default(),
            &AtomicBool::new(true),
            &AtomicUsize::new(0)
        )
        .unwrap_err()
        .contains("previous selection"));
    }
}
