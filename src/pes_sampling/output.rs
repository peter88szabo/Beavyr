use crate::{
    molecule::Molecule,
    trajectory::{TrajectoryFrame, TrajectoryState},
};
use bevy::prelude::Vec3;
use std::{
    io::{BufWriter, Write},
    path::Path,
};
#[derive(Clone, Debug)]
pub struct Sample {
    pub coordinates: Vec<f64>,
    pub distance: Option<f64>,
    pub distance_index: Option<usize>,
    pub sample_at_distance: Option<usize>,
}
#[derive(Clone, Debug)]
pub struct SamplingOutput {
    pub atoms: Vec<String>,
    pub samples: Vec<Sample>,
    pub seed: u64,
    pub rejected: usize,
    pub warnings: Vec<String>,
    pub partner_sizes: Vec<usize>,
}
impl SamplingOutput {
    pub fn sample_comment(&self, i: usize) -> String {
        let s = &self.samples[i];
        let mut text = format!("PES_sample={} seed={} phases=uniform_0_2pi", i+1, self.seed);
        if let Some(d) = s.distance {
            text.push_str(&format!(" COM_distance_A={d:.10} distance_index={} sample_at_distance={} partner_A_atoms={}", s.distance_index.unwrap_or(0), s.sample_at_distance.unwrap_or(0), self.partner_sizes.first().copied().unwrap_or(0)));
        }
        text
    }

    pub fn write_xyz(&self, mut writer: impl Write) -> Result<(), String> {
        for (i, s) in self.samples.iter().enumerate() {
            writeln!(writer, "{}", self.atoms.len()).map_err(|e| e.to_string())?;
            writeln!(writer, "{}", self.sample_comment(i)).map_err(|e| e.to_string())?;
            for (a, p) in self.atoms.iter().zip(s.coordinates.chunks_exact(3)) {
                writeln!(
                    writer,
                    "{a:<3} {:18.10} {:18.10} {:18.10}",
                    p[0], p[1], p[2]
                )
                .map_err(|e| e.to_string())?;
            }
        }
        Ok(())
    }
    pub fn save(&self, path: &Path) -> Result<(), String> {
        let suffix = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or_default();
        let temporary = path.with_file_name(format!(
            ".{}.pes-{}-{suffix}.tmp",
            path.file_name().unwrap_or_default().to_string_lossy(),
            std::process::id()
        ));
        let result = (|| {
            let file = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)
                .map_err(|e| e.to_string())?;
            let mut writer = BufWriter::new(file);
            self.write_xyz(&mut writer)?;
            writer.flush().map_err(|e| e.to_string())?;
            writer.get_ref().sync_all().map_err(|e| e.to_string())?;
            drop(writer);
            std::fs::rename(&temporary, path).map_err(|e| e.to_string())
        })();
        if result.is_err() {
            let _ = std::fs::remove_file(temporary);
        }
        result
    }
    pub fn load(&self, mol: &mut Molecule, traj: &mut TrajectoryState) -> bool {
        let frames = self
            .samples
            .iter()
            .enumerate()
            .map(|(i,s)| TrajectoryFrame { metadata: crate::trajectory_data::Metadata::comment(self.sample_comment(i)),
                atoms: self.atoms.clone(),
                pos: s
                    .coordinates
                    .chunks_exact(3)
                    .map(|p| Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32))
                    .collect(),
            })
            .collect();
        crate::ts_generation::show_frames(frames, mol, traj)
    }
}
