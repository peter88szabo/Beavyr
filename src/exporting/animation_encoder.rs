//! GIF is self-contained; MP4 uses a user-selectable FFmpeg executable, without a shell.
use std::{
    fs::{self, File, OpenOptions},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    process::{Child, ChildStdin, Command, Stdio},
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc, Mutex,
    },
    time::{SystemTime, UNIX_EPOCH},
};
#[path = "gif_encoder.rs"]
mod gif_encoder;
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Gif,
    Mp4,
}
impl Format {
    pub fn ext(self) -> &'static str {
        match self {
            Self::Gif => "gif",
            Self::Mp4 => "mp4",
        }
    }
}
enum Stream {
    Gif(gif_encoder::Gif<BufWriter<File>>),
    Mp4(ChildStdin),
}
pub struct Encoder {
    stream: Option<Stream>,
    pub child: Option<Arc<Mutex<Child>>>,
    pub cancel: Arc<AtomicBool>,
    temporary: PathBuf,
    destination: PathBuf,
    stderr: PathBuf,
    done: bool,
}
impl Encoder {
    pub fn new(
        destination: &Path,
        format: Format,
        width: u16,
        height: u16,
        fps: u16,
        repeat: bool,
        ffmpeg: &str,
    ) -> Result<Self, String> {
        if width < 2 || height < 2 || fps == 0 || fps > 100 || width > 4096 || height > 4096 {
            return Err(
                "Choose dimensions from 2 to 4096 pixels and a frame rate from 1 to 100.".into(),
            );
        }
        if format == Format::Mp4 && (width % 2 != 0 || height % 2 != 0) {
            return Err("MP4 width and height must both be even.".into());
        }
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let temporary = destination.with_file_name(format!(
            ".beavyr-animation-{}-{stamp}.tmp",
            std::process::id()
        ));
        let file = OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&temporary)
            .map_err(|e| e.to_string())?;
        let mut result = Self {
            stream: None,
            child: None,
            cancel: Arc::new(AtomicBool::new(false)),
            stderr: temporary.with_extension("stderr"),
            temporary,
            destination: destination.to_owned(),
            done: false,
        };
        result.stream = Some(match format {
            Format::Gif => Stream::Gif(
                gif_encoder::Gif::new(BufWriter::new(file), width, height, fps, repeat)
                    .map_err(|e| e.to_string())?,
            ),
            Format::Mp4 => {
                drop(file);
                let stderr = OpenOptions::new()
                    .create_new(true)
                    .write(true)
                    .open(&result.stderr)
                    .map_err(|e| e.to_string())?;
                let mut child=Command::new(ffmpeg).args(["-v","error","-nostdin","-f","rawvideo","-pixel_format","rgba","-video_size",&format!("{width}x{height}"),"-framerate",&fps.to_string(),"-i","pipe:0","-an","-c:v","libx264","-threads","4","-pix_fmt","yuv420p","-movflags","+faststart","-f","mp4","-y"]).arg(&result.temporary)
                    .stdin(Stdio::piped()).stdout(Stdio::null()).stderr(stderr).spawn()
                    .map_err(|e|format!("Cannot start FFmpeg at {ffmpeg}: {e}. Install FFmpeg or select its executable path. GIF export works without FFmpeg."))?;
                let stdin = child
                    .stdin
                    .take()
                    .ok_or("FFmpeg input pipe is unavailable")?;
                result.child = Some(Arc::new(Mutex::new(child)));
                Stream::Mp4(stdin)
            }
        });
        Ok(result)
    }
    pub fn frame(&mut self, rgba: &[u8]) -> Result<(), String> {
        match self.stream.as_mut().ok_or("Encoder is closed")? {
            Stream::Gif(gif) => gif.frame(rgba),
            Stream::Mp4(pipe) => pipe.write_all(rgba),
        }
        .map_err(|e| {
            format!(
                "Cannot encode frame: {e}\n{}",
                fs::read_to_string(&self.stderr).unwrap_or_default()
            )
        })
    }
    pub fn finish(&mut self) -> Result<(), String> {
        if let Some(stream) = self.stream.take() {
            match stream {
                Stream::Gif(gif) => gif.finish().map_err(|e| e.to_string())?,
                Stream::Mp4(pipe) => drop(pipe),
            }
        }
        if let Some(child) = &self.child {
            let status = loop {
                let status = child
                    .lock()
                    .map_err(|e| e.to_string())?
                    .try_wait()
                    .map_err(|e| e.to_string())?;
                if let Some(status) = status {
                    break status;
                }
                std::thread::sleep(std::time::Duration::from_millis(20));
            };
            if !status.success() {
                return Err(format!(
                    "FFmpeg failed ({status}):\n{}",
                    fs::read_to_string(&self.stderr).unwrap_or_default()
                ));
            }
        }
        if self.cancel.load(Ordering::SeqCst) {
            return Err("Export cancelled.".into());
        }
        fs::rename(&self.temporary, &self.destination)
            .map_err(|e| format!("Cannot save {}: {e}", self.destination.display()))?;
        self.done = true;
        Ok(())
    }
}
impl Drop for Encoder {
    fn drop(&mut self) {
        self.stream.take();
        if let Some(child) = &self.child {
            if let Ok(mut child) = child.lock() {
                let _ = child.kill();
                let _ = child.wait();
            }
        }
        if !self.done {
            let _ = fs::remove_file(&self.temporary);
        }
        let _ = fs::remove_file(&self.stderr);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn dir() -> PathBuf {
        crate::qchem_interfaces::xtb_optimize::create_xtb_run_dir(
            &std::env::temp_dir().join("beavyr-animation-tests"),
        )
        .unwrap()
    }
    #[test]
    fn cancellation_and_missing_ffmpeg_preserve_existing_file() {
        let dir = dir();
        let path = dir.join("existing.gif");
        fs::write(&path, "previous export").unwrap();
        {
            let mut encoder = Encoder::new(&path, Format::Gif, 2, 2, 20, true, "").unwrap();
            encoder.frame(&[255; 16]).unwrap();
        }
        assert_eq!(fs::read_to_string(&path).unwrap(), "previous export");
        let error = Encoder::new(
            &path,
            Format::Mp4,
            2,
            2,
            20,
            false,
            "/nonexistent/beavyr-test-ffmpeg",
        )
        .err()
        .unwrap();
        assert!(error.contains("FFmpeg") && error.contains("GIF export works"));
        assert_eq!(fs::read_to_string(&path).unwrap(), "previous export");
        assert_eq!(fs::read_dir(&dir).unwrap().count(), 1);
        fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn gif_export_finishes_and_rejects_invalid_dimensions() {
        let dir = dir();
        let path = dir.join("animation.gif");
        assert!(Encoder::new(&path, Format::Mp4, 3, 2, 20, false, "").is_err());
        let mut encoder = Encoder::new(&path, Format::Gif, 2, 2, 20, true, "").unwrap();
        encoder.frame(&[255; 16]).unwrap();
        encoder.frame(&[0; 16]).unwrap();
        encoder.finish().unwrap();
        drop(encoder);
        let bytes = fs::read(&path).unwrap();
        assert!(bytes.starts_with(b"GIF89a"));
        assert_eq!(bytes.last(), Some(&0x3b));
        // Retain a small fixture for an independent Pillow/ImageMagick decode check.
        let evidence = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/animation-test.gif");
        fs::write(evidence, bytes).unwrap();
        assert_eq!(fs::read_dir(&dir).unwrap().count(), 1);
        fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn cancel_before_finalization_never_publishes_the_temporary_file() {
        let dir = dir();
        let path = dir.join("original.gif");
        fs::write(&path, "original").unwrap();
        let mut encoder = Encoder::new(&path, Format::Gif, 2, 2, 20, true, "").unwrap();
        encoder.frame(&[255; 16]).unwrap();
        encoder.cancel.store(true, Ordering::SeqCst);
        assert!(encoder.finish().unwrap_err().contains("cancelled"));
        drop(encoder);
        assert_eq!(fs::read_to_string(&path).unwrap(), "original");
        assert_eq!(fs::read_dir(&dir).unwrap().count(), 1);
        fs::remove_dir_all(dir).unwrap();
    }
}
