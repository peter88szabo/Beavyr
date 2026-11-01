//! Direct executable invocation, captured logs and cooperative process-tree cancellation.
use std::{
    path::Path,
    process::{Command, Stdio},
    sync::atomic::{AtomicBool, Ordering},
    time::Duration,
};
pub fn run(
    mut command: Command,
    directory: &Path,
    stem: &str,
    cancel: &AtomicBool,
) -> Result<String, String> {
    if cancel.load(Ordering::Relaxed) {
        return Err("Cancelled.".into());
    }
    std::fs::create_dir_all(directory).map_err(|e| e.to_string())?;
    let out = directory.join(format!("{stem}.out"));
    let err = directory.join(format!("{stem}.err"));
    command
        .current_dir(directory)
        .stdout(Stdio::from(
            std::fs::File::create(&out).map_err(|e| e.to_string())?,
        ))
        .stderr(Stdio::from(
            std::fs::File::create(&err).map_err(|e| e.to_string())?,
        ));
    #[cfg(unix)]
    {
        use std::os::unix::process::CommandExt;
        command.process_group(0);
    }
    let mut child = command.spawn().map_err(|e| {
        format!(
            "Cannot start calculation: {e}\nRun directory: {}",
            directory.display()
        )
    })?;
    let status = loop {
        if cancel.load(Ordering::Relaxed) {
            #[cfg(unix)]
            {
                let _ = Command::new("kill")
                    .args(["-KILL", "--", &format!("-{}", child.id())])
                    .status();
            }
            #[cfg(windows)]
            {
                let _ = Command::new("taskkill")
                    .args(["/PID", &child.id().to_string(), "/T", "/F"])
                    .status();
            }
            let _ = child.kill();
            let _ = child.wait();
            return Err("Cancelled.".into());
        }
        match child.try_wait() {
            Ok(Some(status)) => break status,
            Ok(None) => std::thread::sleep(Duration::from_millis(50)),
            Err(e) => {
                let _ = child.kill();
                let _ = child.wait();
                return Err(e.to_string());
            }
        }
    };
    let text = std::fs::read_to_string(&out).unwrap_or_default();
    if !status.success() {
        let errors = std::fs::read_to_string(&err).unwrap_or_default();
        let tail = text
            .lines()
            .rev()
            .take(35)
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect::<Vec<_>>()
            .join("\n");
        return Err(format!(
            "Calculation failed ({status}).\n{tail}\n{errors}\nRun directory: {}",
            directory.display()
        ));
    }
    Ok(text)
}
pub fn directory() -> Result<std::path::PathBuf, String> {
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(|e| e.to_string())?
        .as_nanos();
    let p = std::env::temp_dir().join(format!("beavyr-vibronic-{}-{stamp}", std::process::id()));
    std::fs::create_dir_all(&p).map_err(|e| e.to_string())?;
    Ok(p)
}
