//! Save and restore the PES Sampling panel settings (RON), including the
//! loaded Hessians, so a sampling run can be repeated or shared exactly.
use super::state::Settings;

const HEADER: &str = "// Beavyr PES sampling settings, format 1\n";

pub fn to_text(settings: &Settings) -> Result<String, String> {
    let body = ron::ser::to_string_pretty(settings, ron::ser::PrettyConfig::default())
        .map_err(|e| e.to_string())?;
    Ok(format!("{HEADER}{body}\n"))
}

pub fn from_text(text: &str) -> Result<Settings, String> {
    let settings: Settings = ron::from_str(text)
        .map_err(|e| format!("Not a readable PES sampling settings file: {e}"))?;
    super::sampling::check_loaded(&settings)?;
    Ok(settings)
}
