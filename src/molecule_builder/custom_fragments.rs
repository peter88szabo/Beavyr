//! Named, user-owned fragments. Placement goes through the existing builder.
use crate::{
    molecule::{parse_xyz_angstrom, Molecule},
    ui_style::ResponseExt,
};
use bevy::prelude::*;
use bevy_egui::egui;
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeSet,
    fs,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

const PREFIX: &str = "Saved: ";

#[derive(Clone, Serialize, Deserialize)]
pub struct Fragment {
    version: u32,
    pub name: String,
    pub xyz: String,
    #[serde(skip)]
    path: PathBuf,
}

impl Fragment {
    pub fn key(&self) -> String {
        format!("{PREFIX}{}", self.name)
    }

    pub(crate) fn from_selection(
        name: &str,
        mol: &Molecule,
        indices: &[usize],
        anchor: usize,
    ) -> Result<Self, String> {
        let name = name.trim();
        validate_name(name)?;
        if indices.is_empty() || !indices.contains(&anchor) {
            return Err("Choose an attachment atom included in the fragment.".into());
        }
        if mol.atoms.len() != mol.pos.len() || indices.iter().any(|&i| i >= mol.atoms.len()) {
            return Err("The selection no longer belongs to this molecule.".into());
        }
        let mut ordered = vec![anchor];
        ordered.extend(indices.iter().copied().filter(|i| *i != anchor));
        let mut seen = BTreeSet::new();
        let mut xyz = String::new();
        for i in ordered {
            if !seen.insert(i) {
                continue;
            }
            let p = mol.pos[i] - mol.pos[anchor];
            if !p.is_finite()
                || mol.atoms[i].is_empty()
                || mol.atoms[i].chars().any(char::is_whitespace)
            {
                return Err("The fragment contains invalid coordinates or elements.".into());
            }
            xyz.push_str(&format!(
                "{} {:.9} {:.9} {:.9}\n",
                mol.atoms[i], p.x, p.y, p.z
            ));
        }
        Ok(Self {
            version: 1,
            name: name.into(),
            xyz,
            path: PathBuf::new(),
        })
    }

    fn validate(&self) -> Result<(), String> {
        validate_name(&self.name)?;
        let (n, symbols, coords) = parse_xyz_angstrom(&self.xyz);
        if self.version != 1
            || n == 0
            || symbols.len() != n
            || coords.len() != n
            || self.xyz.lines().filter(|s| !s.trim().is_empty()).count() != n
            || coords
                .iter()
                .any(|p| p.len() != 3 || p.iter().any(|x| !x.is_finite()))
        {
            return Err("Invalid saved fragment".into());
        }
        Ok(())
    }
}

fn validate_name(name: &str) -> Result<(), String> {
    if name.trim().is_empty() {
        return Err("Give the fragment a name.".into());
    }
    if name.chars().count() > 80 || name.chars().any(char::is_control) {
        return Err("Use a name of at most 80 characters, without line breaks.".into());
    }
    Ok(())
}

fn directory() -> Option<PathBuf> {
    #[cfg(target_os = "windows")]
    if let Some(base) = std::env::var_os("APPDATA") {
        return Some(PathBuf::from(base).join("Beavyr/fragments"));
    }
    #[cfg(target_os = "macos")]
    if let Some(base) = std::env::var_os("HOME") {
        return Some(PathBuf::from(base).join("Library/Application Support/Beavyr/fragments"));
    }
    std::env::var_os("XDG_CONFIG_HOME")
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|p| PathBuf::from(p).join(".config")))
        .map(|base| base.join("beavyr/fragments"))
}

#[derive(Clone)]
pub struct Library {
    pub entries: Vec<Fragment>,
    directory: Option<PathBuf>,
    pub message: Option<String>,
    removed: Option<Fragment>,
}

impl Default for Library {
    fn default() -> Self {
        Self::at(directory())
    }
}

impl Library {
    fn at(directory: Option<PathBuf>) -> Self {
        let mut library = Self {
            entries: vec![],
            directory,
            message: None,
            removed: None,
        };
        library.refresh();
        library
    }
    pub fn find(&self, key: &str) -> Option<&Fragment> {
        let name = key.strip_prefix(PREFIX)?;
        self.entries.iter().find(|f| f.name == name)
    }
    fn check_name(&self, name: &str) -> Result<(), String> {
        validate_name(name)?;
        if self
            .entries
            .iter()
            .any(|f| f.name.to_lowercase() == name.trim().to_lowercase())
        {
            return Err("This name is already in My fragments. Choose another name.".into());
        }
        Ok(())
    }
    fn refresh(&mut self) {
        let Some(directory) = &self.directory else {
            self.message = Some("No user configuration directory is available.".into());
            return;
        };
        let files = match fs::read_dir(directory) {
            Ok(files) => files,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
                self.entries.clear();
                return;
            }
            Err(e) => {
                self.message = Some(format!("Cannot read My fragments: {e}"));
                return;
            }
        };
        let mut entries = Vec::new();
        let mut failures = 0;
        for file in files {
            let Ok(file) = file else {
                failures += 1;
                continue;
            };
            let path = file.path();
            if path.extension().is_none_or(|e| e != "beavyr-fragment") {
                continue;
            }
            let read = fs::read_to_string(&path)
                .map_err(|e| e.to_string())
                .and_then(|text| ron::from_str::<Fragment>(&text).map_err(|e| e.to_string()))
                .and_then(|mut f| {
                    f.validate()?;
                    f.path = path;
                    Ok(f)
                });
            match read {
                Ok(fragment) => entries.push(fragment),
                Err(_) => failures += 1,
            }
        }
        entries.sort_by_key(|f| f.name.to_lowercase());
        self.entries = entries;
        self.message =
            (failures > 0).then(|| format!("{failures} saved fragment file(s) could not be read."));
    }
    fn save(&mut self, mut fragment: Fragment) -> Result<String, String> {
        self.refresh();
        self.check_name(&fragment.name)?;
        fragment.validate()?;
        let directory = self
            .directory
            .as_ref()
            .ok_or("No fragment directory is available")?;
        fs::create_dir_all(directory).map_err(|e| e.to_string())?;
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_err(|e| e.to_string())?
            .as_nanos();
        fragment.path = directory.join(format!(
            "fragment-{stamp}-{}.beavyr-fragment",
            std::process::id()
        ));
        let text =
            ron::ser::to_string_pretty(&fragment, Default::default()).map_err(|e| e.to_string())?;
        crate::project::atomic_write(&fragment.path, text.as_bytes()).map_err(|e| e.to_string())?;
        let key = fragment.key();
        self.message = Some(format!("Saved {}", fragment.name));
        self.entries.push(fragment);
        self.entries.sort_by_key(|f| f.name.to_lowercase());
        Ok(key)
    }
    fn remove(&mut self, index: usize) -> Result<(), String> {
        let fragment = self.entries.get(index).ok_or("Fragment no longer exists")?;
        fs::rename(
            &fragment.path,
            fragment.path.with_extension("removed-fragment"),
        )
        .map_err(|e| e.to_string())?;
        self.removed = Some(self.entries.remove(index));
        self.message = Some("Removed from My fragments. Undo is available.".into());
        Ok(())
    }
    fn undo_remove(&mut self) -> Result<(), String> {
        let fragment = self.removed.as_ref().ok_or("Nothing to restore")?;
        self.check_name(&fragment.name)?;
        fs::rename(
            fragment.path.with_extension("removed-fragment"),
            &fragment.path,
        )
        .map_err(|e| e.to_string())?;
        self.removed = None;
        self.refresh();
        Ok(())
    }
}

#[derive(Clone)]
pub struct SaveForm {
    name: String,
    whole: bool,
    indices: String,
    anchor: usize,
}
impl Default for SaveForm {
    fn default() -> Self {
        Self {
            name: String::new(),
            whole: true,
            indices: String::new(),
            anchor: 1,
        }
    }
}

fn selection(text: &str, count: usize) -> Result<Vec<usize>, String> {
    let mut indices = BTreeSet::new();
    for token in text
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|s| !s.is_empty())
    {
        let (start, end) = token.split_once('-').unwrap_or((token, token));
        let start = start
            .parse::<usize>()
            .map_err(|_| "Use atom numbers or ranges, such as 1-6, 9.")?;
        let end = end
            .parse::<usize>()
            .map_err(|_| "Use atom numbers or ranges, such as 1-6, 9.")?;
        if start == 0 || end < start || end > count {
            return Err(format!("Atom numbers must be between 1 and {count}."));
        }
        indices.extend((start..=end).map(|n| n - 1));
    }
    if indices.is_empty() {
        return Err("Select atoms to save.".into());
    }
    Ok(indices.into_iter().collect())
}

pub fn show(
    ui: &mut egui::Ui,
    library: &mut Library,
    form: &mut SaveForm,
    chosen: &mut String,
    mol: &Molecule,
    selected: Option<usize>,
    moving: Option<&[usize]>,
    search: &str,
) {
    ui.horizontal(|ui| {
        if ui
            .small_button("↻")
            .on_hover_text("Refresh saved fragments")
            .clicked_once()
        {
            library.refresh();
        }
        if library.removed.is_some() && ui.small_button("Undo removal").clicked_once() {
            if let Err(e) = library.undo_remove() {
                library.message = Some(e);
            }
        }
    });
    let mut remove = None;
    egui::ScrollArea::vertical().id_salt("saved_fragment_buttons").max_height(112.0).show(ui, |ui| {
        ui.horizontal_wrapped(|ui| {
        for (i, fragment) in library.entries.iter().enumerate() {
            if !fragment.name.to_lowercase().contains(&search.trim().to_lowercase()) { continue; }
            let response = ui.add(egui::Button::new(&fragment.name).selected(*chosen == fragment.key()))
                .on_hover_text("Select this fragment, then Place, Add or Replace. Right-click for library actions.");
            if response.clicked_once() { *chosen = fragment.key(); }
            response.context_menu(|ui| {
                if ui.button("Remove from My fragments").clicked_once() { remove = Some(i); ui.close(); }
            });
        }
    });
    });
    if let Some(i) = remove {
        if let Err(e) = library.remove(i) {
            library.message = Some(e);
        }
    }
    egui::CollapsingHeader::new("Save molecule as fragment…").show(ui, |ui| {
        ui.add(egui::TextEdit::singleline(&mut form.name).hint_text("Fragment name (required)").desired_width(f32::INFINITY));
        ui.horizontal_wrapped(|ui| {
            ui.selectable_value(&mut form.whole, true, "Whole molecule");
            ui.selectable_value(&mut form.whole, false, "Selected atoms");
        });
        if !form.whole {
            ui.add(egui::TextEdit::singleline(&mut form.indices).hint_text("Atom numbers: 1-6, 9").desired_width(f32::INFINITY));
            ui.horizontal_wrapped(|ui| {
                if ui.add_enabled(selected.is_some(), egui::Button::new("Add picked atom")).clicked_once() {
                    if let Some(i) = selected { if !form.indices.is_empty() { form.indices.push_str(", "); } form.indices.push_str(&(i + 1).to_string()); }
                }
                if ui.add_enabled(moving.is_some_and(|m| !m.is_empty()), egui::Button::new("Use edited fragment")).clicked_once() {
                    if let Some(moving) = moving { form.indices = moving.iter().map(|i| (i+1).to_string()).collect::<Vec<_>>().join(", "); }
                }
            });
        }
        let indices = if form.whole { Ok((0..mol.atoms.len()).collect()) } else { selection(&form.indices, mol.atoms.len()) };
        ui.horizontal_wrapped(|ui| {
            ui.label("Attachment atom").on_hover_text("This atom connects to the host when the fragment is attached. Include only atoms you want to insert.");
            ui.add(egui::DragValue::new(&mut form.anchor).range(1..=mol.atoms.len().max(1)));
            if let Some(symbol) = form.anchor.checked_sub(1).and_then(|i| mol.atoms.get(i)) { ui.weak(symbol); }
            if ui.add_enabled(selected.is_some(), egui::Button::new("Use picked")).clicked_once() { if let Some(i) = selected { form.anchor = i + 1; } }
        });
        let fragment = library.check_name(&form.name).and_then(|_| indices)
            .and_then(|indices| Fragment::from_selection(&form.name, mol, &indices, form.anchor.saturating_sub(1)));
        let reason = fragment.as_ref().err().cloned().unwrap_or_default();
        if ui.add_enabled(fragment.is_ok(), crate::ui_style::primary("Save to My fragments"))
            .on_disabled_hover_text(&reason).clicked_once() {
            match fragment.and_then(|f| library.save(f)) {
                Ok(key) => { *chosen = key; form.name.clear(); },
                Err(error) => library.message = Some(error),
            }
        }
        if !form.name.trim().is_empty() && !reason.is_empty() { ui.small(reason); }
    });
    if let Some(message) = &library.message {
        ui.small(message);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn molecule() -> Molecule {
        Molecule {
            atoms: vec!["C".into(), "O".into(), "H".into()],
            pos: vec![
                Vec3::new(3.0, 0.0, 0.0),
                Vec3::new(4.4, 0.0, 0.0),
                Vec3::new(5.3, 0.0, 0.0),
            ],
            bonds: vec![],
            hydrogen_bonds: vec![],
        }
    }
    #[test]
    fn named_selection_puts_the_attachment_first_and_preserves_geometry() {
        let mol = molecule();
        let f = Fragment::from_selection("Hydroxyl", &mol, &[1, 2], 1).unwrap();
        let (n, atoms, xyz) = parse_xyz_angstrom(&f.xyz);
        assert_eq!(n, 2);
        assert_eq!(atoms, ["O", "H"]);
        assert_eq!(xyz[0], [0.0, 0.0, 0.0]);
        assert!((xyz[1][0] - 0.9).abs() < 1e-6);
        assert!(Fragment::from_selection(" ", &mol, &[0, 1], 0).is_err());
        assert!(Fragment::from_selection("Named", &mol, &[1, 2], 0).is_err());
        assert_eq!(selection("1-2,2, 3", 3).unwrap(), vec![0, 1, 2]);
        assert!(selection("0, 3", 3).is_err());
        assert!(selection("1-999999999", 3).is_err());
    }
    #[test]
    fn named_buttons_survive_reload_and_removal_can_be_undone() {
        let path =
            std::env::temp_dir().join(format!("beavyr-fragments-test-{}", std::process::id()));
        let mut library = Library::at(Some(path.clone()));
        let fragment = Fragment::from_selection("My alcohol", &molecule(), &[0, 1, 2], 1).unwrap();
        let key = library.save(fragment).unwrap();
        let mut library = Library::at(Some(path.clone()));
        assert_eq!(library.find(&key).unwrap().name, "My alcohol");
        assert!(library.check_name("my alcohol").is_err());
        library.remove(0).unwrap();
        assert!(library.find(&key).is_none());
        assert!(Library::at(Some(path.clone())).entries.is_empty());
        library.undo_remove().unwrap();
        assert!(library.find(&key).is_some());
        fs::remove_dir_all(path).unwrap();
    }

    #[test]
    fn saved_fragment_attachment_uses_the_existing_fragment_editor_targets() {
        use super::super::{
            builder_ui::{commit_fragment_connect, ZMatrixBuilderState},
            zmat2xyz,
        };
        let fragment = Fragment::from_selection("Hydroxyl", &molecule(), &[1, 2], 1).unwrap();
        let mut mol = Molecule {
            atoms: vec!["C".into()],
            pos: vec![Vec3::ZERO],
            bonds: vec![],
            hydrogen_bonds: vec![],
        };
        let mut builder = ZMatrixBuilderState::default();
        builder.custom_fragments = Library::at(None);
        builder.zmat = zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos);
        builder.frag_name = fragment.key();
        builder.custom_fragments.entries.push(fragment);
        builder.selected_index = Some(0);
        commit_fragment_connect(
            &mut builder,
            &mut mol,
            &mut crate::settings::MolSettings::default(),
        );
        assert!(builder.last_error.is_none(), "{:?}", builder.last_error);
        assert_eq!(mol.atoms, ["C", "O", "H"]);
        assert!((mol.pos[1].distance(mol.pos[2]) - 0.9).abs() < 1e-5);
        assert_eq!(builder.placed_editor.moving, Some(vec![1, 2]));
        assert_eq!(&builder.placed_editor.picks[..2], &[0, 1]);
        assert_eq!(builder.placed_editor.session_origin, Some(mol.pos.clone()));
    }
}
