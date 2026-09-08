//! Opening files named on the command line.
//!
//! `beavyr ZZAllyl.xyz` should put the structure on screen before the window
//! is even drawn, and `beavyr water.hess` should have the frequencies analysed
//! by the time you look at it. This module decides what each file is and hands
//! it to the tool that owns it.
//!
//! The deciding is separated from the loading on purpose: [`classify`] is a
//! pure function of a path and the first few kilobytes of the file, so every
//! rule about what counts as what is testable without a running app.
//!
//! See `docs/superpowers/specs/2026-09-08-open-a-file-from-the-command-line-design.md`.

use bevy::prelude::*;
use std::path::{Path, PathBuf};

use crate::events::MoleculeChanged;
use crate::molecule::Molecule;
use crate::settings::MolSettings;
use crate::ui_layout::{Tab, UiLayout};

/// Upper bound on how much of an output file is read to decide what it is.
///
/// Only `.out` and `.log` are ever read at all, and they are read whole up to
/// this size. Reading a prefix does not work: a Gaussian job restarted with
/// `Link1` puts its banner *after* the continuation header, and in
/// `examples/c164-ts1-2-1001.log` the word "Gaussian" first appears 172 kB in.
/// Any fixed peek is therefore a guess about someone else's file layout, and
/// this one guessed wrong on a real file in this repository.
///
/// Reading twelve megabytes once at startup costs a few milliseconds and is
/// invisible next to compiling the shaders.
const MAX_SNIFF_BYTES: u64 = 64 * 1024 * 1024;

/// What a file on the command line turned out to be.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileKind {
    /// XYZ coordinates, one frame or many.
    Xyz,
    /// An ORCA `.hess`, or a Gaussian log with frequencies in it.
    Hessian,
    /// Orbitals and a basis set.
    Molden,
    /// An ORCA Cartesian gradient.
    Engrad,
    /// An output file with excited states in it.
    TdDft,
    /// An output file we recognise as an output file, but which holds nothing
    /// Beavyr can use. Distinct from `Unrecognised` so the user can be told
    /// what was looked for.
    OutputWithNothingUsable,
    /// Not a file Beavyr knows.
    Unrecognised,
}

impl FileKind {
    /// The window to open once a file of this kind has loaded, if any.
    pub fn window(self) -> Option<Tab> {
        match self {
            FileKind::Xyz => Some(Tab::Structure),
            FileKind::Hessian | FileKind::Engrad => Some(Tab::Vibrations),
            FileKind::Molden => Some(Tab::Surface),
            FileKind::TdDft => Some(Tab::UvVis),
            FileKind::OutputWithNothingUsable | FileKind::Unrecognised => None,
        }
    }
}

/// Decide what a file is, from its name and its opening bytes.
///
/// Extension first, because that is what the user thinks in. `.out` and `.log`
/// are written by every program in the field and settle nothing, so those are
/// decided by looking inside.
///
/// `.inp` is deliberately absent: it is an ORCA *input* file, and the Molden
/// dialog's `inp` filter is for the `molden.input` spelling. Guessing between
/// the two would be wrong about half the time.
pub fn classify(path: &Path, head: &str) -> FileKind {
    // Files off a Windows share or an old cluster often arrive shouting.
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();

    if name.ends_with(".molden.input") {
        return FileKind::Molden;
    }

    let ext = path
        .extension()
        .map(|e| e.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();

    match ext.as_str() {
        "xyz" => FileKind::Xyz,
        "hess" => FileKind::Hessian,
        "molden" => FileKind::Molden,
        "engrad" => FileKind::Engrad,
        "out" | "log" => classify_output(head),
        _ => FileKind::Unrecognised,
    }
}

/// The ambiguous half: an `.out` or a `.log` could come from anywhere.
///
/// A Gaussian frequency job and an ORCA excited-state job cannot both match:
/// an ORCA `Freq` output carries no force constants at all -- those live in the
/// companion `.hess` -- so there is no file that is both.
fn classify_output(head: &str) -> FileKind {
    if crate::qchem_interfaces::gaussian_log::is_gaussian_log(head) {
        return FileKind::Hessian;
    }
    if has_excited_states(head) {
        return FileKind::TdDft;
    }
    FileKind::OutputWithNothingUsable
}

/// Whether an output file looks like it contains excited states.
///
/// Only a hint: the real answer is whether the parser succeeds, and the caller
/// reports that. This exists to route the file to the right parser, not to
/// promise the parse will work.
fn has_excited_states(head: &str) -> bool {
    const MARKERS: [&str; 5] = [
        "ABSORPTION SPECTRUM",
        "TD-DFT",
        "TDDFT",
        "EXCITED STATES",
        "CIS-EXCITED STATES",
    ];
    let upper = head.to_ascii_uppercase();
    MARKERS.iter().any(|m| upper.contains(m))
}

/// Read as much of a file as [`classify`] needs to decide what it is.
///
/// Only `.out` and `.log` need reading at all -- every other extension settles
/// the question by itself -- so nothing else is opened. An unreadable file
/// gives an empty string rather than an error: a binary `.gbw` is not a parse
/// failure worth reporting, it is simply not one of ours.
fn sniff(path: &Path) -> String {
    let ext = path
        .extension()
        .map(|e| e.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    if ext != "out" && ext != "log" {
        return String::new();
    }

    use std::io::Read;
    let Ok(file) = std::fs::File::open(path) else {
        return String::new();
    };
    let mut buf = Vec::new();
    if file.take(MAX_SNIFF_BYTES).read_to_end(&mut buf).is_err() {
        return String::new();
    }
    String::from_utf8_lossy(&buf).into_owned()
}

/// What the command line asked for.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct CommandLine {
    pub paths: Vec<PathBuf>,
    pub help: bool,
    /// Arguments that begin with `-` and are not `--help`. Reported rather
    /// than ignored, so a mistyped flag does not look like it worked.
    pub unknown_flags: Vec<String>,
}

/// Split the arguments into paths, `--help`, and anything unrecognised.
///
/// No argument-parsing crate: this is one flag and a list of filenames, and
/// the dependency list is short on purpose.
pub fn parse_args<I: IntoIterator<Item = String>>(args: I) -> CommandLine {
    let mut out = CommandLine::default();
    for arg in args {
        match arg.as_str() {
            "-h" | "--help" => out.help = true,
            _ if arg.starts_with('-') && arg.len() > 1 => out.unknown_flags.push(arg),
            _ => out.paths.push(PathBuf::from(arg)),
        }
    }
    out
}

/// What `beavyr --help` prints.
pub fn help_text() -> String {
    "\
Beavyr -- molecular viewer and quantum-chemistry front-end

    beavyr                    open with an empty viewport
    beavyr FILE...            open with these files loaded
    beavyr --help             this message

Files are recognised by their extension:

    .xyz            a structure; several frames become a trajectory
    .hess           an ORCA Hessian -- frequencies are analysed on opening
    .molden         orbitals and basis set
    .engrad         an ORCA gradient, for the reaction-path projection
    .out, .log      looked inside: a Gaussian frequency job, or an
                    excited-state calculation

Several files may be given at once and each is loaded into its own tool. Where
more than one carries a geometry, the last one named is the structure left on
screen."
        .to_string()
}

/// What happened while loading the command line, for the user to read.
///
/// Startup problems have nowhere obvious to appear -- there is no panel that
/// owns "the file you named does not exist" -- so they are collected here and
/// shown at the top of the Structure panel, which is always present.
#[derive(Resource, Default)]
pub struct StartupLoadReport {
    /// `(message, is_error)`, in the order the files were given.
    pub messages: Vec<(String, bool)>,
}

impl StartupLoadReport {
    fn ok(&mut self, message: String) {
        println!("{message}");
        self.messages.push((message, false));
    }

    fn failed(&mut self, message: String) {
        eprintln!("beavyr: {message}");
        self.messages.push((message, true));
    }
}

/// Resources the loaders need, grouped because Bevy allows a system only
/// sixteen parameters and each loader wants two or three of its own.
type LoadTargets<'w> = (
    ResMut<'w, crate::ui::XyzBuffer>,
    ResMut<'w, crate::trajectory::TrajectoryState>,
    ResMut<'w, crate::orbitals::OrbitalState>,
    ResMut<'w, crate::qchem_interfaces::xtb_freq::XtbFreqPanelState>,
    ResMut<'w, crate::qchem_interfaces::xtb_freq::XtbFrequencyTask>,
    ResMut<'w, crate::uvvis::UvVisState>,
);

/// Load whatever was named on the command line, once, at startup.
///
/// Ordered before the camera is centred, so the view frames the molecule that
/// was loaded rather than an empty scene.
pub fn load_command_line_files(
    mut report: ResMut<StartupLoadReport>,
    mut layout: ResMut<UiLayout>,
    mut mol: ResMut<Molecule>,
    mut settings: ResMut<MolSettings>,
    mut cam: ResMut<crate::camera::OrbitCamera>,
    mut ev_changed: MessageWriter<MoleculeChanged>,
    mut targets: LoadTargets,
) {
    // `--help` never reaches here: `main` answers it before Bevy is built, so
    // the text prints without a window appearing first.
    let cmd = parse_args(std::env::args().skip(1));

    for flag in &cmd.unknown_flags {
        report.failed(format!("unknown option {flag} -- try --help"));
    }

    if cmd.paths.is_empty() {
        return;
    }

    let (xyz_buf, traj, orbitals, freq_panel, freq_task, uvvis) = &mut targets;

    // An .engrad only means something against a Hessian, and the Hessian may
    // be named after it on the command line. Hold gradients back and apply
    // them once everything else has loaded.
    let mut pending_gradients: Vec<PathBuf> = Vec::new();

    for path in &cmd.paths {
        if !path.exists() {
            report.failed(format!("{}: no such file", path.display()));
            continue;
        }

        let kind = classify(path, &sniff(path));
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.display().to_string());

        let loaded = match kind {
            FileKind::Xyz => match crate::ui::load_xyz_from_path(
                path,
                xyz_buf,
                traj,
                &mut mol,
                &mut settings,
                &mut cam,
                &mut ev_changed,
            ) {
                Ok(1) => {
                    report.ok(format!("{name}: {} atoms", mol.atoms.len()));
                    true
                }
                Ok(frames) => {
                    report.ok(format!(
                        "{name}: {frames} frames, {} atoms",
                        mol.atoms.len()
                    ));
                    layout.open[Tab::Trajectory.index()] = true;
                    true
                }
                Err(err) => {
                    report.failed(err);
                    false
                }
            },

            FileKind::Hessian => {
                let ok = crate::qchem_interfaces::xtb_freq::load_hessian_from_path(
                    path, freq_panel, freq_task,
                );
                match freq_task.last_message.clone() {
                    Some(message) if ok => report.ok(message),
                    Some(message) => report.failed(message),
                    None if ok => report.ok(format!("{name}: loaded")),
                    None => report.failed(format!("{name}: could not be read")),
                }
                if ok {
                    // The analysed geometry becomes the structure on screen,
                    // the same as when a Hessian is loaded from the panel.
                    // Taken, not cloned: the Vibrations panel applies this on
                    // the first frame it is drawn, and doing it twice would
                    // re-emit the structure-changed event for no reason. Done
                    // here as well as there because the camera is centred at
                    // startup, before any panel has drawn.
                    if let Some((atoms, pos)) = freq_task.pending_geometry.take() {
                        mol.atoms = atoms;
                        mol.pos = pos;
                        mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
                        ev_changed.write(MoleculeChanged::parse_xyz(true));
                        if let Some(c) = crate::ui::compute_centroid(&mol.pos) {
                            cam.target = c;
                        }
                    }
                }
                ok
            }

            FileKind::Molden => {
                let ok = crate::orbitals::ui::load_molden_from_path(
                    path,
                    orbitals,
                    &mut mol,
                    &mut settings,
                    &mut ev_changed,
                );
                match (&orbitals.warning, ok) {
                    (_, true) => report.ok(format!("{name}: {} atoms", mol.atoms.len())),
                    (Some(warning), false) => report.failed(format!("{name}: {warning}")),
                    (None, false) => report.failed(format!("{name}: could not be read")),
                }
                ok
            }

            FileKind::TdDft => {
                let ok = crate::uvvis::ui::load_tddft_from_path(path, uvvis);
                match (&uvvis.load_error, ok) {
                    (_, true) => report.ok(format!("{name}: excited states loaded")),
                    (Some(err), false) => report.failed(err.clone()),
                    (None, false) => report.failed(format!("{name}: could not be read")),
                }
                ok
            }

            FileKind::Engrad => {
                pending_gradients.push(path.clone());
                // Reported once it has been applied, below.
                false
            }

            FileKind::OutputWithNothingUsable => {
                report.failed(format!(
                    "{name}: an output file with no frequencies and no excited states in it. \
                     An ORCA Freq job keeps its force constants in the companion .hess file, \
                     not in the .out -- load that instead."
                ));
                false
            }

            FileKind::Unrecognised => {
                report.failed(format!(
                    "{name}: not a file Beavyr reads. Try --help for the list."
                ));
                false
            }
        };

        if loaded {
            if let Some(tab) = kind.window() {
                layout.open[tab.index()] = true;
            }
        }
    }

    for path in pending_gradients {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.display().to_string());
        let ok = crate::qchem_interfaces::xtb_freq::load_gradient_from_path(
            &path, freq_panel, freq_task,
        );
        match freq_task.last_message.clone() {
            Some(message) if ok => report.ok(message),
            Some(message) => report.failed(format!("{name}: {message}")),
            None if ok => report.ok(format!("{name}: gradient loaded")),
            None => report.failed(format!("{name}: gradient could not be read")),
        }
        if ok {
            layout.open[Tab::Vibrations.index()] = true;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kind_of(name: &str, head: &str) -> FileKind {
        classify(Path::new(name), head)
    }

    #[test]
    fn each_plain_extension_maps_to_its_tool() {
        assert_eq!(kind_of("ZZAllyl.xyz", ""), FileKind::Xyz);
        assert_eq!(kind_of("water.hess", ""), FileKind::Hessian);
        assert_eq!(kind_of("mo.molden", ""), FileKind::Molden);
        assert_eq!(kind_of("job.engrad", ""), FileKind::Engrad);
    }

    /// Files off a Windows share or an older cluster arrive capitalised.
    #[test]
    fn extensions_are_recognised_whatever_their_case() {
        assert_eq!(kind_of("ZZAllyl.XYZ", ""), FileKind::Xyz);
        assert_eq!(kind_of("WATER.Hess", ""), FileKind::Hessian);
        assert_eq!(kind_of("MO.MOLDEN", ""), FileKind::Molden);
    }

    /// Molden's own spelling of its format, which has two dots in it and so
    /// has an "extension" of `.input`.
    #[test]
    fn the_molden_input_spelling_is_recognised() {
        assert_eq!(kind_of("orbitals.molden.input", ""), FileKind::Molden);
    }

    /// An ORCA input file, which is a different thing entirely and must not be
    /// mistaken for the Molden spelling above.
    #[test]
    fn an_orca_input_file_is_not_claimed() {
        assert_eq!(kind_of("job.inp", ""), FileKind::Unrecognised);
    }

    #[test]
    fn a_gaussian_log_is_read_for_its_frequencies() {
        let head = " Entering Gaussian System, Link 0=g16\n Gaussian 16:  ES64L-G16RevC.01\n";
        assert_eq!(kind_of("freq.log", head), FileKind::Hessian);
        assert_eq!(kind_of("freq.out", head), FileKind::Hessian);
    }

    #[test]
    fn an_output_with_an_absorption_spectrum_goes_to_uv_vis() {
        let head = "* O   R   C   A *\n-----------------\nABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\n";
        assert_eq!(kind_of("exc.out", head), FileKind::TdDft);
    }

    /// The case that matters most: an ORCA Freq output has no force constants
    /// in it at all. Saying so beats opening nothing.
    #[test]
    fn an_output_with_neither_is_reported_rather_than_guessed() {
        let head = "* O   R   C   A *\n\nSCF converged.\nTOTAL SCF ENERGY  -76.4\n";
        assert_eq!(kind_of("scf.out", head), FileKind::OutputWithNothingUsable);
    }

    #[test]
    fn an_unknown_extension_is_not_guessed_at() {
        assert_eq!(kind_of("wavefunction.gbw", ""), FileKind::Unrecognised);
        assert_eq!(kind_of("notes.pdf", ""), FileKind::Unrecognised);
        assert_eq!(kind_of("noextension", ""), FileKind::Unrecognised);
    }

    #[test]
    fn every_loadable_kind_names_a_window_to_open() {
        for kind in [
            FileKind::Xyz,
            FileKind::Hessian,
            FileKind::Molden,
            FileKind::Engrad,
            FileKind::TdDft,
        ] {
            assert!(kind.window().is_some(), "{kind:?} opens no window");
        }
        assert!(FileKind::Unrecognised.window().is_none());
        assert!(FileKind::OutputWithNothingUsable.window().is_none());
    }

    #[test]
    fn arguments_split_into_paths_and_flags() {
        let cmd = parse_args(["a.xyz".to_string(), "b.hess".to_string()]);
        assert_eq!(cmd.paths.len(), 2);
        assert!(!cmd.help);
        assert!(cmd.unknown_flags.is_empty());
    }

    #[test]
    fn help_is_recognised_both_ways() {
        assert!(parse_args(["--help".to_string()]).help);
        assert!(parse_args(["-h".to_string()]).help);
    }

    /// A mistyped flag must not be silently treated as a filename, or the user
    /// gets "no such file: --hepl" and has to work out why.
    #[test]
    fn an_unknown_flag_is_kept_apart_from_the_filenames() {
        let cmd = parse_args(["--hepl".to_string(), "a.xyz".to_string()]);
        assert_eq!(cmd.unknown_flags, vec!["--hepl".to_string()]);
        assert_eq!(cmd.paths, vec![PathBuf::from("a.xyz")]);
    }

    /// A lone "-" is a filename by convention, not a flag.
    #[test]
    fn a_bare_dash_is_treated_as_a_path() {
        let cmd = parse_args(["-".to_string()]);
        assert!(cmd.unknown_flags.is_empty());
        assert_eq!(cmd.paths, vec![PathBuf::from("-")]);
    }

    #[test]
    fn the_help_text_names_every_extension_that_is_accepted() {
        let help = help_text();
        for ext in [".xyz", ".hess", ".molden", ".engrad", ".out", ".log"] {
            assert!(help.contains(ext), "help does not mention {ext}");
        }
    }

    /// The files the user actually has. These are read from disk, so this also
    /// checks the sniffing survives a real ORCA banner.
    #[test]
    fn the_projects_own_example_files_are_recognised() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let cases: [(&str, FileKind); 4] = [
            ("examples/ZZAllyl.xyz", FileKind::Xyz),
            (
                "examples/TS_Gamma-3-6_ZZ-S_11_Compound_1.hess",
                FileKind::Hessian,
            ),
            (
                "examples/Exc_TPSSh-D4-def2-TZVP_SCMwater_16.molden",
                FileKind::Molden,
            ),
            (
                "examples/Exc_root_OnlyFreq_from_optimized.engrad",
                FileKind::Engrad,
            ),
        ];
        for (relative, expected) in cases {
            let path = root.join(relative);
            assert!(path.exists(), "missing example file {relative}");
            assert_eq!(classify(&path, &sniff(&path)), expected, "{relative}");
        }
    }

    #[test]
    fn a_real_orca_excited_state_output_is_recognised() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let path = root.join("examples/AbsorptionSpectrum/Exc_TPSSh-D4-def2-TZVP_SCMwater_16.out");
        assert!(path.exists());
        assert_eq!(classify(&path, &sniff(&path)), FileKind::TdDft);
    }


    /// The file that caught the original mistake. It is a Gaussian job
    /// restarted with `Link1`, so it opens with the continuation header and
    /// does not say "Gaussian" until 172 kB in -- an 8 kB peek called it an
    /// unusable output file.
    #[test]
    fn a_gaussian_log_is_found_even_when_its_banner_is_far_into_the_file() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let path = root.join("examples/c164-ts1-2-1001.log");
        let text = std::fs::read_to_string(&path).unwrap();
        let banner_at = text.find("Gaussian").expect("no banner in this file at all");
        assert!(
            banner_at > 8192,
            "this test is pointless if the banner is near the top (found at {banner_at})"
        );
        assert_eq!(classify(&path, &sniff(&path)), FileKind::Hessian);
    }

    #[test]
    fn a_real_gaussian_log_is_recognised() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let path = root.join("examples/c164-ts1-2-1001.log");
        assert!(path.exists());
        assert_eq!(classify(&path, &sniff(&path)), FileKind::Hessian);
    }
}
