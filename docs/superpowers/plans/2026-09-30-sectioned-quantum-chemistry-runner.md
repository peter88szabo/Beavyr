# Sectioned Quantum Chemistry Runner Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the quantum-chemistry runner's eleven-row form with six collapsed sections that each show their current setting, behind a switch that keeps the old layout reachable.

**Architecture:** A pure section model decides which sections exist, what each header reads and whether each applies, with no dependency on the drawing library so every rule is a tested function. The drawing code reads that model. The existing flat layout moves to a file of its own, unchanged, and a remembered setting chooses between them.

**Tech Stack:** Rust 2021, Bevy 0.19, bevy_egui 0.42. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-09-30-panel-usability-redesign.md`

## Global Constraints

- **Do not commit.** Peter commits his own code. Every task ends by running the tests and leaving the working tree for him. No `git add`, no `git commit`, no branches.
- **Cap all builds and test runs at four CPUs:** prefix with `taskset -c 0-3` and pass `--jobs 4`, per `src/AGENTS.md`.
- **Run tests serially:** `taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1`.
- **No new crate dependencies.** The project keeps five and the README advertises that.
- **No capability is removed.** Everything the panel can do today it must still do.
- **No load-bearing information lives only in a hover tooltip.** Where a user needs a fact in order to act, it goes on screen. Hover may still carry elaboration.
- **Nothing appears or disappears as settings change.** A section that stops applying says so on its header, in place.
- **The new layout is the default**, and the choice is remembered between sessions.

## Review Focus

Five things the spec implies, that no single task's happy path exercises, and that would bite a real user. Each has a test placed in the task that owns the code.

1. **DREIDING has no level of theory.** Its Level of theory header must read that it has none rather than the section disappearing. Test in Task 1.
2. **Three programs have no input file.** xTB, Behemoth and Sparrow are driven by settings; their Input file header must say so rather than vanishing. Test in Task 1.
3. **A job with no excited states must not resize the Job section.** Only TD-DFT adds those controls, and the header must read the same way either way. Test in Task 1.
4. **Switching program must update every header at once.** A stale header showing the previous program's basis is worse than no header. Test in Task 3.
5. **A long path or functional name must not push the expander off the row.** Headers are built from user-supplied strings of unbounded length. Test in Task 1.

---

### Task 1: The section model

A pure description of the panel: which sections there are, what each header reads for the current state, and whether each applies. No egui. This is where every rule lives and where every rule is tested.

**Files:**
- Create: `src/qchem_panel/sections.rs`
- Modify: `src/qchem_panel/mod.rs` (add `pub mod sections;`)

**Interfaces:**
- Consumes: `QcPanelState` from `src/qchem_panel/mod.rs`; `QcProgram`; `crate::qchem_interfaces::job`.
- Produces:
  - `pub enum Section { Program, LevelOfTheory, Job, Molecule, Resources, InputFile }`
  - `pub const ALL: [Section; 6]`
  - `pub fn title(Section) -> &'static str`
  - `pub fn header(Section, &QcPanelState) -> String` — the current setting, shown on the collapsed header
  - `pub fn applies(Section, &QcPanelState) -> bool` — false only where the section has nothing to set
  - `pub fn unavailable_note(Section, &QcPanelState) -> Option<&'static str>` — why, when `applies` is false

- [ ] **Step 1: Write the failing tests**

Create `src/qchem_panel/sections.rs` containing only the tests for now:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::qchem_interfaces::job::{JobType, QcMethod};
    use crate::qchem_interfaces::program::QcProgram;
    use crate::qchem_panel::QcPanelState;

    fn orca() -> QcPanelState {
        let mut state = QcPanelState::default();
        state.set_program(QcProgram::Orca);
        state.method = Some(QcMethod::Dft);
        state.config.functional = "B3LYP".to_string();
        state.config.basis = "def2-SVP".to_string();
        state
    }

    /// Every section names itself, and the list is the six the spec gives.
    #[test]
    fn there_are_six_sections_and_each_is_named() {
        assert_eq!(ALL.len(), 6);
        for section in ALL {
            assert!(!title(section).is_empty(), "{section:?}");
        }
    }

    /// The header carries the current setting, which is what makes the whole
    /// panel readable without opening anything.
    #[test]
    fn the_header_carries_the_current_setting() {
        let state = orca();
        assert!(header(Section::Program, &state).contains("ORCA"));
        let theory = header(Section::LevelOfTheory, &state);
        assert!(theory.contains("B3LYP"), "{theory}");
        assert!(theory.contains("def2-SVP"), "{theory}");
        assert!(header(Section::Job, &state).contains("Optimization"));
    }

    /// Review Focus 1. A force field has no level of theory. The section must
    /// say so rather than disappear, or the form moves under the cursor --
    /// which is the complaint this redesign exists to fix.
    #[test]
    fn a_force_field_says_it_has_no_level_of_theory() {
        let mut state = QcPanelState::default();
        state.set_program(QcProgram::Dreiding);
        assert!(!applies(Section::LevelOfTheory, &state));
        let note = unavailable_note(Section::LevelOfTheory, &state)
            .expect("a reason is given");
        assert!(note.contains("force field"), "{note}");
        // And it still has a header, because it is still on screen.
        assert!(!header(Section::LevelOfTheory, &state).is_empty());
    }

    /// Review Focus 2. Three programs are driven by settings, not by a file.
    #[test]
    fn programs_without_an_input_file_say_so() {
        let mut state = QcPanelState::default();
        for program in [QcProgram::Xtb, QcProgram::Behemoth, QcProgram::SparrowPy] {
            state.set_program(program);
            assert!(!applies(Section::InputFile, &state), "{program:?}");
            assert!(
                unavailable_note(Section::InputFile, &state).is_some(),
                "{program:?} must say why"
            );
        }
        state.set_program(QcProgram::Orca);
        assert!(applies(Section::InputFile, &state));
    }

    /// Review Focus 3. Only TD-DFT adds excited-state controls, and the Job
    /// header must read the same shape either way so the section does not
    /// change size as the job changes.
    #[test]
    fn the_job_header_reads_the_same_shape_for_every_job() {
        let mut state = orca();
        for job in [JobType::Optimize, JobType::TdDft, JobType::SinglePoint] {
            state.job = job;
            let text = header(Section::Job, &state);
            assert_eq!(text, job.label(), "{job:?}");
        }
    }

    /// Review Focus 5. Headers are built from strings the user supplies and
    /// those have no length limit. A header that grows without bound pushes
    /// the expander off the row and the section cannot be opened at all.
    #[test]
    fn a_header_is_bounded_however_long_the_setting_is() {
        let mut state = orca();
        state.config.functional = "x".repeat(400);
        state.config.basis = "y".repeat(400);
        for section in ALL {
            let text = header(section, &state);
            assert!(
                text.chars().count() <= MAX_HEADER_CHARS,
                "{section:?} header is {} chars",
                text.chars().count()
            );
        }
    }

    /// Every section that applies has a header, and every one that does not
    /// has a reason. Neither may be silent.
    #[test]
    fn nothing_is_silent() {
        let mut state = QcPanelState::default();
        for program in QcProgram::ALL {
            state.set_program(program);
            for section in ALL {
                assert!(!header(section, &state).is_empty(), "{program:?} {section:?}");
                if !applies(section, &state) {
                    assert!(
                        unavailable_note(section, &state).is_some(),
                        "{program:?} {section:?} is unavailable without saying why"
                    );
                }
            }
        }
    }
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1 qchem_panel::sections`
Expected: compile errors — `Section`, `ALL`, `title`, `header`, `applies`, `unavailable_note`, `MAX_HEADER_CHARS` are not defined.

- [ ] **Step 3: Write the model**

Prepend to `src/qchem_panel/sections.rs`, above the tests:

```rust
//! What the runner panel is made of, and what each part currently says.
//!
//! No drawing here. Every rule the panel follows -- which sections exist, what
//! a header reads, whether a section applies to the chosen program -- is a
//! plain function with a test, so "a force field has no level of theory" is
//! something a test asserts rather than something only a click reveals.

use super::QcPanelState;
use crate::qchem_interfaces::program::QcProgram;

/// The longest a collapsed header may be.
///
/// Headers are built from strings the user supplies -- a hand-typed functional,
/// a deep path -- and those have no length limit. An unbounded header pushes
/// the expander off the row, and a section that cannot be opened is worse than
/// one with a truncated label.
pub const MAX_HEADER_CHARS: usize = 60;

/// One collapsible part of the panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Section {
    Program,
    LevelOfTheory,
    Job,
    Molecule,
    Resources,
    InputFile,
}

/// The sections, in the order they are drawn: what runs it, at what level of
/// theory, to do what, on what, with which resources, and the file that
/// results.
pub const ALL: [Section; 6] = [
    Section::Program,
    Section::LevelOfTheory,
    Section::Job,
    Section::Molecule,
    Section::Resources,
    Section::InputFile,
];

pub fn title(section: Section) -> &'static str {
    match section {
        Section::Program => "Program",
        Section::LevelOfTheory => "Level of theory",
        Section::Job => "Job",
        Section::Molecule => "Molecule",
        Section::Resources => "Resources",
        Section::InputFile => "Input file",
    }
}

/// Shortens `text` to fit a header, keeping the end.
///
/// The end, not the beginning: the distinctive part of a long path or a long
/// functional specification is at the end, and a row of leading directories
/// tells nobody anything.
fn fit(text: &str) -> String {
    let count = text.chars().count();
    if count <= MAX_HEADER_CHARS {
        return text.to_string();
    }
    let keep = MAX_HEADER_CHARS.saturating_sub(1);
    let skip = count - keep;
    format!("\u{2026}{}", text.chars().skip(skip).collect::<String>())
}

/// What this section's collapsed header reads for the current state.
pub fn header(section: Section, state: &QcPanelState) -> String {
    let text = match section {
        Section::Program => state.program.label().to_string(),
        Section::LevelOfTheory => {
            if !applies(section, state) {
                "none: a force field".to_string()
            } else if state.needs_basis() {
                format!("{}/{}", theory_name(state), state.config.basis)
            } else {
                theory_name(state)
            }
        }
        Section::Job => state.job.label().to_string(),
        Section::Molecule => format!(
            "charge {}, {}",
            state.charge,
            multiplicity_name(state.multiplicity)
        ),
        Section::Resources => format!(
            "{} thread{}, {} MB",
            state.config.nproc,
            if state.config.nproc == 1 { "" } else { "s" },
            state.config.memory_mb
        ),
        Section::InputFile => {
            if !applies(section, state) {
                "none: driven by settings".to_string()
            } else if state.advanced {
                "written by hand".to_string()
            } else {
                "built from the settings above".to_string()
            }
        }
    };
    fit(&text)
}

/// The level of theory as a name, without the basis set.
fn theory_name(state: &QcPanelState) -> String {
    use crate::qchem_interfaces::job::QcMethod;
    match state.method {
        None => "force field".to_string(),
        Some(QcMethod::Dft) => match state.program {
            QcProgram::Orca => {
                crate::qchem_interfaces::orca_method::functional_label(&state.config.functional)
                    .to_string()
            }
            QcProgram::PySCF => {
                crate::qchem_interfaces::pyscf_method::functional_label(&state.config.functional)
                    .to_string()
            }
            QcProgram::Psi4 | QcProgram::Psi4Py => {
                crate::qchem_interfaces::psi4_method::functional_label(&state.config.functional)
                    .to_string()
            }
            _ => state.config.functional.clone(),
        },
        Some(QcMethod::Semiempirical) => {
            crate::qchem_interfaces::sparrow_method::method_label(&state.config.functional)
                .to_string()
        }
        Some(other) => other.label().to_string(),
    }
}

/// A spin state named rather than numbered, because "singlet" is what a
/// chemist reads and "1" is what a program wants.
fn multiplicity_name(multiplicity: i32) -> String {
    match multiplicity {
        1 => "singlet".to_string(),
        2 => "doublet".to_string(),
        3 => "triplet".to_string(),
        4 => "quartet".to_string(),
        5 => "quintet".to_string(),
        other => format!("multiplicity {other}"),
    }
}

/// Whether this section has anything to set.
///
/// False is not "hide it". A section that does not apply stays on screen and
/// says why, because a form whose rows come and go moves under the cursor --
/// which is the complaint this redesign exists to answer.
pub fn applies(section: Section, state: &QcPanelState) -> bool {
    match section {
        Section::LevelOfTheory => state.has_level_of_theory(),
        Section::InputFile => state.has_input_file(),
        // The rest always have something to set.
        Section::Program | Section::Job | Section::Molecule | Section::Resources => true,
    }
}

/// Why a section has nothing to set, when it has nothing to set.
pub fn unavailable_note(section: Section, state: &QcPanelState) -> Option<&'static str> {
    if applies(section, state) {
        return None;
    }
    Some(match section {
        Section::LevelOfTheory => {
            "DREIDING is a generic force field. Its parameters follow from the atom types, \
             which are perceived from the structure, so there is no level of theory to set."
        }
        Section::InputFile => {
            "This program is driven by settings rather than by an input file, so there is \
             nothing to edit or save. xTB and Behemoth take everything on their command \
             lines, Sparrow is driven through SCINE, and DREIDING runs inside Beavyr."
        }
        _ => "Nothing to set here.",
    })
}
```

Then add to `src/qchem_panel/mod.rs`, beside the other module declarations:

```rust
pub mod sections;
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1 qchem_panel::sections`
Expected: PASS, 7 tests.

- [ ] **Step 5: Run the whole suite and leave it for review**

Run: `taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1`
Expected: every test passes; the count rises by 7.

Do not commit. Report the count and stop.

---

### Task 2: Keep the old layout reachable

Move the current flat layout into a file of its own, unchanged, and add the setting that chooses between layouts. After this task both branches still draw the old layout: the switch exists and works, the new layout does not exist yet. That is deliberate — it proves the escape hatch before anything depends on it.

**Files:**
- Create: `src/qchem_panel/ui_classic.rs`
- Modify: `src/qchem_panel/mod.rs` (add `PanelLayout`, the `layout` field, `pub mod ui_classic;`)
- Modify: `src/qchem_panel/ui.rs` (move `body` and its helpers out; call the chosen layout)
- Modify: `src/qchem_interfaces/config.rs` (remember the choice)

**Interfaces:**
- Consumes: everything `src/qchem_panel/ui.rs` already uses.
- Produces:
  - `pub enum PanelLayout { Sections, Classic }` in `src/qchem_panel/mod.rs`
  - `QcPanelState::layout: PanelLayout`
  - `pub fn load_layout() -> PanelLayout` and `pub fn save_layout(PanelLayout)` in `src/qchem_interfaces/config.rs`
  - `pub fn body(ui: &mut egui::Ui, panel: &mut QcPanelState, task: &mut QcRunTask, opt_panel: &mut XtbPanelState, mol: &Molecule)` in `src/qchem_panel/ui_classic.rs` — the current `body`, verbatim, with the same parameter list. Task 3's `ui_sections::body` takes the identical signature, so the caller's two match arms differ only in the module they name.

- [ ] **Step 1: Write the failing test**

Add to `src/qchem_panel/mod.rs`, in its `#[cfg(test)] mod tests` block (create the block if absent):

```rust
#[cfg(test)]
mod tests {
    use super::*;

    /// The new layout is what a user gets unless they have chosen otherwise.
    #[test]
    fn the_new_layout_is_the_default() {
        assert_eq!(QcPanelState::default().layout, PanelLayout::Sections);
    }

    /// The choice survives a restart. Being dropped back into a layout you
    /// rejected is worse than having no switch at all, and the existing
    /// classic/windowed switch already has that wart; this one must not.
    #[test]
    fn the_layout_choice_is_remembered() {
        use crate::qchem_interfaces::config;

        let dir = std::env::temp_dir().join(format!(
            "beavyr_layout_test_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        // SAFETY: the suite runs single-threaded via --test-threads=1 and no
        // other test in this crate touches these variables concurrently.
        unsafe {
            std::env::remove_var("XDG_CONFIG_HOME");
            std::env::set_var("HOME", &dir);
        }

        // Nothing saved yet: the default.
        assert_eq!(config::load_layout(), PanelLayout::Sections);

        config::save_layout(PanelLayout::Classic);
        assert_eq!(config::load_layout(), PanelLayout::Classic);

        config::save_layout(PanelLayout::Sections);
        assert_eq!(config::load_layout(), PanelLayout::Sections);

        let _ = std::fs::remove_dir_all(&dir);
    }
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1 qchem_panel::tests`
Expected: compile errors — `PanelLayout`, the `layout` field, `config::load_layout` and `config::save_layout` do not exist.

- [ ] **Step 3: Add the layout setting**

In `src/qchem_panel/mod.rs`, add beside the other declarations:

```rust
pub mod ui_classic;
```

and above `QcPanelState`:

```rust
/// Which arrangement of the runner panel is drawn.
///
/// The sectioned layout replaced a flat eleven-row form. The old one stays
/// reachable rather than being deleted, because a redesign is a claim about
/// what someone finds usable and that claim should be testable by them, not
/// argued with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PanelLayout {
    /// Six collapsed sections, each showing its current setting.
    Sections,
    /// The original flat form, every control visible at once.
    Classic,
}

impl PanelLayout {
    pub fn label(self) -> &'static str {
        match self {
            PanelLayout::Sections => "sections",
            PanelLayout::Classic => "the old layout",
        }
    }

    /// The other one, for a control that toggles.
    pub fn other(self) -> PanelLayout {
        match self {
            PanelLayout::Sections => PanelLayout::Classic,
            PanelLayout::Classic => PanelLayout::Sections,
        }
    }
}
```

Add the field to `QcPanelState`:

```rust
    /// Which arrangement is drawn. Remembered between sessions.
    pub layout: PanelLayout,
```

and in `Default for QcPanelState`, load it rather than assuming:

```rust
            layout: crate::qchem_interfaces::config::load_layout(),
```

In `src/qchem_interfaces/config.rs`, append:

```rust
/// Where the runner panel's layout choice is remembered.
const LAYOUT_FILE: &str = "qc_panel_layout.txt";

/// The remembered layout, or the sectioned one when nothing is saved.
///
/// Never errors: an unreadable or absent config is "nothing chosen yet", which
/// is the default, not a failure worth interrupting anyone for.
pub fn load_layout() -> crate::qchem_panel::PanelLayout {
    use crate::qchem_panel::PanelLayout;
    let Some(base) = config_path_for(QcProgram::Xtb).and_then(|p| p.parent().map(Path::to_path_buf))
    else {
        return PanelLayout::Sections;
    };
    match std::fs::read_to_string(base.join(LAYOUT_FILE))
        .unwrap_or_default()
        .trim()
    {
        "classic" => PanelLayout::Classic,
        _ => PanelLayout::Sections,
    }
}

/// Best-effort save, like the executable paths: forgetting a layout across
/// restarts is an annoyance, never a reason to interrupt a session.
pub fn save_layout(layout: crate::qchem_panel::PanelLayout) {
    use crate::qchem_panel::PanelLayout;
    let Some(base) = config_path_for(QcProgram::Xtb).and_then(|p| p.parent().map(Path::to_path_buf))
    else {
        return;
    };
    let _ = std::fs::create_dir_all(&base);
    let text = match layout {
        PanelLayout::Classic => "classic",
        PanelLayout::Sections => "sections",
    };
    let _ = std::fs::write(base.join(LAYOUT_FILE), text);
}
```

Add `use std::path::Path;` to the top of `config.rs` if it is not already there.

- [ ] **Step 4: Run the test to verify it passes**

Run: `taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1 qchem_panel::tests`
Expected: PASS, 2 tests.

- [ ] **Step 5: Move the old layout, unchanged**

Cut the whole `fn body(...)` from `src/qchem_panel/ui.rs`, together with every helper only it calls — `program_row`, `charge_and_multiplicity`, `level_of_theory`, `f12_basis_row`, `adopt_method_defaults`, `job_row`, `excited_state_row`, `wavefunction_row`, `resources_row`, `input_row`, `run_row`, `results_row` — and paste them into a new `src/qchem_panel/ui_classic.rs`.

Head that file with:

```rust
//! The runner panel's original layout: every control in one flat column.
//!
//! Kept, unchanged, behind the layout switch. It is here rather than deleted
//! because the sectioned layout that replaced it is a claim about what is
//! usable, and a claim like that should be testable by the person it is made
//! for.
//!
//! Nothing new goes in this file. A capability added to the panel is added to
//! the sectioned layout; if the two ever have to differ, that is the signal
//! that this one has outlived its purpose.
```

Make `body` public: `pub fn body(...)`. Add the imports it needs at the top of
the new file, copied from `ui.rs`. `adopt_method_defaults` and `f12_basis_row`
stay private to this file: they are called only from `level_of_theory`, which
Task 3 reuses whole.

- [ ] **Step 6: Call the chosen layout**

In `src/qchem_panel/ui.rs`, replace the body call inside `qc_panel` with:

```rust
        .show(ctx, |ui| {
            match panel.layout {
                // Not yet written: Task 3 replaces this arm. Until then both
                // choices draw the old layout, which proves the switch works
                // before anything depends on it.
                PanelLayout::Sections => super::ui_classic::body(ui, panel, task, opt_panel, mol),
                PanelLayout::Classic => super::ui_classic::body(ui, panel, task, opt_panel, mol),
            }
            layout_switch(ui, panel);
        });
```

and add, in `ui.rs`:

```rust
/// The control that swaps layouts, at the foot of the panel.
///
/// Worded as what it does rather than as a setting, and placed last so it is
/// out of the way of ordinary use without being hidden.
fn layout_switch(ui: &mut egui::Ui, panel: &mut QcPanelState) {
    ui.separator();
    ui.horizontal(|ui| {
        if ui
            .small_button(format!("use {}", panel.layout.other().label()))
            .clicked()
        {
            panel.layout = panel.layout.other();
            crate::qchem_interfaces::config::save_layout(panel.layout);
        }
    });
}
```

Add `use super::PanelLayout;` to `ui.rs`.

- [ ] **Step 7: Run the whole suite and check the build**

Run: `taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1`
Expected: every test passes.

Run: `taskset -c 0-3 cargo build --jobs 4 2>&1 | grep -c "qchem_panel/"`
Expected: `0` — no warnings in the panel.

Do not commit. Report and stop.

---

### Task 3: Draw the sectioned layout

**Files:**
- Create: `src/qchem_panel/ui_sections.rs`
- Modify: `src/qchem_panel/ui.rs` (point the `Sections` arm at it)
- Modify: `src/qchem_panel/mod.rs` (add `pub mod ui_sections;`)

**Interfaces:**
- Consumes: `sections::{Section, ALL, title, header, applies, unavailable_note}` from Task 1; `ui_classic::{adopt_method_defaults, f12_basis_row}` from Task 2; the existing `method_ui` row helpers.
- Produces: `pub fn body(ui: &mut egui::Ui, panel: &mut QcPanelState, task: &mut QcRunTask, opt_panel: &mut XtbPanelState, mol: &Molecule)` — the same signature as `ui_classic::body`, so the caller's match arms differ only in which module they name.

- [ ] **Step 1: Write the failing test**

Add to `src/qchem_panel/mod.rs`'s test block:

```rust
    /// Review Focus 4. Switching program must update every header at once. A
    /// header still showing the previous program's basis is worse than no
    /// header, because it reads as current.
    #[test]
    fn switching_program_updates_every_header() {
        use crate::qchem_interfaces::job::QcMethod;
        use crate::qchem_interfaces::program::QcProgram;
        use sections::{header, Section};

        let mut state = QcPanelState::default();
        state.set_program(QcProgram::Orca);
        state.method = Some(QcMethod::Dft);
        state.config.functional = "M062X".to_string();
        let before = header(Section::LevelOfTheory, &state);
        assert!(before.contains("M06-2X"), "{before}");

        state.set_program(QcProgram::PySCF);
        let after = header(Section::LevelOfTheory, &state);
        assert!(
            !after.contains("M06-2X"),
            "ORCA's spelling must not survive the switch: {after}"
        );
        assert!(header(Section::Program, &state).contains("PySCF"));
    }
```

- [ ] **Step 2: Run it to verify it fails**

Run: `taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1 switching_program_updates_every_header`
Expected: FAIL — `sections` is in scope from Task 1, so this compiles; it fails only if `set_program` leaves the functional behind. If it already passes, that is fine: record it as a regression guard and continue.

- [ ] **Step 3: Write the sectioned layout**

Create `src/qchem_panel/ui_sections.rs`:

```rust
//! The runner panel as six collapsed sections.
//!
//! Every rule this draws lives in [`super::sections`], which has no drawing
//! dependency and carries the tests. This file only draws.
//!
//! Two properties are the point, and both are easy to lose in a later edit:
//!
//! * the whole setup reads without opening anything, because each header
//!   carries its own current setting;
//! * nothing appears or disappears, because a section that does not apply
//!   stays on screen and says why. A form whose rows come and go moves under
//!   the cursor, which is the complaint this layout answers.

use bevy::prelude::*;
use bevy_egui::egui;

use super::run::QcRunTask;
use super::sections::{self, Section};
use super::QcPanelState;
use crate::molecule::Molecule;
use crate::qchem_interfaces::valence::validate_electronic_state;
use crate::qchem_interfaces::xtb_optimize::{python_environment_row, XtbPanelState};

pub fn body(
    ui: &mut egui::Ui,
    panel: &mut QcPanelState,
    task: &mut QcRunTask,
    opt_panel: &mut XtbPanelState,
    mol: &Molecule,
) {
    let running = task.is_running();

    ui.label(egui::RichText::new(panel.summary_line()).strong());
    ui.separator();

    let mut python_ready = true;
    let mut level_ok = true;

    for section in sections::ALL {
        let available = sections::applies(section, panel);
        let heading = format!(
            "{}  \u{2014}  {}",
            sections::title(section),
            sections::header(section, panel)
        );

        egui::CollapsingHeader::new(heading)
            .id_salt(format!("qc_section_{section:?}"))
            // Only the Program section starts open, so a panel opened fresh
            // shows the whole setup at a glance rather than one expanded wall.
            .default_open(section == Section::Program)
            .show(ui, |ui| {
                if let Some(note) = sections::unavailable_note(section, panel) {
                    ui.label(egui::RichText::new(note).small().weak());
                    return;
                }
                match section {
                    Section::Program => {
                        super::ui_classic::program_row(ui, panel, opt_panel, running);
                        python_ready = python_environment_row(ui, panel.program, running);
                    }
                    Section::LevelOfTheory => {
                        level_ok = super::ui_classic::level_of_theory(ui, panel, running);
                    }
                    Section::Job => {
                        super::ui_classic::job_row(ui, panel, running);
                        if panel.needs_excited_state_options() {
                            super::ui_classic::excited_state_row(ui, panel, running);
                        }
                    }
                    Section::Molecule => {
                        super::ui_classic::charge_and_multiplicity(ui, panel, running);
                        electronic_state(ui, panel, mol);
                    }
                    Section::Resources => {
                        super::ui_classic::resources_row(ui, panel, running);
                    }
                    Section::InputFile => {
                        super::ui_classic::wavefunction_row(ui, panel, running);
                        super::ui_classic::input_row(ui, panel, mol, running);
                    }
                }
                let _ = available;
            });
    }

    let electronic = validate_electronic_state(mol, panel.charge, panel.multiplicity);

    ui.separator();
    super::ui_classic::run_row(
        ui,
        panel,
        task,
        opt_panel,
        mol,
        running,
        python_ready,
        level_ok,
        electronic.is_ok(),
    );
    super::ui_classic::results_row(ui, panel, task);
}

/// The charge and spin state's own verdict, inside the Molecule section.
///
/// On screen rather than in a tooltip: a run refused for an impossible
/// electronic state is something the user has to act on, and a reason they
/// cannot reach is no reason at all.
fn electronic_state(ui: &mut egui::Ui, panel: &QcPanelState, mol: &Molecule) {
    match validate_electronic_state(mol, panel.charge, panel.multiplicity) {
        Err(error) => {
            ui.colored_label(
                egui::Color32::from_rgb(230, 120, 90),
                format!("\u{26a0} {error}"),
            );
        }
        Ok((_, warnings)) if panel.show_warnings => {
            for warning in warnings {
                ui.label(egui::RichText::new(format!("\u{26a0} {}", warning.message)).weak());
            }
        }
        Ok(_) => {}
    }
}
```

In `src/qchem_panel/mod.rs` add `pub mod ui_sections;`.

In `src/qchem_panel/ui_classic.rs`, change these helpers from private to `pub(crate)` so this file can reuse them rather than duplicating: `program_row`, `charge_and_multiplicity`, `level_of_theory`, `job_row`, `excited_state_row`, `wavefunction_row`, `resources_row`, `input_row`, `run_row`, `results_row`.

In `src/qchem_panel/ui.rs`, point the arm at the new module:

```rust
                PanelLayout::Sections => super::ui_sections::body(ui, panel, task, opt_panel, mol),
```

- [ ] **Step 4: Run the whole suite**

Run: `taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1`
Expected: every test passes.

- [ ] **Step 5: Build and check for warnings**

Run: `taskset -c 0-3 cargo build --jobs 4 2>&1 | grep -A4 "qchem_panel/"`
Expected: no output.

Do not commit. Report and stop.

---

### Task 4: Promote the hover-only information

Three facts a user needs in order to act currently live only in hover text. On a mouse that is merely inconvenient; the spec forbids it because tablets come later and retrofitting is expensive. Move them on screen.

**Files:**
- Modify: `src/qchem_panel/ui_classic.rs` (`job_row`, `input_row`)
- Modify: `src/qchem_panel/summary.rs` (`optimization_buttons`)

**Interfaces:**
- Consumes: `crate::qchem_interfaces::job::unsupported_reason`.
- Produces: no new interfaces; behaviour change only.

- [ ] **Step 1: Write the failing test**

Add to `src/qchem_panel/mod.rs`'s test block:

```rust
    /// A job the chosen program cannot do must explain itself on screen.
    ///
    /// The reason currently lives in a hover tooltip. A user who cannot see
    /// why an option is unavailable is left guessing, and a tooltip is
    /// unreachable the moment this panel is used with a finger.
    #[test]
    fn an_unavailable_job_explains_itself_without_hovering() {
        use crate::qchem_interfaces::job::{self, JobType};
        use crate::qchem_interfaces::program::QcProgram;

        // The model already carries the reason; this test pins that it is
        // available to draw on screen rather than only as hover text.
        let reason = job::unsupported_reason(QcProgram::Xtb, JobType::Irc)
            .expect("xTB has no IRC");
        assert!(reason.contains("ORCA") || reason.contains("Psi4"), "{reason}");

        // And that every refusal is short enough to show inline.
        for program in QcProgram::ALL {
            for a_job in JobType::ALL {
                if let Some(text) = job::unsupported_reason(program, a_job) {
                    assert!(
                        text.chars().count() < 200,
                        "{program:?}/{a_job:?} is too long to show inline: {} chars",
                        text.chars().count()
                    );
                }
            }
        }
    }
```

- [ ] **Step 2: Run it to verify it fails or passes**

Run: `taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1 an_unavailable_job_explains_itself`
Expected: PASS if the reasons are already short, FAIL naming any that are too long. If any fail, shorten that reason in `src/qchem_interfaces/job.rs` and rerun.

- [ ] **Step 3: Show the reason under the job selector**

In `src/qchem_panel/ui_classic.rs`, at the end of `job_row`, after the
description label, add:

```rust
    // Why the jobs this program cannot do are unavailable, on screen rather
    // than behind a hover. A user who can see an option but not why it is
    // refused is left guessing, and hover text does not exist on a touch
    // screen, which these panels are meant to reach later.
    //
    // Deduplicated: several jobs share one reason, and repeating it once per
    // job would bury the panel in its own explanations.
    let mut refused: Vec<&'static str> = JobType::ALL
        .into_iter()
        .filter_map(|other| job::unsupported_reason(panel.program, other))
        .collect();
    refused.sort_unstable();
    refused.dedup();
    for reason in refused {
        ui.label(egui::RichText::new(reason).small().weak());
    }
```

`job_row` already has `job` and `JobType` in scope from its own use of them.

- [ ] **Step 4: Show why the optimisation summary is unavailable**

In `src/qchem_panel/summary.rs`, in `optimization_buttons`, the disabled Optimization Summary button carries its reason in `on_disabled_hover_text`. Keep that, and add the same sentence as a visible weak label beneath the row when `has_report` is false:

```rust
    if !has_report {
        ui.label(
            egui::RichText::new(
                "No optimization summary: this program printed no per-cycle convergence table \
                 Beavyr can read. The energy plot and the cycle table above come from the \
                 trajectory instead.",
            )
            .small()
            .weak(),
        );
    }
```

- [ ] **Step 5: Run the whole suite and build**

Run: `taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1`
Expected: every test passes.

Run: `taskset -c 0-3 cargo build --jobs 4 2>&1 | tail -2`
Expected: it finishes.

Do not commit. Report the final test count and stop.

---

## After this plan

The builder is stage 2 and gets its own plan, written once this one has been
used for real. That ordering is deliberate: the builder is a rewrite of a
6254-line file, and what is learned from living with the sectioned runner
should shape it.
