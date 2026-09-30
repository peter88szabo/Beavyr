# SDD ledger — plan: docs/superpowers/plans/2026-09-30-sectioned-quantum-chemistry-runner.md

Setup: Ruling: no git worktree, work in place on master — Peter tests the built
binary from the main tree throughout this session and a worktree would take
that away — cost if wrong: none, nothing is committed so `git checkout` reverts
everything.

Setup: Ruling: no commits at any step — Peter commits his own code, stated
explicitly this session and held in memory; the plan's Global Constraints say
the same — cost if wrong: none, he commits when he is satisfied.

Pre-flight (shared interfaces):
- Task 1 produces sections::{Section, ALL, title, header, applies,
  unavailable_note, MAX_HEADER_CHARS}; Task 3 consumes all of them. Names and
  arities agree between the Produces block and Task 3's call sites. Clean.
- Task 2 produces ui_classic::body(ui, panel, task, opt_panel, mol) and
  PanelLayout; Task 3 produces ui_sections::body with the identical signature
  and Task 2's caller switches on PanelLayout. Signatures agree. Clean.
- Task 2 makes ten row helpers pub(crate); Task 3 calls exactly those ten.
  Checked name by name against ui.rs as it stands. Clean.
- Task 4 consumes job::unsupported_reason, which exists today and is unchanged
  by Tasks 1-3. Clean.

Task 1: complete (no commit per Setup ruling; tests: cargo test -- --test-threads=1 -> 1192 passed )
  7 new tests in qchem_panel::sections, all watched RED before GREEN.

Task 2: RETRACTED completion line — written before the suite output was read,
  and the suite failed: ui.rs's test block referenced helpers that had moved to
  ui_classic.rs. Corrected below.
  3 new tests. Old layout moved verbatim to ui_classic.rs with 11 helpers made
  pub(crate); PanelLayout added and persisted; both switch arms draw classic.
Task 2: Ruling: PanelLayout derives Default with #[default] on Sections rather
  than the plan's hand-written impl — same behaviour, and it lets the test
  assert the default through the trait instead of through QcPanelState, which
  would have needed a config file on disk — cost if wrong: none.
Task 2: Note: sections.rs is dead code until Task 3 consumes it; 10 warnings
  expected to clear there.
Task 2: complete (no commit per Setup ruling; tests: cargo test -- --test-threads=1 -> 1195 passed, 0 failed)
Task 2: Ruling: ui.rs's test block moved to ui_classic.rs along with the code
  it exercises — the plan said to move the helpers but not the tests, and tests
  belong with the helpers they call — cost if wrong: none, they are the same
  tests in a different file.
Task 3: Step 2 test `switching_program_updates_every_header` passed on first run (set_program already resets the functional). Plan anticipated this: recorded as a regression guard, not a RED step.
Task 3: complete (no commit per Peter's standing instruction; files src/qchem_panel/ui_sections.rs new, mod.rs + ui.rs wired). Build clean, zero warnings under qchem_panel/. Tests: taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1 -> 1196 passed; 0 failed; 22 ignored.
Task 4: Step 2 test `an_unavailable_job_explains_itself_without_hovering` passed first run (every refusal already under 200 chars). Recorded as a guard, as the plan allowed.
Task 4: Ruling: also promoted two hovers the plan's steps did not name -- the Energy Plot button's "no per-cycle energies" and input_row's "Rebuilding discards anything typed in the box". Why: the plan listed input_row in its Files block but wrote no step for it, and one visible reason beside two greyed buttons reads as if it explained both. The rebuild warning guards a destructive click. Cost if wrong: three extra weak labels Peter can delete in one line each.
Task 4: Ruling: the job_row and input_row edits land in ui_classic.rs despite that file's "nothing new goes in this file" header. Why: both rows are shared -- ui_sections calls them -- so this is a change to a shared control, not a new capability in the old layout. Cost if wrong: none; the alternative would duplicate the rows.
Task 4: complete (no commit per Peter's standing instruction). Tests: taskset -c 0-3 cargo test --jobs 4 -- --test-threads=1 -> 1197 passed; 0 failed; 22 ignored. Build finishes, zero warnings under qchem_panel/.
Final: review package built from the working tree, not a commit range, because Peter commits his own code and nothing here is committed.
Final: Peter lifted the no-commit rule for this change mid-run ("I allow you to commit the changes if you are done"). Commit after the fix pass.
Final: self-found defect ahead of the review -- ui_sections.rs sets python_ready inside the Program CollapsingHeader body, which egui does not run when collapsed, so Run is enabled for a Python backend the environment cannot import. Fix queued for the single fix pass.
Final: pre-existing observation -- ui_classic::level_of_theory returns true on both exit paths (304, 423) despite documenting a real verdict, so the level_ok gate has never done anything. run.rs:202 catches a missing method at run time, so no user-visible harm.
Final: fixed the collapsed-section gate -- python_ready and level_ok now computed in ui_sections::body outside the CollapsingHeader loop; new pure predicate xtb_optimize::python_environment_ready, which python_environment_row now delegates to so row and gate cannot drift. Test python_readiness_is_answerable_without_drawing added (a guard, not RED->GREEN: the defect is in UI structure and is not unit-testable; the test pins the predicate the fix depends on). Suite 1198 passed, 0 failed, 22 ignored.
