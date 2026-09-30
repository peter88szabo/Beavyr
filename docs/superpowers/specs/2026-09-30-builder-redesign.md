# Molecule builder and fragment editor redesign

**Date:** 2026-09-30
**Status:** implemented on branch `builder-redesign`, awaiting Peter's test
**Supersedes:** the stage-2 sketch in the withdrawn panel-usability spec

## Why

The builder panel is too tall, and three of its parts are the reason: the
Z-matrix table, the bond thresholds, and the fragment editor's three stacked
sub-sections. Separately, atom insertion and fragment insertion are two panels
doing the same job, add and replace hide behind a Mode toggle, and there are
six Undo buttons so pressing one does not reverse what was last done.

**No collapsible sections.** An earlier attempt hid one-line rows behind
collapsible headers; it cost a click and saved no height. Everything visible in
the panel is a control. Height is won by moving whole features into their own
windows, not by folding rows away.

## The panel

```
Molecule Builder
Click atoms to select.   Selected: 3, 7  [Clear]
Insert
  What    [ CH3                    v ]
  Bond 1.540   Angle 109.5   Dihedral 180.0
  [ Add to selection ]  [ Replace selection ]
  [ Remove selected atoms ]
  [ Fragment editor... ] [ Z-matrix... ] [ Bond thresholds... ]
  [ Undo ] [ Redo ]   last: added CH3 at atom 7
```

## The six changes

1. **Z-matrix moves to its own window**, opened by a button. Behaviour
   unchanged; it is a move, not a rewrite.
2. **Atom and fragment insertion merge.** One `What` dropdown over the
   fragment list, with `Single atom...` as its first entry. Choosing it opens a
   small form (element dropdown or typed symbol, Use this / Cancel); afterwards
   the atom is inserted through the same path as any fragment, with the same
   Bond/Angle/Dihedral fields.
3. **Add and Replace are two buttons**, not a Mode toggle. `Replace selection`
   keeps today's Replace meaning: the selected atom is replaced by the fragment.
4. **The fragment editor is one form in one window**: axis, the atoms that will
   move, then rotate / distance / bend as three numbers with a live preview, and
   one Apply. Cancel restores. Its three Undo buttons go away because nothing is
   committed until Apply.
5. **One Undo and one Redo** over a single list of named steps (add atom, add
   fragment, replace, remove, fragment-editor apply, cleanup). A `last: ...`
   line says what Undo will take back before it is pressed.
6. **The top line always says what the next click does**, and lists the current
   selection, which is also drawn in the viewport.

## Not changing

No capability is removed. The viewport, the backends and every other panel are
untouched. Bond thresholds move to a window with no change to their behaviour.

## Reversibility

The work lands on branch `builder-redesign`. There is no in-panel layout
switch: reverting is `git checkout master`.

## Testing

The rules go in plain functions with tests away from the drawing code: which
actions a selection supports, what the next-click line reads, and that undo
followed by redo returns the structure unchanged for every step type.

## What was built

All six changes, one commit each so any of them can be bisected or reverted
on its own:

| Commit | Change |
|---|---|
| `dfb1900` | Z-matrix table into its own window, still fully editable |
| `251c17f` | one Insert section; Add and Replace as two buttons |
| `0cd632e` | fragment editor as one window, one form, live preview |
| `1df470c` | one Undo and Redo over every action |
| `72a4a36` | the line saying what the next click will do |

Bond thresholds moved to their own window as part of the fragment editor
commit; they used to be nested *inside* the fragment editor, which is not
where anyone would look for them.

Suite: 1195 passing, 0 failing, 22 ignored. The only warnings under
`molecule_builder/` are `min_clearance` and `centroid_excluding`, both already
unused before this work.

### Capabilities kept

Nothing was removed. "Add Atom (auto)" is what "Add to selection" does for a
single atom; picking bond, angle and dihedral references by hand is still
there as "Pick references..."; the Z-matrix edits exactly as before.
