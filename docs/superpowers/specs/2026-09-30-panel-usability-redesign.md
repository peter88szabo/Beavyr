# Making the builder and the runner usable

**Date:** 2026-09-30
**Status:** design, awaiting review
**Stages:** the runner first, then the builder, each behind its own switch

## The brief

Two panels are hard to use: the molecule builder and the general
quantum-chemistry runner. Asked what actually goes wrong, the answer was three
things, and one telling omission.

What goes wrong:

* **Too much on screen at once.** Every control is visible whether or not it
  applies, so the window has to be tall to be usable, and rows appearing and
  disappearing move things under the cursor.
* **Losing track of the mode.** In the builder, add-atom, add-fragment and the
  fragment editor all want the same mouse clicks and nothing says which is
  listening.
* **Undo is unreliable.** The builder has a separate Undo per operation rather
  than one history, so pressing one does not reverse what was actually last
  done.

What was *not* chosen: "too many steps for a routine job". So nothing in this
design tries to save clicks. Every change here buys orientation or
reversibility, and where the two conflict, orientation wins.

One requirement given directly: **add fragment and replace fragment must be
separate actions**, not one control with a mode switch.

## The shared principle

Both panels today are flat lists of every control that could ever apply. The
principle replacing that is: **show the state, and only the actions that state
supports.**

For the builder that means the selection is the state. For the runner it means
the current setup is the state. Neither panel should ever hold something the
user cannot see.

## Stage 1: the runner

Six collapsed sections, each showing its current setting on its own header.
Only the section being worked in is open.

```
  Program          ORCA  (~/Programs/orca)        >
  Level of theory  B3LYP/def2-SVP + D4            >
  Job              Optimization                   >
  Molecule         neutral singlet                >
  Resources        4 threads, 4 GB                >
  Input file       built from the settings above  >
  ------------------------------------------------
                  [ Run ]
```

Two properties matter more than the tidiness.

**The whole setup reads without opening anything.** Today the level of theory
is legible only if you scroll to it; here it is on the header. This is what
makes the panel answerable at a glance, which the eleven-row form is not.

**Nothing appears or disappears.** Today a row vanishes when it stops applying
and the form jumps. Here a section that does not apply says so on its own
header, in place. A force field's level-of-theory header reads "none: a force
field", rather than the section being absent.

The sections map onto what already exists, so this is a rearrangement rather
than new behaviour:

| Section | Holds |
|---|---|
| Program | the program, its path or its Python environment |
| Level of theory | method, functional, basis, dispersion, extra keywords |
| Job | the job, and the excited-state options when the job is TD-DFT |
| Molecule | charge, multiplicity, and the valence warnings |
| Resources | memory, threads |
| Input file | the orbital toggle, Save input, and the direct editor |

Run, the status line and the summary button stay below the sections, always
visible, because they are what the panel is for.

## Stage 2: the builder

**Clicking an atom always means select it.** There are no tools and no modes,
so there is no hidden state to lose track of, and the selection is drawn in the
viewport so it is never invisible.

```
  Selection: atoms 3, 7                     [Clear]
  -------------------------------------------------
    [ Add atom here ]        [ Add fragment ]
    [ Replace with fragment ][ Delete these ]
    [ Set distance ]         [ Rotate about 3-7 ]

  Pick a third atom to set an angle.
  -------------------------------------------------
  History                      [ Undo ]  [ Redo ]
    4. replaced atom 12 with OH
    3. added CH3 at atom 7
    2. deleted atom 12
    1. set 3-7 to 1.54 A
```

What is offered follows from what is selected, and nothing else is shown:

| Selected | Offers |
|---|---|
| nothing | add the first atom; load or clear |
| one atom | add here, add fragment, replace with fragment, delete |
| two atoms | the above, plus set distance and rotate about that bond |
| three | plus set angle |
| four | plus set dihedral |

An action that cannot apply is **absent**, not greyed, and a line says what to
pick next to reach more. That line is the replacement for a greyed button: it
tells you how to proceed rather than only that you cannot.

**Add fragment and replace fragment are two buttons.** The current panel makes
them one control with a mode dropdown, which is precisely the hidden state this
design exists to remove.

The fragment editor's rotate, translate and bend become actions on a selection
like any other. They lose their own sub-panel and their own three undo buttons.

### History

Every action is one committed step with a name that can be read. Undo and Redo
walk that list, and the list is visible, so what is about to be undone can be
seen before pressing. That visibility is the cure: undo feels unreliable today
because several separate undos exist and pressing one does not reverse the last
thing done.

This history is in memory and lasts the session. It holds 200 steps, which is
far more than a building session produces and small enough that the memory it
costs is irrelevant; past that the oldest are dropped. Loading a different
structure clears it, because undoing across a structure change would restore
atoms that no longer belong to anything.

**Recent Structures does not change and is never called history.** It stays
what it is: twenty milestones on disk, surviving a restart. The two answer
different questions -- take back that last click, versus go back to what I had
this morning -- and the only real risk is naming them alike, so they are not.

## The switch

Each new panel replaces its old one behind a toggle, and the old layout stays
one click away until it is no longer wanted.

* A setting with two values per panel, defaulting to the new layout.
* A plainly worded control in each panel: "use the old layout".
* The old code stays where it is, called from the other branch of the toggle,
  and is deleted only when the new one has earned it.

The choice is remembered between sessions, through the same config-file
mechanism the program paths already use. The existing classic/windowed layout
switch is *not* remembered today, which is a small wart; this one should be,
because being dropped back into a layout you rejected is worse than the wart.

## What is not changing

No capability is removed from either panel. Everything both can do today, they
can still do; only how it is reached changes.

No wizard anywhere. Clicks are not the problem, and a wizard would buy
orientation at a price that was not asked for.

The viewport, the trajectory player, the Vibrations and UV-Vis panels, and
every backend are untouched.

## A constraint held in reserve

These panels should later work from a tablet. That is explicitly out of scope
here -- this design targets a PC with a mouse -- but one habit is worth keeping
now because it is cheap now and expensive later: **no load-bearing information
lives only in a hover tooltip.**

Today a good deal does: why a job is disabled for a program, which Python
interpreter was searched, why the optimisation summary button is greyed. On a
touch screen there is no hover, so all of it would have to move. Where this
design introduces such text it goes on screen, and where an existing tooltip
carries something the user needs in order to act, it is promoted as the
sections are rebuilt.

## What this costs

The runner is a few hundred lines rearranged. Low risk, and it proves the
sectioned pattern.

The builder is not. `builder_ui.rs` is 6254 lines with two functions of roughly
a thousand lines each, and the picking logic, the three competing modes and the
per-operation undos are woven through all of it. Selection-first is a rewrite
of that file, split into several smaller ones. It is the largest single piece
of work discussed for this project, and it is why the runner goes first and the
old builder stays reachable.

## Testing

The rules each panel follows -- which actions a selection supports, what a
section's header says, whether a program's level of theory applies -- go in
plain functions with tests, away from the drawing code, as
`qchem_interfaces::job` and the method catalogues already do. That is what
makes "with two atoms you can set a distance" something a test can assert
rather than something only a click can reveal.

The history gets tests for the property that makes it trustworthy: that undo
followed by redo returns the structure unchanged, for every action type.
