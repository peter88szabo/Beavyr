# Packaging Beavyr for Linux

Everything here builds one file — `beavyr_0.2.0_amd64.deb` — that you can mail
to a colleague or put on a share. It installs Beavyr the way any other Ubuntu
program installs, and removes the same way.

## For whoever receives the file

**Install** — double-click it, or:

```sh
sudo apt install ./beavyr_0.2.0_amd64.deb
```

Afterwards Beavyr is in the applications menu with its icon, and `beavyr` works
in any terminal.

**Remove** — the Uninstall button in Ubuntu's App Center, or:

```sh
sudo apt remove beavyr
```

That takes every file the package installed. Nothing is left behind, and
nothing was ever written outside the paths listed below.

**Requirements** — Ubuntu 22.04 or newer, 64-bit, and a machine with working
graphics drivers. Beavyr does not need Python, and does not need xTB or ORCA to
start: those are only wanted when you actually run a calculation, and you point
Beavyr at them from inside the program.

## For whoever builds it

```sh
./packaging/build-deb.sh              # builds the binary, then the package
./packaging/build-deb.sh --no-build   # package an already-built release binary
```

The result lands in `target/deb/`. The script needs only `cargo`, `dpkg-deb`
and `strip`, all of which are already on an Ubuntu machine.

## What the package puts on the system

| File | Path |
|---|---|
| the program | `/usr/bin/beavyr` |
| icon | `/usr/share/icons/hicolor/scalable/apps/beavyr.svg` |
| menu entry | `/usr/share/applications/beavyr.desktop` |
| App Center entry | `/usr/share/metainfo/be.kuleuven.Beavyr.metainfo.xml` |
| licence, README | `/usr/share/doc/beavyr/` |

Five files, all owned by the package manager. That is the whole installation.

## The files here

| File | What it is |
|---|---|
| `build-deb.sh` | builds the package |
| `beavyr.svg` | the application icon — the beaver in profile, on a pale tile |
| `beavyr.desktop` | the menu entry: name, icon, category |
| `be.kuleuven.Beavyr.metainfo.xml` | the App Center description |
| `icons/` | the six icon candidates and the page for comparing them |

The design decisions behind all of this — why a `.deb` and not a Flatpak, why
the icon ships as SVG rather than PNG, and why the graphics libraries have to
be listed by hand — are in
`docs/superpowers/specs/2026-09-08-debian-package-design.md`.

---

## A note on the size of this repository

`.git` is 366 MB, against 1.8 MB of source. Almost all of it is binary ORCA
scratch files under `examples/` — `.gbw`, `.cis`, `.densities`, `.mkl` — 355 MB
across 19 committed files, none of which Beavyr can read. There is no parser
for any of those formats.

They are now in `.gitignore`, which stops **new** ones being added. It does not
remove the ones already committed. Two further steps exist, and both are your
call because they have consequences:

**1. Stop tracking them.** The repository stops growing, and the files stay on
your disk untouched:

```sh
git rm --cached $(git ls-files examples | grep -E '\.(gbw|cis|densities|densitiesinfo|mkl|cpcm|cpcm_corr|carthess|bibtex)$')
git commit -m "Stop tracking ORCA scratch files nothing reads"
```

A fresh clone still downloads all 366 MB, because the old versions remain in
the history. Only new clones stop getting the files in their working copy.

**2. Remove them from the history as well.** This is what actually shrinks the
repository, to roughly 11 MB, using `git filter-repo` and a force push. It
rewrites every commit, so anyone else who has cloned the repository has to
delete their copy and clone again. With a private repository and one author
that is usually fine; with collaborators it is not.

**Do not ignore `examples/` as a whole.** 97 of the 628 tests read the text
files in it — geometries, ORCA outputs, Hessians — at run time. Removing the
folder from the repository makes those 97 fail for anyone who clones it. This
was measured, not assumed: hiding the folder and running `cargo test` gives
`531 passed; 97 failed`.
