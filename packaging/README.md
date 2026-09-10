# Packaging Beavyr for Linux

Two packages are built here, one per family of distribution:

| file | for | installs with |
|---|---|---|
| `beavyr_0.2.0_amd64.deb` | Debian, Ubuntu, Mint | `apt` |
| `beavyr-0.2.0-1-x86_64.pkg.tar.zst` | Arch, Manjaro, EndeavourOS | `pacman` |

Either is one file you can mail to a colleague or put on a share, and either
installs and removes the way anything else on that distribution does. A `.deb`
cannot be installed on Arch and the Arch package cannot be installed on Ubuntu,
so hand over the one that matches.

## For whoever receives the .deb

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

## For whoever receives the Arch package

**Install**

```sh
sudo pacman -U ./beavyr-0.2.0-1-x86_64.pkg.tar.zst
```

**Remove**

```sh
sudo pacman -R beavyr
```

**Requirements** — a 64-bit Arch or Arch-derived system with working graphics
drivers. The package requires only `libcap` and `systemd-libs`,
which every desktop install already has. The Vulkan loader, a Vulkan driver
for your GPU and the X11 or Wayland libraries are listed as optional
dependencies because Beavyr opens them by name at run time rather than linking
them; in practice a working desktop has them, but `pacman` will name them if
something is missing.

## For whoever builds it

```sh
./packaging/build-deb.sh              # builds the binary, then the package
./packaging/build-deb.sh --no-build   # package an already-built release binary
```

The result lands in `target/deb/`. The script needs only `cargo`, `dpkg-deb`
and `strip`, all of which are already on an Ubuntu machine.

### The Arch package

```sh
./packaging/build-arch.sh             # builds the binary, then the package
./packaging/build-arch.sh --no-build  # package an already-built release binary
```

The result lands in `target/arch/`. This one has to run **on Arch**, because it
needs `makepkg`, and it has to run as an ordinary user, because `makepkg`
refuses to run as root. From another distribution, use a container:

```sh
docker run --rm -v "$PWD":/src -w /src archlinux:base-devel bash -c '
  pacman -Syu --noconfirm --needed git rust libcap systemd-libs pkgconf
  useradd -m builder && chown -R builder /src
  sudo -u builder ./packaging/build-arch.sh'
```

`PKGBUILD` is an *in-tree* one: it packages the release binary already built
from this repository rather than fetching a source tarball, because building
Bevy under `makepkg` takes about a quarter of an hour and the release workflow
has already built exactly the binary we want to ship. A submission to the AUR
would want a different PKGBUILD, with `source=()` pointing at a release tag and
a real `build()` running `cargo build --release`; the `depends` and
`optdepends` here transfer to it unchanged.

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
| `build-deb.sh` | builds the Debian/Ubuntu package |
| `build-arch.sh` | builds the Arch package |
| `PKGBUILD` | the Arch recipe `build-arch.sh` drives |
| `beavyr.svg` | the application icon — the beaver in profile, on a pale tile |
| `beavyr.desktop` | the menu entry: name, icon, category |
| `be.kuleuven.Beavyr.metainfo.xml` | the App Center description |
| `icons/` | the six icon candidates and the page for comparing them |

The design decisions behind all of this — why native packages and not a
Flatpak, why the icon ships as SVG rather than PNG, and why the graphics
libraries have to be listed by hand — are in
`docs/superpowers/specs/2026-09-08-debian-package-design.md`. The Arch package
follows the same reasoning and the same dependency split; only the package
manager differs.

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
