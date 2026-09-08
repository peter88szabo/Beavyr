# A Debian Package for Beavyr — Design

Date: 2026-09-08
Status: approved for implementation

## Goal

One file — `beavyr_0.1.0_amd64.deb` — that a colleague can be sent by mail or
download from a share, install by double-clicking, and remove again from
Ubuntu's own software list. After installing they get Beavyr in the
applications menu with an icon, and `beavyr` as a command in any terminal.

No Flatpak and no Snap. Both sandbox the program, and Beavyr's whole job is
launching xTB and Behemoth from wherever they happen to live on the user's
machine or their cluster mount. A sandbox turns that into a permanent fight for
no benefit here, where the audience is a research group on Ubuntu.

## Removal is a requirement, not an afterthought

This is the reason to build a real package rather than a script that copies
files. Removal is:

```
sudo apt remove beavyr          # or the App Center's Uninstall button
```

`dpkg` owns every file the package placed, so removal takes all of them and
leaves nothing behind. A shell installer cannot promise that, which is why one
is not being written.

## What goes where

| File | Path |
|---|---|
| binary | `/usr/bin/beavyr` |
| icon | `/usr/share/icons/hicolor/scalable/apps/beavyr.svg` |
| menu entry | `/usr/share/applications/beavyr.desktop` |
| App Center entry | `/usr/share/metainfo/beavyr.metainfo.xml` |
| licence and notes | `/usr/share/doc/beavyr/` |

The binary is installed as lowercase `beavyr`. Cargo builds it as `Beavyr`,
matching the crate name, which is awkward to type at a shell prompt; the
package renames it on the way in and the display name stays capitalised.

### The icon is shipped as SVG, not PNG

`hicolor/scalable/` is the standard place for a vector icon and every current
desktop renders it at whatever size it needs. This avoids shipping eight PNG
renditions, and avoids requiring a rasteriser on the build machine — this one
has neither `rsvg-convert` nor Inkscape installed. If a desktop is ever found
that will not draw it, adding PNG sizes is a change to the build script only.

The icon is candidate 01, the beaver in profile, drawn on a pale rounded tile.
The tile is deliberate: application menus are dark on some desktops and light
on others, and a dark outline on transparency vanishes against the first kind.

## Dependencies

The release binary links only `libc6`, `libgcc-s1` and `libm` — the Rust
standard library is linked statically once `--no-default-features` turns off
Bevy's dynamic linking, so no Bevy shared object is shipped or needed.

Everything to do with graphics and windowing is loaded at run time by name
rather than linked, so `dpkg` cannot discover it from the binary and it has to
be declared by hand:

* **Depends**: `libc6`, `libgcc-s1`, `libx11-6`, `libxkbcommon0` — without a
  window the program cannot start at all.
* **Recommends**: `libvulkan1`, `mesa-vulkan-drivers`, `libwayland-client0` —
  present on any desktop Ubuntu; a headless or unusual machine may legitimately
  lack them, and `Recommends` says so without blocking installation.

Because these are opened by name and not linked, the list is a judgement made
from how Bevy and wgpu behave, not something `dpkg-shlibdeps` derived. It is
written here so the next person knows it was a decision.

## Build

`packaging/build-deb.sh`, which:

1. builds `cargo build --release --no-default-features`;
2. strips the binary — a Bevy release binary carries a large symbol table that
   nothing in normal use needs;
3. lays out a staging tree in the shape above;
4. writes the `control` file with the size and dependencies;
5. calls `dpkg-deb --build --root-owner-group`.

`dpkg-deb` is part of a base Ubuntu install, so the script needs nothing that
is not already on the machine. No `cargo-deb`, which would have to be compiled
first and adds a build dependency for no gain at this size.

## File associations — deliberately left out for now

The menu entry declares no `MimeType`, so double-clicking an `.xyz` file will
not open Beavyr yet. That would only work once Beavyr reads a filename from its
command line, which is specified separately in
`2026-09-08-open-a-file-from-the-command-line-design.md` and is not built. When
it is, the association is three lines in the `.desktop` file plus a
`update-desktop-database` trigger.

Claiming the association first would be worse than not having it: the file
would open a program that ignores it.

## Out of scope

* An apt repository, signing, or automatic updates. The deliverable is a file
  to send to colleagues; hosting is a different question.
* Anything but `amd64`.
* Packaging Behemoth. It is a separate program, and being a command-line tool
  it would have no menu entry or App Center listing anyway.
