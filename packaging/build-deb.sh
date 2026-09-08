#!/usr/bin/env bash
#
# Build beavyr_<version>_amd64.deb — one file to hand to a colleague.
#
# Needs nothing that is not already on an Ubuntu machine: cargo, dpkg-deb and
# strip. Run it from anywhere; it finds the project itself.
#
#   ./packaging/build-deb.sh            build the package
#   ./packaging/build-deb.sh --no-build use the release binary already built
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
VERSION="$(grep -m1 '^version' "$ROOT/Cargo.toml" | cut -d'"' -f2)"
ARCH="amd64"
PKG="beavyr"
STAGE="$ROOT/target/deb/${PKG}_${VERSION}_${ARCH}"
OUT="$ROOT/target/deb/${PKG}_${VERSION}_${ARCH}.deb"

if [[ "${1:-}" != "--no-build" ]]; then
  echo "==> building release binary (this takes a few minutes)"
  ( cd "$ROOT" && cargo build --release --no-default-features )
fi

BIN="$ROOT/target/release/Beavyr"
[[ -x "$BIN" ]] || { echo "no release binary at $BIN" >&2; exit 1; }

echo "==> laying out the package"
rm -rf "$STAGE"
install -d "$STAGE/DEBIAN"
install -d "$STAGE/usr/bin"
install -d "$STAGE/usr/share/applications"
install -d "$STAGE/usr/share/icons/hicolor/scalable/apps"
install -d "$STAGE/usr/share/metainfo"
install -d "$STAGE/usr/share/doc/$PKG"

# Cargo names the binary after the crate, with a capital B. On a command line
# that is awkward, so it goes in as lowercase.
install -m 0755 "$BIN" "$STAGE/usr/bin/beavyr"
strip --strip-unneeded "$STAGE/usr/bin/beavyr" 2>/dev/null || \
  echo "    (strip unavailable — shipping the binary unstripped)"

install -m 0644 "$HERE/beavyr.desktop"      "$STAGE/usr/share/applications/beavyr.desktop"
install -m 0644 "$HERE/beavyr.svg"          "$STAGE/usr/share/icons/hicolor/scalable/apps/beavyr.svg"
install -m 0644 "$HERE/be.kuleuven.Beavyr.metainfo.xml" "$STAGE/usr/share/metainfo/be.kuleuven.Beavyr.metainfo.xml"
install -m 0644 "$ROOT/LICENSE"             "$STAGE/usr/share/doc/$PKG/copyright"
install -m 0644 "$ROOT/README.md"           "$STAGE/usr/share/doc/$PKG/README.md"

INSTALLED_KB="$(du -sk "$STAGE" | cut -f1)"

# The first five Depends were read off the built binary with ldd and dpkg -S,
# not guessed. The X11 pair, and everything in Recommends, is opened by name at
# run time so it cannot be found that way; those are a judgement call. See
# docs/superpowers/specs/2026-09-08-debian-package-design.md.
cat > "$STAGE/DEBIAN/control" <<EOF
Package: $PKG
Version: $VERSION
Section: science
Priority: optional
Architecture: $ARCH
Maintainer: Peter Szabo <peter.szabo@kuleuven.be>
Installed-Size: $INSTALLED_KB
Depends: libc6 (>= 2.35), libgcc-s1, libasound2t64, libcap2, libudev1, libx11-6, libxkbcommon0
Recommends: libvulkan1, mesa-vulkan-drivers, libwayland-client0, libxcursor1, libxrandr2, libxi6
Description: Molecular viewer and quantum-chemistry front-end
 Beavyr displays molecules and trajectories at interactive frame rates, drives
 external quantum-chemistry programs, and reads back what they produce:
 geometries, Hessians, orbitals and excited states.
 .
 It computes normal modes with Eckart projection, thermochemistry with
 Grimme's quasi-RRHO treatment, IR and UV-Vis spectra, and molecular-orbital
 and density isosurfaces. It reads XYZ, Molden, ORCA .hess and .engrad, and
 Gaussian frequency output, and can drive xTB and Behemoth directly.
 .
 To remove it again: sudo apt remove beavyr
EOF

# The icon cache and the desktop database are only refreshed if the tools are
# present; on a machine without them the menu entry still works.
cat > "$STAGE/DEBIAN/postinst" <<'EOF'
#!/bin/sh
set -e
if [ "$1" = "configure" ]; then
    if command -v gtk-update-icon-cache >/dev/null 2>&1; then
        gtk-update-icon-cache -q -f /usr/share/icons/hicolor 2>/dev/null || true
    fi
    if command -v update-desktop-database >/dev/null 2>&1; then
        update-desktop-database -q /usr/share/applications 2>/dev/null || true
    fi
fi
exit 0
EOF

cat > "$STAGE/DEBIAN/postrm" <<'EOF'
#!/bin/sh
set -e
if [ "$1" = "remove" ] || [ "$1" = "purge" ]; then
    if command -v gtk-update-icon-cache >/dev/null 2>&1; then
        gtk-update-icon-cache -q -f /usr/share/icons/hicolor 2>/dev/null || true
    fi
    if command -v update-desktop-database >/dev/null 2>&1; then
        update-desktop-database -q /usr/share/applications 2>/dev/null || true
    fi
fi
exit 0
EOF

chmod 0755 "$STAGE/DEBIAN/postinst" "$STAGE/DEBIAN/postrm"

echo "==> building the .deb"
dpkg-deb --build --root-owner-group "$STAGE" "$OUT" >/dev/null

echo
echo "    $OUT"
echo "    $(du -h "$OUT" | cut -f1)"
echo
echo "    install:  sudo apt install $OUT"
echo "    remove:   sudo apt remove beavyr"
