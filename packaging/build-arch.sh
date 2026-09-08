#!/usr/bin/env bash
#
# Build beavyr-<version>-1-x86_64.pkg.tar.zst — the Arch Linux equivalent of
# build-deb.sh.
#
# Must be run on an Arch system (or in an archlinux container), because it
# needs makepkg. Everything else it needs is in base-devel.
#
#   ./packaging/build-arch.sh            build the package
#   ./packaging/build-arch.sh --no-build use the release binary already built
#
# Note that makepkg refuses to run as root. In a container, create an ordinary
# user with passwordless sudo and run this as them; the release workflow does
# exactly that.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
VERSION="$(grep -m1 '^version' "$ROOT/Cargo.toml" | cut -d'"' -f2)"
OUTDIR="$ROOT/target/arch"

command -v makepkg >/dev/null || {
  echo "makepkg not found: this script builds an Arch package and has to run" >&2
  echo "on Arch. On another distro, use the release workflow or a container:" >&2
  echo "  docker run --rm -v \"\$PWD\":/src archlinux:base-devel /src/packaging/..." >&2
  exit 1
}
[[ "$(id -u)" -ne 0 ]] || {
  echo "makepkg refuses to run as root. Run this as an ordinary user." >&2
  exit 1
}

if [[ "${1:-}" != "--no-build" ]]; then
  echo "==> building release binary (this takes a few minutes)"
  ( cd "$ROOT" && cargo build --release --no-default-features )
fi

BIN="$ROOT/target/release/Beavyr"
[[ -x "$BIN" ]] || { echo "no release binary at $BIN" >&2; exit 1; }

echo "==> stripping the binary"
strip --strip-unneeded "$BIN" 2>/dev/null || \
  echo "    (strip unavailable — shipping the binary unstripped)"

# makepkg runs in packaging/ itself, so that the PKGBUILD's $startdir points
# where it expects: $startdir is packaging/ and $startdir/.. the project root.
# Its working directories are redirected into target/ so the repository is
# left unmodified.
echo "==> building the package"
rm -rf "$OUTDIR"
install -d "$OUTDIR" "$OUTDIR/work"

(
  cd "$HERE"
  PKGDEST="$OUTDIR" BUILDDIR="$OUTDIR/work" SRCDEST="$OUTDIR/work" \
    makepkg --force --nodeps --noconfirm
)

OUT="$(ls "$OUTDIR"/*.pkg.tar.zst)"
echo
echo "    $OUT"
echo "    $(du -h "$OUT" | cut -f1)"
echo
echo "    install:  sudo pacman -U $OUT"
echo "    remove:   sudo pacman -R beavyr"
