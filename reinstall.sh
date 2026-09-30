#!/usr/bin/env bash
#
# Rebuild Beavyr from the current source and reinstall it system-wide.
#
# Builds a fresh .deb with packaging/build-deb.sh, then installs it. Earlier
# this script only reinstalled whatever .deb was left in target/deb, so it
# installed an old build -- or failed outright once that file was gone.
#
# The build is held to four cores so it does not starve anything else running.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERSION="$(grep -m1 '^version' "$ROOT/Cargo.toml" | cut -d'"' -f2)"
DEB="$ROOT/target/deb/beavyr_${VERSION}_amd64.deb"

echo "==> building Beavyr ${VERSION} from the current source (four cores)"
CARGO_BUILD_JOBS=4 taskset -c 0-3 "$ROOT/packaging/build-deb.sh"

[[ -f "$DEB" ]] || { echo "the build did not produce $DEB" >&2; exit 1; }

echo "==> installing $DEB"
sudo apt install --reinstall -y "$DEB"

echo "==> done: run 'beavyr' or open it from the applications menu"
