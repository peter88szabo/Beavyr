#!/usr/bin/env bash
#
# Assemble Beavyr.app and wrap it in a .dmg. Runs on a Mac -- GitHub's macOS
# runner -- because `iconutil` and `hdiutil` are macOS-only.
#
#   build-app.sh <version> <universal-binary> <icons-dir> <output-dir>
#
set -euo pipefail

VERSION="$1"
BINARY="$2"
ICONS="$3"
OUT="$4"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

APP="$OUT/Beavyr.app"
rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"

install -m 0755 "$BINARY" "$APP/Contents/MacOS/beavyr"

# macOS wants its icon as a .icns, which is a bundle of PNG sizes.
ICONSET="$OUT/beavyr.iconset"
rm -rf "$ICONSET"; mkdir -p "$ICONSET"
cp "$ICONS/16.png"   "$ICONSET/icon_16x16.png"
cp "$ICONS/32.png"   "$ICONSET/icon_16x16@2x.png"
cp "$ICONS/32.png"   "$ICONSET/icon_32x32.png"
cp "$ICONS/64.png"   "$ICONSET/icon_32x32@2x.png"
cp "$ICONS/128.png"  "$ICONSET/icon_128x128.png"
cp "$ICONS/256.png"  "$ICONSET/icon_128x128@2x.png"
cp "$ICONS/256.png"  "$ICONSET/icon_256x256.png"
cp "$ICONS/512.png"  "$ICONSET/icon_256x256@2x.png"
cp "$ICONS/512.png"  "$ICONSET/icon_512x512.png"
cp "$ICONS/1024.png" "$ICONSET/icon_512x512@2x.png"
iconutil -c icns "$ICONSET" -o "$APP/Contents/Resources/beavyr.icns"

sed "s/VERSION_PLACEHOLDER/$VERSION/g" "$HERE/Info.plist" > "$APP/Contents/Info.plist"

# Ad-hoc signature. This is NOT Apple notarisation and does not stop Gatekeeper
# warning about an unidentified developer -- that needs a paid Apple developer
# account. What it does do is make the app launchable at all on Apple Silicon,
# where an unsigned binary is refused outright rather than merely warned about.
codesign --force --deep --sign - "$APP"

# A .dmg with a shortcut to /Applications, which is how a Mac user expects to
# install something: drag the left icon onto the right one.
STAGE="$OUT/dmg"
rm -rf "$STAGE"; mkdir -p "$STAGE"
cp -R "$APP" "$STAGE/"
ln -s /Applications "$STAGE/Applications"
cp "$HERE/INSTALL.txt" "$STAGE/HOW TO OPEN THIS.txt"

hdiutil create -volname "Beavyr $VERSION" \
    -srcfolder "$STAGE" -ov -format UDZO \
    "$OUT/beavyr-$VERSION-macos.dmg"

echo "built $OUT/beavyr-$VERSION-macos.dmg"
