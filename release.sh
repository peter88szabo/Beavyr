#!/usr/bin/env bash
# Edit this version, then run ./release.sh (Git and Python 3 are required).
VERSION="0.2.1"
REMOTE="origin"

# The existing .github/workflows/release.yml chooses the operating systems:
# Ubuntu .deb, Arch .pkg.tar.zst, Windows setup.exe, and universal macOS .dmg.
# This script commits changes and pushes a tag; GitHub builds the installers.
# ./release.sh --dry-run previews the release without changing anything.
set -euo pipefail

fail() { printf 'Error: %s\n' "$*" >&2; exit 1; }
case "${1:-}" in
    '') [[ $# == 0 ]] || fail 'Unexpected arguments.'; dry_run=false ;;
    --dry-run) [[ $# == 1 ]] || fail 'Unexpected arguments.'; dry_run=true ;;
    --help|-h)
        printf '%s\n' 'Edit VERSION at the top, then run ./release.sh.' \
            'Use ./release.sh --dry-run to preview without writing or pushing.' \
            'Includes tracked changes and new source, tests, packaging, and CI files.' \
            'Other new files must be staged with git add before running this script.'
        exit 0 ;;
    *) fail 'Usage: ./release.sh [--dry-run|--help]' ;;
esac
command -v git >/dev/null || fail 'Install Git first.'
command -v python3 >/dev/null || fail 'Install Python 3 first.'

# Restrict this shell and its children to at most four available Linux cores.
python3 - <<'PY'
import os
if hasattr(os, 'sched_setaffinity'):
    pid = os.getppid()
    os.sched_setaffinity(pid, sorted(os.sched_getaffinity(pid))[:4])
PY
export CARGO_BUILD_JOBS=4 RAYON_NUM_THREADS=4 OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
# Git's compression and index workers also respect the four-thread limit.
git() { command git -c pack.threads=4 -c index.threads=4 "$@"; }
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
[[ -f Cargo.toml && -f Cargo.lock && -f .github/workflows/release.yml ]] ||
    fail 'Keep release.sh in the Beavyr repository root.'
[[ "$(git rev-parse --show-toplevel)" == "$PWD" ]] || fail 'Not at the repository root.'
branch=$(git symbolic-ref --quiet --short HEAD) || fail 'Check out a branch first.'
tag="v$VERSION"
git check-ref-format "refs/tags/$tag" >/dev/null || fail 'Invalid version tag.'
git remote get-url "$REMOTE" >/dev/null || fail "Remote $REMOTE is missing."
[[ -z "$(git ls-files --unmerged)" ]] || fail 'Resolve merge conflicts first.'
for state in MERGE_HEAD CHERRY_PICK_HEAD REVERT_HEAD rebase-merge rebase-apply; do
    [[ ! -e "$(git rev-parse --git-path "$state")" ]] || fail 'Finish the current Git operation first.'
done

# Validate both files before changing either. Update only the root package's
# version in Cargo.lock; dependency versions and checksums stay untouched.
update_versions() {
    python3 - "$VERSION" "$1" <<'PY'
import pathlib
import re
import sys

version, mode = sys.argv[1:]
if not re.fullmatch(r'(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)', version):
    sys.exit('VERSION must have three numbers, for example 0.2.1.')
if any(int(part) > 65535 for part in version.split('.')):
    sys.exit('Each VERSION component must be <= 65535 for the Windows installer.')
manifest = pathlib.Path('Cargo.toml')
lock = pathlib.Path('Cargo.lock')
manifest_text, lock_text = manifest.read_text(), lock.read_text()
package = re.search(r'(?ms)^\[package\]\s*\n(.*?)(?=^\[|\Z)', manifest_text)
if package is None:
    sys.exit('Cargo.toml has no [package] section.')
name = re.search(r'(?m)^name\s*=\s*"([^"]+)"', package[1])
if name is None:
    sys.exit('Cannot find the package name.')
pattern = r'(?m)^(version\s*=\s*")[^"]+("[^\n]*)$'
new_package, count = re.subn(pattern, lambda m: m[1] + version + m[2], package[1])
if count != 1:
    sys.exit('Expected exactly one package version in Cargo.toml.')
entries = list(re.finditer(r'(?ms)^\[\[package\]\]\s*\n(.*?)(?=^\[\[package\]\]|\Z)', lock_text))
roots = [m for m in entries if re.search(r'(?m)^name = "' + re.escape(name[1]) + r'"$', m[1])
         and not re.search(r'(?m)^source\s*=', m[1])]
if len(roots) != 1:
    sys.exit('Expected exactly one local root package in Cargo.lock.')
root = roots[0]
new_root, count = re.subn(pattern, lambda m: m[1] + version + m[2], root[1])
if count != 1:
    sys.exit('Expected exactly one root package version in Cargo.lock.')
if mode == 'write':
    manifest.write_text(manifest_text[:package.start(1)] + new_package + manifest_text[package.end(1):])
    lock.write_text(lock_text[:root.start(1)] + new_root + lock_text[root.end(1):])
PY
}
update_versions check
if git show-ref --verify --quiet "refs/tags/$tag"; then
    fail "Local tag $tag already exists. Choose a new version, or retry the previous push as printed on failure."
fi
printf 'Release %s from branch %s via %s\n' "$tag" "$branch" "$REMOTE"
printf '%s\n' 'Installers: Ubuntu, Arch Linux, Windows, macOS (Intel + Apple Silicon).' \
    'Will update Cargo.toml/Cargo.lock, commit tracked edits and new project files, and push branch + tag.'
# Do not automatically add unrelated untracked research data or machine-local scripts.
paths=(src tests packaging .github .cargo release.sh)
git status --short --untracked-files=normal
if "$dry_run"; then
    printf '%s\n' 'Dry run complete. No files, commits, tags, or remote refs were changed.'
    exit 0
fi

git var GIT_AUTHOR_IDENT >/dev/null || fail 'Configure your Git name and email first.'
git var GIT_COMMITTER_IDENT >/dev/null || fail 'Configure your Git name and email first.'
# Check credentials, remote tags, and branch ancestry before editing or committing.
remote_refs=$(git ls-remote "$REMOTE" "refs/heads/$branch" "refs/tags/$tag") ||
    fail 'Cannot read the remote. Check your GitHub connection and credentials.'
if [[ "$remote_refs" == *$'\t'"refs/tags/$tag"* ]]; then
    fail "Remote tag $tag already exists; choose a new VERSION."
fi
if [[ "$remote_refs" == *$'\t'"refs/heads/$branch"* ]]; then
    git fetch --no-tags "$REMOTE" "refs/heads/$branch"
    git merge-base --is-ancestor FETCH_HEAD HEAD ||
        fail 'The remote branch has changes you do not have. Integrate them before releasing.'
fi
update_versions write
git add -u -- .
for path in "${paths[@]}"; do
    [[ ! -e "$path" ]] || git add -A -- "$path"
done
if ! git diff --cached --quiet; then
    git commit -m "Release $tag"
fi
git tag -a "$tag" -m "Beavyr $VERSION"
# Atomic push prevents publishing only the branch or only the tag. Never force.
if ! git push --atomic "$REMOTE" "HEAD:refs/heads/$branch" "refs/tags/$tag:refs/tags/$tag"; then
    printf '\nThe local release commit/tag were kept. After fixing the push error, retry:\n' >&2
    printf '  git push --atomic %q %q %q\n' "$REMOTE" "HEAD:refs/heads/$branch" "refs/tags/$tag:refs/tags/$tag" >&2
    exit 1
fi
printf '\nPushed %s. Follow the Release workflow in the GitHub Actions tab.\n' "$tag"
printf '%s\n' 'When all builds succeed, the installers appear on the GitHub Releases page.'
