#!/usr/bin/env bash
# Publish a built site into the gh-pages branch, retrying if another workflow pushed first.
#
#   publish-pages.sh SRC_DIR DEST     DEST is a folder (e.g. lite) or "." for the site root.
#
# A folder DEST is replaced as a whole.  The root is replaced too, except for the folders
# other workflows publish (PRESERVE, default "lite coupled-model").
# Needs GITHUB_TOKEN, GITHUB_REPOSITORY and GITHUB_SHA (set by GitHub Actions).
set -euo pipefail
SRC=$(cd "$1" && pwd); DEST=$2
PRESERVE=${PRESERVE:-"lite coupled-model"}
REPO=${PAGES_REPO:-"https://x-access-token:${GITHUB_TOKEN}@github.com/${GITHUB_REPOSITORY}.git"}   # PAGES_REPO: for local tests
WORK=$(mktemp -d)

for attempt in 1 2 3 4 5; do
  rm -rf "$WORK/site"
  git clone --quiet --depth 1 --branch gh-pages "$REPO" "$WORK/site"
  cd "$WORK/site"
  if [ "$DEST" = "." ]; then
    keep=(-name .git); for d in $PRESERVE; do keep+=(-o -name "$d"); done
    find . -mindepth 1 -maxdepth 1 ! \( "${keep[@]}" \) -exec rm -rf {} +
  else
    rm -rf "$DEST"; mkdir -p "$DEST"
  fi
  cp -a "$SRC"/. "$DEST"/
  touch .nojekyll
  git add -A
  if git diff --cached --quiet; then echo "gh-pages already up to date"; exit 0; fi
  git -c user.name="github-actions[bot]" -c user.email="41898282+github-actions[bot]@users.noreply.github.com" \
      commit --quiet -m "deploy ${DEST} from ${GITHUB_SHA}"
  if git push --quiet origin gh-pages; then echo "published ${DEST} (attempt ${attempt})"; exit 0; fi
  echo "push rejected (another publish landed first); retrying" >&2
  cd "$WORK"; sleep $((attempt * 5))
done
echo "could not publish after 5 attempts" >&2
exit 1
