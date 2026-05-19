#!/bin/bash
set -euo pipefail

test -d .git || git init
test -f .gitmodules || { echo "ERROR: .gitmodules not found" >&2; exit 1; }

url_for_path() {
  local wanted_path="$1"
  git config -f .gitmodules --get-regexp '^submodule\..*\.path$' |
  while read -r key path; do
    if [ "$path" = "$wanted_path" ]; then
      local name="${key#submodule.}"
      name="${name%.path}"
      git config -f .gitmodules --get "submodule.${name}.url"
      return 0
    fi
  done
}

restore_submodule() {
  local sha="$1"
  local path="$2"
  local url

  url="$(url_for_path "$path" || true)"
  test -n "$url" || { echo "ERROR: no url for $path in .gitmodules" >&2; exit 1; }

  echo "==> Restoring $path @ $sha"
  mkdir -p "$(dirname "$path")"

  if [ -d "$path/.git" ] || [ -f "$path/.git" ]; then
    git -C "$path" remote get-url origin >/dev/null 2>&1 || git -C "$path" remote add origin "$url"
    git -C "$path" remote set-url origin "$url"
  else
    rm -rf "$path"
    git clone "$url" "$path"
  fi

  git -C "$path" fetch --tags origin '+refs/heads/*:refs/remotes/origin/*' '+refs/tags/*:refs/tags/*'
  git -C "$path" checkout --detach "$sha"

  git add "$path"
}

restore_submodule 984e3f194862b17916536b5fade40cba6e47a6fe third-party/cereal
restore_submodule eddb0241389718a23a42db6af5f0164b6e0139af third-party/google-benchmark
restore_submodule 52eb8108c5bdec04579160ae17225d66034bd723 third-party/google-test
restore_submodule 83edb60836d87cf1b406e8846b9059c03031e8f5 third-party/gperftools

git add .gitmodules
git submodule absorbgitdirs || true
git submodule status --recursive