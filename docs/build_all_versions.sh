#!/usr/bin/env bash
# Builds the whole documentation site, every version in one folder:
#
#   dev       this checkout (main)
#   vX.Y      each release tag v*, as it was at that release
#   stable    the newest release
#
# plus versions.js (what the version picker in the top bar lists) and a
# root index.html that opens stable (dev before the first release).
#
#   docs/build_all_versions.sh [output folder, default docs/pages-site]
#
# The site is published in one piece (.github/workflows/Documenter.yml
# uploads this folder to GitHub Pages; there is no gh-pages branch — see
# docs/make.jl), so every version is rebuilt each time. A release's pages
# keep their own content, with this checkout's presentation: its theme
# (docs/src/.vitepress) and figure styling (docs/make_figures.jl).
set -euo pipefail

root="$(git rev-parse --show-toplevel)"
site="$(realpath -m "${1:-$root/docs/pages-site}")"
rm -rf "$site"
mkdir -p "$site"

# Documenter decides each build's version folder from these, as on GitHub
# Actions; the token and actor only have to be present.
export GITHUB_ACTIONS=true GITHUB_EVENT_NAME=push GITHUB_REPOSITORY=Richard7987/AsteroidPipeline.jl
export GITHUB_TOKEN="${GITHUB_TOKEN:-unused}" GITHUB_ACTOR="${GITHUB_ACTOR:-local}"

build() {   # build <checkout> <git ref>: compile its docs, copy each version folder into the site
    local dir="$1" ref="$2"
    echo "=== building the documentation of $ref"
    (cd "$dir" && GITHUB_REF="$ref" julia --project=docs -e \
        'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()' &&
     GITHUB_REF="$ref" julia --project=docs docs/make.jl)
    if [ ! -f "$dir/docs/build/bases.txt" ]; then
        echo "=== $ref: no version folder of its own (a newer release replaces it)"
        return
    fi
    local i=1 base
    while read -r base; do
        [ -n "$base" ] || { echo "error: $ref was built without a version folder" >&2; exit 1; }
        rm -rf "${site:?}/$base"
        cp -r "$dir/docs/build/$i" "$site/$base"
        i=$((i + 1))
    done < "$dir/docs/build/bases.txt"
}

build "$root" refs/heads/main

# Only the newest release of each series is built — each vX.Y below 1.0,
# each vX from 1.0 on — the folders DocumenterVitepress keeps (keep=:breaking):
# an older release in the same series would build no folder of its own.
releases="$(git -C "$root" tag -l 'v[0-9]*' --sort=v:refname |
    awk -F'[v.]' '{ series = ($2 == 0) ? $2 "." $3 : $2; last[series] = $0 }
                  END { for (s in last) print last[s] }' | sort -V)"

for tag in $releases; do
    wt="$(mktemp -d)"
    git -C "$root" worktree add -q --detach "$wt" "$tag"
    rm -rf "$wt/docs/src/.vitepress"
    cp -r "$root/docs/src/.vitepress" "$wt/docs/src/.vitepress"
    cp "$root/docs/make_figures.jl" "$wt/docs/make_figures.jl"
    build "$wt" "refs/tags/$tag"
    git -C "$root" worktree remove --force "$wt"
done

# The version picker's list: stable, then releases newest first, then dev.
versions=()
[ -d "$site/stable" ] && versions+=(stable)
while read -r v; do versions+=("$v"); done < <(ls "$site" | grep -E '^v[0-9]' | sort -rV)
versions+=(dev)
{
    echo 'var DOC_VERSIONS = ['
    for v in "${versions[@]}"; do echo "  \"$v\","; done
    echo '];'
    [ -d "$site/stable" ] && echo 'var DOCUMENTER_STABLE = "stable";'
} > "$site/versions.js"

home=dev
[ -d "$site/stable" ] && home=stable
printf '%s\n' "<meta http-equiv=\"refresh\" content=\"0; url=./$home/\"/>" > "$site/index.html"
touch "$site/.nojekyll"
echo "=== site in $site: ${versions[*]}"
