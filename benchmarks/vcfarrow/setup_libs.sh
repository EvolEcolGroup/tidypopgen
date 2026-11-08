#!/usr/bin/env bash
# Install tidypopgen refs into their own libraries with R CMD INSTALL (-O2),
# never from load_all() (which compiles the C++ at -O0).
#   ./setup_libs.sh                    v1.0.0 and vcf_fast (the default)
#   ./setup_libs.sh [label=]ref ...    e.g. v0.4.4=5ab2687c9 v1.0.0 dev
# Each ref goes to lib/<label>/ (label defaults to the ref, with slashes
# turned into "_"). Refs are taken from the git repository this folder is in
# (or TIDYPOPGEN_REPO), as committed: uncommitted changes are not installed.
# For 0.4.4 use 5ab2687c9 (main, CRAN 0.4.4): the v0.4.4 tag is one commit
# before the version bump (DESCRIPTION 0.4.3.9006), with the same R/ and src/.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
repo="${TIDYPOPGEN_REPO:-$(git -C "$here" rev-parse --show-toplevel)}"
mkdir -p "$here/logs"
refs=("$@")
[ ${#refs[@]} -eq 0 ] && refs=(v1.0.0 vcf_fast)
for spec in "${refs[@]}"; do
  if [[ "$spec" == *=* ]]; then
    tag="${spec%%=*}"
    ref="${spec#*=}"
  else
    ref="$spec"
    tag="${ref//\//_}"
  fi
  src="$here/src/$tag"
  lib="$here/lib/$tag"
  rm -rf "$src" "$lib"
  mkdir -p "$src" "$lib"
  git -C "$repo" archive "$ref" | tar -x -C "$src"
  git -C "$repo" rev-parse "$ref" > "$lib/COMMIT"
  echo "== installing $ref ($(cat "$lib/COMMIT")) into $lib"
  R CMD INSTALL --no-docs --no-multiarch --no-test-load --library="$lib" "$src" \
    > "$here/logs/install_$tag.log" 2>&1
  so="$lib/tidypopgen/libs/tidypopgen.so"
  n_o0=$(objdump -d "$so" | grep -c 'call.*vectorI.*ixEm' || true)
  echo "   non-inlined Rcpp indexing calls (0 = no -O0 code): $n_o0"
  Rscript -e "library(tidypopgen, lib.loc = '$lib'); cat('   loaded', as.character(packageVersion('tidypopgen', lib.loc = '$lib')), '\n')"
done
