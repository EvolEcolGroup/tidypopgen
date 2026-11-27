#!/usr/bin/env bash
# Set up the merge/concat benchmark:
# - bcftools and htslib (bgzip, tabix) in a conda env in env/bcftools, made
#   with micromamba (downloaded into env/bin if it is not on the PATH);
# - tidypopgen installed with R CMD INSTALL (-O2, never load_all(), which
#   compiles the C++ at -O0) into lib/, from a git ref (HEAD by default; only
#   committed changes are installed).
#   ./setup.sh [ref]
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
repo="$(git -C "$here" rev-parse --show-toplevel)"
ref="${1:-HEAD}"

mamba_bin="$(command -v micromamba || true)"
if [ -z "$mamba_bin" ]; then
  mkdir -p "$here/env/bin"
  mamba_bin="$here/env/bin/micromamba"
  [ -x "$mamba_bin" ] || curl -fsSL -o "$mamba_bin" \
    https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
  chmod +x "$mamba_bin"
fi
if [ ! -x "$here/env/bcftools/bin/bcftools" ]; then
  MAMBA_ROOT_PREFIX="$here/env/root" "$mamba_bin" create -y -q \
    -p "$here/env/bcftools" -c conda-forge -c bioconda bcftools htslib
fi
"$here/env/bcftools/bin/bcftools" --version | head -1

rm -rf "$here/src" "$here/lib"
mkdir -p "$here/src" "$here/lib"
git -C "$repo" archive "$ref" | tar -x -C "$here/src"
git -C "$repo" rev-parse "$ref" > "$here/lib/COMMIT"
R CMD INSTALL --no-docs --no-multiarch --library="$here/lib" "$here/src" \
  > "$here/lib/install.log" 2>&1
so="$here/lib/tidypopgen/libs/tidypopgen.so"
echo "non-inlined Rcpp indexing calls (0 = no -O0 code):" \
  "$(objdump -d "$so" | grep -c 'call.*vectorI.*ixEm' || true)"
