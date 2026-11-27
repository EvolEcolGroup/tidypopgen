#!/usr/bin/env bash
# Generate the synthetic VCFs (make_data.R), then compress them with bgzip and
# index them with tabix (bcftools merge needs indexed files).
#   ./prepare_data.sh [scale]
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
B="$here/env/bcftools/bin"
rm -rf "$here/data"
Rscript "$here/make_data.R" "$here/data" "${1:-1}"
for f in "$here"/data/*/*.vcf; do
  "$B/bgzip" -@ 8 -f "$f"
  "$B/tabix" -f -p vcf "$f.gz"
done
du -sh "$here"/data/*
