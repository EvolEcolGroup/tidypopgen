#!/usr/bin/env bash
# generate the sgdp and many scenarios, then bgzip all VCFs (niced, one core)
set -euo pipefail
cd /media/andrea/annex/tmp/claude/tidypopgen/merge_vcf
for sc in sgdp many; do
  nice -n 19 taskset -c 3 Rscript scripts/make_merge_data.R data $sc
done
for f in data/*/*.vcf; do
  nice -n 19 taskset -c 3 bgzip -f "$f"
done
du -sh data/*
echo ALL_DONE
