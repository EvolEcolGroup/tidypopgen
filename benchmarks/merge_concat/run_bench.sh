#!/usr/bin/env bash
# ./run_bench.sh <threads> <rep> [scenarios]  (appends to results.tsv)
# scenarios: a regular expression (default: all of concat, concat_shuffled,
# merge_pops and merge_sgdp)
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"

B="$here/env/bcftools/bin"
LIB="$here/lib"
T=$1
REP=$2
SCENARIOS="${3:-.}"
want() { [[ "$1" =~ $SCENARIOS ]]; }
D="$here/data"
OUT="$here/out"
RES="$here/results.tsv"
[ -f "$RES" ] || echo -e "scenario\tthreads\trep\tstep\twall_s\tinner_s\tpeak_rss_mb\tn_indiv\tn_loci\tr_heap_mb" > "$RES"

timed() { # label, command...
  local label=$1; shift
  local log="$OUT/$label.time"
  /usr/bin/time -f "%e %M" -o "$log" "$@" > "$OUT/$label.out" 2> "$OUT/$label.err"
  read -r wall rss < "$log"
  local inner=NA ni=NA nl=NA heap=NA
  if grep -q '^STEP' "$OUT/$label.out"; then
    read -r _ inner ni nl heap < <(grep '^STEP' "$OUT/$label.out")
  fi
  echo -e "$SC\t$T\t$REP\t$label\t$wall\t$inner\t$((rss / 1024))\t$ni\t$nl\t$heap" >> "$RES"
  echo "$SC T=$T rep=$REP $label: ${wall}s, $((rss / 1024)) MB"
}

run_scenario() { # scenario, bcftools subcommand + options, gt mode, files dir
  SC=$1; local bcf_args=$2 mode=$3 dir=$4
  rm -rf "$OUT"; mkdir -p "$OUT"
  local files=("$D/$dir"/*.vcf.gz)
  cat "${files[@]}" > /dev/null # warm the page cache
  # the old route: bcftools, then gen_tibble() on its output
  # shellcheck disable=SC2086
  timed bcftools "$B/bcftools" $bcf_args --threads "$T" -Oz -o "$OUT/combined.vcf.gz" "${files[@]}"
  timed gen_tibble Rscript "$here/gt_step.R" "$LIB" gen_tibble "$T" "$OUT/gen_tibble" "$OUT/combined.vcf.gz"
  # the new route
  timed "$mode" Rscript "$here/gt_step.R" "$LIB" "$mode" "$T" "$OUT/new" "${files[@]}"
  # same gen_tibble? (samples, loci, ploidy, and the bytes of the backing
  # file). Different gen_tibbles stop the benchmark (as does a comparison
  # that fails, e.g. on a missing file, through set -e)
  local md5_old md5_new same_meta=0
  md5_old=$(md5sum < "$OUT/gen_tibble.bk")
  md5_new=$(md5sum < "$OUT/new.bk")
  Rscript -e "a <- readRDS('$OUT/gen_tibble_meta.rds'); b <- readRDS('$OUT/new_meta.rds'); quit(status = if (identical(a, b)) 0 else 3)" \
    || same_meta=$?
  if [ "$same_meta" -ne 0 ] && [ "$same_meta" -ne 3 ]; then
    echo "$SC: could not compare the gen_tibbles" >&2
    return 1
  fi
  if [ "$same_meta" -eq 0 ] && [ "$md5_old" == "$md5_new" ]; then
    echo "$SC: identical gen_tibbles"
  else
    echo "$SC: DIFFERENT gen_tibbles" >&2
    return 1
  fi
}

if want concat; then
  run_scenario concat "concat" concat concat
fi
if want concat_shuffled; then
  # (bcftools concat needs the samples in the same order in all files)
  SC=concat_shuffled; rm -rf "$OUT"; mkdir -p "$OUT"
  cat "$D/concat_shuffled"/*.vcf.gz > /dev/null
  timed concat Rscript "$here/gt_step.R" "$LIB" concat "$T" "$OUT/new" "$D/concat_shuffled"/*.vcf.gz
fi
if want merge_pops; then
  run_scenario merge_pops "merge" merge_intersection merge_pops
fi
if want merge_sgdp; then
  run_scenario merge_sgdp "merge --missing-to-ref" merge_union_ref merge_sgdp
fi
