#!/usr/bin/env bash
# Time gt_merge_vcfs() phases for several versions, alternated A B B A A B.
#   run_timing.sh <scenario> <n_cores> <reps> <version:lib>...
# Each run is a fresh R process (merge_phases.R), pinned with taskset to
# cores 0..n_cores-1 with OMP_NUM_THREADS = n_cores; appends to
# results/phases.tsv (phase times) and results/runs.tsv (wall time, peak
# RSS, load, date, md5 of the .bk).
set -euo pipefail
W=${TPG_WORK:-/media/andrea/Elements1/tmp/claude/tidypopgen}/merge_vcf
SC=$1
T=$2
REPS=$3
shift 3
VERSIONS=("$@")
export TMPDIR=$W/rtmp
export OMP_NUM_THREADS=$T
CORES=${CORES:-"0-$((T - 1))"}
RUNS=$W/results/runs.tsv
[ -f "$RUNS" ] || echo -e "date\tversion\tscenario\tn_cores\trep\twall_s\tpeak_rss_mb\tload1\tbk_md5\tmode\tcycles_u\tinstr_u\tcommit" > "$RUNS"
MODE=${MERGE_MODE:-phases}
run_one() { # version:lib rep
  local v=${1%%:*} lib=${1#*:} rep=$2
  local out=$W/out/$SC/$v
  rm -rf "$out"; mkdir -p "$out"
  local load1
  load1=$(cut -d' ' -f1 /proc/loadavg)
  local extra=""
  [ "$MODE" == "real" ] && extra="real"
  MERGE_REP=$rep perf stat -x, -e cycles:u,instructions:u -o "$out/perf" \
    /usr/bin/time -f "%e %M" -o "$out/time" taskset -c "$CORES" \
    Rscript $W/scripts/merge_phases.R "$lib" "$SC" "$T" "$out/merged" \
    "$W/results/phases_${v}.tsv" $extra > "$out/log" 2>&1 || { cat "$out/log"; exit 1; }
  read -r wall rss < "$out/time"
  local md5
  md5=$(md5sum < "$out/merged.bk" | cut -c1-8)_$(md5sum < "$out/merged.loci.tsv" | cut -c1-8)
  # keep the samples and loci of this run, to compare the versions
  cp "$out/merged_meta.rds" "$W/results/meta/${SC}_${v}_T${T}_r${rep}.rds"
  local cyc ins commit
  cyc=$(grep "cycles:u" "$out/perf" | cut -d, -f1)
  ins=$(grep "instructions:u" "$out/perf" | cut -d, -f1)
  commit=$(cat "$lib/COMMIT" 2>/dev/null || echo "?")
  echo -e "$(date +%FT%T)\t$v\t$SC\t$T\t$rep\t$wall\t$((rss / 1024))\t$load1\t$md5\t$MODE\t$cyc\t$ins\t$commit" >> "$RUNS"
  echo "$v $SC T=$T rep=$rep: ${wall}s $((rss / 1024))MB load=$load1 md5=$md5"
  sleep 5
}
n=${#VERSIONS[@]}
for ((r = 1; r <= REPS; r++)); do
  if (( r % 2 == 1 )); then
    for ((i = 0; i < n; i++)); do run_one "${VERSIONS[$i]}" "$r"; done
  else
    for ((i = n - 1; i >= 0; i--)); do run_one "${VERSIONS[$i]}" "$r"; done
  fi
done
