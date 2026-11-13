#!/bin/bash
# Time gen_tibble() on a .ped file, for several versions of tidypopgen.
#   run_bench.sh <ped> <work_dir> <reps> <label>=<lib|worktree>:<n_cores> ...
# e.g. run_bench.sh data/synthetic.ped /tmp/w 3 old=/tmp/base_lib:1 \
#        new=~/tidypopgen:1 new4=~/tidypopgen:4
# Prints a line per run: label, rep, elapsed (s), peak RSS (MB), peak R heap
# (MB). Check that the machine is idle first (uptime, top).
ped=$1
work=$2
reps=$3
shift 3
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$work"
echo -e "label\trep\telapsed_s\tpeak_rss_mb\tpeak_heap_mb\tload"
for rep in $(seq 1 "$reps"); do
  for spec in "$@"; do
    label=${spec%%=*}
    rest=${spec#*=}
    src=${rest%:*}
    cores=${rest##*:}
    bk="$work/${label}_$rep"
    rm -f "$bk".bk "$bk".rds "$bk".gt "$bk"_meta.rds
    load=$(cut -d' ' -f1 /proc/loadavg)
    out=$( { /usr/bin/time -f "MAXRSS %M" Rscript "$here/bench_step.R" \
      "$src" "$ped" "$bk" "$cores"; } 2>&1 )
    step=$(echo "$out" | grep '^STEP')
    rss=$(echo "$out" | grep '^MAXRSS' | cut -d' ' -f2)
    elapsed=$(echo "$step" | cut -d' ' -f2)
    heap=$(echo "$step" | cut -d' ' -f5)
    echo -e "$label\t$rep\t$elapsed\t$((rss / 1024))\t$heap\t$load"
    rm -f "$bk".bk
  done
done
