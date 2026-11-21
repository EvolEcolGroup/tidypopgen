#!/bin/bash
# Balanced runs of breakdown_full.R: versions as name:lib:n_cores, in the
# order A B C, C B A, A B C; the input files are dropped from the page cache
# before each run (cold reads), and the .bk is deleted after each run
# usage: run_balanced.sh <bed> <names_as_int> <log> <version>...
bed=$1; nai=$2; log=$3; shift 3
vers=("$@"); rev=(); for ((i=${#vers[@]}-1; i>=0; i--)); do rev+=("${vers[i]}"); done
W=$(cd "$(dirname "$0")/.." && pwd)
for v in "${vers[@]}" "${rev[@]}" "${vers[@]}"; do
  IFS=: read name lib nc <<< "$v"
  sync; dd if=$bed iflag=nocache count=0 2>/dev/null; dd if=${bed%.bed}.bim iflag=nocache count=0 2>/dev/null
  echo "== $name $(cat $W/$lib/COMMIT) n_cores $nc $(date +%T) load $(cut -d' ' -f1-3 /proc/loadavg) free $(df -h ~ | awk 'NR==2{print $4}')" | tee -a $log
  OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores TMPDIR=$W/rtmp taskset -c 0-7 \
    Rscript $W/scripts/breakdown_full.R $W/$lib $bed $W/out $nc $nai 2>&1 | grep -E "elapsed|TOTAL" | tee -a $log
  rm -f $W/out/full.*
  sleep 20
done
