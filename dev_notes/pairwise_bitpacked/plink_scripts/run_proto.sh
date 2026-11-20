#!/bin/bash
# Balanced timings of the prototypes against the base build (tidypopgen only).
# Usage: run_proto.sh <threads> <tasks, e.g. "full grouped"> <rounds, e.g. "F R F">
# Versions: base (lib_new, b8255226), p1 (tiles + 64-word chunks), p2 (p1 +
# AVX2 popcount for KING), p2s (p2 with the AVX2 kernel turned off: tiles,
# block layout, scalar popcount).
set -u
P=${TPG_WORK:-/media/andrea/Elements1/tmp/claude/tidypopgen}/plink_king
T=$1
TASKS=$2
shift 2
ROUNDS="$*"
if [ "$T" = 1 ]; then CPUS=2; else CPUS=0-$((T - 1)); fi
export TMPDIR=$P/rtmp OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=1
LOG=$P/out/proto_log.txt
TSV=$P/out/proto_$T.tsv
declare -A LIBS=([base]=${TPG_WORK:-/media/andrea/Elements1/tmp/claude/tidypopgen}/bench_inb/lib_new
  [p1]=$P/build/lib_p1 [p2]=$P/build/lib_p2 [p2s]=$P/build/lib_p2)
VERS="base p1 p2 p2s"
RVERS=$(echo $VERS | tr ' ' '\n' | tac | tr '\n' ' ')
load() { uptime | sed 's/.*average: //; s/,//g' | cut -d' ' -f1; }
echo "# $(date -Is) $(hostname) T=$T CPUS=$CPUS tasks: $TASKS rounds: $ROUNDS; p1=$(git -C $P/wt rev-parse --short king_plink_ideas~1) p2=$(git -C $P/wt rev-parse --short king_plink_ideas)" >> $LOG
for r in $ROUNDS; do
  if [ "$r" = F ]; then list=$VERS; else list=$RVERS; fi
  for task in $TASKS; do
    for v in $list; do
      l=$(load)
      echo "# $(date -Is) $v $task T=$T load=$l" >> $LOG
      if [ $v = p2s ]; then noavx="TIDYPOPGEN_NO_AVX2=1"; else noavx="X_UNUSED=1"; fi
      env $noavx taskset -c $CPUS Rscript $P/scripts/time_tpg.R ${LIBS[$v]} $task $T $v $TSV > /dev/null 2>> $LOG
      sleep 3
    done
  done
done
echo "# $(date -Is) done T=$T" >> $LOG
