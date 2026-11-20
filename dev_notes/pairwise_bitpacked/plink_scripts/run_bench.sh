#!/bin/bash
# Balanced timings of KING-robust on HGDP650: PLINK 2 vs tidypopgen.
# Usage: run_bench.sh <threads> <rounds order, e.g. "F R F">
# Each task runs in its own process; rounds alternate the order of the tasks.
# Results: out/bench_<threads>.tsv (PLINK) and out/tpg_<threads>.tsv
# (tidypopgen), with the load before each run in out/bench_log.txt.
set -u
P=${TPG_WORK:-/media/andrea/Elements1/tmp/claude/tidypopgen}/plink_king
LIB=${TPG_WORK:-/media/andrea/Elements1/tmp/claude/tidypopgen}/bench_inb/lib_new
PL=$P/bin/plink2
DATA=$P/data/hgdp650_pop
T=$1
shift
ROUNDS="$*"
if [ "$T" = 1 ]; then CPUS=2; else CPUS=0-$((T - 1)); fi
export TMPDIR=$P/rtmp OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=1
OUT=$P/out/run_$T
mkdir -p $OUT
LOG=$P/out/bench_log.txt
PLTSV=$P/out/bench_plink_$T.tsv
TPGTSV=$P/out/bench_tpg_$T.tsv
PLOPT="--threads $T --memory 2000"
POPS=$(ls $P/data/keep | sed 's/\.txt$//')

load() { uptime | sed 's/.*average: //; s/,//g' | cut -d' ' -f1; }

# time a PLINK command: label, then the command
tpl() {
  local label=$1
  shift
  local l=$(load)
  local tf=$OUT/time.txt
  /usr/bin/time -f "%e %U %S" -o $tf taskset -c $CPUS "$@" > /dev/null 2>&1
  printf "%s\t%s\t%d\t%s\t%s\n" "$(date +%H:%M:%S)" "$label" "$T" \
    "$(tr ' ' '\t' < $tf)" "$l" >> $PLTSV
  echo "# $(date -Is) $label T=$T load=$l: $(cat $tf)" >> $LOG
  sleep 3
}

pl_bed_full() { tpl pl_bed_full $PL --bfile $DATA --make-king triangle bin $PLOPT --out $OUT/bed_full; }
pl_pgen_full() { tpl pl_pgen_full $PL --pfile $DATA --make-king triangle bin $PLOPT --out $OUT/pgen_full; }
pl_bed_within() { tpl pl_bed_within $PL --bfile $DATA --make-king-table rel-check $PLOPT --out $OUT/bed_within; }
pl_pgen_within() { tpl pl_pgen_within $PL --pfile $DATA --make-king-table rel-check $PLOPT --out $OUT/pgen_within; }
pl_convert() { tpl pl_convert $PL --bfile $DATA --make-pgen $PLOPT --out $OUT/conv; }
pl_load() { tpl pl_load $PL --pfile $DATA --write-snplist $PLOPT --out $OUT/load; }
pl_read() { tpl pl_read $PL --pfile $DATA --geno-counts $PLOPT --out $OUT/read; }
keeploop() {
  for p in $POPS; do
    $PL --pfile $DATA --keep $P/data/keep/$p.txt --make-king triangle bin $PLOPT --out $OUT/keep_$p
  done
}
export -f keeploop
pl_keeploop() { tpl pl_pgen_keeploop env PL=$PL DATA=$DATA PLOPT="$PLOPT" P=$P OUT=$OUT POPS="$POPS" bash -c keeploop; }

ttpg() {
  local task=$1
  local l=$(load)
  echo "# $(date -Is) tpg_$task T=$T load=$l" >> $LOG
  taskset -c $CPUS Rscript $P/scripts/time_tpg.R $LIB $task $T new $TPGTSV > /dev/null 2>> $LOG
  sleep 3
}
tpg_full() { ttpg full; }
tpg_full_kernel() { ttpg full_kernel; }
tpg_grouped() { ttpg grouped; }
tpg_pack() { ttpg pack; }

TASKS="pl_pgen_full tpg_full pl_bed_full tpg_full_kernel pl_pgen_within tpg_grouped pl_bed_within tpg_pack pl_keeploop pl_convert pl_load pl_read"
RTASKS=$(echo $TASKS | tr ' ' '\n' | tac | tr '\n' ' ')

echo "# $(date -Is) $(hostname) T=$T CPUS=$CPUS rounds: $ROUNDS; $($PL --version); tidypopgen $(git -C ${TPG_REPO:-/media/andrea/annex/git/tidypopgen} rev-parse --short HEAD) (lib_new)" >> $LOG
# warm the page cache
cat $DATA.bed $DATA.pgen $DATA.pvar $DATA.bim ${TPG_WORK:-/media/andrea/Elements1/tmp/claude/tidypopgen}/bench_inb/data/hgdp650.bk > /dev/null
for r in $ROUNDS; do
  if [ "$r" = F ]; then list=$TASKS; else list=$RTASKS; fi
  for t in $list; do $t; done
done
echo "# $(date -Is) done T=$T" >> $LOG
