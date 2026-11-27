#!/bin/bash
# reader-only A/B of TGENO: v1 (no huge pages) vs new (huge pages), 1 and 4 cores, interleaved
B=${1:?usage: run_ab.sh <work dir>}
run() { rm -f $B/out_*; sync; Rscript $B/bench_one.R $1 $B/data/aadr_like.tgeno $2 reader $B/out_run; }
for order in "v1 new" "new v1" "v1 new"; do for v in $order; do for c in 1 4; do lib=$B/lib_$v; run $lib $c; done; done; done
rm -f $B/out_*
