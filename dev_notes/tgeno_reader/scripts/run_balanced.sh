#!/bin/bash
# A B C D E, E D C B A, A B C D E
B=${1:?usage: run_balanced.sh <work dir with lib_base, lib_new, data/aadr_like.*>}
run() { rm -f $B/out_*; sync; Rscript $B/bench_one.R $1 $2 $3 gen_tibble $B/out_run; }
A="$B/lib_base $B/data/aadr_like.geno 1"
Bv="$B/lib_new $B/data/aadr_like.geno 1"
C="$B/lib_new $B/data/aadr_like.geno 4"
D="$B/lib_new $B/data/aadr_like.tgeno 1"
E="$B/lib_new $B/data/aadr_like.tgeno 4"
for order in "A Bv C D E" "E D C Bv A" "A Bv C D E"; do for v in $order; do run ${!v}; done; done
rm -f $B/out_*
