#!/usr/bin/env bash
# The cost of the locus sidecar (foo.loci.tsv, written by
# gen_tibble(names_as_int = TRUE) since commit f88a81e8), against the
# baseline without it (88f4ba82, the tip of vcf_fast). Install both first:
#   ./setup_libs.sh base=88f4ba82 sidecar=f88a81e8
#   ./sidecar_bench.sh
# Settings (environment variables):
#   SIZES       numbers of variants (VCFs made by run_all.sh in data/)
#   DEEP_SIZES  sizes at which lookups, restore, update and the compression
#               of the sidecar are measured too
#   REPS        repeats of each read (alternating the order of the versions;
#               the results have 3 more at 2.3 M variants, run by hand with
#               sidecar_bench.py as in run() below, reps 4-6)
#   DEEP_REPS   repeats of restore, update and the cold lookup
#   THREADS     n_cores of gen_tibble() (4)
#   MEMORY_MAX  memory limit of each R process (50G)
#   MAX_BUSY    wait before each measurement until at most this many cores
#               are busy with other work (3.5)
# Each measurement is a fresh R process (sidecar_bench.py runs
# sidecar_step.R), one at a time. Results: results/sidecar_cost.tsv (one row
# per timed step), results/sidecar_compress.tsv (sidecar_compress.cpp, and
# zstd and gzip as reference points) and results/sidecar_chromosomes.tsv
# (size of the planned chromosome sidecar). The outputs of each size are
# deleted once measured.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results work
SIZES=${SIZES:-"2322460 4644921 18579683"}
DEEP_SIZES=${DEEP_SIZES:-"2322460 18579683"}
REPS=${REPS:-3}
DEEP_REPS=${DEEP_REPS:-2}
THREADS=${THREADS:-4}
MEMORY_MAX=${MEMORY_MAX:-50G}
MAX_BUSY=${MAX_BUSY:-3.5}
RESULTS=${RESULTS:-results/sidecar_cost.tsv}
COMP_RESULTS=${COMP_RESULTS:-results/sidecar_compress.tsv}
CHR_RESULTS=${CHR_RESULTS:-results/sidecar_chromosomes.tsv}
log="logs/sidecar_bench_$(date +%Y%m%d_%H%M).log"
exec > >(tee -a "$log") 2>&1

echo "sidecar_bench.sh $(date): SIZES=$SIZES DEEP_SIZES=$DEEP_SIZES" \
  "REPS=$REPS DEEP_REPS=$DEEP_REPS THREADS=$THREADS MEMORY_MAX=$MEMORY_MAX"
uptime
free -g
for v in base sidecar; do
  [ -d "lib/$v/tidypopgen" ] || { echo "lib/$v missing: run setup_libs.sh"; exit 1; }
  so="lib/$v/tidypopgen/libs/tidypopgen.so"
  n_o0=$(objdump -d "$so" | grep -c 'call.*vectorI.*ixEm' || true)
  echo "$v: commit $(cat "lib/$v/COMMIT"), -O0 calls: $n_o0"
  [ "$n_o0" -eq 0 ] || { echo "lib/$v has -O0 code"; exit 1; }
done

# the compression benchmark (standalone C++, see sidecar_compress.cpp)
comp_bin=work/sidecar_compress/sidecar_compress
mkdir -p work/sidecar_compress
g++ -O2 -Wall -fopenmp -o "$comp_bin" sidecar_compress.cpp -lz -ldeflate

run() { # run <label> <n> <action> <rep> [extra args of sidecar_bench.py]
  local label=$1 n=$2 action=$3 rep=$4
  shift 4
  python3 sidecar_bench.py --lib "lib/$label" --label "$label" \
    --vcf "data/discosnp_${n}.vcf" --work "work/sidecar_${label}_${n}" \
    --results "$RESULTS" --action "$action" --rep "$rep" \
    --threads "$THREADS" --memory-max "$MEMORY_MAX" --max-busy "$MAX_BUSY" "$@"
}

# zstd and gzip on the sidecar, as reference points (single-threaded, on
# P-cores; median of 3)
cli_compress() { # cli_compress <n> <sidecar> <out_dir>
  python3 - "$1" "$2" "$3" "$COMP_RESULTS" <<'EOF'
import os, statistics, subprocess, sys, time
n, src, out_dir, results = sys.argv[1:]
tools = [("zstd-1", ["zstd", "-1", "-T1", "-q", "-c"], ["zstd", "-dc"]),
         ("zstd-3", ["zstd", "-3", "-T1", "-q", "-c"], ["zstd", "-dc"]),
         ("gzip-1", ["gzip", "-1", "-c"], ["gzip", "-dc"])]
pin = ["taskset", "-c", "0-7"]
def timed(cmd, stdin, stdout):
    with open(stdin, "rb") as i, open(stdout, "wb") as o:
        t0 = time.time()
        subprocess.run(pin + cmd, stdin=i, stdout=o, check=True)
        return time.time() - t0
with open(results, "a") as res:
    for name, comp, decomp in tools:
        dst = os.path.join(out_dir, f"sidecar.{name}")
        tc = statistics.median(timed(comp, src, dst) for _ in range(3))
        td = statistics.median(timed(decomp, dst, os.devnull)
                               for _ in range(3))
        size_in, size_out = os.path.getsize(src), os.path.getsize(dst)
        row = [n, name, 1, size_in, size_out, f"{size_in / size_out:.3f}",
               f"{tc:.4f}", "NA", f"{td:.4f}"] + ["NA"] * 5 + ["CLI, whole stream"]
        res.write("\t".join(str(x) for x in row) + "\n")
        print("   ", "\t".join(str(x) for x in row), flush=True)
        os.remove(dst)
EOF
}

[ -f "$COMP_RESULTS" ] || printf "n_variants\tcodec\tthreads\tin_bytes\tout_bytes\tratio\tcomp_1t_s\tcomp_nt_s\tdecomp_1t_s\tdecomp_nt_s\tdecomp_scan_1t_s\tlookup_us\tlookup64k_us\tn_blocks\tcheck\n" > "$COMP_RESULTS"
[ -f "$CHR_RESULTS" ] || printf "n_variants\tn_loci\tn_chromosomes\tloci_sidecar_bytes\tchromosome_sidecar_bytes\n" > "$CHR_RESULTS"

for n in $SIZES; do
  vcf="data/discosnp_${n}.vcf"
  [ -f "$vcf" ] || { echo "$vcf missing: run run_all.sh"; exit 1; }
  deep=false
  [[ " $DEEP_SIZES " == *" $n "* ]] && deep=true
  echo "== $n variants (deep: $deep), $(date)"

  # 1. reading with names_as_int = TRUE: baseline and sidecar, alternating
  # which goes first; then names_as_int = FALSE, for reference
  for rep in $(seq 1 "$REPS"); do
    if [ $((rep % 2)) -eq 1 ]; then order="base sidecar"; else order="sidecar base"; fi
    for v in $order; do
      run "$v" "$n" read_int "$rep" --warm-vcf
    done
    run sidecar "$n" read_chr "$rep" --warm-vcf
  done
  sidecar="work/sidecar_sidecar_${n}/gt.loci.tsv"
  ls -l "work/sidecar_base_${n}" "work/sidecar_sidecar_${n}"

  # 4. the planned chromosome sidecar (int_chromosome, chromosome: one line
  # per distinct chromosome of the biallelic SNPs, in order of first
  # appearance, as the parser's ChromosomeTable): its size, from the locus
  # sidecar (these chromosome names are not numbers, so their codes are
  # 1..K)
  awk -F'\t' -v n="$n" -v sz="$(stat -c %s "$sidecar")" '
    NR == 1 { bytes = length("int_chromosome\tchromosome\n"); next }
    { loci++ }
    !($3 in seen) { seen[$3] = 1; k++; bytes += length(k) + length($3) + 2 }
    END { printf "%s\t%d\t%d\t%s\t%d\n", n, loci, k, sz, bytes }
  ' "$sidecar" | tee -a "$CHR_RESULTS"

  if $deep; then
    # 3. lookups (first call builds the index; then warm), with the sidecar
    # in the page cache; then the first lookup with the sidecar evicted
    run sidecar "$n" lookup 1
    for rep in $(seq 1 "$DEEP_REPS"); do
      run sidecar "$n" lookup_cold "$rep" --evict "$sidecar"
    done
    for rep in $(seq 1 "$DEEP_REPS"); do
      run sidecar "$n" restore "$rep"
    done
    # where the time of the restore goes (to results/sidecar_restore_parts.tsv)
    for a in restore_parts restore_unchecked restore restore_unchecked restore; do
      RESULTS=results/sidecar_restore_parts.tsv run sidecar "$n" "$a" 1
    done
    for rep in $(seq 1 "$DEEP_REPS"); do
      if [ $((rep % 2)) -eq 1 ]; then order="base sidecar"; else order="sidecar base"; fi
      for v in $order; do
        run "$v" "$n" update "$rep"
      done
    done

    # 2. compressing the sidecar (on P-cores 0-7: 1 thread, then 4)
    comp_dir="work/sidecar_compress/${n}"
    mkdir -p "$comp_dir"
    python3 -c "import sys; sys.path.insert(0, '.'); import sidecar_bench as s; print('    busy cores:', round(s.wait_idle($MAX_BUSY, 900), 1))"
    taskset -c 0-7 "$comp_bin" "$sidecar" "$comp_dir" 4 3 2000 |
      tail -n +2 | sed "s/^/${n}\t/" | tee -a "$COMP_RESULTS"
    for f in "$comp_dir"/sidecar.*.gz; do
      if gzip -dc "$f" | cmp -s - "$sidecar"; then
        echo "    $f: gzip -dc gives the sidecar back"
      else
        echo "    $f: gzip -dc DIFFERS"
      fi
    done
    cli_compress "$n" "$sidecar" "$comp_dir"
    rm -rf "$comp_dir"
  fi

  # delete the outputs of this size
  rm -f work/sidecar_base_"${n}"/gt* work/sidecar_sidecar_"${n}"/gt*
  rm -rf work/sidecar_base_"${n}"/out work/sidecar_sidecar_"${n}"/out
done
echo "done $(date)"
