#!/usr/bin/env bash
# Full run of the VCFArrow benchmark, for when the machine is free: 1.0.0
# against the branch vcf_fast (install both first with ./setup_libs.sh).
#   ./run_all.sh                      # default plan, several hours
#   SIZES="2322460" ./run_all.sh      # one size only
# Settings (environment variables):
#   SIZES       numbers of variants: 1/8, 1/4 and all of the 18,579,683 of
#               the published dataset
#   OLD, NEW    labels of the libraries in lib/ (see setup_libs.sh)
#   OLD_MAX     largest size run with OLD (all by default; 0.4.4 needs about
#               39 GB at full size, as in the published run)
#   THREADS     threads given to tidypopgen (the comparison used 4 and 15)
#   MEMORY_MAX  memory limit for each step (the comparison used 50 GB)
#   TIMEOUT     time limit for each step, in s (the comparison used 60 min
#               at full size); "scaled" (default) gives 3600 s at full size
#               and proportionally less on smaller sizes, but at least 600 s
#   STEPS       all (17 functions of the comparison + extras), or core
#   OLD_STEPS   steps run with OLD (STEPS by default), e.g. all steps but
#               write_vcf, which is quadratic in the number of chromosomes in
#               1.0.0 and so only times out on this data (> 2 min already at
#               200,000 variants)
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs results
SIZES=${SIZES:-"2322460 4644921 18579683"}
OLD=${OLD:-v1.0.0}
NEW=${NEW:-vcf_fast}
OLD_MAX=${OLD_MAX:-18579683}
THREADS=${THREADS:-4}
MEMORY_MAX=${MEMORY_MAX:-12G}
TIMEOUT=${TIMEOUT:-scaled}
STEPS=${STEPS:-all}
OLD_STEPS=${OLD_STEPS:-$STEPS}
RESULTS=${RESULTS:-results/vcfarrow_bench.tsv}
log="logs/run_all_$(date +%Y%m%d_%H%M).log"
exec > >(tee -a "$log") 2>&1

echo "run_all.sh $(date): SIZES=$SIZES OLD=$OLD NEW=$NEW OLD_MAX=$OLD_MAX" \
  "THREADS=$THREADS MEMORY_MAX=$MEMORY_MAX TIMEOUT=$TIMEOUT STEPS=$STEPS" \
  "OLD_STEPS=$OLD_STEPS"
uptime
free -g
for v in "$OLD" "$NEW"; do
  [ -d "lib/$v/tidypopgen" ] || { echo "lib/$v missing: run setup_libs.sh"; exit 1; }
  echo "$v: commit $(cat "lib/$v/COMMIT")"
done

# the subsample of the real data shipped with VCFArrow, which make_vcf.R
# resamples (downloaded, not committed: it is VCFArrow's data)
mkdir -p data
for f in vaillantii_discosnp_sub.vcf.gz strata; do
  if [ ! -f "data/$f" ]; then
    curl -fsSL -o "data/$f" \
      "https://raw.githubusercontent.com/legalLab/VCFArrow/main/inst/extdata/$f"
  fi
done

for n in $SIZES; do
  vcf="data/discosnp_${n}.vcf"
  if [ ! -f "$vcf" ] || [ ! -f "data/discosnp_${n}.n_variants" ]; then
    echo "== generating $vcf"
    Rscript make_vcf.R "$n" "$vcf" 1
  fi
  for v in "$OLD" "$NEW"; do
    if [ "$v" = "$OLD" ] && [ "$n" -gt "$OLD_MAX" ]; then
      echo "== skipping $v at $n variants (OLD_MAX=$OLD_MAX)"
      continue
    fi
    if [ "$TIMEOUT" = "scaled" ]; then
      timeout_s=$(( 3600 * n / 18579683 ))
      [ "$timeout_s" -lt 600 ] && timeout_s=600
    else
      timeout_s=$TIMEOUT
    fi
    steps=$STEPS
    [ "$v" = "$OLD" ] && steps=$OLD_STEPS
    echo "== $v, $n variants, time limit ${timeout_s} s per step"
    python3 bench.py --lib "lib/$v" --label "$v" --vcf "$vcf" \
      --work "work/${v}_${n}_t${THREADS}" --results "$RESULTS" \
      --threads "$THREADS" --memory-max "$MEMORY_MAX" --timeout "$timeout_s" \
      --steps "$steps"
    Rscript summarise.R "$RESULTS" > results/summary.md
  done
done
echo "done $(date); summary in results/summary.md"
