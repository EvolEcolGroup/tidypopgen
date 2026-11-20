#!/usr/bin/env bash
# Time bcftools merge and plink2 (import + pmerge) on a scenario.
#   tools_timing.sh <scenario> <n_threads> <rep> [steps]
# steps: a regular expression over bcftools_ou, bcftools_oz, plink_import,
# plink_merge (default: all). Appends to results/tools.tsv (wall time, peak
# RSS, user-space cycles and instructions, load, date).
set -euo pipefail
W=${TPG_WORK:-/media/andrea/Elements1/tmp/claude/tidypopgen}/merge_vcf
B=$W/tools/bcftools/bin/bcftools
P=${TPG_WORK:-/media/andrea/Elements1/tmp/claude/tidypopgen}/plink_king/bin/plink2
SC=$1
T=$2
REP=$3
STEPS=${4:-.}
want() { [[ "$1" =~ $STEPS ]]; }
CORES=${CORES:-"0-$((T - 1))"}
# a counter from perf stat -x, summed over the core types of a hybrid CPU
# (the desktop prints one line for cpu_atom and one for cpu_core)
perf_sum() { awk -F, -v e="$1" '$3 ~ e && $1 ~ /^[0-9]+$/ { s += $1; n++ } END { if (n) print s; else print "NA" }' "$2"; }
# count cycles only if perf can (kernel.perf_event_paranoid <= 2; a reboot
# resets it): otherwise time without perf, and record NA for the counters
if perf stat -e cycles:u true > /dev/null 2>&1; then PERF_OK=1; else PERF_OK=0
  echo "perf cannot count cycles: counters recorded as NA" >&2; fi
counted() {  # counted <perf output file> <command...>
  local o=$1; shift
  if [ "$PERF_OK" = 1 ]; then perf stat -x, -e cycles:u,instructions:u -o "$o" "$@"
  else : > "$o"; "$@"; fi
}
RES=$W/results/tools.tsv
[ -f "$RES" ] || echo -e "date\ttool\tversion\tscenario\tthreads\trep\tstep\twall_s\tpeak_rss_mb\tload1\tcycles_u\tinstr_u\tout_mb" > "$RES"
OUT=$W/out/tools/$SC
mkdir -p "$OUT"
FILES=("$W/data/$SC"/*.vcf.gz)
BV=$($B --version | head -1 | cut -d' ' -f2)
PV=$($P --version | cut -d' ' -f2)
# bcftools merge options matching gt_merge_vcfs() (as far as possible):
# separate records for different alleles (-m none), and missing as REF for
# the per-sample files (-0)
if [ "$SC" == "pops" ]; then MOPT=(-m none); else MOPT=(-m none -0); fi
timed() { # tool version step outfile command...
  local tool=$1 ver=$2 step=$3 outf=$4
  shift 4
  local load1
  load1=$(cut -d' ' -f1 /proc/loadavg)
  counted "$OUT/perf" \
    /usr/bin/time -f "%e %M" -o "$OUT/time" "$@" > "$OUT/$step.log" 2>&1
  local wall rss cyc ins size=NA
  read -r wall rss < "$OUT/time"
  cyc=$(perf_sum cycles:u "$OUT/perf")
  ins=$(perf_sum instructions:u "$OUT/perf")
  if [ -n "$outf" ]; then size=$(du -cm $outf | tail -1 | cut -f1); fi
  echo -e "$(date +%FT%T)\t$tool\t$ver\t$SC\t$T\t$REP\t$step\t$wall\t$((rss / 1024))\t$load1\t$cyc\t$ins\t$size" >> "$RES"
  echo "$tool $step $SC T=$T rep=$REP: ${wall}s $((rss / 1024))MB load=$load1"
  sleep 5
}
if want bcftools_ou; then
  timed bcftools "$BV" bcftools_ou "" taskset -c "$CORES" \
    "$B" merge "${MOPT[@]}" --threads "$T" -Ou -o /dev/null "${FILES[@]}"
fi
if want bcftools_oz; then
  rm -f "$OUT/merged.vcf.gz"
  timed bcftools "$BV" bcftools_oz "$OUT/merged.vcf.gz" taskset -c "$CORES" \
    "$B" merge "${MOPT[@]}" --threads "$T" -Oz -o "$OUT/merged.vcf.gz" "${FILES[@]}"
fi
if want plink_import || want plink_merge; then
  PD=$OUT/plink
  if want plink_import; then
    rm -rf "${PD:?}"
  fi
  mkdir -p "$PD"
  # one import per file (as plink2 reads one VCF at a time), with IDs made of
  # the chromosome, position and alleles, so that loci are matched by
  # position and alleles
  # (plink2 refuses single-sample VCFs without homozygous REF genotypes,
  # such as those of the sgdp and many scenarios, unless told otherwise)
  NONVAR=""
  if [ "$SC" != "pops" ]; then NONVAR="--vcf-allow-no-nonvar"; fi
  cat > "$PD/import.sh" <<EOF
#!/usr/bin/env bash
set -e
for f in ${FILES[@]}; do
  b=\$(basename "\$f" .vcf.gz)
  "$P" --vcf "\$f" --make-pgen --out "$PD/\$b" --threads $T --memory 4000 \
    --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 10 $NONVAR --silent
done
EOF
  chmod +x "$PD/import.sh"
  if want plink_import; then
    timed plink2 "$PV" plink_import "$PD/*.pgen $PD/*.pvar" taskset -c "$CORES" "$PD/import.sh"
  fi
  if want plink_merge; then
    ls "$PD"/*.pgen | sed 's/\.pgen$//' > "$PD/list.txt"
    FILL=()
    if [ "$SC" != "pops" ]; then FILL=(--fill-missing-with-ref); fi
    timed plink2 "$PV" plink_merge "$PD/merged.pgen $PD/merged.pvar" taskset -c "$CORES" \
      "$P" --pmerge-list "$PD/list.txt" --make-pgen --out "$PD/merged" \
      --threads "$T" --memory 4000 --indiv-sort none "${FILL[@]}" --silent
  fi
fi
exit 0
