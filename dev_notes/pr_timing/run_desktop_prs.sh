#!/bin/bash
# Unattended desktop (manica01) timings of PR #239 (king_avx2: KING AVX2
# kernel and tiles) and PR #240 (vcf_merge_speed: faster gt_merge_vcfs()),
# against dev, with -O2 installs of each in its own library:
#   - KING-robust on HGDP650 (full matrix and within populations): dev,
#     king_avx2, and king_avx2 with TIDYPOPGEN_NO_AVX2=1 (scalar kernel);
#   - gt_merge_vcfs() on the pops, sgdp and many datasets: dev, vcf_merge_speed;
#   - gen_tibble() on one VCF with 100 samples and one with 1: dev,
#     vcf_merge_speed (its parallel parsing is in the single-file reader).
# 1 thread on core 2, 8 threads on the P-cores (0-7); builds and checks on
# the E-cores (8-23).
#
# It waits until the timings of the tidypopgl session are over (the marker
# file ~/tpgl_prof/TIMING_BUSY is gone) and the machine is idle; nothing,
# not even the builds, runs before. Reuses the data and scripts of the
# 20 Nov job (~/claude_work/tidypopgen/merge_vcf and plink_king, through
# symlinks in $WORK/tpg_work, so that the old results are not touched), and
# its own clone of the repository ($WORK/repo, never ~/git/tidypopgen).
#
#   mkdir -p ~/claude_work/tidypopgen/pr_timing/scripts
#   cp dev_notes/pr_timing/* ~/claude_work/tidypopgen/pr_timing/scripts/
#   systemd-run --user --on-active=60 --unit=tpg_pr_timing \
#     bash ~/claude_work/tidypopgen/pr_timing/scripts/run_desktop_prs.sh
#
# Phases leave markers in $WORK/done (a rerun continues where it stopped);
# no phase starts after NO_START_AFTER; the log is $WORK/status.md; at the
# end (also after a failure) the tables and summary.md are committed to the
# branch split_plan (dev_notes/pr_timing/desktop_<date>/) and the log is
# posted on claude_jobs#5.
set -uo pipefail

# ---- CONFIG ------------------------------------------------------------------
JOB=${JOB:-tpg_pr_timing}
OLD=${OLD:-$HOME/claude_work/tidypopgen}   # the work folder of the 20 Nov job
WORK=${WORK:-$OLD/pr_timing}
S=$WORK/scripts                            # this script and its R helpers
REPO=$WORK/repo
REMOTE=${REMOTE:-git@github.com:EvolEcolGroup/tidypopgen.git}
BRANCH=${BRANCH:-split_plan}               # where the results are committed
RESULTS_IN_REPO=${RESULTS_IN_REPO:-dev_notes/pr_timing/desktop_$(date +%F)}
# the versions: label, git ref
REFS=${REFS:-"dev:origin/dev king:origin/king_avx2 vcf:origin/vcf_merge_speed"}
BUSY_MARKER=${BUSY_MARKER:-$HOME/tpgl_prof/TIMING_BUSY}
ROUNDS=${ROUNDS:-"F R F R F R"}            # balanced order of the KING runs
ISSUE=${ISSUE:-dramanica/claude_jobs/issues/5}
IDLE_LOAD=${IDLE_LOAD:-1.0}
NO_START_AFTER=${NO_START_AFTER:-2026-11-23 07:00}
ATTRIB=${ATTRIB:-"Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"}
export OMP_PROC_BIND=close OMP_PLACES=cores
export TPG_WORK=$WORK/tpg_work            # read by the scripts of 20 Nov
# (only read for vignettes/articles/data/hgdp/hgdp650_id_pop.txt)
export TPG_REPO=${TPG_REPO:-$HOME/git/tidypopgen}
export TMPDIR=$WORK/rtmp
# ------------------------------------------------------------------------------

STATUS=$WORK/status.md; DONE=$WORK/done; LOGS=$WORK/logs; OUT=$WORK/out
mkdir -p "$DONE" "$LOGS" "$OUT" "$TMPDIR"

say() { echo "- $(date '+%a %d %b %H:%M'): $*" | tee -a "$STATUS"; }
load() { cut -d' ' -f1-3 /proc/loadavg; }
lt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a < b) }'; }
too_late() { [ "$(date +%s)" -ge "$(date -d "$NO_START_AFTER" +%s)" ]; }
lib() { echo "$WORK/libs/$1"; }

# idle: the tidypopgl timings are over and the 5-min load is below IDLE_LOAD
# at three checks five minutes apart
wait_idle() {
  local ok=0
  while [ $ok -lt 3 ]; do
    too_late && return 1
    if [ -e "$BUSY_MARKER" ]; then
      ok=0
    elif lt "$(cut -d' ' -f2 /proc/loadavg)" "$IDLE_LOAD"; then
      ok=$((ok + 1))
    else
      ok=0
    fi
    [ $ok -lt 3 ] && sleep 300
  done
  return 0
}
# before each phase: the marker is still absent and the 1-min load below 1
settle() {
  local n=0
  while [ -e "$BUSY_MARKER" ] || ! lt "$(cut -d' ' -f1 /proc/loadavg)" 1.0; do
    sleep 60; n=$((n + 1))
    if [ $n -ge 30 ]; then
      say "busy (load $(load), marker $([ -e "$BUSY_MARKER" ] && echo present || echo absent)) after 30 min: waiting for idle"
      wait_idle || return 1; n=0
    fi
  done
}
phase() {
  local name=$1; shift
  [ -e "$DONE/$name" ] && return 0
  if too_late; then say "not started: $name (after $NO_START_AFTER)"; return 1; fi
  settle || { say "not started: $name (machine busy)"; return 1; }
  say "start $name (load $(load))"
  local t0=$SECONDS
  # (in a subshell, so that set -e in a phase does not reach this script)
  if ("$@") >> "$LOGS/phase_$name.log" 2>&1; then
    touch "$DONE/$name"; say "done $name ($(( (SECONDS - t0) / 60 )) min)"
  else
    say "FAILED $name ($(( (SECONDS - t0) / 60 )) min), see logs/phase_$name.log"
  fi
}

# ---- set-up and checks (untimed, E-cores) ------------------------------------
setup() {
  set -e
  if [ ! -d "$REPO/.git" ]; then git clone -q "$REMOTE" "$REPO"; fi
  git -C "$REPO" fetch -q origin
  git -C "$REPO" checkout -q "$BRANCH"
  git -C "$REPO" reset -q --hard "origin/$BRANCH"
  local pair v ref commit
  for pair in $REFS; do
    v=${pair%%:*}; ref=${pair#*:}
    commit=$(git -C "$REPO" rev-parse --short=8 "$ref")
    rm -rf "$WORK/src/$v" "$(lib "$v")"
    mkdir -p "$WORK/src/$v" "$(lib "$v")"
    git -C "$REPO" archive "$ref" | tar -x -C "$WORK/src/$v"
    echo "Install $v ($ref, $commit)"
    MAKEFLAGS=-j16 taskset -c 8-23 R CMD INSTALL --no-test-load \
      --library="$(lib "$v")" "$WORK/src/$v"
    echo "$commit" > "$(lib "$v")/COMMIT"
    # an -O2 build: no out-of-line Rcpp indexing (which only -O0 leaves)
    if [ "$(objdump -d "$(lib "$v")/tidypopgen/libs/tidypopgen.so" |
      grep -c 'call.*vectorI.*ixEm')" -gt 0 ]; then
      echo "$v: -O0 code in the build"; return 1
    fi
  done
  # the folders the scripts of 20 Nov expect, with their data and scripts
  local T=$TPG_WORK
  mkdir -p "$T/merge_vcf/results/meta" "$T/merge_vcf/out" "$T/merge_vcf/rtmp" \
    "$T/plink_king/build" "$T/plink_king/out" "$T/plink_king/rtmp"
  ln -sfn "$OLD/merge_vcf/data" "$T/merge_vcf/data"
  ln -sfn "$OLD/merge_vcf/scripts" "$T/merge_vcf/scripts"
  ln -sfn "$OLD/plink_king/scripts" "$T/plink_king/scripts"
  ln -sfn "$OLD/plink_king/data" "$T/plink_king/data"
  ln -sfn "$OLD/bench_inb" "$T/bench_inb"
  # run_proto.sh's labels: base = dev, p3 = king_avx2 (p3s: without AVX2)
  ln -sfn "$(lib dev)" "$T/plink_king/build/lib_base"
  ln -sfn "$(lib king)" "$T/plink_king/build/lib_p3"
}

# the KING outputs of king_avx2, with and without AVX2, must be identical to
# dev's (check_proto.R of 20 Nov: KING full, grouped and subset, IBS, allele
# sharing, grouped inbreeding)
check() {
  set -e
  local P=$TPG_WORK/plink_king/scripts
  taskset -c 8-23 Rscript "$P/check_proto.R" "$(lib dev)" "$OUT/check_dev.rds"
  taskset -c 8-23 Rscript "$P/check_proto.R" "$(lib king)" "$OUT/check_king.rds"
  TIDYPOPGEN_NO_AVX2=1 taskset -c 8-23 Rscript "$P/check_proto.R" \
    "$(lib king)" "$OUT/check_king_scalar.rds"
  Rscript -e 'a <- readRDS(commandArgs(TRUE)[1])
    for (f in commandArgs(TRUE)[-1]) {
      b <- readRDS(f)
      ok <- vapply(names(a), function(n) identical(a[[n]], b[[n]]), TRUE)
      cat(basename(f), if (all(ok)) "identical to dev" else
        paste("DIFFERS from dev:", paste(names(a)[!ok], collapse = ", ")), "\n")
    }' "$OUT/check_dev.rds" "$OUT/check_king.rds" "$OUT/check_king_scalar.rds" |
    tee "$OUT/checks.txt"
}

# ---- timed phases (P-cores) ----------------------------------------------------
# KING: run_proto.sh (one process per run, the faster of two calls is used)
king() {  # king <threads>
  local t=$1 n
  n=$(wc -w <<< "$ROUNDS")
  VERS="base p3 p3s" bash "$TPG_WORK/plink_king/scripts/run_proto.sh" "$t" \
    "full grouped" $ROUNDS
  [ "$(wc -l < "$TPG_WORK/plink_king/out/proto_$t.tsv")" -ge $((6 * n)) ]
}
# gt_merge_vcfs(): run_timing.sh (3 runs per version, A B B A A B)
vcf() {  # vcf <threads> <cores>
  local t=$1 c=$2 sc
  (cd "$TPG_WORK/merge_vcf" &&
    for sc in pops sgdp many; do
      CORES=$c MERGE_MODE=real scripts/run_timing.sh "$sc" "$t" 3 \
        "dev:$(lib dev)" "vcf:$(lib vcf)" || exit 1
    done)
}
# gen_tibble() on one VCF of each kind, in a balanced order
gt() {  # gt <threads> <cores>
  local t=$1 c=$2 r f v list
  for r in $ROUNDS; do
    if [ "$r" = F ]; then list="dev vcf"; else list="vcf dev"; fi
    for f in "$OLD/merge_vcf/data/pops/pop1.vcf.gz" \
      "$OLD/merge_vcf/data/sgdp/S001.vcf.gz"; do
      for v in $list; do
        OMP_NUM_THREADS=$t taskset -c "$c" Rscript "$S/time_gen_tibble.R" \
          "$(lib "$v")" "$f" "$t" "$v" "$OUT/gen_tibble.tsv" || return 1
        sleep 3
      done
    done
  done
}

report() {
  trap - EXIT
  say "reporting"
  local dest=$REPO/$RESULTS_IN_REPO pushed=""
  mkdir -p "$dest"
  cp "$TPG_WORK"/merge_vcf/results/runs.tsv "$TPG_WORK"/merge_vcf/results/phases_*.tsv \
    "$TPG_WORK"/plink_king/out/proto_*.tsv "$TPG_WORK"/plink_king/out/proto_log.txt \
    "$OUT"/gen_tibble.tsv "$OUT"/checks.txt "$dest/" 2>/dev/null
  Rscript "$S/summarise_prs.R" "$dest" > "$dest/summary.md" 2>> "$LOGS/report.log"
  cp "$STATUS" "$dest/"
  if (cd "$REPO" && git checkout -q "$BRANCH" && git pull -q --ff-only &&
      git add -A "$RESULTS_IN_REPO" &&
      git commit -qm "[skip ci] Desktop run of $(date +%F): $JOB results

$ATTRIB" -- "$RESULTS_IN_REPO" && git push -q origin HEAD) >> "$LOGS/report.log" 2>&1; then
    pushed="$BRANCH@$(git -C "$REPO" rev-parse --short HEAD)"
  fi
  {
    echo "## $(date +%F): unattended desktop run of $JOB (PRs #239 and #240) finished"
    echo
    echo "Pushed: ${pushed:-nothing (see $LOGS/report.log on manica01)}; results in \`$RESULTS_IN_REPO\`. Not interpreted yet."
    echo
    cat "$dest/summary.md" 2>/dev/null
    echo
    cat "$STATUS"
  } > "$WORK/comment.md"
  gh api "repos/$ISSUE/comments" -F body=@"$WORK/comment.md" >> "$LOGS/report.log" 2>&1 ||
    say "could not post the comment; it is in $WORK/comment.md"
}

trap report EXIT
say "started; waiting for the tidypopgl timings to end ($BUSY_MARKER) and an idle machine (5-min load < $IDLE_LOAD)"
if ! wait_idle; then say "never idle before $NO_START_AFTER: nothing run"; exit 0; fi
say "machine idle (load $(load))"

# ---- PHASES, in order of priority ------------------------------------------------
phase setup setup                # E-cores; ~10 min (three installs)
[ -e "$DONE/setup" ] || { say "set-up failed: nothing timed"; exit 1; }
phase check check                # E-cores; ~5 min
phase king_t1 king 1             # core 2; ~10 min
phase vcf_t1 vcf 1 2             # core 2; ~20 min
phase king_t8 king 8             # P-cores; ~6 min
phase vcf_t8 vcf 8 0-7           # P-cores; ~10 min
phase gt_t8 gt 8 0-7             # P-cores; ~5 min
phase gt_t1 gt 1 2               # core 2; ~8 min
# ------------------------------------------------------------------------------
say "all phases run"
