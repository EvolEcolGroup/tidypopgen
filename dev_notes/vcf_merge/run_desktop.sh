#!/bin/bash
# Unattended desktop (manica01) timings of the gt_merge_vcfs() prototypes
# (vcf_merge.md): base 5e207d78, P1 d3b982d0, P3 bc6856b6, P4 6982994c, and
# bcftools 1.25 and plink2, on the pops, sgdp and many datasets; then KING
# against PLINK 2 and the KING prototypes (pairwise_bitpacked.md): base
# b8255226, p1 72377615, p3 e2692136 (p3s: p3 without AVX2). Run
# setup_desktop.sh and the KING set-up ("KING against PLINK 2" in
# handover.md) first, then schedule this with a systemd user timer:
#
#   cp ~/git/tidypopgen/dev_notes/vcf_merge/run_desktop.sh ~/claude_work/tidypopgen/
#   systemd-run --user --on-calendar='YYYY-MM-DD 01:00' --unit=vcf_merge_king \
#     bash ~/claude_work/tidypopgen/run_desktop.sh
#   loginctl show-user "$USER" -p Linger    # must be yes to fire after logout
#
# Run a copy: bash reads a script as it runs, so a git pull in the checkout
# during the run could change it under the running job.
#
# From the benchmarking skill's template: waits for an idle machine, runs each
# phase once (markers in $WORK/done), starts no phase after NO_START_AFTER,
# logs to $WORK/status.md, and at the end commits and pushes the result tables
# to the checked-out branch and posts the log on claude_jobs#5.
set -uo pipefail
# run_timing.sh counts cycles with perf: needs kernel.perf_event_paranoid <= 2

# ---- CONFIG ------------------------------------------------------------------
JOB=${JOB:-vcf_merge_king}
REPO=${REPO:-$HOME/git/tidypopgen}
export TPG_WORK=${TPG_WORK:-$HOME/claude_work/tidypopgen}
W=$TPG_WORK/merge_vcf
WORK=${WORK:-$W/desktop}
RESULTS_IN_REPO=${RESULTS_IN_REPO:-dev_notes/vcf_merge/desktop_$(date +%Y-%m)}
BRANCH=${BRANCH:-pca_missing}   # results are committed only if REPO is on it
export TPG_REPO=$REPO
K=$TPG_WORK/plink_king
KING_IN_REPO=${KING_IN_REPO:-dev_notes/pairwise_bitpacked/king_desktop_$(date +%F)}
KROUNDS=${KROUNDS:-"F R F R F R"}  # each round runs every task once; R reverses
ISSUE=${ISSUE:-dramanica/claude_jobs/issues/5}
IDLE_LOAD=${IDLE_LOAD:-1.0}
NO_START_AFTER=${NO_START_AFTER:-$(date -d 'tomorrow 07:00' '+%F %H:%M')}
export OMP_PROC_BIND=close OMP_PLACES=cores
ATTRIB=${ATTRIB:-"Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"}
V="base:$W/libs/base p1:$W/libs/p1 p3:$W/libs/p3 p4:$W/libs/p4"
# ------------------------------------------------------------------------------

STATUS=$WORK/status.md; DONE=$WORK/done; LOGS=$WORK/logs
mkdir -p "$DONE" "$LOGS"

say() { echo "- $(date '+%a %d %b %H:%M'): $*" | tee -a "$STATUS"; }
load() { cut -d' ' -f1-3 /proc/loadavg; }
lt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a < b) }'; }
too_late() { [ "$(date +%s)" -ge "$(date -d "$NO_START_AFTER" +%s)" ]; }

wait_idle() {
  local ok=0
  while [ $ok -lt 3 ]; do
    too_late && return 1
    if lt "$(cut -d' ' -f2 /proc/loadavg)" "$IDLE_LOAD"; then ok=$((ok + 1)); else ok=0; fi
    [ $ok -lt 3 ] && sleep 300
  done
  return 0
}
settle() {
  local n=0
  while ! lt "$(cut -d' ' -f1 /proc/loadavg)" 1.0; do
    sleep 60; n=$((n + 1))
    if [ $n -ge 30 ]; then
      say "load still $(load) after 30 min: waiting for idle"
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
  if (cd "$W" && "$@") >> "$LOGS/phase_$name.log" 2>&1; then
    touch "$DONE/$name"; say "done $name ($(( (SECONDS - t0) / 60 )) min)"
  else
    say "FAILED $name ($(( (SECONDS - t0) / 60 )) min), see logs/phase_$name.log"
  fi
}

# the merge timings of one scenario list at one thread count, pinned
vcf() {  # vcf <threads> <cores>
  local t=$1 c=$2 sc
  for sc in pops sgdp many; do
    CORES=$c MERGE_MODE=real scripts/run_timing.sh "$sc" "$t" 3 $V || return 1
  done
}
tools() {  # tools <threads> <cores> <scenario>...
  local t=$1 c=$2 sc r mem
  shift 2
  for sc in "$@"; do
    mem=4000  # plink2's workspace (MiB); the sgdp merge runs out of memory with it
    [ "$sc" = sgdp ] && mem=32000
    for r in 1 2 3; do
      CORES=$c PLINK_MEM=$mem scripts/tools_timing.sh "$sc" "$t" "$r" || return 1
    done
  done
}

# KING: run_bench.sh (PLINK 2 against pairwise_king() of b8255226) and
# run_proto.sh (the prototypes); each fails if rows are missing from its tables
king_bench() {  # king_bench <threads>
  local t=$1 n=$(wc -w <<< "$KROUNDS")
  bash "$K/scripts/run_bench.sh" "$t" $KROUNDS
  [ "$(wc -l < "$K/out/bench_plink_$t.tsv")" -ge $((8 * n)) ] &&
    [ "$(wc -l < "$K/out/bench_tpg_$t.tsv")" -ge $((4 * n)) ] &&
    ! grep -q "exited with non-zero" "$K/out/bench_plink_$t.tsv"
}
king_proto() {  # king_proto <threads>
  local t=$1 n=$(wc -w <<< "$KROUNDS")
  VERS="base p1 p3 p3s" bash "$K/scripts/run_proto.sh" "$t" "full grouped" $KROUNDS
  [ "$(wc -l < "$K/out/proto_$t.tsv")" -ge $((8 * n)) ]
}

report() {
  trap - EXIT
  say "reporting"
  local dest=$REPO/$RESULTS_IN_REPO kdest=$REPO/$KING_IN_REPO pushed="" br
  mkdir -p "$dest" "$kdest"
  cp "$W"/results/runs.tsv "$W"/results/tools.tsv "$W"/results/phases_*.tsv "$dest/" 2>/dev/null
  cp "$STATUS" "$dest/" 2>/dev/null
  cp "$K"/out/bench_*.tsv "$K"/out/proto_*.tsv "$K"/out/*_log.txt "$kdest/" 2>/dev/null
  Rscript "$K/scripts/summarise.R" "$K/out" > "$kdest/summary.txt" 2>&1
  br=$(cd "$REPO" && git rev-parse --abbrev-ref HEAD)
  if [ "$br" != "$BRANCH" ]; then
    echo "$REPO is on $br, not $BRANCH: nothing committed" >> "$LOGS/report.log"
  elif (cd "$REPO" && git pull -q --ff-only &&
      git add -A "$RESULTS_IN_REPO" "$KING_IN_REPO" &&
      git commit -qm "Desktop run of $(date +%F): $JOB results

$ATTRIB" -- "$RESULTS_IN_REPO" "$KING_IN_REPO" && git push -q origin HEAD) >> "$LOGS/report.log" 2>&1; then
    pushed="$br@$(cd "$REPO" && git rev-parse --short HEAD)"
  fi
  {
    echo "## $(date +%F): unattended desktop run of $JOB finished"
    echo
    echo "Pushed: ${pushed:-nothing (see $LOGS/report.log on manica01)}; results in \`$RESULTS_IN_REPO\` and \`$KING_IN_REPO\`. Not interpreted yet."
    echo
    cat "$STATUS"
  } > "$WORK/comment.md"
  gh api "repos/$ISSUE/comments" -F body=@"$WORK/comment.md" >> "$LOGS/report.log" 2>&1 ||
    say "could not post the comment; it is in $WORK/comment.md"
}

trap report EXIT
say "scheduled run started; waiting for the machine to be idle (5-min load < $IDLE_LOAD)"
if ! wait_idle; then say "the machine was never idle before $NO_START_AFTER: nothing run"; exit 0; fi
say "machine idle (load $(load))"

# ---- PHASES, in order of priority ------------------------------------------------
# (1 and 8 threads, and the tools at 8; the laptop has 4-thread runs of 20 Nov)
phase vcf_t1   vcf 1 2                  # ~30-40 min
phase vcf_t8   vcf 8 0-7                # ~20 min
phase tools_t8 tools 8 0-7 pops many    # ~12 min
phase tools_t8_sgdp tools 8 0-7 sgdp    # ~12 min (bcftools ~2.5 min a run)
phase king_bench_t1 king_bench 1        # core 2; ~15 min
phase king_bench_t8 king_bench 8        # P-cores; ~10 min
phase king_proto_t1 king_proto 1        # core 2; ~15 min
phase king_proto_t8 king_proto 8        # P-cores; ~8 min
# ------------------------------------------------------------------------------
say "all phases run"
