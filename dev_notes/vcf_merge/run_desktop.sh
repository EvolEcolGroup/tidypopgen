#!/bin/bash
# Unattended desktop (manica01) timings of the gt_merge_vcfs() prototypes
# (vcf_merge.md): base 5e207d78, P1 d3b982d0, P3 bc6856b6, P4 6982994c, and
# bcftools 1.25 and plink2, on the pops, sgdp and many datasets. Run
# setup_desktop.sh first, then schedule this with a systemd user timer:
#
#   systemd-run --user --on-calendar='YYYY-MM-DD 01:00' --unit=vcf_merge \
#     bash ~/git/tidypopgen/dev_notes/vcf_merge/run_desktop.sh
#   loginctl show-user "$USER" -p Linger    # must be yes to fire after logout
#
# From the benchmarking skill's template: waits for an idle machine, runs each
# phase once (markers in $WORK/done), starts no phase after NO_START_AFTER,
# logs to $WORK/status.md, and at the end commits and pushes the result tables
# to the checked-out branch and posts the log on claude_jobs#5.
set -uo pipefail
# run_timing.sh counts cycles with perf: needs kernel.perf_event_paranoid <= 2

# ---- CONFIG ------------------------------------------------------------------
JOB=${JOB:-vcf_merge}
REPO=${REPO:-$HOME/git/tidypopgen}
export TPG_WORK=${TPG_WORK:-$HOME/claude_work/tidypopgen}
W=$TPG_WORK/merge_vcf
WORK=${WORK:-$W/desktop}
RESULTS_IN_REPO=${RESULTS_IN_REPO:-dev_notes/vcf_merge/desktop_$(date +%Y-%m)}
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
tools() {  # tools <threads> <cores>
  local t=$1 c=$2 sc r
  for sc in pops many sgdp; do
    for r in 1 2 3; do CORES=$c scripts/tools_timing.sh "$sc" "$t" "$r" || return 1; done
  done
}

report() {
  trap - EXIT
  say "reporting"
  local dest=$REPO/$RESULTS_IN_REPO pushed=""
  mkdir -p "$dest"
  cp "$W"/results/runs.tsv "$W"/results/tools.tsv "$W"/results/phases_*.tsv "$dest/" 2>/dev/null
  cp "$STATUS" "$dest/" 2>/dev/null
  if (cd "$REPO" && git add -A "$RESULTS_IN_REPO" &&
      git commit -qm "Desktop run of $(date +%F): $JOB results

$ATTRIB" && git push -q origin HEAD) >> "$LOGS/report.log" 2>&1; then
    pushed="$(cd "$REPO" && git rev-parse --abbrev-ref HEAD)@$(cd "$REPO" && git rev-parse --short HEAD)"
  fi
  {
    echo "## $(date +%F): unattended desktop run of $JOB finished"
    echo
    echo "Pushed: ${pushed:-nothing (see $LOGS/report.log on manica01)}; results in \`$RESULTS_IN_REPO\`. Not interpreted yet."
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
# (the laptop has 4-thread runs of 20 Nov; 1 thread and the tools are missing)
phase vcf_t1   vcf 1 2                  # ~30-40 min
phase vcf_t8   vcf 8 0-7                # ~20 min
phase tools_t8 tools 8 0-7              # ~1 h (bcftools on sgdp ~4 min a run)
# ------------------------------------------------------------------------------
say "all phases run"
