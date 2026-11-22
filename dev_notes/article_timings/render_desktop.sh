#!/bin/bash
# Re-render the timed articles (benchmark_hgdp, big_data, rad_seq, roh_hgdp)
# with an -O2 build on a quiet manica01 (used on 22 Nov 2026 for dev at
# e812c8d3). Set-up, untimed, before the machine is quiet:
#   git clone --branch dev https://github.com/EvolEcolGroup/tidypopgen.git $W/src
#   ln -s ~/git/tidypopgen/vignettes/articles/data $W/src/vignettes/articles/data
#   git -C $W/src archive HEAD | tar -x -C $W/pkg
#   R CMD INSTALL --library=$W/lib $W/pkg   (check: no -O0 code in the .so)
#   git -C $W/src show HEAD:vignettes/articles/benchmark_hgdp.Rmd |
#     sed -n '/^# Running the benchmark on a laptop/,$p' > $W/laptop_section.md
# Run it as a systemd user unit named tpg_pr_timing*, e.g.
#   systemd-run --user --unit=tpg_pr_timing_articles --setenv=W=$W \
#     bash dev_notes/article_timings/render_desktop.sh
# Two passes: the
# first warms the caches and measures the run-to-run noise; the second is the
# one committed. Every knit starts only when the 1-minute load is below 1 at
# three checks 30 s apart. Runs as the systemd user unit
# tpg_pr_timing_articles (the tidypopgl queue holds its phases while a
# tpg_pr_timing* unit is active) and holds ~/tpgl_prof/TIMING_BUSY.
set -u
W=${W:-$HOME/claude_work/tidypopgen/articles_$(date +%F)}
LIB=$W/lib
ART=$W/src/vignettes/articles
MARK=${MARK:-$HOME/tpgl_prof/TIMING_BUSY}
TAG="tidypopgen article timing $$"
exec >> "$W/render.out" 2>&1

# remove the marker only if this job wrote it
cleanup() {
  if grep -qF "$TAG" "$MARK" 2>/dev/null; then rm -f "$MARK"; fi
  echo "finished $(date -Is); marker: $(cat "$MARK" 2>/dev/null || echo none)"
}
trap cleanup EXIT
[ -e "$MARK" ] || echo "$TAG, $(date -Is)" > "$MARK"

cd "$ART" || exit 1
export PLINK=/usr/bin/plink1.9
unset OMP_NUM_THREADS OMP_THREAD_LIMIT OPENBLAS_NUM_THREADS

# wait until the 1-minute load is below 1 at three checks 30 s apart
# (give up after 4 hours)
wait_quiet() {
  local quiet=0 start load
  start=$(date +%s)
  while [ $quiet -lt 3 ]; do
    if [ $(($(date +%s) - start)) -gt 14400 ]; then
      echo "gave up waiting for a quiet machine $(date -Is)"
      exit 2
    fi
    load=$(cut -d' ' -f1 /proc/loadavg)
    if awk "BEGIN{exit !($load < 1.0)}"; then quiet=$((quiet + 1)); else quiet=0; fi
    sleep 30
  done
}

declare -A out=([hgdp_benchmark]=benchmark_hgdp [big_data]=big_data
  [rad_seq]=rad_seq [roh_hgdp]=roh_hgdp)
echo "started $(date -Is) host $(hostname) commit $(git rev-parse --short HEAD)"
for pass in 1 2; do
  mkdir -p "$W/pass$pass/figure"
  for a in hgdp_benchmark big_data rad_seq roh_hgdp; do
    wait_quiet
    echo "=== pass $pass $a start $(date -Is)"
    echo "load before: $(uptime)"
    top -bn1 -o %CPU -w 200 | sed -n '7,11p'
    /usr/bin/time -f "wall %e s, peak RSS %M KB" -o "$W/time_${a}_pass$pass.txt" \
      env R_LIBS="$LIB" Rscript -e "options(keep.source = TRUE); cat('tidypopgen from', find.package('tidypopgen'), as.character(packageVersion('tidypopgen')), '\n'); source('precompile_$a.R', chdir = TRUE)" \
      > "$W/render_${a}_pass$pass.log" 2>&1
    echo "exit $? end $(date -Is); $(cat "$W/time_${a}_pass$pass.txt")"
    echo "load after: $(uptime)"
    grep -n -i 'error' "$W/render_${a}_pass$pass.log" | head -5
    # put back the laptop section of the HGDP benchmark (it needs the laptop)
    if [ "$a" = hgdp_benchmark ]; then
      cat "$W/laptop_section.md" >> "$ART/benchmark_hgdp.Rmd"
      echo "laptop section re-appended"
    fi
    cp "$ART/${out[$a]}.Rmd" "$W/pass$pass/"
    cp -r "$ART/figure/$a" "$W/pass$pass/figure/" 2>/dev/null
  done
done
echo "all done $(date -Is)"
