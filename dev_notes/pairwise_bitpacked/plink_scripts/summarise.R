# Medians (min-max) of the timings in out/*.tsv
# Usage: Rscript summarise.R <out dir>
out <- commandArgs(trailingOnly = TRUE)[1]
fmt <- function(x) {
  sprintf("%.3f (%.3f-%.3f) n=%d", stats::median(x), min(x), max(x),
    length(x))
}
for (f in list.files(out, "^bench_plink_[0-9]+\\.tsv$", full.names = TRUE)) {
  d <- utils::read.delim(f, header = FALSE,
    col.names = c("time", "task", "threads", "wall", "user", "sys", "load"))
  cat("\n#", basename(f), "(wall s; user s)\n")
  for (t in unique(d$task)) {
    x <- d[d$task == t, ]
    cat(sprintf("%-18s T=%d wall %s; user %s; load %s\n", t, x$threads[1],
      fmt(x$wall), fmt(x$user), paste(x$load, collapse = ",")))
  }
}
for (f in list.files(out, "^(bench_tpg|proto)_[0-9]+\\.tsv$",
  full.names = TRUE)) {
  d <- utils::read.delim(f, header = FALSE,
    col.names = c("time", "version", "task", "threads", "t_load", "first",
      "second", "user", "load"))
  cat("\n#", basename(f), "(first call, second call, gt_load; s)\n")
  for (v in unique(d$version)) {
    for (t in unique(d$task)) {
      x <- d[d$version == v & d$task == t, ]
      if (nrow(x) == 0) next
      cat(sprintf("%-5s %-12s T=%d first %s; second %s; load %s\n", v, t,
        x$threads[1], fmt(x$first), fmt(x$second), fmt(x$t_load)))
    }
  }
}
