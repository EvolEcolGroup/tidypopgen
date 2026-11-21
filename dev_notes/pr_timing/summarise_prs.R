# Markdown summary of the desktop timings of PRs #239 and #240
# (run_desktop_prs.sh): medians (min-max) in seconds, and whether the
# versions gave the same outputs.
# Usage: Rscript summarise_prs.R <results dir> > summary.md
d <- commandArgs(trailingOnly = TRUE)[1]
med <- function(x) {
  sprintf("%.2f (%.2f-%.2f)", stats::median(x), min(x), max(x))
}
ratio <- function(base, new) sprintf("%.2fx", stats::median(base) / stats::median(new))

cat("Desktop (", Sys.info()[["nodename"]], "), ", format(Sys.Date()),
  "; medians (min-max) in seconds.\n", sep = "")

# KING (run_proto.sh): base = dev, p3 = king_avx2, p3s = king_avx2 without
# AVX2; each run is the faster of the two calls in its process
king_cols <- c("time", "version", "task", "threads", "t_load", "first",
  "second", "user", "load")
king_files <- list.files(d, "^proto_[0-9]+\\.tsv$", full.names = TRUE)
if (length(king_files) > 0) {
  cat("\n## KING-robust on HGDP650 (PR #239)\n\n")
  cat("| task | threads | dev | king_avx2 (AVX2) | king_avx2 (scalar) | AVX2 / scalar vs dev |\n")
  cat("|---|---|---|---|---|---|\n")
  for (f in king_files) {
    k <- utils::read.delim(f, header = FALSE, col.names = king_cols)
    k$best <- pmin(k$first, k$second)
    for (task in unique(k$task)) {
      x <- function(v) k$best[k$version == v & k$task == task]
      if (length(x("base")) == 0) next
      cat(sprintf("| %s | %d | %s | %s | %s | %s / %s |\n", task,
        k$threads[1], med(x("base")), med(x("p3")), med(x("p3s")),
        ratio(x("base"), x("p3")), ratio(x("base"), x("p3s"))))
    }
  }
}
if (file.exists(file.path(d, "checks.txt"))) {
  cat("\nOutputs (check_proto.R):\n\n")
  cat(paste0("- ", readLines(file.path(d, "checks.txt"))), sep = "\n")
}

# gt_merge_vcfs() (run_timing.sh)
runs_file <- file.path(d, "runs.tsv")
if (file.exists(runs_file)) {
  r <- utils::read.delim(runs_file)
  cat("\n## gt_merge_vcfs() (PR #240)\n\n")
  cat("| scenario | threads | dev | vcf_merge_speed | speed-up | same .bk and sidecar |\n")
  cat("|---|---|---|---|---|---|\n")
  for (t in sort(unique(r$n_cores))) {
    for (sc in c("pops", "sgdp", "many")) {
      x <- r[r$scenario == sc & r$n_cores == t, ]
      if (nrow(x) == 0) next
      a <- x$wall_s[x$version == "dev"]
      b <- x$wall_s[x$version == "vcf"]
      if (length(a) == 0 || length(b) == 0) next
      same <- length(unique(x$bk_md5)) == 1
      cat(sprintf("| %s | %d | %s | %s | %s | %s |\n", sc, t, med(a), med(b),
        ratio(a, b), if (same) "yes" else "NO"))
    }
  }
}

# gen_tibble() (time_gen_tibble.R)
gt_file <- file.path(d, "gen_tibble.tsv")
if (file.exists(gt_file)) {
  g <- utils::read.delim(gt_file)
  cat("\n## gen_tibble() on one VCF (PR #240)\n\n")
  cat("| file | threads | dev | vcf_merge_speed | speed-up | same .bk |\n")
  cat("|---|---|---|---|---|---|\n")
  for (t in sort(unique(g$n_cores))) {
    for (f in unique(g$file)) {
      x <- g[g$file == f & g$n_cores == t, ]
      a <- x$wall_s[x$version == "dev"]
      b <- x$wall_s[x$version == "vcf"]
      if (length(a) == 0 || length(b) == 0) next
      cat(sprintf("| %s | %d | %s | %s | %s | %s |\n", f, t, med(a), med(b),
        ratio(a, b), if (length(unique(x$bk_md5)) == 1) "yes" else "NO"))
    }
  }
}

# the commits and the load
if (file.exists(runs_file)) {
  r <- utils::read.delim(runs_file)
  cat("\nCommits:", paste(unique(paste0(r$version, "=", r$commit)),
    collapse = ", "), "\n")
  cat("Load (1 min) before the merge runs:", med(r$load1), "\n")
}
