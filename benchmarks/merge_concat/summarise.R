# Summarise results.tsv (the best of the repetitions of each step) as a
# markdown table: the old route (bcftools, then gen_tibble() on its output)
# against gt_merge_vcfs() / gt_concat_vcfs().
#   Rscript summarise.R [results.tsv]
args <- commandArgs(trailingOnly = TRUE)
res_file <- if (length(args) > 0) args[1] else "results.tsv"
res <- read.delim(res_file, stringsAsFactors = FALSE)
best <- aggregate(
  cbind(wall_s, peak_rss_mb, r_heap_mb) ~ scenario + threads + step,
  data = res, FUN = min, na.action = na.pass
)
fmt_row <- function(sc, t) {
  b <- best[best$scenario == sc & best$threads == t, ]
  get <- function(step, col) {
    v <- b[[col]][b$step == step]
    if (length(v) == 0) NA else v
  }
  new_step <- setdiff(b$step, c("bcftools", "gen_tibble"))
  old <- get("bcftools", "wall_s") + get("gen_tibble", "wall_s")
  new <- get(new_step, "wall_s")
  sprintf(
    "| %s | %d | %s | %s | %s | %.1f | %s | %s |",
    sc, t,
    ifelse(is.na(get("bcftools", "wall_s")), "-",
      sprintf("%.1f", get("bcftools", "wall_s"))),
    ifelse(is.na(get("gen_tibble", "wall_s")), "-",
      sprintf("%.1f", get("gen_tibble", "wall_s"))),
    ifelse(is.na(old), "-", sprintf("%.1f", old)), new,
    ifelse(is.na(old), "-", sprintf("%.1fx", old / new)),
    sprintf("%.0f / %.0f", get(new_step, "peak_rss_mb"), get(new_step, "r_heap_mb"))
  )
}
cat(
  "| scenario | threads | bcftools (s) | gen_tibble (s) | old route (s) |",
  "new function (s) | speed-up | new: peak RSS / R heap (MB) |\n"
)
cat("|---|---|---|---|---|---|---|---|\n")
for (sc in unique(best$scenario)) {
  for (t in sort(unique(best$threads[best$scenario == sc]), decreasing = TRUE)) {
    cat(fmt_row(sc, t), "\n")
  }
}
