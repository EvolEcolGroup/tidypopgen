# Read a .ped file with gen_tibble(), and report the time it took.
#   Rscript bench_step.R <lib|worktree> <ped> <backingfile> <n_cores>
# <lib> is a library holding an installed tidypopgen; a path to a source
# tree is loaded with pkgload::load_all(compile = FALSE) instead (compile it
# first with pkgbuild::compile_dll(debug = FALSE), never time an -O0 build)
args <- commandArgs(trailingOnly = TRUE)
src <- args[1]
ped <- args[2]
bk <- args[3]
n_cores <- as.integer(args[4])
if (file.exists(file.path(src, "DESCRIPTION"))) {
  suppressMessages(pkgload::load_all(src, compile = FALSE, quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(tidypopgen, lib.loc = src))
}
invisible(gc(reset = TRUE))
t0 <- proc.time()[["elapsed"]]
x <- suppressWarnings(
  gen_tibble(ped, backingfile = bk, n_cores = n_cores, quiet = TRUE)
)
elapsed <- proc.time()[["elapsed"]] - t0
# the R heap (Mb) at its largest, as seen by gc()
heap_mb <- sum(gc()[, 6])
# the loci and individuals, to check that both versions give the same tables
saveRDS(
  list(
    loci = show_loci(x),
    indiv = as.data.frame(x[, names(x) != "genotypes"]),
    bk_md5 = unname(tools::md5sum(gt_get_file_names(x)[2]))
  ),
  paste0(bk, "_meta.rds")
)
cat(sprintf("STEP %.2f %d %d %.0f\n", elapsed, nrow(x), count_loci(x), heap_mb))
