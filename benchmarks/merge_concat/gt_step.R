# Rscript gt_step.R <lib> <mode> <threads> <backingfile> <files...>
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
mode <- args[2]
threads <- as.integer(args[3])
bk <- args[4]
files <- args[-(1:4)]
suppressPackageStartupMessages(library(tidypopgen, lib.loc = lib))
invisible(gc(reset = TRUE))
t0 <- proc.time()[["elapsed"]]
x <- switch(mode,
  gen_tibble = gen_tibble(files,
    backingfile = bk, n_cores = threads, names_as_int = TRUE, quiet = TRUE
  ),
  merge_union_ref = gt_merge_vcfs(files,
    variants = "union", missing_to_ref = TRUE,
    backingfile = bk, n_cores = threads, quiet = TRUE
  ),
  merge_intersection = gt_merge_vcfs(files,
    backingfile = bk, n_cores = threads, quiet = TRUE
  ),
  concat = gt_concat_vcfs(files,
    backingfile = bk, n_cores = threads, quiet = TRUE
  )
)
elapsed <- proc.time()[["elapsed"]] - t0
# the samples and loci, to check that both routes give the same gen_tibble
saveRDS(
  list(id = x$id, loci = show_loci(x)[, -1], ploidy = indiv_ploidy(x)),
  paste0(bk, "_meta.rds")
)
heap_mb <- sum(gc()[, 6])
cat(sprintf("STEP %.2f %d %d %.0f\n", elapsed, nrow(x), count_loci(x), heap_mb))
