args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(tidypopgen, lib.loc = args[1]))
gt <- gt_load("/media/andrea/Elements1/tmp/claude/tidypopgen/bench_inb/data/hgdp650.gt",
  reattach_to = "/media/andrea/Elements1/tmp/claude/tidypopgen/bench_inb/data/hgdp650.rds")
set.seed(42)
sub <- gt[sort(sample(nrow(gt), 301)), ] |>
  select_loci(sort(sample(count_loci(gt), 50001)))
saveRDS(list(grm_sub = pairwise_grm(sub, n_cores = 3, method = "plink"),
  grm_full = pairwise_grm(gt, n_cores = 4, method = "plink")), args[2])
