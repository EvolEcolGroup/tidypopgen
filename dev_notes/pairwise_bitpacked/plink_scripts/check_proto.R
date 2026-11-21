# Outputs of the bit-packed kernel on HGDP650 for one build, to check that a
# prototype gives the same results as the base build.
# Usage: Rscript check_proto.R <lib> <out.rds>
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(tidypopgen, lib.loc = args[1]))
d <- file.path(Sys.getenv("TPG_REPO", "/media/andrea/annex/git/tidypopgen"), "vignettes/articles/data/hgdp")
gt <- gt_load(
  file.path(Sys.getenv("TPG_WORK", "/media/andrea/Elements1/tmp/claude/tidypopgen"), "bench_inb/data/hgdp650.gt"),
  # the .gt still points at the backing files' old place on the annex
  reattach_to = file.path(Sys.getenv("TPG_WORK", "/media/andrea/Elements1/tmp/claude/tidypopgen"), "bench_inb/data/hgdp650.rds")
)
pops <- read.delim(file.path(d, "hgdp650_id_pop.txt"))
gt$population <- pops$population[match(gt$id, pops$Id)]
# a subset with scattered loci and individuals, and a partial last word
set.seed(42)
sub <- gt[sort(sample(nrow(gt), 301)), ] |>
  select_loci(sort(sample(count_loci(gt), 50001)))
res <- list(
  king_full = pairwise_king(gt, as_matrix = TRUE, n_cores = 4),
  king_grouped = pairwise_king(dplyr::group_by(gt, population),
    as_matrix = TRUE, n_cores = 4),
  king_sub = pairwise_king(sub, as_matrix = TRUE, n_cores = 3),
  ibs_sub = pairwise_ibs(sub, as_matrix = TRUE, n_cores = 3),
  ibs_counts_sub = pairwise_ibs(sub, type = "adjusted_counts", n_cores = 2),
  as_sub = pairwise_allele_sharing(sub, as_matrix = TRUE, n_cores = 3),
  inb_grouped = indiv_inbreeding(dplyr::group_by(gt, population))
)
saveRDS(res, args[2])
cat("saved", args[2], "\n")
