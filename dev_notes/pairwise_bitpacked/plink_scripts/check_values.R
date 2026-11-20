# Compare tidypopgen's KING-robust with PLINK 2 on HGDP650 (full matrix, and
# within populations). Usage: Rscript check_values.R <lib> <plink check dir>
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(tidypopgen, lib.loc = args[1]))
chk <- args[2]
d <- "/media/andrea/annex/git/tidypopgen/vignettes/articles/data/hgdp"
gt <- gt_load("/media/andrea/annex/tmp/claude/tidypopgen/bench_inb/data/hgdp650.gt")
pops <- read.delim(file.path(d, "hgdp650_id_pop.txt"))
gt$population <- pops$population[match(gt$id, pops$Id)]

# full matrix
k_tpg <- pairwise_king(gt, as_matrix = TRUE, n_cores = 1)
saveRDS(k_tpg, file.path(chk, "tpg_full.rds"))
ids <- read.delim(file.path(chk, "full.king.id"))
n <- nrow(ids)
stopifnot(identical(ids$IID, gt$id))
# triangle bin: rows 2..n, row i has i - 1 values (lower triangle, by row)
v <- readBin(file.path(chk, "full.king.bin"), "double", n = n * (n - 1) / 2)
k_pl <- matrix(NA_real_, n, n)
k_pl[upper.tri(k_pl)] <- v # column-major upper == row-major lower
k_pl[lower.tri(k_pl)] <- t(k_pl)[lower.tri(k_pl)]
diff_full <- abs(k_pl - unclass(k_tpg))
diag(diff_full) <- NA
cat(sprintf("full: max |diff| = %.3g, mean |diff| = %.3g, pairs = %d\n",
  max(diff_full, na.rm = TRUE), mean(diff_full, na.rm = TRUE),
  sum(!is.na(diff_full)) / 2))
w <- which(diff_full == max(diff_full, na.rm = TRUE), arr.ind = TRUE)[1, ]
cat(sprintf("  worst pair %s %s: plink %.8f tidypopgen %.8f\n",
  gt$id[w[1]], gt$id[w[2]], k_pl[w[1], w[2]], k_tpg[w[1], w[2]]))
cat(sprintf("  diagonal of tidypopgen: %s\n",
  paste(range(diag(k_tpg)), collapse = " - ")))

# within populations
k_grp <- pairwise_king(dplyr::group_by(gt, population), as_matrix = TRUE,
  n_cores = 1)
keys <- dplyr::group_keys(dplyr::group_by(gt, population))$population
kin0 <- read.delim(file.path(chk, "within.kin0"))
kin0$tpg <- NA_real_
for (g in seq_along(keys)) {
  rows <- which(kin0$X.FID1 == keys[g])
  m <- k_grp[[g]]
  kin0$tpg[rows] <- m[cbind(kin0$IID1[rows], kin0$IID2[rows])]
}
dw <- abs(kin0$KINSHIP - kin0$tpg)
cat(sprintf("within: %d pairs (PLINK), %d (tidypopgen), max |diff| = %.3g\n",
  nrow(kin0), sum(vapply(k_grp, function(m) nrow(m) * (nrow(m) - 1) / 2,
    numeric(1))), max(dw)))
cat(sprintf("  PLINK .kin0 prints 6 significant digits; max relative diff = %.3g\n",
  max(dw / abs(kin0$tpg))))
# the grouped matrices agree with the full matrix of tidypopgen
dg <- 0
for (g in seq_along(keys)) {
  ii <- rownames(k_grp[[g]])
  dg <- max(dg, abs(k_grp[[g]] - k_tpg[ii, ii]))
}
cat(sprintf("tidypopgen grouped vs full: max |diff| = %.3g\n", dg))
