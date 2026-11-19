# Test fixture for loci_ld_prune() and loci_ld_clump(): a small simulated
# dataset with linkage disequilibrium and missing genotypes, and the outputs
# of PLINK 1.9 on it (--r2 for nearby pairs of loci, and --indep-pairwise
# with several windows, steps and thresholds).
#
# 60 diploid individuals, two chromosomes of 100 SNPs at irregular distances
# (1 to 20 kb). Each haplotype copies the previous locus with probability
# 0.85, so that nearby loci are in LD. Each individual misses 0 to 50% of its
# genotypes. Locus 10 is monomorphic, loci 20 and 21 are identical (a tie in
# minor allele frequency), and locus 30 is genotyped in 2 individuals only.
#
# Run from the package root, with PLINK 1.9 on the PATH (or in PLINK):
#   Rscript data-raw/validate/create_plink_ld_test_data.R
devtools::load_all(quiet = TRUE)
out_dir <- "tests/testthat/testdata/plink_ld"
plink <- Sys.getenv("PLINK", "plink")

set.seed(7)
n <- 60
m <- 200
haplotypes <- function() {
  h <- matrix(0L, n, m)
  p <- stats::runif(m, 0.05, 0.95)
  h[, 1] <- stats::rbinom(n, 1, p[1])
  for (l in 2:m) {
    copy <- stats::rbinom(n, 1, 0.85)
    h[, l] <- ifelse(copy == 1, h[, l - 1], stats::rbinom(n, 1, p[l]))
  }
  h
}
geno <- haplotypes() + haplotypes()
geno[, 10] <- 0L
geno[, 21] <- geno[, 20]
rate <- stats::runif(n, 0, 0.5)
geno[matrix(stats::runif(n * m), n) < rate] <- NA
geno[-(1:2), 30] <- NA

gt <- gen_tibble(geno,
  indiv_meta = data.frame(id = paste0("ind", seq_len(n)), population = "pop"),
  loci = data.frame(
    name = paste0("snp", seq_len(m)), chromosome = rep(1:2, each = m / 2),
    position = unlist(lapply(1:2, function(k) {
      cumsum(sample(1000:20000, m / 2))
    })),
    genetic_dist = 0, allele_ref = "A", allele_alt = "G"
  ),
  backingfile = tempfile(), quiet = TRUE
)
gt_as_plink(gt, file = file.path(out_dir, "ld.bed"), type = "bed")

bfile <- file.path(out_dir, "ld")
# r^2 of the pairs of loci within 5 loci of each other
system2(plink, c(
  "--bfile", bfile, "--r2", "--ld-window", "5", "--ld-window-kb", "99999",
  "--ld-window-r2", "0", "--allow-no-sex", "--out", bfile
))
# pruning, with windows in loci and in kb
prune <- list(
  w10_s1_r05 = c("10", "1", "0.5"), w20_s5_r02 = c("20", "5", "0.2"),
  w2_s1_r08 = c("2", "1", "0.8"), w50kb_s1_r02 = c("50kb", "1", "0.2")
)
for (name in names(prune)) {
  out <- file.path(out_dir, name)
  system2(plink, c(
    "--bfile", bfile, "--indep-pairwise", prune[[name]], "--allow-no-sex",
    "--out", out
  ))
  unlink(paste0(out, c(".log", ".nosex", ".prune.out")))
}
writeLines(system2(plink, "--version", stdout = TRUE),
  file.path(out_dir, "plink_version.txt")
)
# keep only the files the tests use
unlink(paste0(bfile, c(".log", ".nosex")))
