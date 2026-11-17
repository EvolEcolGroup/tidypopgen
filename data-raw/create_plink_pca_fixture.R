# Test fixture for gt_pca(): a small simulated dataset with very high
# missingness, and the outputs of PLINK 1.9 on it (--make-rel square, --pca).
#
# Three populations of 20 diploid individuals (Balding-Nichols, Fst about
# 0.05), 1,000 SNPs, one of them monomorphic; within each population, blocks
# of individuals have 0, 20, 50, 80 and 90% missing genotypes.
#
# Run from the package root, with PLINK 1.9 on the PATH (or in PLINK):
#   Rscript data-raw/create_plink_pca_fixture.R
devtools::load_all(quiet = TRUE)
out_dir <- "tests/testthat/testdata/plink_pca"
plink <- Sys.getenv("PLINK", "plink")

set.seed(42)
n_pop <- 3
n_ind <- 20
m <- 1000
p_anc <- stats::runif(m, 0.1, 0.9)
f <- 0.05
freqs <- replicate(n_pop, stats::rbeta(m, p_anc * (1 - f) / f, (1 - p_anc) * (1 - f) / f))
geno <- do.call(rbind, lapply(seq_len(n_pop), function(k) {
  t(replicate(n_ind, stats::rbinom(m, 2, freqs[, k])))
}))
geno[, 1] <- 0 # a monomorphic locus
rate <- rep(rep(c(0, 0.2, 0.5, 0.8, 0.9), each = n_ind / 5), n_pop)
geno[matrix(stats::runif(length(geno)), nrow(geno)) < rate] <- NA

gt <- gen_tibble(geno,
  indiv_meta = data.frame(
    id = paste0("ind", seq_len(nrow(geno))),
    population = rep(paste0("pop", seq_len(n_pop)), each = n_ind)
  ),
  loci = data.frame(
    name = paste0("snp", seq_len(m)), chromosome = "1",
    position = seq_len(m) * 1000L, genetic_dist = 0,
    allele_ref = "A", allele_alt = "G"
  ),
  backingfile = tempfile(), quiet = TRUE
)
gt_as_plink(gt, file = file.path(out_dir, "missing.bed"), type = "bed")

bfile <- file.path(out_dir, "missing")
system2(plink, c("--bfile", bfile, "--make-rel", "square", "--out", bfile))
system2(plink, c("--bfile", bfile, "--pca", "5", "--out", bfile))
writeLines(system2(plink, "--version", stdout = TRUE),
  file.path(out_dir, "plink_version.txt")
)
# keep only the files the tests use
unlink(paste0(bfile, c(".log", ".nosex", ".rel.id")))
