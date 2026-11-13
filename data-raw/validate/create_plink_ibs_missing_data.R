# Create test data to validate pairwise_ibs() against PLINK 1.9 --distance,
# with both corrections for missing genotypes
#
# By default, PLINK 1.9 --distance corrects for missing genotypes with weights
# that depend on the allele frequencies (missing_correction = "allele_freq"
# in pairwise_ibs()); with the 'flat-missing' modifier, all loci have the same
# weight (missing_correction = "flat"). This script simulates a small dataset
# with many missing genotypes (individuals and loci with 0 to 40% missing
# genotypes, monomorphic loci, a locus missing in all individuals and one
# missing in all founders), saves it as ibs_missing.bed/.bim/.fam in
# inst/extdata/related, and runs PLINK on it:
# - ibs_missing.dist.bin and ibs_missing.mibs.bin:
#   --distance square bin allele-ct ibs (the weighted correction)
# - ibs_missing_flat.dist.bin and ibs_missing_flat.mibs.bin:
#   --distance square bin allele-ct ibs flat-missing (the flat correction)
# - ibs_missing_nonfounders.dist.bin: --distance square bin allele-ct, with
#   ibs_missing_nonfounders.fam, in which the last 8 individuals have parents,
#   so that the allele frequencies (and so the weights) are computed from the
#   first 16 individuals (the founders) only
# The 'bin' outputs hold the matrices as little-endian doubles, so that they
# can be compared at full precision (the text outputs have 6 significant
# digits). They are used by tests/testthat/test_snp_ibs.R.
#
# It needs PLINK 1.9 (https://www.cog-genomics.org/plink/1.9/): by default
# the `plink` executable on the PATH, or the one given by the environment
# variable PLINK. It can be installed in a conda environment with:
# conda create -n plink19 -c bioconda -c conda-forge plink
# and the script then run with:
# PLINK=$(conda run -n plink19 which plink) Rscript \
#   data-raw/validate/create_plink_ibs_missing_data.R
# (the test data were created with PLINK v1.9.0 64-bit (27 Sep 2026))
#
# This script should be run from the root of the package

devtools::load_all()

plink <- Sys.getenv("PLINK", "plink")
out_dir <- "inst/extdata/related"
version <- system2(plink, "--version", stdout = TRUE)
if (!grepl("^PLINK v1\\.9", version[1])) {
  stop("PLINK 1.9 is needed, but ", plink, " is ", version[1])
}
message("Using ", version[1])

# simulate the genotypes ------------------------------------------------------
set.seed(2026)
n_indiv <- 24
n_loci <- 300
# allele frequencies, with some monomorphic loci (fixed for either allele)
freqs <- stats::runif(n_loci, 0.02, 0.98)
freqs[sample(n_loci, 20)] <- 0
freqs[sample(n_loci, 10)] <- 1
genotypes <- vapply(freqs, function(p) {
  stats::rbinom(n_indiv, 2, p)
}, numeric(n_indiv))
# missing genotypes: individuals and loci have different missing rates (0 to
# 40% and 0 to 30%), so that pairs of individuals share many missing loci
p_missing_indiv <- seq(0, 0.4, length.out = n_indiv)
p_missing_loci <- stats::runif(n_loci, 0, 0.3)
p_missing <- outer(p_missing_indiv, p_missing_loci, function(a, b) {
  1 - (1 - a) * (1 - b)
})
genotypes[stats::runif(length(genotypes)) < p_missing] <- NA
# a locus missing in all individuals, and one missing in all the founders of
# ibs_missing_nonfounders.fam (individuals 1 to 16)
genotypes[, 5] <- NA
genotypes[1:16, 6] <- NA
message(
  "Proportion of missing genotypes: ",
  round(mean(is.na(genotypes)), 3)
)

ibs_gt <- gen_tibble(
  x = genotypes,
  loci = data.frame(
    name = paste0("rs", seq_len(n_loci)),
    chromosome = rep(1:2, each = n_loci / 2),
    position = rep(seq_len(n_loci / 2) * 1000L, 2),
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  ),
  indiv_meta = data.frame(
    id = paste0("ind", seq_len(n_indiv)),
    population = "pop1"
  ),
  backingfile = tempfile(),
  quiet = TRUE
)
bfile <- file.path(out_dir, "ibs_missing")
gt_as_plink(ibs_gt, file = paste0(bfile, ".bed"))

# a .fam in which the last 8 individuals are the offspring of pairs of the
# first 16
fam <- utils::read.table(paste0(bfile, ".fam"))
fam$V3[17:24] <- fam$V2[seq(1, 16, by = 2)]
fam$V4[17:24] <- fam$V2[seq(2, 16, by = 2)]
utils::write.table(fam, paste0(bfile, "_nonfounders.fam"),
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

# run PLINK ------------------------------------------------------------------
run_plink <- function(args, out) {
  res <- system2(plink, c(
    "--bfile", bfile, args, "--allow-no-sex", "--out", out
  ), stdout = TRUE)
  if (!is.null(attr(res, "status"))) {
    stop(paste(res, collapse = "\n"))
  }
  # keep only the matrices
  unlink(Sys.glob(paste0(out, c(".log", ".nosex", "*.id"))))
}
run_plink(
  c("--distance", "square", "bin", "allele-ct", "ibs"),
  file.path(out_dir, "ibs_missing")
)
run_plink(
  c("--distance", "square", "bin", "allele-ct", "ibs", "flat-missing"),
  file.path(out_dir, "ibs_missing_flat")
)
run_plink(
  c(
    "--fam", paste0(bfile, "_nonfounders.fam"),
    "--distance", "square", "bin", "allele-ct"
  ),
  file.path(out_dir, "ibs_missing_nonfounders")
)
