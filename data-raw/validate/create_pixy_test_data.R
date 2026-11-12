# Create test data to validate Tajima's D with missing data against pixy
#
# pixy (Bailey et al. 2025, Molecular Ecology Resources, e14104) accounts for
# missing genotypes by using the number of non-missing allele copies at each
# site as its sample size. tidypopgen implements the same estimator in
# pop_tajimas_d(), windows_pop_tajimas_d() and jackknife_pop_tajimas_d().
# This script simulates a small dataset with ragged missingness (i.e. a
# different number of genotyped individuals at each locus), runs pixy on it,
# and saves both the input and the pixy output in
# tests/testthat/testdata/pixy, where they are used by
# tests/testthat/test_pixy_tajimas_d.R
#
# Note that pixy versions 2.1.0 to 2.2.1 inflated the denominator of Tajima's
# D when the number of observed alleles varied among sites (see
# https://pixy.readthedocs.io/en/latest/interpreting_tajima_d.html), so a
# later version must be used to generate these results.
#
# Create the conda environment with:
# conda create -n pixy -c conda-forge -c bioconda python=3.12 pixy htslib
#
# This script should be run from the root of the package

devtools::load_all()

out_dir <- "tests/testthat/testdata/pixy"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# helper to run a command within the pixy conda environment
run_in_pixy_env <- function(command, args) {
  # use the conda executable of the active installation, if available
  conda_bin <- Sys.getenv("CONDA_EXE", unset = "conda")
  res <- system2(conda_bin,
    args = c("run", "-n", "pixy", command, args),
    stdout = TRUE, stderr = TRUE
  )
  if (!is.null(attr(res, "status")) && attr(res, "status") != 0) {
    stop(paste(res, collapse = "\n"))
  }
  res
}

# record the version of pixy used to generate the data
pixy_version <- run_in_pixy_env("pixy", "--version")
# only keep the line with the version number (skipping the ascii logo)
pixy_version <- grep("version", pixy_version, value = TRUE)
print(pixy_version)
writeLines(pixy_version, file.path(out_dir, "pixy_version.txt"))

###############################################################################
# Simulate genotypes
###############################################################################
set.seed(123)
# three populations of diploid individuals
pop_sizes <- c(popA = 12, popB = 10, popC = 8)
population <- rep(names(pop_sizes), pop_sizes)
n_indiv <- length(population)
# two chromosomes
n_loci_chr <- c(chr1 = 150, chr2 = 100)
chr_length <- c(chr1 = 50000, chr2 = 30000)
n_loci <- sum(n_loci_chr)

# allele frequencies drawn from a U-shaped distribution, with some
# differentiation among populations
anc_freq <- stats::rbeta(n_loci, 0.3, 0.3)
pop_freq <- sapply(names(pop_sizes), function(pop) {
  stats::rbeta(n_loci, anc_freq * 20 + 0.01, (1 - anc_freq) * 20 + 0.01)
})
# draw the genotypes (dosages)
genotypes <- matrix(NA_integer_, nrow = n_indiv, ncol = n_loci)
for (i_ind in seq_len(n_indiv)) {
  genotypes[i_ind, ] <- stats::rbinom(
    n_loci,
    size = 2,
    prob = pop_freq[, population[i_ind]]
  )
}
# ragged missingness: each locus has its own missing rate (up to 50%)
missing_rate <- stats::runif(n_loci, 0, 0.5)
for (i_locus in seq_len(n_loci)) {
  is_missing <- stats::runif(n_indiv) < missing_rate[i_locus]
  genotypes[is_missing, i_locus] <- NA
}
# a few loci where a whole population is missing, or a single individual is
# genotyped
genotypes[population == "popA", 10] <- NA
genotypes[population == "popB", 200] <- NA
genotypes[population == "popC", 50][-1] <- NA
genotypes[population == "popC", 50][1] <- 1

# loci table
positions <- unlist(lapply(names(n_loci_chr), function(chr) {
  sort(sample.int(chr_length[chr], n_loci_chr[chr]))
}))
loci <- data.frame(
  name = paste0("rs", seq_len(n_loci)),
  chromosome = rep(names(n_loci_chr), n_loci_chr),
  position = as.integer(positions),
  genetic_dist = 0,
  allele_ref = "A",
  allele_alt = "T"
)
indiv_meta <- data.frame(
  id = paste0("ind", seq_len(n_indiv)),
  population = population
)
pixy_gt <- gen_tibble(
  x = genotypes,
  loci = loci,
  indiv_meta = indiv_meta,
  backingfile = tempfile(),
  quiet = TRUE
)

###############################################################################
# Write the input files for pixy
###############################################################################
vcf_path <- file.path(out_dir, "pixy_test.vcf")
gt_as_vcf(pixy_gt, file = vcf_path, overwrite = TRUE)
# pixy requires a bgzipped and indexed vcf (we keep the uncompressed vcf for
# the tests, and compress a temporary copy)
tmp_dir <- tempfile()
dir.create(tmp_dir)
file.copy(vcf_path, file.path(tmp_dir, "pixy_test.vcf"))
run_in_pixy_env("bgzip", c("-f", file.path(tmp_dir, "pixy_test.vcf")))
vcf_gz_path <- file.path(tmp_dir, "pixy_test.vcf.gz")
run_in_pixy_env("tabix", c("-p", "vcf", vcf_gz_path))
# populations file (no header: sample id, population)
pop_path <- file.path(out_dir, "pixy_populations.txt")
utils::write.table(indiv_meta,
  file = pop_path, sep = "\t",
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

###############################################################################
# Run pixy
###############################################################################
# windows of 10 kb, and a single window spanning each chromosome
for (window_size in c("10000", "100000")) {
  run_in_pixy_env("pixy", c(
    "--stats", "pi", "watterson_theta", "tajima_d",
    "--vcf", vcf_gz_path,
    "--populations", pop_path,
    "--window_size", window_size,
    "--bypass_invariant_check",
    "--n_cores", "1",
    "--output_folder", tmp_dir,
    "--output_prefix", paste0("pixy_", window_size)
  ))
  # copy the tajima's d output to the test data folder
  file.copy(
    file.path(tmp_dir, paste0("pixy_", window_size, "_tajima_d.txt")),
    out_dir,
    overwrite = TRUE
  )
}
