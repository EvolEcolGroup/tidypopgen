# Create test data to validate windows_indiv_roh() against PLINK 1.9
# (--homozyg)
#
# windows_indiv_roh() implements the sliding-window algorithm of PLINK 1.9
# --homozyg, with the same parameters. This script runs PLINK on:
# - the simulated dataset of inst/extdata/roh (see
#   data-raw/create_roh_dataset.R), with the default parameters and with
#   other values of each of them;
# - a small random dataset with short chromosomes, large gaps and many
#   heterozygous and missing genotypes, to test the edge cases of the
#   algorithm with small windows;
# and saves the parameters of each case (cases.csv), the small dataset and the
# PLINK outputs (<case>.hom, and <case>.hom.summary for some cases) in
# tests/testthat/testdata/plink_roh, where they are used by
# tests/testthat/test_plink_roh.R
#
# It needs PLINK 1.9 (https://www.cog-genomics.org/plink/1.9/): by default
# the `plink` executable on the PATH, or the one given by the environment
# variable PLINK. PLINK 2 has no --homozyg. It can be installed in a conda
# environment with:
# conda create -n plink19 -c bioconda -c conda-forge plink
# and the script then run with:
# PLINK=$(conda run -n plink19 which plink) Rscript \
#   data-raw/validate/create_plink_roh_test_data.R
# (the test data were created with PLINK v1.9.0 64-bit (27 Sep 2026))
#
# This script should be run from the root of the package

devtools::load_all()

plink <- Sys.getenv("PLINK", "plink")
out_dir <- "tests/testthat/testdata/plink_roh"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
version <- system2(plink, "--version", stdout = TRUE)
if (!grepl("^PLINK v1\\.9", version[1])) {
  stop("PLINK 1.9 is needed, but ", plink, " is ", version[1])
}
writeLines(version[1], file.path(out_dir, "plink_version.txt"))

# the small random dataset ---------------------------------------------------
set.seed(2026)
n_per_chrom <- c(4, 9, 30, 60)
n_loci <- sum(n_per_chrom)
chrom <- rep(seq_along(n_per_chrom), n_per_chrom)
# loci 1 to 9 kb apart, with some gaps of 20 to 40 kb
pos <- unlist(lapply(n_per_chrom, function(k) {
  cumsum(sample(c(1000:9000, 20000:40000), k,
    replace = TRUE,
    prob = c(rep(0.95 / 8001, 8001), rep(0.05 / 20001, 20001))
  ))
}))
n_indiv <- 12
# individuals with increasing proportions of heterozygous genotypes, and
# some missing genotypes
p_het <- rep(c(0, 0.02, 0.05, 0.1, 0.2, 0.4), each = 2)
genotypes <- t(vapply(p_het, function(p) {
  sample(c(0, 1, 2, NA), n_loci,
    replace = TRUE,
    prob = c((1 - p - 0.05) / 2, p, (1 - p - 0.05) / 2, 0.05)
  )
}, numeric(n_loci)))
small_gt <- gen_tibble(
  genotypes,
  indiv_meta = data.frame(
    id = paste0("ind", seq_len(n_indiv)),
    population = "pop"
  ),
  loci = data.frame(
    name = paste0("snp", seq_len(n_loci)),
    chromosome = chrom,
    position = as.integer(pos),
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  ),
  backingfile = tempfile(),
  quiet = TRUE
)
gt_as_plink(small_gt, file = file.path(out_dir, "small"))

# the cases ------------------------------------------------------------------
# one row per case, with the dataset (the simulated one, "roh_sim", or the
# small one, "small") and the arguments of windows_indiv_roh() (NA for the
# default)
cases <- data.frame(
  case = c(
    "sim_default", "sim_window", "sim_het", "sim_lengths",
    "small_a", "small_b", "small_c"
  ),
  dataset = c(rep("roh_sim", 4), rep("small", 3)),
  window_snp = c(NA, 20, NA, 30, 5, 3, 8),
  window_het = c(NA, 0, NA, 2, 1, 0, 2),
  window_missing = c(NA, 1, NA, NA, 1, 0, 2),
  window_threshold = c(NA, 0.1, NA, 0.3, 0.05, 0.5, 0.2),
  min_snp = c(NA, 30, 50, NA, 3, 2, 5),
  min_kb = c(NA, 300, 500, 800, 1, 0.5, 10),
  max_kb_per_snp = c(NA, NA, NA, 20, 10, 5, 50),
  max_gap_kb = c(NA, 500, NA, NA, 15, 30, 1000),
  max_het = c(NA, NA, 0, 2, 1, 0, 2),
  subtract_1_from_lengths = c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE)
)
utils::write.csv(cases, file.path(out_dir, "cases.csv"), row.names = FALSE)

# run PLINK on each case -----------------------------------------------------
plink_flags <- c(
  window_snp = "--homozyg-window-snp",
  window_het = "--homozyg-window-het",
  window_missing = "--homozyg-window-missing",
  window_threshold = "--homozyg-window-threshold",
  min_snp = "--homozyg-snp",
  min_kb = "--homozyg-kb",
  max_kb_per_snp = "--homozyg-density",
  max_gap_kb = "--homozyg-gap",
  max_het = "--homozyg-het"
)
for (i in seq_len(nrow(cases))) {
  bfile <- if (cases$dataset[i] == "roh_sim") {
    file.path("inst/extdata/roh/roh_sim")
  } else {
    file.path(out_dir, "small")
  }
  args <- c("--bfile", bfile, "--homozyg")
  if (cases$subtract_1_from_lengths[i]) {
    args <- c(args, "subtract-1-from-lengths")
  }
  for (arg in names(plink_flags)) {
    if (!is.na(cases[[arg]][i])) {
      args <- c(args, plink_flags[[arg]], format(cases[[arg]][i]))
    }
  }
  out <- tempfile()
  status <- system2(plink, c(args, "--out", out))
  if (status != 0) stop("PLINK failed on case ", cases$case[i])
  file.copy(paste0(out, ".hom"),
    file.path(out_dir, paste0(cases$case[i], ".hom")),
    overwrite = TRUE
  )
  # the number of individuals in a run at each locus (AFF + UNAFF), for the
  # default parameters and the small dataset (to keep the files small)
  if (cases$case[i] == "sim_default" || cases$dataset[i] == "small") {
    file.copy(paste0(out, ".hom.summary"),
      file.path(out_dir, paste0(cases$case[i], ".hom.summary")),
      overwrite = TRUE
    )
  }
}
