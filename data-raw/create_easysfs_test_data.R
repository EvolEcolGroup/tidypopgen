# Create test data to validate the site frequency spectra against easySFS
#
# easySFS (https://github.com/isaacovercast/easySFS) computes 1D, 2D and
# multidimensional site frequency spectra from a vcf, dealing with missing
# genotypes by projecting each population down to a fixed number of allele
# copies (with dadi's hypergeometric projection). tidypopgen implements the
# same approach in pop_sfs(), pairwise_pop_sfs(), nwise_pop_sfs() and
# pop_sfs_preview(). This script simulates small datasets with ragged
# missingness (i.e. a different number of genotyped individuals at each
# locus), runs easySFS on them, and saves both the inputs and the easySFS
# outputs in tests/testthat/testdata/easysfs, where they are used by
# tests/testthat/test_easysfs.R
#
# easySFS is run at a fixed commit (EASYSFS_COMMIT below). It needs python 3
# with numpy, scipy and pandas, e.g. in a conda environment created with:
# conda create -n easysfs -c conda-forge python=3.12 numpy scipy pandas
# By default, easySFS is run with `conda run -n easysfs python`; set the
# environment variable EASYSFS_PYTHON to the path of a python executable to
# use it instead.
#
# This script should be run from the root of the package

devtools::load_all()

out_dir <- "tests/testthat/testdata/easysfs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

EASYSFS_COMMIT <- "87d29965f26398fea217a0df520127f39a2368d4"

# clone easySFS at the fixed commit
easysfs_dir <- tempfile("easySFS")
system2("git", c(
  "clone", "--quiet", "https://github.com/isaacovercast/easySFS",
  easysfs_dir
))
system2("git", c("-C", easysfs_dir, "checkout", "--quiet", EASYSFS_COMMIT))
writeLines(
  c(
    "easySFS (https://github.com/isaacovercast/easySFS)",
    paste("commit", EASYSFS_COMMIT)
  ),
  file.path(out_dir, "easysfs_version.txt")
)

# run easySFS with the given arguments, returning its standard output
run_easysfs <- function(args) {
  python <- Sys.getenv("EASYSFS_PYTHON", unset = "")
  script <- file.path(easysfs_dir, "easySFS.py")
  if (python == "") {
    conda_bin <- Sys.getenv("CONDA_EXE", unset = "conda")
    res <- system2(conda_bin,
      args = c("run", "-n", "easysfs", "python", script, args),
      stdout = TRUE, stderr = TRUE
    )
  } else {
    res <- system2(python, args = c(script, args), stdout = TRUE, stderr = TRUE)
  }
  if (!is.null(attr(res, "status")) && attr(res, "status") != 0) {
    stop(paste(res, collapse = "\n"))
  }
  res
}

# write a vcf with genotypes (a matrix of alternate allele counts,
# individuals x loci, NA for missing) of the given ploidy
write_test_vcf <- function(genotypes, ploidy, ids, loci, file) {
  if (ploidy == 2) {
    codes <- c("0/0", "0/1", "1/1")
    missing_code <- "./."
  } else {
    codes <- c("0", "1")
    missing_code <- "."
  }
  gt_strings <- matrix(codes[genotypes + 1], nrow = nrow(genotypes))
  gt_strings[is.na(genotypes)] <- missing_code
  header <- c(
    "##fileformat=VCFv4.2",
    paste0("##contig=<ID=", unique(loci$chromosome), ">"),
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", ids
    ), collapse = "\t")
  )
  body <- vapply(seq_len(nrow(loci)), function(j) {
    paste(c(
      loci$chromosome[j], loci$position[j], loci$name[j], loci$allele_ref[j],
      loci$allele_alt[j], ".", "PASS", ".", "GT", gt_strings[, j]
    ), collapse = "\t")
  }, character(1))
  writeLines(c(header, body), file)
}

# simulate genotypes for populations of the given sizes, with ragged
# missingness
simulate_genotypes <- function(pop_sizes, n_loci, ploidy) {
  population <- rep(names(pop_sizes), pop_sizes)
  n_indiv <- length(population)
  # allele frequencies drawn from a U-shaped distribution, with some
  # differentiation among populations
  anc_freq <- stats::rbeta(n_loci, 0.3, 0.3)
  pop_freq <- sapply(names(pop_sizes), function(pop) {
    stats::rbeta(n_loci, anc_freq * 20 + 0.01, (1 - anc_freq) * 20 + 0.01)
  })
  genotypes <- matrix(NA_integer_, nrow = n_indiv, ncol = n_loci)
  for (i_ind in seq_len(n_indiv)) {
    genotypes[i_ind, ] <- stats::rbinom(
      n_loci,
      size = ploidy,
      prob = pop_freq[, population[i_ind]]
    )
  }
  # ragged missingness: each locus has its own missing rate (up to 50%)
  missing_rate <- stats::runif(n_loci, 0, 0.5)
  for (i_locus in seq_len(n_loci)) {
    is_missing <- stats::runif(n_indiv) < missing_rate[i_locus]
    genotypes[is_missing, i_locus] <- NA
  }
  list(genotypes = genotypes, population = population)
}

# a loci table over 3 chromosomes
simulate_loci <- function(n_loci) {
  chromosome <- rep(c("chr1", "chr2", "chr3"), length.out = n_loci)
  chromosome <- sort(chromosome)
  data.frame(
    name = paste0("rs", seq_len(n_loci)),
    chromosome = chromosome,
    position = as.integer(stats::ave(
      seq_len(n_loci), chromosome,
      FUN = function(x) seq_along(x) * 100
    )),
    allele_ref = "A",
    allele_alt = "T"
  )
}

###############################################################################
# Diploid data: three populations
###############################################################################
set.seed(42)
pop_sizes <- c(popA = 8, popB = 6, popC = 5)
n_loci <- 300
sim <- simulate_genotypes(pop_sizes, n_loci, ploidy = 2)
genotypes <- sim$genotypes
population <- sim$population
# loci where a whole population is missing, where a single individual is
# genotyped in a population, and that are monomorphic
genotypes[population == "popA", 10] <- NA
genotypes[population == "popB", 200] <- NA
genotypes[population == "popC", 50] <- c(1L, rep(NA, 4))
genotypes[, 100] <- 0L
genotypes[, 101] <- 2L
loci <- simulate_loci(n_loci)
ids <- paste0("ind", seq_along(population))

write_test_vcf(genotypes, 2, ids, loci, file.path(out_dir, "easysfs_test.vcf"))
pop_path <- file.path(out_dir, "easysfs_populations.txt")
utils::write.table(data.frame(ids, population),
  file = pop_path, sep = "\t",
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

# preview of the projections
preview <- run_easysfs(c(
  "-i", file.path(out_dir, "easysfs_test.vcf"), "-p", pop_path, "-a", "-f",
  "-y", "--preview"
))
writeLines(preview, file.path(out_dir, "easysfs_preview.txt"))

# spectra with projections smaller than the number of allele copies (16, 12
# and 10), folded and unfolded, with the populations in another order, and
# with the total length of the sequences (to set the monomorphic sites)
run_spectra <- function(out_name, extra_args) {
  tmp_out <- file.path(tempdir(), out_name)
  run_easysfs(c(
    "-i", file.path(out_dir, "easysfs_test.vcf"), "-p", pop_path, "-a",
    "-f", "-y", "-o", tmp_out, "--prefix", "test", extra_args
  ))
  # keep the dadi and fastsimcoal2 outputs
  dest <- file.path(out_dir, out_name)
  unlink(dest, recursive = TRUE)
  dir.create(dest)
  file.copy(file.path(tmp_out, c("dadi", "fastsimcoal2")), dest,
    recursive = TRUE
  )
}
run_spectra("folded", c("--proj", "10,8,6"))
run_spectra("unfolded", c("--proj", "10,8,6", "--unfolded"))
run_spectra("ordered", c("--proj", "6,10,8", "--order", "popC,popA,popB"))
run_spectra("total_length", c("--proj", "10,8,6", "--total-length", "100000"))

###############################################################################
# Haploid data: two populations
###############################################################################
set.seed(43)
pop_sizes_hap <- c(popA = 7, popB = 6)
sim_hap <- simulate_genotypes(pop_sizes_hap, 200, ploidy = 1)
loci_hap <- simulate_loci(200)
ids_hap <- paste0("hap", seq_along(sim_hap$population))
write_test_vcf(sim_hap$genotypes, 1, ids_hap, loci_hap,
  file = file.path(out_dir, "easysfs_haploid.vcf")
)
pop_path_hap <- file.path(out_dir, "easysfs_haploid_populations.txt")
utils::write.table(data.frame(ids_hap, sim_hap$population),
  file = pop_path_hap, sep = "\t",
  quote = FALSE, row.names = FALSE, col.names = FALSE
)
tmp_out <- file.path(tempdir(), "haploid")
run_easysfs(c(
  "-i", file.path(out_dir, "easysfs_haploid.vcf"), "-p", pop_path_hap,
  "-a", "-f", "-y", "-o", tmp_out, "--prefix", "test", "--ploidy", "1",
  "--proj", "5,4"
))
dest <- file.path(out_dir, "haploid")
unlink(dest, recursive = TRUE)
dir.create(dest)
file.copy(file.path(tmp_out, c("dadi", "fastsimcoal2")), dest,
  recursive = TRUE
)
