# delete the backing files of this file when it ends (see helper_cleanup.R)
local_tempdir_cleanup()

# tests of phased gen_tibbles (gen_tibble(phased = TRUE) and the functions
# that have to keep, or use, the phase of the genotypes)

phased_vcf <- system.file("extdata", "phased_example.vcf",
  package = "tidypopgen"
)

# the haplotypes of the biallelic loci of phased_example.vcf (two rows per
# individual), and the corresponding dosages
phased_vcf_haps <- rbind(
  c(0, 1, 0, 0, 1, 0, 1), # ind1
  c(0, 0, 1, 0, 1, 0, 0),
  c(0, 1, 0, 1, 0, 0, 0), # ind2
  c(1, 1, 0, 0, 1, 0, 1),
  c(1, 0, 0, 0, 1, 1, 0), # ind3
  c(0, 0, 0, 1, 0, 1, 0),
  c(1, 0, 1, 0, 0, 0, 1), # ind4
  c(1, 1, 0, 0, 0, 1, 1),
  c(0, NA, 1, 1, 0, NA, 1), # ind5
  c(1, NA, 1, 0, 1, NA, 0)
)
storage.mode(phased_vcf_haps) <- "integer"
dimnames(phased_vcf_haps) <- list(
  paste0(rep(paste0("ind", 1:5), each = 2), c("_1", "_2")),
  paste0("rs", c(1, 2, 4:8))
)

# random haplotypes (two rows per individual) for n_ind individuals and
# n_loci loci, with a proportion p_missing of missing genotypes
random_haplotypes <- function(n_ind, n_loci, p_missing = 0.05) {
  freq <- stats::runif(n_loci, 0.05, 0.95)
  haps <- matrix(
    stats::rbinom(2 * n_ind * n_loci, 1, rep(freq, each = 2 * n_ind)),
    nrow = 2 * n_ind
  )
  # both haplotypes of a missing genotype are missing
  missing <- matrix(stats::runif(n_ind * n_loci) < p_missing, nrow = n_ind)
  haps[missing[rep(seq_len(n_ind), each = 2), , drop = FALSE]] <- NA
  haps
}

random_loci <- function(n_loci) {
  data.frame(
    name = paste0("rs", seq_len(n_loci)),
    chromosome = rep(c("1", "2"), times = c(ceiling(n_loci / 2), n_loci %/% 2)),
    position = as.integer(seq_len(n_loci) * 100),
    genetic_dist = 0,
    allele_ref = "A",
    allele_alt = "T"
  )
}

# dosages (individuals x loci) from haplotypes
haps_to_dosages <- function(haps) {
  n_ind <- nrow(haps) / 2
  haps[seq(1, by = 2, length.out = n_ind), , drop = FALSE] +
    haps[seq(2, by = 2, length.out = n_ind), , drop = FALSE]
}

################################################################################
# reading phased vcf files
################################################################################

test_that("gen_tibble(phased = TRUE) reads a phased vcf", {
  x <- gen_tibble(phased_vcf,
    phased = TRUE,
    backingfile = tempfile("phased_"), quiet = TRUE
  )
  expect_true(is_phased(x))
  expect_true(is_phased(x$genotypes))
  expect_identical(.gt_get_fbm(x)$code256, code_phased)
  # the multiallelic locus (rs3) is dropped
  expect_identical(loci_names(x), paste0("rs", c(1, 2, 4:8)))
  expect_identical(show_haplotypes(x), phased_vcf_haps)
  # the genotypes are the dosages of the alternate allele
  expect_identical(
    show_genotypes(x),
    unname(haps_to_dosages(phased_vcf_haps) * 1)
  )
  # the raw bytes follow the phased encoding (1|0 is byte 4, missing 3)
  raw_bytes <- .gt_get_fbm(x)$copy(code = as.numeric(0:255))[, ]
  expect_identical(raw_bytes[, 1], c(0, 1, 4, 2, 1))
  expect_identical(raw_bytes[5, c(2, 6)], c(3, 3))
})

test_that("gen_tibble(phased = TRUE) reads gzipped vcf files", {
  vcf_gz <- tempfile(fileext = ".vcf.gz")
  con <- gzfile(vcf_gz, "w")
  writeLines(readLines(phased_vcf), con)
  close(con)
  x <- gen_tibble(vcf_gz,
    phased = TRUE, backingfile = tempfile(),
    quiet = TRUE
  )
  expect_true(is_phased(x))
  expect_identical(show_haplotypes(x), phased_vcf_haps)
})

test_that("gen_tibble(phased = TRUE) stores integer names and chromosomes", {
  names_gt <- gen_tibble(phased_vcf,
    phased = TRUE,
    backingfile = tempfile("phased_"), quiet = TRUE
  )
  for (n_cores in c(1, 2)) {
    int_gt <- gen_tibble(phased_vcf,
      phased = TRUE,
      backingfile = tempfile("phased_"), quiet = TRUE, n_cores = n_cores,
      names_as_int = TRUE, chromosomes_as_int = TRUE
    )
    expect_true(is_phased(int_gt))
    expect_identical(show_loci(int_gt)$name, 1:7)
    expect_true(is.integer(show_loci(int_gt)$chromosome))
    expect_identical(
      unname(show_haplotypes(int_gt)),
      unname(phased_vcf_haps)
    )
    # the sidecars give back the names and chromosomes
    expect_identical(
      show_loci(gt_restore_loci(int_gt)),
      show_loci(names_gt)
    )
  }
})

test_that("gen_tibble ignores the phase of a phased vcf by default", {
  x <- gen_tibble(phased_vcf,
    backingfile = tempfile("unphased_"), quiet = TRUE
  )
  expect_false(is_phased(x))
  expect_identical(.gt_get_fbm(x)$code256, bigsnpr::CODE_012)
  expect_identical(
    show_genotypes(x),
    unname(haps_to_dosages(phased_vcf_haps) * 1)
  )
  expect_error(show_haplotypes(x), "not phased")
})

# write a small vcf with the given genotypes (one string per individual, for
# each locus)
write_test_vcf <- function(genotypes, file = tempfile(fileext = ".vcf")) {
  header <- c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", paste0("s", seq_along(genotypes[[1]]))
    ), collapse = "\t")
  )
  body <- vapply(seq_along(genotypes), function(i) {
    paste(c(
      "1", 100 * i, paste0("rs", i), "A", "T", ".", ".", ".", "GT",
      genotypes[[i]]
    ), collapse = "\t")
  }, character(1))
  writeLines(c(header, body), file)
  file
}

# read a vcf with phased = TRUE
read_phased <- function(vcf, ...) {
  gen_tibble(vcf, phased = TRUE, backingfile = tempfile(), quiet = TRUE, ...)
}

test_that("gen_tibble(phased = TRUE) reads missing genotypes", {
  vcf <- write_test_vcf(list(
    c(".", "0|1", "1|1", ".|."),
    c("0|0", "1|0", "0|.", ".|1"),
    c("1|1", "0|0", "1|0", "0|0")
  ))
  haps <- unname(show_haplotypes(read_phased(vcf)))
  expected <- matrix(c(
    NA, 0L, 1L,
    NA, 0L, 1L,
    0L, 1L, 0L,
    1L, 0L, 0L,
    1L, NA, 1L,
    1L, NA, 0L,
    NA, NA, 0L,
    NA, NA, 0L
  ), ncol = 3, byrow = TRUE)
  expect_identical(haps, expected)
})

test_that("gen_tibble(phased = TRUE) refuses any unphased genotype", {
  # heterozygotes, homozygotes and missing genotypes
  for (unphased in c("0/1", "1/0", "0/0", "1/1", "./.")) {
    vcf <- write_test_vcf(list(
      c("0|0", "0|1", "1|1"),
      c("0|1", "1|0", unphased)
    ))
    expect_error(
      read_phased(vcf),
      "individual s3 at locus 1:200 is unphased (with '/')",
      fixed = TRUE
    )
  }
  # the first unphased genotype of the file is named (also when the lines
  # are read in parallel)
  vcf <- write_test_vcf(list(
    c("0|0", "0|1", "0/0"),
    c("1/1", "1|0", "0|1")
  ))
  for (n_cores in c(1, 2)) {
    expect_error(
      read_phased(vcf, n_cores = n_cores),
      "individual s3 at locus 1:100 is unphased"
    )
  }
  # the same file can be read (unphased) by default
  x <- gen_tibble(vcf, backingfile = tempfile(), quiet = TRUE)
  expect_identical(show_genotypes(x), rbind(c(0, 2), c(1, 1), c(0, 1)))
  # no backing file is left behind by the failed attempts
  bk <- tempfile()
  expect_error(gen_tibble(vcf, phased = TRUE, backingfile = bk, quiet = TRUE))
  expect_false(file.exists(paste0(bk, ".bk")))
})

test_that("gen_tibble(phased = TRUE) refuses genotypes that are not diploid", {
  # haploid individual (found when estimating ploidy)
  vcf <- write_test_vcf(list(c("0|0", "1"), c("0|1", "0")))
  expect_error(read_phased(vcf), "only be read for diploid individuals")
  # a haploid genotype after the first locus
  vcf <- write_test_vcf(list(c("0|0", "1|1"), c("0|1", "1")))
  expect_error(
    read_phased(vcf),
    "individual s2 at locus 1:200 is not diploid"
  )
  # a line with fewer genotypes than individuals
  vcf <- write_test_vcf(list(c("0|0", "1|1", "0|1"), c("0|1", "1|0", "1|1")))
  vcf_lines <- readLines(vcf)
  vcf_lines[5] <- sub("\t1|1", "", vcf_lines[5], fixed = TRUE)
  writeLines(vcf_lines, vcf)
  expect_error(
    read_phased(vcf),
    "locus 1:200 has fewer genotypes than individuals"
  )
  # a triploid genotype after the first locus
  vcf <- write_test_vcf(list(c("0|0", "1|1"), c("0|1|1", "1|0")))
  expect_error(
    read_phased(vcf),
    "individual s1 at locus 1:200 is not diploid"
  )
})

test_that("phased = TRUE is only available for vcf files with parser cpp", {
  bed_file <- system.file("extdata", "pop_a.bed", package = "tidypopgen")
  expect_error(
    gen_tibble(bed_file, phased = TRUE, quiet = TRUE),
    "only available for VCF files read with parser = 'cpp'"
  )
  expect_error(
    gen_tibble(phased_vcf, phased = TRUE, parser = "vcfR", quiet = TRUE),
    "only available for VCF files read with parser = 'cpp'"
  )
})

################################################################################
# matrix of haplotypes
################################################################################

test_that("gen_tibble(phased = TRUE) creates a gen_tibble from haplotypes", {
  haps <- matrix(c(
    0, 1, 1,
    1, 1, 0,
    0, 0, NA,
    0, 1, 1
  ), nrow = 4, byrow = TRUE)
  loci <- data.frame(
    name = paste0("rs", 1:3), chromosome = 1, position = c(10, 20, 30),
    genetic_dist = 0, allele_ref = "A", allele_alt = "T"
  )
  indiv <- data.frame(id = c("a", "b"), population = "pop1")
  x <- gen_tibble(haps,
    indiv_meta = indiv, loci = loci, phased = TRUE,
    quiet = TRUE
  )
  expect_true(is_phased(x))
  expected <- haps
  expected[4, 3] <- NA # the genotype is missing
  storage.mode(expected) <- "integer"
  dimnames(expected) <- list(c("a_1", "a_2", "b_1", "b_2"), loci$name)
  expect_identical(show_haplotypes(x), expected)
  expect_identical(show_genotypes(x), rbind(c(1, 2, 1), c(0, 1, NA)))
  # subsets of individuals and loci
  expect_identical(
    show_haplotypes(x, indiv_indices = 2, loci_indices = c(3, 1)),
    expected[3:4, c(3, 1)]
  )
  expect_identical(
    unname(show_haplotypes(x$genotypes)),
    unname(expected)
  )
})

test_that("gen_tibble(phased = TRUE) validates the matrix of haplotypes", {
  loci <- data.frame(
    name = c("rs1", "rs2"), chromosome = 1, position = c(10, 20),
    genetic_dist = 0, allele_ref = "A", allele_alt = "T"
  )
  indiv <- data.frame(id = c("a", "b"))
  phased_matrix <- function(x, ...) {
    gen_tibble(x, indiv, loci, phased = TRUE, ...)
  }
  expect_error(phased_matrix(matrix(0, 3, 2)), "two rows")
  expect_error(phased_matrix(matrix(2, 4, 2)), "0, 1 or NA")
  expect_error(phased_matrix(matrix("0", 4, 2)), "numeric matrix")
  expect_error(phased_matrix(matrix(0, 4, 2), ploidy = 1), "diploid")
  # the metadata are validated as for unphased data
  expect_error(
    phased_matrix(matrix(0, 6, 2)),
    "mismatch between the number of individuals"
  )
})

################################################################################
# standard functions on phased data
################################################################################

set.seed(123)
n_ind <- 40
n_loci <- 300
rand_haps <- random_haplotypes(n_ind, n_loci, p_missing = 0.05)
rand_loci <- random_loci(n_loci)
rand_indiv <- data.frame(
  id = paste0("ind", seq_len(n_ind)),
  population = rep(c("pop1", "pop2", "pop3"), length.out = n_ind)
)
phased_gt <- gen_tibble(rand_haps,
  phased = TRUE,
  indiv_meta = rand_indiv,
  loci = rand_loci, backingfile = tempfile("rand_phased_"), quiet = TRUE
)
unphased_gt <- gen_tibble(haps_to_dosages(rand_haps),
  indiv_meta = rand_indiv,
  loci = rand_loci, backingfile = tempfile("rand_unphased_"), quiet = TRUE
)

test_that("phased and unphased gen_tibbles have the same genotypes", {
  expect_true(is_phased(phased_gt))
  expect_false(is_phased(unphased_gt))
  # there are heterozygotes in both orientations
  raw_bytes <- .gt_get_fbm(phased_gt)$copy(code = as.numeric(0:255))[, ]
  expect_true(all(0:4 %in% raw_bytes))
  expect_identical(show_genotypes(phased_gt), show_genotypes(unphased_gt))
})

test_that("loci and individual statistics are the same with phased data", {
  expect_equal(loci_alt_freq(phased_gt), loci_alt_freq(unphased_gt))
  expect_equal(loci_maf(phased_gt), loci_maf(unphased_gt))
  expect_equal(loci_missingness(phased_gt), loci_missingness(unphased_gt))
  expect_equal(loci_hwe(phased_gt), loci_hwe(unphased_gt))
  expect_equal(loci_pi(phased_gt), loci_pi(unphased_gt))
  expect_equal(indiv_het_obs(phased_gt), indiv_het_obs(unphased_gt))
  expect_equal(indiv_missingness(phased_gt), indiv_missingness(unphased_gt))
  expect_equal(qc_report_loci(phased_gt), qc_report_loci(unphased_gt))
  expect_equal(qc_report_indiv(phased_gt), qc_report_indiv(unphased_gt))
  # grouped
  phased_grp <- phased_gt %>% dplyr::group_by(population)
  unphased_grp <- unphased_gt %>% dplyr::group_by(population)
  expect_equal(loci_alt_freq(phased_grp), loci_alt_freq(unphased_grp))
  expect_equal(
    loci_alt_freq(phased_grp, type = "tidy"),
    loci_alt_freq(unphased_grp, type = "tidy")
  )
  expect_equal(loci_missingness(phased_grp), loci_missingness(unphased_grp))
  expect_equal(loci_hwe(phased_grp), loci_hwe(unphased_grp))
  expect_equal(loci_pi(phased_grp), loci_pi(unphased_grp))
})

test_that("population statistics are the same with phased data", {
  phased_grp <- phased_gt %>% dplyr::group_by(population)
  unphased_grp <- unphased_gt %>% dplyr::group_by(population)
  expect_equal(pop_het_obs(phased_grp), pop_het_obs(unphased_grp))
  expect_equal(pop_het_exp(phased_grp), pop_het_exp(unphased_grp))
  expect_equal(pop_fis(phased_grp), pop_fis(unphased_grp))
  expect_equal(pop_global_stats(phased_grp), pop_global_stats(unphased_grp))
  expect_equal(pop_fst(phased_grp), pop_fst(unphased_grp))
  expect_equal(
    pairwise_pop_fst(phased_grp),
    pairwise_pop_fst(unphased_grp)
  )
  expect_equal(
    pairwise_pop_fst(phased_grp, method = "WC84"),
    pairwise_pop_fst(unphased_grp, method = "WC84")
  )
  expect_equal(pop_tajimas_d(phased_grp), pop_tajimas_d(unphased_grp))
  expect_equal(nwise_pop_pbs(phased_grp), nwise_pop_pbs(unphased_grp))
  expect_equal(
    pop_sfs(phased_grp, n_proj = c(10, 10, 10)),
    pop_sfs(unphased_grp, n_proj = c(10, 10, 10))
  )
  # windowed and jackknife versions
  expect_equal(
    windows_pairwise_pop_fst(phased_grp, window_size = 50, step_size = 25),
    windows_pairwise_pop_fst(unphased_grp, window_size = 50, step_size = 25)
  )
  expect_equal(
    windows_nwise_pop_pbs(phased_grp, window_size = 50, step_size = 25),
    windows_nwise_pop_pbs(unphased_grp, window_size = 50, step_size = 25)
  )
  expect_equal(
    windows_pop_tajimas_d(phased_grp, window_size = 50, step_size = 25),
    windows_pop_tajimas_d(unphased_grp, window_size = 50, step_size = 25)
  )
  expect_equal(
    jackknife_pairwise_pop_fst(phased_grp, block_size = 50),
    jackknife_pairwise_pop_fst(unphased_grp, block_size = 50)
  )
  expect_equal(
    jackknife_pop_tajimas_d(phased_grp, block_size = 50),
    jackknife_pop_tajimas_d(unphased_grp, block_size = 50)
  )
})

test_that("tibbles that lost the gen_tbl class keep the phase", {
  # e.g. from as_tibble(), or pick() and reframe(), whose class the generics
  # restore from the genotypes column (see restore_genotypes_tbl())
  plain <- tibble::as_tibble(phased_gt)
  expect_false(inherits(plain, "gen_tbl"))
  expect_true(is_phased(plain))
  expect_identical(show_haplotypes(plain), show_haplotypes(phased_gt))
  expect_equal(loci_alt_freq(plain), loci_alt_freq(unphased_gt))
  reframe_stats <- function(x) {
    x %>%
      dplyr::group_by(population) %>%
      dplyr::reframe(
        phased = is_phased(dplyr::pick(dplyr::everything())),
        hwe = loci_hwe(dplyr::pick(dplyr::everything()))
      )
  }
  phased_stats <- reframe_stats(phased_gt)
  unphased_stats <- reframe_stats(unphased_gt)
  expect_true(all(phased_stats$phased))
  expect_false(any(unphased_stats$phased))
  expect_equal(phased_stats$hwe, unphased_stats$hwe)
})

test_that("pairwise statistics of individuals are the same with phased data", {
  expect_equal(pairwise_ibs(phased_gt), pairwise_ibs(unphased_gt))
  expect_equal(
    pairwise_ibs(phased_gt, missing_correction = "allele_freq"),
    pairwise_ibs(unphased_gt, missing_correction = "allele_freq")
  )
  ibs_counts_phased <- pairwise_ibs(phased_gt, type = "raw_counts")
  ibs_counts_unphased <- pairwise_ibs(unphased_gt, type = "raw_counts")
  expect_identical(ibs_counts_phased$ibs[], ibs_counts_unphased$ibs[])
  expect_identical(ibs_counts_phased$valid_n[], ibs_counts_unphased$valid_n[])
  expect_equal(pairwise_king(phased_gt), pairwise_king(unphased_gt))
  expect_equal(
    pairwise_allele_sharing(phased_gt),
    pairwise_allele_sharing(unphased_gt)
  )
  expect_equal(pairwise_grm(phased_gt), pairwise_grm(unphased_gt))
  expect_equal(indiv_inbreeding(phased_gt), indiv_inbreeding(unphased_gt))
  # the same as computing the statistics from the dosages in R
  dosages <- show_genotypes(unphased_gt)
  i <- 3
  j <- 7
  valid <- !is.na(dosages[i, ]) & !is.na(dosages[j, ])
  ibs_ij <- sum(2 - abs(dosages[i, valid] - dosages[j, valid])) /
    (2 * sum(valid))
  ibs <- pairwise_ibs(phased_gt, as_matrix = TRUE)
  expect_equal(unname(ibs[i, j]), ibs_ij)
})

test_that("PCA and LD clumping are the same with phased data", {
  # PCA needs complete data (without any missing genotype in the backing
  # file)
  complete_haps <- random_haplotypes(n_ind, 100, p_missing = 0)
  phased_complete <- gen_tibble(complete_haps,
    phased = TRUE,
    indiv_meta = rand_indiv, loci = random_loci(100),
    backingfile = tempfile(), quiet = TRUE
  )
  unphased_complete <- gen_tibble(haps_to_dosages(complete_haps),
    indiv_meta = rand_indiv, loci = random_loci(100),
    backingfile = tempfile(), quiet = TRUE
  )
  expect_true(is_phased(phased_complete))
  pca_phased <- gt_pca_partialSVD(phased_complete, k = 3)
  pca_unphased <- gt_pca_partialSVD(unphased_complete, k = 3)
  expect_equal(abs(pca_phased$u), abs(pca_unphased$u), tolerance = 1e-6)
  expect_equal(pca_phased$d, pca_unphased$d, tolerance = 1e-6)
  # ld clumping
  expect_identical(
    loci_ld_clump(phased_complete, thr_r2 = 0.2),
    loci_ld_clump(unphased_complete, thr_r2 = 0.2)
  )
})

test_that("exports of phased data are the same as of unphased data", {
  # plink files have no phase: the same bed file is written
  bed_phased <- gt_as_plink(phased_gt, file = tempfile(), type = "bed")
  bed_unphased <- gt_as_plink(unphased_gt, file = tempfile(), type = "bed")
  expect_identical(
    unname(tools::md5sum(bed_phased)),
    unname(tools::md5sum(bed_unphased))
  )
  raw_phased <- gt_as_plink(phased_gt, file = tempfile(), type = "raw")
  raw_unphased <- gt_as_plink(unphased_gt, file = tempfile(), type = "raw")
  expect_identical(readLines(raw_phased), readLines(raw_unphased))
  # genind
  skip_if_not_installed("adegenet")
  expect_equal(
    gt_as_genind(phased_gt)@tab,
    gt_as_genind(unphased_gt)@tab
  )
})

test_that("phased gen_tibbles can not be imputed", {
  expect_error(gt_impute_simple(phased_gt), "phased genotypes can not be")
  skip_if_not_installed("xgboost")
  expect_error(gt_impute_xgboost(phased_gt), "phased genotypes can not be")
})

################################################################################
# keeping the phase
################################################################################

test_that("subsetting, saving and loading keep phased data", {
  # subset of individuals and loci (in a different order)
  x_sub <- phased_gt %>%
    dplyr::filter(population == "pop2") %>%
    select_loci(c(10, 3, 200))
  expect_true(is_phased(x_sub))
  rows <- which(rand_indiv$population == "pop2")
  hap_rows <- as.vector(rbind(2 * rows - 1, 2 * rows))
  expect_identical(
    unname(show_haplotypes(x_sub)),
    unname(show_haplotypes(phased_gt)[hap_rows, c(10, 3, 200)])
  )
  # save and load
  gt_file <- gt_save(x_sub, quiet = TRUE)[1]
  x_loaded <- gt_load(gt_file)
  expect_true(is_phased(x_loaded))
  expect_identical(show_haplotypes(x_loaded), show_haplotypes(x_sub))
  # a new backing file with only the individuals and loci in the gen_tibble
  x_new <- gt_update_backingfile(x_sub,
    backingfile = tempfile(),
    quiet = TRUE
  )
  expect_true(is_phased(x_new))
  expect_identical(show_haplotypes(x_new), show_haplotypes(x_sub))
  # reordering loci
  x_shuffled <- phased_gt %>% select_loci(sample(count_loci(phased_gt)))
  x_ordered <- gt_order_loci(x_shuffled,
    use_current_table = FALSE,
    quiet = TRUE
  )
  expect_true(is_phased(x_ordered))
  expect_identical(show_haplotypes(x_ordered), show_haplotypes(phased_gt))
})

test_that("gt_as_vcf writes phased genotypes, which can be read back", {
  vcf <- gt_as_vcf(phased_gt, file = tempfile(fileext = ".vcf"))
  vcf_lines <- readLines(vcf)
  genotypes <- unlist(strsplit(
    vcf_lines[!startsWith(vcf_lines, "#")],
    "\t"
  ))
  # the genotypes are phased (missing genotypes are .|.)
  gts <- genotypes[grepl("^[0-9.][/|][0-9.]$", genotypes)]
  expect_true(all(gts %in% c("0|0", "0|1", "1|0", "1|1", ".|.")))
  expect_true(all(c("0|1", "1|0", ".|.") %in% gts))
  # round trip
  x <- gen_tibble(vcf, phased = TRUE, backingfile = tempfile(), quiet = TRUE)
  expect_true(is_phased(x))
  expect_identical(show_haplotypes(x), show_haplotypes(phased_gt))
  expect_identical(show_loci(x)$name, show_loci(phased_gt)$name)
  expect_identical(show_loci(x)$position, show_loci(phased_gt)$position)
  # the same vcf read without phased = TRUE gives the unphased genotypes
  x_unphased <- gen_tibble(vcf, backingfile = tempfile(), quiet = TRUE)
  expect_false(is_phased(x_unphased))
  expect_identical(
    show_genotypes(x_unphased),
    show_genotypes(unphased_gt)
  )
  # unphased gen_tibbles are still written as unphased
  vcf_unphased <- gt_as_vcf(unphased_gt, file = tempfile(fileext = ".vcf"))
  vcf_lines <- readLines(vcf_unphased)
  vcf_body <- vcf_lines[!startsWith(vcf_lines, "#")]
  expect_false(any(grepl("|", vcf_body, fixed = TRUE)))
})

# two phased gen_tibbles with the same loci, and a version of the second one
# in which the reference and alternate alleles of some loci are swapped (its
# haplotypes then carry the other allele)
set.seed(42)
rb_haps_a <- random_haplotypes(6, 20, p_missing = 0.1)
rb_haps_b <- random_haplotypes(5, 20, p_missing = 0.1)
rb_loci <- random_loci(20)
rb_indiv_a <- data.frame(id = paste0("a", 1:6), population = "a")
rb_indiv_b <- data.frame(id = paste0("b", 1:5), population = "b")
rb_swapped <- c(2, 5, 6, 11, 17)
rb_loci_b_swapped <- rb_loci
rb_loci_b_swapped$allele_ref[rb_swapped] <- rb_loci$allele_alt[rb_swapped]
rb_loci_b_swapped$allele_alt[rb_swapped] <- rb_loci$allele_ref[rb_swapped]
rb_haps_b_swapped <- rb_haps_b
rb_haps_b_swapped[, rb_swapped] <- 1 - rb_haps_b[, rb_swapped]
rb_expected <- rbind(rb_haps_a, rb_haps_b)
storage.mode(rb_expected) <- "integer"

test_that("rbind keeps the phase of two phased gen_tibbles", {
  x_a <- gen_tibble(rb_haps_a,
    phased = TRUE,
    indiv_meta = rb_indiv_a, loci = rb_loci,
    backingfile = tempfile(), quiet = TRUE
  )
  x_b <- gen_tibble(rb_haps_b,
    phased = TRUE,
    indiv_meta = rb_indiv_b, loci = rb_loci,
    backingfile = tempfile(), quiet = TRUE
  )
  merged <- rbind(x_a, x_b, quiet = TRUE, backingfile = tempfile())
  expect_true(is_phased(merged))
  expect_identical(unname(show_haplotypes(merged)), unname(rb_expected))
  expect_identical(show_loci(merged)$allele_ref, rb_loci$allele_ref)
  expect_identical(show_loci(merged)$allele_alt, rb_loci$allele_alt)
  # only the loci in common are kept
  x_b_sub <- x_b %>% select_loci(c(1:10, 15))
  merged_sub <- rbind(x_a, x_b_sub, quiet = TRUE, backingfile = tempfile())
  expect_true(is_phased(merged_sub))
  expect_identical(
    unname(show_haplotypes(merged_sub)),
    unname(rb_expected[, c(1:10, 15)])
  )
})

test_that("rbind does not swap or flip the alleles of phased gen_tibbles", {
  x_a <- gen_tibble(rb_haps_a,
    phased = TRUE,
    indiv_meta = rb_indiv_a, loci = rb_loci,
    backingfile = tempfile(), quiet = TRUE
  )
  x_b <- gen_tibble(rb_haps_b,
    phased = TRUE,
    indiv_meta = rb_indiv_b, loci = rb_loci,
    backingfile = tempfile(), quiet = TRUE
  )
  x_b_swapped <- gen_tibble(rb_haps_b_swapped,
    phased = TRUE,
    indiv_meta = rb_indiv_b, loci = rb_loci_b_swapped,
    backingfile = tempfile(), quiet = TRUE
  )
  # swapped alleles (no backing file is created)
  bk <- tempfile()
  expect_error(
    rbind(x_a, x_b_swapped, quiet = TRUE, backingfile = bk),
    "alleles of 5 loci are swapped between the two phased gen_tibbles"
  )
  expect_false(file.exists(paste0(bk, ".bk")))
  # the loci can be found with rbind_dry_run(), and removed
  report <- rbind_dry_run(x_a, x_b_swapped, quiet = TRUE)
  to_swap <- which(report$target$to_swap)
  expect_equal(to_swap, rb_swapped)
  merged <- rbind(x_a, x_b_swapped %>% select_loci(-dplyr::all_of(to_swap)),
    quiet = TRUE, backingfile = tempfile()
  )
  expect_true(is_phased(merged))
  expect_identical(
    unname(show_haplotypes(merged)),
    unname(rb_expected[, -rb_swapped])
  )
  # strand flipping
  expect_error(
    rbind(x_a, x_b, flip_strand = TRUE, quiet = TRUE, backingfile = bk),
    "strands can not be flipped when merging phased gen_tibbles"
  )
  expect_false(file.exists(paste0(bk, ".bk")))
})

test_that("rbind of phased and unphased gen_tibbles discards the phase", {
  x_a <- gen_tibble(rb_haps_a,
    phased = TRUE,
    indiv_meta = rb_indiv_a, loci = rb_loci,
    backingfile = tempfile(), quiet = TRUE
  )
  x_b_swapped <- gen_tibble(rb_haps_b_swapped,
    phased = TRUE,
    indiv_meta = rb_indiv_b, loci = rb_loci_b_swapped,
    backingfile = tempfile(), quiet = TRUE
  )
  x_a_unphased <- gen_tibble(haps_to_dosages(rb_haps_a),
    indiv_meta = rb_indiv_a, loci = rb_loci,
    backingfile = tempfile(), quiet = TRUE
  )
  x_b_unphased <- gen_tibble(haps_to_dosages(rb_haps_b_swapped),
    indiv_meta = rb_indiv_b, loci = rb_loci_b_swapped,
    backingfile = tempfile(), quiet = TRUE
  )
  # (the merged data are unphased, so alleles are swapped as for any
  # unphased data)
  expected_dosages <- unname(haps_to_dosages(rb_expected) * 1)
  # (the report of rbind is printed, so we capture it)
  utils::capture.output(expect_message(
    merged_ab <- rbind(x_a, x_b_unphased, backingfile = tempfile()),
    "only the ref gen_tibble is phased"
  ))
  expect_false(is_phased(merged_ab))
  expect_identical(.gt_get_fbm(merged_ab)$code256, bigsnpr::CODE_012)
  expect_identical(show_genotypes(merged_ab), expected_dosages)
  # the raw bytes are unphased (no byte 4)
  raw_bytes <- .gt_get_fbm(merged_ab)$copy(code = as.numeric(0:255))[, ]
  expect_true(all(raw_bytes %in% 0:3))
  utils::capture.output(expect_message(
    merged_ba <- rbind(x_a_unphased, x_b_swapped, backingfile = tempfile()),
    "only the target gen_tibble is phased"
  ))
  expect_false(is_phased(merged_ba))
  expect_identical(show_genotypes(merged_ba), expected_dosages)
  raw_bytes <- .gt_get_fbm(merged_ba)$copy(code = as.numeric(0:255))[, ]
  expect_true(all(raw_bytes %in% 0:3))
  # unphased data are merged as before
  merged_unphased <- rbind(x_a_unphased, x_b_unphased,
    quiet = TRUE,
    backingfile = tempfile()
  )
  expect_false(is_phased(merged_unphased))
  expect_identical(show_genotypes(merged_unphased), expected_dosages)
})

test_that("rbind_byte_luts recodes phased bytes", {
  # phased data that stay phased are copied as they are (and never swapped)
  luts <- rbind_byte_luts(phased = TRUE, keep_phase = TRUE)
  expect_identical(luts$copy, 0:255)
  expect_identical(luts$swap, 0:255)
  # without keeping the phase, 1|0 becomes 0|1 (the unphased heterozygote),
  # also after swapping the alleles of both haplotypes
  luts <- rbind_byte_luts(phased = TRUE, keep_phase = FALSE)
  expect_identical(luts$copy[1:6], c(0L, 1L, 2L, 3L, 1L, 5L))
  expect_identical(luts$swap[1:6], c(2L, 1L, 0L, 3L, 1L, 5L))
  # unphased (and imputed) data
  luts <- rbind_byte_luts(phased = FALSE, keep_phase = FALSE)
  expect_identical(luts$copy, 0:255)
  expect_identical(luts$swap[1:8], c(2L, 1L, 0L, 3L, 6L, 5L, 4L, 7L))
})

################################################################################
# merging and concatenating phased vcf files
################################################################################

merge_loci <- data.frame(
  chrom = "1", pos = c(100, 200, 300), id = c("rs1", "rs2", "rs3"),
  ref = "A", alt = "T"
)
merge_geno_a <- rbind(c("0|1", "1|1"), c("1|0", ".|."), c("0|0", "0|1"))
merge_geno_b <- rbind(c("1|0", "0|0"), c("0|1", "1|0"), c("1|1", "0|1"))
merge_vcf_a <- write_merge_vcf(c("a1", "a2"), merge_loci, merge_geno_a)
merge_vcf_b <- write_merge_vcf(c("b1", "b2"), merge_loci, merge_geno_b)

test_that("gt_merge_vcfs reads phased vcf files", {
  merged_gt <- gt_merge_vcfs(c(merge_vcf_a, merge_vcf_b),
    phased = TRUE, backingfile = tempfile(), quiet = TRUE
  )
  # the same data in a single vcf
  single_gt <- gen_tibble(
    write_merge_vcf(
      c("a1", "a2", "b1", "b2"), merge_loci,
      cbind(merge_geno_a, merge_geno_b)
    ),
    phased = TRUE,
    names_as_int = TRUE, backingfile = tempfile(), quiet = TRUE
  )
  expect_true(is_phased(merged_gt))
  expect_identical(show_haplotypes(merged_gt), show_haplotypes(single_gt))
  # by default, the phase is ignored
  unphased_gt <- gt_merge_vcfs(c(merge_vcf_a, merge_vcf_b),
    backingfile = tempfile(), quiet = TRUE
  )
  expect_false(is_phased(unphased_gt))
  # unphased genotypes are an error (leaving no backing file)
  vcf_unphased <- write_merge_vcf(
    c("b1", "b2"), merge_loci,
    sub("|", "/", merge_geno_b, fixed = TRUE)
  )
  bk <- tempfile()
  expect_error(
    gt_merge_vcfs(c(merge_vcf_a, vcf_unphased),
      phased = TRUE, backingfile = bk, quiet = TRUE
    ),
    "individual b1 at locus 1:100 of the VCF"
  )
  expect_false(file.exists(paste0(bk, ".bk")))
  # the alleles of phased genotypes are never swapped
  expect_error(
    gt_merge_vcfs(c(merge_vcf_a, merge_vcf_b),
      phased = TRUE, ref_mismatch = "swap", quiet = TRUE
    ),
    "never swapped"
  )
})

test_that("gt_concat_vcfs reads phased vcf files", {
  # the samples are in a different order in the second file
  vcf_1 <- write_merge_vcf(
    c("a1", "a2"), merge_loci[1:2, ],
    merge_geno_a[1:2, ]
  )
  vcf_2 <- write_merge_vcf(
    c("a2", "a1"), merge_loci[3, ],
    merge_geno_a[3, 2:1, drop = FALSE]
  )
  concat_gt <- gt_concat_vcfs(c(vcf_1, vcf_2),
    phased = TRUE, backingfile = tempfile(), quiet = TRUE
  )
  single_gt <- gen_tibble(merge_vcf_a,
    phased = TRUE,
    names_as_int = TRUE, backingfile = tempfile(), quiet = TRUE
  )
  expect_true(is_phased(concat_gt))
  expect_identical(show_haplotypes(concat_gt), show_haplotypes(single_gt))
})

test_that("is_phased checks its input", {
  expect_error(is_phased(1:3), "should be a gen_tibble")
  expect_error(is_phased(data.frame(a = 1)), "not a gen_tibble")
  expect_identical(code_phased[1:6], c(0, 1, 2, NA, 1, NA))
  expect_identical(
    unname(phased_bytes[1:5, ]),
    matrix(c(0L, 0L, 1L, NA, 1L, 0L, 1L, 1L, NA, 0L), ncol = 2)
  )
})
