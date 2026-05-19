# limit number of threads for tests
data.table::setDTthreads(2)
if (rlang::is_installed("RhpcBLASctl")) {
  RhpcBLASctl::blas_set_num_threads(2)
  RhpcBLASctl::omp_set_num_threads(2)
}

test_that("autoplot when one loci is missing for all indivs", {
  # Create example
  example_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d", "e"),
    population = c("pop1", "pop1", "pop2", "pop2", "pop2")
  )
  example_genotypes <- rbind(
    c(1, NA, 0, 1, 1, 0),
    c(2, NA, 0, 0, NA, 0),
    c(1, NA, 0, 0, 1, 1),
    c(0, NA, 0, 1, 2, 1),
    c(1, NA, NA, 2, 1, 0)
  )
  example_loci <- data.frame(
    name = c("rs1", "rs2", "rs3", "rs4", "x1", "x2"),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
    genetic_dist = c(0, 0, 0, 0, 0, 0),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  example_gt <- gen_tibble(example_genotypes,
    indiv_meta = example_indiv_meta,
    loci = example_loci,
    backingfile = tempfile(),
    quiet = TRUE
  )
  example_gt <- example_gt %>% group_by(population)

  ex_loci_report <- example_gt %>%
    qc_report_loci()

  expect_message(
    {
      plt <- autoplot(ex_loci_report, type = "overview")
    },
    "One or more loci are missing for every individual"
  )

  # check the intersection counts are correct
  # 3 SNPs pass, 1 SNP is over missingness threshold (rs4),
  # and 1 SNP is monomorphic and over missingngess threshold (rs3)
  expect_equal(plt[[1]]$data$count, c(3, 1, 1))

  # set a very high missingness threshold, so that no loci are filtered out
  # check that we still get the message
  expect_message(
    {
      plt2 <- autoplot(ex_loci_report, type = "overview", miss_threshold = 0.5)
    },
    "One or more loci are missing for every individual"
  )
  # check the intersection counts are correct
  # because now 4 SNPs pass all checks and 1 SNP is monomorphic (rs3)
  expect_equal(
    plt2[[1]]$data,
    tibble(
      Missing = c(TRUE, TRUE),
      MAF = c(TRUE, FALSE),
      HWE = c(TRUE, TRUE),
      count = c(4, 1),
      intersection_id = c(1, 2)
    )
  )
})

test_that("loci report autoplot with all SNPs under MAF threshold", {
  example_indiv_meta <- data.frame(
    id = c("a", "b", "c", "d", "e"),
    population = c("pop1", "pop1", "pop1", "pop1", "pop1")
  )

  # All SNPs have a MAF of 0.1, no missingness, and hwe p = 0.5
  example_genotypes <- rbind(
    c(2, 0, 0, 0, 1, 0),
    c(2, 0, 0, 0, 2, 0),
    c(1, 1, 0, 0, 2, 0),
    c(2, 0, 0, 1, 2, 0),
    c(2, 0, 1, 0, 2, 1)
  )
  example_loci <- data.frame(
    name = c("rs1", "rs2", "rs3", "rs4", "x1", "x2"),
    chromosome = c(1, 1, 1, 1, 2, 2),
    position = c(3, 5, 65, 343, 23, 456),
    genetic_dist = c(0, 0, 0, 0, 0, 0),
    allele_ref = c("A", "T", "C", "G", "C", "T"),
    allele_alt = c("T", "C", NA, "C", "G", "A")
  )
  example_gt <- gen_tibble(example_genotypes,
    indiv_meta = example_indiv_meta,
    loci = example_loci,
    backingfile = tempfile(),
    quiet = TRUE
  )

  example_gt <- example_gt %>% group_by(population)

  ex_loci_report <- example_gt %>%
    qc_report_loci()

  # If no SNPs pass, the set size bar for that filter remains empty
  plt2 <- autoplot(ex_loci_report, type = "overview", maf_threshold = 0.1)
  # Expect that MAF is FALSE and 6 loci pass the other thresholds
  expect_equal(
    plt2[[1]]$data,
    tibble(
      Missing = TRUE,
      MAF = FALSE,
      HWE = TRUE,
      count = c(6),
      intersection_id = c(1)
    )
  )
})

test_that("reordering", {
  bed_path <- system.file(
    "extdata/related/families.bed",
    package = "tidypopgen"
  )
  families_bigsnp_path <- bigsnpr::snp_readBed(
    bed_path,
    backingfile = tempfile()
  )
  families <- gen_tibble(
    families_bigsnp_path,
    quiet = TRUE,
    valid_alleles = c("1", "2")
  )
  hwe <- loci_hwe(families)
  hwe_loci <- which(hwe < 0.05)
  families_sub <- families %>% select_loci(all_of(hwe_loci))

  expect_message(
    report <- qc_report_loci(families_sub),
    "This gen_tibble is not grouped"
  )

  # with default values
  total <- count_loci(families_sub)
  n_miss <- length(which(report$missingness < 0.01)) # 19 pass missingness
  n_maf <- length(which(report$maf > 0.05)) # all 30 pass maf
  n_hwe <- length(which(report$hwe_p > 0.01)) # 26 pass hwe

  all_pass <- report %>%
    filter(maf > 0.05, missingness < 0.01, hwe_p > 0.01)
  maf_hwe <- report %>%
    filter(maf > 0.05, hwe_p > 0.01, missingness >= 0.01)
  maf_only <- report %>%
    filter(maf > 0.05, hwe_p <= 0.01, missingness >= 0.01)
  miss_maf <- report %>%
    filter(maf > 0.05, hwe_p <= 0.01, missingness < 0.01)

  plt <- autoplot(report)

  expect_equal(
    plt[[1]]$data,
    tibble(
      Missing = c(TRUE, FALSE, FALSE, TRUE),
      MAF = c(TRUE, TRUE, TRUE, TRUE), # all pass MAF
      HWE = c(TRUE, TRUE, FALSE, FALSE),
      count = c(
        nrow(all_pass), nrow(maf_hwe),
        nrow(maf_only), nrow(miss_maf)
      ),
      intersection_id = c(1, 2, 3, 4)
    )
  )
})
