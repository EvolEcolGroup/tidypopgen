test_that("filter_high_relatedness removes necessary individuals", {
  # Create gentibble
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

  # Create an example relatedness matrix
  families_fbm <- tidypopgen:::.gt_get_fbm(families)
  king_matrix <- snp_king(families_fbm)

  # Find which individuals are over an arbitrary relateness threshold
  res <- filter_high_relatedness(
    king_matrix,
    kings_threshold = 0.2,
    verbose = FALSE
  )

  # individuals 11 and 12 at 0.2294
  # individuals 9 and 10 at 0.2742

  families$population <- c(rep("pop1", 12))
  families <- families %>% group_by(population)
  qc_report_fam <- qc_report_indiv(families, kings_threshold = 0.2)

  # Subset the matrix using res[[1]] - the individuals to keep
  sub_matrix <- king_matrix[
    c(as.numeric(res[[1]])),
    c(as.numeric(res[[1]]))
  ]

  # If we rerun, all individuals should be retained
  # (all TRUE in res2 output [[3]])
  res2 <- filter_high_relatedness(
    sub_matrix,
    kings_threshold = 0.2,
    verbose = FALSE
  )
  expect_true(all(res2[[3]]))
})
