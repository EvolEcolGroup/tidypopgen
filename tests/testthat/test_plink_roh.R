# Compare the runs of homozygosity of windows_indiv_roh() with PLINK 1.9
# --homozyg. The PLINK outputs are created by
# data-raw/validate/create_plink_roh_test_data.R (which also lists the
# cases). The tests are skipped if the data have not been created.

plink_roh_dir <- test_path("testdata", "plink_roh")

test_that("windows_indiv_roh gives the runs of PLINK --homozyg", {
  skip_if_not(
    file.exists(file.path(plink_roh_dir, "plink_version.txt")),
    paste(
      "PLINK test data not created",
      "(see data-raw/validate/create_plink_roh_test_data.R)"
    )
  )
  cases <- utils::read.csv(file.path(plink_roh_dir, "cases.csv"))
  datasets <- list(
    roh_sim = system.file("extdata", "roh", "roh_sim.bed",
      package = "tidypopgen"
    ),
    small = file.path(plink_roh_dir, "small.bed")
  )
  arg_names <- setdiff(names(cases), c("case", "dataset"))
  for (i in seq_len(nrow(cases))) {
    test_gt <- gen_tibble(datasets[[cases$dataset[i]]],
      backingfile = tempfile(), quiet = TRUE
    )
    args <- as.list(cases[i, arg_names])
    args <- args[!is.na(args)]
    roh <- do.call(windows_indiv_roh, c(list(test_gt), args))
    hom <- utils::read.table(
      file.path(plink_roh_dir, paste0(cases$case[i], ".hom")),
      header = TRUE, stringsAsFactors = FALSE
    )
    label <- cases$case[i]
    expect_identical(roh$id, as.character(hom$IID), label = label)
    expect_identical(
      as.character(roh$chromosome), as.character(hom$CHR),
      label = label
    )
    expect_identical(roh$start, hom$POS1, label = label)
    expect_identical(roh$end, hom$POS2, label = label)
    expect_identical(roh$n_loci, hom$NSNP, label = label)
    # PLINK prints the length in kb, and the proportion of heterozygous
    # genotypes, with 3 decimals
    expect_equal(roh$length_bp / 1000, hom$KB, tolerance = 1e-6,
      label = label
    )
    expect_equal(round(roh$n_het / roh$n_loci, 3), hom$PHET,
      tolerance = 1e-6, label = label
    )
  }
})
