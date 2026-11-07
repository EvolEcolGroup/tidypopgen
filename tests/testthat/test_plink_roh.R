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
    # genotypes, with 3 decimals (which is enough to recover the number of
    # heterozygous genotypes of runs of fewer than 1000 loci)
    expect_equal(roh$length_bp / 1000, hom$KB, tolerance = 1e-6,
      label = label
    )
    expect_identical(roh$n_het, as.integer(round(hom$PHET * hom$NSNP)),
      label = label
    )
  }
})

test_that("gt_roh_loci counts the individuals in runs as PLINK", {
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
  summaries <- list.files(plink_roh_dir, pattern = "\\.hom\\.summary$")
  expect_true(length(summaries) > 0)
  for (summary_file in summaries) {
    i <- match(sub("\\.hom\\.summary$", "", summary_file), cases$case)
    test_gt <- gen_tibble(datasets[[cases$dataset[i]]],
      backingfile = tempfile(), quiet = TRUE
    )
    args <- as.list(cases[i, arg_names])
    args <- args[!is.na(args)]
    # grouped, but counted over all individuals
    roh <- do.call(
      windows_indiv_roh,
      c(list(test_gt %>% group_by(population)), args)
    )
    roh_loci <- gt_roh_loci(roh, show_loci(test_gt), by_group = FALSE)
    plink_summary <- utils::read.table(
      file.path(plink_roh_dir, summary_file),
      header = TRUE, stringsAsFactors = FALSE
    )
    expect_identical(roh_loci$name, plink_summary$SNP)
    expect_identical(roh_loci$position, plink_summary$BP)
    expect_identical(roh_loci$n_roh,
      as.integer(plink_summary$AFF + plink_summary$UNAFF),
      label = summary_file
    )
  }
})
